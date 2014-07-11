# -*- coding: utf-8 -*-
"""
Created on 19 Dec 2013

@author: Kimon Tsitsikas

Copyright © 2012-2013 Kimon Tsitsikas, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the
terms  of the GNU General Public License version 2 as published by the Free
Software  Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR  PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Odemis. If not, see http://www.gnu.org/licenses/.
"""

from __future__ import division

from concurrent.futures._base import CancelledError, CANCELLED, FINISHED
import logging
import math
import numpy
from odemis import model
from odemis.dataio import hdf5
import heapq
import os
import time
from . import coordinates, transform


MAX_TRIALS_NUMBER = 2  # Maximum number of scan grid repetitions

def _DoFindOverlay(future, repetitions, dwell_time, max_allowed_diff, escan,
                   ccd, detector):
    """
    Scans a spots grid using the e-beam and captures the CCD image, isolates the 
    spots in the CCD image and finds the coordinates of their centers, matches the 
    coordinates of the spots in the CCD image to those of SEM image and calculates 
    the transformation values from optical to electron image (i.e. ScanGrid->
    DivideInNeighborhoods->FindCenterCoordinates-> ReconstructCoordinates->MatchCoordinates->
    CalculateTransform). In case matching the coordinates is infeasible, it automatically 
    repeats grid scan -and thus all steps until matching- with different parameters.
    future (model.ProgressiveFuture): Progressive future provided by the wrapper
    repetitions (tuple of ints): The number of CL spots are used
    dwell_time (float): Time to scan each spot (in s)
    max_allowed_diff (float): Maximum allowed difference in electron coordinates #m
    escan (model.Emitter): The e-beam scanner
    ccd (model.DigitalCamera): The CCD
    detector (model.Detector): The electron detector
    returns translation (Tuple of 2 floats), 
            scaling (Float), 
            rotation (Float): Transformation parameters
            transform_data : Transformation metadata
    raises:    
            CancelledError if cancelled
            ValueError if procedure failed 
    """
    logging.debug("Starting Overlay...")

    try:
        # Repeat until we can find overlay (matching coordinates is feasible)
        for trial in range(MAX_TRIALS_NUMBER):
            # Grid scan
            if future._find_overlay_state == CANCELLED:
                raise CancelledError()

            # Update progress of the future (it may be the second trial)
            future.set_end_time(time.time() +
                                estimateOverlayTime(future._scanner.dwell_time,
                                                    repetitions))

            # Wait for ScanGrid to finish
            optical_image, electron_coordinates, electron_scale = future._scanner.DoAcquisition()
            if future._find_overlay_state == CANCELLED:
                raise CancelledError()

            # Update remaining time to 6secs (hardcoded estimation)
            future.set_end_time(time.time() + 6)

            # Check if ScanGrid gave one image or list of images
            # If it is a list, follow the "one image per spot" procedure
            logging.debug("Isolating spots...")
            if isinstance(optical_image, list):
                opt_img_shape = optical_image[0].shape
                subimages = []
                subimage_coordinates = []
                for img in optical_image:
                    subspots, subspot_coordinates = coordinates.DivideInNeighborhoods(img, (1, 1), img.shape[0] / 2)
                    subimages.append(subspots[0])
                    subimage_coordinates.append(subspot_coordinates[0])
            else:
                # Distance between spots in the optical image (in optical pixels)
                optical_dist = escan.pixelSize.value[0] * electron_scale[0] / optical_image.metadata[model.MD_PIXEL_SIZE][0]
                opt_img_shape = optical_image.shape

                # Isolate spots
                if future._find_overlay_state == CANCELLED:
                    raise CancelledError()
                subimages, subimage_coordinates = coordinates.DivideInNeighborhoods(optical_image, repetitions, optical_dist)

            if not subimages:
                raise ValueError("Overlay failure")

            # Find the centers of the spots
            if future._find_overlay_state == CANCELLED:
                raise CancelledError()
            logging.debug("Finding spot centers with %d subimages...", len(subimages))
            spot_coordinates = coordinates.FindCenterCoordinates(subimages)

            # Reconstruct the optical coordinates
            if future._find_overlay_state == CANCELLED:
                raise CancelledError()
            optical_coordinates = coordinates.ReconstructCoordinates(subimage_coordinates, spot_coordinates)

            # Check if SEM calibration is correct. If this is not the case
            # generate a warning message and provide the ratio of X/Y scale.
            ratio = _computeGridRatio(optical_coordinates, repetitions)
            if not (0.9 < ratio < 1.1):
                logging.warning("SEM may needs calibration. X/Y ratio is %f.", ratio)
            else:
                logging.info("SEM X/Y ratio is %f.", ratio)

            opt_offset = (opt_img_shape[1] / 2, opt_img_shape[0] / 2)

            optical_coordinates = [(x - opt_offset[0], y - opt_offset[1]) for x, y in optical_coordinates]

            # TODO: Make function for scale calculation
            # Estimate the scale by measuring the distance between the first
            # two spots in optical and electron coordinates. 
            #  * For electrons, it's easy as we've placed them.
            #  * For optical, we pick one spot, and measure the distance to the
            #    closest spot. Should be more precise than previous estimation.
            p1 = optical_coordinates[0]
            def dist_to_p1(p):
                diff = [a - b for a, b in zip(p1, p)]
                return math.hypot(*diff)
            optical_dist = min(dist_to_p1(p) for p in optical_coordinates[1:])
            scale = electron_scale[0] / optical_dist

            # max_allowed_diff in pixels
            max_allowed_diff_px = max_allowed_diff / escan.pixelSize.value[0]

            # Match the electron to optical coordinates
            if future._find_overlay_state == CANCELLED:
                raise CancelledError()

            logging.debug("Matching coordinates...")
            known_ec, known_oc = coordinates.MatchCoordinates(optical_coordinates,
                                                         electron_coordinates,
                                                         scale,
                                                         max_allowed_diff_px)
            if known_ec:
                break
            else:
                if trial < MAX_TRIALS_NUMBER - 1:
                    future._scanner.dwell_time = future._scanner.dwell_time * 1.2 + 0.1
                    logging.warning("Trying with dwell time = %g s...", future._scanner.dwell_time)
        else:
            # Make failure report
            _MakeReport(optical_image, repetitions, escan.magnification.value, escan.pixelSize.value, dwell_time, electron_coordinates)
            raise ValueError("Overlay failure")


        # Calculate transformation parameters
        if future._find_overlay_state == CANCELLED:
            raise CancelledError()

        # We are almost done... about 1 s left
        future.set_end_time(time.time() + 1)

        logging.debug("Calculating transformation...")
        ret = transform.CalculateTransform(known_ec, known_oc)

        if future._find_overlay_state == CANCELLED:
            raise CancelledError()
        logging.debug("Calculating transform metadata...")

        transform_data = _transformMetadata(optical_image, ret, escan, ccd)
        if transform_data == []:
            raise ValueError('Metadata is missing')

        # Also indicate which dwell time eventually worked
        transform_data[model.MD_DWELL_TIME] = dwell_time

        logging.debug("Overlay done.")
        return ret, transform_data
    except Exception as exp:
        logging.debug("Finding overlay failed", exc_info=1)
        raise exp
    finally:
        with future._overlay_lock:
            if future._find_overlay_state == CANCELLED:
                raise CancelledError()
            future._find_overlay_state = FINISHED

def _CancelFindOverlay(future):
    """
    Canceller of _DoFindOverlay task.
    """
    logging.debug("Cancelling overlay...")

    with future._overlay_lock:
        if future._find_overlay_state == FINISHED:
            return False
        future._find_overlay_state = CANCELLED
        future._scanner.CancelAcquisition()
        logging.debug("Overlay cancelled.")

    return True

def _computeGridRatio(coord, shape):
    """
    coord (list of tuple of 2 floats): coordinates
    shape (2 ints): X and Y number of coordinates 
    return (float): ratio X/Y
    """
    x_cors = [i[0] for i in coord]
    y_cors = [i[1] for i in coord]
    x_max_cors = numpy.mean(heapq.nlargest(shape[0], x_cors))
    x_min_cors = numpy.mean(heapq.nsmallest(shape[0], x_cors))
    y_max_cors = numpy.mean(heapq.nlargest(shape[1], y_cors))
    y_min_cors = numpy.mean(heapq.nsmallest(shape[1], y_cors))
    x_scale = x_max_cors - x_min_cors
    y_scale = y_max_cors - y_min_cors
    return x_scale / y_scale

def estimateOverlayTime(dwell_time, repetitions):
    """
    Estimates overlay procedure duration
    """
    return 6 + dwell_time * numpy.prod(repetitions)  # s

def _transformMetadata(optical_image, transformation_values, escan, ccd):
    """
    Returns the transform metadata for the optical image based on the 
    transformation values
    """
    escan_pxs = escan.pixelSize.value
    logging.debug("Ebeam pixel size: %g ", escan_pxs[0])
    ((calc_translation_x, calc_translation_y),
             (calc_scaling_x, calc_scaling_y),
                                calc_rotation) = transformation_values

    # Update scaling
    scale = (escan_pxs[0] * calc_scaling_x,
             escan_pxs[1] * calc_scaling_y)

    transform_md = {model.MD_ROTATION_COR: calc_rotation}

    # X axis is same direction in image and physical referentials
    # Y axis is opposite direction, that's why we don't need a "-"
    position_cor = (-scale[0] * calc_translation_x,
                    scale[1] * calc_translation_y)
    logging.debug("Center shift correction: %s", position_cor)
    transform_md[model.MD_POS_COR] = position_cor
    if isinstance(optical_image, list):
        opt_img_pxs = optical_image[0]
    else:
        opt_img_pxs = optical_image
    try:
        pixel_size = opt_img_pxs.metadata[model.MD_PIXEL_SIZE]
    except KeyError:
        logging.warning("No MD_PIXEL_SIZE data available")
        return transform_md
    pixel_size_cor = (scale[0] / pixel_size[0],
                      scale[1] / pixel_size[1])
    logging.debug("Pixel size correction: %s", pixel_size_cor)
    transform_md[model.MD_PIXEL_SIZE_COR] = pixel_size_cor

    return transform_md


def _MakeReport(optical_image, repetitions, magnification, pixel_size, dwell_time, electron_coordinates):
    """
    Creates failure report in case we cannot match the coordinates.
    optical_image (2d array): Image from CCD
    repetitions (tuple of ints): The number of CL spots are used
    dwell_time (float): Time to scan each spot (in s)
    electron_coordinates (list of tuples): Coordinates of e-beam grid
    """
    path = os.path.join(os.path.expanduser(u"~"), u"odemis-overlay-report",
                        time.strftime(u"%Y%m%d-%H%M%S"))
    os.makedirs(path)
    hdf5.export(os.path.join(path, u"OpticalGrid.h5"), optical_image)
    report = open(os.path.join(path, u"report.txt"), 'w')
    report.write("\n****Overlay Failure Report****\n\n"
                 + "\nSEM magnification:\n" + str(magnification)
                 + "\nSEM pixel size:\n" + str(pixel_size)
                 + "\nGrid size:\n" + str(repetitions)
                 + "\n\nMaximum dwell time used:\n" + str(dwell_time)
                 + "\n\nElectron coordinates of the scanned grid:\n" + str(electron_coordinates)
                 + "\n\nThe optical image of the grid can be seen in OpticalGrid.h5\n\n")
    report.close()

    logging.warning("Failed to find overlay. Please check the failure report in %s.",
                    path)


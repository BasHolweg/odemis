# !/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created on 10 February 2021

@author: Sabrina Rossberger


Copyright © 2021 Sabrina Rossberger, Delmic



"""
from __future__ import division

import logging
import sys
import time
import numpy

from odemis import model
from odemis.acq.align.spot import FindGridSpots
from odemis.util.driver import get_backend_status, BACKEND_RUNNING

std_dark_gain = 'True'
dwell_time = 10e-6

def acquire_MegaField(mppc, dataflow, descanner, multibeam, stage, ccd, beamshift, mm, field_images, dwell_time):
    """
    TODO
    """
    # the location of the MPPC array mapped to the diagnostic camera TODO I don't understand this comment
    mean_spot = (756, 560)  # in pixels; (65.2, 92.6)um

    # setting of the scanner
    multibeam.scanOffset.value = (-0.0935 / 1, 0.0935 / 1)
    multibeam.scanGain.value = (0.0935 / 1, -0.0935 / 1)

    # setting of the descanner
    offset_x = 0.1624
    offset_y = 0.216955
    descanner.scanOffset.value = (offset_x, offset_y)
    descanner.scanGain.value = (offset_x + 0.0082, offset_y - 0.0082)

    # routine to align the multiprobe pattern with the MPPC using the mapping of the MPPC to diagnosic camera
    # this part uploads the standard descan offset
    beamshift.shift.value = (0, 0)
    mppc.dataContent.value = "empty"  # don't return an image
    mppc.filename.value = time.strftime("testing_megafield_id-%Y-%m-%d-%H-%M-%S")
    multibeam.dwellTime.value = 1e-6
    dataflow.subscribe(image_received)
    dataflow.next((1001, 1001))
    time.sleep(dwell_time * mppc.cellCompleteResolution.value[0] * mppc.cellCompleteResolution.value[1] + 0.5)
    dataflow.unsubscribe(image_received)
    ccd_image = ccd.data.get(asap=False)
    spot_coordinates, *_ = FindGridSpots(ccd_image, (8, 8))
    spot_coordinates[:, 1] = ccd_image.shape[1] - spot_coordinates[:, 1]
    shift_descan = numpy.mean(spot_coordinates, axis=0) - mean_spot
    offset_x = offset_x + shift_descan[0] * 0.0024642857142857144  # 1/140 * (3.45/10)
    offset_y = offset_y + shift_descan[1] * 0.0024642857142857144
    descanner.scanOffset.value = (offset_x, offset_y)
    descanner.scanGain.value = (offset_x + 0.0082 / 1, offset_y - 0.0082 / 1)

    # final settings for the acquisition
    multibeam.dwellTime.value = dwell_time
    mppc.overVoltage.value = 1.8
    multibeam.scanDelay.value = (0e-5, 0.0)
    mppc.acqDelay.value = 0.00  # acqDelay >= scanner.scanDelay
    descanner.physicalFlybackTime = 250e-4  # hardcoded+
    mppc.cellCompleteResolution.value = (900, 900)
    multibeam.resolution.value = (800 * 8, 800 * 8)  # change to (900*8, 900*8) for FIELD CORRECTIONS
    mppc.cellTranslation.value = tuple(tuple((0, 0) for i in range(0, mppc.shape[0])) for i in range(0, mppc.shape[1]))

    if std_dark_gain == 'True':  # TODO should this option stay?
        mppc.cellDarkOffset.value = ((32610, 32512, 32631, 32705, 32833, 32587, 32129, 32924),
                                     (32776, 32781, 32609, 32733, 32977, 32764, 32984, 32790),
                                     (32479, 32402, 32931, 32902, 32735, 33041, 32780, 32821),
                                     (33014, 32957, 32472, 32462, 32234, 32571, 32560, 32737),
                                     (32870, 32688, 32587, 32591, 32896, 32751, 32749, 32876),
                                     (32499, 32963, 32297, 32608, 32680, 32585, 32799, 32562),
                                     (32359, 32625, 32648, 32628, 33083, 32741, 32603, 32763),
                                     (32688, 32810, 32393, 32497, 32740, 32280, 32633, 32578))

        mppc.cellDigitalGain.value = ((1.4033838418883835, 1.3815675127292752, 1.408208449247768, 1.1587998586046597,
                                       1.149068548526435, 1.2830011181580079, 1.197738151157632, 1.357734293731753),
                                      (1.2338445250472097, 1.3005552730745185, 1.342377686636198, 1.1638021287242348,
                                       1.1847154644164415, 1.2459790451569392, 1.1276893842388696, 1.3360160349190446),
                                      (1.2003976000630705, 1.1364007415897017, 1.2354973047369953, 1.2014533950701192,
                                       1.2018868539752536, 1.2038216228144334, 1.0947659586205654, 1.3940386115155805),
                                      (1.241508610162757, 1.0892308414346563, 1.0690610417160937, 1.0737925448880836,
                                       1.1035458304227468, 1.2149661227387953, 1.0583906475139082, 1.3283332490759001),
                                      (1.270294902441606, 1.1642958064806992, 1.056405338256749, 1.01059850445321,
                                       1.0667511185564946, 1.2427729708033957, 1.0738405392021235, 1.1870507455532333),
                                      (1.201524752392076, 1.1375566758141369, 1.0834190313448255, 1.0182433195297147,
                                       1.0183692719045843, 1.1671973417743966, 1.0820640322253232, 1.2880083559414344),
                                      (1.2117217304291876, 1.0734889479555605, 1.161479139763327, 1.1959017403008299,
                                       1.0018479817582293, 1.1012583187834029, 1.1227877982571943, 1.3500507703903017),
                                      (1.376763516610481, 1.2815030552541757, 1.2483421772462209, 1.461928114869715,
                                       1.1607017269348425, 1.1559352001248653, 1.2795062327345301, 1.4106657020179465))

    # current overscan parameters
    mppc.cellTranslation.value = (((0, 0), (4, 2), (13, 11), (17, 16), (21, 23), (26, 28), (33, 33), (40, 39)),
                                  ((7, 0), (11, 5), (19, 12), (27, 21), (30, 26), (34, 31), (42, 36), (50, 44)),
                                  ((14, 4), (19, 10), (26, 16), (36, 21), (36, 26), (43, 35), (52, 40), (60, 47)),
                                  ((17, 6), (20, 12), (32, 16), (39, 26), (44, 31), (47, 36), (61, 45), (70, 52)),
                                  ((28, 7), (30, 15), (37, 20), (54, 26), (54, 31), (57, 37), (64, 42), (70, 41)),
                                  ((40, 5), (41, 13), (46, 19), (54, 24), (64, 30), (67, 32), (72, 36), (79, 41)),
                                  ((50, 1), (50, 10), (54, 17), (61, 23), (76, 27), (76, 33), (78, 36), (86, 38)),
                                  ((65, 1), (62, 13), (63, 20), (69, 23), (80, 30), (83, 33), (87, 35), (92, 38)))

    # mppc.cellDarkOffset.value = tuple(tuple(0 for i in range(0, mppc.shape[0])) for i in range(0, mppc.shape[1]))
    # mppc.cellDigitalGain.value = tuple(tuple(1 for i in range(0, mppc.shape[0])) for i in range(0, mppc.shape[1]))

    image_counter = 0  # used to acquire the first field twice due to intensity issues (to be removed)
    # Note: numbering of offloaded images is incorrect then: image 0 equal image 1

    mppc.dataContent.value = "empty"
    mppc.filename.value = time.strftime("testing_megafield_id-%Y-%m-%d-%H-%M-%S")
    dataflow = mppc.data
    dataflow.subscribe(image_received)

    pos = stage.position.value

    overlap = 0.875  # for debugging acquire with overlap of field images
    move_col_x = -(3.2 * 8 * overlap) * 1e-6  # or: shift = 3200 * 4 * 1e-9   # half a field 3200 px, 4nm pixel size

    # center of multiprobe pattern on diagnostic camera TODO code moved?

    # stage position correction
    # read out position stage and position metrology module
    stage_pos_start = (stage.position.value['x'], stage.position.value['y'])
    # mm_pos_start = (mm.position.value['x'], mm.position.value['y'])
    # map MM position to stage position: calculate the offset between the 2 coordinate systems
    # translation_stage_to_mm = numpy.array(stage_pos_start) - numpy.array(mm_pos_start)

    stage_pos_x = stage_pos_start[0]
    stage_pos_y = stage_pos_start[1]

    for row in range(field_images[1]):  # y axis, move stage on y,

        for col in range(field_images[0]):  # x axis, move stage on x

            logging.debug("col (x): {} row (y): {} stage: {}".format(col, row, stage.position.value))

            if image_counter < 1:
                dataflow.next((1000, 1000))  # acquire the first field twice
                time.sleep(
                    dwell_time * mppc.cellCompleteResolution.value[0] * mppc.cellCompleteResolution.value[1] + 0.5)
                time.sleep(0.5)

            ####################################################################################
            # Correct for beam shift.
            ccd_image = ccd.data.get(asap=False)
            spot_coordinates, *_ = FindGridSpots(ccd_image, (8, 8))
            logging.debug("spot_coordinates position: {}".format(numpy.mean(spot_coordinates, axis=0)))
            # Transfer to a coordinate system with the origin in the bottom left.
            spot_coordinates[:, 1] = ccd_image.shape[1] - spot_coordinates[:, 1]
            shift = numpy.mean(spot_coordinates, axis=0) - mean_spot
            # print("mean position: {}".format(numpy.mean(spot_coordinates, axis=0)))
            # print("Shift x: {}, shift y: {}".format(-shift[0] / scan_gain[0], -shift[1] / scan_gain[1]))
            # Compensate for shift with dc beam shift
            logging.debug("beam shift position: {}".format(beamshift.shift.value))
            # convert shift from pixels to um
            shift_um = shift * 3.45e-6 / 40  # pixelsize ueye cam and 40x magnification
            logging.debug("beam shift um: {}".format(shift_um))
            cur_beam_shift_pos = numpy.array(beamshift.shift.value)
            beamshift.shift.value = tuple(cur_beam_shift_pos + shift_um)
            #############################\#################################################

            # acquire field
            dataflow.next((col, row))  # acquire the field image (y, x)
            # #TODO investigate order!! with testcases in wrapper order not checked with filename pattern
            image_counter += 1
            time.sleep(dwell_time * mppc.cellCompleteResolution.value[0] * mppc.cellCompleteResolution.value[1] + 0.5)
            # wait a bit so the image can be acquired TODO implement proper future for image acquisition (will be done in acq manager)

            # TODO only execute when more than 1 col
            ##old: stage move without mm correction ####################################################################
            # Move stage one field in the positive x direction of the scanner.
            # Sample should move on TSF screen to the right if we apply a positive x stage move.
            # (0,0) is top left corner -> move sample to the left, which is negative x stage move.
            # stage.moveRel({"x": move_col_x}).result()  # in meter!

            ##stage move with mm correction ############################################################################
            # # TODO move into separate funciton call
            # # stage pos correction due to inaccuracy in stage pos
            # # move the stage by absolute position based on first field image position -> calculate next stage position
            # # (to avoid accumulation of relative stage movement errors)
            # # Move stage one field in the positive x direction of the scanner.
            # # Sample should move on TSF screen to the right if we apply a positive x stage move.
            # # (0,0) is top left corner -> move sample to the left, which is negative x stage move.
            stage_pos_x = stage_pos_start[0] + move_col_x * (col + 1)
            stage.moveAbs({"x": stage_pos_x}).result()  # in meter!

            # read new stage position
            stage_pos_target = (stage_pos_x, stage_pos_y)
            # read new MM position
            # stage_pos_actual = (mm.position.value['x'], mm.position.value['y'])
            # calc difference from theoretical stage pos and MM (actual stage) pos taking the offset into account
            # pos_corr_x = stage_pos_target[0] - stage_pos_actual[0] - translation_stage_to_mm[
            #     0]  # 2 - 11.5 - (-10) = 0.5
            # pos_corr_y = stage_pos_target[1] - stage_pos_actual[1] - translation_stage_to_mm[1]
            # correct with beam shift (first read, then set as also absolute pos)
            # beamshift_pos_curr = beamshift.shift.value
            # beamshift.shift.value = (beamshift_pos_curr[0] + pos_corr_x, beamshift_pos_curr[1] + pos_corr_y)
            # print("beam shift after column move um: {}".format(beamshift.shift.value))
            ##############################################################################
            time.sleep(0.5)
        # Move stage back to starting position to start a new row
        stage.moveAbs(pos).result()  # TODO integrate into absolute movement below
        # move in y
        move_row_y = (row + 1) * (3.2 * 8 * overlap) * 1e-6
        logging.debug(move_row_y)

        ##old: stage move without mm correction ####################################################################
        # stage.moveRel({'y': move_row_y}).result()   # in meter!

        ##stage move with mm correction ############################################################################
        stage_pos_y = stage_pos_start[1] + move_row_y * (row + 1)
        stage_pos_x = stage_pos_start[0]
        stage.moveAbs({'y': stage_pos_y}).result()  # in meter!
        stage.moveAbs({'x': stage_pos_x}).result()
        # read new stage position
        stage_pos_target = (stage_pos_x, stage_pos_y)
        # read new MM position
        # stage_pos_actual = (mm.position.value['x'], mm.position.value['y'])
        # calc difference from theoretical stage pos and MM (actual stage) pos taking the offset into account
        # pos_corr_x = stage_pos_target[0] - stage_pos_actual[0] - translation_stage_to_mm[0]  # 2 - 11.5 - (-10) = 0.5
        # pos_corr_y = stage_pos_target[1] - stage_pos_actual[1] - translation_stage_to_mm[1]
        # correct with beam shift (first read, then set as also absolute pos)
        logging.debug("moving row")
        # beamshift_pos_curr = beamshift.shift.value
        # beamshift.shift.value = (beamshift_pos_curr[0] + pos_corr_x, beamshift_pos_curr[1] + pos_corr_y)
        # print("beam shift row move um: {}".format(beamshift.shift.value))
        ##############################################################################
        time.sleep(0.5)

    time.sleep(dwell_time * mppc.cellCompleteResolution.value[0] * mppc.cellCompleteResolution.value[
        1] + 1)  # Allow 1.5 seconds per field image to offload.
    dataflow.unsubscribe(image_received)

    # Move stage back to starting position of the stage
    stage.moveAbs(pos).result()


##############################################################################
# stage position correction

# read out position stage and position metrology module
# map MM position to stage position: calculate the offset between the 2 coordinate systems
# move the stage by absolute position based on first field image position -> calulate next stage position
# (to avoid accumulation of relative stage movement errors)
# read new stage position
# read new MM position
# calc difference from theoretical stage pos and MM (actual stage) pos taking the offset into account
# correct with beam shift (first read, then set as also absolute pos)
##############################################################################


def image_received(self, *args):
    """
    Subscriber for test cases which counts the number of times it is notified.
    *args contains the image/data which is received from the subscriber.
    """
    print("image received")


def main(args):
    logging.getLogger().setLevel(logging.DEBUG)

    if get_backend_status() != BACKEND_RUNNING:
        raise ValueError("Backend is not running.")

    ccd = model.getComponent(role="diagnostic-ccd")
    mppc = model.getComponent(role="mppc")
    dataflow = mppc.data
    MirrorDescanner = model.getComponent(role="descanner")
    MultiBeamScanner = model.getComponent(role="multibeam")
    scanner = model.getComponent(role="mb-scanner")
    stage = model.getComponent(role="stage-scan")
    mm = model.getComponent(role="stage-pos")
    beamshift = model.getComponent(role="ebeam-shift")

    # Set initial hardware settings
    scanner.multiBeamMode.value = True  # True selects multibeam mode
    scanner.external.value = True  # Put microscope in external mode
    # scanner.power.value = True  # Switch on beam
    scanner.blanker.value = False  # unblank the beam
    scanner.horizontalFoV.value = 2.2e-05  # corresponds to mag = 5000x in quad view
    # Provide coordinate transform for beam shift.
    beamshift.updateMetadata({model.MD_CALIB: scanner.beamShiftTransformationMatrix.value})

    time.sleep(2)  # wait a bit TODO why?

    # size of mega field
    field_images = (1, 1)  # (x, y) Note: x (horizontal) = col, y (vertical) = row
    dwell_time = 10e-6
    # debugging: (1,3): expect 3 images in y direction (rows), files 0_0_, 1_0_, 2_0_
    # debugging: (2,3): expect 3 images in y direction (rows), 2 images in x dir (columns) files 0_0_, 1_0_, 2_0_, 0_1_, 1_1_, 2_1_

    try:
        acquire_MegaField(mppc, dataflow, MirrorDescanner, MultiBeamScanner, stage, ccd, beamshift, mm, field_images,
                          dwell_time)
    except Exception as exp:
        logging.error("%s", exp, exc_info=True)
    finally:
        dataflow.unsubscribe(image_received)
        scanner.blanker.value = True  # blank the beam
        # scanner.power.value = False  # Switch off beam
    print("Done")


if __name__ == '__main__':
    ret = main(sys.argv)
    exit(ret)

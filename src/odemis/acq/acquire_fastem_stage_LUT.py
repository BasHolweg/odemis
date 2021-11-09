# !/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created on 10 February 2021

@author: Sabrina Rossberger, Thera Pals, Wilco Zuidema


Copyright © 2021 Sabrina Rossberger, Delmic



"""
from __future__ import division

import logging
import os
import sys
import threading
import time
from datetime import datetime

import numpy

from odemis import dataio
from odemis import model
from odemis.acq.align.spot import FindGridSpots
from odemis.util.driver import get_backend_status, BACKEND_RUNNING

std_dark_gain = False
MEAN_SPOT = (776, 568)  # in pixels on DC; (725 * 3.45, 417 * 3.45) um with 3.45um = pixelsize of DC


def mppc2mp(ccd, multibeam, descanner, mppc, dataflow):

    # routine to align the spot grid with the MPPC using the mapping of the MPPC to diagnostic camera

    # mppc.cellTranslation.value = tuple(tuple((0, 0) for i in range(0, mppc.shape[0])) for i in range(0, mppc.shape[1]))
    # mppc.cellDarkOffset.value = tuple(tuple(2**15 for i in range(0, mppc.shape[0])) for i in range(0, mppc.shape[1]))
    # mppc.cellDigitalGain.value = tuple(tuple(1 for i in range(0, mppc.shape[0])) for i in range(0, mppc.shape[1]))

    # setting of the scanner
    multibeam.scanOffset.value = (-0.0935 / 1, 0.0935 / 1)
    multibeam.scanGain.value = (0.0935 / 1, -0.0935 / 1)

    # setting of the descanner
    # good offset positions
    good_offset_x = 0.032043  # a.u.
    good_offset_y = 0.053406
    # bad values for debug:
    # offset_x = 0.02043
    # offset_y = 0.033406

    offset_x = 0.0213
    offset_y = -0.0671
    print("inital descan offset x: {}; inital descan offset y: {}".format(offset_x, offset_y))

    descanner.scanOffset.value = (offset_x, offset_y)
    descanner.scanGain.value = (offset_x + 0.0082, offset_y - 0.0082)

    mppc.filename.value = time.strftime("testing_megafield_id-%Y-%m-%d-%H-%M-%S")
    multibeam.dwellTime.value = 2e-6

    # create folder to store calibration images
    dir_name = "mppc-to-mp-translation_acq"

    user_path = os.path.expanduser("~")
    path_images = os.path.join(user_path, "development", "fastem-calibrations", "images")
    os.makedirs(path_images, exist_ok=True)  # check if directory already existing, if not create it

    # create directory with name as specified in parameter dir_name
    path_calib = os.path.join(path_images, dir_name)
    os.makedirs(path_calib, exist_ok=True)  # check if directory already existing, if not create it

    # create directory for saving the calibration images with date and timestamp as name
    path = os.path.join(path_calib, datetime.now().strftime("%Y%m%d-%H%M%S"))
    os.makedirs(path)

    # upload the settings to the HW via the ASM
    # TODO replace with get call; disadvantage it will be labelled as (0,0), which is not handy for debugging
    dataflow.subscribe(on_field_image)  # no event set yet; only when receiving data
    # clear event that is set in the callback
    image_received.clear()
    dataflow.next((1001, 1001))  # generate data,
    logging.debug("requested calibration image acquired via ASM")
    if not image_received.wait(
            multibeam.dwellTime.value
            * mppc.cellCompleteResolution.value[0]
            * mppc.cellCompleteResolution.value[1]
            + 1):
        # wait until the current status is returned; if status not True (no data received by callback function and
        # thus event not set to true, but False (timeout) -> raise error
        # wait until the current status is returned; if status not True (no data received by callback function and
        # thus event not set to true), but False (timeout) -> raise error
        dataflow.unsubscribe(on_field_image)
        raise TimeoutError("Calibration timed out")
    # blocking; callback function receive data; event is set to True; wait until flag is true
    logging.debug("received calibration image from ASM")
    dataflow.unsubscribe(on_field_image)

    ccd_image = ccd.data.get(asap=False)  # asap=False: wait until new image is acquired (don't read from buffer)

    logging.debug("received diagnostic camera image")

    # save the calibration image
    file_name = "dc_image_before_alignment_acq.tiff"
    file_path = os.path.join(path, file_name)
    dataio.tiff.export(file_path, ccd_image)

    # calculate the difference in the multiprobe position on the diagnostic camera image compared to the "good"
    # (mean_spot) position, where the mp signal is mapped correctly onto the mppc detector
    # TODO: replace (8,8) with shape attribute
    spot_coordinates, *_ = FindGridSpots(ccd_image, (8, 8))
    spot_coordinates[:, 1] = ccd_image.shape[1] - spot_coordinates[:, 1]
    shift_descan = numpy.mean(spot_coordinates, axis=0) - MEAN_SPOT
    print("correction for descan offset: {}".format(shift_descan * (3.45 / 10)))

    offset_x = offset_x + shift_descan[0] * 0.000196022
    offset_y = offset_y + shift_descan[1] * 0.000196022
    descanner.scanOffset.value = (offset_x, offset_y)
    descanner.scanGain.value = (offset_x + 0.0082 / 1, offset_y - 0.0082 / 1)

    print("expected descan offset x: {}; calibrated descan offset y: {}".format(good_offset_x, good_offset_y))
    print("calibrated descan offset x: {}; calibrated descan offset y: {}".format(offset_x, offset_y))

    # upload the calculated offsets to the ASM
    dataflow.get(dataContent="empty")

    # save the calibration image after applying the new settings
    ccd_image = ccd.data.get(asap=False)  # asap=False: wait until new image is acquired (don't read from buffer)

    # save the calibration image
    file_name = "dc_image_after_alignment_acq.tiff"
    file_path = os.path.join(path, file_name)
    dataio.tiff.export(file_path, ccd_image)


def correct_stage_magnetic_field(ccd, beamshift):

    logging.debug("Perform paramagnetic field correction.")
    ccd_image = ccd.data.get(asap=False)  # asap=False: wait until new image is acquired (don't read from buffer)
    spot_coordinates, *_ = FindGridSpots(ccd_image, (8, 8))
    # Transfer to a coordinate system with the origin in the bottom left.
    spot_coordinates[:, 1] = ccd_image.shape[1] - spot_coordinates[:, 1]
    print("spot_coordinates position: {}".format(numpy.mean(spot_coordinates, axis=0)))
    shift = numpy.mean(spot_coordinates, axis=0) - MEAN_SPOT
    # Compensate for shift with dc beam shift
    # convert shift from pixels to um
    shift_um = shift * 3.45e-6 / 40  # pixelsize ueye cam and 40x magnification
    print("beam shift adjustment required due to stage magnetic field (um): {}".format(shift_um))
    cur_beam_shift_pos = numpy.array(beamshift.shift.value)
    print("current beam shift um: {}".format(beamshift.shift.value))
    beamshift.shift.value = (shift_um + cur_beam_shift_pos)
    print("new beam shift um: {}".format(beamshift.shift.value))


def correct_stage_pos(stage, mm, beamshift, row, col , x_stage_pos, y_stage_pos, x_stage_pos_rot, y_stage_pos_rot):

    # stage pos correction due to inaccuracy in stage pos

    ##############################################################################
    # plan for stage position correction

    # read out position stage and position metrology module
    # map MM position to stage position: calculate the offset between the 2 coordinate systems
    # move the stage by absolute position based on first field image position -> calulate next stage position
    # (to avoid accumulation of relative stage movement errors)
    # read new stage position
    # read new MM position
    # calc difference from theoretical stage pos and MM (actual stage) pos taking the offset into account
    # correct with beam shift (first read, then set as also absolute pos)
    ##############################################################################

    x_stage_pos_mm = mm.position.value['x'] + x_stage_pos_rot
    y_stage_pos_mm = mm.position.value['y'] + y_stage_pos_rot

    # read new MM position
    stage_pos_actual = (0, 0)

    # calc difference from theoretical stage pos and MM (actual stage) pos taking the offset into account
    pos_corr_x = stage_pos_actual[0] - x_stage_pos_mm[row, col]
    pos_corr_y = stage_pos_actual[1] - y_stage_pos_mm[row, col]
    tries = 0
    while numpy.abs(pos_corr_x) > 700e-9 or numpy.abs(pos_corr_y) > 700e-9:
        stage.moveAbs({"x": x_stage_pos[row, col]+10e-6}).result()  # in meter!
        time.sleep(0.4)
        stage.moveAbs({"y": y_stage_pos[row, col]+10e-6}).result()  # in meter!
        time.sleep(0.4)
        stage.moveAbs({"x": x_stage_pos[row, col]}).result()  # in meter!
        time.sleep(0.4)
        stage.moveAbs({"y": y_stage_pos[row, col]}).result()  # in meter!
        time.sleep(0.4)
        stage_pos_actual = (mm.position.value['x'], mm.position.value['y'])
        pos_corr_x = stage_pos_actual[0] - x_stage_pos_mm[row, col]
        pos_corr_y = stage_pos_actual[1] - y_stage_pos_mm[row, col]
        print('the position correction required is x:%e y:%e' % (pos_corr_x, pos_corr_y))
        tries = tries + 1
        if tries > 3:
            print("stage is not complying with our orders")
            break

    print('the position correction required is x:%e y:%e' % (pos_corr_x, pos_corr_y))
    # correct with beam shift (first read, then set as also absolute pos)
    beamshift_pos_curr = beamshift.shift.value
    if numpy.abs(beamshift_pos_curr[0] + pos_corr_x) < 0.041e-3 or \
            numpy.abs(beamshift_pos_curr[1] + pos_corr_y) < 0.041e-3:
        beamshift.shift.value = (beamshift_pos_curr[0] + pos_corr_x, beamshift_pos_curr[1] + pos_corr_y)
    else:
        print("The beam shift is out of range, the field cannot be aligned with beam shift")

    print("current beam shift um: {}".format(beamshift.shift.value))


def settings_megafield(multibeam, descanner, mppc, dwell_time):

    # adjust settings for megafield image acquisition

    multibeam.dwellTime.value = dwell_time
    mppc.overVoltage.value = 1.5
    multibeam.scanDelay.value = (0e-5, 0.0)
    mppc.acqDelay.value = 0.00  # acqDelay >= scanner.scanDelay
    descanner.physicalFlybackTime = 250e-4  # hardcoded
    mppc.cellCompleteResolution.value = (900, 900)
    multibeam.resolution.value = (800*8, 800*8)  # change to (900*8, 900*8) for FIELD CORRECTIONS

    if std_dark_gain:
        mppc.cellDarkOffset.value = ((32609, 32512, 32632, 32707, 32833, 32590, 32129, 32915),
                                     (32775, 32782, 32607, 32731, 32974, 32770, 32980, 32787),
                                     (32480, 32403, 32931, 32905, 32734, 33041, 32779, 32821),
                                     (33012, 32958, 32473, 32462, 32235, 32570, 32558, 32735),
                                     (32871, 32687, 32583, 32591, 32898, 32750, 32748, 32875),
                                     (32497, 32963, 32294, 32607, 32681, 32585, 32799, 32562),
                                     (32357, 32625, 32650, 32630, 33080, 32741, 32606, 32766),
                                     (32687, 32809, 32393, 32497, 32740, 32279, 32634, 32579))
        mppc.cellDigitalGain.value = ((1.0310130234815644, 0.98911785399073, 1.0181627423471624, 0.8577040961758531, 0.8880908034347579, 0.9890579976899265, 0.9383996622711986, 1.1113391502317609),
                                      (0.8593799077749432, 0.8986987222165352, 0.9363027653136717, 0.822184009271164, 0.8479924018374657, 0.9021588056776584, 0.837452960039753, 1.058894453119832),
                                      (0.8465728724395954, 0.7903281330555539, 0.8599569977952789, 0.8457909791302854, 0.8256918838727555, 0.8533798498272129, 0.8002763945317064, 1.0731022584055032),
                                      (0.8621492151552735, 0.7531764833868415, 0.7352312175442495, 0.7401340718521829, 0.7597036584238763, 0.844865700403647, 0.7770241838127158, 0.9983141992621981),
                                      (0.8665074854440534, 0.8098175493432067, 0.7311517895733345, 0.6948023507297094, 0.7318181050458212, 0.8681639239369103, 0.7755633808103727, 0.8923224747397783),
                                      (0.8335233462552188, 0.7881700089377985, 0.749662487713841, 0.7158504604573277, 0.701682390230767, 0.7977500467587251, 0.7914494398257048, 0.9763284195350683),
                                      (0.8521377276528598, 0.7643034555179854, 0.7983025624619302, 0.8379968529320432, 0.7114937102240616, 0.7570616654458873, 0.8074662374392703, 1.0288508064914936),
                                      (0.9513008682286824, 0.8964435843490056, 0.8529914655577808, 0.995611934322976, 0.8125292911222118, 0.8228977944874873, 0.9327903752318961, 1.040315455335232))

    # TODO values should be set in calibration (field corrections)
    # debug
    # mppc.cellDarkOffset.value = tuple(tuple(0 for i in range(0, mppc.shape[0])) for i in range(0, mppc.shape[1]))
    # mppc.cellDigitalGain.value = tuple(tuple(1 for i in range(0, mppc.shape[0])) for i in range(0, mppc.shape[1]))

    # debug
    # mppc.cellTranslation.value = tuple(tuple((0, 0) for i in range(0, mppc.shape[0])) for i in range(0, mppc.shape[1]))

    # current overscan parameters
    # mppc.cellTranslation.value = (((0, 0), (4, 2), (13, 11), (17, 16), (21, 23), (26, 28), (33, 33), (40, 39)),
    #                               ((7, 0), (11, 5), (19, 12), (27, 21), (30, 26), (34, 31), (42, 36), (50, 44)),
    #                               ((14, 4), (19, 10), (26, 16), (36, 21), (36, 26), (43, 35), (52, 40), (60, 47)),
    #                               ((17, 6), (20, 12), (32, 16), (39, 26), (44, 31), (47, 36), (61, 45), (70, 52)),
    #                               ((28, 7), (30, 15), (37, 20), (54, 26), (54, 31), (57, 37), (64, 42), (70, 41)),
    #                               ((40, 5), (41, 13), (46, 19), (54, 24), (64, 30), (67, 32), (72, 36), (79, 41)),
    #                               ((50, 1), (50, 10), (54, 17), (61, 23), (76, 27), (76, 33), (78, 36), (86, 38)),
    #                               ((65, 1), (62, 13), (63, 20), (69, 23), (80, 30), (83, 33), (87, 35), (92, 38)))
    # mppc.cellTranslation.value = (((17, 30), (17, 26), (17, 25), (17, 25), (17, 25), (17, 28), (19, 14), (19, 15)), ((15, 18), (16, 22), (17, 19), (17, 19), (17, 19), (17, 21), (19, 14), (19, 14)), ((15, 15), (16, 18), (17, 15), (17, 13), (9, 10), (14, 15), (14, 15), (14, 15)), ((15, 6), (16, 16), (17, 15), (17, 7), (17, 7), (31, 7), (8, 15), (8, 15)), ((15, 6), (8, 7), (23, 3), (23, 3), (23, 1), (23, 1), (23, 0), (23, 0)), ((8, 8), (8, 6), (9, 6), (9, 2), (18, 0), (17, 0), (17, 0), (20, 0)), ((8, 9), (8, 9), (9, 6), (9, 2), (9, 0), (9, 0), (9, 0), (9, 0)), ((8, 7), (8, 6), (7, 6), (2, 2), (3, 0), (2, 0), (0, 0), (0, 0)))


# def acquire_megafield(ccd, stage, multibeam, descanner, mppc, beamshift, mm, field_images, dwell_time):
def acquire_megafield(ccd, stage, multibeam, descanner, mppc, beamshift, field_images, dwell_time):

    # calculate lookup table for stage positions
    overlap = 0.075  # the fractional overlap
    # TODO positions should be calculated on the fly and not via lookuptable
    # create grid of field positions in x and y in meter
    # negative x indices as negative move with stage-bare moves feature to the left and thus moved right on sample
    # (2nd field should be right of first field)
    # transpose y as it is column vecto TODO why?
    # TODO why is there a multiplication at the end with the opposite axis?
    # TODO do we rotate the stage positions in the mp coordinate system into the stage-bare system -confirm!
    x_stage_pos = - numpy.array([numpy.linspace(0, ((field_images[0] - 1) * 3.195e-6 * 8 * (1 - overlap)), field_images[0])] * field_images[1])
    y_stage_pos = + numpy.array([numpy.linspace(0, ((field_images[1] - 1) * 3.195e-6 * 8 * (1 - overlap)), field_images[1])] * field_images[0]).transpose()

    # # rotate the field positions in the multiprobe coordinate system into the stage-bare coordinate system
    # # Note: rotation of axes needs the inverse rotation matrix (left-handed)
    # # TODO this should be handled by the component; is this the stage to mp calibration? use stage-scan????
    # mpp_angle = 1
    # x_stage_pos_rot = numpy.cos(numpy.radians(mpp_angle)) * x_stage_pos + numpy.sin(numpy.radians(mpp_angle)) * y_stage_pos
    # y_stage_pos_rot = numpy.cos(numpy.radians(mpp_angle)) * y_stage_pos - numpy.sin(numpy.radians(mpp_angle)) * x_stage_pos
    #
    # # add an offset (aka current stage-bare positions) to the field positions in the stage-bare coordinate system
    # x_stage_pos = x_stage_pos_rot + stage.position.value['x']
    # y_stage_pos = y_stage_pos_rot + stage.position.value['y']

    # add an offset (aka current stage-scan positions) to the field positions in the stage-scan coordinate system
    x_stage_pos = x_stage_pos + stage.position.value['x']
    y_stage_pos = y_stage_pos + stage.position.value['y']

    # save start position of megafield
    start_pos = stage.position.value

    mppc.dataContent.value = "empty"
    mppc.filename.value = time.strftime("testing_megafield_id-%Y-%m-%d-%H-%M-%S")
    dataflow = mppc.data

    # start recording data
    dataflow.subscribe(on_field_image)

    for row in range(field_images[1]):  # y axis, move stage on y,
        for col in range(field_images[0]):  # x axis, move stage on x

            # Correct for paramagnetic field of the stage using the beam shift.
            if col > 0 or row > 0:  # skip the first image as not needed for it
                correct_stage_magnetic_field(ccd, beamshift)

            # Note: move the stage by absolute position based on first field image position -> calculate next stage
            # position (to avoid accumulation of relative stage movement errors)
            # Sample should move on TSF screen to the right if we apply a positive x stage move.
            # (0,0) is top left corner -> move sample to the left, which is negative x stage move.
            if col > 0 or row > 0:
                logging.debug("moving the stage")
                stage.moveAbs({"x": x_stage_pos[row, col]}).result()  # in meter!
                stage.moveAbs({"y": y_stage_pos[row, col]}).result()  # in meter!
                time.sleep(0.4)  # TODO needed for settling of stage?

                # correct stage position
                # correct_stage_pos(stage, mm, beamshift, row, col, x_stage_pos, y_stage_pos,
                # x_stage_pos_rot, y_stage_pos_rot)

            print("Acquire field image with col (x): {} and row (y): {} at stage position: {}"
                  .format(col, row, stage.position.value))

            # clear event that is set in the callback
            image_received.clear()
            # request next field
            dataflow.next((col, row))   # acquire the field image (y, x)
            # Note: 1 or 2 additional seconds were not enough
            # TODO these extra seconds should be independent of the dwell_time, like flyback etc.
            if not image_received.wait(dwell_time * mppc.cellCompleteResolution.value[0] * mppc.cellCompleteResolution.value[1] + 60):
                # wait returns the current status; if status not True then no data received by callback function and
                # event was not set to true, but False (timeout) -> raise error
                raise TimeoutError("Did not receive field image in time! Timed out during field number row %s "
                                   "and col %s." % (row, col))

            # TODO investigate order!! with testcases in wrapper order not checked with filename pattern

    dataflow.unsubscribe(on_field_image)

    # Move stage back to starting position of the stage
    stage.moveAbs(start_pos).result()


image_received = threading.Event()


def on_field_image(df, da):
    """
    Subscriber for test cases which counts the number of times it is notified.
    *args contains the image/data which is received from the subscriber.
    """
    print("image received")
    image_received.set()


def main(args):

    logging.getLogger().setLevel(logging.DEBUG)

    if get_backend_status() != BACKEND_RUNNING:
        raise ValueError("Backend is not running.")

    ccd = model.getComponent(role="diagnostic-ccd")
    mppc = model.getComponent(role="mppc")
    dataflow = mppc.data
    descanner = model.getComponent(role="descanner")
    multibeam = model.getComponent(role="multibeam")
    scanner = model.getComponent(role="e-beam")
    # stage = model.getComponent(role="stage-bare")  # TODO use stage-scan!!???
    stage = model.getComponent(role="stage-scan")
    # mm = model.getComponent(role="stage-pos")
    beamshift = model.getComponent(role="ebeam-shift")

    # check that the MM is referenced
    # while "x" not in mm.position.value.keys():
    #     logging.debug("Wait a bit for the metrology module to finish referencing. Can take up to 8 minutes.")
    #     time.sleep(1)
    # logging.debug("Metrology module is referenced.")

    scanner.blanker.value = False  # unblank the beam
    scanner.horizontalFoV.value = 2.2e-05  # corresponds to mag = 5000x in quad view
    scanner.multiBeamMode.value = True  # True selects multibeam mode
    scanner.immersion.value = True  # put the microscope in immersion mode if not already
    scanner.external.value = True  # put microscope in external mode

    # Provide coordinate transform for beam shift.
    beamshift.updateMetadata({model.MD_CALIB: scanner.beamShiftTransformationMatrix.value})
    # reset beamshift
    beamshift.shift.value = (0, 0)

    # adjust megafield size
    field_images = (1, 1)   # (x, y) Note: x (horizontal) = col, y (vertical) = row
    dwell_time = 10e-6
    # external storage: <row>_<col>_<zl>.tiff
    # debugging: (1,3): expect 3 images in y direction (rows), files 0_0_, 1_0_, 2_0_
    # debugging: (2,3): expect 3 images in y direction (rows), 2 images in x dir (columns) files 0_0_, 1_0_, 2_0_, 0_1_, 1_1_, 2_1_

    try:
        # clear event that is set in the callback
        image_received.clear()
        # align detector with scanner
        mppc2mp(ccd, multibeam, descanner, mppc, dataflow)
        # adjust settings for megafield image acquisition
        settings_megafield(multibeam, descanner, mppc, dwell_time)
        # acquire the image data
        # acquire_megafield(ccd, stage, multibeam, descanner, mppc, beamshift, mm, field_images, dwell_time)
        acquire_megafield(ccd, stage, multibeam, descanner, mppc, beamshift, field_images, dwell_time)
    except Exception as exp:
        logging.error("%s", exp, exc_info=True)
    finally:
        dataflow.unsubscribe(on_field_image)
        scanner.blanker.value = True  # blank the beam

    print("Done")


if __name__ == '__main__':
    ret = main(sys.argv)
    exit(ret)

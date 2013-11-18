# -*- coding: utf-8 -*-
'''
Created on 17 Jul 2012

@author: Éric Piel

Copyright © 2012 Éric Piel, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the terms 
of the GNU General Public License version 2 as published by the Free Software 
Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
Odemis. If not, see http://www.gnu.org/licenses/.
'''
from __future__ import division

import calendar
from libtiff import TIFF
import logging
import numpy
from numpy.polynomial import polynomial
from odemis import model
import odemis
import re
import time

import libtiff.libtiff_ctypes as T # for the constant names
import xml.etree.ElementTree as ET


#pylint: disable=E1101
# Note about libtiff: it's a pretty ugly library, with 2 different wrappers.
# We use only the C wrapper (not the Python implementation).
# Note concerning the image format: it follows the numpy convention. The first
# dimension is the height, and second one is the width. (This is so because
# in memory the height is the slowest changing dimension, so it is first in C
# order.)
# So an RGB image of W horizontal pixels and H vertical pixels is an
# array of shape (H, W, 3). It is recommended to have the image in memory in C
# order (but that should not matter). For raw data, the convention is to have
# 5 dimensions in the order CTZYX.
# PIL and wxPython have images with the size expressed as (width, height), although
# in memory it corresponds to the same representation.
# User-friendly name
FORMAT = "TIFF"
# list of file-name extensions possible, the first one is the default when saving a file
EXTENSIONS = [".ome.tiff", ".ome.tif", ".tiff", ".tif"]

# We try to make it as much as possible looking like a normal (multi-page) TIFF,
# with as much metadata as possible saved in the known TIFF tags. In addition,
# we ensure it's compatible with OME-TIFF, which support much more metadata, and
# more than 3 dimensions for the data.
# Note that it could be possible to use more TIFF tags, from DNG and LSM for
# example.
# See also MicroManager OME-TIFF format with additional tricks for ImageJ:
# http://micro-manager.org/wiki/Micro-Manager_File_Formats

# TODO: make sure that _all_ the metadata is saved, either in TIFF tags, OME-TIFF,
# or in a separate mechanism.

def _convertToTiffTag(metadata):
    """
    Converts DataArray tags to libtiff tags.
    metadata (dict of tag -> value): the metadata of a DataArray
    returns (dict of tag -> value): the metadata as compatible for libtiff
    """
    # Note on strings: We accept unicode and utf-8 encoded strings, but TIFF
    # only supports ASCII (7-bits) officially. It seems however that utf-8
    # strings are often accepted, so for now we use that. If it turns out to
    # bring too much problem we might need .encode("ascii", "ignore") or
    # unidecode (but they are lossy).

    tiffmd = {}
    # we've got choice between inches and cm... so it's easy
    tiffmd[T.TIFFTAG_RESOLUTIONUNIT] = T.RESUNIT_CENTIMETER
    for key, val in metadata.items():
        if key == model.MD_SW_VERSION:
            tiffmd[T.TIFFTAG_SOFTWARE] = odemis.__shortname__ + " " + val.encode("utf-8")
        elif key == model.MD_HW_NAME:
            tiffmd[T.TIFFTAG_MAKE] = val.encode("utf-8")
        elif key == model.MD_HW_VERSION:
            tiffmd[T.TIFFTAG_MODEL] = val.encode("utf-8")
        elif key == model.MD_ACQ_DATE:
            tiffmd[T.TIFFTAG_DATETIME] = time.strftime("%Y:%m:%d %H:%M:%S", time.gmtime(val))
        elif key == model.MD_PIXEL_SIZE:
            # convert m/px -> px/cm
            try:
                tiffmd[T.TIFFTAG_XRESOLUTION] = (1 / val[0]) / 100
                tiffmd[T.TIFFTAG_YRESOLUTION] = (1 / val[1]) / 100
            except ZeroDivisionError:
                logging.debug("Pixel size tag is incorrect: %r", val)
        elif key == model.MD_POS:
            # convert m -> cm
            # XYPosition doesn't support negative values. So instead, we shift
            # everything by 1 m, which should be enough as samples are typically
            # a few cm big. The most important is that the image positions are
            # correct relatively to each other (for a given sample).
            tiffmd[T.TIFFTAG_XPOSITION] = 100 + val[0] * 100
            tiffmd[T.TIFFTAG_YPOSITION] = 100 + val[1] * 100
        elif key == model.MD_ROTATION:
            # TODO: should use the coarse grain rotation to update Orientation
            # and update rotation information to -45< rot < 45 -> maybe GeoTIFF's ModelTransformationTag?
            # or actually rotate the data?
            logging.info("Metadata tag '%s' skipped when saving TIFF file", key)
        # TODO MD_BPP : the actual bit size of the detector
        # Use SMINSAMPLEVALUE and SMAXSAMPLEVALUE ?
        # N = SPP in the specification, but libtiff duplicates the values
        elif key == model.MD_DESCRIPTION:
            # We don't use description as it's used for OME-TIFF
            tiffmd[T.TIFFTAG_PAGENAME] = val.encode("utf-8")
        # TODO save the brightness and contrast applied by the user?
        # Could use GrayResponseCurve, DotRange, or TransferFunction?
        # TODO save the tint applied by the user? maybe WhitePoint can help
        # TODO save username as "Artist" ? => not gonna fly if the user is "odemis"
        else:
            logging.debug("Metadata tag '%s' skipped when saving TIFF metadata", key)

    return tiffmd

def _GetFieldDefault(tfile, tag, default=None):
    """
    Same as TIFF.GetField(), but if the tag is not defined, return default
    Note: the C libtiff has GetFieldDefaulted() which returns the dafault value
    of the specification, this function is different.
    tag (int or string): tag id or name
    default (value): value to return if the tag is not defined
    """
    ret = tfile.GetField(tag)
    if ret is None:
        return default
    else:
        return ret

# factor for value -> m
resunit_to_m = {T.RESUNIT_INCH: 0.0254, T.RESUNIT_CENTIMETER: 0.01}
def _readTiffTag(tfile):
    """
    Reads the tiff tags of the current page and convert them into metadata
    It tries to do the reverse of _convertToTiffTag(). Support for other
    metadata and other ways to encode metadata is best-effort only.
    tfile (TIFF): the opened tiff file
    return (dict of tag -> value): the metadata of a DataArray
    """
    md = {}

    # scale + position
    resunit = _GetFieldDefault(tfile, T.TIFFTAG_RESOLUTIONUNIT, T.RESUNIT_INCH)
    factor = resunit_to_m.get(resunit, 1) # default to 1

    xres = tfile.GetField(T.TIFFTAG_XRESOLUTION)
    yres = tfile.GetField(T.TIFFTAG_YRESOLUTION)
    if xres is not None and yres is not None:
        try:
            md[model.MD_PIXEL_SIZE] = (factor / xres, factor / yres)
        except ZeroDivisionError:
            pass

    xpos = tfile.GetField(T.TIFFTAG_XPOSITION)
    ypos = tfile.GetField(T.TIFFTAG_YPOSITION)
    if xpos is not None and ypos is not None:
        # -1 m for the shift
        md[model.MD_POS] = (factor * xpos - 1, factor * ypos - 1)

    # informative metadata
    val = tfile.GetField(T.TIFFTAG_PAGENAME)
    if val is not None:
        md[model.MD_DESCRIPTION] = val
    val = tfile.GetField(T.TIFFTAG_SOFTWARE)
    if val is not None:
        md[model.MD_SW_VERSION] = val
    val = tfile.GetField(T.TIFFTAG_MAKE)
    if val is not None:
        md[model.MD_HW_NAME] = val
    val = tfile.GetField(T.TIFFTAG_MODEL)
    if val is not None:
        md[model.MD_HW_VERSION] = val
    val = tfile.GetField(T.TIFFTAG_DATETIME)
    if val is not None:
        try:
            t = calendar.timegm(time.strptime(val, "%Y:%m:%d %H:%M:%S"))
            md[model.MD_ACQ_DATE] = t
        except (OverflowError, ValueError):
            logging.info("Failed to parse date '%s'", val)

    return md

def _isThumbnail(tfile):
    """
    Detects whether the current image if a file is a thumbnail or not
    returns (boolean): True if the image is a thumbnail, False otherwise
    """
    # Best method is to check for the official TIFF flag
    subft = _GetFieldDefault(tfile, T.TIFFTAG_SUBFILETYPE, 0)
    if subft & T.FILETYPE_REDUCEDIMAGE:
        return True

    return False

def _guessModelName(das):
    """
    Detect the model of the Delmic microscope from the type of images acquired
    das (list of DataArrays)
    return (String or None): None if unknown, or the name of the model
    """
    # If any image has MD_IN_WL => fluorescence => SECOM
    # If any image has MD_WL_* => spectra => SPARC
    # If any image has MD_AR_POLE => angular resolved => SPARC
    for da in das:
        md = da.metadata
        if (model.MD_WL_LIST in md or model.MD_WL_POLYNOMIAL in md or
            model.MD_AR_POLE in md):
            return "SPARC"
        elif model.MD_IN_WL in md:
            return "SECOM"

    return None

def _indent(elem, level=0):
    """
    In-place pretty-print formatter
    From http://effbot.org/zone/element-lib.htm#prettyprint
    elem (ElementTree)
    """
    i = u"\n" + level * u"    "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + u"    "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            _indent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

_ROI_NS = "http://www.openmicroscopy.org/Schemas/ROI/2012-06"

def _convertToOMEMD(images):
    """
    Converts DataArray tags to OME-TIFF tags.
    images (list of DataArrays): the images that will be in the TIFF file, in order
      They should have 5 dimensions in this order: CTZYX, with the exception that
      all the first dimensions of size 1 can be skipped.
    returns (string): the XML data as compatible with OME
    Note: the images will be considered from the same detectors only if they are
      consecutive (and the metadata confirms the detector is the same)
    """
    # An OME-TIFF is a TIFF file with one OME-XML metadata embedded. For an
    # overview of OME-XML, see:
    # http://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2012-06/ome.html

    # It is not very clear in OME how to express that it's different acquisitions
    # of the _same_ sample by different instruments. However, our interpretation
    # of the format is the following:
#    + OME
#      + Experiment        # To describe the type of microscopy
#      + Experimenter      # To describe the user
#      + Instrument (*)    # To describe the acquisition technical details for each
#        + Microscope      # set of emitter/detector.
#        + LightSource (*)
#        + Detector (*)
#        + Objective (*)
#        + Filter (*)
#      + Image (*)         # To describe a set of images by the same instrument
#        + Description     # Not sure what to put (Image has "Name" attribute) => simple user note?
#        + AcquisitionDate # time of acquisition of the (first) image
#        + ExperimentRef
#        + ExperimenterRef
#        + InstrumentRef
#        + ImagingEnvironment # To describe the physical conditions (temp...)
#        + Pixels          # technical dimensions of the images (XYZ, T, C)
#          + Channel (*)   # emitter settings for the given channel (light wavelength)
#            + DetectorSettings
#          + Plane (*)     # physical dimensions/position of each images
#          + TiffData (*)  # where to find the data in the tiff file (IFD)
#                          # we explicitly reference each DataArray to avoid
#                          # potential ordering problems of the "Image" elements
#        + ROIRef (*)
#      + ROI (*)
#        . ID
#        . Name
#        + Union
#          + Shape (*)
#            . ID
#            . TheC
#            + Point
#                . X
#                . Y

# For fluorescence microscopy, the schema usage is straight-forwards. For SEM
# images, some metadata might be difficult to fit. For spectrum acquisitions
# (i.e., a C11YX cube), we use the channels to encode each spectrum.
# For AR acquisition, CCDX and CCDY are mapped to X/Y, and AR_POLE is recorded
# as an ROI Point with Name "PolePosition" and coordinates in pixel (like AR_POLE).

    # To create and manipulate the XML, we use the Python ElementTree API.

    # Note: pylibtiff has a small OME support, but it's so terrible that we are
    # much better ignoring it completely
    root = ET.Element('OME', attrib={
            "xmlns": "http://www.openmicroscopy.org/Schemas/OME/2012-06",
            "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
            "xsi:schemaLocation": "http://www.openmicroscopy.org/Schemas/OME/2012-06 http://www.openmicroscopy.org/Schemas/OME/2012-06/ome.xsd",
            })
    com_txt = ("Warning: this comment is an OME-XML metadata block, which "
               "contains crucial dimensional parameters and other important "
               "metadata. Please edit cautiously (if at all), and back up the "
               "original data before doing so. For more information, see the "
               "OME-TIFF web site: http://ome-xml.org/wiki/OmeTiff.")
    root.append(ET.Comment(com_txt))

    # add the microscope
    instr = ET.SubElement(root, "Instrument", attrib={
                                      "ID": "Instrument:0"})
    micro = ET.SubElement(instr, "Microscope", attrib={
                               "Manufacturer": "Delmic",
                                })
    model_name = _guessModelName(images)
    if model_name:
        micro.attrib["Model"] = model_name

    # for each set of images from the same instrument, add them
    groups = _findImageGroups(images)

    # Detectors
    for ifd, g in groups.items():
        did = ifd # our convention: ID is the first IFD
        da0 = g[0]
        if model.MD_HW_NAME in da0.metadata:
            detect = ET.SubElement(instr, "Detector", attrib={
                                      "ID": "Detector:%d" % did,
                                      "Model": da0.metadata[model.MD_HW_NAME]})

    # Objectives
    for ifd, g in groups.items():
        oid = ifd # our convention: ID is the first IFD
        da0 = g[0]
        if model.MD_LENS_MAG in da0.metadata:
            obj = ET.SubElement(instr, "Objective", attrib={
                      "ID": "Objective:%d" % oid,
                      "CalibratedMagnification": "%f" % da0.metadata[model.MD_LENS_MAG]
                      })
            if model.MD_LENS_NAME in da0.metadata:
                obj.attrib["Model"] = da0.metadata[model.MD_LENS_NAME]

    # TODO: filters (with TransmittanceRange)

    rois = {} # dict str->ET.Element (ROI ID -> ROI XML element)
    for ifd, g in groups.items():
        _addImageElement(root, g, ifd, rois)

    # ROIs have to come _after_ images, so add them only now
    root.extend(rois.values())

    # TODO add tag to each image with "Odemis", so that we can find them back
    # easily in a database?

    # make it more readable
    _indent(root)
    ometxt = ('<?xml version="1.0" encoding="UTF-8"?>' +
              ET.tostring(root, encoding="utf-8"))
    return ometxt

def _findElementByID(root, eid, tag=None):
    """
    Find the element with the given ID.
    Note: OME conformant documents cannot have multiple elements with the same
      ID. It is assumed to be correct. 
    root (ET.Element): the root element to start the search
    eid (str): the ID to match
    tag (str or None): the tag of the element to match. If None, any element
      with the given ID will be looked for.
    return (ET.Element): the element with .ID=eid, or None
    """
    for el in root.iter(tag):
        if "ID" in el.attrib and el.attrib["ID"] == eid:
            return el
         
    return None

def _updateMDFromOME(root, das):
    """
    Updates the metadata of DAs according to OME XML
    root (ET.Element): the root (i.e., OME) element of the XML description
    data (list of DataArrays): DataArrays at the same place as the TIFF IFDs
    return None: only the metadata of DA's inside is updated 
    """
    # For each Image in the XML, gorge ourself from all the metadata we can
    # find, and then use it to update the metadata of each IFD referenced.
    for ime in root.findall("Image"):
        md = {}
        try:
            md[model.MD_DESCRIPTION] = ime.attrib["Name"]
        except KeyError: 
            pass

        try:
            md[model.MD_USER_NOTE] = ime.attrib["Description"]
        except KeyError:
            pass

        acq_date = ime.find("AcquisitionDate")
        if acq_date is not None:
            try:
                # the earliest time of all the acquisitions in this image (see DeltaT)
                val = acq_date.text
                md[model.MD_ACQ_DATE] = calendar.timegm(time.strptime(val, "%Y-%m-%dT%H:%M:%S"))
            except (OverflowError, ValueError):
                pass

        objse = ime.find("ObjectiveSettings")
        try:
            obje = _findElementByID(root, objse.attrib["ID"], "Objective")
            mag = obje.attrib["CalibratedMagnification"]
            md[model.MD_LENS_MAG] = float(mag)
        except (AttributeError, KeyError, ValueError):
            pass

        pxe = ime.find("Pixels") # there must be only one per Image
        try:
            psx = float(pxe.attrib["PhysicalSizeX"]) * 1e-6 # µm -> m
            psy = float(pxe.attrib["PhysicalSizeY"]) * 1e-6
            md[model.MD_PIXEL_SIZE] = (psx, psy)
        except (KeyError, ValueError):
            pass

        try:
            md[model.MD_BPP] = int(pxe.attrib["SignificantBits"])
        except (KeyError, ValueError):
            pass

        ctz_2_ifd = _getIFDsFromOME(pxe)

        # Channels are a bit tricky, because apparently they are associated to
        # each C only by the order they are specified.
        wl_list = [] # we'll know it only once all the channels are passed
        chan = 0
        for che in pxe.findall("Channel"):
            mdc = {}
            try:
                mdc[model.MD_DESCRIPTION] = che.attrib["Name"]
            except KeyError:
                pass

            try:
                iwl = float(che.attrib["ExcitationWavelength"]) * 1e-9 # nm -> m
                mdc[model.MD_IN_WL] = (iwl, iwl)
            except (KeyError, ValueError):
                pass

            try:
                owl = float(che.attrib["EmissionWavelength"]) * 1e-9 # nm -> m
                if che.attrib["AcquisitionMode"] == "SpectralImaging":
                    # Spectrum => on the whole data cube
                    wl_list.append(owl)
                else:
                    # Fluorescence
                    mdc[model.MD_OUT_WL] = (owl, owl)
            except (KeyError, ValueError):
                pass

            try:
                hex_str = che.attrib["Color"] # hex string
                hex_str = hex_str[-8:] # almost copy of conversion.hex_to_rgb
                tint = tuple(int(hex_str[i:i + 2], 16) for i in [0, 2, 4, 6])
                mdc[model.MD_USER_TINT] = tint[:3] # only RGB
            except (KeyError, ValueError):
                pass

            # TODO: parse detector

            d_settings = che.find("DetectorSettings")
            if d_settings is not None:
                try:
                    bin_str = d_settings.attrib["Binning"]
                    m = re.match("(?P<b1>\d+)\s*x\s*(?P<b2>\d+)", bin_str)
                    mdc[model.MD_BINNING] = (int(m.group("b1")), int(m.group("b2")))
                except KeyError:
                    pass
                try:
                    mdc[model.MD_GAIN] = float(d_settings.attrib["Gain"])
                except (KeyError, ValueError):
                    pass
                try:
                    ror = float(d_settings.attrib["ReadOutRate"]) # MHz
                    mdc[model.MD_READOUT_TIME] = 1e-6 / ror # s
                except (KeyError, ValueError):
                    pass                

            # update all the IFDs related to this channel
            for ifd in ctz_2_ifd[chan].flat:
                if ifd == -1:
                    continue # no IFD known, it's alright, might be just 3D array
                da = das[ifd]
                if da is None:
                    continue # might be a thumbnail, it's alright
                # First apply the global MD, then per-channel
                da.metadata.update(md)
                da.metadata.update(mdc)

            chan += 1

        # Update metadata of each da, so that they will be merged
        if wl_list:
            if len(wl_list) != chan:
                logging.warning("WL_LIST has length %d, while expected %d",
                                len(wl_list), chan)
            for ifd in ctz_2_ifd.flat:
                if ifd == -1:
                    continue
                da = das[ifd]
                if da is None:
                    continue
                da.metadata.update({model.MD_WL_LIST: wl_list})


        # Plane (= one per CTZ -> IFD)
        for ple in pxe.findall("Plane"):
            mdp = {}
            pos = []
            try:
                for d in "CTZ": # that's our fixed order
                    ds = int(ple.attrib["The%s" % d]) # required tag
                    pos.append(ds)
            except KeyError:
                logging.warning("Failed to parse Plane element, skipping metadata")
                continue
            
            try:
                deltat = float(ple.attrib["DeltaT"]) # s
                mdp[model.MD_ACQ_DATE] = md[model.MD_ACQ_DATE] + deltat
            except (KeyError, ValueError):
                pass
            
            try:
                # FIXME: could actually be the dwell time (if scanned)
                mdp[model.MD_EXP_TIME] = float(ple.attrib["ExposureTime"]) # s
            except (KeyError, ValueError):
                pass

            try:
                # We assume it's in meters, as we write it (but there is no official unit)
                psx = float(pxe.attrib["PositionX"])
                psy = float(pxe.attrib["PositionY"])
                mdp[model.MD_POS] = (psx, psy)
            except (KeyError, ValueError):
                pass

            ifd = ctz_2_ifd[tuple(pos)]
            if ifd == -1:
                continue # no IFD known, it's alright, might be just 3D array
            da = das[ifd]
            if da is None:
                continue # might be a thumbnail, it's alright
            da.metadata.update(mdp)

def _getIFDsFromOME(pixele):
    """
    Return the IFD containing the 2D data for each high dimension of an array.
    Note: this doesn't take into account if the data is 3D.
    pixele (ElementTree): the element to Pixels of an image
    return (numpy.array of int): shape is the shape of the 3 high dimensions CTZ,
     the value is the IFD number of -1 if not specified.
    """
    hdshape = []
    for d in "CTZ": # that's our fixed order
        ds = int(pixele.get("Size%s" % d, "1"))
        hdshape.append(ds)

    imsetn = numpy.empty(hdshape, dtype="int")
    imsetn[:] = -1
    for tfe in pixele.findall("TiffData"):
        ifd = int(tfe.get("IFD", "0"))
        pos = []
        for d in "CTZ": # that's our fixed order
            ds = int(tfe.get("First%s" % d, "0"))
            pos.append(ds)
        # TODO: if no IFD specified, PC should default to all the IFDs
        pc = int(tfe.get("PlaneCount", "1"))

        # TODO: If has UUID tag, it means it's actually in a separate document
        # (cf tiff series) => load it and add it to the data as a new IFD?

        # If PlaneCount is > 1: it's in the same order as DimensionOrder
        # TODO: for now fixed to ZTC, but that should as DimensionOrder
        for i in range(pc):
            imsetn[tuple(pos)] = ifd + i
            if i == (pc - 1): # don't compute next position if it's over
                break
            # compute next position (with modulo)
            pos[-1] += 1
            for d in range(len(pos) - 1, -1, -1):
                if pos[d] >= hdshape[d]:
                    pos[d] = 0
                    pos[d - 1] += 1 # will fail if d = 0, on purpose

    return imsetn

# List of metadata which is allowed to merge (and possibly loose it partially)
WHITELIST_MD_MERGE = frozenset([model.MD_DESCRIPTION, model.MD_FILTER_NAME,
                                model.MD_HW_NAME, model.MD_HW_VERSION,
                                model.MD_SW_VERSION, model.MD_LENS_NAME,
                                model.MD_SENSOR_TEMP, model.MD_ACQ_DATE])
def _canBeMerged(das):
    """
    Check whether multiple DataArrays can be merged into a larger DA without 
      metadata loss.
    Note: this is about merging DataArrays for use in Odemis. For use in 
      OME-TIFF, the conditions are different.
    das (list of DataArrays): all the DataArrays
    return (boolean): True if they can be merged, False otherwise.
    """
    if len(das) <= 1:
        return True

    shape = das[0].shape
    md = das[0].metadata.copy()
    for da in das[1:]:
        # shape must be the same
        if shape != da.shape:
            return False
        # all the important metadata must be the same, or not defined
        for mdk in (set(md.keys()) | set(da.metadata.keys())):
            if (mdk not in WHITELIST_MD_MERGE and
                md.get(mdk) != da.metadata.get(mdk)):
                return False
        
    return True

def _mergeDA(das, hdim_index):
    """
    Merge multiple DataArrays into a higher dimension DataArray.
    das (list of DataArrays): ordered list of DataArrays (can contain more 
      arrays than what is used in the high dimension arrays
    hdim_index (ndarray of int >= 0): an array representing the higher 
      dimensions of the final merged arrays. Each value is the index of the 
      small array in das.
    return (DataArray): the merge of all the DAs. The shape is hdim_index.shape
     + shape of original DataArray. The metadata is the metadata of the first
     DataArray inserted
    """
    fim = das[hdim_index.flat[0]]
    tshape = hdim_index.shape + fim.shape
    imset = numpy.empty(tshape, fim.dtype)
    for hi, i in numpy.ndenumerate(hdim_index):
        imset[hi] = das[i]
    
    return model.DataArray(imset, metadata=fim.metadata)
    

def _foldArraysFromOME(root, das):
    """
    Reorganize DataArrays with more than 2 dimensions according to OME XML
    Note: it expects _updateMDFromOME has been run before and so each array
     has its metadata filled up.
    Note: Officially OME supports only base arrays of 2D. But we also support
     base arrays of 3D if the data is RGB (3rd dimension has length 3).
    root (ET.Element): the root (i.e., OME) element of the XML description
    data (list of DataArrays): DataArrays at the same place as the TIFF IFDs
    return (list of DataArrays): new shorter list of DAs 
    """
    omedas = []

    for ime in root.findall("Image"):
        pxe = ime.find("Pixels") # there must be only one per Image

        # The relation between channel and planes is not very clear. Each channel
        # can have multiple SamplesPerPixel, apparently to indicate they have
        # multiple planes. However, the Interleaved attribute is global for
        # Pixels, and seems to imply that RGB data could be saved as a whole,
        # although OME-TIFF normally only has 2D arrays.
        # So far the understanding is Channel refers to the "Logical channels",
        # and Plane refers to the C dimension.
#        spp = int(pxe.get("Channel/SamplesPerPixel", "1"))

        imsetn = _getIFDsFromOME(pxe)
        # For now we expect RGB as (SPP=3,) SizeC=3, PlaneCount=1, and 1 3D IFD,
        # or as (SPP=3,) SizeC=3, PlaneCount=3 and 3 2D IFDs.

        # Read the complete shape of the dataset
        dshape = list(imsetn.shape)
        for d in "YX":
            ds = int(pxe.get("Size%s" % d, "1"))
            dshape.append(ds)

        # Check if the IFDs are 2D or 3D, based on the first one
        fifd = imsetn[0, 0, 0]
        if fifd == -1:
            raise ValueError("Not all IFDs defined for image %d" % len(omedas) + 1)

        fim = das[fifd]
        if fim == None:
            continue # thumbnail
        is_3d = (len(fim.shape) == 3 and fim.shape[0] > 1)

        # Remove C dim if 3D
        if is_3d:
            if imsetn.shape != fim.shape[0]:
                # 3D data arrays are not officially supported in OME-TIFF anyway
                raise NotImplementedError("Loading of %d channel from images "
                       "with %d channels not supported" % (imsetn.shape, fim.shape[0]))
            imsetn = imsetn[0]

        # Short-circuit for dataset with only one IFD
        if all(d == 1 for d in imsetn.shape):
            omedas.append(fim)
            continue

        if -1 in imsetn:
            raise ValueError("Not all IFDs defined for image %d" % len(omedas) + 1)

        # TODO: Position might also be different. Probably should be grouped
        # by position too, as Odemis doesn't support such case.
        
        # In Odemis, arrays can be merged along C _only_ if they are continuous
        # (like for a spectrum acquisition). If it's like fluorescence, with
        # complex channel metadata, they need to be kept separated.
        # Check for all the C with T=0, Z=0 (should be the same at other indices)
        das_tz0n = list(imsetn[..., 0, 0])
        das_tz0 = [das[i] for i in das_tz0n]
        if not _canBeMerged(das_tz0):
            for sub_imsetn in imsetn:
                # Combine all the IFDs into a (1+)4D array    
                sub_imsetn.shape = (1,) + sub_imsetn.shape
                imset = _mergeDA(das, sub_imsetn)
                omedas.append(imset)
        else:
            # Combine all the IFDs into a 5D array
            imset = _mergeDA(das, imsetn)
            if is_3d:
                # move the C axis back to first position (it's currently TZCYX)
                imset = numpy.rollaxis(imset, -3)
            omedas.append(imset)

    return omedas

def _countNeededIFDs(da):
    """
    return the number of IFD (aka TIFF pages, aka planes) needed for storing the given array
    da (DataArray): can have any dimensions, should be ordered ...CTZYX
    return (int > 1)
    """
    # Storred as a sequence of 2D arrays... excepted if it contains RGB images,
    # then we store RGB images (i.e., 3D arrays).
    rep_hdim = list(da.shape[:-2])
    if da.ndim >= 5 and da.shape[-5] == 3: # RGB
        rep_hdim[-3] = 1
    return numpy.prod(rep_hdim)


def _findImageGroups(das):
    """
    Find groups of images which should be considered part of the same acquisition
    (aka "Image" in OME-XML).
    das (list of DataArray): all the images of the final TIFF file
    returns (dict int -> list of DataArrays):
        IFD (index) of the first DataArray of a group -> "group" of DataArrays
    """
    # We consider images to be part of the same group if they have:
    # * same shape
    # * metadata that show they were acquired by the same instrument
    groups = dict()
    current_ifd = 0
    prev_da = None
    for da in das:
        # check if it can be part of the current group (compare just to the previous DA)
        if (prev_da is None or
            prev_da.shape != da.shape or
            prev_da.metadata.get(model.MD_HW_NAME, None) != da.metadata.get(model.MD_HW_NAME, None) or
            prev_da.metadata.get(model.MD_HW_VERSION, None) != da.metadata.get(model.MD_HW_VERSION, None)
            ):
            # new group
            group_ifd = current_ifd
        groups.setdefault(group_ifd, []).append(da)
        
        # increase ifd by the number of planes
        current_ifd += _countNeededIFDs(da)
        prev_da = da
            
    return groups

def _dtype2OMEtype(dtype):
    """
    Converts a numpy dtype to a OME type
    dtype (numpy.dtype)
    returns (string): OME type
    """
    # OME type is one of : int8, int16, int32, uint8, uint16, uint32, float,
    # bit, double, complex, double-complex
    # Meaning is not explicit, but probably the same as TIFF/C:
    # Float is 4 bytes, while double is 8 bytes
    # So complex is 8 bytes, and double complex == 16 bytes
    if dtype.kind in "ui":
        # int and uint seems to be compatible in general
        return "%s" % dtype
    elif dtype.kind == "f":
        if dtype.itemsize <= 4:
            return "float"
        else:
            return "double"
    elif dtype.kind == "c":
        if dtype.itemsize <= 4:
            return "complex"
        else:
            return "double-complex"
    elif dtype.kind == "b":
        return "bit"
    else:
        raise NotImplementedError("data type %s is not support by OME" % dtype)
    

def _addImageElement(root, das, ifd, rois):
    """
    Add the metadata of a list of DataArray to a OME-XML root element
    root (Element): the root element
    das (list of DataArray): all the images to describe. Each DataArray
     can have up to 5 dimensions in the order CTZYX. IOW, RGB images have C=3. 
    ifd (int): the IFD of the first DataArray
    rois (dict str -> ET element): all ROIs added so, and will be updated as
      needed. 
    Note: the images in das must be added in the final TIFF in the same order
     and contiguously
    """
    assert(len(das) > 0)
    # all image have the same shape?
    assert all([das[0].shape == im.shape for im in das]) 

    idnum = len(root.findall("Image"))
    ime = ET.SubElement(root, "Image", attrib={"ID": "Image:%d" % idnum})

    # compute a common metadata
    globalMD = {}
    for da in das:
        globalMD.update(da.metadata)

    globalAD = None
    if model.MD_ACQ_DATE in globalMD:
        # Need to find the earliest time
        globalAD = min(d.metadata[model.MD_ACQ_DATE]
                        for d in das if model.MD_ACQ_DATE in d.metadata)
        ad = ET.SubElement(ime, "AcquisitionDate")
        ad.text = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(globalAD))

    # find out about the common attribute (Name)
    if model.MD_DESCRIPTION in globalMD:
        ime.attrib["Name"] = globalMD[model.MD_DESCRIPTION]

    # find out about the common sub-elements (time, user note, shape)
    if model.MD_USER_NOTE in globalMD:
        desc = ET.SubElement(ime, "Description")
        desc.text = globalMD[model.MD_USER_NOTE]

    if model.MD_LENS_MAG in globalMD:
        # add ref to Objective
        ose = ET.SubElement(ime, "ObjectiveSettings",
                            attrib={"ID": "Objective:%d" % ifd})


    # Find a dimension along which the DA can be concatenated. That's a
    # dimension which is of size 1. 
    # For now, if there are many possibilities, we pick the first one.
    da0 = das[0]
    
    dshape = das[0].shape 
    if len(dshape) < 5:
        dshape = [1] * (5-len(dshape)) + list(dshape)
    if not 1 in dshape:
        raise ValueError("No dimension found to concatenate images: %s" % dshape)
    concat_axis = dshape.index(1)
    # global shape (same as dshape, but the axis of concatenation is the number of images)
    gshape = list(dshape)
    gshape[concat_axis] = len(das)
    gshape = tuple(gshape)
    
    # Note: it seems officially OME-TIFF doesn't support RGB TIFF (instead,
    # each colour should go in a separate channel). However, that'd defeat
    # the purpose of the thumbnail, and it seems at OMERO handles this
    # not too badly (all the other images get 3 components).
    is_rgb = (dshape[0] == 3)

    # TODO: check that all the DataArrays have the same shape
    pixels = ET.SubElement(ime, "Pixels", attrib={
                              "ID": "Pixels:%d" % idnum,
                              "DimensionOrder": "XYZTC", # we don't have ZT so it doesn't matter
                              "Type": "%s" % _dtype2OMEtype(da0.dtype),
                              "SizeX": "%d" % gshape[4], # numpy shape is reversed
                              "SizeY": "%d" % gshape[3],
                              "SizeZ": "%d" % gshape[2],
                              "SizeT": "%d" % gshape[1],
                              "SizeC": "%d" % gshape[0],
                              })
    # Add optional values
    if model.MD_PIXEL_SIZE in globalMD:
        pxs = globalMD[model.MD_PIXEL_SIZE]
        pixels.attrib["PhysicalSizeX"] = "%f" % (pxs[0] * 1e6) # in µm
        pixels.attrib["PhysicalSizeY"] = "%f" % (pxs[1] * 1e6)
    
    if model.MD_BPP in globalMD:
        bpp = globalMD[model.MD_BPP]
        pixels.attrib["SignificantBits"] = "%d" % bpp # in bits

    # For each DataArray, add a Channel, TiffData, and Plane, but be careful
    # because they all have to be grouped and in this order.

    # For spectrum images, 1 channel per wavelength. As suggested on the
    # OME-devel mailing list, we use EmissionWavelength to store it.
    # http://trac.openmicroscopy.org.uk/ome/ticket/7355 mentions spectrum can
    # be stored as a fake Filter with the CutIn/CutOut wavelengths, but that's
    # not so pretty.
    if model.MD_WL_POLYNOMIAL in globalMD:
        # we store it by computing the actual values of each channel
        pn = globalMD[model.MD_WL_POLYNOMIAL]
        l = gshape[0]
        npn = polynomial.Polynomial(pn, domain=[0, l - 1], window=[0, l - 1])
        wl_list = npn.linspace(l)[1]
#        if len(pn) >= 1:
#            pixels.attrib["WaveStart"] = "%f" % (pn[0] * 1e9) # in nm
#        if len(pn) >= 2:
#            pixels.attrib["WaveIncrement"] = "%f" % (pn[1] * 1e9) # in nm/px
    elif model.MD_WL_LIST in globalMD:
        wl_list = globalMD[model.MD_WL_LIST]
        if len(wl_list) != da.shape[0]:
            logging.warning("WL_LIST metadata has length of %d, while expected "
                            "%d, skipping it", len(wl_list), da.shape[0])
            wl_list = None

#        inc = (lwl[-1] - lwl[0]) / (len(lwl) - 1)
#        pixels.attrib["WaveStart"] = "%f" % (lwl[0] * 1e9) # in nm
#        pixels.attrib["WaveIncrement"] = "%f" % (inc * 1e9) # in nm/px
    else:
        wl_list = None

    subid = 0
    for da in das:
        if is_rgb or len(da.shape) < 5:
            num_chan = 1
        else:
            num_chan = da.shape[0]
        for c in range(num_chan):
            chan = ET.SubElement(pixels, "Channel", attrib={
                                   "ID": "Channel:%d:%d" % (idnum, subid)})
            if is_rgb:
                chan.attrib["SamplesPerPixel"] = "%d" % dshape[0]
                
            # Name can be different for each channel in case of fluroescence
            if model.MD_DESCRIPTION in da.metadata:
                chan.attrib["Name"] = da.metadata[model.MD_DESCRIPTION]
    
            # TODO Fluor attrib for the dye?
            # TODO create a Filter with the cut range?
            if model.MD_IN_WL in da.metadata:
                iwl = da.metadata[model.MD_IN_WL]
                xwl = numpy.mean(iwl) * 1e9 # in nm
                chan.attrib["ExcitationWavelength"] = "%d" % round(xwl)
    
                # if input wavelength range is small, it means we are in epifluoresence
                if abs(iwl[1] - iwl[0]) < 100e-9:
                    chan.attrib["IlluminationType"] = "Epifluorescence"
                    chan.attrib["AcquisitionMode"] = "WideField"
                    chan.attrib["ContrastMethod"] = "Fluorescence"
                else:
                    chan.attrib["IlluminationType"] = "Epifluorescence"
                    chan.attrib["AcquisitionMode"] = "WideField"
                    chan.attrib["ContrastMethod"] = "Brightfield"
    
            if model.MD_OUT_WL in da.metadata:
                owl = da.metadata[model.MD_OUT_WL]
                ewl = numpy.mean(owl) * 1e9 # in nm
                chan.attrib["EmissionWavelength"] = "%d" % round(ewl)
    
            if wl_list is not None and len(wl_list) > 0:
                if model.MD_OUT_WL in da.metadata:
                    logging.warning("DataArray contains both OUT_WL (%s) and "
                                    "incompatible WL_LIST metadata",
                                    da.metadata[model.MD_OUT_WL])
                else:
                    chan.attrib["AcquisitionMode"] = "SpectralImaging"
                    # It should be an int, but that looses too much precision
                    chan.attrib["EmissionWavelength"] = "%f" % (wl_list[c] * 1e9)
                 
            if model.MD_USER_TINT in da.metadata:
                # user tint is 3 tuple int
                # color is hex RGBA (eg: #FFFFFFFF)
                tint = da.metadata[model.MD_USER_TINT]
                if len(tint) == 3:
                    tint = tint + (255,) # need alpha channel
                hex_str = "".join("%.2x" % c for c in tint) # copy of conversion.rgb_to_hex()
                chan.attrib["Color"] = "#%s" % hex_str


            # Add info on detector
            attrib = {}
            if model.MD_BINNING in da.metadata:
                attrib["Binning"] = "%dx%d" % da.metadata[model.MD_BINNING]
            if model.MD_GAIN in da.metadata:
                attrib["Gain"] = "%f" % da.metadata[model.MD_GAIN]
            if model.MD_READOUT_TIME in da.metadata:
                ror = (1 / da.metadata[model.MD_READOUT_TIME]) / 1e6 # MHz
                attrib["ReadOutRate"] = "%f" % ror
            if model.MD_EBEAM_ENERGY in da.metadata:
                # Schema only mentions PMT, but we use it for the e-beam too
                attrib["Voltage"] = "%f" % da.metadata[model.MD_EBEAM_ENERGY] # V
    
            if attrib:
                # detector of the group has the same id as first IFD of the group
                attrib["ID"] = "Detector:%d" % ifd
                ds = ET.SubElement(chan, "DetectorSettings", attrib=attrib)

            subid += 1

    # TiffData Element: describe every single IFD image
    # TODO: could be more compact for DAs of dim > 2, with PlaneCount = first dim > 1?
    subid = 0
    rep_hdim = list(gshape[:-2])
    if is_rgb:
        rep_hdim[0] = 1
    for index in numpy.ndindex(*rep_hdim):
        tde = ET.SubElement(pixels, "TiffData", attrib={
                                "IFD": "%d" % (ifd + subid),
                                "FirstC": "%d" % index[0],
                                "FirstT": "%d" % index[1],
                                "FirstZ": "%d" % index[2],
                                "PlaneCount": "1"
                                })
        subid += 1

    # Plane Element
    subid = 0
    for index in numpy.ndindex(*rep_hdim):
        da = das[index[concat_axis]]
        plane = ET.SubElement(pixels, "Plane", attrib={
                               "TheC": "%d" % index[0],
                               "TheT": "%d" % index[1],
                               "TheZ": "%d" % index[2],
                               })
        if model.MD_ACQ_DATE in da.metadata:
            diff = da.metadata[model.MD_ACQ_DATE] - globalAD
            plane.attrib["DeltaT"] = "%.12f" % diff

        if model.MD_EXP_TIME in da.metadata:
            exp = da.metadata[model.MD_EXP_TIME]
            plane.attrib["ExposureTime"] = "%.12f" % exp
        elif model.MD_DWELL_TIME in da.metadata:
            # typical for scanning techniques => more or less correct
            exp = da.metadata[model.MD_DWELL_TIME] * numpy.prod(da.shape)
            plane.attrib["ExposureTime"] = "%.12f" % exp

        # Note that Position has no official unit, which prevents Tiling to be
        # usable. In one OME-TIFF official example of tiles, they use pixels
        # (and ModuloAlongT "tile")
        if model.MD_POS in da.metadata:
            pos = da.metadata[model.MD_POS]
            plane.attrib["PositionX"] = "%.12f" % pos[0] # any unit is allowed => m
            plane.attrib["PositionY"] = "%.12f" % pos[1]

        subid += 1

    # ROIs (= ROIRefs + new ROI elements)
    # For now, we use them only to store AR_POLE metadata
    for chan, da in enumerate(das):
        # Note: we assume the
        if model.MD_AR_POLE in da.metadata:
            rid = _createPointROI(rois, "PolePosition",
                                  da.metadata[model.MD_AR_POLE],
                                  shp_attrib={"TheC": "%d" % chan})
            ET.SubElement(ime, "ROIRef", attrib={"xmlns": _ROI_NS, "ID": rid})

def _createPointROI(rois, name, p, shp_attrib=None):
    """
    Create a new Point ROI XML element and add it to the dict
    rois (dict str-> ET.Element): list of all current ROIs
    name (str or None): name of the ROI (if None -> no name)
    p (tuple of 2 floats): values to put in X and Y (arbitrary)
    shp_attrib (dict or None): attributes for the shape element
    return id (str): ROI ID for referencing
    """
    shp_attrib = shp_attrib or {}
    # Find an ID not yet used
    for n in range(len(rois), -1, -1):
        rid = "ROI:%d" % n
        if not rid in rois: # very likely
            shapeid = "Shape:%d" % n # for now, assume 1 shape <-> 1 ROI
            break
    else:
        raise IndexError("Couldn't find an available ID")

    # Create ROI/Union/Shape/Point
    roie = ET.Element("ROI", attrib={"xmlns": _ROI_NS, "ID": rid})
    if name is not None:
        roie.attrib["Name"] = name
    unione = ET.SubElement(roie, "Union")
    shapee = ET.SubElement(unione, "Shape", attrib={"ID": shapeid})
    shapee.attrib.update(shp_attrib)
    pointe = ET.SubElement(shapee, "Point", attrib={"X": "%f" % p[0],
                                                    "Y": "%f" % p[1]})

    # Add the element to all the ROIs
    rois[rid] = roie
    return rid


def _saveAsTiffLT(filename, data, thumbnail):
    """
    Saves a DataArray as a TIFF file.
    filename (string): name of the file to save
    data (ndarray): 2D data of int or float
    """
    _saveAsMultiTiffLT(filename, [data], thumbnail)

def _saveAsMultiTiffLT(filename, ldata, thumbnail, compressed=True):
    """
    Saves a list of DataArray as a multiple-page TIFF file.
    filename (string): name of the file to save
    ldata (list of DataArray): list of 2D data of int or float. Should have at least one array
    thumbnail (None or DataArray): see export
    compressed (boolean): whether the file is LZW compressed or not.
    """
    f = TIFF.open(filename, mode='w')

    # According to this page: http://www.openmicroscopy.org/site/support/file-formats/ome-tiff/ome-tiff-data
    # LZW is a good trade-off between compatibility and small size (reduces file
    # size by about 2). => that's why we use it by default
    if compressed:
        compression = "lzw"
    else:
        compression = None

    # OME tags: a XML document in the ImageDescription of the first image
    if thumbnail is not None:
        # OME expects channel as 5th dimension. If thumbnail is RGB as HxWx3,
        # reorganise as 3x1x1xHxW
        if len(thumbnail.shape) == 3:
            OME_thumbnail = numpy.rollaxis(thumbnail, 2)
            OME_thumbnail = OME_thumbnail[:,numpy.newaxis,numpy.newaxis,:,:] #5D
        else:
            OME_thumbnail = thumbnail
        alldata = [OME_thumbnail] + ldata
    else:
        alldata = ldata
    # TODO: reorder the data so that data from the same sensor are together

    ometxt = _convertToOMEMD(alldata)

    if thumbnail is not None:
        # save the thumbnail just as the first image
        # FIXME: Note that this is contrary to the specification which states:
        # "If multiple subfiles are written, the first one must be the
        # full-resolution image." The main problem is that most thumbnailers
        # use the first image as thumbnail. Maybe we should provide our own
        # clever thumbnailer?

        f.SetField(T.TIFFTAG_IMAGEDESCRIPTION, ometxt)
        ometxt = None
        f.SetField(T.TIFFTAG_PAGENAME, "Composited image")
        # Flag for saying it's a thumbnail
        f.SetField(T.TIFFTAG_SUBFILETYPE, T.FILETYPE_REDUCEDIMAGE)

        # Warning, upstream Pylibtiff has a bug: it can only write RGB images are
        # organised as 3xHxW, while normally in numpy, it's HxWx3.
        # Our version is fixed

        # write_rgb makes it clever to detect RGB vs. Greyscale
        f.write_image(thumbnail, compression=compression, write_rgb=True)
        

        # TODO also save it as thumbnail of the image (in limited size)
        # see  http://www.libtiff.org/man/thumbnail.1.html

        # libtiff.py doesn't support yet SubIFD's so it's not going to fly


        # from http://stackoverflow.com/questions/11959617/in-a-tiff-create-a-sub-ifd-with-thumbnail-libtiff
#        //Define the number of sub-IFDs you are going to write
#        //(assuming here that we are only writing one thumbnail for the image):
#        int number_of_sub_IFDs = 1;
#        toff_t sub_IFDs_offsets[1] = { 0UL };
#
#        //set the TIFFTAG_SUBIFD field:
#        if(!TIFFSetField(created_TIFF, TIFFTAG_SUBIFD, number_of_sub_IFDs,
#            sub_IFDs_offsets))
#        {
#            //there was an error setting the field
#        }
#
#        //Write your main image raster data to the TIFF (using whatever means you need,
#        //such as TIFFWriteRawStrip, TIFFWriteEncodedStrip, TIFFWriteEncodedTile, etc.)
#        //...
#
#        //Write your main IFD like so:
#        TIFFWriteDirectory(created_TIFF);
#
#        //Now here is the trick: like the comment in the libtiff source states, the
#        //next n directories written will be sub-IFDs of the main IFD (where n is
#        //number_of_sub_IFDs specified when you set the TIFFTAG_SUBIFD field)
#
#        //Set up your sub-IFD
#        if(!TIFFSetField(created_TIFF, TIFFTAG_SUBFILETYPE, FILETYPE_REDUCEDIMAGE))
#        {
#            //there was an error setting the field
#        }
#
#        //set the rest of the required tags here, as well as any extras you would like
#        //(remember, these refer to the thumbnail, not the main image)
#        //...
#
#        //Write this sub-IFD:
#        TIFFWriteDirectory(created_TIFF);

    for data in ldata:
        # TODO: see if we need to set FILETYPE_PAGE + Page number for each image? data?
        tags = _convertToTiffTag(data.metadata)
        if ometxt: # save OME tags if not yet done
            f.SetField(T.TIFFTAG_IMAGEDESCRIPTION, ometxt)
            ometxt = None

        # for data > 2D: write as a sequence of 2D images or RGB images
        if data.ndim == 5 and data.shape[0] == 3: # RGB
            # Write an RGB image, instead of 3 images along C
            write_rgb = True
            hdim = data.shape[1:3]
            data = numpy.rollaxis(data, 0, -2) # move C axis near YX
        else:
            write_rgb = False
            hdim = data.shape[:-2]
            
        for i in numpy.ndindex(*hdim):
            # Save metadata (before the image)
            for key, val in tags.items():
                f.SetField(key, val)
            f.write_image(data[i], write_rgb=write_rgb, compression=compression)

def _thumbsFromTIFF(filename):
    """
    Read thumbnails from an TIFF file.
    return (list of model.DataArray)
    """
    f = TIFF.open(filename, mode='r')
    # open each image/page as a separate data
    data = []
    f.SetDirectory(0)
    while True:
        if _isThumbnail(f):
            md = _readTiffTag(f) # reads tag of the current image
            image = f.read_image()
            da = model.DataArray(image, metadata=md)
            data.append(da)
        
        # TODO: also check SubIFD for sub directories that might contain 
        # thumbnails
        if f.ReadDirectory() == 0: # reads _next_ directory
            break

    return data

def _reconstructFromOMETIFF(xml, data):
    """
    Update DAs to reflect shape and metadata contained in OME XML
    xml (string): String containing the OME XML declaration
    data (list of model.DataArray): each
    return (list of model.DataArray): new list with the DAs following the OME 
      XML description. Note that DAs are either updated or completely recreated.  
    """
    # Remove "xmlns" which is the default namespace and is appended everywhere
    # It's not beautiful, but the simplest with ET to handle expected namespaces.
    xml = re.sub('xmlns="http://www.openmicroscopy.org/Schemas/OME/....-.."',
                 "", xml, count=1)
    root = ET.fromstring(xml)
    _updateMDFromOME(root, data)
    omedata = _foldArraysFromOME(root, data)

    return omedata

def _dataFromTIFF(filename):
    """
    Read microscopy data from a TIFF file.
    filename (string): path of the file to read
    return (list of model.DataArray)
    """
    f = TIFF.open(filename, mode='r')
    
    # open each image/page as a separate image
    data = []
    for image in f.iter_images():
        # If it's a thumbnail, skip it, but leave the space free to not mess with the IFD number
        if _isThumbnail(f):
            data.append(None)
            continue
        md = _readTiffTag(f) # reads tag of the current image
        da = model.DataArray(image, metadata=md)
        data.append(da)
        
    # If looks like OME TIFF, reconstruct >2D data and add metadata
    # It's OME TIFF, if it has a valid ome-tiff XML in the first T.TIFFTAG_IMAGEDESCRIPTION
    # Warning: we support what we write, not the whole OME-TIFF specification.
    f.SetDirectory(0)
    desc = f.GetField(T.TIFFTAG_IMAGEDESCRIPTION)
    if (desc and ((desc.startswith("<?xml") and "<ome " in desc.lower())
                  or desc[:4].lower() == '<ome')):
        try:
            data = _reconstructFromOMETIFF(desc, data)
        except Exception:
            # fallback to pretend there was no OME XML
            logging.exception("Failed to decode OME XML string: '%s'", desc)

    # Remove all the None (=thumbnails) from the list
    data = [i for i in data if i is not None]
    return data

def export(filename, data, thumbnail=None):
    '''
    Write a TIFF file with the given image and metadata
    filename (string): filename of the file to create (including path)
    data (list of model.DataArray, or model.DataArray): the data to export.
       Metadata is taken directly from the DA object. If it's a list, a multiple
       page file is created. It must have 5 dimensions in this order: Channel, 
       Time, Z, Y, X. However, all the first dimensions of size 1 can be omitted
       (ex: an array of 111YX can be given just as YX, but RGB images are 311YX,
       so must always be 5 dimensions).
    thumbnail (None or numpy.array): Image used as thumbnail for the file. Can be of any
      (reasonable) size. Must be either 2D array (greyscale) or 3D with last
      dimension of length 3 (RGB). If the exporter doesn't support it, it will
      be dropped silently.
    '''
    if isinstance(data, list):
        _saveAsMultiTiffLT(filename, data, thumbnail)
    else:
        # TODO should probably not enforce it: respect duck typing
        assert(isinstance(data, model.DataArray))
        _saveAsTiffLT(filename, data, thumbnail)

def read_data(filename):
    """
    Read an TIFF file and return its content (skipping the thumbnail).
    filename (string): filename of the file to read
    return (list of model.DataArray): the data to import (with the metadata 
     as .metadata). It might be empty.
     Warning: reading back a file just exported might give a smaller number of
     DataArrays! This is because export() tries to aggregate data which seems
     to be from the same acquisition but on different dimensions C, T, Z.
     read_data() cannot separate them back explicitly.
    raises:
        IOError in case the file format is not as expected.
    """
    # TODO: support filename to be a File or Stream (but it seems very difficult
    # to do it without looking at the .filename attribute)
    # see http://pytables.github.io/cookbook/inmemory_hdf5_files.html

    return _dataFromTIFF(filename)

def read_thumbnail(filename):
    """
    Read the thumbnail data of a given TIFF file.
    filename (string): filename of the file to read
    return (list of model.DataArray): the thumbnails attached to the file. If 
     the file contains multiple thumbnails, all of them are returned. If it 
     contains none, an empty list is returned.
    raises:
        IOError in case the file format is not as expected.
    """
    # TODO: support filename to be a File or Stream

    return _thumbsFromTIFF(filename)

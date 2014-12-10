# -*- coding: utf-8 -*-
"""
Created on 17 Aug 2012

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

"""
from __future__ import division

# for listing all the types of file format supported
import importlib
import logging
from odemis.dataio import tiff
import os


# The interface of a "format manager" is as follows:
#  * one module
#  * FORMAT (string): user friendly name of the format
#  * EXTENSIONS (list of strings): possible file-name extensions
#  * export (callable): write model.DataArray into a file
#  * get_data (callable): read a file into model.DataArray
#  if it doesn't support writing, then is has no .export(), and if it doesn't
#  support reading, then it has not get_data().
__all__ = ["tiff", "hdf5", "png"]


def get_available_formats(mode=os.O_RDWR, allowlossy=False):
    """
    Find the available file formats

    :param mode: (os.O_RDONLY, os.O_WRONLY, or os.O_RDWR): whether only list
        formats which can be read, which can be written, or all of them.
    :returns: (dict string -> list of strings): name of each format -> list of
        extensions
    """
    formats = {}
    # Look dynamically which format is available
    for module_name in __all__:
        try:
            exporter = importlib.import_module("." + module_name, "odemis.dataio")
        except Exception: #pylint: disable=W0703
            continue # module cannot be loaded
        if not allowlossy and hasattr(exporter, "LOSSY") and exporter.LOSSY:
            logging.debug("Skipping exporter %s as it is lossy", module_name)
            continue
        if ((mode == os.O_RDONLY and not hasattr(exporter, "read_data")) or
            (mode == os.O_WRONLY and not hasattr(exporter, "export"))):
            continue
        formats[exporter.FORMAT] = exporter.EXTENSIONS

    if not formats:
        logging.error("No file converter found!")

    return formats


def get_converter(fmt):
    """ Return the converter corresponding to a format name

    :param fmt: (string) the format name
    :returns: (module) the converter

    :raises ValueError: in case no exporter can be found

    """

    # Look dynamically which format is available
    for module_name in __all__:
        try:
            converter = importlib.import_module("." + module_name, "odemis.dataio")
        except (ValueError, TypeError, ImportError):
            logging.exception("Import of converter failed for fmt %s", fmt)
            continue  # module cannot be loaded

        if fmt == converter.FORMAT:
            return converter

    raise ValueError("No converter for format %s found" % fmt)


def find_fittest_exporter(filename, default=tiff):
    """
    Find the most fitting exporter according to a filename (actually, its extension)
    filename (string): (path +) filename with extension
    default (dataio. Module): default exporter to pick if no really fitting
      exporter is found
    returns (dataio. Module): the right exporter
    """
    # Find the extension of the file
    basename = os.path.basename(filename).lower()
    if basename == "":
        raise ValueError("Filename should have at least one letter: '%s'" % filename)

    # make sure we pick the format with the longest fitting extention
    best_len = 0
    best_fmt = default
    for module_name in __all__:
        try:
            exporter = importlib.import_module("." + module_name, "odemis.dataio")
        except Exception:  #pylint: disable=W0702
            continue # module cannot be loaded
        for e in exporter.EXTENSIONS:
            if filename.endswith(e) and len(e) > best_len:
                best_len = len(e)
                best_fmt = exporter

    if best_len > 0:
        logging.debug("Determined that '%s' corresponds to %s format",
                      basename, best_fmt.FORMAT)
    return best_fmt

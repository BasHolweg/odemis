# -*- coding: utf-8 -*-
'''
Created on 26 Mar 2012

@author: Éric Piel

Copyright © 2012-2013 Éric Piel, Delmic

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
import logging
import os
import subprocess

# Generic metadata about the package

def _get_version_git():
    """
    Get the version via git
    raises LookupError if no version info found
    """
    # change directory to root
    rootdir = os.path.join(os.path.dirname(__file__), "..", "..") # odemis/src/odemis/../..

    if not os.path.isdir(rootdir) or not os.path.isdir(os.path.join(rootdir, ".git")):
        raise LookupError("Not in a git directory")

    try:
        out = subprocess.check_output(args=["git", "describe", "--tags", "--dirty", "--always"],
                                      cwd=rootdir)
        return out.strip()
    except EnvironmentError:
        raise LookupError("Unable to run git")

def _get_version_setuptools():
    """
    Gets the version via setuptools/pkg_resources
    raises LookupError if no version info found
    """
    import pkg_resources
    try:
        return pkg_resources.get_distribution("odemis").version
    except pkg_resources.DistributionNotFound:
        raise LookupError("Not packaged via setuptools")

def _get_version():
    try:
        return _get_version_git()
    except LookupError:
        # fallback to setuptools (if it's not in git, it should be packaged)
        return _get_version_setuptools()
    except LookupError:
        logging.warning("Unable to find the actual version")
        return "Unknown"

def get_major_version():
    """ This function returns a short version string of the form "vX.X" """
    return _get_version().split("-")[0]

__version__ = _get_version()
__fullname__ = "Open Delmic Microscope Software"
__shortname__ = "Odemis"
__copyright__ = "Copyright © 2012-2015 Delmic"
__authors__ = ["Éric Piel", "Rinze de Laat", "Kimon Tsitsikas"]
__license__ = "GNU General Public License version 2"
__licensetxt__ = (
"""Odemis is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License version 2 as published by the Free Software
Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Odemis. If not, see http://www.gnu.org/licenses/.
""")

# vim:tabstop=4:shiftwidth=4:expandtab:spelllang=en_gb:spell:

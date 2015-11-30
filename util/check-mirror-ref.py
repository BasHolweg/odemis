#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Checks the reference switches of each axis of the mirror.
# Normally, they should be opposite from each other. If not, it's a sign of a
# faulty hardware.
'''
Created on November 2015

@author: Éric Piel

Copyright © 2015 Éric Piel, Delmic

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

import argparse
import logging
import Pyro4
from odemis import model
from odemis.driver import tmcm
import sys
import time


# Standard way to configure the TMCM board handling the Redux stage
REDUX_KWARGS = {
    "port": "/dev/ttyTMCM*",
    "address": 4,
    "axes": ["s", "l"],
    "ustepsize": [5.9e-9, 5.9e-9],  # m/µstep (doesn't really matter here)
    "refproc": "Standard",
    "refswitch": {"s": 0, "l": 0},
}

def check_ref(mirror):
    """
    Check reference switches of both axes of the mirror
    """
    # we rely on the fact that the axes S/L are 0/1 with reference switches
    # powered via output port 0

    # Turn on the reference switches power
    mirror.SetIO(2, 0, 1)

    # Check left/right switch values are opposite
    try:
        error = False
        for name, axis in {"s": 0, "l": 1}.items():
            vr = mirror.GetAxisParam(axis, 10)
            vl = mirror.GetAxisParam(axis, 11)
            if vl == vr:
                logging.error("Reference switches for axis %s both have value %d", name, vr)
                error = True
    finally:
        # Turn on the reference switches power
        mirror.SetIO(2, 0, 0)

    if error:
        raise ValueError("Reference switches didn't appear correct")

def check_via_backend():
    """
    Use the backend to access the mirror actuator component and check the
      reference switches.
    raise:
        CommunicationError if no backend present
        LookupError: backend is present but doesn't have mirror
        IOError: if move failed
    """
    mirror = model.getComponent(role="mirror")
    logging.debug("Using the backend to check the mirror reference switches")
    check_ref(mirror)


def check_direct():
    """
    Try to directly connect to the TMCM board and park the mirror
    """
    mirror = tmcm.TMCLController("Mirror stage", "mirror", **REDUX_KWARGS)
    logging.info("Connected to %s", mirror.hwVersion)

    try:
        check(mirror)
    finally:
        mirror.terminate()


def main(args):
    """
    Handles the command line arguments
    args is the list of arguments passed
    return (int): value to return to the OS as program exit code
    """

    # arguments handling
    parser = argparse.ArgumentParser(prog="check-mirror-ref",
                        description="Attempt to check the mirror reference switches of the SPARCv2")

    parser.add_argument("--log-level", dest="loglev", metavar="<level>", type=int,
                        default=1, help="set verbosity level (0-2, default = 1)")

    options = parser.parse_args(args[1:])

    # Set up logging before everything else
    if options.loglev < 0:
        logging.error("Log-level must be positive.")
        return 127
    loglev_names = (logging.WARNING, logging.INFO, logging.DEBUG)
    loglev = loglev_names[min(len(loglev_names) - 1, options.loglev)]
    logging.getLogger().setLevel(loglev)

    try:
        try:
            check_via_backend()
        except (Pyro4.errors.CommunicationError, IOError, LookupError):
            logging.info("Failed to access the backend, will try directly")
            check_direct()
    except KeyboardInterrupt:
        logging.info("Interrupted before the end of the execution")
        return 1
    except ValueError as exp:
        logging.error("%s", exp)
        return 127
    except IOError as exp:
        logging.error("%s", exp)
        return 129
    except Exception:
        logging.exception("Unexpected error while performing action.")
        return 130

    return 0


if __name__ == '__main__':
    ret = main(sys.argv)
    exit(ret)

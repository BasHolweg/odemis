#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 12 Mar 2012

@author: Éric Piel
Testing class for driver.andorcam3 .

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
import logging
from odemis.driver import andorcam3
import unittest
from unittest.case import skip

from cam_test_abs import VirtualTestCam, VirtualStaticTestCam, \
    VirtualTestSynchronized


logging.getLogger().setLevel(logging.DEBUG)

CLASS = andorcam3.AndorCam3
KWARGS = dict(name="camera", role="ccd", device=0, transpose=[2, -1],
              bitflow_install_dirs="/usr/share/bitflow/")

class StaticTestAndorCam3(VirtualStaticTestCam, unittest.TestCase):
    camera_type = CLASS
    camera_kwargs = KWARGS

# Inheritance order is important for setUp, tearDown
@skip("simple")
class TestAndorCam3(VirtualTestCam, unittest.TestCase):
    """
    Test directly the AndorCam3 class.
    """
    camera_type = CLASS
    camera_kwargs = KWARGS
    
#@skip("simple")
class TestSynchronized(VirtualTestSynchronized, unittest.TestCase):
    """
    Test the synchronizedOn(Event) interface, using the fake SEM
    """
    camera_type = CLASS
    camera_kwargs = KWARGS

if __name__ == '__main__':
    unittest.main()


#from odemis.driver import andorcam3
#import logging
#logging.getLogger().setLevel(logging.DEBUG)
#
#a = andorcam3.AndorCam3("test", "cam", 0)
#a.targetTemperature.value = -15
#a.fanSpeed.value = 0
#rr = a.readoutRate.value
#a.data.get()
#rt = a.GetFloat(u"ReadoutTime")
#res = a.resolution.value
#res[0] * res[1] / rr
#a.data.get()
#a.resolution.value = (128, 128)


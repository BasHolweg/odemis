# -*- coding: utf-8 -*-
'''
Created on 6 Nov 2013

@author: Éric Piel

Copyright © 2013 Éric Piel, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Odemis. If not, see http://www.gnu.org/licenses/.
'''
from __future__ import division

import logging
from odemis.driver import omicronxx
import os
import time
import unittest
from unittest.case import skip


logging.getLogger().setLevel(logging.DEBUG)

CLASS = omicronxx.MultixX

if os.name == "nt":
    PORT = "COM1"
else:
    PORTS = "/dev/ttyOXX*" #"/dev/tty*"

class TestActuator(unittest.TestCase):
    def setUp(self):
        self.dev = CLASS("test", "light", PORTS)

    def tearDown(self):
        self.dev.terminate()

    def test_simple(self):
        # should start off
        self.assertEqual(self.dev.power.value, 0)

        # turn on first source to 10%
        self.dev.power.value = self.dev.power.range[1]
        em = self.dev.emissions.value
        em[0] = 0.1
        self.dev.emissions.value = em
        self.assertGreater(self.dev.emissions.value[0], 0)

    def test_cycle(self):
        """
        Test each emission source for 2 seconds at maximum intensity and then 1s
        at 10%.
        """
        em = self.dev.emissions.value
        em = [0 for v in em]
        self.dev.power.value = self.dev.power.range[1]

        # can fully checked only by looking what the hardware is doing
        print "Starting emission source cycle..."
        for i in range(len(em)):
            print "Turning on wavelength %g" % self.dev.spectra.value[i][2]
            em[i] = 0.1
            self.dev.emissions.value = em
            time.sleep(5)
            self.assertEqual(self.dev.emissions.value, em)
            em[i] = 0
            self.dev.emissions.value = em
            self.assertEqual(self.dev.emissions.value, em)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 22 Nov 2016

@author: Éric Piel

Copyright © 2016 Éric Piel, Delmic

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

import logging
from odemis.driver import pwrcomedi, omicronxx, emitter, rigol, simulated
import time
import unittest

logging.getLogger().setLevel(logging.DEBUG)

DEPENDENCY1_CLASS = pwrcomedi.Light
DEPENDENCY1_KWARGS = {"name": "test1", "role": None,
                 "device": "/dev/comedi0", # Simulator, if comedi_test is loaded
          "channels": [0, 2],
          "spectra": [(615.e-9, 625.e-9, 633.e-9, 640.e-9, 650.e-9),
                      (525.e-9, 540.e-9, 550.e-9, 555.e-9, 560.e-9)],
          # Omicron has power max = 1.4W => need to have at least 30% of that on each source
          "pwr_curve": [{-3: 0, # V -> W
                         3: 1,
                        },
                        {# Missing 0W => 0V -> 0W
                         0.1: 0.1,
                         0.3: 0.2,
                         0.5: 0.4,
                         0.7: 0.8,
                         1: 1.2,
                        }
                        ]
         }
DEPENDENCY2_CLASS = omicronxx.HubxX
DEPENDENCY2_KWARGS = {"name": "test2", "role": None, "port": "/dev/fakehub"}
CONFIG_DG1000Z = {"name": "Rigol Wave Gen", "role": "pc-emitter",
                "host": "fake",
                "port": 5555, "channel": 1,
                "limits": (-10.0, 10.0)
}
KWARGS = {"name": "test", "role": "light"}


class TestMultiplexLight(unittest.TestCase):

    def setUp(self):
        self.dependency1 = DEPENDENCY1_CLASS(**DEPENDENCY1_KWARGS)
        self.dependency2 = DEPENDENCY2_CLASS(**DEPENDENCY2_KWARGS)
        self.dev = emitter.MultiplexLight("test", "light",
                                          dependencies={"c1": self.dependency1, "c2": self.dependency2})

    def tearDown(self):
        self.dev.terminate()
        self.dependency1.terminate()
        self.dependency2.terminate()

    def test_simple(self):
        # should start off
        self.assertEqual(self.dev.power.value, 0)

        self.assertEqual(len(self.dev.emissions.value), len(self.dev.spectra.value))

        # turn on first source to 50%
        self.dev.power.value = self.dev.power.range[1]
        em = self.dev.emissions.value
        em[0] = 0.5
        self.dev.emissions.value = em
        self.assertGreater(self.dev.emissions.value[0], 0)

        # turn on second source to 90%
        self.dev.power.value = self.dev.power.range[1]
        em[0:2] = [0, 0.9]
        self.dev.emissions.value = em
        self.assertGreater(self.dev.emissions.value[1], 0)

    def test_multi(self):
        """
        simultaneous source activation
        """
        self.dev.power.value = self.dev.power.range[1]
        em = [1] * len(self.dev.emissions.value)
        self.dev.emissions.value = em
        # They should all be on
        self.assertTrue(all(e > 0 for e in self.dev.emissions.value))

        # Not all should be at the max, due to clamping
        self.assertTrue(any(e < 1 for e in self.dev.emissions.value))

    def test_cycle(self):
        """
        Test each emission source for 2 seconds at maximum intensity and then 1s
        at 30%.
        """
        em = [0] * len(self.dev.emissions.value)
        self.dev.power.value = self.dev.power.range[1]

        # can fully checked only by looking what the hardware is doing
        logging.info("Starting emission source cycle...")
        for i in range(len(em)):
            logging.info("Turning on wavelength %g", self.dev.spectra.value[i][2])
            em[i] = 1
            self.dev.emissions.value = em
            time.sleep(1)
            self.assertGreater(self.dev.emissions.value, 0)  # Can't check for equality due to clamping

            em[i] = 0.3
            self.dev.emissions.value = em
            time.sleep(1)
            self.assertEqual(self.dev.emissions.value, em)
#             # value so small that it's == 0 for the hardware
#             self.dev.emissions.value[i] = 1e-8
#             em[i] = 0
#             self.assertEqual(self.dev.emissions.value, em)


class TestExtendedLight(unittest.TestCase):
    """
    Tests for extended light
    """

    def setUp(self):
        self.wg = rigol.WaveGenerator(**CONFIG_DG1000Z)    # specify IP of actual device
        self.light = simulated.Light("test", "light")
        CONFIG_EX_LIGHT = {"name": "Test Extended Light", "role": None,
                           "dependencies": {"light": self.light, "clock": self.wg }
                           }
        self.ex_light = emitter.ExtendedLight(**CONFIG_EX_LIGHT)

    def tearDown(self):
        self.wg.terminate() # free up gsocket.
        time.sleep(1.0) # give some time to make sure socket is released.

    def test_power(self):
        '''
        Test 1: If there are emissions, and power > 0, the wave generator power
        should be active (1) and the light power should be the same as ex_light power
        '''
        self.ex_light.power.value = 5
        em = self.ex_light.emissions.value
        em[0] = 0.5
        self.ex_light.emissions.value = em
        self.assertEqual(self.ex_light.power.value, 5)
        self.assertEqual(self.wg.power.value, 1)
        self.assertEqual(self.light.power.value, self.ex_light.power.value)
        '''
        Test 2: If there are no emissions, and power > 0, the wave generator power
        should be off (0) and the light power should be the same as ex_light power
        '''
        self.ex_light.power.value = 5
        em = [0] * len(self.ex_light.emissions.value)
        self.ex_light.emissions.value = em
        self.assertEqual(self.ex_light.power.value, 5)
        self.assertEqual(self.wg.power.value, 0)
        self.assertEqual(self.light.power.value, self.ex_light.power.value)
        '''
        Test 3: If there are emissions, but power = 0, the wave generator power
        should be off (0) and the light power should be the same as ex_light power
        '''
        self.ex_light.power.value = 0
        em = self.ex_light.emissions.value
        em[0] = 0.5
        self.ex_light.emissions.value = em
        self.assertEqual(self.ex_light.power.value, 0)
        self.assertEqual(self.wg.power.value, 0)
        self.assertEqual(self.light.power.value, self.ex_light.power.value)

    def test_period(self):
        self.ex_light.power.value = 0
        self.assertEqual(self.ex_light.power.value, 0)
        for i in range(1000, 10000, 1000):    # specify range of frequencies to increment
            self.ex_light.period.value = 1 / i
            self.assertEqual(self.ex_light.period.value, 1 / i)
            self.ex_light.power.value = 5
            time.sleep(0.1)
            self.ex_light.power.value = 0

if __name__ == "__main__":
    unittest.main()

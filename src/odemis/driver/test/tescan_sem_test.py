#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 12 May 2014

Copyright © 2014 Kimon Tsitsikas, Delmic

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

import Pyro4
import copy
import logging
from odemis import model
from odemis.driver import tescan_sem
import os
import pickle
import threading
import time
import unittest
from odemis.dataio import hdf5
from unittest.case import skip


logging.getLogger().setLevel(logging.DEBUG)

# arguments used for the creation of basic components
CONFIG_SED = {"name": "sed", "role": "sed", "channel": 0}
CONFIG_BSD = {"name": "bsd", "role": "bsd"}
CONFIG_STG = {"name": "stg", "role": "stage"}
CONFIG_CM = {"name": "camera", "role": "camera"}
CONFIG_FOCUS = {"name": "focus", "role": "focus", "axes": ["z"]}
CONFIG_PRESSURE = {"name": "pressure", "role": "pressure"}
CONFIG_SCANNER = {"name": "scanner", "role": "ebeam",
                  "fov_range": [196.e-9, 25586.e-6]}
CONFIG_SEM = {"name": "sem", "role": "sem",
              "children": {"detector": CONFIG_SED, "scanner": CONFIG_SCANNER,
                           "stage": CONFIG_STG, "focus": CONFIG_FOCUS,
                           "camera": CONFIG_CM, "pressure": CONFIG_PRESSURE},
              "host": "192.168.92.91"
              }

class TestSEMStatic(unittest.TestCase):
    """
    Tests which don't need a SEM component ready
    """
    def test_creation(self):
        """
        Doesn't even try to acquire an image, just create and delete components
        """
        sem = tescan_sem.TescanSEM(**CONFIG_SEM)
        self.assertEqual(len(sem.children), 6)
        
        for child in sem.children:
            if child.name == CONFIG_SED["name"]:
                sed = child
            elif child.name == CONFIG_SCANNER["name"]:
                scanner = child
            elif child.name == CONFIG_STG["name"]:
                stage = child
            elif child.name == CONFIG_FOCUS["name"]:
                focus = child
            elif child.name == CONFIG_CM["name"]:
                camera = child
            elif child.name == CONFIG_PRESSURE["name"]:
                pressure = child
        
        self.assertEqual(len(scanner.resolution.value), 2)
        self.assertIsInstance(sed.data, model.DataFlow)
        
        self.assertTrue(sem.selfTest(), "SEM self test failed.")
        sem.terminate()
    
    def test_error(self):
        wrong_config = copy.deepcopy(CONFIG_SEM)
        wrong_config["children"]["scanner"]["channels"] = [1, 1]
        self.assertRaises(Exception, tescan_sem.TescanSEM, **wrong_config)
    
    def test_pickle(self):
        try:
            os.remove("test")
        except OSError:
            pass
        daemon = Pyro4.Daemon(unixsocket="test")
        
        sem = tescan_sem.TescanSEM(daemon=daemon, **CONFIG_SEM)
                
        dump = pickle.dumps(sem, pickle.HIGHEST_PROTOCOL)
#        print "dump size is", len(dump)
        sem_unpickled = pickle.loads(dump)
        self.assertEqual(len(sem_unpickled.children), 6)
        sem.terminate()
    
class TestSEM(unittest.TestCase):
    """
    Tests which can share one SEM device
    """
    @classmethod
    def setUpClass(cls):
        cls.sem = tescan_sem.TescanSEM(**CONFIG_SEM)
        
        for child in cls.sem.children:
            if child.name == CONFIG_SED["name"]:
                cls.sed = child
            elif child.name == CONFIG_SCANNER["name"]:
                cls.scanner = child
            elif child.name == CONFIG_STG["name"]:
                cls.stage = child
            elif child.name == CONFIG_FOCUS["name"]:
                cls.focus = child
            elif child.name == CONFIG_CM["name"]:
                cls.camera = child
            elif child.name == CONFIG_PRESSURE["name"]:
                cls.pressure = child

    @classmethod
    def tearUpClass(cls):
        cls.sem.terminate()
        time.sleep(3)

    def setUp(self):
        # reset resolution and dwellTime
        self.scanner.scale.value = (1, 1)
        self.scanner.resolution.value = (512, 256)
        self.size = self.scanner.resolution.value
        self.scanner.dwellTime.value = self.scanner.dwellTime.range[0]
        self.acq_dates = (set(), set()) # 2 sets of dates, one for each receiver
        self.acq_done = threading.Event()

    def tearUp(self):
#        print gc.get_referrers(self.camera)
#        gc.collect()
        pass

    def assertTupleAlmostEqual(self, first, second, places=None, msg=None, delta=None):
        """
        check two tuples are almost equal (value by value)
        """
        for f, s in zip(first, second):
            self.assertAlmostEqual(f, s, places=places, msg=msg, delta=delta)

    def compute_expected_duration(self):
        dwell = self.scanner.dwellTime.value
        settle = 5.e-6
        size = self.scanner.resolution.value
        return size[0] * size[1] * dwell + size[1] * settle

    def test_acquire(self):
        self.scanner.dwellTime.value = 10e-6  # s
        expected_duration = self.compute_expected_duration()

        start = time.time()
        im = self.sed.data.get()
        hdf5.export("test.h5", model.DataArray(im))
        duration = time.time() - start

        self.assertEqual(im.shape, self.size[::-1])
        self.assertGreaterEqual(duration, expected_duration, "Error execution took %f s, less than exposure time %d." % (duration, expected_duration))
        self.assertIn(model.MD_DWELL_TIME, im.metadata)

    def test_roi(self):
        """
        check that .translation and .scale work
        """
        # First, test simple behaviour on the VA
        # max resolution
        max_res = self.scanner.resolution.range[1]
        self.scanner.scale.value = (1, 1)
        self.scanner.resolution.value = max_res
        self.scanner.translation.value = (-1, 1)  # will be set back to 0,0 as it cannot move
        self.assertEqual(self.scanner.translation.value, (0, 0))

        # scale up
        self.scanner.scale.value = (16, 16)
        exp_res = (max_res[0] // 16, max_res[1] // 16)
        self.assertTupleAlmostEqual(self.scanner.resolution.value, exp_res)
        self.scanner.translation.value = (-1, 1)
        self.assertEqual(self.scanner.translation.value, (0, 0))

        # shift
        exp_res = (max_res[0] // 32, max_res[1] // 32)
        self.scanner.resolution.value = exp_res
        self.scanner.translation.value = (-1, 1)
        self.assertTupleAlmostEqual(self.scanner.resolution.value, exp_res)
        self.assertEqual(self.scanner.translation.value, (-1, 1))

        # change scale to some float
        self.scanner.resolution.value = (max_res[0] // 16, max_res[1] // 16)
        self.scanner.scale.value = (1.5, 2.3)
        exp_res = (max_res[0] // 1.5, max_res[1] // 2.3)
        self.assertTupleAlmostEqual(self.scanner.resolution.value, exp_res)
        self.assertEqual(self.scanner.translation.value, (0, 0))

        self.scanner.scale.value = (1, 1)
        self.assertTupleAlmostEqual(self.scanner.resolution.value, max_res, delta=1.1)
        self.assertEqual(self.scanner.translation.value, (0, 0))

        # Then, check metadata fits with the expectations
        center = (1e3, -2e3)  # m
        # simulate the information on the position (normally from the mdupdater)
        self.scanner.updateMetadata({model.MD_POS: center})

        self.scanner.resolution.value = max_res
        self.scanner.scale.value = (16, 16)
        self.scanner.dwellTime.value = self.scanner.dwellTime.range[0]

        # normal acquisition
        im = self.sed.data.get()
        self.assertEqual(im.shape, self.scanner.resolution.value[-1::-1])
        self.assertTupleAlmostEqual(im.metadata[model.MD_POS], center)

        # shift a bit
        # reduce the size of the image so that we can have translation
        self.scanner.resolution.value = (max_res[0] // 32, max_res[1] // 32)
        self.scanner.translation.value = (-1.26, 10)  # px
        pxs = self.scanner.pixelSize.value
        exp_pos = (center[0] + (-1.26 * pxs[0]),
                   center[1] - (10 * pxs[1]))  # because translation Y is opposite from physical one
        im = self.sed.data.get()
        self.assertEqual(im.shape, self.scanner.resolution.value[-1::-1])
        self.assertTupleAlmostEqual(im.metadata[model.MD_POS], exp_pos)

        # only one point
        self.scanner.resolution.value = (1, 1)
        im = self.sed.data.get()
        self.assertEqual(im.shape, self.scanner.resolution.value[-1::-1])
        self.assertTupleAlmostEqual(im.metadata[model.MD_POS], exp_pos)

    def test_acquire_high_osr(self):
        """
        small resolution, but large osr, to force acquisition not by whole array
        """
        self.scanner.resolution.value = (256, 200)
        self.size = self.scanner.resolution.value
        self.scanner.dwellTime.value = self.scanner.dwellTime.range[0] * 1000
        expected_duration = self.compute_expected_duration()  # about 1 min

        start = time.time()
        im = self.sed.data.get()
        duration = time.time() - start

        self.assertEqual(im.shape, self.size[-1:-3:-1])
        self.assertGreaterEqual(duration, expected_duration, "Error execution took %f s, less than exposure time %d." % (duration, expected_duration))
        self.assertIn(model.MD_DWELL_TIME, im.metadata)

    def test_acquire_flow(self):
        expected_duration = self.compute_expected_duration()

        number = 5
        self.left = number
        self.sed.data.subscribe(self.receive_image)

        self.acq_done.wait(number * (2 + expected_duration * 1.1))  # 2s per image should be more than enough in any case

        self.assertEqual(self.left, 0)

    def test_df_fast_sub_unsub(self):
        """
        Test the dataflow on a very fast cycle subscribing/unsubscribing
        SEMComedi had a bug causing the threads not to start again
        """
        self.scanner.dwellTime.value = self.scanner.dwellTime.range[0]
        number = 10
        expected_duration = self.compute_expected_duration()

        self.left = 10000  # don't unsubscribe automatically

        for i in range(number):
            self.sed.data.subscribe(self.receive_image)
            time.sleep(i)
            self.sed.data.unsubscribe(self.receive_image)

        # now this one should work
        self.sed.data.subscribe(self.receive_image)
        time.sleep(expected_duration * 5)  # make sure we received at least one image
        self.sed.data.unsubscribe(self.receive_image)

        self.assertLessEqual(self.left, 10000 - 1)

    def test_df_alternate_sub_unsub(self):
        """
        Test the dataflow on a quick cycle subscribing/unsubscribing
        Andorcam3 had a real bug causing deadlock in this scenario
        """
        self.scanner.dwellTime.value = 10e-6
        number = 5
        expected_duration = self.compute_expected_duration()

        self.left = 10000 + number  # don't unsubscribe automatically

        for i in range(number):
            self.sed.data.subscribe(self.receive_image)
            time.sleep(expected_duration * 6)  # make sure we received at least one image
            self.sed.data.unsubscribe(self.receive_image)

        # if it has acquired a least 5 pictures we are already happy
        self.assertLessEqual(self.left, 10000)

    def onEvent(self):
        self.events += 1

    def receive_image(self, dataflow, image):
        """
        callback for df of test_acquire_flow()
        """
        self.assertEqual(image.shape, self.size[-1:-3:-1])
        self.assertIn(model.MD_DWELL_TIME, image.metadata)
        self.acq_dates[0].add(image.metadata[model.MD_ACQ_DATE])
#         print "Received an image"
        self.left -= 1
        if self.left <= 0:
            dataflow.unsubscribe(self.receive_image)
            self.acq_done.set()

if __name__ == "__main__":
    unittest.main()

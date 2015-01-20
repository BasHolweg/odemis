# -*- coding: utf-8 -*-
'''
Created on 10 Jan 2014

@author: Kimon Tsitsikas

Copyright © 2014 Kimon Tsitsikas, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Odemis. If not, see http://www.gnu.org/licenses/.
'''
from __future__ import division

import numpy
from odemis import model
from odemis.dataio import hdf5
from odemis.util import polar
import unittest


class TestPolarConversion(unittest.TestCase):
    """
    Test AngleResolved2Polar
    """
    def setUp(self):
        data = hdf5.read_data("ar-example-input.h5")
        mag = 0.4917
        spxs = (13e-6, 13e-6)
        binning = (4, 4)
        # data[0].metadata[model.MD_BASELINE] = 820
        data[0].metadata[model.MD_BINNING] = binning
        data[0].metadata[model.MD_SENSOR_PIXEL_SIZE] = spxs
        data[0].metadata[model.MD_LENS_MAG] = mag
        data[0].metadata[model.MD_AR_POLE] = (141, 255 - 139.449038462)
        mag = data[0].metadata[model.MD_LENS_MAG]
        pxs = (spxs[0] * binning[0] / mag,
               spxs[1] * binning[1] / mag)
        data[0].metadata[model.MD_PIXEL_SIZE] = pxs
        self.data = data

        white_data_512 = model.DataArray(numpy.empty((512, 512), dtype="uint16"))
        white_data_512[...] = 255
        white_mag_512 = 0.4917
        white_spxs_512 = (13e-6, 13e-6)
        white_binning_512 = (2, 2)
        white_data_512.metadata[model.MD_AR_POLE] = (283, 259)
        white_pxs_512 = (white_spxs_512[0] * white_binning_512[0] / white_mag_512,
               white_spxs_512[1] * white_binning_512[1] / white_mag_512)
        white_data_512.metadata[model.MD_PIXEL_SIZE] = white_pxs_512
        self.white_data_512 = white_data_512

        white_data_1024 = model.DataArray(numpy.empty((1024, 1024), dtype="uint16"))
        white_data_1024[...] = 255
        white_mag_1024 = 0.4917
        white_spxs_1024 = (13e-6, 13e-6)
        white_binning_1024 = (2, 2)
        white_data_1024.metadata[model.MD_AR_POLE] = (283, 259)
        white_pxs_1024 = (white_spxs_1024[0] * white_binning_1024[0] / white_mag_1024,
               white_spxs_1024[1] * white_binning_1024[1] / white_mag_1024)
        white_data_1024.metadata[model.MD_PIXEL_SIZE] = white_pxs_1024
        self.white_data_1024 = white_data_1024

        white_data_2500 = model.DataArray(numpy.empty((2560, 2160), dtype="uint16"))
        white_data_2500[...] = 255
        white_mag_2500 = 0.4917
        white_spxs_2500 = (13e-6, 13e-6)
        white_binning_2500 = (2, 2)
        white_data_2500.metadata[model.MD_AR_POLE] = (283, 259)
        # These values makes the computation much harder:
#         white_mag_2500 = 0.53
#         white_spxs_2500 = (6.5e-6, 6.5e-6)
#         white_binning_2500 = (1, 1)
#         white_data_2500.metadata[model.MD_AR_POLE] = (1480, 1129)
        white_pxs_2500 = (white_spxs_2500[0] * white_binning_2500[0] / white_mag_2500,
               white_spxs_2500[1] * white_binning_2500[1] / white_mag_2500)
        white_data_2500.metadata[model.MD_PIXEL_SIZE] = white_pxs_2500
        self.white_data_2500 = white_data_2500

    def test_precomputed(self):
        data = self.data
        C, T, Z, Y, X = data[0].shape
        data[0].shape = Y, X
        result = polar.AngleResolved2Polar(data[0], 201)

        desired_output = hdf5.read_data("desired201x201image.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1e-04)

    def test_uint16_input(self):
        """
        Tests for input of DataArray with uint16 ndarray.
        """
        data = self.data
        C, T, Z, Y, X = data[0].shape
        data[0].shape = Y, X
        result = polar.AngleResolved2Polar(data[0], 201)

        desired_output = hdf5.read_data("desired201x201image.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1e-04)

    def test_int8_input(self):
        """
        Tests for input of DataArray with int8 ndarray.
        """
        data = self.data
        # scipy.misc.bytescale(data)
        data[0] = data[0].astype(numpy.int64)
        data[0] = numpy.right_shift(data[0], 8)
        data[0] = data[0].astype(numpy.int8)
        C, T, Z, Y, X = data[0].shape
        data[0].shape = Y, X
        result = polar.AngleResolved2Polar(data[0], 201)

        desired_output = polar.AngleResolved2Polar(data[0].astype(float), 201)

        numpy.testing.assert_allclose(result, desired_output, rtol=1e-04)

    def test_float_input(self):
        """
        Tests for input of DataArray with float ndarray.
        """
        data = self.data
        data[0] = data[0].astype(numpy.float)
        C, T, Z, Y, X = data[0].shape
        data[0].shape = Y, X
        result = polar.AngleResolved2Polar(data[0], 201)

        desired_output = hdf5.read_data("desired201x201image.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1e-04)

    def test_100x100(self):
        data = self.data
        C, T, Z, Y, X = data[0].shape
        data[0].shape = Y, X

        result = polar.AngleResolved2Polar(data[0], 101)

        desired_output = hdf5.read_data("desired100x100image.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1e-04)

    def test_1000x1000(self):
        data = self.data
        data[0] = data[0].astype(numpy.int64)
        data[0] = numpy.right_shift(data[0], 8)
        data[0] = data[0].astype(numpy.int8)
        C, T, Z, Y, X = data[0].shape
        data[0].shape = Y, X

        result = polar.AngleResolved2Polar(data[0], 1001)

        desired_output = hdf5.read_data("desired1000x1000image.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1)

    def test_2000x2000(self):
        data = self.data
        C, T, Z, Y, X = data[0].shape
        data[0].shape = Y, X

        result = polar.AngleResolved2Polar(data[0], 2001)

        desired_output = hdf5.read_data("desired2000x2000image.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1e-04)

    def test_512x512(self):
        """
        Test for 512x512 white image input
        """
        white_data_512 = self.white_data_512
        Y, X = white_data_512.shape
        result = polar.AngleResolved2Polar(white_data_512, 201)

        desired_output = hdf5.read_data("desired_white_512.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1e-04)

    def test_1024x1024(self):
        """
        Test for 1024x1024 white image input
        """
        white_data_1024 = self.white_data_1024
        Y, X = white_data_1024.shape
        result = polar.AngleResolved2Polar(white_data_1024, 201)

        desired_output = hdf5.read_data("desired_white_1024.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1e-04)

    def test_2560x2160(self):
        """
        Test for 2560x2160 white image input
        """
        white_data_2500 = self.white_data_2500
        Y, X = white_data_2500.shape
        result = polar.AngleResolved2Polar(white_data_2500, 2000, dtype=numpy.float16)

        desired_output = hdf5.read_data("desired_white_2500.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        # FIXME: doesn't seem to pass on 64 bits ?! floating point computation differences?

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1e-04)

    def test_background_substraction_precomputed(self):
        """
        Test clean up before polar conversion
        """
        data = self.data
        C, T, Z, Y, X = data[0].shape
        data[0].shape = Y, X
        clean_data = polar.ARBackgroundSubtract(data[0])
        result = polar.AngleResolved2Polar(clean_data, 201)

        desired_output = hdf5.read_data("substracted_background_image.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1e-04)

    def test_background_substraction_uint16_input(self):
        """
        Tests for input of DataArray with uint16 ndarray.
        """
        data = self.data
        data[0] = data[0].astype(numpy.uint16)
        C, T, Z, Y, X = data[0].shape
        data[0].shape = Y, X
        clean_data = polar.ARBackgroundSubtract(data[0])
        result = polar.AngleResolved2Polar(clean_data, 201)

        desired_output = hdf5.read_data("substracted_background_image.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1e-04)

    def test_background_substraction_int8_input(self):
        """
        Tests for input of DataArray with int8 ndarray.
        """
        data = self.data
        # scipy.misc.bytescale(data)
        data[0] = data[0].astype(numpy.int64)
        data[0] = numpy.right_shift(data[0], 8)
        data[0] = data[0].astype(numpy.int8)
        C, T, Z, Y, X = data[0].shape
        data[0].shape = Y, X
        clean_data = polar.ARBackgroundSubtract(data[0])
        result = polar.AngleResolved2Polar(clean_data, 201)

        desired_output = polar.AngleResolved2Polar(data[0].astype(float), 201)

        numpy.testing.assert_allclose(result, desired_output, rtol=1)

    def test_background_substraction_float_input(self):
        """
        Tests for input of DataArray with float ndarray.
        """
        data = self.data
        data[0] = data[0].astype(numpy.float)
        C, T, Z, Y, X = data[0].shape
        data[0].shape = Y, X
        clean_data = polar.ARBackgroundSubtract(data[0])
        result = polar.AngleResolved2Polar(clean_data, 201)

        desired_output = hdf5.read_data("substracted_background_image.h5")
        C, T, Z, Y, X = desired_output[0].shape
        desired_output[0].shape = Y, X

        numpy.testing.assert_allclose(result, desired_output[0], rtol=1e-04)


if __name__ == "__main__":
#     import sys;sys.argv = ['', 'TestPolarConversionOutput.test_2000x2000']
    unittest.main()
#    suite = unittest.TestLoader().loadTestsFromTestCase(TestPolarConversionOutput)
#    unittest.TextTestRunner(verbosity=2).run(suite)


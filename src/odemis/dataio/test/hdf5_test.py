# -*- coding: utf-8 -*-
'''
Created on 14 Jan 2013

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
from numpy.polynomial import polynomial
from odemis import model
from odemis.dataio import hdf5
from unittest.case import skip
import h5py
import logging
import numpy
import os
import time
import unittest

logging.getLogger().setLevel(logging.DEBUG)

FILENAME = "test" + hdf5.EXTENSIONS[0]
class TestHDF5IO(unittest.TestCase):

    def tearDown(self):
        # clean up
        try:
            os.remove(FILENAME)
        except Exception:
            pass

    # Be careful: numpy's notation means that the pixel coordinates are Y,X,C
    def testExportOnePage(self):
        # create a simple greyscale image
        size = (256, 512) # (width, height)
        dtype = numpy.dtype("uint16")
        data = model.DataArray(numpy.zeros(size[-1:-3:-1], dtype))
        white = (12, 52) # non symmetric position
        # less that 2**15 so that we don't have problem with PIL.getpixel() always returning an signed int
        data[white[-1:-3:-1]] = 124

        # export
        hdf5.export(FILENAME, data)

        # check it's here
        st = os.stat(FILENAME) # this test also that the file is created
        self.assertGreater(st.st_size, 0)

        f = h5py.File(FILENAME, "r")
        # need to transform to a full numpy.array just to remove the dimensions
        im = numpy.array(f["Acquisition0/ImageData/Image"])
        im.shape = im.shape[3:5]
        self.assertEqual(im.shape, data.shape)
        self.assertEqual(im[white[-1:-3:-1]], data[white[-1:-3:-1]])

#    @skip("Doesn't work")
    def testExportMultiPage(self):
        # create a simple greyscale image
        size = (512, 256)
        white = (12, 52) # non symmetric position
        dtype = numpy.dtype("uint8")
        ldata = []
        num = 2
        for i in range(num):
            a = model.DataArray(numpy.zeros(size[-1:-3:-1], dtype))
            a[white[-1:-3:-1]] = 124
            ldata.append(a)

        # export
        hdf5.export(FILENAME, ldata)

        # check it's here
        st = os.stat(FILENAME) # this test also that the file is created
        self.assertGreater(st.st_size, 0)
        f = h5py.File(FILENAME, "r")

        # check the number of channels
        im = numpy.array(f["Acquisition0/ImageData/Image"])
        for i in range(num):
            subim = im[i, 0, 0] # just one channel
            self.assertEqual(subim.shape, size[-1:-3:-1])
            self.assertEqual(subim[white[-1:-3:-1]], 124)

        os.remove(FILENAME)

#    @skip("Doesn't work")
    def testExportThumbnail(self):
        # create 2 simple greyscale images
        size = (512, 256)
        dtype = numpy.dtype("uint16")
        ldata = []
        num = 2
        for i in range(num):
            ldata.append(model.DataArray(numpy.zeros(size[-1:-3:-1], dtype)))

        # thumbnail : small RGB completely red
        tshape = (size[1] // 8, size[0] // 8, 3)
        tdtype = numpy.uint8
        thumbnail = model.DataArray(numpy.zeros(tshape, tdtype))
        thumbnail[:, :, 0] += 255 # red
        blue = (12, 22) # non symmetric position
        thumbnail[blue[-1:-3:-1]] = [0, 0, 255]

        # export
        hdf5.export(FILENAME, ldata, thumbnail)

        # check it's here
        st = os.stat(FILENAME) # this test also that the file is created
        self.assertGreater(st.st_size, 0)
        f = h5py.File(FILENAME, "r")

        # look for the thumbnail
        im = f["Preview/Image"]
        self.assertEqual(im.shape, tshape)
        self.assertEqual(im[0, 0].tolist(), [255, 0, 0])
        self.assertEqual(im[blue[-1:-3:-1]].tolist(), [0, 0, 255])

        # check the number of channels
        im = numpy.array(f["Acquisition0/ImageData/Image"])
        for i in range(num):
            subim = im[i, 0, 0] # just one channel
            self.assertEqual(subim.shape, size[-1::-1])


    def testExportCube(self):
        """
        Check it's possible to export a 3D data (typically: 2D area with full
         spectrum for each point)
        """
        dtype = numpy.dtype("uint16")
        size3d = (512, 256, 220) # X, Y, C
        size = (512, 256)
        metadata3d = {model.MD_SW_VERSION: "1.0-test",
                    model.MD_HW_NAME: "fake spec",
                    model.MD_DESCRIPTION: "test3d",
                    model.MD_ACQ_DATE: time.time(),
                    model.MD_BPP: 12,
                    model.MD_BINNING: (1, 1), # px, px
                    model.MD_PIXEL_SIZE: (1e-6, 2e-5), # m/px
                    model.MD_WL_POLYNOMIAL: [500e-9, 1e-9], # m, m/px: wl polynomial
                    model.MD_POS: (1e-3, -30e-3), # m
                    model.MD_EXP_TIME: 1.2, #s
                    model.MD_IN_WL: (500e-9, 520e-9), #m
                    }
        metadata = {model.MD_SW_VERSION: "1.0-test",
                    model.MD_HW_NAME: u"", # check empty unicode strings
                    model.MD_DESCRIPTION: u"tÉst",
                    model.MD_ACQ_DATE: time.time(),
                    model.MD_BPP: 12,
                    model.MD_BINNING: (1, 2), # px, px
                    model.MD_PIXEL_SIZE: (1e-6, 2e-5), # m/px
                    model.MD_POS: (1e-3, -30e-3), # m
                    model.MD_DWELL_TIME: 1.2, #s
                    model.MD_IN_WL: (500e-9, 520e-9), #m
                    }
        ldata = []
        # 3D data generation (+ metadata): gradient along the wavelength
        data3d = numpy.empty(size3d[-1::-1], dtype=dtype)
        end = 2 ** metadata3d[model.MD_BPP]
        step = end // size3d[2]
        lin = numpy.arange(0, end, step, dtype=dtype)[:size3d[2]]
        lin.shape = (size3d[2], 1, 1) # to be able to copy it on the first dim
        data3d[:] = lin
        # introduce Time and Z dimension to state the 3rd dim is channel
        data3d = data3d[:, numpy.newaxis, numpy.newaxis, :, :]
        ldata.append(model.DataArray(data3d, metadata3d))

        # an additional 2D data, for the sake of it
        ldata.append(model.DataArray(numpy.zeros(size[-1::-1], dtype), metadata))

        # export
        hdf5.export(FILENAME, ldata)

        # check it's here
        st = os.stat(FILENAME) # this test also that the file is created
        self.assertGreater(st.st_size, 0)
        f = h5py.File(FILENAME, "r")

        # check the 3D data
        im = f["Acquisition0/ImageData/Image"]
        self.assertEqual(im[1, 0, 0, 1, 1], step)
        self.assertEqual(im.shape, data3d.shape)
        self.assertEqual(im.attrs["IMAGE_SUBCLASS"], "IMAGE_GRAYSCALE")

        # check basic metadata
        self.assertEqual(im.dims[4].label, "X")
        self.assertEqual(im.dims[0].label, "C")
        # wl polynomial is linear
        cres = im.dims[0][0][()] # first dimension (C), first scale, first and only value
        self.assertAlmostEqual(metadata3d[model.MD_WL_POLYNOMIAL][1], cres)
        coff = f["Acquisition0/ImageData/COffset"][()]
        self.assertAlmostEqual(metadata3d[model.MD_WL_POLYNOMIAL][0], coff)

        # check the 2D data
        im = numpy.array(f["Acquisition1/ImageData/Image"])
        subim = im[0, 0, 0] # just one channel
        self.assertEqual(subim.shape, size[-1::-1])


    def testMetadata(self):
        """
        checks that the metadata is saved with every picture
        """
        size = (512, 256, 1)
        dtype = numpy.dtype("uint16")
        metadata = {model.MD_SW_VERSION: "1.0-test",
                    model.MD_HW_NAME: "fake hw",
                    model.MD_DESCRIPTION: "test",
                    model.MD_ACQ_DATE: time.time(),
                    model.MD_BPP: 12,
                    model.MD_BINNING: (1, 2), # px, px
                    model.MD_PIXEL_SIZE: (1e-6, 2e-5), # m/px
                    model.MD_POS: (1e-3, -30e-3), # m
                    model.MD_EXP_TIME: 1.2, #s
                    model.MD_IN_WL: (500e-9, 520e-9), #m
                    }

        data = model.DataArray(numpy.zeros((size[1], size[0]), dtype), metadata=metadata)

        # thumbnail : small RGB completely red
        tshape = (size[1] // 8, size[0] // 8, 3)
        tdtype = numpy.uint8
        thumbnail = model.DataArray(numpy.zeros(tshape, tdtype))
        thumbnail[:, :, 0] += 255 # red
        blue = (12, 22) # non symmetric position
        thumbnail[blue[-1:-3:-1]] = [0, 0, 255]

        # export
        hdf5.export(FILENAME, data, thumbnail)

        # check it's here
        st = os.stat(FILENAME) # this test also that the file is created
        self.assertGreater(st.st_size, 0)
        f = h5py.File(FILENAME, "r")

        # check format
        im = f["Acquisition0/ImageData/Image"]
        self.assertEqual(im.attrs["IMAGE_SUBCLASS"], "IMAGE_GRAYSCALE")

        # check basic metadata
        self.assertEqual(im.dims[4].label, "X")
        yres = im.dims[3][0][()] # second last dimension (Y), first scale, first and only value
        self.assertAlmostEqual(metadata[model.MD_PIXEL_SIZE][1], yres)
        ypos = f["Acquisition0/ImageData/YOffset"][()]
        self.assertAlmostEqual(metadata[model.MD_POS][1], ypos)

        # Check physical metadata
        desc = f["Acquisition0/PhysicalData/Title"][()]
        self.assertAlmostEqual(metadata[model.MD_DESCRIPTION], desc)

        iwl = f["Acquisition0/PhysicalData/ExcitationWavelength"][()] # m
        self.assertTrue((metadata[model.MD_IN_WL][0] <= iwl and
                         iwl <= metadata[model.MD_IN_WL][1]))

        expt = f["Acquisition0/PhysicalData/IntegrationTime"][()] # s
        self.assertAlmostEqual(metadata[model.MD_EXP_TIME], expt)


    def testExportRead(self):
        """
        Checks that we can read back an image and a thumbnail
        """
        # create 2 simple greyscale images
        sizes = [(512, 256), (500, 400)] # different sizes to ensure different acquisitions
        dtype = numpy.dtype("uint16")
        white = (12, 52) # non symmetric position
        ldata = []
        num = 2
        # TODO: check support for combining channels when same data shape
        for i in range(num):
            a = model.DataArray(numpy.zeros(sizes[i][-1:-3:-1], dtype))
            a[white[-1:-3:-1]] = 1027
            ldata.append(a)

        # thumbnail : small RGB completely red
        tshape = (sizes[0][1] // 8, sizes[0][0] // 8, 3)
        tdtype = numpy.uint8
        thumbnail = model.DataArray(numpy.zeros(tshape, tdtype))
        thumbnail[:, :, 0] += 255 # red
        blue = (12, 22) # non symmetric position
        thumbnail[blue[-1:-3:-1]] = [0, 0, 255]

        # export
        hdf5.export(FILENAME, ldata, thumbnail)

        # check it's here
        st = os.stat(FILENAME) # this test also that the file is created
        self.assertGreater(st.st_size, 0)

        # check data
        rdata = hdf5.read_data(FILENAME)
        self.assertEqual(len(rdata), num)

        for i, im in enumerate(rdata):
            subim = im[0, 0, 0] # remove C,T,Z dimensions
            self.assertEqual(subim.shape, sizes[i][-1::-1])
            self.assertEqual(subim[white[-1:-3:-1]], ldata[i][white[-1:-3:-1]])

        # check thumbnail
        rthumbs = hdf5.read_thumbnail(FILENAME)
        self.assertEqual(len(rthumbs), 1)
        im = rthumbs[0]
        self.assertEqual(im.shape, tshape)
        self.assertEqual(im[0, 0].tolist(), [255, 0, 0])
        self.assertEqual(im[blue[-1:-3:-1]].tolist(), [0, 0, 255])


    def testReadMDSpec(self):
        """
        Checks that we can read back the metadata of an image
        """
        metadata = [{model.MD_SW_VERSION: "1.0-test",
                     model.MD_HW_NAME: "fake hw",
                     model.MD_DESCRIPTION: "test",
                     model.MD_ACQ_DATE: time.time(),
                     model.MD_BPP: 12,
                     model.MD_BINNING: (1, 2), # px, px
                     model.MD_PIXEL_SIZE: (1e-6, 2e-5), # m/px
                     model.MD_POS: (1e-3, -30e-3), # m
                     model.MD_EXP_TIME: 1.2, # s
                     model.MD_LENS_MAG: 1200, # ratio
                    },
                    {model.MD_SW_VERSION: "1.0-test",
                     model.MD_HW_NAME: "fake spec",
                     model.MD_DESCRIPTION: "test3d",
                     model.MD_ACQ_DATE: time.time(),
                     model.MD_BPP: 12,
                     model.MD_BINNING: (1, 1), # px, px
                     model.MD_PIXEL_SIZE: (1e-6, 2e-5), # m/px
                     model.MD_WL_POLYNOMIAL: [500e-9, 1e-9], # m, m/px: wl polynomial
                     model.MD_POS: (1e-3, -30e-3), # m
                     model.MD_EXP_TIME: 1.2, # s
                    },
                    ]
        # create 2 simple greyscale images
        sizes = [(512, 256), (500, 400, 1, 1, 220)] # different sizes to ensure different acquisitions
        dtype = numpy.dtype("uint8")
        ldata = []
        for i, s in enumerate(sizes):
            a = model.DataArray(numpy.zeros(s[::-1], dtype), metadata[i])
            ldata.append(a)

        # thumbnail : small RGB completely red
        tshape = (sizes[0][1] // 8, sizes[0][0] // 8, 3)
        tdtype = numpy.uint8
        thumbnail = model.DataArray(numpy.zeros(tshape, tdtype))
        thumbnail[:, :, 1] += 255 # green

        # export
        hdf5.export(FILENAME, ldata, thumbnail)

        # check it's here
        st = os.stat(FILENAME) # this test also that the file is created
        self.assertGreater(st.st_size, 0)

        # check data
        rdata = hdf5.read_data(FILENAME)
        self.assertEqual(len(rdata), len(ldata))

        for i, im in enumerate(rdata):
            md = metadata[i]
            self.assertEqual(im.metadata[model.MD_DESCRIPTION], md[model.MD_DESCRIPTION])
            self.assertEqual(im.metadata[model.MD_POS], md[model.MD_POS])
            self.assertEqual(im.metadata[model.MD_PIXEL_SIZE], md[model.MD_PIXEL_SIZE])
            self.assertEqual(im.metadata[model.MD_ACQ_DATE], md[model.MD_ACQ_DATE])
            if model.MD_LENS_MAG in md:
                self.assertEqual(im.metadata[model.MD_LENS_MAG], md[model.MD_LENS_MAG])

            # None of the images are using light => no MD_IN_WL
            self.assertFalse(model.MD_IN_WL in im.metadata,
                             "Reporting excitation wavelength while there is none")

            if model.MD_WL_POLYNOMIAL in md:
                pn = md[model.MD_WL_POLYNOMIAL]
                # 2 formats possible
                if model.MD_WL_LIST in im.metadata:
                    l = ldata[i].shape[2]
                    npn = polynomial.Polynomial(pn,
                                    domain=[0, l - 1],
                                    window=[0, l - 1])
                    wl = npn.linspace(l)
                    self.assertEqual(im.metadata[model.MD_WL_LIST], wl)
                else:
                    self.assertEqual(im.metadata[model.MD_WL_POLYNOMIAL], pn)

        # check thumbnail
        rthumbs = hdf5.read_thumbnail(FILENAME)
        self.assertEqual(len(rthumbs), 1)
        im = rthumbs[0]
        self.assertEqual(im.shape, tshape)
        self.assertEqual(im[0, 0].tolist(), [0, 255, 0])

# TODO: test compatibility with Hyperspy for loading spectra (exported by Odemis)





    def testReadMDFluo(self):
        """
        Checks that we can read back the metadata of a fluoresence image
        The OME-TIFF file will contain just one big array, but three arrays 
        should be read back with the right data.
        """
        metadata = [{model.MD_SW_VERSION: "1.0-test",
                     model.MD_HW_NAME: "fake hw",
                     model.MD_DESCRIPTION: "brightfield",
                     model.MD_ACQ_DATE: time.time(),
                     model.MD_BPP: 12,
                     model.MD_BINNING: (1, 1), # px, px
                     model.MD_PIXEL_SIZE: (1e-6, 1e-6), # m/px
                     model.MD_POS: (13.7e-3, -30e-3), # m
                     model.MD_EXP_TIME: 1.2, # s
                     model.MD_IN_WL: (400e-9, 630e-9), # m
                     model.MD_OUT_WL: (400e-9, 630e-9), # m
                    },
                    {model.MD_SW_VERSION: "1.0-test",
                     model.MD_HW_NAME: "fake hw",
                     model.MD_DESCRIPTION: "blue dye",
                     model.MD_ACQ_DATE: time.time() + 10,
                     model.MD_BPP: 12,
                     model.MD_BINNING: (1, 1), # px, px
                     model.MD_PIXEL_SIZE: (1e-6, 1e-6), # m/px
                     model.MD_POS: (13.7e-3, -30e-3), # m
                     model.MD_EXP_TIME: 1.2, # s
                     model.MD_IN_WL: (500e-9, 520e-9), # m
                     model.MD_OUT_WL: (600e-9, 630e-9), # m
                    },
                    {model.MD_SW_VERSION: "1.0-test",
                     model.MD_HW_NAME: "fake hw",
                     model.MD_DESCRIPTION: "green dye",
                     model.MD_ACQ_DATE: time.time() + 20,
                     model.MD_BPP: 12,
                     model.MD_BINNING: (1, 1), # px, px
                     model.MD_PIXEL_SIZE: (1e-6, 1e-6), # m/px
                     model.MD_POS: (13.7e-3, -30e-3), # m
                     model.MD_EXP_TIME: 1, # s
                     model.MD_IN_WL: (600e-9, 620e-9), # m
                     model.MD_OUT_WL: (620e-9, 650e-9), # m
                    },
                    ]
        # create 3 greyscale images of same size
        size = (512, 256)
        dtype = numpy.dtype("uint16")
        ldata = []
        for i, md in enumerate(metadata):
            a = model.DataArray(numpy.zeros(size[::-1], dtype), md)
            a[i, i] = i # "watermark" it
            ldata.append(a)

        # thumbnail : small RGB completely red
        tshape = (size[1] // 8, size[0] // 8, 3)
        tdtype = numpy.uint8
        thumbnail = model.DataArray(numpy.zeros(tshape, tdtype))
        thumbnail[:, :, 1] += 255 # green

        # export
        hdf5.export(FILENAME, ldata, thumbnail)

        # check it's here
        st = os.stat(FILENAME) # this test also that the file is created
        self.assertGreater(st.st_size, 0)

        # check data
        rdata = hdf5.read_data(FILENAME)
        self.assertEqual(len(rdata), len(ldata))

        # TODO: rdata and ldata don't have to be in the same order
        for i, im in enumerate(rdata):
            md = metadata[i]
            self.assertEqual(im.metadata[model.MD_DESCRIPTION], md[model.MD_DESCRIPTION])
            self.assertAlmostEqual(im.metadata[model.MD_POS][0], md[model.MD_POS][0])
            self.assertAlmostEqual(im.metadata[model.MD_POS][1], md[model.MD_POS][1])
            self.assertAlmostEqual(im.metadata[model.MD_PIXEL_SIZE][0], md[model.MD_PIXEL_SIZE][0])
            self.assertAlmostEqual(im.metadata[model.MD_PIXEL_SIZE][1], md[model.MD_PIXEL_SIZE][1])

            iwl = im.metadata[model.MD_IN_WL] # nm
            self.assertTrue((md[model.MD_IN_WL][0] <= iwl[0] and
                             iwl[1] <= md[model.MD_IN_WL][1]))

            owl = im.metadata[model.MD_OUT_WL] # nm
            self.assertTrue((md[model.MD_OUT_WL][0] <= owl[0] and
                             owl[1] <= md[model.MD_OUT_WL][1]))

            # SVI HDF5 only records one acq time per T dimension
            # so only check for the first channel
            if i == 0:
                self.assertAlmostEqual(im.metadata[model.MD_ACQ_DATE], md[model.MD_ACQ_DATE], delta=1)

            # SVI HDF5 doesn't this metadata:
#            self.assertEqual(im.metadata[model.MD_BPP], md[model.MD_BPP])
#            self.assertEqual(im.metadata[model.MD_BINNING], md[model.MD_BINNING])
            self.assertEqual(im.metadata[model.MD_EXP_TIME], md[model.MD_EXP_TIME])


        # check thumbnail
        rthumbs = hdf5.read_thumbnail(FILENAME)
        self.assertEqual(len(rthumbs), 1)
        im = rthumbs[0]
        self.assertEqual(im.shape, tshape)
        self.assertEqual(im[0, 0].tolist(), [0, 255, 0])


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

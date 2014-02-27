# -*- coding: utf-8 -*-
"""
Created on 3 Jan 2014

@author: kimon

Copyright © 2013-2014 Éric Piel & Kimon Tsitsikas, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the
terms  of the GNU General Public License version 2 as published by the Free
Software  Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR  PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Odemis. If not, see http://www.gnu.org/licenses/.
"""

from __future__ import division

import logging
import numpy
import math

from numpy import arange
from numpy import fft

def CalculateDrift(previous_img, current_img, precision=1):
    """
    Given two images, it calculates the drift in x and y axis. It first computes
    the cross-correlation of the two images and then locates the peak. The coordinates 
    of the peak of the cross-correlation define the shift vector between the two images. 
    previous_img (numpy.array): 2d array with the previous frame
    current_img (numpy.array): 2d array with the last frame, must be of same 
      shape as previous_img
    precision (1<=int): Calculate drift within 1/precision of a pixel
    returns (tuple of floats): Drift in pixels
    """
    if precision < 1:
        raise ValueError("Precision cannot be less than 1, got %s." % (precision,))
    assert previous_img.shape == current_img.shape

    previous_fft = fft.fft2(previous_img)
    current_fft = fft.fft2(current_img)
    (m, n) = previous_fft.shape

    if precision == 1:
        # Cross-correlation computation
        CC = fft.ifft2(previous_fft * current_fft.conj())

        # Locate the peak
        ACC = abs(CC)
        loc1 = ACC.argmax(0)
        max1 = ACC[(loc1, range(ACC.shape[1]))]
        loc2 = max1.argmax(0)
#        max1, loc1 = abs(CC).max(0), abs(CC).argmax(0)
#        max2, loc2 = abs(max1).max(0), abs(max1).argmax(0)

        rloc = loc1[loc2]
        cloc = loc2

        # Calculate shift from the peak
        md2 = numpy.fix(m / 2)
        nd2 = numpy.fix(n / 2)
        if rloc > md2:
            row_shift = rloc - m
        else:
            row_shift = rloc

        if cloc > nd2:
            col_shift = cloc - n
        else:
            col_shift = cloc

    else:
        mlarge, nlarge = m * 2, n * 2

        # Upsample by factor of 2 to obtain initial estimation and
        # embed Fourier data in a 2x larger array
        CC = numpy.zeros((mlarge, nlarge), dtype=numpy.complex)
        CC[m - numpy.fix(m / 2):m + 1 + numpy.fix((m - 1) / 2),
           n - numpy.fix(n / 2):n + 1 + numpy.fix((n - 1) / 2)] = \
				fft.fftshift(previous_fft) * fft.fftshift(current_fft).conj()

        # Cross-correlation computation
        CC = fft.ifft2(fft.ifftshift(CC))

        # Locate the peak
        ACC = abs(CC)
        loc1 = ACC.argmax(0)
        max1 = ACC[(loc1, range(ACC.shape[1]))]
        loc2 = max1.argmax(0)

        rloc = loc1[loc2]
        cloc = loc2
       
        # Calculate shift in previous pixel grid from the position of the peak
        (m, n) = CC.shape
        md2 = numpy.fix(m / 2)
        nd2 = numpy.fix(n / 2)

        if rloc > md2:
            row_shift = rloc - m
        else:
            row_shift = rloc

        if cloc > nd2:
            col_shift = cloc - n
        else:
            col_shift = cloc

        row_shift = row_shift / 2
        col_shift = col_shift / 2

        # DFT computation
        # Initial shift estimation in upsampled grid
        row_shift = numpy.round(row_shift * precision) / precision
        col_shift = numpy.round(col_shift * precision) / precision
        dft_shift = numpy.fix(numpy.ceil(precision * 1.5) / 2)  # Center of output at dft_shift+1

        # Matrix multiply DFT around the current shift estimation
        CC = (_UpsampledDFT(current_fft * previous_fft.conj(),
                            numpy.ceil(precision * 1.5),
                            numpy.ceil(precision * 1.5),
                            precision,
                            dft_shift - row_shift * precision,
                            dft_shift - col_shift * precision)
              ) / (md2 * nd2 * (precision ** 2))
        # was .conj(), but as we just need the abs(), it's not needed

        # Locate maximum and map back to original pixel grid
        ACC = abs(CC)
        loc1 = ACC.argmax(0)
        max1 = ACC[(loc1, range(ACC.shape[1]))]
        loc2 = max1.argmax(0)

        rloc = loc1[loc2]
        cloc = loc2

        rloc -= dft_shift
        cloc -= dft_shift

        row_shift += rloc / precision
        col_shift += cloc / precision

        if md2 == 1:
            row_shift = 0
        if nd2 == 1:
            col_shift = 0

    return row_shift, col_shift


def _UpsampledDFT(data, nor, noc, precision=1, roff=0, coff=0):
    """
    Upsampled DFT by matrix multiplies. 
    data (numpy.array): 2d array 
    nor, noc (ints): Number of pixels in the output upsampled DFT, in units
    of upsampled pixels
    precision (int): Calculate drift within 1/precision of a pixel
    roff, coff (ints): Row and column offsets, allow to shift the output array
                    to a region of interest on the DFT 
    returns (tuple of floats): Drift in pixels
    """
    z = 1j  # imaginary unit
    nr, nc = data.shape

    # Compute kernels and obtain DFT by matrix products
    kernc = numpy.power(math.e, (-z * 2 * math.pi / (nc * precision)) * ((fft.ifftshift((arange(0, nc)))[:, None]).transpose() \
                                - numpy.floor(nc / 2)) * (arange(0, noc) - coff)[:, None])

    kernr = numpy.power(math.e, (-z * 2 * math.pi / (nr * precision)) * ((fft.ifftshift(arange(0, nr)))[:, None] \
                                - numpy.floor(nr / 2)) * ((arange(0, nor)[:, None]).transpose() - roff))

    return numpy.dot(numpy.dot((kernr.transpose()), data), kernc.transpose())

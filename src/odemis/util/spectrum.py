# -*- coding: utf-8 -*-

"""
Created on 28 Feb 2014

@author: Éric Piel

Copyright © 2014 Éric Piel, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License version 2 as published by the Free
Software Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Odemis. If not, see http://www.gnu.org/licenses/.

"""

# Various function to handle spectral data
# Note that the spectrum is normally contained on the C dimension, which is
# by convention the first of the 5 dimensions of a DataArray (CTZYX).

from __future__ import division

import logging
from numpy.polynomial import polynomial
from odemis import model


def get_wavelength_per_pixel(da):
    """
    Computes the wavelength for each pixel along the C dimension

    :param da: (model.DataArray of shape C...): the DataArray with metadata
        either MD_WL_POLYNOMIAL or MD_WL_LIST
    :return: (list of float of length C): the wavelength (in m) for each pixel
        in C
    :raises:
        AttributeError: if no metadata is present
        KeyError: if no metadata is available
        ValueError: if the metadata doesn't provide enough information
    """

    if not hasattr(da, 'metadata'):
        raise AttributeError("No metadata found in data array")

    # MD_WL_LIST has priority
    if model.MD_WL_LIST in da.metadata:
        wl = da.metadata[model.MD_WL_LIST]
        if len(wl) == da.shape[0]:
            return wl
        else:
            logging.warning("MD_WL_LIST is not the same length as the data")

    if model.MD_WL_POLYNOMIAL in da.metadata:
        pn = da.metadata[model.MD_WL_POLYNOMIAL]
        pn = polynomial.polytrim(pn)
        if len(pn) >= 2:
            npn = polynomial.Polynomial(pn,  #pylint: disable=E1101
                                        domain=[0, da.shape[0] - 1],
                                        window=[0, da.shape[0] - 1])
            return npn.linspace(da.shape[0])[1]
        else:
            # a polynomial of 0 or 1 value is useless
            raise ValueError("Wavelength polynomial has only %d degree"
                             % len(pn))

    raise KeyError("No MD_WL_* metadata available")


def coefficients_to_dataarray(coef):
    """
    Convert a spectrum efficiency coefficient array to a DataArray as expected
    by odemis.acq.calibration
    coef (numpy array of shape (N,2)): first column is the wavelength in nm, second
      column is the coefficient (float > 0)
    returns (DataArray of shape (N,1,1,1,1)): the same data with the wavelength
     encoded in the metadata, as understood by the rest of Odemis.
    """
    n = coef.shape[0]

    # Create the content of the DataArray directly from the second column
    da = model.DataArray(coef[:, 1])
    da.shape += (1, 1, 1, 1) # add another 4 dims

    # metadata from the first column, converting from nm to m
    wl_list = (1e-9 * coef[:, 0]).tolist()
    da.metadata[model.MD_WL_LIST] = wl_list

    return da

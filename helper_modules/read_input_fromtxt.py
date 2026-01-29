#!/usr/bin/env python3

"""
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de

"""

import hoomd
from hoomd import *
from hoomd import sph
import numpy as np
import itertools
import gsd.hoomd
import os
import array 


def sanity_check_input(input_dict):
    """
    

    Parameters
    ----------
    input_dict : dict
        all input data from get_input_data_from_file().

    Returns
    -------
    None.

    """

    fname = input_dict['rawfilename']

    # check if file exists
    if os.path.exists(fname) == False:
        raise ValueError(f'File {fname} does not exist!')

    # check kernel
    list_of_kernels = ['WendlandC2', 'WendlandC4', 'WendlandC6', 'Quintic', 'CubicSpline']
    kernel = input_dict['kernel']
    if kernel not in list_of_kernels:
        raise ValueError(f'Kernel {kernel} not a member of available kernels.')

    # check domain sizes
    nx = input_dict['nx']
    ny = input_dict['ny']
    nz = input_dict['nz']

    if nx <= 0 or ny <= 0 or nz <= 0:
        raise ValueError('Domain sizes can not be <= 0.')

    # check fluid parameters
    fdensity = input_dict['fdensity']
    if fdensity <= 0:
        raise ValueError('Density has to be > 0.')

    fviscosity = input_dict['fviscosity']
    if fviscosity <= 0:
        raise ValueError('Viscosity has to be > 0.')

    d_handle = input_dict['delete_flag']
    if d_handle != 0 and d_handle != 1:
        raise ValueError('Flag on deleting solids not set properly.')

    porosity = input_dict['porosity']
    if porosity < 0.0 or porosity > 1.0:
        raise ValueError('Porosity has to be 0 <= phi <= 1.')



def get_input_data_from_file(inputfile):
    """
    
    Parameters
    ----------
    inputfile : file
        .txt file with all input parameter.

    Returns
    -------
    dictonary.

    """
    parameter_dict = {}

    strs = np.genfromtxt(inputfile, str, skip_header = 2, skip_footer = 7, usecols = (0, ))
    flts = np.genfromtxt(inputfile, np.float64, skip_header = 4, usecols = (0, ))


    parameter_dict.update({'rawfilename' : str(strs[0]) })
    parameter_dict.update({'kernel'      : str(strs[1]) })

    parameter_dict.update({'vsize'       : np.float64(flts[0]) })
    parameter_dict.update({'nx'          : np.int32(flts[1]) })
    parameter_dict.update({'ny'          : np.int32(flts[2]) })
    parameter_dict.update({'nz'          : np.int32(flts[3]) })
    parameter_dict.update({'fdensity'    : np.float64(flts[4]) })
    parameter_dict.update({'fviscosity'  : np.float64(flts[5]) })
    parameter_dict.update({'delete_flag' : np.int32(flts[6]) })
    parameter_dict.update({'porosity'    : np.float64(flts[7]) })

    sanity_check_input(parameter_dict)

    return parameter_dict

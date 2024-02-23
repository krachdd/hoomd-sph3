#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""
# ----- HEADER -----------------------------------------------
import hoomd
from hoomd import *
from hoomd import sph
from hoomd.sph import _sph
import numpy as np
import math
# import itertools
from datetime import datetime
import export_gsd2vtu, delete_solids_initial_timestep 
import sph_info, sph_helper, read_input_fromtxt
import os, sys
from optparse import OptionParser
import array

import gsd.hoomd
# -----------------------------------------------------------

# get stuff from input file
infile = str(sys.argv[1])
params = read_input_fromtxt.get_input_data_from_file(infile)
nx, ny, nz = np.int32(params['nx']), np.int32(params['ny']), np.int32(params['nz']) 
voxelsize           = np.float64(params['vsize'])               # [ m ]

rawfile = params['rawfilename']
rawfile = open(rawfile,'rb')
    
tids = np.fromfile(rawfile, dtype= np.uint8, count=nx*ny*nz)
tids = np.reshape(tids, (nx, ny, nz), order = 'F')
    
rawfile.close()

porosity = np.sum(tids)/(nx * ny * nz)
print(f'Porosity: {porosity}')




for i in range(11):
    nodes = i
    alls = []
    for j in range(i):
        # print(f'{i} and {j}')
        alls.append(tids)
    if i == 0: continue
    d = np.concatenate(alls, axis = 2)
    print(d.shape)

    fn = f'lowporous_sym_{d.shape[0]}_{d.shape[1]}_{d.shape[2]}_vs_1e-3_nodes{i}.raw'
    d = d.flatten('F')
    d = np.array(d, dtype=np.bool_)
    file = open(fn, mode = 'w')
    d.tofile(file)
    file.close()
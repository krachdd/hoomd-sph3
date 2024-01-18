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
import export_gsd2vtu 
import read_input_fromtxt
import delete_solids_initial_timestep
import sys, os

import gsd.hoomd
# ------------------------------------------------------------

filename = str(sys.argv[1])

device = hoomd.device.CPU(notice_level=2)
# device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)

filename = str(sys.argv[1])

if device.communicator.rank == 0:
    print(f'{os.path.basename(__file__)}: input file: {filename} ')

sim.create_state_from_gsd(filename = filename)

# Print the number of particles on each rank.
with sim.state.cpu_local_snapshot as snap:
    N = len(snap.particles.auxiliary3)
    for i in range(N):
        print(snap.particles.auxiliary3[i])
    print(f'{N} particles on rank {device.communicator.rank}')


if device.communicator.rank == 0:
    export_gsd2vtu.export_all(filename)
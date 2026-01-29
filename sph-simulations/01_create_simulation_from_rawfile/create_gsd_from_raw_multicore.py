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
import delete_solids_initial_timestep
import read_input_fromtxt


# device = hoomd.device.CPU(notice_level=2)
device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)

# get stuff from input file
infile = str(sys.argv[1])
params = read_input_fromtxt.get_input_data_from_file(infile)
if device.communicator.rank == 0:
    print(params)

# Fluid and particle properties
voxelsize  = np.float64(params['vsize'])
DX   = voxelsize
V    = DX * DX * DX
RHO0 = np.float64(params['fdensity'])
MU   = np.float64(params['fviscosity'])
M    = RHO0 * V

# get kernel properties
KERNEL  = params['kernel']
H       = hoomd.sph.kernel.OptimalH[KERNEL]*DX       # m
RCUT    = hoomd.sph.kernel.Kappa[KERNEL]*H           # m

# get simulation box sizes etc.
NX, NY, NZ = np.int32(params['nx']), np.int32(params['ny']), np.int32(params['nz']) 
LX, LY, LZ = NX*voxelsize, NY*voxelsize, NZ*voxelsize
# box dimensions
box_Lx, box_Ly, box_Lz = LX, LY, LZ  

# Number of Particles
N_particles = NX * NY * NZ 

# get Type id data from raw file
rawfile = params['rawfilename']
rawf_handle = open(rawfile, 'rb')
tids = array.array("B")
tids.fromfile(rawf_handle, NX * NY * NZ)
tids = np.array(tids, dtype = np.uint8)
tids = tids.reshape((NZ, NY, NX))
tids = tids.flatten(order = 'F')
rawf_handle.close()
porosity = np.sum(tids)/(NX * NY * NZ)



# define meshgrid and add properties
x, y, z = np.meshgrid(*(np.linspace(-box_Lx / 2, box_Lx / 2, NX, endpoint=False),),
                      *(np.linspace(-box_Ly / 2, box_Ly / 2, NY, endpoint=False),),
                      *(np.linspace(-box_Lz / 2, box_Lz / 2, NZ, endpoint=False),))

positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
masses     = np.ones((positions.shape[0]), dtype = np.float32) * M
slengths   = np.ones((positions.shape[0]), dtype = np.float32) * H
density    = np.ones((positions.shape[0]), dtype = np.float32) * RHO0
# dpes       = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)


# print(positions)

snapshot = hoomd.Snapshot(device.communicator)
if snapshot.communicator.rank == 0:
    snapshot.particles.N = N_particles
    snapshot.configuration.box = [box_Lx, box_Ly, box_Lz] + [0, 0, 0]
    snapshot.particles.position[:] = positions
    snapshot.particles.typeid[:] = tids
    snapshot.particles.types = ['F', 'S']
    snapshot.particles.velocity[:] = velocities
    snapshot.particles.mass[:] = masses
    snapshot.particles.slength[:] = slengths
    snapshot.particles.density[:] = density



# # create Snapshot 
# snapshot = hoomd.Snapshot(device.communicator)
# snapshot.configuration.box = [box_Lx, box_Ly, box_Lz] + [0, 0, 0]
# snapshot.particles.N = N_particles
# snapshot.particles.position[:] = positions
# snapshot.particles.typeid[:] = tids
# snapshot.particles.types = ['F', 'S']
# snapshot.particles.velocity[:] = velocities
# snapshot.particles.mass[:] = masses
# snapshot.particles.slength[:] = slengths
# snapshot.particles.dpe[:] = dpes

sim.create_state_from_snapshot(snapshot)

deletesolid_flag = params['delete_flag']
if deletesolid_flag == 1:
    if device.communicator.rank == 0:
        print(f'Delete solid particles')
    sim, ndel_particles = delete_solids_initial_timestep.delete_solids(sim, device, KERNEL, 0.000001, MU, DX, RHO0)
    N_particles = N_particles - ndel_particles


init_filename = rawfile.replace('.raw', '_init.gsd')

hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

print(f'Filename: {init_filename}, Number of particles: {N_particles}')

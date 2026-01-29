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
# ----- HEADER -----------------------------------------------
import hoomd
from hoomd import *
from hoomd import sph
from hoomd.sph import _sph
import numpy as np
from mpi4py import MPI
import math
# import itertools
from datetime import datetime
# import export_pgsd2vtu 
import read_input_fromtxt
import delete_solids_initial_timestep
import sys, os
import array
import logging

import pgsd.hoomd
# ------------------------------------------------------------

# logging.basicConfig(filename='parallel_gsd.log', encoding='utf-8', level=logging.DEBUG)

# device = hoomd.device.CPU(notice_level=2)
device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)

new_comm = MPI.COMM_WORLD
rank = new_comm.Get_rank()
size = new_comm.Get_size()
# new_comm = device.communicator
# rank = device.communicator.rank
# size = device.communicator.num_ranks
offset = 0


# get stuff from input file
infile = str(sys.argv[1])
params = read_input_fromtxt.get_input_data_from_file(infile)
if rank == 0:
    print(params)

# Fluid and particle properties
voxelsize        = np.float64(params['vsize'])       
dx               = voxelsize
specific_volume  = dx * dx * dx
rho0             = np.float64(params['fdensity'])
viscosity        = np.float64(params['fviscosity'])
mass             = rho0 * specific_volume

# get kernel properties
kernel  = params['kernel']
slength = hoomd.sph.kernel.OptimalH[kernel]*dx       # [ m ]
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength     # [ m ]



# get simulation box sizes etc.
nx, ny, nz = np.int32(params['nx']), np.int32(params['ny']), np.int32(params['nz']) 
lx, ly, lz = nx*voxelsize, ny*voxelsize, nz*voxelsize
# box dimensions
box_lx, box_ly, box_lz = lx, ly, lz 

# Number of Particles
n_particles = nx * ny* nz

# get Type id data from raw file
rawfile = params['rawfilename']
rawf_handle = open(rawfile, 'rb')
tids = array.array("B")
tids.fromfile(rawf_handle, nx * ny * nz)
tids = np.array(tids, dtype = np.uint8)
tids = tids.reshape((nz, ny, nx))
tids = tids.flatten(order = 'F')
rawf_handle.close()
porosity = 1.0 - np.sum(tids)/(nx * ny * nz)



# get rank specific parameters
if n_particles != tids.shape[0]:
    raise ValueError('Number of particles does not match the number of Type IDs.')

# Compute number of particles per rank
n_particles_rank = math.floor( n_particles/size )

if rank < n_particles%size:
    n_particles_rank = n_particles_rank + 1
for i in range(rank):
    offset = offset + math.floor( n_particles/size )
    if rank < n_particles%size:
        offset = offset+1


# define meshgrid and add properties
x, y, z = np.meshgrid(*(np.linspace(-box_lx / 2, box_lx / 2, nx, endpoint=True),),
                      *(np.linspace(-box_ly / 2, box_ly / 2, ny, endpoint=True),),
                      *(np.linspace(-box_lz / 2, box_lz / 2, nz, endpoint=True),))

x = x.ravel()
y = y.ravel()
z = z.ravel()

x = x[offset:offset+n_particles_rank]
y = y[offset:offset+n_particles_rank]
z = z[offset:offset+n_particles_rank]

positions = np.array((x.ravel(), y.ravel(), z.ravel())).T
tids      = tids[offset:offset+n_particles_rank]

velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
masses     = np.ones((positions.shape[0]), dtype = np.float32) * mass
slengths   = np.ones((positions.shape[0]), dtype = np.float32) * slength
densities  = np.ones((positions.shape[0]), dtype = np.float32) * rho0

part_dist = np.zeros(size)
a = n_particles_rank*np.ones(1)
new_comm.Allgather([a, MPI.FLOAT], [part_dist, MPI.FLOAT])

if rank == 0:
    print(f'[INFO]: Particle Distribution on {size} ranks: {part_dist}')

# # create Snapshot 
# # snapshot = hoomd.Snapshot(device.communicator)
snapshot = pgsd.hoomd.Frame(size)
snapshot.part_dist = part_dist
snapshot.configuration.box     = [box_lx, box_ly, box_lz] + [0, 0, 0]
snapshot.particles.N           = n_particles_rank
snapshot.particles.position    = positions
snapshot.particles.typeid      = tids
snapshot.particles.types       = ['F', 'S']
snapshot.particles.velocity    = velocities
snapshot.particles.mass        = masses
snapshot.particles.slength     = slengths
snapshot.particles.density     = densities



sim.create_state_from_snapshot(snapshot)

# # deletesolid_flag = params['delete_flag']
# # if deletesolid_flag == 1:
# #     print(f'Delete solid particles')
# #     sim, ndel_particles = delete_solids_initial_timestep.delete_solids(sim, device, kernel, 0.000001, viscosity, dx, rho0)
# #     n_particles = n_particles - ndel_particles

init_filename = rawfile.replace('.raw', '_init.gsd')
hoomd.write.PGSD.write(state = sim.state, mode = 'wb', filename = init_filename)

# print(f'Filename: {init_filename}, Number of particles: {n_particles}')

# if device.communicator.rank == 0:
#     export_pgsd2vtu.export_spf(init_filename)

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
import math
import itertools
from datetime import datetime

# import export_gsd2vtu 
# import read_input_fromtxt
# import delete_solids_initial_timestep

import gsd.hoomd
# ------------------------------------------------------------



device = hoomd.device.CPU(notice_level=2)
# device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)


# Fluid and particle properties
num_length          = 400
lref                = 1.2
voxelsize           = lref/num_length
dx                  = voxelsize
specific_volume     = dx * dx * dx
rho0                = 1000.0
mass                = rho0 * specific_volume
refvel              = 10
Re                  = 10
viscosity           = rho0 * lref * refvel/Re



# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx                  # [ m ]
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength                # [ m ]

# particles per Kernel Radius
part_rcut  = math.ceil(rcut/dx) 
part_depth = math.ceil(1.5 * hoomd.sph.kernel.Kappa[kernel] * rcut/dx) 


# get simulation box sizes etc.
nx, ny, nz = num_length + (2*part_rcut), num_length + (2*part_rcut) , part_depth
lx, ly, lz = nx*voxelsize, ny*voxelsize, nz*voxelsize
# box dimensions
box_lx, box_ly, box_lz = lx, ly, lz  

print(lx, ly, lz)

print(part_rcut)

# Number of Particles
n_particles = nx * ny * nz 

# define meshgrid and add properties
x, y, z = np.meshgrid(*(np.linspace(-box_lx / 2, box_lx / 2, nx, endpoint=True),),
                      *(np.linspace(-box_ly / 2, box_ly / 2, ny, endpoint=True),),
                      *(np.linspace(-box_lz / 2, box_lz / 2, nz, endpoint=True),))

positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
masses     = np.ones((positions.shape[0]), dtype = np.float32) * mass
slengths   = np.ones((positions.shape[0]), dtype = np.float32) * slength
densities  = np.ones((positions.shape[0]), dtype = np.float32) * rho0

# # # create Snapshot 
snapshot = gsd.hoomd.Snapshot()
snapshot.configuration.box     = [box_lx, box_ly, box_lz] + [0, 0, 0]
snapshot.particles.N           = n_particles
snapshot.particles.position    = positions
snapshot.particles.typeid      = [0] * n_particles
snapshot.particles.types       = ['F', 'S']
snapshot.particles.velocity    = velocities
snapshot.particles.mass        = masses
snapshot.particles.slength     = slengths
snapshot.particles.density     = densities


x    = snapshot.particles.position[:]
tid  = snapshot.particles.typeid[:]
vels = snapshot.particles.velocity[:]
for i in range(len(x)):
    xi,yi,zi  = x[i][0], x[i][1], x[i][2]
    tid[i]    = 0
    # solid walls 
    if ( xi < -0.5 * lref or xi > 0.5 * lref or yi < -0.5 * lref or yi > 0.5 * lref):
        tid[i] = 1
    if ( yi >= 0.5 * lref):
        vels[i][0] = refvel

snapshot.particles.typeid[:]     = tid
snapshot.particles.velocity[:]   = vels

sim.create_state_from_snapshot(snapshot)

init_filename = f'liddrivencavity_{nx}_{ny}_{nz}_vs_{voxelsize}_init.gsd'
# hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

with gsd.hoomd.open(name = init_filename, mode = 'w') as f:
    f.append(snapshot)

# if device.communicator.rank == 0:
#     export_gsd2vtu.export_spf(init_filename)

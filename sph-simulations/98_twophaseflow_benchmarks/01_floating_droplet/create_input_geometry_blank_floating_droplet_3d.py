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


device = hoomd.device.CPU(notice_level=2)
# device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)

# Fluid and particle properties
num_length          = int(sys.argv[1])                #  int
lref                = 0.001                           # [ m ]
radius              = 0.5 * lref                      # [ m ]
voxelsize           = lref/num_length                 # [ m ]
dx                  = voxelsize                       # [ m ]
specific_volume     = dx * dx * dx                    # [ m**3 ]
rho01               = 1000.0                          # [ kg/m**3 ]
rho02               = 1000.0                          # [ kg/m**3 ]
mass1               = rho01 * specific_volume         # [ kg ]
mass2               = rho02 * specific_volume         # [ kg ]
fx                  = 0.0                             # [ m/s ]
viscosity1          = 0.001                           # [ Pa s ]
viscosity2          = 0.001                           # [ Pa s ]
backbressure        = 0.3                             # [ - ]
sigma               = 0.01                            # [ Pa/m**2 ]
contact_angle       = 60                              # [ ° ]

# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx        # [ m ]
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength      # [ m ]

# particles per Kernel Radius
part_rcut  = math.ceil(rcut/dx) 
part_depth = math.ceil(2.5 * hoomd.sph.kernel.Kappa[kernel] * rcut/dx) 

# get simulation box sizes etc.
nx, ny, nz = int(num_length), int(num_length), int(num_length)
lx, ly, lz = float(nx) * voxelsize, float(ny) * voxelsize, float(nz) * voxelsize
# box dimensions
box_lx, box_ly, box_lz = lx, ly, lz

# Number of Particles
n_particles = nx * ny * nz 

# define meshgrid and add properties
# define meshgrid and add properties
x, y, z = np.meshgrid(*(np.linspace(-box_lx / 2 + (dx/2), box_lx / 2 - (dx/2), nx, endpoint=True),),
                      *(np.linspace(-box_ly / 2 + (dx/2), box_ly / 2 - (dx/2), ny, endpoint=True),),
                      *(np.linspace(-box_lz / 2 + (dx/2), box_lz / 2 - (dx/2), nz, endpoint=True),))

positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
masses     = np.ones((positions.shape[0]), dtype = np.float32) * mass1
slengths   = np.ones((positions.shape[0]), dtype = np.float32) * slength
densities  = np.ones((positions.shape[0]), dtype = np.float32) * rho01

# # # create Snapshot 
snapshot = gsd.hoomd.Frame()
snapshot.configuration.box     = [box_lx, box_ly, box_lz] + [0, 0, 0]
snapshot.particles.N           = n_particles
snapshot.particles.position    = positions
snapshot.particles.typeid      = [0] * n_particles
snapshot.particles.types       = ['W', 'N']
snapshot.particles.velocity    = velocities
snapshot.particles.mass        = masses
snapshot.particles.slength     = slengths
snapshot.particles.density     = densities

x    = snapshot.particles.position[:]
tid  = snapshot.particles.typeid[:]
rhos = snapshot.particles.density[:]
mass = snapshot.particles.mass[:]

for i in range(len(x)):
    xi,yi,zi  = x[i][0], x[i][1], x[i][2]
    mass[i]   = mass1
    rhos[i]   = rho01
    # Non wetting fluid
    if (xi > -0.25 * lref and xi < 0.25 * lref and 
        yi > -0.25 * lref and yi < 0.25 * lref and 
        zi > -0.25 * lref and zi < 0.25 *lref):
        mass[i]   = mass1
        rhos[i]   = rho01
        tid[i]    = 1
    # # solid walls
    # if ( xi < -0.5 * lref or xi > 0.5 * lref or 
    #      yi < -0.5 * lref or yi > 0.5 * lref or 
    #      zi < -0.5 * lref or zi > 0.5 * lref):
    #     tid[i] = 2

snapshot.particles.typeid[:]     = tid
snapshot.particles.density[:]    = rhos
snapshot.particles.mass[:]       = mass

sim.create_state_from_snapshot(snapshot)

init_filename = f'blank_droplet_3d_{nx}_{ny}_{nz}_vs_{voxelsize}_init.gsd'
# hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

with gsd.hoomd.open(name = init_filename, mode = 'w') as f:
    f.append(snapshot)

# if device.communicator.rank == 0:
#     export_gsd2vtu.export_spf(init_filename)



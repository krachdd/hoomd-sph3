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

import gsd.hoomd
# ------------------------------------------------------------

device = hoomd.device.CPU(notice_level=2)
# device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)

# Option parser
parser = OptionParser()
parser.add_option("-n","--resolution"     ,type=int,   dest="resolution"  ,default=100)
parser.add_option("-S","--initgsd"        ,type=str,   dest="initgsd"     ,default=None)
parser.add_option("-i","--steps"          ,type=int,   dest="steps"       ,default=20000)
parser.add_option("-R","--reynolds"       ,type=float, dest="reynolds"    ,default=1)
(options, args) = parser.parse_args()

# Fluid and particle properties
num_length          = options.resolution                        # [ - ]
lref                = 0.01                                       # [ m ]
voxelsize           = lref/float(2*num_length)                    # [ m ]
dx                  = voxelsize                                 # [ m ]
specific_volume     = dx * dx * dx                              # [ m^3 ]
rho0                = 1000.0                                       # [ kg/m^3 ]
mass1               = rho0 * specific_volume                    # [ kg ]
mass2               = rho0 * specific_volume * 1.05                    # [ kg ]
fx                  = 1.0                                       # [ m/s ]
viscosity           = 1e-03   # [ Pa s ]

# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx                  # [ m ]
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength                # [ m ]

# particles per Kernel Radius
part_rcut  = math.ceil(rcut/dx) 
part_depth = math.ceil(1.5 * hoomd.sph.kernel.Kappa[kernel] * rcut/dx) 

# get simulation box sizes etc.
nx, ny, nz = int(2*num_length + (2*part_rcut)), int(0.5*num_length), int(part_depth)
lx, ly, lz = float(nx) * voxelsize, float(ny) * voxelsize, float(nz) * voxelsize
# box dimensions
box_lx, box_ly, box_lz = lx, ly, lz

# Number of Particles
n_particles = nx * ny * nz 

# define meshgrid and add properties
x, y, z = np.meshgrid(*(np.linspace(-box_lx / 2 + (dx/2), box_lx / 2 - (dx/2), nx, endpoint=True),),
                      *(np.linspace(-box_ly / 2 + (dx/2), box_ly / 2 - (dx/2), ny, endpoint=True),),
                      *(np.linspace(-box_lz / 2 + (dx/2), box_lz / 2 - (dx/2), nz, endpoint=True),))

positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
masses     = np.ones((positions.shape[0]), dtype = np.float32) * mass1
slengths   = np.ones((positions.shape[0]), dtype = np.float32) * slength
densities  = np.ones((positions.shape[0]), dtype = np.float32) * rho0

# create Snapshot 
snapshot = gsd.hoomd.Frame()
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
masses = snapshot.particles.mass[:]
for i in range(len(x)):
    xi,yi,zi  = x[i][0], x[i][1], x[i][2]
    tid[i]    = 0
    if ( xi > 0.3 * lref and xi <= 0.5 * lref):
        tid[i] = 0
        masses[i] = mass2
    # solid walls 
    if ( xi < -0.5 * lref or xi > 0.5 * lref):
        tid[i] = 1
        masses[i] = mass1

snapshot.particles.typeid[:]     = tid
snapshot.particles.velocity[:]   = vels
snapshot.particles.mass[:]       = masses

sim.create_state_from_snapshot(snapshot)

init_filename = f'dd_{nx}_{ny}_{nz}_vs_{voxelsize}_re_{options.reynolds}_init.gsd'
# hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

with gsd.hoomd.open(name = init_filename, mode = 'w') as f:
    f.append(snapshot)

if device.communicator.rank == 0:
    export_gsd2vtu.export_spf(init_filename)
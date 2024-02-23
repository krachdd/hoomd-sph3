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
import sys, os
import array

import gsd.hoomd
# ------------------------------------------------------------

# device = hoomd.device.CPU(notice_level=2)
device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)

# get stuff from input file
infile = str(sys.argv[1])
params = read_input_fromtxt.get_input_data_from_file(infile)
if device.communicator.rank == 0:
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


# define meshgrid and add properties
x, y, z = np.meshgrid(*(np.linspace(-box_lx / 2, box_lx / 2, nx, endpoint=True),),
                      *(np.linspace(-box_ly / 2, box_ly / 2, ny, endpoint=True),),
                      *(np.linspace(-box_lz / 2, box_lz / 2, nz, endpoint=True),))

positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
masses     = np.ones((positions.shape[0]), dtype = np.float32) * mass
slengths   = np.ones((positions.shape[0]), dtype = np.float32) * slength
densities  = np.ones((positions.shape[0]), dtype = np.float32) * rho0

# create Snapshot 
# snapshot = hoomd.Snapshot(device.communicator)
snapshot = gsd.hoomd.Frame()
snapshot.configuration.box     = [box_lx, box_ly, box_lz] + [0, 0, 0]
snapshot.particles.N           = n_particles
snapshot.particles.position    = positions
snapshot.particles.typeid      = tids
snapshot.particles.types       = ['F', 'S']
snapshot.particles.velocity    = velocities
snapshot.particles.mass        = masses
snapshot.particles.slength     = slengths
snapshot.particles.density     = densities


sim.create_state_from_snapshot(snapshot)

# deletesolid_flag = params['delete_flag']
# if deletesolid_flag == 1:
#     print(f'Delete solid particles')
#     sim, ndel_particles = delete_solids_initial_timestep.delete_solids(sim, device, kernel, 0.000001, viscosity, dx, rho0)
#     n_particles = n_particles - ndel_particles

init_filename = rawfile.replace('.raw', '_init.gsd')
hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

print(f'Filename: {init_filename}, Number of particles: {n_particles}')

# if device.communicator.rank == 0:
#     export_gsd2vtu.export_spf(init_filename)
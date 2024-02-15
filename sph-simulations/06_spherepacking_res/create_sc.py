#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""

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

import gsd.hoomd

device = hoomd.device.CPU(notice_level=2)
# device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)

# get stuff from input file
infile = str(sys.argv[1])
params = read_input_fromtxt.get_input_data_from_file(infile)
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
# add densities
# for i in range(len(dpes)): dpes[i][0] = RHO0

# create Snapshot 
snapshot = gsd.hoomd.Snapshot()
snapshot.configuration.box = [box_Lx, box_Ly, box_Lz] + [0, 0, 0]
snapshot.particles.N = N_particles
snapshot.particles.position = positions
snapshot.particles.typeid = tids
snapshot.particles.types = ['F', 'S']
snapshot.particles.velocity = velocities
snapshot.particles.mass = masses
snapshot.particles.slength = slengths
snapshot.particles.density = density

sim.create_state_from_snapshot(snapshot)

deletesolid_flag = params['delete_flag']
if deletesolid_flag == 1:
    print(f'Delete solid particles')
    sim, ndel_particles = delete_solids_initial_timestep.delete_solids(sim, device, KERNEL, 0.000001, MU, DX, RHO0)
    N_particles = N_particles - ndel_particles

init_filename = rawfile.replace('.raw', f'vs_{voxelsize}_init.gsd')
# hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

with gsd.hoomd.open(name = init_filename, mode = 'w') as f:
    f.append(snapshot)

print(f'Filename: {init_filename}, Number of particles: {N_particles}')

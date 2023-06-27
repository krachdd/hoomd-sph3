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
num_length          = int(sys.argv[1])
lref                = 0.001               # [m]
voxelsize           = lref/num_length
dx                  = voxelsize
specific_volume     = dx * dx * dx
rho0                = 1000.0
mass                = rho0 * specific_volume
print('mass:', mass)
print('dx:', dx)
# refvel              = 10
fx                  = 0.1                # [m/s]
viscosity           = 0.01               # [Pa s]


# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx       # m
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength     # m
print('slength:', slength)
print('rcut:', rcut)


# particles per Kernel Radius
part_rcut  = math.ceil(rcut/dx) 
# part_depth = math.ceil(2.5 * hoomd.sph.kernel.Kappa[kernel] * rcut/dx) 

# get simulation box sizes etc.
nx, ny, nz = int(num_length), int(num_length + (3*part_rcut)), int(num_length + (3*part_rcut))
lx, ly, lz = float(nx) * voxelsize, float(ny) * voxelsize, float(nz) * voxelsize
# box dimensions
box_lx, box_ly, box_lz = lx, ly, lz

# Number of Particles
n_particles = nx * ny * nz 

# get simulation box sizes etc.
nxin, nyin, nzin = int(num_length), int(num_length), int(num_length)
lxin, lyin, lzin = float(nxin) * voxelsize, float(nyin) * voxelsize, float(nzin) * voxelsize
n_initialparticles = nxin * nyin * nzin
# box dimensions
box_lxin, box_lyin, box_lzin = lxin, lyin, lzin

# define meshgrid and add properties
x, y, z = np.meshgrid(*(np.linspace(-box_lx / 2 + dx/2, box_lx / 2 - dx/2, nx, endpoint=True),),
                      *(np.linspace(-box_ly / 2 + dx/2, box_ly / 2 - dx/2, ny, endpoint=True),),
                      *(np.linspace(-box_lz / 2 + dx/2, box_lz / 2 - dx/2, nz, endpoint=True),))

positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

# # define meshgrid and add properties
# xin, yin, zin = np.meshgrid(*(np.linspace(-box_lxin / 2  + dx/2, box_lxin / 2 - dx/2 - dx/2, nxin, endpoint=True),),
#                       *(np.linspace(-box_lyin / 2  + dx/2, box_lyin / 2 - dx/2, nyin, endpoint=True),),
#                       *(np.linspace(-box_lzin / 2  + dx/2, box_lzin / 2 - dx/2, nzin, endpoint=True),))

# define meshgrid and add properties
xin, yin, zin = np.meshgrid(*(np.linspace(-box_lxin / 2 + dx/2, box_lxin / 2 - dx/2, nxin, endpoint=True),),
                      *(np.linspace(-box_lyin / 2 + dx/2, box_lyin / 2 - dx/2, nyin, endpoint=True),),
                      *(np.linspace(-box_lzin / 2 + dx/2, box_lzin / 2 - dx/2, nzin, endpoint=True),))

positionsin = np.array((xin.ravel(), yin.ravel(), zin.ravel())).T

velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
masses     = np.ones((positions.shape[0]), dtype = np.float32) * mass
slengths   = np.ones((positions.shape[0]), dtype = np.float32) * slength
densities  = np.ones((positions.shape[0]), dtype = np.float32) * rho0

velocitiesin = np.zeros((positionsin.shape[0], positionsin.shape[1]), dtype = np.float32)
massesin     = np.ones((positionsin.shape[0]), dtype = np.float32) * mass
slengthsin   = np.ones((positionsin.shape[0]), dtype = np.float32) * slength
densitiesin  = np.ones((positionsin.shape[0]), dtype = np.float32) * rho0

#n_particles.append(n_initialparticles)
positions = np.concatenate((positions,positionsin))
velocities = np.concatenate((velocities,velocitiesin))
masses = np.concatenate((masses,massesin))
slengths = np.concatenate((slengths,slengthsin))
densities = np.concatenate((densities,densitiesin))



# # # create Snapshot 
snapshot = gsd.hoomd.Frame()
snapshot.configuration.box     = [box_lx, box_ly, box_lz] + [0, 0, 0]
snapshot.particles.N           = n_particles  + n_initialparticles
snapshot.particles.position    = positions
snapshot.particles.typeid      = [0] * (n_particles + n_initialparticles)
snapshot.particles.types       = ['F', 'S', 'I']
snapshot.particles.velocity    = velocities
snapshot.particles.mass        = masses
snapshot.particles.slength     = slengths
snapshot.particles.density     = densities


x    = snapshot.particles.position[:n_particles]
tid  = snapshot.particles.typeid[:]
vels = snapshot.particles.velocity[:]
for i in range(len(x)):
    xi,yi,zi  = x[i][0], x[i][1], x[i][2]
    tid[i]    = 0
    # solid walls 
    if ( yi < -0.5 * lref or yi > 0.5 * lref):
        tid[i] = 1
    if ( zi < -0.5 * lref or zi > 0.5 * lref):
        tid[i] = 1

xin    = snapshot.particles.position[n_particles:]
tidin  = snapshot.particles.typeid[n_particles:]
for i in range(len(xin)):
    xi,yi,zi  = xin[i][0], xin[i][1], xin[i][2]
    tidin[i]    = 2

snapshot.particles.typeid[:n_particles]     = tid
snapshot.particles.typeid[n_particles:]     = tidin
snapshot.particles.velocity[:]   = vels

sim.create_state_from_snapshot(snapshot)

init_filename = f'channel_flow_{nx}_{ny}_{nz}_vs_{voxelsize}_init.gsd'
# hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

with gsd.hoomd.open(name = init_filename, mode = 'wb') as f:
    f.append(snapshot)

if device.communicator.rank == 0:
    export_gsd2vtu.export_spf(init_filename)
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
import os, sys

import gsd.hoomd
# ------------------------------------------------------------

device = hoomd.device.CPU(notice_level=2)
# device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)

# Fluid and particle properties
num_length          = int(sys.argv[1])

# Inner diameter and radius of cylinder
innerD = 0.1                              # [m]
innerR = 0.5*innerD                       # [m]

# Mix features
fw  = innerR*0.1                          # [m]
fh  = innerR*0.9                          # [m]

# Angular velocity in rad/s and 1/s
angvel_s = 1.0                           # 1/s
angvel   = angvel_s*(2*np.pi)            # rad/s

# Characteristic length and velocity
lref = innerR - fh                        # [m]
refvel = angvel*innerD

# Discretization parameters
voxelsize           = innerD/float(num_length)
dx                  = voxelsize
specific_volume     = dx * dx * dx
rho0                = 1000.0              # [kg / m^3]
mass                = rho0 * specific_volume
viscosity           = 1.0            # [Pa s]

# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx       # m
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength     # m

# particles per Kernel Radius
part_rcut  = math.ceil(rcut/dx) 

# outer diameter of cylinder
outerD = innerD + 2 * part_rcut * voxelsize                # m
outerR = 0.5*outerD

# water free surface
frsurf = 0.5*outerR      # m

# get simulation box sizes etc.
nx, ny, nz = int(num_length + (3*part_rcut)), int(num_length + (3*part_rcut)), int(num_length)
lx, ly, lz = float(nx) * voxelsize, float(ny) * voxelsize, float(nz) * voxelsize
# box dimensions
box_lx, box_ly, box_lz = lx, ly, lz

# Number of Particles
n_particles = nx * ny * nz 

# define meshgrid and add properties
x, y, z = np.meshgrid(*(np.linspace(-box_lx / 2 + dx/2, box_lx / 2 - dx/2, nx, endpoint=True),),
                      *(np.linspace(-box_ly / 2 + dx/2, box_ly / 2 - dx/2, ny, endpoint=True),),
                      *(np.linspace(-box_lz / 2 + dx/2, box_lz / 2 - dx/2, nz, endpoint=True),))

positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
masses     = np.ones((positions.shape[0]), dtype = np.float32) * mass
slengths   = np.ones((positions.shape[0]), dtype = np.float32) * slength
densities  = np.ones((positions.shape[0]), dtype = np.float32) * rho0

# create Snapshot 
snapshot = gsd.hoomd.Frame()
snapshot.configuration.box     = [box_lx, box_ly, box_lz] + [0, 0, 0]
snapshot.particles.N           = n_particles
snapshot.particles.position    = positions
snapshot.particles.typeid      = [0] * n_particles
snapshot.particles.types       = ['F', 'S', 'R', 'D']
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
    if (np.sqrt((xi)**2 + (yi)**2) > innerR):
        tid[i] = 1
    # inner mixer
    if (xi**2/fh**2 + yi**2/(0.4*fh)**2 < 1):
        #print('HelloAgainMixer')
        tid[i] = 2
    # # fluid particles to be deleted
    # if ( tid[i] == 0 and yi > frsurf ):
    #     tid[i] = 3
    # delete redundant solid particles
    if (np.sqrt((xi)**2 + (yi)**2) > outerR):
        tid[i] = 3

snapshot.particles.typeid[:]     = tid
snapshot.particles.velocity[:]   = vels

sim.create_state_from_snapshot(snapshot)

# Identify solid particles with zero charge and delete them ( redundant solid particles )
tags    = []
deleted = 0

with sim.state.cpu_local_snapshot as snap:
    for i in range(len(snap.particles.position)):
        # print(data.particles.mass[i])
        if snap.particles.typeid[i] == 3:
            tags.append(snap.particles.tag[i])
            #print('tag:',snap.particles.tag[i])
            # print(f'Rank: {device.communicator.rank} -> Delete Particle {data.particles.tag[i]}')
            deleted += 1

# if device.communicator.rank == 0:
for t in tags:
    # print(f'Rank: {device.communicator.rank} --> Remove particle {t} of {deleted}')
    sim.state.removeParticle(t)

init_filename = f'mixer_{nx}_{ny}_{nz}_vs_{voxelsize}_init.gsd'
hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

# with gsd.hoomd.open(name = init_filename, mode = 'wb') as f:
#     f.append(snapshot)

if device.communicator.rank == 0:
    export_gsd2vtu.export_spf(init_filename)
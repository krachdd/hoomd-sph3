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
import delete_initialfluids_initial_timestep
import sys, os

import gsd.hoomd

# ------------------------------------------------------------


device = hoomd.device.CPU(notice_level=2)
# device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)

# Fluid and particle properties
D                   = 100
lref                = 0.1               # [m]
voxelsize           = lref/D
dx                  = voxelsize
specific_volume     = dx * dx * dx
rho0                = 1000.0
mass                = rho0 * specific_volume
radius              = 0.02
fx                  = 1.5e-07                # [m/s]
viscosity           = 0.001               # [Pa s] 


# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx       # m
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength     # m

# particles per Kernel Radius
part_rcut  = math.ceil(rcut/dx) 
part_depth = math.ceil(2.5 * hoomd.sph.kernel.Kappa[kernel] * rcut/dx) 

# get simulation box sizes etc.
nx, ny, nz = int(D), int(D), int(part_depth)
lx, ly, lz = float(nx) * voxelsize, float(ny) * voxelsize, float(nz) * voxelsize
# box dimensions
box_lx, box_ly, box_lz = lx, ly, lz



# Number of Particles
n_particles = nx * ny * nz 
n_initialparticles = nx * ny * nz

print('Number of particles:', 2*n_particles)

# define meshgrid and add properties
x, y, z = np.meshgrid(*(np.linspace(-box_lx / 2 + dx/2, box_lx / 2 - dx/2, nx, endpoint=True),),
                      *(np.linspace(-box_ly / 2 + dx/2, box_ly / 2 - dx/2, ny, endpoint=True),),
                      *(np.linspace(-box_lz / 2 + dx/2, box_lz / 2 - dx/2, nz, endpoint=True),))

positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

# define meshgrid and add properties
xin, yin, zin = np.meshgrid(*(np.linspace(-box_lx / 2 + dx/2, box_lx / 2 - dx/2, nx, endpoint=True),),
                      *(np.linspace(-box_ly / 2 + dx/2, box_ly / 2 - dx/2, ny, endpoint=True),),
                      *(np.linspace(-box_lz / 2 + dx/2, box_lz / 2 - dx/2, nz, endpoint=True),))

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
    # if ( yi < -2.5 * lref or yi > 2.5 * lref):
    #     tid[i] = 1
    centerx = 0.0 
    centery = 0.0
    distance = np.sqrt((yi - centery)**2 + (xi - centerx)**2)
    if (distance < radius):
        tid[i] = 1

xin    = snapshot.particles.position[n_particles:]
tidin  = snapshot.particles.typeid[n_particles:]
for i in range(len(xin)):
    xi,yi,zi  = xin[i][0], xin[i][1], xin[i][2]
    tidin[i]    = 2

snapshot.particles.typeid[:n_particles]     = tid
snapshot.particles.typeid[n_particles:]     = tidin
snapshot.particles.velocity[:]   = vels

#print(snapshot.particles.typeid[n_particles:])

sim.create_state_from_snapshot(snapshot)

N_particles = n_particles + n_initialparticles

# deletesolid_flag = params['delete_flag']
# if deletesolid_flag == 1:
#     print(f'Delete solid particles')
#     sim, ndel_particles = delete_solids_initial_timestep.delete_solids(sim, device, KERNEL, 0.000001, MU, DX, RHO0)
#     N_particles = N_particles - ndel_particles

if device.communicator.rank == 0:
    print(f'Delete initial fluid particles')
    sim, ndel_particles = delete_initialfluids_initial_timestep.delete_initialfluids(sim, device, kernel, 0.000001, viscosity, dx, rho0)
    N_particles = N_particles - ndel_particles

# if device.communicator.rank == 0:
#     print(f'Delete solid particles')
#     sim, ndel_particles = delete_solids_initial_timestep.delete_solids(sim, device, kernel, 0.000001, viscosity, dx, rho0)
#     N_particles = N_particles - ndel_particles


init_filename = f'cylinder_body_force_{nx}_{ny}_{nz}_vs_{voxelsize}_init.gsd'
hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

# with gsd.hoomd.open(name = init_filename, mode = 'wb') as f:
#     f.append(snapshot)

print(f'Filename: {init_filename}, Number of particles: {N_particles}')

if device.communicator.rank == 0:
    export_gsd2vtu.export_spf(init_filename)
#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: drostan, daniel.rostan@mib.uni-stuttgart.de
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

# Discretization parameters
num_length          = int(sys.argv[1])
lref                = 0.01125              # m
voxelsize           = lref/num_length
dx                  = voxelsize
specific_volume     = dx * dx * dx
rho0                = 1920.0              # [kg / m^3]
mass                = rho0 * specific_volume
viscosity0          = 0.5               # [Pa s]
shearstress0        = 40.0
m_model             = 50

# Inner diameter and radius
innerD = 0.2667                              # [m]
innerR = 0.5*innerD
# Outer diameter and radius
outerD = 0.2892                              # [m]
outerR = 0.5*outerD

# Angular velocity in rad/s and 1/s
angvel_s = 1.3333                           # 1/s
angvel   = angvel_s*(2*np.pi)            # rad/s

# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx       # m
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength     # m

# types 
types = ['F', 'S', 'R', 'D']
#types = ['F', 'S', 'R']


# particles per Kernel Radius
part_rcut  = math.ceil(rcut/dx) 
print('part_rcut: ',part_rcut)

L_all = 0.664
L = 0.40
L_bar = 0.133
L_barbar = 0.131
alpha = (2/3)*np.pi
Ly = L_all + 4*rcut
Lx = outerD + 3*rcut
Lz = outerD + 3*rcut
print('Lx: ',Lx)
print('Ly: ',Ly)
print('Lz: ',Lz)
n_x = Lx/dx
n_y = Ly/dx
n_z = Lz/dx
print('n_x: ',n_x)
print('n_y: ',n_y)
print('n_z: ',n_z)

Rs = 0.05
# alpha_half = 120.0*0.5
alpha_half = 0.5*alpha
L_tri = innerR/np.tan(alpha_half)
print('L_tri: ', L_tri)

part_x = math.ceil(Lx/dx) 
part_y = math.ceil(Ly/dx) 
part_z = math.ceil(Lz/dx) 


print(part_x)
print(part_y)
print(part_z)
print(dx)


# Lx = innerD/voxelsize + 3*part_rcut
# Ly = lref/voxelsize + 3*part_rcut
# Lz = innerD/voxelsize + 3*part_rcut
# outerR = innerR + 2.0*part_rcut * voxelsize


# get simulation box sizes etc.
# nx, ny, nz = int(num_length + (3*part_rcut)), int(num_length + (3*part_rcut)), int(0.2*num_length + (3*part_rcut))
nx, ny, nz = int(part_x), int(part_y), int(part_z)
lx, ly, lz = float(nx) * voxelsize, float(ny) * voxelsize, float(nz) * voxelsize
# box dimensions
box_lx, box_ly, box_lz = lx, ly, lz

# Number of Particles
n_particles = nx * ny * nz 
print('n_particles', n_particles)

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
snapshot.particles.types       = types
snapshot.particles.velocity    = velocities
snapshot.particles.mass        = masses
snapshot.particles.slength     = slengths
snapshot.particles.density     = densities

print('test:',-(0.5*L_all-L_bar))
print('test2:', innerR/np.tan(alpha_half))

x    = snapshot.particles.position[:]
tid  = snapshot.particles.typeid[:]
vels = snapshot.particles.velocity[:]
for i in range(len(x)):
    xi,yi,zi  = x[i][0], x[i][1], x[i][2]
    tid[i]    = 0
    # solid outside boundary walls 
    if (np.sqrt((xi)**2 + (zi)**2) > outerR):
        tid[i] = 1
    # rotating middle 
    if (np.sqrt((xi)**2 + (zi)**2) < innerR):
        tid[i] = 2
    # fluid L_barbar
    if (yi < 0.5 * L_all) and (yi > (0.5*L_all - L_barbar)) and (np.sqrt((xi)**2 + (zi)**2) > Rs):
        tid[i] = 0
    # fluid L_bar (0.13335 = innerR)
    if (yi > -0.5 * L_all) and (yi < -(0.5*L_all - L_barbar)) and (np.sqrt((xi)**2 + (zi)**2) > (0.13335-np.tan(alpha_half)*(-yi-(0.5*L_all-L_bar)))):
        tid[i] = 0
    # lower plate
    if (yi < -0.5 * L_all):
        tid[i] = 1
    # solid outside boundary walls 
    if (np.sqrt((xi)**2 + (zi)**2) > outerR):
        tid[i] = 1
    # particles to delete free surface
    if (np.sqrt((xi)**2 + (zi)**2) < outerR) and (yi > 0.5*L_all) and (np.sqrt((xi)**2 + (zi)**2) > Rs):
        tid[i] = 3
    # particles to delete L_barbar
    # if (yi < 0.5 * L_all) and (yi > (0.5*L_all - L_barbar)) and (np.sqrt((xi)**2 + (zi)**2) < (Rs-2*rcut)):
    if (yi < 0.5 * L_all) and (yi > 0.0) and (np.sqrt((xi)**2 + (zi)**2) < (Rs-2*rcut)):
        tid[i] = 3
    # particles to delete L
    if (yi < (0.5*L_all - L_barbar - 2*rcut)) and (yi > (-0.5*L_all + L_bar + 2*rcut)) and (np.sqrt((xi)**2 + (zi)**2) < (innerR-2*rcut)):
        tid[i] = 3
    # particles to delete L_bar
    if (yi < (-0.5*L_all + L_bar + 2*rcut)) and (yi > (-0.5*L_all + L_bar + 2*rcut - innerR/np.tan(alpha_half))) and (np.sqrt((xi)**2 + (zi)**2) < (0.13335-rcut-np.tan(alpha_half)*(-yi-(0.5*L_all-L_bar-2*rcut)))):
        tid[i] = 3
    # if (yi < (-0.5*L_all + L_bar + 2*rcut)) and (yi > (-0.5*L_all + L_bar + 2*rcut - innerR/np.tan(alpha_half))) and (np.sqrt((xi)**2 + (zi)**2) < (np.tan(alpha_half)*(-yi-(0.5*L_all-L_bar-2*rcut)))):
    #     tid[i] = 3
    # if (yi < (-0.5*L_all + L_bar + 2*rcut)) and (yi > (-0.5*L_all + L_bar + 2*rcut + 0.5*innerR*np.tan(alpha_half))) and (np.sqrt((xi)**2 + (zi)**2) < Rs):
    #     tid[i] = 0
    # if (yi > -0.5 * L_all) and (yi < -(0.5*L_all - L_barbar)) and (np.sqrt((xi)**2 + (zi)**2) > (0.13335-np.tan(alpha_half)*(-yi-(0.5*L_all-L_bar)))):
    #     tid[i] = 0
    # # particles to delete L_bar
    # if (yi < 0.5 * L_all) and (yi > (0.5*L_all - L_barbar)) and (np.sqrt((xi)**2 + (zi)**2) < (Rs-2*rcut)):
    #     tid[i] = 3
    # particles to delete out of cylinder
    if (np.sqrt((xi)**2 + (zi)**2) > outerR + 2*rcut):
        tid[i] = 3

 #(0.1335-(np.tan(alpha_half)*((-1)*yi-(0.5*L_all-L_bar)))
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

print('deleted: ', deleted)

# print(f'Delete solid particles')
# sim, ndel_particles = delete_solids_initial_timestep.delete_solids_rigid(sim, device, kernel, 0.000001, viscosity0, dx, rho0)
# n_particles = n_particles - ndel_particles

init_filename = f'pp1_rheometer_closed_ramp_{nx}_{ny}_{nz}_vs_{voxelsize}_init.gsd'
hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

if device.communicator.rank == 0:
    export_gsd2vtu.export_spf(init_filename)
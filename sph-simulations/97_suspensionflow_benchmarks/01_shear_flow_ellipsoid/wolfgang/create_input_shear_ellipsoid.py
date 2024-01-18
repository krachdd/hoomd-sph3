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
lref                = 0.01               # [m]
voxelsize           = lref/float(num_length)
dx                  = voxelsize
specific_volume     = dx * dx * dx
rho0                = 1000.0
rho1                = 1000.0
mass                = rho0 * specific_volume
#fx                  = 0.1                # [m/s]
viscosity           = 0.01               # [Pa s]
lidvel              = 0.001


# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx       # m
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength     # m

# particles per Kernel Radius
part_rcut  = math.ceil(rcut/dx) 
part_depth = math.ceil(2.5 * hoomd.sph.kernel.Kappa[kernel] * rcut/dx) 

# get simulation box sizes etc.
nx, ny, nz = int(num_length), int(num_length + (2*part_rcut)), int(num_length)
lx, ly, lz = float(nx) * voxelsize, float(ny) * voxelsize, float(nz) * voxelsize
# box dimensions
box_lx, box_ly, box_lz = lx, ly, lz

# Number of Particles
n_particles = nx * ny * nz 

# ## Position and properites sphere
# rs1  = lref*0.2
# xs1  = 0.0*lref
# ys1  = 0.0*lref
# zs1  = 0.0*lref
# rhos = rho1

# Data ellipsoid 1
# a1   = lref*0.17
# b1   = lref*0.2
# c1   = lref*0.1
a1   = lref*0.1
b1   = lref*0.33
c1   = lref*0.1
xs1  = 0.0*lref
ys1  = 0.0*lref
zs1  = 0.0*lref
rhos = rho1


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
snapshot.particles.types       = ['F', 'S', 'O']
snapshot.particles.velocity    = velocities
snapshot.particles.mass        = masses
snapshot.particles.slength     = slengths
snapshot.particles.density     = densities



x    = snapshot.particles.position[:]
tid  = snapshot.particles.typeid[:]
vels = snapshot.particles.velocity[:]
dens = snapshot.particles.density[:]
m    = snapshot.particles.mass[:]
for i in range(len(x)):
    xi,yi,zi  = x[i][0], x[i][1], x[i][2]
    tid[i]    = 0
    # solid walls 
    if ( yi < -0.5 * lref or yi > 0.5 * lref):
        tid[i] = 1
    if (yi > 0.5 * lref):
        vels[i][0] = lidvel
    if (yi < -0.5 * lref):
        vels[i][0] = -lidvel
    # # Suspended Object
    # if ( (xi-xs1)**2 + (yi-ys1)**2 + (zi-zs1)**2 < rs1**2 ):
    #     tid[i] = 2
    #     dens[i]= rhos
    #     m[i]   = specific_volume*rhos
    # Suspended Ellipsoid 1
    if ( ((xi-xs1)**2/a1**2) + ((yi-ys1)**2/b1**2) + ((zi-zs1)**2/c1**2) <= 1.0 ):
        tid[i]    = 2
        dens[i]   = rhos
        m[i]      = specific_volume*rhos



snapshot.particles.typeid[:]     = tid
snapshot.particles.velocity[:]   = vels
snapshot.particles.density[:]    = dens
snapshot.particles.mass[:]       = m



sim.create_state_from_snapshot(snapshot)

init_filename = f'shear_ellispoid_wolfgang_{a1}_{b1}_{c1}_{nx}_{ny}_{nz}_vs_{voxelsize}_init.gsd'
# hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

with gsd.hoomd.open(name = init_filename, mode = 'wb') as f:
    f.append(snapshot)

if device.communicator.rank == 0:
    export_gsd2vtu.export_spf(init_filename)
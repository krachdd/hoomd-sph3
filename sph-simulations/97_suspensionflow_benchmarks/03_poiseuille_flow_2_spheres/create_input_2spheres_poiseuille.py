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

# Fluid and particle properties
# Res old file 60
num_length          = int(sys.argv[1])
lref                = 0.01               # [m]
voxelsize           = lref/float(num_length)
dx                  = voxelsize
specific_volume     = dx * dx * dx
rho0                = 1000.0
rho1                = 1000.0
mass                = rho0 * specific_volume
fx                  = 1.0                # [m/s]
viscosity           = 0.01               # [Pa s]
lidvel              = 0.0
re                  = 10.0
kc                  = 1.0
dc                  = 1.0 

print('dx: ', dx)

# Reference velocity
if lidvel == 0.0:
    viscosity   = 0.01                  # Pa s
    refvel = fx*lref*lref*0.25/(viscosity/rho0)
else:
    refvel = lidvel
    viscosity   = (lref*lidvel*rho0)/re       # Pa s  re = rho0*lidvel*lref/viscosity 

# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx       # m
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength     # m

print('slength:', slength)
print('rcut: ', rcut)

# particles per Kernel Radius
part_rcut  = math.ceil(rcut/dx) 
part_depth = math.ceil(2.5 * hoomd.sph.kernel.Kappa[kernel] * rcut/dx) 

print('part_rcut:', part_rcut)
print('part_depth:', part_depth)


# get simulation box sizes etc.
nx, ny, nz = int(num_length), int(num_length + (2*part_rcut)), int(num_length)
lx, ly, lz = float(nx) * voxelsize, float(ny) * voxelsize, float(nz) * voxelsize
print('lx: ', lx)
print('ly: ', ly)
print('lz: ', lz)

# box dimensions
box_lx, box_ly, box_lz = lx, ly, lz

# Number of Particles
n_particles = nx * ny * nz 


# Data sphere 1
rs1  = lref*0.1
xs1  = 0.0*lref
ys1  = -0.4*lref
zs1  = 0.0*lref
rho_s1 = rho1

# Data sphere 2
rs2  = lref*0.1
xs2  = 0.0*lref
ys2  = 0.4*lref
zs2  = 0.0*lref
rho_s2 = rho1

# Data ellipsoid 1
a1   = lref*0.17
b1   = lref*0.1
c1   = lref*0.2
xs3  = 0.3*lref
ys3  = -0.2*lref
zs3  = 0.05*lref
rho_e1 = rho1

# Data ellipsoid 1
a2   = lref*0.17
b2   = lref*0.1
c2   = lref*0.2
xs4  = -0.3*lref
ys4  = 0.2*lref
zs4  = -0.05*lref
rho_e2 = rho1

# ## Position and properites sphere
# rs1  = lref*0.2
# xs1  = 0.0*lref
# ys1  = 0.0*lref
# zs1  = 0.0*lref
# rhos = rho1

# # Data ellipsoid 1
# # a1   = lref*0.17
# # b1   = lref*0.2
# # c1   = lref*0.1
# a1   = lref*0.1
# b1   = lref*0.33
# c1   = lref*0.1
# xs1  = 0.0*lref
# ys1  = 0.0*lref
# zs1  = 0.0*lref
# rhos = rho1


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
snapshot.particles.types       = ['F', 'S', 'A', 'B']
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
    # Upper Wall
    if ( yi > 0.5 * lref ):
        tid[i] = 1
        vels[i][0] = lidvel
    # Lower Wall
    if ( yi < -0.5 * lref    ):
        tid[i]  = 1
        #vels[i][0] = - lidvel
    # Suspended Objects
    # Suspended Sphere 1
    if ( (xi-xs1)**2 + (yi-ys1)**2 + (zi-zs1)**2 < rs1**2 ):
        tid[i]    = 2
        dens[i] = rho_s1
        m[i]      = specific_volume*rho_s1
    # Suspended Sphere 2
    if ( (xi-xs2)**2 + (yi-ys2)**2 + (zi-zs2)**2 < rs2**2 ):
        tid[i]    = 3
        dens[i] = rho_s2
        m[i]      = specific_volume*rho_s2
    # # Suspended Ellipsoid 1
    # if ( ((xi-xs3)**2/a1**2) + ((yi-ys3)**2/b1**2) + ((zi-zs3)**2/c1**2) <= 1.0 ):
    #     tid[i]    = 2
    #     dens[i] = rho_e1
    #     m[i]      = specific_volume*rho_e1
    # # Suspended Ellipsoid 2
    # if ( ((xi-xs4)**2/a2**2) + ((yi-ys4)**2/b2**2) + ((zi-zs4)**2/c2**2) <= 1.0 ):
    #     tid[i]    = 3
    #     dens[i] = rho_e2
    #     m[i]      = specific_volume*rho_e2   
    # # solid walls 
    # if ( yi < -0.5 * lref or yi > 0.5 * lref):
    #     tid[i] = 1
    # if (yi > 0.5 * lref):
    #     vels[i][0] = lidvel
    # if (yi < -0.5 * lref):
    #     vels[i][0] = -lidvel
    # # # Suspended Object
    # # if ( (xi-xs1)**2 + (yi-ys1)**2 + (zi-zs1)**2 < rs1**2 ):
    # #     tid[i] = 2
    # #     dens[i]= rhos
    # #     m[i]   = specific_volume*rhos
    # # Suspended Ellipsoid 1
    # if ( ((xi-xs1)**2/a1**2) + ((yi-ys1)**2/b1**2) + ((zi-zs1)**2/c1**2) <= 1.0 ):
    #     tid[i]    = 2
    #     dens[i]   = rhos
    #     m[i]      = specific_volume*rhos



snapshot.particles.typeid[:]     = tid
snapshot.particles.velocity[:]   = vels
snapshot.particles.density[:]    = dens
snapshot.particles.mass[:]       = m



sim.create_state_from_snapshot(snapshot)

init_filename = f'2spheres_poiseuille_{nx}_{ny}_{nz}_vs_{voxelsize}_init.gsd'
#hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

with gsd.hoomd.open(name = init_filename, mode = 'wb') as f:
    f.append(snapshot)

if device.communicator.rank == 0:
    export_gsd2vtu.export_spf(init_filename)
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

import gsd.hoomd
# ------------------------------------------------------------

resolutions = [20, 30, 50, 100]

for i in range(len(resolutions)):

    device = hoomd.device.CPU(notice_level=2)
    # device = hoomd.device.CPU(notice_level=10)
    sim = hoomd.Simulation(device=device)

    # Fluid and particle properties
    num_length          = resolutions[i]
    lref                = 0.1               # [m]
    # radius              = 0.5 * lref
    voxelsize           = lref/num_length
    dx                  = voxelsize
    specific_volume     = dx * dx * dx
    rho0                = 1000.0
    mass                = rho0 * specific_volume
    refvel              = 10
    Re                  = 1
    viscosity           = rho0 * lref * refvel/Re


    # get kernel properties
    kernel  = 'WendlandC4'
    slength = hoomd.sph.kernel.OptimalH[kernel]*dx       # m
    rcut    = hoomd.sph.kernel.Kappa[kernel]*slength     # m

    # particles per Kernel Radius
    part_rcut  = math.ceil(rcut/dx) 
    part_depth = math.ceil(2.5 * hoomd.sph.kernel.Kappa[kernel] * rcut/dx) 

    # get simulation box sizes etc.
    nx, ny, nz = num_length, num_length + (2*part_rcut), part_depth
    lx, ly, lz = nx*voxelsize, ny*voxelsize, nz*voxelsize
    # box dimensions
    box_lx, box_ly, box_lz = lx, ly, lz

    # Number of Particles
    n_particles = nx * ny * nz 

    # define meshgrid and add properties
    x, y, z = np.meshgrid(*(np.linspace(-box_lx / 2, box_lx / 2, nx, endpoint=True),),
                          *(np.linspace(-box_ly / 2, box_ly / 2, ny, endpoint=True),),
                          *(np.linspace(-box_lz / 2, box_lz / 2, nz, endpoint=True),))

    positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

    velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
    masses     = np.ones((positions.shape[0]), dtype = np.float32) * mass
    slengths   = np.ones((positions.shape[0]), dtype = np.float32) * slength
    densities  = np.ones((positions.shape[0]), dtype = np.float32) * rho0

    # # # create Snapshot 
    snapshot = gsd.hoomd.Snapshot()
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
    for i in range(len(x)):
        xi,yi,zi  = x[i][0], x[i][1], x[i][2]
        tid[i]    = 0
        # solid walls 
        if ( yi < -0.5 * lref or yi > 0.5 * lref):
            tid[i] = 1
        if (yi > 0.5 * lref):
            vels[i][0] = 0.1

    snapshot.particles.typeid[:]     = tid
    snapshot.particles.velocity[:]   = vels

    sim.create_state_from_snapshot(snapshot)

    init_filename = f'couette_flow_{nx}_{ny}_{nz}_vs_{voxelsize}_init.gsd'
    # hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = init_filename)

    with gsd.hoomd.open(name = init_filename, mode = 'wb') as f:
        f.append(snapshot)

    if device.communicator.rank == 0:
        export_gsd2vtu.export_spf(init_filename)
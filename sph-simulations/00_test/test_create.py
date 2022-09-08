#!/usr/bin/env python3
import hoomd
from hoomd import *
from hoomd import sph
import numpy as np
import itertools
import gsd.hoomd



device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)


#  ________________________________
#               ---> FX             0.001 m
#  ________________________________  
#

# System sizes
LREF = 0.001                    # m

LX = LREF*2
LY = LREF*0.5
LZ = LREF*0.5

# Parameters
KERNEL  = 'WendlandC4'
NL      = 10                       # INT
FX      = 0.1                      # m/s^2

DX      = LREF/NL                  # m
V       = DX*DX*DX                 # m^3

RHO0 = 1000.0                      # kg / m^3
M    = RHO0*V                      # kg
DRHO = 0.05                        # %
MU   = 0.01                        # Pa s

UREF = FX*LREF*LREF*0.25/(MU/RHO0)

H       = hoomd.sph.kernel.OptimalH[KERNEL]*DX       # m
RCUT    = hoomd.sph.kernel.Kappa[KERNEL]*H           # m


Nx = int(LX/DX)                 # particles per box direction
Ny = int((LY + 3 * RCUT)/DX)    # particles per box direction
Nz = int((LZ + 3 * RCUT)/DX)    # particles per box direction
N_particles = Nx * Ny * Nz      # Number of Particles

box_Lx = LREF*2    # box dimension
box_Ly = LREF*0.5  # box dimension
box_Lz = LREF*0.5  # box dimension

x, y, z = np.meshgrid(*(np.linspace(-box_Lx / 2, box_Lx / 2, Nx, endpoint=False),),
                      *(np.linspace(-box_Ly / 2, box_Ly / 2, Ny, endpoint=False),),
                      *(np.linspace(-box_Lz / 2, box_Lz / 2, Nz, endpoint=False),))

positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
masses     = np.ones((positions.shape[0]), dtype = np.float32) * M
slengths   = np.ones((positions.shape[0]), dtype = np.float32) * H
dpes       = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)

snapshot = hoomd.Snapshot(device.communicator)
snapshot.configuration.box = [box_Lx, box_Ly, box_Lz] + [0, 0, 0]
snapshot.particles.N = N_particles
snapshot.particles.position[:] = positions
snapshot.particles.typeid[:] = [0] * N_particles
snapshot.particles.types = ['F', 'S']
snapshot.particles.velocity[:] = velocities
snapshot.particles.mass[:] = masses
snapshot.particles.slength[:] = slengths
snapshot.particles.dpe[:] = dpes

x   = snapshot.particles.position[:]
dpe = snapshot.particles.dpe[:]
tid = snapshot.particles.typeid[:]

for i in range(len(x)):
    xi,yi,zi  = x[i][0], x[i][1], x[i][2]
    dpe[i][0] = RHO0
    tid[i]    = 0
    if ( yi < -0.4*LY or yi > 0.4*LY ):
        tid[i] = 1

snapshot.particles.dpe[:]      = dpe
snapshot.particles.typeid[:]   = tid

sim.create_state_from_snapshot(snapshot)

hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = "test_parallelplates.gsd")
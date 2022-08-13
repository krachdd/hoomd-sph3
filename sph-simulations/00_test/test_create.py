#!/home/david/anaconda3/envs/sph3/bin/python
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


# Nx = int(LX/DX)
# Ny = int((LY + 3 * RCUT)/DX)
# Nz = int((LZ + 3 * RCUT)/DX)
# N_particles = Nx * Ny * Nz

# xdim = np.linspace(0, LX, int(LX/DX), endpoint=False)
# ydim = np.linspace(0, (LY+3*RCUT), int((LY+3*RCUT)/DX), endpoint=False)
# zdim = np.linspace(0, (LZ+3*RCUT), int((LZ+3*RCUT)/DX), endpoint=False)



# xdim = np.linspace(-LX/2, LX/2, int(LX/DX), endpoint=True)
# ydim = np.linspace(-(LY/2+3*RCUT), (LY/2+3*RCUT), int((LY+3*RCUT)/DX), endpoint=True)
# zdim = np.linspace(-(LZ/2+3*RCUT), (LZ/2+3*RCUT), int((LZ+3*RCUT)/DX), endpoint=True)

# Create a simple cubic configuration of particles
N = 5  # particles per box direction
box_L = 20  # box dimension
N_particles = N**3

x, y, z = np.meshgrid(*(np.linspace(-box_L / 2, box_L / 2, N, endpoint=False),)
                      * 3)

# x, y, z = np.meshgrid(xdim, ydim, zdim)
positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

# print("xmax, min, {0}, {1}".format(np.max(positions[:, 0]), np.min(positions[:, 0])))
# print("ymax, min, {0}, {1}".format(np.max(positions[:, 1]), np.min(positions[:, 1])))
# print("zmax, min, {0}, {1}".format(np.max(positions[:, 2]), np.min(positions[:, 2])))

# print(LX, LY, LZ, RCUT, LY + 3* RCUT, LZ + 3* RCUT)

velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
# velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
masses     = np.ones((positions.shape[0]), dtype = np.float32) * M
slengths   = np.ones((positions.shape[0]), dtype = np.float32) * H
dpes       = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)

snapshot = hoomd.Snapshot(device.communicator)
snapshot.configuration.box = [box_L] * 3 + [0, 0, 0]
# snapshot.configuration.box = [LX*1.00001, (LY + 3 * RCUT)*1.00001, (LZ + 3 * RCUT)* 1.00001, 0, 0, 0]
snapshot.particles.N = N_particles
snapshot.particles.position[:] = positions
snapshot.particles.typeid[:] = [0] * N_particles
snapshot.particles.types = ['F', 'S']
snapshot.particles.velocity[:] = velocities
snapshot.particles.mass[:] = masses
snapshot.particles.slength[:] = slengths
snapshot.particles.dpe[:] = dpes

# snapshot = gsd.hoomd.Snapshot()



# # if device.communicator.rank == 0:
# # m   = snapshot.particles.mass[:]
# # v   = snapshot.particles.velocity[:]
# x   = snapshot.particles.position[:]
# # h   = snapshot.particles.slength[:]
# dpe = snapshot.particles.dpe[:]
# tid = snapshot.particles.typeid[:]

# for i in range(len(x)):
#     xi,yi,zi  = x[i][0], x[i][1], x[i][2]
#     dpe[i][0] = RHO0
#     tid[i]    = 0
#     if ( yi < -0.5*LY or yi > 0.5*LY ):
#         tid[i] = 1

# # snapshot.particles.mass[:]     = m
# # snapshot.particles.slength[:]  = h
# snapshot.particles.dpe[:]      = dpe
# snapshot.particles.typeid[:]   = tid

print(snapshot.configuration.box)
# print(LX, LY, LZ, RCUT, positions)

sim.create_state_from_snapshot(snapshot)

hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = "test_parallelplates.gsd")


# with gsd.hoomd.open(name='lattice.gsd', mode='xb') as f:
#     f.append(snapshot)


# # Print the domain decomposition.
# domain_decomposition = sim.state.domain_decomposition
# if device.communicator.rank == 0:
#     print("Domain decomposition: {0}".format(domain_decomposition))

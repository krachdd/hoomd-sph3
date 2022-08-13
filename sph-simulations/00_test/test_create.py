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
NL      = 20                       # INT
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


Nx = int(LX/DX)
Ny = int((LY + 3 * RCUT)/DX)
Nz = int((LZ + 3 * RCUT)/DX)
N_particles = Nx * Ny * Nz

xdim = np.linspace(-LX/2, LX/2, int(LX/DX), endpoint=True)
ydim = np.linspace(-(LY/2+3*RCUT), (LY/2+3*RCUT), int((LY+3*RCUT)/DX), endpoint=True)
zdim = np.linspace(-(LZ/2+3*RCUT), (LZ/2+3*RCUT), int((LZ+3*RCUT)/DX), endpoint=True)

x, y, z = np.meshgrid(xdim, ydim, zdim)
positions = np.array((x.ravel(), y.ravel(), z.ravel())).T

velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
# velocities = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)
masses     = np.ones((positions.shape[0]), dtype = np.float32) * M
slengths   = np.ones((positions.shape[0]), dtype = np.float32) * H
dpes       = np.zeros((positions.shape[0], positions.shape[1]), dtype = np.float32)

snapshot = hoomd.Snapshot(device.communicator)
snapshot.particles.N = N_particles
snapshot.particles.position[:] = positions
snapshot.particles.typeid[:] = [0] * N_particles
snapshot.particles.types = ['F', 'S']
snapshot.particles.velocity[:] = velocities
snapshot.particles.mass[:] = masses
snapshot.particles.slength[:] = slengths
snapshot.particles.dpe[:] = dpes

# snapshot = gsd.hoomd.Snapshot()



# if device.communicator.rank == 0:
m   = snapshot.particles.mass[:]
v   = snapshot.particles.velocity[:]
x   = snapshot.particles.position[:]
h   = snapshot.particles.slength[:]
dpe = snapshot.particles.dpe[:]
tid = snapshot.particles.typeid[:]

for i in range(len(x)):
    xi,yi,zi  = x[i][0], x[i][1], x[i][2]
    dpe[i][0] = RHO0
    tid[i]    = 0
    if ( yi < -0.5*LY or yi > 0.5*LY ):
        tid[i] = 1

sim.create_state_from_snapshot(snapshot)

hoomd.write.GSD.write(state = sim.state, mode = 'wb', filename = "test_parallelplates.gsd")


# with gsd.hoomd.open(name='lattice.gsd', mode='xb') as f:
#     f.append(snapshot)


# Print the domain decomposition.
domain_decomposition = sim.state.domain_decomposition
if device.communicator.rank == 0:
    print("Domain decomposition: {0}".format(domain_decomposition))

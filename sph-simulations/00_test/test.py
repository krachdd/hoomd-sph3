#!/home/david/anaconda3/envs/sph3/bin/python
import hoomd
from hoomd import *
from hoomd import sph
import numpy as np
import itertools
import gsd.hoomd



device = hoomd.device.CPU()
sim = hoomd.Simulation(device=device)



# sim.operations.integrator = hoomd.sph.Integrator()


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


print("H = {0}".format(H))
print("RCUT = {0}".format(RCUT))
print("int(RCUT/DX) = {0}".format(int(RCUT/DX)))

print("int(LX/DX) = {0}".format(int(LX/DX)))
print("LX = {0}".format(LX))
print("int(LY/DX) = {0}".format(int(LY/DX)))
print("int(LZ/DX) = {0}".format(int(LZ/DX)))

xdim = np.linspace(-LX/2, LX/2, int(LX/DX))
ydim = np.linspace(-(LY/2+RCUT), (LY/2+RCUT), int((LY+RCUT)/DX))
zdim = np.linspace(-(LZ/2+RCUT), (LZ/2+RCUT), int((LZ+RCUT)/DX))

position = list(itertools.product(xdim, ydim, zdim, repeat=1))
N_particles = len(position)

snapshot = gsd.hoomd.Snapshot()
snapshot.particles.N = N_particles
snapshot.particles.position = position[:]

if device.communicator.rank == 0:
    # m   = snapshot.particles.mass[:]
    # v   = snapshot.particles.velocity[:]
    x   = snapshot.particles.position[:]
    # h   = snapshot.particles.slength[:]
    # dpe = snapshot.particles.dpe[:]
    # tid = snapshot.particles.typeid[:]
    # Set initial conditions
    # snapshot.particles.types = ['F','S']






# Print the domain decomposition.
# domain_decomposition = sim.state.domain_decomposition
# if device.communicator.rank == 0:
#     print(domain_decomposition)

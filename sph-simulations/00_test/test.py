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


# print("H = {0}".format(H))
# print("RCUT = {0}".format(RCUT))
# print("int(RCUT/DX) = {0}".format(int(RCUT/DX)))

# print("int(LX/DX) = {0}".format(int(LX/DX)))
# print("LX = {0}".format(LX))
# print("int(LY/DX) = {0}".format(int(LY/DX)))
# print("int(LZ/DX) = {0}".format(int(LZ/DX)))

xdim = np.linspace(-LX/2, LX/2, int(LX/DX), endpoint=True)
ydim = np.linspace(-(LY/2+RCUT), (LY/2+RCUT), int((LY+RCUT)/DX), endpoint=True)
zdim = np.linspace(-(LZ/2+RCUT), (LZ/2+RCUT), int((LZ+RCUT)/DX), endpoint=True)

position = list(itertools.product(xdim, ydim, zdim, repeat=1))



N_particles = len(position)


velocity = []
for i in range(N_particles): velocity.append((np.float64(0.), np.float64(0.), np.float64(0.)))
dpe = np.copy(velocity)
mass = []
for i in range(N_particles): mass.append(np.float64(0.))
slength = np.copy(mass)
typeid = []
for i in range(N_particles): typeid.append(np.int32(0))


# print(type(velocity[0][0]))


snapshot = gsd.hoomd.Snapshot()
snapshot.particles.N = N_particles
snapshot.particles.position = position[:]
# snapshot.particles.velocity = velocity[:]
# # snapshot.particles.mass = mass[:]
# # snapshot.particles.slength = slength[:]
# # snapshot.particles.dpe = dpe[:]
# # snapshot.particles.typeid = typeid[:]
# # snapshot.particles.types = ['F','S']

# if device.communicator.rank == 0:
#     m   = snapshot.particles.mass[:]
#     v   = snapshot.particles.velocity[:]
#     x   = snapshot.particles.position[:]
#     # h   = snapshot.particles.slength[:]
#     # dpe = snapshot.particles.dpe[:]
#     # tid = snapshot.particles.typeid[:]
#     # Set initial conditions
#     # snapshot.particles.types = ['F','S']


#     for i in range(len(x)):
#         xi,yi,zi  = x[i][0], x[i][1], x[i][2]
#         # h[i]      = H
#         # dpe[i][0] = RHO0
#         m[i]      = M
#         v[i][0]   = 0.0
#         v[i][1]   = 0.0
#         v[i][2]   = 0.0
#         # tid[i]    = 0
#         # if ( yi < -0.5*LY or yi > 0.5*LY ):
#         #     tid[i] = 1

#     snapshot.particles.velocity[:] = v
#     snapshot.particles.mass[:]     = m
#     # snapshot.particles.slength[:]  = h
#     # snapshot.particles.dpe[:]      = dpe
#     snapshot.particles.typeid[:]   = tid


with gsd.hoomd.open(name='parallelplates.gsd', mode='xb+') as f:
    f.append(snapshot)


# Print the domain decomposition.
# domain_decomposition = sim.state.domain_decomposition
# if device.communicator.rank == 0:
#     print(domain_decomposition)

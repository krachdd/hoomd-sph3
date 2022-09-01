#!/home/david/anaconda3/envs/sph3/bin/python
import hoomd
from hoomd import *
from hoomd import sph
from hoomd.sph import _sph
import numpy as np
import itertools
import gsd.hoomd


# -----------------------------------------------------------------------------------
# KEEP SAME AS IN TEST CREATE 
# -----------------------------------------------------------------------------------


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

# ------------------------------------------------------------------------------------



# device = hoomd.device.CPU(notice_level=10)
device = hoomd.device.CPU(notice_level=2)
sim = hoomd.Simulation(device=device)
# if device.communicator.rank == 0:
#     print("hoomd.version.mpi_enabled: {0}".format(hoomd.version.mpi_enabled))


filename = "test_parallelplates.gsd"
dumpname = "run_{0}".format(filename)

sim.create_state_from_gsd(filename = filename)

# Print the domain decomposition.
domain_decomposition = sim.state.domain_decomposition
if device.communicator.rank == 0:
    print(domain_decomposition)

# Print the location of the split planes.
split_fractions = sim.state.domain_decomposition_split_fractions
if device.communicator.rank == 0:
    print(split_fractions)

# Print the number of particles on each rank.
with sim.state.cpu_local_snapshot as snap:
    N = len(snap.particles.position)
    print(f'{N} particles on rank {device.communicator.rank}')

# Kernel
KERNEL  = 'WendlandC4'
H       = hoomd.sph.kernel.OptimalH[KERNEL]*DX       # m
RCUT    = hoomd.sph.kernel.Kappa[KERNEL]*H           # m
Kernel = hoomd.sph.kernel.Kernels[KERNEL]()
Kappa  = Kernel.Kappa()

print(Kernel.__dict__)

# print(dir(hoomd.sph._sph))
# print(hoomd.sph.kernel.Kernels['WendlandC4']
# print(dir(hoomd.sph.kernel.WendlandC4))



# # Neighbor list
NList = hoomd.nsearch.nlist.Cell(buffer = RCUT*0.05, rebuild_check_delay = 1, kappa = Kappa)

print(NList.__dict__)

# print(NList._getattr_param('kappa'))

# # print(help(NList._setattr_param))
# print(NList._setattr_param('kappa', 2.1))



# print(NList._getattr_param('kappa'))

# Kernel.setNeighborList(NList)

# # Equation of State
EOS = hoomd.sph.eos.Tait()
EOS.set_params(RHO0,0.05)

# # Define groups
groupFLUID  = hoomd.filter.Type(['F']) # is zero
groupSOLID  = hoomd.filter.Type(['S']) # is one

# print(type(groupSOLID))

with sim.state.cpu_local_snapshot as data:
    print(f'{np.count_nonzero(data.particles.typeid == 0)} fluid particles on rank {device.communicator.rank}')
    print(f'{np.count_nonzero(data.particles.typeid == 1)} solid particles on rank {device.communicator.rank}')

# # Set up SPH solver
# model = hoomd.sph.models.SinglePhaseFlow(Kernel,EOS,NList,groupFLUID,groupSOLID)

integrator = hoomd.sph.Integrator(dt=0.005)

VelocityVerlet = hoomd.sph.methods.VelocityVerlet(filter=groupFLUID)
integrator.methods.append(VelocityVerlet)

model = hoomd.sph.sphmodel.SinglePhaseFlow(kernel = Kernel,
                                           eos    = EOS,
                                           nlist  = NList,
                                           fluidgroup = groupFLUID,
                                           solidgroup = groupSOLID,
                                           #mu     = 0.01
                                           )

integrator.forces.append(model)



print(model)
# model.densitymethod = 'SUMMATION'
# model.set_params(MU)
# model.activateArtificialViscosity(0.2,0.0)
# model.deactivateShepardRenormalization()
# model.deactivateDensityDiffusion()
# model.setBodyAcceleration(FX,0.0,0.0,1000)
# dt = model.compute_dt(LREF,UREF,DRHO)



print("integrator Forces: {0}".format(integrator.forces[:]))
print("integrator Methods: {0}".format(integrator.methods[:]))



hoomd.write.GSD.write(filename = dumpname, state = sim.state, mode = 'wb')
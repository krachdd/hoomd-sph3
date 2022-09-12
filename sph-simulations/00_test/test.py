#!/usr/bin/env python3
import hoomd
from hoomd import *
from hoomd import sph
from hoomd.sph import _sph
import numpy as np
import itertools
from datetime import datetime
import export_gsd2vtu 
# from mpi4py import MPI
# import gsd.hoomd


# -----------------------------------------------------------------------------------
# KEEP SAME AS IN TEST CREATE 
# -----------------------------------------------------------------------------------


# System sizes
LREF = 0.001                    # m

LX = LREF*2
LY = LREF*2
LZ = LREF*2

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

# ------------------------------------------------------------------------------------



device = hoomd.device.CPU(notice_level = 10)
# device = hoomd.device.CPU(notice_level = 7)
# device = hoomd.device.CPU(notice_level = 2)
sim = hoomd.Simulation(device = device)
# if device.communicator.rank == 0:
#     print("hoomd.version.mpi_enabled: {0}".format(hoomd.version.mpi_enabled))


filename = "test_tube.gsd"
dumpname = "run_{0}".format(filename)

sim.create_state_from_gsd(filename = filename)



# Print the domain decomposition.
domain_decomposition = sim.state.domain_decomposition
if device.communicator.rank == 0:
    print(f'Domain Decomposition: {domain_decomposition}')

# Print the location of the split planes.
split_fractions = sim.state.domain_decomposition_split_fractions
if device.communicator.rank == 0:
    print(f'Locations of SplitPlanes: {split_fractions}')

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

# print(Kernel.__dict__)

# print(dir(hoomd.sph._sph))
# print(hoomd.sph.kernel.Kernels['WendlandC4']
# print(dir(hoomd.sph.kernel.WendlandC4))

# # Neighbor list
# NList = hoomd.nsearch.nlist.Cell(buffer = RCUT*0.05, rebuild_check_delay = 1, kappa = Kappa)


# print(NList.__dict__)
# print(NList.loggables)

NList = hoomd.nsearch.nlist.Stencil(buffer = RCUT*0.05, cell_width=1.5, rebuild_check_delay = 1, kappa = Kappa)
# print(NList.__dict__)
# print(NList._getattr_param('kappa'))
# # print(help(NList._setattr_param))
# print(NList._setattr_param('kappa', 2.1))
# print(NList._getattr_param('kappa'))
# Kernel.setNeighborList(NList)



# # Equation of State
EOS = hoomd.sph.eos.Tait()
EOS.set_params(RHO0,0.05)

# Define groups/filters
filterFLUID  = hoomd.filter.Type(['F']) # is zero
filterSOLID  = hoomd.filter.Type(['S']) # is one

with sim.state.cpu_local_snapshot as data:
    print(f'{np.count_nonzero(data.particles.typeid == 0)} fluid particles on rank {device.communicator.rank}')
    print(f'{np.count_nonzero(data.particles.typeid == 1)} solid particles on rank {device.communicator.rank}')

# # Set up SPH solver
model = hoomd.sph.sphmodel.SinglePhaseFlow(kernel = Kernel,
                                           eos    = EOS,
                                           nlist  = NList,
                                           fluidgroup_filter = filterFLUID,
                                           solidgroup_filter = filterSOLID)

print("SetModelParameter")

model.mu = 0.01
model.densitymethod = 'SUMMATION'
model.gx = FX
model.damp = 1000
model.artificialviscosity = True
model.alpha = 0.2
model.beta = 0.0
# Access the local snapshot, this is not ideal! 

print(model.__dict__)

with sim.state.cpu_local_snapshot as snapshot:
    model.max_sl = np.max(snapshot.particles.slength[:])

dt = model.compute_dt(LREF,UREF,DRHO)

integrator = hoomd.sph.Integrator(dt=dt)

VelocityVerlet = hoomd.sph.methods.VelocityVerlet(filter=filterFLUID)

integrator.methods.append(VelocityVerlet)
integrator.forces.append(model)

compute_filter_all = hoomd.filter.All()
spf_properties = hoomd.sph.compute.SinglePhaseFlowBasicProperties(compute_filter_all)
sim.operations.computes.append(spf_properties)

if device.communicator.rank == 0:
    print(f'Computed Time step: {dt}')
    print("Integrator Forces: {0}".format(integrator.forces[:]))
    print("Integrator Methods: {0}".format(integrator.methods[:]))
    print("Simulation Computes: {0}".format(sim.operations.computes[:]))




gsd_trigger = hoomd.trigger.Periodic(100)
gsd_writer = hoomd.write.GSD(filename=dumpname,
                             trigger=gsd_trigger,
                             mode='wb' #,
                             # filter=hoomd.filter.Null()
                             )
sim.operations.writers.append(gsd_writer)



# hoomd.write.GSD.write(filename = dumpname, state = sim.state, mode = 'wb')
log_trigger = hoomd.trigger.Periodic(1)
logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps', 'walltime'])
logger.add(spf_properties, quantities=['kinetic_energy'])
table = hoomd.write.Table(trigger=log_trigger,
                          logger=logger)
# file = open('log.txt', mode='x', newline='\n')
# table_file = hoomd.write.Table(output=file,
#                                trigger=hoomd.trigger.Periodic(period=5000),
#                                logger=logger)
# sim.operations.writers.append(table_file)
sim.operations.writers.append(table)

sim.operations.integrator = integrator

# print(model.loggables)
# print(sim.loggables)
# print(spf_properties.loggables)

if device.communicator.rank == 0:
    print("Starting Run at {0}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

sim.run(1001)

if device.communicator.rank == 0:
    export_gsd2vtu.export_basic(dumpname)

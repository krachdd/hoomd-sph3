#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""

import hoomd
from hoomd import *
from hoomd import sph
from hoomd.sph import _sph
import numpy as np
# import itertools
import datetime as datet
from datetime import datetime
import export_gsd2vtu 
# from mpi4py import MPI
# import gsd.hoomd

# -----------------------------------------------------------------------------------
# DEFINE SOME LOGGING CLASSES
# -----------------------------------------------------------------------------------

# class Status():

#     def __init__(self, sim):
#         self.sim = sim

#     @property
#     def seconds_remaining(self):
#         try:
#             return (self.sim.final_timestep - self.sim.timestep) / self.sim.tps
#         except ZeroDivisionError:
#             return 0

#     @property
#     def etr(self):
#         return str(datetime.timedelta(seconds=self.seconds_remaining))



# -----------------------------------------------------------------------------------
# KEEP SAME AS IN TEST CREATE 
# -----------------------------------------------------------------------------------


# System sizes
LREF = 0.002                    # m

LX = LREF
LY = LREF
LZ = LREF

# Parameters
KERNEL  = 'WendlandC4'
NL      = 20                       # INT
FX      = 0.1                      # m/s^2

DX      = LREF/NL                  # m
V       = DX*DX*DX                 # m^3

RHO0 = 1000.0                      # kg / m^3
M    = RHO0*V                      # kg
MU   = 0.01                        # Pa s

UREF = FX*LREF*LREF*0.25/(MU/RHO0)

densitymethod = 'CONTINUITY'
# densitymethod = 'SUMMATION'

steps = 400001

# ------------------------------------------------------------------------------------



# device = hoomd.device.CPU(notice_level = 10)# , msg_file = 'log.txt')
# device = hoomd.device.CPU(notice_level = 7, msg_file = 'log.txt')
device = hoomd.device.CPU(notice_level = 2)# , msg_file = 'log.txt')
sim = hoomd.Simulation(device = device)
# if device.communicator.rank == 0:
#     print("hoomd.version.mpi_enabled: {0}".format(hoomd.version.mpi_enabled))


# filename = "test_tube.gsd"
# filename = "test_tube_343000.gsd" # NL = 60
# filename = "test_tube_125000.gsd" # NL = 40
filename = "test_tube_27000.gsd" # NL = 20
dumpname = "run_bf{0}_{1}".format(FX ,filename)

sim.create_state_from_gsd(filename = filename)

    

if SHOW_DECOMP_INFO:
    sph_info.print_decomp_info(sim, device)
# Kernel
KERNEL  = 'WendlandC4'
H       = hoomd.sph.kernel.OptimalH[KERNEL]*DX       # m
RCUT    = hoomd.sph.kernel.Kappa[KERNEL]*H           # m
Kernel = hoomd.sph.kernel.Kernels[KERNEL]()
Kappa  = Kernel.Kappa()

# print(f'H: {H}, RCUT: {RCUT}, Kappa: {Kappa}')


# print(dir(hoomd.sph._sph))
# print(hoomd.sph.kernel.Kernels['WendlandC4']
# print(dir(hoomd.sph.kernel.WendlandC4))

# # Neighbor list
NList = hoomd.nsearch.nlist.Cell(buffer = RCUT*0.05, rebuild_check_delay = 1, kappa = Kappa)

# NList = hoomd.nsearch.nlist.Stencil(buffer = RCUT*0.05, cell_width=1.5, rebuild_check_delay = 1, kappa = Kappa)
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
filterAll    = hoomd.filter.All()

with sim.state.cpu_local_snapshot as data:
    print(f'{np.count_nonzero(data.particles.typeid == 0)} fluid particles on rank {device.communicator.rank}')
    print(f'{np.count_nonzero(data.particles.typeid == 1)} solid particles on rank {device.communicator.rank}')

# # Set up SPH solver
model = hoomd.sph.sphmodel.SinglePhaseFlow(kernel = Kernel,
                                           eos    = EOS,
                                           nlist  = NList,
                                           fluidgroup_filter = filterFLUID,
                                           solidgroup_filter = filterSOLID)
if device.communicator.rank == 0:
    print("SetModelParameter on all ranks")

model.mu = MU
model.densitymethod = densitymethod
model.gx = FX
model.damp = 1000
# model.artificialviscosity = True
model.artificialviscosity = True 
model.alpha = 0.2
model.beta = 0.0
model.densitydiffusion = False
model.shepardrenormanlization = False 
# Access the local snapshot, this is not ideal! 

with sim.state.cpu_local_snapshot as snapshot:
    model.max_sl = np.max(snapshot.particles.slength[:])

dt = model.compute_dt(LREF, UREF, DX, DRHO)

integrator = hoomd.sph.Integrator(dt=dt)

# VelocityVerlet = hoomd.sph.methods.VelocityVerlet(filter=filterFLUID, densitymethod = densitymethod)
VelocityVerlet = hoomd.sph.methods.VelocityVerletBasic(filter=filterFLUID, densitymethod = densitymethod)

integrator.methods.append(VelocityVerlet)
integrator.forces.append(model)

compute_filter_all = hoomd.filter.All()
compute_filter_fluid = hoomd.filter.Type(['F'])
spf_properties = hoomd.sph.compute.SinglePhaseFlowBasicProperties(compute_filter_fluid)
sim.operations.computes.append(spf_properties)

if device.communicator.rank == 0:
    print(f'Computed Time step: {dt}')
    print("Integrator Forces: {0}".format(integrator.forces[:]))
    print("Integrator Methods: {0}".format(integrator.methods[:]))
    print("Simulation Computes: {0}".format(sim.operations.computes[:]))




gsd_trigger = hoomd.trigger.Periodic(100)
gsd_writer = hoomd.write.GSD(filename=dumpname,
                             trigger=gsd_trigger,
                             mode='wb',
                             dynamic = ['property', 'momentum']
                             )
sim.operations.writers.append(gsd_writer)

# hoomd.write.GSD.write(filename = dumpname, state = sim.state, mode = 'wb')
log_trigger = hoomd.trigger.Periodic(100)
logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps', 'walltime', 'timestep_size'])
# logger[('custom', 'ETA')] = (lambda: str(datet.timedelta(seconds=((sim.final_timestep - sim.timestep) / (1e-9 + sim.tps) ))), 'scalar')
logger[('custom', 'RE')] = (lambda: RHO0 * spf_properties.abs_velocity * LREF / (MU * spf_properties.num_particles), 'scalar')

logger.add(spf_properties, quantities=['abs_velocity', 'num_particles', 'fluid_vel_x_sum', 'mean_density'])
table = hoomd.write.Table(trigger=log_trigger, 
                          logger=logger, max_header_len = 10)
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

sim.run(steps, write_at_start=True)

if device.communicator.rank == 0:
    export_gsd2vtu.export_spf(dumpname)

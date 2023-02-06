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
# import itertools
from datetime import datetime
import export_gsd2vtu 
import read_input_fromtxt

# ------------------------------------------------------------




device = hoomd.device.CPU(notice_level=2)
# device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)




# System sizes
LREF = 0.01                    # m
radius = 0.35 * LREF

LX = LREF
LY = LREF
LZ = LREF

# Parameters
KERNEL  = 'WendlandC4'
NL      = 100                       # INT
# FX      = 0.1                      # m/s^2

DX      = LREF/NL                  # m
V       = DX*DX*DX                 # m^3

RHO0 = 1000.0                      # kg / m^3
M    = RHO0*V                      # kg
DRHO = 0.01                        # %
MU   = 0.01                        # Pa s


filename = 'Sphere_in_box_866085_init.gsd'

FX    = np.float64(sys.argv[1])

dt_string = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
logname  = filename.replace('_init.gsd', '')
logname  = f'{logname}_run_{FX}_{dt_string}.log'
dumpname = filename.replace('_init.gsd', '')
dumpname = f'{dumpname}_run_{FX}.gsd'

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


voxelsize = DX
# define model parameters
densitymethod = 'CONTINUITY'
steps = 100001
FX    = np.float64(sys.argv[1])
if device.communicator.rank == 0:
    print(f'Body Force: {FX}')

dumpname = f'bf_{FX}_{dumpname}'
if device.communicator.rank == 0:
    print(f'Dumpname: {dumpname}')


# get kernel properties
H       = hoomd.sph.kernel.OptimalH[KERNEL]*DX       # m
RCUT    = hoomd.sph.kernel.Kappa[KERNEL]*H           # m
Kernel = hoomd.sph.kernel.Kernels[KERNEL]()
Kappa  = Kernel.Kappa()

# Neighbor list
NList = hoomd.nsearch.nlist.Cell(buffer = RCUT*0.05, rebuild_check_delay = 1, kappa = Kappa)

# Equation of State
EOS = hoomd.sph.eos.Linear()
EOS.set_params(RHO0,0.05)

# Define groups/filters
filterFLUID  = hoomd.filter.Type(['F']) # is zero
filterSOLID  = hoomd.filter.Type(['S']) # is one
filterAll    = hoomd.filter.All()

with sim.state.cpu_local_snapshot as snap:
    print(f'{np.count_nonzero(snap.particles.typeid == 0)} fluid particles on rank {device.communicator.rank}')
    print(f'{np.count_nonzero(snap.particles.typeid == 1)} solid particles on rank {device.communicator.rank}')

# Set up SPH solver
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
model.damp = 5000
# model.artificialviscosity = True
model.artificialviscosity = True 
model.alpha = 0.2
model.beta = 0.0
model.densitydiffusion = False
# model.ddiff = 0.1
model.shepardrenormanlization = False 
model.densitydiffusion = False




# Access the local snapshot, this is not ideal! 
# with sim.state.cpu_local_snapshot as snap:
#     model.max_sl = np.max(snap.particles.slength[:])
    # fluid_particle_ratio = np.count_nonzero(snap.particles.typeid[:] == 0)/(len(snap.particles.mass[:]))


maximum_smoothing_length = 0.0
# Call get_snapshot on all ranks.
snapshot = sim.state.get_snapshot()
# Access particle data on rank 0 only.
if snapshot.communicator.rank == 0:
    maximum_smoothing_length = np.max(snapshot.particles.slength)
    # total_number_fluid_particles = snapshot.particles.N
maximum_smoothing_length = device.communicator.bcast_double(maximum_smoothing_length)
model.max_sl = maximum_smoothing_length

UREF = FX*LREF*LREF*0.25/(MU/RHO0)
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




gsd_trigger = hoomd.trigger.Periodic(500)
gsd_writer = hoomd.write.GSD(filename=dumpname,
                             trigger=gsd_trigger,
                             mode='wb',
                             dynamic = ['property', 'momentum']
                             )
sim.operations.writers.append(gsd_writer)



# hoomd.write.GSD.write(filename = dumpname, state = sim.state, mode = 'wb')
log_trigger = hoomd.trigger.Periodic(100)
logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps', 'walltime'])
logger.add(spf_properties, quantities=['abs_velocity', 'num_particles', 'fluid_vel_x_sum', 'mean_density'])
logger[('custom', 'RE')] = (lambda: RHO0 * spf_properties.abs_velocity * LREF / (MU), 'scalar')
# logger[('custom', 'k_1[1e-9]')] = (lambda: (MU / (RHO0 * FX)) * (spf_properties.abs_velocity) * porosity *1.0e9, 'scalar')
table = hoomd.write.Table(trigger=log_trigger, 
                          logger=logger, max_header_len = 10)

# file = open('log.txt', mode='x', newline='\n')
# table_file = hoomd.write.Table(output=file,
#                                trigger=hoomd.trigger.Periodic(period=5000),
#                                logger=logger)
# sim.operations.writers.append(table_file)
sim.operations.writers.append(table)

file = open(logname, mode='w+', newline='\n')
table_file = hoomd.write.Table(output=file,
                               trigger=log_trigger,
                               logger=logger, max_header_len = 10)
sim.operations.writers.append(table_file)

sim.operations.integrator = integrator

# print(model.loggables)
# print(sim.loggables)
# print(spf_properties.loggables)

if device.communicator.rank == 0:
    print("Starting Run at {0}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

sim.run(steps, write_at_start=True)

# if device.communicator.rank == 0:
#     export_gsd2vtu.export_spf(dumpname)
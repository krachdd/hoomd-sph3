#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: drostan, daniel.rostan@mib.uni-stuttgart.de
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
import sys, os

import gsd.hoomd
# ------------------------------------------------------------


device = hoomd.device.CPU(notice_level=2)
# device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)

filename = str(sys.argv[2])

if device.communicator.rank == 0:
    print(f'{os.path.basename(__file__)}: input file: {filename} ')

dt_string = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
logname  = filename.replace('_init.gsd', '')
logname  = f'{logname}_run.log'
dumpname = filename.replace('_init.gsd', '')
dumpname = f'{dumpname}_run.gsd'

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


# Fluid and particle properties
num_length          = int(sys.argv[1])

# Inner diameterand radius of cylinder
innerD = 0.1                              # [m]
innerR = 0.5*innerD                       # [m]

# Mix features
fw  = innerR*0.1                          # [m]
fh  = innerR*0.9                          # [m]

# Angular velocity in rad/s and 1/s
angvel_s = 1.0                           # 1/s
angvel   = angvel_s*(2*np.pi)            # rad/s

# Characteristic length and velocity
lref = innerR - fh                       # [m]
refvel = angvel*innerD

# Discretization parameters
voxelsize           = innerD/float(num_length)
dx                  = voxelsize
specific_volume     = dx * dx * dx
rho0                = 971.0              # [kg / m^3]
#rho0                = 1000.0              # [kg / m^3]
mass                = rho0 * specific_volume
viscosity           = 0.00971            # [Pa s]
#viscosity           = 0.1            # [Pa s]
#gz                  = -0.981             # [m/s**2]
gz = 0.0
backp              = 0.15
#backp              = 0.01

# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx       # m
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength     # m

# define model parameters
densitymethod = 'SUMMATION'
steps = int(sys.argv[3])

drho = 0.01                        # %

kernel_obj = hoomd.sph.kernel.Kernels[kernel]()
kappa      = kernel_obj.Kappa()

# Neighbor list
nlist = hoomd.nsearch.nlist.Cell(buffer = rcut*0.05, rebuild_check_delay = 1, kappa = kappa)

# Equation of State
eos = hoomd.sph.eos.Tait()
eos.set_params(rho0,backp)

# Define groups/filters
filterfluid  = hoomd.filter.Type(['F']) # is zero
filtersolid  = hoomd.filter.Type(['S','R']) # is one
filterrigid  = hoomd.filter.Type(['R']) # is two
filterall    = hoomd.filter.All()

with sim.state.cpu_local_snapshot as snap:
    print(f'{np.count_nonzero(snap.particles.typeid == 0)} fluid particles on rank {device.communicator.rank}')
    print(f'{np.count_nonzero(snap.particles.typeid == 1)} solid particles on rank {device.communicator.rank}')
    print(f'{np.count_nonzero(snap.particles.typeid == 2)} solid particles on rank {device.communicator.rank}')


# Set up SPH solver
model = hoomd.sph.sphmodel.SinglePhaseFlow(kernel = kernel_obj,
                                           eos    = eos,
                                           nlist  = nlist,
                                           fluidgroup_filter = filterfluid,
                                           solidgroup_filter = filtersolid,
                                           )
if device.communicator.rank == 0:
    print("SetModelParameter on all ranks")

model.mu = viscosity
model.densitymethod = densitymethod
model.gz = gz
model.damp = 1000
# model.artificialviscosity = True
model.artificialviscosity = True 
model.alpha = 0.2
model.beta = 0.0
model.densitydiffusion = False
model.shepardrenormanlization = False

maximum_smoothing_length = 0.0
# Call get_snapshot on all ranks.
snapshot = sim.state.get_snapshot()
# Access particle data on rank 0 only.
if snapshot.communicator.rank == 0:
    maximum_smoothing_length = np.max(snapshot.particles.slength)

maximum_smoothing_length = device.communicator.bcast_double(maximum_smoothing_length)
model.max_sl = maximum_smoothing_length

# compute dt
dt = model.compute_dt(lref, refvel, dx, drho)


integrator = hoomd.sph.Integrator(dt=dt)

VelocityVerlet = hoomd.sph.methods.VelocityVerlet(filter=filterFLUID, densitymethod = densitymethod)
#velocityverlet = hoomd.sph.methods.VelocityVerletBasic(filter=filterfluid, densitymethod = densitymethod)
# rigidbodyintegrator = hoomd.sph.methods.RigidBodyIntegrator(filter=filterrigid, transvel = [0.0,0.0,0.0],
#                                          rotvel   = variant.Ramp(0.0, angvel_s, 0, 1000),
#                                          pivotpnt = [0.0,0.0,0.0],
#                                          rotaxis  = [0.0,0.0,1.0])

rigidbodyintegrator = hoomd.sph.methods.RigidBodyIntegrator(filter=filterrigid, transvel_x = 0.0, transvel_y = 0.0, transvel_z = 0.0,
                                          rotvel   = variant.Ramp(0.0, angvel_s, 0, 1000),
                                          pivotpnt = [0.0,0.0,0.0],
                                          rotaxis  = [0.0,0.0,1.0])

integrator.methods.append(velocityverlet)
integrator.methods.append(rigidbodyintegrator)
integrator.forces.append(model)

compute_filter_all = hoomd.filter.All()
compute_filter_fluid = hoomd.filter.Type(['F'])
spf_properties = hoomd.sph.compute.SinglePhaseFlowBasicProperties(compute_filter_fluid)
sim.operations.computes.append(spf_properties)

if device.communicator.rank == 0:
    print(f'Computed Time step: {dt}')
    print(f'Integrator Forces: {integrator.forces[:]}')
    print(f'Integrator Methods: {integrator.methods[:]}')
    print(f'Simulation Computes: {sim.operations.computes[:]}')

gsd_trigger = hoomd.trigger.Periodic(100)
gsd_writer = hoomd.write.GSD(filename=dumpname,
                             trigger=gsd_trigger,
                             mode='wb',
                             dynamic = ['property', 'momentum']
                             )
sim.operations.writers.append(gsd_writer)

log_trigger = hoomd.trigger.Periodic(100)
logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps', 'walltime'])
logger.add(spf_properties, quantities=['abs_velocity', 'num_particles', 'fluid_vel_x_sum', 'mean_density'])

table = hoomd.write.Table(trigger=log_trigger, 
                          logger=logger, max_header_len = 10)
sim.operations.writers.append(table)

file = open(logname, mode='w+', newline='\n')
table_file = hoomd.write.Table(output=file,
                               trigger=log_trigger,
                               logger=logger, max_header_len = 10)
sim.operations.writers.append(table_file)

sim.operations.integrator = integrator

if device.communicator.rank == 0:
    print(f'Starting Run at {dt_string}')

sim.run(steps, write_at_start=True)

if device.communicator.rank == 0:
    export_gsd2vtu.export_spf(dumpname)
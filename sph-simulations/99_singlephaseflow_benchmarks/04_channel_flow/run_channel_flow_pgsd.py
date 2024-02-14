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
import math
# import itertools
from datetime import datetime
import export_gsd2vtu, delete_solids_initial_timestep 
import sph_info, sph_helper, read_input_fromtxt
import sys, os
from mpi4py import MPI

import pgsd.hoomd
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

sim.create_state_from_pgsd(filename = filename)

if SHOW_DECOMP_INFO:
    sph_info.print_decomp_info(sim, device)

# Fluid and particle properties
SHOW_PROC_PART_INFO = False
SHOW_DECOMP_INFO    = False
num_length          = int(sys.argv[1])                          # [ - ]
lref                = 0.001                                     # [ m ]
voxelsize           = lref/float(num_length)                    # [ m ]
dx                  = voxelsize                                 # [ m ]
specific_volume     = dx * dx * dx                              # [ m^3 ]
rho0                = 1000.0                                    # [ kg/m^3 ]
mass                = rho0 * specific_volume                    # [ kg ]
fx                  = 0.1                                       # [ m/s^2 ]
viscosity           = 0.01                                      # [ Pa s ]
drho                = 0.01                                      # [ % ]
backpress           = 0.01                                      # [ - ]
lidvel              = 0.01                                      # [ m/s ]
refvel              = lidvel * lref**2 * 0.25 / (viscosity/rho0)# [ m/s ]

# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx                  # [ m ]
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength                # [ m ]

# define model parameters
densitymethod = 'SUMMATION'
steps = int(sys.argv[3])

kernel_obj = hoomd.sph.kernel.Kernels[kernel]()
kappa      = kernel_obj.Kappa()

# Neighbor list
nlist = hoomd.nsearch.nlist.Cell(buffer = rcut*0.05, rebuild_check_delay = 1, kappa = kappa)

# Equation of State
eos = hoomd.sph.eos.Tait()
eos.set_params( rho0, backpress )

# Define groups/filters
filterfluid  = hoomd.filter.Type(['F']) # is zero
filtersolid  = hoomd.filter.Type(['S']) # is one
filterall    = hoomd.filter.All()

with sim.state.cpu_local_snapshot as snap:
    print(f'{np.count_nonzero(snap.particles.typeid == 0)} fluid particles on rank {device.communicator.rank}')
    print(f'{np.count_nonzero(snap.particles.typeid == 1)} solid particles on rank {device.communicator.rank}')

# Set up SPH solver
model = hoomd.sph.sphmodel.SinglePhaseFlow(kernel = kernel_obj,
                                           eos    = eos,
                                           nlist  = nlist,
                                           fluidgroup_filter = filterfluid,
                                           solidgroup_filter = filtersolid)
if device.communicator.rank == 0:
    print("SetModelParameter on all ranks")

model.mu = viscosity
model.densitymethod = densitymethod
model.gx = fx
model.damp = 1000
model.artificialviscosity = True 
model.alpha = 0.2
model.beta = 0.0
model.densitydiffusion = False
model.shepardrenormanlization = False

maximum_smoothing_length = sph_helper.set_max_sl(sim, device, snapshot, model)

c, c_condition = model.compute_speedofsound(LREF = lref, UREF = refvel, 
                                            DX = dx, DRHO = drho, H = maximum_smoothing_length, 
                                            MU = viscosity, RHO0 = rho0)

if device.communicator.rank == 0:
    print(f'Speed of sound [m/s]: {c}, Used: {c_condition}')

# compute dt
dt, dt_condition = model.compute_dt(LREF = lref, UREF = refvel, 
                                          DX = dx, DRHO = drho, H = maximum_smoothing_length, 
                                          MU = viscosity, RHO0 = rho0)

if device.communicator.rank == 0:
    print(f'Timestep size [s]: {dt}, Used: {dt_condition}')

integrator = hoomd.sph.Integrator(dt=dt)

velocityverlet = hoomd.sph.methods.VelocityVerletBasic(filter=filterfluid, densitymethod = densitymethod)

integrator.methods.append(velocityverlet)
integrator.forces.append(model)

compute_filter_all = hoomd.filter.All()
compute_filter_fluid = hoomd.filter.Type(['F'])
spf_properties = hoomd.sph.compute.SinglePhaseFlowBasicProperties(compute_filter_fluid)
sim.operations.computes.append(spf_properties)

if device.communicator.rank == 0:
    print(f'Integrator Forces: {integrator.forces[:]}')
    print(f'Integrator Methods: {integrator.methods[:]}')
    print(f'Simulation Computes: {sim.operations.computes[:]}')

pgsd_trigger = hoomd.trigger.Periodic(100)
pgsd_writer = hoomd.write.PGSD(filename=dumpname,
                             trigger=pgsd_trigger,
                             mode='wb',
                             dynamic = ['property', 'momentum']
                             )
sim.operations.writers.append(pgsd_writer)

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
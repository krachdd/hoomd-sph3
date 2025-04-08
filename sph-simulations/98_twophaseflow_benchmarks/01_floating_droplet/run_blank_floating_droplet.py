#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""
# ----- HEADER -----------------------------------------------
import sys, os
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
from optparse import OptionParser

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

# Fluid and particle properties
SHOW_PROC_PART_INFO = False
SHOW_DECOMP_INFO    = False
num_length          = int(sys.argv[1])                #  int
lref                = 0.001                           # [ m ]
radius              = 0.5 * lref                      # [ m ]
voxelsize           = lref/num_length                 # [ m ]
dx                  = voxelsize                       # [ m ]
specific_volume     = dx * dx * dx                    # [ m**3 ]
rho01               = 1000.0                          # [ kg/m**3 ]
rho02               = 1000.0                          # [ kg/m**3 ]
mass1               = rho01 * specific_volume         # [ kg ]
mass2               = rho02 * specific_volume         # [ kg ]
fx                  = 0.0                             # [ m/s ]
viscosity1          = 0.001                           # [ Pa s ]
viscosity2          = 0.001                           # [ Pa s ]
backpress           = 0.01                             # [ - ]
sigma               = 0.01                            # [ Pa/m**2 ]
contact_angle       = 60                              # [ ° ]
refvel              = fx
drho                = 0.01                            # [ - ] %
steps               = 1001                           # [ - ]

# get kernel properties
kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel]*dx        # [ m ]
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength      # [ m ]
kernel_obj = hoomd.sph.kernel.Kernels[kernel]()
kappa      = kernel_obj.Kappa()

# define model parameters
densitymethod = 'SUMMATION'
colorgradientmethod = 'DENSITYRATIO'

if SHOW_DECOMP_INFO:
    sph_info.print_decomp_info(sim, device)

# Neighbor list
nlist = hoomd.nsearch.nlist.Cell(buffer = rcut*0.05, rebuild_check_delay = 1, kappa = kappa)

# Equation of State
eos1 = hoomd.sph.eos.Tait()
eos2 = hoomd.sph.eos.Tait()
eos1.set_params( rho01, backpress )
eos2.set_params( rho02, backpress )

# Define groups/filters
filterfluidW = hoomd.filter.Type(['W']) # is zero
filterfluidN = hoomd.filter.Type(['N']) # is one
# filtersolid  = hoomd.filter.Type(['S']) # is two
filterall    = hoomd.filter.All()

if SHOW_PROC_PART_INFO:
    with sim.state.cpu_local_snapshot as snap:
        print(f'{np.count_nonzero(snap.particles.typeid == 0)} fluid particles on rank {device.communicator.rank}')
        print(f'{np.count_nonzero(snap.particles.typeid == 1)} solid particles on rank {device.communicator.rank}')

# Set up SPH solver
model = hoomd.sph.sphmodel.TwoPhaseFlow(kernel = kernel_obj,
                                        eos1   = eos1,
                                        eos2   = eos2,
                                        nlist  = nlist,
                                        fluidgroup1_filter = filterfluidW,
                                        fluidgroup2_filter = filterfluidN,
                                        #solidgroup_filter = filtersolid)
                                        densitymethod = densitymethod,
                                        colorgradientmethod = colorgradientmethod)

model.mu1 = viscosity1
model.mu2 = viscosity2
model.sigma12 = sigma
model.omega = contact_angle
model.densitymethod = densitymethod
model.gx = fx
model.damp = 1000
model.artificialviscosity = True 
model.alpha = 0.4
model.beta = 0.2
model.densitydiffusion = False
model.shepardrenormanlization = False
model.fickian_shifting = False

maximum_smoothing_length = sph_helper.set_max_sl(sim, device, snapshot, model)

c1, c1_condition, c2, c2_condition = model.compute_speedofsound(LREF = lref, UREF = refvel, 
                                            DX = dx, DRHO = drho, H = maximum_smoothing_length, 
                                            MU1 = viscosity1, MU2 = viscosity2, 
                                            RHO01 = rho01, RHO02 = rho02,
                                            SIGMA12 = sigma)

if device.communicator.rank == 0:
    print(f'Fluid1: Speed of sound [m/s]: {c1}, Used: {c1_condition}')
    print(f'Fluid2: Speed of sound [m/s]: {c2}, Used: {c2_condition}')

# sph_helper.update_min_c0_tpf(device, model, c1, c2, mode = 'plain', lref = lref, uref = refvel, bforce = fx, cfactor = 10.0)

# compute dt
dt, dt_condition = model.compute_dt(LREF = lref, UREF = refvel, 
                                          DX = dx, DRHO = drho, H = maximum_smoothing_length, 
                                          MU1 = viscosity1, MU2 = viscosity2, 
                                          RHO01 = rho01, RHO02 = rho02,
                                          SIGMA12 = sigma)

if device.communicator.rank == 0:
    print(f'Timestep size [s]: {dt}, Used: {dt_condition}')

integrator = hoomd.sph.Integrator(dt=dt)

# VelocityVerlet = hoomd.sph.methods.VelocityVerlet(filter=filterFLUID, densitymethod = densitymethod)
velocityverletW = hoomd.sph.methods.VelocityVerletBasic(filter=filterfluidW, densitymethod = densitymethod)
velocityverletN = hoomd.sph.methods.VelocityVerletBasic(filter=filterfluidN, densitymethod = densitymethod)

integrator.methods.append(velocityverletW)
integrator.methods.append(velocityverletN)
integrator.forces.append(model)

# compute_filter_all = hoomd.filter.All()
# compute_filter_fluid = hoomd.filter.Type(['F'])
# spf_properties = hoomd.sph.compute.SinglePhaseFlowBasicProperties(compute_filter_fluid)
# sim.operations.computes.append(spf_properties)

if device.communicator.rank == 0:
    print(f'Computed Time step: {dt}')
    print(f'Integrator Forces: {integrator.forces[:]}')
    print(f'Integrator Methods: {integrator.methods[:]}')
    print(f'Simulation Computes: {sim.operations.computes[:]}')

gsd_trigger = hoomd.trigger.Periodic(1)
gsd_writer = hoomd.write.GSD(filename=dumpname,
                             trigger=gsd_trigger,
                             mode='wb',
                             dynamic = ['property', 'momentum']
                             )
sim.operations.writers.append(gsd_writer)

log_trigger = hoomd.trigger.Periodic(1)
logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps', 'walltime'])
# logger.add(spf_properties, quantities=['abs_velocity', 'num_particles', 'fluid_vel_x_sum', 'mean_density'])

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
    export_gsd2vtu.export_tpf(dumpname)
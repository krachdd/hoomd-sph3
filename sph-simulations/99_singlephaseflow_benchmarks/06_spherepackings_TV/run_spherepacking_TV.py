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
import export_gsd2vtu, delete_solids_initial_timestep 
import sph_info, sph_helper, read_input_fromtxt

# ------------------------------------------------------------

device = hoomd.device.CPU(notice_level=2)
# device = hoomd.device.CPU(notice_level=10)
sim = hoomd.Simulation(device=device)

# get stuff from input file
infile = str(sys.argv[1])
params = read_input_fromtxt.get_input_data_from_file(infile)

FX    = 0.01

rawfilename = params['rawfilename']
filename = rawfilename.replace('.raw', '_init.gsd')
dt_string = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
logname  = filename.replace('_init.gsd', '')
logname  = f'{logname}_run_{FX}_{dt_string}_TV.log'
dumpname = filename.replace('_init.gsd', '')
dumpname = f'{dumpname}_run_{FX}_TV.gsd'

sim.create_state_from_gsd(filename = filename)



# Define necessary parameters
# Fluid and particle properties
SHOW_PROC_PART_INFO = False
SHOW_DECOMP_INFO    = False
if SHOW_DECOMP_INFO:
    sph_info.print_decomp_info(sim, device)
lref                = 0.001                                     # [ m ]
voxelsize           = np.float64(params['vsize'])               # [ m ]
dx                  = voxelsize                                 # [ m ]
specific_volume     = dx * dx * dx                              # [ m^3 ]
rho0                = np.float64(params['fdensity'])            # [ kg/m^3 ]
mass                = rho0 * specific_volume                    # [ kg ]
fx                  = 0.1                                       # [ m/s^2 ]
viscosity           = np.float64(params['fviscosity'])          # [ Pa s ]
drho                = 0.01                                      # [ % ]
backpress           = 0.01                                      # [ - ]
refvel              = fx * lref**2 * 0.25 / (viscosity/rho0)    # [ m/s ]
porosity            = np.float64(params['porosity'])            # [ - ]

# get kernel properties
kernel  = params['kernel']
slength = hoomd.sph.kernel.OptimalH[kernel]*dx                  # [ m ]
rcut    = hoomd.sph.kernel.Kappa[kernel]*slength                # [ m ]

kernel_obj = hoomd.sph.kernel.Kernels[kernel]()
kappa      = kernel_obj.Kappa()

# get simulation box sizes etc.
NX, NY, NZ = np.int32(params['nx']), np.int32(params['ny']), np.int32(params['nz']) 
lref = NX * voxelsize
Re = 0.01
refvel = (Re * viscosity)/(rho0 * lref)

# define model parameters
densitymethod = 'SUMMATION'
steps = 30001

# Neighbor list
NList = hoomd.nsearch.nlist.Cell(buffer = rcut*0.05, rebuild_check_delay = 1, kappa = kappa)

# Equation of State
eos = hoomd.sph.eos.Tait()
eos.set_params( rho0, backpress )

# Define groups/filters
filterFLUID  = hoomd.filter.Type(['F']) # is zero
filterSOLID  = hoomd.filter.Type(['S']) # is one
filterAll    = hoomd.filter.All()

with sim.state.cpu_local_snapshot as snap:
    print(f'{np.count_nonzero(snap.particles.typeid == 0)} fluid particles on rank {device.communicator.rank}')
    print(f'{np.count_nonzero(snap.particles.typeid == 1)} solid particles on rank {device.communicator.rank}')

# Set up SPH solver
model = hoomd.sph.sphmodel.SinglePhaseFlowTV(kernel = kernel_obj,
                                           eos    = eos,
                                           nlist  = NList,
                                           fluidgroup_filter = filterFLUID,
                                           solidgroup_filter = filterSOLID,
                                           densitymethod = densitymethod)
if device.communicator.rank == 0:
    print("SetModelParameter on all ranks")

model.mu = viscosity
model.densitymethod = densitymethod
model.gx = fx
model.damp = 5000
model.artificialviscosity = True 
model.alpha = 0.2
model.beta = 0.2
model.densitydiffusion = False
model.shepardrenormanlization = False 

maximum_smoothing_length = sph_helper.set_max_sl(sim, device, model)

c, c_condition = model.compute_speedofsound(LREF = lref, UREF = refvel, 
                                            DX = dx, DRHO = drho, H = maximum_smoothing_length, 
                                            MU = viscosity, RHO0 = rho0)

if device.communicator.rank == 0:
    print(f'Speed of sound [m/s]: {c}, Used: {c_condition}')

# cfactor = 10
# if c < cfactor * refvel:
#     model.set_speedofsound(cfactor * refvel)
#     if device.communicator.rank == 0:
#         print(f'Increase Speed of Sound to adami condition: {cfactor} * revel: {model.get_speedofsound()}')

# compute dt
dt, dt_condition = model.compute_dt(LREF = lref, UREF = refvel, 
                                          DX = dx, DRHO = drho, H = maximum_smoothing_length, 
                                          MU = viscosity, RHO0 = rho0)

if device.communicator.rank == 0:
    print(f'Timestep size [s]: {dt}, Used: {dt_condition}')

integrator = hoomd.sph.Integrator(dt=dt)

# VelocityVerlet = hoomd.sph.methods.VelocityVerlet(filter=filterFLUID, densitymethod = densitymethod)
kdktv = hoomd.sph.methods.KickDriftKickTV(filter=filterFLUID, densitymethod = densitymethod)

integrator.methods.append(kdktv)
integrator.forces.append(model)

compute_filter_all = hoomd.filter.All()
compute_filter_fluid = hoomd.filter.Type(['F'])
spf_properties = hoomd.sph.compute.SinglePhaseFlowBasicProperties(compute_filter_fluid)
sim.operations.computes.append(spf_properties)

if device.communicator.rank == 0:
    print(f'Integrator Forces: {integrator.forces[:]}')
    print(f'Integrator Methods: {integrator.methods[:]}')
    print(f'Simulation Computes: {sim.operations.computes[:]}')




gsd_trigger = hoomd.trigger.Periodic(3000)
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
logger.add(spf_properties, quantities=['abs_velocity', 'num_particles', 'fluid_vel_x_sum', 'mean_density', 'e_kin_fluid'])
logger[('custom', 'RE')] = (lambda: rho0 * spf_properties.abs_velocity * lref / (viscosity), 'scalar')
logger[('custom', 'k_1[1e-9]')] = (lambda: (viscosity / (rho0 * fx)) * (spf_properties.abs_velocity) * porosity *1.0e9, 'scalar')
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

# if device.communicator.rank == 0:
#     export_gsd2vtu.export_spf(dumpname)
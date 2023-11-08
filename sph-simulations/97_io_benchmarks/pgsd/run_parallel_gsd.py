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

# get stuff from input file
infile = str(sys.argv[1])
params = read_input_fromtxt.get_input_data_from_file(infile)
if device.communicator.rank == 0:
    print(params)

FX    = 0.001

rawfilename = params['rawfilename']
filename = rawfilename.replace('.raw', '_init.gsd')
dt_string = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
logname  = filename.replace('_init.gsd', '')
logname  = f'{logname}_run_gsd_{dt_string}.log'
dumpname = filename.replace('_init.gsd', '')
dumpname = f'{dumpname}_run_gsd.gsd'

# sim.create_state_from_pgsd(filename = filename)
# MPI.COMM_WORLD.Barrier()

# # Print the domain decomposition.
# domain_decomposition = sim.state.domain_decomposition
# if device.communicator.rank == 0:
#     print(f'Domain Decomposition: {domain_decomposition}')

# # Print the location of the split planes.
# split_fractions = sim.state.domain_decomposition_split_fractions
# if device.communicator.rank == 0:
#     print(f'Locations of SplitPlanes: {split_fractions}')

# # Print the number of particles on each rank.
# with sim.state.cpu_local_snapshot as snap:
#     N = len(snap.particles.position)
#     print(f'{N} particles on rank {device.communicator.rank}')

# # Define necessary parameters
# # Fluid and particle properties
# voxelsize  = np.float64(params['vsize'])
# DX   = voxelsize
# V    = DX * DX * DX
# RHO0 = np.float64(params['fdensity'])
# MU   = np.float64(params['fviscosity'])
# M    = RHO0 * V
# porosity = np.float64(params['porosity'])
# # get simulation box sizes etc.
# NX, NY, NZ = np.int32(params['nx']), np.int32(params['ny']), np.int32(params['nz']) 
# LREF = NX * voxelsize
# # define model parameters
# densitymethod = 'CONTINUITY'
# steps = 31

# DRHO = 0.01                        # %

# # get kernel properties
# KERNEL  = params['kernel']
# H       = hoomd.sph.kernel.OptimalH[KERNEL]*DX       # m
# RCUT    = hoomd.sph.kernel.Kappa[KERNEL]*H           # m
# Kernel = hoomd.sph.kernel.Kernels[KERNEL]()
# Kappa  = Kernel.Kappa()

# # Neighbor list
# NList = hoomd.nsearch.nlist.Cell(buffer = RCUT*0.05, rebuild_check_delay = 1, kappa = Kappa)

# # Equation of State
# EOS = hoomd.sph.eos.Tait()
# EOS.set_params(RHO0,0.1)

# # Define groups/filters
# filterFLUID  = hoomd.filter.Type(['F']) # is zero
# filterSOLID  = hoomd.filter.Type(['S']) # is one
# filterAll    = hoomd.filter.All()

# with sim.state.cpu_local_snapshot as snap:
#     print(f'{np.count_nonzero(snap.particles.typeid == 0)} fluid particles on rank {device.communicator.rank}')
#     print(f'{np.count_nonzero(snap.particles.typeid == 1)} solid particles on rank {device.communicator.rank}')

# # Set up SPH solver
# model = hoomd.sph.sphmodel.SinglePhaseFlow(kernel = Kernel,
#                                            eos    = EOS,
#                                            nlist  = NList,
#                                            fluidgroup_filter = filterFLUID,
#                                            solidgroup_filter = filterSOLID)
# if device.communicator.rank == 0:
#     print("SetModelParameter on all ranks")

# model.mu = MU
# model.densitymethod = densitymethod
# model.gx = FX
# model.damp = 5000
# # model.artificialviscosity = True
# model.artificialviscosity = True 
# model.alpha = 0.2
# model.beta = 0.0
# model.densitydiffusion = False
# model.shepardrenormanlization = False 




# # Access the local snapshot, this is not ideal! 
# # with sim.state.cpu_local_snapshot as snap:
# #     model.max_sl = np.max(snap.particles.slength[:])
#     # fluid_particle_ratio = np.count_nonzero(snap.particles.typeid[:] == 0)/(len(snap.particles.mass[:]))


# maximum_smoothing_length = 0.0
# # Call get_snapshot on all ranks.
# snapshot = sim.state.get_snapshot()
# # Access particle data on rank 0 only.
# if snapshot.communicator.rank == 0:
#     maximum_smoothing_length = np.max(snapshot.particles.slength)
#     # total_number_fluid_particles = snapshot.particles.N
# maximum_smoothing_length = device.communicator.bcast_double(maximum_smoothing_length)
# model.max_sl = maximum_smoothing_length


# reference_length = NX * DX
# # reference_length = NX * DX
# # Compute reference_velocity via Reynolds number definition
# Re = 0.01
# reference_velocity = (Re * MU)/(RHO0 * reference_length)
# # Compute reference_velocity via permeability
# # reference_velocity = porosity * 
# # mydict['pestimate']*((mydict['lref']**2)/(8*options.mu))*options.rho0*options.fz

# dt = model.compute_dt(reference_length, FX, DX, DRHO)

# integrator = hoomd.sph.Integrator(dt=dt)

# # VelocityVerlet = hoomd.sph.methods.VelocityVerlet(filter=filterFLUID, densitymethod = densitymethod)
# VelocityVerlet = hoomd.sph.methods.VelocityVerletBasic(filter=filterFLUID, densitymethod = densitymethod)

# integrator.methods.append(VelocityVerlet)
# integrator.forces.append(model)

# compute_filter_all = hoomd.filter.All()
# compute_filter_fluid = hoomd.filter.Type(['F'])
# spf_properties = hoomd.sph.compute.SinglePhaseFlowBasicProperties(compute_filter_fluid)
# sim.operations.computes.append(spf_properties)

# if device.communicator.rank == 0:
#     print(f'Computed Time step: {dt}')
#     print(f'Integrator Forces: {integrator.forces[:]}')
#     print(f'Integrator Methods: {integrator.methods[:]}')
#     print(f'Simulation Computes: {sim.operations.computes[:]}')




# gsd_trigger = hoomd.trigger.Periodic(10)
# gsd_writer = hoomd.write.GSD(filename=dumpname,
#                              trigger=gsd_trigger,
#                              mode='wb',
#                              dynamic = ['property', 'momentum']
#                              )
# sim.operations.writers.append(gsd_writer)



# # hoomd.write.GSD.write(filename = dumpname, state = sim.state, mode = 'wb')
# log_trigger = hoomd.trigger.Periodic(1)
# logger = hoomd.logging.Logger(categories=['scalar', 'string'])
# logger.add(sim, quantities=['timestep', 'tps', 'walltime'])
# logger.add(spf_properties, quantities=['abs_velocity', 'num_particles', 'fluid_vel_x_sum', 'mean_density'])
# logger[('custom', 'RE')] = (lambda: RHO0 * spf_properties.abs_velocity * LREF / (MU), 'scalar')
# logger[('custom', 'k_1[1e-9]')] = (lambda: (MU / (RHO0 * FX)) * (spf_properties.abs_velocity) * porosity *1.0e9, 'scalar')
# table = hoomd.write.Table(trigger=log_trigger, 
#                           logger=logger, max_header_len = 10)
# sim.operations.writers.append(table)

# file = open(logname, mode='w+', newline='\n')
# table_file = hoomd.write.Table(output=file,
#                                trigger=log_trigger,
#                                logger=logger, max_header_len = 10)
# sim.operations.writers.append(table_file)

# sim.operations.integrator = integrator

# if device.communicator.rank == 0:
#     print(f'Starting Run at {dt_string}')

# sim.run(steps, write_at_start=True)

# if device.communicator.rank == 0:
#     export_gsd2vtu.export_spf(dumpname)
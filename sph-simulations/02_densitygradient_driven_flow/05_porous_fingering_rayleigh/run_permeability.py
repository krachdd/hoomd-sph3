#!/usr/bin/env python3
"""
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de

Darcy permeability pre-run for the disordered-pack fingering benchmark.

Measures K of the disordered cylinder pack on the fully periodic cell from
    create_input_geometry.py <num_diameter> [seed] --perm
via body-force-driven Stokes flow in x:  K = mu * phi * <u> / (rho * fx),
averaged over the last 25% of the run.  Kozeny-Carman with the mean grain
diameter serves as the order-of-magnitude reference.

The measured K defines the Rayleigh number of run_fingering_ra.py:
    Ra = drho * g * K * H / (phi * mu * D).

Usage:
    python3 run_permeability.py <init_gsd> [fx] [steps] [damp]
      init_gsd : porous_ra_perm_*_init.gsd
      fx       : body force in x [m/s^2]  (default 1e-4, Darcy regime)
      steps    : simulation steps         (default 20001)
      damp     : body-force ramp [steps]  (default 2000)
"""

import sys, os
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules', 'gsd2vtu'))

import hoomd
from hoomd import sph
import numpy as np
from datetime import datetime
import gsd.hoomd
import sph_helper

try:
    import export_gsd2vtu
    HAS_VTU = True
except ImportError:
    HAS_VTU = False

# ─── Device & simulation ─────────────────────────────────────────────────────
device = hoomd.device.CPU(notice_level=2)
sim    = hoomd.Simulation(device=device)

filename = str(sys.argv[1])
fx       = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4
steps    = int(sys.argv[3])   if len(sys.argv) > 3 else 20001
damp     = int(sys.argv[4])   if len(sys.argv) > 4 else 2000

dt_string = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
tag       = f'fx{fx:.4e}'
logname   = filename.replace('_init.gsd', f'_{tag}_run.log')
dumpname  = filename.replace('_init.gsd', f'_{tag}_run.gsd')

sim.create_state_from_gsd(filename=filename)

# ─── Physical parameters (must match create_input_geometry.py) ───────────────
d_mean = 0.4e-3        # mean grain diameter       [m]

rho0      = 1000.0     # rest density              [kg/m³]
viscosity = 0.001      # dynamic viscosity         [Pa·s]
drho      = 0.01       # allowed density variation [–]
backpress = 0.01       # background pressure coeff [–]

phi_val   = 0.0
vsize_val = 0.0
snapshot  = sim.state.get_snapshot()
if snapshot.communicator.rank == 0:
    n_total   = float(len(snapshot.particles.typeid))
    n_fluid   = float(np.sum(snapshot.particles.typeid == 0))
    phi_val   = n_fluid / n_total
    vsize_val = float(np.max(snapshot.particles.slength)) \
                / hoomd.sph.kernel.OptimalH['WendlandC4']

phi = device.communicator.bcast_double(phi_val)
dx  = device.communicator.bcast_double(vsize_val)

# Kozeny–Carman reference (order of magnitude for a disordered pack)
K_KC = d_mean**2 * phi**3 / (180.0 * (1.0 - phi)**2)

U_D_ref = K_KC * rho0 * fx / viscosity
refvel  = max(U_D_ref / phi, 1e-10)
lref    = d_mean                       # pore scale

if device.communicator.rank == 0:
    print(f'Disordered-pack permeability run:  fx = {fx:.4e} m/s²')
    print(f'  dx = {dx:.4e} m,  d_mean = {d_mean:.1e} m,  φ = {phi:.4f}')
    print(f'  K_KC = {K_KC:.4e} m²  (Kozeny–Carman reference)')
    print(f'  U_D prediction = {U_D_ref:.4e} m/s')

# ─── Kernel ──────────────────────────────────────────────────────────────────
kernel     = 'WendlandC4'
slength    = hoomd.sph.kernel.OptimalH[kernel] * dx
rcut       = hoomd.sph.kernel.Kappa[kernel] * slength
kernel_obj = hoomd.sph.kernel.Kernels[kernel]()
kappa      = kernel_obj.Kappa()

nlist = hoomd.nsearch.nlist.Cell(buffer=rcut * 0.05,
                                  rebuild_check_delay=1, kappa=kappa)

eos = hoomd.sph.eos.Tait()
eos.set_params(rho0, backpress)

filterfluid = hoomd.filter.Type(['F'])
filtersolid = hoomd.filter.Type(['S'])

model = hoomd.sph.sphmodel.SinglePhaseFlow(
    kernel=kernel_obj, eos=eos, nlist=nlist,
    fluidgroup_filter=filterfluid, solidgroup_filter=filtersolid,
    densitymethod='SUMMATION')

model.mu                  = viscosity
model.gx                  = fx
model.damp                = damp
model.artificialviscosity = True
model.alpha               = 0.2
model.beta                = 0.0
model.densitydiffusion    = False

# ─── Speed of sound & timestep ───────────────────────────────────────────────
maximum_smoothing_length = sph_helper.set_max_sl(sim, device, model)

c, cond = model.compute_speedofsound(
    LREF=lref, UREF=refvel, DX=dx, DRHO=drho,
    H=maximum_smoothing_length, MU=viscosity, RHO0=rho0)

if device.communicator.rank == 0:
    print(f'Speed of sound: {c:.4f} m/s  ({cond})')

dt, dt_cond = model.compute_dt(
    LREF=lref, UREF=refvel, DX=dx, DRHO=drho,
    H=maximum_smoothing_length, MU=viscosity, RHO0=rho0)

if device.communicator.rank == 0:
    print(f'Timestep: {dt:.3e} s  ({dt_cond})')

integrator = hoomd.sph.Integrator(dt=dt)
vvb = hoomd.sph.methods.VelocityVerletBasic(filter=filterfluid,
                                              densitymethod='SUMMATION')
integrator.methods.append(vvb)
integrator.forces.append(model)
sim.operations.integrator = integrator

spf_properties = hoomd.sph.compute.SinglePhaseFlowBasicProperties(filterfluid)
sim.operations.computes.append(spf_properties)

# ─── Output ──────────────────────────────────────────────────────────────────
try:
    os.remove(dumpname)
except (FileNotFoundError, OSError):
    pass
device.communicator.barrier()

gsd_period = max(1, steps // 20)
gsd_writer = hoomd.write.GSD(filename=dumpname,
                              trigger=hoomd.trigger.Periodic(gsd_period),
                              mode='wb',
                              dynamic=['property', 'momentum'])
sim.operations.writers.append(gsd_writer)

logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps', 'walltime'])
logger.add(spf_properties,
           quantities=['abs_velocity', 'mean_density', 'e_kin_fluid'])
logger[('custom', 'k_m2')] = (
    lambda: viscosity * phi * spf_properties.abs_velocity / (rho0 * fx),
    'scalar')

log_period = max(1, steps // 200)
table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(log_period),
                          logger=logger, max_header_len=10)
sim.operations.writers.append(table)

log_file = open(logname, mode='w+', newline='\n')
table_file = hoomd.write.Table(output=log_file,
                                trigger=hoomd.trigger.Periodic(log_period),
                                logger=logger, max_header_len=10)
sim.operations.writers.append(table_file)

# ─── Run ─────────────────────────────────────────────────────────────────────
if device.communicator.rank == 0:
    print(f'Starting permeability run  fx={fx:.4e}  at {dt_string}')
sim.run(steps, write_at_start=True)
gsd_writer.flush()

# ─── Post-processing: steady-state K ─────────────────────────────────────────
if device.communicator.rank == 0:
    log_file.seek(0)
    raw = log_file.read()
    log_file.close()

    lines = [l for l in raw.splitlines() if l and not l.startswith('#')]
    header_idx = next((i for i, l in enumerate(lines) if 'k_m2' in l), None)

    if header_idx is not None:
        header = lines[header_idx].split()
        data_lines = [l for l in lines[header_idx + 1:]
                      if l.strip() and (l.lstrip()[0].isdigit() or l.lstrip()[0] == '-')]
        n_tail = max(1, len(data_lines) // 4)
        tail   = data_lines[-n_tail:]
        col_k  = next(i for i, h in enumerate(header) if 'k_m2' in h)
        k_ss   = float(np.mean([float(r.split()[col_k]) for r in tail]))

        print(f'\n── Disordered-pack permeability summary  (fx = {fx:.4e} m/s²) ──')
        print(f'  Steady-state K = {k_ss:.4e} m²   (K/K_KC = {k_ss/K_KC:.3f})')
        print(f'  φ = {phi:.4f}')
        print(f'\n  Pass this K to run_fingering_ra.py:')
        print(f'    python3 run_fingering_ra.py <pack_init.gsd> {k_ss:.4e} <D_m2s> [steps]')

        summary_file = filename.replace('_init.gsd', '_permeability_summary.dat')
        write_header = not os.path.exists(summary_file)
        with open(summary_file, 'a') as sf:
            if write_header:
                sf.write('# fx[m/s2]  K[m2]  K_KC[m2]  phi\n')
            sf.write(f'{fx:.6e}  {k_ss:.6e}  {K_KC:.6e}  {phi:.6f}\n')
    else:
        print('Warning: could not parse log file for post-processing.')

if HAS_VTU and device.communicator.rank == 0:
    export_gsd2vtu.export_spf(dumpname)

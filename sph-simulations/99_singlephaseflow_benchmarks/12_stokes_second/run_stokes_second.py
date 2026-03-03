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

Stokes' second problem (oscillating plate) — WCSPH benchmark.

BENCHMARK DESCRIPTION
---------------------
A fluid layer of height $H = l_\mathrm{ref}$ above a plate oscillating at
$v_\mathrm{wall}(t) = U_0 \cos(\omega t)$.  The upper plate is stationary.
$H \gg \delta$ so that the upper plate has negligible influence.

Analytical steady-state solution (Stokes oscillating boundary layer):
    $v_x(y', t) = U_0 \exp(-y'/\delta) \cos(\omega t - y'/\delta)$
    where $y' = y + H/2$ (distance from oscillating plate)
    and $\delta = \sqrt{2\nu/\omega}$ (Stokes layer thickness)

This benchmark tests:
  - Correct transient viscous diffusion from a moving boundary
  - Adami 2012 no-slip boundary with time-varying wall velocity

The oscillating plate velocity is updated every time step via a Python loop.
The run is split into batches of one period each; the final velocity profile
(at the phase omega*t = 0) is compared to the analytical amplitude profile.

Usage:
    python3 run_stokes_second.py <num_length> <init_gsd_file> [n_periods]
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

num_length = int(sys.argv[1])
filename   = str(sys.argv[2])
n_periods  = int(sys.argv[3]) if len(sys.argv) > 3 else 20

dt_string = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
logname   = filename.replace('_init.gsd', '_runWC.log')
dumpname  = filename.replace('_init.gsd', '_runWC.gsd')

sim.create_state_from_gsd(filename=filename, domain_decomposition=(None, None, 1))

# ─── Physical parameters ─────────────────────────────────────────────────────
lref       = 0.001          # fluid layer height H                  [m]
dx         = lref / num_length
rho0       = 1000.0         # rest density                          [kg/m³]
viscosity  = 0.01           # dynamic viscosity                     [Pa·s]
U0         = 0.001          # plate velocity amplitude              [m/s]
omega      = 2.0 * np.pi * 100.0  # angular frequency (100 Hz)     [rad/s]
drho       = 0.01
backpress  = 0.01
nu         = viscosity / rho0
delta      = np.sqrt(2.0 * nu / omega)  # $\delta = \sqrt{2\nu/\omega}$ [m]
T_period   = 2.0 * np.pi / omega       # $T = 2\pi/\omega$ [s]

if device.communicator.rank == 0:
    print(f'Stokes 2nd: nu={nu:.4e}, omega={omega:.2f} rad/s')
    print(f'  Stokes layer delta = {delta:.4e} m  (vs dx = {dx:.4e} m)')
    print(f'  Ratio H/delta = {lref/delta:.2f}  (should be >> 1)')

# ─── Kernel ──────────────────────────────────────────────────────────────────
kernel     = 'WendlandC4'
slength    = hoomd.sph.kernel.OptimalH[kernel] * dx
rcut       = hoomd.sph.kernel.Kappa[kernel] * slength
kernel_obj = hoomd.sph.kernel.Kernels[kernel]()
kappa      = kernel_obj.Kappa()

nlist = hoomd.nsearch.nlist.Cell(buffer=rcut * 0.05,
                                  rebuild_check_delay=1, kappa=kappa)

# ─── EOS ─────────────────────────────────────────────────────────────────────
eos = hoomd.sph.eos.Tait()
eos.set_params(rho0, backpress)

# ─── Filters ─────────────────────────────────────────────────────────────────
filterfluid = hoomd.filter.Type(['F'])
filtersolid = hoomd.filter.Type(['S'])

# ─── SPH model ───────────────────────────────────────────────────────────────
model = hoomd.sph.sphmodel.SinglePhaseFlow(
    kernel=kernel_obj, eos=eos, nlist=nlist,
    fluidgroup_filter=filterfluid, solidgroup_filter=filtersolid,
    densitymethod='SUMMATION')

model.mu               = viscosity
model.gx               = 0.0
model.damp             = 0
model.artificialviscosity = True
model.alpha            = 0.1
model.beta             = 0.0
model.densitydiffusion = False

# ─── Speed of sound & timestep ───────────────────────────────────────────────
maximum_smoothing_length = sph_helper.set_max_sl(sim, device, model)

c, cond = model.compute_speedofsound(
    LREF=lref, UREF=U0, DX=dx, DRHO=drho,
    H=maximum_smoothing_length, MU=viscosity, RHO0=rho0)

if device.communicator.rank == 0:
    print(f'Speed of sound: {c:.4f} m/s  ({cond})')

dt, dt_cond = model.compute_dt(
    LREF=lref, UREF=U0, DX=dx, DRHO=drho,
    H=maximum_smoothing_length, MU=viscosity, RHO0=rho0)

# Override: need enough time resolution to capture oscillation (< T/50)
dt_osc = T_period / 50.0
if dt_osc < dt:
    dt = dt_osc
if device.communicator.rank == 0:
    print(f'Timestep: {dt:.3e} s  ({dt_cond}, capped to T/50)')

steps_per_period = max(1, int(np.round(T_period / dt)))
dt_actual = T_period / steps_per_period  # re-synchronise to exact period

# ─── Integrator ──────────────────────────────────────────────────────────────
integrator = hoomd.sph.Integrator(dt=dt_actual)
vvb = hoomd.sph.methods.VelocityVerletBasic(filter=filterfluid,
                                              densitymethod='SUMMATION')
integrator.methods.append(vvb)
integrator.forces.append(model)
sim.operations.integrator = integrator

# ─── Output ──────────────────────────────────────────────────────────────────
gsd_writer = hoomd.write.GSD(filename=dumpname,
                              trigger=hoomd.trigger.Periodic(steps_per_period // 10),
                              mode='wb',
                              dynamic=['property', 'momentum'])
sim.operations.writers.append(gsd_writer)

logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps', 'walltime'])
compute_fluid = hoomd.sph.compute.SinglePhaseFlowBasicProperties(filter=filterfluid)
sim.operations.computes.append(compute_fluid)
logger.add(compute_fluid, quantities=['e_kin_fluid', 'mean_density'])
table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(steps_per_period // 10),
                          logger=logger, max_header_len=10)
sim.operations.writers.append(table)

log_file = open(logname, mode='w+', newline='\n')
table_file = hoomd.write.Table(output=log_file,
                                trigger=hoomd.trigger.Periodic(steps_per_period // 10),
                                logger=logger, max_header_len=10)
sim.operations.writers.append(table_file)

# ─── Run: oscillate bottom plate, one period at a time ───────────────────────
if device.communicator.rank == 0:
    print(f'Starting Stokes 2nd run: {n_periods} periods × {steps_per_period} steps/period'
          f' at {dt_string}')

sim.run(0, write_at_start=True)

for period_i in range(n_periods):
    # Update bottom plate velocity at start of each period step-by-step
    for step_i in range(steps_per_period):
        t = (period_i * steps_per_period + step_i) * dt_actual
        v_wall = U0 * np.cos(omega * t)
        with sim.state.cpu_local_snapshot as snap:
            # Recompute mask each step: local particle count changes as fluid
            # particles migrate between MPI ranks after sim.run(1).
            pos_local = np.array(snap.particles.position)
            tid_local = np.array(snap.particles.typeid)
            bot_mask  = (tid_local == 1) & (pos_local[:, 1] < -0.5 * lref)
            snap.particles.velocity[bot_mask, 0] = np.float32(v_wall)
        sim.run(1)

gsd_writer.flush()

# ─── Post-processing: Stokes layer amplitude profile ─────────────────────────
if device.communicator.rank == 0:
    with gsd.hoomd.open(dumpname, 'r') as traj:
        snap = traj[-1]
        pos = snap.particles.position
        tid = snap.particles.typeid
        vel = snap.particles.velocity

    # Extract fluid-only velocity at omega*t = 0 (cosine maximum)
    # Last period ends at t = n_periods * T_period → phase = 0
    fluid   = (tid == 0)
    y_f     = pos[fluid, 1]
    vx_f    = vel[fluid, 0]
    yp_f    = y_f + 0.5 * lref  # distance from oscillating plate

    # Analytical amplitude (velocity envelope at $\omega t = 0$): $A(y') = U_0 \exp(-y'/\delta)$
    inside  = yp_f < lref * 0.9  # exclude near-top region (influenced by upper plate)
    yp_ev   = yp_f[inside]
    vx_ev   = vx_f[inside]
    va_ev   = U0 * np.exp(-yp_ev / delta)

    L2_err = np.sqrt(np.mean((vx_ev - va_ev)**2)) / U0 * 100.0

    print(f'\n── Stokes 2nd check (after {n_periods} periods, step {snap.configuration.step}) ──')
    print(f'  Stokes layer delta = {delta:.4e} m')
    print(f'  L₂ velocity error vs analytical amplitude = {L2_err:.2f} %')
    print(f'  (Compare: delta/dx = {delta/dx:.2f} — at least 3 points per delta needed)')

if HAS_VTU and device.communicator.rank == 0:
    export_gsd2vtu.export_spf(dumpname)

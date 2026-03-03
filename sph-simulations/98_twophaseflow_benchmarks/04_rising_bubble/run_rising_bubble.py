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

Rising bubble — WCSPH run script.

BENCHMARK DESCRIPTION
---------------------
A light spherical gas bubble ($\rho_2 = 100$ kg/m³) rises through a heavier liquid
($\rho_1 = 1000$ kg/m³) under gravity.  Surface tension ($\sigma = 0.05$ N/m) keeps the
bubble coherent.  Solid walls bound y; x and z are periodic.

Key parameters:
  Eötvös   $Eo = (\rho_1 - \rho_2) \, g \, D^2 / \sigma \approx 28.2$
  Morton   $Mo = g \, \mu_1^4 \, (\rho_1 - \rho_2) / (\rho_1^2 \, \sigma^3) \approx 706$
  Density ratio  $\rho_1 / \rho_2 = 10$
  Viscosity ratio  $\mu_1 / \mu_2 = 10$

The terminal rise velocity and bubble shape are the primary validation
quantities.  Compare against experimental data (e.g. Grace 1973 diagram) or
reference numerical solutions (e.g. Adami et al. 2010).

Usage:
    python3 run_rising_bubble.py <num_length> <init_gsd_file>
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

dt_string = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
logname   = filename.replace('_init.gsd', '_run.log')
dumpname  = filename.replace('_init.gsd', '_run.gsd')

sim.create_state_from_gsd(filename=filename)

# ─── Physical parameters ─────────────────────────────────────────────────────
lref       = 0.001        # reference length (domain width)   [m]
R_bub      = 0.2 * lref   # bubble radius                     [m]
dx         = lref / num_length
rho01      = 1000.0       # rest density liquid 'W'           [kg/m³]
rho02      = 100.0        # rest density gas    'N'           [kg/m³]
viscosity1 = 0.01         # dynamic viscosity liquid          [Pa·s]
viscosity2 = 0.001        # dynamic viscosity gas             [Pa·s]
sigma      = 5e-5         # surface tension                   [N/m]
                          # σ = 5e-5 N/m gives Eo = 28.2 with D = 0.4 mm, g = 9.81 m/s²
                          # (the previous value of 0.05 N/m gave Eo = 0.028, so buoyancy
                          #  was overwhelmed by CSF parasitic currents and the bubble did
                          #  not rise)
gy         = -9.81        # gravitational acceleration        [m/s²]
backpress  = 0.01         # background pressure coeff         [–]
drho       = 0.01         # allowed density variation         [–]
# Steps to reach terminal velocity:
#   τ_relax = ρ₂ 2R² / (9 μ₁) ≈ 89 µs → 5τ ≈ 444 µs ≈ 700 steps (high-res).
#   For a clean U_T measurement the last third of the run must be at steady rise;
#   50 001 steps → ~31.6 ms physical time (high-res, dx=25 µm) → bubble rises
#   ≈ 0.18 mm ≈ 0.45 D, which gives a well-resolved terminal-velocity window.
steps      = int(sys.argv[3]) if len(sys.argv) > 3 else 20001  # simulation steps

# Estimate reference velocity from buoyancy (Hadamard-Rybczynski upper bound):
# $U_\mathrm{ref} \sim \sqrt{(\rho_1 - \rho_2) \, g \, R / \rho_1}$  (rough estimate)
g_mag  = abs(gy)
drho_b = rho01 - rho02
U_ref  = np.sqrt(drho_b * g_mag * R_bub / rho01)  # rough estimate

# ─── Kernel ──────────────────────────────────────────────────────────────────
kernel     = 'WendlandC4'
slength    = hoomd.sph.kernel.OptimalH[kernel] * dx
rcut       = hoomd.sph.kernel.Kappa[kernel] * slength
kernel_obj = hoomd.sph.kernel.Kernels[kernel]()
kappa      = kernel_obj.Kappa()

nlist = hoomd.nsearch.nlist.Cell(buffer=rcut * 0.05,
                                  rebuild_check_delay=1, kappa=kappa)

# ─── EOS (Tait) ──────────────────────────────────────────────────────────────
eos1 = hoomd.sph.eos.Tait()
eos2 = hoomd.sph.eos.Tait()
eos1.set_params(rho01, backpress)
eos2.set_params(rho02, backpress)

# ─── Particle filters ────────────────────────────────────────────────────────
filterfluidW = hoomd.filter.Type(['W'])   # liquid (phase 1)
filterfluidN = hoomd.filter.Type(['N'])   # gas bubble (phase 2)
filtersolid  = hoomd.filter.Type(['S'])   # solid walls (top + bottom)

# ─── TwoPhaseFlow model ──────────────────────────────────────────────────────
model = hoomd.sph.sphmodel.TwoPhaseFlow(
    kernel=kernel_obj,
    eos1=eos1, eos2=eos2,
    nlist=nlist,
    fluidgroup1_filter=filterfluidW,
    fluidgroup2_filter=filterfluidN,
    solidgroup_filter=filtersolid,
    densitymethod='SUMMATION',
    colorgradientmethod='DENSITYRATIO',
)

model.mu1              = viscosity1
model.mu2              = viscosity2
model.sigma12          = sigma
model.omega            = 90          # neutral wetting at solid walls
model.gy               = gy          # gravity in y direction
model.damp             = 1000
model.artificialviscosity = True
model.alpha            = 0.2
model.beta             = 0.0
model.densitydiffusion = False

# ─── Speed of sound & timestep ───────────────────────────────────────────────
maximum_smoothing_length = sph_helper.set_max_sl(sim, device, model)

c1, cond1, c2, cond2 = model.compute_speedofsound(
    LREF=lref, UREF=U_ref, DX=dx, DRHO=drho,
    H=maximum_smoothing_length,
    MU1=viscosity1, MU2=viscosity2,
    RHO01=rho01, RHO02=rho02, SIGMA12=sigma)

if device.communicator.rank == 0:
    print(f'Phase W speed of sound: {c1:.4f} m/s  ({cond1})')
    print(f'Phase N speed of sound: {c2:.4f} m/s  ({cond2})')

sph_helper.update_min_c0_tpf(device, model, c1, c2,
                              mode='plain', lref=lref, uref=U_ref, cfactor=10.0)

dt, dt_cond = model.compute_dt(
    LREF=lref, UREF=U_ref, DX=dx, DRHO=drho,
    H=maximum_smoothing_length,
    MU1=viscosity1, MU2=viscosity2,
    RHO01=rho01, RHO02=rho02, SIGMA12=sigma)

if device.communicator.rank == 0:
    print(f'Timestep: {dt:.3e} s  ({dt_cond})')
    D   = 2 * R_bub
    Eo  = drho_b * g_mag * D**2 / sigma
    print(f'Eötvös number Eo = {Eo:.1f}')
    print(f'Density ratio    = {rho01/rho02:.0f}:1')

# ─── Integrator ──────────────────────────────────────────────────────────────
integrator = hoomd.sph.Integrator(dt=dt)
vvbW = hoomd.sph.methods.VelocityVerletBasic(filter=filterfluidW,
                                              densitymethod='SUMMATION')
vvbN = hoomd.sph.methods.VelocityVerletBasic(filter=filterfluidN,
                                              densitymethod='SUMMATION')
integrator.methods.append(vvbW)
integrator.methods.append(vvbN)
integrator.forces.append(model)
sim.operations.integrator = integrator

# ─── Output ──────────────────────────────────────────────────────────────────
gsd_writer = hoomd.write.GSD(filename=dumpname,
                              trigger=hoomd.trigger.Periodic(200),
                              mode='wb',
                              dynamic=['property', 'momentum'])
sim.operations.writers.append(gsd_writer)

logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps', 'walltime'])
compute_W = hoomd.sph.compute.SinglePhaseFlowBasicProperties(filter=filterfluidW)
compute_N = hoomd.sph.compute.SinglePhaseFlowBasicProperties(filter=filterfluidN)
sim.operations.computes.append(compute_W)
sim.operations.computes.append(compute_N)
logger.add(compute_W, quantities=['e_kin_fluid', 'mean_density'])
logger.add(compute_N, quantities=['e_kin_fluid', 'mean_density'])
table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(200), logger=logger,
                          max_header_len=10)
sim.operations.writers.append(table)

log_file = open(logname, mode='w+', newline='\n')
table_file = hoomd.write.Table(output=log_file,
                                trigger=hoomd.trigger.Periodic(200),
                                logger=logger, max_header_len=10)
sim.operations.writers.append(table_file)

# ─── Run ─────────────────────────────────────────────────────────────────────
if device.communicator.rank == 0:
    print(f'Starting WCSPH rising-bubble run at {dt_string}')
sim.run(steps, write_at_start=True)
gsd_writer.flush()

# ─── Post-processing: bubble centroid trajectory ──────────────────────────────
if device.communicator.rank == 0:
    y_cg = []
    t_vals = []
    with gsd.hoomd.open(dumpname, 'r') as traj:
        for snap in traj:
            tid = snap.particles.typeid
            pos = snap.particles.position
            bub = (tid == 1)
            if bub.any():
                y_cg.append(float(np.mean(pos[bub, 1])))
                t_vals.append(snap.configuration.step * dt)
        snap_last = traj[-1]

    y_cg   = np.array(y_cg)
    t_vals = np.array(t_vals)
    if len(t_vals) >= 2:
        # Estimate terminal velocity from last third of the run
        n_last = max(2, len(t_vals) // 3)
        U_T    = np.mean(np.gradient(y_cg[-n_last:], t_vals[-n_last:]))
    else:
        U_T = float('nan')

    print(f'\n── Rising bubble summary (last frame, step {snap_last.configuration.step}) ──')
    print(f'  Bubble centroid y = {y_cg[-1]*1e3:.3f} mm')
    print(f'  Estimated U_T     = {U_T:.4f} m/s')

if HAS_VTU and device.communicator.rank == 0:
    export_gsd2vtu.export_tpf(dumpname)

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

Capillary rise — WCSPH TwoPhaseFlow run script, GPU variant.

Identical physics and post-processing to run_capillary_rise.py; only the
execution device differs (hoomd.device.GPU, single rank).

BENCHMARK DESCRIPTION
---------------------
A wetting liquid ('W') fills the lower region of a 3-D square-box domain that
contains a vertical round capillary tube (solid ring 'S').  Gas ('N') fills
the remaining space above the initial flat liquid–gas interface.

Gravity acts in the $-y$ direction.  Surface tension $\sigma$ and the prescribed
contact angle $\theta$ (via model.omega) drive capillary rise ($\theta < 90°$) or
capillary depression ($\theta > 90°$) of the meniscus inside the tube.

VALIDATION — Jurin's law
  $h_\mathrm{Jurin} = 2 \sigma \cos(\theta) / (\rho_1 g R_\mathrm{cap})$

The measured quantity is:
  h_meas = y_meniscus_inside $-$ y_reservoir_surface_outside
         $\approx$ h_Jurin  at steady state (independent of reservoir size)

PHYSICAL PARAMETERS
  Read from capillary_rise_params.txt for the selected case_id.
  Defaults: $\sigma$=0.01 N/m, $\rho_1$=1000 kg/m$^3$, $\rho_2$=100 kg/m$^3$,
            $\mu_1$=0.1 $\mathrm{Pa \cdot s}$, $\mu_2$=0.001 $\mathrm{Pa \cdot s}$, g=9.81 m/s$^2$.

Usage:
    mpirun -np 1 python3 run_capillary_rise_gpu.py <num_length> <init_gsd> [case_id] [steps]
    (submit with caprise_gpu.sh on wolfgang; SPH_DEVICE=cpu runs the same script on the CPU)
      num_length : particles across R_cap (must match the geometry file)
      init_gsd   : GSD file from create_capillary_geometry.py
      case_id    : row index in capillary_rise_params.txt  (default: 0)
      steps      : simulation steps  (default: 50001)
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
# GPU variant: the whole two-phase model (forces, normals, CSF, particle
# shifting) runs on the device; the integrator still runs on the host.
# SPH_DEVICE=cpu switches back to the CPU for a one-to-one comparison run.
if os.environ.get('SPH_DEVICE', 'gpu').lower() == 'cpu':
    device = hoomd.device.CPU(notice_level=2)
else:
    device = hoomd.device.GPU(notice_level=2)
sim    = hoomd.Simulation(device=device)

num_length = int(sys.argv[1])
filename   = str(sys.argv[2])
case_id    = int(sys.argv[3])   if len(sys.argv) > 3 else 0
steps      = int(sys.argv[4])   if len(sys.argv) > 4 else 50001

# ─── Read parameter set from params file ─────────────────────────────────────
params_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'capillary_rise_params.txt')

def read_params(fname, cid):
    with open(fname) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # strip inline comments
            line = line.split('#')[0].strip()
            if not line:
                continue
            parts = line.split()
            if int(parts[0]) == cid:
                return dict(
                    theta_deg = float(parts[1]),
                    sigma     = float(parts[2]),
                    rho1      = float(parts[3]),
                    rho2      = float(parts[4]),
                    mu1       = float(parts[5]),
                    mu2       = float(parts[6]),
                )
    raise ValueError(f'Case {cid} not found in {fname}')

params    = read_params(params_file, case_id)
theta_deg = params['theta_deg']
sigma     = params['sigma']
rho01     = params['rho1']
rho02     = params['rho2']
mu1       = params['mu1']
mu2       = params['mu2']
theta_rad = np.radians(theta_deg)

# ─── File names ───────────────────────────────────────────────────────────────
dt_string = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
base      = filename.replace('_init.gsd', '')
_tag = ('_gpu_K' + os.environ.get('SPH_NORMAL_RELAX_ITERS', '1')
        + ('_cip' if os.environ.get('SPH_CIP', '0') == '1' else '')
        + ('_noav' if os.environ.get('SPH_AV', '1') == '0' else '')
        + (('_a' + os.environ['SPH_AV_ALPHA']) if 'SPH_AV_ALPHA' in os.environ else '')
        + (('_cf' + os.environ['SPH_CFACTOR']) if 'SPH_CFACTOR' in os.environ else '')
        + '_n%d_run' % steps)
logname   = f'{base}_case{case_id}_theta{int(theta_deg):03d}' + _tag + '.log'
dumpname  = f'{base}_case{case_id}_theta{int(theta_deg):03d}' + _tag + '.gsd'

sim.create_state_from_gsd(filename=filename)

# ─── Physical parameters ─────────────────────────────────────────────────────
R_cap      = 0.001                      # inner tube radius = reference length [m]
dx         = R_cap / num_length          # particle spacing                     [m]
gy         = -9.81                       # gravitational acceleration            [m/s²]
backpress  = 0.01                        # background pressure coefficient       [–]
drho       = 0.01                        # allowed density variation             [–]
n_wall     = 2                           # wall thickness (matches geometry)
R_inner    = R_cap
R_outer    = R_cap + n_wall * dx

# Analytical Jurin height
h_Jurin = 2.0 * sigma * np.cos(theta_rad) / (rho01 * abs(gy) * R_inner)

# Reference velocity: capillary velocity scale = $\sigma / (\rho_1 R)$
# (upper bound on initial meniscus speed)
U_cap  = sigma / (rho01 * R_cap)          # capillary velocity scale  [m/s]
U_ref  = max(U_cap, 1e-8)                 # ensure nonzero refvel

if device.communicator.rank == 0:
    print(f'─── Capillary Rise  case={case_id}  θ={theta_deg:.0f}°  ───')
    print(f'  σ={sigma:.4f} N/m,  ρ₁={rho01:.0f} kg/m³,  ρ₂={rho02:.0f} kg/m³')
    print(f'  μ₁={mu1:.4f} Pa·s,  μ₂={mu2:.4f} Pa·s')
    print(f'  h_Jurin = 2σcos(θ)/(ρ₁gR) = {h_Jurin*1e3:+.3f} mm  '
          f'({"+rise" if h_Jurin > 0 else ("neutral" if h_Jurin == 0 else "depression")})')
    print(f'  U_cap = σ/(ρ₁R) = {U_cap:.4f} m/s')

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
filterfluidW = hoomd.filter.Type(['W'])   # liquid  (phase 1)
filterfluidN = hoomd.filter.Type(['N'])   # gas     (phase 2)
filtersolid  = hoomd.filter.Type(['S'])   # solid walls + tube ring

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

model.mu1              = mu1
model.mu2              = mu2
model.sigma12          = sigma
model.omega            = theta_deg      # contact angle [°] at solid–liquid interface
model.gy               = gy             # gravity in −y direction
model.damp             = 2000           # body-force ramp steps (does not affect steady state)
# SPH_AV=0 switches Monaghan artificial viscosity off (creeping flow: the
# physical viscosity is sufficient; AV adds ~rho*alpha*c*h/8 ~ 0.2 Pa s here).
model.artificialviscosity = os.environ.get('SPH_AV', '1') != '0'
model.alpha            = float(os.environ.get('SPH_AV_ALPHA', '0.2'))   # Monaghan AV linear coefficient
model.beta             = 0.0
model.densitydiffusion = False
# Normal-relaxation sweeps per step (curvature bottleneck study); K=1 = original.
normal_relax_iters = int(os.environ.get('SPH_NORMAL_RELAX_ITERS', '1'))
# Consistent interface pressure (Hu & Adams 2009) on cross-phase pairs.
use_cip = os.environ.get('SPH_CIP', '0') == '1'
model.consistent_interface_pressure = use_cip

# ─── Speed of sound & timestep ───────────────────────────────────────────────
maximum_smoothing_length = sph_helper.set_max_sl(sim, device, model)

c1, cond1, c2, cond2 = model.compute_speedofsound(
    LREF=R_cap, UREF=U_ref, DX=dx, DRHO=drho,
    H=maximum_smoothing_length,
    MU1=mu1, MU2=mu2,
    RHO01=rho01, RHO02=rho02, SIGMA12=sigma)

if device.communicator.rank == 0:
    print(f'  Phase W speed of sound: {c1:.4f} m/s  ({cond1})')
    print(f'  Phase N speed of sound: {c2:.4f} m/s  ({cond2})')

cfactor_used = float(os.environ.get('SPH_CFACTOR', '10.0'))
sph_helper.update_min_c0_tpf(device, model, c1, c2,
                              mode='plain', lref=R_cap, uref=U_ref,
                              cfactor=cfactor_used)   # sound-speed safety factor

dt, dt_cond = model.compute_dt(
    LREF=R_cap, UREF=U_ref, DX=dx, DRHO=drho,
    H=maximum_smoothing_length,
    MU1=mu1, MU2=mu2,
    RHO01=rho01, RHO02=rho02, SIGMA12=sigma)

if device.communicator.rank == 0:
    t_sim = steps * dt
    t_LW  = h_Jurin**2 * 2.0 * mu1 / max(R_inner * sigma * abs(np.cos(theta_rad)), 1e-20)
    print(f'  Timestep: {dt:.3e} s  ({dt_cond})')
    print(f'  Simulated time: {t_sim:.3e} s')
    print(f'  Lucas-Washburn equilib. time (est.): {t_LW:.3e} s')
    if t_sim < 3 * t_LW and abs(np.cos(theta_rad)) > 0.1:
        print(f'  WARNING: simulated time may be insufficient for full equilibration.'
              f'  Consider increasing steps.')

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
if normal_relax_iters > 1:
    sim.run(0)                                   # attach, then configure the C++ model
    model.activateNormalRelaxation(normal_relax_iters)
    if device.communicator.rank == 0:
        print(f'  Normal relaxation: {normal_relax_iters} sweeps per step')
if use_cip and device.communicator.rank == 0:
    print('  Consistent interface pressure: ON')

# ─── Output ──────────────────────────────────────────────────────────────────
gsd_period = max(1, steps // 200)
log_period = max(1, steps // 1000)

# Remove any stale output GSD left by a previous crashed run.
# ALL ranks attempt the removal so each node's NFS metadata cache is flushed;
# the first rank to run succeeds, the rest get FileNotFoundError (ignored).
try:
    os.remove(dumpname)
except (FileNotFoundError, OSError):
    pass
device.communicator.barrier()

gsd_writer = hoomd.write.GSD(filename=dumpname,
                              trigger=hoomd.trigger.Periodic(gsd_period),
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

table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(log_period * 5),
                          logger=logger, max_header_len=10)
sim.operations.writers.append(table)

log_file  = open(logname, mode='w+', newline='\n')
table_log = hoomd.write.Table(output=log_file,
                              trigger=hoomd.trigger.Periodic(log_period),
                              logger=logger, max_header_len=10)
sim.operations.writers.append(table_log)

# ─── Run ─────────────────────────────────────────────────────────────────────
if device.communicator.rank == 0:
    print(f'Starting capillary rise run (θ={theta_deg:.0f}°, case {case_id}) at {dt_string} on {device.__class__.__name__}')
sim.run(steps, write_at_start=True)
gsd_writer.flush()

# ─── Post-processing: measure equilibrium rise height ────────────────────────
if device.communicator.rank == 0:
    # Work on the last GSD frame
    with gsd.hoomd.open(dumpname, 'r') as traj:
        snap = traj[-1]

    pos = snap.particles.position
    tid = snap.particles.typeid
    step_last = snap.configuration.step

    x_all = pos[:, 0]
    y_all = pos[:, 1]
    z_all = pos[:, 2]
    r_all = np.sqrt(x_all**2 + z_all**2)

    # W particles strictly inside the tube (avoid tube-wall neighbours)
    mask_W_inside = (tid == 0) & (r_all < R_inner - 0.5 * dx)

    # W particles in the reservoir (outside the tube outer wall + gap)
    mask_W_outside = (tid == 0) & (r_all > R_outer + 0.5 * dx)

    n_W_inside  = int(np.sum(mask_W_inside))
    n_W_outside = int(np.sum(mask_W_outside))

    if n_W_inside > 4 and n_W_outside > 4:
        # Interface positions: 97th percentile of y for W particles
        # (robust against a few outlier particles near the solid walls)
        y_meniscus  = float(np.percentile(pos[mask_W_inside,  1], 97))
        y_reservoir = float(np.percentile(pos[mask_W_outside, 1], 97))
        h_meas      = y_meniscus - y_reservoir

        rel_err = ((h_meas - h_Jurin) / max(abs(h_Jurin), dx)) * 100.0

        print(f'\n── Capillary rise result (case {case_id}, θ={theta_deg:.0f}°, '
              f'step {step_last}) ──')
        print(f'  h_Jurin (analytical)  = {h_Jurin*1e3:+7.3f} mm')
        print(f'  y_meniscus  (inside)  = {y_meniscus*1e3:+7.3f} mm  '
              f'(from {n_W_inside} W particles inside tube)')
        print(f'  y_reservoir (outside) = {y_reservoir*1e3:+7.3f} mm  '
              f'(from {n_W_outside} W particles in reservoir)')
        print(f'  h_meas = meniscus − reservoir = {h_meas*1e3:+7.3f} mm')
        if abs(h_Jurin) > dx:
            print(f'  Relative error = {rel_err:+.1f} %')
        else:
            print(f'  (θ≈90°: Jurin height ≈ 0, absolute residual = {abs(h_meas)*1e3:.4f} mm)')

        # Write one-line summary to a results file
        summary_file = os.path.join(os.path.dirname(os.path.abspath(filename)),
                                    'capillary_rise_summary_gpu.dat')
        write_header = not os.path.isfile(summary_file)
        with open(summary_file, 'a') as sf:
            if write_header:
                sf.write('# case_id  theta_deg  h_Jurin[m]  h_meas[m]  err_pct  sigma  mu1  K_relax  cip  steps  av  alpha  cfactor\n')
            sf.write(f'{case_id:8d}  {theta_deg:9.1f}  {h_Jurin:10.4e}  '
                     f'{h_meas:10.4e}  {rel_err:8.2f}  {sigma:6.4f}  {mu1:.4f}  {normal_relax_iters:d}  {int(use_cip):d}  {steps:d}  {int(model.artificialviscosity):d}  {model.alpha:g}  {cfactor_used:g}\n')
        print(f'  Summary appended to: {summary_file}')
    else:
        print(f'  WARNING: too few W particles inside ({n_W_inside}) or '
              f'outside ({n_W_outside}) tube for measurement.')

if HAS_VTU and device.communicator.rank == 0:
    export_gsd2vtu.export_tpf(dumpname)

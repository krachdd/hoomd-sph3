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

Sessile Droplet Shear + Snap-back — three-phase WCSPH run script.

BENCHMARK DESCRIPTION
---------------------
The simulation runs in three consecutive phases without restart:

Phase 1 — Relaxation
  A cubic patch of fluid 'W' (droplet) sits on the stationary bottom solid
  wall, surrounded by fluid 'N' (ambient).  Surface tension (σ = 0.01 N/m)
  and the prescribed contact angle (θ = 90°) reshape the cube into a sessile
  hemispherical droplet.  Both solid walls are stationary.  Gravity (gy)
  acts in the −y direction; with equal densities (ρ₁ = ρ₂) the Bond number
  Bo ≈ 0.054, so the equilibrium shape remains close to a hemisphere.

Phase 2 — Shear
  After relaxation the top solid wall is set to velocity U_wall in the
  x direction (Couette-like driving).  The resulting shear flow deforms the
  sessile droplet.  Contact-angle hysteresis (θ_rec=75°, θ_adv=105°) pins
  the contact line, so the droplet tilts but does not slide.  The extent of
  deformation is characterised by the capillary number Ca = μ U_wall / σ = 0.001.

Phase 3 — Snap-back
  The top wall velocity is reset to zero.  With hysteresis pinning active the
  droplet elastically recovers toward its relaxed position.  The residual
  x-displacement after snap-back (relative to the shear-induced drift) is used
  as a recovery metric.

Domain:
  x : 4·lref    periodic
  y : 2·lref     solid walls  (Adami 2012 no-slip BC, 3 layers each side)
  z : 2·lref     periodic

Key dimensionless numbers:
  Bond number  Bo = ρ g R² / σ  ≈ 1.0    (gravity and surface tension comparable)
  Capillary    Ca = μ U_wall / σ = 0.001  (surface-tension dominated)
  Reynolds     Re = ρ U_wall lref / μ      = 10

Step budget (dt ≈ 1.25 µs):
  Capillary time  τ_cap = μ R / σ  ≈  80 steps
  Relaxation default (20001 steps) ≈ 250 τ_cap  → well-converged sessile shape
  Shear strain    1/γ̇ = H/U_wall  ≈ 160 000 steps
  Shear default  (50001 steps)     ≈ 0.3 shear strains → observable deformation
  Snap-back default (50001 steps)  ≈ same duration as shear phase

Usage:
    python3 run_sessile_droplet_snapback.py <num_length> <init_gsd_file> \\
            [steps_relax] [steps_shear] [steps_snapback] [gsd_period]

    num_length    : resolution used to create the init file (e.g. 20)
    init_gsd_file : path to the *_init.gsd produced by create_input_geometry.py
    steps_relax   : relaxation steps        (default 20001)
    steps_shear   : shear steps             (default 50001)
    steps_snapback: snap-back steps         (default 50001)
    gsd_period    : GSD write interval      (default 200)
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

num_length    = int(sys.argv[1])
filename      = str(sys.argv[2])
steps_relax   = int(sys.argv[3]) if len(sys.argv) > 3 else 20001
steps_shear   = int(sys.argv[4]) if len(sys.argv) > 4 else 50001
steps_snapback = int(sys.argv[5]) if len(sys.argv) > 5 else 50001
gsd_period    = int(sys.argv[6]) if len(sys.argv) > 6 else 200

dt_string = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
logname   = filename.replace('_init.gsd', '_run.log')
dumpname  = filename.replace('_init.gsd', '_run.gsd')

sim.create_state_from_gsd(filename=filename)

# ─── Physical parameters ─────────────────────────────────────────────────────
lref       = 0.001        # reference length                              [m]
dx         = lref / num_length
rho01      = 1000.0       # rest density phase W (droplet)                [kg/m³]
rho02      = 1000.0       # rest density phase N (ambient)                [kg/m³]
viscosity1 = 0.001        # dynamic viscosity phase W                     [Pa·s]
viscosity2 = 0.001        # dynamic viscosity phase N                     [Pa·s]
sigma      = 0.01         # surface tension                               [N/m]
theta_eq   = 90           # equilibrium contact angle with solid wall     [°]
theta_adv  = 105.0        # advancing contact angle for hysteresis        [°]
theta_rec  =  75.0        # receding  contact angle for hysteresis        [°]
gy         = -9.81        # gravitational acceleration (−y direction)     [m/s²]
backpress  = 0.01         # background pressure coefficient               [–]
drho       = 0.01         # allowed density variation                     [–]
U_wall     = 0.01         # top wall velocity in shear phase              [m/s]

H_flu          = 2 * lref                                   # fluid channel height  [m]
R_drop_target  = H_flu / 2                                  # = lref               [m]
a_cube         = (2 * np.pi * R_drop_target**3 / 3) ** (1/3)  # ≈ 1.28·lref       [m]

# Derived dimensionless numbers (for information)
n_cube = int(round(a_cube / dx))
R_hemi = (3 * n_cube**3 * dx**3 / (2 * np.pi)) ** (1 / 3)
Bo = rho01 * abs(gy) * R_hemi**2 / sigma
Ca = viscosity1 * U_wall / sigma
Re = rho01 * U_wall * lref / viscosity1

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
filterfluidW = hoomd.filter.Type(['W'])   # droplet phase
filterfluidN = hoomd.filter.Type(['N'])   # ambient phase
filtersolid  = hoomd.filter.Type(['S'])   # solid walls (bottom + top)

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
model.omega            = theta_eq     # equilibrium / nominal contact angle [°]
model.omega_adv        = theta_adv    # advancing contact angle for hysteresis [°]
model.omega_rec        = theta_rec    # receding  contact angle for hysteresis [°]
model.hysteresis       = True         # enable contact-angle hysteresis
model.gx               = 0.0
model.gy               = gy
model.damp             = 1000
model.artificialviscosity = True
model.alpha            = 0.2
model.beta             = 0.0
model.densitydiffusion = False

# ─── Speed of sound & timestep ───────────────────────────────────────────────
maximum_smoothing_length = sph_helper.set_max_sl(sim, device, model)

c1, cond1, c2, cond2 = model.compute_speedofsound(
    LREF=lref, UREF=U_wall, DX=dx, DRHO=drho,
    H=maximum_smoothing_length,
    MU1=viscosity1, MU2=viscosity2,
    RHO01=rho01, RHO02=rho02, SIGMA12=sigma)

if device.communicator.rank == 0:
    print(f'Phase W speed of sound: {c1:.4f} m/s  ({cond1})')
    print(f'Phase N speed of sound: {c2:.4f} m/s  ({cond2})')

sph_helper.update_min_c0_tpf(device, model, c1, c2,
                              mode='plain', lref=lref, uref=U_wall, cfactor=10.0)

dt, dt_cond = model.compute_dt(
    LREF=lref, UREF=U_wall, DX=dx, DRHO=drho,
    H=maximum_smoothing_length,
    MU1=viscosity1, MU2=viscosity2,
    RHO01=rho01, RHO02=rho02, SIGMA12=sigma)

if device.communicator.rank == 0:
    print(f'Timestep: {dt:.3e} s  ({dt_cond})')
    print(f'Bond number  Bo = {Bo:.4f}')
    print(f'Capillary    Ca = {Ca:.4f}')
    print(f'Reynolds     Re = {Re:.1f}')
    print(f'R_sessile (hemisphere) ≈ {R_hemi*1e3:.3f} mm')

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
table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(200), logger=logger,
                          max_header_len=10)
sim.operations.writers.append(table)

log_file = open(logname, mode='w+', newline='\n')
table_file = hoomd.write.Table(output=log_file,
                                trigger=hoomd.trigger.Periodic(200),
                                logger=logger, max_header_len=10)
sim.operations.writers.append(table_file)

# ─── Phase 1: Relaxation (stationary walls) ───────────────────────────────────
if device.communicator.rank == 0:
    print(f'\n── Phase 1: Relaxation  ({steps_relax} steps) ─────────────────────')
    print(f'   Both solid walls stationary.  Surface tension morphs the cube')
    print(f'   into a sessile hemispherical droplet (θ = {theta_eq}°).')
    print(f'   Starting at {dt_string}')

sim.run(steps_relax, write_at_start=True)
gsd_writer.flush()

# ─── Measure droplet shape after relaxation ───────────────────────────────────
cx_r = cy_r = cz_r = 0.0
if device.communicator.rank == 0:
    with gsd.hoomd.open(dumpname, 'r') as traj:
        snap_r = traj[-1]

    pos_r = snap_r.particles.position
    tid_r = snap_r.particles.typeid
    drop_r = (tid_r == 0)   # 'W' — droplet

    if drop_r.any():
        p = pos_r[drop_r]
        cx_r = float(np.mean(p[:, 0]))
        cy_r = float(np.mean(p[:, 1]))
        cz_r = float(np.mean(p[:, 2]))
        dy_r = float(p[:, 1].max() - p[:, 1].min())   # droplet height
        dx_r = float(p[:, 0].max() - p[:, 0].min())   # droplet extent in x
        dz_r = float(p[:, 2].max() - p[:, 2].min())   # droplet extent in z
        print(f'\n── Relaxation complete at step {snap_r.configuration.step} ──')
        print(f'   Droplet centroid   : ({cx_r*1e3:.4f}, {cy_r*1e3:.4f}, {cz_r*1e3:.4f}) mm')
        print(f'   Droplet height (y) : {dy_r*1e3:.4f} mm')
        print(f'   Droplet extent  x  : {dx_r*1e3:.4f} mm')
        print(f'   Droplet extent  z  : {dz_r*1e3:.4f} mm')
    else:
        print('   Warning: no droplet (W) particles found in final relaxation frame.')

# ─── Phase transition: activate top wall shear ───────────────────────────────
if device.communicator.rank == 0:
    print(f'\n── Phase transition: setting top wall velocity to {U_wall} m/s ──')

with sim.state.cpu_local_snapshot as snap:
    pos_l = np.asarray(snap.particles.position)
    tid_l = np.asarray(snap.particles.typeid)
    vel_l = np.asarray(snap.particles.velocity)
    # Top solid wall: type 'S' (index 2) and y >= H_flu/2
    top_mask = (tid_l == 2) & (pos_l[:, 1] >= H_flu / 2)
    vel_l[top_mask, 0] = U_wall   # set vx = U_wall; vy, vz remain zero

if device.communicator.rank == 0:
    print(f'   Top wall velocity updated.')

# ─── Phase 2: Shear (top wall moving at U_wall) ──────────────────────────────
if device.communicator.rank == 0:
    print(f'\n── Phase 2: Shear  ({steps_shear} steps) ──────────────────────────')
    print(f'   Top wall velocity U_wall = {U_wall} m/s.')
    print(f'   Capillary number Ca = {Ca:.4f}.')

sim.run(steps_shear)
gsd_writer.flush()

# ─── Post-processing: droplet shape after shear ───────────────────────────────
cx_s = cy_s = cz_s = 0.0
if device.communicator.rank == 0:
    with gsd.hoomd.open(dumpname, 'r') as traj:
        snap_s = traj[-1]

    pos_s = snap_s.particles.position
    tid_s = snap_s.particles.typeid
    drop_s = (tid_s == 0)   # 'W' — droplet

    if drop_s.any():
        p = pos_s[drop_s]
        cx_s = float(np.mean(p[:, 0]))
        cy_s = float(np.mean(p[:, 1]))
        cz_s = float(np.mean(p[:, 2]))
        dy_s = float(p[:, 1].max() - p[:, 1].min())
        dx_s = float(p[:, 0].max() - p[:, 0].min())
        dz_s = float(p[:, 2].max() - p[:, 2].min())
        # Deformation parameter D = (L - B) / (L + B)
        # L = longest axis, B = shortest axis (in x-y plane)
        L = max(dx_s, dy_s)
        B = min(dx_s, dy_s)
        D_param = (L - B) / (L + B) if (L + B) > 0 else 0.0
        delta_x_shear = cx_s - cx_r
        print(f'\n── Shear complete at step {snap_s.configuration.step} ──')
        print(f'   Droplet centroid   : ({cx_s*1e3:.4f}, {cy_s*1e3:.4f}, {cz_s*1e3:.4f}) mm')
        print(f'   Droplet height (y) : {dy_s*1e3:.4f} mm')
        print(f'   Droplet extent  x  : {dx_s*1e3:.4f} mm')
        print(f'   Droplet extent  z  : {dz_s*1e3:.4f} mm')
        print(f'   Deformation D      : {D_param:.4f}  (0 = undeformed, 1 = fully elongated)')
        print(f'   Centroid x-drift   : Δx = {delta_x_shear*1e3:.4f} mm')
    else:
        print('   Warning: no droplet (W) particles found in final shear frame.')

# ─── Phase transition: deactivate top wall (snap-back) ───────────────────────
if device.communicator.rank == 0:
    print(f'\n── Phase transition: resetting top wall velocity to 0 m/s ──')

with sim.state.cpu_local_snapshot as snap:
    pos_l = np.asarray(snap.particles.position)
    tid_l = np.asarray(snap.particles.typeid)
    vel_l = np.asarray(snap.particles.velocity)
    top_mask = (tid_l == 2) & (pos_l[:, 1] >= H_flu / 2)
    vel_l[top_mask, 0] = 0.0

if device.communicator.rank == 0:
    print(f'   Top wall velocity reset to zero.')

# ─── Phase 3: Snap-back (top wall stationary again) ──────────────────────────
if device.communicator.rank == 0:
    print(f'\n── Phase 3: Snap-back  ({steps_snapback} steps) ────────────────────')
    print(f'   Top wall stopped.  Hysteresis pinning allows elastic recovery.')

sim.run(steps_snapback)
gsd_writer.flush()

# ─── Post-processing: droplet shape after snap-back ───────────────────────────
if device.communicator.rank == 0:
    with gsd.hoomd.open(dumpname, 'r') as traj:
        snap_b = traj[-1]

    pos_b = snap_b.particles.position
    tid_b = snap_b.particles.typeid
    drop_b = (tid_b == 0)   # 'W' — droplet

    if drop_b.any():
        p = pos_b[drop_b]
        cx_b = float(np.mean(p[:, 0]))
        cy_b = float(np.mean(p[:, 1]))
        cz_b = float(np.mean(p[:, 2]))
        dy_b = float(p[:, 1].max() - p[:, 1].min())
        dx_b = float(p[:, 0].max() - p[:, 0].min())
        dz_b = float(p[:, 2].max() - p[:, 2].min())
        L = max(dx_b, dy_b)
        B = min(dx_b, dy_b)
        D_param_b = (L - B) / (L + B) if (L + B) > 0 else 0.0
        delta_x_shear    = cx_s - cx_r     # drift during shear
        delta_x_recovery = cx_b - cx_r     # residual after snap-back
        print(f'\n── Snap-back complete at step {snap_b.configuration.step} ──')
        print(f'   Droplet centroid   : ({cx_b*1e3:.4f}, {cy_b*1e3:.4f}, {cz_b*1e3:.4f}) mm')
        print(f'   Droplet height (y) : {dy_b*1e3:.4f} mm')
        print(f'   Droplet extent  x  : {dx_b*1e3:.4f} mm')
        print(f'   Droplet extent  z  : {dz_b*1e3:.4f} mm')
        print(f'   Deformation D      : {D_param_b:.4f}  (0 = undeformed, 1 = fully elongated)')
        print(f'   Centroid x-drift (shear)    : Δx_shear    = {delta_x_shear*1e3:.4f} mm')
        print(f'   Centroid x-drift (residual) : Δx_recovery = {delta_x_recovery*1e3:.4f} mm')
        if abs(delta_x_shear) > 0 and abs(delta_x_recovery) < abs(delta_x_shear) * 0.2:
            verdict = "Droplet fully recovered (hysteresis pinning effective)"
        else:
            verdict = "Partial recovery"
        print(f'   Verdict            : {verdict}')
    else:
        print('   Warning: no droplet (W) particles found in final snap-back frame.')

    print(f'\nOutput GSD : {dumpname}')
    print(f'Output log : {logname}')

if HAS_VTU and device.communicator.rank == 0:
    export_gsd2vtu.export_tpf(dumpname)

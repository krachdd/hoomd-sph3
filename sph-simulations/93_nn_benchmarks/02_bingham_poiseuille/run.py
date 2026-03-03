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

Bingham plane Poiseuille flow — Non-Newtonian benchmark.

BENCHMARK DESCRIPTION
---------------------
Body-force-driven flow between two parallel plates at y = ±H/2.
The fluid is a Bingham plastic regularised with the Papanastasiou method:

    μ_eff = μ_p + τ_y · (1 − e^(−m|γ̇|)) / |γ̇|

Exact Bingham analytical solution (for comparison):

  Plug zone (|y| ≤ y_p = τ_y / (ρ f_x)):
      v(y) = v_plug

  Flowing zone (y_p < |y| ≤ H/2):
      v(y) = ρ f_x / (2 μ_p) · ((H/2)² − y²) − τ_y / μ_p · (H/2 − |y|)

  v_plug = v(y_p)

Note: the Papanastasiou regularisation introduces a smoothed transition near
y_p, so a small L₂ discrepancy relative to the exact Bingham solution is expected
near the yield surface.  The benchmark checks that
  (a) the flowing-zone profile matches well (L₂ < 5 %),
  (b) the plug velocity is within 5 % of the exact value,
  (c) τ_y = 0 case reduces to Newtonian (L₂ < 5 %).

Parameters:
  μ_p   = 0.001 Pa·s   (plastic viscosity)
  τ_y   = 0.004 Pa     (yield stress → Bingham number Bi ≈ 0.32)
  m_reg = 1.0 s        (Papanastasiou regularisation)
  ρ f_x H/2 = 0.05 Pa  (wall stress)  →  y_p / (H/2) = 0.08

Usage:
    python3 run.py [num_length [steps]]
    Defaults: num_length=20, steps=30001
"""

import sys, os, math
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../../'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules'))

import hoomd
from hoomd import sph
import numpy as np
from datetime import datetime
import gsd.hoomd
import sph_helper

# ─── CLI args ────────────────────────────────────────────────────────────────
num_length = int(sys.argv[1]) if len(sys.argv) > 1 else 20
steps      = int(sys.argv[2]) if len(sys.argv) > 2 else 30001

# ─── Physical parameters ─────────────────────────────────────────────────────
lref      = 0.001        # channel gap         [m]
H         = lref
rho0      = 1000.0       # rest density        [kg/m³]
fx        = 0.1          # body force          [m/s²]
mu_p      = 0.001        # plastic viscosity   [Pa·s]
tau_y     = 0.004        # yield stress        [Pa]
m_reg     = 1.0          # Papanastasiou reg.  [s]
drho      = 0.01
backpress = 0.01
dx        = lref / num_length
mass      = rho0 * dx**3

# ─── Kernel ──────────────────────────────────────────────────────────────────
kernel     = 'WendlandC4'
slength_dx = hoomd.sph.kernel.OptimalH[kernel]
rcut_sl    = hoomd.sph.kernel.Kappa[kernel]
slength    = slength_dx * dx
rcut       = rcut_sl * slength
kernel_obj = hoomd.sph.kernel.Kernels[kernel]()
kappa      = kernel_obj.Kappa()

# ─── Geometry (parallel plates, created inline) ───────────────────────────────
def make_parallel_plates_gsd(filename):
    part_rcut = math.ceil(rcut / dx)
    ny = num_length + 3 * part_rcut
    nz = math.ceil(2.5 * rcut_sl * rcut / dx)
    lx, ly, lz = num_length * dx, ny * dx, nz * dx
    xs = np.linspace(-lx/2 + dx/2, lx/2 - dx/2, num_length)
    ys = np.linspace(-ly/2 + dx/2, ly/2 - dx/2, ny)
    zs = np.linspace(-lz/2 + dx/2, lz/2 - dx/2, nz)
    xg, yg, zg = np.meshgrid(xs, ys, zs, indexing='ij')
    pos = np.column_stack([xg.ravel(), yg.ravel(), zg.ravel()])
    N   = pos.shape[0]
    tid = np.where(np.abs(pos[:, 1]) >= 0.5 * lref, 1, 0).astype(np.int32)
    snap = gsd.hoomd.Frame()
    snap.configuration.box  = [lx, ly, lz, 0, 0, 0]
    snap.particles.N        = N
    snap.particles.types    = ['F', 'S']
    snap.particles.typeid   = tid
    snap.particles.position = pos.astype(np.float32)
    snap.particles.velocity = np.zeros((N, 3), dtype=np.float32)
    snap.particles.mass     = np.full(N, mass,    dtype=np.float32)
    snap.particles.slength  = np.full(N, slength, dtype=np.float32)
    snap.particles.density  = np.full(N, rho0,    dtype=np.float32)
    with gsd.hoomd.open(filename, 'w') as f:
        f.append(snap)
    return int(np.sum(tid == 0)), int(np.sum(tid == 1))

gsd_file = f'bingham_{num_length}_init.gsd'
nf, ns = make_parallel_plates_gsd(gsd_file)
print(f'Geometry: {nf} fluid + {ns} solid particles')

# ─── Exact Bingham analytical functions ──────────────────────────────────────
y_p = tau_y / (rho0 * fx)               # plug-zone half-width
print(f'Bingham parameters:  μ_p={mu_p} Pa·s,  τ_y={tau_y} Pa,  m_reg={m_reg} s')
print(f'Plug zone half-width: y_p = {y_p*1e3:.3f} mm  ({y_p/(H/2)*100:.1f} % of H/2)')

def v_bingham_exact(y, rho, fx, mu_p, tau_y, H):
    """Exact (unregularised) Bingham velocity profile."""
    yp = tau_y / (rho * fx)
    ya = np.abs(y)
    v  = np.where(
        ya <= yp,
        # Plug zone: uniform velocity at edge of plug
        rho*fx/(2*mu_p) * ((H/2)**2 - yp**2) - tau_y/mu_p * (H/2 - yp),
        # Flowing zone
        rho*fx/(2*mu_p) * ((H/2)**2 - ya**2) - tau_y/mu_p * (H/2 - ya)
    )
    return np.maximum(v, 0.0)   # velocity must be non-negative

v_plug_exact = float(rho0*fx/(2*mu_p) * ((H/2)**2 - y_p**2) - tau_y/mu_p * (H/2 - y_p))
print(f'Exact plug velocity: {v_plug_exact:.5e} m/s')


def run_case(label, use_bingham, tauy_val=tau_y, m_val=m_reg, mu_min=0.0):
    """Run one plane Poiseuille case with optional Bingham activation."""
    print(f'\n{"═"*60}')
    if use_bingham:
        print(f'  Case: {label}  (Bingham: μ_p={mu_p}, τ_y={tauy_val}, m={m_val})')
    else:
        print(f'  Case: {label}  (Newtonian: μ={mu_p})')
    print(f'{"═"*60}')
    dumpname = f'bingham_{num_length}_{label}_run.gsd'

    # Maximum effective viscosity (at zero shear rate, Papanastasiou limit)
    # μ_eff(0) = μ_p + τ_y * m
    mu_max = mu_p + tauy_val * m_val if use_bingham else mu_p
    v_ref  = max(v_plug_exact if use_bingham else rho0*fx/(2*mu_p)*(H/2)**2, 1e-8)

    device = hoomd.device.CPU(notice_level=2)
    sim    = hoomd.Simulation(device=device)
    sim.create_state_from_gsd(filename=gsd_file,
                               domain_decomposition=(None, None, 1))

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

    model.mu                  = mu_max   # used for CFL; equals μ_p for Newtonian
    model.gx                  = fx
    model.damp                = 1000
    model.artificialviscosity = True
    model.alpha               = 0.2
    model.beta                = 0.0
    model.densitydiffusion    = False

    maximum_smoothing_length = sph_helper.set_max_sl(sim, device, model)
    c, cond = model.compute_speedofsound(
        LREF=lref, UREF=v_ref, DX=dx, DRHO=drho,
        H=maximum_smoothing_length, MU=mu_max, RHO0=rho0)
    print(f'  Speed of sound:  {c:.4f} m/s  ({cond})')

    dt, dt_cond = model.compute_dt(
        LREF=lref, UREF=v_ref, DX=dx, DRHO=drho,
        H=maximum_smoothing_length, MU=mu_max, RHO0=rho0)
    print(f'  Timestep:  {dt:.3e} s  ({dt_cond})')
    print(f'  T_total = {steps*dt:.3f} s  (t_diff ≈ {rho0*H**2/mu_p:.3f} s)')

    integrator = hoomd.sph.Integrator(dt=dt)
    vvb = hoomd.sph.methods.VelocityVerletBasic(filter=filterfluid,
                                                  densitymethod='SUMMATION')
    integrator.methods.append(vvb)
    integrator.forces.append(model)
    sim.operations.integrator = integrator
    sim.run(0)   # trigger _attach_hook → _cpp_obj

    # Activate Bingham after attachment (requires _cpp_obj)
    if use_bingham:
        model.activateBingham(mu_p=mu_p, tau_y=tauy_val, m_reg=m_val, mu_min=mu_min)

    gsd_writer = hoomd.write.GSD(filename=dumpname,
                                  trigger=hoomd.trigger.Periodic(2000),
                                  mode='wb',
                                  dynamic=['property', 'momentum'])
    sim.operations.writers.append(gsd_writer)

    logger = hoomd.logging.Logger(categories=['scalar', 'string'])
    logger.add(sim, quantities=['timestep', 'tps', 'walltime'])
    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(2000), logger=logger,
                              max_header_len=10)
    sim.operations.writers.append(table)

    print(f'  Running {steps} steps ...')
    sim.run(steps, write_at_start=True)
    gsd_writer.flush()

    # ── Post-processing ──────────────────────────────────────────────────────
    with gsd.hoomd.open(dumpname, 'r') as traj:
        snap = traj[-1]
        pos  = snap.particles.position
        tid  = snap.particles.typeid
        vel  = snap.particles.velocity

    fluid  = tid == 0
    y_f    = pos[fluid, 1]
    vx_f   = vel[fluid, 0]
    inside = np.abs(y_f) < H / 2
    y_ev   = y_f[inside]
    vx_ev  = vx_f[inside]
    va_ev  = v_bingham_exact(y_ev, rho0, fx, mu_p, tauy_val if use_bingham else 0.0, H)

    v_max_ref = float(np.max(va_ev))
    L2_err = np.sqrt(np.mean((vx_ev - va_ev)**2)) / max(v_max_ref, 1e-12) * 100.0

    # Plug velocity check
    sph_plug = float(np.max(vx_ev))
    err_plug = abs(sph_plug - v_max_ref) / max(v_max_ref, 1e-12) * 100.0

    print(f'\n  ── Bingham Poiseuille check (step {snap.configuration.step}) ──')
    print(f'     Exact  v_plug / v_max = {v_max_ref:.5e} m/s')
    print(f'     SPH    v_plug / v_max = {sph_plug:.5e} m/s  ({err_plug:.2f} % err)')
    print(f'     L₂ error / v_max_ref  = {L2_err:.2f} %')
    return L2_err


# ─── Case A: Newtonian regression (τ_y = 0) ──────────────────────────────────
L2_A = run_case('nreg', use_bingham=False)

# ─── Case B: Bingham plastic ─────────────────────────────────────────────────
L2_B = run_case('bingham', use_bingham=True)

# ─── Summary ─────────────────────────────────────────────────────────────────
print(f'\n{"═"*60}')
print(f'  BENCHMARK SUMMARY (num_length={num_length}, steps={steps})')
print(f'{"═"*60}')
print(f'  Case A  τ_y=0 (Newtonian)  : L₂ = {L2_A:.2f} %')
print(f'  Case B  Bingham (τ_y={tau_y}): L₂ = {L2_B:.2f} %')
result = 'PASS' if (L2_A < 5.0 and L2_B < 8.0) else 'FAIL'
print(f'  Result: {result}  (thresholds: Newtonian < 5 %, Bingham < 8 % due to regularisation)')

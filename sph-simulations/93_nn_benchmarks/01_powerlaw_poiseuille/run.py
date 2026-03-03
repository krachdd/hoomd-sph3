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

Power-Law flow — Non-Newtonian benchmark.

BENCHMARK DESCRIPTION
---------------------
Two sub-cases with the Power-Law constitutive model  μ_eff = K · |γ̇|^(n−1):

  Case A — Plane Poiseuille, n=1, K=0.001 Pa·s  (Newtonian regression)
    Analytical: v(y) = (ρ f_x / K)^(1/n) · n/(n+1) · ((H/2)^(1+1/n) − |y|^(1+1/n))
    For n=1, K=μ this reduces to the standard Newtonian result.

  Case B — Plane Couette, n=2, K=0.001 Pa·s²  (shear-thickening, uniform γ̇)
    In Couette flow the shear rate is constant: γ̇ = U_wall / H.
    μ_eff = K · γ̇^(n−1) = K · (U_wall/H) = 0.01 Pa·s throughout.
    Analytical: v(y) = U_wall · (y + H/2) / H  (linear, same shape as Newtonian).
    The key check is that the profile is linear and L₂ < 5 %.

Usage:
    python3 run.py [num_length [steps]]
    Defaults: num_length=20, steps=10001
"""

import sys, os, math
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../../'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules'))

import hoomd
from hoomd import sph
import numpy as np
import gsd.hoomd
import sph_helper

# ─── CLI args ────────────────────────────────────────────────────────────────
num_length = int(sys.argv[1]) if len(sys.argv) > 1 else 20
steps      = int(sys.argv[2]) if len(sys.argv) > 2 else 10001

# ─── Fixed physical parameters ────────────────────────────────────────────────
lref      = 0.001        # channel gap         [m]
H         = lref
rho0      = 1000.0       # rest density        [kg/m³]
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

# ─── Geometry helpers ────────────────────────────────────────────────────────
def _base_lattice(nx, dx, slength, rcut, rcut_sl, lref):
    """Return (positions, ny, nz) for a plane-channel geometry."""
    part_rcut = math.ceil(rcut / dx)
    ny = nx + 3 * part_rcut
    nz = math.ceil(2.5 * rcut_sl * rcut / dx)
    lx, ly, lz = nx*dx, ny*dx, nz*dx
    xs = np.linspace(-lx/2 + dx/2, lx/2 - dx/2, nx)
    ys = np.linspace(-ly/2 + dx/2, ly/2 - dx/2, ny)
    zs = np.linspace(-lz/2 + dx/2, lz/2 - dx/2, nz)
    xg, yg, zg = np.meshgrid(xs, ys, zs, indexing='ij')
    pos = np.column_stack([xg.ravel(), yg.ravel(), zg.ravel()])
    return pos, lx, ly, lz

def make_poiseuille_gsd(filename):
    """Parallel plates, no-slip, body-force driven."""
    pos, lx, ly, lz = _base_lattice(num_length, dx, slength, rcut, rcut_sl, lref)
    N   = pos.shape[0]
    tid = np.where(np.abs(pos[:, 1]) >= 0.5*lref, 1, 0).astype(np.int32)
    vel = np.zeros((N, 3), dtype=np.float32)
    snap = gsd.hoomd.Frame()
    snap.configuration.box  = [lx, ly, lz, 0, 0, 0]
    snap.particles.N, snap.particles.types, snap.particles.typeid = N, ['F','S'], tid
    snap.particles.position = pos.astype(np.float32)
    snap.particles.velocity = vel
    snap.particles.mass     = np.full(N, mass,    dtype=np.float32)
    snap.particles.slength  = np.full(N, slength, dtype=np.float32)
    snap.particles.density  = np.full(N, rho0,    dtype=np.float32)
    with gsd.hoomd.open(filename, 'w') as f: f.append(snap)
    return int(np.sum(tid==0)), int(np.sum(tid==1))

def make_couette_gsd(filename, U_wall):
    """Parallel plates; top wall moves at U_wall, bottom stationary."""
    pos, lx, ly, lz = _base_lattice(num_length, dx, slength, rcut, rcut_sl, lref)
    N   = pos.shape[0]
    y   = pos[:, 1]
    top = y >  0.5 * lref
    bot = y < -0.5 * lref
    tid = np.where(top | bot, 1, 0).astype(np.int32)
    vel = np.zeros((N, 3), dtype=np.float32)
    vel[top, 0] = U_wall
    snap = gsd.hoomd.Frame()
    snap.configuration.box  = [lx, ly, lz, 0, 0, 0]
    snap.particles.N, snap.particles.types, snap.particles.typeid = N, ['F','S'], tid
    snap.particles.position = pos.astype(np.float32)
    snap.particles.velocity = vel
    snap.particles.mass     = np.full(N, mass,    dtype=np.float32)
    snap.particles.slength  = np.full(N, slength, dtype=np.float32)
    snap.particles.density  = np.full(N, rho0,    dtype=np.float32)
    with gsd.hoomd.open(filename, 'w') as f: f.append(snap)
    return int(np.sum(tid==0)), int(np.sum(tid==1))

# ─── Analytical solutions ─────────────────────────────────────────────────────
def v_powerlaw_poiseuille(y, K, n, rho, fx, H):
    return (rho*fx/K)**(1.0/n) * (n/(n+1)) * ((H/2)**(1+1/n) - np.abs(y)**(1+1/n))

def v_couette_linear(y, U_wall, H):
    return U_wall * (y + H/2) / H

# ─── Case runner ─────────────────────────────────────────────────────────────
def run_powerlaw(label, gsd_file, K, n, mu_rep, v_ref, fx_val, U_wall,
                 va_fn, mu_min=0.0):
    """Run one Power-Law case and return L₂ error."""
    print(f'\n{"═"*60}')
    print(f'  Case: {label}  K={K}, n={n}')
    print(f'{"═"*60}')
    dumpname = f'{label}_run.gsd'

    device = hoomd.device.CPU(notice_level=2)
    sim    = hoomd.Simulation(device=device)
    sim.create_state_from_gsd(filename=gsd_file,
                               domain_decomposition=(None, None, 1))

    nlist = hoomd.nsearch.nlist.Cell(buffer=rcut*0.05, rebuild_check_delay=1,
                                     kappa=kappa)
    eos = hoomd.sph.eos.Tait(); eos.set_params(rho0, backpress)
    ff  = hoomd.filter.Type(['F'])
    fs  = hoomd.filter.Type(['S'])

    model = hoomd.sph.sphmodel.SinglePhaseFlow(
        kernel=kernel_obj, eos=eos, nlist=nlist,
        fluidgroup_filter=ff, solidgroup_filter=fs,
        densitymethod='SUMMATION')

    model.mu                  = mu_rep
    model.gx                  = fx_val
    model.damp                = 1000
    model.artificialviscosity = True
    model.alpha               = 0.2
    model.beta                = 0.0
    model.densitydiffusion    = False

    max_sl = sph_helper.set_max_sl(sim, device, model)
    c, cc  = model.compute_speedofsound(LREF=lref, UREF=v_ref, DX=dx,
                                         DRHO=drho, H=max_sl, MU=mu_rep, RHO0=rho0)
    print(f'  c₀ = {c:.4f} m/s  ({cc}),  v_ref = {v_ref:.4e} m/s')
    dt, dc = model.compute_dt(LREF=lref, UREF=v_ref, DX=dx,
                               DRHO=drho, H=max_sl, MU=mu_rep, RHO0=rho0)
    print(f'  dt = {dt:.3e} s  ({dc}),  T_total = {steps*dt:.3f} s')

    integrator = hoomd.sph.Integrator(dt=dt)
    vvb = hoomd.sph.methods.VelocityVerletBasic(filter=ff, densitymethod='SUMMATION')
    integrator.methods.append(vvb)
    integrator.forces.append(model)
    sim.operations.integrator = integrator
    sim.run(0)                        # trigger _attach_hook → _cpp_obj

    model.activatePowerLaw(K=K, n=n, mu_min=mu_min)

    gw = hoomd.write.GSD(filename=dumpname,
                          trigger=hoomd.trigger.Periodic(max(steps//10, 1)),
                          mode='wb', dynamic=['property', 'momentum'])
    sim.operations.writers.append(gw)
    log = hoomd.logging.Logger(categories=['scalar', 'string'])
    log.add(sim, quantities=['timestep', 'tps', 'walltime'])
    sim.operations.writers.append(
        hoomd.write.Table(trigger=hoomd.trigger.Periodic(max(steps//10, 1)),
                          logger=log, max_header_len=10))

    print(f'  Running {steps} steps ...')
    sim.run(steps, write_at_start=True)
    gw.flush()

    with gsd.hoomd.open(dumpname, 'r') as traj:
        snap = traj[-1]
        pos, tid, vel = snap.particles.position, snap.particles.typeid, snap.particles.velocity

    fluid  = (tid == 0)
    y_f    = pos[fluid, 1]
    vx_f   = vel[fluid, 0]
    inside = np.abs(y_f) < H/2
    y_ev   = y_f[inside]; vx_ev = vx_f[inside]
    va_ev  = va_fn(y_ev)
    v_ref2 = max(float(np.max(np.abs(va_ev))), 1e-12)
    L2     = np.sqrt(np.mean((vx_ev - va_ev)**2)) / v_ref2 * 100.0

    print(f'\n  ── check (step {snap.configuration.step}) ──')
    print(f'     analytical v_max = {v_ref2:.5e} m/s')
    print(f'     SPH        v_max = {float(np.max(vx_ev)):.5e} m/s')
    print(f'     L₂ / v_max       = {L2:.2f} %')
    return L2


# ═══════════════════════════════════════════════════════════════════════
# Case A: n=1 Poiseuille  (Newtonian regression)
# ═══════════════════════════════════════════════════════════════════════
fx_A   = 0.1
K_A    = 0.001    # Pa·s  (same as mu for n=1)
n_A    = 1.0
# For n=1 the wall effective viscosity = K; v_max = f_x/(2*nu)*(H/2)^2
v_max_A = (rho0*fx_A/K_A)**(1/n_A) * (n_A/(n_A+1)) * (H/2)**(1+1/n_A)
gsd_A = f'poiseuille_{num_length}_init.gsd'
nf, ns = make_poiseuille_gsd(gsd_A)
print(f'[A] Poiseuille geometry: {nf} fluid + {ns} solid')

L2_A = run_powerlaw(
    label   = 'caseA_poiseuille_n1',
    gsd_file= gsd_A,
    K=K_A, n=n_A,
    mu_rep  = K_A,
    v_ref   = v_max_A,
    fx_val  = fx_A,
    U_wall  = 0.0,
    va_fn   = lambda y: v_powerlaw_poiseuille(y, K_A, n_A, rho0, fx_A, H),
)

# ═══════════════════════════════════════════════════════════════════════
# Case B: n=2 Couette  (shear-thickening, uniform shear rate)
# ═══════════════════════════════════════════════════════════════════════
U_wall_B = 0.01   # [m/s]
K_B      = 0.001  # Pa·s²  (n=2)
n_B      = 2.0
# Uniform shear rate: γ̇ = U_wall/H  →  μ_eff = K*(U_wall/H)^(n-1) = 0.01 Pa·s
gamma_B  = U_wall_B / H
mu_eff_B = K_B * gamma_B**(n_B - 1.0)
gsd_B = f'couette_{num_length}_init.gsd'
nf, ns = make_couette_gsd(gsd_B, U_wall_B)
print(f'\n[B] Couette geometry:   {nf} fluid + {ns} solid')
print(f'    γ̇ = {gamma_B:.1f} s⁻¹,  μ_eff = {mu_eff_B:.4f} Pa·s')

L2_B = run_powerlaw(
    label   = 'caseB_couette_n2',
    gsd_file= gsd_B,
    K=K_B, n=n_B,
    mu_rep  = mu_eff_B,
    v_ref   = U_wall_B,
    fx_val  = 0.0,
    U_wall  = U_wall_B,
    va_fn   = lambda y: v_couette_linear(y, U_wall_B, H),
)

# ─── Summary ─────────────────────────────────────────────────────────────────
print(f'\n{"═"*60}')
print(f'  BENCHMARK SUMMARY (num_length={num_length}, steps={steps})')
print(f'{"═"*60}')
print(f'  Case A  Poiseuille n=1 (Newtonian regression) : L₂ = {L2_A:.2f} %')
print(f'  Case B  Couette    n=2 (shear-thickening)     : L₂ = {L2_B:.2f} %')
result = 'PASS' if (L2_A < 5.0 and L2_B < 5.0) else 'FAIL'
print(f'  Result: {result}  (threshold: L₂ < 5 %)')

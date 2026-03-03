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

Two-phase Couette flow with Power-Law lower fluid — Non-Newtonian benchmark.

BENCHMARK DESCRIPTION
---------------------
Two immiscible fluid layers sheared between a stationary bottom wall
(y = −H/2) and a top wall moving at U_wall (y = +H/2).  No gravity.

  Layer W (lower, y ∈ [−H/2, 0]): Power Law  K=0.004 Pa·s^n, n varies
  Layer N (upper, y ∈ [0,   H/2]): Newtonian  μ₂=0.01 Pa·s

For steady-state Couette flow the shear stress is uniform across both layers.
With piecewise linear velocity profile:

    v_W(y) = γ̇_W · (y + H/2)    for y ∈ [−H/2, 0]
    v_N(y) = v_i + γ̇_N · y      for y ∈ [0, H/2]

where τ is found by solving (numerically for n ≠ 1):

    (τ/K)^(1/n) · H/2  +  τ/μ₂ · H/2  =  U_wall

and  v_i = (τ/K)^(1/n) · H/2.

Two sub-cases are run:
  Case A  n=1 (Newtonian regression):
    v_i = U_wall · μ₂/(K+μ₂) = 0.01 × 0.01/0.014 = 7.14 mm/s
    Verify that activatePowerLaw1(K, n=1) gives the same result as Newtonian.

  Case B  n=2 (shear-thickening W layer):
    μ_eff_W = K · γ̇_W,  which is HIGHER than K → slower W-layer flow
    v_i < Newtonian v_i  (shear-thickening resists the top wall)
    Convergence is fast because μ_eff_W >> K → short diffusion time.

Parameters:  K=0.004, μ₂=0.01, H=0.001, U_wall=0.01  (all SI)

Usage:
    python3 run.py [num_length [steps]]
    Defaults: num_length=20, steps=10001
"""

import sys, os
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../../'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules'))

import hoomd
from hoomd import sph
import numpy as np
from scipy.optimize import brentq
import gsd.hoomd
import sph_helper

# ─── CLI args ────────────────────────────────────────────────────────────────
num_length = int(sys.argv[1]) if len(sys.argv) > 1 else 20
steps      = int(sys.argv[2]) if len(sys.argv) > 2 else 10001

# ─── Physical parameters ─────────────────────────────────────────────────────
lref      = 0.001
H         = lref
U_wall    = 0.01          # top wall velocity         [m/s]
rho01     = 1000.0        # phase W rest density      [kg/m³]
rho02     = 1000.0        # phase N rest density      [kg/m³]
K_W       = 0.004         # Power Law consistency     [Pa·s^n]
mu_N      = 0.01          # Newtonian viscosity       [Pa·s]
sigma     = 0.0
backpress = 0.01
drho      = 0.01
dx        = lref / num_length
n_solid   = 3
mass      = rho01 * dx**3

# ─── Kernel ──────────────────────────────────────────────────────────────────
kernel     = 'WendlandC4'
slength    = hoomd.sph.kernel.OptimalH[kernel] * dx
rcut       = hoomd.sph.kernel.Kappa[kernel] * slength
kernel_obj = hoomd.sph.kernel.Kernels[kernel]()
kappa      = kernel_obj.Kappa()

# ─── Geometry ─────────────────────────────────────────────────────────────────
def make_couette_gsd(filename):
    nx = num_length; ny = num_length + 2*n_solid; nz = num_length
    lx, ly, lz = nx*dx, ny*dx, nz*dx
    y_bot, y_top = -(H/2 + n_solid*dx), H/2 + n_solid*dx
    x_arr = np.linspace(-lx/2+dx/2, lx/2-dx/2, nx)
    y_arr = np.linspace(y_bot+dx/2,  y_top-dx/2, ny)
    z_arr = np.linspace(-lz/2+dx/2, lz/2-dx/2, nz)
    xg, yg, zg = np.meshgrid(x_arr, y_arr, z_arr, indexing='ij')
    pos = np.stack((xg.ravel(), yg.ravel(), zg.ravel()), axis=1).astype(np.float32)
    N   = len(pos); y = pos[:, 1]
    tid = np.zeros(N, dtype=np.int32)
    vel = np.zeros((N, 3), dtype=np.float32)
    tid[(y >= 0.0) & (y < H/2)] = 1    # 'N' upper fluid
    sol_bot = y < -H/2;  sol_top = y >= H/2
    tid[sol_bot | sol_top] = 2          # 'S' solid
    vel[sol_top, 0] = U_wall
    snap = gsd.hoomd.Frame()
    snap.configuration.box  = [lx, ly, lz, 0, 0, 0]
    snap.particles.N        = N
    snap.particles.types    = ['W', 'N', 'S']
    snap.particles.typeid   = tid
    snap.particles.position = pos
    snap.particles.velocity = vel
    snap.particles.mass     = np.full(N, mass, dtype=np.float32)
    snap.particles.slength  = np.full(N, slength, dtype=np.float32)
    snap.particles.density  = np.full(N, rho01, dtype=np.float32)
    with gsd.hoomd.open(filename, 'w') as f: f.append(snap)
    return int(np.sum(tid==0)), int(np.sum(tid==1)), int(np.sum(tid==2))

gsd_file = f'couette_{num_length}_init.gsd'
nW, nN, nS = make_couette_gsd(gsd_file)
print(f'Geometry: {nW} W + {nN} N + {nS} S particles')

# ─── Analytical solution ──────────────────────────────────────────────────────
def solve_couette(K, n, mu2, U_wall, H):
    H1 = H/2
    if abs(n-1.0) < 1e-10:
        tau = U_wall / (H1*(1.0/K + 1.0/mu2))
    else:
        def res(tau):
            return (tau/K)**(1.0/n)*H1 + tau/mu2*H1 - U_wall
        tau = brentq(res, 1e-20, max(K,mu2)*U_wall/H1*100, xtol=1e-15)
    gdW = (tau/K)**(1.0/n)
    gdN = tau/mu2
    v_i = gdW * H1
    return tau, v_i, gdW, gdN

def v_piecewise(y, gdW, gdN, v_i, H):
    return np.where(y < 0.0, gdW*(y + H/2), v_i + gdN*y)


# ─── Runner ──────────────────────────────────────────────────────────────────
def run_case(label, n_W):
    print(f'\n{"═"*60}')
    print(f'  Case: {label}  K_W={K_W}, n_W={n_W}, μ_N={mu_N}')
    print(f'{"═"*60}')
    dumpname = f'{label}_run.gsd'

    tau, v_i, gdW, gdN = solve_couette(K_W, n_W, mu_N, U_wall, H)
    mu_eff_W = K_W * gdW**(n_W-1.0)
    print(f'  Analytical τ   = {tau:.5e} Pa')
    print(f'  Analytical v_i = {v_i:.5e} m/s')
    print(f'  μ_eff_W        = {mu_eff_W:.5e} Pa·s')

    mu_rep1 = max(K_W, mu_eff_W)   # conservative upper bound for CFL

    device = hoomd.device.CPU(notice_level=2)
    sim    = hoomd.Simulation(device=device)
    sim.create_state_from_gsd(filename=gsd_file)

    nlist = hoomd.nsearch.nlist.Cell(buffer=rcut*0.05, rebuild_check_delay=1,
                                     kappa=kappa)
    eos1 = hoomd.sph.eos.Tait(); eos1.set_params(rho01, backpress)
    eos2 = hoomd.sph.eos.Tait(); eos2.set_params(rho02, backpress)

    fW = hoomd.filter.Type(['W'])
    fN = hoomd.filter.Type(['N'])
    fS = hoomd.filter.Type(['S'])

    model = hoomd.sph.sphmodel.TwoPhaseFlow(
        kernel=kernel_obj, eos1=eos1, eos2=eos2, nlist=nlist,
        fluidgroup1_filter=fW, fluidgroup2_filter=fN, solidgroup_filter=fS,
        densitymethod='SUMMATION', colorgradientmethod='DENSITYRATIO')

    model.mu1 = mu_rep1;  model.mu2 = mu_N
    model.sigma12 = sigma; model.omega = 90
    model.gx = 0.0;        model.damp = 1000
    model.artificialviscosity = True
    model.alpha = 0.2;     model.beta = 0.0
    model.densitydiffusion = False

    max_sl = sph_helper.set_max_sl(sim, device, model)
    c1, cc1, c2, cc2 = model.compute_speedofsound(
        LREF=lref, UREF=U_wall, DX=dx, DRHO=drho, H=max_sl,
        MU1=mu_rep1, MU2=mu_N, RHO01=rho01, RHO02=rho02, SIGMA12=sigma)
    print(f'  c₁={c1:.4f} m/s ({cc1}),  c₂={c2:.4f} m/s ({cc2})')

    sph_helper.update_min_c0_tpf(device, model, c1, c2,
                                  mode='plain', lref=lref, uref=U_wall, cfactor=10.0)

    dt, dc = model.compute_dt(
        LREF=lref, UREF=U_wall, DX=dx, DRHO=drho, H=max_sl,
        MU1=mu_rep1, MU2=mu_N, RHO01=rho01, RHO02=rho02, SIGMA12=sigma)
    t_D_N = rho02*H**2/mu_N
    print(f'  dt={dt:.3e} s ({dc}),  T_total={steps*dt:.4f} s,  t_D_N={t_D_N:.4f} s')

    integrator = hoomd.sph.Integrator(dt=dt)
    vvbW = hoomd.sph.methods.VelocityVerletBasic(filter=fW, densitymethod='SUMMATION')
    vvbN = hoomd.sph.methods.VelocityVerletBasic(filter=fN, densitymethod='SUMMATION')
    integrator.methods.append(vvbW); integrator.methods.append(vvbN)
    integrator.forces.append(model)
    sim.operations.integrator = integrator
    sim.run(0)   # trigger _attach_hook → _cpp_obj

    model.activatePowerLaw1(K=K_W, n=n_W, mu_min=0.0)

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

    fluid  = (tid == 0) | (tid == 1)
    y_f    = pos[fluid, 1];  vx_f = vel[fluid, 0]
    in_gap = np.abs(y_f) < H/2
    y_ev   = y_f[in_gap];   vx_ev = vx_f[in_gap]
    va_ev  = v_piecewise(y_ev, gdW, gdN, v_i, H)

    near = np.abs(y_ev) < 2*dx
    v_i_sph = float(np.mean(vx_ev[near])) if near.any() else float('nan')

    L2 = np.sqrt(np.mean((vx_ev - va_ev)**2)) / U_wall * 100.0
    err_vi = abs(v_i_sph - v_i) / U_wall * 100.0

    print(f'\n  ── Two-phase Couette check (step {snap.configuration.step}) ──')
    print(f'     Analytical v_i  = {v_i:.5e} m/s')
    print(f'     SPH        v_i ≈ {v_i_sph:.5e} m/s  ({err_vi:.2f} % err of U_wall)')
    print(f'     L₂ / U_wall     = {L2:.2f} %')
    return L2


# ─── Case A: Newtonian regression (n=1) ───────────────────────────────────────
L2_A = run_case('caseA_n1', n_W=1.0)

# ─── Case B: Shear-thickening (n=2) ───────────────────────────────────────────
L2_B = run_case('caseB_n2', n_W=2.0)

# ─── Summary ─────────────────────────────────────────────────────────────────
print(f'\n{"═"*60}')
print(f'  BENCHMARK SUMMARY (num_length={num_length}, steps={steps})')
print(f'{"═"*60}')
print(f'  Case A  n=1 (Newtonian regression) : L₂ = {L2_A:.2f} %')
print(f'  Case B  n=2 (shear-thickening W)   : L₂ = {L2_B:.2f} %')
result = 'PASS' if (L2_A < 5.0 and L2_B < 5.0) else 'FAIL'
print(f'  Result: {result}  (threshold: L₂ < 5 %)')

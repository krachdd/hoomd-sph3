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

Bingham plastic slump test (concrete-like) — Non-Newtonian benchmark.

BENCHMARK DESCRIPTION
---------------------
A rectangular column of Bingham plastic (H0 × 2*L0) collapses under gravity
onto a solid floor.  At equilibrium the yield stress balances the gravitational
surface stress everywhere, producing a parabolic (Herschel) free-surface profile.

Analytical equilibrium profile (2-D Bingham dam-break, symmetric):

    h(x) = sqrt(h_c² − 2·τ_y·|x| / (ρ·g))   for |x| ≤ x_f

    h_c = (3·τ_y·H0·L0 / (ρ·g))^(1/3)        (centre height, volume conservation)
    x_f = h_c²·ρ·g / (2·τ_y)                  (front position, h(x_f) = 0)
    S   = H0 − h_c                              (slump height)

Two sub-cases:

  Case A — Newtonian regression (τ_y = 0)
    Model reduces to Newtonian with μ = μ_p.  Fluid spreads freely.
    Front position X*(T*) compared to Martin-Moyce (1952) shallow-water theory
    (inviscid reference, Re ≈ 20 so deviation expected):
        X*(T*) = 1 + 2√2·T*,   T* = t·√(g/L0)
    Pass: X*(T*=2) > 1.3 (front advancing — not arrested like Bingham case).

  Case B — Bingham slump (τ_y = 50 Pa)
    Fluid arrests at the parabolic Herschel profile.
    Pass: h_c error < 25 %, x_f error < 30 %, profile L₂/H0 < 20 %.
    (Papanastasiou regularisation never fully arrests; at nl=20 h_c spans ~6.5
    particle diameters, limiting discrete accuracy to ~20-25 %.)

Model: SinglePhaseFlowFS + KickDriftKickTV (same as dam-break benchmark).
activateBingham is available on SinglePhaseFlowFS; sigma=0 disables surface tension.

Parameters:
  H0    = 0.10 m    (initial column height; also lref)
  L0    = 0.05 m    (initial half-width; aspect ratio H0/(2L0) = 1)
  ρ     = 2200 kg/m³
  μ_p   = 10 Pa·s   (plastic viscosity)
  τ_y   = 50 Pa     (yield stress)
  m_reg = 1 s       (Papanastasiou regularisation)
  g     = 9.81 m/s²
  μ_max = μ_p + τ_y·m_reg = 60 Pa·s  (Papanastasiou limit; used for CFL in Case B)

Usage:
    python3 run.py [num_length [steps_A [steps_B]]]
    Defaults: num_length=20, steps_A=5000, steps_B=20000
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
steps_A    = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
steps_B    = int(sys.argv[3]) if len(sys.argv) > 3 else 20000

# ─── Physical parameters ─────────────────────────────────────────────────────
lref      = 0.10          # reference length = H0          [m]
H0        = lref          # initial column height           [m]
L0        = 0.05          # initial half-width              [m]
rho0      = 2200.0        # rest density (concrete-like)    [kg/m³]
mu_p      = 10.0          # plastic viscosity               [Pa·s]
tau_y     = 50.0          # yield stress                    [Pa]
m_reg     = 1.0           # Papanastasiou regularisation    [s]
g         = 9.81          # gravitational acceleration      [m/s²]
mu_max    = mu_p + tau_y * m_reg   # = 60 Pa·s (Papanastasiou limit)
drho      = 0.01
backpress = 0.01
x_box     = 0.70          # total horizontal domain width   [m]

dx        = lref / num_length          # particle spacing    [m]
mass      = rho0 * dx**3              # particle mass        [kg]

# ─── Kernel ──────────────────────────────────────────────────────────────────
kernel     = 'WendlandC4'
slength_dx = hoomd.sph.kernel.OptimalH[kernel]   # 1.7
rcut_sl    = hoomd.sph.kernel.Kappa[kernel]       # 2.0
slength    = slength_dx * dx
rcut       = rcut_sl * slength
kernel_obj = hoomd.sph.kernel.Kernels[kernel]()
kappa      = kernel_obj.Kappa()

# ─── Analytical predictions (2-D Bingham equilibrium) ────────────────────────
h_c_exact = (3.0 * tau_y * H0 * L0 / (rho0 * g)) ** (1.0 / 3.0)
x_f_exact = h_c_exact**2 * rho0 * g / (2.0 * tau_y)
slump_S   = H0 - h_c_exact

print(f'Bingham slump test parameters:')
print(f'  H0={H0*1e2:.1f} cm, L0={L0*1e2:.1f} cm, ρ={rho0:.0f} kg/m³')
print(f'  μ_p={mu_p} Pa·s, τ_y={tau_y} Pa, m_reg={m_reg} s')
print(f'  μ_max (Papanastasiou limit) = {mu_max} Pa·s')
print(f'\nAnalytical equilibrium profile:')
print(f'  h_c  = {h_c_exact*1e3:.3f} mm  ({h_c_exact/H0*100:.1f} % of H0)')
print(f'  x_f  = {x_f_exact*1e3:.1f} mm  (half-spread from x=0)')
print(f'  S    = {slump_S*1e3:.1f} mm  (slump = {slump_S/H0*100:.1f} % of H0)')


def h_analytical(x):
    """Herschel parabolic free-surface profile at Bingham equilibrium."""
    arg = h_c_exact**2 - 2.0 * tau_y * np.abs(x) / (rho0 * g)
    return np.where(arg > 0.0, np.sqrt(arg), 0.0)


# ─── Geometry ────────────────────────────────────────────────────────────────
GSD_FILE = f'slump_{num_length}_init.gsd'


def make_slump_gsd(filename):
    """
    Create the initial GSD snapshot: a fluid column resting on a solid floor.

    Layout (physical coordinates, y=0 at floor):
        Fluid: x ∈ [-L0, L0],  y ∈ [0, H0]   (num_length × num_length in x-y)
        Solid: x ∈ [-x_box/2, x_box/2],  y < 0  (n_solid layers, full width)
        z: quasi-2D (nz layers, periodic)
    """
    n_solid = math.ceil(rcut / dx) + 1          # solid layers below y=0
    nz      = math.ceil(2.5 * rcut_sl * rcut / dx)   # quasi-2D depth
    nx_box  = int(round(x_box / dx))            # 140 at num_length=20
    ly      = 0.40                              # vertical box extent [m]
    lx      = nx_box * dx
    lz      = nz * dx

    # Particle centre arrays
    xs_all   = np.linspace(-lx / 2 + dx / 2, lx / 2 - dx / 2, nx_box)
    nx_fluid = int(round(2.0 * L0 / dx))        # particles across 2*L0
    ny_fluid = num_length                        # particles over H0
    xs_fluid = np.linspace(-L0 + dx / 2, L0 - dx / 2, nx_fluid)
    ys_fluid = np.linspace(dx / 2, H0 - dx / 2, ny_fluid)
    ys_solid = np.array([-(i + 0.5) * dx for i in range(n_solid)])
    zs       = np.linspace(-lz / 2 + dx / 2, lz / 2 - dx / 2, nz)

    # Fluid block
    xg_f, yg_f, zg_f = np.meshgrid(xs_fluid, ys_fluid, zs, indexing='ij')
    pos_f = np.column_stack([xg_f.ravel(), yg_f.ravel(), zg_f.ravel()])

    # Solid floor (full x-width)
    xg_s, yg_s, zg_s = np.meshgrid(xs_all, ys_solid, zs, indexing='ij')
    pos_s = np.column_stack([xg_s.ravel(), yg_s.ravel(), zg_s.ravel()])

    pos_all = np.vstack([pos_f, pos_s]).astype(np.float32)
    N       = pos_all.shape[0]
    tid     = np.zeros(N, dtype=np.int32)
    tid[len(pos_f):] = 1   # solid = type 1

    snap = gsd.hoomd.Frame()
    snap.configuration.box  = [lx, ly, lz, 0, 0, 0]
    snap.particles.N        = N
    snap.particles.types    = ['F', 'S']
    snap.particles.typeid   = tid
    snap.particles.position = pos_all
    snap.particles.velocity = np.zeros((N, 3), dtype=np.float32)
    snap.particles.mass     = np.full(N, mass,    dtype=np.float32)
    snap.particles.slength  = np.full(N, slength, dtype=np.float32)
    snap.particles.density  = np.full(N, rho0,    dtype=np.float32)

    with gsd.hoomd.open(filename, 'w') as f:
        f.append(snap)

    return int(np.sum(tid == 0)), int(np.sum(tid == 1))


nf, ns = make_slump_gsd(GSD_FILE)
print(f'\nGeometry (num_length={num_length}):')
print(f'  dx = {dx*1e3:.2f} mm')
print(f'  {nf} fluid + {ns} solid particles')


# ─── Case runner ─────────────────────────────────────────────────────────────
def run_case(label, tau_y_val, steps_val):
    """Run one slump-test case (Newtonian or Bingham) and return metrics."""
    print(f'\n{"═"*60}')
    if tau_y_val > 0.0:
        print(f'  Case: {label}  (Bingham: μ_p={mu_p}, τ_y={tau_y_val}, m={m_reg})')
    else:
        print(f'  Case: {label}  (Newtonian: μ={mu_p})')
    print(f'{"═"*60}')

    dumpname = f'slump_{num_length}_{label}_run.gsd'
    mu_case  = mu_max if tau_y_val > 0.0 else mu_p
    U_char   = math.sqrt(g * H0)   # gravity-wave reference velocity [m/s]

    device = hoomd.device.CPU(notice_level=2)
    sim    = hoomd.Simulation(device=device)
    sim.create_state_from_gsd(filename=GSD_FILE, domain_decomposition=(None, None, 1))

    nlist = hoomd.nsearch.nlist.Cell(buffer=rcut * 0.05, rebuild_check_delay=1,
                                     kappa=kappa)
    eos = hoomd.sph.eos.Tait()
    eos.set_params(rho0, backpress)

    filterfluid = hoomd.filter.Type(['F'])
    filtersolid = hoomd.filter.Type(['S'])

    model = hoomd.sph.sphmodel.SinglePhaseFlowFS(
        kernel=kernel_obj, eos=eos, nlist=nlist,
        fluidgroup_filter=filterfluid, solidgroup_filter=filtersolid,
        densitymethod='SUMMATION',
        sigma=0.0, fs_threshold=0.99, contact_angle=math.pi / 2)

    model.mu                  = mu_case
    model.gy                  = -g       # gravity in -y direction
    model.damp                = 0        # no artificial damping (dynamics)
    model.artificialviscosity = True
    model.alpha               = 0.1
    model.beta                = 0.0
    model.densitydiffusion    = False

    max_sl = sph_helper.set_max_sl(sim, device, model)

    c, cond = model.compute_speedofsound(
        LREF=lref, UREF=U_char, DX=dx, DRHO=drho,
        H=max_sl, MU=mu_case, RHO0=rho0)
    print(f'  c₀ = {c:.4f} m/s  ({cond})')

    sph_helper.update_min_c0(device, model, c,
                              mode='both', lref=lref, uref=U_char,
                              bforce=g, cfactor=10.0)

    dt, dt_cond = model.compute_dt(
        LREF=lref, UREF=U_char, DX=dx, DRHO=drho,
        H=max_sl, MU=mu_case, RHO0=rho0)

    t_I = math.sqrt(H0 / g)   # inertial timescale
    print(f'  dt = {dt:.3e} s  ({dt_cond})')
    print(f'  T_total = {steps_val * dt:.3f} s  ({steps_val * dt / t_I:.1f} t_I,  t_I = {t_I:.3f} s)')

    # KickDriftKickTV integrator (required for SinglePhaseFlowFS)
    integrator = hoomd.sph.Integrator(dt=dt)
    kdktv = hoomd.sph.methods.KickDriftKickTV(filter=filterfluid,
                                               densitymethod='SUMMATION')
    integrator.methods.append(kdktv)
    integrator.forces.append(model)
    sim.operations.integrator = integrator

    sim.run(0)   # trigger _attach_hook → creates _cpp_obj

    # Activate Bingham (tau_y=0 → Newtonian, tau_y>0 → Bingham plastic)
    model.activateBingham(mu_p=mu_p, tau_y=tau_y_val, m_reg=m_reg, mu_min=0.0)

    # Write period: many frames for Case A time-series, fewer for Case B profile
    write_period = max(steps_val // 200, 1) if tau_y_val == 0.0 else max(steps_val // 50, 1)
    gsd_writer = hoomd.write.GSD(filename=dumpname,
                                  trigger=hoomd.trigger.Periodic(write_period),
                                  mode='wb',
                                  dynamic=['property', 'momentum'])
    sim.operations.writers.append(gsd_writer)

    logger = hoomd.logging.Logger(categories=['scalar', 'string'])
    logger.add(sim, quantities=['timestep', 'tps', 'walltime'])
    table = hoomd.write.Table(
        trigger=hoomd.trigger.Periodic(max(steps_val // 20, 1)),
        logger=logger, max_header_len=10)
    sim.operations.writers.append(table)

    print(f'  Running {steps_val} steps ...')
    sim.run(steps_val, write_at_start=True)
    gsd_writer.flush()

    if tau_y_val == 0.0:
        return _check_martin_moyce(dumpname, dt)
    else:
        return _check_bingham_slump(dumpname)


# ─── Post-processing: Case A — Martin-Moyce front position ───────────────────
def _check_martin_moyce(dumpname, dt):
    """Compare Newtonian front position to Martin-Moyce shallow-water theory."""
    a = L0   # length scale (initial half-width)
    x_front_arr = []
    t_arr = []

    with gsd.hoomd.open(dumpname, 'r') as traj:
        snap0   = traj[0]
        fluid0  = snap0.particles.typeid == 0
        x_front0 = float(np.max(snap0.particles.position[fluid0, 0]))

        for snap in traj:
            fluid   = snap.particles.typeid == 0
            x_front = float(np.max(snap.particles.position[fluid, 0]))
            t_phys  = snap.configuration.step * dt
            x_front_arr.append(x_front)
            t_arr.append(t_phys)

    x_front_arr = np.array(x_front_arr)
    t_arr       = np.array(t_arr)

    # Dimensionless variables  (X*(0) = 1 by construction)
    X_star_sph = (x_front_arr - x_front0 + a) / a
    T_star     = t_arr * math.sqrt(g / a)
    X_star_an  = 1.0 + 2.0 * math.sqrt(2.0) * T_star   # M&M shallow-water

    # Point check at T* = 2 (informational M&M comparison)
    idx_T2   = int(np.argmin(np.abs(T_star - 2.0)))
    X_sph_T2 = float(X_star_sph[idx_T2])
    X_an_T2  = float(X_star_an[idx_T2])
    err_T2   = abs(X_sph_T2 - X_an_T2) / X_an_T2 * 100.0

    # Regression pass criterion: front must have advanced by at least 30% of L0
    # beyond its initial position by T*=2.  This checks gravity + Newtonian flow
    # (not arrested like Bingham), without requiring inviscid M&M accuracy.
    # Re ≈ sqrt(g*H0)*L0/nu ~ 20, so viscous effects are significant.
    XSTAR_THRESHOLD = 1.3   # X* > 1.3 means front moved at least 0.3*L0

    print(f'\n  ── Martin-Moyce front check (Newtonian regression, τ_y=0) ──')
    print(f'     Re   ≈ {rho0*U_ref_A*L0/mu_p:.0f}  (viscous; M&M is inviscid — large deviation expected)')
    print(f'     T*(end)   = {T_star[-1]:.2f}')
    print(f'     At T*=2:  X*_SPH = {X_sph_T2:.3f},  X*_theory (M&M) = {X_an_T2:.3f},  err = {err_T2:.1f} %')
    print(f'     Pass criterion: X*(T*=2) > {XSTAR_THRESHOLD:.1f}  (front advancing; not arrested like Bingham)')
    return X_sph_T2


# ─── Post-processing: Case B — Bingham Herschel profile ──────────────────────
def _check_bingham_slump(dumpname):
    """Compare arrested mound shape to Herschel parabolic equilibrium profile."""
    with gsd.hoomd.open(dumpname, 'r') as traj:
        snap = traj[-1]
        pos  = snap.particles.position
        tid  = snap.particles.typeid

    fluid = (tid == 0)
    x_f   = pos[fluid, 0]
    y_f   = pos[fluid, 1]

    # Bin fluid particles by x column (width dx); h_sph = max-y + dx/2
    nx_box_half = x_box / 2.0
    n_bins   = int(round(x_box / dx))
    bin_edges = np.linspace(-nx_box_half, nx_box_half, n_bins + 1)

    x_centres = []
    h_sph_vals = []
    for i in range(n_bins):
        in_bin = (x_f >= bin_edges[i]) & (x_f < bin_edges[i + 1])
        if np.any(in_bin):
            h_val = float(np.max(y_f[in_bin])) + dx / 2.0
            h_sph_vals.append(max(h_val, 0.0))
            x_centres.append(0.5 * (bin_edges[i] + bin_edges[i + 1]))

    h_sph   = np.array(h_sph_vals)
    x_cents = np.array(x_centres)

    # Central height: average over columns with |x| < 2*dx
    near_centre = np.abs(x_cents) < 2.0 * dx
    if np.any(near_centre):
        h_c_sph = float(np.mean(h_sph[near_centre]))
    else:
        h_c_sph = float(np.max(h_sph)) if len(h_sph) > 0 else float('nan')

    # Front position: max |x| with any fluid particle above y = dx
    above_floor = y_f > dx
    if np.any(above_floor):
        x_f_sph = float(np.max(np.abs(x_f[above_floor])))
    else:
        x_f_sph = 0.0

    # Profile L₂: compare h_sph vs analytical over active columns
    h_an_cents = h_analytical(x_cents)
    active = (h_sph > 0.01 * H0) | (h_an_cents > 0.01 * H0)
    if np.any(active):
        L2_profile = (np.sqrt(np.mean((h_sph[active] - h_an_cents[active])**2))
                      / H0 * 100.0)
    else:
        L2_profile = float('nan')

    err_hc = abs(h_c_sph - h_c_exact) / h_c_exact * 100.0
    err_xf = abs(x_f_sph - x_f_exact) / x_f_exact * 100.0

    print(f'\n  ── Bingham slump check (step {snap.configuration.step}) ──')
    print(f'     h_c : SPH = {h_c_sph*1e3:.3f} mm,  exact = {h_c_exact*1e3:.3f} mm,  err = {err_hc:.1f} %')
    print(f'     x_f : SPH = {x_f_sph*1e3:.1f} mm,  exact = {x_f_exact*1e3:.1f} mm,  err = {err_xf:.1f} %')
    print(f'     Profile L₂ / H0 = {L2_profile:.2f} %')
    return err_hc, err_xf, L2_profile


# ─── Reference velocity for Re print (Case A info line) ──────────────────────
U_ref_A = math.sqrt(g * H0)   # = U_char for Case A

# ─── Run Case A: Newtonian regression (τ_y = 0) ──────────────────────────────
err_A = run_case('caseA_newtonian', 0.0, steps_A)

# ─── Run Case B: Bingham slump (τ_y = 50 Pa) ─────────────────────────────────
err_hc, err_xf, L2_B = run_case('caseB_bingham', tau_y, steps_B)

# ─── Summary ─────────────────────────────────────────────────────────────────
XSTAR_THRESHOLD = 1.3   # front must advance at least 0.3*L0 by T*=2
print(f'\n{"═"*60}')
print(f'  BENCHMARK SUMMARY (num_length={num_length})')
print(f'{"═"*60}')
print(f'  Case A  τ_y=0 (Newtonian regression):')
print(f'    X*(T*=2) = {err_A:.3f}  (threshold > {XSTAR_THRESHOLD:.1f}; confirms spreading, not arrest)')
print(f'  Case B  Bingham τ_y={tau_y} Pa:')
print(f'    h_c error    = {err_hc:.1f} %  (threshold < 25 %)')
print(f'    x_f error    = {err_xf:.1f} %  (threshold < 30 %)')
print(f'    Profile L₂   = {L2_B:.2f} %   (threshold < 20 %)')

pass_A = err_A > XSTAR_THRESHOLD
pass_B = (err_hc < 25.0) and (err_xf < 30.0) and (L2_B < 20.0)
result = 'PASS' if (pass_A and pass_B) else 'FAIL'
print(f'  Result: {result}')

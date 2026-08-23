#!/usr/bin/env python3
r"""
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

Pore-scale density-driven fingering — single-mode brine/water benchmark
(SinglePhaseFlowGDGD, miscible, Boussinesq buoyancy, no-flux grain BC).

BENCHMARK DESCRIPTION
---------------------
A layer of brine (concentration C = 1, heavier by $\Delta\rho = \rho_0 |\beta_s| \Delta C$)
rests on top of water (C = 0) inside a square array of cylindrical grains
(see create_input_geometry.py).  The interface is perturbed with a single
cosine mode spanning the periodic width:

    $y_\mathrm{int}(x) = \delta_0 \cos(2\pi x / L_x)$

Gravity destabilises the layering; the finger grows through the pore space.
On the Darcy scale, linear stability of a sharp interface between two fluids
in a porous medium (equal viscosities) predicts exponential growth
$A(t) = \delta_0 e^{\sigma t}$ with

    $\sigma = \frac{K \Delta\rho g \, k}{\varphi(\mu_1+\mu_2)} - D k^2
            = \frac{U_b k}{2\varphi} - D k^2$,
    $U_b = K \Delta\rho g / \mu$,   $k = 2\pi/L_x$

(the factor 2 arises because the pressure perturbation must drive flow in
BOTH layers; buoyancy velocity $U_b$ from the measured permeability K of
run_permeability.py; the $-Dk^2$ term is the diffusive stabilisation of the
miscible interface).  The pore-scale simulation resolves individual grains
($d/\lambda = 0.1$), so agreement within a few tens of percent is the
expected outcome — the benchmark reports the measured/theoretical ratio.
Shake-down result at 5 dx per grain diameter: $\sigma_{meas}/\sigma_{theory}
\approx 0.95$ once the ramp/transient window is excluded.

The light-over-heavy configuration (--stable) must remain stable: the mode
amplitude may only decay.

CONVENTIONS
-----------
Boussinesq buoyancy (SinglePhaseFlowGDGD.cc): $\Delta F_b = m g \,(-\beta_s (C - C_\mathrm{ref}))$.
With $g_y < 0$ and $\beta_s < 0$: C = 1 gets a extra downward force -> heavier. Here
$\beta_s = -0.05$, $\Delta C = 1$, $C_\mathrm{ref} = 0.5$  =>  $\Delta\rho/\rho_0 = 5\%$ (Boussinesq-valid).

Scalar wall BC = 'noflux': grains are impermeable to solute (zero-flux),
they do NOT act as concentration sources (unlike the default Dirichlet BC).

PARAMETER FEASIBILITY (real brine)
----------------------------------
Real NaCl brine has D ~ 1.6e-9 m^2/s (Sc ~ 600), unreachable at pore-scale
SPH resolutions (grid Peclet >> 1).  The benchmark therefore scales D up to
5e-7 m^2/s and compensates via geometry/viscosity so that the governing
dimensionless numbers (density contrast 5 %, U_b*dx/D < 1, sigma > 0) stay in
the physically relevant regime.  VRD mode (boussinesq=False) exists for larger
density contrasts but is not yet validated.

Usage:
    python3 run_fingering_single_mode.py <init_gsd> <K_m2> [steps] [--stable]

    init_gsd : porous_finger_*_init.gsd from create_input_geometry.py
    K_m2     : measured permeability [m^2] from run_permeability.py
               (fallback: Drummond-Tahir estimate, pass 'DT')
    steps    : simulation steps (default: enough for sigma*t = 4)
    --stable : swap layers (light brine analog on bottom) -> decay check
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

args     = [a for a in sys.argv[1:] if a != '--stable']
stable   = '--stable' in sys.argv[1:]
filename = str(args[0])

dt_string = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
suffix    = '_stable_run' if stable else '_run'
logname   = filename.replace('_init.gsd', f'{suffix}.log')
dumpname  = filename.replace('_init.gsd', f'{suffix}.gsd')

sim.create_state_from_gsd(filename=filename, domain_decomposition=(None, None, 1))

# ─── Physical parameters ─────────────────────────────────────────────────────
# Geometry constants (must match create_input_geometry.py)
a_grain = 0.8e-3            # grain lattice spacing            [m]
R_grain = 0.2e-3            # grain radius                     [m]
d_grain = 2.0 * R_grain
throat  = a_grain - d_grain # pore-throat width                [m]

rho0      = 1000.0          # rest density                     [kg/m³]
viscosity = 0.001           # dynamic viscosity                [Pa·s]
gy        = -9.81           # gravity                          [m/s²]
beta_s    = -0.05           # solutal expansion coefficient    [–]  (<0: C=1 heavy)
kappa_s   = 5.0e-7          # solute diffusivity D             [m²/s]
C_hi      = 1.0             # brine concentration
C_lo      = 0.0             # water concentration
C_ref     = 0.5 * (C_hi + C_lo)
DeltaC    = C_hi - C_lo
drho_rel  = 0.01            # allowed (acoustic) density variation [–]
backpress = 0.01            # background pressure coefficient      [–]
damp      = 500             # body-force ramp [steps] — buoyancy ramps consistently

Delta_rho = rho0 * abs(beta_s) * DeltaC   # 50 kg/m³ → Δρ/ρ0 = 5 %

# Permeability: measured (run_permeability.py) or Drummond–Tahir fallback
c_solid = np.pi * R_grain**2 / a_grain**2
K_DT = R_grain**2 / (8.0 * c_solid) * (
    -np.log(c_solid) - 1.476 + 2.0 * c_solid
    - 1.774 * c_solid**2 + 4.076 * c_solid**3)
K = K_DT if args[1] == 'DT' else float(args[1])

# ─── Box & derived quantities ────────────────────────────────────────────────
box = sim.state.box
Lx  = box.Lx
phi_val   = 0.0
vsize_val = 0.0
ly_fluid_val = 0.0
snapshot = sim.state.get_snapshot()
if snapshot.communicator.rank == 0:
    vsize_val = float(np.max(snapshot.particles.slength)) \
                / hoomd.sph.kernel.OptimalH['WendlandC4']
    y   = snapshot.particles.position[:, 1]
    tid = snapshot.particles.typeid
    # fluid region ends where the wall layers begin
    ly_fluid_val = 2.0 * float(np.max(np.abs(y[tid == 0]))) + vsize_val
    in_lat = np.abs(y) < ly_fluid_val / 2
    phi_val = float(np.sum((tid == 0) & in_lat)) / float(np.sum(in_lat))

phi      = device.communicator.bcast_double(phi_val)
dx       = device.communicator.bcast_double(vsize_val)
ly_fluid = device.communicator.bcast_double(ly_fluid_val)

delta0 = 0.5 * R_grain      # initial interface perturbation amplitude [m]
k_wave = 2.0 * np.pi / Lx

# Darcy-scale linear theory (two-fluid, equal viscosities: factor 1/2)
U_b       = K * Delta_rho * abs(gy) / viscosity              # buoyancy velocity  [m/s]
sigma_lin = U_b * k_wave / (2.0 * phi) - kappa_s * k_wave**2 # growth rate        [1/s]
refvel    = max(U_b / phi, 1e-10)                            # interstitial scale [m/s]
Ra_like   = U_b / (2.0 * phi * kappa_s * k_wave)             # >1 → unstable

steps = int(args[2]) if len(args) > 2 else \
    (int(4.0 / max(sigma_lin, 1e-3) / 1e-5) if sigma_lin > 0 else 50001)

if device.communicator.rank == 0:
    Pe_grid = U_b * dx / kappa_s
    print(f'Porous fingering  ({"STABLE (decay check)" if stable else "unstable single mode"})')
    print(f'  Lx = {Lx*1e3:.2f} mm,  fluid height = {ly_fluid*1e3:.2f} mm,  '
          f'dx = {dx*1e6:.1f} µm,  φ = {phi:.4f}')
    print(f'  K = {K:.4e} m²  (Drummond–Tahir: {K_DT:.4e} m²)')
    print(f'  Δρ = {Delta_rho:.1f} kg/m³ ({abs(beta_s)*DeltaC*100:.1f} %),  D = {kappa_s:.1e} m²/s')
    print(f'  U_b = {U_b:.4e} m/s,  grid Péclet U_b·dx/D = {Pe_grid:.2f}')
    print(f'  σ_theory = U_b·k/(2φ) − D·k² = {sigma_lin:.3f} s⁻¹  '
          f'(U_b·k/(2φ) / D·k² = {Ra_like:.1f})')
    if sigma_lin > 0:
        print(f'  τ = 1/σ = {1.0/sigma_lin*1e3:.1f} ms')

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

# ─── SinglePhaseFlowGDGD model ───────────────────────────────────────────────
model = hoomd.sph.sphmodel.SinglePhaseFlowGDGD(
    kernel=kernel_obj, eos=eos, nlist=nlist,
    fluidgroup_filter=filterfluid, solidgroup_filter=filtersolid,
    densitymethod='SUMMATION',
    kappa_s=kappa_s,
    beta_s=beta_s,
    scalar_ref=C_ref,
    boussinesq=True,
    scalar_wall_bc='noflux',
)

model.mu                  = viscosity
model.gx                  = 0.0
model.gy                  = gy
model.gz                  = 0.0
model.damp                = damp
model.artificialviscosity = True
model.alpha               = 0.1
model.beta                = 0.0
model.densitydiffusion    = False

# ─── Speed of sound & timestep ───────────────────────────────────────────────
maximum_smoothing_length = sph_helper.set_max_sl(sim, device, model)

# LREF = fluid height: hydrostatic pressure builds over the FULL stratified
# column, so the gravity-wave condition c² ≥ g·LREF/Δρ must use the fluid
# height — a pore-scale LREF lets the column compress by (H/LREF)·Δρ
# (verified: 8.7 % density overshoot with LREF = throat width).
c, cond = model.compute_speedofsound(
    LREF=ly_fluid, UREF=refvel, DX=dx, DRHO=drho_rel,
    H=maximum_smoothing_length, MU=viscosity, RHO0=rho0)

eos.set_speedofsound(c)

dt, dt_cond = model.compute_dt(
    LREF=ly_fluid, UREF=refvel, DX=dx, DRHO=drho_rel,
    H=maximum_smoothing_length, MU=viscosity, RHO0=rho0)

if device.communicator.rank == 0:
    print(f'Speed of sound: {c:.4f} m/s  ({cond})')
    print(f'Timestep: {dt:.3e} s  ({dt_cond})')
    if sigma_lin > 0:
        print(f'  steps per e-fold: {int(1.0/(sigma_lin*dt))},  '
              f'running {steps} steps = σ·t = {sigma_lin*steps*dt:.2f}')

# ─── Integrator ──────────────────────────────────────────────────────────────
integrator = hoomd.sph.Integrator(dt=dt)
vvb = hoomd.sph.methods.VelocityVerletBasic(filter=filterfluid,
                                              densitymethod='SUMMATION')
integrator.methods.append(vvb)
integrator.forces.append(model)
sim.operations.integrator = integrator

# ─── Initialise concentration field C (aux4.x) ───────────────────────────────
# Unstable: brine (C = 1, heavy) above the perturbed interface, water below.
# Stable (--stable): layers swapped.  Solids get C_ref (excluded from
# diffusion by the no-flux BC, but keeps the Boussinesq reference clean).
with sim.state.cpu_local_snapshot as snap:
    pos  = snap.particles.position[:]
    tid  = snap.particles.typeid[:]
    aux4 = snap.particles.auxiliary4

    x_p   = np.array(pos[:, 0])
    y_p   = np.array(pos[:, 1])
    y_int = delta0 * np.cos(k_wave * x_p)

    is_fluid = np.array(tid) == 0
    is_upper = is_fluid & (y_p >  y_int)
    is_lower = is_fluid & (y_p <= y_int)

    aux4[:, 0] = C_ref
    if stable:
        aux4[is_upper, 0] = C_lo   # light on top → stable
        aux4[is_lower, 0] = C_hi
    else:
        aux4[is_upper, 0] = C_hi   # brine on top → unstable
        aux4[is_lower, 0] = C_lo

# ─── Output ──────────────────────────────────────────────────────────────────
try:
    os.remove(dumpname)
except (FileNotFoundError, OSError):
    pass
device.communicator.barrier()

gsd_period = max(1, steps // 100)   # ≈ 100 frames
gsd_writer = hoomd.write.GSD(filename=dumpname,
                              trigger=hoomd.trigger.Periodic(gsd_period),
                              mode='wb',
                              dynamic=['property', 'momentum'])
sim.operations.writers.append(gsd_writer)

logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps', 'walltime'])
compute_fluid = hoomd.sph.compute.SinglePhaseFlowBasicProperties(filter=filterfluid)
sim.operations.computes.append(compute_fluid)
logger.add(compute_fluid, quantities=['e_kin_fluid', 'mean_density'])
log_period = max(1, steps // 200)
table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(log_period),
                          logger=logger, max_header_len=10)
sim.operations.writers.append(table)

log_file   = open(logname, mode='w+', newline='\n')
table_file = hoomd.write.Table(output=log_file,
                                trigger=hoomd.trigger.Periodic(log_period),
                                logger=logger, max_header_len=10)
sim.operations.writers.append(table_file)

# ─── Run ─────────────────────────────────────────────────────────────────────
if device.communicator.rank == 0:
    print(f'Starting porous fingering run at {dt_string}')
    print(f'  δ₀ = {delta0*1e6:.1f} µm = {delta0/dx:.1f} dx,  damp = {damp} steps')

sim.run(steps, write_at_start=True)
gsd_writer.flush()

# ─── Post-processing: mode amplitude & growth rate vs Darcy theory ───────────
# The interface height h(x, t) is extracted per x-bin as the y-level where the
# bin-averaged concentration profile crosses C_ref; the single-mode amplitude
# is its cosine Fourier coefficient  A(t) = |(2/Lx) ∫ h(x) cos(kx) dx|.
if device.communicator.rank == 0:
    nbin_x = 32
    xedges = np.linspace(-Lx / 2, Lx / 2, nbin_x + 1)
    xmid   = 0.5 * (xedges[1:] + xedges[:-1])
    nbin_y = 80
    yedges = np.linspace(-ly_fluid / 2, ly_fluid / 2, nbin_y + 1)
    ymid   = 0.5 * (yedges[1:] + yedges[:-1])

    def mode_amplitude(pos_s, C_s):
        """Fourier amplitude of the C = C_ref contour of the binned C field."""
        h = np.full(nbin_x, np.nan)
        ix = np.clip(np.digitize(pos_s[:, 0], xedges) - 1, 0, nbin_x - 1)
        iy = np.clip(np.digitize(pos_s[:, 1], yedges) - 1, 0, nbin_y - 1)
        csum = np.zeros((nbin_x, nbin_y))
        cnt  = np.zeros((nbin_x, nbin_y))
        np.add.at(csum, (ix, iy), C_s)
        np.add.at(cnt,  (ix, iy), 1.0)
        with np.errstate(invalid='ignore'):
            cbar = csum / cnt
        for b in range(nbin_x):
            prof = cbar[b]
            ok   = np.isfinite(prof)
            if ok.sum() < 4:
                continue
            yv, cv = ymid[ok], prof[ok]
            # topmost crossing of C_ref moving downward (brine sits on top)
            cross = np.where(np.diff(np.sign(cv - C_ref)) != 0)[0]
            if len(cross) == 0:
                continue
            j = cross[np.argmin(np.abs(yv[cross]))]   # crossing nearest y = 0
            y1, y2, c1, c2 = yv[j], yv[j + 1], cv[j], cv[j + 1]
            h[b] = y1 + (C_ref - c1) * (y2 - y1) / (c2 - c1)
        ok = np.isfinite(h)
        if ok.sum() < nbin_x // 2:
            return np.nan
        hm = h[ok] - np.mean(h[ok])
        return abs(2.0 * np.mean(hm * np.cos(k_wave * xmid[ok])))

    frames_t, frames_A, frames_mass = [], [], []
    with gsd.hoomd.open(dumpname, 'r') as traj:
        for snap in traj:
            t_step = snap.configuration.step
            tid_s  = snap.particles.typeid
            fluid  = tid_s == 0
            pos_s  = snap.particles.position[fluid]
            C_s    = snap.particles.auxiliary4[fluid, 0]
            m_s    = snap.particles.mass[fluid]
            frames_t.append(t_step * dt)
            frames_A.append(mode_amplitude(pos_s, C_s))
            frames_mass.append(float(np.sum(m_s * C_s)))

    frames_t    = np.array(frames_t)
    frames_A    = np.array(frames_A)
    frames_mass = np.array(frames_mass)

    print(f'\n── Porous fingering check ──')
    print(f'  σ_theory = {sigma_lin:.3f} s⁻¹,  δ₀ = {delta0*1e6:.1f} µm')
    print(f'  {"t [ms]":>9}  {"A [µm]":>9}  {"A/δ₀":>7}  {"σ·t":>6}')
    for i in range(len(frames_t)):
        print(f'  {frames_t[i]*1e3:>9.2f}  {frames_A[i]*1e6:>9.2f}'
              f'  {frames_A[i]/delta0:>7.3f}  {max(sigma_lin,0)*frames_t[i]:>6.2f}')

    # solute conservation (no-flux BC + Lagrangian advection)
    drift = abs(frames_mass[-1] - frames_mass[0]) / frames_mass[0]
    print(f'\n  Solute content Σ m·C drift: {drift*100:.3f} %  (pass: < 1 %)')

    if stable:
        A_ok = frames_A[np.isfinite(frames_A)]
        print(f'  STABLE config: A start {A_ok[0]*1e6:.1f} µm → end {A_ok[-1]*1e6:.1f} µm '
              f'(pass: no growth)')
    else:
        # fit sigma over the established exponential window: skip the body-
        # force ramp AND the flow-establishment transient (σ·t < 0.3, during
        # which the Darcy velocity field builds up and diffusion dominates
        # the contour position), cap at 0.2/k (linear-regime limit)
        t_ramp  = damp * dt
        t_start = max(2.0 * t_ramp, 0.3 / max(sigma_lin, 1e-6))
        valid = (np.isfinite(frames_A) & (frames_t > t_start)
                 & (frames_A < 0.2 / k_wave))
        if valid.sum() >= 3:
            p = np.polyfit(frames_t[valid], np.log(frames_A[valid]), 1)
            sigma_meas = p[0]
            print(f'  σ_measured = {sigma_meas:.3f} s⁻¹  '
                  f'(fit window t > {t_start*1e3:.1f} ms)'
                  f'  →  σ_meas/σ_theory = {sigma_meas/sigma_lin:.2f}'
                  f'  (pass: 0.7 – 1.3, d/λ = {d_grain*k_wave/(2*np.pi):.2f} '
                  f'scale separation is marginal)')
        else:
            print('  (not enough frames in the exponential window for a σ fit —'
                  ' run longer or dump more frames)')

if HAS_VTU and device.communicator.rank == 0:
    export_gsd2vtu.export_gdgd(dumpname)

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

Multi-mode pore-scale density fingering at controlled Rayleigh number
(brine over water in a disordered cylinder pack, SinglePhaseFlowGDGD).

BENCHMARK DESCRIPTION
---------------------
An initially FLAT brine/water interface (no imposed perturbation — fingers are
seeded by the grain disorder) in the pack from create_input_geometry.py.
The governing dimensionless group is the (solutal) Rayleigh number

    $Ra = \Delta\rho \, g \, K \, H / (\varphi \, \mu \, D)$

with the measured permeability K (run_permeability.py), fluid height H, and
solute diffusivity D chosen per run.  Sweep D (see run_all_ra.sh) to sample
e.g. Ra ~= 200 / 1000 / 5000.

MEASURED QUANTITIES (written to <run>.npz + printed)
----------------------------------------------------
  onset time   : departure of interface-slice concentration variance from the
                 pure-diffusion (erf-profile) solution
  wavelength   : dominant mode of the transverse FFT of C(x) in a slice one
                 mean grain diameter below the initial interface
  mixing width : W(t) = y(<C> = 0.99) - y(<C> = 0.01) of the transverse-
                 averaged profile <C>(y, t); diffusive W ~ sqrt(t) at low Ra,
                 convective W ~ t once fingers develop (expected for Ra >~ 1000)

Expected trends: onset time and dominant wavelength decrease with Ra; the
--stable configuration (light layer on top) stays purely diffusive.

APPLICATION CONTEXT
-------------------
Convective dissolution of CO2 in brine aquifers, seawater intrusion, and
density-driven mixing in underground H2 storage (see also the two-phase
benchmark 01_twophaseflow_benchmarks/08_h2_bubble_brine_shear).  Real brine
(D ~ 1.6e-9 m^2/s, Sc ~ 600) is out of reach at pore-scale resolution: D is
scaled up while Ra is held at the target value; the density contrast stays at
5% (Boussinesq-valid). VRD mode exists for larger contrasts (not yet validated).

Usage:
    python3 run_fingering_ra.py <init_gsd> <K_m2> <D_m2s> [steps] [--stable] [--mu MU]

    init_gsd : porous_ra_pack_*_init.gsd (or *_small_*)
    K_m2     : measured permeability [m^2] from run_permeability.py
    D_m2s    : solute diffusivity [m^2/s] (sets Ra)
    --mu MU  : dynamic viscosity [Pa s] (default 0.001)

RESOLUTION CONSTRAINT: the grid Peclet is tied to Ra by the identity
    Pe_grid = U_b dx / D = Ra * phi * dx / H     (independent of mu and D!)
so the maximum resolvable Rayleigh number at a given resolution is
    Ra_max ~ 4 H / (phi dx).
Full pack at 10 dx/diameter (H = 12 mm, dx = 40 um): Ra_max ~ 1600;
Ra = 5000 needs 20 dx/diameter or a taller domain.
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

argv   = sys.argv[1:]
stable = '--stable' in argv
viscosity = 0.001
if '--mu' in argv:
    viscosity = float(argv[argv.index('--mu') + 1])
    del argv[argv.index('--mu') + 1]
args = [a for a in argv if not a.startswith('--')]

filename = str(args[0])
K        = float(args[1])
kappa_s  = float(args[2])          # solute diffusivity D [m²/s]

dt_string = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
suffix    = f'_D{kappa_s:.2e}' + ('_stable' if stable else '') + '_run'
logname   = filename.replace('_init.gsd', f'{suffix}.log')
dumpname  = filename.replace('_init.gsd', f'{suffix}.gsd')

sim.create_state_from_gsd(filename=filename, domain_decomposition=(None, None, 1))

# ─── Physical parameters ─────────────────────────────────────────────────────
d_mean    = 0.4e-3          # mean grain diameter (matches create script) [m]
rho0      = 1000.0          # rest density                     [kg/m³]
gy        = -9.81           # gravity                          [m/s²]
beta_s    = -0.05           # solutal expansion coefficient    [–]  (<0: C=1 heavy)
C_hi, C_lo = 1.0, 0.0
C_ref     = 0.5 * (C_hi + C_lo)
DeltaC    = C_hi - C_lo
drho_rel  = 0.01
backpress = 0.01
damp      = 500

Delta_rho = rho0 * abs(beta_s) * DeltaC   # 50 kg/m³

# ─── Box & derived quantities ────────────────────────────────────────────────
box = sim.state.box
Lx  = box.Lx
phi_val = vsize_val = ly_fluid_val = 0.0
snapshot = sim.state.get_snapshot()
if snapshot.communicator.rank == 0:
    vsize_val = float(np.max(snapshot.particles.slength)) \
                / hoomd.sph.kernel.OptimalH['WendlandC4']
    y   = snapshot.particles.position[:, 1]
    tid = snapshot.particles.typeid
    ly_fluid_val = 2.0 * float(np.max(np.abs(y[tid == 0]))) + vsize_val
    in_lat = np.abs(y) < ly_fluid_val / 2
    phi_val = float(np.sum((tid == 0) & in_lat)) / float(np.sum(in_lat))

phi      = device.communicator.bcast_double(phi_val)
dx       = device.communicator.bcast_double(vsize_val)
ly_fluid = device.communicator.bcast_double(ly_fluid_val)
H        = ly_fluid

U_b = K * Delta_rho * abs(gy) / viscosity          # buoyancy velocity [m/s]
Ra  = Delta_rho * abs(gy) * K * H / (phi * viscosity * kappa_s)
refvel = max(U_b / phi, 1e-10)

# characteristic convective time scale (Darcy): t_c = phi H / U_b
t_conv = phi * H / max(U_b, 1e-12)
steps  = int(args[3]) if len(args) > 3 else 50001

if device.communicator.rank == 0:
    Pe_grid = U_b * dx / kappa_s
    print(f'Porous Ra fingering  ({"STABLE (decay check)" if stable else "multi-mode"})')
    print(f'  Lx = {Lx*1e3:.2f} mm,  H = {H*1e3:.2f} mm,  dx = {dx*1e6:.1f} µm,  φ = {phi:.4f}')
    print(f'  K = {K:.4e} m²,  D = {kappa_s:.2e} m²/s,  μ = {viscosity:.4f} Pa·s')
    print(f'  Δρ = {Delta_rho:.1f} kg/m³,  U_b = {U_b:.4e} m/s')
    print(f'  Ra = ΔρgKH/(φμD) = {Ra:.0f}')
    print(f'  grid Péclet U_b·dx/D = {Pe_grid:.2f}'
          + ('   ** WARNING: > 4, raise --mu or D **' if Pe_grid > 4 else ''))
    print(f'  convective time t_c = φH/U_b = {t_conv:.3f} s')

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

model = hoomd.sph.sphmodel.SinglePhaseFlowGDGD(
    kernel=kernel_obj, eos=eos, nlist=nlist,
    fluidgroup_filter=filterfluid, solidgroup_filter=filtersolid,
    densitymethod='SUMMATION',
    kappa_s=kappa_s, beta_s=beta_s, scalar_ref=C_ref,
    boussinesq=True, scalar_wall_bc='noflux')

model.mu                  = viscosity
model.gx, model.gy, model.gz = 0.0, gy, 0.0
model.damp                = damp
model.artificialviscosity = True
model.alpha               = 0.1
model.beta                = 0.0
model.densitydiffusion    = False

# ─── Speed of sound & timestep ───────────────────────────────────────────────
# LREF = fluid height: the hydrostatic pressure of the stratified column sets
# the required stiffness (a pore-scale LREF lets the column compress by
# (H/LREF)·Δρ — verified as a several-percent density overshoot).
maximum_smoothing_length = sph_helper.set_max_sl(sim, device, model)

c, cond = model.compute_speedofsound(
    LREF=H, UREF=refvel, DX=dx, DRHO=drho_rel,
    H=maximum_smoothing_length, MU=viscosity, RHO0=rho0)
eos.set_speedofsound(c)

dt, dt_cond = model.compute_dt(
    LREF=H, UREF=refvel, DX=dx, DRHO=drho_rel,
    H=maximum_smoothing_length, MU=viscosity, RHO0=rho0)

if device.communicator.rank == 0:
    print(f'Speed of sound: {c:.4f} m/s  ({cond})')
    print(f'Timestep: {dt:.3e} s  ({dt_cond})')
    print(f'  simulated time: {steps*dt:.3f} s = {steps*dt/t_conv:.2f} t_c')

integrator = hoomd.sph.Integrator(dt=dt)
vvb = hoomd.sph.methods.VelocityVerletBasic(filter=filterfluid,
                                              densitymethod='SUMMATION')
integrator.methods.append(vvb)
integrator.forces.append(model)
sim.operations.integrator = integrator

# ─── Initialise concentration: flat layered interface at y = 0 ───────────────
with sim.state.cpu_local_snapshot as snap:
    pos  = snap.particles.position[:]
    tid  = snap.particles.typeid[:]
    aux4 = snap.particles.auxiliary4

    y_p      = np.array(pos[:, 1])
    is_fluid = np.array(tid) == 0
    is_upper = is_fluid & (y_p > 0.0)
    is_lower = is_fluid & (y_p <= 0.0)

    aux4[:, 0] = C_ref
    if stable:
        aux4[is_upper, 0] = C_lo
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

gsd_period = max(1, steps // 100)
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
    print(f'Starting Ra fingering run at {dt_string}  (Ra = {Ra:.0f})')
sim.run(steps, write_at_start=True)
gsd_writer.flush()

# ─── Post-processing ─────────────────────────────────────────────────────────
if device.communicator.rank == 0:
    nbin_x = 64
    xedges = np.linspace(-Lx / 2, Lx / 2, nbin_x + 1)
    nbin_y = 120
    yedges = np.linspace(-H / 2, H / 2, nbin_y + 1)
    ymid   = 0.5 * (yedges[1:] + yedges[:-1])

    y_slice = -d_mean          # slice one mean grain diameter below interface
    slice_halfwidth = d_mean / 2

    ts, Wmix, var_slice, spectra, profiles = [], [], [], [], []
    kx = np.fft.rfftfreq(nbin_x, d=Lx / nbin_x) * 2.0 * np.pi   # [rad/m]

    with gsd.hoomd.open(dumpname, 'r') as traj:
        for snap in traj:
            t = snap.configuration.step * dt
            tid_s = snap.particles.typeid
            fluid = tid_s == 0
            pos_s = snap.particles.position[fluid]
            C_s   = snap.particles.auxiliary4[fluid, 0]

            # transverse-averaged profile <C>(y)
            iy = np.clip(np.digitize(pos_s[:, 1], yedges) - 1, 0, nbin_y - 1)
            csum = np.bincount(iy, weights=C_s, minlength=nbin_y)
            cnt  = np.bincount(iy, minlength=nbin_y).astype(float)
            with np.errstate(invalid='ignore'):
                cbar = csum / cnt

            # mixing-zone width W(t): span where 0.01 < <C> < 0.99
            ok = np.isfinite(cbar)
            mixed = ok & (cbar > 0.01) & (cbar < 0.99)
            W = float(ymid[mixed].max() - ymid[mixed].min()) if mixed.any() else 0.0

            # slice below the interface: C(x) → variance + FFT spectrum
            in_slice = np.abs(pos_s[:, 1] - y_slice) < slice_halfwidth
            ix = np.clip(np.digitize(pos_s[in_slice, 0], xedges) - 1, 0, nbin_x - 1)
            cs = np.bincount(ix, weights=C_s[in_slice], minlength=nbin_x)
            cn = np.bincount(ix, minlength=nbin_x).astype(float)
            with np.errstate(invalid='ignore'):
                cx = np.where(cn > 0, cs / cn, 0.0)
            var = float(np.var(cx))
            spec = np.abs(np.fft.rfft(cx - cx.mean())) / nbin_x

            ts.append(t)
            Wmix.append(W)
            var_slice.append(var)
            spectra.append(spec)
            profiles.append(cbar)

    ts        = np.array(ts)
    Wmix      = np.array(Wmix)
    var_slice = np.array(var_slice)
    spectra   = np.array(spectra)
    profiles  = np.array(profiles)

    # onset time: slice variance exceeds a small threshold (pure diffusion
    # keeps the slice at C ~ 0 until the diffusive front arrives:
    # t_diff = (d_mean)^2 / (4 D) — fingers arriving earlier signal convection)
    thresh = 1e-3
    onset_idx = np.argmax(var_slice > thresh) if np.any(var_slice > thresh) else -1
    t_onset = ts[onset_idx] if onset_idx >= 0 else np.nan
    t_diff  = d_mean**2 / (4.0 * kappa_s)

    # dominant wavelength at the last frame with meaningful spectrum
    lam_dom = np.nan
    for spec in spectra[::-1]:
        if spec[1:].max() > 1e-3:
            k_dom = kx[1 + int(np.argmax(spec[1:]))]
            lam_dom = 2.0 * np.pi / k_dom
            break

    print(f'\n── Ra fingering summary  (Ra = {Ra:.0f}) ──')
    print(f'  {"t [s]":>8}  {"W_mix [mm]":>11}  {"var_slice":>10}')
    for i in range(len(ts)):
        print(f'  {ts[i]:>8.4f}  {Wmix[i]*1e3:>11.3f}  {var_slice[i]:>10.4e}')
    print(f'\n  onset time (slice variance > {thresh}):  t = {t_onset:.4f} s'
          f'   (pure diffusion needs t_diff ≈ {t_diff:.4f} s)')
    print(f'  dominant finger wavelength: λ = {lam_dom*1e3:.2f} mm'
          if np.isfinite(lam_dom) else '  no finger signal in slice spectrum')
    print(f'  mixing width: W(start) = {Wmix[0]*1e3:.2f} mm → W(end) = {Wmix[-1]*1e3:.2f} mm')
    if stable:
        print(f'  STABLE config: expect W ~ sqrt(t) (diffusive) and no onset before t_diff')

    npzname = dumpname.replace('.gsd', '.npz')
    np.savez(npzname, t=ts, W=Wmix, var_slice=var_slice, spectra=spectra,
             kx=kx, profiles=profiles, ymid=ymid, Ra=Ra, K=K, D=kappa_s,
             phi=phi, H=H, U_b=U_b, dt=dt)
    print(f'  time series written to {npzname}')

if HAS_VTU and device.communicator.rank == 0:
    export_gsd2vtu.export_gdgd(dumpname)

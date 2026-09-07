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

Bingham slug mobilization benchmark (06) — single run.

Runs one (tau_y, g_x) point: a Bingham slug in a Newtonian carrier inside a
periodic tube, body-force driven. The slug's mean axial velocity is sampled at
mid-run and at the end; an ARRESTED slug plateaus at the Papanastasiou creep
floor, a MOBILIZED slug is still accelerating on the mu_p viscous time scale
(rho r^2/mu_p >> t_sim), so the growth ratio v_end/v_mid is a calibration-free
classifier. Appends to "mobilization_results.txt":

    tau_y  g  v_mid  v_end  u_creep_pred  u_newt_pred

  u_creep_pred = rho g (L/l) r^2 / (8 mu_max)  — regularized creep floor
                 (the periodic pressure field concentrates the column drive
                 onto the slug in the arrested limit: G_slug = rho g L/l)
  u_newt_pred  = rho g r^2 / (8 mu_p)          — Newtonian Poiseuille scale

Usage:
    python3 run_mobilization.py <init.gsd> <tau_y> <g> [steps=5000]
"""

import sys, os
import numpy as np

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules'))

import hoomd
from hoomd import sph
import gsd.hoomd
import sph_helper

device = hoomd.device.CPU(notice_level=1)
sim = hoomd.Simulation(device=device)

initfile = sys.argv[1]
tau_y    = float(sys.argv[2])
gx       = float(sys.argv[3])
steps    = int(sys.argv[4]) if len(sys.argv) > 4 else 5000

sim.create_state_from_gsd(filename=initfile)

# ─── Physical parameters (must match create_slug_tube.py) ───────────────────
dx    = 5.0e-4
r_in  = 8 * dx
L     = sim.state.box.Lx
l_s   = 16 * dx
rho0  = 1000.0
mu1   = 0.1           # carrier (Newtonian)             [Pa s]
mu_p  = 0.1           # slug plastic viscosity          [Pa s]
drho  = 0.01
U_ref = 0.02          # velocity scale for regularization & c

gdot_c  = U_ref / r_in
m_reg   = 3.0 / gdot_c if tau_y > 0 else 0.0
mu_max2 = mu_p + tau_y * m_reg

g_crit  = 2.0 * tau_y * l_s / (rho0 * r_in * L)
u_creep = rho0 * gx * (L / l_s) * r_in**2 / (8.0 * mu_max2)
u_newt  = rho0 * gx * r_in**2 / (8.0 * mu_p)

if device.communicator.rank == 0:
    print(f"── Slug mobilization: tau_y={tau_y} Pa  g={gx} m/s²  "
          f"(g/g_crit={gx/g_crit if g_crit>0 else float('inf'):.2f}) ──")

# ─── Kernel / nlist ─────────────────────────────────────────────────────────
kernel     = "WendlandC4"
rcut       = sph.kernel.Kappa[kernel] * sph.kernel.OptimalH[kernel] * dx
kernel_obj = sph.kernel.Kernels[kernel]()
nlist = hoomd.nsearch.nlist.Cell(buffer=rcut*0.05, rebuild_check_delay=1,
                                 kappa=kernel_obj.Kappa())

eos1 = sph.eos.Tait(); eos1.set_params(rho0, 0.01)
eos2 = sph.eos.Tait(); eos2.set_params(rho0, 0.01)

fA = hoomd.filter.Type(["A"])
fB = hoomd.filter.Type(["B"])
fS = hoomd.filter.Type(["S"])

model = sph.sphmodel.TwoPhaseFlow(
    kernel=kernel_obj, eos1=eos1, eos2=eos2, nlist=nlist,
    fluidgroup1_filter=fA, fluidgroup2_filter=fB, solidgroup_filter=fS,
    densitymethod="SUMMATION", colorgradientmethod="DENSITYRATIO")

# mu2 carries the Papanastasiou-regularized viscosity ceiling so compute_dt
# respects the explicit stability limit (cf. benchmarks 02 and 05).
model.mu1     = mu1
model.mu2     = mu_max2
model.sigma12 = 0.0     # no interfacial tension: isolate the yield criterion
model.omega   = 90.0
model.gx      = gx
model.damp    = 500
model.artificialviscosity = True
model.alpha   = 0.1
model.beta    = 0.0
model.densitydiffusion = False

max_sl = sph_helper.set_max_sl(sim, device, model)
c1, cond1, c2, cond2 = model.compute_speedofsound(
    LREF=L, UREF=U_ref, DX=dx, DRHO=drho, H=max_sl,
    MU1=mu1, MU2=mu_max2, RHO01=rho0, RHO02=rho0, SIGMA12=0.0)
sph_helper.update_min_c0_tpf(device, model, c1, c2,
                             mode="plain", lref=L, uref=U_ref, cfactor=1.5)
dt, dt_cond = model.compute_dt(
    LREF=L, UREF=U_ref, DX=dx, DRHO=drho, H=max_sl,
    MU1=mu1, MU2=mu_max2, RHO01=rho0, RHO02=rho0, SIGMA12=0.0)

if device.communicator.rank == 0:
    print(f"   dt={dt:.3e} s ({dt_cond})  t_sim={steps*dt:.4f} s  "
          f"mu_max2={mu_max2:.3f} Pa·s")

integrator = sph.Integrator(dt=dt)
integrator.methods.append(sph.methods.VelocityVerletBasic(filter=fA, densitymethod="SUMMATION"))
integrator.methods.append(sph.methods.VelocityVerletBasic(filter=fB, densitymethod="SUMMATION"))
integrator.forces.append(model)
sim.operations.integrator = integrator

sim.run(0)
if tau_y > 0:
    model._cpp_obj.activateBingham2(mu_p, tau_y, m_reg, 0.0)

def slug_velocity():
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        B = snap.particles.typeid == 1
        return float(snap.particles.velocity[B, 0].mean())
    return 0.0

sim.run(steps // 2)
v_mid = slug_velocity()
sim.run(steps - steps // 2)
v_end = slug_velocity()

if sim.device.communicator.rank == 0:
    line = f"{tau_y:.4f} {gx:.6f} {v_mid:.6e} {v_end:.6e} {u_creep:.6e} {u_newt:.6e}"
    with open("mobilization_results.txt", "a") as f:
        f.write(line + "\n")
    growth = v_end / v_mid if v_mid > 0 else float("inf")
    print(f"   <v_x>_slug: mid={v_mid:.4e}  end={v_end:.4e}  growth={growth:.2f}  "
          f"(creep {u_creep:.3e}, Newtonian {u_newt:.3e})")

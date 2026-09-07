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

Viscoplastic invasion benchmark (05) — run script.

BENCHMARK DESCRIPTION
---------------------
A Newtonian phase A invades a sphere packing saturated with phase B (Newtonian
baseline for tau_y = 0, Bingham defender for tau_y > 0). Body-force driven
along $+x$, density-matched phases, neutral wetting ($\\theta = 90°$) so that
the two cases differ only in defender rheology.

The Papanastasiou regularization parameter is chosen as
$m = 3 / \\dot{\\gamma}_c$ with $\\dot{\\gamma}_c = U_\\mathrm{ref} / r_c$, and the
regularized viscosity ceiling $\\mu_\\mathrm{max} = \\mu_p + \\tau_y m$ is passed to
compute_dt so the explicit viscous stability limit is respected (same approach
as benchmark 02_bingham_poiseuille).

Reported dimensionless groups:
    Ca = mu_1 U_ref / sigma            (capillary number)
    Bn = tau_y r_c / (mu_p U_ref)      (Bingham number)
    Y  = tau_y r_c / sigma             (yield-capillary number)

Usage:
    mpirun -np <n> python3 run_invasion.py <init.gsd> <tau_y> [t_end] [max_steps] [a_vox]
      tau_y     : defender yield stress [Pa]; 0 = Newtonian baseline
      t_end     : simulated physical time [s]     (default 0.030)
      max_steps : hard cap on time steps          (default 30000)
      a_vox     : lattice constant in dx units, must match create_packing.py
                  (default 16; used only for the r_c estimate)
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

initfile  = sys.argv[1]
tau_y     = float(sys.argv[2])
t_end     = float(sys.argv[3]) if len(sys.argv) > 3 else 0.030   # [s]
max_steps = int(sys.argv[4])   if len(sys.argv) > 4 else 30000
a_vox     = int(sys.argv[5])   if len(sys.argv) > 5 else 16

sim.create_state_from_gsd(filename=initfile)

# ─── Physical parameters (geometry read from the GSD box) ───────────────────
dx     = 2.0e-4                                # must match create_packing.py
L      = sim.state.box.Lx
a      = a_vox * dx
R_frac = 0.45
r_c    = (np.sqrt(2.0)/2.0 - R_frac) * a       # pore-channel radius
rho0   = 1000.0
mu1    = 0.1          # invader (Newtonian)              [Pa s]
mu_p   = 0.1          # defender plastic/base viscosity  [Pa s]
sigma  = 0.001        # interfacial tension              [N/m]
omega  = 90.0         # neutral wetting (isolate rheology)
gx     = 59.0         # body-force drive along +x        [m/s²]
drho   = 0.01
U_ref  = 0.08         # expected pore velocity scale     [m/s]

# Papanastasiou regularization: m chosen so m * gamma_dot_c ~ 3
gdot_c = U_ref / r_c
m_reg  = 3.0 / gdot_c if tau_y > 0 else 0.0
mu_max2 = mu_p + tau_y * m_reg               # regularized viscosity ceiling

Bn = tau_y * r_c / (mu_p * U_ref) if tau_y > 0 else 0.0
Y  = tau_y * r_c / sigma if tau_y > 0 else 0.0
Ca = mu1 * U_ref / sigma

label    = f"tauy{tau_y:g}"
dumpname = initfile.replace("_init.gsd", f"_{label}_run.gsd")
logname  = initfile.replace("_init.gsd", f"_{label}_run.log")
dtname   = initfile.replace("_init.gsd", f"_{label}_dt.txt")

if device.communicator.rank == 0:
    print(f"── Viscoplastic invasion: tau_y={tau_y} Pa  Bn={Bn:.2f}  Y={Y:.2f}  Ca={Ca:.1f} ──")
    print(f"   m_reg={m_reg:.4f} s  mu_max2={mu_max2:.3f} Pa·s  gdot_c={gdot_c:.1f} 1/s")

# ─── Kernel / nlist ─────────────────────────────────────────────────────────
kernel     = "WendlandC4"
rcut       = sph.kernel.Kappa[kernel] * sph.kernel.OptimalH[kernel] * dx
kernel_obj = sph.kernel.Kernels[kernel]()
nlist = hoomd.nsearch.nlist.Cell(buffer=rcut*0.05, rebuild_check_delay=1,
                                 kappa=kernel_obj.Kappa())

# ─── EOS ────────────────────────────────────────────────────────────────────
eos1 = sph.eos.Tait(); eos1.set_params(rho0, 0.01)
eos2 = sph.eos.Tait(); eos2.set_params(rho0, 0.01)

fA = hoomd.filter.Type(["A"])
fB = hoomd.filter.Type(["B"])
fS = hoomd.filter.Type(["S"])

model = sph.sphmodel.TwoPhaseFlow(
    kernel=kernel_obj, eos1=eos1, eos2=eos2, nlist=nlist,
    fluidgroup1_filter=fA, fluidgroup2_filter=fB, solidgroup_filter=fS,
    densitymethod="SUMMATION", colorgradientmethod="DENSITYRATIO")

# Store the regularized viscosity ceiling as mu2 so compute_dt respects the
# Papanastasiou-regularized stability limit. Physically harmless: once
# activateBingham2 is on, phase 2 uses mu_p from the Bingham parameters.
model.mu1     = mu1
model.mu2     = mu_max2
model.sigma12 = sigma
model.omega   = omega
model.gx      = gx
model.damp    = 300
model.artificialviscosity = True
model.alpha   = 0.1
model.beta    = 0.0
model.densitydiffusion = False

# ─── Speed of sound & timestep ──────────────────────────────────────────────
max_sl = sph_helper.set_max_sl(sim, device, model)

c1, cond1, c2, cond2 = model.compute_speedofsound(
    LREF=L, UREF=U_ref, DX=dx, DRHO=drho, H=max_sl,
    MU1=mu1, MU2=mu_max2, RHO01=rho0, RHO02=rho0, SIGMA12=sigma)
# modest 1.5x margin: the gravity-wave condition already yields Ma ~ 0.01
sph_helper.update_min_c0_tpf(device, model, c1, c2,
                             mode="plain", lref=L, uref=U_ref, cfactor=1.5)
dt, dt_cond = model.compute_dt(
    LREF=L, UREF=U_ref, DX=dx, DRHO=drho, H=max_sl,
    MU1=mu1, MU2=mu_max2, RHO01=rho0, RHO02=rho0, SIGMA12=sigma)

steps = min(int(np.ceil(t_end / dt)), max_steps)
if device.communicator.rank == 0:
    print(f"   c1={c1:.3f} ({cond1})  c2={c2:.3f} ({cond2})")
    print(f"   dt={dt:.3e} s ({dt_cond})  steps={steps}  t_sim={steps*dt:.4f} s")
    with open(dtname, "w") as f:
        f.write(f"{dt:.10e}\n")

# ─── Integrator ─────────────────────────────────────────────────────────────
integrator = sph.Integrator(dt=dt)
integrator.methods.append(sph.methods.VelocityVerletBasic(filter=fA, densitymethod="SUMMATION"))
integrator.methods.append(sph.methods.VelocityVerletBasic(filter=fB, densitymethod="SUMMATION"))
integrator.forces.append(model)
sim.operations.integrator = integrator

# Attach so the C++ object exists, then switch phase B to Bingham rheology.
sim.run(0)
if tau_y > 0:
    model._cpp_obj.activateBingham2(mu_p, tau_y, m_reg, 0.0)
    if device.communicator.rank == 0:
        print(f"   Phase B: Bingham activated (mu_p={mu_p}, tau_y={tau_y}, m={m_reg:.4f})")

# ─── Output ─────────────────────────────────────────────────────────────────
# Remove any stale output GSD left by a previous crashed run.
try:
    os.remove(dumpname)
except OSError:
    pass
device.communicator.barrier()

gsd_period = max(1, steps // 20)
gsd_writer = hoomd.write.GSD(filename=dumpname,
                             trigger=hoomd.trigger.Periodic(gsd_period),
                             mode="wb", dynamic=["property", "momentum"])
sim.operations.writers.append(gsd_writer)

logger = hoomd.logging.Logger(categories=["scalar", "string"])
logger.add(sim, quantities=["timestep", "tps", "walltime"])
propsA = sph.compute.SinglePhaseFlowBasicProperties(filter=fA)
sim.operations.computes.append(propsA)
logger.add(propsA, quantities=["e_kin_fluid", "abs_velocity", "mean_density"])
log_file = open(logname, mode="w", newline="\n")
sim.operations.writers.append(hoomd.write.Table(
    output=log_file, trigger=hoomd.trigger.Periodic(50),
    logger=logger, max_header_len=10))

# ─── Run ────────────────────────────────────────────────────────────────────
sim.run(steps, write_at_start=True)
gsd_writer.flush()

if device.communicator.rank == 0:
    with gsd.hoomd.open(dumpname, "r") as traj:
        s0, s1 = traj[0], traj[-1]
    for tag, s in (("first", s0), ("last", s1)):
        A = s.particles.typeid == 0
        xA = s.particles.position[A, 0]
        vA = s.particles.velocity[A]
        print(f"   [{tag}] step={s.configuration.step}  "
              f"front_x={np.percentile(xA, 99)*1e3:+.3f} mm  "
              f"<|v_A|>={np.linalg.norm(vA, axis=1).mean():.4f} m/s")
    print(f"   done -> {dumpname}")

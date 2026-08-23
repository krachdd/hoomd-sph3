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

Pseudo-3D BODY-FORCE-driven invasion (WP2, pressure-controlled complement).

The piston protocol imposes the flux, so true arrest cannot occur there; this
runner imposes the driving (a uniform body force g_x on BOTH fluids, the
pressure-gradient analogue) on the piston-free reservoir domain from
create_piston_domain.py --no-piston. Purpose: bracket the arrest boundary in
the low-Ca / high-Y corner, where the observable is the invader flux
(mean phase-A velocity): flowing runs settle at a finite plateau, arrested
runs decay toward the regularization creep floor.

Arrest estimate for a slot throat of half-width r_t:
    g_crit ~ tau_y / (rho * r_t)

Usage:
    python3 run_bodyforce_invasion.py <init.gsd> <tau_y> --g GX [--sigma 0.01]
                                      [--steps N] [--res 15] [--vcap 1.0]
                                      [--damp 1000]
"""

import sys, os, argparse
import numpy as np

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules'))

import hoomd
from hoomd import sph
import gsd.hoomd
import sph_helper

ap = argparse.ArgumentParser()
ap.add_argument("init")
ap.add_argument("tau_y", type=float)
ap.add_argument("--g", type=float, required=True, help="body force along +x [m/s^2]")
ap.add_argument("--sigma", type=float, default=0.01)
ap.add_argument("--steps", type=int, default=30000)
ap.add_argument("--res", type=int, default=15)
ap.add_argument("--vcap", type=float, default=1.0,
                help="velocity limiter cap [m/s] (no U_p scale here)")
ap.add_argument("--damp", type=int, default=1000,
                help="body-force ramp-up steps (model.damp)")
ap.add_argument("--m_gd", type=float, default=3.0,
                help="regularization sharpness m*gdot_c; dt_Fourier ~ 1/(mu_p + tau_y*m)")
args = ap.parse_args()

device = hoomd.device.CPU(notice_level=1)

# Fractured-MPI-world guard: if hoomd is not linked against the MPI that
# launched this job, MPI_Init falls back to a singleton world and EVERY
# process believes it is rank 0 of 1 -- they all then race on the same
# output files (observed on the cluster as "GSD: Not a GSD file" while the
# file exists). Abort immediately with a clear message instead.
_launched = max(int(os.environ.get("OMPI_COMM_WORLD_SIZE", "1")),
                int(os.environ.get("PMI_SIZE", "1")),
                int(os.environ.get("SLURM_STEP_NUM_TASKS", "1")))
if _launched > 1 and device.communicator.num_ranks == 1:
    raise RuntimeError(
        f"MPI world fractured: the launcher started {_launched} processes but "
        "hoomd sees num_ranks == 1. hoomd is not built against the MPI that "
        "launched this job (module mismatch). Rebuild hoomd with the loaded "
        "MPI module, or launch with the matching mpirun/srun.")

sim = hoomd.Simulation(device=device)

initfile, tau_y, gx = args.init, args.tau_y, args.g
sim.create_state_from_gsd(filename=initfile)

# ─── Physical parameters (must match create_piston_domain.py) ───────────────
res     = args.res
dx      = 2.0e-4 * 8.0 / res
r_t     = (res / 2.0) * dx  # half of the res-dx throat = 0.8 mm at any res
rho0    = 1000.0
mu1     = 0.1               # invader
mu_p    = 0.1               # defender plastic/base viscosity
sigma   = args.sigma
omega   = 90.0
drho    = 0.01
U_ref   = 0.1               # velocity headroom for c/dt (flow scale unknown a priori)

gdot_c  = U_ref / r_t
m_reg   = args.m_gd / gdot_c if tau_y > 0 else 0.0
mu_max2 = mu_p + tau_y * m_reg

g_crit = tau_y / (rho0 * r_t) if tau_y > 0 else 0.0
Y      = tau_y * r_t / sigma

label    = f"tauy{tau_y:g}_g{gx:g}_sig{sigma:g}_m{args.m_gd:g}"
dumpname = initfile.replace("_init.gsd", f"_{label}_run.gsd")
logname  = initfile.replace("_init.gsd", f"_{label}_run.log")

if device.communicator.rank == 0:
    print(f"── Body-force invasion: tau_y={tau_y} Pa  g={gx} m/s²  "
          f"Y={Y:.1f}  g/g_crit={gx/g_crit if g_crit else float('inf'):.2f} ──")

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
fS = hoomd.filter.Type(["S", "P"])   # P declared but empty in this domain

model = sph.sphmodel.TwoPhaseFlow(
    kernel=kernel_obj, eos1=eos1, eos2=eos2, nlist=nlist,
    fluidgroup1_filter=fA, fluidgroup2_filter=fB, solidgroup_filter=fS,
    densitymethod="SUMMATION", colorgradientmethod="DENSITYRATIO")

model.mu1     = mu1
model.mu2     = mu_max2     # regularized ceiling for compute_dt
model.sigma12 = sigma
model.omega   = omega
model.artificialviscosity = True
model.alpha   = 0.1
model.beta    = 0.0
model.densitydiffusion = False
model.gx      = gx
model.damp    = args.damp   # gentle body-force ramp-up

max_sl = sph_helper.set_max_sl(sim, device, model)
c1, cond1, c2, cond2 = model.compute_speedofsound(
    LREF=sim.state.box.Lx, UREF=U_ref, DX=dx, DRHO=drho, H=max_sl,
    MU1=mu1, MU2=mu_max2, RHO01=rho0, RHO02=rho0, SIGMA12=sigma)

# Pressure-based sound-speed calibration: the body force integrates to a
# pressure difference rho*g*Lx around the periodic loop; require
# rho0 c^2 >= p_drive / drho so neither phase compresses under the drive.
p_drive = rho0 * gx * sim.state.box.Lx
c_press = np.sqrt(p_drive / (rho0 * drho))
c_target = max(1.5 * c1, c_press)
model.set_speedofsound(c_target, c_target)
if device.communicator.rank == 0:
    print(f"   p_drive = rho*g*Lx = {p_drive:.0f} Pa -> c = {c_target:.2f} m/s")
dt, dt_cond = model.compute_dt(
    LREF=sim.state.box.Lx, UREF=U_ref, DX=dx, DRHO=drho, H=max_sl,
    MU1=mu1, MU2=mu_max2, RHO01=rho0, RHO02=rho0, SIGMA12=sigma)

steps = args.steps
if device.communicator.rank == 0:
    print(f"   dt={dt:.3e} s ({dt_cond})  steps={steps}  t_sim={steps*dt:.4f} s")

integrator = sph.Integrator(dt=dt)
integrator.methods.append(sph.methods.VelocityVerletBasic(filter=fA, densitymethod="SUMMATION"))
integrator.methods.append(sph.methods.VelocityVerletBasic(filter=fB, densitymethod="SUMMATION"))
integrator.forces.append(model)
sim.operations.integrator = integrator


class VelocityLimiter(hoomd.custom.Action):
    """Cap fluid particle speeds at v_cap (see run_piston_invasion.py)."""

    def __init__(self, v_cap, fluid_typeids=(0, 1)):
        super().__init__()
        self._v_cap = float(v_cap)
        self._fluid_typeids = np.asarray(fluid_typeids)

    def act(self, timestep):
        with self._state.cpu_local_snapshot as snap:
            v = snap.particles.velocity
            sp = np.linalg.norm(v, axis=1)
            mask = (sp > self._v_cap) & np.isin(snap.particles.typeid,
                                                self._fluid_typeids)
            if np.any(mask):
                v[mask] *= (self._v_cap / sp[mask])[:, None]


sim.operations.updaters.append(
    hoomd.update.CustomUpdater(trigger=hoomd.trigger.Periodic(1),
                               action=VelocityLimiter(args.vcap)))

# switch defender to Bingham after attach; same stabilizer stack as the
# piston runner (delta+ shifting with interface projection, Riemann pair
# dissipation) — see run_piston_invasion.py for the full rationale.
sim.run(0)
if tau_y > 0:
    model._cpp_obj.activateBingham2(mu_p, tau_y, m_reg, 0.0)
model._cpp_obj.activateParticleShifting(0.1, 0.2, 4, True)
model._cpp_obj.activateRiemannDissipation(1.0)

# ─── Output ─────────────────────────────────────────────────────────────────
try:
    os.remove(dumpname)
except OSError:
    pass
device.communicator.barrier()

gsd_period = max(1, steps // 30)
gsd_writer = hoomd.write.GSD(filename=dumpname,
                             trigger=hoomd.trigger.Periodic(gsd_period),
                             mode="wb", dynamic=["property", "momentum"])
sim.operations.writers.append(gsd_writer)

logger = hoomd.logging.Logger(categories=["scalar", "string"])
logger.add(sim, quantities=["timestep", "tps", "walltime"])
propsA = sph.compute.SinglePhaseFlowBasicProperties(filter=fA)
propsB = sph.compute.SinglePhaseFlowBasicProperties(filter=fB)
sim.operations.computes.append(propsA)
sim.operations.computes.append(propsB)
logger.add(propsA, quantities=["abs_velocity", "mean_density"])
logger.add(propsB, quantities=["abs_velocity"])
log_file = open(logname, mode="w", newline="\n")
sim.operations.writers.append(hoomd.write.Table(
    output=log_file, trigger=hoomd.trigger.Periodic(max(1, steps // 60)),
    logger=logger, max_header_len=10))

# ─── Run ────────────────────────────────────────────────────────────────────
snap0 = sim.state.get_snapshot()
front0 = None
if snap0.communicator.rank == 0:
    A = snap0.particles.typeid == 0
    front0 = float(np.percentile(snap0.particles.position[A, 0], 99))

sim.run(steps)
gsd_writer.flush()

snap = sim.state.get_snapshot()
if snap.communicator.rank == 0:
    A = snap.particles.typeid == 0
    front = float(np.percentile(snap.particles.position[A, 0], 99))
    uA = propsA.abs_velocity
    # classifier hint (final verdict belongs to postprocessing on the log):
    # flowing if the phase-A flux is well above the regularization creep floor
    u_creep = rho0 * gx * r_t**2 / (8.0 * mu_max2)
    verdict = "FLOWING" if uA > 2.0 * u_creep else "ARRESTED(creep-level)"
    print(f"   invader front x = {front*1e3:+.2f} mm (moved "
          f"{(front-front0)*1e3:+.3f} mm)")
    print(f"   mean |u_A| = {uA:.4e} m/s  creep floor = {u_creep:.4e} m/s "
          f"-> {verdict}")
    print(f"   done -> {dumpname}")

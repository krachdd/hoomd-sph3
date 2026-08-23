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

Pseudo-3D piston-driven invasion (WP2 variant) — run script.

A solid piston (type P, prescribed velocity U_p) pushes the Newtonian invader
(A) from its reservoir through the cylinder matrix saturated with the defender
(B; Newtonian for tau_y = 0, Bingham otherwise). Velocity-controlled
displacement: Ca is imposed exactly; the yield stress appears as the piston
force, logged via ComputeSolidProperties on type P.

The piston moves kinematically: its GSD velocity is U_p x^ (so the Adami
no-slip BC sees a moving wall) and a CustomUpdater advances its positions by
U_p*dt each step. The run stops when the piston has traveled 80% of reservoir
A (or at max_steps).

NOTE: the kinematic mover writes through state.cpu_local_snapshot and is
validated single-rank; for MPI runs verify piston-particle migration first.

Usage:
    python3 run_piston_invasion.py <init.gsd> <tau_y> [--U_p 0.02] [--sigma 0.001]
                                   [--steps N] [--res 8] [--resA_mm 12.0]
                                   [--ramp 2000]
      --steps   : hard step budget (the run still stops earlier when the piston
                  has traveled 0.8 * resA_mm)
      --ramp    : soft-start ramp length in steps; U_p rises linearly from 0 so
                  the yield peak in the force trace is not buried in the
                  water-hammer transient
"""

import sys, os, argparse
import numpy as np

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules'))

import hoomd
from hoomd import sph
import gsd.hoomd
import sph_helper

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

ap = argparse.ArgumentParser()
ap.add_argument("init")
ap.add_argument("tau_y", type=float)
ap.add_argument("--U_p", type=float, default=0.02)
ap.add_argument("--sigma", type=float, default=0.001)
ap.add_argument("--steps", type=int, default=30000)
ap.add_argument("--res", type=int, default=8)
ap.add_argument("--resA_mm", type=float, default=12.0)
ap.add_argument("--ramp", type=int, default=2000)
ap.add_argument("--m_gd", type=float, default=3.0,
                help="regularization sharpness m*gdot_c; dt_Fourier ~ 1/(mu_p + tau_y*m)")
args = ap.parse_args()

initfile  = args.init
tau_y     = args.tau_y
U_p       = args.U_p
max_steps = args.steps
res       = args.res

sim.create_state_from_gsd(filename=initfile)

# ─── Physical parameters (must match create_piston_domain.py) ───────────────
dx      = 2.0e-4 * 8.0 / res
r_t     = (res / 2.0) * dx  # half of the res-dx throat = 0.8 mm at any res
rho0    = 1000.0
mu1     = 0.1               # invader
mu_p    = 0.1               # defender plastic/base viscosity
sigma   = args.sigma
omega   = 90.0
drho    = 0.01
U_ref   = 4.0 * U_p         # headroom for pore-scale acceleration

gdot_c  = U_ref / r_t
m_reg   = args.m_gd / gdot_c if tau_y > 0 else 0.0
mu_max2 = mu_p + tau_y * m_reg

Ca = mu1 * U_p / sigma
Bn = tau_y * r_t / (mu_p * U_p) if tau_y > 0 else 0.0

label    = f"tauy{tau_y:g}_Up{U_p:g}_sig{sigma:g}_m{args.m_gd:g}"
dumpname = initfile.replace("_init.gsd", f"_{label}_run.gsd")
logname  = initfile.replace("_init.gsd", f"_{label}_run.log")

if device.communicator.rank == 0:
    print(f"── Piston invasion: tau_y={tau_y} Pa  U_p={U_p} m/s  "
          f"Ca={Ca:.2f}  Bn={Bn:.2f} ──")

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
fS = hoomd.filter.Type(["S", "P"])   # matrix grains + piston are boundaries
fP = hoomd.filter.Type(["P"])

model = sph.sphmodel.TwoPhaseFlow(
    kernel=kernel_obj, eos1=eos1, eos2=eos2, nlist=nlist,
    fluidgroup1_filter=fA, fluidgroup2_filter=fB, solidgroup_filter=fS,
    densitymethod="SUMMATION", colorgradientmethod="DENSITYRATIO")

model.mu1     = mu1
model.mu2     = mu_max2     # regularized ceiling for compute_dt (cf. run_invasion.py)
model.sigma12 = sigma
model.omega   = omega
model.artificialviscosity = True
model.alpha   = 0.1
model.beta    = 0.0
model.densitydiffusion = False

max_sl = sph_helper.set_max_sl(sim, device, model)
c1, cond1, c2, cond2 = model.compute_speedofsound(
    LREF=sim.state.box.Lx, UREF=U_ref, DX=dx, DRHO=drho, H=max_sl,
    MU1=mu1, MU2=mu_max2, RHO01=rho0, RHO02=rho0, SIGMA12=sigma)

# Pressure-based sound-speed calibration: the piston must transmit the full
# Bingham driving pressure across the matrix without compressing the invader
# reservoir. Estimate p_drive as yield threshold + yielded viscous Poiseuille
# contribution across the matrix and require rho0 c^2 >= p_drive / drho.
L_mat   = 96 * res / 8.0 * dx                 # matrix length
phi     = 0.65                                # matrix porosity
p_yield = 2.0 * tau_y * L_mat / r_t
p_visc  = 8.0 * mu_p * (U_p / phi) * L_mat / r_t**2
p_drive = p_yield + p_visc
c_press = np.sqrt(p_drive / (rho0 * drho))
c_target = max(1.5 * c1, c_press)
model.set_speedofsound(c_target, c_target)
if device.communicator.rank == 0:
    print(f"   p_drive estimate = {p_drive:.0f} Pa (yield {p_yield:.0f} + "
          f"viscous {p_visc:.0f}) -> c = {c_target:.2f} m/s")
dt, dt_cond = model.compute_dt(
    LREF=sim.state.box.Lx, UREF=U_ref, DX=dx, DRHO=drho, H=max_sl,
    MU1=mu1, MU2=mu_max2, RHO01=rho0, RHO02=rho0, SIGMA12=sigma)

# stop when the piston has traveled 80% of reservoir A
travel_max = 0.8 * args.resA_mm * 1e-3
steps = min(int(np.ceil(travel_max / (U_p * dt))), max_steps)
if device.communicator.rank == 0:
    print(f"   dt={dt:.3e} s ({dt_cond})  steps={steps}  "
          f"piston travel={steps*dt*U_p*1e3:.2f} mm  t_sim={steps*dt:.4f} s")

integrator = sph.Integrator(dt=dt)
integrator.methods.append(sph.methods.VelocityVerletBasic(filter=fA, densitymethod="SUMMATION"))
integrator.methods.append(sph.methods.VelocityVerletBasic(filter=fB, densitymethod="SUMMATION"))
integrator.forces.append(model)
sim.operations.integrator = integrator


class PistonDrive(hoomd.custom.Action):
    """Advance piston (type P) positions kinematically by u(t)*dt per step.

    Soft start: u(t) = U_p * min(1, timestep/ramp_steps). During the ramp both
    the kinematic displacement AND the piston particles' velocity field are
    updated (the Adami no-slip BC reads the velocity), so the fluid is
    accelerated gently and the yield peak in the force trace is not buried in
    the acoustic water-hammer transient of an impulsive start.
    """

    def __init__(self, u_p, dt, ramp_steps=0):
        super().__init__()
        self._u_p = u_p
        self._dt = dt
        self._ramp = max(int(ramp_steps), 0)

    def act(self, timestep):
        f = min(1.0, timestep / self._ramp) if self._ramp > 0 else 1.0
        u = self._u_p * f
        with self._state.cpu_local_snapshot as snap:
            piston = snap.particles.typeid == 3
            snap.particles.position[piston, 0] += u * self._dt
            # always (re)write the velocity: makes the run's U_p authoritative
            # over whatever U_p the init GSD was generated with
            snap.particles.velocity[piston, 0] = u


piston_drive = PistonDrive(U_p, dt, ramp_steps=args.ramp)
sim.operations.updaters.append(
    hoomd.update.CustomUpdater(trigger=hoomd.trigger.Periodic(1),
                               action=piston_drive))


class VelocityLimiter(hoomd.custom.Action):
    """Cap fluid particle speeds at v_cap.

    When the invasion front pinches a thin defender film against a grain, the
    trapped particle is ejected tangentially ("watermelon-seed" squirt). The
    squirts are transient and mostly self-resolve, but under a fixed dt one
    occasionally feeds back into a runaway (observed deterministically at
    ~17.4k steps, first grain column, mirrored at +-y) and leaves the box.
    Healthy flow tops out at ~0.2-0.5 m/s here; the cap only ever touches
    numerically pathological particles.
    """

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


v_cap = 25.0 * U_p
sim.operations.updaters.append(
    hoomd.update.CustomUpdater(trigger=hoomd.trigger.Periodic(1),
                               action=VelocityLimiter(v_cap)))

# switch defender to Bingham after attach
sim.run(0)
if tau_y > 0:
    model._cpp_obj.activateBingham2(mu_p, tau_y, m_reg, 0.0)
model.computeSolidForces()          # needed for the piston force observable
# delta+-SPH particle shifting (Sun et al. 2017): prevents the pairing/wall-
# penetration instability observed when the yield-resisting defender is
# squeezed against grain surfaces in the throats. Interface-normal shift is
# projected out; modest amplitude following the implementation guidance.
model._cpp_obj.activateParticleShifting(0.1, 0.2, 4, True)
# Riemann-type pair dissipation (Zhang, Hu & Adams 2017): damps the relative
# approach velocity of particle pairs. Both res-8 attempts (with and without
# shifting) died at the SAME step ~17.4k when the invasion front reaches the
# third grain column: a local pressure spike ejects a particle through a
# grain. Riemann dissipation targets exactly that acoustic spike without the
# bulk-viscosity pollution Monaghan AV would add to the Bingham rheology.
model._cpp_obj.activateRiemannDissipation(1.0)

# ─── Output ─────────────────────────────────────────────────────────────────
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
piston_props = sph.compute.SolidProperties(filter=fP)
sim.operations.computes.append(piston_props)
logger.add(piston_props, quantities=["total_drag_x"])
propsA = sph.compute.SinglePhaseFlowBasicProperties(filter=fA)
sim.operations.computes.append(propsA)
logger.add(propsA, quantities=["abs_velocity", "mean_density"])
log_file = open(logname, mode="w", newline="\n")
sim.operations.writers.append(hoomd.write.Table(
    output=log_file, trigger=hoomd.trigger.Periodic(50),
    logger=logger, max_header_len=10))

# ─── Run ────────────────────────────────────────────────────────────────────
x_piston_start = None
snap0 = sim.state.get_snapshot()
if snap0.communicator.rank == 0:
    P = snap0.particles.typeid == 3
    x_piston_start = float(snap0.particles.position[P, 0].mean())

sim.run(steps)
gsd_writer.flush()

snap = sim.state.get_snapshot()
if snap.communicator.rank == 0:
    P = snap.particles.typeid == 3
    A = snap.particles.typeid == 0
    x_piston = float(snap.particles.position[P, 0].mean())
    travel = x_piston - x_piston_start
    front = float(np.percentile(snap.particles.position[A, 0], 99))
    print(f"   piston travel = {travel*1e3:.3f} mm "
          f"(kinematic prediction {steps*dt*U_p*1e3:.3f} mm)")
    print(f"   invader front x = {front*1e3:+.2f} mm  "
          f"piston force x = {piston_props.total_drag_x:.4e} N")
    print(f"   done -> {dumpname}")

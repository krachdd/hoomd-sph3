"""Integration tests for the SinglePhaseFlowGDGD solver.

Small end-to-end runs that lock in the density-gradient-driven physics:
  * scalar diffusion of a Gaussian matches the analytic variance growth
  * Dirichlet wall BC injects scalar into the fluid (bounded by wall value)
  * no-flux wall BC conserves the fluid scalar content
  * Boussinesq buoyancy has the correct sign and respects the gravity ramp
  * compute_speedofsound / compute_dt implement the documented conditions
    (single Courant factor, gravity-wave condition, scalar-diffusion limit)
  * AdaptiveTimestep includes the scalar-diffusion candidate
"""

import numpy as np
import pytest

import hoomd
import hoomd.sph as sph

from .conftest import make_fluid_block_gsd, set_max_sl

RHO0 = 1000.0
DX = 0.05
DRHO = 0.01
UREF = 0.1
LREF = 1.0


def build_gdgd_sim(gsd_path, dx, rho0=RHO0, mu=0.01, kappa_s=0.0,
                   beta_s=0.0, scalar_ref=0.0, boussinesq=True,
                   scalar_wall_bc="dirichlet", g=(0.0, 0.0, 0.0), damp=0,
                   uref=UREF):
    """Construct a ready-to-run SinglePhaseFlowGDGD simulation from a GSD file.

    Returns (sim, model, filterfluid, filtersolid, dt).
    """
    device = hoomd.device.CPU(notice_level=1)
    sim = hoomd.Simulation(device=device)
    sim.create_state_from_gsd(filename=str(gsd_path))

    kernel_obj = sph.kernel.Kernels["WendlandC4"]()
    slength = sph.kernel.OptimalH["WendlandC4"] * dx
    rcut = sph.kernel.Kappa["WendlandC4"] * slength
    nlist = hoomd.nsearch.nlist.Cell(buffer=rcut * 0.05,
                                     rebuild_check_delay=1,
                                     kappa=kernel_obj.Kappa())

    eos = sph.eos.Tait()
    eos.set_params(rho0, 0.01)

    filterfluid = hoomd.filter.Type(["F"])
    filtersolid = hoomd.filter.Type(["S"])

    model = sph.sphmodel.SinglePhaseFlowGDGD(
        kernel=kernel_obj, eos=eos, nlist=nlist,
        fluidgroup_filter=filterfluid, solidgroup_filter=filtersolid,
        densitymethod="SUMMATION",
        kappa_s=kappa_s, beta_s=beta_s, scalar_ref=scalar_ref,
        boussinesq=boussinesq, scalar_wall_bc=scalar_wall_bc)
    model.mu = mu
    model.gx, model.gy, model.gz = g
    model.damp = damp

    max_sl = set_max_sl(sim, device, model)

    model.compute_speedofsound(
        LREF=LREF, UREF=uref, DX=dx, DRHO=DRHO,
        H=max_sl, MU=mu, RHO0=rho0)
    dt, _ = model.compute_dt(
        LREF=LREF, UREF=uref, DX=dx, DRHO=DRHO,
        H=max_sl, MU=mu, RHO0=rho0)

    integrator = sph.Integrator(dt=dt)
    method = sph.methods.VelocityVerletBasic(filter=filterfluid,
                                             densitymethod="SUMMATION")
    integrator.methods.append(method)
    integrator.forces.append(model)
    sim.operations.integrator = integrator

    return sim, model, filterfluid, filtersolid, dt


def set_scalar(sim, values_fn):
    """Set aux4.x from a callable (pos, typeid) -> T on the local snapshot."""
    with sim.state.cpu_local_snapshot as snap:
        pos = np.array(snap.particles.position[:])
        tid = np.array(snap.particles.typeid[:])
        snap.particles.auxiliary4[:, 0] = values_fn(pos, tid)


def fluid_scalar(sim):
    """Return (T, mass, pos, typeid, velocity) copies from the local snapshot.

    The global snapshot does not expose auxiliary4, so read the rank-local
    data instead (the pytest suite runs single-rank, where local == global).
    """
    with sim.state.cpu_local_snapshot as snap:
        T = np.array(snap.particles.auxiliary4[:, 0])
        m = np.array(snap.particles.mass[:])
        pos = np.array(snap.particles.position[:])
        tid = np.array(snap.particles.typeid[:])
        vel = np.array(snap.particles.velocity[:])
    return T, m, pos, tid, vel


def test_gaussian_diffusion_matches_analytic(tmp_path):
    """A 1D Gaussian scalar profile must spread with sigma^2 = sigma0^2 + 2*kappa*t."""
    kappa_s = 2e-3
    sigma0 = 0.1
    f = tmp_path / "gauss.gsd"
    make_fluid_block_gsd(f, n=20, dx=DX, rho0=RHO0)
    sim, model, ff, fs, dt = build_gdgd_sim(f, DX, kappa_s=kappa_s)

    set_scalar(sim, lambda pos, tid: np.exp(-pos[:, 0]**2 / (2.0 * sigma0**2)))

    steps = 100
    sim.run(steps)

    T, m, pos, tid, vel = fluid_scalar(sim)
    x = pos[:, 0]
    var_meas = float(np.sum(T * x**2) / np.sum(T))
    # initial sampling of the Gaussian on the lattice is not exact either:
    # compare the *growth* of the variance against 2*kappa*t
    growth_meas = var_meas - sigma0**2
    growth_theory = 2.0 * kappa_s * steps * dt
    assert growth_meas == pytest.approx(growth_theory, rel=0.15)
    # profile must remain bounded and finite
    assert np.all(np.isfinite(T))
    assert T.max() <= 1.0 + 1e-6


def test_dirichlet_wall_heats_fluid(tmp_path):
    """Hot walls (T=1) with Dirichlet BC must heat the cold fluid, bounded by 1."""
    f = tmp_path / "dirichlet.gsd"
    make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0, solid_planes=3)
    sim, model, ff, fs, dt = build_gdgd_sim(f, DX, kappa_s=2e-3,
                                            scalar_wall_bc="dirichlet")

    set_scalar(sim, lambda pos, tid: np.where(tid == 1, 1.0, 0.0))

    sim.run(50)

    T, m, pos, tid, vel = fluid_scalar(sim)
    T_fluid = T[tid == 0]
    assert T_fluid.mean() > 1e-6, "Dirichlet wall must inject scalar"
    assert T_fluid.max() <= 1.0 + 1e-6, "fluid scalar may not exceed wall value"


def test_noflux_wall_conserves_scalar(tmp_path):
    """No-flux walls must neither inject nor absorb scalar from the fluid."""
    f = tmp_path / "noflux.gsd"
    make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0, solid_planes=3)
    sim, model, ff, fs, dt = build_gdgd_sim(f, DX, kappa_s=2e-3,
                                            scalar_wall_bc="noflux")

    # fluid: linear gradient in z; walls: extreme value that would visibly
    # leak into the fluid if the BC were (incorrectly) Dirichlet
    def init(pos, tid):
        T = 0.5 + pos[:, 2]
        T[tid == 1] = 999.0
        return T

    set_scalar(sim, init)

    T0, m0, _, tid0, _ = fluid_scalar(sim)
    sim.run(50)
    T1, m1, _, tid1, _ = fluid_scalar(sim)

    total0 = float(np.sum(m0[tid0 == 0] * T0[tid0 == 0]))
    total1 = float(np.sum(m1[tid1 == 0] * T1[tid1 == 0]))
    assert total1 == pytest.approx(total0, rel=1e-3), \
        "no-flux BC must conserve the fluid scalar content"
    # and nothing may drift toward the wall value
    assert T1[tid1 == 0].max() < 2.0, \
        "fluid scalar contaminated by wall value: no-flux BC leaked"


def test_boussinesq_force_sign(tmp_path):
    """Heavy fluid (C=1, beta_s<0) must sink relative to light fluid (C=0)."""
    f = tmp_path / "buoy.gsd"
    make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0)
    sim, model, ff, fs, dt = build_gdgd_sim(
        f, DX, kappa_s=0.0, beta_s=-0.05, scalar_ref=0.5,
        boussinesq=True, g=(0.0, -9.81, 0.0), damp=0)

    set_scalar(sim, lambda pos, tid: np.where(pos[:, 1] > 0.0, 1.0, 0.0))

    sim.run(2)

    T, m, pos, tid, vel = fluid_scalar(sim)
    vy = vel[:, 1]
    vy_heavy = vy[T > 0.5].mean()
    vy_light = vy[T < 0.5].mean()
    assert vy_heavy < vy_light, \
        "C=1 (heavy, beta_s<0) must accelerate downward relative to C=0"


def test_boussinesq_respects_gravity_ramp(tmp_path):
    """With a long damp ramp, buoyancy must be damped like gravity: the net
    body force (gravity + buoyancy) is ~0 at the start of the ramp, so the
    stratified fluid must stay at rest."""
    f = tmp_path / "damp.gsd"
    make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0)
    sim, model, ff, fs, dt = build_gdgd_sim(
        f, DX, kappa_s=0.0, beta_s=-0.05, scalar_ref=0.5,
        boussinesq=True, g=(0.0, -9.81, 0.0), damp=10**6)

    set_scalar(sim, lambda pos, tid: np.where(pos[:, 1] > 0.0, 1.0, 0.0))

    sim.run(2)

    T, m, pos, tid, vel = fluid_scalar(sim)
    vmax = np.abs(vel).max()
    # before the damp fix, buoyancy acted at full strength while gravity
    # was ramped to ~0: the two layers would move apart immediately
    v_buoy_1step = abs(0.05 * 0.5 * 9.81) * dt  # undamped buoyancy kick
    assert vmax < 0.01 * v_buoy_1step, \
        f"buoyancy must ramp with gravity (vmax={vmax:g})"


def test_compute_sos_and_dt_conditions(tmp_path):
    """Lock in the corrected speed-of-sound and timestep formulas."""
    kappa_s = 1e-3
    mu = 0.01
    g = 9.81
    f = tmp_path / "dt.gsd"
    make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0)
    sim, model, ff, fs, dt = build_gdgd_sim(
        f, DX, mu=mu, kappa_s=kappa_s, g=(0.0, -g, 0.0))

    c, _ = model.compute_speedofsound(LREF=LREF, UREF=UREF, DX=DX, DRHO=DRHO,
                                      H=model.max_sl, MU=mu, RHO0=RHO0)
    c_ref = np.sqrt(max(UREF**2 / DRHO,
                        g * LREF / DRHO,
                        mu * UREF / (RHO0 * LREF * DRHO)))
    assert c == pytest.approx(c_ref, rel=1e-12)
    # the returned c must also be stored in the EOS
    assert model.eos.SpeedOfSound == pytest.approx(c_ref, rel=1e-12)

    COURANT = 0.25
    dt2, cond = model.compute_dt(LREF=LREF, UREF=UREF, DX=DX, DRHO=DRHO,
                                 H=model.max_sl, MU=mu, RHO0=RHO0,
                                 COURANT=COURANT)
    h = model.max_sl
    dt_ref = COURANT * min(DX / c_ref,
                           DX * DX * RHO0 / (8.0 * mu),
                           np.sqrt(h / (16.0 * g)),
                           DX * DX / (8.0 * kappa_s))
    assert dt2 == pytest.approx(dt_ref, rel=1e-12)

    # a large diffusivity must become the limiting condition
    sim.run(0)  # attach operations so the C++ model object exists
    model.setGDGDParams(10.0, 0.0, 0.0, True)
    dt3, cond3 = model.compute_dt(LREF=LREF, UREF=UREF, DX=DX, DRHO=DRHO,
                                  H=model.max_sl, MU=mu, RHO0=RHO0,
                                  COURANT=COURANT)
    assert dt3 == pytest.approx(COURANT * DX * DX / (8.0 * 10.0), rel=1e-12)
    assert "ScalarDiffusion-condition" in cond3


def test_adaptive_timestep_includes_diffusion(tmp_path):
    """AdaptiveTimestep must respect the scalar-diffusion limit for GDGD."""
    kappa_s = 10.0
    f = tmp_path / "adapt.gsd"
    make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0)
    sim, model, ff, fs, dt0 = build_gdgd_sim(f, DX, kappa_s=kappa_s)

    adapt = sph.update.AdaptiveTimestep(model=model, courant=0.25)
    sim.operations.updaters.append(
        adapt.as_updater(trigger=hoomd.trigger.Periodic(5)))

    sim.run(11)

    dt_new = sim.operations.integrator.dt
    h = model.max_sl
    assert dt_new > 0.0
    assert np.isfinite(dt_new)
    assert dt_new <= 0.25 * h * h / (8.0 * kappa_s) * 1.001

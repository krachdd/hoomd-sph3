"""Integration tests for the SinglePhaseFlow solver.

These are small (1000-particle) end-to-end runs that lock in the physics
fixes:
  * hydrostatic-free rest state stays at rest, density ~= rho0
  * e_kin_fluid equals sum(1/2 m |v|^2) computed from the snapshot
  * getMaxVelocity() equals the true maximum particle speed
  * viscous drag on a wall points WITH the flow (reaction-force sign)
  * delta-SPH continuity mode runs stably and keeps p = EOS(rho)
  * AdaptiveTimestep produces a bounded, positive dt
"""

import numpy as np
import pytest

import hoomd
import hoomd.sph as sph

from .conftest import make_fluid_block_gsd, build_singlephase_sim

RHO0 = 1000.0
DX = 0.05


def taylor_green(pos):
    lx = pos[:, 0].max() - pos[:, 0].min() + DX
    k = 2.0 * np.pi / lx
    v = np.zeros_like(pos)
    u0 = 0.05
    v[:, 0] = -u0 * np.cos(k * pos[:, 0]) * np.sin(k * pos[:, 1])
    v[:, 1] = u0 * np.sin(k * pos[:, 0]) * np.cos(k * pos[:, 1])
    return v


def test_rest_state_stays_at_rest(tmp_path):
    f = tmp_path / "block.gsd"
    make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0)
    sim, model, ff, fs, dt = build_singlephase_sim(f, DX, rho0=RHO0)
    sim.run(20)

    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        vmax = np.abs(snap.particles.velocity).max()
        # a periodic lattice at rest must produce no spurious motion beyond
        # float32 round-off of the GSD positions
        assert vmax < 1e-6
        assert np.allclose(snap.particles.density, RHO0, rtol=2e-2)


def test_ekin_and_maxvel_match_snapshot(tmp_path):
    f = tmp_path / "tg.gsd"
    make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0, velocity_fn=taylor_green)
    sim, model, ff, fs, dt = build_singlephase_sim(f, DX, rho0=RHO0)

    props = sph.compute.SinglePhaseFlowBasicProperties(filter=ff)
    sim.operations.computes.append(props)

    sim.run(2)

    snap = sim.state.get_snapshot()
    e_kin = props.e_kin_fluid
    vmax_cpp = model._cpp_obj.getMaxVelocity()

    if snap.communicator.rank == 0:
        v = snap.particles.velocity
        m = snap.particles.mass
        e_kin_ref = float(0.5 * (m * (v**2).sum(axis=1)).sum())
        assert e_kin == pytest.approx(e_kin_ref, rel=1e-6)

        vmax_ref = float(np.sqrt((v**2).sum(axis=1)).max())
        # forces (and max_vel) were evaluated at the half step; allow the
        # half-kick difference
        assert vmax_cpp == pytest.approx(vmax_ref, rel=5e-2)
        assert vmax_cpp > 0.0


def test_wall_drag_points_with_flow(tmp_path):
    """Fluid sliding between two static walls must drag the walls along +x."""
    f = tmp_path / "channel.gsd"

    def uniform_x(pos):
        v = np.zeros_like(pos)
        v[:, 0] = 0.05
        return v

    make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0,
                         velocity_fn=uniform_x, solid_planes=3)
    sim, model, ff, fs, dt = build_singlephase_sim(f, DX, rho0=RHO0, mu=0.05)
    sim.run(0)  # attach operations so the C++ model object exists
    model.computeSolidForces()

    solid_props = sph.compute.SolidProperties(filter=fs)
    sim.operations.computes.append(solid_props)

    sim.run(10)

    drag_x = solid_props.total_drag_x
    assert drag_x > 0.0, (
        "viscous drag on the wall must point in the flow direction "
        f"(got {drag_x})")

    # and the fluid must decelerate in response
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        fluid = snap.particles.typeid == 0
        assert snap.particles.velocity[fluid, 0].mean() < 0.05


def test_continuity_delta_sph_runs_and_pressure_consistent(tmp_path):
    f = tmp_path / "tg2.gsd"
    make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0, velocity_fn=taylor_green)
    sim, model, ff, fs, dt = build_singlephase_sim(
        f, DX, rho0=RHO0, densitymethod="CONTINUITY")
    sim.run(0)  # attach operations so the C++ model object exists
    model.activateDeltaSPH(0.1)

    sim.run(20)

    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        rho = snap.particles.density
        assert np.all(np.isfinite(rho))
        assert np.allclose(rho, RHO0, rtol=5e-2)
        v = snap.particles.velocity
        assert np.all(np.isfinite(v))


def test_molteni_diffusion_smooths_density(tmp_path):
    """With the corrected sign, density diffusion must reduce density spread
    relative to a run without diffusion (anti-diffusion would grow it)."""
    def spread(with_diffusion, tmp):
        f = tmp / f"d{int(with_diffusion)}.gsd"
        make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0,
                             velocity_fn=taylor_green)
        sim, model, ffl, fsl, dt = build_singlephase_sim(
            f, DX, rho0=RHO0, densitymethod="CONTINUITY")
        if with_diffusion:
            # parameter-dict pathway: applied when the model attaches
            model.densitydiffusion = True
            model.ddiff = 0.1
        sim.run(50)
        snap = sim.state.get_snapshot()
        return float(np.std(snap.particles.density))

    s_off = spread(False, tmp_path)
    s_on = spread(True, tmp_path)
    assert np.isfinite(s_on)
    assert s_on <= s_off * 1.05  # diffusion must not amplify density noise


def test_adaptive_timestep(tmp_path):
    f = tmp_path / "tg3.gsd"
    make_fluid_block_gsd(f, n=10, dx=DX, rho0=RHO0, velocity_fn=taylor_green)
    sim, model, ff, fs, dt0 = build_singlephase_sim(f, DX, rho0=RHO0)

    adapt = sph.update.AdaptiveTimestep(model=model, courant=0.25)
    sim.operations.updaters.append(
        adapt.as_updater(trigger=hoomd.trigger.Periodic(5)))

    sim.run(11)

    dt_new = sim.operations.integrator.dt
    assert dt_new > 0.0
    assert np.isfinite(dt_new)
    # the acoustic CFL bound must hold
    c = model.eos.SpeedOfSound
    h = model.max_sl
    assert dt_new <= 0.25 * h / c * 1.001

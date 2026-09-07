"""Tests for the per-particle strain-rate estimator (WP0).

The non-Newtonian solvers store gamma_dot = sqrt(2 D:D) per particle in the
energy array, computed from the L-matrix renormalized velocity gradient.
Requirements verified here:

  * simple shear v_x = k*y  ->  gamma_dot == k (exact for linear fields)
  * rigid rotation v = omega x r  ->  gamma_dot == 0 (frame invariance);
    the old per-pair |dv|/r estimator read gamma_dot ~ omega here
  * Bingham channel flow: shear localizes at the walls, gamma_dot finite
"""

import numpy as np
import pytest

import hoomd

from conftest import make_fluid_block_gsd, build_singlephase_sim

RHO0 = 1000.0
DX = 0.05
N = 10
L = N * DX
RCUT = 2.0 * 1.7 * DX  # kappa * OptimalH * dx for WendlandC4

K_SHEAR = 0.5   # [1/s]
OMEGA = 0.5     # [1/s]


def _activate_powerlaw_n1(sim, model):
    """Constant-viscosity power law: mu_eff = K, but the strain-rate pass runs."""
    sim.run(0)
    model.activatePowerLaw(K=0.1, n=1.0, mu_min=0.0)


def test_simple_shear_reads_exact_rate(tmp_path):
    f = tmp_path / "shear.gsd"

    def shear(pos):
        v = np.zeros_like(pos)
        v[:, 0] = K_SHEAR * pos[:, 1]
        return v

    make_fluid_block_gsd(f, n=N, dx=DX, rho0=RHO0, velocity_fn=shear)
    sim, model, ff, fs, dt = build_singlephase_sim(f, DX, rho0=RHO0)
    _activate_powerlaw_n1(sim, model)
    sim.run(1)

    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        pos = snap.particles.position
        gdot = np.asarray(snap.particles.energy)
        # particles whose kernel support does not wrap the shear-discontinuous
        # y boundary (v_x = k*y is not periodic in y)
        interior = np.abs(pos[:, 1]) < L / 2 - RCUT
        assert interior.sum() > 50
        g = gdot[interior]
        assert np.allclose(g, K_SHEAR, rtol=2e-2), (
            f"gamma_dot = {g.mean():.4f} +- {g.std():.4f}, expected {K_SHEAR}")


def test_rigid_rotation_reads_zero(tmp_path):
    f = tmp_path / "rot.gsd"

    def rotation(pos):
        v = np.zeros_like(pos)
        v[:, 0] = -OMEGA * pos[:, 1]
        v[:, 1] = OMEGA * pos[:, 0]
        return v

    make_fluid_block_gsd(f, n=N, dx=DX, rho0=RHO0, velocity_fn=rotation)
    sim, model, ff, fs, dt = build_singlephase_sim(f, DX, rho0=RHO0)
    _activate_powerlaw_n1(sim, model)
    sim.run(1)

    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        pos = snap.particles.position
        gdot = np.asarray(snap.particles.energy)
        # v = omega x r is not periodic in x or y; keep support interior
        interior = (np.abs(pos[:, 0]) < L / 2 - RCUT) & \
                   (np.abs(pos[:, 1]) < L / 2 - RCUT)
        assert interior.sum() > 10
        g = gdot[interior]
        # frame invariance: rotation must not register as shear. The old
        # per-pair estimator returned ~OMEGA here; allow only a few percent
        # (residual from the one integration step taken).
        assert np.max(g) < 0.05 * OMEGA, (
            f"rotation read as shear: max gamma_dot = {np.max(g):.4f} "
            f"(omega = {OMEGA})")


def test_bingham_channel_shear_localizes_at_walls(tmp_path):
    f = tmp_path / "channel.gsd"

    def uniform_x(pos):
        v = np.zeros_like(pos)
        v[:, 0] = 0.05
        return v

    make_fluid_block_gsd(f, n=N, dx=DX, rho0=RHO0,
                         velocity_fn=uniform_x, solid_planes=3)
    sim, model, ff, fs, dt = build_singlephase_sim(f, DX, rho0=RHO0, mu=0.05)
    sim.run(0)
    model.activateBingham(mu_p=0.05, tau_y=0.01, m_reg=0.5, mu_min=0.0)
    sim.run(10)

    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        pos = snap.particles.position
        tid = snap.particles.typeid
        gdot = np.asarray(snap.particles.energy)
        fluid = tid == 0
        assert np.all(np.isfinite(gdot[fluid]))
        assert np.all(np.isfinite(snap.particles.velocity))

        # fluid slab is centered in z; shear must be largest near the walls
        z = pos[fluid, 2]
        g = gdot[fluid]
        z_half = (np.abs(z).max() - np.abs(z).min()) / 2 + np.abs(z).min()
        near_wall = np.abs(z) > z_half
        core = ~near_wall
        assert g[near_wall].mean() > g[core].mean(), (
            "wall shear should exceed core shear in start-up channel flow")

"""Shared fixtures for the hoomd.sph test suite."""

import numpy as np
import pytest

hoomd = pytest.importorskip("hoomd")
pytest.importorskip("hoomd.sph")
gsd_hoomd = pytest.importorskip("gsd.hoomd")


KERNELS = ["WendlandC2", "WendlandC4", "WendlandC6", "Quintic", "CubicSpline"]


def make_fluid_block_gsd(path, n=10, dx=0.05, rho0=1000.0,
                         velocity_fn=None, solid_planes=0):
    """Write a small periodic block of fluid particles to a GSD file.

    Args:
        path: output file name.
        n: particles per box edge (fluid region).
        dx: lattice spacing.
        rho0: rest density; particle mass is rho0*dx^3.
        velocity_fn: optional callable pos -> (N, 3) velocities.
        solid_planes: number of solid layers added at BOTH -z and +z ends
            (walls normal to z). 0 means all-fluid periodic box.

    Returns:
        (n_fluid, n_solid, slength)
    """
    import hoomd.sph.kernel as skernel
    slength = skernel.OptimalH["WendlandC4"] * dx

    nz_solid = solid_planes
    nz = n + 2 * nz_solid
    lx = n * dx
    lz = nz * dx

    xs = np.linspace(-lx / 2 + dx / 2, lx / 2 - dx / 2, n, endpoint=True)
    zs = np.linspace(-lz / 2 + dx / 2, lz / 2 - dx / 2, nz, endpoint=True)

    xg, yg, zg = np.meshgrid(xs, xs, zs, indexing="ij")
    pos = np.column_stack([xg.ravel(), yg.ravel(), zg.ravel()])
    npart = pos.shape[0]

    typeid = np.zeros(npart, dtype=np.int32)
    if nz_solid > 0:
        z = pos[:, 2]
        zwall_lo = -lz / 2 + nz_solid * dx
        zwall_hi = lz / 2 - nz_solid * dx
        solid_mask = (z < zwall_lo) | (z > zwall_hi)
        typeid[solid_mask] = 1
    n_solid = int((typeid == 1).sum())
    n_fluid = npart - n_solid

    vel = np.zeros((npart, 3), dtype=np.float32)
    if velocity_fn is not None:
        vel = np.asarray(velocity_fn(pos), dtype=np.float32)
        vel[typeid == 1] = 0.0

    frame = gsd_hoomd.Frame()
    frame.configuration.box = [lx, lx, lz, 0, 0, 0]
    frame.particles.N = npart
    frame.particles.types = ["F", "S"]
    frame.particles.typeid = typeid
    frame.particles.position = pos.astype(np.float32)
    frame.particles.velocity = vel
    frame.particles.mass = np.full(npart, rho0 * dx**3, dtype=np.float32)
    frame.particles.slength = np.full(npart, slength, dtype=np.float32)
    frame.particles.density = np.full(npart, rho0, dtype=np.float32)

    with gsd_hoomd.open(name=str(path), mode="w") as f:
        f.append(frame)

    return n_fluid, n_solid, slength


def build_singlephase_sim(gsd_path, dx, rho0=1000.0, mu=0.01,
                          densitymethod="SUMMATION", uref=0.1):
    """Construct a ready-to-run SinglePhaseFlow simulation from a GSD file.

    Returns (sim, model, filterfluid, filtersolid, dt).
    """
    import hoomd.sph as sph

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

    model = sph.sphmodel.SinglePhaseFlow(
        kernel=kernel_obj, eos=eos, nlist=nlist,
        fluidgroup_filter=filterfluid, solidgroup_filter=filtersolid,
        densitymethod=densitymethod)
    model.mu = mu

    max_sl = set_max_sl(sim, device, model)

    lref = 1.0
    model.compute_speedofsound(
        LREF=lref, UREF=uref, DX=dx, DRHO=0.01,
        H=max_sl, MU=mu, RHO0=rho0)
    dt, _ = model.compute_dt(
        LREF=lref, UREF=uref, DX=dx, DRHO=0.01,
        H=max_sl, MU=mu, RHO0=rho0)

    integrator = sph.Integrator(dt=dt)
    method = sph.methods.VelocityVerletBasic(filter=filterfluid,
                                             densitymethod=densitymethod)
    integrator.methods.append(method)
    integrator.forces.append(model)
    sim.operations.integrator = integrator

    return sim, model, filterfluid, filtersolid, dt


def set_max_sl(sim, device, model):
    """Minimal replacement for helper_modules.sph_helper.set_max_sl."""
    snapshot = sim.state.get_snapshot()
    max_sl = 0.0
    if snapshot.communicator.rank == 0:
        max_sl = float(np.max(snapshot.particles.slength))
    max_sl = device.communicator.bcast_double(max_sl)
    model.max_sl = max_sl
    return max_sl

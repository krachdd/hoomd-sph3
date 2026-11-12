==========
hoomd-sph3
==========

**hoomd-sph3** is a smoothed particle hydrodynamics (SPH) solver built on
`HOOMD-blue <https://hoomd-blue.readthedocs.io/>`_ 7.2.2. It solves the weakly
compressible Navier-Stokes equations for single-phase, two-phase, free-surface,
non-Newtonian, and density-gradient-driven flows. Simulations are set up and
controlled from Python and run on CPUs with MPI domain decomposition and OpenMP
threading.

The SPH solver lives in `hoomd.sph`, the neighbor search in `hoomd.nsearch`.
The HOOMD-blue core (`hoomd`) provides the simulation state, particle
filters, triggers, logging, and file output. The molecular dynamics, hard
particle Monte Carlo, and MPCD packages of HOOMD-blue are not part of this
build.

Resources
=========

- `GitHub repository <https://github.com/krachdd/hoomd-sph3>`_:
  Source code, benchmark cases, and issue tracker.
- :doc:`Compilation guide </building>`:
  Instructions for compiling **hoomd-sph3** and its I/O libraries.
- :doc:`Example simulations </examples>`:
  The benchmark cases in ``sph-simulations/``.
- `HOOMD-blue documentation <https://hoomd-blue.readthedocs.io/>`_:
  The upstream documentation of the core simulation engine.

Example script
==============

Poiseuille flow through a tube. The initial configuration ``tube_init.gsd``
holds fluid (type ``F``) and wall (type ``S``) particles on a lattice with
spacing ``dx``; ``create_input_geometries_tube.py`` in
``sph-simulations/00_singlephaseflow_benchmarks/01_poiseuille`` creates it.

.. code:: python

    import hoomd
    import hoomd.sph

    dx = 1e-4             # particle spacing [m]
    rho0 = 1000.0         # rest density [kg/m^3]
    viscosity = 0.01      # dynamic viscosity [Pa s]
    fx = 0.1              # body force [m/s^2]

    sim = hoomd.Simulation(device=hoomd.device.CPU())
    sim.create_state_from_gsd(filename='tube_init.gsd')

    kernel = hoomd.sph.kernel.WendlandC4()
    nlist = hoomd.nsearch.nlist.Cell(buffer=0.1 * dx,
                                     rebuild_check_delay=1,
                                     kappa=kernel.Kappa())
    eos = hoomd.sph.eos.Tait()
    eos.set_params(rho0, 0.01)

    model = hoomd.sph.sphmodel.SinglePhaseFlow(
        kernel=kernel,
        eos=eos,
        nlist=nlist,
        fluidgroup_filter=hoomd.filter.Type(['F']),
        solidgroup_filter=hoomd.filter.Type(['S']),
        densitymethod='SUMMATION')
    model.mu = viscosity
    model.gx = fx
    # Largest smoothing length in the initial configuration; it sets the
    # neighbor list cutoff (kappa * max_sl).
    h = hoomd.sph.kernel.OptimalH['WendlandC4'] * dx
    model.max_sl = h

    # Set the speed of sound and the time step from the flow parameters.
    lref = 1e-3
    uref = fx * lref**2 / (4 * viscosity / rho0)
    model.compute_speedofsound(LREF=lref, UREF=uref, DX=dx, DRHO=0.01, H=h,
                               MU=viscosity, RHO0=rho0)
    dt, _ = model.compute_dt(LREF=lref, UREF=uref, DX=dx, DRHO=0.01, H=h,
                             MU=viscosity, RHO0=rho0)

    integrator = hoomd.sph.Integrator(dt=dt)
    integrator.methods.append(hoomd.sph.methods.VelocityVerletBasic(
        filter=hoomd.filter.Type(['F']), densitymethod='SUMMATION'))
    integrator.forces.append(model)
    sim.operations.integrator = integrator

    sim.operations.writers.append(hoomd.write.GSD(
        filename='tube_run.gsd', trigger=hoomd.trigger.Periodic(100),
        mode='wb', dynamic=['property', 'momentum']))

    sim.run(10000)

The benchmark scripts in ``sph-simulations/`` read the maximum smoothing length
from the initial configuration with ``helper_modules/sph_helper.py`` and convert the output to VTU files for ParaView with
``helper_modules/gsd2vtu``. See :doc:`examples`.

.. rubric::
    Contents

.. toctree::
    :maxdepth: 2
    :caption: Guides

    getting-started

.. toctree::
   :maxdepth: 1
   :caption: Python API

   hoomd/module-sph
   hoomd/module-nsearch
   module-hoomd

.. toctree::
   :maxdepth: 2
   :caption: Reference

   documentation
   changes
   developers

   open-source

   indices

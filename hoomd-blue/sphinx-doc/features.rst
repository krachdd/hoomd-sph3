Features
========

SPH models
----------

`hoomd.sph.sphmodel` implements weakly compressible SPH models:

* Single-phase flow (`hoomd.sph.sphmodel.SinglePhaseFlow`) and its
  transport-velocity variant (`hoomd.sph.sphmodel.SinglePhaseFlowTV`).
* Two-phase flow of immiscible fluids with surface tension
  (`hoomd.sph.sphmodel.TwoPhaseFlow`, `hoomd.sph.sphmodel.TwoPhaseFlowTV`).
* Free-surface flow with free-surface detection, contact angles, and
  continuum surface force surface tension
  (`hoomd.sph.sphmodel.SinglePhaseFlowFS`).
* Density-gradient-driven flow with a transported scalar field and buoyancy
  (`hoomd.sph.sphmodel.SinglePhaseFlowGDGD`).
* Non-Newtonian rheology: power-law, Carreau, Bingham, and Herschel-Bulkley
  viscosity models.
* Density by summation or by the continuity equation, artificial viscosity,
  density diffusion (including delta-SPH), and Shepard renormalization.
* Solid walls represented by particles, with no-slip boundary conditions.

Building blocks
---------------

* Smoothing kernels (`hoomd.sph.kernel`): Wendland C2, C4, C6, quintic, and
  cubic spline.
* Equations of state (`hoomd.sph.eos`): Tait and linear.
* Integration methods (`hoomd.sph.methods`): velocity Verlet and
  kick-drift-kick with transport-velocity correction.
* Adaptive time stepping (`hoomd.sph.update.AdaptiveTimestep`) from the
  acoustic, viscous, and body-force conditions.
* Flow properties for logging (`hoomd.sph.compute`), such as kinetic energy,
  mean density, and forces on solids.
* Neighbor lists (`hoomd.nsearch.nlist`): cell, stencil, and tree.

Fluid and solid particles are selected by particle type with
`hoomd.filter` objects.

Parallel execution
------------------

hoomd-sph3 runs on CPUs. Use MPI to decompose the domain over several
processes, OpenMP threads to parallelize the SPH particle loops within a
process, or both. Set the number of threads with ``OMP_NUM_THREADS``. With
hybrid MPI and OpenMP, release the core binding of the MPI ranks (for example
``mpirun --bind-to none``) so that the threads of a rank run on separate
cores.

Results are not bitwise reproducible across different thread counts because
floating point sums are reordered. GPU execution is not implemented for the
SPH solver.

Input and output
----------------

The SPH per-particle fields (smoothing length, density, pressure, energy, and
four auxiliary vectors) are part of the simulation state and of the GSD files:

* `hoomd.write.GSD` writes trajectories through MPI rank 0 in single or
  double precision. Read them with ``gsd-sph`` (``gsd.hoomd``), the GSD 5
  Python package extended with the SPH fields.
* `hoomd.write.PGSD` writes trajectories with parallel MPI-IO, every rank
  writes its own particles. Read them with ``pgsd.hoomd`` or ``gsd.hoomd``.
* `hoomd.Simulation.create_state_from_gsd` and
  `hoomd.Simulation.create_state_from_pgsd`
  read the initial configuration.
* ``helper_modules/gsd2vtu`` converts GSD files to VTU files for ParaView.

HOOMD-blue core
---------------

hoomd-sph3 uses the core of HOOMD-blue 7.2.2: the Python simulation interface
(`hoomd.Simulation`, `hoomd.State`), particle filters, triggers, variants,
logging (`hoomd.logging`), custom Python actions (`hoomd.custom`), and the
writers in `hoomd.write`. See the
`HOOMD-blue documentation <https://hoomd-blue.readthedocs.io/>`_ for the core
concepts.

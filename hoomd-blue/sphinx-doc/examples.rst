Example simulations
===================

``sph-simulations/`` in the repository contains benchmark cases for every SPH
model. Each case has a script that creates the initial configuration and one or
more scripts that run the simulation, for example::

    $ cd sph-simulations/00_singlephaseflow_benchmarks/01_poiseuille
    $ python3 create_input_geometries_tube.py 20
    $ mpirun -n 4 python3 run_tube.py 20 poiseuille_flow_20_32_32_vs_5e-05_init.gsd 10001

The argument of the create script sets the resolution (particles per
reference length). ``run.sh`` in each case directory lists the commands used for
the reference runs. The run scripts import modules from ``helper_modules/``,
so source ``addsphpath.sh`` first (see :doc:`building`).

Single-phase flow
-----------------

``00_singlephaseflow_benchmarks``:

.. list-table::
   :header-rows: 1

   * - Case
     - Description
   * - ``01_poiseuille``
     - Poiseuille flow in a tube
   * - ``02_parallel_plates``
     - Poiseuille flow between parallel plates
   * - ``03_couette_flow``
     - Couette flow
   * - ``04_channel_flow``
     - Pressure-driven channel flow
   * - ``05_liddriven_cavity``
     - Lid-driven cavity
   * - ``06_spherepackings``
     - Flow through sphere packings
   * - ``07_adami_cylinder1``, ``08_adami_cylinder2``
     - Flow around cylinders (Adami et al.)
   * - ``09_drag_cylinder``
     - Drag on a cylinder
   * - ``10_hydrostatic``
     - Hydrostatic pressure
   * - ``11_taylor_green``
     - Taylor-Green vortex decay
   * - ``12_stokes_second``
     - Stokes' second problem (oscillating plate)
   * - ``13_bcc_spherepacking_permeability``
     - Permeability of a BCC sphere packing

Two-phase flow
--------------

``01_twophaseflow_benchmarks``:

.. list-table::
   :header-rows: 1

   * - Case
     - Description
   * - ``01_layered_couette``
     - Layered two-fluid Couette flow
   * - ``02_static_droplet``
     - Static droplet, Laplace pressure jump
   * - ``03_capillary_rise``
     - Capillary rise
   * - ``04_rising_bubble``
     - Rising bubble
   * - ``05_rayleigh_taylor``
     - Rayleigh-Taylor instability
   * - ``06_sessile_droplet_shear``
     - Sessile droplet under shear
   * - ``07_sessile_droplet_snapback``
     - Sessile droplet snap-back
   * - ``08_h2_bubble_brine_shear``
     - Hydrogen bubble in brine under shear

Density-gradient-driven flow
----------------------------

``02_densitygradient_driven_flow``:

.. list-table::
   :header-rows: 1

   * - Case
     - Description
   * - ``01_thermal_couette``
     - Thermal Couette flow
   * - ``02_heated_cavity``
     - Differentially heated cavity
   * - ``03_boussinesq_rayleigh_taylor``
     - Rayleigh-Taylor instability, Boussinesq approximation

Non-Newtonian flow
------------------

``03_nonnewtonian_benchmarks``:

.. list-table::
   :header-rows: 1

   * - Case
     - Description
   * - ``01_powerlaw_poiseuille``
     - Power-law fluid, Poiseuille flow
   * - ``02_bingham_poiseuille``
     - Bingham fluid, Poiseuille flow
   * - ``03_twophase_powerlaw_couette``
     - Two-phase power-law Couette flow
   * - ``04_bingham_slump_test``
     - Bingham slump test
   * - ``05_viscoplastic_invasion``
     - Viscoplastic invasion
   * - ``06_bingham_slug_mobilization``
     - Mobilization of a Bingham slug

Free-surface flow
-----------------

``04_free_surface_flow`` (`hoomd.sph.sphmodel.SinglePhaseFlowFS`):

.. list-table::
   :header-rows: 1

   * - Case
     - Description
   * - ``01_hydrostatic_fs``
     - Hydrostatic column with a free surface
   * - ``02_dam_break``
     - Dam break
   * - ``03_sessile_droplet``
     - Sessile droplet on a substrate

Input and output
----------------

``05_io_benchmarks`` measures the read and write performance of
`hoomd.write.GSD` (``01_gsd``) and `hoomd.write.PGSD` (``02_pgsd``) on domains
created by ``00_create_domains``.

Post-processing
---------------

``helper_modules/gsd2vtu`` converts the GSD output to VTU and PVD files for
ParaView. Most run scripts call it at the end of the run.

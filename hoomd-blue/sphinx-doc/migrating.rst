Migrating to the latest version
===============================

Migrating from hoomd-sph3 based on HOOMD-blue 5.2
-------------------------------------------------

hoomd-sph3 now builds on HOOMD-blue 7.2.2, gsd-sph 5.0.1, and pgsd-sph with the
GSD 2.1 file format. The SPH solver produces the same results as before;
these changes affect scripts and post-processing:

* Rebuild all three parts (see :doc:`building`). The ``gsd`` and ``pgsd``
  Python packages must match the new ``hoomd`` build.
* `hoomd.write.GSD` writes GSD schema 2.0 files. Read them with gsd-sph 5.0 or
  newer; older ``gsd`` versions cannot read them. ``pgsd.hoomd`` reads both
  schema versions.
* `hoomd.write.GSD` has a ``precision`` argument. Set ``precision='double'``
  to store positions, velocities, and the SPH fields (density, pressure,
  energy, smoothing length, auxiliary vectors) in double precision. The
  default is ``'single'``, as before.
* New GSD and PGSD files use file format version 2.1. Opening an older file
  for writing (for example with mode ``'r+'`` or ``'a'``) upgrades its header
  to 2.1; open files with mode ``'r'`` to leave them unchanged.
* Neither ``hoomd``, ``gsd``, nor ``pgsd`` installs a ``SIGTERM`` handler
  anymore. Install your own handler if a job script needs a clean shutdown
  on ``SIGTERM``.
* ``gsd`` no longer flushes implicitly. Call ``flush()`` or close the file
  before reading a file that is still open for writing.
* ``gsd.hoomd.Frame`` has no ``state`` chunks in the file. ``Frame.state``
  remains available as in-memory storage.

Changes in HOOMD-blue
---------------------

The following notes from HOOMD-blue apply to the core features that
hoomd-sph3 uses. See the
`HOOMD-blue migration guide <https://hoomd-blue.readthedocs.io/en/latest/migrating.html>`_
for older versions.

Migrating to HOOMD-blue 7
^^^^^^^^^^^^^^^^^^^^^^^^^

Breaking changes to existing functionalities
""""""""""""""""""""""""""""""""""""""""""""

* HOOMD-blue 7 writes GSD schema 2.0 files. You may need to update any tools that read GSD
  to their latest version to ensure compatibility.

Migrating to HOOMD-blue 6
^^^^^^^^^^^^^^^^^^^^^^^^^

Breaking changes to existing functionalities
""""""""""""""""""""""""""""""""""""""""""""

Check your scripts for anywhere you set type parameters with two tuples. For example:
``lj.r_cut[('A', 'B'), ('C', 'D')] = ...``. In HOOMD-blue 5.x and prior versions, this
would take the product and set the A-C, A-D, B-C, and B-D parameters. In HOOMD-blue
6 and later, this will set the A-B and C-D parameters. Identify what parameters  you
intend to set and adjust your script accordingly. Many users expect the new behavior
because it is consistent with setting 3 (or more) tuples:
``lj.r_cut[('A', 'B'), ('C', 'D'), ('E', 'F')] = ...`` sets parameters for A-B, C-D, and E-F.

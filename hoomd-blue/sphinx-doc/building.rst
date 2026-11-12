Building from source
====================

hoomd-sph3 is built from source together with its two I/O libraries:

* **pgsd-sph**: MPI-IO version of the GSD file layer, used by
  `hoomd.write.PGSD` and the ``pgsd`` Python package.
* **gsd-sph**: the GSD 5 Python package (``gsd``) extended with the SPH
  particle fields, used for reading and post-processing.

Both are git submodules in ``dependencies/``.

Prerequisites
-------------

- C++17 compiler
- CMake
- MPI (tested with OpenMPI)
- Python 3.10 or newer with NumPy, Cython, and mpi4py
- pybind11 2.12 or newer
- Eigen
- cereal
- OpenMP (optional): threads the SPH particle loops. The build stays
  serial when CMake does not find OpenMP.

``requirements.txt`` in the repository root lists the Python packages, including
those for post-processing (``vtk``, ``pyevtk``) and for building this
documentation.

Get the source
--------------

Clone the repository with its submodules::

    $ git clone --recurse-submodules https://github.com/krachdd/hoomd-sph3.git
    $ cd hoomd-sph3

Build
-----

``buildall.sh`` links the pgsd C sources into ``hoomd-blue/hoomd/extern``
(``link_pgsd_module.sh``) and builds pgsd-sph, gsd-sph, and HOOMD-blue with
the SPH components, in this order, from scratch::

    $ bash buildall.sh

``fast_buildall.sh`` rebuilds without deleting the build directories.

.. note::

    ``link_pgsd_module.sh`` runs ``git checkout main`` in all submodules.
    When you work on other submodule branches, create the links manually and
    build the three parts with CMake as shown below.

The steps of ``buildall.sh`` are::

    $ ln -s $PWD/dependencies/pgsd-sph/pgsd/pgsd/pgsd.c hoomd-blue/hoomd/extern/pgsd.c
    $ ln -s $PWD/dependencies/pgsd-sph/pgsd/pgsd/pgsd.h hoomd-blue/hoomd/extern/pgsd.h
    $ CC=mpicc CXX=mpicxx cmake -B dependencies/pgsd-sph/pgsd/build -S dependencies/pgsd-sph/pgsd
    $ cmake --build dependencies/pgsd-sph/pgsd/build
    $ cmake -B dependencies/gsd-sph/gsd/build -S dependencies/gsd-sph/gsd
    $ cmake --build dependencies/gsd-sph/gsd/build
    $ cmake -B hoomd-blue/build -S hoomd-blue
    $ cmake --build hoomd-blue/build -j 4

CMake options for the HOOMD-blue part:

- ``ENABLE_MPI``: build with MPI support (default: ``on``).
- ``CMAKE_BUILD_TYPE``: ``Release`` (default), ``RelWithDebInfo``, or ``Debug``.
- ``Python_EXECUTABLE``: the Python interpreter to build against.
- ``BUILD_TESTING``: build the C++ unit tests (default: ``off``).

The molecular dynamics, HPMC, MPCD, and metal packages of HOOMD-blue are not
included and their ``BUILD_*`` options have no effect. GPU builds
(``ENABLE_GPU``) are not supported by the SPH components.

Set up the environment
----------------------

The packages are used from their build directories. Source ``addsphpath.sh``
from the repository root to add them and ``helper_modules/`` to
``PYTHONPATH``::

    $ source addsphpath.sh
    $ python3 -c "import hoomd, hoomd.sph; print(hoomd.version.version)"

Run the tests
-------------

Run the SPH unit tests from their directory in the build tree (the test
modules import ``conftest`` directly)::

    $ cd hoomd-blue/build/hoomd/sph/pytest
    $ python3 -m pytest

Run the HOOMD-blue core tests from the build directory::

    $ cd hoomd-blue/build
    $ python3 -m pytest hoomd/pytest

Core tests for features removed in hoomd-sph3 (bonds, angles, and other
topology, MD integrators) fail.

Run the gsd-sph and pgsd-sph tests from their build directories::

    $ python3 -m pytest dependencies/gsd-sph/gsd/build/gsd/test
    $ python3 -m pytest dependencies/pgsd-sph/pgsd/build/pgsd/test

The pgsd-sph tests write files on 1, 2, and 4 MPI ranks with ``mpirun``.

Build the documentation
-----------------------

Install the packages listed in ``hoomd-blue/sphinx-doc/requirements.in`` and
build the HTML documentation::

    $ cd hoomd-blue/sphinx-doc
    $ python3 generate-toctree.py
    $ sphinx-build -b html . _build/html

``generate-toctree.py`` creates the API pages for new modules and classes. Run
it after adding a class to the ``__all__`` list of a module.

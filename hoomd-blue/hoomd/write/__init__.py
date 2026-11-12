# Copyright (c) 2009-2026 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Writers write the state of the simulation, logger quantities, or calculated
results to output files or streams:

* `GSD` saves the simulation trajectory to a file, including the SPH
  per-particle fields.
* `PGSD` saves the trajectory with parallel MPI-IO: every rank writes its own
  particles. Read the files with the ``pgsd`` or ``gsd`` Python packages.
* `Burst` provides a sliding window of a simulation trajectory wrote out at
  `Burst.dump` for use in selective high frequency trajectory data.
* Combine `GSD` with a `hoomd.logging.Logger` to save system properties or
  per-particle calculated results.
* Use `HDF5Log` to store logged data in HDF5 resizable datasets.
* Use `Table` to display the status of the simulation periodically to standard
  out.
* Implement custom output formats with `CustomWriter`.

Writers do not modify the system state.

Tip:
    `OVITO <https://www.ovito.org/>`_ has native support for GSD files,
    including logged per-particle array quantities and particle shapes.
"""

from hoomd.write.custom_writer import CustomWriter
from hoomd.write.gsd import GSD
from hoomd.write.pgsd import PGSD
from hoomd.write.gsd_burst import Burst
# from hoomd.write.dcd import DCD
from hoomd.write.table import Table
from hoomd.write.hdf5 import HDF5Log

__all__ = [
    # "DCD",
    "GSD",
    "PGSD",
    "Burst",
    "CustomWriter",
    "HDF5Log",
    "Table",
]

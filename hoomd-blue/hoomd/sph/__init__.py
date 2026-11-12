# Copyright (c) 2025-2026 David Krach, Daniel Rostan.
# All rights reserved. BSD-3-Clause license, see repository LICENSE.
#
# maintainer: dkrach, david.krach@mib.uni-stuttgart.de

"""Smoothed particle hydrodynamics (SPH).

`hoomd.sph` solves the weakly compressible Navier-Stokes equations with the
SPH method. Fluid and solid (boundary) particles are distinguished by particle
type and selected with `hoomd.filter` objects.

To run an SPH simulation, assign a `hoomd.sph.Integrator` to
`hoomd.Simulation.operations` and give it:

* one or more integration methods from `hoomd.sph.methods` that advance the
  fluid particles,
* an SPH model from `hoomd.sph.sphmodel` (a `hoomd.sph.force.Force`) that
  computes density, pressure, and the momentum equation. A model needs a
  smoothing kernel (`hoomd.sph.kernel`), an equation of state
  (`hoomd.sph.eos`), and a neighbor list (`hoomd.nsearch.nlist`).

Use `hoomd.sph.compute` to log flow properties and
`hoomd.sph.update.AdaptiveTimestep` to adapt the time step during a run.
SPH-specific per-particle fields (density, pressure, energy, smoothing length,
auxiliary vectors) are stored in GSD files written by `hoomd.write.GSD` and in
parallel by `hoomd.write.PGSD`.
"""

from hoomd.sph import integrate
from hoomd.sph.integrate import Integrator
from hoomd.sph import force
from hoomd.sph import eos
from hoomd.sph import sphmodel
from hoomd.sph import constrain
from hoomd.sph import kernel
from hoomd.sph import methods
from hoomd.sph import compute
from hoomd.sph import update
from hoomd.sph.half_step_hook import HalfStepHook


__all__ = [
    "Integrator",
    "force",
    "eos",
    "sphmodel",
    "constrain",
    "kernel",
    "methods",
    "compute",
    "update",
    "HalfStepHook",
]
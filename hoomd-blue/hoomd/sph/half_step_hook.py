# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Perform user defined computations during the half-step of a \
`hoomd.sph.Integrator`.

HalfStepHook can be subclassed to define custom operations at the middle of
each integration step. Examples of use cases include evaluating collective
variables or biasing the simulation.
"""

from hoomd.sph import _sph


class HalfStepHook(_sph.HalfStepHook):
    """HalfStepHook base class.

    HalfStepHook provides an interface to perform computations during the
    half-step of a hoomd.sph.Integrator.
    """

    def update(self, timestep):
        """Called during the half-step of a `hoomd.sph.Integrator`.

        This method should provide the implementation of any computation that
        the user wants to execute at each timestep in the middle of the
        integration routine.
        """
        raise TypeError(
            "Use a hoomd.sph.HalfStepHook derived class implementing the "
            "corresponding update method.")
# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Constraint force classes apply forces and the resulting virial to particles that
force classes apply forces and the resulting virial to particles that
enforce specific constraints on the positions of the particles. The constraint
is satisfied at all times, so there is no potential energy associated with the
constraint.

Each constraint removes a number of degrees of freedom from the system.
`hoomd.md.compute.ThermodynamicQuantities` accounts for these lost degrees of
freedom when computing kinetic temperature and pressure. See
`hoomd.State.update_group_dof` for details on when the degrees of freedom for a
group are calculated.

Warning:
    Do not apply multiple constraint class instances to the same particle. Each
    instance solves for its constraints independently.
"""
import hoomd.sph
from hoomd.sph import _sph
from hoomd.data.parameterdicts import ParameterDict, TypeParameterDict
from hoomd.data.typeparam import TypeParameter
from hoomd.data.typeconverter import OnlyIf, to_type_converter
from hoomd.sph.force import Force
import hoomd
# from hoomd.operation import _HOOMDBaseObject



class Constraint(Force):
    """Constraint force base class.

    `Constraint` is the base class for all constraint forces.

    Warning:
        This class should not be instantiated by users. The class can be used
        for `isinstance` or `issubclass` checks.
    """



    def _attach_hook(self):
        """Create the c++ mirror class."""
        if isinstance(self._simulation.device, hoomd.device.CPU):
            cpp_cls = getattr(_sph, self._cpp_class_name)
        else:
            cpp_cls = getattr(_sph, self._cpp_class_name + "GPU")

        self._cpp_obj = cpp_cls(self._simulation.state._cpp_sys_def)

        # super()._attach()


class Distance(Constraint):
    """Constrain pairwise particle distances.

    Args:
        tolerance (float): Relative tolerance for constraint violation warnings.

    `Distance` applies forces between particles that constrain the distances
    between particles to specific values. The algorithm implemented is described
    in:

    1. M. Yoneya, H. J. C. Berendsen, and K. Hirasawa, "A Non-Iterative
       Matrix Method for Constraint Molecular Dynamics Simulations," Molecular
       Simulation, vol. 13, no. 6, pp. 395--405, 1994.
    2. M. Yoneya, "A Generalized Non-iterative Matrix Method for Constraint
       Molecular Dynamics Simulations," Journal of Computational Physics,
       vol. 172, no. 1, pp. 188--197, Sep. 2001.

    Each distance constraint takes the form:

    .. math::

        \\chi_{ij}(r) = \\mathrm{minimum\\_image}(\\vec{r}_j - \\vec{r}_i)^2
            - d_{ij}^2 = 0

    Where :math:`i` and :math:`j` are the the particle tags in the
    ``constraint_group`` and :math:`d_{ij}` is the constraint distance as given
    by the `system state <hoomd.State>`_.

    The method sets the second derivative of the Lagrange multipliers with
    respect to time to zero, such that both the distance constraints and their
    time derivatives are conserved within the accuracy of the Velocity Verlet
    scheme (:math:`O(\\delta t^2)`. It solves the corresponding linear system of
    equations to determine the force. The constraints are satisfied at :math:`t
    + 2 \\delta t`, so the scheme is self-correcting and avoids drifts.

    Add an instance of `Distance` to the integrator constraints list
    `hoomd.md.Integrator.constraints` to apply the force during the simulation.

    Warning:
        In MPI simulations, it is an error when molecules defined by constraints
        extend over more than half the local domain size because all particles
        connected through constraints will be communicated between ranks as
        ghost particles.

    Note:
        `tolerance` sets the tolerance to detect constraint violations and
        issue a warning message. It does not influence the computation of the
        constraint force.

    Attributes:
        tolerance (float): Relative tolerance for constraint violation warnings.
    """

    _cpp_class_name = "ForceDistanceConstraint"

    def __init__(self, tolerance=1e-3):
        self._param_dict.update(ParameterDict(tolerance=float(tolerance)))



__all__ = [
    "Constraint",
    "Distance",
]
# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Apply forces to particles."""

from abc import abstractmethod

import hoomd
import hoomd.sph
from hoomd import _hoomd
from hoomd.sph import _sph
from hoomd.operation import Compute
from hoomd.logging import log
from hoomd.data.typeparam import TypeParameter
from hoomd.data.typeconverter import OnlyTypes
from hoomd.data.parameterdicts import ParameterDict, TypeParameterDict
from hoomd.filter import ParticleFilter
# from hoomd.md.manifold import Manifold
import numpy


class Force(Compute):
    r"""Defines a force for molecular dynamics simulations.

    `Force` is the base class for all molecular dynamics forces and provides
    common methods.

    A `Force` class computes the force and torque on each particle in the
    simulation state :math:`\vec{F}_i` and :math:`\vec{\tau}_i`. With a few
    exceptions (noted in the documentation of the specific force classes),
    `Force` subclasses also compute the contribution to the system's potential
    energy :math:`U` and the the virial tensor :math:`W`. `Force` breaks the
    computation of the total system :math:`U` and :math:`W` into per-particle
    and additional terms as detailed in the documentation for each specific
    `Force` subclass.

    .. math::

        U & = U_\mathrm{additional} + \sum_{i=0}^{N_\mathrm{particles}-1} U_i \\
        W & = W_\mathrm{additional} + \sum_{i=0}^{N_\mathrm{particles}-1} W_i

    `Force` represents virial tensors as six element arrays listing the
    components of the tensor in this order:

    .. math::

        (W^{xx}, W^{xy}, W^{xz}, W^{yy}, W^{yz}, W^{zz}).

    The components of the virial tensor for a force on a single particle are:

    .. math::

        W^{kl}_i = F^k \cdot r_i^l

    where the superscripts select the x,y, and z components of the vectors.
    To properly account for periodic boundary conditions, pairwise interactions
    evaluate the virial:

    .. math::

        W^{kl}_i = \frac{1}{2} \sum_j F^k_{ij} \cdot
        \mathrm{minimum\_image}(\vec{r}_j - \vec{r}_i)^l

    Warning:
        This class should not be instantiated by users. The class can be used
        for `isinstance` or `issubclass` checks.
    """

    def __init__(self):
        self._in_context_manager = False

    @log(requires_run=True)
    def energy(self):
        """float: The potential energy :math:`U` of the system from this force \
        :math:`[\\mathrm{energy}]`."""
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.calcEnergySum()

    @log(category="particle", requires_run=True)
    def energies(self):
        """(*N_particles*, ) `numpy.ndarray` of ``float``: Energy \
        contribution :math:`U_i` from each particle :math:`[\\mathrm{energy}]`.

        Attention:
            In MPI parallel execution, the array is available on rank 0 only.
            `energies` is `None` on ranks >= 1.
        """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.getEnergies()

    @log(requires_run=True)
    def additional_energy(self):
        """float: Additional energy term :math:`U_\\mathrm{additional}` \
        :math:`[\\mathrm{energy}]`."""
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.getExternalEnergy()

    @log(category="particle", requires_run=True)
    def forces(self):
        """(*N_particles*, 3) `numpy.ndarray` of ``float``: The \
        force :math:`\\vec{F}_i` applied to each particle \
        :math:`[\\mathrm{force}]`.

        Attention:
            In MPI parallel execution, the array is available on rank 0 only.
            `forces` is `None` on ranks >= 1.
        """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.getForces()

    # @log(category="particle", requires_run=True)
    # def torques(self):
    #     """(*N_particles*, 3) `numpy.ndarray` of ``float``: The torque \
    #     :math:`\\vec{\\tau}_i` applied to each particle \
    #     :math:`[\\mathrm{force} \\cdot \\mathrm{length}]`.

    #     Attention:
    #         In MPI parallel execution, the array is available on rank 0 only.
    #         `torques` is `None` on ranks >= 1.
    #     """
    #     self._cpp_obj.compute(self._simulation.timestep)
    #     return self._cpp_obj.getTorques()

    # @log(category="particle", requires_run=True)
    # def virials(self):
    #     """(*N_particles*, 6) `numpy.ndarray` of ``float``: Virial tensor \
    #     contribution :math:`W_i` from each particle :math:`[\\mathrm{energy}]`.

    #     Attention:
    #         To improve performance `Force` objects only compute virials when
    #         needed. When not computed, `virials` is `None`. Virials are computed
    #         on every step when using a `md.methods.NPT` or `md.methods.NPH`
    #         integrator, on steps where a writer is triggered (such as
    #         `write.GSD` which may log pressure or virials), or when
    #         `Simulation.always_compute_pressure` is `True`.

    #     Attention:
    #         In MPI parallel execution, the array is available on rank 0 only.
    #         `virials` is `None` on ranks >= 1.
    #     """
    #     self._cpp_obj.compute(self._simulation.timestep)
    #     return self._cpp_obj.getVirials()

    # @log(category="sequence", requires_run=True)
    # def additional_virial(self):
    #     """(1, 6) `numpy.ndarray` of ``float``: Additional virial tensor \
    #     term :math:`W_\\mathrm{additional}` :math:`[\\mathrm{energy}]`."""
    #     self._cpp_obj.compute(self._simulation.timestep)
    #     virial = []
    #     for i in range(6):
    #         virial.append(self._cpp_obj.getExternalVirial(i))
    #     return numpy.array(virial, dtype=numpy.float64)

    @property
    def cpu_local_force_arrays(self):
        """hoomd.sph.data.ForceLocalAccess: Expose force arrays on the CPU.

        Provides direct access to the force, potential energy, torque, and
        virial data of the particles in the system on the cpu through a context
        manager. All data is MPI rank-local.

        The `hoomd.sph.data.ForceLocalAccess` object returned by this property
        has four arrays through which one can modify the force data:

        Note:
            The local arrays are read only for built-in forces. Use `Custom` to
            implement custom forces.

        Examples::

            with self.cpu_local_force_arrays as arrays:
                arrays.force[:] = ...
                arrays.potential_energy[:] = ...
                arrays.torque[:] = ...
                arrays.virial[:] = ...
        """
        if self._in_context_manager:
            raise RuntimeError("Cannot enter cpu_local_force_arrays context "
                               "manager inside another local_force_arrays "
                               "context manager")
        if not self._attached:
            raise hoomd.error.DataAccessError("cpu_local_force_arrays")
        return hoomd.sph.data.ForceLocalAccess(self, self._simulation.state)

    @property
    def gpu_local_force_arrays(self):
        """hoomd.sph.data.ForceLocalAccessGPU: Expose force arrays on the GPU.

        Provides direct access to the force, potential energy, torque, and
        virial data of the particles in the system on the gpu through a context
        manager. All data is MPI rank-local.

        The `hoomd.sph.data.ForceLocalAccessGPU` object returned by this property
        has four arrays through which one can modify the force data:

        Note:
            The local arrays are read only for built-in forces. Use `Custom` to
            implement custom forces.

        Examples::

            with self.gpu_local_force_arrays as arrays:
                arrays.force[:] = ...
                arrays.potential_energy[:] = ...
                arrays.torque[:] = ...
                arrays.virial[:] = ...

        Note:
            GPU local force data is not available if the chosen device for the
            simulation is `hoomd.device.CPU`.
        """
        if not isinstance(self._simulation.device, hoomd.device.GPU):
            raise RuntimeError(
                "Cannot access gpu_local_force_arrays without a GPU device")
        if self._in_context_manager:
            raise RuntimeError(
                "Cannot enter gpu_local_force_arrays context manager inside "
                "another local_force_arrays context manager")
        if not self._attached:
            raise hoomd.error.DataAccessError("gpu_local_force_arrays")
        return hoomd.sph.data.ForceLocalAccessGPU(self, self._simulation.state)


class Custom(Force):
    """Custom forces implemented in python.

    Derive a custom force class from `Custom`, and override the `set_forces`
    method to compute forces on particles. Users have direct, zero-copy access
    to the C++ managed buffers via either the `cpu_local_force_arrays` or
    `gpu_local_force_arrays` property. Choose the property that corresponds to
    the device you wish to alter the data on. In addition to zero-copy access to
    force buffers, custom forces have access to the local snapshot API via the
    ``_state.cpu_local_snapshot`` or the ``_state.gpu_local_snapshot`` property.

    See Also:
      See the documentation in `hoomd.State` for more information on the local
      snapshot API.

    Examples::

        class MyCustomForce(hoomd.force.Custom):
            def __init__(self):
                super().__init__(aniso=True)

            def set_forces(self, timestep):
                with self.cpu_local_force_arrays as arrays:
                    arrays.force[:] = -5
                    arrays.torque[:] = 3
                    arrays.potential_energy[:] = 27
                    arrays.virial[:] = np.arange(6)[None, :]

    In addition, since data is MPI rank-local, there may be ghost particle data
    associated with each rank. To access this read-only ghost data, access the
    property name with either the prefix ``ghost_`` of the suffix
    ``_with_ghost``.

    Note:
        Pass ``aniso=True`` to the `sph.force.Custom` constructor if your custom
        force produces non-zero torques on particles.

    Examples::

        class MyCustomForce(hoomd.force.Custom):
            def __init__(self):
                super().__init__()

            def set_forces(self, timestep):
                with self.cpu_local_force_arrays as arrays:
                    # access only the ghost particle forces
                    ghost_force_data = arrays.ghost_force

                    # access torque data on this rank and ghost torque data
                    torque_data = arrays.torque_with_ghost

    Note:
        When accessing the local force arrays, always use a context manager.

    Note:
        The shape of the exposed arrays cannot change while in the context
        manager.

    Note:
        All force data buffers are MPI rank local, so in simulations with MPI,
        only the data for a single rank is available.

    Note:
        Access to the force buffers is constant (O(1)) time.

    .. versionchanged:: 3.1.0
        `Custom` zeros the force, torque, energy, and virial arrays before
        calling the user-provided `set_forces`.
    """

    def __init__(self): # , aniso=False):
        super().__init__()
        # self._aniso = aniso

        self._state = None  # to be set on attaching

    def _attach_hook(self):
        self._state = self._simulation.state
        self._cpp_obj = _sph.CustomForceCompute(self._state._cpp_sys_def,
                                               self.set_forces) # , self._aniso)

    @abstractmethod
    def set_forces(self, timestep):
        """Set the forces in the simulation loop.

        Args:
            timestep (int): The current timestep in the simulation.
        """
        pass




__all__ = [
    "Custom",
    "Force",
]
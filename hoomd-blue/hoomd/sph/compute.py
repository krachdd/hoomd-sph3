"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""

"""Compute properties of mechanical System .

The SPH compute classes compute instantaneous properties of the simulation state
and provide results as loggable quantities for use with `hoomd.logging.Logger`
or by direct access via the Python API.
"""

from hoomd.sph import _sph
from hoomd.operation import Compute
from hoomd.data.parameterdicts import ParameterDict
from hoomd.logging import log
import hoomd

class SinglePhaseFlowBasicProperties(Compute):
    """Compute mechanical properties of a subset of the system.

    Args:
        filter (`hoomd.filter`): Particle filter to compute thermodynamic
            properties for.

    `SinglePhaseFlowBasicProperties` acts on a subset of particles in the system and
    calculates mechanical properties of those particles. Add a
    `SinglePhaseFlowBasicProperties` instance to a logger to save these quantities to a
    file, see `hoomd.logging.Logger` for more details.
    
    --------------------------------------------------------------------------------
    I keep this info here for later use of rigid bodies
    Note:
        For compatibility with `hoomd.md.constrain.Rigid`,
        `ThermodynamicQuantities` performs all sums
        :math:`\\sum_{i \\in \\mathrm{filter}}` over free particles and
        rigid body centers - ignoring constituent particles to avoid double
        counting.
    .---------------------------------------------------------------------------------

    Examples::

        f = filter.Type('A')
        compute.ThermodynamicQuantities(filter=f)
    """

    def __init__(self, filter):
        super().__init__()
        self._filter = filter

    def _attach_hook(self):
        if isinstance(self._simulation.device, hoomd.device.CPU):
            spfbasic_cls = _sph.ComputeSPFBasicProperties
        else:
            spfbasic_cls = _sph.ComputeSPFBasicPropertiesGPU
        group = self._simulation.state._get_group(self._filter)
        self._cpp_obj = spfbasic_cls(self._simulation.state._cpp_sys_def, group)
        # super()._attach()

    @log(requires_run=True)
    def abs_velocity(self):
        """Absolute velocity (norm of the vector) of the subset """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.abs_velocity

    @log(requires_run=True)
    def num_particles(self):
        """Number of particles :math:`N` in the subset."""
        return self._cpp_obj.num_particles

    @log(requires_run=True)
    def volume(self):
        """Volume :math:`V` of the simulation box (area in 2D) \
        :math:`[\\mathrm{length}^{D}]`."""
        return self._cpp_obj.volume

    @log(requires_run=True)
    def fluid_vel_x_sum(self):
        """Sum of Fluid Particle velocity in xdir """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.fluid_vel_x_sum

    @log(requires_run=True)
    def fluid_vel_y_sum(self):
        """Sum of Fluid Particle velocity in xdir """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.fluid_vel_y_sum

    @log(requires_run=True)
    def fluid_vel_z_sum(self):
        """Sum of Fluid Particle velocity in xdir """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.fluid_vel_z_sum

    @log(requires_run=True)
    def mean_density(self):
        """Sum of Fluid Particle velocity in xdir """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.mean_density

class SinglePhaseFlowNNBasicProperties(Compute):
    """Compute mechanical properties of a subset of the system.

    Args:
        filter (`hoomd.filter`): Particle filter to compute thermodynamic
            properties for.

    `SinglePhaseFlowBasicProperties` acts on a subset of particles in the system and
    calculates mechanical properties of those particles. Add a
    `SinglePhaseFlowBasicProperties` instance to a logger to save these quantities to a
    file, see `hoomd.logging.Logger` for more details.
    
    --------------------------------------------------------------------------------
    I keep this info here for later use of rigid bodies
    Note:
        For compatibility with `hoomd.md.constrain.Rigid`,
        `ThermodynamicQuantities` performs all sums
        :math:`\\sum_{i \\in \\mathrm{filter}}` over free particles and
        rigid body centers - ignoring constituent particles to avoid double
        counting.
    .---------------------------------------------------------------------------------

    Examples::

        f = filter.Type('A')
        compute.ThermodynamicQuantities(filter=f)
    """

    def __init__(self, filter):
        super().__init__()
        self._filter = filter

    def _attach_hook(self):
        if isinstance(self._simulation.device, hoomd.device.CPU):
            spfbasic_cls = _sph.ComputeSPFNNBasicProperties
        else:
            spfbasic_cls = _sph.ComputeSPFNNBasicPropertiesGPU
        group = self._simulation.state._get_group(self._filter)
        self._cpp_obj = spfbasic_cls(self._simulation.state._cpp_sys_def, group)
        # super()._attach()

    @log(requires_run=True)
    def abs_velocity(self):
        """Absolute velocity (norm of the vector) of the subset """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.abs_velocity

    @log(requires_run=True)
    def num_particles(self):
        """Number of particles :math:`N` in the subset."""
        return self._cpp_obj.num_particles

    @log(requires_run=True)
    def volume(self):
        """Volume :math:`V` of the simulation box (area in 2D) \
        :math:`[\\mathrm{length}^{D}]`."""
        return self._cpp_obj.volume

    @log(requires_run=True)
    def fluid_vel_x_sum(self):
        """Sum of Fluid Particle velocity in xdir """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.fluid_vel_x_sum

    @log(requires_run=True)
    def fluid_vel_y_sum(self):
        """Sum of Fluid Particle velocity in xdir """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.fluid_vel_y_sum

    @log(requires_run=True)
    def fluid_vel_z_sum(self):
        """Sum of Fluid Particle velocity in xdir """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.fluid_vel_z_sum

    @log(requires_run=True)
    def mean_density(self):
        """Sum of Fluid Particle velocity in xdir """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.mean_density

    @log(requires_run=True)
    def mean_viscosity(self):
        """Sum of Fluid Particle velocity in xdir """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.mean_viscosity

    @log(requires_run=True)
    def max_viscosity(self):
        """Sum of Fluid Particle velocity in xdir """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.max_viscosity

    @log(requires_run=True)
    def max_shearrate(self):
        """Sum of Fluid Particle velocity in xdir """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.max_shearrate

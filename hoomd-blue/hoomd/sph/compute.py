"""
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de
"""

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
import inspect

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
    def e_kin_fluid(self):
        """Kinetic energy of the subset """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.e_kin_fluid

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



class SolidProperties(Compute):
    """Compute mechanical properties of a subset of the system.

    Args:
        filter (`hoomd.filter`): Particle filter to compute thermodynamic
            properties for.

    `SolidProperties` acts on a subset of particles in the system and
    calculates mechanical properties of those particles. Add a
    `SolidProperties` instance to a logger to save these quantities to a
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
            spfbasic_cls = _sph.ComputeSolidProperties
        else:
            spfbasic_cls = _sph.ComputeSolidPropertiesGPU
        group = self._simulation.state._get_group(self._filter)
        self._cpp_obj = spfbasic_cls(self._simulation.state._cpp_sys_def, group)
        # super()._attach()

    @log(requires_run=True)
    def total_drag_x(self):
        """Absolute velocity (norm of the vector) of the subset """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.total_drag_x

    @log(requires_run=True)
    def total_drag_y(self):
        """Absolute velocity (norm of the vector) of the subset """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.total_drag_y

    @log(requires_run=True)
    def total_drag_z(self):
        """Absolute velocity (norm of the vector) of the subset """
        self._cpp_obj.compute(self._simulation.timestep)
        return self._cpp_obj.total_drag_z

    @log(requires_run=True)
    def num_particles(self):
        """Number of particles :math:`N` in the subset."""
        return self._cpp_obj.num_particles

    @log(requires_run=True)
    def volume(self):
        """Volume :math:`V` of the simulation box (area in 2D) \
        :math:`[\\mathrm{length}^{D}]`."""
        return self._cpp_obj.volume

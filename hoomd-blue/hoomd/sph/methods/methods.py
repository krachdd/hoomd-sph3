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

SPH integration methods.

"""

from hoomd.sph import _sph
import hoomd
from hoomd.operation import AutotunedObject
from hoomd.data.parameterdicts import ParameterDict, TypeParameterDict
from hoomd.data.typeparam import TypeParameter
from hoomd.data.typeconverter import OnlyTypes, OnlyIf, to_type_converter
from hoomd.filter import ParticleFilter
from hoomd.variant import Variant
from collections.abc import Sequence
import inspect


class Method(AutotunedObject):
    """Base class integration method.

    Provides common methods for all subclasses.

    Note:
        Users should use the subclasses and not instantiate `Method` directly.
    """

    __doc__ = (
        inspect.cleandoc(__doc__)
        + "\n"
        + inspect.cleandoc(AutotunedObject._doc_inherited)
    )


    def _attach_hook(self):
        self._simulation.state.update_group_dof()

    def _detach_hook(self):
        if self._simulation is not None:
            self._simulation.state.update_group_dof()



class VelocityVerlet(Method):
    r"""Velocity-Verlet SPH integration method.

    Args:
        filter (hoomd.filter.filter_like): Subset of particles on which to
            apply this method.
        densitymethod (str): Density computation strategy,
            ``'SUMMATION'`` or ``'CONTINUITY'``.

    Integrates translational degrees of freedom using the Velocity-Verlet
    (leapfrog) scheme:

    .. math::

        \mathbf{v}_i^{n+1/2} &= \mathbf{v}_i^n
                               + \frac{\Delta t}{2} \mathbf{a}_i^n \\
        \mathbf{r}_i^{n+1}   &= \mathbf{r}_i^n
                               + \Delta t\, \mathbf{v}_i^{n+1/2} \\
        \mathbf{v}_i^{n+1}   &= \mathbf{v}_i^{n+1/2}
                               + \frac{\Delta t}{2} \mathbf{a}_i^{n+1}

    Attributes:
        filter (hoomd.filter.filter_like): Subset of particles on which to
            apply this method.
    """

    DENSITYMETHODS = {'SUMMATION':_sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      'CONTINUITY':_sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {'HARMONICAVERAGE':_sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    def __init__(self, filter, densitymethod):
        # store metadata
        param_dict = ParameterDict(filter=ParticleFilter,)
        param_dict.update(dict(filter=filter, densitymethod=densitymethod))

        # set defaults
        self._param_dict.update(param_dict)

        self.str_densitymethod = self._param_dict._dict["densitymethod"]

        if self.str_densitymethod == str('SUMMATION'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str('CONTINUITY'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

    def _attach_hook(self):
        sim = self._simulation
        # initialize the reflected c++ class
        if isinstance(sim.device, hoomd.device.CPU) or not hasattr(_sph, 'VelocityVerletGPU'):
            self._cpp_obj = _sph.VelocityVerlet(sim.state._cpp_sys_def,
                                           sim.state._get_group(self.filter))
        else:
            self._cpp_obj = _sph.VelocityVerletGPU(sim.state._cpp_sys_def,
                                              sim.state._get_group(self.filter))

        # Reload density and viscosity methods from __dict__
        self.str_densitymethod = self._param_dict._dict["densitymethod"]
        if self.str_densitymethod == str('SUMMATION'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str('CONTINUITY'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        self.setdensitymethod(self.str_densitymethod)

        # Attach param_dict and typeparam_dict
        super()._attach_hook()

    # @property
    def densitymethod(self):
        # Invert key mapping
        invD = dict((v,k) for k, v in self.DENSITYMETHODS.items())
        return invD[self._cpp_obj.getDensityMethod()]

    # @densitymethod.setter
    def setdensitymethod(self, method):
        if method not in self.DENSITYMETHODS:
            raise ValueError("Undefined DensityMethod.")
        self._cpp_obj.setDensityMethod(self.DENSITYMETHODS[method])

        



class VelocityVerletBasic(Method):
    r"""Basic Velocity-Verlet SPH integration method (single-phase only).

    Args:
        filter (hoomd.filter.filter_like): Subset of particles on which to
            apply this method.
        densitymethod (str): Density computation strategy,
            ``'SUMMATION'`` or ``'CONTINUITY'``.

    Simplified Velocity-Verlet integrator without transport-velocity
    corrections; equivalent to `VelocityVerlet` but without the TV
    half-step for the position update.

    Attributes:
        filter (hoomd.filter.filter_like): Subset of particles on which to
            apply this method.
    """

    DENSITYMETHODS = {'SUMMATION':_sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      'CONTINUITY':_sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {'HARMONICAVERAGE':_sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    def __init__(self, filter, densitymethod):
        # store metadata
        param_dict = ParameterDict(filter=ParticleFilter,)
        param_dict.update(dict(filter=filter, densitymethod=densitymethod))

        # set defaults
        self._param_dict.update(param_dict)

        self.str_densitymethod = self._param_dict._dict["densitymethod"]

        if self.str_densitymethod == str('SUMMATION'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str('CONTINUITY'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

    def _attach_hook(self):
        sim = self._simulation
        # initialize the reflected c++ class
        if isinstance(sim.device, hoomd.device.CPU) or not hasattr(_sph, 'VelocityVerletBasicGPU'):
            self._cpp_obj = _sph.VelocityVerletBasic(sim.state._cpp_sys_def,
                                           sim.state._get_group(self.filter))
        else:
            self._cpp_obj = _sph.VelocityVerletBasicGPU(sim.state._cpp_sys_def,
                                              sim.state._get_group(self.filter))

        # Reload density and viscosity methods from __dict__
        self.str_densitymethod = self._param_dict._dict["densitymethod"]
        if self.str_densitymethod == str('SUMMATION'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str('CONTINUITY'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        self.setdensitymethod(self.str_densitymethod)

        # Attach param_dict and typeparam_dict
        super()._attach_hook()

    # @property
    def densitymethod(self):
        # Invert key mapping
        invD = dict((v,k) for k, v in self.DENSITYMETHODS.items())
        return invD[self._cpp_obj.getDensityMethod()]

    # @densitymethod.setter
    def setdensitymethod(self, method):
        if method not in self.DENSITYMETHODS:
            raise ValueError("Undefined DensityMethod.")
        self._cpp_obj.setDensityMethod(self.DENSITYMETHODS[method])

        

class KickDriftKickTV(Method):
    r"""Kick-Drift-Kick integration method with transport-velocity (TV) correction.

    Args:
        filter (hoomd.filter.filter_like): Subset of particles on which to
            apply this method.
        densitymethod (str): Density computation strategy,
            ``'SUMMATION'`` or ``'CONTINUITY'``.
        vlimit (bool): Enable velocity limiter.  Default ``False``.
        vlimit_val (float): Maximum allowed velocity magnitude.
        xlimit (bool): Enable position limiter.  Default ``False``.
        xlimit_val (float): Maximum allowed position displacement per step.

    Implements the transport-velocity formulation of Adami et al. (2013).
    Particle positions are advanced along the smooth transport velocity
    :math:`\tilde{\mathbf{v}}_i` rather than the physical velocity
    :math:`\mathbf{v}_i`:

    .. math::

        \tilde{\mathbf{v}}_i^{n+1/2} &= \mathbf{v}_i^n
                                       + \frac{\Delta t}{2}
                                         \mathbf{a}_i^n \\
        \mathbf{r}_i^{n+1}           &= \mathbf{r}_i^n
                                       + \Delta t\,
                                         \tilde{\mathbf{v}}_i^{n+1/2} \\
        \mathbf{v}_i^{n+1}           &= \mathbf{v}_i^{n+1/2}
                                       + \frac{\Delta t}{2}
                                         \mathbf{a}_i^{n+1}

    Attributes:
        filter (hoomd.filter.filter_like): Subset of particles on which to
            apply this method.
    """

    DENSITYMETHODS = {'SUMMATION':_sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      'CONTINUITY':_sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {'HARMONICAVERAGE':_sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    def __init__(self, filter, densitymethod, vlimit = False, vlimit_val = 0.0, xlimit = False, xlimit_val = 0.0):
        # store metadata
        param_dict = ParameterDict(filter=ParticleFilter,)
        param_dict.update(dict(filter=filter, densitymethod=densitymethod, 
                               vlimit = vlimit, vlimit_val = vlimit_val,
                               xlimit = xlimit, xlimit_val = xlimit_val))

        # set defaults
        self._param_dict.update(param_dict)

        self.str_densitymethod = self._param_dict._dict["densitymethod"]
        self.mvlimit = self._param_dict._dict["vlimit"]
        self.mxlimit = self._param_dict._dict["xlimit"]
        self.mvlimit_val = self._param_dict._dict["vlimit_val"]
        self.mxlimit_val = self._param_dict._dict["xlimit_val"]

        if self.str_densitymethod == str('SUMMATION'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str('CONTINUITY'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

    def _attach_hook(self):
        sim = self._simulation
        # initialize the reflected c++ class
        if isinstance(sim.device, hoomd.device.CPU):
            self._cpp_obj = _sph.KickDriftKickTV(sim.state._cpp_sys_def,
                                           sim.state._get_group(self.filter))
        else:
            self._cpp_obj = _sph.KickDriftKickTVGPU(sim.state._cpp_sys_def,
                                              sim.state._get_group(self.filter))

        # Reload density and viscosity methods from __dict__
        self.str_densitymethod = self._param_dict._dict["densitymethod"]
        if self.str_densitymethod == str('SUMMATION'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str('CONTINUITY'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        self.setdensitymethod(self.str_densitymethod)

        if self.vlimit == True:
            self._cpp_obj.setvLimit(self.vlimit_val)

        if self.xlimit == True:
            self._cpp_obj.setxLimit(self.xlimit_val)

        # Attach param_dict and typeparam_dict
        super()._attach_hook()

    # @property
    def densitymethod(self):
        # Invert key mapping
        invD = dict((v,k) for k, v in self.DENSITYMETHODS.items())
        return invD[self._cpp_obj.getDensityMethod()]

    # @densitymethod.setter
    def setdensitymethod(self, method):
        if method not in self.DENSITYMETHODS:
            raise ValueError("Undefined DensityMethod.")
        self._cpp_obj.setDensityMethod(self.DENSITYMETHODS[method])

    # @property
    def getvlimit(self):
        return self._cpp_obj.getvLimit()

    # @densitymethod.setter
    def setvLimit(self, limit_val):
        if limit_val > 0:
            self._cpp_obj.setvLimit(limit_val)
        else:
            raise ValueError("vlimit_val must be positive.")

    # @property
    def getxlimit(self):
        return self._cpp_obj.getxLimit()

    # @densitymethod.setter
    def setxLimit(self, xlimit_val):
        if xlimit_val > 0:
            self._cpp_obj.setxLimit(xlimit_val)
        else:
            raise ValueError("xlimit_val must be positive.")



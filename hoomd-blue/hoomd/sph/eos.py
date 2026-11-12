# Copyright (c) 2025-2026 David Krach, Daniel Rostan.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software without
#    specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# maintainer: dkrach, david.krach@mib.uni-stuttgart.de

"""Equations of state.

An equation of state relates the particle density to the pressure in the weakly
compressible SPH formulation. Set the rest density and background pressure with
``set_params(rho0, backpressure)`` and pass the object to a model in
`hoomd.sph.sphmodel`. The models set the speed of sound
(see ``compute_speedofsound``).
"""

import hoomd 
import hoomd.sph
from hoomd import _hoomd 
from hoomd.sph import _sph
from hoomd.operation import _HOOMDBaseObject
import numpy

class _StateEquation(_HOOMDBaseObject):
    r"""
    Constructs the equation of state meta class
    """

    def __init__(self, name = None):
        self._in_context_manager = False

        self.SpeedOfSound = 0;
        self.BackgroundPressure = 0;
        self.TransportVelocityPressure = 0;
        self.RestDensity = 0;

        self.cpp_stateequation = None;

        # Allow kernel class to store a name.
        if name is None:
            self.name = "";
        else:
            self.name="_" + name;

    def check_initialization(self):
        # check that we have been initialized properly
        if self.cpp_stateequation is None:
            # hoomd.context.msg.error('Bug in hoomd_script: cpp_stateequation not set, please report\n');
            raise RuntimeError("Bug in hoomd_script: cpp_stateequation not set, please report\n");

    def set_params(self,rho0,bp,tvp = 0):
        self.check_initialization();
        self.RestDensity        = rho0.item() if isinstance(rho0, numpy.generic) else rho0
        self.BackgroundPressure = bp.item()   if isinstance(bp, numpy.generic)   else bp
        self.TransportVelocityPressure = tvp.item()   if isinstance(tvp, numpy.generic)   else tvp
        self.cpp_stateequation.setParams(self.RestDensity,0.1,self.BackgroundPressure,self.TransportVelocityPressure)

    def set_speedofsound(self,c):
        self.check_initialization();
        self.SpeedOfSound       = c.item()    if isinstance(c, numpy.generic)    else c
        self.cpp_stateequation.setParams(self.RestDensity,self.SpeedOfSound,self.BackgroundPressure,self.TransportVelocityPressure)

    def pressure(self, rho):
        self.check_initialization()
        mrho = rho.item() if isinstance(rho, numpy.generic) else rho
        return self.cpp_stateequation.Pressure(mrho)


class Tait(_StateEquation):
    R"""Tait (weakly-compressible) equation of state.

    .. math::

        p(\rho) = \frac{\rho_0 c_0^2}{\gamma}
                  \left[\left(\frac{\rho}{\rho_0}\right)^\gamma - 1\right]
                  + p_b

    where :math:`\rho_0` is the rest density, :math:`c_0` is the reference
    speed of sound, :math:`\gamma = 7` (water-like exponent), and :math:`p_b`
    is the background pressure.
    """
    def __init__(self):
        _StateEquation.__init__(self, "Tait")
        self.cpp_stateequation = _sph.Tait()

class Linear(_StateEquation):
    R"""Linear (acoustic) equation of state.

    .. math::

        p(\rho) = c_0^2 (\rho - \rho_0) + p_b

    where :math:`\rho_0` is the rest density, :math:`c_0` is the reference
    speed of sound, and :math:`p_b` is the background pressure.
    """
    def __init__(self):
        _StateEquation.__init__(self, "Linear")
        self.cpp_stateequation = _sph.Linear()

__all__ = [
    "Tait",
    "Linear",
]
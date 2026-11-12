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

"""SPH flow models.

An SPH model is a `hoomd.sph.force.Force` that evaluates the SPH density,
pressure, and momentum equations for the fluid particles and imposes boundary
conditions at the solid particles. Append one model to
`hoomd.sph.Integrator.forces`.

* `SinglePhaseFlow`: weakly compressible single-phase flow.
* `SinglePhaseFlowTV`: single-phase flow with the transport-velocity
  formulation.
* `SinglePhaseFlowGDGD`: single-phase flow with a transported scalar field that
  drives buoyancy (density-gradient-driven flow).
* `SinglePhaseFlowFS`: single-phase flow with a free surface.
* `TwoPhaseFlow`: two immiscible fluids with surface tension.
* `TwoPhaseFlowTV`: two-phase flow with the transport-velocity formulation.
"""

from .sphmodel import (SPHModel, SinglePhaseFlow, SinglePhaseFlowTV,
                        SinglePhaseFlowGDGD, SinglePhaseFlowFS,
                        TwoPhaseFlow, TwoPhaseFlowTV)

__all__ = [
    "SPHModel",
    "SinglePhaseFlow",
    "SinglePhaseFlowTV",
    "SinglePhaseFlowGDGD",
    "SinglePhaseFlowFS",
    "TwoPhaseFlow",
    "TwoPhaseFlowTV",
]

/* ---------------------------------------------------------
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
----------------------------------------------------------*/


#ifndef __EVALUATION_METHOD_DEFINITON_H__
#define __EVALUATION_METHOD_DEFINITON_H__

#include "hoomd/HOOMDMath.h"
/*! \file EvaluationMethodDefinition.h
    \brief Data structures to define Evaluation methods for density, viscosity etc
    */

namespace hoomd
    {
namespace sph
    {
    //! Enum for various density evaluation approaches
    enum DensityMethod
    {
        DENSITYSUMMATION,    //!< Summation approach
        DENSITYCONTINUITY,    //!< Continuity approach
    };

    //! Enum for various viscosity evaluation approaches
    enum ViscosityMethod
    {
        HARMONICAVERAGE, //!< Viscosity operator based on inter-particle averaged shear stress
    };

    //! Enum for various Colorgradient evaluation approaches
    enum ColorGradientMethod
    {
        DENSITYRATIO, //!< Method to compute TPF Color gradient
        NUMBERDENSITY, //!< Method to compute TPF Color gradient
    };

    } // end namespace sph
    } // end namespace hoomd

#endif

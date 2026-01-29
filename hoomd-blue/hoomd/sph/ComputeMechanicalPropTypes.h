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


#ifndef __COMPUTE_MECHANICALPROP_TYPES_H__
#define __COMPUTE_MECHANICALPROP_TYPES_H__

#include "hoomd/HOOMDMath.h"
/*! \file ComputeMechanicalPropTypes.h
    \brief Data structures common to both CPU and GPU implementations of ComputeMechanicalProp
    */

namespace hoomd
    {
namespace sph
    {
//! Enum for indexing the GPUArray of computed values
struct singlephaseflow_logger_index
    {
    //! The enum
    enum Enum
        {
        sum_fluid_velocity_x=0, //!< Index for the sum of fluid x-velocity in the GPUArray
        sum_fluid_velocity_y,   //!< Index for the sum of fluid y-velocity in the GPUArray
        sum_fluid_velocity_z,   //!< Index for the sum of fluid z-velocity in the GPUArray
        sum_fluid_density,      //!< Index for the mean fluid particle density
        abs_velocity,           //!< Index for the overall absolute velocity per particle
        e_kin_fluid,            //!< Index for the overall kinetic energy of the system
        num_quantities // final element to count number of quantities
        };
    };

struct solidphase_logger_index
    {
    //! The enum
    enum Enum
        {
        total_drag_x,           //!< Index for the overall total drag on the solid phase
        total_drag_y,           //!< Index for the overall total drag on the solid phase
        total_drag_z,           //!< Index for the overall total drag on the solid phase
        num_quantities // final element to count number of quantities
        };
    };

    } // end namespace sph
    } // end namespace hoomd

#endif
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

/*! \file VelocityVerletGPU.cuh
    \brief Declarations of GPU kernel wrapper functions for VelocityVerletGPU.
*/

#ifndef __VELOCITY_VERLET_GPU_CUH__
#define __VELOCITY_VERLET_GPU_CUH__

#include "hip/hip_runtime.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"

namespace hoomd
{
namespace sph
{
namespace kernel
{

/*! First half-step of velocity-Verlet on the GPU.
 *
 * Updates density, pressure, velocity, and position by half a time-step,
 * then wraps particle positions back into the periodic box.
 */
hipError_t gpu_vv_step1(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    Scalar4*              d_pos,
    Scalar4*              d_vel,
    const Scalar3*        d_accel,
    Scalar*               d_density,
    Scalar*               d_pressure,
    const Scalar3*        d_dpedt,
    int3*                 d_image,
    BoxDim                box,
    Scalar                deltaT,
    unsigned int          block_size);

/*! Second half-step of velocity-Verlet on the GPU.
 *
 * Computes acceleration from net force, updates dpedt from net rate-of-change,
 * completes the density/pressure/position/velocity integration.
 */
hipError_t gpu_vv_step2(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    Scalar4*              d_pos,
    Scalar4*              d_vel,
    Scalar3*              d_accel,
    Scalar*               d_density,
    Scalar*               d_pressure,
    Scalar3*              d_dpedt,
    const Scalar4*        d_net_force,
    const Scalar4*        d_net_ratedpe,
    Scalar                deltaT,
    unsigned int          block_size);

} // namespace kernel
} // namespace sph
} // namespace hoomd

#endif // __VELOCITY_VERLET_GPU_CUH__

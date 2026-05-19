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

/*! \file VelocityVerletGPU.cu
    \brief GPU kernels for velocity-Verlet SPH time integration.
*/

#include "VelocityVerletGPU.cuh"
#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"

namespace hoomd
{
namespace sph
{
namespace kernel
{

// =========================================================================
// Step 1 kernel: density/pressure/velocity/position half-step + box wrap
// =========================================================================

__global__ void gpu_vv_step1_kernel(
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
    Scalar                deltaT)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size)
        return;

    unsigned int j = d_index_array[group_idx];
    Scalar half_dt = Scalar(0.5) * deltaT;

    // density and pressure: first half-step
    d_density[j]  += half_dt * d_dpedt[j].x;
    d_pressure[j] += half_dt * d_dpedt[j].y;

    // velocity: first half-step
    Scalar3 acc = d_accel[j];
    Scalar4 vel = d_vel[j];
    vel.x += half_dt * acc.x;
    vel.y += half_dt * acc.y;
    vel.z += half_dt * acc.z;
    d_vel[j] = vel;

    // position: first half-step (using updated half-step velocity)
    Scalar4 pos = d_pos[j];
    pos.x += half_dt * vel.x;
    pos.y += half_dt * vel.y;
    pos.z += half_dt * vel.z;

    // wrap back into periodic box
    int3 img = d_image[j];
    box.wrap(pos, img);
    d_pos[j]   = pos;
    d_image[j] = img;
    }

// =========================================================================
// Step 2 kernel: acceleration update + density/pressure/position/velocity
// =========================================================================

__global__ void gpu_vv_step2_kernel(
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
    Scalar                deltaT)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size)
        return;

    unsigned int j = d_index_array[group_idx];
    Scalar half_dt = Scalar(0.5) * deltaT;

    // compute acceleration from net force (mass stored in vel.w)
    Scalar minv = Scalar(1.0) / d_vel[j].w;
    Scalar4 nf  = d_net_force[j];
    Scalar3 acc;
    acc.x = nf.x * minv;
    acc.y = nf.y * minv;
    acc.z = nf.z * minv;
    d_accel[j] = acc;

    // update dpedt from net rate-of-change (energy term ignored, matches CPU)
    Scalar4 nr = d_net_ratedpe[j];
    Scalar3 dpedt;
    dpedt.x = nr.x;
    dpedt.y = nr.y;
    dpedt.z = Scalar(0.0);
    d_dpedt[j] = dpedt;

    // density and pressure: second half-step
    d_density[j]  += half_dt * dpedt.x;
    d_pressure[j] += half_dt * dpedt.y;

    // position: second half-step (half-step velocity before vel update)
    Scalar4 vel = d_vel[j];
    Scalar4 pos = d_pos[j];
    pos.x += half_dt * vel.x;
    pos.y += half_dt * vel.y;
    pos.z += half_dt * vel.z;
    d_pos[j] = pos;

    // velocity: second half-step
    vel.x += half_dt * acc.x;
    vel.y += half_dt * acc.y;
    vel.z += half_dt * acc.z;
    d_vel[j] = vel;
    }

// =========================================================================
// Kernel wrappers
// =========================================================================

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
    unsigned int          block_size)
    {
    unsigned int grid = (group_size + block_size - 1) / block_size;
    hipLaunchKernelGGL(gpu_vv_step1_kernel,
                       dim3(grid), dim3(block_size), 0, 0,
                       group_size, d_index_array,
                       d_pos, d_vel, d_accel,
                       d_density, d_pressure, d_dpedt,
                       d_image, box, deltaT);
    return hipSuccess;
    }

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
    unsigned int          block_size)
    {
    unsigned int grid = (group_size + block_size - 1) / block_size;
    hipLaunchKernelGGL(gpu_vv_step2_kernel,
                       dim3(grid), dim3(block_size), 0, 0,
                       group_size, d_index_array,
                       d_pos, d_vel, d_accel,
                       d_density, d_pressure, d_dpedt,
                       d_net_force, d_net_ratedpe, deltaT);
    return hipSuccess;
    }

} // namespace kernel
} // namespace sph
} // namespace hoomd

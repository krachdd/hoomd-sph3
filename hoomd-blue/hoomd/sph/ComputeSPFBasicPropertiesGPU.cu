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
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*! \file ComputeSPFBasicPropertiesGPU.cu
    \brief GPU reduction for SinglePhaseFlow basic properties.

    Strategy
    --------
    One __global__ kernel maps each particle in the group onto a thread.
    Each thread computes a local SPFProps contribution.  Shared-memory
    reduction within each block collapses to a single per-block value, which
    is written to a device-side block-result buffer.  A second, small kernel
    (one block) then reduces the block-result buffer to the final single
    SPFProps written to d_out.

    This avoids the CUB dependency while still being fully on-device.
    Only the 6-double (48-byte) result crosses PCIe per log step instead
    of the full ~64 MB particle arrays.
*/

#include "ComputeSPFBasicPropertiesGPU.cuh"
#include "hoomd/HOOMDMath.h"

namespace hoomd
{
namespace sph
{
namespace kernel
{

// =========================================================================
// Pass 1: per-block partial reduction
// =========================================================================

__global__ void gpu_spf_partial_reduce_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    SPFProps*             d_block_out)
    {
    extern __shared__ double smem[];   // 6 doubles per thread

    unsigned int tid  = threadIdx.x;
    unsigned int gidx = blockIdx.x * blockDim.x + tid;

    // Each thread initialises its shared slot
    double* s_vx      = smem + 0 * blockDim.x;
    double* s_vy      = smem + 1 * blockDim.x;
    double* s_vz      = smem + 2 * blockDim.x;
    double* s_abs_v   = smem + 3 * blockDim.x;
    double* s_density = smem + 4 * blockDim.x;
    double* s_ekin    = smem + 5 * blockDim.x;

    double vx = 0.0, vy = 0.0, vz = 0.0, abs_v = 0.0, dens = 0.0, ekin = 0.0;

    if (gidx < group_size)
        {
        unsigned int i = d_index_array[gidx];
        Scalar4 vel    = d_vel[i];          // .w = mass
        double  mass   = (double)vel.w;
        vx   = (double)vel.x;
        vy   = (double)vel.y;
        vz   = (double)vel.z;
        abs_v  = sqrt(vx*vx + vy*vy + vz*vz);
        dens   = (double)d_density[i];
        ekin   = 0.5 * mass * abs_v * abs_v;
        }

    s_vx[tid]      = vx;
    s_vy[tid]      = vy;
    s_vz[tid]      = vz;
    s_abs_v[tid]   = abs_v;
    s_density[tid] = dens;
    s_ekin[tid]    = ekin;
    __syncthreads();

    // Tree reduction in shared memory
    for (unsigned int stride = blockDim.x >> 1; stride > 0; stride >>= 1)
        {
        if (tid < stride)
            {
            s_vx[tid]      += s_vx[tid + stride];
            s_vy[tid]      += s_vy[tid + stride];
            s_vz[tid]      += s_vz[tid + stride];
            s_abs_v[tid]   += s_abs_v[tid + stride];
            s_density[tid] += s_density[tid + stride];
            s_ekin[tid]    += s_ekin[tid + stride];
            }
        __syncthreads();
        }

    if (tid == 0)
        {
        d_block_out[blockIdx.x].vx      = s_vx[0];
        d_block_out[blockIdx.x].vy      = s_vy[0];
        d_block_out[blockIdx.x].vz      = s_vz[0];
        d_block_out[blockIdx.x].abs_v   = s_abs_v[0];
        d_block_out[blockIdx.x].density = s_density[0];
        d_block_out[blockIdx.x].e_kin   = s_ekin[0];
        }
    }

// =========================================================================
// Pass 2: reduce the per-block partial sums.
// Fixed block size of 256 (smem = 6*256*8 = 12 KB, always fits).
// Grid-stride loop accumulates any number of pass-1 blocks.
// =========================================================================

#define PASS2_BS 256

__global__ void gpu_spf_final_reduce_kernel(
    unsigned int  n_blocks,
    const SPFProps* d_block_in,
    SPFProps*       d_out)
    {
    extern __shared__ double smem2[];

    unsigned int tid = threadIdx.x;

    double* s_vx      = smem2 + 0 * PASS2_BS;
    double* s_vy      = smem2 + 1 * PASS2_BS;
    double* s_vz      = smem2 + 2 * PASS2_BS;
    double* s_abs_v   = smem2 + 3 * PASS2_BS;
    double* s_density = smem2 + 4 * PASS2_BS;
    double* s_ekin    = smem2 + 5 * PASS2_BS;

    double vx = 0.0, vy = 0.0, vz = 0.0, abs_v = 0.0, dens = 0.0, ekin = 0.0;

    for (unsigned int i = tid; i < n_blocks; i += PASS2_BS)
        {
        vx   += d_block_in[i].vx;
        vy   += d_block_in[i].vy;
        vz   += d_block_in[i].vz;
        abs_v  += d_block_in[i].abs_v;
        dens   += d_block_in[i].density;
        ekin   += d_block_in[i].e_kin;
        }

    s_vx[tid]      = vx;
    s_vy[tid]      = vy;
    s_vz[tid]      = vz;
    s_abs_v[tid]   = abs_v;
    s_density[tid] = dens;
    s_ekin[tid]    = ekin;
    __syncthreads();

    for (unsigned int stride = PASS2_BS >> 1; stride > 0; stride >>= 1)
        {
        if (tid < stride)
            {
            s_vx[tid]      += s_vx[tid + stride];
            s_vy[tid]      += s_vy[tid + stride];
            s_vz[tid]      += s_vz[tid + stride];
            s_abs_v[tid]   += s_abs_v[tid + stride];
            s_density[tid] += s_density[tid + stride];
            s_ekin[tid]    += s_ekin[tid + stride];
            }
        __syncthreads();
        }

    if (tid == 0)
        {
        d_out->vx      = s_vx[0];
        d_out->vy      = s_vy[0];
        d_out->vz      = s_vz[0];
        d_out->abs_v   = s_abs_v[0];
        d_out->density = s_density[0];
        d_out->e_kin   = s_ekin[0];
        }
    }

// =========================================================================
// Host-side launcher
// =========================================================================

void gpu_spf_basic_props_reduce(
    unsigned int          group_size,
    unsigned int          n_blocks,
    const unsigned int*   d_index_array,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    SPFProps*             d_block_out,
    SPFProps*             d_out)
    {
    if (group_size == 0) return;

    // Pass 1: 256 threads per block, smem = 6*256*8 = 12 KB.
    gpu_spf_partial_reduce_kernel<<<n_blocks, 256, 6 * 256 * sizeof(double)>>>(
        group_size, d_index_array, d_vel, d_density, d_block_out);

    // Pass 2: fixed PASS2_BS=256 threads with grid-stride loop, smem = 12 KB.
    gpu_spf_final_reduce_kernel<<<1, PASS2_BS, 6 * PASS2_BS * sizeof(double)>>>(
        n_blocks, d_block_out, d_out);
    }

} // namespace kernel
} // namespace sph
} // namespace hoomd

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

/*! \file ComputeSPFBasicPropertiesGPU.cuh
    \brief Device-side declarations for the GPU SPF basic-properties reduction.

    All six sums (vx, vy, vz, |v|, density, e_kin) are computed with a single
    CUB DeviceReduce pass over the particle-group index array, writing results
    into a small device-side output buffer.  No device->host transfer of particle
    arrays is required; only the 6-scalar result buffer crosses PCIe.
*/

#pragma once

#include "hoomd/HOOMDMath.h"

namespace hoomd
{
namespace sph
{
namespace kernel
{

//! POD carrying the six accumulated sums from one reduce pass.
struct SPFProps
    {
    double vx;       //!< sum of velocity x-components
    double vy;       //!< sum of velocity y-components
    double vz;       //!< sum of velocity z-components
    double abs_v;    //!< sum of |v| per particle
    double density;  //!< sum of densities
    double e_kin;    //!< sum of 0.5*m*|v|^2
    };

/*! Launch the two-pass block reduction.
    \param group_size    Number of particles in the group
    \param n_blocks      Number of pass-1 blocks (= ceil(group_size/256))
    \param d_index_array Device pointer to the group index array (size >= group_size)
    \param d_vel         Device pointer to particle velocities (Scalar4; .w = mass)
    \param d_density     Device pointer to particle densities (Scalar)
    \param d_block_out   Caller-allocated device workspace (size >= n_blocks SPFProps)
    \param d_out         Device pointer to output SPFProps (exactly 1 element)
*/
void gpu_spf_basic_props_reduce(
    unsigned int          group_size,
    unsigned int          n_blocks,
    const unsigned int*   d_index_array,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    SPFProps*             d_block_out,
    SPFProps*             d_out);

} // namespace kernel
} // namespace sph
} // namespace hoomd

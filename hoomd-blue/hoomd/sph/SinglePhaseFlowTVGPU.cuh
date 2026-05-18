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

/*! \file SinglePhaseFlowTVGPU.cuh
    \brief Declares GPU kernel wrapper for SinglePhaseFlowTVGPU.

    Only the TV-specific force computation kernel is declared here.  The shared
    ndensity / pressure / noslip / solid-forces kernels are in SinglePhaseFlowGPU.cuh
    and reused without modification.
*/

#ifndef __SINGLE_PHASE_FLOW_TV_GPU_CUH__
#define __SINGLE_PHASE_FLOW_TV_GPU_CUH__

#include "hip/hip_runtime.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"

#include "SPHDeviceFunctions.cuh"

namespace hoomd
{
namespace sph
{
namespace kernel
{

/*! Compute transport-velocity SPH forces on fluid particles.
 *
 *  Extension of gpu_sph_forcecomputation that adds:
 *    - Per-particle artificial-stress tensor A_ij (Adami 2013)
 *    - Background pressure contribution written to d_bpc (aux2)
 *  Both use the transport velocity d_tv (aux3).
 *
 * \param d_tv          Transport velocity array (aux3)
 * \param d_bpc         [out] Background pressure contribution to TV (aux2)
 * \param Pb            Transport-velocity background pressure scalar
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_tv_forcecomputation(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar3*        d_tv,
    const Scalar*         d_h,
    Scalar4*              d_force,
    Scalar4*              d_ratedpe,
    Scalar3*              d_bpc,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    uint32_t*             d_max_vel_bits,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHEOSDevParams       eos,
    SPHNNViscParams       nn,
    Scalar                Pb,
    int                   density_method,
    int                   artificial_viscosity,
    Scalar                avalpha,
    Scalar                avbeta,
    int                   tensil_correction,
    Scalar                tensil_eps_pos,
    Scalar                tensil_eps_neg,
    int                   density_diffusion,
    Scalar                ddiff,
    unsigned int          block_size);

} // namespace kernel
} // namespace sph
} // namespace hoomd

#endif // __SINGLE_PHASE_FLOW_TV_GPU_CUH__

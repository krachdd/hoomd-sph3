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

/*! \file TwoPhaseFlowCSFGPU.cuh
    \brief GPU kernel wrappers for the two-phase interface machinery:
           colour gradients (normals), normal smoothing, prescribed-contact-angle
           wall blend, curvature-form CSF surface force and the optional
           pairwise wall-adhesion force.

    Mirrors TwoPhaseFlow::compute_colorgradients() and
    TwoPhaseFlow::compute_surfaceforce() in TwoPhaseFlow.cc pass by pass.
    Only the smoothing-kernel type is a template parameter (the EOS does not
    enter), so 5 instantiations cover all two-phase solvers.
*/

#ifndef __TWO_PHASE_FLOW_CSF_GPU_CUH__
#define __TWO_PHASE_FLOW_CSF_GPU_CUH__

#include "hip/hip_runtime.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"

#include "SPHDeviceFunctions.cuh"

namespace hoomd
{
namespace sph
{

/*! POD parameters of the interface machinery. */
struct SPHCSFParams
    {
    int    cg_method;      //!< 0 = DENSITYRATIO, 1 = NUMBERDENSITY
    int    blend_active;   //!< 1 if the prescribed-contact-angle wall blend is applied (omega != 90 && sigma12 > 0)
    Scalar cos_omega;      //!< cos(omega)
    Scalar sin_omega;      //!< sin(omega)
    Scalar sigma12;        //!< fluid-fluid surface tension
    int    adhesion_active;//!< 1 if the pairwise wall-adhesion force is applied
    Scalar beta_adh;       //!< wall-adhesion calibration coefficient (SPH_BETA_ADH)
    };

namespace kernel
{

//! Raw colour gradients for ALL local particles (aux2 = solid normal, aux3 = fluid normal)
template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_colorgradient_raw(
    unsigned int          N_local,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    Scalar3*              d_sn,
    Scalar3*              d_fn,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHCSFParams          cp,
    unsigned int          block_size);

//! Shepard-smoothed fluid normals of the fluid group -> d_fn_smooth (Adami 2010)
template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_normal_smooth(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    const Scalar3*        d_fn,
    Scalar3*              d_fn_smooth,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    unsigned int          block_size);

//! Scatter d_src[i] -> d_dst[i] for the members of a group
hipError_t gpu_sph_2pf_copy_group(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar3*        d_src,
    Scalar3*              d_dst,
    unsigned int          block_size);

//! Blended prescribed-contact-angle wall correction of the fluid normals (in place)
template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_wall_blend(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    const Scalar3*        d_sn,
    Scalar3*              d_fn,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHCSFParams          cp,
    unsigned int          block_size);

//! CSF pass 1: unit normals + reliability flag for all local+ghost slots
hipError_t gpu_sph_2pf_csf_pass1(
    unsigned int          N_total,
    const Scalar4*        d_pos,
    const Scalar*         d_h,
    const Scalar3*        d_fn,
    Scalar3*              d_nhat,
    unsigned int*         d_rel,
    const unsigned int*   d_type_property_map,
    SPHKernelDevParams    kp,
    unsigned int          block_size);

//! CSF pass 2: Morris-corrected curvature and surface force (aux4) for the fluid group
template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_csf_pass2(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    const Scalar3*        d_fn,
    const Scalar3*        d_nhat,
    const unsigned int*   d_rel,
    Scalar3*              d_sf,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHCSFParams          cp,
    unsigned int          block_size);

//! Optional pairwise wall-adhesion force added to aux4 (SPH_BETA_ADH > 0)
template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_wall_adhesion(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    Scalar3*              d_sf,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHCSFParams          cp,
    unsigned int          block_size);

} // namespace kernel
} // namespace sph
} // namespace hoomd

#endif // __TWO_PHASE_FLOW_CSF_GPU_CUH__

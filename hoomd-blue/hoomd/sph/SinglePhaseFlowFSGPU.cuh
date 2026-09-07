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

/*! \file SinglePhaseFlowFSGPU.cuh
    \brief Declarations of GPU kernel wrapper functions for SinglePhaseFlowFSGPU.

    Four kernels mirror the four FS-specific CPU stages in SinglePhaseFlowFS:
      1. gpu_sph_fs_detect_freesurface  — kernel-completeness + outward normal
      2. gpu_sph_fs_compute_curvature   — mean curvature from normal divergence
      3. gpu_sph_fs_pressure_clamp      — clamp tensile surface pressure to 0
      4. gpu_sph_fs_forcecomputation    — TV forces + CSF + Young's wetting force

    Shared ndensity / pressure / noslip / solid-forces kernels are reused from
    SinglePhaseFlowGPU.cuh without modification.
*/

#ifndef __SINGLE_PHASE_FLOW_FS_GPU_CUH__
#define __SINGLE_PHASE_FLOW_FS_GPU_CUH__

#include "hip/hip_runtime.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"

#include "SPHDeviceFunctions.cuh"

namespace hoomd
{
namespace sph
{

// POD parameter struct for free-surface physics — usable on host and device.
struct SPHFSParams
    {
    Scalar sigma;         //!< Surface tension coefficient [N/m]
    Scalar fs_threshold;  //!< lambda threshold for surface-particle detection
    Scalar rho0_ref;      //!< Reference density for Shepard-sum volumes
    Scalar cos_ca;        //!< cos(contact_angle)
    Scalar sin_ca;        //!< sin(contact_angle)
    int    apply_ca;      //!< 1 if |contact_angle - pi/2| > 0.01 (contact-angle correction active)
    };

namespace kernel
{

/*! Compute kernel-completeness lambda and outward free-surface normal.
 *
 *  Writes:
 *    d_aux4[i].x <- lambda   (Shepard sum; fluid neighbours + self contribution)
 *    d_aux4[i].y <- 0        (kappa filled later by gpu_sph_fs_compute_curvature)
 *    d_aux4[i].z <- gnorm    (|grad_lambda|; 0 for bulk particles)
 *    d_fs_n[i]   <- n_fs     (unit outward normal; {0,0,0} for bulk particles)
 *
 *  Contact-angle correction (Huber 2016) is applied when fs.apply_ca != 0.
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_fs_detect_freesurface(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    Scalar3*              d_fs_n,       //!< [out] aux2: outward fs normals
    Scalar3*              d_aux4,       //!< [out] aux4: lambda / kappa / gnorm
    const Scalar*         d_h,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHFSParams           fs,
    unsigned int          block_size);

/*! Compute mean curvature kappa_i for surface particles.
 *
 *  Writes d_aux4[i].y <- kappa_i.  Skipped (set to 0) for bulk particles and
 *  contact-line particles (has solid neighbour).
 *
 *  \pre  gpu_sph_fs_detect_freesurface has been called.  In MPI runs,
 *        ghost-particle normals (aux2) must be synced before calling this.
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_fs_compute_curvature(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar3*        d_fs_n,       //!< [in]  aux2: outward normals
    Scalar3*              d_aux4,       //!< [r/w] aux4: kappa written to .y
    const Scalar*         d_h,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    Scalar                fs_threshold,
    unsigned int          block_size);

/*! Clamp tensile pressure to zero at free-surface particles.
 *
 *  Non-template: no kernel or EOS dependency.
 */
hipError_t gpu_sph_fs_pressure_clamp(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    Scalar*               d_pressure,
    const Scalar3*        d_aux4,       //!< [in] aux4: lambda in .x
    Scalar                fs_threshold,
    unsigned int          block_size);

/*! TV pair forces + CSF surface tension + Young's wetting force.
 *
 *  On entry d_bpc (aux2) holds the outward fs normals written by
 *  detect_freesurface.  Each thread reads its own normal, then accumulates
 *  the background pressure contribution (BPC) into the same slot.
 *
 * \param d_aux4    [in]  lambda(.x) / kappa(.y) / gnorm(.z) per particle
 * \param d_bpc     [r/w] aux2: read fs normal, overwrite with BPC
 * \param Pb        Transport-velocity background pressure scalar
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_fs_forcecomputation(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar3*        d_tv,
    const Scalar3*        d_aux4,       //!< lambda / kappa / gnorm (read)
    const Scalar*         d_h,
    Scalar4*              d_force,
    Scalar4*              d_ratedpe,
    Scalar3*              d_bpc,        //!< [r/w] aux2: in=fs_normal, out=bpc
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
    SPHFSParams           fs,
    int                   density_method,
    int                   artificial_viscosity,
    Scalar                avalpha,
    Scalar                avbeta,
    int                   density_diffusion,
    Scalar                ddiff,
    unsigned int          block_size);

} // namespace kernel
} // namespace sph
} // namespace hoomd

#endif // __SINGLE_PHASE_FLOW_FS_GPU_CUH__

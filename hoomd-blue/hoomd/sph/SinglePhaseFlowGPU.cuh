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

/*! \file SinglePhaseFlowGPU.cuh
    \brief Declarations of GPU kernel wrapper functions for SinglePhaseFlowGPU.
*/

#ifndef __SINGLE_PHASE_FLOW_GPU_CUH__
#define __SINGLE_PHASE_FLOW_GPU_CUH__

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

/*! Compute number density for fluid particles.
 *
 * \param group_size   Number of fluid particles
 * \param d_index_array Fluid group index → global particle index
 * \param d_pos        Position (xyz) + type (w)
 * \param d_vel        Velocity (xyz) + mass (w)
 * \param d_density    [out] Density array
 * \param d_h          Per-particle smoothing length (ignored when const_slength)
 * \param d_n_neigh    Number of neighbours per particle
 * \param d_nlist      Neighbour list (flat)
 * \param d_head_list  Head pointer per particle into d_nlist
 * \param box          Simulation box (for minImage)
 * \param kp           Kernel parameters
 * \param block_size   CUDA/HIP block size
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_ndensity(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    Scalar*               d_density,
    const Scalar*         d_h,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    BoxDim                box,
    SPHKernelDevParams    kp,
    unsigned int          block_size);

/*! Compute EOS pressure for fluid particles.
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_pressure(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar*         d_density,
    Scalar*               d_pressure,
    SPHEOSDevParams       eos,
    unsigned int          block_size);

/*! Compute fictitious solid-particle velocity and pressure (no-slip BC).
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_noslip(
    unsigned int          solid_group_size,
    const unsigned int*   d_solid_index,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    Scalar*               d_density,
    Scalar*               d_pressure,
    Scalar3*              d_vf,
    const Scalar3*        d_accel,
    const Scalar*         d_h,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHEOSDevParams       eos,
    Scalar3               bodyforce,
    unsigned int          block_size);

/*! Compute SPH forces on fluid particles.
 *
 * \param d_max_vel_bits  [in/out] Bit-cast uint for atomicMax max-velocity reduction.
 *                                 Caller must zero before the first call; read after last call.
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_forcecomputation(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar*         d_h,
    Scalar4*              d_force,
    Scalar4*              d_ratedpe,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    uint32_t*             d_max_vel_bits,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHEOSDevParams       eos,
    SPHNNViscParams       nn,
    int                   density_method,      // 0=SUMMATION 1=CONTINUITY
    int                   artificial_viscosity,// bool as int
    Scalar                avalpha,
    Scalar                avbeta,
    int                   tensil_correction,   // bool as int
    Scalar                tensil_eps_pos,
    Scalar                tensil_eps_neg,
    int                   density_diffusion,   // bool as int
    Scalar                ddiff,
    unsigned int          block_size);

/*! Compute SPH forces on solid particles (reaction force to fluid pressure/viscosity).
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_solid_forces(
    unsigned int          solid_group_size,
    const unsigned int*   d_solid_index,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar*         d_h,
    Scalar4*              d_force,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHNNViscParams       nn,
    int                   density_method,
    unsigned int          block_size);

} // namespace kernel
} // namespace sph
} // namespace hoomd

#endif // __SINGLE_PHASE_FLOW_GPU_CUH__

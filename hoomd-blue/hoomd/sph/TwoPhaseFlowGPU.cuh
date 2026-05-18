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

/*! \file TwoPhaseFlowGPU.cuh
    \brief Declarations of GPU kernel wrapper functions for TwoPhaseFlowGPU.

    Two new kernels handle the two-phase-specific physics:
      1. gpu_sph_2pf_forcecomputation  — pressure + viscosity + surface force
      2. gpu_sph_2pf_solid_forces      — solid particle reaction forces

    Shared ndensity / pressure / noslip kernels are reused from
    SinglePhaseFlowGPU.cuh with EOS1 parameters.
*/

#ifndef __TWO_PHASE_FLOW_GPU_CUH__
#define __TWO_PHASE_FLOW_GPU_CUH__

#include "hip/hip_runtime.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"

#include "SPHDeviceFunctions.cuh"

namespace hoomd
{
namespace sph
{

/*! POD parameter struct for two-phase force kernel parameters. */
struct SPHTwoPhaseParams
    {
    Scalar cmax;          //!< Maximum speed of sound (for Monaghan AV)
    Scalar avalpha;       //!< Monaghan AV: linear coefficient
    Scalar avbeta;        //!< Monaghan AV: quadratic coefficient
    Scalar riemann_beta;  //!< Riemann dissipation scaling coefficient
    Scalar ddiff;         //!< Molteni–Colagrossi density diffusion coefficient
    Scalar gvec_x;        //!< Gravity x-component (for CIP)
    Scalar gvec_y;        //!< Gravity y-component (for CIP)
    Scalar gvec_z;        //!< Gravity z-component (for CIP)
    int density_method;          //!< 0 = DENSITYSUMMATION, 1 = DENSITYCONTINUITY
    int artificial_viscosity;    //!< 1 if Monaghan AV is active
    int riemann_dissipation;     //!< 1 if Riemann dissipation is active
    int cip;                     //!< 1 if consistent interface pressure is active
    int density_diffusion;       //!< 1 if density diffusion is active
    };

namespace kernel
{

/*! Two-phase flow force computation kernel.
 *
 *  Computes pressure + viscous forces for each fluid particle using
 *  per-phase EOS and viscosity parameters.  Surface force density is
 *  read directly from aux4 (pre-computed on CPU by compute_surfaceforce).
 *
 *  Features:
 *   - Consistent interface pressure (Hu & Adams 2009)
 *   - Monaghan (1992) artificial viscosity OR Riemann dissipation (Zhang 2017)
 *   - Molteni–Colagrossi density diffusion (two-phase corrected)
 *   - Non-Newtonian viscosity for both phases
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
hipError_t gpu_sph_2pf_forcecomputation(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,          //!< aux1: fictitious solid velocities
    const Scalar3*        d_sf,          //!< aux4: surface force density (pre-computed)
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
    SPHEOSDevParams       eos1,
    SPHEOSDevParams       eos2,
    SPHNNViscParams       nn1,
    SPHNNViscParams       nn2,
    SPHTwoPhaseParams     fparams,
    unsigned int          block_size);

/*! Two-phase solid particle reaction forces.
 *
 *  Computes pressure and viscous forces on solid particles from
 *  neighbouring fluid particles.  Viscosity is selected per-neighbour
 *  based on whether it is fluid 1 or fluid 2.
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
hipError_t gpu_sph_2pf_solid_forces(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
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
    SPHNNViscParams       nn1,
    SPHNNViscParams       nn2,
    int                   density_method,
    unsigned int          block_size);

} // namespace kernel
} // namespace sph
} // namespace hoomd

#endif // __TWO_PHASE_FLOW_GPU_CUH__

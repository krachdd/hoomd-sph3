/* ---------------------------------------------------------
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.
Redistribution and use permitted under BSD 3-Clause License.
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*! \file SinglePhaseFlowGDGDGPU.cuh
    \brief Declares GPU kernel wrapper for SinglePhaseFlowGDGDGPU (GDGD force computation).
*/

#ifndef __SINGLE_PHASE_FLOW_GDGD_GPU_CUH__
#define __SINGLE_PHASE_FLOW_GDGD_GPU_CUH__

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

/*! Compute GDGD (gradient-density-gradient-driven) SPH forces.
 *
 *  Extends the base force kernel with:
 *    - Per-particle VRD rest density from aux4.x (scalar T)
 *    - On-the-fly VRD pressure (SUMMATION + VRD mode)
 *    - Scalar diffusion rate into rdpe.z
 *    - Boussinesq buoyancy correction
 *
 * \param d_aux4         Scalar field T (aux4, only .x component used)
 * \param kappa_s        Scalar diffusivity
 * \param beta_s         Expansion coefficient
 * \param scalar_ref     Reference scalar value
 * \param rho0           Global reference density
 * \param boussinesq     1 = Boussinesq mode, 0 = VRD mode
 * \param bodyforce      Body acceleration vector (for Boussinesq buoyancy)
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_gdgd_forcecomputation(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar3*        d_aux4,
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
    Scalar                kappa_s,
    Scalar                beta_s,
    Scalar                scalar_ref,
    Scalar                rho0,
    int                   boussinesq,
    Scalar3               bodyforce,
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

#endif // __SINGLE_PHASE_FLOW_GDGD_GPU_CUH__

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

/*! \file TwoPhaseFlowGPU.cu
    \brief GPU kernels for two-phase SPH force computation.

    Implements two GPU kernels:
      1. gpu_sph_2pf_forcecomputation — pressure + viscosity + surface force for fluid particles
         (one templated kernel; the "fast" wrapper fixes the viscosity models at compile time)
      2. gpu_sph_2pf_solid_forces     — exact reaction of the fluid pair forces on solid particles

    The pair physics lives in TwoPhaseFlowGPUPair.cuh, shared with the
    transport-velocity variant, and mirrors TwoPhaseFlow.cc term by term.
*/

#include "TwoPhaseFlowGPU.cuh"
#include "TwoPhaseFlowGPUPair.cuh"
#include "SPHDeviceFunctions.cuh"

namespace hoomd
{
namespace sph
{
namespace kernel
{

// =========================================================================
// Kernel 1: Two-phase force computation (fluid particles)
//
// Uses tpp=4 threads per particle: each group of 4 threads cooperates on
// one particle's neighbour loop (thread t processes neighbours t,t+4,t+8,...),
// then the partial force accumulations are reduced via warp shuffles.
// =========================================================================

static constexpr unsigned int SPH_2PF_TPP = 4; // threads per particle

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_,
         bool FAST, NonNewtonianModel NM1_, NonNewtonianModel NM2_>
__launch_bounds__(384, 1)
__global__ void gpu_sph_2pf_forcecomputation_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar3*        d_sf,
    const Scalar*         d_h,
    const Scalar*         d_gdot,
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
    SPHTwoPhaseParams     fparams)
    {
    constexpr unsigned int tpp = SPH_2PF_TPP;
    const unsigned int ppb  = blockDim.x / tpp;          // particles per block
    const unsigned int pid  = blockIdx.x * ppb + threadIdx.x / tpp;
    const unsigned int lane = threadIdx.x % tpp;          // thread index within particle
    const bool active       = (pid < group_size);

    unsigned int i = 0;
    SPH2PFParticle I;
    size_t       myHead = 0;
    unsigned int size   = 0;

    Scalar3 fi   = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    Scalar  drho = Scalar(0);

    if (active)
        {
        i = d_index_array[pid];
        I = sph_2pf_load_particle<FAST, NM1_, NM2_>(i, d_pos, d_vel, d_density, d_pressure,
                                                    d_h, d_gdot, d_type_property_map,
                                                    kp, eos1, eos2, nn1, nn2, fparams);
        myHead = d_head_list[i];
        size   = d_n_neigh[i];

        // Neighbour loop — thread 'lane' processes neighbours lane, lane+tpp, ...
        for (unsigned int j = lane; j < size; j += tpp)
            {
            unsigned int k = d_nlist[myHead + j];
            sph_2pf_pair<KT_, FAST, NM1_, NM2_>(I, k, d_pos, d_vel, d_density, d_pressure,
                                                d_vf, d_h, d_gdot, d_type_property_map,
                                                box, kp, eos1, eos2, nn1, nn2, fparams,
                                                fi, drho);
            }
        }

    // Warp-level reduction across the tpp threads of each particle group
    // (inactive threads contribute zeros).
    #pragma unroll
    for (int d = tpp / 2; d > 0; d >>= 1)
        {
        fi.x += __shfl_xor_sync(0xffffffff, fi.x, d);
        fi.y += __shfl_xor_sync(0xffffffff, fi.y, d);
        fi.z += __shfl_xor_sync(0xffffffff, fi.z, d);
        drho += __shfl_xor_sync(0xffffffff, drho, d);
        }

    if (active && lane == 0)
        {
        // Surface force density pre-computed in aux4 (CPU compute_surfaceforce)
        Scalar3 sf = d_sf[i];
        fi.x += sf.x;
        fi.y += sf.y;
        fi.z += sf.z;

        d_force[i]   = make_scalar4(fi.x, fi.y, fi.z, Scalar(0));
        // Only drho/dt is provided; pressure is re-evaluated from the EOS
        // every step in computeForces() (as on the CPU), so dp/dt stays 0.
        d_ratedpe[i] = make_scalar4(drho, Scalar(0), Scalar(0), Scalar(0));

        // Max-velocity reduction for the adaptive timestep
        float vi_total = float(sqrt(dot(I.vel, I.vel)));
        atomicMax(d_max_vel_bits, __float_as_uint(vi_total));
        }
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_,
         bool FAST, NonNewtonianModel NM1_, NonNewtonianModel NM2_>
static hipError_t launch_2pf_forcecomputation(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar3*        d_sf,
    const Scalar*         d_h,
    const Scalar*         d_gdot,
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
    unsigned int          block_size)
    {
    if (group_size == 0) return hipSuccess;
    constexpr unsigned int tpp = SPH_2PF_TPP;
    // Clamp block_size to the kernel's register-imposed maximum and to a
    // multiple of tpp (any warp multiple satisfies the latter).
    hipFuncAttributes attr;
    hipFuncGetAttributes(&attr,
        (const void*)(gpu_sph_2pf_forcecomputation_kernel<KT_, SET1_, SET2_, FAST, NM1_, NM2_>));
    block_size = min(block_size, (unsigned int)attr.maxThreadsPerBlock);
    block_size = (block_size / tpp) * tpp;
    if (block_size == 0) block_size = tpp;

    unsigned int ppb = block_size / tpp;
    dim3 grid((group_size + ppb - 1) / ppb, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_forcecomputation_kernel<KT_, SET1_, SET2_, FAST, NM1_, NM2_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array,
                       d_pos, d_vel, d_density, d_pressure, d_vf, d_sf, d_h, d_gdot,
                       d_force, d_ratedpe,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map,
                       d_max_vel_bits, box, kp, eos1, eos2, nn1, nn2, fparams);
    return hipSuccess;
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
hipError_t gpu_sph_2pf_forcecomputation(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar3*        d_sf,
    const Scalar*         d_h,
    const Scalar*         d_gdot,
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
    unsigned int          block_size)
    {
    return launch_2pf_forcecomputation<KT_, SET1_, SET2_, false, NEWTONIAN, NEWTONIAN>(
        group_size, d_index_array, d_pos, d_vel, d_density, d_pressure, d_vf, d_sf, d_h, d_gdot,
        d_force, d_ratedpe, d_n_neigh, d_nlist, d_head_list, d_type_property_map,
        d_max_vel_bits, box, kp, eos1, eos2, nn1, nn2, fparams, block_size);
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_,
         NonNewtonianModel NMI_, NonNewtonianModel NMJ_>
hipError_t gpu_sph_2pf_forcecomputation_fast(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar3*        d_sf,
    const Scalar*         d_h,
    const Scalar*         d_gdot,
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
    unsigned int          block_size)
    {
    return launch_2pf_forcecomputation<KT_, SET1_, SET2_, true, NMI_, NMJ_>(
        group_size, d_index_array, d_pos, d_vel, d_density, d_pressure, d_vf, d_sf, d_h, d_gdot,
        d_force, d_ratedpe, d_n_neigh, d_nlist, d_head_list, d_type_property_map,
        d_max_vel_bits, box, kp, eos1, eos2, nn1, nn2, fparams, block_size);
    }

// =========================================================================
// Kernel 2: Two-phase solid particle reaction forces
//
// Mirrors TwoPhaseFlow::compute_solid_forces(): for each solid particle,
// accumulate the exact reaction of the pair forces the fluid loop applied
// to its fluid neighbours.  The pair expressions are identical to the fluid
// loop (symmetric in i<->j with dx, dv flipping sign), so Newton's third
// law holds without extra sign flips or mass-ratio scaling.
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
__global__ void gpu_sph_2pf_solid_forces_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar*         d_h,
    const Scalar*         d_gdot,
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
    int                   nn_active)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;

    unsigned int i = d_index_array[group_idx];

    Scalar4 posi = d_pos[i];
    Scalar3 pi   = make_scalar3(posi.x, posi.y, posi.z);
    // Fictitious (Adami) velocity: the fluid loop used dv = v_f - v~_s
    Scalar3 vfi  = d_vf[i];
    Scalar3 vi   = make_scalar3(vfi.x, vfi.y, vfi.z);
    Scalar  mi   = d_vel[i].w;
    Scalar  Pi   = d_pressure[i];
    Scalar  rhoi = d_density[i];
    Scalar  Vi   = mi / rhoi;
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];

    Scalar3 fi = make_scalar3(Scalar(0), Scalar(0), Scalar(0));

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];

    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];

        Scalar4 posk = d_pos[k];
        if (sph_checksolid(d_type_property_map, posk.w)) continue;
        bool j_isfluid1 = sph_checkfluid1(d_type_property_map, posk.w);

        Scalar3 dx  = box.minImage(make_scalar3(pi.x - posk.x, pi.y - posk.y, pi.z - posk.z));
        Scalar  rsq = dot(dx, dx);
        if (kp.const_slength && rsq > kp.rcutsq) continue;

        Scalar4 velk = d_vel[k];
        Scalar3 vj   = make_scalar3(velk.x, velk.y, velk.z);
        Scalar  mj   = velk.w;
        Scalar  rhoj = d_density[k];
        Scalar  Vj   = mj / rhoj;
        Scalar  Pj   = d_pressure[k];

        Scalar3 dv = make_scalar3(vi.x - vj.x, vi.y - vj.y, vi.z - vj.z);
        Scalar  r  = sqrt(rsq);

        Scalar meanh  = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = (r > Scalar(1e-8) * meanh) ? dwdr / r : Scalar(0);

        // Pressure reaction: same expression as the fluid loop
        Scalar temp0;
        if (density_method == 0)
            temp0 = (Vi * Vi + Vj * Vj) * ((rhoj * Pi + rhoi * Pj) / (rhoi + rhoj));
        else
            temp0 = mi * mj * (Pi + Pj) / (rhoi * rhoj);
        fi.x -= temp0 * dwdr_r * dx.x;
        fi.y -= temp0 * dwdr_r * dx.y;
        fi.z -= temp0 * dwdr_r * dx.z;

        // Viscous reaction: viscosity of the fluid neighbour at its own
        // per-particle shear rate (solid pairs use mu_eff_j on both sides)
        Scalar gdot_j = nn_active ? d_gdot[k] : Scalar(0);
        Scalar mu_eff_j = sph_2pf_mu_eff<false, NEWTONIAN, NEWTONIAN>(j_isfluid1, gdot_j, nn1, nn2);
        temp0 = mu_eff_j * (Vi * Vi + Vj * Vj) * dwdr_r;
        fi.x += temp0 * dv.x;
        fi.y += temp0 * dv.y;
        fi.z += temp0 * dv.z;
        }

    // Accumulate into the global force (fluid forces were written by forcecomputation)
    d_force[i].x += fi.x;
    d_force[i].y += fi.y;
    d_force[i].z += fi.z;
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
hipError_t gpu_sph_2pf_solid_forces(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar*         d_h,
    const Scalar*         d_gdot,
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
    int                   nn_active,
    unsigned int          block_size)
    {
    if (group_size == 0) return hipSuccess;
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_solid_forces_kernel<KT_, SET1_, SET2_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array,
                       d_pos, d_vel, d_density, d_pressure, d_vf, d_h, d_gdot,
                       d_force,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map,
                       box, kp, nn1, nn2, density_method, nn_active);
    return hipSuccess;
    }

// =========================================================================
// Explicit template instantiations
// =========================================================================

#define INST_2PF_GPU(KT, SET1, SET2) \
    template hipError_t gpu_sph_2pf_forcecomputation<KT, SET1, SET2>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, const Scalar*, const Scalar*, \
        const Scalar3*, const Scalar3*, const Scalar*, const Scalar*, \
        Scalar4*, Scalar4*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, uint32_t*, BoxDim, \
        SPHKernelDevParams, SPHEOSDevParams, SPHEOSDevParams, \
        SPHNNViscParams, SPHNNViscParams, SPHTwoPhaseParams, unsigned int); \
    template hipError_t gpu_sph_2pf_solid_forces<KT, SET1, SET2>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, const Scalar*, const Scalar*, \
        const Scalar3*, const Scalar*, const Scalar*, \
        Scalar4*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, BoxDim, SPHKernelDevParams, \
        SPHNNViscParams, SPHNNViscParams, int, int, unsigned int);

INST_2PF_GPU(wendlandc2, tait,   tait)
INST_2PF_GPU(wendlandc2, tait,   linear)
INST_2PF_GPU(wendlandc2, linear, tait)
INST_2PF_GPU(wendlandc2, linear, linear)
INST_2PF_GPU(wendlandc4, tait,   tait)
INST_2PF_GPU(wendlandc4, tait,   linear)
INST_2PF_GPU(wendlandc4, linear, tait)
INST_2PF_GPU(wendlandc4, linear, linear)
INST_2PF_GPU(wendlandc6, tait,   tait)
INST_2PF_GPU(wendlandc6, tait,   linear)
INST_2PF_GPU(wendlandc6, linear, tait)
INST_2PF_GPU(wendlandc6, linear, linear)
INST_2PF_GPU(quintic,    tait,   tait)
INST_2PF_GPU(quintic,    tait,   linear)
INST_2PF_GPU(quintic,    linear, tait)
INST_2PF_GPU(quintic,    linear, linear)
INST_2PF_GPU(cubicspline,tait,   tait)
INST_2PF_GPU(cubicspline,tait,   linear)
INST_2PF_GPU(cubicspline,linear, tait)
INST_2PF_GPU(cubicspline,linear, linear)

#undef INST_2PF_GPU

// ── Fast-path (viscosity-model-templated) instantiations ─────────────────
// NEWTONIAN/POWERLAW combinations; other models take the runtime kernel.

#define INST_2PF_FAST(KT, SET1, SET2, NM1, NM2) \
    template hipError_t gpu_sph_2pf_forcecomputation_fast<KT, SET1, SET2, NM1, NM2>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, const Scalar*, const Scalar*, \
        const Scalar3*, const Scalar3*, const Scalar*, const Scalar*, \
        Scalar4*, Scalar4*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, uint32_t*, BoxDim, \
        SPHKernelDevParams, SPHEOSDevParams, SPHEOSDevParams, \
        SPHNNViscParams, SPHNNViscParams, SPHTwoPhaseParams, unsigned int);

#define INST_2PF_FAST_KT(KT, SET1, SET2) \
    INST_2PF_FAST(KT, SET1, SET2, NEWTONIAN, NEWTONIAN) \
    INST_2PF_FAST(KT, SET1, SET2, POWERLAW,  NEWTONIAN) \
    INST_2PF_FAST(KT, SET1, SET2, NEWTONIAN, POWERLAW ) \
    INST_2PF_FAST(KT, SET1, SET2, POWERLAW,  POWERLAW )

INST_2PF_FAST_KT(wendlandc2, tait,   tait)
INST_2PF_FAST_KT(wendlandc2, tait,   linear)
INST_2PF_FAST_KT(wendlandc2, linear, tait)
INST_2PF_FAST_KT(wendlandc2, linear, linear)
INST_2PF_FAST_KT(wendlandc4, tait,   tait)
INST_2PF_FAST_KT(wendlandc4, tait,   linear)
INST_2PF_FAST_KT(wendlandc4, linear, tait)
INST_2PF_FAST_KT(wendlandc4, linear, linear)
INST_2PF_FAST_KT(wendlandc6, tait,   tait)
INST_2PF_FAST_KT(wendlandc6, tait,   linear)
INST_2PF_FAST_KT(wendlandc6, linear, tait)
INST_2PF_FAST_KT(wendlandc6, linear, linear)
INST_2PF_FAST_KT(quintic,    tait,   tait)
INST_2PF_FAST_KT(quintic,    tait,   linear)
INST_2PF_FAST_KT(quintic,    linear, tait)
INST_2PF_FAST_KT(quintic,    linear, linear)
INST_2PF_FAST_KT(cubicspline,tait,   tait)
INST_2PF_FAST_KT(cubicspline,tait,   linear)
INST_2PF_FAST_KT(cubicspline,linear, tait)
INST_2PF_FAST_KT(cubicspline,linear, linear)

#undef INST_2PF_FAST_KT
#undef INST_2PF_FAST

} // namespace kernel
} // namespace sph
} // namespace hoomd

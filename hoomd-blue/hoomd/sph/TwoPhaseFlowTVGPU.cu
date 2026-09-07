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

/*! \file TwoPhaseFlowTVGPU.cu
    \brief GPU kernel for two-phase transport-velocity SPH force computation.

    Extends the shared two-phase pair interaction (TwoPhaseFlowGPUPair.cuh,
    mirroring TwoPhaseFlow::forcecomputation) with the transport-velocity
    terms of TwoPhaseFlowTV::forcecomputation:
      - Artificial-stress tensor (Adami 2013) for tensile instability suppression
      - Background-pressure contribution (BPC) written to aux2
    Both TV terms are applied to fluid-fluid pairs only (not fluid-solid).
*/

#include "TwoPhaseFlowTVGPU.cuh"
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
// Kernel: Two-phase transport-velocity force computation
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
__global__ void gpu_sph_2pf_tv_forcecomputation_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar3*        d_tv,
    const Scalar3*        d_sf,
    const Scalar*         d_h,
    const Scalar*         d_gdot,
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
    SPHEOSDevParams       eos1,
    SPHEOSDevParams       eos2,
    SPHNNViscParams       nn1,
    SPHNNViscParams       nn2,
    SPHTwoPhaseParams     fparams,
    Scalar                Pb1,
    Scalar                Pb2)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;

    unsigned int i = d_index_array[group_idx];

    SPH2PFParticle I = sph_2pf_load_particle<false, NEWTONIAN, NEWTONIAN>(
        i, d_pos, d_vel, d_density, d_pressure, d_h, d_gdot, d_type_property_map,
        kp, eos1, eos2, nn1, nn2, fparams);

    // Transport velocity and per-phase background pressure
    Scalar3 tvi = d_tv[i];
    Scalar  Pbi = I.isfluid1 ? Pb1 : Pb2;

    // Artificial stress tensor for particle i: A = rho_i * v_i (x) (tv_i - v_i)
    const Scalar3 vi = I.vel;
    Scalar A11i = I.rho * vi.x * (tvi.x - vi.x);
    Scalar A12i = I.rho * vi.x * (tvi.y - vi.y);
    Scalar A13i = I.rho * vi.x * (tvi.z - vi.z);
    Scalar A21i = I.rho * vi.y * (tvi.x - vi.x);
    Scalar A22i = I.rho * vi.y * (tvi.y - vi.y);
    Scalar A23i = I.rho * vi.y * (tvi.z - vi.z);
    Scalar A31i = I.rho * vi.z * (tvi.x - vi.x);
    Scalar A32i = I.rho * vi.z * (tvi.y - vi.y);
    Scalar A33i = I.rho * vi.z * (tvi.z - vi.z);

    Scalar3 fi   = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    Scalar  drho = Scalar(0);
    Scalar3 bpc  = make_scalar3(Scalar(0), Scalar(0), Scalar(0));

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];

    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];

        // Pressure + dissipation + viscous force and continuity rate
        SPH2PFPairOut o = sph_2pf_pair<KT_, false, NEWTONIAN, NEWTONIAN>(
            I, k, d_pos, d_vel, d_density, d_pressure, d_vf, d_h, d_gdot,
            d_type_property_map, box, kp, eos1, eos2, nn1, nn2, fparams, fi, drho);
        if (!o.valid || o.j_issolid) continue;

        // ── Transport-velocity terms (fluid-fluid pairs only) ─────────────
        Scalar3 tvk = d_tv[k];
        Scalar4 velk = d_vel[k];
        Scalar3 vk = make_scalar3(velk.x, velk.y, velk.z);
        Scalar A11k = o.rhoj * vk.x * (tvk.x - vk.x);
        Scalar A12k = o.rhoj * vk.x * (tvk.y - vk.y);
        Scalar A13k = o.rhoj * vk.x * (tvk.z - vk.z);
        Scalar A21k = o.rhoj * vk.y * (tvk.x - vk.x);
        Scalar A22k = o.rhoj * vk.y * (tvk.y - vk.y);
        Scalar A23k = o.rhoj * vk.y * (tvk.z - vk.z);
        Scalar A31k = o.rhoj * vk.z * (tvk.x - vk.x);
        Scalar A32k = o.rhoj * vk.z * (tvk.y - vk.y);
        Scalar A33k = o.rhoj * vk.z * (tvk.z - vk.z);

        Scalar tv_temp = Scalar(0.5) * o.vijsqr * o.dwdr_r;
        Scalar A1ij = (A11i+A11k)*o.dx.x + (A12i+A12k)*o.dx.y + (A13i+A13k)*o.dx.z;
        Scalar A2ij = (A21i+A21k)*o.dx.x + (A22i+A22k)*o.dx.y + (A23i+A23k)*o.dx.z;
        Scalar A3ij = (A31i+A31k)*o.dx.x + (A32i+A32k)*o.dx.y + (A33i+A33k)*o.dx.z;
        fi.x += tv_temp * A1ij;
        fi.y += tv_temp * A2ij;
        fi.z += tv_temp * A3ij;

        // Background-pressure contribution (aux2)
        Scalar bcoef = o.vijsqr * Pbi / I.m * o.dwdr_r;
        bpc.x -= bcoef * o.dx.x;
        bpc.y -= bcoef * o.dx.y;
        bpc.z -= bcoef * o.dx.z;
        } // end neighbour loop

    // Surface force (pre-computed on CPU by compute_surfaceforce)
    Scalar3 sf = d_sf[i];
    fi.x += sf.x;
    fi.y += sf.y;
    fi.z += sf.z;

    d_force[i]   = make_scalar4(fi.x, fi.y, fi.z, Scalar(0));
    // Only drho/dt is provided; pressure is re-evaluated from the EOS every
    // step in computeForces() (as on the CPU), so dp/dt stays 0.
    d_ratedpe[i] = make_scalar4(drho, Scalar(0), Scalar(0), Scalar(0));
    d_bpc[i]     = bpc;

    // Max-velocity reduction for the adaptive timestep
    float vi_total = float(sqrt(dot(vi, vi)));
    atomicMax(d_max_vel_bits, __float_as_uint(vi_total));
    }

// =========================================================================
// Wrapper function
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
hipError_t gpu_sph_2pf_tv_forcecomputation(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar3*        d_tv,
    const Scalar3*        d_sf,
    const Scalar*         d_h,
    const Scalar*         d_gdot,
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
    SPHEOSDevParams       eos1,
    SPHEOSDevParams       eos2,
    SPHNNViscParams       nn1,
    SPHNNViscParams       nn2,
    SPHTwoPhaseParams     fparams,
    Scalar                Pb1,
    Scalar                Pb2,
    unsigned int          block_size)
    {
    if (group_size == 0) return hipSuccess;
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_tv_forcecomputation_kernel<KT_, SET1_, SET2_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array,
                       d_pos, d_vel, d_density, d_pressure,
                       d_vf, d_tv, d_sf, d_h, d_gdot,
                       d_force, d_ratedpe, d_bpc,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map,
                       d_max_vel_bits, box, kp, eos1, eos2, nn1, nn2, fparams, Pb1, Pb2);
    return hipSuccess;
    }

// =========================================================================
// Explicit template instantiations (5 x 2 x 2 = 20)
// =========================================================================

#define INST_2PF_TV_GPU(KT, SET1, SET2) \
    template hipError_t gpu_sph_2pf_tv_forcecomputation<KT, SET1, SET2>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, const Scalar*, const Scalar*, \
        const Scalar3*, const Scalar3*, const Scalar3*, const Scalar*, const Scalar*, \
        Scalar4*, Scalar4*, Scalar3*, \
        const unsigned int*, const unsigned int*, const size_t*, const unsigned int*, \
        uint32_t*, BoxDim, SPHKernelDevParams, \
        SPHEOSDevParams, SPHEOSDevParams, SPHNNViscParams, SPHNNViscParams, \
        SPHTwoPhaseParams, Scalar, Scalar, unsigned int);

INST_2PF_TV_GPU(wendlandc2, tait,   tait)
INST_2PF_TV_GPU(wendlandc2, tait,   linear)
INST_2PF_TV_GPU(wendlandc2, linear, tait)
INST_2PF_TV_GPU(wendlandc2, linear, linear)
INST_2PF_TV_GPU(wendlandc4, tait,   tait)
INST_2PF_TV_GPU(wendlandc4, tait,   linear)
INST_2PF_TV_GPU(wendlandc4, linear, tait)
INST_2PF_TV_GPU(wendlandc4, linear, linear)
INST_2PF_TV_GPU(wendlandc6, tait,   tait)
INST_2PF_TV_GPU(wendlandc6, tait,   linear)
INST_2PF_TV_GPU(wendlandc6, linear, tait)
INST_2PF_TV_GPU(wendlandc6, linear, linear)
INST_2PF_TV_GPU(quintic,    tait,   tait)
INST_2PF_TV_GPU(quintic,    tait,   linear)
INST_2PF_TV_GPU(quintic,    linear, tait)
INST_2PF_TV_GPU(quintic,    linear, linear)
INST_2PF_TV_GPU(cubicspline,tait,   tait)
INST_2PF_TV_GPU(cubicspline,tait,   linear)
INST_2PF_TV_GPU(cubicspline,linear, tait)
INST_2PF_TV_GPU(cubicspline,linear, linear)

#undef INST_2PF_TV_GPU

} // namespace kernel
} // namespace sph
} // namespace hoomd

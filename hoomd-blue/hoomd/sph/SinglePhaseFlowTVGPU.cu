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

/*! \file SinglePhaseFlowTVGPU.cu
    \brief GPU kernel for transport-velocity SPH force computation.

    Extends the base SPH force kernel with:
      - Artificial stress tensor A_ij (Adami 2013 transport velocity correction)
      - Background pressure contribution written to d_bpc (aux2)
    Both terms use the transport velocity from d_tv (aux3).
*/

#include "SinglePhaseFlowTVGPU.cuh"
#include "SPHDeviceFunctions.cuh"

namespace hoomd
{
namespace sph
{
namespace kernel
{

template<SmoothingKernelType KT_, StateEquationType SET_>
__global__ void gpu_sph_tv_forcecomputation_kernel(
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
    Scalar                ddiff)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;

    unsigned int i = d_index_array[group_idx];

    Scalar3 pi   = make_scalar3(d_pos[i].x, d_pos[i].y, d_pos[i].z);
    Scalar3 vi   = make_scalar3(d_vel[i].x, d_vel[i].y, d_vel[i].z);
    Scalar  mi   = d_vel[i].w;
    Scalar  Pi   = d_pressure[i];
    Scalar  rhoi = d_density[i];
    Scalar  Vi   = mi / rhoi;
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];

    Scalar3 tvi  = d_tv[i]; // transport velocity of particle i

    // Artificial stress tensor components for particle i: A_ij = rho_i * v_i * (tv_i - v_i)
    Scalar A11i = rhoi * vi.x * (tvi.x - vi.x);
    Scalar A12i = rhoi * vi.x * (tvi.y - vi.y);
    Scalar A13i = rhoi * vi.x * (tvi.z - vi.z);
    Scalar A21i = rhoi * vi.y * (tvi.x - vi.x);
    Scalar A22i = rhoi * vi.y * (tvi.y - vi.y);
    Scalar A23i = rhoi * vi.y * (tvi.z - vi.z);
    Scalar A31i = rhoi * vi.z * (tvi.x - vi.x);
    Scalar A32i = rhoi * vi.z * (tvi.y - vi.y);
    Scalar A33i = rhoi * vi.z * (tvi.z - vi.z);

    Scalar4 fi   = make_scalar4(Scalar(0), Scalar(0), Scalar(0), Scalar(0));
    Scalar4 rdpe = make_scalar4(Scalar(0), Scalar(0), Scalar(0), Scalar(0));
    Scalar3 bpc  = make_scalar3(Scalar(0), Scalar(0), Scalar(0));

    // Pre-compute tensile correction reference value (constant h only)
    Scalar tensil_od_wdeltap = Scalar(0);
    if (tensil_correction && kp.const_slength)
        {
        Scalar wdeltap = sph_wij<KT_>(kp.alpha, kp.ch, kp.ch / Scalar(1.5));
        tensil_od_wdeltap = (wdeltap > Scalar(0)) ? (Scalar(1) / wdeltap) : Scalar(0);
        }

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];

    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];

        bool issolid = sph_checksolid(d_type_property_map, d_pos[k].w);

        Scalar3 pj  = make_scalar3(d_pos[k].x, d_pos[k].y, d_pos[k].z);
        Scalar3 dx  = box.minImage(make_scalar3(pi.x - pj.x, pi.y - pj.y, pi.z - pj.z));
        Scalar  rsq = dot(dx, dx);

        if (kp.const_slength && rsq > kp.rcutsq) continue;

        Scalar mj   = d_vel[k].w;
        Scalar rhoj = d_density[k];
        Scalar Vj   = mj / rhoj;
        Scalar Pj   = d_pressure[k];

        // Physical velocity of j (for fictitious solid: use vf; for fluid: vel)
        Scalar3 vj;
        if (issolid)
            vj = make_scalar3(d_vf[k].x, d_vf[k].y, d_vf[k].z);
        else
            vj = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);

        Scalar3 tvj = d_tv[k]; // transport velocity of neighbour k

        // Artificial stress tensor for particle k
        Scalar A11j = rhoj * vj.x * (tvj.x - vj.x);
        Scalar A12j = rhoj * vj.x * (tvj.y - vj.y);
        Scalar A13j = rhoj * vj.x * (tvj.z - vj.z);
        Scalar A21j = rhoj * vj.y * (tvj.x - vj.x);
        Scalar A22j = rhoj * vj.y * (tvj.y - vj.y);
        Scalar A23j = rhoj * vj.y * (tvj.z - vj.z);
        Scalar A31j = rhoj * vj.z * (tvj.x - vj.x);
        Scalar A32j = rhoj * vj.z * (tvj.y - vj.y);
        Scalar A33j = rhoj * vj.z * (tvj.z - vj.z);

        Scalar3 dv  = make_scalar3(vi.x - vj.x, vi.y - vj.y, vi.z - vj.z);
        Scalar  r   = sqrtf(rsq);

        Scalar meanh  = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar epssqr = Scalar(0.01) * meanh * meanh;

        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = dwdr / (r + epssqr);

        Scalar vijsqr = Vi * Vi + Vj * Vj;

        // TV pressure formula: (rho_j*P_i + rho_i*P_j) / (rho_i + rho_j)
        Scalar P_tv = (rhoj * Pi + rhoi * Pj) / (rhoi + rhoj);

        // Artificial viscosity (Monaghan 1983)
        Scalar avc = Scalar(0);
        if (artificial_viscosity && !issolid)
            {
            Scalar dotdvdx = dot(dv, dx);
            if (dotdvdx < Scalar(0))
                {
                Scalar muij    = meanh * dotdvdx / (rsq + epssqr);
                Scalar meanrho = Scalar(0.5) * (rhoi + rhoj);
                avc = (-avalpha * eos.c * muij + avbeta * muij * muij) / meanrho;
                }
            }

        // Pressure force
        fi.x -= vijsqr * (P_tv + avc) * dwdr_r * dx.x;
        fi.y -= vijsqr * (P_tv + avc) * dwdr_r * dx.y;
        fi.z -= vijsqr * (P_tv + avc) * dwdr_r * dx.z;

        // Tensile instability correction (Monaghan 1994, fluid-fluid only)
        if (tensil_correction && !issolid)
            {
            Scalar od_wdeltap;
            if (kp.const_slength)
                od_wdeltap = tensil_od_wdeltap;
            else
                {
                Scalar wdeltap = sph_wij<KT_>(kp.alpha, meanh, meanh / Scalar(1.5));
                od_wdeltap = (wdeltap > Scalar(0)) ? (Scalar(1) / wdeltap) : Scalar(0);
                }
            Scalar wij_val = sph_wij<KT_>(kp.alpha, meanh, r);
            Scalar fab = wij_val * od_wdeltap;
            fab = fab * fab * fab * fab;
            Scalar ti = (Pi / (rhoi * rhoi)) * (Pi > Scalar(0) ?  tensil_eps_pos : -tensil_eps_neg);
            Scalar tj = (Pj / (rhoj * rhoj)) * (Pj > Scalar(0) ?  tensil_eps_pos : -tensil_eps_neg);
            Scalar tc = mi * mj * fab * (ti + tj);
            fi.x -= tc * dwdr_r * dx.x;
            fi.y -= tc * dwdr_r * dx.y;
            fi.z -= tc * dwdr_r * dx.z;
            }

        // Viscous force
        {
        Scalar dvnorm    = sqrtf(dot(dv, dv));
        Scalar gamma_dot = dvnorm / (r + sqrtf(epssqr));
        Scalar mu_eff    = sph_nn_viscosity(nn.mu, gamma_dot, nn.model,
                                            nn.K, nn.n, nn.mu0, nn.muinf,
                                            nn.lambda_NN, nn.tauy, nn.m_reg, nn.mu_min);
        Scalar temp0 = mu_eff * vijsqr * dwdr_r;
        fi.x += temp0 * dv.x;
        fi.y += temp0 * dv.y;
        fi.z += temp0 * dv.z;
        }

        // Artificial stress tensor contribution (TV transport correction)
        {
        Scalar temp0 = Scalar(0.5) * vijsqr * dwdr_r;
        Scalar A1ij = (A11i + A11j) * dx.x + (A12i + A12j) * dx.y + (A13i + A13j) * dx.z;
        Scalar A2ij = (A21i + A21j) * dx.x + (A22i + A22j) * dx.y + (A23i + A23j) * dx.z;
        Scalar A3ij = (A31i + A31j) * dx.x + (A32i + A32j) * dx.y + (A33i + A33j) * dx.z;
        fi.x += temp0 * A1ij;
        fi.y += temp0 * A2ij;
        fi.z += temp0 * A3ij;
        }

        // Background pressure contribution (for transport velocity update in integrator)
        bpc.x -= vijsqr * (Pb / mi) * dwdr_r * dx.x;
        bpc.y -= vijsqr * (Pb / mi) * dwdr_r * dx.y;
        bpc.z -= vijsqr * (Pb / mi) * dwdr_r * dx.z;

        // Density rate (DENSITYCONTINUITY) — use physical advection velocity
        if (density_method == 1)
            {
            // For solid neighbours: use physical velocity (not fictitious)
            Scalar3 vj_adv;
            if (issolid)
                vj_adv = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);
            else
                vj_adv = vj;
            Scalar3 dv_adv = make_scalar3(vi.x - vj_adv.x, vi.y - vj_adv.y, vi.z - vj_adv.z);
            rdpe.x += rhoi * Vj * dot(dv_adv, dwdr_r * dx);

            if (!issolid && density_diffusion)
                rdpe.x -= (Scalar(2) * ddiff * meanh * eos.c * mj
                           * (rhoi / rhoj - Scalar(1))
                           * dot(dx, dwdr_r * dx)) / (rsq + epssqr);
            }
        }

    // dp/dt chain rule (DENSITYCONTINUITY)
    if (density_method == 1)
        rdpe.y = sph_dpressuredrho<SET_>(eos.rho0, eos.c, rhoi) * rdpe.x;

    d_force[i]   = fi;
    d_ratedpe[i] = rdpe;
    d_bpc[i]     = bpc;

    // Max-velocity reduction
    Scalar vi_total = sqrtf(vi.x * vi.x + vi.y * vi.y + vi.z * vi.z);
    atomicMax(d_max_vel_bits, __float_as_uint(vi_total));
    }

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
    unsigned int          block_size)
    {
    if (group_size == 0) return hipSuccess;
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_tv_forcecomputation_kernel<KT_, SET_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array,
                       d_pos, d_vel, d_density, d_pressure, d_vf, d_tv, d_h,
                       d_force, d_ratedpe, d_bpc,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map,
                       d_max_vel_bits, box, kp, eos, nn, Pb,
                       density_method, artificial_viscosity, avalpha, avbeta,
                       tensil_correction, tensil_eps_pos, tensil_eps_neg,
                       density_diffusion, ddiff);
    return hipSuccess;
    }

// =========================================================================
// Explicit template instantiations
// =========================================================================

#define INST_TV_GPU(KT, SET) \
    template hipError_t gpu_sph_tv_forcecomputation<KT, SET>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, const Scalar*, const Scalar*, \
        const Scalar3*, const Scalar3*, const Scalar*, \
        Scalar4*, Scalar4*, Scalar3*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, uint32_t*, \
        BoxDim, SPHKernelDevParams, SPHEOSDevParams, SPHNNViscParams, \
        Scalar, int, int, Scalar, Scalar, int, Scalar, Scalar, int, Scalar, \
        unsigned int);

INST_TV_GPU(wendlandc2, tait)
INST_TV_GPU(wendlandc2, linear)
INST_TV_GPU(wendlandc4, tait)
INST_TV_GPU(wendlandc4, linear)
INST_TV_GPU(wendlandc6, tait)
INST_TV_GPU(wendlandc6, linear)
INST_TV_GPU(quintic,    tait)
INST_TV_GPU(quintic,    linear)
INST_TV_GPU(cubicspline,tait)
INST_TV_GPU(cubicspline,linear)

#undef INST_TV_GPU

} // namespace kernel
} // namespace sph
} // namespace hoomd

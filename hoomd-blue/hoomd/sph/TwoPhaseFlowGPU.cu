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
      2. gpu_sph_2pf_solid_forces     — solid particle reaction forces
*/

#include "TwoPhaseFlowGPU.cuh"
#include "SPHDeviceFunctions.cuh"

namespace hoomd
{
namespace sph
{
namespace kernel
{

// =========================================================================
// Kernel 1: Two-phase force computation (fluid particles)
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
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
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;

    unsigned int i = d_index_array[group_idx];

    // ── Particle i properties ──────────────────────────────────────────────
    bool i_isfluid1 = sph_checkfluid1(d_type_property_map, d_pos[i].w);
    SPHEOSDevParams eosi = i_isfluid1 ? eos1 : eos2;
    SPHNNViscParams nni  = i_isfluid1 ? nn1  : nn2;

    Scalar3 pi   = make_scalar3(d_pos[i].x, d_pos[i].y, d_pos[i].z);
    Scalar3 vi   = make_scalar3(d_vel[i].x, d_vel[i].y, d_vel[i].z);
    Scalar  mi   = d_vel[i].w;
    Scalar  Pi   = d_pressure[i];
    Scalar  rhoi = d_density[i];
    Scalar  Vi   = mi / rhoi;
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];

    Scalar4 fi   = make_scalar4(Scalar(0), Scalar(0), Scalar(0), Scalar(0));
    Scalar4 rdpe = make_scalar4(Scalar(0), Scalar(0), Scalar(0), Scalar(0));

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];

    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];

        bool issolid  = sph_checksolid(d_type_property_map, d_pos[k].w);
        bool j_isfluid1 = sph_checkfluid1(d_type_property_map, d_pos[k].w);

        Scalar3 pj  = make_scalar3(d_pos[k].x, d_pos[k].y, d_pos[k].z);
        Scalar3 dx  = box.minImage(make_scalar3(pi.x - pj.x, pi.y - pj.y, pi.z - pj.z));
        Scalar  rsq = dot(dx, dx);

        if (kp.const_slength && rsq > kp.rcutsq) continue;

        Scalar mj  = d_vel[k].w;
        if (mj < Scalar(0)) continue; // removed solid particle

        Scalar3 vj;
        if (issolid)
            vj = make_scalar3(d_vf[k].x, d_vf[k].y, d_vf[k].z);
        else
            vj = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);

        Scalar rhoj = d_density[k];
        Scalar Vj   = mj / rhoj;
        Scalar Pj   = d_pressure[k];

        SPHEOSDevParams eosj = issolid ? eosi : (j_isfluid1 ? eos1 : eos2);
        SPHNNViscParams nnj  = issolid ? nni  : (j_isfluid1 ? nn1  : nn2);

        Scalar3 dv = make_scalar3(vi.x - vj.x, vi.y - vj.y, vi.z - vj.z);
        Scalar r   = sqrtf(rsq);

        Scalar meanh  = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar epssqr = Scalar(0.01) * meanh * meanh;

        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = dwdr / (r + sqrtf(epssqr));

        // ── Pressure force ─────────────────────────────────────────────────
        Scalar prefactor;
        Scalar p_bar;
        if (fparams.density_method == 0) // DENSITYSUMMATION
            {
            prefactor = Vi * Vi + Vj * Vj;
            // Consistent interface pressure for cross-phase fluid pairs
            if (fparams.cip && !issolid && (i_isfluid1 != j_isfluid1))
                {
                Scalar3 gvec = make_scalar3(fparams.gvec_x, fparams.gvec_y, fparams.gvec_z);
                p_bar = (eosj.rho0 * Pi + eosi.rho0 * Pj
                         + eosi.rho0 * eosj.rho0 * dot(gvec, dx))
                        / (eosi.rho0 + eosj.rho0);
                }
            else
                p_bar = (rhoj * Pi + rhoi * Pj) / (rhoi + rhoj);
            }
        else // DENSITYCONTINUITY
            {
            prefactor = mi * mj;
            p_bar = (Pi + Pj) / (rhoi * rhoj);
            }

        // ── Momentum dissipation (fluid–fluid pairs only) ──────────────────
        Scalar avc = Scalar(0);
        if (!issolid)
            {
            Scalar dotdvdx = dot(dv, dx);
            // [A] Monaghan artificial viscosity (Monaghan 1992)
            if (fparams.artificial_viscosity && dotdvdx < Scalar(0))
                {
                Scalar muij    = meanh * dotdvdx / (rsq + epssqr);
                Scalar meanrho = Scalar(0.5) * (rhoi + rhoj);
                avc = (-fparams.avalpha * fparams.cmax * muij
                       + fparams.avbeta * muij * muij) / meanrho;
                }
            // [B] Riemann dissipation (Zhang, Hu & Adams 2017)
            else if (fparams.riemann_dissipation && dotdvdx < Scalar(0))
                {
                Scalar r_eps = r + sqrtf(epssqr);
                Scalar uij   = dotdvdx / r_eps;
                Scalar Zi    = rhoi * eosi.c;
                Scalar Zj    = rhoj * eosj.c;
                Scalar Zstar = (Zi * Zj) / (Zi + Zj);
                Scalar meanrho = Scalar(0.5) * (rhoi + rhoj);
                avc = -fparams.riemann_beta * Zstar * uij / meanrho;
                }
            }

        fi.x -= prefactor * (p_bar + avc) * dwdr_r * dx.x;
        fi.y -= prefactor * (p_bar + avc) * dwdr_r * dx.y;
        fi.z -= prefactor * (p_bar + avc) * dwdr_r * dx.z;

        // ── Viscous force (harmonic mean viscosity) ────────────────────────
        {
        Scalar dvnorm    = sqrtf(dot(dv, dv));
        Scalar gamma_dot = dvnorm / (r + sqrtf(epssqr));
        Scalar mu_eff_i  = sph_nn_viscosity(nni.mu, gamma_dot, nni.model,
                                             nni.K, nni.n, nni.mu0, nni.muinf,
                                             nni.lambda_NN, nni.tauy, nni.m_reg, nni.mu_min);
        Scalar mu_eff_j  = sph_nn_viscosity(nnj.mu, gamma_dot, nnj.model,
                                             nnj.K, nnj.n, nnj.mu0, nnj.muinf,
                                             nnj.lambda_NN, nnj.tauy, nnj.m_reg, nnj.mu_min);
        Scalar denom = mu_eff_i + mu_eff_j;
        Scalar mu_harm = (denom > Scalar(0))
                         ? Scalar(2) * mu_eff_i * mu_eff_j / denom
                         : Scalar(0);
        Scalar tv = mu_harm * prefactor * dwdr_r;
        fi.x += tv * dv.x;
        fi.y += tv * dv.y;
        fi.z += tv * dv.z;
        }

        // ── Density continuity rate (DENSITYCONTINUITY only) ──────────────
        if (fparams.density_method == 1)
            {
            Scalar3 vj_adv;
            if (issolid)
                vj_adv = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);
            else
                vj_adv = vj;
            Scalar3 dv_adv = make_scalar3(vi.x - vj_adv.x, vi.y - vj_adv.y, vi.z - vj_adv.z);
            rdpe.x += rhoi * Vj * dot(dv_adv, dwdr_r * dx);

            // Molteni density diffusion (two-phase corrected: rest-density normalised)
            if (!issolid && fparams.density_diffusion)
                {
                Scalar drive = rhoi / eosi.rho0 - rhoj / eosj.rho0;
                rdpe.x -= (Scalar(2) * fparams.ddiff * meanh * fparams.cmax * mj
                           * drive * dot(dx, dwdr_r * dx)) / (rsq + epssqr);
                }
            }
        }

    // dp/dt via chain rule (DENSITYCONTINUITY only)
    if (fparams.density_method == 1)
        {
        Scalar dpdrho_i;
        if (i_isfluid1)
            dpdrho_i = sph_dpressuredrho<SET1_>(eos1.rho0, eos1.c, rhoi);
        else
            dpdrho_i = sph_dpressuredrho<SET2_>(eos2.rho0, eos2.c, rhoi);
        rdpe.y = dpdrho_i * rdpe.x;
        }

    // Add pre-computed surface force from aux4
    fi.x += d_sf[i].x;
    fi.y += d_sf[i].y;
    fi.z += d_sf[i].z;

    d_force[i]   = fi;
    d_ratedpe[i] = rdpe;

    // Max-velocity reduction
    Scalar vi_total = sqrtf(vi.x * vi.x + vi.y * vi.y + vi.z * vi.z);
    atomicMax(d_max_vel_bits, __float_as_uint(vi_total));
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
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_forcecomputation_kernel<KT_, SET1_, SET2_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array,
                       d_pos, d_vel, d_density, d_pressure, d_vf, d_sf, d_h,
                       d_force, d_ratedpe,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map,
                       d_max_vel_bits, box, kp, eos1, eos2, nn1, nn2, fparams);
    return hipSuccess;
    }

// =========================================================================
// Kernel 2: Two-phase solid particle reaction forces
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
__global__ void gpu_sph_2pf_solid_forces_kernel(
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
    int                   density_method)
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

    Scalar4 fi = make_scalar4(Scalar(0), Scalar(0), Scalar(0), Scalar(0));

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];

    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];

        if (sph_checksolid(d_type_property_map, d_pos[k].w)) continue;

        bool j_isfluid1 = sph_checkfluid1(d_type_property_map, d_pos[k].w);
        SPHNNViscParams nnj = j_isfluid1 ? nn1 : nn2;

        Scalar3 pj  = make_scalar3(d_pos[k].x, d_pos[k].y, d_pos[k].z);
        Scalar3 dx  = box.minImage(make_scalar3(pi.x - pj.x, pi.y - pj.y, pi.z - pj.z));
        Scalar  rsq = dot(dx, dx);

        if (kp.const_slength && rsq > kp.rcutsq) continue;

        Scalar3 vj  = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);
        Scalar  mj  = d_vel[k].w;
        Scalar  rhoj = d_density[k];
        Scalar  Vj   = mj / rhoj;
        Scalar  Pj   = d_pressure[k];

        Scalar3 dv = make_scalar3(vi.x - vj.x, vi.y - vj.y, vi.z - vj.z);
        Scalar r   = sqrtf(rsq);

        Scalar meanh  = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar epssqr = Scalar(0.01) * meanh * meanh;

        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = dwdr / (r + sqrtf(epssqr));

        Scalar pressure_term;
        if (density_method == 0)
            pressure_term = -(Vi * Vi + Vj * Vj) * (rhoj * Pi + rhoi * Pj) / (rhoi + rhoj);
        else
            pressure_term = -mi * mj * (Pi + Pj) / (rhoi * rhoj);

        Scalar scale = mj / mi;
        fi.x -= scale * pressure_term * dwdr_r * dx.x;
        fi.y -= scale * pressure_term * dwdr_r * dx.y;
        fi.z -= scale * pressure_term * dwdr_r * dx.z;

        {
        Scalar dvnorm    = sqrtf(dot(dv, dv));
        Scalar gamma_dot = dvnorm / (r + sqrtf(epssqr));
        Scalar mu_eff_j  = sph_nn_viscosity(nnj.mu, gamma_dot, nnj.model,
                                             nnj.K, nnj.n, nnj.mu0, nnj.muinf,
                                             nnj.lambda_NN, nnj.tauy, nnj.m_reg, nnj.mu_min);
        Scalar viscous_term;
        if (density_method == 0)
            viscous_term = mu_eff_j * (Vi * Vi + Vj * Vj) * dwdr_r;
        else
            viscous_term = mu_eff_j * mi * mj * dwdr_r;
        fi.x -= scale * viscous_term * dv.x;
        fi.y -= scale * viscous_term * dv.y;
        fi.z -= scale * viscous_term * dv.z;
        }
        }

    // Accumulate into global force (readwrite: forcecomputation already wrote fluid forces)
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
    unsigned int          block_size)
    {
    if (group_size == 0) return hipSuccess;
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_solid_forces_kernel<KT_, SET1_, SET2_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array,
                       d_pos, d_vel, d_density, d_pressure, d_h,
                       d_force,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map,
                       box, kp, nn1, nn2, density_method);
    return hipSuccess;
    }

// =========================================================================
// Explicit template instantiations
// =========================================================================

#define INST_2PF_GPU(KT, SET1, SET2) \
    template hipError_t gpu_sph_2pf_forcecomputation<KT, SET1, SET2>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, const Scalar*, const Scalar*, \
        const Scalar3*, const Scalar3*, const Scalar*, \
        Scalar4*, Scalar4*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, uint32_t*, BoxDim, \
        SPHKernelDevParams, SPHEOSDevParams, SPHEOSDevParams, \
        SPHNNViscParams, SPHNNViscParams, SPHTwoPhaseParams, unsigned int); \
    template hipError_t gpu_sph_2pf_solid_forces<KT, SET1, SET2>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, const Scalar*, const Scalar*, const Scalar*, \
        Scalar4*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, BoxDim, SPHKernelDevParams, \
        SPHNNViscParams, SPHNNViscParams, int, unsigned int);

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

} // namespace kernel
} // namespace sph
} // namespace hoomd

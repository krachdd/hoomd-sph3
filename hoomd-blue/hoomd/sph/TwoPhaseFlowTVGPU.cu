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

    Extends the two-phase force kernel with:
      - Artificial-stress tensor (Adami 2013) for tensile instability suppression
      - Background-pressure contribution (BPC) written to aux2
    Both TV terms are applied to fluid-fluid pairs only (not fluid-solid).
*/

#include "TwoPhaseFlowTVGPU.cuh"
#include "TwoPhaseFlowGPU.cuh"
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

    Scalar3 pi   = make_scalar3(d_pos[i].x, d_pos[i].y, d_pos[i].z);
    Scalar3 vi   = make_scalar3(d_vel[i].x, d_vel[i].y, d_vel[i].z);
    Scalar  mi   = d_vel[i].w;
    Scalar  Pi   = d_pressure[i];
    Scalar  rhoi = d_density[i];
    Scalar  Vi   = mi / rhoi;
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];

    bool i_isfluid1 = sph_checkfluid1(d_type_property_map, d_pos[i].w);
    SPHNNViscParams nni  = i_isfluid1 ? nn1 : nn2;
    SPHEOSDevParams eosi = i_isfluid1 ? eos1 : eos2;

    // Transport velocity and per-phase background pressure
    Scalar3 tvi = d_tv[i];
    Scalar  Pbi = i_isfluid1 ? Pb1 : Pb2;

    // Artificial stress tensor for particle i: A = rho_i * v_i ⊗ (tv_i - v_i)
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

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];

    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];

        bool j_issolid  = sph_checksolid(d_type_property_map, d_pos[k].w);
        bool j_isfluid1 = sph_checkfluid1(d_type_property_map, d_pos[k].w);

        SPHNNViscParams nnj  = j_isfluid1 ? nn1 : nn2;
        SPHEOSDevParams eosj = j_isfluid1 ? eos1 : eos2;
        if (j_issolid) { nnj = nni; eosj = eosi; }

        Scalar3 pj  = make_scalar3(d_pos[k].x, d_pos[k].y, d_pos[k].z);
        Scalar3 dx  = box.minImage(make_scalar3(pi.x - pj.x, pi.y - pj.y, pi.z - pj.z));
        Scalar  rsq = dot(dx, dx);

        if (kp.const_slength && rsq > kp.rcutsq) continue;

        // Neighbour velocity: fictitious for solid, physical for fluid
        Scalar3 vj;
        if (j_issolid)
            vj = make_scalar3(d_vf[k].x, d_vf[k].y, d_vf[k].z);
        else
            vj = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);

        Scalar  mj   = d_vel[k].w;
        Scalar  rhoj = d_density[k];
        Scalar  Vj   = mj / rhoj;
        Scalar  Pj   = d_pressure[k];

        Scalar3 dv = make_scalar3(vi.x - vj.x, vi.y - vj.y, vi.z - vj.z);
        Scalar  r  = sqrtf(rsq);

        Scalar meanh  = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar epssqr = Scalar(0.01) * meanh * meanh;
        Scalar eps    = sqrtf(epssqr);

        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = dwdr / (r + eps);

        // ── Pressure force ────────────────────────────────────────────────
        Scalar prefactor, pterm;
        if (fparams.density_method == 0) // DENSITYSUMMATION
            {
            if (fparams.cip && !j_issolid && (i_isfluid1 != j_isfluid1))
                pterm = (eosj.rho0*Pi + eosi.rho0*Pj
                         + eosi.rho0*eosj.rho0*(fparams.gvec_x*dx.x
                                                +fparams.gvec_y*dx.y
                                                +fparams.gvec_z*dx.z))
                        / (eosi.rho0 + eosj.rho0);
            else
                pterm = (rhoj*Pi + rhoi*Pj) / (rhoi + rhoj);
            prefactor = Vi*Vi + Vj*Vj;
            }
        else // DENSITYCONTINUITY
            {
            pterm    = (Pi + Pj) / (rhoi * rhoj);
            prefactor = mi * mj;
            }

        // ── Momentum dissipation (fluid-fluid only) ───────────────────────
        Scalar avc = Scalar(0);
        if (!j_issolid)
            {
            Scalar dotdvdx = dot(dv, dx);
            if (fparams.artificial_viscosity && dotdvdx < Scalar(0))
                {
                Scalar muij    = meanh * dotdvdx / (rsq + epssqr);
                Scalar meanrho = Scalar(0.5) * (rhoi + rhoj);
                avc = (-fparams.avalpha * fparams.cmax * muij + fparams.avbeta * muij * muij) / meanrho;
                }
            else if (fparams.riemann_dissipation && dotdvdx < Scalar(0))
                {
                Scalar uij     = dotdvdx / (r + eps);
                Scalar Zi      = rhoi * eosi.c;
                Scalar Zj      = rhoj * eosj.c;
                Scalar Zstar   = (Zi * Zj) / (Zi + Zj);
                Scalar meanrho = Scalar(0.5) * (rhoi + rhoj);
                avc = -fparams.riemann_beta * Zstar * uij / meanrho;
                }
            }

        fi.x -= prefactor * (pterm + avc) * dwdr_r * dx.x;
        fi.y -= prefactor * (pterm + avc) * dwdr_r * dx.y;
        fi.z -= prefactor * (pterm + avc) * dwdr_r * dx.z;

        // ── Viscosity ─────────────────────────────────────────────────────
        {
        Scalar dvnorm    = sqrtf(dot(dv, dv));
        Scalar gamma_dot = dvnorm / (r + eps);
        Scalar mu_eff_i = sph_nn_viscosity(nni.mu, gamma_dot, nni.model,
                                            nni.K, nni.n, nni.mu0, nni.muinf,
                                            nni.lambda_NN, nni.tauy, nni.m_reg, nni.mu_min);
        Scalar mu_eff_j;
        if (j_issolid)
            mu_eff_j = mu_eff_i;
        else
            mu_eff_j = sph_nn_viscosity(nnj.mu, gamma_dot, nnj.model,
                                         nnj.K, nnj.n, nnj.mu0, nnj.muinf,
                                         nnj.lambda_NN, nnj.tauy, nnj.m_reg, nnj.mu_min);
        Scalar mu_harm = Scalar(2) * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
        Scalar visc_fac = mu_harm * (Vi*Vi + Vj*Vj) * dwdr_r;
        fi.x += visc_fac * dv.x;
        fi.y += visc_fac * dv.y;
        fi.z += visc_fac * dv.z;
        }

        // ── Transport-velocity terms (fluid-fluid pairs only) ─────────────
        if (!j_issolid)
            {
            Scalar3 tvj = d_tv[k];
            Scalar3 vkphys = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);

            Scalar A11j = rhoj * vkphys.x * (tvj.x - vkphys.x);
            Scalar A12j = rhoj * vkphys.x * (tvj.y - vkphys.y);
            Scalar A13j = rhoj * vkphys.x * (tvj.z - vkphys.z);
            Scalar A21j = rhoj * vkphys.y * (tvj.x - vkphys.x);
            Scalar A22j = rhoj * vkphys.y * (tvj.y - vkphys.y);
            Scalar A23j = rhoj * vkphys.y * (tvj.z - vkphys.z);
            Scalar A31j = rhoj * vkphys.z * (tvj.x - vkphys.x);
            Scalar A32j = rhoj * vkphys.z * (tvj.y - vkphys.y);
            Scalar A33j = rhoj * vkphys.z * (tvj.z - vkphys.z);

            Scalar vijsqr  = Vi*Vi + Vj*Vj;
            Scalar tv_temp = Scalar(0.5) * vijsqr * dwdr_r;
            Scalar A1ij = (A11i+A11j)*dx.x + (A12i+A12j)*dx.y + (A13i+A13j)*dx.z;
            Scalar A2ij = (A21i+A21j)*dx.x + (A22i+A22j)*dx.y + (A23i+A23j)*dx.z;
            Scalar A3ij = (A31i+A31j)*dx.x + (A32i+A32j)*dx.y + (A33i+A33j)*dx.z;
            fi.x += tv_temp * A1ij;
            fi.y += tv_temp * A2ij;
            fi.z += tv_temp * A3ij;

            // BPC contribution
            bpc.x -= vijsqr * Pbi / mi * dwdr_r * dx.x;
            bpc.y -= vijsqr * Pbi / mi * dwdr_r * dx.y;
            bpc.z -= vijsqr * Pbi / mi * dwdr_r * dx.z;
            }

        // ── Density continuity rate ───────────────────────────────────────
        if (fparams.density_method == 1) // DENSITYCONTINUITY
            {
            Scalar3 dv_cont = dv;
            if (j_issolid)
                {
                // Use physical velocity for continuity with solids
                Scalar3 vj_phys = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);
                dv_cont = make_scalar3(vi.x - vj_phys.x, vi.y - vj_phys.y, vi.z - vj_phys.z);
                }
            rdpe.x += rhoi * Vj * dot(dv_cont, dwdr_r * dx);

            if (!j_issolid && fparams.density_diffusion)
                {
                Scalar drive = rhoi / eosi.rho0 - rhoj / eosj.rho0;
                rdpe.x -= (Scalar(2)*fparams.ddiff*meanh*fparams.cmax*mj*drive
                           * dot(dx, dwdr_r*dx)) / (rsq + epssqr);
                }
            }
        } // end neighbour loop

    // dp/dt via chain rule
    if (fparams.density_method == 1)
        {
        Scalar dpdrho_i = i_isfluid1 ? sph_dpressuredrho<SET1_>(eosi.rho0, eosi.c, rhoi)
                                     : sph_dpressuredrho<SET2_>(eosi.rho0, eosi.c, rhoi);
        rdpe.y = dpdrho_i * rdpe.x;
        }

    // Surface force (pre-computed on CPU by compute_surfaceforce)
    fi.x += d_sf[i].x;
    fi.y += d_sf[i].y;
    fi.z += d_sf[i].z;

    d_force[i]   = fi;
    d_ratedpe[i] = rdpe;
    d_bpc[i]     = bpc;

    // Max-velocity reduction
    Scalar vi_total = sqrtf(vi.x*vi.x + vi.y*vi.y + vi.z*vi.z);
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
                       d_vf, d_tv, d_sf, d_h,
                       d_force, d_ratedpe, d_bpc,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map,
                       d_max_vel_bits, box, kp, eos1, eos2, nn1, nn2, fparams, Pb1, Pb2);
    return hipSuccess;
    }

// =========================================================================
// Explicit template instantiations (5 × 2 × 2 = 20)
// =========================================================================

#define INST_2PF_TV_GPU(KT, SET1, SET2) \
    template hipError_t gpu_sph_2pf_tv_forcecomputation<KT, SET1, SET2>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, const Scalar*, const Scalar*, \
        const Scalar3*, const Scalar3*, const Scalar3*, const Scalar*, \
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

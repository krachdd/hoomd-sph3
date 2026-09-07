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

/*! \file TwoPhaseFlowGPUPair.cuh
    \brief Shared device code for the two-phase pair interaction.

    Single source of truth for the fluid-particle pair force of the GPU
    two-phase solvers (TwoPhaseFlowGPU and TwoPhaseFlowTVGPU).  It mirrors
    TwoPhaseFlow<KT_,SET1_,SET2_>::forcecomputation() in TwoPhaseFlow.cc
    term by term:

      - pressure force, symmetric volume form, with optional consistent
        interface pressure (Hu & Adams 2009) for cross-phase pairs
      - Monaghan (1992) artificial viscosity  diss = m_i m_j Pi_ij
        OR Riemann dissipation (Zhang 2017)   diss = (V_i^2+V_j^2) p_d
      - harmonic-mean viscous force, always (V_i^2+V_j^2) weighted,
        non-Newtonian viscosity evaluated at the PER-PARTICLE shear rate
        stored in the energy array by compute_strain_rate()
      - continuity-equation density rate with Molteni-Colagrossi diffusion
        (rest-density normalised drive term)

    Any change to the CPU pair terms must be replicated here — see
    the CPU-vs-GPU consistency check used when this file was written.
*/

#ifndef __TWO_PHASE_FLOW_GPU_PAIR_CUH__
#define __TWO_PHASE_FLOW_GPU_PAIR_CUH__

#include "hip/hip_runtime.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"

#include "SPHDeviceFunctions.cuh"
#include "TwoPhaseFlowGPU.cuh"

namespace hoomd
{
namespace sph
{
namespace kernel
{

//! State of the central fluid particle i, loaded once per particle
struct SPH2PFParticle
    {
    Scalar3 pos;      //!< position
    Scalar3 vel;      //!< physical velocity
    Scalar  m;        //!< mass
    Scalar  P;        //!< pressure
    Scalar  rho;      //!< density
    Scalar  V;        //!< volume m/rho
    Scalar  h;        //!< smoothing length (kp.ch when constant)
    Scalar  rho0;     //!< rest density of the particle's phase
    Scalar  c;        //!< speed of sound of the particle's phase
    Scalar  mu_eff;   //!< effective viscosity at the particle's own shear rate
    bool    isfluid1; //!< true if the particle belongs to fluid 1
    };

//! Pair quantities handed back to the caller (needed by the TV extension)
struct SPH2PFPairOut
    {
    bool    valid;     //!< false if the pair was skipped (outside cutoff)
    bool    j_issolid; //!< neighbour is a solid (wall) particle
    Scalar3 dx;        //!< r_i - r_j (minimum image)
    Scalar  dwdr_r;    //!< (dW/dr)/r
    Scalar  vijsqr;    //!< V_i^2 + V_j^2
    Scalar  rhoj;      //!< neighbour density
    Scalar  Vj;        //!< neighbour volume
    Scalar  mj;        //!< neighbour mass
    };

/*! Effective viscosity of a phase at a given shear rate.
 *  FAST=true resolves the model at compile time (NM1_ for fluid 1, NM2_ for
 *  fluid 2); FAST=false uses the runtime switch in sph_nn_viscosity.
 */
template<bool FAST, NonNewtonianModel NM1_, NonNewtonianModel NM2_>
__device__ __forceinline__ Scalar sph_2pf_mu_eff(bool isfluid1, Scalar gdot,
                                                 const SPHNNViscParams& nn1,
                                                 const SPHNNViscParams& nn2)
    {
    if constexpr (FAST)
        {
        return isfluid1 ? sph_nn_viscosity_t<NM1_>(nn1, gdot)
                        : sph_nn_viscosity_t<NM2_>(nn2, gdot);
        }
    else
        {
        const SPHNNViscParams& nn = isfluid1 ? nn1 : nn2;
        return sph_nn_viscosity(nn.mu, gdot, nn.model, nn.K, nn.n, nn.mu0, nn.muinf,
                                nn.lambda_NN, nn.tauy, nn.m_reg, nn.mu_min);
        }
    }

/*! Load the per-particle state of fluid particle i. */
template<bool FAST, NonNewtonianModel NM1_, NonNewtonianModel NM2_>
__device__ __forceinline__ SPH2PFParticle sph_2pf_load_particle(
    unsigned int          i,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar*         d_h,
    const Scalar*         d_gdot,
    const unsigned int*   d_type_property_map,
    const SPHKernelDevParams& kp,
    const SPHEOSDevParams& eos1,
    const SPHEOSDevParams& eos2,
    const SPHNNViscParams& nn1,
    const SPHNNViscParams& nn2,
    const SPHTwoPhaseParams& fp)
    {
    SPH2PFParticle I;
    Scalar4 posi = d_pos[i];
    Scalar4 veli = d_vel[i];
    I.isfluid1 = sph_checkfluid1(d_type_property_map, posi.w);
    I.pos  = make_scalar3(posi.x, posi.y, posi.z);
    I.vel  = make_scalar3(veli.x, veli.y, veli.z);
    I.m    = veli.w;
    I.P    = d_pressure[i];
    I.rho  = d_density[i];
    I.V    = I.m / I.rho;
    I.h    = kp.const_slength ? kp.ch : d_h[i];
    I.rho0 = I.isfluid1 ? eos1.rho0 : eos2.rho0;
    I.c    = I.isfluid1 ? eos1.c    : eos2.c;
    // Per-particle shear rate from compute_strain_rate() (energy array);
    // zero when no non-Newtonian model is active (then ignored anyway).
    Scalar gdot_i = fp.nn_active ? d_gdot[i] : Scalar(0);
    I.mu_eff = sph_2pf_mu_eff<FAST, NM1_, NM2_>(I.isfluid1, gdot_i, nn1, nn2);
    return I;
    }

/*! Pair interaction of fluid particle i with neighbour k.
 *
 *  Accumulates the pressure + dissipation + viscous force into \a fi and the
 *  continuity density rate into \a drho.  Returns pair geometry so the TV
 *  kernel can add its fluid-fluid-only terms on top.
 */
template<SmoothingKernelType KT_, bool FAST, NonNewtonianModel NM1_, NonNewtonianModel NM2_>
__device__ __forceinline__ SPH2PFPairOut sph_2pf_pair(
    const SPH2PFParticle& I,
    unsigned int          k,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar*         d_h,
    const Scalar*         d_gdot,
    const unsigned int*   d_type_property_map,
    const BoxDim&         box,
    const SPHKernelDevParams& kp,
    const SPHEOSDevParams& eos1,
    const SPHEOSDevParams& eos2,
    const SPHNNViscParams& nn1,
    const SPHNNViscParams& nn2,
    const SPHTwoPhaseParams& fp,
    Scalar3&              fi,
    Scalar&               drho)
    {
    SPH2PFPairOut o;
    o.valid = false;

    Scalar4 posk = d_pos[k];
    Scalar4 velk = d_vel[k];

    o.j_issolid     = sph_checksolid (d_type_property_map, posk.w);
    bool j_isfluid1 = sph_checkfluid1(d_type_property_map, posk.w);

    // Neighbour phase properties; a solid neighbour takes the properties of i
    Scalar rho0j = o.j_issolid ? I.rho0 : (j_isfluid1 ? eos1.rho0 : eos2.rho0);
    Scalar cj    = o.j_issolid ? I.c    : (j_isfluid1 ? eos1.c    : eos2.c);

    Scalar3 dx  = box.minImage(make_scalar3(I.pos.x - posk.x, I.pos.y - posk.y, I.pos.z - posk.z));
    Scalar  rsq = dot(dx, dx);
    if (kp.const_slength && rsq > kp.rcutsq) return o;

    // Neighbour velocity: fictitious (Adami) velocity for solids
    Scalar3 vj;
    if (o.j_issolid)
        { Scalar3 vfk = d_vf[k]; vj = make_scalar3(vfk.x, vfk.y, vfk.z); }
    else
        vj = make_scalar3(velk.x, velk.y, velk.z);

    Scalar mj   = velk.w;
    Scalar rhoj = d_density[k];
    Scalar Vj   = mj / rhoj;
    Scalar Pj   = d_pressure[k];

    Scalar3 dv = make_scalar3(I.vel.x - vj.x, I.vel.y - vj.y, I.vel.z - vj.z);
    Scalar  r  = sqrt(rsq);

    Scalar meanh  = kp.const_slength ? kp.ch : Scalar(0.5) * (I.h + d_h[k]);
    Scalar eps    = Scalar(0.1) * meanh;
    Scalar epssqr = eps * eps;

    Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
    Scalar dwdr_r = (r > Scalar(1e-8) * meanh) ? dwdr / r : Scalar(0);

    // ── Pressure force (symmetric volume form) ─────────────────────────
    Scalar prefactor, pbar;
    if (fp.density_method == 0) // DENSITYSUMMATION
        {
        if (fp.cip && !o.j_issolid && (I.isfluid1 != j_isfluid1))
            {
            Scalar gdotdx = fp.gvec_x * dx.x + fp.gvec_y * dx.y + fp.gvec_z * dx.z;
            pbar = (rho0j * I.P + I.rho0 * Pj + I.rho0 * rho0j * gdotdx) / (I.rho0 + rho0j);
            }
        else
            pbar = (rhoj * I.P + I.rho * Pj) / (I.rho + rhoj);
        prefactor = I.V * I.V + Vj * Vj;
        }
    else // DENSITYCONTINUITY
        {
        pbar      = (I.P + Pj) / (I.rho * rhoj);
        prefactor = I.m * mj;
        }

    // ── Momentum dissipation (fluid-fluid pairs only) ──────────────────
    // [A] Monaghan AV: diss = m_i m_j Pi_ij     [B] Riemann: diss = (V_i^2+V_j^2) p_d
    Scalar diss = Scalar(0);
    if (fp.artificial_viscosity && !o.j_issolid)
        {
        Scalar dotdvdx = dot(dv, dx);
        if (dotdvdx < Scalar(0))
            {
            Scalar muij    = meanh * dotdvdx / (rsq + epssqr);
            Scalar meanrho = Scalar(0.5) * (I.rho + rhoj);
            diss = I.m * mj * (-fp.avalpha * fp.cmax * muij + fp.avbeta * muij * muij) / meanrho;
            }
        }
    else if (fp.riemann_dissipation && !o.j_issolid)
        {
        Scalar dotdvdx = dot(dv, dx);
        if (dotdvdx < Scalar(0))
            {
            Scalar uij   = dotdvdx / (r + eps);
            Scalar Zi    = I.rho * I.c;
            Scalar Zj    = rhoj * cj;
            Scalar Zstar = (Zi * Zj) / (Zi + Zj);
            Scalar pd    = -fp.riemann_beta * Zstar * uij;
            diss = (I.V * I.V + Vj * Vj) * pd;
            }
        }

    Scalar pcoef = (prefactor * pbar + diss) * dwdr_r;
    fi.x -= pcoef * dx.x;
    fi.y -= pcoef * dx.y;
    fi.z -= pcoef * dx.z;

    // ── Viscous force (harmonic-mean viscosity, per-particle shear rate) ──
    {
    Scalar mu_eff_j;
    if (o.j_issolid)
        mu_eff_j = I.mu_eff;
    else
        {
        Scalar gdot_j = fp.nn_active ? d_gdot[k] : Scalar(0);
        mu_eff_j = sph_2pf_mu_eff<FAST, NM1_, NM2_>(j_isfluid1, gdot_j, nn1, nn2);
        }
    Scalar denom   = I.mu_eff + mu_eff_j;
    Scalar mu_harm = (denom > Scalar(0)) ? Scalar(2) * I.mu_eff * mu_eff_j / denom : Scalar(0);
    Scalar vcoef   = mu_harm * (I.V * I.V + Vj * Vj) * dwdr_r;
    fi.x += vcoef * dv.x;
    fi.y += vcoef * dv.y;
    fi.z += vcoef * dv.z;
    }

    // ── Continuity-equation density rate (DENSITYCONTINUITY only) ──────
    if (fp.density_method == 1)
        {
        Scalar3 dv_adv = dv;
        if (o.j_issolid)
            {
            // physical advection velocity of the solid, not the fictitious one
            dv_adv = make_scalar3(I.vel.x - velk.x, I.vel.y - velk.y, I.vel.z - velk.z);
            }
        drho += I.rho * Vj * dot(dv_adv, dwdr_r * dx);

        // Molteni-Colagrossi diffusion, rest-density normalised drive term
        if (!o.j_issolid && fp.density_diffusion)
            drho += Scalar(2) * fp.ddiff * meanh * fp.cmax * (mj / rhoj) * I.rho0
                    * (I.rho / I.rho0 - rhoj / rho0j) * dwdr_r;
        }

    o.valid  = true;
    o.dx     = dx;
    o.dwdr_r = dwdr_r;
    o.vijsqr = I.V * I.V + Vj * Vj;
    o.rhoj   = rhoj;
    o.Vj     = Vj;
    o.mj     = mj;
    return o;
    }

} // namespace kernel
} // namespace sph
} // namespace hoomd

#endif // __TWO_PHASE_FLOW_GPU_PAIR_CUH__

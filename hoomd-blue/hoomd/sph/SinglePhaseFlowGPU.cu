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

/*! \file SinglePhaseFlowGPU.cu
    \brief GPU kernels for single-phase SPH force computation.

    Five kernels mirror the five CPU compute stages in SinglePhaseFlow:
      1. gpu_sph_ndensity_kernel      — number density (DENSITYSUMMATION)
      2. gpu_sph_pressure_kernel      — EOS pressure
      3. gpu_sph_noslip_kernel        — fictitious solid properties (Adami 2012)
      4. gpu_sph_forcecomputation_kernel — main force loop
      5. gpu_sph_solid_forces_kernel  — reaction forces on solid particles

    Each kernel is launched with one thread per particle, serial inner
    neighbour loop.  All five are wrapped as template<KT_, SET_> hipError_t
    functions with explicit instantiations for all 10 type combinations.
*/

#include "SinglePhaseFlowGPU.cuh"
#include "SPHDeviceFunctions.cuh"

namespace hoomd
{
namespace sph
{
namespace kernel
{

// =========================================================================
// Kernel 1: Number density
// =========================================================================

template<SmoothingKernelType KT_>
__global__ void gpu_sph_ndensity_kernel(
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
    SPHKernelDevParams    kp)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;

    unsigned int i = d_index_array[group_idx];

    Scalar3 pi = make_scalar3(d_pos[i].x, d_pos[i].y, d_pos[i].z);
    Scalar  mi = d_vel[i].w;
    Scalar  hi = kp.const_slength ? kp.ch : d_h[i];

    Scalar ni = sph_w0<KT_>(kp.alpha, kp.self_density, hi);

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];

    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];

        Scalar3 pj = make_scalar3(d_pos[k].x, d_pos[k].y, d_pos[k].z);
        Scalar3 dx = box.minImage(make_scalar3(pi.x - pj.x, pi.y - pj.y, pi.z - pj.z));
        Scalar  rsq = dot(dx, dx);

        if (kp.const_slength && rsq > kp.rcutsq) continue;

        Scalar r    = sqrtf(rsq);
        Scalar meanh = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        ni += sph_wij<KT_>(kp.alpha, meanh, r);
        }

    d_density[i] = ni * mi;
    }

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
    unsigned int          block_size)
    {
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_ndensity_kernel<KT_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array, d_pos, d_vel, d_density,
                       d_h, d_n_neigh, d_nlist, d_head_list, box, kp);
    return hipSuccess;
    }

// =========================================================================
// Kernel 2: EOS pressure
// =========================================================================

template<StateEquationType SET_>
__global__ void gpu_sph_pressure_kernel(
    unsigned int        group_size,
    const unsigned int* d_index_array,
    const Scalar*       d_density,
    Scalar*             d_pressure,
    SPHEOSDevParams     eos)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;

    unsigned int i = d_index_array[group_idx];
    d_pressure[i] = sph_pressure<SET_>(eos.rho0, eos.c, eos.bp, d_density[i]);
    }

template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_pressure(
    unsigned int        group_size,
    const unsigned int* d_index_array,
    const Scalar*       d_density,
    Scalar*             d_pressure,
    SPHEOSDevParams     eos,
    unsigned int        block_size)
    {
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_pressure_kernel<SET_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array, d_density, d_pressure, eos);
    return hipSuccess;
    }

// =========================================================================
// Kernel 3: No-slip fictitious solid properties (Adami 2012)
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
__global__ void gpu_sph_noslip_kernel(
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
    Scalar3               bodyforce)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= solid_group_size) return;

    unsigned int i = d_solid_index[group_idx];

    Scalar3 pi  = make_scalar3(d_pos[i].x, d_pos[i].y, d_pos[i].z);
    Scalar3 vi  = make_scalar3(d_vel[i].x, d_vel[i].y, d_vel[i].z);
    Scalar  hi  = kp.const_slength ? kp.ch : d_h[i];

    // Read acceleration (NaN-safe)
    Scalar3 accel_i = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    Scalar3 ai = d_accel[i];
    if (ai.x == ai.x && ai.y == ai.y && ai.z == ai.z)
        accel_i = ai;

    Scalar3 uf_c0     = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    Scalar  pf_c0     = Scalar(0);
    Scalar3 ph_c0     = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    Scalar  wij_c0    = Scalar(0);
    unsigned int fluid_neighbors = 0;

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];

    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];

        if (sph_checksolid(d_type_property_map, d_pos[k].w)) continue;

        Scalar3 pj = make_scalar3(d_pos[k].x, d_pos[k].y, d_pos[k].z);
        Scalar3 dx = box.minImage(make_scalar3(pi.x - pj.x, pi.y - pj.y, pi.z - pj.z));
        Scalar  rsq = dot(dx, dx);

        if (kp.const_slength && rsq > kp.rcutsq) continue;

        Scalar r     = sqrtf(rsq);
        Scalar meanh = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar wij   = sph_wij<KT_>(kp.alpha, meanh, r);

        Scalar3 vj = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);
        Scalar  Pj = d_pressure[k];

        uf_c0.x += vj.x * wij;
        uf_c0.y += vj.y * wij;
        uf_c0.z += vj.z * wij;

        pf_c0 += Pj * wij;

        ph_c0.x += d_density[k] * dx.x * wij;
        ph_c0.y += d_density[k] * dx.y * wij;
        ph_c0.z += d_density[k] * dx.z * wij;

        wij_c0 += wij;
        fluid_neighbors++;
        }

    if (fluid_neighbors > 0 && wij_c0 > Scalar(0))
        {
        Scalar inv_w = Scalar(1.0) / wij_c0;

        d_vf[i].x = Scalar(2.0) * vi.x - inv_w * uf_c0.x;
        d_vf[i].y = Scalar(2.0) * vi.y - inv_w * uf_c0.y;
        d_vf[i].z = Scalar(2.0) * vi.z - inv_w * uf_c0.z;

        Scalar3 hp_factor = make_scalar3(bodyforce.x - accel_i.x,
                                         bodyforce.y - accel_i.y,
                                         bodyforce.z - accel_i.z);

        ph_c0.x *= inv_w;
        ph_c0.y *= inv_w;
        ph_c0.z *= inv_w;

        Scalar P_s = inv_w * pf_c0 + dot(hp_factor, ph_c0);

        if (P_s < Scalar(0))
            {
            d_pressure[i] = eos.bp;
            d_density[i]  = eos.rho0;
            }
        else
            {
            d_pressure[i] = P_s;
            d_density[i]  = sph_density_from_p<SET_>(eos.rho0, eos.c, eos.bp, P_s);
            }
        }
    else
        {
        d_vf[i].x = Scalar(0); d_vf[i].y = Scalar(0); d_vf[i].z = Scalar(0);
        d_pressure[i] = eos.bp;
        d_density[i]  = eos.rho0;
        }
    }

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
    unsigned int          block_size)
    {
    if (solid_group_size == 0) return hipSuccess;
    dim3 grid((solid_group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_noslip_kernel<KT_, SET_>),
                       grid, threads, 0, 0,
                       solid_group_size, d_solid_index,
                       d_pos, d_vel, d_density, d_pressure, d_vf, d_accel,
                       d_h, d_n_neigh, d_nlist, d_head_list,
                       d_type_property_map, box, kp, eos, bodyforce);
    return hipSuccess;
    }

// =========================================================================
// Kernel 4: Main force computation
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
__global__ void gpu_sph_forcecomputation_kernel(
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

    Scalar3 pi  = make_scalar3(d_pos[i].x, d_pos[i].y, d_pos[i].z);
    Scalar3 vi  = make_scalar3(d_vel[i].x, d_vel[i].y, d_vel[i].z);
    Scalar  mi  = d_vel[i].w;
    Scalar  Pi  = d_pressure[i];
    Scalar  rhoi = d_density[i];
    Scalar  Vi   = mi / rhoi;
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];

    Scalar4 fi   = make_scalar4(Scalar(0), Scalar(0), Scalar(0), Scalar(0));
    Scalar4 rdpe = make_scalar4(Scalar(0), Scalar(0), Scalar(0), Scalar(0));

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

        Scalar3 pj = make_scalar3(d_pos[k].x, d_pos[k].y, d_pos[k].z);
        Scalar3 dx = box.minImage(make_scalar3(pi.x - pj.x, pi.y - pj.y, pi.z - pj.z));
        Scalar  rsq = dot(dx, dx);

        if (kp.const_slength && rsq > kp.rcutsq) continue;

        Scalar mj   = d_vel[k].w;
        Scalar rhoj = d_density[k];
        Scalar Vj   = mj / rhoj;
        Scalar Pj   = d_pressure[k];

        Scalar3 vj;
        if (issolid)
            vj = make_scalar3(d_vf[k].x, d_vf[k].y, d_vf[k].z);
        else
            vj = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);

        Scalar3 dv = make_scalar3(vi.x - vj.x, vi.y - vj.y, vi.z - vj.z);
        Scalar  r  = sqrtf(rsq);

        Scalar meanh  = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar epssqr = Scalar(0.01) * meanh * meanh;

        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = dwdr / (r + epssqr);

        // Pressure force
        Scalar temp0;
        if (density_method == 0) // DENSITYSUMMATION
            temp0 = (Vi * Vi + Vj * Vj) * ((rhoj * Pi + rhoi * Pj) / (rhoi + rhoj));
        else // DENSITYCONTINUITY
            temp0 = mi * mj * (Pi + Pj) / (rhoi * rhoj);

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
                avc *= (density_method == 0) ? (Vi * Vi + Vj * Vj) : (mi * mj);
                }
            }

        fi.x -= (temp0 + avc) * dwdr_r * dx.x;
        fi.y -= (temp0 + avc) * dwdr_r * dx.y;
        fi.z -= (temp0 + avc) * dwdr_r * dx.z;

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
        temp0 = mu_eff * (Vi * Vi + Vj * Vj) * dwdr_r;
        fi.x += temp0 * dv.x;
        fi.y += temp0 * dv.y;
        fi.z += temp0 * dv.z;
        }

        // Density rate (DENSITYCONTINUITY)
        if (density_method == 1)
            {
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

    // Max-velocity reduction via atomicMax on IEEE 754 bit-cast
    Scalar vi_total = sqrtf(vi.x * vi.x + vi.y * vi.y + vi.z * vi.z);
    atomicMax(d_max_vel_bits, __float_as_uint(vi_total));
    }

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
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_forcecomputation_kernel<KT_, SET_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array,
                       d_pos, d_vel, d_density, d_pressure, d_vf, d_h,
                       d_force, d_ratedpe,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map,
                       d_max_vel_bits, box, kp, eos, nn,
                       density_method, artificial_viscosity, avalpha, avbeta,
                       tensil_correction, tensil_eps_pos, tensil_eps_neg,
                       density_diffusion, ddiff);
    return hipSuccess;
    }

// =========================================================================
// Kernel 5: Solid particle forces (reaction)
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
__global__ void gpu_sph_solid_forces_kernel(
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
    int                   density_method)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= solid_group_size) return;

    unsigned int i = d_solid_index[group_idx];

    Scalar3 pi  = make_scalar3(d_pos[i].x, d_pos[i].y, d_pos[i].z);
    Scalar3 vi  = make_scalar3(d_vel[i].x, d_vel[i].y, d_vel[i].z);
    Scalar  mi  = d_vel[i].w;
    Scalar  Pi  = d_pressure[i];
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

        Scalar3 pj = make_scalar3(d_pos[k].x, d_pos[k].y, d_pos[k].z);
        Scalar3 dx = box.minImage(make_scalar3(pi.x - pj.x, pi.y - pj.y, pi.z - pj.z));
        Scalar  rsq = dot(dx, dx);

        if (kp.const_slength && rsq > kp.rcutsq) continue;

        Scalar  mj   = d_vel[k].w;
        Scalar  rhoj = d_density[k];
        Scalar  Vj   = mj / rhoj;
        Scalar  Pj   = d_pressure[k];

        Scalar3 vj  = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);
        Scalar3 dv  = make_scalar3(vi.x - vj.x, vi.y - vj.y, vi.z - vj.z);
        Scalar  r   = sqrtf(rsq);

        Scalar meanh  = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar epssqr = Scalar(0.01) * meanh * meanh;

        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = dwdr / (r + epssqr);

        Scalar temp0;
        if (density_method == 0)
            temp0 = -(Vi * Vi + Vj * Vj) * ((rhoj * Pi + rhoi * Pj) / (rhoi + rhoj));
        else
            temp0 = -mi * mj * (Pi + Pj) / (rhoi * rhoj);

        Scalar ratio = mj / mi;
        fi.x -= ratio * temp0 * dwdr_r * dx.x;
        fi.y -= ratio * temp0 * dwdr_r * dx.y;
        fi.z -= ratio * temp0 * dwdr_r * dx.z;

        {
        Scalar dvnorm    = sqrtf(dot(dv, dv));
        Scalar gamma_dot = dvnorm / (r + sqrtf(epssqr));
        Scalar mu_eff    = sph_nn_viscosity(nn.mu, gamma_dot, nn.model,
                                            nn.K, nn.n, nn.mu0, nn.muinf,
                                            nn.lambda_NN, nn.tauy, nn.m_reg, nn.mu_min);
        temp0 = mu_eff * (Vi * Vi + Vj * Vj) * dwdr_r;
        fi.x -= ratio * temp0 * dv.x;
        fi.y -= ratio * temp0 * dv.y;
        fi.z -= ratio * temp0 * dv.z;
        }
        }

    d_force[i] = fi;
    }

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
    unsigned int          block_size)
    {
    if (solid_group_size == 0) return hipSuccess;
    dim3 grid((solid_group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_solid_forces_kernel<KT_, SET_>),
                       grid, threads, 0, 0,
                       solid_group_size, d_solid_index,
                       d_pos, d_vel, d_density, d_pressure, d_h,
                       d_force, d_n_neigh, d_nlist, d_head_list,
                       d_type_property_map, box, kp, nn, density_method);
    return hipSuccess;
    }

// =========================================================================
// Explicit template instantiations for all 10 <KT_, SET_> combinations
// =========================================================================

#define INSTANTIATE_SPH_GPU(KT, SET) \
    template hipError_t gpu_sph_ndensity<KT, SET>( \
        unsigned int, const unsigned int*, const Scalar4*, const Scalar4*, \
        Scalar*, const Scalar*, const unsigned int*, const unsigned int*, \
        const size_t*, BoxDim, SPHKernelDevParams, unsigned int); \
    template hipError_t gpu_sph_pressure<KT, SET>( \
        unsigned int, const unsigned int*, const Scalar*, Scalar*, \
        SPHEOSDevParams, unsigned int); \
    template hipError_t gpu_sph_noslip<KT, SET>( \
        unsigned int, const unsigned int*, const Scalar4*, const Scalar4*, \
        Scalar*, Scalar*, Scalar3*, const Scalar3*, const Scalar*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, BoxDim, SPHKernelDevParams, SPHEOSDevParams, \
        Scalar3, unsigned int); \
    template hipError_t gpu_sph_forcecomputation<KT, SET>( \
        unsigned int, const unsigned int*, const Scalar4*, const Scalar4*, \
        const Scalar*, const Scalar*, const Scalar3*, const Scalar*, \
        Scalar4*, Scalar4*, const unsigned int*, const unsigned int*, \
        const size_t*, const unsigned int*, uint32_t*, BoxDim, \
        SPHKernelDevParams, SPHEOSDevParams, SPHNNViscParams, \
        int, int, Scalar, Scalar, int, Scalar, Scalar, int, Scalar, unsigned int); \
    template hipError_t gpu_sph_solid_forces<KT, SET>( \
        unsigned int, const unsigned int*, const Scalar4*, const Scalar4*, \
        const Scalar*, const Scalar*, const Scalar*, Scalar4*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, BoxDim, SPHKernelDevParams, SPHNNViscParams, \
        int, unsigned int);

INSTANTIATE_SPH_GPU(wendlandc2, tait)
INSTANTIATE_SPH_GPU(wendlandc2, linear)
INSTANTIATE_SPH_GPU(wendlandc4, tait)
INSTANTIATE_SPH_GPU(wendlandc4, linear)
INSTANTIATE_SPH_GPU(wendlandc6, tait)
INSTANTIATE_SPH_GPU(wendlandc6, linear)
INSTANTIATE_SPH_GPU(quintic,    tait)
INSTANTIATE_SPH_GPU(quintic,    linear)
INSTANTIATE_SPH_GPU(cubicspline,tait)
INSTANTIATE_SPH_GPU(cubicspline,linear)

#undef INSTANTIATE_SPH_GPU

} // namespace kernel
} // namespace sph
} // namespace hoomd

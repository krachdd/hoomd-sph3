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

/*! \file TwoPhaseFlowCSFGPU.cu
    \brief GPU kernels for the two-phase interface machinery.

    Each kernel is one pass of TwoPhaseFlow::compute_colorgradients() /
    TwoPhaseFlow::compute_surfaceforce() (TwoPhaseFlow.cc), term by term:

      colorgradient_raw : raw colour gradients, all local particles
      normal_smooth     : Shepard-smoothed fluid normals (Adami et al. 2010)
      copy_group        : write the smoothed normals back
      wall_blend        : blended prescribed-contact-angle correction
      csf_pass1         : unit normals + reliability flag (local + ghosts)
      csf_pass2         : Morris-corrected curvature and CSF force
      wall_adhesion     : optional pairwise wetting force (SPH_BETA_ADH)
*/

#include "TwoPhaseFlowCSFGPU.cuh"
#include "SPHDeviceFunctions.cuh"

namespace hoomd
{
namespace sph
{
namespace kernel
{

static constexpr Scalar CSF_EPS_NORM = Scalar(1e-6);

__device__ __forceinline__ Scalar3 sph_scalar3(const Scalar4& v)
    {
    return make_scalar3(v.x, v.y, v.z);
    }

// =========================================================================
// Raw colour gradients (all local particles)
// =========================================================================

template<SmoothingKernelType KT_>
__global__ void gpu_sph_2pf_colorgradient_raw_kernel(
    unsigned int          N_local,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    Scalar3*              d_sn,
    Scalar3*              d_fn,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHCSFParams          cp)
    {
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N_local) return;

    Scalar4 posi = d_pos[i];
    Scalar3 pi   = sph_scalar3(posi);
    Scalar  mi   = d_vel[i].w;
    Scalar  rhoi = d_density[i];
    Scalar  Vi   = mi / rhoi;
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];

    bool i_issolid  = sph_checksolid (d_type_property_map, posi.w);
    bool i_isfluid1 = sph_checkfluid1(d_type_property_map, posi.w);
    bool i_isfluid2 = sph_checkfluid2(d_type_property_map, posi.w);

    Scalar3 sn = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    Scalar3 fn = make_scalar3(Scalar(0), Scalar(0), Scalar(0));

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];
    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];
        Scalar4 posk = d_pos[k];

        bool j_issolid  = sph_checksolid (d_type_property_map, posk.w);
        bool j_isfluid1 = sph_checkfluid1(d_type_property_map, posk.w);
        bool j_isfluid2 = sph_checkfluid2(d_type_property_map, posk.w);

        // same phase: no colour gradient
        if ((i_issolid && j_issolid) || (i_isfluid1 && j_isfluid1) || (i_isfluid2 && j_isfluid2))
            continue;

        Scalar3 dx  = box.minImage(make_scalar3(pi.x - posk.x, pi.y - posk.y, pi.z - posk.z));
        Scalar  rsq = dot(dx, dx);
        if (kp.const_slength && rsq > kp.rcutsq) continue;

        Scalar mj   = d_vel[k].w;
        Scalar rhoj = d_density[k];
        Scalar Vj   = mj / rhoj;
        Scalar r    = sqrt(rsq);

        Scalar meanh  = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = (r > Scalar(1e-8) * meanh) ? dwdr / r : Scalar(0);

        Scalar temp0;
        if (cp.cg_method == 0) // DENSITYRATIO (Adami)
            {
            temp0 = rhoi / (rhoi + rhoj) * (Vi * Vi + Vj * Vj) / Vi;
            // solid-fluid bands split 50/50 (wetting fix 2026-08-22)
            if (i_issolid || j_issolid)
                temp0 = Scalar(0.5) * (Vi * Vi + Vj * Vj) / Vi;
            }
        else // NUMBERDENSITY
            temp0 = Vj * Vj / Vi;

        Scalar c = temp0 * dwdr_r;
        if (i_issolid)
            {
            // per-phase solid-fluid gradients: fluid 1 -> sn, fluid 2 -> fn
            if (j_isfluid1)      { sn.x += c * dx.x; sn.y += c * dx.y; sn.z += c * dx.z; }
            else if (j_isfluid2) { fn.x += c * dx.x; fn.y += c * dx.y; fn.z += c * dx.z; }
            }
        else if (j_issolid)
            { sn.x += c * dx.x; sn.y += c * dx.y; sn.z += c * dx.z; }
        else
            { fn.x += c * dx.x; fn.y += c * dx.y; fn.z += c * dx.z; }
        }

    // orientation: solid -> fluid, fluid 1 -> fluid 2
    if (i_issolid)
        {
        sn = make_scalar3(-sn.x, -sn.y, -sn.z);
        fn = make_scalar3(-fn.x, -fn.y, -fn.z);
        }
    if (i_isfluid1)
        fn = make_scalar3(-fn.x, -fn.y, -fn.z);

    d_sn[i] = sn;
    d_fn[i] = fn;
    }

template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_colorgradient_raw(
    unsigned int N_local, const Scalar4* d_pos, const Scalar4* d_vel, const Scalar* d_density,
    const Scalar* d_h, Scalar3* d_sn, Scalar3* d_fn, const unsigned int* d_n_neigh,
    const unsigned int* d_nlist, const size_t* d_head_list, const unsigned int* d_type_property_map,
    BoxDim box, SPHKernelDevParams kp, SPHCSFParams cp, unsigned int block_size)
    {
    if (N_local == 0) return hipSuccess;
    block_size = sph_clamp_block_size((const void*)(gpu_sph_2pf_colorgradient_raw_kernel<KT_>), block_size);
    dim3 grid((N_local + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_colorgradient_raw_kernel<KT_>), grid, threads, 0, 0,
                       N_local, d_pos, d_vel, d_density, d_h, d_sn, d_fn,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map, box, kp, cp);
    return hipSuccess;
    }

// =========================================================================
// Shepard smoothing of the fluid normals (Adami et al. 2010)
// =========================================================================

template<SmoothingKernelType KT_>
__global__ void gpu_sph_2pf_normal_smooth_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    const Scalar3*        d_fn,
    Scalar3*              d_fn_smooth,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;
    unsigned int i = d_index_array[group_idx];

    Scalar3 fni = d_fn[i];
    Scalar norm_i = sqrt(dot(fni, fni));
    if (norm_i < CSF_EPS_NORM)
        {
        d_fn_smooth[i] = fni;
        return;
        }

    Scalar3 pi   = sph_scalar3(d_pos[i]);
    Scalar  Vi   = d_vel[i].w / d_density[i];
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];
    Scalar  w0_i = sph_w0<KT_>(kp.alpha, kp.self_density, hi);

    Scalar3 acc = make_scalar3(Vi * fni.x * w0_i, Vi * fni.y * w0_i, Vi * fni.z * w0_i);
    Scalar  w_acc = Vi * w0_i;

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];
    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];
        Scalar4 posk = d_pos[k];
        if (sph_checksolid(d_type_property_map, posk.w)) continue;   // fluid neighbours only

        Scalar3 dx  = box.minImage(make_scalar3(pi.x - posk.x, pi.y - posk.y, pi.z - posk.z));
        Scalar  rsq = dot(dx, dx);
        if (kp.const_slength && rsq > kp.rcutsq) continue;

        Scalar r     = sqrt(rsq);
        Scalar meanh = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar wij   = sph_wij<KT_>(kp.alpha, meanh, r);
        Scalar Vk    = d_vel[k].w / d_density[k];
        Scalar3 fnk  = d_fn[k];
        acc.x += Vk * fnk.x * wij;
        acc.y += Vk * fnk.y * wij;
        acc.z += Vk * fnk.z * wij;
        w_acc += Vk * wij;
        }

    if (w_acc > CSF_EPS_NORM)
        {
        Scalar inv_w = Scalar(1.0) / w_acc;
        d_fn_smooth[i] = make_scalar3(acc.x * inv_w, acc.y * inv_w, acc.z * inv_w);
        }
    else
        d_fn_smooth[i] = fni;
    }

template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_normal_smooth(
    unsigned int group_size, const unsigned int* d_index_array, const Scalar4* d_pos,
    const Scalar4* d_vel, const Scalar* d_density, const Scalar* d_h, const Scalar3* d_fn,
    Scalar3* d_fn_smooth, const unsigned int* d_n_neigh, const unsigned int* d_nlist,
    const size_t* d_head_list, const unsigned int* d_type_property_map, BoxDim box,
    SPHKernelDevParams kp, unsigned int block_size)
    {
    if (group_size == 0) return hipSuccess;
    block_size = sph_clamp_block_size((const void*)(gpu_sph_2pf_normal_smooth_kernel<KT_>), block_size);
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_normal_smooth_kernel<KT_>), grid, threads, 0, 0,
                       group_size, d_index_array, d_pos, d_vel, d_density, d_h, d_fn, d_fn_smooth,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map, box, kp);
    return hipSuccess;
    }

// =========================================================================
// Group scatter copy
// =========================================================================

__global__ void gpu_sph_2pf_copy_group_kernel(
    unsigned int group_size, const unsigned int* d_index_array,
    const Scalar3* d_src, Scalar3* d_dst)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;
    unsigned int i = d_index_array[group_idx];
    d_dst[i] = d_src[i];
    }

hipError_t gpu_sph_2pf_copy_group(
    unsigned int group_size, const unsigned int* d_index_array,
    const Scalar3* d_src, Scalar3* d_dst, unsigned int block_size)
    {
    if (group_size == 0) return hipSuccess;
    block_size = sph_clamp_block_size((const void*)gpu_sph_2pf_copy_group_kernel, block_size);
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL(gpu_sph_2pf_copy_group_kernel, grid, threads, 0, 0,
                       group_size, d_index_array, d_src, d_dst);
    return hipSuccess;
    }

// =========================================================================
// Blended prescribed-contact-angle wall correction (in place on fn)
// =========================================================================

template<SmoothingKernelType KT_>
__global__ void gpu_sph_2pf_wall_blend_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    const Scalar3*        d_sn,
    Scalar3*              d_fn,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHCSFParams          cp)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;
    unsigned int i = d_index_array[group_idx];

    Scalar3 sni = d_sn[i];
    Scalar normsn = sqrt(dot(sni, sni));
    if (normsn < CSF_EPS_NORM) return;
    Scalar3 fni = d_fn[i];
    Scalar normfn = sqrt(dot(fni, fni));
    if (normfn < CSF_EPS_NORM) return;

    // wall proximity weight from the solid kernel sum phi = sum_solid V_j W_ij
    Scalar3 pi = sph_scalar3(d_pos[i]);
    Scalar  hi = kp.const_slength ? kp.ch : d_h[i];
    Scalar phi_s = Scalar(0);
    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];
    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];
        Scalar4 posk = d_pos[k];
        if (!sph_checksolid(d_type_property_map, posk.w)) continue;
        Scalar3 dx  = box.minImage(make_scalar3(pi.x - posk.x, pi.y - posk.y, pi.z - posk.z));
        Scalar  rsq = dot(dx, dx);
        if (kp.const_slength && rsq > kp.rcutsq) continue;
        Scalar meanh = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        phi_s += (d_vel[k].w / d_density[k]) * sph_wij<KT_>(kp.alpha, meanh, sqrt(rsq));
        }
    Scalar wblend = Scalar(3) * phi_s;
    if (wblend > Scalar(1)) wblend = Scalar(1);
    if (wblend <= Scalar(0)) return;

    // nw = -sn/|sn| (into the wall); target n = -cos(omega) nw + sin(omega) t_hat
    Scalar3 nw = make_scalar3(-sni.x / normsn, -sni.y / normsn, -sni.z / normsn);
    Scalar fdotn = dot(fni, nw);
    Scalar3 t = make_scalar3(fni.x - fdotn * nw.x, fni.y - fdotn * nw.y, fni.z - fdotn * nw.z);
    Scalar normt = sqrt(dot(t, t));
    if (normt < Scalar(1e-3) * normfn) return;
    Scalar it = Scalar(1) / normt;
    Scalar3 ntgt = make_scalar3(-cp.cos_omega * nw.x + cp.sin_omega * t.x * it,
                                -cp.cos_omega * nw.y + cp.sin_omega * t.y * it,
                                -cp.cos_omega * nw.z + cp.sin_omega * t.z * it);
    Scalar wb1 = Scalar(1) - wblend;
    Scalar3 nblend = make_scalar3(wblend * ntgt.x + wb1 * fni.x / normfn,
                                  wblend * ntgt.y + wb1 * fni.y / normfn,
                                  wblend * ntgt.z + wb1 * fni.z / normfn);
    Scalar nb = sqrt(dot(nblend, nblend));
    if (nb < Scalar(1e-6)) return;
    d_fn[i] = make_scalar3(normfn * nblend.x / nb, normfn * nblend.y / nb, normfn * nblend.z / nb);
    }

template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_wall_blend(
    unsigned int group_size, const unsigned int* d_index_array, const Scalar4* d_pos,
    const Scalar4* d_vel, const Scalar* d_density, const Scalar* d_h, const Scalar3* d_sn,
    Scalar3* d_fn, const unsigned int* d_n_neigh, const unsigned int* d_nlist,
    const size_t* d_head_list, const unsigned int* d_type_property_map, BoxDim box,
    SPHKernelDevParams kp, SPHCSFParams cp, unsigned int block_size)
    {
    if (group_size == 0) return hipSuccess;
    block_size = sph_clamp_block_size((const void*)(gpu_sph_2pf_wall_blend_kernel<KT_>), block_size);
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_wall_blend_kernel<KT_>), grid, threads, 0, 0,
                       group_size, d_index_array, d_pos, d_vel, d_density, d_h, d_sn, d_fn,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map, box, kp, cp);
    return hipSuccess;
    }

// =========================================================================
// CSF pass 1: unit normals + reliability (local + ghost slots)
// =========================================================================

__global__ void gpu_sph_2pf_csf_pass1_kernel(
    unsigned int N_total, const Scalar4* d_pos, const Scalar* d_h, const Scalar3* d_fn,
    Scalar3* d_nhat, unsigned int* d_rel, const unsigned int* d_type_property_map,
    SPHKernelDevParams kp)
    {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N_total) return;
    d_nhat[idx] = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    d_rel[idx]  = 0;
    if (sph_checksolid(d_type_property_map, d_pos[idx].w)) return;
    Scalar3 fni = d_fn[idx];
    Scalar nf = sqrt(dot(fni, fni));
    Scalar hi = kp.const_slength ? kp.ch : d_h[idx];
    if (nf * hi < Scalar(0.01)) return;   // unreliable normal
    d_nhat[idx] = make_scalar3(fni.x / nf, fni.y / nf, fni.z / nf);
    d_rel[idx]  = 1;
    }

hipError_t gpu_sph_2pf_csf_pass1(
    unsigned int N_total, const Scalar4* d_pos, const Scalar* d_h, const Scalar3* d_fn,
    Scalar3* d_nhat, unsigned int* d_rel, const unsigned int* d_type_property_map,
    SPHKernelDevParams kp, unsigned int block_size)
    {
    if (N_total == 0) return hipSuccess;
    block_size = sph_clamp_block_size((const void*)gpu_sph_2pf_csf_pass1_kernel, block_size);
    dim3 grid((N_total + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL(gpu_sph_2pf_csf_pass1_kernel, grid, threads, 0, 0,
                       N_total, d_pos, d_h, d_fn, d_nhat, d_rel, d_type_property_map, kp);
    return hipSuccess;
    }

// =========================================================================
// CSF pass 2: curvature-form surface force
// =========================================================================

template<SmoothingKernelType KT_>
__global__ void gpu_sph_2pf_csf_pass2_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    const Scalar3*        d_fn,
    const Scalar3*        d_nhat,
    const unsigned int*   d_rel,
    Scalar3*              d_sf,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHCSFParams          cp)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;
    unsigned int i = d_index_array[group_idx];
    if (!d_rel[i]) return;

    Scalar3 fni  = d_fn[i];
    Scalar  nfi  = sqrt(dot(fni, fni));
    Scalar3 pi   = sph_scalar3(d_pos[i]);
    Scalar  Vi   = d_vel[i].w / d_density[i];
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];
    Scalar3 ni   = d_nhat[i];

    Scalar num = Scalar(0);   // sum (n_k - n_i) . gradW_ik V_k
    Scalar den = Scalar(0);   // sum W_ik V_k (support completeness)

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];
    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];
        if (!d_rel[k]) continue;
        Scalar4 posk = d_pos[k];
        Scalar3 dx  = box.minImage(make_scalar3(pi.x - posk.x, pi.y - posk.y, pi.z - posk.z));
        Scalar  rsq = dot(dx, dx);
        if (kp.const_slength && rsq > kp.rcutsq) continue;
        Scalar r     = sqrt(rsq);
        Scalar meanh = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        if (r < Scalar(1e-8) * meanh) continue;
        Scalar Vk     = d_vel[k].w / d_density[k];
        Scalar dwdr_r = sph_dwijdr<KT_>(kp.alpha, meanh, r) / r;
        Scalar3 nk    = d_nhat[k];
        num += Vk * dwdr_r * ((nk.x - ni.x) * dx.x + (nk.y - ni.y) * dx.y + (nk.z - ni.z) * dx.z);
        den += Vk * sph_wij<KT_>(kp.alpha, meanh, r);
        }
    den += Vi * sph_w0<KT_>(kp.alpha, kp.self_density, hi);   // self-contribution
    if (den < Scalar(0.1)) return;                             // too little reliable support

    Scalar kappa = -num / den;
    Scalar coef  = cp.sigma12 * kappa * nfi * Vi;
    Scalar3 sf = d_sf[i];
    d_sf[i] = make_scalar3(sf.x + coef * ni.x, sf.y + coef * ni.y, sf.z + coef * ni.z);
    }

template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_csf_pass2(
    unsigned int group_size, const unsigned int* d_index_array, const Scalar4* d_pos,
    const Scalar4* d_vel, const Scalar* d_density, const Scalar* d_h, const Scalar3* d_fn,
    const Scalar3* d_nhat, const unsigned int* d_rel, Scalar3* d_sf,
    const unsigned int* d_n_neigh, const unsigned int* d_nlist, const size_t* d_head_list,
    BoxDim box, SPHKernelDevParams kp, SPHCSFParams cp, unsigned int block_size)
    {
    if (group_size == 0) return hipSuccess;
    block_size = sph_clamp_block_size((const void*)(gpu_sph_2pf_csf_pass2_kernel<KT_>), block_size);
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_csf_pass2_kernel<KT_>), grid, threads, 0, 0,
                       group_size, d_index_array, d_pos, d_vel, d_density, d_h, d_fn, d_nhat, d_rel,
                       d_sf, d_n_neigh, d_nlist, d_head_list, box, kp, cp);
    return hipSuccess;
    }

// =========================================================================
// Optional pairwise wall-adhesion force (Tartakovsky-Meakin type)
// =========================================================================

template<SmoothingKernelType KT_>
__global__ void gpu_sph_2pf_wall_adhesion_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    Scalar3*              d_sf,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHCSFParams          cp)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;
    unsigned int i = d_index_array[group_idx];

    Scalar4 posi = d_pos[i];
    bool i_f1 = sph_checkfluid1(d_type_property_map, posi.w);
    bool i_f2 = sph_checkfluid2(d_type_property_map, posi.w);
    if (!i_f1 && !i_f2) return;

    Scalar3 pi = sph_scalar3(posi);
    Scalar  Vi = d_vel[i].w / d_density[i];
    Scalar  hi = kp.const_slength ? kp.ch : d_h[i];
    Scalar3 fadh = make_scalar3(Scalar(0), Scalar(0), Scalar(0));

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];
    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];
        Scalar4 posk = d_pos[k];
        if (!sph_checksolid(d_type_property_map, posk.w)) continue;
        Scalar3 dx  = box.minImage(make_scalar3(pi.x - posk.x, pi.y - posk.y, pi.z - posk.z));
        Scalar  rsq = dot(dx, dx);
        if (kp.const_slength && rsq > kp.rcutsq) continue;
        Scalar meanh = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar r = sqrt(rsq);
        if (r < Scalar(1e-8) * meanh) continue;
        Scalar wij = sph_wij<KT_>(kp.alpha, meanh, r);
        Scalar Vk  = d_vel[k].w / d_density[k];
        // differential form: favoured phase attracted, the other repelled
        Scalar s_a  = cp.beta_adh * cp.sigma12 / (meanh * meanh) * Scalar(0.5)
                      * (i_f1 ? cp.cos_omega : -cp.cos_omega);
        Scalar fmag = -s_a * Vi * Vk * wij / r;
        fadh.x += fmag * dx.x;
        fadh.y += fmag * dx.y;
        fadh.z += fmag * dx.z;
        }
    Scalar3 sf = d_sf[i];
    d_sf[i] = make_scalar3(sf.x + fadh.x, sf.y + fadh.y, sf.z + fadh.z);
    }

template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_wall_adhesion(
    unsigned int group_size, const unsigned int* d_index_array, const Scalar4* d_pos,
    const Scalar4* d_vel, const Scalar* d_density, const Scalar* d_h, Scalar3* d_sf,
    const unsigned int* d_n_neigh, const unsigned int* d_nlist, const size_t* d_head_list,
    const unsigned int* d_type_property_map, BoxDim box, SPHKernelDevParams kp,
    SPHCSFParams cp, unsigned int block_size)
    {
    if (group_size == 0) return hipSuccess;
    block_size = sph_clamp_block_size((const void*)(gpu_sph_2pf_wall_adhesion_kernel<KT_>), block_size);
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_wall_adhesion_kernel<KT_>), grid, threads, 0, 0,
                       group_size, d_index_array, d_pos, d_vel, d_density, d_h, d_sf,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map, box, kp, cp);
    return hipSuccess;
    }


// =========================================================================
// delta+-SPH particle shifting (Sun et al. 2017) -- mirrors
// TwoPhaseFlow::compute_particle_shift() pass by pass
// =========================================================================

template<SmoothingKernelType KT_>
__global__ void gpu_sph_2pf_shift_pass1_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    const Scalar3*        d_fn,
    Scalar3*              d_shift,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHShiftParams        sp)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;
    unsigned int i = d_index_array[group_idx];

    const Scalar eps      = Scalar(1e-10);
    const Scalar eps_norm = Scalar(1e-6);

    Scalar4 posi = d_pos[i];
    Scalar4 veli = d_vel[i];
    Scalar3 pi   = sph_scalar3(posi);
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];
    // W_ref: kernel at the approximate initial spacing dp ~ 0.5 h
    Scalar w_ref = sph_wij<KT_>(kp.alpha, hi, Scalar(0.5) * hi);
    if (w_ref < eps) w_ref = eps;

    Scalar3 grad_sum = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];
    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];
        // solid dummies participate (they close the kernel support at walls);
        // guard against marked-removed solids with zero density
        Scalar mk   = d_vel[k].w;
        Scalar rhok = d_density[k];
        if (rhok < Scalar(1e-12)) continue;
        Scalar hk   = kp.const_slength ? kp.ch : d_h[k];
        Scalar Vk   = mk / rhok;
        Scalar4 posk = d_pos[k];
        Scalar3 dx  = box.minImage(make_scalar3(pi.x - posk.x, pi.y - posk.y, pi.z - posk.z));
        Scalar  rsq = dot(dx, dx);
        if (rsq > kp.rcutsq) continue;
        Scalar r      = sqrt(rsq);
        Scalar meanh  = Scalar(0.5) * (hi + hk);
        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar wij    = sph_wij<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = (r > Scalar(1e-8) * meanh) ? dwdr / r : Scalar(0);
        // enhancement factor [1 + R (W_ij/W_ref)^n]
        Scalar ratio = wij / w_ref;
        Scalar Rpow  = Scalar(1);
        for (unsigned int q = 0; q < sp.n; q++) Rpow *= ratio;
        Scalar enhance = Scalar(1) + sp.R * Rpow;
        Scalar c = enhance * Vk * dwdr_r;
        grad_sum.x += c * dx.x;
        grad_sum.y += c * dx.y;
        grad_sum.z += c * dx.z;
        }

    // delta r_i = -A Ma_i (2 h_i)^2 sum_j [...] V_j grad W_ij
    Scalar ci    = sph_checkfluid1(d_type_property_map, posi.w) ? sp.c1 : sp.c2;
    Scalar vmagi = sqrt(veli.x * veli.x + veli.y * veli.y + veli.z * veli.z);
    Scalar coeff = -sp.A * (vmagi / ci) * Scalar(4.0) * hi * hi;
    Scalar3 dr = make_scalar3(coeff * grad_sum.x, coeff * grad_sum.y, coeff * grad_sum.z);

    // interface condition: remove the component normal to the fluid-fluid interface
    if (sp.interface_condition)
        {
        Scalar3 fn_i = d_fn[i];
        Scalar fn_mag = sqrt(dot(fn_i, fn_i));
        if (fn_mag > eps_norm)
            {
            Scalar inv = Scalar(1) / fn_mag;
            Scalar3 nh = make_scalar3(fn_i.x * inv, fn_i.y * inv, fn_i.z * inv);
            Scalar dr_n = dot(dr, nh);
            dr.x -= dr_n * nh.x; dr.y -= dr_n * nh.y; dr.z -= dr_n * nh.z;
            }
        }

    // NaN-safe magnitude cap at 0.1 h_i
    Scalar drmag2 = dot(dr, dr);
    const Scalar drcap = Scalar(0.1) * hi;
    if (!isfinite(drmag2))
        dr = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    else if (drmag2 > drcap * drcap)
        {
        Scalar rescale = drcap / sqrt(drmag2);
        dr.x *= rescale; dr.y *= rescale; dr.z *= rescale;
        }
    d_shift[i] = dr;
    }

template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_shift_pass1(
    unsigned int group_size, const unsigned int* d_index_array, const Scalar4* d_pos,
    const Scalar4* d_vel, const Scalar* d_density, const Scalar* d_h, const Scalar3* d_fn,
    Scalar3* d_shift, const unsigned int* d_n_neigh, const unsigned int* d_nlist,
    const size_t* d_head_list, const unsigned int* d_type_property_map, BoxDim box,
    SPHKernelDevParams kp, SPHShiftParams sp, unsigned int block_size)
    {
    if (group_size == 0) return hipSuccess;
    block_size = sph_clamp_block_size((const void*)(gpu_sph_2pf_shift_pass1_kernel<KT_>), block_size);
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_shift_pass1_kernel<KT_>), grid, threads, 0, 0,
                       group_size, d_index_array, d_pos, d_vel, d_density, d_h, d_fn, d_shift,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map, box, kp, sp);
    return hipSuccess;
    }

template<SmoothingKernelType KT_>
__global__ void gpu_sph_2pf_shift_pass2_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_h,
    const Scalar3*        d_shift,
    Scalar*               d_drho,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    BoxDim                box,
    SPHKernelDevParams    kp)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;
    unsigned int i = d_index_array[group_idx];

    Scalar3 pi   = sph_scalar3(d_pos[i]);
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];
    Scalar  rhoi = d_density[i];
    Scalar3 dri  = d_shift[i];
    Scalar delta_rho = Scalar(0);

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];
    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];
        Scalar mk   = d_vel[k].w;
        Scalar rhok = d_density[k];
        if (rhok < Scalar(1e-12)) continue;
        Scalar hk   = kp.const_slength ? kp.ch : d_h[k];
        Scalar Vk   = mk / rhok;
        Scalar4 posk = d_pos[k];
        Scalar3 dx  = box.minImage(make_scalar3(pi.x - posk.x, pi.y - posk.y, pi.z - posk.z));
        Scalar  rsq = dot(dx, dx);
        if (rsq > kp.rcutsq) continue;
        Scalar r      = sqrt(rsq);
        Scalar meanh  = Scalar(0.5) * (hi + hk);
        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = (r > Scalar(1e-8) * meanh) ? dwdr / r : Scalar(0);
        // ghost / solid slots carry a zero shift in d_shift
        Scalar3 drk = d_shift[k];
        Scalar3 ddr = make_scalar3(dri.x - drk.x, dri.y - drk.y, dri.z - drk.z);
        delta_rho += rhoi * Vk * dwdr_r * dot(ddr, dx);
        }
    d_drho[i] = delta_rho;
    }

template<SmoothingKernelType KT_>
hipError_t gpu_sph_2pf_shift_pass2(
    unsigned int group_size, const unsigned int* d_index_array, const Scalar4* d_pos,
    const Scalar4* d_vel, const Scalar* d_density, const Scalar* d_h, const Scalar3* d_shift,
    Scalar* d_drho, const unsigned int* d_n_neigh, const unsigned int* d_nlist,
    const size_t* d_head_list, BoxDim box, SPHKernelDevParams kp, unsigned int block_size)
    {
    if (group_size == 0) return hipSuccess;
    block_size = sph_clamp_block_size((const void*)(gpu_sph_2pf_shift_pass2_kernel<KT_>), block_size);
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_2pf_shift_pass2_kernel<KT_>), grid, threads, 0, 0,
                       group_size, d_index_array, d_pos, d_vel, d_density, d_h, d_shift, d_drho,
                       d_n_neigh, d_nlist, d_head_list, box, kp);
    return hipSuccess;
    }

__global__ void gpu_sph_2pf_shift_pass3_kernel(
    unsigned int group_size, const unsigned int* d_index_array, Scalar4* d_pos, int3* d_image,
    Scalar* d_density, const Scalar3* d_shift, const Scalar* d_drho, BoxDim local_box)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;
    unsigned int i = d_index_array[group_idx];
    Scalar4 p = d_pos[i];
    Scalar3 dr = d_shift[i];
    p.x += dr.x; p.y += dr.y; p.z += dr.z;
    int3 img = d_image[i];
    // local box: periodic only along non-decomposed directions (same as the integrators)
    local_box.wrap(p, img);
    d_pos[i]   = p;
    d_image[i] = img;
    if (d_drho) d_density[i] += d_drho[i];
    }

hipError_t gpu_sph_2pf_shift_pass3(
    unsigned int group_size, const unsigned int* d_index_array, Scalar4* d_pos, int3* d_image,
    Scalar* d_density, const Scalar3* d_shift, const Scalar* d_drho, BoxDim local_box,
    unsigned int block_size)
    {
    if (group_size == 0) return hipSuccess;
    block_size = sph_clamp_block_size((const void*)gpu_sph_2pf_shift_pass3_kernel, block_size);
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL(gpu_sph_2pf_shift_pass3_kernel, grid, threads, 0, 0,
                       group_size, d_index_array, d_pos, d_image, d_density, d_shift, d_drho, local_box);
    return hipSuccess;
    }

// =========================================================================
// Explicit instantiations (kernel type only)
// =========================================================================

#define INST_CSF_GPU(KT) \
    template hipError_t gpu_sph_2pf_colorgradient_raw<KT>( \
        unsigned int, const Scalar4*, const Scalar4*, const Scalar*, const Scalar*, \
        Scalar3*, Scalar3*, const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, BoxDim, SPHKernelDevParams, SPHCSFParams, unsigned int); \
    template hipError_t gpu_sph_2pf_normal_smooth<KT>( \
        unsigned int, const unsigned int*, const Scalar4*, const Scalar4*, const Scalar*, \
        const Scalar*, const Scalar3*, Scalar3*, const unsigned int*, const unsigned int*, \
        const size_t*, const unsigned int*, BoxDim, SPHKernelDevParams, unsigned int); \
    template hipError_t gpu_sph_2pf_wall_blend<KT>( \
        unsigned int, const unsigned int*, const Scalar4*, const Scalar4*, const Scalar*, \
        const Scalar*, const Scalar3*, Scalar3*, const unsigned int*, const unsigned int*, \
        const size_t*, const unsigned int*, BoxDim, SPHKernelDevParams, SPHCSFParams, unsigned int); \
    template hipError_t gpu_sph_2pf_csf_pass2<KT>( \
        unsigned int, const unsigned int*, const Scalar4*, const Scalar4*, const Scalar*, \
        const Scalar*, const Scalar3*, const Scalar3*, const unsigned int*, Scalar3*, \
        const unsigned int*, const unsigned int*, const size_t*, BoxDim, SPHKernelDevParams, \
        SPHCSFParams, unsigned int); \
    template hipError_t gpu_sph_2pf_wall_adhesion<KT>( \
        unsigned int, const unsigned int*, const Scalar4*, const Scalar4*, const Scalar*, \
        const Scalar*, Scalar3*, const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, BoxDim, SPHKernelDevParams, SPHCSFParams, unsigned int); \
    template hipError_t gpu_sph_2pf_shift_pass1<KT>( \
        unsigned int, const unsigned int*, const Scalar4*, const Scalar4*, const Scalar*, \
        const Scalar*, const Scalar3*, Scalar3*, const unsigned int*, const unsigned int*, \
        const size_t*, const unsigned int*, BoxDim, SPHKernelDevParams, SPHShiftParams, unsigned int); \
    template hipError_t gpu_sph_2pf_shift_pass2<KT>( \
        unsigned int, const unsigned int*, const Scalar4*, const Scalar4*, const Scalar*, \
        const Scalar*, const Scalar3*, Scalar*, const unsigned int*, const unsigned int*, \
        const size_t*, BoxDim, SPHKernelDevParams, unsigned int);

INST_CSF_GPU(wendlandc2)
INST_CSF_GPU(wendlandc4)
INST_CSF_GPU(wendlandc6)
INST_CSF_GPU(quintic)
INST_CSF_GPU(cubicspline)

#undef INST_CSF_GPU

} // namespace kernel
} // namespace sph
} // namespace hoomd

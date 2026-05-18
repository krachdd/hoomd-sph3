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

/*! \file SinglePhaseFlowFSGPU.cu
    \brief GPU kernels for the free-surface SPH force pipeline.

    Implements the four FS-specific compute stages as GPU kernels:
      1. gpu_sph_fs_detect_freesurface_kernel  — lambda + outward normals
      2. gpu_sph_fs_compute_curvature_kernel   — mean curvature
      3. gpu_sph_fs_pressure_clamp_kernel      — surface pressure floor
      4. gpu_sph_fs_forcecomputation_kernel    — TV forces + CSF + Young's
*/

#include "SinglePhaseFlowFSGPU.cuh"
#include "SPHDeviceFunctions.cuh"

namespace hoomd
{
namespace sph
{
namespace kernel
{

// =========================================================================
// Kernel 1: Detect free surface (lambda + outward normal)
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
__global__ void gpu_sph_fs_detect_freesurface_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    Scalar3*              d_fs_n,
    Scalar3*              d_aux4,
    const Scalar*         d_h,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHFSParams           fs)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;

    unsigned int i = d_index_array[group_idx];

    Scalar3 pi = make_scalar3(d_pos[i].x, d_pos[i].y, d_pos[i].z);
    Scalar  mi = d_vel[i].w;
    Scalar  hi = kp.const_slength ? kp.ch : d_h[i];
    Scalar  Vi = mi / fs.rho0_ref;

    // Self contribution W(0): same as sph_w0 (poly(0) * normfactor)
    Scalar lambda = Vi * sph_wij<KT_>(kp.alpha, hi, Scalar(0));

    Scalar3 grad_lambda = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    Scalar3 n_wall_acc  = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    bool    has_solid_nbr = false;

    const Scalar eps_sq = Scalar(1e-12);

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
        Scalar Vj   = mj / fs.rho0_ref;

        Scalar r     = sqrtf(rsq);
        Scalar meanh = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar epssqr = Scalar(0.01) * meanh * meanh;

        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = dwdr / (r + epssqr);

        if (!issolid)
            {
            Scalar wij = sph_wij<KT_>(kp.alpha, meanh, r);
            lambda        += Vj * wij;
            grad_lambda.x += Vj * dwdr_r * dx.x;
            grad_lambda.y += Vj * dwdr_r * dx.y;
            grad_lambda.z += Vj * dwdr_r * dx.z;
            }
        else
            {
            if (fs.apply_ca)
                {
                n_wall_acc.x += Vj * dwdr_r * dx.x;
                n_wall_acc.y += Vj * dwdr_r * dx.y;
                n_wall_acc.z += Vj * dwdr_r * dx.z;
                has_solid_nbr = true;
                }
            }
        }

    d_aux4[i].x = lambda;
    d_aux4[i].y = Scalar(0);  // kappa filled by compute_curvature
    d_aux4[i].z = Scalar(0);  // gnorm; set below for surface particles

    Scalar gnorm_sq = dot(grad_lambda, grad_lambda);

    if (lambda < fs.fs_threshold && gnorm_sq > eps_sq)
        {
        Scalar gnorm = sqrtf(gnorm_sq);
        d_aux4[i].z = gnorm;

        // Outward normal: -grad_lambda / gnorm  (grad_lambda points inward toward bulk)
        Scalar3 n_fs;
        n_fs.x = -grad_lambda.x / gnorm;
        n_fs.y = -grad_lambda.y / gnorm;
        n_fs.z = -grad_lambda.z / gnorm;

        // Contact-angle correction at the triple line (Huber 2016).
        if (has_solid_nbr)
            {
            Scalar nw_sq = dot(n_wall_acc, n_wall_acc);
            if (nw_sq > eps_sq)
                {
                Scalar nw_norm = sqrtf(nw_sq);
                // n_wall_acc = sum(Vj grad W) for solid neighbours: points FROM fluid INTO solid.
                // Negate to get wall normal pointing INTO the fluid.
                Scalar3 n_w;
                n_w.x = -n_wall_acc.x / nw_norm;
                n_w.y = -n_wall_acc.y / nw_norm;
                n_w.z = -n_wall_acc.z / nw_norm;

                // Tangential component of n_fs in the wall plane.
                Scalar n_dot_nw = dot(n_fs, n_w);
                Scalar3 n_t;
                n_t.x = n_fs.x - n_dot_nw * n_w.x;
                n_t.y = n_fs.y - n_dot_nw * n_w.y;
                n_t.z = n_fs.z - n_dot_nw * n_w.z;

                Scalar nt_sq = dot(n_t, n_t);
                if (nt_sq > eps_sq)
                    {
                    Scalar nt_norm = sqrtf(nt_sq);
                    n_t.x /= nt_norm;
                    n_t.y /= nt_norm;
                    n_t.z /= nt_norm;

                    // n_corrected = sin(theta) * t_wall + cos(theta) * n_wall
                    n_fs.x = fs.sin_ca * n_t.x + fs.cos_ca * n_w.x;
                    n_fs.y = fs.sin_ca * n_t.y + fs.cos_ca * n_w.y;
                    n_fs.z = fs.sin_ca * n_t.z + fs.cos_ca * n_w.z;

                    // Renormalise after blending.
                    Scalar nfs_norm = sqrtf(dot(n_fs, n_fs));
                    if (nfs_norm > eps_sq)
                        {
                        n_fs.x /= nfs_norm;
                        n_fs.y /= nfs_norm;
                        n_fs.z /= nfs_norm;
                        }
                    }
                }
            }

        d_fs_n[i] = n_fs;
        }
    else
        {
        // Bulk particle: zero normal (sentinel for subsequent kernels).
        d_fs_n[i] = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
        }
    }

template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_fs_detect_freesurface(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    Scalar3*              d_fs_n,
    Scalar3*              d_aux4,
    const Scalar*         d_h,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    SPHFSParams           fs,
    unsigned int          block_size)
    {
    if (group_size == 0) return hipSuccess;
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_fs_detect_freesurface_kernel<KT_, SET_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array, d_pos, d_vel,
                       d_fs_n, d_aux4, d_h,
                       d_n_neigh, d_nlist, d_head_list,
                       d_type_property_map, box, kp, fs);
    return hipSuccess;
    }

// =========================================================================
// Kernel 2: Compute curvature
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
__global__ void gpu_sph_fs_compute_curvature_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar3*        d_fs_n,
    Scalar3*              d_aux4,
    const Scalar*         d_h,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    Scalar                fs_threshold)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;

    unsigned int i = d_index_array[group_idx];

    // Skip bulk particles.
    Scalar lambda_i = d_aux4[i].x;
    if (lambda_i >= fs_threshold)
        {
        d_aux4[i].y = Scalar(0);
        return;
        }

    Scalar3 n_i = d_fs_n[i];
    const Scalar eps_sq = Scalar(1e-12);
    if (dot(n_i, n_i) < eps_sq)
        {
        d_aux4[i].y = Scalar(0);
        return;
        }

    Scalar3 pi  = make_scalar3(d_pos[i].x, d_pos[i].y, d_pos[i].z);
    Scalar  mi  = d_vel[i].w;
    Scalar  rhoi = d_density[i];
    Scalar  Vi   = mi / rhoi;
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];

    Scalar  kappa_nj = Scalar(0);
    Scalar3 grad_lam = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    bool    has_solid = false;

    size_t       myHead = d_head_list[i];
    unsigned int size   = d_n_neigh[i];

    for (unsigned int j = 0; j < size; j++)
        {
        unsigned int k = d_nlist[myHead + j];

        if (sph_checksolid(d_type_property_map, d_pos[k].w))
            {
            has_solid = true;
            continue;
            }

        Scalar3 pj  = make_scalar3(d_pos[k].x, d_pos[k].y, d_pos[k].z);
        Scalar3 dx  = box.minImage(make_scalar3(pi.x - pj.x, pi.y - pj.y, pi.z - pj.z));
        Scalar  rsq = dot(dx, dx);

        if (kp.const_slength && rsq > kp.rcutsq) continue;

        Scalar mj   = d_vel[k].w;
        Scalar rhoj = d_density[k];
        Scalar Vj   = mj / rhoj;

        Scalar3 n_j = d_fs_n[k];

        Scalar r     = sqrtf(rsq);
        Scalar meanh = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar epssqr = Scalar(0.01) * meanh * meanh;

        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = dwdr / (r + epssqr);

        Scalar3 gradW;
        gradW.x = dwdr_r * dx.x;
        gradW.y = dwdr_r * dx.y;
        gradW.z = dwdr_r * dx.z;

        kappa_nj   += Vj * dot(n_j, gradW);
        grad_lam.x += Vj * gradW.x;
        grad_lam.y += Vj * gradW.y;
        grad_lam.z += Vj * gradW.z;
        }

    // Contact-line particles (solid neighbour detected): set kappa = 0.
    // The CA-corrected normal from detect_freesurface is still propagated
    // to neighbouring surface particles via their kappa_nj accumulation.
    if (has_solid)
        {
        d_aux4[i].y = Scalar(0);
        return;
        }

    // kappa_i = (sum Vj n_j . gradW  -  n_i . sum Vj gradW) / V_i
    Scalar kappa = kappa_nj - dot(n_i, grad_lam);
    d_aux4[i].y = kappa / Vi;
    }

template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_fs_compute_curvature(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar3*        d_fs_n,
    Scalar3*              d_aux4,
    const Scalar*         d_h,
    const unsigned int*   d_n_neigh,
    const unsigned int*   d_nlist,
    const size_t*         d_head_list,
    const unsigned int*   d_type_property_map,
    BoxDim                box,
    SPHKernelDevParams    kp,
    Scalar                fs_threshold,
    unsigned int          block_size)
    {
    if (group_size == 0) return hipSuccess;
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_fs_compute_curvature_kernel<KT_, SET_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array, d_pos, d_vel, d_density,
                       d_fs_n, d_aux4, d_h,
                       d_n_neigh, d_nlist, d_head_list,
                       d_type_property_map, box, kp, fs_threshold);
    return hipSuccess;
    }

// =========================================================================
// Kernel 3: Pressure clamp at free surface (non-template)
// =========================================================================

__global__ void gpu_sph_fs_pressure_clamp_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    Scalar*               d_pressure,
    const Scalar3*        d_aux4,
    Scalar                fs_threshold)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;

    unsigned int i = d_index_array[group_idx];
    if (d_aux4[i].x < fs_threshold && d_pressure[i] < Scalar(0))
        d_pressure[i] = Scalar(0);
    }

hipError_t gpu_sph_fs_pressure_clamp(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    Scalar*               d_pressure,
    const Scalar3*        d_aux4,
    Scalar                fs_threshold,
    unsigned int          block_size)
    {
    if (group_size == 0) return hipSuccess;
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL(gpu_sph_fs_pressure_clamp_kernel,
                       grid, threads, 0, 0,
                       group_size, d_index_array, d_pressure, d_aux4, fs_threshold);
    return hipSuccess;
    }

// =========================================================================
// Kernel 4: Force computation — TV + CSF surface tension + Young's wetting
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
__global__ void gpu_sph_fs_forcecomputation_kernel(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar3*        d_tv,
    const Scalar3*        d_aux4,
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
    SPHFSParams           fs,
    int                   density_method,
    int                   artificial_viscosity,
    Scalar                avalpha,
    Scalar                avbeta,
    int                   density_diffusion,
    Scalar                ddiff)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size) return;

    unsigned int i = d_index_array[group_idx];

    // ── Read free-surface normal from aux2 BEFORE zeroing for BPC ─────────
    // detect_freesurface stored the outward unit normal in aux2.
    // We save it here, then overwrite aux2 with the BPC accumulation.
    Scalar3 n_fs_i = d_bpc[i];

    Scalar lambda_i = d_aux4[i].x;
    Scalar kappa_i  = d_aux4[i].y;
    const Scalar eps_nfs = Scalar(1e-12);
    bool is_surface = (lambda_i < fs.fs_threshold)
                      && (dot(n_fs_i, n_fs_i) > eps_nfs);

    // ── Particle i properties ──────────────────────────────────────────────
    Scalar3 pi   = make_scalar3(d_pos[i].x, d_pos[i].y, d_pos[i].z);
    Scalar3 vi   = make_scalar3(d_vel[i].x, d_vel[i].y, d_vel[i].z);
    Scalar  mi   = d_vel[i].w;
    Scalar  Pi   = d_pressure[i];
    Scalar  rhoi = d_density[i];
    Scalar  Vi   = mi / rhoi;
    Scalar  hi   = kp.const_slength ? kp.ch : d_h[i];

    Scalar3 tvi  = d_tv[i];

    // Artificial stress tensor for transport velocity correction.
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

    // Accumulators for Young's wetting force (contact-line surface particles).
    Scalar3 grad_lam_fy   = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    Scalar3 n_wall_acc_fy = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
    bool    has_solid_fy  = false;

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

        Scalar mj  = d_vel[k].w;

        // Skip solid particles marked for removal (mass set to -999 by
        // mark_solid_particles_toremove). Their huge Vj would blow up forces.
        if (mj < Scalar(0)) continue;

        Scalar3 vj;
        if (issolid)
            vj = make_scalar3(d_vf[k].x, d_vf[k].y, d_vf[k].z);
        else
            vj = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);

        Scalar rhoj = d_density[k];
        Scalar Vj   = mj / rhoj;
        Scalar Pj   = d_pressure[k];

        Scalar3 tvj = d_tv[k];

        // Artificial stress tensor for neighbour k.
        Scalar A11j = rhoj * vj.x * (tvj.x - vj.x);
        Scalar A12j = rhoj * vj.x * (tvj.y - vj.y);
        Scalar A13j = rhoj * vj.x * (tvj.z - vj.z);
        Scalar A21j = rhoj * vj.y * (tvj.x - vj.x);
        Scalar A22j = rhoj * vj.y * (tvj.y - vj.y);
        Scalar A23j = rhoj * vj.y * (tvj.z - vj.z);
        Scalar A31j = rhoj * vj.z * (tvj.x - vj.x);
        Scalar A32j = rhoj * vj.z * (tvj.y - vj.y);
        Scalar A33j = rhoj * vj.z * (tvj.z - vj.z);

        Scalar3 dv = make_scalar3(vi.x - vj.x, vi.y - vj.y, vi.z - vj.z);
        Scalar  r  = sqrtf(rsq);

        Scalar meanh  = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar epssqr = Scalar(0.01) * meanh * meanh;

        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = dwdr / (r + epssqr);

        Scalar vijsqr = Vi * Vi + Vj * Vj;

        // Accumulate Young's wetting force terms (done for every pair).
        if (!issolid)
            {
            grad_lam_fy.x += Vj * dwdr_r * dx.x;
            grad_lam_fy.y += Vj * dwdr_r * dx.y;
            grad_lam_fy.z += Vj * dwdr_r * dx.z;
            }
        else
            {
            n_wall_acc_fy.x += Vj * dwdr_r * dx.x;
            n_wall_acc_fy.y += Vj * dwdr_r * dx.y;
            n_wall_acc_fy.z += Vj * dwdr_r * dx.z;
            has_solid_fy = true;
            }

        // TV pressure (Adami 2013): P_bar = (rho_j * P_i + rho_i * P_j) / (rho_i + rho_j)
        Scalar temp0 = (rhoj * Pi + rhoi * Pj) / (rhoi + rhoj);

        // Adami pressure floor for solids below the fluid particle.
        // Prevents contact-line particles from being pulled into the floor.
        if (issolid && dx.y > Scalar(0) && temp0 < Scalar(200.0))
            temp0 = Scalar(200.0);

        // Artificial viscosity (Monaghan 1983): fluid-fluid only.
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

        // Pressure force.
        fi.x -= vijsqr * (temp0 + avc) * dwdr_r * dx.x;
        fi.y -= vijsqr * (temp0 + avc) * dwdr_r * dx.y;
        fi.z -= vijsqr * (temp0 + avc) * dwdr_r * dx.z;

        // Viscous force.
        {
        Scalar dvnorm    = sqrtf(dot(dv, dv));
        Scalar gamma_dot = dvnorm / (r + sqrtf(epssqr));
        Scalar mu_eff    = sph_nn_viscosity(nn.mu, gamma_dot, nn.model,
                                            nn.K, nn.n, nn.mu0, nn.muinf,
                                            nn.lambda_NN, nn.tauy, nn.m_reg, nn.mu_min);
        Scalar tv = mu_eff * vijsqr * dwdr_r;
        fi.x += tv * dv.x;
        fi.y += tv * dv.y;
        fi.z += tv * dv.z;
        }

        // Artificial stress tensor contribution (TV correction).
        {
        Scalar tv = Scalar(0.5) * vijsqr * dwdr_r;
        Scalar A1ij = (A11i + A11j) * dx.x + (A12i + A12j) * dx.y + (A13i + A13j) * dx.z;
        Scalar A2ij = (A21i + A21j) * dx.x + (A22i + A22j) * dx.y + (A23i + A23j) * dx.z;
        Scalar A3ij = (A31i + A31j) * dx.x + (A32i + A32j) * dx.y + (A33i + A33j) * dx.z;
        fi.x += tv * A1ij;
        fi.y += tv * A2ij;
        fi.z += tv * A3ij;
        }

        // Background pressure contribution to transport velocity (aux2 / BPC).
        bpc.x -= vijsqr * (Pb / mi) * dwdr_r * dx.x;
        bpc.y -= vijsqr * (Pb / mi) * dwdr_r * dx.y;
        bpc.z -= vijsqr * (Pb / mi) * dwdr_r * dx.z;

        // Density continuity rate (DENSITYCONTINUITY only).
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

    // dp/dt chain rule (DENSITYCONTINUITY).
    if (density_method == 1)
        rdpe.y = sph_dpressuredrho<SET_>(eos.rho0, eos.c, rhoi) * rdpe.x;

    // ── Young's-equation wetting force ────────────────────────────────────
    // Applied at contact-line free-surface particles (has_solid && is_surface
    // && sigma > 0 && contact angle != 90 degrees).
    if (has_solid_fy && is_surface && fs.sigma > Scalar(0) && fs.apply_ca)
        {
        const Scalar fy_eps = Scalar(1e-12);
        Scalar gl_sq = dot(grad_lam_fy, grad_lam_fy);
        Scalar nw_sq = dot(n_wall_acc_fy, n_wall_acc_fy);
        if (gl_sq > fy_eps && nw_sq > fy_eps)
            {
            // Natural (uncorrected) outward normal.
            Scalar gl_norm = sqrtf(gl_sq);
            Scalar3 n_nat;
            n_nat.x = -grad_lam_fy.x / gl_norm;
            n_nat.y = -grad_lam_fy.y / gl_norm;
            n_nat.z = -grad_lam_fy.z / gl_norm;

            // Wall outward normal (pointing into fluid).
            Scalar nw_norm = sqrtf(nw_sq);
            Scalar3 n_w;
            n_w.x = -n_wall_acc_fy.x / nw_norm;
            n_w.y = -n_wall_acc_fy.y / nw_norm;
            n_w.z = -n_wall_acc_fy.z / nw_norm;

            // Snap n_w to the nearest axis unit vector to remove asymmetric
            // solid-coverage artefacts that would drive floor penetration.
            Scalar ax = fabsf(n_w.x);
            Scalar ay = fabsf(n_w.y);
            Scalar az = fabsf(n_w.z);
            Scalar3 n_w_ax = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
            if (ax >= ay && ax >= az)
                n_w_ax.x = (n_w.x > Scalar(0)) ? Scalar(1) : Scalar(-1);
            else if (ay >= ax && ay >= az)
                n_w_ax.y = (n_w.y > Scalar(0)) ? Scalar(1) : Scalar(-1);
            else
                n_w_ax.z = (n_w.z > Scalar(0)) ? Scalar(1) : Scalar(-1);

            // Actual contact-angle cosine (snapped wall normal).
            Scalar cos_theta_act = dot(n_nat, n_w_ax);

            // Wall-tangential spreading direction t (lies exactly in wall plane).
            Scalar3 n_t;
            n_t.x = n_nat.x - cos_theta_act * n_w_ax.x;
            n_t.y = n_nat.y - cos_theta_act * n_w_ax.y;
            n_t.z = n_nat.z - cos_theta_act * n_w_ax.z;
            Scalar nt_sq = dot(n_t, n_t);
            if (nt_sq > fy_eps)
                {
                Scalar nt_norm = sqrtf(nt_sq);
                n_t.x /= nt_norm;
                n_t.y /= nt_norm;
                n_t.z /= nt_norm;

                // Contact-line arc length per particle ≈ Vi^(1/3)
                Scalar Vi_cbrt = cbrtf(Vi);

                // F_Y = sigma * (cos_ca + cos_theta_act) * Vi^(1/3) * t_hat
                // Equilibrium when cos_theta_act = -cos_ca -> theta_actual = theta_eq
                Scalar F_Y_mag = fs.sigma * (fs.cos_ca + cos_theta_act) * Vi_cbrt;
                fi.x += F_Y_mag * n_t.x;
                fi.y += F_Y_mag * n_t.y;
                fi.z += F_Y_mag * n_t.z;
                }
            }
        }

    // ── CSF surface tension force ──────────────────────────────────────────
    // F = -sigma * kappa_i * (mi/rhoi)^2 * gnorm_i * n_fs_i
    // (Hu & Adams 2006 corrected CSF formulation)
    if (is_surface && fs.sigma > Scalar(0))
        {
        Scalar gnorm_i = d_aux4[i].z;
        Scalar surf_coeff = -fs.sigma * kappa_i * (mi / rhoi) * (mi / rhoi) * gnorm_i;
        fi.x += surf_coeff * n_fs_i.x;
        fi.y += surf_coeff * n_fs_i.y;
        fi.z += surf_coeff * n_fs_i.z;
        }

    d_force[i]   = fi;
    d_ratedpe[i] = rdpe;
    d_bpc[i]     = bpc;

    // Max-velocity reduction via atomicMax on IEEE 754 bit-cast.
    Scalar vi_total = sqrtf(vi.x * vi.x + vi.y * vi.y + vi.z * vi.z);
    atomicMax(d_max_vel_bits, __float_as_uint(vi_total));
    }

template<SmoothingKernelType KT_, StateEquationType SET_>
hipError_t gpu_sph_fs_forcecomputation(
    unsigned int          group_size,
    const unsigned int*   d_index_array,
    const Scalar4*        d_pos,
    const Scalar4*        d_vel,
    const Scalar*         d_density,
    const Scalar*         d_pressure,
    const Scalar3*        d_vf,
    const Scalar3*        d_tv,
    const Scalar3*        d_aux4,
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
    SPHFSParams           fs,
    int                   density_method,
    int                   artificial_viscosity,
    Scalar                avalpha,
    Scalar                avbeta,
    int                   density_diffusion,
    Scalar                ddiff,
    unsigned int          block_size)
    {
    if (group_size == 0) return hipSuccess;
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_fs_forcecomputation_kernel<KT_, SET_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array,
                       d_pos, d_vel, d_density, d_pressure, d_vf, d_tv,
                       d_aux4, d_h,
                       d_force, d_ratedpe, d_bpc,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map,
                       d_max_vel_bits, box, kp, eos, nn, Pb, fs,
                       density_method, artificial_viscosity, avalpha, avbeta,
                       density_diffusion, ddiff);
    return hipSuccess;
    }

// =========================================================================
// Explicit template instantiations
// =========================================================================

#define INST_FS_GPU(KT, SET) \
    template hipError_t gpu_sph_fs_detect_freesurface<KT, SET>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, Scalar3*, Scalar3*, const Scalar*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, BoxDim, SPHKernelDevParams, SPHFSParams, unsigned int); \
    template hipError_t gpu_sph_fs_compute_curvature<KT, SET>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, const Scalar*, \
        const Scalar3*, Scalar3*, const Scalar*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, BoxDim, SPHKernelDevParams, Scalar, unsigned int); \
    template hipError_t gpu_sph_fs_forcecomputation<KT, SET>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, const Scalar*, const Scalar*, \
        const Scalar3*, const Scalar3*, const Scalar3*, const Scalar*, \
        Scalar4*, Scalar4*, Scalar3*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, uint32_t*, BoxDim, \
        SPHKernelDevParams, SPHEOSDevParams, SPHNNViscParams, \
        Scalar, SPHFSParams, int, int, Scalar, Scalar, int, Scalar, unsigned int);

INST_FS_GPU(wendlandc2, tait)
INST_FS_GPU(wendlandc2, linear)
INST_FS_GPU(wendlandc4, tait)
INST_FS_GPU(wendlandc4, linear)
INST_FS_GPU(wendlandc6, tait)
INST_FS_GPU(wendlandc6, linear)
INST_FS_GPU(quintic,    tait)
INST_FS_GPU(quintic,    linear)
INST_FS_GPU(cubicspline,tait)
INST_FS_GPU(cubicspline,linear)

#undef INST_FS_GPU

} // namespace kernel
} // namespace sph
} // namespace hoomd

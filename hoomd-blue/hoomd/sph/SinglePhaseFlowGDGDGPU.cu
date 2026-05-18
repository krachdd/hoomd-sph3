/* ---------------------------------------------------------
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.
Redistribution and use permitted under BSD 3-Clause License.
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*! \file SinglePhaseFlowGDGDGPU.cu
    \brief GPU kernel for gradient-density-gradient-driven (GDGD) SPH force computation.
*/

#include "SinglePhaseFlowGDGDGPU.cuh"
#include "SPHDeviceFunctions.cuh"

namespace hoomd
{
namespace sph
{
namespace kernel
{

template<SmoothingKernelType KT_, StateEquationType SET_>
__global__ void gpu_sph_gdgd_forcecomputation_kernel(
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
    Scalar  Ti   = d_aux4[i].x;  // scalar field (T or concentration)

    // Per-particle VRD rest density
    Scalar rho0_i = boussinesq ? rho0 : rho0 * (Scalar(1) - beta_s * (Ti - scalar_ref));

    // On-the-fly VRD pressure for particle i (SUMMATION + VRD only)
    Scalar Pi_use = Pi;
    if (!boussinesq && density_method == 0) // DENSITYSUMMATION
        Pi_use = sph_pressure_vrd<SET_>(eos.c, eos.bp, rho0_i, rhoi);

    Scalar4 fi   = make_scalar4(Scalar(0), Scalar(0), Scalar(0), Scalar(0));
    Scalar4 rdpe = make_scalar4(Scalar(0), Scalar(0), Scalar(0), Scalar(0));

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
        Scalar Tj   = d_aux4[k].x;

        Scalar3 vj;
        if (issolid)
            vj = make_scalar3(d_vf[k].x, d_vf[k].y, d_vf[k].z);
        else
            vj = make_scalar3(d_vel[k].x, d_vel[k].y, d_vel[k].z);

        // VRD pressure for neighbour (SUMMATION + VRD + fluid only)
        Scalar Pj_use = d_pressure[k];
        if (!boussinesq && density_method == 0 && !issolid)
            {
            Scalar rho0_k = rho0 * (Scalar(1) - beta_s * (Tj - scalar_ref));
            Pj_use = sph_pressure_vrd<SET_>(eos.c, eos.bp, rho0_k, rhoj);
            }

        Scalar3 dv  = make_scalar3(vi.x - vj.x, vi.y - vj.y, vi.z - vj.z);
        Scalar  r   = sqrtf(rsq);

        Scalar meanh  = kp.const_slength ? kp.ch : Scalar(0.5) * (hi + d_h[k]);
        Scalar epssqr = Scalar(0.01) * meanh * meanh;

        Scalar dwdr   = sph_dwijdr<KT_>(kp.alpha, meanh, r);
        Scalar dwdr_r = dwdr / (r + epssqr);

        Scalar vijsqr = Vi * Vi + Vj * Vj;

        // Pressure force
        Scalar temp0;
        if (density_method == 0) // DENSITYSUMMATION
            temp0 = vijsqr * (rhoj * Pi_use + rhoi * Pj_use) / (rhoi + rhoj);
        else // DENSITYCONTINUITY
            temp0 = mi * mj * (Pi_use + Pj_use) / (rhoi * rhoj);

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
                avc *= (density_method == 0) ? vijsqr : (mi * mj);
                }
            }

        fi.x -= (temp0 + avc) * dwdr_r * dx.x;
        fi.y -= (temp0 + avc) * dwdr_r * dx.y;
        fi.z -= (temp0 + avc) * dwdr_r * dx.z;

        // Viscous force
        {
        Scalar dvnorm    = sqrtf(dot(dv, dv));
        Scalar gamma_dot = dvnorm / (r + sqrtf(epssqr));
        Scalar mu_eff    = sph_nn_viscosity(nn.mu, gamma_dot, nn.model,
                                            nn.K, nn.n, nn.mu0, nn.muinf,
                                            nn.lambda_NN, nn.tauy, nn.m_reg, nn.mu_min);
        Scalar temp1 = mu_eff * vijsqr * dwdr_r;
        fi.x += temp1 * dv.x;
        fi.y += temp1 * dv.y;
        fi.z += temp1 * dv.z;
        }

        // Scalar diffusion (Morris-Fox-Zhu 1997): dT/dt += (kappa_s/Vi) * vijsqr * (Ti-Tj) * dwdr_r
        rdpe.z += kappa_s / Vi * vijsqr * (Ti - Tj) * dwdr_r;

        // Density continuity rate (DENSITYCONTINUITY)
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
        {
        Scalar dpdrho = boussinesq
                        ? sph_dpressuredrho<SET_>(eos.rho0, eos.c, rhoi)
                        : sph_dpressure_vrd_drho<SET_>(eos.c, rho0_i, rhoi);
        rdpe.y = dpdrho * rdpe.x;
        }

    // Boussinesq buoyancy correction: dF_b = mi * g * (-beta_s * (Ti - T_ref))
    if (boussinesq)
        {
        Scalar buoy = -beta_s * (Ti - scalar_ref);
        fi.x += mi * bodyforce.x * buoy;
        fi.y += mi * bodyforce.y * buoy;
        fi.z += mi * bodyforce.z * buoy;
        }

    d_force[i]   = fi;
    d_ratedpe[i] = rdpe;

    Scalar vi_total = sqrtf(vi.x * vi.x + vi.y * vi.y + vi.z * vi.z);
    atomicMax(d_max_vel_bits, __float_as_uint(vi_total));
    }

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
    unsigned int          block_size)
    {
    if (group_size == 0) return hipSuccess;
    dim3 grid((group_size + block_size - 1) / block_size, 1, 1);
    dim3 threads(block_size, 1, 1);
    hipLaunchKernelGGL((gpu_sph_gdgd_forcecomputation_kernel<KT_, SET_>),
                       grid, threads, 0, 0,
                       group_size, d_index_array,
                       d_pos, d_vel, d_density, d_pressure, d_vf, d_aux4, d_h,
                       d_force, d_ratedpe,
                       d_n_neigh, d_nlist, d_head_list, d_type_property_map,
                       d_max_vel_bits, box, kp, eos, nn,
                       kappa_s, beta_s, scalar_ref, rho0, boussinesq, bodyforce,
                       density_method, artificial_viscosity, avalpha, avbeta,
                       density_diffusion, ddiff);
    return hipSuccess;
    }

// =========================================================================
// Explicit template instantiations
// =========================================================================

#define INST_GDGD_GPU(KT, SET) \
    template hipError_t gpu_sph_gdgd_forcecomputation<KT, SET>( \
        unsigned int, const unsigned int*, \
        const Scalar4*, const Scalar4*, const Scalar*, const Scalar*, \
        const Scalar3*, const Scalar3*, const Scalar*, \
        Scalar4*, Scalar4*, \
        const unsigned int*, const unsigned int*, const size_t*, \
        const unsigned int*, uint32_t*, \
        BoxDim, SPHKernelDevParams, SPHEOSDevParams, SPHNNViscParams, \
        Scalar, Scalar, Scalar, Scalar, int, Scalar3, \
        int, int, Scalar, Scalar, int, Scalar, unsigned int);

INST_GDGD_GPU(wendlandc2, tait)
INST_GDGD_GPU(wendlandc2, linear)
INST_GDGD_GPU(wendlandc4, tait)
INST_GDGD_GPU(wendlandc4, linear)
INST_GDGD_GPU(wendlandc6, tait)
INST_GDGD_GPU(wendlandc6, linear)
INST_GDGD_GPU(quintic,    tait)
INST_GDGD_GPU(quintic,    linear)
INST_GDGD_GPU(cubicspline,tait)
INST_GDGD_GPU(cubicspline,linear)

#undef INST_GDGD_GPU

} // namespace kernel
} // namespace sph
} // namespace hoomd

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

/*! \file TwoPhaseFlowCSFGPU.h
    \brief Host-side driver of the GPU interface machinery (colour gradients,
           normal smoothing, wall blend, CSF surface force).

    Shared by TwoPhaseFlowGPU and TwoPhaseFlowTVGPU (which do not share a GPU
    base class), so the host orchestration lives here once and both classes'
    compute_colorgradients()/compute_surfaceforce() overrides just call it.
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __TWO_PHASE_FLOW_CSF_GPU_H__
#define __TWO_PHASE_FLOW_CSF_GPU_H__

#include "hoomd/GPUArray.h"
#include "hoomd/ExecutionConfiguration.h"
#include "hoomd/ParticleData.h"
#include "hoomd/ParticleGroup.h"
#include "hoomd/nsearch/NeighborList.h"
#include "TwoPhaseFlowCSFGPU.cuh"
#include "EvaluationMethodDefinition.h"

#include <cmath>
#include <cstdlib>
#include <memory>

namespace hoomd
{
namespace sph
{

/*! Device scratch buffers of the interface machinery, resized on demand. */
struct CSFGPUBuffers
    {
    GPUArray<Scalar3>      fn_smooth; //!< smoothed fluid normals (N_total)
    GPUArray<Scalar3>      nhat;      //!< unit normals of reliable slots (N_total)
    GPUArray<unsigned int> rel;       //!< reliability flags (N_total)
    GPUArray<Scalar3>      shift;     //!< particle shift vectors (N_total, ghosts zero)
    GPUArray<Scalar>       drho;      //!< ALE density correction (N_total)

    CSFGPUBuffers(std::shared_ptr<const ExecutionConfiguration> exec_conf)
        : fn_smooth(1, exec_conf), nhat(1, exec_conf), rel(1, exec_conf),
          shift(1, exec_conf), drho(1, exec_conf) { }

    void ensure(unsigned int n)
        {
        if (fn_smooth.getNumElements() < n) fn_smooth.resize(n);
        if (nhat.getNumElements()      < n) nhat.resize(n);
        if (rel.getNumElements()       < n) rel.resize(n);
        if (shift.getNumElements()     < n) shift.resize(n);
        if (drho.getNumElements()      < n) drho.resize(n);
        }
    };

/*! Free-function equivalent of the CHECK_CUDA_ERROR() macro (which needs this->m_exec_conf). */
inline void csf_check_cuda_error(const std::shared_ptr<const ExecutionConfiguration>& exec_conf)
    {
    hipError_t err_sync = hipPeekAtLastError();
    exec_conf->handleHIPError(err_sync, __FILE__, __LINE__);
    exec_conf->setDevice();
    hipError_t err_async = hipDeviceSynchronize();
    exec_conf->handleHIPError(err_async, __FILE__, __LINE__);
    }

/*! Build the interface parameters from the model state (mirrors the
 *  conditions in TwoPhaseFlow::compute_colorgradients/compute_surfaceforce).
 */
inline SPHCSFParams make_csfparams(ColorGradientMethod cg_method, Scalar omega, Scalar sigma12)
    {
    SPHCSFParams cp;
    cp.cg_method    = (cg_method == NUMBERDENSITY) ? 1 : 0;
    cp.blend_active = (omega != Scalar(90) && sigma12 > Scalar(0)) ? 1 : 0;
    cp.cos_omega    = Scalar(cos(omega * (M_PI / Scalar(180))));
    cp.sin_omega    = Scalar(sin(omega * (M_PI / Scalar(180))));
    cp.sigma12      = sigma12;
    // Calibration override SPH_BETA_ADH (default 0 = off since the
    // curvature-form CSF, 2026-08-23), same as the CPU path.
    static const Scalar beta_adh = getenv("SPH_BETA_ADH")
        ? Scalar(atof(getenv("SPH_BETA_ADH"))) : Scalar(0.0);
    cp.beta_adh        = beta_adh;
    cp.adhesion_active = (cp.blend_active && beta_adh != Scalar(0)) ? 1 : 0;
    return cp;
    }

/*! GPU equivalent of TwoPhaseFlow::compute_colorgradients():
 *  aux2 = solid normals, aux3 = smoothed + wall-blended fluid normals.
 */
template<SmoothingKernelType KT_>
void gpu_csf_compute_colorgradients(std::shared_ptr<ParticleData>          pdata,
                                    std::shared_ptr<nsearch::NeighborList> nlist,
                                    std::shared_ptr<ParticleGroup>         fluidgroup,
                                    const GPUArray<unsigned int>&          type_property_map,
                                    const SPHKernelDevParams&              kp,
                                    const SPHCSFParams&                    cp,
                                    CSFGPUBuffers&                         buf,
                                    unsigned int                           block_size,
                                    std::shared_ptr<const ExecutionConfiguration> exec_conf)
    {
    const bool error_check = exec_conf->isCUDAErrorCheckingEnabled();
    const BoxDim box = pdata->getGlobalBox();
    const unsigned int N_local = pdata->getN();
    const unsigned int N_total = N_local + pdata->getNGhosts();
    const unsigned int group_size = fluidgroup->getNumMembers();
    buf.ensure(N_total);

    ArrayHandle<Scalar3> d_sn(pdata->getAuxiliaries2(), access_location::device, access_mode::overwrite);
    ArrayHandle<Scalar3> d_fn(pdata->getAuxiliaries3(), access_location::device, access_mode::overwrite);
    ArrayHandle<Scalar4> d_pos(pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_vel(pdata->getVelocities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_density(pdata->getDensities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h(pdata->getSlengths(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_n_neigh(nlist->getNNeighArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(nlist->getNListArray(), access_location::device, access_mode::read);
    ArrayHandle<size_t>  d_head_list(nlist->getHeadList(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_index(fluidgroup->getIndexArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_type_map(type_property_map, access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_fn_smooth(buf.fn_smooth, access_location::device, access_mode::overwrite);

    // zero (incl. ghost slots, filled later by the ghost exchange)
    hipMemset(d_sn.data, 0, sizeof(Scalar3) * pdata->getAuxiliaries2().getNumElements());
    hipMemset(d_fn.data, 0, sizeof(Scalar3) * pdata->getAuxiliaries3().getNumElements());

    kernel::gpu_sph_2pf_colorgradient_raw<KT_>(
        N_local, d_pos.data, d_vel.data, d_density.data, d_h.data, d_sn.data, d_fn.data,
        d_n_neigh.data, d_nlist.data, d_head_list.data, d_type_map.data, box, kp, cp, block_size);
    if (error_check) csf_check_cuda_error(exec_conf);

    kernel::gpu_sph_2pf_normal_smooth<KT_>(
        group_size, d_index.data, d_pos.data, d_vel.data, d_density.data, d_h.data,
        d_fn.data, d_fn_smooth.data,
        d_n_neigh.data, d_nlist.data, d_head_list.data, d_type_map.data, box, kp, block_size);
    if (error_check) csf_check_cuda_error(exec_conf);

    kernel::gpu_sph_2pf_copy_group(group_size, d_index.data, d_fn_smooth.data, d_fn.data, block_size);
    if (error_check) csf_check_cuda_error(exec_conf);

    if (cp.blend_active)
        {
        kernel::gpu_sph_2pf_wall_blend<KT_>(
            group_size, d_index.data, d_pos.data, d_vel.data, d_density.data, d_h.data,
            d_sn.data, d_fn.data,
            d_n_neigh.data, d_nlist.data, d_head_list.data, d_type_map.data, box, kp, cp, block_size);
        if (error_check) csf_check_cuda_error(exec_conf);
        }
    }

/*! GPU equivalent of TwoPhaseFlow::compute_surfaceforce(): aux4 = surface
 *  force density (curvature-form CSF + optional wall adhesion).
 *  \pre aux2/aux3 hold ghost-synced normals.
 */
template<SmoothingKernelType KT_>
void gpu_csf_compute_surfaceforce(std::shared_ptr<ParticleData>          pdata,
                                  std::shared_ptr<nsearch::NeighborList> nlist,
                                  std::shared_ptr<ParticleGroup>         fluidgroup,
                                  const GPUArray<unsigned int>&          type_property_map,
                                  const SPHKernelDevParams&              kp,
                                  const SPHCSFParams&                    cp,
                                  CSFGPUBuffers&                         buf,
                                  unsigned int                           block_size,
                                  std::shared_ptr<const ExecutionConfiguration> exec_conf)
    {
    const bool error_check = exec_conf->isCUDAErrorCheckingEnabled();
    const BoxDim box = pdata->getGlobalBox();
    const unsigned int N_total = pdata->getN() + pdata->getNGhosts();
    const unsigned int group_size = fluidgroup->getNumMembers();
    buf.ensure(N_total);

    ArrayHandle<Scalar3> d_sf(pdata->getAuxiliaries4(), access_location::device, access_mode::overwrite);
    ArrayHandle<Scalar3> d_fn(pdata->getAuxiliaries3(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_pos(pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_vel(pdata->getVelocities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_density(pdata->getDensities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h(pdata->getSlengths(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_n_neigh(nlist->getNNeighArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(nlist->getNListArray(), access_location::device, access_mode::read);
    ArrayHandle<size_t>  d_head_list(nlist->getHeadList(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_index(fluidgroup->getIndexArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_type_map(type_property_map, access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_nhat(buf.nhat, access_location::device, access_mode::overwrite);
    ArrayHandle<unsigned int> d_rel(buf.rel, access_location::device, access_mode::overwrite);

    hipMemset(d_sf.data, 0, sizeof(Scalar3) * pdata->getAuxiliaries4().getNumElements());

    kernel::gpu_sph_2pf_csf_pass1(
        N_total, d_pos.data, d_h.data, d_fn.data, d_nhat.data, d_rel.data, d_type_map.data,
        kp, block_size);
    if (error_check) csf_check_cuda_error(exec_conf);

    kernel::gpu_sph_2pf_csf_pass2<KT_>(
        group_size, d_index.data, d_pos.data, d_vel.data, d_density.data, d_h.data,
        d_fn.data, d_nhat.data, d_rel.data, d_sf.data,
        d_n_neigh.data, d_nlist.data, d_head_list.data, box, kp, cp, block_size);
    if (error_check) csf_check_cuda_error(exec_conf);

    if (cp.adhesion_active)
        {
        kernel::gpu_sph_2pf_wall_adhesion<KT_>(
            group_size, d_index.data, d_pos.data, d_vel.data, d_density.data, d_h.data, d_sf.data,
            d_n_neigh.data, d_nlist.data, d_head_list.data, d_type_map.data, box, kp, cp, block_size);
        if (error_check) csf_check_cuda_error(exec_conf);
        }
    }

/*! GPU equivalent of TwoPhaseFlow::compute_particle_shift() (delta+-SPH,
 *  Sun et al. 2017): pass 1 shift vectors, pass 2 ALE density correction
 *  (DENSITYCONTINUITY, evaluated from the pre-shift densities), pass 3 apply
 *  and wrap with the local box.
 *  \pre aux3 holds ghost-synced fluid-fluid normals.
 */
template<SmoothingKernelType KT_>
void gpu_csf_compute_particle_shift(std::shared_ptr<ParticleData>          pdata,
                                    std::shared_ptr<nsearch::NeighborList> nlist,
                                    std::shared_ptr<ParticleGroup>         fluidgroup,
                                    const GPUArray<unsigned int>&          type_property_map,
                                    const SPHKernelDevParams&              kp,
                                    const SPHShiftParams&                  sp,
                                    bool                                   density_continuity,
                                    CSFGPUBuffers&                         buf,
                                    unsigned int                           block_size,
                                    std::shared_ptr<const ExecutionConfiguration> exec_conf)
    {
    const bool error_check = exec_conf->isCUDAErrorCheckingEnabled();
    const BoxDim box       = pdata->getGlobalBox();
    const BoxDim local_box = pdata->getBox();
    const unsigned int N_total = pdata->getN() + pdata->getNGhosts();
    const unsigned int group_size = fluidgroup->getNumMembers();
    buf.ensure(N_total);

    { // pass 1 + 2: read-only particle data
    ArrayHandle<Scalar4> d_pos(pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_vel(pdata->getVelocities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_density(pdata->getDensities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h(pdata->getSlengths(), access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_fn(pdata->getAuxiliaries3(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_n_neigh(nlist->getNNeighArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(nlist->getNListArray(), access_location::device, access_mode::read);
    ArrayHandle<size_t>  d_head_list(nlist->getHeadList(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_index(fluidgroup->getIndexArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_type_map(type_property_map, access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_shift(buf.shift, access_location::device, access_mode::overwrite);
    ArrayHandle<Scalar>  d_drho(buf.drho, access_location::device, access_mode::overwrite);

    hipMemset(d_shift.data, 0, sizeof(Scalar3) * buf.shift.getNumElements());

    kernel::gpu_sph_2pf_shift_pass1<KT_>(
        group_size, d_index.data, d_pos.data, d_vel.data, d_density.data, d_h.data, d_fn.data,
        d_shift.data, d_n_neigh.data, d_nlist.data, d_head_list.data, d_type_map.data,
        box, kp, sp, block_size);
    if (error_check) csf_check_cuda_error(exec_conf);

    if (density_continuity)
        {
        kernel::gpu_sph_2pf_shift_pass2<KT_>(
            group_size, d_index.data, d_pos.data, d_vel.data, d_density.data, d_h.data,
            d_shift.data, d_drho.data, d_n_neigh.data, d_nlist.data, d_head_list.data,
            box, kp, block_size);
        if (error_check) csf_check_cuda_error(exec_conf);
        }
    }

    { // pass 3: apply
    ArrayHandle<Scalar4> d_pos(pdata->getPositions(), access_location::device, access_mode::readwrite);
    ArrayHandle<int3>    d_image(pdata->getImages(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar>  d_density(pdata->getDensities(), access_location::device, access_mode::readwrite);
    ArrayHandle<unsigned int> d_index(fluidgroup->getIndexArray(), access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_shift(buf.shift, access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_drho(buf.drho, access_location::device, access_mode::read);
    kernel::gpu_sph_2pf_shift_pass3(
        group_size, d_index.data, d_pos.data, d_image.data, d_density.data, d_shift.data,
        density_continuity ? d_drho.data : nullptr, local_box, block_size);
    if (error_check) csf_check_cuda_error(exec_conf);
    }
    }

} // namespace sph
} // namespace hoomd

#endif // __TWO_PHASE_FLOW_CSF_GPU_H__

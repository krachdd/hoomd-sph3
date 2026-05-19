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

/*! \file TwoPhaseFlowGPU.cc
    \brief Host-side implementation of TwoPhaseFlowGPU.
*/

#include "TwoPhaseFlowGPU.h"
#include "SinglePhaseFlowGPU.cuh"
#include "TwoPhaseFlowGPU.cuh"

#include <cmath>
#include <stdexcept>
#include <pybind11/pybind11.h>

using namespace std;

namespace hoomd
{
namespace sph
{

// =========================================================================
// Constructor
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
TwoPhaseFlowGPU<KT_, SET1_, SET2_>::TwoPhaseFlowGPU(
    std::shared_ptr<SystemDefinition>      sysdef,
    std::shared_ptr<SmoothingKernel<KT_>>  skernel,
    std::shared_ptr<StateEquation<SET1_>>  equationofstate1,
    std::shared_ptr<StateEquation<SET2_>>  equationofstate2,
    std::shared_ptr<nsearch::NeighborList> nlist,
    std::shared_ptr<ParticleGroup>         fluidgroup1,
    std::shared_ptr<ParticleGroup>         fluidgroup2,
    std::shared_ptr<ParticleGroup>         solidgroup,
    DensityMethod       mdensitymethod,
    ViscosityMethod     mviscositymethod,
    ColorGradientMethod mcolorgradientmethod)
    : TwoPhaseFlow<KT_, SET1_, SET2_>(sysdef, skernel, equationofstate1, equationofstate2,
                                      nlist, fluidgroup1, fluidgroup2, solidgroup,
                                      mdensitymethod, mviscositymethod, mcolorgradientmethod),
      m_max_vel_bits(1, this->m_exec_conf),
      m_max_vel_bits_host(nullptr)
    {
    hipHostMalloc(&m_max_vel_bits_host, sizeof(uint32_t), hipHostMallocDefault);
    *m_max_vel_bits_host = 0;
    if (!this->m_exec_conf->isCUDAEnabled())
        {
        this->m_exec_conf->msg->error()
            << "Creating a TwoPhaseFlowGPU without a GPU in the execution configuration"
            << endl;
        throw std::runtime_error("Error initializing TwoPhaseFlowGPU");
        }

    m_tuner_ndensity  = std::make_shared<Autotuner<1>>(std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(this->m_exec_conf)}, this->m_exec_conf, "sph_2pf_ndensity");
    m_tuner_pressure  = std::make_shared<Autotuner<1>>(std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(this->m_exec_conf)}, this->m_exec_conf, "sph_2pf_pressure");
    m_tuner_noslip    = std::make_shared<Autotuner<1>>(std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(this->m_exec_conf)}, this->m_exec_conf, "sph_2pf_noslip");
    // Force kernel uses tpp=4 threads per particle; block_size is clamped dynamically
    // inside the wrapper via hipFuncGetAttributes.  The Autotuner scans all warp multiples
    // up to 512; the wrapper will silently reduce any over-limit value.
    m_tuner_force     = std::make_shared<Autotuner<1>>(std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(this->m_exec_conf)}, this->m_exec_conf, "sph_2pf_force", 5, false, [](const std::array<unsigned int, 1>& p) { return p[0] <= 512; });
    m_tuner_solidforce= std::make_shared<Autotuner<1>>(std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(this->m_exec_conf)}, this->m_exec_conf, "sph_2pf_solidforce");

    this->m_autotuners.push_back(m_tuner_ndensity);
    this->m_autotuners.push_back(m_tuner_pressure);
    this->m_autotuners.push_back(m_tuner_noslip);
    this->m_autotuners.push_back(m_tuner_force);
    this->m_autotuners.push_back(m_tuner_solidforce);
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
TwoPhaseFlowGPU<KT_, SET1_, SET2_>::~TwoPhaseFlowGPU()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying TwoPhaseFlowGPU" << endl;
    if (m_max_vel_bits_host) hipHostFree(m_max_vel_bits_host);
    }

// =========================================================================
// Helper: build parameter structs
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
SPHKernelDevParams TwoPhaseFlowGPU<KT_, SET1_, SET2_>::make_kparams() const
    {
    SPHKernelDevParams kp;
    kp.alpha         = this->m_skernel->getAlpha();
    kp.self_density  = this->m_skernel->getSelfDensity();
    kp.ch            = this->m_ch;
    kp.rcutsq        = this->m_rcutsq;
    kp.const_slength = this->m_const_slength ? 1 : 0;
    kp.kappa         = this->m_kappa;
    return kp;
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
SPHEOSDevParams TwoPhaseFlowGPU<KT_, SET1_, SET2_>::make_eos1params() const
    {
    SPHEOSDevParams eos;
    eos.rho0 = this->m_eos1->getRestDensity();
    eos.c    = this->m_eos1->getSpeedOfSound();
    eos.bp   = this->m_eos1->getBackgroundPressure();
    return eos;
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
SPHEOSDevParams TwoPhaseFlowGPU<KT_, SET1_, SET2_>::make_eos2params() const
    {
    SPHEOSDevParams eos;
    eos.rho0 = this->m_eos2->getRestDensity();
    eos.c    = this->m_eos2->getSpeedOfSound();
    eos.bp   = this->m_eos2->getBackgroundPressure();
    return eos;
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
SPHNNViscParams TwoPhaseFlowGPU<KT_, SET1_, SET2_>::make_nn1params() const
    {
    SPHNNViscParams nn;
    nn.mu        = this->m_mu1;
    nn.K         = this->m_nn_K1;
    nn.n         = this->m_nn_n1;
    nn.mu0       = this->m_nn_mu0_1;
    nn.muinf     = this->m_nn_muinf_1;
    nn.lambda_NN = this->m_nn_lambda1;
    nn.tauy      = this->m_nn_tauy1;
    nn.m_reg     = this->m_nn_m1;
    nn.mu_min    = this->m_nn_mu_min1;
    nn.model     = static_cast<int>(this->m_nn_model1);
    return nn;
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
SPHNNViscParams TwoPhaseFlowGPU<KT_, SET1_, SET2_>::make_nn2params() const
    {
    SPHNNViscParams nn;
    nn.mu        = this->m_mu2;
    nn.K         = this->m_nn_K2;
    nn.n         = this->m_nn_n2;
    nn.mu0       = this->m_nn_mu0_2;
    nn.muinf     = this->m_nn_muinf_2;
    nn.lambda_NN = this->m_nn_lambda2;
    nn.tauy      = this->m_nn_tauy2;
    nn.m_reg     = this->m_nn_m2;
    nn.mu_min    = this->m_nn_mu_min2;
    nn.model     = static_cast<int>(this->m_nn_model2);
    return nn;
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
SPHTwoPhaseParams TwoPhaseFlowGPU<KT_, SET1_, SET2_>::make_2pfparams(uint64_t timestep)
    {
    SPHTwoPhaseParams fp;
    fp.cmax                = this->m_cmax;
    fp.avalpha             = this->m_avalpha;
    fp.avbeta              = this->m_avbeta;
    fp.riemann_beta        = this->m_riemann_beta;
    fp.ddiff               = this->m_ddiff;
    Scalar3 gvec           = this->getAcceleration(timestep);
    fp.gvec_x              = gvec.x;
    fp.gvec_y              = gvec.y;
    fp.gvec_z              = gvec.z;
    fp.density_method      = static_cast<int>(this->m_density_method);
    fp.artificial_viscosity= this->m_artificial_viscosity ? 1 : 0;
    fp.riemann_dissipation = this->m_riemann_dissipation  ? 1 : 0;
    fp.cip                 = this->m_consistent_interface_pressure ? 1 : 0;
    fp.density_diffusion   = this->m_density_diffusion ? 1 : 0;
    return fp;
    }

// =========================================================================
// compute_ndensity (shared GPU kernel — uses fluid union group)
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::compute_ndensity(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "TwoPhaseFlowGPU::compute_ndensity" << endl;

    const BoxDim box = this->m_pdata->getGlobalBox();
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    if (group_size == 0) return;

    SPHKernelDevParams kp = make_kparams();

    { // device scope
    ArrayHandle<Scalar>  d_density(this->m_pdata->getDensities(),
                                    access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_pos(this->m_pdata->getPositions(),
                                access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_vel(this->m_pdata->getVelocities(),
                                access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h(this->m_pdata->getSlengths(),
                              access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(),
                                         access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(this->m_nlist->getNListArray(),
                                       access_location::device, access_mode::read);
    ArrayHandle<size_t>  d_head_list(this->m_nlist->getHeadList(),
                                      access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_index_array(this->m_fluidgroup->getIndexArray(),
                                             access_location::device, access_mode::read);

    m_tuner_ndensity->begin();
    kernel::gpu_sph_ndensity<KT_, SET1_>(
        group_size, d_index_array.data,
        d_pos.data, d_vel.data, d_density.data, d_h.data,
        d_n_neigh.data, d_nlist.data, d_head_list.data,
        box, kp,
        m_tuner_ndensity->getParam()[0]);
    if (this->m_exec_conf->isCUDAErrorCheckingEnabled()) CHECK_CUDA_ERROR();
    m_tuner_ndensity->end();
    }
    }

// =========================================================================
// compute_pressure (per-phase EOS — fluidgroup1 with eos1, fluidgroup2 with eos2)
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::compute_pressure(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "TwoPhaseFlowGPU::compute_pressure" << endl;

    SPHEOSDevParams eos1 = make_eos1params();
    SPHEOSDevParams eos2 = make_eos2params();

    // Phase 1
    {
    unsigned int group_size1 = this->m_fluidgroup1->getNumMembers();
    if (group_size1 > 0)
        {
        ArrayHandle<Scalar>  d_density (this->m_pdata->getDensities(),
                                         access_location::device, access_mode::read);
        ArrayHandle<Scalar>  d_pressure(this->m_pdata->getPressures(),
                                         access_location::device, access_mode::readwrite);
        ArrayHandle<unsigned int> d_index_array(this->m_fluidgroup1->getIndexArray(),
                                                 access_location::device, access_mode::read);

        m_tuner_pressure->begin();
        kernel::gpu_sph_pressure<KT_, SET1_>(
            group_size1, d_index_array.data,
            d_density.data, d_pressure.data, eos1,
            m_tuner_pressure->getParam()[0]);
        if (this->m_exec_conf->isCUDAErrorCheckingEnabled()) CHECK_CUDA_ERROR();
        m_tuner_pressure->end();
        }
    }

    // Phase 2
    {
    unsigned int group_size2 = this->m_fluidgroup2->getNumMembers();
    if (group_size2 > 0)
        {
        ArrayHandle<Scalar>  d_density (this->m_pdata->getDensities(),
                                         access_location::device, access_mode::read);
        ArrayHandle<Scalar>  d_pressure(this->m_pdata->getPressures(),
                                         access_location::device, access_mode::readwrite);
        ArrayHandle<unsigned int> d_index_array(this->m_fluidgroup2->getIndexArray(),
                                                 access_location::device, access_mode::read);

        m_tuner_pressure->begin();
        kernel::gpu_sph_pressure<KT_, SET2_>(
            group_size2, d_index_array.data,
            d_density.data, d_pressure.data, eos2,
            m_tuner_pressure->getParam()[0]);
        if (this->m_exec_conf->isCUDAErrorCheckingEnabled()) CHECK_CUDA_ERROR();
        m_tuner_pressure->end();
        }
    }
    }

// =========================================================================
// compute_noslip (shared GPU kernel with eos1 parameters)
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::compute_noslip(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "TwoPhaseFlowGPU::compute_noslip" << endl;

    unsigned int solid_group_size = this->m_solidgroup->getNumMembers();
    if (solid_group_size == 0) return;

    const BoxDim box = this->m_pdata->getGlobalBox();
    SPHKernelDevParams kp  = make_kparams();
    SPHEOSDevParams    eos = make_eos1params();

    Scalar3 bodyforce = this->getAcceleration(timestep);

    { // device scope
    ArrayHandle<Scalar>  d_density(this->m_pdata->getDensities(),
                                    access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar>  d_pressure(this->m_pdata->getPressures(),
                                     access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3> d_vf(this->m_pdata->getAuxiliaries1(),
                               access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_pos(this->m_pdata->getPositions(),
                                access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_vel(this->m_pdata->getVelocities(),
                                access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_accel(this->m_pdata->getAccelerations(),
                                  access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h(this->m_pdata->getSlengths(),
                              access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(),
                                         access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(this->m_nlist->getNListArray(),
                                       access_location::device, access_mode::read);
    ArrayHandle<size_t>  d_head_list(this->m_nlist->getHeadList(),
                                      access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_solid_index(this->m_solidgroup->getIndexArray(),
                                             access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_type_map(this->m_type_property_map,
                                          access_location::device, access_mode::read);

    m_tuner_noslip->begin();
    kernel::gpu_sph_noslip<KT_, SET1_>(
        solid_group_size, d_solid_index.data,
        d_pos.data, d_vel.data, d_density.data, d_pressure.data,
        d_vf.data, d_accel.data, d_h.data,
        d_n_neigh.data, d_nlist.data, d_head_list.data,
        d_type_map.data, box, kp, eos, bodyforce,
        m_tuner_noslip->getParam()[0]);
    if (this->m_exec_conf->isCUDAErrorCheckingEnabled()) CHECK_CUDA_ERROR();
    m_tuner_noslip->end();
    }
    }

// =========================================================================
// compute_colorgradients — zero Aux2/Aux3 on device when sigma=0
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::compute_colorgradients(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "TwoPhaseFlowGPU::compute_colorgradients" << endl;

    if (this->m_sigma12 == Scalar(0.0) &&
        this->m_sigma01  == Scalar(0.0) &&
        this->m_sigma02  == Scalar(0.0))
        {
        // No surface tension: zero Aux2 (surface normal) and Aux3 (curvature normal) on device.
        ArrayHandle<Scalar3> d_sn(this->m_pdata->getAuxiliaries2(),
                                   access_location::device, access_mode::overwrite);
        ArrayHandle<Scalar3> d_fn(this->m_pdata->getAuxiliaries3(),
                                   access_location::device, access_mode::overwrite);
        hipMemset(d_sn.data, 0, sizeof(Scalar3) * this->m_pdata->getAuxiliaries2().getNumElements());
        hipMemset(d_fn.data, 0, sizeof(Scalar3) * this->m_pdata->getAuxiliaries3().getNumElements());
        }
    else
        {
        TwoPhaseFlow<KT_, SET1_, SET2_>::compute_colorgradients(timestep);
        }
    }

// =========================================================================
// compute_surfaceforce — zero Aux4 on device when sigma=0
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::compute_surfaceforce(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "TwoPhaseFlowGPU::compute_surfaceforce" << endl;

    if (this->m_sigma12 == Scalar(0.0) &&
        this->m_sigma01  == Scalar(0.0) &&
        this->m_sigma02  == Scalar(0.0))
        {
        // No surface tension: zero Aux4 (surface force density) on device.
        ArrayHandle<Scalar3> d_sf(this->m_pdata->getAuxiliaries4(),
                                   access_location::device, access_mode::overwrite);
        hipMemset(d_sf.data, 0, sizeof(Scalar3) * this->m_pdata->getAuxiliaries4().getNumElements());
        }
    else
        {
        TwoPhaseFlow<KT_, SET1_, SET2_>::compute_surfaceforce(timestep);
        }
    }

// =========================================================================
// forcecomputation — two-phase pressure + viscosity + surface force
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::forcecomputation(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "TwoPhaseFlowGPU::forcecomputation" << endl;

    const BoxDim box        = this->m_pdata->getGlobalBox();
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    if (group_size == 0) return;

    SPHKernelDevParams  kp     = make_kparams();
    SPHEOSDevParams     eos1   = make_eos1params();
    SPHEOSDevParams     eos2   = make_eos2params();
    SPHNNViscParams     nn1    = make_nn1params();
    SPHNNViscParams     nn2    = make_nn2params();
    SPHTwoPhaseParams   fparams= make_2pfparams(timestep);

    { // device scope
    ArrayHandle<Scalar4> d_force  (this->m_force,   access_location::device, access_mode::overwrite);
    ArrayHandle<Scalar4> d_ratedpe(this->m_ratedpe, access_location::device, access_mode::overwrite);

    hipMemset(d_force.data,   0, sizeof(Scalar4) * this->m_force.getNumElements());
    hipMemset(d_ratedpe.data, 0, sizeof(Scalar4) * this->m_ratedpe.getNumElements());

    ArrayHandle<Scalar4> d_pos    (this->m_pdata->getPositions(),
                                    access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_vel    (this->m_pdata->getVelocities(),
                                    access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_density(this->m_pdata->getDensities(),
                                    access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_pressure(this->m_pdata->getPressures(),
                                     access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_vf     (this->m_pdata->getAuxiliaries1(),
                                    access_location::device, access_mode::read);
    // aux4 holds surface force density pre-computed by CPU compute_surfaceforce()
    ArrayHandle<Scalar3> d_sf     (this->m_pdata->getAuxiliaries4(),
                                    access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h      (this->m_pdata->getSlengths(),
                                    access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(),
                                         access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist  (this->m_nlist->getNListArray(),
                                         access_location::device, access_mode::read);
    ArrayHandle<size_t>  d_head_list   (this->m_nlist->getHeadList(),
                                         access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_index_array(this->m_fluidgroup->getIndexArray(),
                                             access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_type_map(this->m_type_property_map,
                                          access_location::device, access_mode::read);
    ArrayHandle<uint32_t> d_max_vel_bits(m_max_vel_bits,
                                          access_location::device, access_mode::overwrite);
    hipMemset(d_max_vel_bits.data, 0, sizeof(uint32_t));

    // Dispatch to the viscosity-model-templated fast kernel when both models
    // are NEWTONIAN or POWERLAW; fall back to the runtime-switch kernel for
    // CARREAU, BINGHAM, and HERSCHELBULKLEY.
    using NM = NonNewtonianModel;
    const NM nm1 = static_cast<NM>(nn1.model);
    const NM nm2 = static_cast<NM>(nn2.model);
    const bool use_fast = (nm1 == NM::NEWTONIAN || nm1 == NM::POWERLAW) &&
                          (nm2 == NM::NEWTONIAN || nm2 == NM::POWERLAW);

#define _FAST_CALL(M1, M2) \
    kernel::gpu_sph_2pf_forcecomputation_fast<KT_, SET1_, SET2_, M1, M2>( \
        group_size, d_index_array.data, \
        d_pos.data, d_vel.data, d_density.data, d_pressure.data, \
        d_vf.data, d_sf.data, d_h.data, \
        d_force.data, d_ratedpe.data, \
        d_n_neigh.data, d_nlist.data, d_head_list.data, \
        d_type_map.data, d_max_vel_bits.data, \
        box, kp, eos1, eos2, nn1, nn2, fparams, \
        m_tuner_force->getParam()[0])

    m_tuner_force->begin();
    if (use_fast)
        {
        if      (nm1 == NM::NEWTONIAN && nm2 == NM::NEWTONIAN) _FAST_CALL(NM::NEWTONIAN, NM::NEWTONIAN);
        else if (nm1 == NM::POWERLAW  && nm2 == NM::NEWTONIAN) _FAST_CALL(NM::POWERLAW,  NM::NEWTONIAN);
        else if (nm1 == NM::NEWTONIAN && nm2 == NM::POWERLAW ) _FAST_CALL(NM::NEWTONIAN, NM::POWERLAW );
        else                                                    _FAST_CALL(NM::POWERLAW,  NM::POWERLAW );
        }
    else
        {
        kernel::gpu_sph_2pf_forcecomputation<KT_, SET1_, SET2_>(
            group_size, d_index_array.data,
            d_pos.data, d_vel.data, d_density.data, d_pressure.data,
            d_vf.data, d_sf.data, d_h.data,
            d_force.data, d_ratedpe.data,
            d_n_neigh.data, d_nlist.data, d_head_list.data,
            d_type_map.data, d_max_vel_bits.data,
            box, kp, eos1, eos2, nn1, nn2, fparams,
            m_tuner_force->getParam()[0]);
        }
#undef _FAST_CALL

    if (this->m_exec_conf->isCUDAErrorCheckingEnabled()) CHECK_CUDA_ERROR();
    m_tuner_force->end();

    // Async readback: copy 4 bytes to pinned host buffer without stalling the GPU.
    // We read the value written by the PREVIOUS step (1-step lag is fine for CFL).
    hipMemcpyAsync(m_max_vel_bits_host, d_max_vel_bits.data,
                   sizeof(uint32_t), hipMemcpyDeviceToHost, 0);
    } // end device scope

    // Use value from previous step's async copy (already complete by now).
    {
    float vel_float;
    memcpy(&vel_float, m_max_vel_bits_host, sizeof(float));
    this->m_timestep_list[5] = static_cast<double>(vel_float);
    }

    this->applyBodyForce(timestep, this->m_fluidgroup);
    }

// =========================================================================
// compute_solid_forces — two-phase solid reaction forces
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::compute_solid_forces(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "TwoPhaseFlowGPU::compute_solid_forces" << endl;

    unsigned int solid_group_size = this->m_solidgroup->getNumMembers();
    if (solid_group_size == 0) return;

    const BoxDim box = this->m_pdata->getGlobalBox();
    SPHKernelDevParams kp  = make_kparams();
    SPHNNViscParams    nn1 = make_nn1params();
    SPHNNViscParams    nn2 = make_nn2params();

    { // device scope
    ArrayHandle<Scalar4> d_force  (this->m_force,
                                    access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_pos    (this->m_pdata->getPositions(),
                                    access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_vel    (this->m_pdata->getVelocities(),
                                    access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_density(this->m_pdata->getDensities(),
                                    access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_pressure(this->m_pdata->getPressures(),
                                     access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h      (this->m_pdata->getSlengths(),
                                    access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(),
                                         access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist  (this->m_nlist->getNListArray(),
                                         access_location::device, access_mode::read);
    ArrayHandle<size_t>  d_head_list   (this->m_nlist->getHeadList(),
                                         access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_solid_index(this->m_solidgroup->getIndexArray(),
                                             access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_type_map(this->m_type_property_map,
                                          access_location::device, access_mode::read);

    m_tuner_solidforce->begin();
    kernel::gpu_sph_2pf_solid_forces<KT_, SET1_, SET2_>(
        solid_group_size, d_solid_index.data,
        d_pos.data, d_vel.data, d_density.data, d_pressure.data, d_h.data,
        d_force.data,
        d_n_neigh.data, d_nlist.data, d_head_list.data,
        d_type_map.data,
        box, kp, nn1, nn2,
        static_cast<int>(this->m_density_method),
        m_tuner_solidforce->getParam()[0]);
    if (this->m_exec_conf->isCUDAErrorCheckingEnabled()) CHECK_CUDA_ERROR();
    m_tuner_solidforce->end();
    }
    }

// =========================================================================
// pybind11 export
// =========================================================================

namespace detail
{
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void export_TwoPhaseFlowGPU(pybind11::module& m, std::string name)
    {
    pybind11::class_<TwoPhaseFlowGPU<KT_, SET1_, SET2_>,
                     TwoPhaseFlow<KT_, SET1_, SET2_>,
                     std::shared_ptr<TwoPhaseFlowGPU<KT_, SET1_, SET2_>>>(m, name.c_str())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<SmoothingKernel<KT_>>,
                            std::shared_ptr<StateEquation<SET1_>>,
                            std::shared_ptr<StateEquation<SET2_>>,
                            std::shared_ptr<nsearch::NeighborList>,
                            std::shared_ptr<ParticleGroup>,
                            std::shared_ptr<ParticleGroup>,
                            std::shared_ptr<ParticleGroup>,
                            DensityMethod,
                            ViscosityMethod,
                            ColorGradientMethod>());
    }
} // namespace detail

// =========================================================================
// Explicit template instantiations (5 kernels × 2 × 2 = 20)
// =========================================================================

#define INST_2PF_GPU(KT, SET1, SET2) \
    template class PYBIND11_EXPORT TwoPhaseFlowGPU<KT, SET1, SET2>; \
    namespace detail { \
    template void export_TwoPhaseFlowGPU<KT, SET1, SET2>(pybind11::module& m, std::string name); \
    }

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

} // namespace sph
} // namespace hoomd

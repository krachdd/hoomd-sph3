/* ---------------------------------------------------------
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.
Redistribution and use permitted under BSD 3-Clause License.
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*! \file SinglePhaseFlowGDGDGPU.cc
    \brief Host-side implementation of SinglePhaseFlowGDGDGPU.
*/

#include "SinglePhaseFlowGDGDGPU.h"
#include "SinglePhaseFlowGPU.cuh"
#include "SinglePhaseFlowGDGDGPU.cuh"

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

template<SmoothingKernelType KT_, StateEquationType SET_>
SinglePhaseFlowGDGDGPU<KT_, SET_>::SinglePhaseFlowGDGDGPU(
    std::shared_ptr<SystemDefinition>      sysdef,
    std::shared_ptr<SmoothingKernel<KT_>>  skernel,
    std::shared_ptr<StateEquation<SET_>>   equationofstate,
    std::shared_ptr<nsearch::NeighborList> nlist,
    std::shared_ptr<ParticleGroup>         fluidgroup,
    std::shared_ptr<ParticleGroup>         solidgroup,
    DensityMethod   mdensitymethod,
    ViscosityMethod mviscositymethod)
    : SinglePhaseFlowGDGD<KT_, SET_>(sysdef, skernel, equationofstate, nlist,
                                      fluidgroup, solidgroup,
                                      mdensitymethod, mviscositymethod),
      m_max_vel_bits(1, this->m_exec_conf)
    {
    if (!this->m_exec_conf->isCUDAEnabled())
        {
        this->m_exec_conf->msg->error()
            << "Creating a SinglePhaseFlowGDGDGPU without a GPU in the execution configuration"
            << endl;
        throw std::runtime_error("Error initializing SinglePhaseFlowGDGDGPU");
        }

    m_tuner_ndensity  = std::make_shared<Autotuner<1>>(std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(this->m_exec_conf)}, this->m_exec_conf, "sph_gdgd_ndensity");
    m_tuner_pressure  = std::make_shared<Autotuner<1>>(std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(this->m_exec_conf)}, this->m_exec_conf, "sph_gdgd_pressure");
    m_tuner_noslip    = std::make_shared<Autotuner<1>>(std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(this->m_exec_conf)}, this->m_exec_conf, "sph_gdgd_noslip");
    m_tuner_force     = std::make_shared<Autotuner<1>>(std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(this->m_exec_conf)}, this->m_exec_conf, "sph_gdgd_force");
    m_tuner_solidforce= std::make_shared<Autotuner<1>>(std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(this->m_exec_conf)}, this->m_exec_conf, "sph_gdgd_solidforce");

    this->m_autotuners.push_back(m_tuner_ndensity);
    this->m_autotuners.push_back(m_tuner_pressure);
    this->m_autotuners.push_back(m_tuner_noslip);
    this->m_autotuners.push_back(m_tuner_force);
    this->m_autotuners.push_back(m_tuner_solidforce);
    }

template<SmoothingKernelType KT_, StateEquationType SET_>
SinglePhaseFlowGDGDGPU<KT_, SET_>::~SinglePhaseFlowGDGDGPU()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying SinglePhaseFlowGDGDGPU" << endl;
    }

// =========================================================================
// Helper: build parameter structs
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
SPHKernelDevParams SinglePhaseFlowGDGDGPU<KT_, SET_>::make_kparams() const
    {
    SPHKernelDevParams kp;
    kp.alpha        = this->m_skernel->getAlpha();
    kp.self_density = this->m_skernel->getSelfDensity();
    kp.ch           = this->m_ch;
    kp.rcutsq       = this->m_rcutsq;
    kp.const_slength= this->m_const_slength ? 1 : 0;
    kp.kappa        = this->m_kappa;
    return kp;
    }

template<SmoothingKernelType KT_, StateEquationType SET_>
SPHEOSDevParams SinglePhaseFlowGDGDGPU<KT_, SET_>::make_eosparams() const
    {
    SPHEOSDevParams eos;
    eos.rho0 = this->m_eos->getRestDensity();
    eos.c    = this->m_eos->getSpeedOfSound();
    eos.bp   = this->m_eos->getBackgroundPressure();
    return eos;
    }

template<SmoothingKernelType KT_, StateEquationType SET_>
SPHNNViscParams SinglePhaseFlowGDGDGPU<KT_, SET_>::make_nnparams() const
    {
    SPHNNViscParams nn;
    nn.mu       = this->m_mu;
    nn.K        = this->m_nn_K;
    nn.n        = this->m_nn_n;
    nn.mu0      = this->m_nn_mu0;
    nn.muinf    = this->m_nn_muinf;
    nn.lambda_NN= this->m_nn_lambda;
    nn.tauy     = this->m_nn_tauy;
    nn.m_reg    = this->m_nn_m;
    nn.mu_min   = this->m_nn_mu_min;
    nn.model    = static_cast<int>(this->m_nn_model);
    return nn;
    }

// =========================================================================
// compute_ndensity
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowGDGDGPU<KT_, SET_>::compute_ndensity(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "SinglePhaseFlowGDGDGPU::compute_ndensity" << endl;

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
    kernel::gpu_sph_ndensity<KT_, SET_>(
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
// compute_pressure
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowGDGDGPU<KT_, SET_>::compute_pressure(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "SinglePhaseFlowGDGDGPU::compute_pressure" << endl;

    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    if (group_size == 0) return;

    SPHEOSDevParams eos = make_eosparams();

    { // device scope
    ArrayHandle<Scalar>  d_density (this->m_pdata->getDensities(),
                                     access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_pressure(this->m_pdata->getPressures(),
                                     access_location::device, access_mode::readwrite);
    ArrayHandle<unsigned int> d_index_array(this->m_fluidgroup->getIndexArray(),
                                             access_location::device, access_mode::read);

    m_tuner_pressure->begin();
    kernel::gpu_sph_pressure<KT_, SET_>(
        group_size, d_index_array.data,
        d_density.data, d_pressure.data, eos,
        m_tuner_pressure->getParam()[0]);
    if (this->m_exec_conf->isCUDAErrorCheckingEnabled()) CHECK_CUDA_ERROR();
    m_tuner_pressure->end();
    }
    }

// =========================================================================
// compute_noslip
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowGDGDGPU<KT_, SET_>::compute_noslip(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "SinglePhaseFlowGDGDGPU::compute_noslip" << endl;

    unsigned int solid_group_size = this->m_solidgroup->getNumMembers();
    if (solid_group_size == 0) return;

    const BoxDim box = this->m_pdata->getGlobalBox();
    SPHKernelDevParams kp  = make_kparams();
    SPHEOSDevParams    eos = make_eosparams();

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
    kernel::gpu_sph_noslip<KT_, SET_>(
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
// forcecomputation (GDGD-specific GPU kernel)
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowGDGDGPU<KT_, SET_>::forcecomputation(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "SinglePhaseFlowGDGDGPU::forcecomputation" << endl;

    const BoxDim box        = this->m_pdata->getGlobalBox();
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    if (group_size == 0) return;

    SPHKernelDevParams kp  = make_kparams();
    SPHEOSDevParams    eos = make_eosparams();
    SPHNNViscParams    nn  = make_nnparams();

    { // device scope
    ArrayHandle<Scalar4> d_force  (this->m_force,    access_location::device, access_mode::overwrite);
    ArrayHandle<Scalar4> d_ratedpe(this->m_ratedpe,  access_location::device, access_mode::overwrite);

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
    ArrayHandle<Scalar3> d_aux4   (this->m_pdata->getAuxiliaries4(),
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

    m_tuner_force->begin();
    kernel::gpu_sph_gdgd_forcecomputation<KT_, SET_>(
        group_size, d_index_array.data,
        d_pos.data, d_vel.data, d_density.data, d_pressure.data,
        d_vf.data, d_aux4.data, d_h.data,
        d_force.data, d_ratedpe.data,
        d_n_neigh.data, d_nlist.data, d_head_list.data,
        d_type_map.data, d_max_vel_bits.data,
        box, kp, eos, nn,
        this->m_kappa_s, this->m_beta_s, this->m_scalar_ref,
        this->m_rho0,
        this->m_boussinesq ? 1 : 0,
        this->m_bodyforce,
        static_cast<int>(this->m_density_method),
        this->m_artificial_viscosity ? 1 : 0,
        this->m_avalpha, this->m_avbeta,
        this->m_density_diffusion  ? 1 : 0,
        this->m_ddiff,
        m_tuner_force->getParam()[0]);
    if (this->m_exec_conf->isCUDAErrorCheckingEnabled()) CHECK_CUDA_ERROR();
    m_tuner_force->end();
    } // end device scope

    // Read back max velocity for adaptive timestep
    {
    ArrayHandle<uint32_t> h_max_vel_bits(m_max_vel_bits,
                                          access_location::host, access_mode::read);
    uint32_t bits = h_max_vel_bits.data[0];
    float vel_float;
    memcpy(&vel_float, &bits, sizeof(float));
    this->m_timestep_list[5] = static_cast<double>(vel_float);
    }

    this->applyBodyForce(timestep, this->m_fluidgroup);
    }

// =========================================================================
// compute_solid_forces
// =========================================================================

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowGDGDGPU<KT_, SET_>::compute_solid_forces(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "SinglePhaseFlowGDGDGPU::compute_solid_forces" << endl;

    unsigned int solid_group_size = this->m_solidgroup->getNumMembers();
    if (solid_group_size == 0) return;

    const BoxDim box = this->m_pdata->getGlobalBox();
    SPHKernelDevParams kp  = make_kparams();
    SPHNNViscParams    nn  = make_nnparams();

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
    kernel::gpu_sph_solid_forces<KT_, SET_>(
        solid_group_size, d_solid_index.data,
        d_pos.data, d_vel.data, d_density.data, d_pressure.data, d_h.data,
        d_force.data,
        d_n_neigh.data, d_nlist.data, d_head_list.data,
        d_type_map.data,
        box, kp, nn,
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
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlowGDGDGPU(pybind11::module& m, std::string name)
    {
    pybind11::class_<SinglePhaseFlowGDGDGPU<KT_, SET_>,
                     SinglePhaseFlowGDGD<KT_, SET_>,
                     std::shared_ptr<SinglePhaseFlowGDGDGPU<KT_, SET_>>>(m, name.c_str())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<SmoothingKernel<KT_>>,
                            std::shared_ptr<StateEquation<SET_>>,
                            std::shared_ptr<nsearch::NeighborList>,
                            std::shared_ptr<ParticleGroup>,
                            std::shared_ptr<ParticleGroup>,
                            DensityMethod,
                            ViscosityMethod>());
    }
} // namespace detail

// =========================================================================
// Explicit template instantiations
// =========================================================================

#define INST_GDGD_GPU(KT, SET) \
    template class PYBIND11_EXPORT SinglePhaseFlowGDGDGPU<KT, SET>; \
    namespace detail { \
    template void export_SinglePhaseFlowGDGDGPU<KT, SET>(pybind11::module& m, std::string name); \
    }

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

} // namespace sph
} // namespace hoomd

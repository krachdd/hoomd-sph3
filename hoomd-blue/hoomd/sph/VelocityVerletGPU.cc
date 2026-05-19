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

/*! \file VelocityVerletGPU.cc
    \brief Host-side implementation of VelocityVerletGPU.
*/

#include "VelocityVerletGPU.h"
#include "VelocityVerletGPU.cuh"

#include <stdexcept>
#include <pybind11/pybind11.h>

using namespace std;

namespace hoomd
{
namespace sph
{

VelocityVerletGPU::VelocityVerletGPU(std::shared_ptr<SystemDefinition> sysdef,
                                     std::shared_ptr<ParticleGroup>    group)
    : VelocityVerlet(sysdef, group)
    {
    m_exec_conf->msg->notice(5) << "Constructing VelocityVerletGPU" << endl;

    if (!m_exec_conf->isCUDAEnabled())
        {
        m_exec_conf->msg->error()
            << "Creating a VelocityVerletGPU without a GPU in the execution configuration"
            << endl;
        throw std::runtime_error("Error initializing VelocityVerletGPU");
        }

    m_tuner_step1 = std::make_shared<Autotuner<1>>(
        std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(m_exec_conf)},
        m_exec_conf, "vv_step1");
    m_tuner_step2 = std::make_shared<Autotuner<1>>(
        std::vector<std::vector<unsigned int>>{AutotunerBase::makeBlockSizeRange(m_exec_conf)},
        m_exec_conf, "vv_step2");

    m_autotuners.push_back(m_tuner_step1);
    m_autotuners.push_back(m_tuner_step2);
    }

VelocityVerletGPU::~VelocityVerletGPU()
    {
    m_exec_conf->msg->notice(5) << "Destroying VelocityVerletGPU" << endl;
    }

void VelocityVerletGPU::integrateStepOne(uint64_t timestep)
    {
    m_exec_conf->msg->notice(9) << "VelocityVerletGPU: integrateStepOne" << endl;

    unsigned int group_size = m_group->getNumMembers();
    if (group_size == 0)
        return;

    const BoxDim& box = m_pdata->getBox();

    ArrayHandle<Scalar4>      d_pos    (m_pdata->getPositions(),    access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4>      d_vel    (m_pdata->getVelocities(),   access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3>      d_accel  (m_pdata->getAccelerations(),access_location::device, access_mode::read);
    ArrayHandle<Scalar>       d_density(m_pdata->getDensities(),    access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar>       d_pressure(m_pdata->getPressures(),   access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3>      d_dpedt  (m_pdata->getDPEdts(),       access_location::device, access_mode::read);
    ArrayHandle<int3>         d_image  (m_pdata->getImages(),       access_location::device, access_mode::readwrite);
    ArrayHandle<unsigned int> d_index  (m_group->getIndexArray(),   access_location::device, access_mode::read);

    m_tuner_step1->begin();
    kernel::gpu_vv_step1(
        group_size, d_index.data,
        d_pos.data, d_vel.data, d_accel.data,
        d_density.data, d_pressure.data, d_dpedt.data,
        d_image.data, box, m_deltaT,
        m_tuner_step1->getParam()[0]);
    if (m_exec_conf->isCUDAErrorCheckingEnabled()) CHECK_CUDA_ERROR();
    m_tuner_step1->end();
    }

void VelocityVerletGPU::integrateStepTwo(uint64_t timestep)
    {
    m_exec_conf->msg->notice(9) << "VelocityVerletGPU: integrateStepTwo" << endl;

    unsigned int group_size = m_group->getNumMembers();
    if (group_size == 0)
        return;

    ArrayHandle<Scalar4>      d_pos      (m_pdata->getPositions(),        access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4>      d_vel      (m_pdata->getVelocities(),        access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3>      d_accel    (m_pdata->getAccelerations(),     access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar>       d_density  (m_pdata->getDensities(),         access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar>       d_pressure (m_pdata->getPressures(),         access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3>      d_dpedt    (m_pdata->getDPEdts(),            access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4>      d_net_force(m_pdata->getNetForce(),          access_location::device, access_mode::read);
    ArrayHandle<Scalar4>      d_net_ratedpe(m_pdata->getNetRateDPEArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_index    (m_group->getIndexArray(),        access_location::device, access_mode::read);

    m_tuner_step2->begin();
    kernel::gpu_vv_step2(
        group_size, d_index.data,
        d_pos.data, d_vel.data, d_accel.data,
        d_density.data, d_pressure.data, d_dpedt.data,
        d_net_force.data, d_net_ratedpe.data,
        m_deltaT, m_tuner_step2->getParam()[0]);
    if (m_exec_conf->isCUDAErrorCheckingEnabled()) CHECK_CUDA_ERROR();
    m_tuner_step2->end();
    }

namespace detail
{
void export_VelocityVerletGPU(pybind11::module& m)
    {
    pybind11::class_<VelocityVerletGPU, VelocityVerlet,
                     std::shared_ptr<VelocityVerletGPU>>(m, "VelocityVerletGPU")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<ParticleGroup>>())
        .def("getDensityMethod", &VelocityVerletGPU::getDensityMethod)
        .def("setDensityMethod", &VelocityVerletGPU::setDensityMethod)
        .def_property("limit",      &VelocityVerletGPU::getLimit,     &VelocityVerletGPU::setLimit)
        .def_property("zero_force", &VelocityVerletGPU::getZeroForce, &VelocityVerletGPU::setZeroForce);
    }
} // namespace detail

} // namespace sph
} // namespace hoomd

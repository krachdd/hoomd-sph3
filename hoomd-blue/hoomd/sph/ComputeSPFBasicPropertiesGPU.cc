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
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*! \file ComputeSPFBasicPropertiesGPU.cc
    \brief Host-side implementation of ComputeSPFBasicPropertiesGPU.

    computeProperties() dispatches the two-pass GPU reduction and then
    copies only the 48-byte SPFProps result to the CPU-side m_properties
    GPUArray so that all inherited getter methods work without change.
*/

#include "ComputeSPFBasicPropertiesGPU.h"

#include <iostream>
#include <stdexcept>
#include <pybind11/pybind11.h>

using namespace std;

namespace hoomd
{
namespace sph
{

// =========================================================================
// Constructor / Destructor
// =========================================================================

ComputeSPFBasicPropertiesGPU::ComputeSPFBasicPropertiesGPU(
    std::shared_ptr<SystemDefinition> sysdef,
    std::shared_ptr<ParticleGroup>    group)
    : ComputeSPFBasicProperties(sysdef, group),
      m_d_props(1, m_exec_conf),
      m_d_block_props(0, m_exec_conf)
    {
    if (!m_exec_conf->isCUDAEnabled())
        {
        m_exec_conf->msg->error()
            << "ComputeSPFBasicPropertiesGPU created without a GPU device" << endl;
        throw std::runtime_error("Error initialising ComputeSPFBasicPropertiesGPU");
        }
    m_exec_conf->msg->notice(5) << "Constructing ComputeSPFBasicPropertiesGPU" << endl;
    }

ComputeSPFBasicPropertiesGPU::~ComputeSPFBasicPropertiesGPU()
    {
    m_exec_conf->msg->notice(5) << "Destroying ComputeSPFBasicPropertiesGPU" << endl;
    }

// =========================================================================
// computeProperties — GPU reduction
// =========================================================================

void ComputeSPFBasicPropertiesGPU::computeProperties()
    {
    if (m_group->getNumMembersGlobal() == 0) return;

    const unsigned int group_size = m_group->getNumMembers();

    // Ensure temporary block buffer is large enough for n_blocks = ceil(group_size/256).
    const unsigned int n_blocks = (group_size + 255) / 256;
    if (m_d_block_props.getNumElements() < n_blocks)
        m_d_block_props.resize(n_blocks);

    { // device scope — all pointers stay on device
    ArrayHandle<unsigned int>      d_index_array(m_group->getIndexArray(),
                                                  access_location::device, access_mode::read);
    ArrayHandle<Scalar4>           d_vel        (m_pdata->getVelocities(),
                                                  access_location::device, access_mode::read);
    ArrayHandle<Scalar>            d_density    (m_pdata->getDensities(),
                                                  access_location::device, access_mode::read);
    ArrayHandle<kernel::SPFProps>  d_block      (m_d_block_props,
                                                  access_location::device, access_mode::overwrite);
    ArrayHandle<kernel::SPFProps>  d_out        (m_d_props,
                                                  access_location::device, access_mode::overwrite);

    kernel::gpu_spf_basic_props_reduce(
        group_size,
        n_blocks,
        d_index_array.data,
        d_vel.data,
        d_density.data,
        d_block.data,
        d_out.data);
    if (m_exec_conf->isCUDAErrorCheckingEnabled()) CHECK_CUDA_ERROR();
    } // close device scope — ArrayHandles released

    // Transfer the 48-byte result to the host-side m_properties GPUArray.
    // This is the *only* device->host copy: 6 doubles, not ~64 MB of particle data.
    {
    ArrayHandle<kernel::SPFProps> h_d_props(m_d_props,
                                             access_location::host, access_mode::read);
    ArrayHandle<Scalar> h_props(m_properties,
                                 access_location::host, access_mode::overwrite);

    const kernel::SPFProps& r = h_d_props.data[0];
    h_props.data[singlephaseflow_logger_index::sum_fluid_velocity_x] = Scalar(r.vx);
    h_props.data[singlephaseflow_logger_index::sum_fluid_velocity_y] = Scalar(r.vy);
    h_props.data[singlephaseflow_logger_index::sum_fluid_velocity_z] = Scalar(r.vz);
    h_props.data[singlephaseflow_logger_index::sum_fluid_density]    = Scalar(r.density);
    h_props.data[singlephaseflow_logger_index::abs_velocity]         = Scalar(r.abs_v);
    h_props.data[singlephaseflow_logger_index::e_kin_fluid]          = Scalar(r.e_kin);
    }

#ifdef ENABLE_MPI
    m_properties_reduced = !m_pdata->getDomainDecomposition();
#endif
    }

// =========================================================================
// Python export
// =========================================================================

namespace detail
{
void export_ComputeSPFMechanicalPropertiesGPU(pybind11::module& m)
    {
    pybind11::class_<ComputeSPFBasicPropertiesGPU,
                     ComputeSPFBasicProperties,
                     std::shared_ptr<ComputeSPFBasicPropertiesGPU>>(
        m, "ComputeSPFBasicPropertiesGPU")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<ParticleGroup>>())
        .def_property_readonly("num_particles",   &ComputeSPFBasicPropertiesGPU::getNumParticles)
        .def_property_readonly("abs_velocity",    &ComputeSPFBasicPropertiesGPU::getAbsoluteVelocity)
        .def_property_readonly("e_kin_fluid",     &ComputeSPFBasicPropertiesGPU::getEkinFluid)
        .def_property_readonly("fluid_vel_x_sum", &ComputeSPFBasicPropertiesGPU::getSumFluidXVelocity)
        .def_property_readonly("fluid_vel_y_sum", &ComputeSPFBasicPropertiesGPU::getSumFluidYVelocity)
        .def_property_readonly("fluid_vel_z_sum", &ComputeSPFBasicPropertiesGPU::getSumFluidZVelocity)
        .def_property_readonly("mean_density",    &ComputeSPFBasicPropertiesGPU::getMeanFluidDensity)
        .def_property_readonly("volume",          &ComputeSPFBasicPropertiesGPU::getVolume);
    }
} // namespace detail

} // namespace sph
} // namespace hoomd

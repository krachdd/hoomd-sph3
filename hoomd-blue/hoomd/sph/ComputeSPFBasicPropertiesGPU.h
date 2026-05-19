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

/*! \file ComputeSPFBasicPropertiesGPU.h
    \brief GPU-accelerated version of ComputeSPFBasicProperties.

    Overrides computeProperties() to run a fully on-device two-pass
    tree-reduction via gpu_spf_basic_props_reduce().  The only host<->device
    transfer is 48 bytes (one SPFProps struct) at the moment a Python property
    is accessed, rather than copying the full particle-velocity and density
    arrays (~64 MB) on every log step.
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __COMPUTE_SPF_BASIC_PROPERTIES_GPU_H__
#define __COMPUTE_SPF_BASIC_PROPERTIES_GPU_H__

#include "ComputeSPFBasicProperties.h"
#include "ComputeSPFBasicPropertiesGPU.cuh"
#include "hoomd/GPUArray.h"

#include <pybind11/pybind11.h>

namespace hoomd
{
namespace sph
{

/*! GPU-side reduction for SinglePhaseFlow basic mechanical properties.

    Inherits all getter methods from ComputeSPFBasicProperties; only
    computeProperties() is replaced with a GPU kernel dispatch.
*/
class PYBIND11_EXPORT ComputeSPFBasicPropertiesGPU : public ComputeSPFBasicProperties
    {
    public:
    ComputeSPFBasicPropertiesGPU(std::shared_ptr<SystemDefinition> sysdef,
                                  std::shared_ptr<ParticleGroup>    group);

    virtual ~ComputeSPFBasicPropertiesGPU();

    protected:
    //! GPU-side reduction; writes into m_d_props on the device.
    virtual void computeProperties() override;

    private:
    /// Device-side output buffer: one SPFProps (6 doubles = 48 bytes).
    GPUArray<kernel::SPFProps> m_d_props;

    /// Temporary per-block partial sums (resized lazily to fit n_blocks).
    GPUArray<kernel::SPFProps> m_d_block_props;
    };

namespace detail
{
void export_ComputeSPFMechanicalPropertiesGPU(pybind11::module& m);
} // namespace detail

} // namespace sph
} // namespace hoomd

#endif // __COMPUTE_SPF_BASIC_PROPERTIES_GPU_H__

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

/*! \file VelocityVerletGPU.h
    \brief GPU-accelerated velocity-Verlet SPH integrator.

    Inherits VelocityVerlet and overrides integrateStepOne / integrateStepTwo
    with GPU kernel launches that keep all particle data on the device.
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __VelocityVerletGPU_H__
#define __VelocityVerletGPU_H__

#include "VelocityVerlet.h"
#include "VelocityVerletGPU.cuh"
#include "hoomd/Autotuner.h"

#include <pybind11/pybind11.h>

namespace hoomd
{
namespace sph
{

class PYBIND11_EXPORT VelocityVerletGPU : public VelocityVerlet
    {
    public:
        VelocityVerletGPU(std::shared_ptr<SystemDefinition> sysdef,
                          std::shared_ptr<ParticleGroup>    group);
        virtual ~VelocityVerletGPU();

        virtual void integrateStepOne(uint64_t timestep) override;
        virtual void integrateStepTwo(uint64_t timestep) override;

    private:
        std::shared_ptr<Autotuner<1>> m_tuner_step1;
        std::shared_ptr<Autotuner<1>> m_tuner_step2;
    };

namespace detail
{
void export_VelocityVerletGPU(pybind11::module& m);
} // namespace detail

} // namespace sph
} // namespace hoomd

#endif // __VelocityVerletGPU_H__

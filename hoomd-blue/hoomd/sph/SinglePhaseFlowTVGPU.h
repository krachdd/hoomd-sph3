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

/*! \file SinglePhaseFlowTVGPU.h
    \brief GPU-accelerated transport-velocity SPH solver.

    Inherits SinglePhaseFlowTV<KT_, SET_> and overrides the virtual compute
    methods.  ndensity / pressure / noslip / solid-forces use the shared GPU
    kernels from SinglePhaseFlowGPU.cuh; forcecomputation uses the TV-specific
    kernel from SinglePhaseFlowTVGPU.cuh that adds the artificial stress tensor
    and background pressure contribution.
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __SinglePhaseFlowTVGPU_H__
#define __SinglePhaseFlowTVGPU_H__

#include "SinglePhaseFlowTV.h"
#include "SinglePhaseFlowGPU.cuh"
#include "SinglePhaseFlowTVGPU.cuh"
#include "hoomd/Autotuner.h"
#include "hoomd/GPUArray.h"

#include <pybind11/pybind11.h>

namespace hoomd
{
namespace sph
{

template<SmoothingKernelType KT_, StateEquationType SET_>
class PYBIND11_EXPORT SinglePhaseFlowTVGPU : public SinglePhaseFlowTV<KT_, SET_>
    {
    public:
        SinglePhaseFlowTVGPU(std::shared_ptr<SystemDefinition>      sysdef,
                             std::shared_ptr<SmoothingKernel<KT_>>  skernel,
                             std::shared_ptr<StateEquation<SET_>>   equationofstate,
                             std::shared_ptr<nsearch::NeighborList> nlist,
                             std::shared_ptr<ParticleGroup>         fluidgroup,
                             std::shared_ptr<ParticleGroup>         solidgroup,
                             DensityMethod   mdensitymethod  = DENSITYSUMMATION,
                             ViscosityMethod mviscositymethod = HARMONICAVERAGE);

        virtual ~SinglePhaseFlowTVGPU();

    protected:
        std::shared_ptr<Autotuner<1>> m_tuner_ndensity;
        std::shared_ptr<Autotuner<1>> m_tuner_pressure;
        std::shared_ptr<Autotuner<1>> m_tuner_noslip;
        std::shared_ptr<Autotuner<1>> m_tuner_force;
        std::shared_ptr<Autotuner<1>> m_tuner_solidforce;

        GPUArray<uint32_t> m_max_vel_bits;

        virtual void compute_ndensity(uint64_t timestep) override;
        virtual void compute_pressure(uint64_t timestep) override;
        virtual void compute_noslip(uint64_t timestep)   override;
        virtual void forcecomputation(uint64_t timestep) override;
        virtual void compute_solid_forces(uint64_t timestep) override;

    private:
        SPHKernelDevParams make_kparams()  const;
        SPHEOSDevParams    make_eosparams() const;
        SPHNNViscParams    make_nnparams()  const;
    };

namespace detail
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlowTVGPU(pybind11::module& m, std::string name);
} // namespace detail

} // namespace sph
} // namespace hoomd

#endif // __SinglePhaseFlowTVGPU_H__

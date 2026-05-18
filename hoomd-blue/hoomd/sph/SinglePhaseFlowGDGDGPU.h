/* ---------------------------------------------------------
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.
Redistribution and use permitted under BSD 3-Clause License.
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*! \file SinglePhaseFlowGDGDGPU.h
    \brief GPU-accelerated GDGD SPH solver.

    Inherits SinglePhaseFlowGDGD<KT_, SET_> and overrides the five virtual
    compute methods.  ndensity / pressure / noslip / solid-forces reuse the
    shared GPU kernels from SinglePhaseFlowGPU.cuh; forcecomputation uses the
    GDGD-specific kernel from SinglePhaseFlowGDGDGPU.cuh.
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __SinglePhaseFlowGDGDGPU_H__
#define __SinglePhaseFlowGDGDGPU_H__

#include "SinglePhaseFlowGDGD.h"
#include "SinglePhaseFlowGPU.cuh"
#include "SinglePhaseFlowGDGDGPU.cuh"
#include "hoomd/Autotuner.h"
#include "hoomd/GPUArray.h"

#include <pybind11/pybind11.h>

namespace hoomd
{
namespace sph
{

template<SmoothingKernelType KT_, StateEquationType SET_>
class PYBIND11_EXPORT SinglePhaseFlowGDGDGPU : public SinglePhaseFlowGDGD<KT_, SET_>
    {
    public:
        SinglePhaseFlowGDGDGPU(std::shared_ptr<SystemDefinition>      sysdef,
                               std::shared_ptr<SmoothingKernel<KT_>>  skernel,
                               std::shared_ptr<StateEquation<SET_>>   equationofstate,
                               std::shared_ptr<nsearch::NeighborList> nlist,
                               std::shared_ptr<ParticleGroup>         fluidgroup,
                               std::shared_ptr<ParticleGroup>         solidgroup,
                               DensityMethod   mdensitymethod  = DENSITYSUMMATION,
                               ViscosityMethod mviscositymethod = HARMONICAVERAGE);

        virtual ~SinglePhaseFlowGDGDGPU();

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
void export_SinglePhaseFlowGDGDGPU(pybind11::module& m, std::string name);
} // namespace detail

} // namespace sph
} // namespace hoomd

#endif // __SinglePhaseFlowGDGDGPU_H__

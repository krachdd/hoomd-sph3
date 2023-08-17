// SuspensionFlowFactory.h
#pragma once

#include "SmoothingKernel.h"
#include "StateEquations.h"

#include "SuspensionFlow.h"  // Include the SuspensionFlow class header


namespace hoomd 
{
namespace sph
{

class SuspensionFlowFactory {
public:
    template <SmoothingKernelType KT_, StateEquationType SET_>
    static std::shared_ptr<SuspensionFlow<KT_, SET_>> CreateSuspensionFlow(SmoothingKernelType kt,
                                           StateEquationType set) {
        // Create and return an instance of SuspensionFlow with the provided constructor arguments
        return std::make_shared<SuspensionFlow<KT_, SET_>>(
            /* constructor arguments */);
    }
};
}
}
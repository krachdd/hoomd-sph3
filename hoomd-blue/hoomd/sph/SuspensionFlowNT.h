// SuspensionFlowNT.h
#pragma once


//#include "SuspensionFlow.h"  // Include the SuspensionFlow class header


namespace hoomd 
{
namespace sph
{

class SuspensionFlowNT {
    public:

    /// Get the number of aggregates (global)
    unsigned int getNAggregatesGlobal() const

    virtual void setDeltaT(Scalar dt);

    // /// Validate rigid body constituent particles. The method purposely does not check
    // /// positions or orientation.
    virtual void validateRigidBodies();

    virtual void updateCompositeParticles(uint64_t timestep);

    protected:

    //     Scalar m_deltaT; //!< timestep size (required for some types of non-conservative forces)
    // virtual void someCommonOperation() = 0;
};

} // end namespace hoomd
} // end namespace sph
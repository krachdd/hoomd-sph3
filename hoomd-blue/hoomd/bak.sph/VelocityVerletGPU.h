// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "VelocityVerlet.h"

//#ifndef __VELOCITY_VERLET_H__
//#define __VELOCITY_VERLET_H__

/*! \file VelocityVerletGPU.h
 *    \brief Declares an VelocityVerletGPU integrator for performing two-step integration on SPH particles
 */

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include "hoomd/Autotuner.h"

//! Integrates part of the system forward in two steps in the NVE ensemble on the GPU
/*! Implements velocity-verlet NVE integration through the IntegrationMethodTwoStep interface, runs on the GPU
 * 
 *    \ingroup updaters
 */
class VelocityVerletGPU : public VelocityVerlet
{
public:
    //! Constructs the integration method and associates it with the system
    VelocityVerletGPU(boost::shared_ptr<SystemDefinition> sysdef, boost::shared_ptr<ParticleGroup> group, bool skip_restart);
    virtual ~VelocityVerletGPU() {};
    
    //! Performs the first step of the integration
    virtual void integrateStepOne(unsigned int timestep);
    
    //! Performs the second step of the integration
    virtual void integrateStepTwo(unsigned int timestep);
    
    //! Set autotuner parameters
    /*! \param enable Enable/disable autotuning
     *            \param period period (approximate) in time steps when returning occurs
     */
    virtual void setAutotunerParams(bool enable, unsigned int period)
    {
        VelocityVerlet::setAutotunerParams(enable, period);
        m_tuner_one->setPeriod(period);
        m_tuner_one->setEnabled(enable);
        m_tuner_two->setPeriod(period);
        m_tuner_two->setEnabled(enable);
    }
    
private:
    boost::scoped_ptr<Autotuner> m_tuner_one; //!< Autotuner for block size (step one kernel)
    boost::scoped_ptr<Autotuner> m_tuner_two; //!< Autotuner for block size (step two kernel)
};

//! Exports the VelocityVerletGPU class to python
void export_VelocityVerletGPU();

//#endif // #ifndef __VELOCITY_VERLET_H__

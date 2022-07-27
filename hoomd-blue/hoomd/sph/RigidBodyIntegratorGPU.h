// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "RigidBodyIntegrator.h"

//#ifndef __RIGID_BODY_INTEGRATOR_H__
//#define __RIGID_BODY_INTEGRATOR_H__

/*! \file RigidBodyIntegratorGPU.h
 *    \brief Declares the RigidBodyIntegratorGPU class
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
class RigidBodyIntegratorGPU : public RigidBodyIntegrator
{
public:
    //! Constructs the integration method and associates it with the system
    RigidBodyIntegratorGPU(boost::shared_ptr<SystemDefinition> sysdef,
                           boost::shared_ptr<ParticleGroup> group,
                           boost::shared_ptr<Variant> transvel_x,
                           boost::shared_ptr<Variant> transvel_y,
                           boost::shared_ptr<Variant> transvel_z,
                           boost::shared_ptr<Variant> rotatvel,
                           Scalar pivotpnt_x,
                           Scalar pivotpnt_y,
                           Scalar pivotpnt_z,
                           Scalar rotaxis_x,
                           Scalar rotaxis_y,
                           Scalar rotaxis_z,
                           bool skip_restart=false);
    //! Destructor
    virtual ~RigidBodyIntegratorGPU();
    
    //! Performs the first step of the integration
    virtual void integrateStepOne(unsigned int timestep);
    
    //! Performs the second step of the integration
    virtual void integrateStepTwo(unsigned int timestep);
    
    virtual void setAutotunerParams(bool enable, unsigned int period)
    {
        RigidBodyIntegrator::setAutotunerParams(enable, period);
        m_tuner_one->setPeriod(period);
        m_tuner_one->setEnabled(enable);
        m_tuner_two->setPeriod(period);
        m_tuner_two->setEnabled(enable);
    }
    
private:
    boost::scoped_ptr<Autotuner> m_tuner_one; //!< Autotuner for block size (step one kernel)
    boost::scoped_ptr<Autotuner> m_tuner_two; //!< Autotuner for block size (step two kernel)
};

//! Exports the RigidBodyIntegratorGPU class to python
void export_RigidBodyIntegratorGPU();

//#endif // #ifndef __RIGID_BODY_INTEGRATOR_H__


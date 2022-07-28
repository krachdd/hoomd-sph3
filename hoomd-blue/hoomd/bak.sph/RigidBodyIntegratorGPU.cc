// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "RigidBodyIntegratorGPU.h"
#include "RigidBodyIntegratorGPU.cuh"

using namespace std;

/*! \file RigidBodyIntegratorGPU.h
 *    \brief Contains code for the RigidBodyIntegratorGPU class
 */

/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
 *    \param group The group of particles this integration method is to work on
 */
RigidBodyIntegratorGPU::RigidBodyIntegratorGPU(boost::shared_ptr<SystemDefinition> sysdef,
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
                                               bool skip_restart)
: RigidBodyIntegrator(sysdef, group, transvel_x, transvel_y, transvel_z,
                      rotatvel, pivotpnt_x, pivotpnt_y, pivotpnt_z, rotaxis_x,
                      rotaxis_y, rotaxis_z, skip_restart)
{
    // only one GPU is supported
    if (!m_exec_conf->isCUDAEnabled())
    {
        m_exec_conf->msg->error() << "Creating a RigidBodyIntegratorGPU when CUDA is disabled" << endl;
        throw std::runtime_error("Error initializing RigidBodyIntegratorGPU");
    }
    
    m_exec_conf->msg->notice(5) << "Constructing RigidBodyIntegrator GPU" << endl;
    
    // initialize autotuner
    std::vector<unsigned int> valid_params;
    for (unsigned int block_size = 32; block_size <= 1024; block_size += 32)
        valid_params.push_back(block_size);
    
    m_tuner_one.reset(new Autotuner(valid_params, 5, 100000, "RigidBodyIntegrator_step_one", this->m_exec_conf));
    m_tuner_two.reset(new Autotuner(valid_params, 5, 100000, "RigidBodyIntegrator_step_two", this->m_exec_conf));
}

RigidBodyIntegratorGPU::~RigidBodyIntegratorGPU()
{
    m_exec_conf->msg->notice(5) << "Destroying RigidBodyIntegratorGPU" << endl;
}

/*! \param timestep Current time step
 *    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2 per the velocity verlet
 *          method.
 */
void RigidBodyIntegratorGPU::integrateStepOne(unsigned int timestep)
{
    unsigned int group_size = m_group->getNumMembers();
    
    // profile this step
    if (m_prof)
        m_prof->push("RigidBodyIntegrator Integrate Step 1");
    
    // Evaluate variant expressions
    Scalar transvel_x    = m_transvel_x->getValue(timestep);
    Scalar transvel_y    = m_transvel_y->getValue(timestep);
    Scalar transvel_z    = m_transvel_z->getValue(timestep);
    Scalar rotatvel      = m_rotatvel->getValue(timestep);
    Scalar transvel_x_dt = m_transvel_x->getRate(timestep);
    Scalar transvel_y_dt = m_transvel_y->getRate(timestep);
    Scalar transvel_z_dt = m_transvel_z->getRate(timestep);
    Scalar rotatvel_dt   = m_rotatvel->getRate(timestep);
    
    // access all the needed data
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3> d_accel(m_pdata->getAccelerations(), access_location::device, access_mode::readwrite);
    
    BoxDim box = m_pdata->getBox();
    ArrayHandle< unsigned int > d_index_array(m_group->getIndexArray(), access_location::device, access_mode::read);
    
    // perform the update on the GPU
    m_tuner_one->begin();
    gpu_RigidBody_step_one(d_pos.data,
                           d_vel.data,
                           d_accel.data,
                           d_index_array.data,
                           group_size,
                           box,
                           m_deltaT,
                           m_pivotpnt,
                           m_rotaxis,
                           transvel_x,
                           transvel_y,
                           transvel_z,
                           rotatvel,
                           transvel_x_dt,
                           transvel_y_dt,
                           transvel_z_dt,
                           rotatvel_dt,
                           m_tuner_one->getParam());
    m_tuner_one->end();
    
    if(m_exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();
    
    
    // done profiling
    if (m_prof)
        m_prof->pop(m_exec_conf);
}

/*! \param timestep Current time step
 *    \post particle velocities are moved forward to timestep+1 on the GPU
 */
void RigidBodyIntegratorGPU::integrateStepTwo(unsigned int timestep)
{
    unsigned int group_size = m_group->getNumMembers();
    
    // profile this step
    if (m_prof)
        m_prof->push("RigidBodyIntegrator Integrate GPU Step 2");
    
    ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::readwrite);
    ArrayHandle< unsigned int > d_index_array(m_group->getIndexArray(), access_location::device, access_mode::read);
    
    m_tuner_two->begin();
    //perform the update on the GPU
    gpu_RigidBody_step_two(d_pos.data,
                           d_vel.data,
                           d_index_array.data,
                           group_size,
                           m_deltaT,
                           m_tuner_one->getParam());
    m_tuner_two->end();
    
    if(m_exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();
    
    
    // done profiling
    if (m_prof)
        m_prof->pop(m_exec_conf);
}

void export_RigidBodyIntegratorGPU()
{
    class_<RigidBodyIntegratorGPU, boost::shared_ptr<RigidBodyIntegratorGPU>, bases<RigidBodyIntegrator>, boost::noncopyable> 
         ("RigidBodyIntegratorGPU", init< boost::shared_ptr<SystemDefinition>,
         boost::shared_ptr<ParticleGroup>,
         boost::shared_ptr<Variant>,
         boost::shared_ptr<Variant>,
         boost::shared_ptr<Variant>,
         boost::shared_ptr<Variant>,
         Scalar,
         Scalar,
         Scalar,
         Scalar,
         Scalar,
         Scalar,
         bool
         >())
    .def("setRotationSpeed",&RigidBodyIntegratorGPU::setRotationSpeed)
    .def("setPivotPoint",&RigidBodyIntegratorGPU::setPivotPoint)
    .def("setRotationAxis",&RigidBodyIntegratorGPU::setRotationAxis)
    .def("setTranslationalVelocity",&RigidBodyIntegratorGPU::setTranslationalVelocity)
    ;
}

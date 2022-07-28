// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon


#include "VelocityVerletGPU.h"
#include "VelocityVerletGPU.cuh"

using namespace std;

/*! \file VelocityVerletGPU.h
 *    \brief Contains code for the VelocityVerletGPU class
 */

/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
 *    \param group The group of particles this integration method is to work on
 */
VelocityVerletGPU::VelocityVerletGPU(boost::shared_ptr<SystemDefinition> sysdef,
                                     boost::shared_ptr<ParticleGroup> group, 
                                     bool skip_restart)
: VelocityVerlet(sysdef, group, skip_restart)
{
    // only one GPU is supported
    if (!m_exec_conf->isCUDAEnabled())
    {
        m_exec_conf->msg->error() << "Creating a VelocityVerletGPU when CUDA is disabled" << endl;
        throw std::runtime_error("Error initializing VelocityVerletGPU");
    }
    
    // initialize autotuner
    std::vector<unsigned int> valid_params;
    for (unsigned int block_size = 32; block_size <= 1024; block_size += 32)
        valid_params.push_back(block_size);
    
    m_tuner_one.reset(new Autotuner(valid_params, 5, 100000, "VelocityVerlet_step_one", this->m_exec_conf));
    m_tuner_two.reset(new Autotuner(valid_params, 5, 100000, "VelocityVerlet_step_two", this->m_exec_conf));
}

/*! \param timestep Current time step
 *    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2 per the velocity verlet
 *          method.
 */
void VelocityVerletGPU::integrateStepOne(unsigned int timestep)
{
    unsigned int group_size = m_group->getNumMembers();
    
    // profile this step
    if (m_prof)
        m_prof->push(m_exec_conf, "SPH Integrate Step 1");
    
    // access all the needed data
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3> d_accel(m_pdata->getAccelerations(), access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_dpe(m_pdata->getDPEs(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3> d_dpedt(m_pdata->getDPEdts(), access_location::device, access_mode::read);
    ArrayHandle<int3> d_image(m_pdata->getImages(), access_location::device, access_mode::readwrite);
    
    BoxDim box = m_pdata->getBox();
    ArrayHandle< unsigned int > d_index_array(m_group->getIndexArray(), access_location::device, access_mode::read);
    
    // perform the update on the GPU
    m_tuner_one->begin();
    gpu_VelocityVerlet_step_one(d_pos.data,
                                d_vel.data,
                                d_accel.data,
                                d_image.data,
                                d_dpe.data,
                                d_dpedt.data,
                                d_index_array.data,
                                group_size,
                                box,
                                m_deltaT,
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
void VelocityVerletGPU::integrateStepTwo(unsigned int timestep)
{
    unsigned int group_size = m_group->getNumMembers();
    
    const GPUArray< Scalar4 >& net_force = m_pdata->getNetForce();
    
    // profile this step
    if (m_prof)
        m_prof->push(m_exec_conf, "VelocityVerlet step 2");
    
    ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3> d_accel(m_pdata->getAccelerations(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3> d_dpe(m_pdata->getDPEs(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3> d_dpedt(m_pdata->getDPEdts(), access_location::device, access_mode::readwrite);
    
    ArrayHandle<Scalar4> d_net_force(net_force, access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_net_ratedpe(m_pdata->getNetRateDPEArray(), access_location::device, access_mode::read);
    ArrayHandle< unsigned int > d_index_array(m_group->getIndexArray(), access_location::device, access_mode::read);

    m_tuner_two->begin();
    // perform the update on the GPU
    gpu_VelocityVerlet_step_two(d_pos.data,
                                d_vel.data,
                                d_accel.data,
                                d_dpe.data,
                                d_dpedt.data,
                                d_index_array.data,
                                group_size,
                                d_net_force.data,
                                d_net_ratedpe.data,
                                m_deltaT,
                                m_tuner_two->getParam()
    );
    m_tuner_two->end();
    
    if(m_exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();
    
    // done profiling
    if (m_prof)
        m_prof->pop(m_exec_conf);
}

void export_VelocityVerletGPU()
{
    class_<VelocityVerletGPU, boost::shared_ptr<VelocityVerletGPU>, bases<VelocityVerlet>, boost::noncopyable>
    ("VelocityVerletGPU", init< boost::shared_ptr<SystemDefinition>, boost::shared_ptr<ParticleGroup>, bool >())
    ;
}

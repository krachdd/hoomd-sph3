// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "SuspendedObjectIntegratorGPU.h"
#include "SuspendedObjectIntegratorGPU.cuh"
#include "GeneralGPUFunctions.cuh"


#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#include "hoomd/HOOMDMPI.h"
#endif

using namespace std;

/*! \file SuspendedObjectIntegratorGPU.h
 *    \brief Contains code for the SuspendedObjectIntegratorGPU class
 */

/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
 *    \param group The group of particles this integration method is to work on
 *    \param skip_restart Skip initialization of the restart information
 */
SuspendedObjectIntegratorGPU::SuspendedObjectIntegratorGPU(boost::shared_ptr<SystemDefinition> sysdef,
                                                           boost::shared_ptr<ParticleGroup> group,
                                                           bool skip_restart)

: SuspendedObjectIntegrator(sysdef, group, skip_restart)
{
    if (!m_exec_conf->isCUDAEnabled())
    {
        m_exec_conf->msg->error() << "Creating a SuspendedObjectIntegratorGPU with CUDA disabled" << endl;
        throw std::runtime_error("Error initializing SuspendedObjectIntegratorGPU");
    }
    
    m_exec_conf->msg->notice(5) << "Constructing SuspendedObjectIntegratorGPU" << endl;
    
    m_reduction_block_size = 512;
    
    // this breaks memory scaling (calculate memory requirements from global group size)
    // unless we reallocate memory with every change of the maximum particle number
    m_num_blocks = m_group->getNumMembersGlobal() / m_reduction_block_size + 1;
    GPUArray< Scalar > scratch(m_num_blocks, m_exec_conf);
    m_scratch.swap(scratch);
    
    //GPUArray< Scalar> temperature(1, m_exec_conf);
    //m_temperature.swap(temperature);
}

SuspendedObjectIntegratorGPU::~SuspendedObjectIntegratorGPU()
{
    m_exec_conf->msg->notice(5) << "Destroying SuspendedObjectIntegratorGPU" << endl;
}

//TODO: spaeter wenn -mode=gpu Compute... routines

/*! \param timestep Current time step
 *    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2 per the Nose-Hoover
 *     thermostat and Anderson barostat
 */
void SuspendedObjectIntegratorGPU::integrateStepOne(unsigned int timestep)
{
    if (m_group->getNumMembersGlobal() == 0)
    {
        m_exec_conf->msg->error() << "integrate.sph(): Integration group empty." << std::endl;
        throw std::runtime_error("Error during SPH integration.");
    }
    
    unsigned int group_size = m_group->getNumMembers();
    
    // profile this step
    if (m_prof)
        m_prof->push("Suspended Object Integrate Step 1");
    
    // Local copy of the simulation box
    const BoxDim& box = m_pdata->getGlobalBox();
    
    m_translationvel += Scalar(1.0/2.0)*m_translationaccel*m_deltaT;
    m_angularvel += Scalar(1.0/2.0)*m_angularaccel*m_deltaT;
    vec3<Scalar> angularvel(m_angularvel.x, m_angularvel.y, m_angularvel.z);
    
    // Array handles
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::readwrite);
    ArrayHandle< unsigned int > d_index_array(m_group->getIndexArray(), access_location::device, access_mode::read);
    
    gpu_SuspendedOIntegrator_step_one(d_pos.data,
                                      d_vel.data,
                                      d_index_array.data,
                                      group_size,
                                      box,
                                      m_centerofmass,
                                      angularvel,
                                      m_translationvel, m_deltaT);
    
    #ifdef ENABLE_MPI
    MPI_Barrier(m_exec_conf->getMPICommunicator());
    #endif
    
    // Recompute center of mass
    ComputeCenterOfMass(false);
    
    // done profiling
    if (m_prof)
        m_prof->pop();
    
}

/*! \param timestep Current time step
 *    \post particle velocities are moved forward to timestep+1
 */
void SuspendedObjectIntegratorGPU::integrateStepTwo(unsigned int timestep)
{
    unsigned int group_size = m_group->getNumMembers();
    
    //const GPUArray< Scalar4 >& net_force = m_pdata->getNetForce();
    
    // profile this step
    if (m_prof)
        m_prof->push("Suspended Object Integrate Step 2");
    
    
    // Local copy of the simulation box
    const BoxDim& box = m_pdata->getGlobalBox();
    
    // Get array handles
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_net_force(m_pdata->getNetForce(), access_location::device, access_mode::read);
    ArrayHandle< unsigned int > d_index_array(m_group->getIndexArray(), access_location::device, access_mode::read);
    
    // Initialize total force and torque vectors
    GPUArray<Scalar> totalforcetorque(6, m_exec_conf);
    
    {
        ArrayHandle<Scalar> d_totalforcetorque(totalforcetorque, access_location::device, access_mode::readwrite);
        // Loop over group index
        gpu_SuspendedOIntegrator_step_two_a(d_pos.data,
                                            d_net_force.data,
                                            d_index_array.data,
                                            group_size,
                                            box,
                                            m_centerofmass,
                                            d_totalforcetorque.data);
    }
    // TODO gpu cpu sync
    ArrayHandle<Scalar> h_totalforcetorque(totalforcetorque, access_location::host, access_mode::readwrite);
    
    #ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(MPI_IN_PLACE, h_totalforcetorque.data, 6, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
    #endif
    
    // Update slave particle velocities
    ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::overwrite);
    
    // Compute translational acceleration
    m_translationaccel.x = h_totalforcetorque.data[0]/m_totalmass;
    m_translationaccel.y = h_totalforcetorque.data[1]/m_totalmass;
    m_translationaccel.z = h_totalforcetorque.data[2]/m_totalmass;
    
    // Compute angular acceleration
    m_angularaccel.x = m_invJ[0]*h_totalforcetorque.data[3]+m_invJ[1]*h_totalforcetorque.data[4]+m_invJ[2]*h_totalforcetorque.data[5];
    m_angularaccel.y = m_invJ[3]*h_totalforcetorque.data[3]+m_invJ[4]*h_totalforcetorque.data[4]+m_invJ[5]*h_totalforcetorque.data[5];
    m_angularaccel.z = m_invJ[6]*h_totalforcetorque.data[3]+m_invJ[7]*h_totalforcetorque.data[4]+m_invJ[8]*h_totalforcetorque.data[5];
    
    // Integrate Object variables
    m_translationvel += Scalar(1.0/2.0)*m_translationaccel*m_deltaT;
    m_angularvel += Scalar(1.0/2.0)*m_angularaccel*m_deltaT;
    vec3<Scalar> angularvel(m_angularvel.x, m_angularvel.y, m_angularvel.z);
    
    gpu_SuspendedOIntegrator_step_two_b(d_vel.data,
                                        d_pos.data,
                                        d_index_array.data,
                                        group_size,
                                        box,
                                        m_centerofmass,
                                        m_translationvel,
                                        angularvel);
    
    // done profiling
    if (m_prof)
        m_prof->pop();
}



void SuspendedObjectIntegratorGPU::ComputeTotalMass(bool print = false)
    {
    // Initialize mass to zero
    m_totalmass = 0;
    // intermediate result for gpu reduction
    GPUArray<Scalar> total_mass(1, m_exec_conf);
    // no need for total_mass.memclear() as the memory is set to zero during the allcoation above

    { // scope to limit lifetime of d_total_mass as we have to read back the result to the cpu
        unsigned int group_size = m_group->getNumMembers();
        ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::read);
        ArrayHandle<unsigned int> d_indices(m_group->getIndexArray(), access_location::device, access_mode::read);
        ArrayHandle<Scalar> d_total_mass(total_mass, access_location::device, access_mode::readwrite); // rw as we use it as accumulator

        gpu_sum_x(d_vel.data, group_size, d_total_mass.data, d_indices.data);
    }
    // copy back the result to cpu
    // TODO cpu-gpu sync! we are constraint to the interface... can we avoid this somehow?
    ArrayHandle<Scalar> h_total_mass(total_mass, access_location::host, access_mode::read);
    m_totalmass = h_total_mass.data[0];


    // TODO keep this???? this was in the cpu version. Not sure if this is needed for gpu too
#ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(MPI_IN_PLACE, &m_totalmass, 1, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
#endif

    if ( print )
        m_exec_conf->msg->notice(5) << "SuspendedObject Total Mass: " << m_totalmass << endl;
    }

void SuspendedObjectIntegratorGPU::ComputeCenterOfMass(bool print = false)
    {
    // profile this step
    if (m_prof)
        m_prof->push("Suspended Object Compute Center of Mass");

    // Read-Handle to velocity mass and position array
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);

    // grab the box dimensions
    BoxDim globalbox = m_pdata->getGlobalBox();

    // Limit values
    Scalar3 lo = globalbox.getLo();
    Scalar3 Ld = globalbox.getL();

    // Average position on unit circle
    GPUArray<Scalar> angles(6, m_exec_conf); // TODO cache this?
    // allocation sets memory to zeros

    unsigned int group_size = m_group->getNumMembers();
    {
        ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
        ArrayHandle<unsigned int> d_indices(m_group->getIndexArray(), access_location::device, access_mode::read);
        ArrayHandle<Scalar> d_angles(angles, access_location::device, access_mode::readwrite); // rw as we use it as accumulator
        
        gpu_SuspendedOIntegrator_compute_center_of_mass(d_pos.data, d_indices.data, group_size, d_angles.data, lo, Ld);
    }

    // Count total number of particles in group
    unsigned int totalgroupN = group_size;
    // TODO cpu-gpu sync
    ArrayHandle<Scalar> h_angles(angles, access_location::host, access_mode::read);
    // move data out of h_angles as we will modify it.
    Scalar angles_total[6];

#ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(&h_angles.data[0], &angles_total[0], 6, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
    MPI_Allreduce(MPI_IN_PLACE, &totalgroupN, 1, MPI_UNSIGNED, MPI_SUM, this->m_exec_conf->getMPICommunicator());
#else
    memcpy(angles_total, h_angles.data, 6*sizeof(Scalar));
#endif

    // Averaging
    angles_total[0] /= totalgroupN;
    angles_total[1] /= totalgroupN;
    angles_total[2] /= totalgroupN;
    angles_total[3] /= totalgroupN;
    angles_total[4] /= totalgroupN;
    angles_total[5] /= totalgroupN;

    // Set attribute
    m_centerofmass = make_scalar3(0.0, 0.0, 0.0);
    if ( angles_total[1] != 0 && angles_total[0] != 0 )
        m_centerofmass.x = lo.x + ( atan2(-angles_total[1], -angles_total[0]) + Scalar(M_PI) ) * ( Ld.x / Scalar(M_TWOPI) );
    if ( angles_total[3] != 0 && angles_total[2] != 0 )
        m_centerofmass.y = lo.y + ( atan2(-angles_total[3], -angles_total[2]) + Scalar(M_PI) ) * ( Ld.y / Scalar(M_TWOPI) );
    if ( angles_total[5] != 0 && angles_total[4] != 0 )
        m_centerofmass.z = lo.z + ( atan2(-angles_total[5], -angles_total[4]) + Scalar(M_PI) ) * ( Ld.z / Scalar(M_TWOPI) );

    if ( print )
        {
        m_exec_conf->msg->notice(5) << "SuspendedObject consists of " << totalgroupN << " slave particles" << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject center of Mass x: " << m_centerofmass.x << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject center of mass y: " << m_centerofmass.y << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject center of mass z: " << m_centerofmass.z << endl;
        }

    // done profiling
    if (m_prof)
        m_prof->pop();
    }

void SuspendedObjectIntegratorGPU::ComputeTranslationVelocity(bool print = false)
    {
    // Initialize center of mass to zero
    GPUArray<Scalar> translationvel(3, m_exec_conf); // TODO cache
    // set to zero during allocation

    {
        // Read-Handle to velocity mass array
        ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::read);
        ArrayHandle<Scalar> d_translationvel(translationvel, access_location::device, access_mode::readwrite);
        ArrayHandle<unsigned int> d_indices(m_group->getIndexArray(), access_location::device, access_mode::read);
        unsigned int group_size = m_group->getNumMembers();
        
        gpu_SuspendedOIntegrator_compute_translation_velocity(d_vel.data, m_totalmass, d_indices.data, group_size, d_translationvel.data);
    }
    // TODO gpu cpu sync
    ArrayHandle<Scalar> h_translationvel(translationvel, access_location::host, access_mode::readwrite);


#ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(MPI_IN_PLACE, &h_translationvel.data[0], 3, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
#endif

    // Set attribute
    m_translationvel.x = h_translationvel.data[0];
    m_translationvel.y = h_translationvel.data[1];
    m_translationvel.z = h_translationvel.data[2];

    if ( print )
        {
        m_exec_conf->msg->notice(5) << "SuspendedObject translation velocity x: " << m_translationvel.x << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject translation velocity y: " << m_translationvel.y << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject translation velocity z: " << m_translationvel.z << endl;
        }
    }

void SuspendedObjectIntegratorGPU::ComputeMomentOfInertia(bool print = false)
    {
    // Local copy of the simulation box
    const BoxDim& box = m_pdata->getGlobalBox();

    // Initialize moment of inertia tensor to zero
    GPUArray<Scalar> J(9, m_exec_conf);

    {
        // Read-Handle to velocity mass and position array
        ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
        ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::read);
        ArrayHandle<unsigned int> d_indices(m_group->getIndexArray(), access_location::device, access_mode::read);
        ArrayHandle<Scalar> d_J(J, access_location::device, access_mode::readwrite);
        unsigned int group_size = m_group->getNumMembers();

        gpu_SuspendedOIntegrator_compute_moment_of_intertia(
            d_pos.data, d_vel.data, 
            d_indices.data, group_size, 
            m_centerofmass, box, 
            d_J.data
        );
    }
    // TODO cpu gpu sync
    ArrayHandle<Scalar> h_J(J, access_location::host, access_mode::readwrite);

#ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(MPI_IN_PLACE, &h_J.data[0], 9, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
#endif

    // Compute inverse of moment of inertia tensor and store as class attribute
    Scalar detJ = +h_J.data[0]*(h_J.data[4]*h_J.data[8]-h_J.data[7]*h_J.data[5])
                  -h_J.data[1]*(h_J.data[3]*h_J.data[8]-h_J.data[5]*h_J.data[6])
                  +h_J.data[2]*(h_J.data[3]*h_J.data[7]-h_J.data[4]*h_J.data[6]);
    Scalar invdetJ = 1/detJ;

    m_invJ[0] =   (h_J.data[4]*h_J.data[8]-h_J.data[7]*h_J.data[5])*invdetJ;
    m_invJ[1] =  -(h_J.data[1]*h_J.data[8]-h_J.data[2]*h_J.data[7])*invdetJ;
    m_invJ[2] =   (h_J.data[1]*h_J.data[5]-h_J.data[2]*h_J.data[4])*invdetJ;
    m_invJ[3] =  -(h_J.data[3]*h_J.data[8]-h_J.data[5]*h_J.data[6])*invdetJ;
    m_invJ[4] =   (h_J.data[0]*h_J.data[8]-h_J.data[2]*h_J.data[6])*invdetJ;
    m_invJ[5] =  -(h_J.data[0]*h_J.data[5]-h_J.data[3]*h_J.data[2])*invdetJ;
    m_invJ[6] =   (h_J.data[3]*h_J.data[7]-h_J.data[6]*h_J.data[4])*invdetJ;
    m_invJ[7] =  -(h_J.data[0]*h_J.data[7]-h_J.data[6]*h_J.data[1])*invdetJ;
    m_invJ[8] =   (h_J.data[0]*h_J.data[4]-h_J.data[3]*h_J.data[1])*invdetJ;

    if ( print )
        {
        m_exec_conf->msg->notice(5) << "SuspendedObject moment of inertia: (0,0):" << h_J.data[0] << " (0,1):" << h_J.data[1] << " (0,2):" << h_J.data[2] << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject moment of inertia: (1,0):" << h_J.data[3] << " (1,1):" << h_J.data[4] << " (1,2):" << h_J.data[5] << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject moment of inertia: (2,0):" << h_J.data[6] << " (2,1):" << h_J.data[7] << " (2,2):" << h_J.data[8] << endl;
        }
    }

void SuspendedObjectIntegratorGPU::ComputeAngularVelocity(bool print = false)
    {
    // Local copy of the simulation box
    const BoxDim& box = m_pdata->getGlobalBox();

    // Initialize center of mass to zero
    GPUArray<Scalar> angmomentum(3, m_exec_conf);

    {
        // Read-Handle to velocity mass and position array
        ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
        ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::read);
        ArrayHandle<Scalar> d_angmomentum(angmomentum, access_location::device, access_mode::readwrite);
        ArrayHandle<unsigned int> d_indices(m_group->getIndexArray(), access_location::device, access_mode::read);
        unsigned int group_size = m_group->getNumMembers();
        
        gpu_suspendedOIntegrator_compute_angular_velocity(
            d_pos.data, d_vel.data, 
            m_centerofmass, m_translationvel, box,
            d_indices.data, group_size, 
            d_angmomentum.data
        );
    }
    // TODO gpu cpu sync
    ArrayHandle<Scalar> h_angmomentum(angmomentum, access_location::host, access_mode::readwrite);

#ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(MPI_IN_PLACE, &h_angmomentum.data[0], 3, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
#endif

    // Compute angular velocity using inverted moment of inertia tensor
    m_angularvel.x = m_invJ[0]*h_angmomentum.data[0]+m_invJ[1]*h_angmomentum.data[1]+m_invJ[2]*h_angmomentum.data[2];
    m_angularvel.y = m_invJ[3]*h_angmomentum.data[0]+m_invJ[4]*h_angmomentum.data[1]+m_invJ[5]*h_angmomentum.data[2];
    m_angularvel.z = m_invJ[6]*h_angmomentum.data[0]+m_invJ[7]*h_angmomentum.data[1]+m_invJ[8]*h_angmomentum.data[2];

    if ( print )
        {
        m_exec_conf->msg->notice(5) << "SuspendedObject initial angular velocity x: " << m_angularvel.x << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject initial angular velocity y: " << m_angularvel.y << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject initial angular velocity z: " << m_angularvel.z << endl;
        }
    }


void export_SuspendedObjectIntegratorGPU()
{
    class_<SuspendedObjectIntegratorGPU, boost::shared_ptr<SuspendedObjectIntegratorGPU>, bases<SuspendedObjectIntegrator>, boost::noncopyable>
    ("SuspendedObjectIntegratorGPU", init< boost::shared_ptr<SystemDefinition>, boost::shared_ptr<ParticleGroup>, bool >())
    ;
}

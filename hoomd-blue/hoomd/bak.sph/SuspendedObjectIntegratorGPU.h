// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "SuspendedObjectIntegrator.h"
#include "hoomd/Variant.h"

//#ifndef __SUSPENDED_OBJECT_INTEGRATOR_H__
//#define __SUSPENDED_OBJECT_INTEGRATOR_H__

/*! \file SuspendedObjectIntegratorGPU.h
 *    \brief Declares an Velocity-Verlet integrator for suspended objects.
 *           For every time step.
 *              1. Sum up all forces acting on boundary particles
 *                 -> These forces have been computed in QINSSPForceCompute assuming force symmetry.
 *              2. Sum up total force and total torque acting on suspended rigid object.
 *                 -> MPI_ALLREDUCE will be necessary here.
 *              3. Integrate object center of mass, translation velocity and angular velocity in time
 *                 -> Probably it is cheapest to perform this computation on all processors.
 *                    rather than broadcasting the result to all slave processors.
 *                 -> Alternatively, recompute center of mass after advecting boundary particles.
 *              4. Compute velocity of boundary particles and integrate their position in time.
 *           At initialization time:
 *              a) Compute total mass (remains constant during simulation)
 *              b) Compute moment of inertia tensor and inverse (remains constant during simulation)
 *              c) Compute initial translational velocity (integrated in time)
 *              d) Compute initial angular velocity (integrated in time)
 *              e) Compute initial center of mass (integrated in time)
 *              -> MPI_ALLREDUCE operations will be necessary here.
 *            For more information on rigid body dynamics see:
 *            https://en.wikipedia.org/wiki/Rigid_body_dynamics
 */

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif


//! Integrates part of the system forward in two steps in the NPT ensemble
/*! This is a version of SuspendedObjectIntegrator that runs on the GPU.
 * 
 *    \ingroup updaters
 */
class SuspendedObjectIntegratorGPU : public SuspendedObjectIntegrator
{
public:
    //! Constructs the integration method and associates it with the system
    SuspendedObjectIntegratorGPU(boost::shared_ptr<SystemDefinition> sysdef,
                                 boost::shared_ptr<ParticleGroup> group,
                                 bool skip_restart=false);
    
    
    virtual ~SuspendedObjectIntegratorGPU();
    
    //! Performs the first step of the integration
    virtual void integrateStepOne(unsigned int timestep);
    
    //! Performs the second step of the integration
    virtual void integrateStepTwo(unsigned int timestep);
    
protected:
    GPUArray<Scalar> m_scratch;     //!< Scratch space for reduction of squared velocities
    GPUArray<Scalar> m_temperature; //!< Stores temperature after reduction step

    // TODO can we replace these by GPUArray and reduce gpu-cpu-gpu roundtrips?
    //      with MPI: no
    //      w/o  MPI: m_totalmass OK
    //                m_centerofmass minor postprocessing
    //                m_translationvel OK
    //                m_invJ minor postprocessing 
    //                m_angularvel minor postprocessing 
    Scalar  m_totalmass; //!< Total mass of object
    Scalar3 m_centerofmass; //!< Center of mass vector
    Scalar3 m_translationvel; //!< Translational velocity
    //Scalar3 m_translationaccel; //!< Translational acceleration
    //Scalar3 m_angularaccel; //!< Angular acceleration
    Scalar  m_invJ[9]; //!< Inverted moment of inertia tensor
    Scalar3 m_angularvel; //!< Angular velocity

    // TODO how should this be used?
    unsigned int m_num_blocks;             //!< Number of blocks participating in the reduction
    unsigned int m_reduction_block_size;   //!< Block size executed

    // TODO make these virtual
    //! Compute total mass of suspended object
    void ComputeTotalMass(bool print);

    //! Compute center of mass of suspended object
    void ComputeCenterOfMass(bool print);

    //! Compute translational velocity of suspended object
    void ComputeTranslationVelocity(bool print);

    //! Compute moment of inertia tensor of suspended object
    void ComputeMomentOfInertia(bool print);

    //! Compute initial angular velocity of suspended object
    void ComputeAngularVelocity(bool print);
};

//! Exports the SuspendedObjectIntegratorGPU class to python
void export_SuspendedObjectIntegratorGPU();

//#endif // #ifndef __SUSPENDED_OBJECT_INTEGRATOR_H__

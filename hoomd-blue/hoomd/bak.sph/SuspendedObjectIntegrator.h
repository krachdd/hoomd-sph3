// Copyright (c) 2009-2022 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

#include "SPHIntegrationMethodTwoStep.h"
#include <hoomd/HOOMDMath.h>
#include <hoomd/VectorMath.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#ifndef __SUSPENDED_OBJECT_INTEGRATOR_H__
#define __SUSPENDED_OBJECT_INTEGRATOR_H__

/*! \file SuspendedObjectIntegrator.h
    \brief Declares an Velocity-Verlet integrator for suspended objects.
           For every time step.
              1. Sum up all forces acting on boundary particles
                 -> These forces have been computed in QINSSPForceCompute assuming force symmetry.
              2. Sum up total force and total torque acting on suspended rigid object.
                 -> MPI_ALLREDUCE will be necessary here.
              3. Integrate object center of mass, translation velocity and angular velocity in time
                 -> Probably it is cheapest to perform this computation on all processors.
                    rather than broadcasting the result to all slave processors.
                 -> Alternatively, recompute center of mass after advecting boundary particles.
              4. Compute velocity of boundary particles and integrate their position in time.
           At initialization time:
              a) Compute total mass (remains constant during simulation)
              b) Compute moment of inertia tensor and inverse (remains constant during simulation)
              c) Compute initial translational velocity (integrated in time)
              d) Compute initial angular velocity (integrated in time)
              e) Compute initial center of mass (integrated in time)
              -> MPI_ALLREDUCE operations will be necessary here.
            For more information on rigid body dynamics see:
            https://en.wikipedia.org/wiki/Rigid_body_dynamics
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_TWOPI
#define M_TWOPI 6.28318530717958647692
#endif

//! Integrates the a suspended object forward one step
/*! Implements time integration for suspended particles through the IntegrationMethodTwoStep interface
    \ingroup updaters
*/

namespace hoomd
{
namespace sph
{
class SuspendedObjectIntegrator : public SPHIntegrationMethodTwoStep
    {
    public:
        //! Constructor
        SuspendedObjectIntegrator(std::shared_ptr<SystemDefinition> sysdef,
                                  std::shared_ptr<ParticleGroup> group,
                                  bool skip_restart=false);

        //! Destructor
        virtual ~SuspendedObjectIntegrator();

        //! Prepare for the run
        virtual void integrateStepOne(unsigned int timestep);

        //! Prepare for the run
        virtual void integrateStepTwo(unsigned int timestep);

    protected:
        Scalar  m_totalmass; //!< Total mass of object
        Scalar3 m_centerofmass; //!< Center of mass vector
        Scalar3 m_translationvel; //!< Translational velocity
        Scalar3 m_translationaccel; //!< Translational acceleration
        Scalar3 m_angularaccel; //!< Angular acceleration
        Scalar  m_invJ[9]; //!< Inverted moment of inertia tensor
        Scalar3 m_angularvel; //!< Angular velocity

        // TODO make these virtual?
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

namespace detail 
{
//! Exports the SuspendedObjectIntegrator class to python
void export_SuspendedObjectIntegrator(pybind11::module& m);

} // end namespace detail
} // end namespace sph 
} // end namespace hoomd 

#endif // #ifndef __SUSPENDED_OBJECT_INTEGRATOR_H__

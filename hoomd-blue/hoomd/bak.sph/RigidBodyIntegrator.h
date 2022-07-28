// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

#include "SPHIntegrationMethodTwoStep.h"
#include <hoomd/HOOMDMath.h>
#include <hoomd/VectorMath.h>
#include <hoomd/Variant.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#ifndef __RIGID_BODY_INTEGRATOR_H__
#define __RIGID_BODY_INTEGRATOR_H__

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

//! Integrates a rigid body with prescribed motion forward one step
/*! Implements time integration through the IntegrationMethodTwoStep interface
    \ingroup updaters
*/
namespace hoomd 
{
namespace sph
{
class RigidBodyIntegrator : public SPHIntegrationMethodTwoStep
    {
    public:
        //! Constructor
        RigidBodyIntegrator(std::shared_ptr<SystemDefinition> sysdef,
                            std::shared_ptr<ParticleGroup> group,
                            std::shared_ptr<Variant> transvel_x,
                            std::shared_ptr<Variant> transvel_y,
                            std::shared_ptr<Variant> transvel_z,
                            std::shared_ptr<Variant> rotatvel,
                            Scalar pivotpnt_x,
                            Scalar pivotpnt_y,
                            Scalar pivotpnt_z,
                            Scalar rotaxis_x,
                            Scalar rotaxis_y,
                            Scalar rotaxis_z,
                            bool skip_restart=false);

        //! Destructor
        virtual ~RigidBodyIntegrator();

        //! Update the angular velocity
        /*! \param angvel New angular velocity to set
        */
        virtual void setRotationSpeed(std::shared_ptr<Variant> rotatvel)
            {
            m_rotatvel = rotatvel;
            }

        //! Update the rotational pivot point
        /*! \param pivotpoint_x New pivot point x-coordinate
          ! \param pivotpoint_y New pivot point y-coordinate
          ! \param pivotpoint_z New pivot point z-coordinate
        */
        virtual void setPivotPoint(Scalar pivotpnt_x,
                                   Scalar pivotpnt_y,
                                   Scalar pivotpnt_z)
            {
            m_pivotpnt = make_scalar3(pivotpnt_x,pivotpnt_y,pivotpnt_z);
            }

        //! Update the rotational axis
        /*! \param rotaxis_x New rotation axis x-component
          ! \param rotaxis_y New rotation axis y-component
          ! \param rotaxis_z New rotation axis z-component
        */
        virtual void setRotationAxis(Scalar rotaxis_x,
                                     Scalar rotaxis_y,
                                     Scalar rotaxis_z)
            {
            vec3<Scalar> rotaxis(rotaxis_x,rotaxis_y,rotaxis_z);
            rotaxis   = rotaxis/Scalar(sqrt(dot(rotaxis,rotaxis)));
            m_rotaxis = vec_to_scalar3(rotaxis);
            }


        //! Update the angular velocity
        /*! \param angvel New angular velocity to set
        */
        virtual void setTranslationalVelocity(std::shared_ptr<Variant> transvel_x,
                                              std::shared_ptr<Variant> transvel_y,
                                              std::shared_ptr<Variant> transvel_z)
            {
            m_transvel_x = transvel_x;
            m_transvel_y = transvel_y;
            m_transvel_z = transvel_z;
            }

        //! Prepare for the run
        virtual void integrateStepOne(unsigned int timestep);

        //! Prepare for the run
        virtual void integrateStepTwo(unsigned int timestep);

    protected:
        Scalar3 m_pivotpnt; //!< Rotational pivot point
        Scalar3 m_rotaxis; //!< Rotational axis
        std::shared_ptr<Variant> m_transvel_x; //!< Translational velocity x-component
        std::shared_ptr<Variant> m_transvel_y; //!< Translational velocity y-component
        std::shared_ptr<Variant> m_transvel_z; //!< Translational velocity z-component
        std::shared_ptr<Variant> m_rotatvel; //!< Angular velocity
    };

namespace detail 
{

//! Exports the SuspendedObjectIntegrator class to python
void export_RigidBodyIntegrator(pybind11::module& m);

}   // end namespace detail 
}   // end namespace sph 
}   // end namespace hoomd 

#endif // #ifndef __RIGID_BODY_INTEGRATOR_H__

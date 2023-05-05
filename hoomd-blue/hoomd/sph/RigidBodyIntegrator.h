/* ---------------------------------------------------------
maintainer: drostan, daniel.rostan@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SPHIntegrationMethodTwoStep.h"
#include "EvaluationMethodDefinition.h"
#include <hoomd/HOOSPHMath.h>
#include <hoomd/VectorMath.h>
#include <hoomd/Variant.h>

#ifndef __RIGID_BODY_INTEGRATOR_H__
#define __RIGID_BODY_INTEGRATOR_H__

/*! \file RigidBodyIntegrator.h
    \brief Declares the RigidBodyIntegrator class
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace sph
    {
//! Integrates part of the system forward in two steps in the NVE ensemble
/*! Implements velocity-verlet NVE integration through the IntegrationMethodTwoStep interface

    \ingroup updaters
*/
class PYBIND11_EXPORT RigidBodyIntegrator : public SPHIntegrationMethodTwoStep
    {
    public:
    //! Constructs the integration method and associates it with the system
    VelocityVerlet(std::shared_ptr<SystemDefinition> sysdef, 
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
                   Scalar rotaxis_z,);
    virtual ~RigidBodyIntegrator();

    //! Update the angular velocity
    /*! \param angvel New angular velocity to set
    */
    virtual void setRotationSpeed(st::shared_ptr<Variant> rotatvel)
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
    virtual void setTranslationalVelocity(boost::shared_ptr<Variant> transvel_x,
                                          boost::shared_ptr<Variant> transvel_y,
                                          boost::shared_ptr<Variant> transvel_z)
        {
        m_transvel_x = transvel_x;
        m_transvel_y = transvel_y;
        m_transvel_z = transvel_z;
        }

    //! Performs the first step of the integration
    virtual void integrateStepOne(uint64_t timestep);

    //! Performs the second step of the integration
    virtual void integrateStepTwo(uint64_t timestep);

    // //! Getter and Setter methods for density method
    // DensityMethod getDensityMethod()
    //     {
    //     return m_density_method;
    //     }
    // void setDensityMethod(DensityMethod densitymethod)
    //     {
    //     m_density_method = densitymethod;
    //     m_densitymethod_set = true;
    //     }

    protected:
    // bool m_limit;       //!< True if we should limit the distance a particle moves in one step
    // Scalar m_limit_val; //!< The maximum distance a particle is to move in one step
    // bool m_zero_force;  //!< True if the integration step should ignore computed forces
    // bool m_densitymethod_set; //!< True if method was set
    // DensityMethod m_density_method; //!< Density approach to use
    Scalar3 m_pivotpnt; //!< Rotational pivot point
    Scalar3 m_rotaxis; //!< Rotational axis
    std::shared_ptr<Variant> m_transvel_x; //!< Translational velocity x-component
    std::shared_ptr<Variant> m_transvel_y; //!< Translational velocity y-component
    std::shared_ptr<Variant> m_transvel_z; //!< Translational velocity z-component
    std::shared_ptr<Variant> m_rotatvel; //!< Angular velocity
    };

    namespace detail 
{
void export_RigidBodyIntegrator(pybind11::module& m);

} // end namespace detail
    } // end namespace sph
    } // end namespace hoomd

#endif // #ifndef __RIGID_BODY_INTEGRATOR_H__


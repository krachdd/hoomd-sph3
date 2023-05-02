/* ---------------------------------------------------------
maintainer: nkijanksi, nadine.kijanksi@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SPHIntegrationMethodTwoStep.h"
#include "EvaluationMethodDefinition.h"

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

#include <pybind11/pybind11.h>

// #ifndef M_PI
// #define M_PI 3.14159265358979323846
// #endif

// #ifndef M_TWOPI
// #define M_TWOPI 6.28318530717958647692
// #endif

namespace hoomd
    {
namespace sph
    {
//! Integrates part of the system forward in two steps in the NVE ensemble
/*! Implements velocity-verlet NVE integration through the IntegrationMethodTwoStep interface

    \ingroup updaters
*/
class PYBIND11_EXPORT SuspendedObjectIntegrator : public SPHIntegrationMethodTwoStep
    {
    public:
    //! Constructs the integration method and associates it with the system
    SuspendedObjectIntegrator(std::shared_ptr<SystemDefinition> sysdef, std::shared_ptr<ParticleGroup> group);
    virtual ~SuspendedObjectIntegrator();

    //! Get the movement limit
    pybind11::object getLimit();

    //! Sets the movement limit
    void setLimit(pybind11::object limit);

    //! Get zero force
    bool getZeroForce();

    //! Sets the zero force option
    /*! \param zero_force Set to true to specify that the integration with a zero net force on each
       of the particles in the group
    */
    void setZeroForce(bool zero_force);

    //! Performs the first step of the integration
    virtual void integrateStepOne(uint64_t timestep);

    //! Performs the second step of the integration
    virtual void integrateStepTwo(uint64_t timestep);

    //! Getter and Setter methods for density method
    DensityMethod getDensityMethod()
        {
        return m_density_method;
        }
    void setDensityMethod(DensityMethod densitymethod)
        {
        m_density_method = densitymethod;
        m_densitymethod_set = true;
        }

    protected:
    bool m_limit;       //!< True if we should limit the distance a particle moves in one step
    Scalar m_limit_val; //!< The maximum distance a particle is to move in one step
    bool m_zero_force;  //!< True if the integration step should ignore computed forces
    bool m_densitymethod_set; //!< True if method was set
    DensityMethod m_density_method; //!< Density approach to use

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
void export_SuspendedObjectIntegrator(pybind11::module& m);

} // end namespace detail
    } // end namespace sph
    } // end namespace hoomd

#endif // #ifndef __SUSPENDED_OBJECT_INTEGRATOR_H__

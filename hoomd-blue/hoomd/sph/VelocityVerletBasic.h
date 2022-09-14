/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SPHIntegrationMethodTwoStep.h"

#ifndef __VELOCITY_VERLET_BASIC_H__
#define __VELOCITY_VERLET_BASIC_H__

/*! \file VelocityVerletBasic.h
    \brief Declares the VelocityVerletBasic class
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
class PYBIND11_EXPORT VelocityVerletBasic : public SPHIntegrationMethodTwoStep
    {
    public:
    //! Constructs the integration method and associates it with the system
    VelocityVerletBasic(std::shared_ptr<SystemDefinition> sysdef, std::shared_ptr<ParticleGroup> group);
    virtual ~VelocityVerletBasic();

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

    protected:
    bool m_limit;       //!< True if we should limit the distance a particle moves in one step
    Scalar m_limit_val; //!< The maximum distance a particle is to move in one step
    bool m_zero_force;  //!< True if the integration step should ignore computed forces
    };

    namespace detail 
{
void export_VelocityVerletBasic(pybind11::module& m);

} // end namespace detail
    } // end namespace sph
    } // end namespace hoomd

#endif // #ifndef __VELOCITY_VERLET_H__


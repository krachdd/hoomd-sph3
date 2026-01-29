/* ---------------------------------------------------------
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SPHIntegrationMethodTwoStep.h"
#include "EvaluationMethodDefinition.h"

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

    };

    namespace detail 
{
void export_VelocityVerletBasic(pybind11::module& m);

} // end namespace detail
    } // end namespace sph
    } // end namespace hoomd

#endif // #ifndef __VELOCITY_VERLET_H__


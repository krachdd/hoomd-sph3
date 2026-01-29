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

#include "hoomd/ForceCompute.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/ParticleGroup.h"
#include "hoomd/VectorMath.h"
#include <memory>

/*! \file ConstantForceCompute.h
    \brief Declares a class for computing constant forces and torques
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>

#ifndef __CONSTANTFORCECOMPUTE_H__
#define __CONSTANTFORCECOMPUTE_H__

namespace hoomd
    {
namespace sph
    {
//! Adds an constant force to a number of particles
/*! \ingroup computes
 */
class PYBIND11_EXPORT ConstantForceCompute : public ForceCompute
    {
    public:
    //! Constructs the compute
    ConstantForceCompute(std::shared_ptr<SystemDefinition> sysdef,
                         std::shared_ptr<ParticleGroup> group);

    //! Destructor
    ~ConstantForceCompute();

    /** Sets constant force vector for a given particle type
        @param typ Particle type to set constant force vector
        @param v The constant force vector value to set (a 3-tuple)
    */
    void setConstantForce(const std::string& type_name, pybind11::tuple v);

    /// Gets constant force vector for a given particle type
    pybind11::tuple getConstantForce(const std::string& type_name);

    /** Sets constant torque vector for a given particle type
        @param typ Particle type to set constant torque vector
        @param v The constant torque vector value to set (a 3-tuple)
    */
    // void setConstantTorque(const std::string& type_name, pybind11::tuple v);

    /// Gets constant torque vector for a given particle type
    // pybind11::tuple getConstantTorque(const std::string& type_name);

    std::shared_ptr<ParticleGroup>& getGroup()
        {
        return m_group;
        }

    protected:
    //! Actually compute the forces
    virtual void computeForces(uint64_t timestep);

    //! Set forces for particles
    virtual void setForces();

    std::shared_ptr<ParticleGroup> m_group; //!< Group of particles on which this force is applied
    GPUVector<Scalar3>
        m_constant_force; //! constant force unit vectors and magnitudes for each particle type

    GPUVector<Scalar3>
        m_constant_torque; //! constant torque unit vectors and magnitudes for each particle type

    bool m_parameters_updated; //!< True if forces need to be rearranged
    };

    } // end namespace sph
    } // end namespace hoomd
#endif

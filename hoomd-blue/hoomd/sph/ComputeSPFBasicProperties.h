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

#include "ComputeMechanicalPropTypes.h"
#include "hoomd/Compute.h"
#include "hoomd/ParticleGroup.h"

#include <limits>
#include <memory>

/*! \file ComputeSPFBasicProperties.h
    \brief Declares a class for computing quantities
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>

#ifndef __COMPUTE_SPF_BASIC_PROPERTIES_H__
#define __COMPUTE_SPF_BASIC_PROPERTIES_H__

namespace hoomd
    {
namespace sph
    {

//! Computes properties of a group of particles
/*! ComputeSPFMechanical properties calculates instantaneous properties and provides them in Python.
    All computed values are stored in a GPUArray so that they can be accessed on the GPU without
   intermediate copies. Use the enum values in singlephaseflow_logger_index to index the array and extract the
   properties of interest. Convenience functions are provided for accessing the values on the CPU.
   Certain properties, like ndof and num_particles are always known and there is no need for them to
   be accessible via the GPUArray.

    Computed quantities available in the GPUArray:
     - ergaenzen

    Values available all the time
     - number of particles in the group

    \ingroup computes
*/

class PYBIND11_EXPORT ComputeSPFBasicProperties : public Compute
    {
    public:
    //! Constructs the compute
    ComputeSPFBasicProperties(std::shared_ptr<SystemDefinition> sysdef, std::shared_ptr<ParticleGroup> group);

    //! Destructor
    virtual ~ComputeSPFBasicProperties();

    //! Compute the temperature
    virtual void compute(uint64_t timestep);

    //! Returns the total kinetic energy last computed by compute()
    /*! \returns Instantaneous total kinetic energy of the system
     */
    Scalar getAbsoluteVelocity()
        {
#ifdef ENABLE_MPI
        if (!m_properties_reduced)
            reduceProperties();
#endif
        const unsigned int num_global = m_group->getNumMembersGlobal();
            { // GPU Array Scope
            // return only translational component if the flags are not valid
            ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
            return h_properties.data[singlephaseflow_logger_index::abs_velocity]/num_global;
            } // End GPU Array Scope
        }

    //! Returns the total kinetic energy last computed by compute()
    /*! \returns Instantaneous total kinetic energy of the system
     */
    Scalar getEkinFluid()
        {
#ifdef ENABLE_MPI
        if (!m_properties_reduced)
            reduceProperties();
#endif

        // return only translational component if the flags are not valid
        ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
        return h_properties.data[singlephaseflow_logger_index::e_kin_fluid];
        }

    //! Returns the sum of particle fluid velocity in xdir last computed by compute()
    /*! \returns Instantaneous sum of particle fluid velocity in xdir 
     */

    Scalar getSumFluidXVelocity()
        {
#ifdef ENABLE_MPI
        if (!m_properties_reduced)
            reduceProperties();
#endif

        // return only translational component if the flags are not valid
        ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
        return h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_x];
        }

    //! Returns the sum of particle fluid velocity in ydir last computed by compute()
    /*! \returns Instantaneous sum of particle fluid velocity in ydir 
     */

    Scalar getSumFluidYVelocity()
        {
#ifdef ENABLE_MPI
        if (!m_properties_reduced)
            reduceProperties();
#endif

        // return only translational component if the flags are not valid
        ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
        return h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_y];
        }

    //! Returns the sum of particle fluid velocity in zdir last computed by compute()
    /*! \returns Instantaneous sum of particle fluid velocity in zdir 
     */

    Scalar getSumFluidZVelocity()
        {
#ifdef ENABLE_MPI
        if (!m_properties_reduced)
            reduceProperties();
#endif

        // return only translational component if the flags are not valid
        ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
        return h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_z];
        }


    //! Returns the mean fluid partcile density last computed by compute()
    /*! \returns Instantaneous mean fluid particle density 
     */

    Scalar getMeanFluidDensity()
        {
#ifdef ENABLE_MPI
        if (!m_properties_reduced)
            reduceProperties();
#endif
        const unsigned int num_global = m_group->getNumMembersGlobal();
            { // GPU Array Scope 
            // return only translational component if the flags are not valid
            ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
            return h_properties.data[singlephaseflow_logger_index::sum_fluid_density]/num_global;
            } // End GPU Array Scope
        }

    unsigned int getNumParticles()
        {
        return m_group->getNumMembersGlobal();
        }

    //! Get the gpu array of properties
    const GPUArray<Scalar>& getProperties()
        {
#ifdef ENABLE_MPI
        if (!m_properties_reduced)
            reduceProperties();
#endif

        return m_properties;
        }

    /// Get the box volume (or area in 2D)
    const Scalar getVolume()
        {
        bool two_d = m_sysdef->getNDimensions() == 2;
        return m_sysdef->getParticleData()->getGlobalBox().getVolume(two_d);
        }

    protected:
    std::shared_ptr<ParticleGroup> m_group; //!< Group to compute properties for
    GPUArray<Scalar> m_properties;       //!< Stores the computed properties

    /// Store the particle data flags used during the last computation
    PDataFlags m_computed_flags;

    //! Does the actual computation
    virtual void computeProperties();

#ifdef ENABLE_MPI
    bool m_properties_reduced; //!< True if properties have been reduced across MPI

    //! Reduce properties over MPI
    virtual void reduceProperties();
#endif
    };

    } // end namespace sph
    } // end namespace hoomd

#endif
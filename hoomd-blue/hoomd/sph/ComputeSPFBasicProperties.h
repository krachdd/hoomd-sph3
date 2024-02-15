/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "ComputeMechanicalPropTypes.h"
#include "hoomd/Compute.h"
#include "hoomd/GlobalArray.h"
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
    All computed values are stored in a GlobalArray so that they can be accessed on the GPU without
   intermediate copies. Use the enum values in singlephaseflow_logger_index to index the array and extract the
   properties of interest. Convenience functions are provided for accessing the values on the CPU.
   Certain properties, like ndof and num_particles are always known and there is no need for them to
   be accessible via the GlobalArray.

    Computed quantities available in the GlobalArray:
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

        // return only translational component if the flags are not valid
        ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
        return h_properties.data[singlephaseflow_logger_index::abs_velocity]/m_group->getNumMembersGlobal();
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

        // return only translational component if the flags are not valid
        ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
        return h_properties.data[singlephaseflow_logger_index::sum_fluid_density]/m_group->getNumMembersGlobal();
        }

    unsigned int getNumParticles()
        {
        return m_group->getNumMembersGlobal();
        }

    //! Get the gpu array of properties
    const GlobalArray<Scalar>& getProperties()
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
    GlobalArray<Scalar> m_properties;       //!< Stores the computed properties

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
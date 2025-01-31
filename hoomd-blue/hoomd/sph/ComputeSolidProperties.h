/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "ComputeMechanicalPropTypes.h"
#include "hoomd/Compute.h"
#include "hoomd/ParticleGroup.h"

#include <limits>
#include <memory>

/*! \file ComputeSolidProperties.h
    \brief Declares a class for computing quantities
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>

#ifndef __COMPUTE_SOLIDPHASE_PROPERTIES_H__
#define __COMPUTE_SOLIDPHASE_PROPERTIES_H__

namespace hoomd
    {
namespace sph
    {

//! Computes properties of a group of particles
/*! ComputeSPFMechanical properties calculates instantaneous properties and provides them in Python.
    All computed values are stored in a GPUArray so that they can be accessed on the GPU without
   intermediate copies. Use the enum values in solidphase_logger_index to index the array and extract the
   properties of interest. Convenience functions are provided for accessing the values on the CPU.
   Certain properties, like ndof and num_particles are always known and there is no need for them to
   be accessible via the GPUArray.

    Computed quantities available in the GPUArray:
     - ergaenzen

    Values available all the time
     - number of particles in the group

    \ingroup computes
*/

class PYBIND11_EXPORT ComputeSolidProperties : public Compute
    {
    public:
    //! Constructs the compute
    ComputeSolidProperties(std::shared_ptr<SystemDefinition> sysdef, std::shared_ptr<ParticleGroup> group);

    //! Destructor
    virtual ~ComputeSolidProperties();

    //! Compute the temperature
    virtual void compute(uint64_t timestep);


    //! Returns the the total drag force on solid phase last computed by compute()
    /*! \returns Instantaneous total drag force on the solid phase 
     */
    Scalar getTotalDragX()
        {
#ifdef ENABLE_MPI
        if (!m_properties_reduced)
            reduceProperties();
#endif

        ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
        return h_properties.data[solidphase_logger_index::total_drag_x];
        }

    //! Returns the the total drag force on solid phase last computed by compute()
    /*! \returns Instantaneous total drag force on the solid phase 
     */
    Scalar getTotalDragY()
        {
#ifdef ENABLE_MPI
        if (!m_properties_reduced)
            reduceProperties();
#endif

        ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
        return h_properties.data[solidphase_logger_index::total_drag_y];
        }

    //! Returns the the total drag force on solid phase last computed by compute()
    /*! \returns Instantaneous total drag force on the solid phase 
     */
    Scalar getTotalDragZ()
        {
#ifdef ENABLE_MPI
        if (!m_properties_reduced)
            reduceProperties();
#endif

        ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
        return h_properties.data[solidphase_logger_index::total_drag_z];
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
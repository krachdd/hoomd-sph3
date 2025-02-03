/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*! \file ComputeSolidProperties.cc
    \brief Contains code for the ComputeSolidProperties class
*/

#include "ComputeSolidProperties.h"
#include "hoomd/VectorMath.h"

#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#include "hoomd/HOOMDMPI.h"
#endif

#include <iostream>
using namespace std;

namespace hoomd
    {
namespace sph
    {
/*! \param sysdef System for which to compute  properties
    \param group Subset of the system over which properties are calculated
*/
ComputeSolidProperties::ComputeSolidProperties(std::shared_ptr<SystemDefinition> sysdef,
                             std::shared_ptr<ParticleGroup> group)
    : Compute(sysdef), m_group(group)
    {
    m_exec_conf->msg->notice(5) << "Constructing ComputeSolidProperties" << endl;

    assert(m_pdata);
    GPUArray<Scalar> properties(solidphase_logger_index::num_quantities, m_exec_conf);
    m_properties.swap(properties);


    m_computed_flags.reset();

#ifdef ENABLE_MPI
    m_properties_reduced = true;
#endif
    }

ComputeSolidProperties::~ComputeSolidProperties()
    {
    m_exec_conf->msg->notice(5) << "Destroying ComputeSolidProperties" << endl;
    }

/*! Calls computeProperties if the properties need updating
    \param timestep Current time step of the simulation
*/
void ComputeSolidProperties::compute(uint64_t timestep)
    {
    Compute::compute(timestep);
    if (shouldCompute(timestep))
        {
        computeProperties();
        m_computed_flags = m_pdata->getFlags();
        }
    }


/*! Computes all properties of the system in one fell swoop.
 */
void ComputeSolidProperties::computeProperties()
    {
    // just drop out if the group is an empty group
    if (m_group->getNumMembersGlobal() == 0)
        return;

    unsigned int group_size = m_group->getNumMembers();

    assert(m_pdata);

    // PDataFlags flags = m_pdata->getFlags();

        { // GPU Array scope
        // access the particle data
        ArrayHandle<Scalar4> h_net_force(m_pdata->getNetForce(), access_location::host, access_mode::read); // that is the net force of ts -1 

        double total_drag_x = 0.0;
        double total_drag_y = 0.0;
        double total_drag_z = 0.0;
        
        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            // Read particle index
            unsigned int j = m_group->getMemberIndex(group_idx);

            // Sum up drag forces
            total_drag_x += h_net_force.data[j].x;
            total_drag_y += h_net_force.data[j].y;
            total_drag_z += h_net_force.data[j].z;
            }

        ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::overwrite);
        h_properties.data[solidphase_logger_index::total_drag_x]  = Scalar(total_drag_x);
        h_properties.data[solidphase_logger_index::total_drag_y]  = Scalar(total_drag_y);
        h_properties.data[solidphase_logger_index::total_drag_z]  = Scalar(total_drag_z);
        } // end GPU Array scope
#ifdef ENABLE_MPI
    // in MPI, reduce extensive quantities only when they're needed
    m_properties_reduced = !m_pdata->getDomainDecomposition();
#endif // ENABLE_MPI
}


#ifdef ENABLE_MPI
void ComputeSolidProperties::reduceProperties()
    {
    if (m_properties_reduced)
        {
        return;
        }

        { // GPU Array Scope
        // reduce properties
        ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::readwrite);
        MPI_Allreduce(MPI_IN_PLACE,
                      h_properties.data,
                      solidphase_logger_index::num_quantities,
                      MPI_HOOMD_SCALAR,
                      MPI_SUM,
                      m_exec_conf->getMPICommunicator());

        m_properties_reduced = true;
        } // end GPU Array Scope
    }
#endif

namespace detail
    {
void export_ComputeSolidProperties(pybind11::module& m)
    {
    pybind11::class_<ComputeSolidProperties, Compute, std::shared_ptr<ComputeSolidProperties>>(m, "ComputeSolidProperties")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<ParticleGroup>>())
        .def_property_readonly("num_particles", &ComputeSolidProperties::getNumParticles)
        .def_property_readonly("total_drag_x", &ComputeSolidProperties::getTotalDragX)
        .def_property_readonly("total_drag_y", &ComputeSolidProperties::getTotalDragY)
        .def_property_readonly("total_drag_z", &ComputeSolidProperties::getTotalDragZ)
        .def_property_readonly("volume", &ComputeSolidProperties::getVolume);
    }

    } // end namespace detail
    } // end namespace sph
    } // end namespace hoomd
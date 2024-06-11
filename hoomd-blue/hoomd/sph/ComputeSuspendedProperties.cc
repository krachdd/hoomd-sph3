/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*! \file ComputeSuspendedProperties.cc
    \brief Contains code for the ComputeSuspendedProperties class
*/

#include "ComputeSuspendedProperties.h"
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
ComputeSuspendedProperties::ComputeSuspendedProperties(std::shared_ptr<SystemDefinition> sysdef,
                             std::shared_ptr<ParticleGroup> group)
    : Compute(sysdef), m_group(group)
    {
    m_exec_conf->msg->notice(5) << "Constructing ComputeSuspendedProperties" << endl;

    assert(m_pdata);
    GlobalArray<Scalar> properties(solidphase_logger_index::num_quantities, m_exec_conf);
    m_properties.swap(properties);
    TAG_ALLOCATION(m_properties);

#if defined(ENABLE_HIP) && defined(__HIP_PLATFORM_NVCC__)
    if (m_exec_conf->isCUDAEnabled() && m_exec_conf->allConcurrentManagedAccess())
        {
        // store in host memory for faster access from CPU
        cudaMemAdvise(m_properties.get(),
                      m_properties.getNumElements() * sizeof(Scalar),
                      cudaMemAdviseSetPreferredLocation,
                      cudaCpuDeviceId);
        CHECK_CUDA_ERROR();
        }
#endif

    m_computed_flags.reset();

#ifdef ENABLE_MPI
    m_properties_reduced = true;
#endif
    }

ComputeSuspendedProperties::~ComputeSuspendedProperties()
    {
    m_exec_conf->msg->notice(5) << "Destroying ComputeSuspendedProperties" << endl;
    }

/*! Calls computeProperties if the properties need updating
    \param timestep Current time step of the simulation
*/
void ComputeSuspendedProperties::compute(uint64_t timestep)
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
void ComputeSuspendedProperties::computeProperties()
    {
    // just drop out if the group is an empty group
    if (m_group->getNumMembersGlobal() == 0)
        return;

    unsigned int group_size = m_group->getNumMembers();

    assert(m_pdata);

    // PDataFlags flags = m_pdata->getFlags();

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

#ifdef ENABLE_MPI
    // in MPI, reduce extensive quantities only when they're needed
    m_properties_reduced = !m_pdata->getDomainDecomposition();
#endif // ENABLE_MPI
}


#ifdef ENABLE_MPI
void ComputeSuspendedProperties::reduceProperties()
    {
    if (m_properties_reduced)
        return;

    // reduce properties
    ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::readwrite);
    MPI_Allreduce(MPI_IN_PLACE,
                  h_properties.data,
                  solidphase_logger_index::num_quantities,
                  MPI_HOOMD_SCALAR,
                  MPI_SUM,
                  m_exec_conf->getMPICommunicator());

    m_properties_reduced = true;
    }
#endif

namespace detail
    {
void export_ComputeSuspendedProperties(pybind11::module& m)
    {
    pybind11::class_<ComputeSuspendedProperties, Compute, std::shared_ptr<ComputeSuspendedProperties>>(m, "ComputeSuspendedProperties")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<ParticleGroup>>())
        .def_property_readonly("num_particles", &ComputeSuspendedProperties::getNumParticles)
        .def_property_readonly("total_drag_x", &ComputeSuspendedProperties::getTotalDragX)
        .def_property_readonly("total_drag_y", &ComputeSuspendedProperties::getTotalDragY)
        .def_property_readonly("total_drag_z", &ComputeSuspendedProperties::getTotalDragZ)
        .def_property_readonly("volume", &ComputeSuspendedProperties::getVolume);
    }

    } // end namespace detail
    } // end namespace md
    } // end namespace hoomd
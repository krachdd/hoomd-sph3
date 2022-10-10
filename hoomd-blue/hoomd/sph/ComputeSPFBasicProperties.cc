/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*! \file ComputeSPFBasicProperties.cc
    \brief Contains code for the ComputeSPFBasicProperties class
*/

#include "ComputeSPFBasicProperties.h"
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
ComputeSPFBasicProperties::ComputeSPFBasicProperties(std::shared_ptr<SystemDefinition> sysdef,
                             std::shared_ptr<ParticleGroup> group)
    : Compute(sysdef), m_group(group)
    {
    m_exec_conf->msg->notice(5) << "Constructing ComputeSPFBasicProperties" << endl;

    assert(m_pdata);
    GlobalArray<Scalar> properties(singlephaseflow_logger_index::num_quantities, m_exec_conf);
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

ComputeSPFBasicProperties::~ComputeSPFBasicProperties()
    {
    m_exec_conf->msg->notice(5) << "Destroying ComputeSPFBasicProperties" << endl;
    }

/*! Calls computeProperties if the properties need updating
    \param timestep Current time step of the simulation
*/
void ComputeSPFBasicProperties::compute(uint64_t timestep)
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
void ComputeSPFBasicProperties::computeProperties()
    {
    // just drop out if the group is an empty group
    if (m_group->getNumMembersGlobal() == 0)
        return;

    unsigned int group_size = m_group->getNumMembers();

    assert(m_pdata);

    // PDataFlags flags = m_pdata->getFlags();

    // access the particle data
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar> h_density(m_pdata->getDensities(), access_location::host, access_mode::read);

    // ArrayHandle<unsigned int> h_body(m_pdata->getBodies(),
    //                                  access_location::host,
    //                                  access_mode::read);
    // ArrayHandle<unsigned int> h_tag(m_pdata->getTags(), access_location::host, access_mode::read);
    // tag evtl needed for group wise logging

    // if (flags[pdata_flag::kinetic_energy]){
    // TO DO, add flags and use them in the future
    // }

    double fluid_vel_x_sum  = 0.0;
    double fluid_vel_y_sum  = 0.0;
    double fluid_vel_z_sum  = 0.0;
    double sum_density      = 0.0;
    // double fluid_prtl = 0;
    double kinetic_energy = 0.0;
    // double adaptive_tstep = 0.0;
    
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int j = m_group->getMemberIndex(group_idx);


        // Sum velocities
        fluid_vel_x_sum += h_vel.data[j].x;
        fluid_vel_y_sum += h_vel.data[j].y;
        fluid_vel_z_sum += h_vel.data[j].z;
        kinetic_energy  += abs(sqrt(pow(h_vel.data[j].x,2)+pow(h_vel.data[j].y,2)+pow(h_vel.data[j].z,2)));
        sum_density     += h_density.data[j];
    }

    ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::overwrite);
    h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_x]  = Scalar(fluid_vel_x_sum);
    h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_y]  = Scalar(fluid_vel_y_sum);
    h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_z]  = Scalar(fluid_vel_z_sum);
    h_properties.data[singlephaseflow_logger_index::sum_fluid_density]    = Scalar(sum_density);
    // h_properties.data[singlephaseflow_logger_index::total_fluid_particles] = Scalar(fluid_prtl);
    h_properties.data[singlephaseflow_logger_index::kinetic_energy]        = Scalar(kinetic_energy);
    // h_properties.data[singlephaseflow_logger_index::dt_adapt] = Scalar(adaptive_tstep);

#ifdef ENABLE_MPI
    // in MPI, reduce extensive quantities only when they're needed
    m_properties_reduced = !m_pdata->getDomainDecomposition();
#endif // ENABLE_MPI
}


#ifdef ENABLE_MPI
void ComputeSPFBasicProperties::reduceProperties()
    {
    if (m_properties_reduced)
        return;

    // reduce properties
    ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::readwrite);
    MPI_Allreduce(MPI_IN_PLACE,
                  h_properties.data,
                  singlephaseflow_logger_index::num_quantities,
                  MPI_HOOMD_SCALAR,
                  MPI_SUM,
                  m_exec_conf->getMPICommunicator());

    m_properties_reduced = true;
    }
#endif

namespace detail
    {
void export_ComputeSPFMechanicalProperties(pybind11::module& m)
    {
    pybind11::class_<ComputeSPFBasicProperties, Compute, std::shared_ptr<ComputeSPFBasicProperties>>(m, "ComputeSPFBasicProperties")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<ParticleGroup>>())
        .def_property_readonly("num_particles", &ComputeSPFBasicProperties::getNumParticles)
        .def_property_readonly("kinetic_energy", &ComputeSPFBasicProperties::getKineticEnergy)
        .def_property_readonly("fluid_vel_x_sum", &ComputeSPFBasicProperties::getSumFluidXVelocity)
        .def_property_readonly("fluid_vel_y_sum", &ComputeSPFBasicProperties::getSumFluidYVelocity)
        .def_property_readonly("fluid_vel_z_sum", &ComputeSPFBasicProperties::getSumFluidZVelocity)
        .def_property_readonly("mean_density", &ComputeSPFBasicProperties::getMeanFluidDensity)
        .def_property_readonly("volume", &ComputeSPFBasicProperties::getVolume);
    }

    } // end namespace detail
    } // end namespace md
    } // end namespace hoomd
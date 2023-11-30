/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*! \file ComputeSPFNNBasicProperties.cc
    \brief Contains code for the ComputeSPFNNBasicProperties class
*/

#include "ComputeSPFNNBasicProperties.h"
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
ComputeSPFNNBasicProperties::ComputeSPFNNBasicProperties(std::shared_ptr<SystemDefinition> sysdef,
                             std::shared_ptr<ParticleGroup> group)
    : Compute(sysdef), m_group(group)
    {
    m_exec_conf->msg->notice(5) << "Constructing ComputeSPFNNBasicProperties" << endl;

    assert(m_pdata);
    GlobalArray<Scalar> properties(singlephaseflownn_logger_index::num_quantities, m_exec_conf);
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

ComputeSPFNNBasicProperties::~ComputeSPFNNBasicProperties()
    {
    m_exec_conf->msg->notice(5) << "Destroying ComputeSPFNNBasicProperties" << endl;
    }

/*! Calls computeProperties if the properties need updating
    \param timestep Current time step of the simulation
*/
void ComputeSPFNNBasicProperties::compute(uint64_t timestep)
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
void ComputeSPFNNBasicProperties::computeProperties()
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
    ArrayHandle<Scalar3> h_nn(m_pdata->getAuxiliaries3(), access_location::host,access_mode::readwrite);

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
    double abs_velocity = 0.0;
    // double adaptive_tstep = 0.0;
    double mean_viscosity = 0.0;
    double max_viscosity = 0.0;
    double max_shearrate = 0.0;
    
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int j = m_group->getMemberIndex(group_idx);


        // Sum velocities
        fluid_vel_x_sum += h_vel.data[j].x;
        fluid_vel_y_sum += h_vel.data[j].y;
        fluid_vel_z_sum += h_vel.data[j].z;
        abs_velocity  += sqrt(h_vel.data[j].x * h_vel.data[j].x + h_vel.data[j].y * h_vel.data[j].y + h_vel.data[j].z * h_vel.data[j].z);
        sum_density     += h_density.data[j];
        mean_viscosity += h_nn.data[j].y;

        if (h_nn.data[j].y > max_viscosity)
        {
            max_viscosity = h_nn.data[j].y;
        }

        if (h_nn.data[j].x > max_shearrate)
        {
            max_shearrate = h_nn.data[j].x;
        }

    }

    ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::overwrite);
    h_properties.data[singlephaseflownn_logger_index::sum_fluid_velocity_x]  = Scalar(fluid_vel_x_sum);
    h_properties.data[singlephaseflownn_logger_index::sum_fluid_velocity_y]  = Scalar(fluid_vel_y_sum);
    h_properties.data[singlephaseflownn_logger_index::sum_fluid_velocity_z]  = Scalar(fluid_vel_z_sum);
    h_properties.data[singlephaseflownn_logger_index::sum_fluid_density]    = Scalar(sum_density);
    // h_properties.data[singlephaseflownn_logger_index::total_fluid_particles] = Scalar(fluid_prtl);
    h_properties.data[singlephaseflownn_logger_index::abs_velocity]        = Scalar(abs_velocity);
    // h_properties.data[singlephaseflownn_logger_index::dt_adapt] = Scalar(adaptive_tstep);
    h_properties.data[singlephaseflownn_logger_index::mean_viscosity]        = Scalar(mean_viscosity);
    h_properties.data[singlephaseflownn_logger_index::max_viscosity]        = Scalar(max_viscosity);
    h_properties.data[singlephaseflownn_logger_index::max_shearrate]        = Scalar(max_shearrate);



#ifdef ENABLE_MPI
    // in MPI, reduce extensive quantities only when they're needed
    m_properties_reduced = !m_pdata->getDomainDecomposition();
#endif // ENABLE_MPI
}


#ifdef ENABLE_MPI
void ComputeSPFNNBasicProperties::reduceProperties()
    {
    if (m_properties_reduced)
        return;

    // reduce properties
    ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::readwrite);
    MPI_Allreduce(MPI_IN_PLACE,
                  h_properties.data,
                  singlephaseflownn_logger_index::num_quantities,
                  MPI_HOOMD_SCALAR,
                  MPI_SUM,
                  m_exec_conf->getMPICommunicator());

    m_properties_reduced = true;
    }
#endif

namespace detail
    {
void export_ComputeSPFNNMechanicalProperties(pybind11::module& m)
    {
    pybind11::class_<ComputeSPFNNBasicProperties, Compute, std::shared_ptr<ComputeSPFNNBasicProperties>>(m, "ComputeSPFNNBasicProperties")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<ParticleGroup>>())
        .def_property_readonly("num_particles", &ComputeSPFNNBasicProperties::getNumParticles)
        .def_property_readonly("abs_velocity", &ComputeSPFNNBasicProperties::getAbsoluteVelocity)
        .def_property_readonly("fluid_vel_x_sum", &ComputeSPFNNBasicProperties::getSumFluidXVelocity)
        .def_property_readonly("fluid_vel_y_sum", &ComputeSPFNNBasicProperties::getSumFluidYVelocity)
        .def_property_readonly("fluid_vel_z_sum", &ComputeSPFNNBasicProperties::getSumFluidZVelocity)
        .def_property_readonly("mean_density", &ComputeSPFNNBasicProperties::getMeanFluidDensity)
        .def_property_readonly("mean_viscosity", &ComputeSPFNNBasicProperties::getMeanViscosity)
        .def_property_readonly("max_viscosity", &ComputeSPFNNBasicProperties::getMaxViscosity)
        .def_property_readonly("max_shearrate", &ComputeSPFNNBasicProperties::getMaxShearrate)
        .def_property_readonly("volume", &ComputeSPFNNBasicProperties::getVolume);

    }

    } // end namespace detail
    } // end namespace md
    } // end namespace hoomd
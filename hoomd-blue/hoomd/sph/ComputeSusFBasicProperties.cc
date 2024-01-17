/* ---------------------------------------------------------
maintainer: drostan, daniel.rostan@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*! \file ComputeSusFBasicProperties.cc
    \brief Contains code for the ComputeSusFBasicProperties class
*/

#include "ComputeSusFBasicProperties.h"
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
ComputeSusFBasicProperties::ComputeSusFBasicProperties(std::shared_ptr<SystemDefinition> sysdef,
                             std::shared_ptr<ParticleGroup> group)
    : Compute(sysdef), m_group(group)
    {
    m_exec_conf->msg->notice(5) << "Constructing ComputeSusFBasicProperties" << endl;

    assert(m_pdata);
    GlobalArray<Scalar> properties(suspensionflow_logger_index::num_quantities, m_exec_conf);
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

ComputeSusFBasicProperties::~ComputeSusFBasicProperties()
    {
    m_exec_conf->msg->notice(5) << "Destroying ComputeSusFBasicProperties" << endl;
    }

/*! Calls computeProperties if the properties need updating
    \param timestep Current time step of the simulation
*/
void ComputeSusFBasicProperties::compute(uint64_t timestep)
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
void ComputeSusFBasicProperties::computeProperties()
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
    ArrayHandle<Scalar3> h_av(m_pdata->getAuxiliaries3(), access_location::host,access_mode::read);
    ArrayHandle<Scalar3> h_tv(m_pdata->getAuxiliaries4(), access_location::host,access_mode::read);

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
    double angularvel_x = 0.0;
    double angularvel_y = 0.0;
    double angularvel_z = 0.0;
    double translationvel_x = 0.0;
    double translationvel_y = 0.0;
    double translationvel_z = 0.0;

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
        angularvel_x = h_av.data[j].x;
        angularvel_y = h_av.data[j].y;
        angularvel_z = h_av.data[j].z;
        translationvel_x = h_tv.data[j].x;
        translationvel_y = h_tv.data[j].y;
        translationvel_z = h_tv.data[j].z;

    }

    ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::overwrite);
    h_properties.data[suspensionflow_logger_index::sum_fluid_velocity_x]  = Scalar(fluid_vel_x_sum);
    h_properties.data[suspensionflow_logger_index::sum_fluid_velocity_y]  = Scalar(fluid_vel_y_sum);
    h_properties.data[suspensionflow_logger_index::sum_fluid_velocity_z]  = Scalar(fluid_vel_z_sum);
    h_properties.data[suspensionflow_logger_index::sum_fluid_density]    = Scalar(sum_density);
    // h_properties.data[suspensionflow_logger_index::total_fluid_particles] = Scalar(fluid_prtl);
    h_properties.data[suspensionflow_logger_index::angularvel_x]        = Scalar(angularvel_x);
    h_properties.data[suspensionflow_logger_index::angularvel_y]        = Scalar(angularvel_y);
    h_properties.data[suspensionflow_logger_index::angularvel_z]        = Scalar(angularvel_z);
    h_properties.data[suspensionflow_logger_index::translationvel_x]    = Scalar(translationvel_x);
    h_properties.data[suspensionflow_logger_index::translationvel_y]    = Scalar(translationvel_y);
    h_properties.data[suspensionflow_logger_index::translationvel_z]    = Scalar(translationvel_z);


    // h_properties.data[suspensionflow_logger_index::dt_adapt] = Scalar(adaptive_tstep);

#ifdef ENABLE_MPI
    // in MPI, reduce extensive quantities only when they're needed
    m_properties_reduced = !m_pdata->getDomainDecomposition();
#endif // ENABLE_MPI
}


#ifdef ENABLE_MPI
void ComputeSusFBasicProperties::reduceProperties()
    {
    if (m_properties_reduced)
        return;

    // reduce properties
    ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::readwrite);
    MPI_Allreduce(MPI_IN_PLACE,
                  h_properties.data,
                  suspensionflow_logger_index::num_quantities,
                  MPI_HOOMD_SCALAR,
                  MPI_SUM,
                  m_exec_conf->getMPICommunicator());

    m_properties_reduced = true;
    }
#endif

namespace detail
    {
void export_ComputeSusFMechanicalProperties(pybind11::module& m)
    {
    pybind11::class_<ComputeSusFBasicProperties, Compute, std::shared_ptr<ComputeSusFBasicProperties>>(m, "ComputeSusFBasicProperties")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<ParticleGroup>>())
        .def_property_readonly("num_particles", &ComputeSusFBasicProperties::getNumParticles)
        .def_property_readonly("abs_velocity", &ComputeSusFBasicProperties::getAbsoluteVelocity)
        .def_property_readonly("fluid_vel_x_sum", &ComputeSusFBasicProperties::getSumFluidXVelocity)
        .def_property_readonly("fluid_vel_y_sum", &ComputeSusFBasicProperties::getSumFluidYVelocity)
        .def_property_readonly("fluid_vel_z_sum", &ComputeSusFBasicProperties::getSumFluidZVelocity)
        .def_property_readonly("mean_density", &ComputeSusFBasicProperties::getMeanFluidDensity)
        .def_property_readonly("volume", &ComputeSusFBasicProperties::getVolume)
        .def_property_readonly("angularvel_x", &ComputeSusFBasicProperties::getAngularXVelocity)
        .def_property_readonly("angularvel_y", &ComputeSusFBasicProperties::getAngularYVelocity)
        .def_property_readonly("angularvel_z", &ComputeSusFBasicProperties::getAngularZVelocity) 
        .def_property_readonly("translationvel_x", &ComputeSusFBasicProperties::getTranslationalXVelocity)
        .def_property_readonly("translationvel_y", &ComputeSusFBasicProperties::getTranslationalYVelocity)
        .def_property_readonly("translationvel_z", &ComputeSusFBasicProperties::getTranslationalZVelocity);      
    }

    } // end namespace detail
    } // end namespace md
    } // end namespace hoomd
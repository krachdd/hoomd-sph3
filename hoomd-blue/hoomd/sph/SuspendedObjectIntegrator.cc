/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SuspendedObjectIntegrator.h"

#include "hoomd/VectorMath.h"
#include <vector>


#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#endif

using namespace std;

/*! \file SuspendedObjectIntegrator.h
    \brief Contains code for the SuspendedObjectIntegrator class
*/

namespace hoomd
    {
namespace sph
    {
/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
    \param group The group of particles this integration method is to work on
*/
SuspendedObjectIntegrator::SuspendedObjectIntegrator(std::shared_ptr<SystemDefinition> sysdef,
                       std::shared_ptr<ParticleGroup> group)
    : SPHIntegrationMethodTwoStep(sysdef, group), m_limit(false), m_limit_val(1.0), m_zero_force(false)
    {
    m_exec_conf->msg->notice(5) << "Constructing SuspendedObjectIntegrator" << endl;
    m_densitymethod_set = false;

    // Compute total mass of suspended object
    ComputeTotalMass(true);

    // Compute initial center of mass
    ComputeCenterOfMass(true);

    // Compute initial translation velocity
    ComputeTranslationVelocity(true);

    // Compute moment of inertia tensor
    ComputeMomentOfInertia(true);

    // Compute angular velocity
    ComputeAngularVelocity(true);

    // Initialize translational and angular accelerations to zero
    m_translationaccel = make_scalar3(0.0, 0.0, 0.0);
    m_angularaccel = make_scalar3(0.0, 0.0, 0.0);  
    }

SuspendedObjectIntegrator::~SuspendedObjectIntegrator()
    {
    m_exec_conf->msg->notice(5) << "Destroying SuspendedObjectIntegrator" << endl;
    }

/*! \param limit Distance to limit particle movement each time step

    Once the limit is set, future calls to update() will never move a particle
    a distance larger than the limit in a single time step
*/

pybind11::object SuspendedObjectIntegrator::getLimit()
    {
    pybind11::object result;
    if (m_limit)
        {
        result = pybind11::cast(m_limit_val);
        }
    else
        {
        result = pybind11::none();
        }
    return result;
    }

void SuspendedObjectIntegrator::setLimit(pybind11::object limit)
    {
    if (limit.is_none())
        {
        m_limit = false;
        }
    else
        {
        m_limit = true;
        m_limit_val = pybind11::cast<Scalar>(limit);
        }
    }

bool SuspendedObjectIntegrator::getZeroForce()
    {
    return m_zero_force;
    }

void SuspendedObjectIntegrator::setZeroForce(bool zero_force)
    {
    m_zero_force = zero_force;
    }

void SuspendedObjectIntegrator::ComputeTotalMass(bool print = false)
    {
    // Initialize mass to zero
    m_totalmass = 0;

    // Read-Handle to velocity mass array
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::read);

    // Loop over group index
    unsigned int group_size = m_group->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);
        m_totalmass += h_vel.data[j].w;
        }

#ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(MPI_IN_PLACE, &m_totalmass, 1, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
#endif

    if ( print )
        m_exec_conf->msg->notice(5) << "SuspendedObject Total Mass: " << m_totalmass << endl;
    }

void SuspendedObjectIntegrator::ComputeCenterOfMass(bool print = false)
    {
    // Read-Handle to velocity mass and position array
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);

    // grab the box dimensions
    BoxDim globalbox = m_pdata->getGlobalBox();

    // Limit values
    Scalar3 lo = globalbox.getLo();
    Scalar3 Ld = globalbox.getL();

    // Average position on unit circle
    Scalar angles[6] = {0,0,0,0,0,0};

    // Loop over group index
    unsigned int group_size  = m_group->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);
        Scalar3 pos = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);

        // Compute mapped angles
        Scalar theta_x = ( ( pos.x - lo.x ) / Ld.x ) * Scalar(M_TWOPI);
        Scalar theta_y = ( ( pos.y - lo.y ) / Ld.y ) * Scalar(M_TWOPI);
        Scalar theta_z = ( ( pos.z - lo.z ) / Ld.z ) * Scalar(M_TWOPI);

        // Add contribution to mass averaged points on unit circle
        angles[0] += ( Ld.x / Scalar(M_TWOPI) ) * cos(theta_x);
        angles[1] += ( Ld.x / Scalar(M_TWOPI) ) * sin(theta_x);
        angles[2] += ( Ld.y / Scalar(M_TWOPI) ) * cos(theta_y);
        angles[3] += ( Ld.y / Scalar(M_TWOPI) ) * sin(theta_y);
        angles[4] += ( Ld.z / Scalar(M_TWOPI) ) * cos(theta_z);
        angles[5] += ( Ld.z / Scalar(M_TWOPI) ) * sin(theta_z);
        }

    // Count total number of particles in group
    unsigned int totalgroupN = group_size;

#ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(MPI_IN_PLACE, &angles[0], 6, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
    MPI_Allreduce(MPI_IN_PLACE, &totalgroupN, 1, MPI_UNSIGNED, MPI_SUM, this->m_exec_conf->getMPICommunicator());
#endif

    // Averaging
    angles[0] /= totalgroupN;
    angles[1] /= totalgroupN;
    angles[2] /= totalgroupN;
    angles[3] /= totalgroupN;
    angles[4] /= totalgroupN;
    angles[5] /= totalgroupN;

    // Set attribute
    m_centerofmass = make_scalar3(0.0, 0.0, 0.0);
    if ( angles[1] != 0 && angles[0] != 0 )
        m_centerofmass.x = lo.x + ( atan2(-angles[1], -angles[0]) + Scalar(M_PI) ) * ( Ld.x / Scalar(M_TWOPI) );
    if ( angles[3] != 0 && angles[2] != 0 )
        m_centerofmass.y = lo.y + ( atan2(-angles[3], -angles[2]) + Scalar(M_PI) ) * ( Ld.y / Scalar(M_TWOPI) );
    if ( angles[5] != 0 && angles[4] != 0 )
        m_centerofmass.z = lo.z + ( atan2(-angles[5], -angles[4]) + Scalar(M_PI) ) * ( Ld.z / Scalar(M_TWOPI) );

    if ( print )
        {
        m_exec_conf->msg->notice(5) << "SuspendedObject consists of " << totalgroupN << " slave particles" << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject center of Mass x: " << m_centerofmass.x << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject center of mass y: " << m_centerofmass.y << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject center of mass z: " << m_centerofmass.z << endl;
        }


    // std::cout << "m_centerofmass.x: " << m_centerofmass.x << std::endl;
    // std::cout << "m_centerofmass.y: " << m_centerofmass.y << std::endl;
    // std::cout << "m_centerofmass.z: " << m_centerofmass.z << std::endl;

    }

void SuspendedObjectIntegrator::ComputeTranslationVelocity(bool print = false)
    {
    // Initialize center of mass to zero
    Scalar translationvel[3] = {0,0,0};

    // Read-Handle to velocity mass array
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::read);

    // Loop over group index
    unsigned int group_size = m_group->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);
        Scalar3 vel = make_scalar3(h_vel.data[j].x, h_vel.data[j].y, h_vel.data[j].z);
        Scalar m_M  = h_vel.data[j].w/m_totalmass;
        translationvel[0] += m_M*vel.x;
        translationvel[1] += m_M*vel.y;
        translationvel[2] += m_M*vel.z;
        }

#ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(MPI_IN_PLACE, &translationvel[0], 3, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
#endif

    // Set attribute
    m_translationvel.x = translationvel[0];
    m_translationvel.y = translationvel[1];
    m_translationvel.z = translationvel[2];

    if ( print )
        {
        m_exec_conf->msg->notice(5) << "SuspendedObject translation velocity x: " << m_translationvel.x << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject translation velocity y: " << m_translationvel.y << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject translation velocity z: " << m_translationvel.z << endl;
        }
    }

void SuspendedObjectIntegrator::ComputeMomentOfInertia(bool print = false)
    {
    // Local copy of the simulation box
    const BoxDim& box = m_pdata->getGlobalBox();

    // Initialize moment of inertia tensor to zero
    Scalar J[9] = {0,0,0,0,0,0,0,0,0};

    // Read-Handle to velocity mass and position array
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::read);

    // Loop over group index
    unsigned int group_size = m_group->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // Read particle position and mass
        Scalar3 pos = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);
        Scalar mass = h_vel.data[j].w;

        // Relative position to center of mass
        Scalar3 rdiff = pos - m_centerofmass;

        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);

        Scalar dotrdiff = dot(rdiff,rdiff);

        // Add contribution to moment of inertia tensor
        J[0] += mass*(dotrdiff-rdiff.x*rdiff.x); // 00
        J[1] += mass*(        -rdiff.x*rdiff.y); // 01
        J[2] += mass*(        -rdiff.x*rdiff.z); // 02
        J[3] += mass*(        -rdiff.y*rdiff.x); // 10
        J[4] += mass*(dotrdiff-rdiff.y*rdiff.y); // 11
        J[5] += mass*(        -rdiff.y*rdiff.z); // 12
        J[6] += mass*(        -rdiff.z*rdiff.x); // 20
        J[7] += mass*(        -rdiff.z*rdiff.y); // 21
        J[8] += mass*(dotrdiff-rdiff.z*rdiff.z); // 22
        }

#ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(MPI_IN_PLACE, &J[0], 9, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
#endif

    // Compute inverse of moment of inertia tensor and store as class attribute
    Scalar detJ = +J[0]*(J[4]*J[8]-J[7]*J[5])
                  -J[1]*(J[3]*J[8]-J[5]*J[6])
                  +J[2]*(J[3]*J[7]-J[4]*J[6]);
    Scalar invdetJ = 1/detJ;

    m_invJ[0] =   (J[4]*J[8]-J[7]*J[5])*invdetJ;
    m_invJ[1] =  -(J[1]*J[8]-J[2]*J[7])*invdetJ;
    m_invJ[2] =   (J[1]*J[5]-J[2]*J[4])*invdetJ;
    m_invJ[3] =  -(J[3]*J[8]-J[5]*J[6])*invdetJ;
    m_invJ[4] =   (J[0]*J[8]-J[2]*J[6])*invdetJ;
    m_invJ[5] =  -(J[0]*J[5]-J[3]*J[2])*invdetJ;
    m_invJ[6] =   (J[3]*J[7]-J[6]*J[4])*invdetJ;
    m_invJ[7] =  -(J[0]*J[7]-J[6]*J[1])*invdetJ;
    m_invJ[8] =   (J[0]*J[4]-J[3]*J[1])*invdetJ;

    if ( print )
        {
        m_exec_conf->msg->notice(5) << "SuspendedObject moment of inertia: (0,0):" << J[0] << " (0,1):" << J[1] << " (0,2):" << J[2] << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject moment of inertia: (1,0):" << J[3] << " (1,1):" << J[4] << " (1,2):" << J[5] << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject moment of inertia: (2,0):" << J[6] << " (2,1):" << J[7] << " (2,2):" << J[8] << endl;
        }
    }

void SuspendedObjectIntegrator::ComputeAngularVelocity(bool print = false)
    {
    // Local copy of the simulation box
    const BoxDim& box = m_pdata->getGlobalBox();

    // Initialize center of mass to zero
    Scalar angmomentum[3] = {0,0,0};

    // Read-Handle to velocity mass and position array
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::read);

    // Loop over group index
    unsigned int group_size = m_group->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // Read particle position and mass
        Scalar3 pos = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);
        Scalar3 vel = make_scalar3(h_vel.data[j].x, h_vel.data[j].y, h_vel.data[j].z);
        Scalar mass = h_vel.data[j].w;

        // Relative position to center of mass
        Scalar3 rdiff = pos - m_centerofmass;

        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);

        // Relative position to center of mass
        vec3<Scalar> rdiffm(rdiff.x*mass,rdiff.y*mass,rdiff.z*mass);
        vec3<Scalar> vdiff(vel.x-m_translationvel.x,vel.y-m_translationvel.y,vel.z-m_translationvel.z);

        // Angular momentum due to particle j
        Scalar3 angmomentum_j = vec_to_scalar3(cross(rdiffm,vdiff));

        // Compute total angular momentum
        angmomentum[0] += angmomentum_j.x;
        angmomentum[1] += angmomentum_j.y;
        angmomentum[2] += angmomentum_j.z;
        }

#ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(MPI_IN_PLACE, &angmomentum[0], 3, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
#endif

    // Compute angular velocity using inverted moment of inertia tensor
    m_angularvel.x = m_invJ[0]*angmomentum[0]+m_invJ[1]*angmomentum[1]+m_invJ[2]*angmomentum[2];
    m_angularvel.y = m_invJ[3]*angmomentum[0]+m_invJ[4]*angmomentum[1]+m_invJ[5]*angmomentum[2];
    m_angularvel.z = m_invJ[6]*angmomentum[0]+m_invJ[7]*angmomentum[1]+m_invJ[8]*angmomentum[2];

    if ( print )
        {
        m_exec_conf->msg->notice(5) << "SuspendedObject initial angular velocity x: " << m_angularvel.x << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject initial angular velocity y: " << m_angularvel.y << endl;
        m_exec_conf->msg->notice(5) << "SuspendedObject initial angular velocity z: " << m_angularvel.z << endl;
        }
    }

/*! \param timestep Current time step
    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2 per the
   velocity verlet method.
*/
void SuspendedObjectIntegrator::integrateStepOne(uint64_t timestep)
    {

    m_exec_conf->msg->notice(9) << "SuspendedObjectIntegrator: Integrate Step one" << endl;

    // Local copy of the simulation box
    const BoxDim& box = m_pdata->getGlobalBox();

    m_translationvel += Scalar(1.0/2.0)*m_translationaccel*m_deltaT;
    m_angularvel += Scalar(1.0/2.0)*m_angularaccel*m_deltaT;
    vec3<Scalar> angularvel(m_angularvel.x, m_angularvel.y, m_angularvel.z);

    // Array handles
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::overwrite);

    // Loop over group index
    unsigned int group_size = m_group->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // Read particle position and compute relative position
        Scalar3 pos = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);

        // Relative position to center of mass
        Scalar3 rdiff = pos - m_centerofmass;

        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);

        // Relative position to center of mass
        vec3<Scalar> rdiffv(rdiff.x, rdiff.y, rdiff.z);

        // Rotational contribution
        vec3<Scalar> rot = cross(angularvel,rdiffv);

        // Update slave particle velocities
        h_vel.data[j].x = m_translationvel.x + rot.x;
        h_vel.data[j].y = m_translationvel.y + rot.y;
        h_vel.data[j].z = m_translationvel.z + rot.z;

        // Update slave particle positions
        h_pos.data[j].x += h_vel.data[j].x*m_deltaT;
        h_pos.data[j].y += h_vel.data[j].y*m_deltaT;
        h_pos.data[j].z += h_vel.data[j].z*m_deltaT;
        }

#ifdef ENABLE_MPI
    MPI_Barrier(m_exec_conf->getMPICommunicator());
#endif

    // Recompute center of mass
    ComputeCenterOfMass(timestep);

    }

/*! \param timestep Current time step
    \post particle velocities are moved forward to timestep+1
*/
void SuspendedObjectIntegrator::integrateStepTwo(uint64_t timestep)
    {
    m_exec_conf->msg->notice(9) << "SuspendedObjectIntegrator: Integrate Step two" << endl;

    // Local copy of the simulation box
    const BoxDim& box = m_pdata->getGlobalBox();

    // Get array handles
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_net_force(m_pdata->getNetForce(), access_location::host, access_mode::read);

    // Initialize total force and torque vectors
    Scalar totalforcetorque[6] = {0,0,0,0,0,0};

    // Loop over group index
    unsigned int group_size = m_group->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // Net force acting on slave particle
        vec3<Scalar> net_force_j(h_net_force.data[j].x, h_net_force.data[j].y, h_net_force.data[j].z);

        // Read particle position and compute relative position
        Scalar3 pos = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);

        // Relative position to center of mass
        Scalar3 rdiff = pos - m_centerofmass;

        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);

        // Relative position to center of mass
        vec3<Scalar> rdiffv(rdiff.x, rdiff.y, rdiff.z);

        // Torque contribution
        Scalar3 torque_j = vec_to_scalar3(cross(rdiffv,net_force_j));

        // Compute total force
        totalforcetorque[0] += net_force_j.x;
        totalforcetorque[1] += net_force_j.y;
        totalforcetorque[2] += net_force_j.z;

        // Compute total torque
        totalforcetorque[3] += torque_j.x;
        totalforcetorque[4] += torque_j.y;
        totalforcetorque[5] += torque_j.z;
        }

#ifdef ENABLE_MPI
    // Reduce on all processors
    MPI_Allreduce(MPI_IN_PLACE, &totalforcetorque[0], 6, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
#endif

    // Compute translational acceleration
    m_translationaccel.x = totalforcetorque[0]/m_totalmass;
    m_translationaccel.y = totalforcetorque[1]/m_totalmass;
    m_translationaccel.z = totalforcetorque[2]/m_totalmass;

    // Compute angular acceleration
    m_angularaccel.x = m_invJ[0]*totalforcetorque[3]+m_invJ[1]*totalforcetorque[4]+m_invJ[2]*totalforcetorque[5];
    m_angularaccel.y = m_invJ[3]*totalforcetorque[3]+m_invJ[4]*totalforcetorque[4]+m_invJ[5]*totalforcetorque[5];
    m_angularaccel.z = m_invJ[6]*totalforcetorque[3]+m_invJ[7]*totalforcetorque[4]+m_invJ[8]*totalforcetorque[5];

    // Integrate Object variables
    m_translationvel += Scalar(1.0/2.0)*m_translationaccel*m_deltaT;
    m_angularvel += Scalar(1.0/2.0)*m_angularaccel*m_deltaT;
    vec3<Scalar> angularvel(m_angularvel.x, m_angularvel.y, m_angularvel.z);

    // Update slave particle velocities
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::overwrite);

    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // Read particle position and compute relative position
        Scalar3 pos = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);

        // Relative position to center of mass
        Scalar3 rdiff = pos - m_centerofmass;

        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);

        // Relative position to center of mass
        vec3<Scalar> rdiffv(rdiff.x, rdiff.y, rdiff.z);

        // Rotational contribution
        vec3<Scalar> rot = cross(angularvel,rdiffv);

        // Update slave particle velocities
        h_vel.data[j].x = m_translationvel.x + rot.x;
        h_vel.data[j].y = m_translationvel.y + rot.y;
        h_vel.data[j].z = m_translationvel.z + rot.z;
        }
    }

namespace detail
    {
void export_SuspendedObjectIntegrator(pybind11::module& m)
    {
    pybind11::class_<SuspendedObjectIntegrator, SPHIntegrationMethodTwoStep, std::shared_ptr<SuspendedObjectIntegrator>>(
        m,
        "SuspendedObjectIntegrator")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<ParticleGroup>>())
        .def("getDensityMethod", &SuspendedObjectIntegrator::getDensityMethod)
        .def("setDensityMethod", &SuspendedObjectIntegrator::setDensityMethod)
        .def("getCOMX", &SuspendedObjectIntegrator::getCOMX)
        .def("getCOMY", &SuspendedObjectIntegrator::getCOMY)
        .def("getCOMZ", &SuspendedObjectIntegrator::getCOMZ)
        .def("getTranslationX", &SuspendedObjectIntegrator::getTranslationX)
        .def("getTranslationY", &SuspendedObjectIntegrator::getTranslationY)
        .def("getTranslationZ", &SuspendedObjectIntegrator::getTranslationZ)
        .def("getRotationX", &SuspendedObjectIntegrator::getRotationX)
        .def("getRotationY", &SuspendedObjectIntegrator::getRotationY)
        .def("getRotationZ", &SuspendedObjectIntegrator::getRotationZ)
        .def_property("limit", &SuspendedObjectIntegrator::getLimit, &SuspendedObjectIntegrator::setLimit)
        .def_property("zero_force", &SuspendedObjectIntegrator::getZeroForce, &SuspendedObjectIntegrator::setZeroForce);
    }
    } // end namespace detail
    } // end namespace sph
    } // end namespace hoomd


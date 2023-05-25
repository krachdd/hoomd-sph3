/* ---------------------------------------------------------
maintainer: rsivanesapillai, -
----------------------------------------------------------*/

/*
The Velocity Verlet algorithm originally used in Rakulans old code
*/

#include "RigidBodyIntegrator.h"
#include "hoomd/VectorMath.h"
#include <vector>

using namespace std;

/*! \file RigidBodyIntegrator.h
    \brief Contains code for the RigidBodyIntegrator class
*/

namespace hoomd
    {
namespace sph
    {
/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
    \param group The group of particles this integration method is to work on
*/
RigidBodyIntegrator::RigidBodyIntegrator(std::shared_ptr<SystemDefinition> sysdef,
                                         std::shared_ptr<ParticleGroup> group,
                                         std::shared_ptr<Variant> transvel_x,
                                         std::shared_ptr<Variant> transvel_y,
                                         std::shared_ptr<Variant> transvel_z,
                                         std::shared_ptr<Variant> rotatvel,
                                         Scalar pivotpnt_x,
                                         Scalar pivotpnt_y,
                                         Scalar pivotpnt_z,
                                         Scalar rotaxis_x,
                                         Scalar rotaxis_y,
                                         Scalar rotaxis_z)
    : SPHIntegrationMethodTwoStep(sysdef, group), m_transvel_x(transvel_x), m_transvel_y(transvel_y), m_transvel_z(transvel_z), m_rotatvel(rotatvel)
    {
    m_exec_conf->msg->notice(5) << "Constructing RigidBodyIntegrator" << endl;

    // Set Pivot Point
    m_pivotpnt = make_scalar3(pivotpnt_x,pivotpnt_y,pivotpnt_z);

    // Set rotation axis
    vec3<Scalar> rotaxis(rotaxis_x,rotaxis_y,rotaxis_z);
    rotaxis   = rotaxis/Scalar(sqrt(dot(rotaxis,rotaxis)));
    m_rotaxis = vec_to_scalar3(rotaxis);

    }

RigidBodyIntegrator::~RigidBodyIntegrator()
    {
    m_exec_conf->msg->notice(5) << "Destroying RigidBodyIntegrator" << endl;
    }


/*! \param timestep Current time step
    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2.
*/
void RigidBodyIntegrator::integrateStepOne(uint64_t timestep)
    {
    unsigned int group_size = m_group->getNumMembers();
    //std::cout << "group_size_rigid:" << group_size << std::endl;

    m_exec_conf->msg->notice(9) << "RigidBodyIntegrator: Integrate Step one" << endl;

    // Local copy of the simulation box
    const BoxDim& box = m_pdata->getGlobalBox();

    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::overwrite);
    ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::overwrite);
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);

    // Evaluate variant expressions
    Scalar transvel_x    = (*m_transvel_x)(timestep).x;
    Scalar transvel_y    = (*m_transvel_y)(timestep).x;
    Scalar transvel_z    = (*m_transvel_z)(timestep).x;
    Scalar rotatvel      = (*m_rotatvel)(timestep).x;
    Scalar transvel_x_dt = (*m_transvel_x)(timestep).y;
    Scalar transvel_y_dt = (*m_transvel_y)(timestep).y;
    Scalar transvel_z_dt = (*m_transvel_z)(timestep).y;
    Scalar rotatvel_dt   = (*m_rotatvel)(timestep).y;
    vec3<Scalar> angularvel(rotatvel*m_rotaxis.x, rotatvel*m_rotaxis.y, rotatvel*m_rotaxis.z);
    vec3<Scalar> angularaccel(rotatvel_dt*m_rotaxis.x, rotatvel_dt*m_rotaxis.y, rotatvel_dt*m_rotaxis.z);

    // std::cout << "transvel_x:" << transvel_x << std::endl;
    // std::cout << "rotatvel:" << rotatvel << std::endl;
    // std::cout << "transvel_x_dt:" << transvel_x_dt << std::endl;
    // std::cout << "rotatvel_dt:" << rotatvel_dt << std::endl;

    // std::cout << "angularvel_z:" << angularvel.z << std::endl;
    // std::cout << "angularaccel_z:" << angularaccel.z << std::endl;

    // perform the first half step of the integration
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // Read particle position and compute relative position
        Scalar3 pos = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);

        // Relative position to rotational pivot point
        Scalar3 rdiff = pos-m_pivotpnt;

        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);

        // Relative position to center of mass
        vec3<Scalar> rdiffv(rdiff.x, rdiff.y, rdiff.z);

        // Angular velocity
        vec3<Scalar> rot = cross(angularvel,rdiffv);

        // Angular acceleration
        vec3<Scalar> rotaccel = cross(angularaccel,rdiffv);

        // Centrifugal acceleration
        vec3<Scalar> centrifugal = cross(angularvel,rot);

        // Update slave particle velocities
        h_vel.data[j].x = transvel_x + rot.x;
        h_vel.data[j].y = transvel_y + rot.y;
        h_vel.data[j].z = transvel_z + rot.z;

        // Update acceleration of slave particles
        h_accel.data[j].x = transvel_x_dt + rotaccel.x + centrifugal.x;
        h_accel.data[j].y = transvel_y_dt + rotaccel.y + centrifugal.y;
        h_accel.data[j].z = transvel_z_dt + rotaccel.z + centrifugal.z;

        // Update slave particle positions
        h_pos.data[j].x += Scalar(1.0/2.0)*h_vel.data[j].x*m_deltaT;
        h_pos.data[j].y += Scalar(1.0/2.0)*h_vel.data[j].y*m_deltaT;
        h_pos.data[j].z += Scalar(1.0/2.0)*h_vel.data[j].z*m_deltaT;
        }
    }

/*! \param timestep Current time step
    \post particle velocities are moved forward to timestep+1
*/
void RigidBodyIntegrator::integrateStepTwo(uint64_t timestep)
    {
    m_exec_conf->msg->notice(9) << "RigidBodyIntegrator: Integrate Step two" << endl;
    unsigned int group_size = m_group->getNumMembers();

    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);

    // perform the first half step of the integration
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // Update slave particle positions
        h_pos.data[j].x += Scalar(1.0/2.0)*h_vel.data[j].x*m_deltaT;
        h_pos.data[j].y += Scalar(1.0/2.0)*h_vel.data[j].y*m_deltaT;
        h_pos.data[j].z += Scalar(1.0/2.0)*h_vel.data[j].z*m_deltaT;

        }
    }

namespace detail
    {
void export_RigidBodyIntegrator(pybind11::module& m)
    {
    pybind11::class_<RigidBodyIntegrator, SPHIntegrationMethodTwoStep, std::shared_ptr<RigidBodyIntegrator>>(
        m,
        "RigidBodyIntegrator")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, 
                            std::shared_ptr<ParticleGroup>,
                            std::shared_ptr<Variant>,
                            std::shared_ptr<Variant>,
                            std::shared_ptr<Variant>,
                            std::shared_ptr<Variant>,
                            Scalar,
                            Scalar,
                            Scalar,
                            Scalar,
                            Scalar,
                            Scalar>())
        .def("setRotationSpeed",&RigidBodyIntegrator::setRotationSpeed)
        .def("setPivotPoint",&RigidBodyIntegrator::setPivotPoint)
        .def("setRotationAxis",&RigidBodyIntegrator::setRotationAxis)
        .def("setTranslationalVelocity",&RigidBodyIntegrator::setTranslationalVelocity);
    }
    } // end namespace detail
    } // end namespace sph
    } // end namespace hoomd


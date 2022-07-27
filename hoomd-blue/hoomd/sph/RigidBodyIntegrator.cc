// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

/*! \file RigidBodyIntegrator.cc
    \brief Defines the RigidBodyIntegrator class
*/

#include "RigidBodyIntegrator.h"
#include "hoomd/VectorMath.h"

// #include <boost/python.hpp>
// using namespace boost::python;

// #include <boost/bind.hpp>
// using namespace boost;

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#endif

using namespace std;

namespace hoomd 
{
namespace sph
{
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
                                         Scalar rotaxis_z,
                                         bool skip_restart)
    : SPHIntegrationMethodTwoStep(sysdef, group),
    m_transvel_x(transvel_x), m_transvel_y(transvel_y), m_transvel_z(transvel_z), m_rotatvel(rotatvel)
    {
    m_exec_conf->msg->notice(5) << "Constructing RigidBodyIntegrator" << endl;

    // Set Pivot Point
    m_pivotpnt = make_scalar3(pivotpnt_x,pivotpnt_y,pivotpnt_z);

    // Set rotation axis
    vec3<Scalar> rotaxis(rotaxis_x,rotaxis_y,rotaxis_z);
    rotaxis   = rotaxis/Scalar(sqrt(dot(rotaxis,rotaxis)));
    m_rotaxis = vec_to_scalar3(rotaxis);

    if (!skip_restart)
        {
        // set a named, but otherwise blank set of integrator variables
        IntegratorVariables v = getIntegratorVariables();

        if (!restartInfoTestValid(v, "nve", 0))
            {
            v.type = "nve";
            v.variable.resize(0);
            setValidRestart(false);
            }
        else
            setValidRestart(true);

        setIntegratorVariables(v);
        }
    }

RigidBodyIntegrator::~RigidBodyIntegrator()
    {
    m_exec_conf->msg->notice(5) << "Destroying RigidBodyIntegrator" << endl;
    }

/*! \param timestep Current time step
    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2.
*/
void RigidBodyIntegrator::integrateStepOne(unsigned int timestep)
    {

    // profile this step
    if (m_prof)
        m_prof->push("RigidBodyIntegrator Integrate Step 1");

    // Local copy of the simulation box
    const BoxDim& box = m_pdata->getGlobalBox();

    // Evaluate variant expressions
    Scalar transvel_x    = m_transvel_x->getValue(timestep);
    Scalar transvel_y    = m_transvel_y->getValue(timestep);
    Scalar transvel_z    = m_transvel_z->getValue(timestep);
    Scalar rotatvel      = m_rotatvel->getValue(timestep);
    Scalar transvel_x_dt = m_transvel_x->getRate(timestep);
    Scalar transvel_y_dt = m_transvel_y->getRate(timestep);
    Scalar transvel_z_dt = m_transvel_z->getRate(timestep);
    Scalar rotatvel_dt   = m_rotatvel->getRate(timestep);
    vec3<Scalar> angularvel(rotatvel*m_rotaxis.x, rotatvel*m_rotaxis.y, rotatvel*m_rotaxis.z);
    vec3<Scalar> angularaccel(rotatvel_dt*m_rotaxis.x, rotatvel_dt*m_rotaxis.y, rotatvel_dt*m_rotaxis.z);

    // Array handles
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::overwrite);
    ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::overwrite);

    // Loop over group index
    unsigned int group_size = m_group->getNumMembers();
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

    // done profiling
    if (m_prof)
        m_prof->pop();
    }

/*! \param timestep Current time step
    \post particle velocities are moved forward to timestep+1
*/
void RigidBodyIntegrator::integrateStepTwo(unsigned int timestep)
    {
    // profile this step
    if (m_prof)
        m_prof->push("RigidBodyIntegrator Integrate Step 2");

    // Update slave particle velocities
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::overwrite);

    // Loop over group index
    unsigned int group_size = m_group->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // Update slave particle positions
        h_pos.data[j].x += Scalar(1.0/2.0)*h_vel.data[j].x*m_deltaT;
        h_pos.data[j].y += Scalar(1.0/2.0)*h_vel.data[j].y*m_deltaT;
        h_pos.data[j].z += Scalar(1.0/2.0)*h_vel.data[j].z*m_deltaT;

        }

    // done profiling
    if (m_prof)
        m_prof->pop();
    }

namespace detail
{
void export_RigidBodyIntegrator(pybind11::module& m)
    {
    pybind11::class_<RigidBodyIntegrator, std::shared_ptr<RigidBodyIntegrator>>(m, "RigidBodyIntegrator")
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
                            Scalar,
                            bool
                            >())
        .def("setRotationSpeed",&RigidBodyIntegrator::setRotationSpeed)
        .def("setPivotPoint",&RigidBodyIntegrator::setPivotPoint)
        .def("setRotationAxis",&RigidBodyIntegrator::setRotationAxis)
        .def("setTranslationalVelocity",&RigidBodyIntegrator::setTranslationalVelocity)
        ;
    }


} // end namespace detail 
} // end namespace sph 
} // end namespace hoomd 
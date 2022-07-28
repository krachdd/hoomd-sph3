// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

/*! \file VelocityVerlet.cc
    \brief Defines the VelocityVerlet class
*/

#include "VelocityVerlet.h"
#include "hoomd/VectorMath.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#endif

using namespace std;

namespace hoomd
{
namespace sph 
{
/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
    \param group The group of particles this integration method is to work on
    \param skip_restart Skip initialization of the restart information
*/
VelocityVerlet::VelocityVerlet(std::shared_ptr<SystemDefinition> sysdef,
                               std::shared_ptr<ParticleGroup> group,
                               bool skip_restart)
    : SPHIntegrationMethodTwoStep(sysdef, group)
    {
    m_exec_conf->msg->notice(5) << "Constructing VelocityVerlet" << endl;

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

VelocityVerlet::~VelocityVerlet()
    {
    m_exec_conf->msg->notice(5) << "Destroying VelocityVerlet" << endl;
    }

/*! \param timestep Current time step
    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2 per the velocity verlet
          method.
*/
void VelocityVerlet::integrateStepOne(unsigned int timestep)
    {
    // profile this step
    if (m_prof)
        m_prof->push("SPH Integrate Step 1");

    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_dpe(m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_dpedt(m_pdata->getDPEdts(), access_location::host, access_mode::read);

    // perform the first half step of velocity verlet
    unsigned int group_size = m_group->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // dpe(t+deltaT/2) = dpe(t) + (1/2)*dpedt(t)*deltaT
        h_dpe.data[j].x += Scalar(1.0/2.0)*h_dpedt.data[j].x*m_deltaT;
        h_dpe.data[j].y += Scalar(1.0/2.0)*h_dpedt.data[j].y*m_deltaT;
        h_dpe.data[j].z += Scalar(1.0/2.0)*h_dpedt.data[j].z*m_deltaT;

        // v(t+deltaT/2) = v(t) + (1/2)*a(t)*deltaT
        h_vel.data[j].x += Scalar(1.0/2.0)*h_accel.data[j].x*m_deltaT;
        h_vel.data[j].y += Scalar(1.0/2.0)*h_accel.data[j].y*m_deltaT;
        h_vel.data[j].z += Scalar(1.0/2.0)*h_accel.data[j].z*m_deltaT;

        // r(t+deltaT/2) = r(t) + v(t+deltaT/2)*deltaT/2
        h_pos.data[j].x += Scalar(1.0/2.0)*h_vel.data[j].x*m_deltaT;
        h_pos.data[j].y += Scalar(1.0/2.0)*h_vel.data[j].y*m_deltaT;
        h_pos.data[j].z += Scalar(1.0/2.0)*h_vel.data[j].z*m_deltaT;
        }

    // particles may have been moved slightly outside the box by the above steps, wrap them back into place
    const BoxDim& box = m_pdata->getBox();

    ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::readwrite);

    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);
        box.wrap(h_pos.data[j], h_image.data[j]);
        }

    // done profiling
    if (m_prof)
        m_prof->pop();
    }

/*! \param timestep Current time step
    \post particle velocities are moved forward to timestep+1
*/
void VelocityVerlet::integrateStepTwo(unsigned int timestep)
    {
    // profile this step
    if (m_prof)
        m_prof->push("SPH Integrate Step 2");

    // fields to integrate
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_dpe(m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_dpedt(m_pdata->getDPEdts(), access_location::host, access_mode::readwrite);

    // force array
    ArrayHandle<Scalar4> h_net_force(m_pdata->getNetForce(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_net_ratedpe(m_pdata->getNetRateDPEArray(), access_location::host, access_mode::read);

    unsigned int group_size = m_group->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // first, calculate acceleration from the net force
        Scalar minv = Scalar(1.0) / h_vel.data[j].w;
        h_accel.data[j].x = h_net_force.data[j].x*minv;
        h_accel.data[j].y = h_net_force.data[j].y*minv;
        h_accel.data[j].z = h_net_force.data[j].z*minv;

        h_dpedt.data[j].x = h_net_ratedpe.data[j].x;
        h_dpedt.data[j].y = h_net_ratedpe.data[j].y;
        h_dpedt.data[j].z = h_net_ratedpe.data[j].z;

        // dpe(t+deltaT) = dpe(t+deltaT/2) + 1/2 * dpedt(t+deltaT)*deltaT
        h_dpe.data[j].x += Scalar(1.0/2.0)*h_dpedt.data[j].x*m_deltaT;
        h_dpe.data[j].y += Scalar(1.0/2.0)*h_dpedt.data[j].y*m_deltaT;
        h_dpe.data[j].z += Scalar(1.0/2.0)*h_dpedt.data[j].z*m_deltaT;

        // r(t+deltaT) = r(t+deltaT/2) + v(t+deltaT/2)*deltaT/2
        h_pos.data[j].x += Scalar(1.0/2.0)*h_vel.data[j].x*m_deltaT;
        h_pos.data[j].y += Scalar(1.0/2.0)*h_vel.data[j].y*m_deltaT;
        h_pos.data[j].z += Scalar(1.0/2.0)*h_vel.data[j].z*m_deltaT;

        // v(t+deltaT) = v(t+deltaT/2) + 1/2 * a(t+deltaT)*deltaT
        h_vel.data[j].x += Scalar(1.0/2.0)*h_accel.data[j].x*m_deltaT;
        h_vel.data[j].y += Scalar(1.0/2.0)*h_accel.data[j].y*m_deltaT;
        h_vel.data[j].z += Scalar(1.0/2.0)*h_accel.data[j].z*m_deltaT;

        }

    // done profiling
    if (m_prof)
        m_prof->pop();
    }
namespace detail 
{
void export_VelocityVerlet(pybind11::module& m)
    {
    pybind::class_<VelocityVerlet, std::shared_ptr<VelocityVerlet>(m, "VelocityVerlet")
    .def(pybind11::init< std::shared_ptr<SystemDefinition>, 
                         std::shared_ptr<ParticleGroup>, bool >())
        ;
    }
} // end namespace detail 
} // end namespace sph 
} // end namespace hoomd 
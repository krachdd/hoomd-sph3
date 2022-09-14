/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

/*
The Velocity Verlet algorithm originally used in Rakulans old code
*/

#include "VelocityVerlet.h"
#include "hoomd/VectorMath.h"
#include <vector>

using namespace std;

/*! \file VelocityVerlet.h
    \brief Contains code for the VelocityVerlet class
*/

namespace hoomd
    {
namespace sph
    {
/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
    \param group The group of particles this integration method is to work on
*/
VelocityVerlet::VelocityVerlet(std::shared_ptr<SystemDefinition> sysdef,
                       std::shared_ptr<ParticleGroup> group)
    : SPHIntegrationMethodTwoStep(sysdef, group), m_limit(false), m_limit_val(1.0), m_zero_force(false)
    {
    m_exec_conf->msg->notice(5) << "Constructing VelocityVerlet" << endl;
    m_densitymethod_set = false;
    }

VelocityVerlet::~VelocityVerlet()
    {
    m_exec_conf->msg->notice(5) << "Destroying VelocityVerlet" << endl;
    }

/*! \param limit Distance to limit particle movement each time step

    Once the limit is set, future calls to update() will never move a particle
    a distance larger than the limit in a single time step
*/

pybind11::object VelocityVerlet::getLimit()
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

void VelocityVerlet::setLimit(pybind11::object limit)
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

bool VelocityVerlet::getZeroForce()
    {
    return m_zero_force;
    }

void VelocityVerlet::setZeroForce(bool zero_force)
    {
    m_zero_force = zero_force;
    }

/*! \param timestep Current time step
    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2 per the
   velocity verlet method.
*/
void VelocityVerlet::integrateStepOne(uint64_t timestep)
    {
    unsigned int group_size = m_group->getNumMembers();

    m_exec_conf->msg->notice(9) << "VelocityVerlet: Integrate Step one" << endl;

    if (m_densitymethod_set == true){
        if ( m_density_method == DENSITYSUMMATION ){
            std::cout << "Using DENSITYSUMMATION in Verlet" << std::endl;
        }
        else if (m_density_method == DENSITYCONTINUITY ){
            std::cout << "Using DENSITYCONTINUITY in Verlet" << std::endl;
        }
    }


    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_dpe(m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_dpedt(m_pdata->getDPEdts(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);

    // perform the first half step of velocity verlet
    // r(t+deltaT) = r(t) + v(t)*deltaT + (1/2)a(t)*deltaT^2
    // v(t+deltaT/2) = v(t) + (1/2)a*deltaT
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);
        // if (m_zero_force)
        //     h_accel.data[j].x = h_accel.data[j].y = h_accel.data[j].z = 0.0;

        // Original HOOMD Velocity Verlet Two Step NVE
        // Scalar dx = h_vel.data[j].x * m_deltaT + Scalar(1.0 / 2.0) * h_accel.data[j].x * m_deltaT * m_deltaT;
        // Scalar dy = h_vel.data[j].y * m_deltaT + Scalar(1.0 / 2.0) * h_accel.data[j].y * m_deltaT * m_deltaT;
        // Scalar dz = h_vel.data[j].z * m_deltaT + Scalar(1.0 / 2.0) * h_accel.data[j].z * m_deltaT * m_deltaT;

        // limit the movement of the particles
        // if (m_limit)
        //     {
        //     Scalar len = sqrt(dx * dx + dy * dy + dz * dz);
        //     if (len > m_limit_val)
        //         {
        //         dx = dx / len * m_limit_val;
        //         dy = dy / len * m_limit_val;
        //         dz = dz / len * m_limit_val;
        //         }
        //     }

        // h_pos.data[j].x += dx;
        // h_pos.data[j].y += dy;
        // h_pos.data[j].z += dz;

        // David 
        // dpe(t+deltaT/2) = dpe(t) + (1/2)*dpedt(t)*deltaT
        h_dpe.data[j].x += Scalar(1.0/2.0)*h_dpedt.data[j].x*m_deltaT;
        h_dpe.data[j].y += Scalar(1.0/2.0)*h_dpedt.data[j].y*m_deltaT;
        h_dpe.data[j].z += Scalar(1.0/2.0)*h_dpedt.data[j].z*m_deltaT;

        // std::cout << "step 1 mass " << h_vel.data[j].w << std::endl;
        // std::cout << "step 1 accelx " << h_accel.data[j].x << std::endl;

        // Original HOOMD Velocity Verlet Two Step NVE
        h_vel.data[j].x += Scalar(1.0 / 2.0) * h_accel.data[j].x * m_deltaT;
        h_vel.data[j].y += Scalar(1.0 / 2.0) * h_accel.data[j].y * m_deltaT;
        h_vel.data[j].z += Scalar(1.0 / 2.0) * h_accel.data[j].z * m_deltaT;

        // David 
        // r(t+deltaT/2) = r(t) + v(t+deltaT/2)*deltaT/2
        h_pos.data[j].x += Scalar(1.0/2.0)*h_vel.data[j].x*m_deltaT;
        h_pos.data[j].y += Scalar(1.0/2.0)*h_vel.data[j].y*m_deltaT;
        h_pos.data[j].z += Scalar(1.0/2.0)*h_vel.data[j].z*m_deltaT;
        }

    // particles may have been moved slightly outside the box by the above steps, wrap them back
    // into place
    const BoxDim& box = m_pdata->getBox();

    ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::readwrite);

    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);
        box.wrap(h_pos.data[j], h_image.data[j]);
        }
    }

/*! \param timestep Current time step
    \post particle velocities are moved forward to timestep+1
*/
void VelocityVerlet::integrateStepTwo(uint64_t timestep)
    {
    m_exec_conf->msg->notice(9) << "VelocityVerlet: Integrate Step two" << endl;
    unsigned int group_size = m_group->getNumMembers();

    const GlobalArray<Scalar4>& net_force = m_pdata->getNetForce();

    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_dpe(m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_dpedt(m_pdata->getDPEdts(), access_location::host, access_mode::readwrite);

    ArrayHandle<Scalar4> h_net_force(net_force, access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_net_ratedpe(m_pdata->getNetRateDPEArray(), access_location::host, access_mode::read);

    // v(t+deltaT) = v(t+deltaT/2) + 1/2 * a(t+deltaT)*deltaT
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // Original HOOMD Velocity Verlet Two Step NVE
        // if (m_zero_force)
        //     {
        //     h_accel.data[j].x = h_accel.data[j].y = h_accel.data[j].z = 0.0;
        //     }
        // else
        //     {
        //     // first, calculate acceleration from the net force
        //     Scalar minv = Scalar(1.0) / h_vel.data[j].w;
        //     h_accel.data[j].x = h_net_force.data[j].x * minv;
        //     h_accel.data[j].y = h_net_force.data[j].y * minv;
        //     h_accel.data[j].z = h_net_force.data[j].z * minv;
        //     }

        // David 
        // first, calculate acceleration from the net force
        Scalar minv = Scalar(1.0) / h_vel.data[j].w;
        h_accel.data[j].x = h_net_force.data[j].x * minv;
        h_accel.data[j].y = h_net_force.data[j].y * minv;
        h_accel.data[j].z = h_net_force.data[j].z * minv;

        // actually not necessary to compute the next 6 lines if m_densitymethod == DENSITYSUMMATION and j_isfluid 
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

        // Original HOOMD Velocity Verlet Two Step NVE
        // then, update the velocity
        h_vel.data[j].x += Scalar(1.0 / 2.0) * h_accel.data[j].x * m_deltaT;
        h_vel.data[j].y += Scalar(1.0 / 2.0) * h_accel.data[j].y * m_deltaT;
        h_vel.data[j].z += Scalar(1.0 / 2.0) * h_accel.data[j].z * m_deltaT;

        // std::cout << "step 2 vel " << h_vel.data[j].x << std::endl;


        // // limit the movement of the particles
        // if (m_limit)
        //     {
        //     Scalar vel = sqrt(h_vel.data[j].x * h_vel.data[j].x + h_vel.data[j].y * h_vel.data[j].y
        //                       + h_vel.data[j].z * h_vel.data[j].z);
        //     if ((vel * m_deltaT) > m_limit_val)
        //         {
        //         h_vel.data[j].x = h_vel.data[j].x / vel * m_limit_val / m_deltaT;
        //         h_vel.data[j].y = h_vel.data[j].y / vel * m_limit_val / m_deltaT;
        //         h_vel.data[j].z = h_vel.data[j].z / vel * m_limit_val / m_deltaT;
        //         }
        //     }
        }
    }

namespace detail
    {
void export_VelocityVerlet(pybind11::module& m)
    {
    pybind11::class_<VelocityVerlet, SPHIntegrationMethodTwoStep, std::shared_ptr<VelocityVerlet>>(
        m,
        "VelocityVerlet")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<ParticleGroup>>())
        .def("getDensityMethod", &VelocityVerlet::getDensityMethod)
        .def("setDensityMethod", &VelocityVerlet::setDensityMethod)
        .def_property("limit", &VelocityVerlet::getLimit, &VelocityVerlet::setLimit)
        .def_property("zero_force", &VelocityVerlet::getZeroForce, &VelocityVerlet::setZeroForce);
    }
    } // end namespace detail
    } // end namespace sph
    } // end namespace hoomd


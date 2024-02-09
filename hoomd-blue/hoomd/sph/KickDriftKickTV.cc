/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "KickDriftKickTV.h"
#include "hoomd/VectorMath.h"
#include <vector>

using namespace std;

/*! \file KickDriftKickTV.h
    \brief Contains code for the KickDriftKickTV class
*/

namespace hoomd
    {
namespace sph
    {
/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
    \param group The group of particles this integration method is to work on
*/
KickDriftKickTV::KickDriftKickTV(std::shared_ptr<SystemDefinition> sysdef,
                       std::shared_ptr<ParticleGroup> group)
    : SPHIntegrationMethodTwoStep(sysdef, group), m_vlimit(false), m_vlimit_val(0.0), m_xlimit(false), m_xlimit_val(0.0), m_zero_force(false)
    {
    m_exec_conf->msg->notice(5) << "Constructing KickDriftKickTV" << endl;
    m_densitymethod_set = false;
    }

KickDriftKickTV::~KickDriftKickTV()
    {
    m_exec_conf->msg->notice(5) << "Destroying KickDriftKickTV" << endl;
    }

/*! \param limit Distance to limit particle movement each time step

    Once the limit is set, future calls to update() will never move a particle
    a distance larger than the limit in a single time step
*/

Scalar KickDriftKickTV::getxLimit()
    {
    Scalar result;
    if (m_xlimit)
        {
        result = m_xlimit_val;
        }
    else
        {
        result = 0.0;
        }
    return result;
    }

Scalar KickDriftKickTV::getvLimit()
    {
    Scalar result;
    if (m_vlimit)
        {
        result = m_vlimit_val;
        }
    else
        {
        result = 0.0;
        }
    return result;
    }

void KickDriftKickTV::setxLimit(Scalar xlimit)
    {
    if (xlimit <= 0.0)
        {
        m_xlimit = false;
        }
    else
        {
        m_xlimit = true;
        m_xlimit_val = xlimit;
        }
    }

void KickDriftKickTV::setvLimit(Scalar vlimit)
    {
    if (vlimit <= 0.0)
        {
        m_vlimit = false;
        }
    else
        {
        m_vlimit = true;
        m_vlimit_val = vlimit;
        }
    }

bool KickDriftKickTV::getZeroForce()
    {
    return m_zero_force;
    }

void KickDriftKickTV::setZeroForce(bool zero_force)
    {
    m_zero_force = zero_force;
    }

/*! \param timestep Current time step
    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2 per the
   velocity verlet method.
*/
void KickDriftKickTV::integrateStepOne(uint64_t timestep)
    {
    unsigned int group_size = m_group->getNumMembers();

    m_exec_conf->msg->notice(9) << "KickDriftKickTV: Integrate Step one" << endl;

    // if (m_densitymethod_set == true){
    //     if ( m_density_method == DENSITYSUMMATION ){
    //         std::cout << "Using DENSITYSUMMATION in Verlet" << std::endl;
    //     }
    //     else if (m_density_method == DENSITYCONTINUITY ){
    //         std::cout << "Using DENSITYCONTINUITY in Verlet" << std::endl;
    //     }
    // }


    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar>  h_density(m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar>  h_pressure(m_pdata->getPressures(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_dpedt(m_pdata->getDPEdts(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_bpc(this->m_pdata->getAuxiliaries2(), access_location::host,access_mode::read); // background pressure contribution to tv
    ArrayHandle<Scalar3> h_tv(this->m_pdata->getAuxiliaries3(), access_location::host,access_mode::readwrite); // transport velocity of the particle tv

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

        // // limit the movement of the particles
        // if (m_xlimit)
        //     {
        //     Scalar len = sqrt(dx * dx + dy * dy + dz * dz);
        //     if (len > m_xlimit_val)
        //         {
        //         dx = dx / len * m_xlimit_val;
        //         dy = dy / len * m_xlimit_val;
        //         dz = dz / len * m_xlimit_val;
        //         }
        //     }

        // h_pos.data[j].x += dx;
        // h_pos.data[j].y += dy;
        // h_pos.data[j].z += dz;

        // David 
        // dpe(t+deltaT/2) = dpe(t) + (1/2)*dpedt(t)*deltaT
        h_density.data[j] += Scalar(1.0/2.0)*h_dpedt.data[j].x*m_deltaT;
        h_pressure.data[j] += Scalar(1.0/2.0)*h_dpedt.data[j].y*m_deltaT;
        // DK: Energy change can be ignored
        // h_dpe.data[j].z += Scalar(1.0/2.0)*h_dpedt.data[j].z*m_deltaT;

        // Original HOOMD Velocity Verlet Two Step NVE
        h_vel.data[j].x += Scalar(1.0 / 2.0) * h_accel.data[j].x * m_deltaT;
        h_vel.data[j].y += Scalar(1.0 / 2.0) * h_accel.data[j].y * m_deltaT;
        h_vel.data[j].z += Scalar(1.0 / 2.0) * h_accel.data[j].z * m_deltaT;

        // David 
        // r(t+deltaT) = r(t) + v(t+deltaT/2)*deltaT
        // h_pos.data[j].x += h_vel.data[j].x*m_deltaT;
        // h_pos.data[j].y += h_vel.data[j].y*m_deltaT;
        // h_pos.data[j].z += h_vel.data[j].z*m_deltaT;

        // Advection Velocity
        h_tv.data[j].x = h_vel.data[j].x + Scalar(1.0 / 2.0) * h_bpc.data[j].x * m_deltaT; 
        h_tv.data[j].y = h_vel.data[j].y + Scalar(1.0 / 2.0) * h_bpc.data[j].y * m_deltaT; 
        h_tv.data[j].z = h_vel.data[j].z + Scalar(1.0 / 2.0) * h_bpc.data[j].z * m_deltaT; 


        Scalar dx = h_tv.data[j].x * m_deltaT;
        Scalar dy = h_tv.data[j].y * m_deltaT;
        Scalar dz = h_tv.data[j].z * m_deltaT;

        // limit the movement of the particles
        if (m_xlimit)
            {
            if (sqrt(dx*dx) > m_xlimit_val)
                {
                dx = m_xlimit_val;
                }
            if (sqrt(dy*dy) > m_xlimit_val)
                {
                dy = m_xlimit_val;
                }
            if (sqrt(dz*dz) > m_xlimit_val)
                {
                dz = m_xlimit_val;
                }
            }

        // Update position with transport veloicity
        // r(t+deltaT) = r(t) + v(t+deltaT/2)*deltaT
        h_pos.data[j].x += dx;
        h_pos.data[j].y += dy;
        h_pos.data[j].z += dz;

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
void KickDriftKickTV::integrateStepTwo(uint64_t timestep)
    {
    m_exec_conf->msg->notice(9) << "KickDriftKickTV: Integrate Step two" << endl;
    unsigned int group_size = m_group->getNumMembers();

    const GlobalArray<Scalar4>& net_force = m_pdata->getNetForce();

    // ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_density(m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_pressure(m_pdata->getPressures(), access_location::host, access_mode::readwrite);
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
        // h_dpedt.data[j].z = h_net_ratedpe.data[j].z;

        // dpe(t+deltaT) = dpe(t+deltaT/2) + 1/2 * dpedt(t+deltaT)*deltaT
        h_density.data[j] += Scalar(1.0/2.0)*h_dpedt.data[j].x*m_deltaT;
        h_pressure.data[j] += Scalar(1.0/2.0)*h_dpedt.data[j].y*m_deltaT;
        // h_dpe.data[j].z += Scalar(1.0/2.0)*h_dpedt.data[j].z*m_deltaT;


        // not done in original VV Algorithm
        // // r(t+deltaT) = r(t+deltaT/2) + v(t+deltaT/2)*deltaT/2
        // h_pos.data[j].x += Scalar(1.0/2.0)*h_vel.data[j].x*m_deltaT;
        // h_pos.data[j].y += Scalar(1.0/2.0)*h_vel.data[j].y*m_deltaT;
        // h_pos.data[j].z += Scalar(1.0/2.0)*h_vel.data[j].z*m_deltaT;

        // Original HOOMD Velocity Verlet Two Step NVE
        // then, update the velocity
        // v(t+deltaT) = v(t+deltaT/2) + 1/2 * a(t+deltaT) deltaT
        h_vel.data[j].x += Scalar(1.0 / 2.0) * h_accel.data[j].x * m_deltaT;
        h_vel.data[j].y += Scalar(1.0 / 2.0) * h_accel.data[j].y * m_deltaT;
        h_vel.data[j].z += Scalar(1.0 / 2.0) * h_accel.data[j].z * m_deltaT;

        // limit the movement of the particles
        if (m_vlimit)
            {
            if ( sqrt( h_vel.data[j].x * h_vel.data[j].x ) > m_vlimit_val )
                {
                h_vel.data[j].x = m_vlimit_val;
                }
            if ( sqrt( h_vel.data[j].y * h_vel.data[j].y ) > m_vlimit_val )
                {
                h_vel.data[j].y = m_vlimit_val;
                }
            if ( sqrt( h_vel.data[j].z * h_vel.data[j].z ) > m_vlimit_val )
                {
                h_vel.data[j].z = m_vlimit_val;
                }
            }
        }
    }

namespace detail
    {
void export_KickDriftKickTV(pybind11::module& m)
    {
    pybind11::class_<KickDriftKickTV, SPHIntegrationMethodTwoStep, std::shared_ptr<KickDriftKickTV>>(
        m,
        "KickDriftKickTV")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<ParticleGroup>>())
        .def("getDensityMethod", &KickDriftKickTV::getDensityMethod)
        .def("setDensityMethod", &KickDriftKickTV::setDensityMethod)
        .def("getvLimit", &KickDriftKickTV::getvLimit)
        .def("getxLimit", &KickDriftKickTV::getxLimit)
        .def("setvLimit", &KickDriftKickTV::setvLimit)
        .def("setxLimit", &KickDriftKickTV::setxLimit)
        // .def_property("limit", &KickDriftKickTV::getLimit, &KickDriftKickTV::setLimit)
        .def_property("zero_force", &KickDriftKickTV::getZeroForce, &KickDriftKickTV::setZeroForce);
    }
    } // end namespace detail
    } // end namespace sph
    } // end namespace hoomd


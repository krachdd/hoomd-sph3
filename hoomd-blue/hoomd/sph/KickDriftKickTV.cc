/* ---------------------------------------------------------
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
    const unsigned int group_size = m_group->getNumMembers();
    const BoxDim& box = m_pdata->getBox();

    m_exec_conf->msg->notice(9) << "KickDriftKickTV: Integrate Step one" << endl;

        { // GPU Array Scope
        ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar>  h_density(m_pdata->getDensities(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar>  h_pressure(m_pdata->getPressures(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::read);
        ArrayHandle<Scalar3> h_dpedt(m_pdata->getDPEdts(), access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar3> h_bpc(this->m_pdata->getAuxiliaries2(), access_location::host,access_mode::read); // background pressure contribution to tv
        ArrayHandle<Scalar3> h_tv(this->m_pdata->getAuxiliaries3(), access_location::host,access_mode::readwrite); // transport velocity of the particle tv
        // aux4.x stores the scalar field T; h_dpedt.z carries dT/dt (SinglePhaseFlowGDGD only).
        ArrayHandle<Scalar3> h_aux4(this->m_pdata->getAuxiliaries4(), access_location::host, access_mode::readwrite);

        // perform the first half step of velocity verlet
        // r(t+deltaT) = r(t) + v(t)*deltaT + (1/2)a(t)*deltaT^2
        // v(t+deltaT/2) = v(t) + (1/2)a*deltaT
        Scalar temp0 = 0.0;
        // Acquire the group index array once: getMemberIndex() acquires an
        // ArrayHandle per call, which is not thread-safe inside the parallel loop
        ArrayHandle<unsigned int> h_members_omp1(m_group->getIndexArray(), access_location::host, access_mode::read);
        #pragma omp parallel for firstprivate(temp0)
        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            unsigned int j = h_members_omp1.data[group_idx];

            // dpe(t+deltaT/2) = dpe(t) + (1/2)*dpedt(t)*deltaT
            temp0 = Scalar(1.0 / 2.0) * m_deltaT;
            h_density.data[j]  += temp0 * h_dpedt.data[j].x;
            h_pressure.data[j] += temp0 * h_dpedt.data[j].y;
            // Scalar field half-step: T(t+dt/2) = T(t) + (1/2)*dT/dt*dt
            // h_dpedt.data[j].z holds dT/dt set by SinglePhaseFlowGDGD (zero for all other solvers).
            h_aux4.data[j].x   += temp0 * h_dpedt.data[j].z;

            // Original HOOMD Velocity Verlet Two Step NVE
            h_vel.data[j].x += temp0 * h_accel.data[j].x;
            h_vel.data[j].y += temp0 * h_accel.data[j].y;
            h_vel.data[j].z += temp0 * h_accel.data[j].z;

            // Advection Velocity
            // temp0 = Scalar(1.0 / 2.0) * 1.0/h_vel.data[j].w * m_deltaT;
            // temp0 = Scalar(1.0 / 2.0) * m_deltaT;
            h_tv.data[j].x = h_vel.data[j].x + temp0 * h_bpc.data[j].x; 
            h_tv.data[j].y = h_vel.data[j].y + temp0 * h_bpc.data[j].y; 
            h_tv.data[j].z = h_vel.data[j].z + temp0 * h_bpc.data[j].z; 


            Scalar dx = h_tv.data[j].x * m_deltaT;
            Scalar dy = h_tv.data[j].y * m_deltaT;
            Scalar dz = h_tv.data[j].z * m_deltaT;

            // limit the movement of the particles
            // Note: use !(fabs(d) <= limit) rather than fabs(d) > limit so that
            // NaN displacements are caught (NaN comparisons always return false).
            if (m_xlimit)
                {
                if (!(fabs(dx) <= m_xlimit_val))
                    dx = (dx > Scalar(0)) ? m_xlimit_val : -m_xlimit_val;
                if (!(fabs(dy) <= m_xlimit_val))
                    dy = (dy > Scalar(0)) ? m_xlimit_val : -m_xlimit_val;
                if (!(fabs(dz) <= m_xlimit_val))
                    dz = (dz > Scalar(0)) ? m_xlimit_val : -m_xlimit_val;
                }

            // Update position with transport veloicity
            // r(t+deltaT) = r(t) + v(t+deltaT/2)*deltaT
            h_pos.data[j].x += dx;
            h_pos.data[j].y += dy;
            h_pos.data[j].z += dz;

            }

        // particles may have been moved slightly outside the box by the above steps, wrap them back
        // into place
        ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::readwrite);

        #pragma omp parallel for
        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            unsigned int j = h_members_omp1.data[group_idx];
            box.wrap(h_pos.data[j], h_image.data[j]);
            }
        } // End GPU Array Scope
    }

/*! \param timestep Current time step
    \post particle velocities are moved forward to timestep+1
*/
void KickDriftKickTV::integrateStepTwo(uint64_t timestep)
    {
    m_exec_conf->msg->notice(9) << "KickDriftKickTV: Integrate Step two" << endl;
    const unsigned int group_size = m_group->getNumMembers();

    const GPUArray<Scalar4>& net_force = m_pdata->getNetForce();

        { // GPU Array Scope
        ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar> h_density(m_pdata->getDensities(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar> h_pressure(m_pdata->getPressures(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar3> h_dpedt(m_pdata->getDPEdts(), access_location::host, access_mode::readwrite);
        // aux4.x stores the scalar field T (temperature / concentration).
        // Advanced here with dT/dt set by SinglePhaseFlowGDGD; zero for all other solvers.
        ArrayHandle<Scalar3> h_aux4(m_pdata->getAuxiliaries4(), access_location::host, access_mode::readwrite);

        ArrayHandle<Scalar4> h_net_force(net_force, access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_net_ratedpe(m_pdata->getNetRateDPEArray(), access_location::host, access_mode::read);

        Scalar minv = Scalar(0);
        Scalar temp0 = Scalar(0);

        // v(t+deltaT) = v(t+deltaT/2) + 1/2 * a(t+deltaT)*deltaT
        // Acquire the group index array once: getMemberIndex() acquires an
        // ArrayHandle per call, which is not thread-safe inside the parallel loop
        ArrayHandle<unsigned int> h_members_omp3(m_group->getIndexArray(), access_location::host, access_mode::read);
        #pragma omp parallel for firstprivate(minv, temp0)
        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            unsigned int j = h_members_omp3.data[group_idx];

            // first, calculate acceleration from the net force
            minv = Scalar(1.0) / h_vel.data[j].w;
            h_accel.data[j].x = h_net_force.data[j].x * minv;
            h_accel.data[j].y = h_net_force.data[j].y * minv;
            h_accel.data[j].z = h_net_force.data[j].z * minv;

            // Copy net rates from the accumulated force/rate arrays into the per-particle
            // rate array (used for the half-step on the NEXT timestep's step 1).
            // The .z component carries dT/dt when SinglePhaseFlowGDGD is active; zero otherwise.
            h_dpedt.data[j].x = h_net_ratedpe.data[j].x;
            h_dpedt.data[j].y = h_net_ratedpe.data[j].y;
            h_dpedt.data[j].z = h_net_ratedpe.data[j].z; // dT/dt (zero for non-GDGD solvers)

            temp0 = Scalar(1.0/2.0) * m_deltaT;
            // dpe(t+deltaT) = dpe(t+deltaT/2) + 1/2 * dpedt(t+deltaT)*deltaT
            h_density.data[j]  += temp0 * h_dpedt.data[j].x;
            h_pressure.data[j] += temp0 * h_dpedt.data[j].y;
            // Second scalar half-step: T(t+dt) = T(t+dt/2) + (1/2)*dT/dt*dt
            h_aux4.data[j].x   += temp0 * h_dpedt.data[j].z;

            // v(t+deltaT) = v(t+deltaT/2) + 1/2 * a(t+deltaT) deltaT
            h_vel.data[j].x += temp0 * h_accel.data[j].x;
            h_vel.data[j].y += temp0 * h_accel.data[j].y;
            h_vel.data[j].z += temp0 * h_accel.data[j].z;

            // limit the velocity of the particles
            // Note: use !(fabs(v) <= limit) rather than fabs(v) > limit so that
            // NaN velocities are caught (NaN comparisons always return false, so
            // fabs(NaN) > limit is false and NaN would bypass the clamp).
            if (m_vlimit)
                {
                if (!(fabs(h_vel.data[j].x) <= m_vlimit_val))
                    h_vel.data[j].x = (h_vel.data[j].x > Scalar(0)) ? m_vlimit_val : -m_vlimit_val;
                if (!(fabs(h_vel.data[j].y) <= m_vlimit_val))
                    h_vel.data[j].y = (h_vel.data[j].y > Scalar(0)) ? m_vlimit_val : -m_vlimit_val;
                if (!(fabs(h_vel.data[j].z) <= m_vlimit_val))
                    h_vel.data[j].z = (h_vel.data[j].z > Scalar(0)) ? m_vlimit_val : -m_vlimit_val;
                }
            }
        } // End GPU Array Scope
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
        .def_property("zero_force", &KickDriftKickTV::getZeroForce, &KickDriftKickTV::setZeroForce);
    }
    } // end namespace detail
    } // end namespace sph
    } // end namespace hoomd


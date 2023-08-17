/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SolidAggregateIntegrator.h"
#include "hoomd/VectorMath.h"
#include <vector>

using namespace std;

/*! \file SolidAggregateIntegrator.h
    \brief Contains code for the SolidAggregateIntegrator class
*/

namespace hoomd
    {
namespace sph
    {
/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
    \param group The group of particles this integration method is to work on
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
SolidAggregateIntegrator<KT_, SET_>::SolidAggregateIntegrator(std::shared_ptr<SystemDefinition> sysdef,
                       std::shared_ptr<ParticleGroup> group, std::shared_ptr<SmoothingKernel<KT_> > skernel,
                        std::shared_ptr<StateEquation<SET_> > equationofstate,)
    : SPHIntegrationMethodTwoStep<KT_, SET_>(sysdef, group, skernel, equationofstate), m_limit(false), m_limit_val(1.0), m_zero_force(false)
    {
    m_exec_conf->msg->notice(5) << "Constructing SolidAggregateIntegrator" << endl;
    m_densitymethod_set = false;
    }

template<SmoothingKernelType KT_,StateEquationType SET_>
SolidAggregateIntegrator<KT_, SET_>::~SolidAggregateIntegrator()
    {
    m_exec_conf->msg->notice(5) << "Destroying SolidAggregateIntegrator" << endl;
    }

/*! \param limit Distance to limit particle movement each time step

    Once the limit is set, future calls to update() will never move a particle
    a distance larger than the limit in a single time step
*/

template<SmoothingKernelType KT_,StateEquationType SET_>
pybind11::object SolidAggregateIntegrator<KT_, SET_>::getLimit()
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

template<SmoothingKernelType KT_,StateEquationType SET_>
void SolidAggregateIntegrator<KT_, SET_>::setLimit(pybind11::object limit)
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

template<SmoothingKernelType KT_,StateEquationType SET_>
bool SolidAggregateIntegrator<KT_, SET_>::getZeroForce()
    {
    return m_zero_force;
    }

template<SmoothingKernelType KT_,StateEquationType SET_>
void SolidAggregateIntegrator<KT_, SET_>::setZeroForce(bool zero_force)
    {
    m_zero_force = zero_force;
    }

/*! \param timestep Current time step
    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2 per the
   velocity verlet method.
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
void SolidAggregateIntegrator<KT_, SET_>::integrateStepOne(uint64_t timestep)
    {
    unsigned int group_size = m_group->getNumMembers();

    m_exec_conf->msg->notice(9) << "SolidAggregateIntegrator: Integrate Step one" << endl;

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

    ArrayHandle<unsigned int> h_body(this->m_pdata->getBodies(),
                                     access_location::host,
                                     access_mode::read);
    ArrayHandle<unsigned int> h_rtag(this->m_pdata->getRTags(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_tag(this->m_pdata->getTags(), access_location::host, access_mode::read);

    // access body positions and orientations
    ArrayHandle<Scalar3> h_body_pos(m_body_pos, access_location::host, access_mode::read);
    // ArrayHandle<Scalar4> h_body_orientation(m_body_orientation,
    //                                         access_location::host,
    //                                         access_mode::read);
    ArrayHandle<unsigned int> h_body_len(m_body_len, access_location::host, access_mode::read);


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
        h_pos.data[j].x += h_vel.data[j].x*m_deltaT;
        h_pos.data[j].y += h_vel.data[j].y*m_deltaT;
        h_pos.data[j].z += h_vel.data[j].z*m_deltaT;
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
template<SmoothingKernelType KT_,StateEquationType SET_>
void SolidAggregateIntegrator<KT_, SET_>::integrateStepTwo(uint64_t timestep)
    {
    m_exec_conf->msg->notice(9) << "SolidAggregateIntegrator: Integrate Step two" << endl;
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
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SolidAggregateIntegrator(pybind11::module& m)
    {
    pybind11::class_<SolidAggregateIntegrator<KT_, SET_>, SPHIntegrationMethodTwoStep<KT_, SET_>, std::shared_ptr<SolidAggregateIntegrator<KT_, SET_>>>(
        m,
        "SolidAggregateIntegrator")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<ParticleGroup>>())
        .def("getDensityMethod", &SolidAggregateIntegrator<KT_, SET_>::getDensityMethod)
        .def("setDensityMethod", &SolidAggregateIntegrator<KT_, SET_>::setDensityMethod)
        .def_property("limit", &SolidAggregateIntegrator<KT_, SET_>::getLimit, &SolidAggregateIntegrator::setLimit)
        .def_property("zero_force", &SolidAggregateIntegrator<KT_, SET_>::getZeroForce, &SolidAggregateIntegrator::setZeroForce);
    }
    } // end namespace detail
    } // end namespace sph
    } // end namespace hoomd


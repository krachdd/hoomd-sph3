/* ---------------------------------------------------------
maintainer: nkijanksi, nadine.kijanksi@mib.uni-stuttgart.de
---------------------------------------------------------*/

/*
The Suspendend Object Integrator algorithm originally used in Nadines code
*/

#include "SuspendedObjectIntegrator.h"
#include "hoomd/VectorMath.h"
#include <vector>

#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#include "hoomd/HOOMDMPI.h"
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
// noch nicht sicher ob hier gebraucht

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
    // profile this step
    if (m_prof)
        m_prof->push("Suspended Object Compute Center of Mass");

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
    MPI_Allreduce(MPI_IN_PLACE, &angles[0], 6, MPI_HOOSPH_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
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

    // done profiling
    if (m_prof)
        m_prof->pop();
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
    MPI_Allreduce(MPI_IN_PLACE, &translationvel[0], 3, MPI_HOOSPH_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
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
    MPI_Allreduce(MPI_IN_PLACE, &J[0], 9, MPI_HOOSPH_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator());
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


/*! \param timestep Current time step
    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2 per the
   velocity verlet method.
*/
void VelocityVerlet::integrateStepOne(uint64_t timestep)
    {
    unsigned int group_size = m_group->getNumMembers();

    m_exec_conf->msg->notice(9) << "VelocityVerlet: Integrate Step one" << endl;

    // if (m_densitymethod_set == true){
    //     if ( m_density_method == DENSITYSUMMATION ){
    //         std::cout << "Using DENSITYSUMMATION in Verlet" << std::endl;
    //     }
    //     else if (m_density_method == DENSITYCONTINUITY ){
    //         std::cout << "Using DENSITYCONTINUITY in Verlet" << std::endl;
    //     }
    // }


    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
    // ArrayHandle<Scalar3> h_dpe(m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_density(m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_pressure(m_pdata->getPressures(), access_location::host, access_mode::readwrite);
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
        h_density.data[j] += Scalar(1.0/2.0)*h_dpedt.data[j].x*m_deltaT;
        h_pressure.data[j] += Scalar(1.0/2.0)*h_dpedt.data[j].y*m_deltaT;
        // DK: Energy change can be ignored
        // h_dpe.data[j].z += Scalar(1.0/2.0)*h_dpedt.data[j].z*m_deltaT;

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
    ArrayHandle<Scalar>  h_density(m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar>  h_pressure(m_pdata->getPressures(), access_location::host, access_mode::readwrite);
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
        // DK: Energy change can be ignored
        // h_dpedt.data[j].z = h_net_ratedpe.data[j].z;

        // dpe(t+deltaT) = dpe(t+deltaT/2) + 1/2 * dpedt(t+deltaT)*deltaT
        h_density.data[j] += Scalar(1.0/2.0)*h_dpedt.data[j].x*m_deltaT;
        h_pressure.data[j] += Scalar(1.0/2.0)*h_dpedt.data[j].y*m_deltaT;
        // h_dpe.data[j].z += Scalar(1.0/2.0)*h_dpedt.data[j].z*m_deltaT;

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


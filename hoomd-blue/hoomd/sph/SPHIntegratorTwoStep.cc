/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SPHIntegratorTwoStep.h"

#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#endif

#include <pybind11/stl_bind.h>
PYBIND11_MAKE_OPAQUE(std::vector<std::shared_ptr<hoomd::sph::SPHIntegrationMethodTwoStep>>);

using namespace std;

namespace hoomd
    {
namespace sph
    {
SPHIntegratorTwoStep::SPHIntegratorTwoStep(std::shared_ptr<SystemDefinition> sysdef, Scalar deltaT)
    : Integrator(sysdef, deltaT), m_prepared(false)
    {
    m_exec_conf->msg->notice(5) << "Constructing SPHIntegratorTwoStep" << endl;

// #ifdef ENABLE_MPI
//     if (m_sysdef->isDomainDecomposed())
//         {
//         m_comm->getComputeCallbackSignal()
//             .connect<SPHIntegratorTwoStep, &SPHIntegratorTwoStep::updateRigidBodies>(this);
//         }
// #endif
    }

SPHIntegratorTwoStep::~SPHIntegratorTwoStep()
    {
    m_exec_conf->msg->notice(5) << "Destroying SPHIntegratorTwoStep" << endl;

// #ifdef ENABLE_MPI
//     if (m_sysdef->isDomainDecomposed())
//         {
//         m_comm->getComputeCallbackSignal()
//             .disconnect<SPHIntegratorTwoStep, &SPHIntegratorTwoStep::updateRigidBodies>(this);
//         }
// #endif
    }

/*! \param timestep Current time step of the simulation
    \post All integration methods in m_methods are applied in order to move the system state
    variables forward to \a timestep+1.
    \post Internally, all forces present in the m_forces std::vector are evaluated at \a timestep+1
*/
void SPHIntegratorTwoStep::update(uint64_t timestep)
    {
    Integrator::update(timestep);
    // ensure that prepRun() has been called
    assert(m_prepared);

    // perform the first step of the integration on all groups
    for (auto& method : m_methods)
        {
        // deltaT should probably be passed as an argument, but that would require modifying many
        // files. Work around this by calling setDeltaT every timestep.
        // method->setAnisotropic(m_integrate_rotational_dof);
        method->setDeltaT(m_deltaT);
        method->integrateStepOne(timestep);
        }

#ifdef ENABLE_MPI
    if (m_sysdef->isDomainDecomposed())
        {
        // Update the rigid body consituent particles before communicating so that any such
        // particles that move from one domain to another are migrated.
        // updateRigidBodies(timestep + 1);

        // perform all necessary communication steps. This ensures
        // a) that particles have migrated to the correct domains
        // b) that forces are calculated correctly, if ghost atom positions are updated every time
        // step
        m_comm->communicate(timestep + 1);

        // Communicator uses a compute callback to trigger updateRigidBodies again and ensure that
        // all ghost constituent particle positions are set in accordance with any just communicated
        // ghost and/or migrated rigid body centers.
        }
    // else
#endif
        // {
        // // Update rigid body constituent particles in serial simulations.
        // updateRigidBodies(timestep + 1);
        // }


    // compute the net force on all particles
#ifdef ENABLE_HIP
    if (m_exec_conf->isCUDAEnabled())
        computeNetForceGPU(timestep + 1);
    else
#endif
        computeNetForce(timestep + 1);

    // Call HalfStep hook
    if (m_half_step_hook)
        {
        m_half_step_hook->update(timestep + 1);
        }

    // perform the second step of the integration on all groups
    for (auto& method : m_methods)
        {
        method->integrateStepTwo(timestep);
        // method->includeRATTLEForce(timestep + 1);
        }

    /* NOTE: For composite particles, it is assumed that positions and orientations are not updated
       in the second step.

       Otherwise we would have to update ghost positions for central particles
       here in order to update the constituent particles.

       TODO: check this assumptions holds for all integrators
     */

    #ifdef ENABLE_MPI
    if (m_comm)
        {
        // perform all necessary communication steps. This ensures
        // a) that particles have migrated to the correct domains
        // b) that forces are calculated correctly, if ghost atom positions are updated every time step

        // also updates rigid bodies after ghost updating
        m_comm->communicate(timestep+1);
        }
    #endif
    }

/*! \param deltaT new deltaT to set
    \post \a deltaT is also set on all contained integration methods
*/
void SPHIntegratorTwoStep::setDeltaT(Scalar deltaT)
    {
    Integrator::setDeltaT(deltaT);

    // set deltaT on all methods already added
    for (auto& method : m_methods)
        {
        method->setDeltaT(deltaT);
        }
    // if (m_rigid_bodies)
    //     {
    //     m_rigid_bodies->setDeltaT(deltaT);
    //     }
    }

/*! \param group Group over which to count degrees of freedom.

    SPHIntegratorTwoStep totals up the degrees of freedom that each integration method provide to the
    group.

    When the user has only one momentum conserving integration method applied to the all group,
    getNDOF subtracts n_dimensions degrees of freedom from the system to account for the pinned
    center of mass. When the query group is not the group of all particles, spread these these
    removed DOF proportionately so that the results given by one ComputeThermo on the all group are
    consitent with the average of many ComputeThermo's on disjoint subset groups.
*/
Scalar SPHIntegratorTwoStep::getTranslationalDOF(std::shared_ptr<ParticleGroup> group)
    {
    Scalar periodic_dof_removed = 0;

    unsigned int N_filter = group->getNumMembersGlobal();
    unsigned int N_particles = m_pdata->getNGlobal();

    // When using rigid bodies, adjust the number of particles to the number of rigid centers and
    // free particles. The constituent particles are in the system, but not part of the equations
    // of motion.
    // if (m_rigid_bodies)
    //     {
    //     m_rigid_bodies->validateRigidBodies();
    //     N_particles
    //         = m_rigid_bodies->getNMoleculesGlobal() + m_rigid_bodies->getNFreeParticlesGlobal();
    //     N_filter = group->getNCentralAndFreeGlobal();
    //     }
    // proportionately remove n_dimensions DOF when there is only one momentum conserving
    // integration method
    if (m_methods.size() == 1 && m_methods[0]->isMomentumConserving()
        && m_methods[0]->getGroup()->getNumMembersGlobal() == N_particles)
        {
        periodic_dof_removed
            = Scalar(m_sysdef->getNDimensions()) * (Scalar(N_filter) / Scalar(N_particles));
        }

    // loop through all methods and add up the number of DOF They apply to the group
    Scalar total = 0;
    for (auto& method : m_methods)
        {
        total += method->getTranslationalDOF(group);
        }

    return total - periodic_dof_removed - getNDOFRemoved(group);
    }

/*! \param group Group over which to count degrees of freedom.
    SPHIntegratorTwoStep totals up the rotational degrees of freedom that each integration method
   provide to the group.
*/
// Scalar SPHIntegratorTwoStep::getRotationalDOF(std::shared_ptr<ParticleGroup> group)
//     {
//     double res = 0;

//     if (m_integrate_rotational_dof)
//         {
//         for (auto& method : m_methods)
//             {
//             res += method->getRotationalDOF(group);
//             }
//         }

//     return res;
//     }

/*!  \param integrate_rotational_dofs true to integrate orientations, false to not
 */
// void SPHIntegratorTwoStep::setIntegrateRotationalDOF(bool integrate_rotational_dof)
//     {
//     m_integrate_rotational_dof = integrate_rotational_dof;
//     }

// const bool SPHIntegratorTwoStep::getIntegrateRotationalDOF()
//     {
//     return m_integrate_rotational_dof;
//     }

/*! Compute accelerations if needed for the first step.
    If acceleration is available in the restart file, then just call computeNetForce so that
    net_force and net_virial are available in Python. This solves ticket #393
*/
void SPHIntegratorTwoStep::prepRun(uint64_t timestep)
    {
    Integrator::prepRun(timestep);

#ifdef ENABLE_MPI
    if (m_sysdef->isDomainDecomposed())
        {
        // force particle migration and ghost exchange
        m_comm->forceMigrate();

        // perform communication
        m_comm->communicate(timestep);
        }
    // else
#endif
        // if (m_rigid_bodies)
        // {
        // m_rigid_bodies->validateRigidBodies();
        // updateRigidBodies(timestep);
        // }

    // compute the net force on all particles
#ifdef ENABLE_HIP
    if (m_exec_conf->isCUDAEnabled())
        computeNetForceGPU(timestep);
    else
#endif
    computeNetForce(timestep);

    // accelerations only need to be calculated if the accelerations have not yet been set
    if (!m_pdata->isAccelSet())
        {
        computeAccelerations(timestep);
        m_pdata->notifyAccelSet();
        }

    // for (auto& method : m_methods)
    //     method->includeRATTLEForce(timestep);

    m_prepared = true;
    }

/*! Return the combined flags of all integration methods.
 */
PDataFlags SPHIntegratorTwoStep::getRequestedPDataFlags()
    {
    PDataFlags flags;

    // loop through all methods
    for (auto& method : m_methods)
        {
        // or them all together
        flags |= method->getRequestedPDataFlags();
        }

    return flags;
    }

// //! Updates the rigid body constituent particles
// void SPHIntegratorTwoStep::updateRigidBodies(uint64_t timestep)
//     {
//     // update the composite particle positions of any rigid bodies
//     if (m_rigid_bodies)
//         {
//         m_rigid_bodies->updateCompositeParticles(timestep);
//         }
//     }


void SPHIntegratorTwoStep::startAutotuning()
    {
    Integrator::startAutotuning();

    // Start autotuning in all methods.
    for (auto& method : m_methods)
        method->startAutotuning();
    }

/// Check if autotuning is complete.
bool SPHIntegratorTwoStep::isAutotuningComplete()
    {
    bool result = Integrator::isAutotuningComplete();
    for (auto& method : m_methods)
        {
        result = result && method->isAutotuningComplete();
        }
    return result;
    }

/// helper function to compute net force
void SPHIntegratorTwoStep::computeNetForce(uint64_t timestep)
    {
    // if (m_rigid_bodies)
    //     {
    //     m_rigid_bodies->validateRigidBodies();
    //     m_constraint_forces.push_back(m_rigid_bodies);
    //     }
    Integrator::computeNetForce(timestep);
    // if (m_rigid_bodies)
    //     {
    //     m_constraint_forces.pop_back();
    //     }
    }

#ifdef ENABLE_HIP
/// helper function to compute net force/virial on the GPU
void SPHIntegratorTwoStep::computeNetForceGPU(uint64_t timestep)
    {
    // if (m_rigid_bodies)
    //     {
    //     m_rigid_bodies->validateRigidBodies();
    //     m_constraint_forces.push_back(m_rigid_bodies);
    //     }
    Integrator::computeNetForceGPU(timestep);
    // if (m_rigid_bodies)
    //     {
    //     m_constraint_forces.pop_back();
    //     }
    }
#endif

#ifdef ENABLE_MPI
/// helper function to determine the ghost communication flags
CommFlags SPHIntegratorTwoStep::determineFlags(uint64_t timestep)
    {
    auto flags = Integrator::determineFlags(timestep);
    // if (m_rigid_bodies)
    //     {
    //     flags |= m_rigid_bodies->getRequestedCommFlags(timestep);
    //     }
    return flags;
    }
#endif

void SPHIntegratorTwoStep::validateGroups()
    {
    // Check that methods have valid groups.
    size_t group_size = 0;
    for (auto& method : m_methods)
        {
        method->validateGroup();
        group_size += method->getGroup()->getNumMembersGlobal();
        }

    // Check that methods have non-overlapping groups.
    if (m_methods.size() <= 1)
        {
        return;
        }
    auto group_union
        = ParticleGroup::groupUnion(m_methods[0]->getGroup(), m_methods[1]->getGroup());
    for (size_t i = 2; i < m_methods.size(); i++)
        {
        group_union = ParticleGroup::groupUnion(m_methods[i]->getGroup(), group_union);
        }
    if (group_size != group_union->getNumMembersGlobal())
        {
        throw std::runtime_error("Error: the provided groups overlap.");
        }
    }


namespace detail
    {
void export_SPHIntegratorTwoStep(pybind11::module& m)
    {
    pybind11::bind_vector<std::vector<std::shared_ptr<SPHIntegrationMethodTwoStep>>>(
        m,
        "IntegrationMethodList");

    pybind11::class_<SPHIntegratorTwoStep, Integrator, std::shared_ptr<SPHIntegratorTwoStep>>(
        m,
        "SPHIntegratorTwoStep")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, Scalar>())
        .def_property_readonly("methods", &SPHIntegratorTwoStep::getIntegrationMethods)
        // .def_property("rigid", &SPHIntegratorTwoStep::getRigid, &SPHIntegratorTwoStep::setRigid)
        // .def_property("integrate_rotational_dof",
                      // &SPHIntegratorTwoStep::getIntegrateRotationalDOF,
                      // &SPHIntegratorTwoStep::setIntegrateRotationalDOF);
        .def_property("half_step_hook",
                      &SPHIntegratorTwoStep::getHalfStepHook,
                      &SPHIntegratorTwoStep::setHalfStepHook)
        .def("validate_groups", &SPHIntegratorTwoStep::validateGroups);
    }

    } // end namespace detail
    } // end namespace sph
    } // end namespace hoomd

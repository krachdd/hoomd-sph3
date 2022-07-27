// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

#include "hoomd/SystemDefinition.h"
#include "hoomd/ParticleGroup.h"
// #include "hoomd/Profiler.h"


#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#ifndef __SPH_INTEGRATION_METHOD_TWO_STEP_H__
#define __SPH_INTEGRATION_METHOD_TWO_STEP_H__

#ifdef ENABLE_MPI
//! Forward declaration
class Communicator;
#endif

/*! \file SPHIntegrationMethodTwoStep.h
    \brief Declares a base class for all two-step integration methods
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

//! Integrates part of the system forward in two steps
/*! \b See hoomd/nsearch/IntegrationMethodTwoStep.h for documentation.
*/
namespace hoomd
{
namespace sph
{
class SPHIntegrationMethodTwoStep 
    {
    public:
        //! Constructs the integration method and associates it with the system
        SPHIntegrationMethodTwoStep(std::shared_ptr<SystemDefinition> sysdef,
                                 std::shared_ptr<ParticleGroup> group);
        virtual ~SPHIntegrationMethodTwoStep() {};

        //! Abstract method that performs the first step of the integration
        /*! \param timestep Current time step
        */
        virtual void integrateStepOne(unsigned int timestep) {}

        //! Abstract method that performs the second step of the integration
        /*! \param timestep Current time step
        */
        virtual void integrateStepTwo(unsigned int timestep)
            {
            }

        //! Sets the profiler for the integration method to use
        // void setProfiler(std::shared_ptr<Profiler> prof);

        //! Set autotuner parameters
        /*! \param enable Enable/disable autotuning
            \param period period (approximate) in time steps when returning occurs

            Derived classes should override this to set the parameters of their autotuners.
        */
        virtual void setAutotunerParams(bool enable, unsigned int period)
            {
            }

        //! Returns a list of log quantities this compute calculates
        /*! The base class implementation just returns an empty vector. Derived classes should override
            this behavior and return a list of quantities that they log.

            See Logger for more information on what this is about.
        */
        virtual std::vector< std::string > getProvidedLogQuantities()
            {
            return std::vector< std::string >();
            }

        //! Calculates the requested log value and returns it
        /*! \param quantity Name of the log quantity to get
            \param timestep Current time step of the simulation
            \param my_quantity_flag Returns true if this method tracks this quantity

            The base class just returns 0. Derived classes should override this behavior and return
            the calculated value for the given quantity. Only quantities listed in
            the return value getProvidedLogQuantities() will be requested from
            getLogValue().

            See Logger for more information on what this is about.
        */
        virtual Scalar getLogValue(const std::string& quantity, unsigned int timestep,  bool &my_quantity_flag)
            {
            return Scalar(0.0);
            }

        //! Change the timestep
        void setDeltaT(Scalar deltaT);

        //! Access the group
        std::shared_ptr<ParticleGroup> getGroup() { return m_group; }

        //! Get whether this restart was valid
        bool isValidRestart() { return m_valid_restart; }

        //! Get the number of degrees of freedom granted to a given group
        virtual unsigned int getNDOF(std::shared_ptr<ParticleGroup> query_group);

        //! Get needed pdata flags
        /*! Not all fields in ParticleData are computed by default. When derived classes need one of these optional
            fields, they must return the requested fields in getRequestedPDataFlags().
        */
        virtual PDataFlags getRequestedPDataFlags()
            {
            return PDataFlags(0);
            }

        //! Validate that all members in the particle group are valid (throw an exception if they are not)
        virtual void validateGroup();

#ifdef ENABLE_MPI
        //! Set the communicator to use
        /*! \param comm MPI communication class
         */
        void setCommunicator(std::shared_ptr<Communicator> comm)
            {
            assert(comm);
            m_comm = comm;
            }
#endif

    protected:
        const std::shared_ptr<SystemDefinition> m_sysdef; //!< The system definition this method is associated with
        const std::shared_ptr<ParticleGroup> m_group;     //!< The group of particles this method works on
        const std::shared_ptr<ParticleData> m_pdata;      //!< The particle data this method is associated with
        // std::shared_ptr<Profiler> m_prof;                 //!< The profiler this method is to use
        std::shared_ptr<const ExecutionConfiguration> m_exec_conf; //!< Stored shared ptr to the execution configuration
        Scalar m_deltaT;                                    //!< The time step

        //! helper function to get the integrator variables from the particle data
        const IntegratorVariables& getIntegratorVariables()
            {
            return m_sysdef->getIntegratorData()->getIntegratorVariables(m_integrator_id);
            }

        //! helper function to store the integrator variables in the particle data
        void setIntegratorVariables(const IntegratorVariables& variables)
            {
            m_sysdef->getIntegratorData()->setIntegratorVariables(m_integrator_id, variables);
            }

        //! helper function to check if the restart information (if applicable) is useable
        bool restartInfoTestValid(IntegratorVariables& v, std::string type, unsigned int nvariables);

        //! Set whether this restart is valid
        void setValidRestart(bool b) { m_valid_restart = b; }

    protected:
#ifdef ENABLE_MPI
        std::shared_ptr<Communicator> m_comm;             //!< The communicator to use for MPI
#endif
    private:
        unsigned int m_integrator_id;                       //!< Registered integrator id to access the state variables
        bool m_valid_restart;                               //!< True if the restart info was valid when loading
    };

namespace detail
{
//! Exports the SPHIntegrationMethodTwoStep class to python
void export_SPHIntegrationMethodTwoStep(pybind11::module& m);

} // end namespace detail 
} // end namespace sph 
} // end namespace hoomd 

#endif // #ifndef __SPH_INTEGRATION_METHOD_TWO_STEP_H__

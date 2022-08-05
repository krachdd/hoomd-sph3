// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

#include "hoomd/Compute.h"
#include "hoomd/Index1D.h"
#include "hoomd/ParticleGroup.h"
#include "hoomd/ForceCompute.h"
#include "hoomd/Integrator.h"
#include "hoomd/nsearch/NeighborList.h"

#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#include "hoomd/HOOMDMPI.h"
#endif

#include <memory>
#include <vector>
#include <stdexcept>
#include <utility>
#include <set>
#include <algorithm>


#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include "SmoothingKernel.h"
#include "StateEquations.h"
#include "SPHBaseClass.h"
#include "SolidFluidTypeBit.h"


/*! \file SinglePhaseFlow.h
    \brief Contains code for the Quasi-incompressible Navier-Stokes solver
          for Single-phase flow
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __SinglePhaseFlow_H__
#define __SinglePhaseFlow_H__


namespace hoomd 
{
namespace sph
{
//! Enum for indexing the GPUArray of computed values
struct singlephaseflow_logger_index
    {
    enum Enum
        {
        sum_fluid_velocity_x=0, //!< Index for the sum of fluid x-velocity in the GPUArray
        sum_fluid_velocity_y,   //!< Index for the sum of fluid y-velocity in the GPUArray
        sum_fluid_velocity_z,   //!< Index for the sum of fluid z-velocity in the GPUArray
        kinetic_energy,         //!< Index for the overall kinetic energy of the system
        total_fluid_particles,  //!< Total number of fluid particles
        // dt_adapt,               //!< Adaptive timestep size
        num_quantities // final element to count number of quantities
        };
    };


//! Computes SinglePhaseFlow forces on each particle
/*!
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
class PYBIND11_EXPORT SinglePhaseFlow : public SPHBaseClass<KT_, SET_>
    {
    public:

        //! Constructor
        SinglePhaseFlow(std::shared_ptr<SystemDefinition> sysdef,
                        std::shared_ptr<SmoothingKernel<KT_> > skernel,
                        std::shared_ptr<StateEquation<SET_> > equationofstate,
                        std::shared_ptr<nsearch::NeighborList> nlist,
                        std::shared_ptr<ParticleGroup> fluidgroup,
                        std::shared_ptr<ParticleGroup> solidgroup,
                        DensityMethod   mdensitymethod=DENSITYSUMMATION,
                        ViscosityMethod mviscositymethod=HARMONICAVERAGE);

        //! Destructor
        virtual ~SinglePhaseFlow();


        //! Returns sum of fluid x-velocity last computed by compute()
        /*! \returns Instantaneous translational kinetic energy of the system
        */
        Scalar GetSumFluidXVelocity()
            {
            #ifdef ENABLE_MPI
            if (!m_properties_reduced) reduceProperties();
            #endif
            ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
            return h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_x];
            }

        //! Returns sum of fluid y-velocity last computed by compute()
        /*! \returns Sum of fluid y-velocity
        */
        Scalar GetSumFluidYVelocity()
            {
            #ifdef ENABLE_MPI
            if (!m_properties_reduced) reduceProperties();
            #endif
            ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
            return h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_y];
            }

        //! Returns sum of fluid z-velocity last computed by compute()
        /*! \returns Sum of fluid z-velocity
        */
        Scalar GetSumFluidZVelocity()
            {
            #ifdef ENABLE_MPI
            if (!m_properties_reduced) reduceProperties();
            #endif
            ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
            return h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_z];
            }

        Scalar GetFluidParticleNum()
            {
            #ifdef ENABLE_MPI
            if (!m_properties_reduced) reduceProperties();
            #endif
            ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
            return h_properties.data[singlephaseflow_logger_index::total_fluid_particles];
            }
        Scalar GetKineticEnergy()
            {
            #ifdef ENABLE_MPI
            if (!m_properties_reduced) reduceProperties();
            #endif
            ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
            return h_properties.data[singlephaseflow_logger_index::kinetic_energy];
            }
        // Scalar GetAdaptTimestep()
        //     {
        //     #ifdef ENABLE_MPI
        //     if (!m_properties_reduced) reduceProperties();
        //     #endif
        //     ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::read);
        //     return h_properties.data[singlephaseflow_logger_index::dt_adapt];
        //     }

        //! Get the GPU array of properties
        const GPUArray<Scalar>& getProperties()
            {
            #ifdef ENABLE_MPI
            if (!m_properties_reduced) reduceProperties();
            #endif

            return m_properties;
            }

        //! Returns a list of log quantities this compute calculates
        virtual std::vector< std::string > getProvidedLogQuantities();

        //! Calculates the requested log value and returns it
        virtual Scalar getLogValue(const std::string& quantity, unsigned int timestep);

        //! Returns a list of log quantities this compute calculates
        virtual std::vector<double> getProvidedTimestepQuantities(unsigned int timestep);

        /*! Set the parameters
         * \param mu Dynamic viscosity
         */
        virtual void setParams(Scalar mu);

        //! Getter and Setter methods for density method
        DensityMethod getDensityMethod()
            {
            return m_density_method;
            }
        void setDensityMethod(DensityMethod densitymethod)
            {
            m_density_method = densitymethod;
            }

        //! Getter and Setter methods for viscosity method
        ViscosityMethod getViscosityMethod()
            {
            return m_viscosity_method;
            }
        void setViscosityMethod(ViscosityMethod viscositymethod)
            {
            m_viscosity_method = viscositymethod;
            }

        // Set constant smoothing length option to true for faster computation
        void setConstSmoothingLength(Scalar h)
            {
            m_const_slength = true;
            m_ch = h;
            m_rcut = m_kappa * m_ch;
            m_rcutsq = m_rcut * m_rcut;
            }

        /*! Set compute solid forces option to true. This is necessary if suspended object
         *  are present or if solid drag forces are to be evaluated.
         */
        void computeSolidForces()
            {
            m_compute_solid_forces = true;
            }

        /*! Turn Monaghan type artificial viscosity option on.
         * \param alpha Volumetric diffusion coefficient for artificial viscosity operator
         * \param beta Shock diffusion coefficient for artificial viscosity operator
         */
        void activateArtificialViscosity(Scalar alpha, Scalar beta)
            {
            m_artificial_viscosity = true;
            m_avalpha = alpha;
            m_avbeta = beta;
            }

        /*! Turn Monaghan type artificial viscosity option off.
         */
        void deactivateArtificialViscosity()
            {
            m_artificial_viscosity = false;
            }

        /*! Turn Molteni type density diffusion option on.
         * \param ddiff Diffusion coefficient for artificial density diffusion operator
         */
        void activateDensityDiffusion(Scalar ddiff)
            {
            m_density_diffusion = true;
            m_ddiff = ddiff;
            }

        /*! Turn Molteni type density diffusion off.
         */
        void deactivateDensityDiffusion()
            {
            m_density_diffusion = false;
            }

        /*! Turn Shepard type density reinitialization on
         * \param shepardfreq Number of timesteps the renormalization is to be applied
         */
        void activateShepardRenormalization(unsigned int shepardfreq);

        /*! Turn Shepard type density reinitialization off.
         */
        void deactivateShepardRenormalization()
            {
            m_shepard_renormalization = false;
            }

        //! Computes forces
        void computeForces(unsigned int timestep);

    #ifdef ENABLE_MPI
        /// The system's communicator.
        std::shared_ptr<Communicator> m_comm;
    #endif

    #ifdef ENABLE_MPI
        //! Get requested ghost communication flags
        virtual CommFlags getRequestedCommFlags(unsigned int timestep)
            {
            // Request communication of all field required during ForceCompute
            CommFlags flags(0);
            flags[comm_flag::net_force] = 0;
            flags[comm_flag::position] = 1; // Stores position and type
            flags[comm_flag::velocity] = 1; // Stores velocity and mass
            flags[comm_flag::dpe] = 1; // Stores density and pressure
            flags[comm_flag::auxiliary1] = 1; // Stores fictitious velocity
            // Add flags requested by base class
            flags |= ForceCompute::getRequestedCommFlags(timestep);
            return flags;
            }
    #endif

        //! Returns true because we compute dpe array content
        virtual bool ComputesDPE()
            {
            return true;
            }

    protected:

        // Shared pointers
        std::shared_ptr<ParticleGroup> m_fluidgroup; //!< Group of fluid particles
        std::shared_ptr<ParticleGroup> m_solidgroup; //!< Group of fluid particles

        std::shared_ptr<nsearch::NeighborList> m_nlist; //!< the neighborlist to use for the computation

        // Model parameters
        Scalar m_ch; //!< Smoothing length to use if constant for all particles
        Scalar m_rcut; //!< Cut-off length to use if constant for all particles
        Scalar m_rcutsq; //!< Square cut-off length to use if constant for all particles
        DensityMethod m_density_method; //!< Density approach to use
        ViscosityMethod m_viscosity_method; //!< Viscosity approach to use

        // Physical variables
        Scalar m_rho0; //!< Rest density (Read from equation of state class)
        Scalar m_c; //!< Speed of sound (Read from equation of state class)
        Scalar m_kappa; //!< Kernel scaling factor (Read from kernel class)
        Scalar m_mu; //!< Viscosity ( Must be set by user )
        Scalar m_avalpha; //!< Volumetric diffusion coefficient for artificial viscosity operator
        Scalar m_avbeta; //!< Shock diffusion coefficient for artificial viscosity operator
        Scalar m_ddiff; //!< Diffusion coefficient for Molteni type density diffusion
        unsigned int m_shepardfreq; //!< Time step frequency for Shepard reinitialization

        // Auxiliary variables
        std::vector<unsigned int> m_fluidtypes; //!< Fluid type numbers
        std::vector<unsigned int> m_solidtypes; //!< Solid type numbers
        GPUArray<unsigned int> m_type_property_map; //!< to check if a particle type is solid or fluid

        // Flags
        bool m_const_slength; //!< True if using constant smoothing length
        bool m_compute_solid_forces; //!< Set to true if forces acting on solid particle are to be computed
        bool m_artificial_viscosity; //!< Set to true if Monaghan type artificial viscosity is to be used
        bool m_density_diffusion; //!< Set to true if Molteni type density diffusion is to be used
        bool m_shepard_renormalization; //!< Set to true if Shepard type density reinitialization is to be used
        bool m_params_set; //!< True if parameters are set

        // Log parameters
        std::vector<std::string> m_logname_list;  //!< Cache all generated logged quantities names
        GPUArray<Scalar> m_properties; //!< Stores the computed properties
        unsigned int m_log_computed_last_timestep; //!< Last time step where log quantities were computed

        // Timestep parameters
        std::vector<double> m_timestep_list = std::vector<double>(7);  //!< Cache all generated timestep quantities names

        //! Computes log properties
        void computeProperties();

#ifdef ENABLE_MPI
        bool m_properties_reduced;      //!< True if properties have been reduced across MPI
        //! Reduce properties over MPI
        void reduceProperties();
#endif

        /*! Helper function to compute particle number density
         * \post For fluid particles, compute number density. For solid particles,
                 compute fluid normalization constant.
         */
        void compute_ndensity(unsigned int timestep);

        /*! Helper function to compute particle pressures
         *  \post Pressure of fluid particle computed
         */
        void compute_pressure(unsigned int timestep);

        /*! Helper function to compute fictitious solid particle properties (pressures and velocities)
        * \pre Ghost particle number densities (i.e. dpe array) must be up-to-date
        * \pre Solid normalization constant \sum_j w_ij must be computed and stored in dpe array
        * \post Fictitious particle properties are computed and stored in aux1 array
        */
        void compute_noslip(unsigned int timestep);

        /*! Helper function to apply Shepard density filter
        * \post Fluid particle densities are recomputed based on the Shepard renormalization
        */
        void renormalize_density(unsigned int timestep);

        /*! Helper function where the actual force computation takes place
         * \pre Number densities and fictitious solid particle properties must be up-to-date
         * \post h_force stores forces acting on fluid particles and .w component stores rate of change of density
         */
        void forcecomputation(unsigned int timestep);

        /*! Helper function to set communication flags and update ghosts densities
        * \param timestep The time step
        * \post Ghost particle dpe array is up-to-date
        */
        void update_ghost_dpe(unsigned int timestep);

        /*! Helper function to set communication flags and update ghosts auxiliary array 1
        * \param timestep The time step
        * \post Ghost particle auxiliary array 1 is up-to-date
        */
        void update_ghost_aux1(unsigned int timestep);

    private:

    };


namespace detail 
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlow(pybind11::module& m, std::string name);
}

// namespace detail
// {
// //! Exports the SinglePhaseFlow class to python
// // template<SmoothingKernelType KT_,StateEquationType SET_>
// // void export_SinglePhaseFlow_templ();
// // template<SmoothingKernelType KT_,StateEquationType SET_>
// void export_SinglePhaseFlow(pybind11::module& m);

// } // end namespace detail
} // end namespace sph
} // end namespace hoomd

#endif // __SinglePhaseFlow_H__

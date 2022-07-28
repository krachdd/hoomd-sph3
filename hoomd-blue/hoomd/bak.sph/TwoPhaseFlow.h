// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

#include "hoomd/Compute.h"
#include "hoomd/Index1D.h"
#include "hoomd/ParticleGroup.h"
#include "hoomd/ForceCompute.h"
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
#include <pybind11/numpy.h>

#include "SmoothingKernel.h"
#include "StateEquations.h"
#include "SPHBaseClass.h"
#include "SolidFluidTypeBit.h"

/*! \file TwoPhaseFlow.h
    \brief Contains code for the quasi-incompressible two phase flow including capillarity.
           This model takes as input two fluid groups (1,2) and one solid groups.
           For each fluid group there must be a viscosity parameter and eos instance passed as parameters.
           Furthermore, fluid interfacial tension as well as solid contact angle are parameters.
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __TwoPhaseFlow_H__
#define __TwoPhaseFlow_H__

//! Computes TwoPhaseFlow forces on each particle
/*!
*/        
namespace hoomd 
{
namespace sph
{

template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
class TwoPhaseFlow : public SPHBaseClass<KT_,SET1_>
    {
    public:


        //! Constructor
        TwoPhaseFlow(std::shared_ptr<SystemDefinition> sysdef,
                     std::shared_ptr<SmoothingKernel<KT_> > skernel,
                     std::shared_ptr<StateEquation<SET1_> > equationofstate1,
                     std::shared_ptr<StateEquation<SET2_> > equationofstate2,
                     std::shared_ptr<nsearch::NeighborList> nlist,
                     std::shared_ptr<ParticleGroup> fluidgroup1,
                     std::shared_ptr<ParticleGroup> fluidgroup2,
                     std::shared_ptr<ParticleGroup> solidgroup,
                     DensityMethod   mdensitymethod=DENSITYSUMMATION,
                     ViscosityMethod mviscositymethod=HARMONICAVERAGE);

        //! Destructor
        virtual ~TwoPhaseFlow();

        /*! Set the parameters
         * \param mu1 Dynamic viscosity of fluid 1
         * \param mu2 Dynamic viscosity of fluid 2
         * \param sigma12 Fluid interfacial tension
         * \param omega Solid - Fluid 1 contact angle
         */
        virtual void setParams(Scalar mu1, Scalar mu2, Scalar sigma12, Scalar omega);

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
        void activateShepardRenormalization(unsigned int shepardfreq)
            {
            m_shepard_renormalization = true;
            m_shepardfreq = shepardfreq;
            }

        /*! Turn Shepard type density reinitialization off.
         */
        void deactivateShepardRenormalization()
            {
            m_shepard_renormalization = false;
            }

        //! Computes forces
        void computeForces(unsigned int timestep);

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
            flags[comm_flag::auxiliary2] = 1; // Stores solid normal vector field
            flags[comm_flag::auxiliary3] = 1; // Stores fluid interfacial normal vector field
            flags[comm_flag::auxiliary4] = 1; // Stores surface force density
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

        // Log parameters
        std::vector<std::string> m_logname_list;  //!< Cache all generated logged quantities names
        GPUArray<Scalar> m_properties;  //!< Stores the computed properties

        // Shared pointers
        std::shared_ptr<ParticleGroup> m_fluidgroup1; //!< Group of fluid particles 1
        std::shared_ptr<ParticleGroup> m_fluidgroup2; //!< Group of fluid particles 2
        std::shared_ptr<ParticleGroup> m_fluidgroup; //!< Group of fluid particles (union of fluid 1 + 2)
        std::shared_ptr<ParticleGroup> m_solidgroup; //!< Group of fluid particles
        std::shared_ptr<StateEquation<SET1_> > m_eos1; //!< The equation of state class for fluid phase 1
        std::shared_ptr<StateEquation<SET2_> > m_eos2; //!< The equation of state class for fluid phase 2
        GPUArray<unsigned int> m_type_property_map; //!< to check if a particle type is solid or fluid

        // Model parameters
        Scalar m_ch; //!< Smoothing length to use if constant for all particles
        Scalar m_rcut; //!< Cut-off length to use if constant for all particles
        Scalar m_rcutsq; //!< Square cut-off length to use if constant for all particles
        DensityMethod m_density_method; //!< Density approach to use
        ViscosityMethod m_viscosity_method; //!< Viscosity approach to use

        // Physical variables
        Scalar m_rho01; //!< Rest density of fluid 1 (Read from equation of state class)
        Scalar m_rho02; //!< Rest density of fluid 2 (Read from equation of state class)
        Scalar m_c1; //!< Speed of sound of fluid 1 (Read from equation of state class)
        Scalar m_c2; //!< Speed of sound of fluid 2 (Read from equation of state class)
        Scalar m_cmax; //!< Maximum speed of sound
        Scalar m_kappa; //!< Kernel scaling factor (Read from kernel class)

        Scalar m_mu1; //!< Viscosity ( Must be set by user )
        Scalar m_mu2; //!< Viscosity ( Must be set by user )
        Scalar m_sigma12; //!< Interfacial tension between fluid phases ( Must be set by user )
        Scalar m_sigma01; //!< Interfacial tension between solid phase and fluid phase 1 ( Computed form input )
        Scalar m_sigma02; //!< Interfacial tension between solid phase and fluid phase 2 ( Computed form input )
        Scalar m_omega; //!< Contact angle ( Must be set by user )

        Scalar m_avalpha; //!< Volumetric diffusion coefficient for artificial viscosity operator
        Scalar m_avbeta; //!< Shock diffusion coefficient for artificial viscosity operator
        Scalar m_ddiff; //!< Diffusion coefficient for Molteni type density diffusion
        unsigned int m_shepardfreq; //!< Time step frequency for Shepard reinitialization

        // Auxiliary variables
        std::vector<unsigned int> m_fluidtypes1; //!< Fluid 1 type numbers
        std::vector<unsigned int> m_fluidtypes2; //!< Fluid 2 type numbers
        std::vector<unsigned int> m_fluidtypes; //!< Fluid type numbers
        std::vector<unsigned int> m_solidtypes; //!< Solid type numbers

        // Flags
        bool m_const_slength; //!< True if using constant smoothing length
        bool m_compute_solid_forces; //!< Set to true if forces acting on solid particle are to be computed
        bool m_artificial_viscosity; //!< Set to true if Monaghan type artificial viscosity is to be used
        bool m_density_diffusion; //!< Set to true if Molteni type density diffusion is to be used
        bool m_shepard_renormalization; //!< Set to true if Shepard type density reinitialization is to be used
        bool m_params_set; //!< True if parameters are set

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

        /*! Helper function to compute solid-fluid and fluid-fluid color gradient vectors
        * \post Solid color gradient vectors are stored in aux2
        * \post Fluid color gradient vectors are stored in aux3
        */
        void compute_colorgradients(unsigned int timestep);

        /*! Helper function to compute interfacial surface force field
        * \pre Normal vector field have been computed and communicated
        * \post Surface force density vectors are stored in aux4
        */
        void compute_surfaceforce(unsigned int timestep);

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

        /*! Helper function to set communication flags and update ghosts auxiliary arrays
        * \param timestep The time step
        * \post Ghost particle auxiliary array 1 is up-to-date
        * \post Ghost particle auxiliary array 2 is up-to-date
        * \post Ghost particle auxiliary array 3 is up-to-date
        */
        void update_ghost_aux123(unsigned int timestep);

        /*! Helper function to set communication flags and update ghosts auxiliary array 4
        * \param timestep The time step
        * \post Ghost particle auxiliary array 4 is up-to-date
        */
        void update_ghost_aux4(unsigned int timestep);

    private:

    };

namespace detail
{
//! Exports the TwoPhaseFlow class to python
// template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
// void export_TwoPhaseFlow_templ();

template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void export_TwoPhaseFlow(pybind11::module& m);

} // end namespace detail 
} // end namespace sph 
} // end namespace hoomd 

#endif // __TwoPhaseFlow_H__

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

#include "EvaluationMethodDefinition.h"

/*! \file TwoPhaseFlow.h
    \brief Contains code for the Quasi-incompressible Navier-Stokes solver
          for Two-phase flow
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __TwoPhaseFlow_H__
#define __TwoPhaseFlow_H__

namespace hoomd 
{
namespace sph
{

//! Computes TwoPhaseFlow forces on each particle
/*!
*/

template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
class PYBIND11_EXPORT TwoPhaseFlow : public SPHBaseClass<KT_, SET1_>
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
                        ViscosityMethod mviscositymethod=HARMONICAVERAGE,
                        ColorGradientMethod mcolorgradientmethod=DENSITYRATIO);

        //! Destructor
        virtual ~TwoPhaseFlow();

        //! Set the rcut for a single type pair
        virtual void setRcut(unsigned int typ1, unsigned int typ2, Scalar rcut);

        /// Set the rcut for a single type pair using a tuple of strings
        virtual void setRCutPython(pybind11::tuple types, Scalar r_cut);

        /// Validate that types are within Ntypes
        void validateTypes(unsigned int typ1, unsigned int typ2, std::string action);

        /*! Set the parameters
         * \param mu1 Dynamic viscosity
         * \param mu2 Dynamic viscosity
         * \param sigma12 Fluid interfacial tension
         * \param omega Solid - Fluid 1 contact angle
         */
        virtual void setParams(Scalar mu1, Scalar mu2, Scalar sigma12, Scalar omega);

        virtual void setHysteresis(Scalar omega_rec, Scalar omega_adv);

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

        //! Getter and Setter methods for viscosity method
        ColorGradientMethod getColorGradientMethod()
            {
            return m_colorgradient_method;
            }
        void setColorGradientMethod(ColorGradientMethod colorgradientmethod)
            {
            m_colorgradient_method = colorgradientmethod;
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

        /*! Turn Monaghan artificial viscosity (AV) on.
         *
         *  Adds a viscous dissipation term to the pressure force for approaching
         *  fluid particle pairs (\f$v_{ij} \cdot r_{ij} < 0\f$) only.  The Monaghan 1992 form:
         *
         *    \f$\Pi_{ij} = (-\alpha c_\mathrm{max} \mu_{ij} + \beta \mu_{ij}^2) / \bar{\rho}_{ij}\f$
         *
         *  where  \f$\mu_{ij} = \bar{h} (v_{ij} \cdot r_{ij}) / (r_{ij}^2 + \eta^2)\f$, \f$\eta = 0.1\bar{h}\f$,
         *         \f$c_\mathrm{max}\f$ is the global maximum speed of sound (from setParams),
         *         \f$\bar{\rho}_{ij} = (\rho_i + \rho_j) / 2\f$.
         *
         *  The linear term (\f$\alpha\f$) diffuses velocity divergence; the quadratic term (\f$\beta\f$)
         *  is a Rankine–Hugoniot correction active only in shock-like conditions.
         *  For most weakly-compressible flows \f$\beta = 0\f$ and \f$\alpha \in [0.01, 0.1]\f$ suffices.
         *  Applied only to fluid–fluid pairs.  Mutually exclusive with Riemann
         *  dissipation (activateRiemannDissipation).
         *
         *  Reference: Monaghan (1992) Annu. Rev. Astron. Astrophys. 30, 543–574.
         *
         * \param alpha  Linear (volumetric) diffusion coefficient
         * \param beta   Quadratic (shock) diffusion coefficient
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

        // ── Non-Newtonian rheology for fluid 1 ────────────────────────────────
        void activatePowerLaw1(Scalar K, Scalar n, Scalar mu_min = Scalar(0))
            {
            m_nn_model1  = POWERLAW;
            m_nn_K1      = K;
            m_nn_n1      = n;
            m_nn_mu_min1 = mu_min;
            }
        void activateCarreau1(Scalar mu0, Scalar muinf, Scalar lambda_NN, Scalar n)
            {
            m_nn_model1   = CARREAU;
            m_nn_mu0_1    = mu0;
            m_nn_muinf_1  = muinf;
            m_nn_lambda1  = lambda_NN;
            m_nn_n1       = n;
            }
        void activateBingham1(Scalar mu_p, Scalar tauy, Scalar m_reg, Scalar mu_min = Scalar(0))
            {
            m_nn_model1  = BINGHAM;
            m_nn_K1      = mu_p;
            m_nn_tauy1   = tauy;
            m_nn_m1      = m_reg;
            m_nn_mu_min1 = mu_min;
            }
        void activateHerschelBulkley1(Scalar K, Scalar n, Scalar tauy, Scalar m_reg, Scalar mu_min = Scalar(0))
            {
            m_nn_model1  = HERSCHELBULKLEY;
            m_nn_K1      = K;
            m_nn_n1      = n;
            m_nn_tauy1   = tauy;
            m_nn_m1      = m_reg;
            m_nn_mu_min1 = mu_min;
            }
        void deactivateNonNewtonian1()
            {
            m_nn_model1 = NEWTONIAN;
            }

        // ── Non-Newtonian rheology for fluid 2 ────────────────────────────────
        void activatePowerLaw2(Scalar K, Scalar n, Scalar mu_min = Scalar(0))
            {
            m_nn_model2  = POWERLAW;
            m_nn_K2      = K;
            m_nn_n2      = n;
            m_nn_mu_min2 = mu_min;
            }
        void activateCarreau2(Scalar mu0, Scalar muinf, Scalar lambda_NN, Scalar n)
            {
            m_nn_model2   = CARREAU;
            m_nn_mu0_2    = mu0;
            m_nn_muinf_2  = muinf;
            m_nn_lambda2  = lambda_NN;
            m_nn_n2       = n;
            }
        void activateBingham2(Scalar mu_p, Scalar tauy, Scalar m_reg, Scalar mu_min = Scalar(0))
            {
            m_nn_model2  = BINGHAM;
            m_nn_K2      = mu_p;
            m_nn_tauy2   = tauy;
            m_nn_m2      = m_reg;
            m_nn_mu_min2 = mu_min;
            }
        void activateHerschelBulkley2(Scalar K, Scalar n, Scalar tauy, Scalar m_reg, Scalar mu_min = Scalar(0))
            {
            m_nn_model2  = HERSCHELBULKLEY;
            m_nn_K2      = K;
            m_nn_n2      = n;
            m_nn_tauy2   = tauy;
            m_nn_m2      = m_reg;
            m_nn_mu_min2 = mu_min;
            }
        void deactivateNonNewtonian2()
            {
            m_nn_model2 = NEWTONIAN;
            }

        /*! Turn consistent interface pressure (CIP) on.
         *
         *  For DENSITYSUMMATION and cross-phase fluid pairs (fluid 1 \f$\leftrightarrow\f$ fluid 2),
         *  replaces the standard actual-density-weighted pressure average
         *
         *    \f$\bar{p}_{ij} = (\rho_j p_i + \rho_i p_j) / (\rho_i + \rho_j)\f$            [standard]
         *
         *  with a rest-density-weighted form plus a hydrostatic correction:
         *
         *    \f$\bar{p}_{ij} = (\rho_{0j} p_i + \rho_{0i} p_j + \rho_{0i} \rho_{0j} (g \cdot r_{ij})) / (\rho_{0i} + \rho_{0j})\f$
         *
         *  where \f$g\f$ = getAcceleration(timestep) is the body-force acceleration and
         *  \f$r_{ij} = r_i - r_j\f$.
         *
         *  Physical motivation: in hydrostatic equilibrium the pressure gradient is
         *  \f$\nabla p = \rho g\f$, which is discontinuous at a density-ratio interface.  The
         *  standard density-weighted average introduces a systematic error proportional
         *  to \f$(\rho_1 - \rho_2)\f$ that drives spurious "parasitic" interfacial velocities.
         *  Rest-density weighting removes this bias; the \f$g \cdot r_{ij}\f$ term ensures the
         *  discretised pressure gradient exactly reproduces gravity for a static column,
         *  regardless of density ratio.  Especially beneficial for water/air (\f$\rho\f$ ratio
         *  ~ 1000:1) cases.
         *
         *  Same-phase pairs, solid-boundary interactions, and DENSITYCONTINUITY are
         *  unaffected.  May be combined with Riemann dissipation.
         *
         *  Reference: Hu & Adams (2009) J. Comput. Phys. 228(20), 7518–7530.
         *             Adami, Hu & Adams (2012) J. Comput. Phys. 231(21), 7057–7075.
         */
        void activateConsistentInterfacePressure()
            {
            m_consistent_interface_pressure = true;
            }

        /*! Turn consistent interface pressure off.
         */
        void deactivateConsistentInterfacePressure()
            {
            m_consistent_interface_pressure = false;
            }

        /*! Turn Riemann-based dissipation on.
         *
         *  Replaces Monaghan AV with a physically-motivated dissipation derived from
         *  the linearised inter-particle Riemann problem.  Uses the harmonic mean of
         *  acoustic impedances \f$Z = \rho c\f$ to handle density and speed-of-sound contrasts
         *  across the interface naturally:
         *
         *    \f$Z^*_{ij}  =  Z_i \cdot Z_j / (Z_i + Z_j)\f$,   \f$Z = \rho c\f$        [harmonic impedance]
         *    \f$u_{ij}   =  (v_i - v_j) \cdot (r_i - r_j) / (|r_{ij}| + \eta)\f$  [signed radial vel.]
         *    \f$\mathrm{avc}    =  -\beta_R \cdot Z^*_{ij} \cdot u_{ij}^- / \bar{\rho}_{ij}\f$               (\f$u^- = \min(u, 0)\f$)
         *
         *  where \f$\bar{\rho}_{ij} = (\rho_i + \rho_j) / 2\f$.  Applied only to approaching fluid–fluid
         *  pairs (\f$v_{ij} \cdot r_{ij} < 0\f$).  Mutually exclusive with Monaghan AV.
         *
         *  Advantages over Monaghan AV for two-phase flows:
         *    - No \f$\alpha\f$/\f$\beta\f$ parameters to tune per problem (only one global \f$\beta_R \approx 1\f$).
         *    - Impedance mismatch at the interface is handled automatically:
         *      \f$Z^* \rightarrow Z_\mathrm{light}/2\f$ when \f$Z_\mathrm{heavy} \gg Z_\mathrm{light}\f$ (e.g. water/air).
         *    - Low dissipation in smooth regions; increases only near shocks/interfaces.
         *
         *  Reference: Zhang, Hu & Adams (2017) J. Comput. Phys. 340, 439–455.
         *
         * \param beta  Dissipation scaling coefficient (default 1.0; reduce towards
         *              0.5 for smoother flows if over-damping is observed)
         */
        void activateRiemannDissipation(Scalar beta = Scalar(1.0))
            {
            m_riemann_dissipation = true;
            m_riemann_beta = beta;
            }

        /*! Turn Riemann-based dissipation off.
         */
        void deactivateRiemannDissipation()
            {
            m_riemann_dissipation = false;
            }

        /*! Turn Molteni–Colagrossi density diffusion on (DENSITYCONTINUITY only).
         *
         *  Adds a diffusive correction to \f$\mathrm{d}\rho/\mathrm{d}t\f$ to smooth density oscillations:
         *
         *    \f$\mathrm{d}\rho_i/\mathrm{d}t \mathrel{+}= -2\delta\bar{h} c_\mathrm{max} m_j (\rho_i/\rho_{0i} - \rho_j/\rho_{0j}) (r_{ij} \cdot \nabla W_{ij}) / (r^2+\eta^2)\f$
         *
         *  The drive term \f$(\rho_i/\rho_{0i} - \rho_j/\rho_{0j})\f$ is the rest-density-normalised form
         *  rather than the original \f$(\rho_i/\rho_j - 1)\f$.  The original is non-zero at
         *  equilibrium when the two phases have different rest densities (\f$\rho_{01} \neq \rho_{02}\f$),
         *  causing unphysical density drift across the interface.  The normalised form
         *  equals zero at equilibrium for both single- and two-phase flows, correcting
         *  stratified-flow artefacts.
         *
         *  Reference: Molteni & Colagrossi (2009) Comput. Phys. Commun. 180, 861–872.
         *             Two-phase correction: see also Grenier et al. (2013).
         *
         * \param ddiff  Diffusion coefficient \f$\delta\f$ (typically 0.1; range 0.05–0.2)
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

        /*! Turn periodic density resummation on (DENSITYCONTINUITY only)
         * \param densreinitfreq Frequency (in timesteps) of density reinitialization
         */
        void activateDensityReinitialization(unsigned int densreinitfreq);

        /*! Turn periodic density resummation off.
         */
        void deactivateDensityReinitialization()
            {
            m_density_reinitialization = false;
            }

        //! Returns a list of log quantities this compute calculates
        virtual std::vector<double> getProvidedTimestepQuantities(uint64_t timestep);

        /*! Turn Fickian shifting based on particle concentration on
         * \param Used in Computation of CSF
         */
        void activateFickianShifting()
            {
            m_fickian_shifting = true;
            }

        /*! Turn Shepard type density reinitialization off.
         */
        void deactivateFickianShifting()
            {
            m_fickian_shifting = false;
            }

        /*! Activate \f$\delta^+\f$-SPH particle shifting (Sun et al. 2017, Comput. Fluids).
         * Shifts fluid particle positions each step to maintain regularity and
         * prevent clustering. The interface-normal component is projected out so
         * that particles cannot cross the fluid-fluid interface.
         * \param A  Amplitude (default 0.2). Start small (~0.05) for first tests.
         * \param R  Enhancement coefficient (default 0.2, Sun et al. recommended).
         * \param n  Enhancement exponent (default 4, Sun et al. recommended).
         * \param interface_condition  If true, remove normal-shift at interface (recommended).
         */
        void activateParticleShifting(Scalar A = Scalar(0.2),
                                       Scalar R = Scalar(0.2),
                                       unsigned int n = 4,
                                       bool interface_condition = true);

        /*! Turn \f$\delta^+\f$-SPH particle shifting off.
         */
        void deactivateParticleShifting()
            {
            m_particle_shifting = false;
            }

        //! Computes forces
        void computeForces(uint64_t timestep);

    #ifdef ENABLE_MPI
        //! Get requested ghost communication flags
        virtual CommFlags getRequestedCommFlags(uint64_t timestep)
            {
            // Request communication of all field required during ForceCompute
            CommFlags flags(0);
            flags[comm_flag::net_force] = 0;
            flags[comm_flag::position] = 1; // Stores position and type
            flags[comm_flag::velocity] = 1; // Stores velocity and mass
            flags[comm_flag::density] = 1; // Stores density 
            flags[comm_flag::pressure] = 1; // Stores pressure
            flags[comm_flag::energy] = 0; // Stores energy/ Partcile Concentration gradient
            flags[comm_flag::auxiliary1] = 1; // Stores fictitious velocity
            flags[comm_flag::auxiliary2] = 1; // Stores solid normal vector field
            flags[comm_flag::auxiliary3] = 1; // Stores fluid interfacial normal vector
            flags[comm_flag::auxiliary4] = 1; // Stores surface force density
            flags[comm_flag::slength] = 1; // Stores smoothing length TODO is this needed
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
    #ifdef ENABLE_MPI
        std::shared_ptr<Communicator> m_comm; //!< The system's MPI communicator
    #endif

        // Shared pointers
        std::shared_ptr<ParticleGroup> m_fluidgroup; //!< Group of fluid particles (union of fluid 1 + 2)
        std::shared_ptr<ParticleGroup> m_fluidgroup1; //!< Group of fluid particles
        std::shared_ptr<ParticleGroup> m_fluidgroup2; //!< Group of fluid particles
        std::shared_ptr<ParticleGroup> m_solidgroup; //!< Group of solid particles
        std::shared_ptr<StateEquation<SET1_>> m_eos1; //!< The equation of state class for fluid phase 1
        std::shared_ptr<StateEquation<SET2_>> m_eos2; //!< The equation of state class for fluid phase 2

        /// r_cut (not squared) given to the neighbor list
        std::shared_ptr<GPUArray<Scalar>> m_r_cut_nlist;

        // Index for rcut pair info -> nlist
        Index2D m_typpair_idx;        //!< Helper class for indexing per type pair arrays

        // Model parameters
        Scalar m_ch; //!< Smoothing length to use if constant for all particles
        Scalar m_rcut; //!< Cut-off length to use if constant for all particles
        Scalar m_rcutsq; //!< Square cut-off length to use if constant for all particles
        DensityMethod m_density_method; //!< Density approach to use
        ViscosityMethod m_viscosity_method; //!< Viscosity approach to use
        ColorGradientMethod m_colorgradient_method; //!< Colorgradient approach to use

        // Physical variables
        Scalar m_rho01; //!< Rest density (Read from equation of state class)
        Scalar m_rho02; //!< Rest density (Read from equation of state class)
        Scalar m_c1; //!< Speed of sound (Read from equation of state class)
        Scalar m_c2; //!< Speed of sound (Read from equation of state class)
        Scalar m_cmax; //!< Maximum Speed of sound 
        Scalar m_kappa; //!< Kernel scaling factor (Read from kernel class)
        Scalar m_mu1; //!< Viscosity ( Must be set by user )
        Scalar m_mu2; //!< Viscosity ( Must be set by user )
        Scalar m_sigma12; //!< Interfacial tension between fluid phases ( Must be set by user )
        Scalar m_sigma01; //!< Interfacial tension between solid phase and fluid phase1 ( Computed from input )
        Scalar m_sigma02; //!< Interfacial tension between solid phase and fluid phase2 ( Computed from input )
        Scalar m_omega;     //!< Contact angle ( Must be set by user )
        Scalar m_omega_adv; //!< Advancing contact angle for hysteresis [deg]
        Scalar m_omega_rec; //!< Receding  contact angle for hysteresis [deg]
        bool   m_hysteresis; //!< True if contact-angle hysteresis is active

        Scalar m_avalpha; //!< Monaghan AV: linear (volumetric) diffusion coefficient \f$\alpha\f$
        Scalar m_avbeta;  //!< Monaghan AV: quadratic (shock) diffusion coefficient \f$\beta\f$
        Scalar m_riemann_beta; //!< Riemann dissipation: scaling coefficient \f$\beta_R\f$ (Zhang et al. 2017)
        Scalar m_ddiff; //!< Diffusion coefficient for Molteni type density diffusion
        unsigned int m_shepardfreq; //!< Time step frequency for Shepard reinitialization

        // Non-Newtonian rheology for fluid 1
        NonNewtonianModel m_nn_model1;
        Scalar m_nn_K1, m_nn_n1;
        Scalar m_nn_mu0_1, m_nn_muinf_1, m_nn_lambda1;
        Scalar m_nn_tauy1, m_nn_m1, m_nn_mu_min1;

        // Non-Newtonian rheology for fluid 2
        NonNewtonianModel m_nn_model2;
        Scalar m_nn_K2, m_nn_n2;
        Scalar m_nn_mu0_2, m_nn_muinf_2, m_nn_lambda2;
        Scalar m_nn_tauy2, m_nn_m2, m_nn_mu_min2;

        // Auxiliary variables
        std::vector<unsigned int> m_fluidtypes1; //!< Fluid 1 type numbers
        std::vector<unsigned int> m_fluidtypes2; //!< Fluid 2 type numbers
        std::vector<unsigned int> m_fluidtypes; //!< Fluid type numbers
        std::vector<unsigned int> m_solidtypes; //!< Solid type numbers
        GPUArray<unsigned int> m_type_property_map; //!< to check if a particle type is solid or fluid

        // Flags
        bool m_const_slength; //!< True if using constant smoothing length
        bool m_compute_solid_forces; //!< Set to true if forces acting on solid particle are to be computed
        bool m_artificial_viscosity; //!< True if Monaghan (1992) AV is active (mutually exclusive with Riemann)
        bool m_riemann_dissipation;  //!< True if Riemann-based dissipation is active (Zhang et al. 2017; mutually exclusive with Monaghan AV)
        bool m_consistent_interface_pressure; //!< True if Hu & Adams (2009) rest-density-weighted cross-phase pressure is active
        bool m_density_diffusion; //!< Set to true if Molteni type density diffusion is to be used
        bool m_shepard_renormalization; //!< Set to true if Shepard type density reinitialization is to be used
        bool m_params_set; //!< True if parameters are set
        bool m_solid_removed; //!< True if solid Particles have been marked to remove
        bool m_fickian_shifting; //!< True if Fickian Particle Shifting is activated
        bool m_pressure_initialized; //!< True once pressure has been EOS-initialized (DENSITYCONTINUITY only)
        bool m_density_reinitialization; //!< True if periodic density resummation is activated
        unsigned int m_densityreinitfreq; //!< Frequency for density reinitialization

        // Particle shifting (\f$\delta^+\f$-SPH, Sun et al. 2017)
        bool m_particle_shifting;            //!< True if \f$\delta^+\f$-SPH particle shifting is activated
        Scalar m_shift_A;                    //!< Shifting amplitude A (default 0.2)
        Scalar m_shift_R;                    //!< Enhancement coefficient R (default 0.2)
        unsigned int m_shift_n;              //!< Enhancement exponent n (default 4)
        bool m_shift_interface_condition;    //!< True to project out interface-normal shift component

        // Log parameters
        uint64_t m_log_computed_last_timestep; //!< Last time step where log quantities were computed

        // Timestep parameters
        std::vector<double> m_timestep_list = std::vector<double>(8);  //!< Cache all generated timestep quantities names


        void mark_solid_particles_toremove(uint64_t timestep);

        /*! Helper function to compute particle number density
         * \post For fluid particles, compute number density. For solid particles,
                 compute fluid normalization constant.
         */
        virtual void compute_ndensity(uint64_t timestep);

        void compute_particlenumberdensity(uint64_t timestep);

        /*! Helper function to compute particle pressures
         *  \post Pressure of fluid particle computed
         */
        virtual void compute_pressure(uint64_t timestep);

        /*! Helper function to compute fictitious solid particle properties (pressures and velocities)
        * \pre Ghost particle number densities (i.e. density array) must be up-to-date
        * \pre Solid normalization constant \sum_j w_ij must be computed and stored in density array
        * \post Fictitious particle properties are computed and stored in aux1 array
        */
        virtual void compute_noslip(uint64_t timestep);

        /*! Helper function to compute particle concentration gradient for Fickian shifting
         * within the CSF computation
         * It overwrites h_pressure and the dot product of the gradient is stored in h_energy
         * Paper: Lind et al. 2012
        */
        void compute_particle_concentration_gradient(uint64_t timestep);

        /*! Compute and apply \f$\delta^+\f$-SPH particle position shifts (Sun et al. 2017).
         * \pre compute_colorgradients() + update_ghost_aux123() must precede this call
         *      (aux3 must hold up-to-date fluid-fluid interface normals).
         * \post Fluid particle positions updated by \f$\delta r_i\f$.
         * \post For DENSITYCONTINUITY: h_density corrected by ALE remapping term.
         */
        void compute_particle_shift(uint64_t timestep);

        /*! Helper function to apply Shepard density filter
        * \post Fluid particle densities are recomputed based on the Shepard renormalization
        */
        void renormalize_density(uint64_t timestep);

        /*! Helper function to compute solid-fluid and fluid-fluid color gradient vectors
        * \post Solid color gradient vectors are stored in aux2
        * \post Fluid color gradient vectors are stored in aux3
        */
        virtual void compute_colorgradients(uint64_t timestep);

        /*! Helper function to compute interfacial surface force field
        * \pre Normal vector field have been computed and communicated
        * \post Surface force density vectors are stored in aux4
        */
        virtual void compute_surfaceforce(uint64_t timestep);

        /*! Helper function where the actual force computation takes place
         * \pre Number densities and fictitious solid particle properties must be up-to-date
         * \post h_force stores forces acting on fluid particles and .w component stores rate of change of density
         */
        virtual void forcecomputation(uint64_t timestep);

        /*! Helper function that computes the fluid-induced forces on solid particles
         */
        virtual void compute_solid_forces(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts densities
        * \param timestep The time step
        * \post Ghost particle density array is up-to-date
        */
        void update_ghost_density(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts densities and pressures
        * \param timestep The time step
        * \post Ghost particle density and pressue array is up-to-date
        */
        void update_ghost_density_pressure(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts auxiliary arrays
        * \param timestep The time step
        * \post Ghost particle auxiliary array 1 is up-to-date
        * \post Ghost particle auxiliary array 2 is up-to-date
        * \post Ghost particle auxiliary array 3 is up-to-date
        */
        void update_ghost_aux123(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts auxiliary array 4
        * \param timestep The time step
        * \post Ghost particle auxiliary array 4 is up-to-date
        */
        void update_ghost_aux4(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts auxiliary array 4
        * \param timestep The time step
        * \post Ghost particle density, pressure and energy is up-to-date
        */
        void update_ghost_density_pressure_energy(uint64_t timestep);

    private:

    };


namespace detail 
{
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void export_TwoPhaseFlow(pybind11::module& m, std::string name);

} // end namespace detail
} // end namespace sph
} // end namespace hoomd

#endif // __TwoPhaseFlow_H__

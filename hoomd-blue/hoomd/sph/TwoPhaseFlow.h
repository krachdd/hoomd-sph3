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
         *  fluid particle pairs (v_ij · r_ij < 0) only.  The Monaghan 1992 form:
         *
         *    Π_ij = (−α c_max μ_ij + β μ_ij²) / ρ̄_ij
         *
         *  where  μ_ij = h̄ (v_ij · r_ij) / (r_ij² + η²),  η = 0.1 h̄,
         *         c_max is the global maximum speed of sound (from setParams),
         *         ρ̄_ij = (ρ_i + ρ_j) / 2.
         *
         *  The linear term (α) diffuses velocity divergence; the quadratic term (β)
         *  is a Rankine–Hugoniot correction active only in shock-like conditions.
         *  For most weakly-compressible flows β = 0 and α ∈ [0.01, 0.1] suffices.
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
         *  For DENSITYSUMMATION and cross-phase fluid pairs (fluid 1 ↔ fluid 2),
         *  replaces the standard actual-density-weighted pressure average
         *
         *    p̄_ij = (ρ_j · p_i  +  ρ_i · p_j) / (ρ_i + ρ_j)            [standard]
         *
         *  with a rest-density-weighted form plus a hydrostatic correction:
         *
         *    p̄_ij = (ρ₀ⱼ · p_i  +  ρ₀ᵢ · p_j  +  ρ₀ᵢ ρ₀ⱼ (g · r_ij)) / (ρ₀ᵢ + ρ₀ⱼ)
         *
         *  where g = getAcceleration(timestep) is the body-force acceleration and
         *  r_ij = r_i − r_j.
         *
         *  Physical motivation: in hydrostatic equilibrium the pressure gradient is
         *  ∇p = ρ g, which is discontinuous at a density-ratio interface.  The
         *  standard density-weighted average introduces a systematic error proportional
         *  to (ρ₁ − ρ₂) that drives spurious "parasitic" interfacial velocities.
         *  Rest-density weighting removes this bias; the g · r_ij term ensures the
         *  discretised pressure gradient exactly reproduces gravity for a static column,
         *  regardless of density ratio.  Especially beneficial for water/air (ρ ratio
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
         *  acoustic impedances Z = ρ c to handle density and speed-of-sound contrasts
         *  across the interface naturally:
         *
         *    Z*_ij  =  Z_i · Z_j / (Z_i + Z_j),   Z = ρ c        [harmonic impedance]
         *    u_ij   =  (v_i − v_j) · (r_i − r_j) / (|r_ij| + η)  [signed radial vel.]
         *    avc    =  −β_R · Z*_ij · u_ij⁻ / ρ̄_ij               (u⁻ = min(u, 0))
         *
         *  where ρ̄_ij = (ρ_i + ρ_j) / 2.  Applied only to approaching fluid–fluid
         *  pairs (v_ij · r_ij < 0).  Mutually exclusive with Monaghan AV.
         *
         *  Advantages over Monaghan AV for two-phase flows:
         *    - No α/β parameters to tune per problem (only one global β_R ≈ 1).
         *    - Impedance mismatch at the interface is handled automatically:
         *      Z* → Z_light/2 when Z_heavy >> Z_light (e.g. water/air).
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
         *  Adds a diffusive correction to dρ/dt to smooth density oscillations:
         *
         *    dρ_i/dt += −2 δ h̄ c_max m_j (ρ_i/ρ₀ᵢ − ρ_j/ρ₀ⱼ) (r_ij · ∇W_ij) / (r²+η²)
         *
         *  The drive term (ρ_i/ρ₀ᵢ − ρ_j/ρ₀ⱼ) is the rest-density-normalised form
         *  rather than the original (ρ_i/ρ_j − 1).  The original is non-zero at
         *  equilibrium when the two phases have different rest densities (ρ₀₁ ≠ ρ₀₂),
         *  causing unphysical density drift across the interface.  The normalised form
         *  equals zero at equilibrium for both single- and two-phase flows, correcting
         *  stratified-flow artefacts.
         *
         *  Reference: Molteni & Colagrossi (2009) Comput. Phys. Commun. 180, 861–872.
         *             Two-phase correction: see also Grenier et al. (2013).
         *
         * \param ddiff  Diffusion coefficient δ (typically 0.1; range 0.05–0.2)
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

        /*! Activate δ⁺-SPH particle shifting (Sun et al. 2017, Comput. Fluids).
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

        /*! Turn δ⁺-SPH particle shifting off.
         */
        void deactivateParticleShifting()
            {
            m_particle_shifting = false;
            }

        //! Computes forces
        void computeForces(uint64_t timestep);
    
    #ifdef ENABLE_MPI
        /// The system's communicator.
        std::shared_ptr<Communicator> m_comm;
    #endif

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

        Scalar m_avalpha; //!< Monaghan AV: linear (volumetric) diffusion coefficient α
        Scalar m_avbeta;  //!< Monaghan AV: quadratic (shock) diffusion coefficient β
        Scalar m_riemann_beta; //!< Riemann dissipation: scaling coefficient β_R (Zhang et al. 2017)
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

        // Particle shifting (δ⁺-SPH, Sun et al. 2017)
        bool m_particle_shifting;            //!< True if δ⁺-SPH particle shifting is activated
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
        void compute_ndensity(uint64_t timestep);
        
        void compute_particlenumberdensity(uint64_t timestep);

        /*! Helper function to compute particle pressures
         *  \post Pressure of fluid particle computed
         */
        void compute_pressure(uint64_t timestep);

        /*! Helper function to compute fictitious solid particle properties (pressures and velocities)
        * \pre Ghost particle number densities (i.e. density array) must be up-to-date
        * \pre Solid normalization constant \sum_j w_ij must be computed and stored in density array
        * \post Fictitious particle properties are computed and stored in aux1 array
        */
        void compute_noslip(uint64_t timestep);

        /*! Helper function to compute particle concentration gradient for Fickian shifting
         * within the CSF computation
         * It overwrites h_pressure and the dot product of the gradient is stored in h_energy
         * Paper: Lind et al. 2012
        */
        void compute_particle_concentration_gradient(uint64_t timestep);

        /*! Compute and apply δ⁺-SPH particle position shifts (Sun et al. 2017).
         * \pre compute_colorgradients() + update_ghost_aux123() must precede this call
         *      (aux3 must hold up-to-date fluid-fluid interface normals).
         * \post Fluid particle positions updated by δr_i.
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
        void compute_colorgradients(uint64_t timestep);

        /*! Helper function to compute interfacial surface force field
        * \pre Normal vector field have been computed and communicated
        * \post Surface force density vectors are stored in aux4
        */
        void compute_surfaceforce(uint64_t timestep);

        /*! Helper function where the actual force computation takes place
         * \pre Number densities and fictitious solid particle properties must be up-to-date
         * \post h_force stores forces acting on fluid particles and .w component stores rate of change of density
         */
        void forcecomputation(uint64_t timestep);

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

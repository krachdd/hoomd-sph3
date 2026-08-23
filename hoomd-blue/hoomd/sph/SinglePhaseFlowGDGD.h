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
#include "SinglePhaseFlow.h"
#include "SolidFluidTypeBit.h"

#include "EvaluationMethodDefinition.h"


/*! \file SinglePhaseFlowGDGD.h
    \brief General Density Gradient Driven Flow (GDGDF) solver.

    Extends SinglePhaseFlow with a transported scalar field T (temperature or
    concentration) that modifies the local rest density via one of two models:

    1. Variable Reference Density (VRD) — buoyancy emerges implicitly from
       pressure gradients driven by per-particle rest densities:
           \f$\rho_{0,i} = \rho_{0,\mathrm{ref}} \cdot (1 - \beta \cdot (T_i - T_\mathrm{ref}))\f$
       For DENSITYSUMMATION, VRD pressures are computed on-the-fly in the pair
       loop.  For DENSITYCONTINUITY, the VRD derivative \f$\partial P/\partial\rho|_{\rho_{0,i}}\f$ is used
       in the \f$\mathrm{d}p/\mathrm{d}t\f$ chain rule.

    2. Boussinesq approximation — EOS uses the global \f$\rho_0\f$; per-particle
       buoyancy is applied as an explicit body force correction:
           \f$\Delta F_b = m_i \cdot g \cdot (-\beta \cdot (T_i - T_\mathrm{ref}))\f$
       so the total body force per particle is \f$m_i \cdot g \cdot (1 - \beta \cdot (T_i - T_\mathrm{ref}))\f$.

    In both modes a scalar diffusion term
        \f$\mathrm{d}T_i/\mathrm{d}t \mathrel{+}= (\kappa_s / V_i) \cdot (V_i^2 + V_j^2) \cdot (T_i - T_j) \cdot \mathrm{d}W/\mathrm{d}r / r\f$
    is accumulated into h_ratedpe.z and time-marched by the integrator
    (VelocityVerletBasic or KickDriftKickTV) via the h_dpedt.z / h_aux4.x
    storage convention:
        aux4.x  = T (scalar field value)
        dpedt.z = \f$\mathrm{d}T/\mathrm{d}t\f$ (rate, copied from net_ratedpe.z in integrateStepTwo)
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __SinglePhaseFlowGDGD_H__
#define __SinglePhaseFlowGDGD_H__


namespace hoomd
{
namespace sph
{

/*! \brief SPH solver for general density-gradient-driven flows.
 *
 *  Inherits the full SinglePhaseFlow pipeline (density, pressure, no-slip
 *  boundary, pair forces, body force).  Overrides forcecomputation() to add
 *  scalar transport and optional VRD/Boussinesq buoyancy; overrides
 *  computeForces() only to insert the MPI ghost update for aux4 at the
 *  correct point in the timestep.
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
class PYBIND11_EXPORT SinglePhaseFlowGDGD : public SinglePhaseFlow<KT_, SET_>
    {
    public:

        //! Constructor — same arguments as SinglePhaseFlow.
        SinglePhaseFlowGDGD(std::shared_ptr<SystemDefinition> sysdef,
                            std::shared_ptr<SmoothingKernel<KT_>> skernel,
                            std::shared_ptr<StateEquation<SET_>> equationofstate,
                            std::shared_ptr<nsearch::NeighborList> nlist,
                            std::shared_ptr<ParticleGroup> fluidgroup,
                            std::shared_ptr<ParticleGroup> solidgroup,
                            DensityMethod   mdensitymethod  = DENSITYSUMMATION,
                            ViscosityMethod mviscositymethod = HARMONICAVERAGE);

        //! Destructor
        virtual ~SinglePhaseFlowGDGD();

        /*! Set GDGD-specific parameters.
         *
         * \param kappa_s    Scalar diffusivity [\f$\mathrm{m}^2/\mathrm{s}\f$] — thermal conductivity / (\f$\rho \cdot c_p\f$)
         *                   for heat transfer, or mass diffusivity for species transport.
         * \param beta_s     Thermal / solutal expansion coefficient [1/K or 1/(unit T)].
         * \param scalar_ref Reference scalar value T_ref (or c_ref).
         * \param boussinesq If true, use Boussinesq approximation (explicit buoyancy force
         *                   correction, global \f$\rho_0\f$ in EOS).  If false, use Variable
         *                   Reference Density (buoyancy via per-particle EOS).
         * \param scalar_noflux If true, solid neighbours are excluded from the scalar
         *                   diffusion sum (zero-flux / adiabatic wall BC).  If false
         *                   (default), solid particles participate with their prescribed
         *                   aux4.x value (Dirichlet wall BC).  The skip-solid zero-flux
         *                   is consistent to leading order only (O(h): kernel support is
         *                   truncated at the wall); a ghost-mirror extrapolation would be
         *                   O(h^2) and is a possible future upgrade.
         */
        void setGDGDParams(Scalar kappa_s, Scalar beta_s, Scalar scalar_ref, bool boussinesq,
                           bool scalar_noflux = false);

        Scalar getKappaS()      { return m_kappa_s; }
        Scalar getBetaS()       { return m_beta_s; }
        Scalar getScalarRef()   { return m_scalar_ref; }
        bool   getBoussinesq()  { return m_boussinesq; }
        bool   getScalarNoflux(){ return m_scalar_noflux; }

        /*! Compute forces pipeline.
         *
         *  Adds an MPI ghost update of aux4 (scalar T) before delegating to
         *  SinglePhaseFlow::computeForces(), which calls our virtual
         *  forcecomputation() override via dynamic dispatch.
         *
         *  Rationale: integrateStepOne advances local T to t+dt/2; ghost T
         *  values must be synchronised before the pair loop reads neighbour T.
         */
        void computeForces(uint64_t timestep);

    #ifdef ENABLE_MPI
        /// The system's communicator (initialised in constructor).
        std::shared_ptr<Communicator> m_comm;

        /*! Ghost communication flags.
         *
         *  Adds auxiliary4 (scalar T field) on top of the flags the base class
         *  would request.  aux4 must be communicated to ghost particles so that
         *  neighbour scalar values are available in the pair force loop.
         */
        virtual CommFlags getRequestedCommFlags(uint64_t timestep)
            {
            CommFlags flags(0);
            flags[comm_flag::net_force]  = 0;
            flags[comm_flag::position]   = 1;  // position + type
            flags[comm_flag::velocity]   = 1;  // velocity + mass
            flags[comm_flag::density]    = 1;
            flags[comm_flag::pressure]   = 1;
            flags[comm_flag::energy]     = 0;
            flags[comm_flag::auxiliary1] = 1;  // fictitious solid velocity
            flags[comm_flag::auxiliary4] = 1;  // scalar field T (temperature / concentration)
            flags[comm_flag::slength]    = 1;
            flags |= ForceCompute::getRequestedCommFlags(timestep);
            return flags;
            }
    #endif

        //! Returns true: this class accumulates dpe array (dT/dt in .z component).
        virtual bool ComputesDPE()
            {
            return true;
            }

    protected:

        Scalar m_kappa_s;         //!< Scalar diffusivity [\f$\mathrm{m}^2/\mathrm{s}\f$]
        Scalar m_beta_s;          //!< Expansion coefficient [1/K or 1/concentration unit]
        Scalar m_scalar_ref;      //!< Reference scalar value (T_ref or c_ref)
        bool   m_boussinesq;      //!< True = Boussinesq; false = Variable Reference Density
        bool   m_scalar_noflux;   //!< True = zero-flux scalar wall BC (skip solid neighbours in diffusion); false = Dirichlet via solid aux4.x
        bool   m_gdgd_params_set; //!< True once setGDGDParams() has been called

        /*! Force computation: standard SPH pair forces + scalar diffusion + buoyancy.
         *
         *  Extension of SinglePhaseFlow::forcecomputation() with:
         *    - Scalar diffusion rate \f$\mathrm{d}T/\mathrm{d}t\f$ accumulated into h_ratedpe.data[i].z
         *      using the Morris-Fox-Zhu (1997) Laplacian SPH operator.
         *    - VRD mode (m_boussinesq == false):
         *        DENSITYSUMMATION  — pair pressures recomputed on-the-fly using
         *                           PressureVRD(\f$\rho_i\f$, \f$\rho_{0,i}(T_i)\f$) per particle.
         *        DENSITYCONTINUITY — \f$\mathrm{d}p/\mathrm{d}t\f$ chain rule uses dPressureVRDdDensity.
         *    - Boussinesq mode (m_boussinesq == true):
         *        Standard EOS pressures used; per-particle buoyancy correction
         *        \f$\Delta F_b = m_i \cdot g \cdot (-\beta \cdot (T_i - T_\mathrm{ref}))\f$ added to h_force.
         *
         *  Scalar T is read from aux4.x; \f$\mathrm{d}T/\mathrm{d}t\f$ is written to h_ratedpe.data[i].z.
         *  The integrator copies net_ratedpe.z \f$\rightarrow\f$ dpedt.z (step 2) and uses
         *  dpedt.z to advance aux4.x in both half-steps.
         */
        virtual void forcecomputation(uint64_t timestep);

        /*! VRD-aware pressure evaluation.
         *
         *  In VRD mode (m_boussinesq == false) the fluid pressure is evaluated
         *  with the per-particle rest density
         *  \f$\rho_{0,i} = \rho_0 (1 - \beta (T_i - T_\mathrm{ref}))\f$ so that the
         *  pressure written to h_pressure (used by the continuity-mode pair
         *  loop and the Adami wall interpolation) carries the buoyancy signal.
         *  In Boussinesq mode this falls back to the base class implementation.
         */
        virtual void compute_pressure(uint64_t timestep);

    #ifdef ENABLE_MPI
        /*! Update ghost particles for the scalar field (aux4).
         *
         *  Called at the beginning of computeForces() so that neighbour scalar
         *  T values in the pair loop reflect the t+dt/2 values set by the
         *  integrator's first half-step.
         */
        void update_ghost_aux4(uint64_t timestep);
    #endif

    private:

    };


namespace detail
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlowGDGD(pybind11::module& m, std::string name);

} // end namespace detail
} // end namespace sph
} // end namespace hoomd

#endif // __SinglePhaseFlowGDGD_H__

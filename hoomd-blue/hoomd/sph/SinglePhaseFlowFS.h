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
#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include "SmoothingKernel.h"
#include "StateEquations.h"
#include "SPHBaseClass.h"
#include "SinglePhaseFlow.h"
#include "SinglePhaseFlowTV.h"
#include "SolidFluidTypeBit.h"

#include "EvaluationMethodDefinition.h"


/*! \file SinglePhaseFlowFS.h
    \brief Free-surface SPH solver for single-phase flows (thin films, waves, jets).

    Extends SinglePhaseFlowTV with:

    1. Free-surface detection via the Shepard kernel-completeness ratio
           \f$\lambda_i = V_i W_0(h) + \sum_{j \neq i} V_j W(r_{ij}, h)\f$
       Particles with \f$\lambda_i <\f$ fs_threshold are flagged as free-surface particles.
       Their outward unit normal is \f$\hat{n}_i = -\nabla\lambda_i / |\nabla\lambda_i|\f$.
       Storage: aux2 (normals, temporarily), aux4.x (lambda), aux4.y (curvature).

    2. Contact-angle enforcement (Huber et al. 2016) at the triple line:
       near solid walls the free-surface normal is blended with the wall normal:
           \f$\hat{n}_\mathrm{corrected} = \sin(\theta_\mathrm{eq})\,\hat{t}_w + \cos(\theta_\mathrm{eq})\,\hat{n}_w\f$
       where \f$\hat{n}_w\f$ is the wall inward normal and \f$\hat{t}_w\f$ the tangential component
       of the gas-side free-surface normal.

    3. Curvature estimation:
           \f$\kappa_i = (1/V_i) \sum_j V_j (\hat{n}_j - \hat{n}_i) \cdot \nabla W_{ij}\f$
       Only computed for surface particles; stored in aux4.y.

    4. Free-surface pressure clamping:
       Prevents unphysical tensile pressure at the surface (\f$P \leftarrow \max(0, P)\f$).

    5. Continuum surface force (CSF) in forcecomputation:
           \f$F_{\sigma,i} = -\sigma \cdot \kappa_i \cdot \hat{n}_i \cdot (m_i/\rho_i)\f$
       The minus sign follows from CSF conventions: for a convex surface
       (\f$\kappa > 0\f$, \f$\hat{n}\f$ outward), the surface tension pulls inward.

    Aux-array usage summary
    -----------------------
    aux1 : fictitious solid velocity (inherited, Adami 2012)
    aux2 : free-surface outward normal (detect_freesurface \f$\rightarrow\f$ compute_curvature),
           then overwritten by BPC in forcecomputation (as in TV, for KickDriftKickTV)
    aux3 : transport velocity (inherited from TV)
    aux4 : .x = \f$\lambda\f$ (kernel completeness), .y = \f$\kappa\f$ (curvature), .z unused

    References
    ----------
    Marrone et al. (2010), Comput. Fluids — free-surface detection
    Huber et al.  (2016), Int. J. Numer. Meth. Fluids — contact-angle BC
    Colagrossi & Landrini (2003), J. Comput. Phys. — CSF surface tension for SPH
    Adami et al. (2013), J. Comput. Phys. — transport velocity
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __SinglePhaseFlowFS_H__
#define __SinglePhaseFlowFS_H__


namespace hoomd
{
namespace sph
{

/*!  SPH solver with free-surface detection, contact-angle BC, and surface tension.
 *
 *   Inherits the full SinglePhaseFlowTV pipeline (density, pressure, no-slip BC,
 *   TV pair forces, body force, KickDriftKickTV integrator support).
 *   Overrides computeForces() to insert the free-surface pipeline and
 *   overrides forcecomputation() to add the CSF surface tension force.
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
class PYBIND11_EXPORT SinglePhaseFlowFS : public SinglePhaseFlowTV<KT_, SET_>
    {
    public:

        //! Constructor — same signature as SinglePhaseFlowTV.
        SinglePhaseFlowFS(std::shared_ptr<SystemDefinition> sysdef,
                          std::shared_ptr<SmoothingKernel<KT_>> skernel,
                          std::shared_ptr<StateEquation<SET_>> equationofstate,
                          std::shared_ptr<nsearch::NeighborList> nlist,
                          std::shared_ptr<ParticleGroup> fluidgroup,
                          std::shared_ptr<ParticleGroup> solidgroup,
                          DensityMethod   mdensitymethod  = DENSITYSUMMATION,
                          ViscosityMethod mviscositymethod = HARMONICAVERAGE);

        //! Destructor
        virtual ~SinglePhaseFlowFS();

        /*! Set free-surface parameters.
         *
         * \param sigma          Surface tension coefficient \f$\sigma\f$ [N/m].  Set to 0
         *                       to disable surface tension while keeping
         *                       free-surface detection active.
         * \param fs_threshold   Kernel-completeness threshold \f$\lambda \in (0,1)\f$.
         *                       Particles with \f$\lambda <\f$ fs_threshold are treated as
         *                       free-surface particles (typical value: 0.75).
         * \param contact_angle  Equilibrium contact angle \f$\theta\f$ [rad] measured
         *                       inside the liquid from the solid wall.
         *                       \f$\pi/2\f$ = neutral wetting (no correction applied),
         *                       0   = complete wetting,
         *                       \f$\pi\f$   = complete non-wetting.
         */
        void setFSParams(Scalar sigma, Scalar fs_threshold, Scalar contact_angle);

        Scalar getSigma()        { return m_sigma; }
        Scalar getFSThreshold()  { return m_fs_threshold; }
        Scalar getContactAngle() { return m_contact_angle; }

        //! Override: insert free-surface pipeline before forcecomputation
        void computeForces(uint64_t timestep);

    #ifdef ENABLE_MPI
        //! Ghost-communication flags: TV flags + aux4 (lambda / curvature)
        virtual CommFlags getRequestedCommFlags(uint64_t timestep)
            {
            CommFlags flags(0);
            flags[comm_flag::net_force]  = 0;
            flags[comm_flag::position]   = 1;
            flags[comm_flag::velocity]   = 1;
            flags[comm_flag::density]    = 1;
            flags[comm_flag::pressure]   = 1;
            flags[comm_flag::energy]     = 0;
            flags[comm_flag::auxiliary1] = 1;  // fictitious solid velocity
            flags[comm_flag::auxiliary2] = 1;  // BPC / fs normals (shared slot)
            flags[comm_flag::auxiliary3] = 1;  // transport velocity
            flags[comm_flag::auxiliary4] = 1;  // lambda + curvature
            flags[comm_flag::slength]    = 1;
            flags |= ForceCompute::getRequestedCommFlags(timestep);
            return flags;
            }
    #endif

        virtual bool ComputesDPE() { return true; }

    protected:

    #ifdef ENABLE_MPI
        /// MPI communicator (shadowed from TV so FS can access it directly)
        std::shared_ptr<Communicator> m_comm;
    #endif

        // ── Free-surface parameters ───────────────────────────────────────────
        Scalar m_sigma;          //!< Surface tension coefficient \f$\sigma\f$ [N/m]
        Scalar m_fs_threshold;   //!< \f$\lambda\f$ threshold for surface detection
        Scalar m_contact_angle;  //!< Equilibrium contact angle [rad]

        // ── Step 1: detect free surface ───────────────────────────────────────
        /*! Compute kernel-completeness \f$\lambda_i\f$ and outward free-surface normals \f$\hat{n}_i\f$.
         *
         *  Writes:
         *    aux4.x \f$\leftarrow \lambda_i\f$  (Shepard sum including self contribution)
         *    aux2   \f$\leftarrow \hat{n}_i\f$  (unit outward normal; {0,0,0} for bulk particles)
         *
         *  When the particle has solid neighbours and \f$|\theta - \pi/2| > 0.01\f$ rad,
         *  the contact-angle correction is applied to \f$\hat{n}_i\f$.
         *
         *  \pre  Density and smoothing lengths are up to date.
         *  \post aux2 and aux4.x are ready for compute_curvature().
         */
        virtual void detect_freesurface(uint64_t timestep);

        // ── Step 2: compute curvature ─────────────────────────────────────────
        /*! Compute mean curvature \f$\kappa_i = (1/V_i) \sum_j V_j (\hat{n}_j - \hat{n}_i) \cdot \nabla W_{ij}\f$.
         *
         *  Only executed for surface particles (\f$\lambda_i <\f$ m_fs_threshold).
         *  Writes:
         *    aux4.y \f$\leftarrow \kappa_i\f$
         *
         *  \pre  detect_freesurface() has been called (aux2 holds \f$\hat{n}_i\f$).
         *        In MPI runs, update_ghost_aux24() must be called first so that
         *        ghost-particle normals are available.
         *  \post aux4.y is ready for forcecomputation().
         */
        virtual void compute_curvature(uint64_t timestep);

        // ── Step 3: free-surface pressure clamping ────────────────────────────
        /*! Clamp pressure of free-surface particles to \f$P \geq 0\f$.
         *
         *  Prevents unphysical tensile pressures that arise from the kernel
         *  truncation near the surface.  Applied only to particles with
         *  \f$\lambda_i <\f$ m_fs_threshold.
         */
        virtual void apply_freesurface_pressure(uint64_t timestep);

        // ── Step 4: force loop with CSF surface tension ───────────────────────
        /*! TV force loop extended with CSF surface tension.
         *
         *  At the start of each particle's outer loop the free-surface normal
         *  is saved from aux2 (written by detect_freesurface), then aux2 is
         *  zeroed for BPC accumulation exactly as in SinglePhaseFlowTV.
         *  After the neighbour loop the CSF surface tension force is added:
         *      \f$F_{\sigma,i} \mathrel{+}= -\sigma \cdot \kappa_i \cdot \hat{n}_i \cdot (m_i/\rho_i)\f$
         *
         *  \pre  detect_freesurface(), compute_curvature(), and
         *        apply_freesurface_pressure() have been called.
         */
        void forcecomputation(uint64_t timestep);

    #ifdef ENABLE_MPI
        /*! Sync aux2 (fs normals) and aux4 (\f$\lambda\f$ + \f$\kappa\f$) to ghost particles.
         *
         *  Called between detect_freesurface() and compute_curvature() so that
         *  neighbour normals from other MPI ranks are available in the curvature
         *  loop.
         */
        void update_ghost_aux24(uint64_t timestep);
    #endif

    private:

    };


namespace detail
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlowFS(pybind11::module& m, std::string name);

} // end namespace detail
} // end namespace sph
} // end namespace hoomd

#endif // __SinglePhaseFlowFS_H__

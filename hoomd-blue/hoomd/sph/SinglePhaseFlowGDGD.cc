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


#include "SinglePhaseFlowGDGD.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

using namespace std;

namespace hoomd
{
namespace sph
{

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

template<SmoothingKernelType KT_, StateEquationType SET_>
SinglePhaseFlowGDGD<KT_, SET_>::SinglePhaseFlowGDGD(
        std::shared_ptr<SystemDefinition> sysdef,
        std::shared_ptr<SmoothingKernel<KT_>> skernel,
        std::shared_ptr<StateEquation<SET_>> equationofstate,
        std::shared_ptr<nsearch::NeighborList> nlist,
        std::shared_ptr<ParticleGroup> fluidgroup,
        std::shared_ptr<ParticleGroup> solidgroup,
        DensityMethod   mdensitymethod,
        ViscosityMethod mviscositymethod)
    : SinglePhaseFlow<KT_, SET_>(sysdef, skernel, equationofstate, nlist,
                                  fluidgroup, solidgroup,
                                  mdensitymethod, mviscositymethod)
    {
    this->m_exec_conf->msg->notice(5) << "Constructing SinglePhaseFlowGDGD" << std::endl;

    // Initialise all inherited flags to safe defaults (mirrors SinglePhaseFlowTV)
    this->m_const_slength          = false;
    this->m_params_set             = false;
    this->m_compute_solid_forces   = false;
    this->m_artificial_viscosity   = false;
    this->m_density_diffusion      = false;
    this->m_shepard_renormalization= false;
    this->m_density_reinitialization = false;
    this->m_ch     = Scalar(0.0);
    this->m_rcut   = Scalar(0.0);
    this->m_rcutsq = Scalar(0.0);
    this->m_avalpha = Scalar(0.0);
    this->m_avbeta  = Scalar(0.0);
    this->m_ddiff   = Scalar(0.0);
    this->m_shepardfreq         = 1;
    this->m_densityreinitfreq   = 1;

    this->m_solid_removed        = false;
    this->m_pressure_initialized = false;

    // Sanity checks
    assert(this->m_pdata);
    assert(this->m_nlist);
    assert(this->m_skernel);
    assert(this->m_eos);

    // Build fluid / solid type vectors
    this->constructTypeVectors(fluidgroup, &this->m_fluidtypes);
    this->constructTypeVectors(solidgroup,  &this->m_solidtypes);

    // Initialise the type-property map (FLUID / SOLID bits)
    unsigned int num_types = this->m_sysdef->getParticleData()->getNTypes();
    this->m_type_property_map = GPUArray<unsigned int>(num_types, this->m_exec_conf);
        { // GPU Array Scope
        ArrayHandle<unsigned int> h_type_property_map(
            this->m_type_property_map, access_location::host, access_mode::overwrite);
        fill_n(h_type_property_map.data, num_types, SolidFluidTypeBit::NONE);
        for (unsigned int i = 0; i < this->m_fluidtypes.size(); i++)
            h_type_property_map.data[this->m_fluidtypes[i]] |= SolidFluidTypeBit::FLUID;
        for (unsigned int i = 0; i < this->m_solidtypes.size(); i++)
            h_type_property_map.data[this->m_solidtypes[i]] |= SolidFluidTypeBit::SOLID;
        } // End GPU Array Scope

    // Store integration methods
    this->m_density_method   = mdensitymethod;
    this->m_viscosity_method = mviscositymethod;

    // Cache frequently needed kernel / EOS scalars
    this->m_rho0  = equationofstate->getRestDensity();
    this->m_c     = equationofstate->getSpeedOfSound();
    this->m_kappa = skernel->getKernelKappa();

    // Register cut-off matrix with the neighbour list
    this->m_r_cut_nlist = std::make_shared<GPUArray<Scalar>>(
        this->m_typpair_idx.getNumElements(), this->m_exec_conf);
    this->m_nlist->addRCutMatrix(this->m_r_cut_nlist);

    // GDGD-specific parameters — defaults (must call setGDGDParams before run)
    m_kappa_s         = Scalar(0.0);
    m_beta_s          = Scalar(0.0);
    m_scalar_ref      = Scalar(0.0);
    m_boussinesq      = false;
    m_gdgd_params_set = false;

#ifdef ENABLE_MPI
    if (this->m_sysdef->isDomainDecomposed())
        {
        auto comm_weak = this->m_sysdef->getCommunicator();
        assert(comm_weak.lock());
        m_comm = comm_weak.lock();
        }
#endif
    }


// ─────────────────────────────────────────────────────────────────────────────
// Destructor
// ─────────────────────────────────────────────────────────────────────────────

template<SmoothingKernelType KT_, StateEquationType SET_>
SinglePhaseFlowGDGD<KT_, SET_>::~SinglePhaseFlowGDGD()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying SinglePhaseFlowGDGD" << std::endl;
    }


// ─────────────────────────────────────────────────────────────────────────────
// Parameter setter
// ─────────────────────────────────────────────────────────────────────────────

/*! \param kappa_s    Scalar diffusivity [m²/s]
    \param beta_s    Expansion coefficient [1/K or 1/concentration unit]
    \param scalar_ref Reference scalar value (T_ref or c_ref)
    \param boussinesq If true, Boussinesq approximation; if false, VRD approach
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowGDGD<KT_, SET_>::setGDGDParams(
        Scalar kappa_s, Scalar beta_s, Scalar scalar_ref, bool boussinesq)
    {
    m_kappa_s         = kappa_s;
    m_beta_s          = beta_s;
    m_scalar_ref      = scalar_ref;
    m_boussinesq      = boussinesq;
    m_gdgd_params_set = true;
    }


// ─────────────────────────────────────────────────────────────────────────────
// MPI helper: update ghost particles for the scalar field (aux4)
// ─────────────────────────────────────────────────────────────────────────────

#ifdef ENABLE_MPI
template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowGDGD<KT_, SET_>::update_ghost_aux4(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7)
        << "Computing SinglePhaseFlowGDGD::update_ghost_aux4" << std::endl;

    if (this->m_comm)
        {
        // Communicate only aux4 (scalar T field); all other flags zero so that
        // expensive density / pressure ghost updates are not duplicated here.
        CommFlags flags(0);
        flags[comm_flag::auxiliary4] = 1;  // scalar field T
        this->m_comm->setFlags(flags);
        this->m_comm->beginUpdateGhosts(timestep);
        this->m_comm->finishUpdateGhosts(timestep);
        }
    }
#endif // ENABLE_MPI


// ─────────────────────────────────────────────────────────────────────────────
// computeForces — inserts ghost aux4 update then delegates to base pipeline
// ─────────────────────────────────────────────────────────────────────────────

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowGDGD<KT_, SET_>::computeForces(uint64_t timestep)
    {
#ifdef ENABLE_MPI
    // After integrateStepOne has advanced local T to t+dt/2, propagate the
    // updated scalar T values to ghost particles so the pair loop reads
    // consistent T(t+dt/2) for all neighbours.
    update_ghost_aux4(timestep);
#endif

    // Delegate to SinglePhaseFlow::computeForces(), which executes the full
    // pipeline (nlist, density, pressure, noslip, ...) and calls our virtual
    // forcecomputation() override via dynamic dispatch.
    SinglePhaseFlow<KT_, SET_>::computeForces(timestep);
    }


// ─────────────────────────────────────────────────────────────────────────────
// forcecomputation — extended pair-force loop with scalar transport
// ─────────────────────────────────────────────────────────────────────────────

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowGDGD<KT_, SET_>::forcecomputation(uint64_t timestep)
    {
    if (this->m_density_method == DENSITYSUMMATION)
        this->m_exec_conf->msg->notice(7)
            << "Computing SinglePhaseFlowGDGD::Forces (SUMMATION)" << std::endl;
    else if (this->m_density_method == DENSITYCONTINUITY)
        this->m_exec_conf->msg->notice(7)
            << "Computing SinglePhaseFlowGDGD::Forces (CONTINUITY)" << std::endl;

    const BoxDim& box       = this->m_pdata->getGlobalBox();
    const unsigned int group_size = this->m_fluidgroup->getNumMembers();

        { // GPU Array Scope
        // ── Output arrays (zeroed before accumulation) ────────────────────────
        ArrayHandle<Scalar4> h_force  (this->m_force,   access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar4> h_ratedpe(this->m_ratedpe, access_location::host, access_mode::readwrite);

        assert(h_force.data);
        assert(h_ratedpe.data);

        memset((void*)h_force.data,   0, sizeof(Scalar4) * this->m_force.getNumElements());
        memset((void*)h_ratedpe.data, 0, sizeof(Scalar4) * this->m_ratedpe.getNumElements());

        // ── Input particle data ───────────────────────────────────────────────
        ArrayHandle<Scalar4> h_pos     (this->m_pdata->getPositions(),    access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(),   access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_density (this->m_pdata->getDensities(),    access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_pressure(this->m_pdata->getPressures(),    access_location::host, access_mode::read);
        // aux1.x/y/z: fictitious solid-boundary velocity (set by compute_noslip)
        ArrayHandle<Scalar3> h_vf      (this->m_pdata->getAuxiliaries1(), access_location::host, access_mode::read);
        // aux4.x: scalar field T (temperature / concentration).
        // Ghost values updated by update_ghost_aux4() at the start of computeForces().
        ArrayHandle<Scalar3> h_aux4    (this->m_pdata->getAuxiliaries4(), access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_h       (this->m_pdata->getSlengths(),     access_location::host, access_mode::read);

        // ── Neighbour list ────────────────────────────────────────────────────
        ArrayHandle<unsigned int> h_n_neigh  (this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist_arr(this->m_nlist->getNListArray(),  access_location::host, access_mode::read);
        ArrayHandle<size_t>       h_head_list(this->m_nlist->getHeadList(),    access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

        assert(h_pos.data != NULL);

        unsigned int size;
        size_t myHead;
        double max_vel = 0.0;

        // ── Fluid-particle outer loop ─────────────────────────────────────────
        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

            // Position of particle i
            Scalar3 pi;
            pi.x = h_pos.data[i].x;
            pi.y = h_pos.data[i].y;
            pi.z = h_pos.data[i].z;

            // Velocity and mass of particle i
            Scalar3 vi;
            vi.x = h_velocity.data[i].x;
            vi.y = h_velocity.data[i].y;
            vi.z = h_velocity.data[i].z;
            Scalar mi = h_velocity.data[i].w;

            // Pre-computed EOS quantities for particle i
            Scalar Pi   = h_pressure.data[i];
            Scalar rhoi = h_density.data[i];
            Scalar Vi   = mi / rhoi;

            // Scalar field value of particle i (temperature or concentration)
            Scalar Ti = h_aux4.data[i].x;

            // Per-particle local rest density (VRD) or global rest density (Boussinesq).
            // ρ₀_i = ρ₀_ref · (1 − β · (T_i − T_ref))
            // Used in: (a) VRD pressure computation on-the-fly for SUMMATION,
            //          (b) VRD dp/dρ chain rule for CONTINUITY,
            //          (c) Boussinesq: set to global ρ₀ (unchanged EOS).
            Scalar rho0_i = m_boussinesq
                            ? this->m_rho0
                            : this->m_rho0 * (Scalar(1.0) - m_beta_s * (Ti - m_scalar_ref));

            // On-the-fly VRD pressure for particle i (DENSITYSUMMATION + VRD only).
            // For DENSITYCONTINUITY, h_pressure already holds the time-integrated
            // VRD-consistent pressure — do not overwrite it.
            Scalar Pi_use = Pi;
            if (!m_boussinesq && this->m_density_method == DENSITYSUMMATION)
                Pi_use = this->m_eos->PressureVRD(rhoi, rho0_i);

            // Adaptive timestep tracking
            Scalar vi_total = sqrt(vi.x*vi.x + vi.y*vi.y + vi.z*vi.z);
            if (i == 0) { max_vel = vi_total; }
            else if (vi_total > max_vel) { max_vel = vi_total; }

            // ── Neighbour inner loop ──────────────────────────────────────────
            myHead = h_head_list.data[i];
            size   = (unsigned int)h_n_neigh.data[i];

            for (unsigned int j = 0; j < size; j++)
                {
                // Index of neighbour particle (real or ghost)
                unsigned int k = h_nlist_arr.data[myHead + j];
                assert(k < this->m_pdata->getN() + this->m_pdata->getNGhosts());

                // Neighbour position and type
                Scalar3 pj;
                pj.x = h_pos.data[k].x;
                pj.y = h_pos.data[k].y;
                pj.z = h_pos.data[k].z;

                bool issolid = checksolid(h_type_property_map.data, h_pos.data[k].w);

                // Distance vector with periodic boundary conditions
                Scalar3 dx;
                dx.x = pi.x - pj.x;
                dx.y = pi.y - pj.y;
                dx.z = pi.z - pj.z;
                dx = box.minImage(dx);

                Scalar rsq = dot(dx, dx);

                // Skip neighbours outside the fixed cut-off (if constant h is used)
                if (this->m_const_slength && rsq > this->m_rcutsq)
                    continue;

                // Neighbour velocity: fictitious (Adami 2012) for solid, physical for fluid
                Scalar3 vj = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
                Scalar  mj = h_velocity.data[k].w;
                if (issolid)
                    {
                    vj.x = h_vf.data[k].x;
                    vj.y = h_vf.data[k].y;
                    vj.z = h_vf.data[k].z;
                    }
                else
                    {
                    vj.x = h_velocity.data[k].x;
                    vj.y = h_velocity.data[k].y;
                    vj.z = h_velocity.data[k].z;
                    }

                Scalar rhoj = h_density.data[k];
                Scalar Vj   = mj / rhoj;

                // Scalar field of neighbour k.
                // For solid particles, aux4.x may encode a Dirichlet wall temperature
                // boundary condition (set directly from Python before each run step).
                Scalar Tj = h_aux4.data[k].x;

                // VRD pressure for neighbour k (DENSITYSUMMATION + VRD + fluid only).
                // For DENSITYCONTINUITY: h_pressure holds the time-integrated pressure.
                // For solid neighbours: Adami-interpolated pressure is used unchanged.
                Scalar Pj_use = h_pressure.data[k];
                if (!m_boussinesq && this->m_density_method == DENSITYSUMMATION && !issolid)
                    {
                    // Per-particle rest density of neighbour k
                    Scalar rho0_k = this->m_rho0
                                    * (Scalar(1.0) - m_beta_s * (Tj - m_scalar_ref));
                    Pj_use = this->m_eos->PressureVRD(rhoj, rho0_k);
                    }

                // Velocity difference (physical velocities for both fluid and solid here)
                Scalar3 dv;
                dv.x = vi.x - vj.x;
                dv.y = vi.y - vj.y;
                dv.z = vi.z - vj.z;

                // Kernel derivative quantities
                Scalar r      = sqrt(rsq);
                Scalar meanh  = this->m_const_slength
                                ? this->m_ch
                                : Scalar(0.5) * (h_h.data[i] + h_h.data[k]);
                Scalar epssqr = Scalar(0.01) * meanh * meanh;
                Scalar dwdr   = this->m_skernel->dwijdr(meanh, r);
                Scalar dwdr_r = dwdr / (r + epssqr);

                // Symmetric volume factor (Vi² + Vj²) — used for pressure (SUMMATION),
                // viscosity, and scalar diffusion.
                Scalar vijsqr = Vi*Vi + Vj*Vj;

                // ── Pressure force ────────────────────────────────────────────
                // DENSITYSUMMATION: symmetric Adami 2013 formulation.
                //   temp0 = (V_i²+V_j²) · (ρ_j·P_i + ρ_i·P_j) / (ρ_i+ρ_j)
                //   F_p   = −temp0 · dW/dr / r · dx
                // DENSITYCONTINUITY: weakly-compressible pressure interaction.
                //   temp0 = m_i·m_j · (P_i + P_j) / (ρ_i·ρ_j)
                //   (vijsqr factor is absorbed into mi·mj normalisation)
                Scalar temp0;
                if (this->m_density_method == DENSITYSUMMATION)
                    temp0 = vijsqr * (rhoj * Pi_use + rhoi * Pj_use) / (rhoi + rhoj);
                else // DENSITYCONTINUITY
                    temp0 = mi * mj * (Pi_use + Pj_use) / (rhoi * rhoj);

                // ── Artificial viscosity (optional, Monaghan 1983) ─────────────
                Scalar avc = 0.0;
                if (this->m_artificial_viscosity && !issolid)
                    {
                    Scalar dotdvdx = dot(dv, dx);
                    if (dotdvdx < Scalar(0))
                        {
                        Scalar muij    = meanh * dotdvdx / (rsq + epssqr);
                        Scalar meanrho = Scalar(0.5) * (rhoi + rhoj);
                        avc = (-this->m_avalpha * this->m_c * muij
                               + this->m_avbeta * muij * muij) / meanrho;
                        // Scale avc consistently with the pressure scaling
                        if (this->m_density_method == DENSITYSUMMATION)
                            avc *= vijsqr;
                        else // DENSITYCONTINUITY
                            avc *= mi * mj;
                        }
                    }

                // Accumulate pressure + AV force
                h_force.data[i].x -= (temp0 + avc) * dwdr_r * dx.x;
                h_force.data[i].y -= (temp0 + avc) * dwdr_r * dx.y;
                h_force.data[i].z -= (temp0 + avc) * dwdr_r * dx.z;

                // ── Viscous force ─────────────────────────────────────────────
                // F_visc = μ_eff · (V_i²+V_j²) · dW/dr/r · (v_i − v_j)
                // (same formula for both density methods — uses vijsqr)
                {
                Scalar dvnorm    = sqrt(dot(dv, dv));
                Scalar gamma_dot = dvnorm / (r + sqrt(epssqr));
                Scalar mu_eff    = computeNNViscosity(this->m_mu, gamma_dot, this->m_nn_model,
                                                      this->m_nn_K, this->m_nn_n,
                                                      this->m_nn_mu0, this->m_nn_muinf,
                                                      this->m_nn_lambda, this->m_nn_tauy,
                                                      this->m_nn_m, this->m_nn_mu_min);
                temp0 = mu_eff * vijsqr * dwdr_r;
                }
                h_force.data[i].x += temp0 * dv.x;
                h_force.data[i].y += temp0 * dv.y;
                h_force.data[i].z += temp0 * dv.z;

                // ── Scalar diffusion ──────────────────────────────────────────
                // Morris-Fox-Zhu (1997) SPH Laplacian operator applied to T:
                //   dT_i/dt += (κ_s / V_i) · (V_i²+V_j²) · (T_i−T_j) · dW/dr / r
                // Units: [m²/s]/[m³] · [m⁶] · [K] · [m⁻⁵] = [K/s]  ✓
                // Applied for all neighbours including solid particles, enabling
                // Dirichlet wall temperature BCs via prescribed aux4.x on solid particles.
                h_ratedpe.data[i].z += m_kappa_s / Vi * vijsqr * (Ti - Tj) * dwdr_r;

                // ── Density continuity rate (DENSITYCONTINUITY only) ──────────
                if (this->m_density_method == DENSITYCONTINUITY)
                    {
                    if (issolid)
                        {
                        // Use physical velocity for the continuity equation
                        // (not the fictitious Adami-interpolated velocity)
                        vj.x = h_velocity.data[k].x;
                        vj.y = h_velocity.data[k].y;
                        vj.z = h_velocity.data[k].z;
                        dv.x = vi.x - vj.x;
                        dv.y = vi.y - vj.y;
                        dv.z = vi.z - vj.z;
                        }

                    // dρ_i/dt += ρ_i · V_j · (v_i − v_j) · ∇W_{ij}
                    h_ratedpe.data[i].x += rhoi * Vj * dot(dv, dwdr_r * dx);

                    // Density diffusion (Molteni & Colagrossi 2009), if enabled
                    if (!issolid && this->m_density_diffusion)
                        h_ratedpe.data[i].x -=
                            (Scalar(2) * this->m_ddiff * meanh * this->m_c * mj
                             * (rhoi / rhoj - Scalar(1)) * dot(dx, dwdr_r * dx))
                            / (rsq + epssqr);
                    }

                } // end neighbour loop

            // ── Post-pair per-particle updates ────────────────────────────────

            // dp/dt chain rule (DENSITYCONTINUITY only).
            // In VRD mode: ∂P/∂ρ is evaluated at the per-particle rest density ρ₀_i.
            // In Boussinesq mode: global ρ₀ is used (standard single-phase EOS).
            if (this->m_density_method == DENSITYCONTINUITY)
                {
                Scalar dpdrho = m_boussinesq
                                ? this->m_eos->dPressuredDensity(rhoi)
                                : this->m_eos->dPressureVRDdDensity(rhoi, rho0_i);
                h_ratedpe.data[i].y = dpdrho * h_ratedpe.data[i].x;
                }

            // ── Boussinesq buoyancy correction ────────────────────────────────
            // Standard gravity F_g = m·g is applied uniformly by applyBodyForce().
            // Boussinesq adds the density-deviation-driven buoyancy part:
            //   ΔF_b = m_i · g · (−β · (T_i − T_ref))
            // so the net body force per particle becomes m·g·(1 − β·(T_i − T_ref)).
            // In VRD mode this correction is NOT applied: buoyancy emerges naturally
            // from the variable-reference-density pressure gradient.
            if (m_boussinesq)
                {
                Scalar buoy_factor = -m_beta_s * (Ti - m_scalar_ref);
                h_force.data[i].x += mi * this->m_bodyforce.x * buoy_factor;
                h_force.data[i].y += mi * this->m_bodyforce.y * buoy_factor;
                h_force.data[i].z += mi * this->m_bodyforce.z * buoy_factor;
                }

            } // end fluid-particle outer loop

        this->m_timestep_list[5] = max_vel;

        } // End GPU Array Scope

    // Apply uniform gravity / body force to all fluid particles.
    // In Boussinesq mode this contributes F_g = m·g; the correction ΔF_b added
    // inside the loop makes the total m·g·(1 − β·(T−T_ref)).
    // In VRD mode this contributes F_g = m·g; buoyancy is implicit in pressure.
    this->applyBodyForce(timestep, this->m_fluidgroup);
    }


// ─────────────────────────────────────────────────────────────────────────────
// pybind11 export
// ─────────────────────────────────────────────────────────────────────────────

namespace detail
{

template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlowGDGD(pybind11::module& m, std::string name)
    {
    pybind11::class_<SinglePhaseFlowGDGD<KT_, SET_>,
                     SPHBaseClass<KT_, SET_>,
                     std::shared_ptr<SinglePhaseFlowGDGD<KT_, SET_>>>(m, name.c_str())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<SmoothingKernel<KT_>>,
                            std::shared_ptr<StateEquation<SET_>>,
                            std::shared_ptr<nsearch::NeighborList>,
                            std::shared_ptr<ParticleGroup>,
                            std::shared_ptr<ParticleGroup>,
                            DensityMethod,
                            ViscosityMethod>())
        .def("setParams",          &SinglePhaseFlowGDGD<KT_, SET_>::setParams)
        .def("setGDGDParams",      &SinglePhaseFlowGDGD<KT_, SET_>::setGDGDParams)
        .def("getKappaS",          &SinglePhaseFlowGDGD<KT_, SET_>::getKappaS)
        .def("getBetaS",           &SinglePhaseFlowGDGD<KT_, SET_>::getBetaS)
        .def("getScalarRef",       &SinglePhaseFlowGDGD<KT_, SET_>::getScalarRef)
        .def("getBoussinesq",      &SinglePhaseFlowGDGD<KT_, SET_>::getBoussinesq)
        .def("getDensityMethod",   &SinglePhaseFlowGDGD<KT_, SET_>::getDensityMethod)
        .def("setDensityMethod",   &SinglePhaseFlowGDGD<KT_, SET_>::setDensityMethod)
        .def("getViscosityMethod", &SinglePhaseFlowGDGD<KT_, SET_>::getViscosityMethod)
        .def("setViscosityMethod", &SinglePhaseFlowGDGD<KT_, SET_>::setViscosityMethod)
        .def("setConstSmoothingLength",        &SinglePhaseFlowGDGD<KT_, SET_>::setConstSmoothingLength)
        .def("computeSolidForces",             &SinglePhaseFlowGDGD<KT_, SET_>::computeSolidForces)
        .def("activateArtificialViscosity",    &SinglePhaseFlowGDGD<KT_, SET_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity",  &SinglePhaseFlowGDGD<KT_, SET_>::deactivateArtificialViscosity)
        .def("activateDensityDiffusion",       &SinglePhaseFlowGDGD<KT_, SET_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion",     &SinglePhaseFlowGDGD<KT_, SET_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &SinglePhaseFlowGDGD<KT_, SET_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization",&SinglePhaseFlowGDGD<KT_, SET_>::deactivateShepardRenormalization)
        .def("activateDensityReinitialization",&SinglePhaseFlowGDGD<KT_, SET_>::activateDensityReinitialization)
        .def("deactivateDensityReinitialization",&SinglePhaseFlowGDGD<KT_, SET_>::deactivateDensityReinitialization)
        .def("activatePowerLaw",               &SinglePhaseFlowGDGD<KT_, SET_>::activatePowerLaw,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("mu_min") = Scalar(0))
        .def("activateCarreau",                &SinglePhaseFlowGDGD<KT_, SET_>::activateCarreau)
        .def("activateBingham",                &SinglePhaseFlowGDGD<KT_, SET_>::activateBingham,
             pybind11::arg("mu_p"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("activateHerschelBulkley",        &SinglePhaseFlowGDGD<KT_, SET_>::activateHerschelBulkley,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("deactivateNonNewtonian",         &SinglePhaseFlowGDGD<KT_, SET_>::deactivateNonNewtonian)
        .def("setAcceleration",  &SPHBaseClass<KT_, SET_>::setAcceleration)
        .def("setRCut",          &SinglePhaseFlowGDGD<KT_, SET_>::setRCutPython)
        ;
    }

} // end namespace detail


// ─────────────────────────────────────────────────────────────────────────────
// Explicit template instantiations (5 kernels × 2 EOS = 10 variants)
// ─────────────────────────────────────────────────────────────────────────────

template class PYBIND11_EXPORT SinglePhaseFlowGDGD<wendlandc2, linear>;
template class PYBIND11_EXPORT SinglePhaseFlowGDGD<wendlandc2, tait>;
template class PYBIND11_EXPORT SinglePhaseFlowGDGD<wendlandc4, linear>;
template class PYBIND11_EXPORT SinglePhaseFlowGDGD<wendlandc4, tait>;
template class PYBIND11_EXPORT SinglePhaseFlowGDGD<wendlandc6, linear>;
template class PYBIND11_EXPORT SinglePhaseFlowGDGD<wendlandc6, tait>;
template class PYBIND11_EXPORT SinglePhaseFlowGDGD<quintic,    linear>;
template class PYBIND11_EXPORT SinglePhaseFlowGDGD<quintic,    tait>;
template class PYBIND11_EXPORT SinglePhaseFlowGDGD<cubicspline,linear>;
template class PYBIND11_EXPORT SinglePhaseFlowGDGD<cubicspline,tait>;


namespace detail
{
    template void export_SinglePhaseFlowGDGD<wendlandc2, linear>(pybind11::module& m, std::string name = "SinglePFGDGD_WC2_L");
    template void export_SinglePhaseFlowGDGD<wendlandc2, tait>  (pybind11::module& m, std::string name = "SinglePFGDGD_WC2_T");
    template void export_SinglePhaseFlowGDGD<wendlandc4, linear>(pybind11::module& m, std::string name = "SinglePFGDGD_WC4_L");
    template void export_SinglePhaseFlowGDGD<wendlandc4, tait>  (pybind11::module& m, std::string name = "SinglePFGDGD_WC4_T");
    template void export_SinglePhaseFlowGDGD<wendlandc6, linear>(pybind11::module& m, std::string name = "SinglePFGDGD_WC6_L");
    template void export_SinglePhaseFlowGDGD<wendlandc6, tait>  (pybind11::module& m, std::string name = "SinglePFGDGD_WC6_T");
    template void export_SinglePhaseFlowGDGD<quintic,    linear>(pybind11::module& m, std::string name = "SinglePFGDGD_Q_L");
    template void export_SinglePhaseFlowGDGD<quintic,    tait>  (pybind11::module& m, std::string name = "SinglePFGDGD_Q_T");
    template void export_SinglePhaseFlowGDGD<cubicspline,linear>(pybind11::module& m, std::string name = "SinglePFGDGD_CS_L");
    template void export_SinglePhaseFlowGDGD<cubicspline,tait>  (pybind11::module& m, std::string name = "SinglePFGDGD_CS_T");
} // end namespace detail

} // end namespace sph
} // end namespace hoomd

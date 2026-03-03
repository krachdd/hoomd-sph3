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

#include "TwoPhaseFlowTV.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

using namespace std;

namespace hoomd
{
namespace sph
{

/*! Constructor — delegates entirely to TwoPhaseFlow base constructor.
 *  m_tv_buf is sized lazily on the first computeForces() call.
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
TwoPhaseFlowTV<KT_, SET1_, SET2_>::TwoPhaseFlowTV(
                                 std::shared_ptr<SystemDefinition> sysdef,
                                 std::shared_ptr<SmoothingKernel<KT_> > skernel,
                                 std::shared_ptr<StateEquation<SET1_> > equationofstate1,
                                 std::shared_ptr<StateEquation<SET2_> > equationofstate2,
                                 std::shared_ptr<nsearch::NeighborList> nlist,
                                 std::shared_ptr<ParticleGroup> fluidgroup1,
                                 std::shared_ptr<ParticleGroup> fluidgroup2,
                                 std::shared_ptr<ParticleGroup> solidgroup,
                                 DensityMethod mdensitymethod,
                                 ViscosityMethod mviscositymethod,
                                 ColorGradientMethod mcolorgradientmethod)
    : TwoPhaseFlow<KT_, SET1_, SET2_>(sysdef, skernel, equationofstate1, equationofstate2,
                                      nlist, fluidgroup1, fluidgroup2, solidgroup,
                                      mdensitymethod, mviscositymethod, mcolorgradientmethod)
    {
    this->m_exec_conf->msg->notice(5) << "Constructing TwoPhaseFlowTV" << std::endl;
    // m_tv_buf is sized lazily in computeForces()
    }

/*! Destructor
*/
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
TwoPhaseFlowTV<KT_, SET1_, SET2_>::~TwoPhaseFlowTV()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying TwoPhaseFlowTV" << std::endl;
    }


/*! forcecomputation — modified copy of TwoPhaseFlow::forcecomputation().
 *
 *  In addition to the standard two-phase pressure / viscosity / surface-force
 *  terms, this function adds:
 *
 *  (a) Artificial-stress force  (Adami et al. 2013, Eq. 11):
 *        F_i^AS += (1/2) Σ_j (Vi²+Vj²) (A_i + A_j)·∇W_ij
 *        A_k = ρ_k v_k ⊗ (tv_k − v_k)   (rank-2 tensor)
 *      Penalises deviation of the transport velocity from the physical velocity.
 *      Counters the tensile instability (particle clustering in tension regions).
 *
 *  (b) Background-pressure contribution (BPC, written to aux2):
 *        bpc_i += −Σ_j (Vi²+Vj²) P_bg / m_i · (∂W/∂r)/r · r_ij
 *      P_bg = EOS::getTransportVelocityPressure() (uniform positive constant).
 *      KickDriftKickTV reads bpc_i on the next half-step to advect particles
 *      with the transport velocity instead of the physical velocity.
 *
 *  All inherited options remain active here:
 *    - Consistent interface pressure (Hu & Adams 2009)
 *    - Riemann-based dissipation (Zhang, Hu & Adams 2017)
 *    - Corrected Molteni–Colagrossi density diffusion
 *
 *  Pre-conditions set by computeForces() before this is called:
 *    aux1 = fictitious solid velocities   (from compute_noslip)
 *    aux2 = zero                          (zeroed in computeForces restore block)
 *    aux3 = transport velocity            (restored from m_tv_buf)
 *    aux4 = surface force density         (from compute_surfaceforce)
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlowTV<KT_, SET1_, SET2_>::forcecomputation(uint64_t timestep)
    {

    if ( this->m_density_method == DENSITYSUMMATION )
        this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlowTV::Forces using SUMMATION approach " << this->m_density_method << endl;
    else if ( this->m_density_method == DENSITYCONTINUITY )
        this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlowTV::Forces using CONTINUITY approach " << this->m_density_method << endl;

    // Grab handles for particle data
    ArrayHandle<Scalar4> h_force(this->m_force, access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_ratedpe(this->m_ratedpe, access_location::host, access_mode::readwrite);

    // Check input data
    assert(h_force.data);
    assert(h_ratedpe.data);

    // Zero data before force calculation
    memset((void*)h_force.data, 0, sizeof(Scalar4)*this->m_force.getNumElements());
    memset((void*)h_ratedpe.data, 0, sizeof(Scalar4)*this->m_ratedpe.getNumElements());

    // Access the particle data
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_density(this->m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar>  h_pressure(this->m_pdata->getPressures(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_vf(this->m_pdata->getAuxiliaries1(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_bpc(this->m_pdata->getAuxiliaries2(), access_location::host, access_mode::readwrite); // BPC (write)
    ArrayHandle<Scalar3> h_tv(this->m_pdata->getAuxiliaries3(), access_location::host, access_mode::read);       // TV (read)
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_sf(this->m_pdata->getAuxiliaries4(), access_location::host, access_mode::read);

    // Access the neighbor list
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Check input data
    assert(h_pos.data != NULL);

    unsigned int size;
    size_t myHead;

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Local variable to store things
    Scalar temp0 = 0;

    // Body force vector for consistent interface pressure correction (Hu & Adams 2009)
    const Scalar3 gvec = this->m_consistent_interface_pressure
                         ? this->getAcceleration(timestep)
                         : make_scalar3(Scalar(0), Scalar(0), Scalar(0));

    // Maximum velocity for adaptive timestep
    double max_vel = 0.0;

    // For each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

        // Access position, velocity, mass and type
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;

        Scalar3 vi;
        vi.x = h_velocity.data[i].x;
        vi.y = h_velocity.data[i].y;
        vi.z = h_velocity.data[i].z;
        Scalar mi = h_velocity.data[i].w;

        // Read particle i pressure
        Scalar Pi = h_pressure.data[i];

        // Read particle i density and volume
        Scalar rhoi = h_density.data[i];
        Scalar Vi   = mi / rhoi;

        // Read particle i type, viscosity, speed of sound and rest density
        bool i_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[i].w);
        Scalar mui   = i_isfluid1 ? this->m_mu1 : this->m_mu2;
        Scalar rho0i = i_isfluid1 ? this->m_rho01 : this->m_rho02;
        Scalar ci    = i_isfluid1 ? this->m_c1 : this->m_c2;

        // Properties needed for adaptive timestep
        Scalar vi_total = sqrt((vi.x * vi.x) + (vi.y * vi.y) + (vi.z * vi.z));
        if (i == 0) { max_vel = vi_total; }
        else if (vi_total > max_vel) { max_vel = vi_total; }

        // --- Transport velocity terms: particle i ---
        Scalar3 tvi = h_tv.data[i];

        // Background pressure for TV (phase-specific)
        Scalar P_tv_i = i_isfluid1 ? this->m_eos1->getTransportVelocityPressure()
                                   : this->m_eos2->getTransportVelocityPressure();

        // Artificial stress tensor A_i = rho_i * v_i ⊗ (tv_i - v_i)
        Scalar A11i = rhoi * vi.x * (tvi.x - vi.x);
        Scalar A12i = rhoi * vi.x * (tvi.y - vi.y);
        Scalar A13i = rhoi * vi.x * (tvi.z - vi.z);
        Scalar A21i = rhoi * vi.y * (tvi.x - vi.x);
        Scalar A22i = rhoi * vi.y * (tvi.y - vi.y);
        Scalar A23i = rhoi * vi.y * (tvi.z - vi.z);
        Scalar A31i = rhoi * vi.z * (tvi.x - vi.x);
        Scalar A32i = rhoi * vi.z * (tvi.y - vi.y);
        Scalar A33i = rhoi * vi.z * (tvi.z - vi.z);

        // Zero BPC accumulator for particle i
        h_bpc.data[i] = make_scalar3(Scalar(0), Scalar(0), Scalar(0));

        // Loop over all neighbours of particle i
        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {
            unsigned int k = h_nlist.data[myHead + j];

            // Sanity check
            assert(k < this->m_pdata->getN() + this->m_pdata->getNGhosts());

            // Access neighbour position
            Scalar3 pj;
            pj.x = h_pos.data[k].x;
            pj.y = h_pos.data[k].y;
            pj.z = h_pos.data[k].z;

            // Determine neighbour type
            bool j_issolid  = checksolid(h_type_property_map.data, h_pos.data[k].w);
            bool j_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[k].w);

            // Read particle j viscosity, speed of sound and rest density
            Scalar muj   = j_isfluid1 ? this->m_mu1 : this->m_mu2;
            Scalar rho0j = j_isfluid1 ? this->m_rho01 : this->m_rho02;
            Scalar cj    = j_isfluid1 ? this->m_c1 : this->m_c2;
            // For solid neighbours, use particle i properties
            muj   = j_issolid ? mui : muj;
            rho0j = j_issolid ? rho0i : rho0j;
            cj    = j_issolid ? ci : cj;

            // Compute distance vector
            Scalar3 dx;
            dx.x = pi.x - pj.x;
            dx.y = pi.y - pj.y;
            dx.z = pi.z - pj.z;

            // Apply periodic boundary conditions
            dx = box.minImage(dx);

            // Calculate squared distance
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, skip
            if ( this->m_const_slength && rsq > this->m_rcutsq )
                continue;

            // Access neighbour velocity; fictitious for solid, physical for fluid
            Scalar3 vj = make_scalar3(0.0, 0.0, 0.0);
            Scalar mj  = h_velocity.data[k].w;
            if ( j_issolid )
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

            // Read neighbour pressure
            Scalar Pj = h_pressure.data[k];

            // Compute velocity difference
            Scalar3 dv;
            dv.x = vi.x - vj.x;
            dv.y = vi.y - vj.y;
            dv.z = vi.z - vj.z;

            // Calculate absolute distance
            Scalar r = sqrt(rsq);

            // Mean smoothing length and denominator modifier
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar eps    = Scalar(0.1)*meanh;
            Scalar epssqr = eps*eps;

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh, r);
            Scalar dwdr_r = dwdr / (r + eps);

            // ── Inter-particle pressure force ────────────────────────────────────
            // Symmetric volume formulation (Adami et al. 2013):
            //   F_i^p = −Σ_j (Vi² + Vj²) · p̄_ij · (∂W/∂r / r) · r_ij
            //
            // DENSITYSUMMATION — density-weighted average pressure:
            //   p̄_ij = (ρ_j·p_i + ρ_i·p_j) / (ρ_i + ρ_j)
            //
            //   With consistent interface pressure (CIP, Hu & Adams 2009), cross-phase
            //   pairs use rest-density weighting + hydrostatic correction:
            //     p̄_ij = (ρ₀ⱼ·p_i + ρ₀ᵢ·p_j + ρ₀ᵢ ρ₀ⱼ (g·r_ij)) / (ρ₀ᵢ + ρ₀ⱼ)
            //   Suppresses parasitic currents at large-density-ratio interfaces.
            //
            // DENSITYCONTINUITY — mass-flux-consistent form:
            //   prefactor = m_i m_j;   p̄_ij = (p_i + p_j) / (ρ_i ρ_j)
            Scalar prefactor = 0.0;
            if ( this->m_density_method == DENSITYSUMMATION )
                {
                // Consistent interface pressure (Hu & Adams 2009): rest-density weighting
                // + hydrostatic correction for cross-phase pairs; standard formula otherwise
                if ( this->m_consistent_interface_pressure && !j_issolid && (i_isfluid1 != j_isfluid1) )
                    temp0 = (rho0j*Pi + rho0i*Pj + rho0i*rho0j*dot(gvec, dx)) / (rho0i + rho0j);
                else
                    temp0 = (rhoj*Pi + rhoi*Pj) / (rhoi + rhoj);
                prefactor = Vi*Vi + Vj*Vj;
                }
            else if ( this->m_density_method == DENSITYCONTINUITY )
                {
                temp0 = (Pi + Pj) / (rhoi*rhoj);
                prefactor = mi * mj;
                }

            // ── Momentum dissipation (fluid–fluid pairs only) ────────────────────
            // Exactly one branch is active at a time (else-if).
            //
            // [A] Monaghan artificial viscosity (Monaghan 1992):
            //     Π_ij = (−α c_max μ_ij + β μ_ij²) / ρ̄_ij
            //     μ_ij = h̄ (v_ij · r_ij) / (r_ij² + η²)
            //   Activated via activateArtificialViscosity(alpha, beta).
            //
            // [B] Riemann-based dissipation (Zhang, Hu & Adams 2017):
            //     Z*_ij = Z_i Z_j / (Z_i + Z_j),  Z = ρ c   [harmonic mean impedance]
            //     u_ij  = (v_ij · r_ij) / (|r_ij| + η)       [signed radial velocity]
            //     avc   = −β_R · Z*_ij · u_ij⁻ / ρ̄_ij        (only if v_ij·r_ij < 0)
            //   Activated via activateRiemannDissipation(beta).
            Scalar avc = 0.0;
            // [A] Monaghan AV — Monaghan (1992) Annu. Rev. Astron. Astrophys. 30, 543–574
            if ( this->m_artificial_viscosity && !j_issolid )
                {
                Scalar dotdvdx = dot(dv, dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar muij    = meanh*dotdvdx / (rsq + epssqr);
                    Scalar meanrho = Scalar(0.5)*(rhoi + rhoj);
                    avc = (-this->m_avalpha*this->m_cmax*muij + this->m_avbeta*muij*muij) / meanrho;
                    }
                }
            // [B] Riemann dissipation — Zhang, Hu & Adams (2017) J. Comput. Phys. 340, 439–455
            else if ( this->m_riemann_dissipation && !j_issolid )
                {
                Scalar dotdvdx = dot(dv, dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar uij   = dotdvdx / (r + eps);
                    Scalar Zi    = rhoi * ci;
                    Scalar Zj    = rhoj * cj;
                    Scalar Zstar = (Zi * Zj) / (Zi + Zj);
                    Scalar meanrho = Scalar(0.5) * (rhoi + rhoj);
                    avc = -this->m_riemann_beta * Zstar * uij / meanrho;
                    }
                }

            // Add pressure + dissipation force contribution to fluid particle
            h_force.data[i].x -= prefactor * (temp0 + avc) * dwdr_r * dx.x;
            h_force.data[i].y -= prefactor * (temp0 + avc) * dwdr_r * dx.y;
            h_force.data[i].z -= prefactor * (temp0 + avc) * dwdr_r * dx.z;

            // Evaluate viscous interaction forces
            {
            Scalar dvnorm    = sqrt(dot(dv, dv));
            Scalar gamma_dot = dvnorm / (r + eps);
            NonNewtonianModel nn_model_i = i_isfluid1 ? this->m_nn_model1 : this->m_nn_model2;
            Scalar mu_eff_i = computeNNViscosity(mui, gamma_dot, nn_model_i,
                i_isfluid1 ? this->m_nn_K1 : this->m_nn_K2,
                i_isfluid1 ? this->m_nn_n1 : this->m_nn_n2,
                i_isfluid1 ? this->m_nn_mu0_1 : this->m_nn_mu0_2,
                i_isfluid1 ? this->m_nn_muinf_1 : this->m_nn_muinf_2,
                i_isfluid1 ? this->m_nn_lambda1 : this->m_nn_lambda2,
                i_isfluid1 ? this->m_nn_tauy1 : this->m_nn_tauy2,
                i_isfluid1 ? this->m_nn_m1 : this->m_nn_m2,
                i_isfluid1 ? this->m_nn_mu_min1 : this->m_nn_mu_min2);
            Scalar mu_eff_j;
            if (j_issolid)
                mu_eff_j = mu_eff_i;
            else
                {
                NonNewtonianModel nn_model_j = j_isfluid1 ? this->m_nn_model1 : this->m_nn_model2;
                mu_eff_j = computeNNViscosity(muj, gamma_dot, nn_model_j,
                    j_isfluid1 ? this->m_nn_K1 : this->m_nn_K2,
                    j_isfluid1 ? this->m_nn_n1 : this->m_nn_n2,
                    j_isfluid1 ? this->m_nn_mu0_1 : this->m_nn_mu0_2,
                    j_isfluid1 ? this->m_nn_muinf_1 : this->m_nn_muinf_2,
                    j_isfluid1 ? this->m_nn_lambda1 : this->m_nn_lambda2,
                    j_isfluid1 ? this->m_nn_tauy1 : this->m_nn_tauy2,
                    j_isfluid1 ? this->m_nn_m1 : this->m_nn_m2,
                    j_isfluid1 ? this->m_nn_mu_min1 : this->m_nn_mu_min2);
                }
            Scalar mu_harm = Scalar(2) * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
            temp0 = mu_harm * (Vi*Vi+Vj*Vj) * dwdr_r;
            }
            h_force.data[i].x += temp0 * dv.x;
            h_force.data[i].y += temp0 * dv.y;
            h_force.data[i].z += temp0 * dv.z;

            // ── Transport-velocity terms (fluid–fluid pairs only) ───────────────
            // Applied only between two fluid particles; solid boundaries use the
            // physical velocity (no fictitious TV) so the j_issolid guard is
            // required for both the artificial-stress and BPC contributions.
            if ( !j_issolid )
                {
                Scalar3 tvk = h_tv.data[k];
                Scalar3 vk  = make_scalar3(h_velocity.data[k].x,
                                           h_velocity.data[k].y,
                                           h_velocity.data[k].z);

                // Artificial-stress tensor for particle k:
                //   A_k = ρ_k v_k ⊗ (tv_k − v_k)   [rank-2 outer product]
                // A is non-zero only where tv ≠ v, i.e. at locations with
                // local velocity divergence that drives particle clustering.
                Scalar A11k = rhoj * vk.x * (tvk.x - vk.x);
                Scalar A12k = rhoj * vk.x * (tvk.y - vk.y);
                Scalar A13k = rhoj * vk.x * (tvk.z - vk.z);
                Scalar A21k = rhoj * vk.y * (tvk.x - vk.x);
                Scalar A22k = rhoj * vk.y * (tvk.y - vk.y);
                Scalar A23k = rhoj * vk.y * (tvk.z - vk.z);
                Scalar A31k = rhoj * vk.z * (tvk.x - vk.x);
                Scalar A32k = rhoj * vk.z * (tvk.y - vk.y);
                Scalar A33k = rhoj * vk.z * (tvk.z - vk.z);

                Scalar vijsqr = Vi*Vi + Vj*Vj;

                // Artificial-stress force (Adami, Hu & Adams 2013, Eq. 11):
                //   F_i^AS += (1/2) Σ_j (Vi²+Vj²) (A_i + A_j)·∇W_ij
                // (A_i+A_j) is contracted with ∇W_ij = (∂W/∂r / r) · r_ij.
                // Factor 1/2 absorbed into tv_temp = 0.5*(Vi²+Vj²)*(∂W/∂r / r).
                // The force opposes particle clustering in tension regions
                // (tensile instability), restoring a uniform particle distribution.
                Scalar tv_temp = Scalar(0.5) * vijsqr * dwdr_r;
                Scalar A1ij = (A11i+A11k)*dx.x + (A12i+A12k)*dx.y + (A13i+A13k)*dx.z;
                Scalar A2ij = (A21i+A21k)*dx.x + (A22i+A22k)*dx.y + (A23i+A23k)*dx.z;
                Scalar A3ij = (A31i+A31k)*dx.x + (A32i+A32k)*dx.y + (A33i+A33k)*dx.z;
                h_force.data[i].x += tv_temp * A1ij;
                h_force.data[i].y += tv_temp * A2ij;
                h_force.data[i].z += tv_temp * A3ij;

                // Background-pressure contribution (BPC), accumulated in aux2:
                //   bpc_i += −Σ_j (Vi²+Vj²) P_bg / m_i · (∂W/∂r / r) · r_ij
                // P_bg = EOS::getTransportVelocityPressure() is a positive constant
                // (typically ~ ρ₀ c²) that keeps all particle pressures positive so
                // particles can be advected along the smooth transport velocity field.
                // KickDriftKickTV reads aux2 on the next half-step to compute the
                // TV advection increment; it must be zeroed before this call (done
                // by the restore block in computeForces()).
                h_bpc.data[i].x -= vijsqr * P_tv_i / mi * dwdr_r * dx.x;
                h_bpc.data[i].y -= vijsqr * P_tv_i / mi * dwdr_r * dx.y;
                h_bpc.data[i].z -= vijsqr * P_tv_i / mi * dwdr_r * dx.z;
                }

            // Evaluate rate of change of density if CONTINUITY approach is used
            if ( this->m_density_method == DENSITYCONTINUITY )
                {
                if ( j_issolid )
                    {
                    // Use physical advection velocity for solid neighbours
                    vj.x = h_velocity.data[k].x;
                    vj.y = h_velocity.data[k].y;
                    vj.z = h_velocity.data[k].z;

                    dv.x = vi.x - vj.x;
                    dv.y = vi.y - vj.y;
                    dv.z = vi.z - vj.z;
                    }

                h_ratedpe.data[i].x += rhoi*Vj*dot(dv, dwdr_r*dx);

                // Molteni–Colagrossi density diffusion (fluid–fluid pairs only).
                // Ref: Molteni & Colagrossi (2009) Comput. Phys. Commun. 180, 861–872.
                //
                // Drive term is (ρ_i/ρ₀ᵢ − ρ_j/ρ₀ⱼ) — rest-density normalised.
                // The original term (ρ_i/ρ_j − 1) is non-zero at equilibrium when
                // ρ₀₁ ≠ ρ₀₂ (different-phase rest densities), generating unphysical
                // density drift across the interface in stratified-flow setups.
                // The normalised form equals zero at equilibrium for both phases.
                if ( !j_issolid && this->m_density_diffusion )
                    h_ratedpe.data[i].x -= (Scalar(2)*this->m_ddiff*meanh*this->m_cmax*mj*(rhoi/rho0i-rhoj/rho0j)*dot(dx,dwdr_r*dx))/(rsq+epssqr);
                }

            } // Closing Neighbour Loop

        // dp/dt = (dp/dρ) * dρ/dt via chain rule (DENSITYCONTINUITY only)
        if ( this->m_density_method == DENSITYCONTINUITY )
            {
            Scalar dpdrho_i = i_isfluid1 ? this->m_eos1->dPressuredDensity(rhoi)
                                         : this->m_eos2->dPressuredDensity(rhoi);
            h_ratedpe.data[i].y = dpdrho_i * h_ratedpe.data[i].x;
            }

        // Add surface force
        h_force.data[i].x += h_sf.data[i].x;
        h_force.data[i].y += h_sf.data[i].y;
        h_force.data[i].z += h_sf.data[i].z;

        } // Closing Fluid Particle Loop

    this->m_timestep_list[5] = max_vel;
    // Add volumetric force (gravity)
    this->applyBodyForce(timestep, this->m_fluidgroup);

    } // end forcecomputation


/*! computeForces — full override of TwoPhaseFlow::computeForces() with
 *  TV save/restore around compute_colorgradients/compute_surfaceforce, and
 *  a post-forcecomputation ghost update for BPC and TV.
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlowTV<KT_, SET1_, SET2_>::computeForces(uint64_t timestep)
    {

    // Start by updating the neighbor list
    this->m_nlist->compute(timestep);

    // Parameters must be set before run()
    if (!this->m_params_set)
        {
        this->m_exec_conf->msg->error() << "sph.models.TwoPhaseFlowTV requires parameters to be set before run()"
            << std::endl;
        throw std::runtime_error("Error computing TwoPhaseFlowTV forces");
        }

    // Mark solid particles for removal once at timestep 0
    if (!this->m_solid_removed)
        {
        this->m_nlist->forceUpdate();
        this->m_nlist->compute(timestep);
        this->mark_solid_particles_toremove(timestep);
        this->m_solid_removed = true;
        }

    // Shepard density renormalization
    if ( this->m_shepard_renormalization && timestep % this->m_shepardfreq == 0 )
        {
        this->renormalize_density(timestep);
        this->compute_pressure(timestep);
        this->m_pressure_initialized = true;
#ifdef ENABLE_MPI
        this->update_ghost_density_pressure(timestep);
#endif
        }

    // Periodic density reinitialization from summation (DENSITYCONTINUITY only)
    if ( this->m_density_reinitialization && timestep % this->m_densityreinitfreq == 0 )
        {
        this->compute_ndensity(timestep);
        this->compute_pressure(timestep);
        this->m_pressure_initialized = true;
#ifdef ENABLE_MPI
        this->update_ghost_density_pressure(timestep);
#endif
        }

    // Fickian shifting: update concentration gradient
    if ( this->m_fickian_shifting )
        {
        this->compute_particle_concentration_gradient(timestep);
#ifdef ENABLE_MPI
        this->update_ghost_density_pressure_energy(timestep);
#endif
        }

    // Density and pressure computation
    if ( this->m_density_method == DENSITYSUMMATION )
        {
        this->compute_ndensity(timestep);
        this->compute_pressure(timestep);
        }
    else // DENSITYCONTINUITY
        {
        if ( !this->m_pressure_initialized )
            {
            this->compute_pressure(timestep);
            this->m_pressure_initialized = true;
            }
        }

#ifdef ENABLE_MPI
    this->update_ghost_density_pressure_energy(timestep);
#endif

    // Compute fictitious solid particle properties (no-slip)
    this->compute_noslip(timestep);

#ifdef ENABLE_MPI
    this->update_ghost_density_pressure(timestep);
#endif

    // ── Save transport velocity (aux3) before compute_colorgradients() ──
    // compute_colorgradients() will overwrite aux3 with fluid-fluid normals.
    // We save both local and ghost particles so the restore is complete.
    {
    unsigned int n_total = this->m_pdata->getN() + this->m_pdata->getNGhosts();
    m_tv_buf.resize(n_total);
    ArrayHandle<Scalar3> h_tv_save(this->m_pdata->getAuxiliaries3(),
                                   access_location::host, access_mode::read);
    memcpy(m_tv_buf.data(), h_tv_save.data, n_total * sizeof(Scalar3));
    }

    // Compute colour gradients (writes normals into aux2 and aux3)
    this->compute_colorgradients(timestep);

#ifdef ENABLE_MPI
    // Communicate normals to ghost particles
    this->update_ghost_aux123(timestep);
#endif

    // δ⁺-SPH particle shifting (Sun et al. 2017).
    // Uses aux3 (fluid-fluid normals) which are correctly in place here.
    if ( this->m_particle_shifting )
        {
        this->compute_particle_shift(timestep);
        this->m_nlist->forceUpdate();
        this->m_nlist->compute(timestep);
        }

    // Compute surface force from normals in aux2/aux3 → aux4
    this->compute_surfaceforce(timestep);

#ifdef ENABLE_MPI
    this->update_ghost_aux4(timestep);
#endif

    // ── Restore TV to aux3; zero aux2 to receive BPC ──
    {
    unsigned int n_total = this->m_pdata->getN() + this->m_pdata->getNGhosts();
    ArrayHandle<Scalar3> h_tv_rest(this->m_pdata->getAuxiliaries3(),
                                   access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_bpc_zero(this->m_pdata->getAuxiliaries2(),
                                    access_location::host, access_mode::readwrite);
    memcpy(h_tv_rest.data, m_tv_buf.data(), n_total * sizeof(Scalar3));
    memset(h_bpc_zero.data, 0, n_total * sizeof(Scalar3));
    }

    // Force computation: reads aux1 (fict. vel), aux3 (TV), aux4 (surface force);
    // writes aux2 (BPC) for each local fluid particle.
    forcecomputation(timestep);

#ifdef ENABLE_MPI
    // Communicate aux1 (fict. vel), aux2 (BPC), aux3 (TV) to ghost particles
    // so that KickDriftKickTV can read them on all ranks.
    this->update_ghost_aux123(timestep);
#endif

    if ( this->m_compute_solid_forces )
        {
        this->compute_solid_forces(timestep);
        }

    } // End computeForces


namespace detail
{

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void export_TwoPhaseFlowTV(pybind11::module& m, std::string name)
{
    pybind11::class_<TwoPhaseFlowTV<KT_, SET1_, SET2_>,
                     TwoPhaseFlow<KT_, SET1_, SET2_>,
                     std::shared_ptr<TwoPhaseFlowTV<KT_, SET1_, SET2_>>>(m, name.c_str())
        .def(pybind11::init< std::shared_ptr<SystemDefinition>,
                             std::shared_ptr<SmoothingKernel<KT_> >,
                             std::shared_ptr<StateEquation<SET1_> >,
                             std::shared_ptr<StateEquation<SET2_> >,
                             std::shared_ptr<nsearch::NeighborList>,
                             std::shared_ptr<ParticleGroup>,
                             std::shared_ptr<ParticleGroup>,
                             std::shared_ptr<ParticleGroup>,
                             DensityMethod,
                             ViscosityMethod,
                             ColorGradientMethod >())
        .def("setParams", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::setParams)
        .def("getDensityMethod", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::getDensityMethod)
        .def("setDensityMethod", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::setDensityMethod)
        .def("getViscosityMethod", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::getViscosityMethod)
        .def("setViscosityMethod", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::setViscosityMethod)
        .def("getColorGradientMethod", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::getColorGradientMethod)
        .def("setColorGradientMethod", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::setColorGradientMethod)
        .def("setConstSmoothingLength", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::setConstSmoothingLength)
        .def("computeSolidForces", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::computeSolidForces)
        .def("activateArtificialViscosity", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::deactivateArtificialViscosity)
        .def("activateConsistentInterfacePressure", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateConsistentInterfacePressure)
        .def("deactivateConsistentInterfacePressure", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::deactivateConsistentInterfacePressure)
        .def("activateRiemannDissipation", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateRiemannDissipation,
             pybind11::arg("beta") = Scalar(1.0))
        .def("deactivateRiemannDissipation", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::deactivateRiemannDissipation)
        .def("activateDensityDiffusion", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::deactivateShepardRenormalization)
        .def("activateDensityReinitialization", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateDensityReinitialization)
        .def("deactivateDensityReinitialization", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::deactivateDensityReinitialization)
        .def("activateFickianShifting", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateFickianShifting)
        .def("deactivateFickianShifting", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::deactivateFickianShifting)
        .def("activateParticleShifting", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateParticleShifting,
             pybind11::arg("A")                   = Scalar(0.2),
             pybind11::arg("R")                   = Scalar(0.2),
             pybind11::arg("n")                   = 4,
             pybind11::arg("interface_condition") = true)
        .def("deactivateParticleShifting", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::deactivateParticleShifting)
        .def("getProvidedTimestepQuantities", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::getProvidedTimestepQuantities)
        .def("activatePowerLaw1", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activatePowerLaw1,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("mu_min") = Scalar(0))
        .def("activateCarreau1", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateCarreau1)
        .def("activateBingham1", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateBingham1,
             pybind11::arg("mu_p"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("activateHerschelBulkley1", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateHerschelBulkley1,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("deactivateNonNewtonian1", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::deactivateNonNewtonian1)
        .def("activatePowerLaw2", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activatePowerLaw2,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("mu_min") = Scalar(0))
        .def("activateCarreau2", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateCarreau2)
        .def("activateBingham2", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateBingham2,
             pybind11::arg("mu_p"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("activateHerschelBulkley2", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::activateHerschelBulkley2,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("deactivateNonNewtonian2", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::deactivateNonNewtonian2)
        .def("setAcceleration", &SPHBaseClass<KT_, SET1_>::setAcceleration)
        .def("setRCut", &TwoPhaseFlowTV<KT_, SET1_, SET2_>::setRCutPython)
        ;
}

} // end namespace detail

//! Explicit template instantiations
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc2, linear, linear>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc2, linear, tait>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc2, tait,   linear>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc2, tait,   tait>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc4, linear, linear>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc4, linear, tait>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc4, tait,   linear>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc4, tait,   tait>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc6, linear, linear>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc6, linear, tait>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc6, tait,   linear>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<wendlandc6, tait,   tait>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<quintic,    linear, linear>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<quintic,    linear, tait>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<quintic,    tait,   linear>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<quintic,    tait,   tait>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<cubicspline, linear, linear>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<cubicspline, linear, tait>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<cubicspline, tait,   linear>;
template class PYBIND11_EXPORT TwoPhaseFlowTV<cubicspline, tait,   tait>;

namespace detail
{

template void export_TwoPhaseFlowTV<wendlandc2, linear, linear>(pybind11::module& m, std::string name = "TwoPFTV_WC2_LL");
template void export_TwoPhaseFlowTV<wendlandc2, linear, tait  >(pybind11::module& m, std::string name = "TwoPFTV_WC2_LT");
template void export_TwoPhaseFlowTV<wendlandc2, tait,   linear>(pybind11::module& m, std::string name = "TwoPFTV_WC2_TL");
template void export_TwoPhaseFlowTV<wendlandc2, tait,   tait  >(pybind11::module& m, std::string name = "TwoPFTV_WC2_TT");
template void export_TwoPhaseFlowTV<wendlandc4, linear, linear>(pybind11::module& m, std::string name = "TwoPFTV_WC4_LL");
template void export_TwoPhaseFlowTV<wendlandc4, linear, tait  >(pybind11::module& m, std::string name = "TwoPFTV_WC4_LT");
template void export_TwoPhaseFlowTV<wendlandc4, tait,   linear>(pybind11::module& m, std::string name = "TwoPFTV_WC4_TL");
template void export_TwoPhaseFlowTV<wendlandc4, tait,   tait  >(pybind11::module& m, std::string name = "TwoPFTV_WC4_TT");
template void export_TwoPhaseFlowTV<wendlandc6, linear, linear>(pybind11::module& m, std::string name = "TwoPFTV_WC6_LL");
template void export_TwoPhaseFlowTV<wendlandc6, linear, tait  >(pybind11::module& m, std::string name = "TwoPFTV_WC6_LT");
template void export_TwoPhaseFlowTV<wendlandc6, tait,   linear>(pybind11::module& m, std::string name = "TwoPFTV_WC6_TL");
template void export_TwoPhaseFlowTV<wendlandc6, tait,   tait  >(pybind11::module& m, std::string name = "TwoPFTV_WC6_TT");
template void export_TwoPhaseFlowTV<quintic,    linear, linear>(pybind11::module& m, std::string name = "TwoPFTV_Q_LL");
template void export_TwoPhaseFlowTV<quintic,    linear, tait  >(pybind11::module& m, std::string name = "TwoPFTV_Q_LT");
template void export_TwoPhaseFlowTV<quintic,    tait,   linear>(pybind11::module& m, std::string name = "TwoPFTV_Q_TL");
template void export_TwoPhaseFlowTV<quintic,    tait,   tait  >(pybind11::module& m, std::string name = "TwoPFTV_Q_TT");
template void export_TwoPhaseFlowTV<cubicspline, linear, linear>(pybind11::module& m, std::string name = "TwoPFTV_CS_LL");
template void export_TwoPhaseFlowTV<cubicspline, linear, tait  >(pybind11::module& m, std::string name = "TwoPFTV_CS_LT");
template void export_TwoPhaseFlowTV<cubicspline, tait,   linear>(pybind11::module& m, std::string name = "TwoPFTV_CS_TL");
template void export_TwoPhaseFlowTV<cubicspline, tait,   tait  >(pybind11::module& m, std::string name = "TwoPFTV_CS_TT");

} // end namespace detail
} // end namespace sph
} // end namespace hoomd

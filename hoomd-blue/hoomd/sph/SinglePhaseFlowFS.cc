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


#include "SinglePhaseFlowFS.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include <cmath>

using namespace std;

namespace hoomd
{
namespace sph
{

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

template<SmoothingKernelType KT_, StateEquationType SET_>
SinglePhaseFlowFS<KT_, SET_>::SinglePhaseFlowFS(
        std::shared_ptr<SystemDefinition> sysdef,
        std::shared_ptr<SmoothingKernel<KT_>> skernel,
        std::shared_ptr<StateEquation<SET_>> equationofstate,
        std::shared_ptr<nsearch::NeighborList> nlist,
        std::shared_ptr<ParticleGroup> fluidgroup,
        std::shared_ptr<ParticleGroup> solidgroup,
        DensityMethod   mdensitymethod,
        ViscosityMethod mviscositymethod)
    : SinglePhaseFlowTV<KT_, SET_>(sysdef, skernel, equationofstate, nlist,
                                    fluidgroup, solidgroup,
                                    mdensitymethod, mviscositymethod)
    {
    this->m_exec_conf->msg->notice(5) << "Constructing SinglePhaseFlowFS" << std::endl;

    // Initialise inherited flags to safe defaults (same as SinglePhaseFlowTV).
    this->m_const_slength           = false;
    this->m_params_set              = false;
    this->m_compute_solid_forces    = false;
    this->m_artificial_viscosity    = false;
    this->m_density_diffusion       = false;
    this->m_shepard_renormalization = false;
    this->m_density_reinitialization= false;
    this->m_ch     = Scalar(0.0);
    this->m_rcut   = Scalar(0.0);
    this->m_rcutsq = Scalar(0.0);
    this->m_avalpha = Scalar(0.0);
    this->m_avbeta  = Scalar(0.0);
    this->m_ddiff   = Scalar(0.0);
    this->m_shepardfreq       = 1;
    this->m_densityreinitfreq = 1;

    this->m_solid_removed        = false;
    this->m_pressure_initialized = false;

    // Sanity checks
    assert(this->m_pdata);
    assert(this->m_nlist);
    assert(this->m_skernel);
    assert(this->m_eos);

    // Build fluid / solid type vectors and type-property map.
    this->constructTypeVectors(fluidgroup, &this->m_fluidtypes);
    this->constructTypeVectors(solidgroup,  &this->m_solidtypes);

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

    // Store density/viscosity methods and cache EOS/kernel scalars.
    this->m_density_method   = mdensitymethod;
    this->m_viscosity_method = mviscositymethod;

    this->m_rho0  = equationofstate->getRestDensity();
    this->m_c     = equationofstate->getSpeedOfSound();
    this->m_kappa = skernel->getKernelKappa();

    // Register cut-off matrix with the neighbour list.
    this->m_r_cut_nlist = std::make_shared<GPUArray<Scalar>>(
        this->m_typpair_idx.getNumElements(), this->m_exec_conf);
    this->m_nlist->addRCutMatrix(this->m_r_cut_nlist);

    // FS-specific defaults (safe: no surface tension, neutral wetting).
    m_sigma         = Scalar(0.0);
    m_fs_threshold  = Scalar(0.75);
    m_contact_angle = Scalar(M_PI_2);   // 90° = neutral wetting

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
SinglePhaseFlowFS<KT_, SET_>::~SinglePhaseFlowFS()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying SinglePhaseFlowFS" << std::endl;
    }


// ─────────────────────────────────────────────────────────────────────────────
// Parameter setter
// ─────────────────────────────────────────────────────────────────────────────

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowFS<KT_, SET_>::setFSParams(
        Scalar sigma, Scalar fs_threshold, Scalar contact_angle)
    {
    if (sigma < Scalar(0))
        {
        this->m_exec_conf->msg->error()
            << "SinglePhaseFlowFS: sigma must be >= 0" << std::endl;
        throw std::invalid_argument("sigma must be >= 0");
        }
    if (fs_threshold <= Scalar(0) || fs_threshold >= Scalar(1))
        {
        this->m_exec_conf->msg->error()
            << "SinglePhaseFlowFS: fs_threshold must be in (0, 1)" << std::endl;
        throw std::invalid_argument("fs_threshold must be in (0, 1)");
        }
    m_sigma         = sigma;
    m_fs_threshold  = fs_threshold;
    m_contact_angle = contact_angle;
    }


// ─────────────────────────────────────────────────────────────────────────────
// MPI helper: sync fs normals (aux2) and lambda/curvature (aux4) to ghosts
// ─────────────────────────────────────────────────────────────────────────────

#ifdef ENABLE_MPI
template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowFS<KT_, SET_>::update_ghost_aux24(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7)
        << "Computing SinglePhaseFlowFS::update_ghost_aux24" << std::endl;

    if (this->m_comm)
        {
        CommFlags flags(0);
        flags[comm_flag::auxiliary2] = 1;  // fs outward normals
        flags[comm_flag::auxiliary4] = 1;  // lambda (x) + curvature (y)
        this->m_comm->setFlags(flags);
        this->m_comm->beginUpdateGhosts(timestep);
        this->m_comm->finishUpdateGhosts(timestep);
        }
    }
#endif // ENABLE_MPI


// ─────────────────────────────────────────────────────────────────────────────
// Step 1: detect_freesurface
// ─────────────────────────────────────────────────────────────────────────────

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowFS<KT_, SET_>::detect_freesurface(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7)
        << "Computing SinglePhaseFlowFS::detect_freesurface" << std::endl;

    const BoxDim& box = this->m_pdata->getGlobalBox();
    const unsigned int group_size = this->m_fluidgroup->getNumMembers();

    const Scalar cos_ca = std::cos(m_contact_angle);
    const Scalar sin_ca = std::sin(m_contact_angle);
    // Contact-angle correction is only applied when the angle deviates from 90°.
    const bool apply_ca = (std::fabs(m_contact_angle - Scalar(M_PI_2)) > Scalar(0.01));

    const Scalar eps_sq = Scalar(1e-12);

        { // GPU Array Scope
        ArrayHandle<Scalar4> h_pos     (this->m_pdata->getPositions(),    access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(),   access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_density (this->m_pdata->getDensities(),    access_location::host, access_mode::read);
        ArrayHandle<Scalar3> h_fs_n    (this->m_pdata->getAuxiliaries2(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar3> h_aux4    (this->m_pdata->getAuxiliaries4(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar>  h_h       (this->m_pdata->getSlengths(),     access_location::host, access_mode::read);

        ArrayHandle<unsigned int> h_n_neigh  (this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist_arr(this->m_nlist->getNListArray(),  access_location::host, access_mode::read);
        ArrayHandle<size_t>       h_head_list(this->m_nlist->getHeadList(),    access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

        // Use reference density $\rho_0$ (not the computed density) for the Shepard
        // sum volumes.  Computing $V_i = m_i / \rho_i$ with the *summation* density
        // creates a circular feedback: surface particles have low $\rho_i$, so their
        // $V_i$ is large, which inflates $\lambda$ above the threshold and prevents
        // free-surface detection.  Using $V = m/\rho_0$ gives the physically correct
        // completeness measure: $\lambda = \rho_\mathrm{fluid\,only} / \rho_0$.
        const Scalar rho0_ref = this->m_eos->getRestDensity();

        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

            Scalar3 pi;
            pi.x = h_pos.data[i].x;
            pi.y = h_pos.data[i].y;
            pi.z = h_pos.data[i].z;

            Scalar mi   = h_velocity.data[i].w;
            Scalar Vi   = mi / rho0_ref;
            Scalar hi   = this->m_const_slength ? this->m_ch : h_h.data[i];

            // Self contribution to Shepard sum: $V_i W(0, h_i)$
            Scalar W0 = this->m_skernel->wij(hi, Scalar(0));
            Scalar lambda = Vi * W0;

            // Gradient of lambda (from neighbours only; self-gradient is zero).
            Scalar3 grad_lambda = make_scalar3(Scalar(0), Scalar(0), Scalar(0));

            // Wall normal accumulator for contact-angle correction.
            Scalar3 n_wall_acc = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
            bool    has_solid_nbr = false;

            size_t       myHead = h_head_list.data[i];
            unsigned int size   = (unsigned int)h_n_neigh.data[i];

            for (unsigned int j = 0; j < size; j++)
                {
                unsigned int k = h_nlist_arr.data[myHead + j];
                assert(k < this->m_pdata->getN() + this->m_pdata->getNGhosts());

                bool issolid = checksolid(h_type_property_map.data, h_pos.data[k].w);

                Scalar3 pj;
                pj.x = h_pos.data[k].x;
                pj.y = h_pos.data[k].y;
                pj.z = h_pos.data[k].z;

                Scalar3 dx;
                dx.x = pi.x - pj.x;
                dx.y = pi.y - pj.y;
                dx.z = pi.z - pj.z;
                dx = box.minImage(dx);

                Scalar rsq = dot(dx, dx);
                if (this->m_const_slength && rsq > this->m_rcutsq) continue;

                Scalar mj   = h_velocity.data[k].w;
                Scalar Vj   = mj / rho0_ref;

                Scalar r     = sqrt(rsq);
                Scalar meanh = this->m_const_slength ? this->m_ch
                                                     : Scalar(0.5)*(hi + h_h.data[k]);
                Scalar epssqr = Scalar(0.01)*meanh*meanh;

                Scalar dwdr   = this->m_skernel->dwijdr(meanh, r);
                Scalar dwdr_r = dwdr / (r + epssqr);

                if (!issolid)
                    {
                    // Fluid neighbour: contribute to Shepard sum and its gradient.
                    Scalar wij = this->m_skernel->wij(meanh, r);
                    lambda        += Vj * wij;
                    // $\nabla_i W_{ij} = (\partial W/\partial r)/r \cdot \mathbf{r}_{ij}$ (points from $j$ toward $i$)
                    grad_lambda.x += Vj * dwdr_r * dx.x;
                    grad_lambda.y += Vj * dwdr_r * dx.y;
                    grad_lambda.z += Vj * dwdr_r * dx.z;
                    }
                else
                    {
                    // Solid neighbour: accumulate wall normal for contact-angle BC.
                    if (apply_ca)
                        {
                        n_wall_acc.x += Vj * dwdr_r * dx.x;
                        n_wall_acc.y += Vj * dwdr_r * dx.y;
                        n_wall_acc.z += Vj * dwdr_r * dx.z;
                        has_solid_nbr = true;
                        }
                    }
                }

            // Store lambda (already includes self contribution V_i * W_0).
            h_aux4.data[i].x = lambda;
            h_aux4.data[i].y = Scalar(0);  // kappa filled later by compute_curvature

            // Determine outward free-surface normal.
            Scalar gnorm_sq = dot(grad_lambda, grad_lambda);

            if (lambda < m_fs_threshold && gnorm_sq > eps_sq)
                {
                Scalar gnorm = sqrt(gnorm_sq);
                // Outward normal: $-\nabla\lambda/|\nabla\lambda|$ ($\nabla\lambda$ points inward toward bulk)
                Scalar3 n_fs;
                n_fs.x = -grad_lambda.x / gnorm;
                n_fs.y = -grad_lambda.y / gnorm;
                n_fs.z = -grad_lambda.z / gnorm;

                // Contact-angle correction at the triple line.
                if (has_solid_nbr)
                    {
                    Scalar nw_sq = dot(n_wall_acc, n_wall_acc);
                    if (nw_sq > eps_sq)
                        {
                        Scalar nw_norm = sqrt(nw_sq);
                        Scalar3 n_w;
                        n_w.x = n_wall_acc.x / nw_norm;
                        n_w.y = n_wall_acc.y / nw_norm;
                        n_w.z = n_wall_acc.z / nw_norm;

                        // Tangential component of n_fs in the wall plane.
                        Scalar n_dot_nw = dot(n_fs, n_w);
                        Scalar3 n_t;
                        n_t.x = n_fs.x - n_dot_nw * n_w.x;
                        n_t.y = n_fs.y - n_dot_nw * n_w.y;
                        n_t.z = n_fs.z - n_dot_nw * n_w.z;

                        Scalar nt_sq = dot(n_t, n_t);
                        if (nt_sq > eps_sq)
                            {
                            Scalar nt_norm = sqrt(nt_sq);
                            n_t.x /= nt_norm;
                            n_t.y /= nt_norm;
                            n_t.z /= nt_norm;

                            // $\hat{n}_\mathrm{corrected} = \sin\theta\,\hat{t}_w + \cos\theta\,\hat{n}_w$
                            n_fs.x = sin_ca * n_t.x + cos_ca * n_w.x;
                            n_fs.y = sin_ca * n_t.y + cos_ca * n_w.y;
                            n_fs.z = sin_ca * n_t.z + cos_ca * n_w.z;

                            // Renormalise after blending.
                            Scalar nfs_norm = sqrt(dot(n_fs, n_fs));
                            if (nfs_norm > eps_sq)
                                {
                                n_fs.x /= nfs_norm;
                                n_fs.y /= nfs_norm;
                                n_fs.z /= nfs_norm;
                                }
                            }
                        }
                    }

                h_fs_n.data[i].x = n_fs.x;
                h_fs_n.data[i].y = n_fs.y;
                h_fs_n.data[i].z = n_fs.z;
                }
            else
                {
                // Bulk particle: zero normal (acts as sentinel in later loops).
                h_fs_n.data[i].x = Scalar(0);
                h_fs_n.data[i].y = Scalar(0);
                h_fs_n.data[i].z = Scalar(0);
                }
            }
        } // End GPU Array Scope
    }


// ─────────────────────────────────────────────────────────────────────────────
// Step 2: compute_curvature
// ─────────────────────────────────────────────────────────────────────────────

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowFS<KT_, SET_>::compute_curvature(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7)
        << "Computing SinglePhaseFlowFS::compute_curvature" << std::endl;

    const BoxDim& box = this->m_pdata->getGlobalBox();
    const unsigned int group_size = this->m_fluidgroup->getNumMembers();

        { // GPU Array Scope
        ArrayHandle<Scalar4> h_pos     (this->m_pdata->getPositions(),    access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(),   access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_density (this->m_pdata->getDensities(),    access_location::host, access_mode::read);
        // aux2 holds the free-surface normals written by detect_freesurface().
        // In MPI runs they have been synced to ghost slots by update_ghost_aux24().
        ArrayHandle<Scalar3> h_fs_n    (this->m_pdata->getAuxiliaries2(), access_location::host, access_mode::read);
        ArrayHandle<Scalar3> h_aux4    (this->m_pdata->getAuxiliaries4(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar>  h_h       (this->m_pdata->getSlengths(),     access_location::host, access_mode::read);

        ArrayHandle<unsigned int> h_n_neigh  (this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist_arr(this->m_nlist->getNListArray(),  access_location::host, access_mode::read);
        ArrayHandle<size_t>       h_head_list(this->m_nlist->getHeadList(),    access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

        const Scalar eps_sq = Scalar(1e-12);

        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

            // Skip bulk particles.
            Scalar lambda_i = h_aux4.data[i].x;
            if (lambda_i >= m_fs_threshold)
                {
                h_aux4.data[i].y = Scalar(0);
                continue;
                }

            Scalar3 n_i = h_fs_n.data[i];
            if (dot(n_i, n_i) < eps_sq)
                {
                h_aux4.data[i].y = Scalar(0);
                continue;
                }

            Scalar3 pi;
            pi.x = h_pos.data[i].x;
            pi.y = h_pos.data[i].y;
            pi.z = h_pos.data[i].z;

            Scalar mi   = h_velocity.data[i].w;
            Scalar rhoi = h_density.data[i];
            Scalar Vi   = mi / rhoi;
            Scalar hi   = this->m_const_slength ? this->m_ch : h_h.data[i];

            Scalar kappa = Scalar(0);

            size_t       myHead = h_head_list.data[i];
            unsigned int size   = (unsigned int)h_n_neigh.data[i];

            for (unsigned int j = 0; j < size; j++)
                {
                unsigned int k = h_nlist_arr.data[myHead + j];

                // Curvature is a fluid-fluid quantity; skip solid particles.
                if (checksolid(h_type_property_map.data, h_pos.data[k].w)) continue;

                Scalar3 pj;
                pj.x = h_pos.data[k].x;
                pj.y = h_pos.data[k].y;
                pj.z = h_pos.data[k].z;

                Scalar3 dx;
                dx.x = pi.x - pj.x;
                dx.y = pi.y - pj.y;
                dx.z = pi.z - pj.z;
                dx = box.minImage(dx);

                Scalar rsq = dot(dx, dx);
                if (this->m_const_slength && rsq > this->m_rcutsq) continue;

                Scalar mj   = h_velocity.data[k].w;
                Scalar rhoj = h_density.data[k];
                Scalar Vj   = mj / rhoj;

                // Normal of neighbour $k$; zero for bulk neighbours, which is
                // correct: their contribution becomes $-V_j \hat{n}_i \cdot \nabla W_{ij}$.
                Scalar3 n_j = h_fs_n.data[k];

                Scalar r     = sqrt(rsq);
                Scalar meanh = this->m_const_slength ? this->m_ch
                                                     : Scalar(0.5)*(hi + h_h.data[k]);
                Scalar epssqr = Scalar(0.01)*meanh*meanh;
                Scalar dwdr   = this->m_skernel->dwijdr(meanh, r);
                Scalar dwdr_r = dwdr / (r + epssqr);

                // $\kappa_i \mathrel{+}= (1/V_i) \sum_j V_j (\hat{n}_j - \hat{n}_i) \cdot \nabla W_{ij}$
                // $\nabla_i W_{ij} = (\partial W/\partial r)/r \cdot \mathbf{r}_{ij}$ (pointing from $j$ toward $i$)
                Scalar3 dn;
                dn.x = n_j.x - n_i.x;
                dn.y = n_j.y - n_i.y;
                dn.z = n_j.z - n_i.z;

                kappa += Vj * dot(dn, dwdr_r * dx);
                }

            // Divide by V_i to get the divergence.
            h_aux4.data[i].y = kappa / Vi;
            }
        } // End GPU Array Scope
    }


// ─────────────────────────────────────────────────────────────────────────────
// Step 3: apply_freesurface_pressure
// ─────────────────────────────────────────────────────────────────────────────

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowFS<KT_, SET_>::apply_freesurface_pressure(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7)
        << "Computing SinglePhaseFlowFS::apply_freesurface_pressure" << std::endl;

    const unsigned int group_size = this->m_fluidgroup->getNumMembers();

        { // GPU Array Scope
        ArrayHandle<Scalar>  h_pressure(this->m_pdata->getPressures(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar3> h_aux4    (this->m_pdata->getAuxiliaries4(), access_location::host, access_mode::read);

        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

            if (h_aux4.data[i].x < m_fs_threshold)
                {
                // Clamp tensile pressure to zero at the free surface.
                if (h_pressure.data[i] < Scalar(0))
                    h_pressure.data[i] = Scalar(0);
                }
            }
        } // End GPU Array Scope
    }


// ─────────────────────────────────────────────────────────────────────────────
// Step 4: forcecomputation — TV forces + CSF surface tension
// ─────────────────────────────────────────────────────────────────────────────

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowFS<KT_, SET_>::forcecomputation(uint64_t timestep)
    {
    if (this->m_density_method == DENSITYSUMMATION)
        this->m_exec_conf->msg->notice(7)
            << "Computing SinglePhaseFlowFS::forcecomputation (SUMMATION)" << std::endl;
    else
        this->m_exec_conf->msg->notice(7)
            << "Computing SinglePhaseFlowFS::forcecomputation (CONTINUITY)" << std::endl;

    const BoxDim& box = this->m_pdata->getGlobalBox();
    const unsigned int group_size = this->m_fluidgroup->getNumMembers();

        { // GPU Array Scope
        ArrayHandle<Scalar4> h_force  (this->m_force,   access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar4> h_ratedpe(this->m_ratedpe, access_location::host, access_mode::readwrite);

        assert(h_force.data);
        assert(h_ratedpe.data);

        memset((void*)h_force.data,   0, sizeof(Scalar4)*this->m_force.getNumElements());
        memset((void*)h_ratedpe.data, 0, sizeof(Scalar4)*this->m_ratedpe.getNumElements());

        ArrayHandle<Scalar4> h_pos     (this->m_pdata->getPositions(),    access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(),   access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_density (this->m_pdata->getDensities(),    access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar>  h_pressure(this->m_pdata->getPressures(),    access_location::host, access_mode::readwrite);
        // aux1: fictitious solid velocity (Adami 2012)
        ArrayHandle<Scalar3> h_vf      (this->m_pdata->getAuxiliaries1(), access_location::host, access_mode::read);
        // aux2: free-surface normals (from detect_freesurface).
        //       Saved per-particle before being overwritten with BPC.
        ArrayHandle<Scalar3> h_bpc     (this->m_pdata->getAuxiliaries2(), access_location::host, access_mode::readwrite);
        // aux3: transport velocity
        ArrayHandle<Scalar3> h_tv      (this->m_pdata->getAuxiliaries3(), access_location::host, access_mode::read);
        // aux4: x=lambda, y=curvature
        ArrayHandle<Scalar3> h_aux4    (this->m_pdata->getAuxiliaries4(), access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_h       (this->m_pdata->getSlengths(),     access_location::host, access_mode::read);

        ArrayHandle<unsigned int> h_n_neigh  (this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist_arr(this->m_nlist->getNListArray(),  access_location::host, access_mode::read);
        ArrayHandle<size_t>       h_head_list(this->m_nlist->getHeadList(),    access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

        assert(h_pos.data != NULL);

        Scalar  temp0  = Scalar(0);
        double  max_vel = 0.0;

        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

            // ── Save free-surface normal from aux2 BEFORE zeroing for BPC ────
            // detect_freesurface() stored the outward unit normal in aux2.
            // We read it here, then overwrite aux2 with the BPC accumulation.
            Scalar3 n_fs_i;
            n_fs_i.x = h_bpc.data[i].x;
            n_fs_i.y = h_bpc.data[i].y;
            n_fs_i.z = h_bpc.data[i].z;

            Scalar lambda_i = h_aux4.data[i].x;
            Scalar kappa_i  = h_aux4.data[i].y;
            bool   is_surface = (lambda_i < m_fs_threshold)
                                && (dot(n_fs_i, n_fs_i) > Scalar(1e-12));

            // Zero BPC for this particle (overwrites the fs normal in aux2).
            h_bpc.data[i].x = Scalar(0);
            h_bpc.data[i].y = Scalar(0);
            h_bpc.data[i].z = Scalar(0);

            // ── Read particle i properties ────────────────────────────────────
            Scalar3 pi;
            pi.x = h_pos.data[i].x;
            pi.y = h_pos.data[i].y;
            pi.z = h_pos.data[i].z;

            Scalar3 vi;
            vi.x = h_velocity.data[i].x;
            vi.y = h_velocity.data[i].y;
            vi.z = h_velocity.data[i].z;
            Scalar mi = h_velocity.data[i].w;

            Scalar3 tvi;
            tvi.x = h_tv.data[i].x;
            tvi.y = h_tv.data[i].y;
            tvi.z = h_tv.data[i].z;

            Scalar Pi   = h_pressure.data[i];
            Scalar rhoi = h_density.data[i];
            Scalar Vi   = mi / rhoi;

            Scalar vi_total = sqrt(vi.x*vi.x + vi.y*vi.y + vi.z*vi.z);
            if (i == 0) { max_vel = vi_total; }
            else if (vi_total > max_vel) { max_vel = vi_total; }

            // Artificial stress tensor A_i for transport velocity correction.
            Scalar A11i = rhoi * vi.x * (tvi.x - vi.x);
            Scalar A12i = rhoi * vi.x * (tvi.y - vi.y);
            Scalar A13i = rhoi * vi.x * (tvi.z - vi.z);
            Scalar A21i = rhoi * vi.y * (tvi.x - vi.x);
            Scalar A22i = rhoi * vi.y * (tvi.y - vi.y);
            Scalar A23i = rhoi * vi.y * (tvi.z - vi.z);
            Scalar A31i = rhoi * vi.z * (tvi.x - vi.x);
            Scalar A32i = rhoi * vi.z * (tvi.y - vi.y);
            Scalar A33i = rhoi * vi.z * (tvi.z - vi.z);

            // ── Neighbour loop ────────────────────────────────────────────────
            size_t       myHead = h_head_list.data[i];
            unsigned int size   = (unsigned int)h_n_neigh.data[i];

            for (unsigned int j = 0; j < size; j++)
                {
                unsigned int k = h_nlist_arr.data[myHead + j];
                assert(k < this->m_pdata->getN() + this->m_pdata->getNGhosts());

                bool issolid = checksolid(h_type_property_map.data, h_pos.data[k].w);

                Scalar3 pj;
                pj.x = h_pos.data[k].x;
                pj.y = h_pos.data[k].y;
                pj.z = h_pos.data[k].z;

                Scalar3 dx;
                dx.x = pi.x - pj.x;
                dx.y = pi.y - pj.y;
                dx.z = pi.z - pj.z;
                dx = box.minImage(dx);

                Scalar rsq = dot(dx, dx);
                if (this->m_const_slength && rsq > this->m_rcutsq) continue;

                // Velocity: fictitious (Adami 2012) for solid, physical for fluid.
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

                Scalar3 tvj;
                tvj.x = h_tv.data[k].x;
                tvj.y = h_tv.data[k].y;
                tvj.z = h_tv.data[k].z;

                // Artificial stress tensor A_j.
                Scalar A11j = rhoj * vj.x * (tvj.x - vj.x);
                Scalar A12j = rhoj * vj.x * (tvj.y - vj.y);
                Scalar A13j = rhoj * vj.x * (tvj.z - vj.z);
                Scalar A21j = rhoj * vj.y * (tvj.x - vj.x);
                Scalar A22j = rhoj * vj.y * (tvj.y - vj.y);
                Scalar A23j = rhoj * vj.y * (tvj.z - vj.z);
                Scalar A31j = rhoj * vj.z * (tvj.x - vj.x);
                Scalar A32j = rhoj * vj.z * (tvj.y - vj.y);
                Scalar A33j = rhoj * vj.z * (tvj.z - vj.z);

                Scalar Pj = h_pressure.data[k];

                Scalar3 dv;
                dv.x = vi.x - vj.x;
                dv.y = vi.y - vj.y;
                dv.z = vi.z - vj.z;

                Scalar r     = sqrt(rsq);
                Scalar meanh = this->m_const_slength ? this->m_ch
                                                     : Scalar(0.5)*(h_h.data[i] + h_h.data[k]);
                Scalar epssqr = Scalar(0.01)*meanh*meanh;
                Scalar dwdr   = this->m_skernel->dwijdr(meanh, r);
                Scalar dwdr_r = dwdr / (r + epssqr);

                // Pressure interaction (Adami 2013 TV): $\bar{p}_{ij} = (\rho_j p_i + \rho_i p_j)/(\rho_i+\rho_j)$
                temp0 = (rhoj*Pi + rhoi*Pj) / (rhoi + rhoj);

                // Artificial viscosity (Monaghan 1983): $\Pi_{ij} = (-\alpha c \mu_{ij} + \beta \mu_{ij}^2)/\bar{\rho}_{ij}$
                Scalar avc = Scalar(0);
                if (this->m_artificial_viscosity && !issolid)
                    {
                    Scalar dotdvdx = dot(dv, dx);
                    if (dotdvdx < Scalar(0))
                        {
                        Scalar muij    = meanh * dotdvdx / (rsq + epssqr);
                        Scalar meanrho = Scalar(0.5)*(rhoi + rhoj);
                        avc = (-this->m_avalpha*this->m_c*muij
                               + this->m_avbeta*muij*muij) / meanrho;
                        }
                    }

                Scalar vijsqr = Vi*Vi + Vj*Vj;

                // Pressure force.
                h_force.data[i].x -= vijsqr*(temp0 + avc)*dwdr_r*dx.x;
                h_force.data[i].y -= vijsqr*(temp0 + avc)*dwdr_r*dx.y;
                h_force.data[i].z -= vijsqr*(temp0 + avc)*dwdr_r*dx.z;

                // Viscous force.
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

                // Artificial stress correction (TV).
                temp0 = Scalar(0.5)*vijsqr*dwdr_r;
                Scalar A1ij = (A11i+A11j)*dx.x + (A12i+A12j)*dx.y + (A13i+A13j)*dx.z;
                Scalar A2ij = (A21i+A21j)*dx.x + (A22i+A22j)*dx.y + (A23i+A23j)*dx.z;
                Scalar A3ij = (A31i+A31j)*dx.x + (A32i+A32j)*dx.y + (A33i+A33j)*dx.z;
                h_force.data[i].x += temp0 * A1ij;
                h_force.data[i].y += temp0 * A2ij;
                h_force.data[i].z += temp0 * A3ij;

                // Background pressure contribution to transport velocity (aux2).
                h_bpc.data[i].x -= vijsqr * this->m_eos->getTransportVelocityPressure()/mi * dwdr_r * dx.x;
                h_bpc.data[i].y -= vijsqr * this->m_eos->getTransportVelocityPressure()/mi * dwdr_r * dx.y;
                h_bpc.data[i].z -= vijsqr * this->m_eos->getTransportVelocityPressure()/mi * dwdr_r * dx.z;

                // Density continuity rate (DENSITYCONTINUITY only).
                if (this->m_density_method == DENSITYCONTINUITY)
                    {
                    if (issolid)
                        {
                        // Use physical velocity for the continuity equation.
                        vj.x = h_velocity.data[k].x;
                        vj.y = h_velocity.data[k].y;
                        vj.z = h_velocity.data[k].z;
                        dv.x = vi.x - vj.x;
                        dv.y = vi.y - vj.y;
                        dv.z = vi.z - vj.z;
                        }
                    h_ratedpe.data[i].x += rhoi * Vj * dot(dv, dwdr_r*dx);

                    if (!issolid && this->m_density_diffusion)
                        h_ratedpe.data[i].x -=
                            (Scalar(2)*this->m_ddiff*meanh*this->m_c*mj
                             * (rhoi/rhoj - Scalar(1)) * dot(dx, dwdr_r*dx))
                            / (rsq + epssqr);
                    }

                } // End neighbour loop

            // dp/dt chain rule (DENSITYCONTINUITY).
            if (this->m_density_method == DENSITYCONTINUITY)
                h_ratedpe.data[i].y =
                    this->m_eos->dPressuredDensity(rhoi) * h_ratedpe.data[i].x;

            // ── CSF surface tension force ─────────────────────────────────────
            // F_{σ,i} = −σ · κ_i · n̂_i · (m_i / ρ_i)
            // Sign: for a convex surface (κ > 0, n̂ outward) the force points
            // inward, compressing the surface — consistent with CSF conventions.
            if (is_surface && m_sigma > Scalar(0))
                {
                Scalar surf_coeff = -m_sigma * kappa_i * (mi / rhoi);
                h_force.data[i].x += surf_coeff * n_fs_i.x;
                h_force.data[i].y += surf_coeff * n_fs_i.y;
                h_force.data[i].z += surf_coeff * n_fs_i.z;
                }

            } // End fluid-particle outer loop

        this->m_timestep_list[5] = max_vel;

        } // End GPU Array Scope

    // Apply gravitational / external body force to fluid particles.
    this->applyBodyForce(timestep, this->m_fluidgroup);
    }


// ─────────────────────────────────────────────────────────────────────────────
// computeForces — full pipeline with free-surface extensions
// ─────────────────────────────────────────────────────────────────────────────

template<SmoothingKernelType KT_, StateEquationType SET_>
void SinglePhaseFlowFS<KT_, SET_>::computeForces(uint64_t timestep)
    {
    // Update the neighbour list.
    this->m_nlist->compute(timestep);

    if (!this->m_params_set)
        {
        this->m_exec_conf->msg->error()
            << "sph.models.SinglePhaseFlowFS requires parameters to be set before run()"
            << std::endl;
        throw std::runtime_error("Error computing SinglePhaseFlowFS forces");
        }

    // Remove solid particles from the nlist on the first call (timestep 0).
    if (!this->m_solid_removed)
        {
        this->m_nlist->forceUpdate();
        this->m_nlist->compute(timestep);
        this->mark_solid_particles_toremove(timestep);
        this->m_solid_removed = true;
        }

    // Shepard density renormalization (if enabled).
    if (this->m_shepard_renormalization && timestep % this->m_shepardfreq == 0)
        {
        this->renormalize_density(timestep);
        this->compute_pressure(timestep);
        this->m_pressure_initialized = true;
#ifdef ENABLE_MPI
        this->update_ghost_density(timestep);
#endif
        }

    // Periodic density reinitialization from summation (DENSITYCONTINUITY only).
    if (this->m_density_reinitialization && timestep % this->m_densityreinitfreq == 0)
        {
        if (this->m_density_method == DENSITYSUMMATION)
            this->m_exec_conf->msg->error()
                << "SinglePhaseFlowFS: density reinitialization only available with CONTINUITY" << std::endl;
        this->compute_ndensity(timestep);
        this->compute_pressure(timestep);
        this->m_pressure_initialized = true;
        }

    // Density + pressure update.
    if (this->m_density_method == DENSITYSUMMATION)
        {
        this->compute_ndensity(timestep);
        this->compute_pressure(timestep);
        }
    else // DENSITYCONTINUITY
        {
        if (!this->m_pressure_initialized)
            {
            this->compute_pressure(timestep);
            this->m_pressure_initialized = true;
            }
        }

#ifdef ENABLE_MPI
    this->update_ghost_density_pressure(timestep);
#endif

    // Compute fictitious solid velocities (Adami 2012 no-slip BC).
    this->compute_noslip(timestep);

#ifdef ENABLE_MPI
    this->update_ghost_aux1(timestep);
#endif

    // ── Free-surface pipeline ─────────────────────────────────────────────────
    detect_freesurface(timestep);

#ifdef ENABLE_MPI
    // Sync fs normals (aux2) and lambda/kappa (aux4) to ghost particles so
    // compute_curvature() can read neighbour normals across MPI boundaries.
    update_ghost_aux24(timestep);
#endif

    compute_curvature(timestep);
    apply_freesurface_pressure(timestep);
    // ── End free-surface pipeline ─────────────────────────────────────────────

    // TV force loop with CSF surface tension.
    forcecomputation(timestep);

#ifdef ENABLE_MPI
    // Sync aux1 (solid BC), aux2 (BPC for KickDriftKickTV), aux3 (TV) to ghosts.
    this->SinglePhaseFlowTV<KT_, SET_>::update_ghost_aux123(timestep);
#endif
    }


// ─────────────────────────────────────────────────────────────────────────────
// pybind11 export
// ─────────────────────────────────────────────────────────────────────────────

namespace detail
{

template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlowFS(pybind11::module& m, std::string name)
    {
    pybind11::class_<SinglePhaseFlowFS<KT_, SET_>,
                     SPHBaseClass<KT_, SET_>,
                     std::shared_ptr<SinglePhaseFlowFS<KT_, SET_>>>(m, name.c_str())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<SmoothingKernel<KT_>>,
                            std::shared_ptr<StateEquation<SET_>>,
                            std::shared_ptr<nsearch::NeighborList>,
                            std::shared_ptr<ParticleGroup>,
                            std::shared_ptr<ParticleGroup>,
                            DensityMethod,
                            ViscosityMethod>())
        .def("setParams",                      &SinglePhaseFlowFS<KT_, SET_>::setParams)
        .def("setFSParams",                    &SinglePhaseFlowFS<KT_, SET_>::setFSParams)
        .def("getSigma",                       &SinglePhaseFlowFS<KT_, SET_>::getSigma)
        .def("getFSThreshold",                 &SinglePhaseFlowFS<KT_, SET_>::getFSThreshold)
        .def("getContactAngle",                &SinglePhaseFlowFS<KT_, SET_>::getContactAngle)
        .def("getDensityMethod",               &SinglePhaseFlowFS<KT_, SET_>::getDensityMethod)
        .def("setDensityMethod",               &SinglePhaseFlowFS<KT_, SET_>::setDensityMethod)
        .def("getViscosityMethod",             &SinglePhaseFlowFS<KT_, SET_>::getViscosityMethod)
        .def("setViscosityMethod",             &SinglePhaseFlowFS<KT_, SET_>::setViscosityMethod)
        .def("setConstSmoothingLength",        &SinglePhaseFlowFS<KT_, SET_>::setConstSmoothingLength)
        .def("computeSolidForces",             &SinglePhaseFlowFS<KT_, SET_>::computeSolidForces)
        .def("activateArtificialViscosity",    &SinglePhaseFlowFS<KT_, SET_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity",  &SinglePhaseFlowFS<KT_, SET_>::deactivateArtificialViscosity)
        .def("activateDensityDiffusion",       &SinglePhaseFlowFS<KT_, SET_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion",     &SinglePhaseFlowFS<KT_, SET_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &SinglePhaseFlowFS<KT_, SET_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization",&SinglePhaseFlowFS<KT_, SET_>::deactivateShepardRenormalization)
        .def("activateDensityReinitialization",&SinglePhaseFlowFS<KT_, SET_>::activateDensityReinitialization)
        .def("deactivateDensityReinitialization",&SinglePhaseFlowFS<KT_, SET_>::deactivateDensityReinitialization)
        .def("activatePowerLaw",               &SinglePhaseFlowFS<KT_, SET_>::activatePowerLaw,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("mu_min") = Scalar(0))
        .def("activateCarreau",                &SinglePhaseFlowFS<KT_, SET_>::activateCarreau)
        .def("activateBingham",                &SinglePhaseFlowFS<KT_, SET_>::activateBingham,
             pybind11::arg("mu_p"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("activateHerschelBulkley",        &SinglePhaseFlowFS<KT_, SET_>::activateHerschelBulkley,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("deactivateNonNewtonian",         &SinglePhaseFlowFS<KT_, SET_>::deactivateNonNewtonian)
        .def("setAcceleration",                &SPHBaseClass<KT_, SET_>::setAcceleration)
        .def("setRCut",                        &SinglePhaseFlowFS<KT_, SET_>::setRCutPython)
        ;
    }

} // end namespace detail


// ─────────────────────────────────────────────────────────────────────────────
// Explicit template instantiations (5 kernels × 2 EOS = 10 variants)
// ─────────────────────────────────────────────────────────────────────────────

template class PYBIND11_EXPORT SinglePhaseFlowFS<wendlandc2, linear>;
template class PYBIND11_EXPORT SinglePhaseFlowFS<wendlandc2, tait>;
template class PYBIND11_EXPORT SinglePhaseFlowFS<wendlandc4, linear>;
template class PYBIND11_EXPORT SinglePhaseFlowFS<wendlandc4, tait>;
template class PYBIND11_EXPORT SinglePhaseFlowFS<wendlandc6, linear>;
template class PYBIND11_EXPORT SinglePhaseFlowFS<wendlandc6, tait>;
template class PYBIND11_EXPORT SinglePhaseFlowFS<quintic,    linear>;
template class PYBIND11_EXPORT SinglePhaseFlowFS<quintic,    tait>;
template class PYBIND11_EXPORT SinglePhaseFlowFS<cubicspline,linear>;
template class PYBIND11_EXPORT SinglePhaseFlowFS<cubicspline,tait>;


namespace detail
{
    template void export_SinglePhaseFlowFS<wendlandc2, linear>(pybind11::module& m, std::string name = "SinglePFFS_WC2_L");
    template void export_SinglePhaseFlowFS<wendlandc2, tait>  (pybind11::module& m, std::string name = "SinglePFFS_WC2_T");
    template void export_SinglePhaseFlowFS<wendlandc4, linear>(pybind11::module& m, std::string name = "SinglePFFS_WC4_L");
    template void export_SinglePhaseFlowFS<wendlandc4, tait>  (pybind11::module& m, std::string name = "SinglePFFS_WC4_T");
    template void export_SinglePhaseFlowFS<wendlandc6, linear>(pybind11::module& m, std::string name = "SinglePFFS_WC6_L");
    template void export_SinglePhaseFlowFS<wendlandc6, tait>  (pybind11::module& m, std::string name = "SinglePFFS_WC6_T");
    template void export_SinglePhaseFlowFS<quintic,    linear>(pybind11::module& m, std::string name = "SinglePFFS_Q_L");
    template void export_SinglePhaseFlowFS<quintic,    tait>  (pybind11::module& m, std::string name = "SinglePFFS_Q_T");
    template void export_SinglePhaseFlowFS<cubicspline,linear>(pybind11::module& m, std::string name = "SinglePFFS_CS_L");
    template void export_SinglePhaseFlowFS<cubicspline,tait>  (pybind11::module& m, std::string name = "SinglePFFS_CS_T");
} // end namespace detail

} // end namespace sph
} // end namespace hoomd

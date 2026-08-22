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

#include "TwoPhaseFlow.h"

#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

using namespace std;

namespace hoomd 
{
namespace sph
{
/*! Constructor
*/
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
TwoPhaseFlow<KT_, SET1_, SET2_>::TwoPhaseFlow(std::shared_ptr<SystemDefinition> sysdef,
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
    : SPHBaseClass<KT_, SET1_>(sysdef, skernel, equationofstate1, nlist), m_fluidgroup1(fluidgroup1), m_fluidgroup2(fluidgroup2), 
        m_solidgroup(solidgroup), m_eos1(equationofstate1), m_eos2(equationofstate2), m_typpair_idx(this->m_pdata->getNTypes())
    {
        this->m_exec_conf->msg->notice(5) << "Constructing TwoPhaseFlow" << std::endl;

        // Set private attributes to default values
        m_const_slength = false;
        m_params_set = false;
        m_compute_solid_forces = false;
        m_artificial_viscosity = false;
        m_riemann_dissipation  = false;
        m_riemann_beta         = Scalar(1.0);
        m_consistent_interface_pressure = false;
        m_density_diffusion = false;
        m_shepard_renormalization = false;
        m_fickian_shifting = false;
        m_pressure_initialized = false;
        m_density_reinitialization = false;
        m_densityreinitfreq = 1;
        m_particle_shifting       = false;
        m_shift_A                 = Scalar(0.2);
        m_shift_R                 = Scalar(0.2);
        m_shift_n                 = 4;
        m_shift_interface_condition = true;
        m_ch = Scalar(0.0);
        m_rcut = Scalar(0.0);
        m_rcutsq = Scalar(0.0);
        m_avalpha = Scalar(0.0);
        m_avbeta = Scalar(0.0);
        m_ddiff = Scalar(0.0);
        m_shepardfreq = 1;

        m_omega_adv = Scalar(180);
        m_omega_rec = Scalar(0);
        m_hysteresis = false;
        m_nn_model1  = NEWTONIAN;
        m_nn_K1      = Scalar(0.0);
        m_nn_n1      = Scalar(1.0);
        m_nn_mu0_1   = Scalar(0.0);
        m_nn_muinf_1 = Scalar(0.0);
        m_nn_lambda1 = Scalar(0.0);
        m_nn_tauy1   = Scalar(0.0);
        m_nn_m1      = Scalar(0.0);
        m_nn_mu_min1 = Scalar(0.0);

        m_nn_model2  = NEWTONIAN;
        m_nn_K2      = Scalar(0.0);
        m_nn_n2      = Scalar(1.0);
        m_nn_mu0_2   = Scalar(0.0);
        m_nn_muinf_2 = Scalar(0.0);
        m_nn_lambda2 = Scalar(0.0);
        m_nn_tauy2   = Scalar(0.0);
        m_nn_m2      = Scalar(0.0);
        m_nn_mu_min2 = Scalar(0.0);

        m_solid_removed = false;

        // Sanity checks
        assert(this->m_pdata);
        assert(this->m_nlist);
        assert(this->m_skernel);
        assert(this->m_eos1);
        assert(this->m_eos2);

        // If $c_1 \ne c_2$, back-pressures differ; apply $\max(b_1, b_2)$ to both phases
        Scalar bp1 = this->m_eos1->getBackgroundPressure();
        Scalar bp2 = this->m_eos2->getBackgroundPressure();
        if ( bp1 > bp2 )
            this->m_eos2->setBackPressure(bp1);
        if ( bp2 > bp1 )
            this->m_eos1->setBackPressure(bp2);

        // Create new fluid ParticleGroup by forming union of fluid 1 and 2
        this->m_fluidgroup = ParticleGroup::groupUnion(fluidgroup1, fluidgroup2);

        // Contruct type vectors
        this->constructTypeVectors(fluidgroup1,&m_fluidtypes1);
        this->constructTypeVectors(fluidgroup2,&m_fluidtypes2);
        this->constructTypeVectors(this->m_fluidgroup,&m_fluidtypes);
        this->constructTypeVectors(solidgroup,&m_solidtypes);

        // all particle groups are based on the same particle data
        unsigned int num_types = this->m_sysdef->getParticleData()->getNTypes();

        m_type_property_map = GPUArray<unsigned int>(num_types, this->m_exec_conf);
        {
            ArrayHandle<unsigned int> h_type_property_map(m_type_property_map, access_location::host, access_mode::overwrite);
            fill_n(h_type_property_map.data, num_types, SolidFluidTypeBit::NONE);
            // no need to parallelize this as there should only be a few particle types
            for (unsigned int i = 0; i < m_fluidtypes1.size(); i++) {
                h_type_property_map.data[m_fluidtypes1[i]] |= SolidFluidTypeBit::FLUID | SolidFluidTypeBit::FLUID1;
            }
            for (unsigned int i = 0; i < m_fluidtypes2.size(); i++) {
                h_type_property_map.data[m_fluidtypes2[i]] |= SolidFluidTypeBit::FLUID | SolidFluidTypeBit::FLUID2;
            }
            for (unsigned int i = 0; i < m_solidtypes.size(); i++) {
                h_type_property_map.data[m_solidtypes[i]] |= SolidFluidTypeBit::SOLID;
            }
        }

        // Set simulations methods
        m_density_method = mdensitymethod;
        m_viscosity_method = mviscositymethod;
        m_colorgradient_method = mcolorgradientmethod;

        // Get necessary variables from kernel and EOS classes
        m_rho01  = equationofstate1->getRestDensity();
        m_rho02  = equationofstate2->getRestDensity();
        m_c1     = equationofstate1->getSpeedOfSound();
        m_c2     = equationofstate2->getSpeedOfSound();
        m_cmax   = this->m_c1 > this->m_c2 ? this->m_c1 : this->m_c2;
        m_kappa = skernel->getKernelKappa();

        m_r_cut_nlist = std::make_shared<GPUArray<Scalar>>(m_typpair_idx.getNumElements(), this->m_exec_conf);
        this->m_nlist->addRCutMatrix(m_r_cut_nlist);

#ifdef ENABLE_MPI
    if (this->m_sysdef->isDomainDecomposed())
        {
        auto comm_weak = this->m_sysdef->getCommunicator();
        assert(comm_weak.lock());
        m_comm = comm_weak.lock();
        }
#endif

}

/*! Destructor
*/
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
TwoPhaseFlow<KT_, SET1_, SET2_>::~TwoPhaseFlow()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying TwoPhaseFlow" << std::endl;
    }


template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::activateShepardRenormalization(unsigned int shepardfreq)
    {
        if (shepardfreq <= 0)
            {
                this->m_exec_conf->msg->error() << "sph.models.TwoPhaseFlow: Shepard density reinitialization period has to be a positive real number" << std::endl;
                throw std::runtime_error("Error initializing TwoPhaseFlow.");
            }
        m_shepard_renormalization = true;
        m_shepardfreq = shepardfreq;
    }


/*! \post Set model parameters
 */

template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::setParams(Scalar mu1, Scalar mu2, Scalar sigma12, Scalar omega)
    {
    this->m_exec_conf->msg->notice(7) << "Setting TwoPhaseFlow parameters" << std::endl;

    this->m_mu1 = mu1;
    this->m_mu2 = mu2;
    if (this->m_mu1 <= 0 || this->m_mu2 <= 0)
         {
         this->m_exec_conf->msg->error() << "sph.models.TwoPhaseFlow: Dynamic viscosity has to be a positive real number" << std::endl;
         throw std::runtime_error("Error initializing TwoPhaseFlow.");
         }

    this->m_sigma12 = sigma12;
    if (this->m_sigma12 < 0)
         {
         this->m_exec_conf->msg->error() << "sph.models.TwoPhaseFlow: Fluid interfacial tension has to be a positive real number" << std::endl;
         throw std::runtime_error("Error initializing TwoPhaseFlow.");
         }

    this->m_omega = omega;
    if (this->m_omega <= 0 || this->m_omega > Scalar(180))
         {
         this->m_exec_conf->msg->error() << "sph.models.TwoPhaseFlow: Contact angle has to be between 0 and 180 degree" << std::endl;
         throw std::runtime_error("Error initializing TwoPhaseFlow.");
         }

    // Young's equation: $\sigma_{s1} - \sigma_{s2} = \sigma_{12} \cos\theta$
    if ( this->m_omega == Scalar(90) )
        {
        this->m_sigma01 = 0.0;
        this->m_sigma02 = 0.0;
        }
    else if ( this->m_omega < Scalar(90) )
        {
        this->m_sigma01 = this->m_sigma12 * cos( this->m_omega * ( M_PI / Scalar(180) ) );
        this->m_sigma02 = 0.0;
        }
    else if ( this->m_omega > Scalar(90) )
        {
        this->m_sigma01 = 0.0;
        this->m_sigma02 = this->m_sigma12 * cos( (Scalar(180)-m_omega) * ( M_PI / Scalar(180) ) );
        }

    this->m_params_set = true;
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::setHysteresis(Scalar omega_rec, Scalar omega_adv)
    {
    if (omega_rec < 0 || omega_rec > 180 || omega_adv < 0 || omega_adv > 180)
        throw std::runtime_error("Hysteresis angles must be in [0,180] deg.");
    if (omega_rec >= omega_adv)
        throw std::runtime_error("omega_rec must be < omega_adv.");
    m_omega_rec  = omega_rec;
    m_omega_adv  = omega_adv;
    m_hysteresis = true;
    }


template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::update_ghost_density_pressure(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Update Ghost density, pressure" << std::endl;

#ifdef ENABLE_MPI
    if (this->m_comm)
        {
        CommFlags flags(0);
        flags[comm_flag::tag] = 0;
        flags[comm_flag::position] = 0;
        flags[comm_flag::velocity] = 0;
        // flags[comm_flag::dpe] = 1;
        flags[comm_flag::density] = 1;
        flags[comm_flag::pressure] = 1;
        flags[comm_flag::energy] = 0;
        flags[comm_flag::auxiliary1] = 0;
        flags[comm_flag::auxiliary2] = 0;
        flags[comm_flag::auxiliary3] = 0;
        flags[comm_flag::auxiliary4] = 0;
        flags[comm_flag::body] = 0;
        flags[comm_flag::image] = 0;
        flags[comm_flag::net_force] = 0;
        flags[comm_flag::net_ratedpe] = 0;
        this->m_comm->setFlags(flags);
        this->m_comm->beginUpdateGhosts(timestep);
        this->m_comm->finishUpdateGhosts(timestep);
        }
#endif
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::update_ghost_density_pressure_energy(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Update Ghost density, pressure and energy" << std::endl;

#ifdef ENABLE_MPI
    if (this->m_comm)
        {
        CommFlags flags(0);
        flags[comm_flag::tag] = 0;
        flags[comm_flag::position] = 0;
        flags[comm_flag::velocity] = 0;
        flags[comm_flag::density] = 1;
        flags[comm_flag::pressure] = 1;
        flags[comm_flag::energy] = 1; // L2 norm of the color gradient
        flags[comm_flag::auxiliary1] = 0;
        flags[comm_flag::auxiliary2] = 0;
        flags[comm_flag::auxiliary3] = 0;
        flags[comm_flag::auxiliary4] = 0;
        flags[comm_flag::body] = 0;
        flags[comm_flag::image] = 0;
        flags[comm_flag::net_force] = 0;
        flags[comm_flag::net_ratedpe] = 0;
        this->m_comm->setFlags(flags);
        this->m_comm->beginUpdateGhosts(timestep);
        this->m_comm->finishUpdateGhosts(timestep);
        }
#endif
    }



/*! Compute the per-particle scalar shear rate gamma_dot = sqrt(2 D:D).

    Same construction as SinglePhaseFlow::compute_strain_rate: L-matrix
    renormalized velocity gradient over ALL fluid neighbors (the velocity field
    is continuous across the fluid-fluid interface, so both phases contribute),
    with solid neighbors entering through their fictitious (Adami) velocities.
    Result stored in the energy array (mutually exclusive with Fickian
    shifting, which stores |grad C|^2 there — enforced in computeForces).
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_strain_rate(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Strain rate" << std::endl;

    const BoxDim& box = this->m_pdata->getGlobalBox();
    const unsigned int group_size = this->m_fluidgroup->getNumMembers();

        { // GPU Array Scope
        ArrayHandle<Scalar>  h_energy(this->m_pdata->getEnergies(), access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
        ArrayHandle<Scalar3> h_vf(this->m_pdata->getAuxiliaries1(), access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_density(this->m_pdata->getDensities(), access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

        ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
        ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

        // Acquire the group index array once: getMemberIndex() acquires an
        // ArrayHandle per call, which is not thread-safe inside the parallel loop
        ArrayHandle<unsigned int> h_members_omp1(this->m_fluidgroup->getIndexArray(), access_location::host, access_mode::read);
        #pragma omp parallel for
        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            unsigned int i = h_members_omp1.data[group_idx];

            Scalar3 pi;
            pi.x = h_pos.data[i].x;
            pi.y = h_pos.data[i].y;
            pi.z = h_pos.data[i].z;
            Scalar3 vi;
            vi.x = h_velocity.data[i].x;
            vi.y = h_velocity.data[i].y;
            vi.z = h_velocity.data[i].z;

            Scalar A[9] = {0,0,0, 0,0,0, 0,0,0};
            Scalar bx[3] = {0,0,0};
            Scalar by[3] = {0,0,0};
            Scalar bz[3] = {0,0,0};

            size_t myHead = h_head_list.data[i];
            unsigned int size = (unsigned int)h_n_neigh.data[i];
            for (unsigned int j = 0; j < size; j++)
                {
                unsigned int k = h_nlist.data[myHead + j];

                Scalar mk = h_velocity.data[k].w;
                // Skip solid particles marked for removal (mass = -999)
                if (mk < Scalar(0))
                    continue;

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
                if ( this->m_const_slength && rsq > this->m_rcutsq )
                    continue;

                Scalar3 vj;
                if ( checksolid(h_type_property_map.data, h_pos.data[k].w) )
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

                Scalar r = sqrt(rsq);
                Scalar meanh = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
                Scalar dwdr   = this->m_skernel->dwijdr(meanh, r);
                Scalar dwdr_r = (r > Scalar(1e-8)*meanh) ? dwdr/r : Scalar(0);

                Scalar Vk = mk / h_density.data[k];

                Scalar3 gradW;
                gradW.x = dwdr_r * dx.x;
                gradW.y = dwdr_r * dx.y;
                gradW.z = dwdr_r * dx.z;

                Scalar c = -Vk;
                A[0] += c*gradW.x*dx.x; A[1] += c*gradW.x*dx.y; A[2] += c*gradW.x*dx.z;
                A[3] += c*gradW.y*dx.x; A[4] += c*gradW.y*dx.y; A[5] += c*gradW.y*dx.z;
                A[6] += c*gradW.z*dx.x; A[7] += c*gradW.z*dx.y; A[8] += c*gradW.z*dx.z;

                Scalar3 dv;
                dv.x = vj.x - vi.x;
                dv.y = vj.y - vi.y;
                dv.z = vj.z - vi.z;

                bx[0] += Vk*dv.x*gradW.x; bx[1] += Vk*dv.x*gradW.y; bx[2] += Vk*dv.x*gradW.z;
                by[0] += Vk*dv.y*gradW.x; by[1] += Vk*dv.y*gradW.y; by[2] += Vk*dv.y*gradW.z;
                bz[0] += Vk*dv.z*gradW.x; bz[1] += Vk*dv.z*gradW.y; bz[2] += Vk*dv.z*gradW.z;
                }

            Scalar det = A[0]*(A[4]*A[8]-A[5]*A[7])
                       - A[1]*(A[3]*A[8]-A[5]*A[6])
                       + A[2]*(A[3]*A[7]-A[4]*A[6]);

            Scalar M[9];
            if ( fabs(det) > Scalar(0.01) )
                {
                Scalar invdet = Scalar(1.0)/det;
                Scalar Ainv[9];
                Ainv[0] = invdet*(A[4]*A[8]-A[5]*A[7]);
                Ainv[1] = invdet*(A[2]*A[7]-A[1]*A[8]);
                Ainv[2] = invdet*(A[1]*A[5]-A[2]*A[4]);
                Ainv[3] = invdet*(A[5]*A[6]-A[3]*A[8]);
                Ainv[4] = invdet*(A[0]*A[8]-A[2]*A[6]);
                Ainv[5] = invdet*(A[2]*A[3]-A[0]*A[5]);
                Ainv[6] = invdet*(A[3]*A[7]-A[4]*A[6]);
                Ainv[7] = invdet*(A[1]*A[6]-A[0]*A[7]);
                Ainv[8] = invdet*(A[0]*A[4]-A[1]*A[3]);
                for (int a = 0; a < 3; a++)
                    {
                    const Scalar* b = (a == 0) ? bx : ((a == 1) ? by : bz);
                    M[3*a+0] = Ainv[0]*b[0] + Ainv[1]*b[1] + Ainv[2]*b[2];
                    M[3*a+1] = Ainv[3]*b[0] + Ainv[4]*b[1] + Ainv[5]*b[2];
                    M[3*a+2] = Ainv[6]*b[0] + Ainv[7]*b[1] + Ainv[8]*b[2];
                    }
                }
            else
                {
                M[0]=bx[0]; M[1]=bx[1]; M[2]=bx[2];
                M[3]=by[0]; M[4]=by[1]; M[5]=by[2];
                M[6]=bz[0]; M[7]=bz[1]; M[8]=bz[2];
                }

            Scalar Dxx = M[0];
            Scalar Dyy = M[4];
            Scalar Dzz = M[8];
            Scalar Dxy = Scalar(0.5)*(M[1]+M[3]);
            Scalar Dxz = Scalar(0.5)*(M[2]+M[6]);
            Scalar Dyz = Scalar(0.5)*(M[5]+M[7]);
            Scalar DD = Dxx*Dxx + Dyy*Dyy + Dzz*Dzz
                      + Scalar(2)*(Dxy*Dxy + Dxz*Dxz + Dyz*Dyz);

            h_energy.data[i] = sqrt(Scalar(2)*DD);
            } // End fluid particle loop
        } // End GPU Array Scope
    } // End compute_strain_rate

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::update_ghost_density(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Update ghost density" << std::endl;

#ifdef ENABLE_MPI
    if (this->m_comm)
        {
        CommFlags flags(0);
        flags[comm_flag::tag] = 0;
        flags[comm_flag::position] = 0;
        flags[comm_flag::velocity] = 0;
        flags[comm_flag::density] = 1;
        flags[comm_flag::pressure] = 0;
        flags[comm_flag::energy] = 0;
        flags[comm_flag::auxiliary1] = 0;
        flags[comm_flag::auxiliary2] = 0;
        flags[comm_flag::auxiliary3] = 0;
        flags[comm_flag::auxiliary4] = 0;
        flags[comm_flag::body] = 0;
        flags[comm_flag::image] = 0;
        flags[comm_flag::net_force] = 0;
        flags[comm_flag::net_ratedpe] = 0;
        this->m_comm->setFlags(flags);
        this->m_comm->beginUpdateGhosts(timestep);
        this->m_comm->finishUpdateGhosts(timestep);
        }
#endif
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::update_ghost_aux123(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Update Ghost density, pressure, aux1-3" << std::endl;

#ifdef ENABLE_MPI
    if (this->m_comm)
        {
        CommFlags flags(0);
        flags[comm_flag::tag] = 0;
        flags[comm_flag::position] = 0;
        flags[comm_flag::velocity] = 0;
        flags[comm_flag::density] = 1;
        flags[comm_flag::pressure] = 1;
        flags[comm_flag::energy] = 0;
        flags[comm_flag::auxiliary1] = 1; // ficticious velocity 
        flags[comm_flag::auxiliary2] = 1;
        flags[comm_flag::auxiliary3] = 1;
        flags[comm_flag::auxiliary4] = 0;
        flags[comm_flag::body] = 0;
        flags[comm_flag::image] = 0;
        flags[comm_flag::net_force] = 0;
        flags[comm_flag::net_ratedpe] = 0;
        this->m_comm->setFlags(flags);
        this->m_comm->beginUpdateGhosts(timestep);
        this->m_comm->finishUpdateGhosts(timestep);
        }
#endif
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::update_ghost_aux4(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Update Ghost aux4" << std::endl;

#ifdef ENABLE_MPI
    if (this->m_comm)
        {
        CommFlags flags(0);
        flags[comm_flag::tag] = 0;
        flags[comm_flag::position] = 0;
        flags[comm_flag::velocity] = 0;
        flags[comm_flag::density] = 0;
        flags[comm_flag::pressure] = 0;
        flags[comm_flag::energy] = 0;
        flags[comm_flag::auxiliary1] = 0;
        flags[comm_flag::auxiliary2] = 0;
        flags[comm_flag::auxiliary3] = 0;
        flags[comm_flag::auxiliary4] = 1;
        flags[comm_flag::body] = 0;
        flags[comm_flag::image] = 0;
        flags[comm_flag::net_force] = 0;
        flags[comm_flag::net_ratedpe] = 0;
        this->m_comm->setFlags(flags);
        this->m_comm->beginUpdateGhosts(timestep);
        this->m_comm->finishUpdateGhosts(timestep);
        }
#endif
    }


template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::validateTypes(unsigned int typ1,
                                             unsigned int typ2,
                                             std::string action)
    {
    auto n_types = this->m_pdata->getNTypes();
    if (typ1 >= n_types || typ2 >= n_types)
        {
        throw std::runtime_error("Error in" + action + " for pair potential. Invalid type");
        }
    }


/*! \param typ1 First type index in the pair
    \param typ2 Second type index in the pair
    \param rcut Cutoff radius to set
    \note When setting the value for (\a typ1, \a typ2), the parameter for (\a typ2, \a typ1) is
   automatically set.
*/
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::setRcut(unsigned int typ1, unsigned int typ2, Scalar rcut)
    {
    validateTypes(typ1, typ2, "setting r_cut");
        {
        // store r_cut unmodified for so the neighbor list knows what particles to include
        ArrayHandle<Scalar> h_r_cut_nlist(*m_r_cut_nlist,
                                          access_location::host,
                                          access_mode::readwrite);
        h_r_cut_nlist.data[m_typpair_idx(typ1, typ2)] = rcut;
        h_r_cut_nlist.data[m_typpair_idx(typ2, typ1)] = rcut;
        }

    // notify the neighbor list that we have changed r_cut values
    this->m_nlist->notifyRCutMatrixChange();
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::setRCutPython(pybind11::tuple types, Scalar r_cut)
    {
    auto typ1 = this->m_pdata->getTypeByName(types[0].cast<std::string>());
    auto typ2 = this->m_pdata->getTypeByName(types[1].cast<std::string>());
    setRcut(typ1, typ2, r_cut);
    }


/*! Mark solid particles to remove
    set mass of a particle to -999.0
 */

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::mark_solid_particles_toremove(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Mark solid Particles to remove at timestep " << timestep << std::endl;

    const unsigned int group_size = this->m_solidgroup->getNumMembers();
    unsigned int size;
    size_t myHead;
    { // GPU Array Scope
    // Grab handles for particle and neighbor data
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::readwrite);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // For all solid particles
    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp2(this->m_solidgroup->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for private(size, myHead)
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = h_members_omp2.data[group_idx];

        // check if solid particle has any fluid neighbor
        bool solid_w_fluid_neigh = false;
        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {
            unsigned int k = h_nlist.data[myHead + j];
            if ( checkfluid(h_type_property_map.data, h_pos.data[k].w) )
                {
                solid_w_fluid_neigh = true;
                break;
                }
            }
        if ( !(solid_w_fluid_neigh) )
            {
            // Solid particles which do not have fluid neighbors are marked
            // using mass=-999 so that they can be deleted during simulation
            h_velocity.data[i].w = Scalar(-999.0);
            }

        } // End solid particle loop
    } // End GPU Array Scope

    } // End mark solid particles to remove


/*! Perform particle concentration gradient
 * This method computes and stores the particle concentration gradient
 * in h_energy[i] to be reused in the computation of the Surface Force.
 * We overwrite h_pressure druing that, which is uncritical, since it is computed afterwards 
 * in compute_pressure, purely on the density of the particle 
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_particle_concentration_gradient(uint64_t timestep)
{
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Number Density" << std::endl;

    // Grab handles for particle data
    ArrayHandle<Scalar> h_density(this->m_pdata->getDensities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar> h_energy(this->m_pdata->getEnergies(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_pressure(this->m_pdata->getPressures(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Zero data before force calculation
    memset((void*)h_pressure.data,0,sizeof(Scalar)*this->m_pdata->getPressures().getNumElements());

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    unsigned int size;
    size_t myHead;

    // Precompute self-density for homogeneous smoothing lengths
    // Scalar w0 = this->m_skernel->w0(m_ch);

    // Particle loop to compute the particle concentration
    // For each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp3(this->m_fluidgroup->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for private(size, myHead)
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int i = h_members_omp3.data[group_idx];
        
        // set temp variable to zero 

        // Access the particle's position
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;

        // Scalar Ci;
        // Ci = w0;

        // Loop over all of the neighbors of this particle
        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];

        for (unsigned int j = 0; j < size; j++)
        {
            // Index of neighbor
            unsigned int k = h_nlist.data[myHead + j];

            // Access neighbor position
            Scalar3 pj;
            pj.x = h_pos.data[k].x;
            pj.y = h_pos.data[k].y;
            pj.z = h_pos.data[k].z;

            // Compute distance vector
            // Scalar3 dx = pj - pi;
            Scalar3 dx;
            dx.x = pi.x - pj.x;
            dx.y = pi.y - pj.y;
            dx.z = pi.z - pj.z;

            Scalar mj   = h_velocity.data[k].w;
            Scalar rhoj = h_density.data[k];

            // Apply periodic boundary conditions
            dx = box.minImage(dx);

            // Calculate squared distance
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, continue with next neighbor in loop
            if ( this->m_const_slength && rsq > m_rcutsq )
                continue;

            // Calculate distance
            Scalar r = sqrt(rsq);
            // $\sum_j (m_j/\rho_j) W_{ij}$ stored in h_pressure
            h_pressure.data[i] += (mj/rhoj)*(this->m_const_slength ? this->m_skernel->wij(m_ch,r) : this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r));

        } // End neighbour loop

    } // End fluid group loop

    Scalar3 gradCi = make_scalar3(0, 0, 0);
    Scalar  temp0 = Scalar(0);
    // Second loop to compute the actual gradient, stored in h_energy
    group_size = this->m_fluidgroup->getNumMembers();
    #pragma omp parallel for private(size, myHead, gradCi) firstprivate(temp0)
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int i = h_members_omp3.data[group_idx];
        
        // set temp variable to zero 
        gradCi.x = 0.0;
        gradCi.y = 0.0;
        gradCi.z = 0.0;

        // Access the particle's position
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;

        Scalar Ci = h_pressure.data[i];

        // Loop over all of the neighbors of this particle
        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];

        for (unsigned int j = 0; j < size; j++)
        {
            // Index of neighbor
            unsigned int k = h_nlist.data[myHead + j];

            // Access neighbor position
            Scalar3 pj;
            pj.x = h_pos.data[k].x;
            pj.y = h_pos.data[k].y;
            pj.z = h_pos.data[k].z;

            Scalar Cj = h_pressure.data[k];
            
            // Compute distance vector
            // Scalar3 dx = pj - pi;
            Scalar3 dx;
            dx.x = pi.x - pj.x;
            dx.y = pi.y - pj.y;
            dx.z = pi.z - pj.z;

            Scalar mj   = h_velocity.data[k].w;
            Scalar rhoj = h_density.data[k];

            // Apply periodic boundary conditions
            dx = box.minImage(dx);

            // Calculate squared distance
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, continue with next neighbor in loop
            if ( this->m_const_slength && rsq > m_rcutsq )
                continue;

            // Calculate distance
            Scalar r = sqrt(rsq);

            // Mean smoothing length
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = (r > Scalar(1e-8)*meanh) ? dwdr/r : Scalar(0);
            
            temp0 = ( Cj - Ci ) * ( mj/rhoj ); 

            gradCi.x += temp0 * dwdr_r * dx.x;
            gradCi.y += temp0 * dwdr_r * dx.y;
            gradCi.z += temp0 * dwdr_r * dx.z;

        } // End neighbour loop

        // Compute the actual squared L2 norm of the partcile contentration gradient

        h_energy.data[i] = dot( gradCi, gradCi );

    } // End fluid group loop

} // End compute particle concentration gradient


/*! Perform number density computation
 * This method computes and stores
     - the number density based mass density ( rho_i = m_i * \sum w_ij ) for fluid particles
       if the SUMMATION approach is being used.
     - the zeroth order normalization constant for solid particles
   in the density Array.
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_ndensity(uint64_t timestep)
{
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Number Density" << std::endl;

    // Grab handles for particle data
    ArrayHandle<Scalar> h_density(this->m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    unsigned int size;
    size_t myHead;
    Scalar ni; // \sum_j w_{ij}

    // Precompute self-density for constant smoothing length (avoids per-particle w0 call)
    Scalar w0 = m_const_slength ? this->m_skernel->w0(m_ch) : Scalar(0.0);

    // Particle loop
    // For each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp5(this->m_fluidgroup->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for private(size, myHead, ni)
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int i = h_members_omp5.data[group_idx];

        // Self-density contribution: use per-particle h when smoothing length is variable
        ni = m_const_slength ? w0 : this->m_skernel->w0(h_h.data[i]);

        // Access the particle's position
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;

        // Loop over all of the neighbors of this particle
        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];

        for (unsigned int j = 0; j < size; j++)
        {
            // Index of neighbor
            unsigned int k = h_nlist.data[myHead + j];

            // Access neighbor position
            Scalar3 pj;
            pj.x = h_pos.data[k].x;
            pj.y = h_pos.data[k].y;
            pj.z = h_pos.data[k].z;

            // Compute distance vector
            // Scalar3 dx = pj - pi;
            Scalar3 dx;
            dx.x = pi.x - pj.x;
            dx.y = pi.y - pj.y;
            dx.z = pi.z - pj.z;

            // Apply periodic boundary conditions
            dx = box.minImage(dx);

            // Calculate squared distance
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, continue with next neighbor in loop
            if ( this->m_const_slength && rsq > m_rcutsq )
                continue;

            // Calculate distance
            Scalar r = sqrt(rsq);

            ni += this->m_const_slength ? this->m_skernel->wij(m_ch,r) : this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);

        } // End neighbour loop

        // Compute mass density from number density if particle i is a fluid particle
        // rho_i = m_i * \sum_j wij
        h_density.data[i] = ni * h_velocity.data[i].w;

    } // End fluid group loop

} // End compute number density


/*! Perform pressure computation
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_pressure(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Pressure" << std::endl;

    // Define ArrayHandles
    ArrayHandle<Scalar> h_density(this->m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_pressure(this->m_pdata->getPressures(), access_location::host, access_mode::readwrite);

    // For each fluid particle of fluidgroup1 
    unsigned int group_size = this->m_fluidgroup1->getNumMembers();
    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp6(this->m_fluidgroup1->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int i = h_members_omp6.data[group_idx];
        // Evaluate pressure
        h_pressure.data[i] = this->m_eos1->Pressure(h_density.data[i]);
    
    } // End fluid group 1 loop

    // For each fluid particle of fluidgroup2 
    group_size = this->m_fluidgroup2->getNumMembers();
    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp7(this->m_fluidgroup2->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int i = h_members_omp7.data[group_idx];
        // Evaluate pressure
        h_pressure.data[i] = this->m_eos2->Pressure(h_density.data[i]);
    
    } // End fluid group 2 loop

} // End compute pressure



template<SmoothingKernelType KT_,StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_noslip(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::NoSlip NoPenetration" << std::endl;

    // Grab handles for particle and neighbor data
    ArrayHandle<Scalar3> h_vf(this->m_pdata->getAuxiliaries1(), access_location::host,access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_density(this->m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar>  h_pressure(this->m_pdata->getPressures(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_accel(this->m_pdata->getAccelerations(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    unsigned int size;
    size_t myHead;

    // For all solid particles
    unsigned int group_size = this->m_solidgroup->getNumMembers();
    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp8(this->m_solidgroup->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for private(size, myHead)
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = h_members_omp8.data[group_idx];

        // Access the particle's position, velocity, mass and type
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;

        // Read acceleration of solid particle i; default to zero if NaN
        Scalar3 accel_i = make_scalar3(0,0,0);
        if ( h_accel.data[i].x == h_accel.data[i].x &&
             h_accel.data[i].y == h_accel.data[i].y &&
             h_accel.data[i].z == h_accel.data[i].z )
            {
            accel_i.x = h_accel.data[i].x;
            accel_i.y = h_accel.data[i].y;
            accel_i.z = h_accel.data[i].z;
            }

        // Initialize fictitious solid velocity vector
        Scalar3 uf_c0 = make_scalar3(0, 0, 0);

        // Initialize fictitious solid pressure scalar
        Scalar pf_c0= Scalar(0);

        // Initialize hydrostatic pressure contribution
        Scalar3 ph_c0 = make_scalar3(0, 0, 0);

        // Initialize reziprocal solid particle wise zeroth order normalisation constant 
        Scalar wij_c0 = Scalar(0);

        // Loop over all of the neighbors of this particle
        // Count fluid neighbors before setting solid particle properties
        unsigned int fluidneighbors = 0;

        // Skip neighbor loop if this solid particle does not have fluid neighbors
        bool solid_w_fluid_neigh = false;
        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {
            unsigned int k = h_nlist.data[myHead + j];
            if ( checkfluid1(h_type_property_map.data, h_pos.data[k].w) ||
                 checkfluid2(h_type_property_map.data, h_pos.data[k].w))
                {
                solid_w_fluid_neigh = true;
                break;
                }
            }
        if ( !(solid_w_fluid_neigh) )
            {
            // Set fictitious solid velocity to zero
            h_vf.data[i].x = 0;
            h_vf.data[i].y = 0;
            h_vf.data[i].z = 0;
            // If no fluid neighbors are present,
            // Set pressure to background pressure
            h_pressure.data[i] = this->m_eos1->getBackgroundPressure();
            // Density to rest density
            h_density.data[i] = this->m_rho01;

            continue;
            }

        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];
        // loop over all neighbours of the solid particle
        // effectivly, only fluid particles contribute to properties of the solid

        for (unsigned int j = 0; j < size; j++)
            {
            // Index of neighbor (MEM TRANSFER: 1 scalar)
            unsigned int k = h_nlist.data[myHead + j];

            // Sanity check
            assert(k < this->m_pdata->getN() + this->m_pdata->getNGhosts());

            // If neighbor particle is solid, continue with next element in loop
            // i.e. interpolations only apply to fluid particles
            if ( checksolid(h_type_property_map.data, h_pos.data[k].w) )
                continue;
            else
                fluidneighbors += 1;

            // Access neighbor position
            Scalar3 pj;
            pj.x = h_pos.data[k].x;
            pj.y = h_pos.data[k].y;
            pj.z = h_pos.data[k].z;

            // Compute distance vector (FLOPS: 3)
            // in this case i is the solid particle, j its fluid neighbour
            Scalar3 dx;
            dx.x = pi.x - pj.x;
            dx.y = pi.y - pj.y;
            dx.z = pi.z - pj.z;

            // Apply periodic boundary conditions (FLOPS: 9)
            dx = box.minImage(dx);

            // Calculate squared distance (FLOPS: 5)
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, skip this loop
            if ( this->m_const_slength && rsq > this->m_rcutsq )
                continue;

            // Access neighbor velocity and mass
            Scalar3 vj;
            vj.x = h_velocity.data[k].x;
            vj.y = h_velocity.data[k].y;
            vj.z = h_velocity.data[k].z;

            // Read particle k pressure
            Scalar Pj = h_pressure.data[k];

            // Calculate absolute and normalized distance
            Scalar r = sqrt(rsq);

            // Evaluate kernel function
            Scalar wij = this->m_const_slength ? this->m_skernel->wij(m_ch,r) : this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);

            // Add contribution to solid fictitious velocity
            uf_c0.x += vj.x*wij;
            uf_c0.y += vj.y*wij;
            uf_c0.z += vj.z*wij;

            // Add contribution to solid fictitious pressure
            pf_c0 += Pj*wij;

            // Add contribution to hydrostatic pressure term
            // this also includes a direction (included in dx)
            // h_density is the density of the fluid and therefore a real density
            ph_c0.x += h_density.data[k] * dx.x * wij;
            ph_c0.y += h_density.data[k] * dx.y * wij;
            ph_c0.z += h_density.data[k] * dx.z * wij;

            wij_c0 += wij;

            } // End neighbor loop

        // Store fictitious solid particle velocity
        if (fluidneighbors > 0 && wij_c0 > 0 )
            {
            Scalar norm_constant = 1./wij_c0;
            // Set fictitious velocity
            h_vf.data[i].x = 2.0 * h_velocity.data[i].x - norm_constant * uf_c0.x;
            h_vf.data[i].y = 2.0 * h_velocity.data[i].y - norm_constant * uf_c0.y;
            h_vf.data[i].z = 2.0 * h_velocity.data[i].z - norm_constant * uf_c0.z;
            // compute fictitious pressure
            // TODO: There is an addition necessary if the acceleration of the solid 
            // phase is not constant, since there is no function that is updating it
            // see ISSUE # 23
            Scalar3 bodyforce = this->getAcceleration(timestep);
            Scalar3 hp_factor;
            hp_factor.x = bodyforce.x - accel_i.x;
            hp_factor.y = bodyforce.y - accel_i.y;
            hp_factor.z = bodyforce.z - accel_i.z;

            ph_c0.x *= norm_constant;
            ph_c0.y *= norm_constant;
            ph_c0.z *= norm_constant;

            h_pressure.data[i] = norm_constant * pf_c0 + dot(hp_factor , ph_c0);
            // Compute solid densities by inverting equation of state
            // Here: overwrite the normalisation constant
                        // If interpolated solid pressure is negative, set to background pressure
            if ( h_pressure.data[i] < 0 )
                {
                // Set pressure to background pressure
                h_pressure.data[i] = this->m_eos1->getBackgroundPressure();
                // Set Density to rest density
                h_density.data[i] = this->m_rho01;
                }
            else 
                {
                // Compute solid densities by inverting equation of state
                h_density.data[i] = this->m_eos1->Density(h_pressure.data[i]);
                }
            }
        else
            {
            // Set fictitious solid velocity to zero
            h_vf.data[i].x = 0.0;
            h_vf.data[i].y = 0.0;
            h_vf.data[i].z = 0.0;

            // If no fluid neighbors are present,
            // Set pressure to background pressure
            h_pressure.data[i] = this->m_eos1->getBackgroundPressure();
            // Density to rest density
            h_density.data[i] = this->m_rho01;
            }

        } // End solid particle loop

    } // End compute noslip computation


// TODO : THIS has still to be checked
/*! Perform Shepard density renormalization
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::renormalize_density(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Density renormalization" << std::endl;

    // Grab handles for particle data
    ArrayHandle<Scalar>  h_density(this->m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);

    auto tmp_density = this->m_pdata->getDensities();
    ArrayHandle<Scalar> h_density_old(tmp_density, access_location::host, access_mode::read);


    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Precompute self-density for homogeneous smoothing lengths
    Scalar w0 = this->m_skernel->w0(this->m_ch);

    unsigned int size;
    size_t myHead;


    // Particle loop
    // For each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp9(this->m_fluidgroup->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for private(size, myHead)
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = h_members_omp9.data[group_idx];

        // Access the particle's position
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;
        Scalar mi = h_velocity.data[i].w;
        Scalar rhoi = h_density.data[i];

        // First compute renormalization factor
        // Initialize with self density of kernel
        Scalar normalization = this->m_const_slength ? w0 : this->m_skernel->w0(h_h.data[i]);
        normalization = normalization * ( mi / rhoi );

        // Loop over all of the neighbors of this particle
        // and compute normalization constant normwij = \sum_j wij*Vj
        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {
                // Index of neighbor
                unsigned int k = h_nlist.data[myHead + j];

                // Access neighbor position
                Scalar3 pj;
                pj.x = h_pos.data[k].x;
                pj.y = h_pos.data[k].y;
                pj.z = h_pos.data[k].z;

                // Compute distance vector
                // Scalar3 dx = pj - pi;
                Scalar3 dx;
                dx.x = pj.x - pi.x;
                dx.y = pj.y - pi.y;
                dx.z = pj.z - pi.z;

                // Apply periodic boundary conditions
                dx = box.minImage(dx);

                // Calculate squared distance
                Scalar rsq = dot(dx, dx);

                // If particle distance is too large, continue with next neighbor in loop
                if ( this->m_const_slength && rsq > this->m_rcutsq )
                    continue;

                // Calculate distance
                Scalar r = sqrt(rsq);

                // Add contribution to renormalization
                Scalar Vj =  h_velocity.data[k].w / h_density_old.data[k] ;
                normalization += this->m_const_slength ? Vj*this->m_skernel->wij(m_ch,r) : Vj*this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);

            } // End of neighbor loop

        normalization = Scalar(1.0)/normalization;

        // Initialize density with normalized kernel self density
        h_density.data[i] = this->m_const_slength ? w0*(mi*normalization): this->m_skernel->w0(h_h.data[i])*(mi*normalization);

        // Loop over all of the neighbors of this particle
        // and compute renormalied density rho_i = \sum_j wij*mj / normwij
        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {
            // Index of neighbor
            unsigned int k = h_nlist.data[myHead + j];

            // Access neighbor position
            Scalar3 pj;
            pj.x = h_pos.data[k].x;
            pj.y = h_pos.data[k].y;
            pj.z = h_pos.data[k].z;

            // Compute distance vector
            // Scalar3 dx = pj - pi;
            Scalar3 dx;
            dx.x = pj.x - pi.x;
            dx.y = pj.y - pi.y;
            dx.z = pj.z - pi.z;

            // Apply periodic boundary conditions
            dx = box.minImage(dx);

            // Calculate squared distance
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, continue with next neighbor in loop
            if ( this->m_const_slength && rsq > this->m_rcutsq )
                continue;

            // Calculate distance
            Scalar r = sqrt(rsq);

            // Add contribution to normalized density interpolation
            Scalar factor =  h_velocity.data[k].w * normalization ;
            h_density.data[i] += this->m_const_slength ? factor*this->m_skernel->wij(m_ch,r) : factor*this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);
            }
        } // End of particle loop
    } // End renormalize density



/*! Compute interfacial color gradients
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_colorgradients(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Normals/ColorGradient" << std::endl;

    // Grab handles for particle and neighbor data
    ArrayHandle<Scalar3> h_sn(this->m_pdata->getAuxiliaries2(), access_location::host,access_mode::readwrite);
    ArrayHandle<Scalar3> h_fn(this->m_pdata->getAuxiliaries3(), access_location::host,access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_density(this->m_pdata->getDensities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_energy(this->m_pdata->getEnergies(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Zero data before calculation
    memset((void*)h_sn.data,0,sizeof(Scalar3)*this->m_pdata->getAuxiliaries2().getNumElements());
    memset((void*)h_fn.data,0,sizeof(Scalar3)*this->m_pdata->getAuxiliaries3().getNumElements());

    // Particle loop
    #pragma omp parallel for
    for (unsigned int i = 0; i < this->m_pdata->getN(); i++)
        {
        // Access the particle's position, mass and type
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;
        Scalar mi = h_velocity.data[i].w;

        // Read particle i density and volume
        Scalar rhoi = h_density.data[i];
        Scalar Vi   = mi / rhoi;

        // Detect particle i type
        bool i_issolid = checksolid(h_type_property_map.data, h_pos.data[i].w);
        bool i_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[i].w);
        bool i_isfluid2 = checkfluid2(h_type_property_map.data, h_pos.data[i].w);

        // Loop over all of the neighbors of this particle
        size_t myHead = h_head_list.data[i];
        const unsigned int size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {
            // Index of neighbor (MEM TRANSFER: 1 scalar)
            unsigned int k = h_nlist.data[myHead + j];

            // Sanity check
            assert(k < this->m_pdata->getN() + this->m_pdata->getNGhosts());

            // Access neighbor position
            Scalar3 pj;
            pj.x = h_pos.data[k].x;
            pj.y = h_pos.data[k].y;
            pj.z = h_pos.data[k].z;

            // Determine neighbor type
            bool j_issolid  = checksolid(h_type_property_map.data, h_pos.data[k].w);
            bool j_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[k].w);
            bool j_isfluid2 = checkfluid2(h_type_property_map.data, h_pos.data[k].w);

            // Skip color gradient computation if both particles belong to same phase
            if (    ( i_issolid  && j_issolid  ) 
                 || ( i_isfluid1 && j_isfluid1 ) 
                 || ( i_isfluid2 && j_isfluid2 ) )
                continue;

            // Compute distance vector (FLOPS: 3)
            Scalar3 dx = pi - pj;

            // Apply periodic boundary conditions (FLOPS: 9)
            dx = box.minImage(dx);

            // Calculate squared distance (FLOPS: 5)
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, skip this loop
            if ( this->m_const_slength && rsq > this->m_rcutsq )
                continue;

            // Access neighbor mass and density
            Scalar mj   = h_velocity.data[k].w;
            Scalar rhoj = h_density.data[k];
            Scalar Vj   = mj / rhoj;

            // Calculate absolute and normalized distance
            Scalar r = sqrt(rsq);

            // Mean smoothing length
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = (r > Scalar(1e-8)*meanh) ? dwdr/r : Scalar(0);


            Scalar temp0 = 0.0;

            if ( m_colorgradient_method == DENSITYRATIO )
            {
                // Adami type color gradient, also implemented in PySPH
                temp0 = rhoi/( rhoi + rhoj ) * (Vi*Vi + Vj*Vj)/Vi;
            }

            else if ( m_colorgradient_method == NUMBERDENSITY )
            {
                temp0 = (Vj*Vj/Vi);
            }
            else {
                throw std::runtime_error("Error: No valid ColorGradientMethod given.");
            }

            // If either on of the particle is a solid, interface must be solid-fluid
            if ( i_issolid || j_issolid )
            {
                h_sn.data[i].x += temp0*dwdr_r*dx.x;
                h_sn.data[i].y += temp0*dwdr_r*dx.y;
                h_sn.data[i].z += temp0*dwdr_r*dx.z;
            }
            // Otherwise, interface must be fluid-fluid
            else
            {
                h_fn.data[i].x += temp0*dwdr_r*dx.x;
                h_fn.data[i].y += temp0*dwdr_r*dx.y;
                h_fn.data[i].z += temp0*dwdr_r*dx.z;
            }

            } // Closing Neighbor Loop

        // Make sure that color gradients point from solid to fluid
        // and from fluid 1 to fluid 2 (affects sign of normals)
        if ( i_issolid )
            {
                h_sn.data[i].x = -h_sn.data[i].x;
                h_sn.data[i].y = -h_sn.data[i].y;
                h_sn.data[i].z = -h_sn.data[i].z;
            }
        if ( i_isfluid1 )
            {
                h_fn.data[i].x = -h_fn.data[i].x;
                h_fn.data[i].y = -h_fn.data[i].y;
                h_fn.data[i].z = -h_fn.data[i].z;
            }

        } // End of particle loop

    // Normal smoothing pass (Adami et al. 2010):
    // Smooth raw color gradients by Shepard-renormalized weighted average of neighbor normals.
    // This reduces parasitic currents at the fluid-fluid interface significantly.
    const Scalar eps_norm = Scalar(1e-6);

    // Fluid-fluid normals (h_fn stored in aux3)
    unsigned int fluid_size = this->m_fluidgroup->getNumMembers();
    // Allocate temporary buffer for smoothed fluid normals
    std::vector<Scalar3> fn_smooth(this->m_pdata->getN() + this->m_pdata->getNGhosts(),
                                   make_scalar3(0.0, 0.0, 0.0));

    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp10(this->m_fluidgroup->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for
    for (unsigned int group_idx = 0; group_idx < fluid_size; group_idx++)
        {
        unsigned int i = h_members_omp10.data[group_idx];

        Scalar norm_i = sqrt(dot(h_fn.data[i], h_fn.data[i]));
        if ( norm_i < eps_norm )
            {
            fn_smooth[i] = h_fn.data[i];
            continue;
            }

        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;
        Scalar mi = h_velocity.data[i].w;
        Scalar rhoi = h_density.data[i];

        Scalar3 acc_fn = make_scalar3(0.0, 0.0, 0.0);
        Scalar  w_acc  = Scalar(0.0);

        // Self-contribution
        Scalar Vi = mi / rhoi;
        Scalar w0_i = this->m_const_slength ? this->m_skernel->w0(m_ch) : this->m_skernel->w0(h_h.data[i]);
        acc_fn.x += Vi * h_fn.data[i].x * w0_i;
        acc_fn.y += Vi * h_fn.data[i].y * w0_i;
        acc_fn.z += Vi * h_fn.data[i].z * w0_i;
        w_acc    += Vi * w0_i;

        size_t myHead_s = h_head_list.data[i];
        unsigned int sz = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < sz; j++)
            {
            unsigned int k = h_nlist.data[myHead_s + j];

            // Only smooth over fluid neighbors
            if ( checksolid(h_type_property_map.data, h_pos.data[k].w) )
                continue;

            Scalar3 pj;
            pj.x = h_pos.data[k].x;
            pj.y = h_pos.data[k].y;
            pj.z = h_pos.data[k].z;

            Scalar3 dx_s;
            dx_s.x = pi.x - pj.x;
            dx_s.y = pi.y - pj.y;
            dx_s.z = pi.z - pj.z;
            dx_s = box.minImage(dx_s);

            Scalar rsq_s = dot(dx_s, dx_s);
            if ( this->m_const_slength && rsq_s > this->m_rcutsq )
                continue;

            Scalar r_s = sqrt(rsq_s);
            Scalar meanh_s = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar wij_s = this->m_skernel->wij(meanh_s, r_s);

            Scalar mj   = h_velocity.data[k].w;
            Scalar rhoj = h_density.data[k];
            Scalar Vk   = mj / rhoj;

            acc_fn.x += Vk * h_fn.data[k].x * wij_s;
            acc_fn.y += Vk * h_fn.data[k].y * wij_s;
            acc_fn.z += Vk * h_fn.data[k].z * wij_s;
            w_acc    += Vk * wij_s;
            }

        if ( w_acc > eps_norm )
            {
            Scalar inv_w = Scalar(1.0) / w_acc;
            fn_smooth[i].x = acc_fn.x * inv_w;
            fn_smooth[i].y = acc_fn.y * inv_w;
            fn_smooth[i].z = acc_fn.z * inv_w;
            }
        else
            fn_smooth[i] = h_fn.data[i];
        }

    // Write back smoothed fluid normals
    #pragma omp parallel for
    for (unsigned int group_idx = 0; group_idx < fluid_size; group_idx++)
        {
        unsigned int i = h_members_omp10.data[group_idx];
        h_fn.data[i] = fn_smooth[i];
        }

    } // End compute colorgradients



/*! Compute surface force vectors
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_surfaceforce(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::SurfaceForce" << std::endl;

    // Grab handles for particle and neighbor data
    ArrayHandle<Scalar3> h_sf(this->m_pdata->getAuxiliaries4(), access_location::host,access_mode::readwrite);
    ArrayHandle<Scalar3> h_sn(this->m_pdata->getAuxiliaries2(), access_location::host,access_mode::read);
    ArrayHandle<Scalar3> h_fn(this->m_pdata->getAuxiliaries3(), access_location::host,access_mode::read);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_density(this->m_pdata->getDensities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_energy(this->m_pdata->getEnergies(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);


    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Zero data before calculation
    memset((void*)h_sf.data,0,sizeof(Scalar3)*this->m_pdata->getAuxiliaries4().getNumElements());

    // for each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp12(this->m_fluidgroup->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = h_members_omp12.data[group_idx];

        // Access the particle's position and type
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;
        bool i_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[i].w);
        bool i_isfluid2 = checkfluid2(h_type_property_map.data, h_pos.data[i].w);

        // Check if there is any fluid particle near the current particle, if not continue
        // This makes sure that only particle near a fluid interface experience an interfacial force.
        // In other words, fluid particles only near solid interfaces are omitted.
        bool nearfluidinterface = false;

        // Loop over all of the neighbors of this particle
        size_t myHead = h_head_list.data[i];
        unsigned int size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {
            // Index of neighbor (MEM TRANSFER: 1 scalar)
            unsigned int k = h_nlist.data[myHead + j];
            assert(k < this->m_pdata->getN() + this->m_pdata->getNGhosts());
            bool j_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[k].w);
            bool j_isfluid2 = checkfluid2(h_type_property_map.data, h_pos.data[k].w);
            // Near fluid interface if i is fluid1 and j is fluid2, or i is fluid2 and j is fluid1
            if ( (i_isfluid1 && j_isfluid2) || (i_isfluid2 && j_isfluid1) )
                {
                    nearfluidinterface = true;
                    break;
                }
            }
        // Wetting fix (2026-08-21): when a solid-fluid tension is active
        // (omega != 90), particles in the wall film also carry interfacial
        // stress (sigma01/sigma02, Huber-style). The original fluid-fluid-only
        // gate zeroed their stress, so the contact-line divergence saw
        // inconsistent neighbor stresses and the Young traction along the
        // wall vanished -- observed as zero capillary rise at any omega.
        bool nearsolid = false;
        if ( m_sigma01 > 0.0 || m_sigma02 > 0.0 )
            {
            Scalar3 sni_gate;
            sni_gate.x = h_sn.data[i].x;
            sni_gate.y = h_sn.data[i].y;
            sni_gate.z = h_sn.data[i].z;
            nearsolid = dot(sni_gate, sni_gate) > Scalar(0);
            }
        if ( !nearfluidinterface && !nearsolid )
            continue;

        // Access the particle's mass
        Scalar mi = h_velocity.data[i].w;

        // Read particle i density and volume
        Scalar rhoi = h_density.data[i];
        Scalar Vi   = mi / rhoi;

        // Read particle i color gradients
        Scalar3 sni;
        sni.x = h_sn.data[i].x;
        sni.y = h_sn.data[i].y;
        sni.z = h_sn.data[i].z;
        Scalar normsni = sqrt(dot(sni,sni));
        Scalar3 fni;
        fni.x = h_fn.data[i].x;
        fni.y = h_fn.data[i].y;
        fni.z = h_fn.data[i].z;
        Scalar normfni = sqrt(dot(fni,fni));


        // Evaluate particle i interfacial stress tensor
        Scalar istress[6] = {0};
        Scalar temp0 = 0.0;
        Scalar temp1 = 0.0;
        // Get particle Concentration gradient (Shifting)
        // Spactial dimension d = 3
        if ( m_fickian_shifting )
            {
            temp1 = 1./3. * h_energy.data[i];
            }
        else 
            {
            temp1 = 1./3. * normfni * normfni;
            }

        // normal vectors point from solid to fluid and from fluid 1
        // to fluid 2
        // if Fluid1 or Fluid2 that has neighbors of other fluid phase 
        if ( this->m_sigma12 > 0.0 && normfni > 0.0 )
        {
            temp0 = this->m_sigma12/normfni;
            istress[0] += temp0 * ( temp1 - fni.x * fni.x); // xx
            istress[1] += temp0 * ( temp1 - fni.y * fni.y); // yy
            istress[2] += temp0 * ( temp1 - fni.z * fni.z); // zz
            istress[3] -= temp0 * ( fni.x * fni.y);         // xy yx
            istress[4] -= temp0 * ( fni.x * fni.z);         // xz zx
            istress[5] -= temp0 * ( fni.y * fni.z);         // yz zy
        }

        if ( !m_fickian_shifting )
        {
            temp1 = 1./3. * normsni * normsni;
        }

        // --- hysteresis block for particle i ---
        Scalar sigma01_i = this->m_sigma01;
        Scalar sigma02_i = this->m_sigma02;
        if (m_hysteresis && normsni > 0.0 && normfni > 0.0)
        {
            Scalar cos_local = dot(fni, sni) / (normfni * normsni);
            cos_local = fmax(Scalar(-1), fmin(Scalar(1), cos_local));
            Scalar theta_local = acos(cos_local) * (Scalar(180) / M_PI);
            Scalar omega_eff = fmax(m_omega_rec, fmin(m_omega_adv, theta_local));
            if      (omega_eff == Scalar(90)) { sigma01_i = 0; sigma02_i = 0; }
            else if (omega_eff <  Scalar(90)) { sigma01_i = this->m_sigma12 * cos(omega_eff*(M_PI/180)); sigma02_i = 0; }
            else                              { sigma01_i = 0; sigma02_i = this->m_sigma12 * cos((180-omega_eff)*(M_PI/180)); }
        }

        // Fluid phase 1 - Solid interface
        if ( i_isfluid1 && sigma01_i > 0.0 && normsni > 0.0 )
        {
            temp0 = sigma01_i/normsni;
            istress[0] += temp0 * ( temp1 - sni.x * sni.x); // xx
            istress[1] += temp0 * ( temp1 - sni.y * sni.y); // yy
            istress[2] += temp0 * ( temp1 - sni.z * sni.z); // zz
            istress[3] -= temp0 * ( sni.x * sni.y);         // xy yx
            istress[4] -= temp0 * ( sni.x * sni.z);         // xz zx
            istress[5] -= temp0 * ( sni.y * sni.z);         // yz zy
        }

        // Fluid phase 2 - Solid interface
        if ( i_isfluid2 && sigma02_i > 0.0 && normsni > 0.0 )
        {
            temp0 = sigma02_i/normsni;
            istress[0] += temp0 * ( temp1 - sni.x * sni.x); // xx
            istress[1] += temp0 * ( temp1 - sni.y * sni.y); // yy
            istress[2] += temp0 * ( temp1 - sni.z * sni.z); // zz
            istress[3] -= temp0 * ( sni.x * sni.y);         // xy yx
            istress[4] -= temp0 * ( sni.x * sni.z);         // xz zx
            istress[5] -= temp0 * ( sni.y * sni.z);         // yz zy
        }

        // Loop over all of the neighbors of this particle
        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {

            // Index of neighbor (MEM TRANSFER: 1 scalar)
            unsigned int k = h_nlist.data[myHead + j];

            // Sanity check
            assert(k < this->m_pdata->getN() + this->m_pdata->getNGhosts());

            // Access neighbor position
            Scalar3 pj;
            pj.x = h_pos.data[k].x;
            pj.y = h_pos.data[k].y;
            pj.z = h_pos.data[k].z;

            // Determine neighbor type
            bool j_issolid  = checksolid(h_type_property_map.data, h_pos.data[k].w);
            bool j_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[k].w);
            bool j_isfluid2 = checkfluid2(h_type_property_map.data, h_pos.data[k].w);

            // Compute normalized color gradients
            Scalar3 snj;
            snj.x = h_sn.data[k].x;
            snj.y = h_sn.data[k].y;
            snj.z = h_sn.data[k].z;
            Scalar normsnj = sqrt(dot(snj,snj));
            Scalar3 fnj;
            fnj.x = h_fn.data[k].x;
            fnj.y = h_fn.data[k].y;
            fnj.z = h_fn.data[k].z;
            Scalar normfnj = sqrt(dot(fnj,fnj));

            // Compute distance vector (FLOPS: 3)
            Scalar3 dx = pi - pj;

            // Apply periodic boundary conditions (FLOPS: 9)
            dx = box.minImage(dx);

            // Calculate squared distance (FLOPS: 5)
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, skip this loop
            if ( this->m_const_slength && rsq > this->m_rcutsq )
                continue;

            // Calculate absolute and normalized distance
            Scalar r = sqrt(rsq);

            // Access neighbor mass and density
            Scalar mj   = h_velocity.data[k].w;
            Scalar rhoj = h_density.data[k];
            Scalar Vj   = mj / rhoj;

            // Mean smoothing length
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = (r > Scalar(1e-8)*meanh) ? dwdr/r : Scalar(0);

            // temp0 = 0.0;
            // temp1 = 0.0;
            // Get particle Concentration gradient (Shifting)
            // Spactial dimension d = 3
            if ( m_fickian_shifting )
                {
                temp1 = 1./3. * h_energy.data[k];
                }
            else 
                {
                temp1 = 1./3. * normfnj * normfnj;
                }

            // Evaluate particle i interfacial stress tensor
            Scalar jstress[6] = {0};
            // normal vectors point from solid to fluid and from fluid 1
            // to fluid 2
            // if Fluid1 or Fluid2 that has neighbors of other fluid phase 
            if ( !(j_issolid) && this->m_sigma12 > 0.0 && normfnj > 0.0 )
            {
                temp0 = this->m_sigma12/normfnj;
                jstress[0] += temp0 * ( temp1 - fnj.x * fnj.x); // xx
                jstress[1] += temp0 * ( temp1 - fnj.y * fnj.y); // yy
                jstress[2] += temp0 * ( temp1 - fnj.z * fnj.z); // zz
                jstress[3] -= temp0 * ( fnj.x * fnj.y);         // xy yx
                jstress[4] -= temp0 * ( fnj.x * fnj.z);         // xz zx
                jstress[5] -= temp0 * ( fnj.y * fnj.z);         // yz zy
            }

            if ( !m_fickian_shifting )
            {
                temp1 = 1./3. * normsnj * normsnj;
            }

            // --- hysteresis block for particle j ---
            Scalar sigma01_j = this->m_sigma01;
            Scalar sigma02_j = this->m_sigma02;
            if (m_hysteresis && normsnj > 0.0 && normfnj > 0.0)
            {
                Scalar cos_local_j = dot(fnj, snj) / (normfnj * normsnj);
                cos_local_j = fmax(Scalar(-1), fmin(Scalar(1), cos_local_j));
                Scalar theta_local_j = acos(cos_local_j) * (Scalar(180) / M_PI);
                Scalar omega_eff_j = fmax(m_omega_rec, fmin(m_omega_adv, theta_local_j));
                if      (omega_eff_j == Scalar(90)) { sigma01_j = 0; sigma02_j = 0; }
                else if (omega_eff_j <  Scalar(90)) { sigma01_j = this->m_sigma12 * cos(omega_eff_j*(M_PI/180)); sigma02_j = 0; }
                else                                { sigma01_j = 0; sigma02_j = this->m_sigma12 * cos((180-omega_eff_j)*(M_PI/180)); }
            }

            // Fluid phase 1 - Solid interface
            if ( j_isfluid1 && sigma01_j > 0.0 && normsnj > 0.0 )
            {
                temp0 = sigma01_j/normsnj;
                jstress[0] += temp0 * ( temp1 - snj.x * snj.x); // xx
                jstress[1] += temp0 * ( temp1 - snj.y * snj.y); // yy
                jstress[2] += temp0 * ( temp1 - snj.z * snj.z); // zz
                jstress[3] -= temp0 * ( snj.x * snj.y);         // xy yx
                jstress[4] -= temp0 * ( snj.x * snj.z);         // xz zx
                jstress[5] -= temp0 * ( snj.y * snj.z);         // yz zy
            }

            // Fluid phase 2 - Solid interface
            if ( j_isfluid2 && sigma02_j > 0.0 && normsnj > 0.0 )
            {
                temp0 = sigma02_j/normsnj;
                jstress[0] += temp0 * ( temp1 - snj.x * snj.x); // xx
                jstress[1] += temp0 * ( temp1 - snj.y * snj.y); // yy
                jstress[2] += temp0 * ( temp1 - snj.z * snj.z); // zz
                jstress[3] -= temp0 * ( snj.x * snj.y);         // xy yx
                jstress[4] -= temp0 * ( snj.x * snj.z);         // xz zx
                jstress[5] -= temp0 * ( snj.y * snj.z);         // yz zy
            }

            // Add contribution to surface force (volume-squared formulation, anti-symmetric)
            h_sf.data[i].x += dwdr_r*dx.x*(Vi*Vi*istress[0]+Vj*Vj*jstress[0])+
                              dwdr_r*dx.y*(Vi*Vi*istress[3]+Vj*Vj*jstress[3])+
                              dwdr_r*dx.z*(Vi*Vi*istress[4]+Vj*Vj*jstress[4]);
            h_sf.data[i].y += dwdr_r*dx.x*(Vi*Vi*istress[3]+Vj*Vj*jstress[3])+
                              dwdr_r*dx.y*(Vi*Vi*istress[1]+Vj*Vj*jstress[1])+
                              dwdr_r*dx.z*(Vi*Vi*istress[5]+Vj*Vj*jstress[5]);
            h_sf.data[i].z += dwdr_r*dx.x*(Vi*Vi*istress[4]+Vj*Vj*jstress[4])+
                              dwdr_r*dx.y*(Vi*Vi*istress[5]+Vj*Vj*jstress[5])+
                              dwdr_r*dx.z*(Vi*Vi*istress[2]+Vj*Vj*jstress[2]);




            } // End of neighbor loop

        // Set component normal to solid surface at solid interface to zero

        } // Closing Fluid Particle Loop


    } // End compute surface force



template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::forcecomputation(uint64_t timestep)
    {

    if ( m_density_method == DENSITYSUMMATION )
        this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Forces using SUMMATION approach " << m_density_method << endl;
    else if ( m_density_method == DENSITYCONTINUITY )
        this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Forces using CONTINUITY approach " << m_density_method << endl;

    { // Begin GPU Array Scope
    // Grab handles for particle data
    // Access mode overwrite implies that data does not need to be read in
    ArrayHandle<Scalar4> h_force(this->m_force,access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_ratedpe(this->m_ratedpe,access_location::host, access_mode::readwrite);

    // Check input data, can be omitted if need be
    assert(h_force.data);
    assert(h_ratedpe.data);

    // Zero data before force calculation
    memset((void*)h_force.data,0,sizeof(Scalar4)*this->m_force.getNumElements());
    memset((void*)h_ratedpe.data,0,sizeof(Scalar4)*this->m_ratedpe.getNumElements());

    // access the particle data
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_density(this->m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar>  h_pressure(this->m_pdata->getPressures(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_vf(this->m_pdata->getAuxiliaries1(), access_location::host,access_mode::read);
    ArrayHandle<Scalar>  h_gdot(this->m_pdata->getEnergies(), access_location::host, access_mode::read); // per-particle shear rate (non-Newtonian)
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_sf(this->m_pdata->getAuxiliaries4(), access_location::host,access_mode::read);

    // access the neighbor list
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

    // Body-force acceleration for consistent interface pressure (Hu & Adams 2009).
    // Fetched once before the loops so the per-pair dot(gvec, dx) is cheap.
    // When CIP is disabled gvec = 0 and the conditional branch is never taken,
    // so there is no run-time cost in the common case.
    const Scalar3 gvec = m_consistent_interface_pressure
                         ? this->getAcceleration(timestep)
                         : make_scalar3(Scalar(0), Scalar(0), Scalar(0));

    // maximum velocity variable for adaptive timestep
    double max_vel = 0.0;

    // for each fluid particle
    unsigned int group_size = m_fluidgroup->getNumMembers();
    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp13(m_fluidgroup->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for private(size, myHead) firstprivate(temp0) reduction(max:max_vel)
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = h_members_omp13.data[group_idx];

        // Access the particle's position, velocity, mass and type
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
        // bool i_isfluid2 = checkfluid2(h_type_property_map.data, h_pos.data[i].w);
        Scalar mui   = i_isfluid1 ? this->m_mu1 : this->m_mu2;
        Scalar rho0i = i_isfluid1 ? this->m_rho01 : this->m_rho02;
        Scalar ci    = i_isfluid1 ? this->m_c1 : this->m_c2;

        // Properties needed for adaptive timestep
        // Total velocity of particle
        Scalar vi_total = sqrt((vi.x * vi.x) + (vi.y * vi.y) + (vi.z * vi.z));
        if (vi_total > max_vel) { max_vel = vi_total; }

        // Loop over all of the neighbors of this particle
        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {
            // Index of neighbor (MEM TRANSFER: 1 scalar)
            unsigned int k = h_nlist.data[myHead + j];

            // Sanity check
            assert(k < this->m_pdata->getN() + this->m_pdata->getNGhosts());

            // Access neighbor position
            Scalar3 pj;
            pj.x = h_pos.data[k].x;
            pj.y = h_pos.data[k].y;
            pj.z = h_pos.data[k].z;

            // Determine neighbor type
            bool j_issolid  = checksolid(h_type_property_map.data, h_pos.data[k].w);
            bool j_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[k].w);

            // Read particle j viscosity, speed of sound and rest density
            Scalar muj   = j_isfluid1 ? this->m_mu1 : this->m_mu2;
            Scalar rho0j = j_isfluid1 ? this->m_rho01 : this->m_rho02;
            Scalar cj    = j_isfluid1 ? this->m_c1 : this->m_c2;
            // If particle j is solid, set parameters equal to those of particle i
            muj   = j_issolid ? mui : muj;
            rho0j = j_issolid ? rho0i : rho0j;
            cj    = j_issolid ? ci : cj;

            // Compute distance vector (FLOPS: 3)
            // Scalar3 dx = pi - pj;
            Scalar3 dx;
            dx.x = pi.x - pj.x;
            dx.y = pi.y - pj.y;
            dx.z = pi.z - pj.z;

            // Apply periodic boundary conditions (FLOPS: 9)
            dx = box.minImage(dx);

            // Calculate squared distance (FLOPS: 5)
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, skip this loop
            if ( this->m_const_slength && rsq > this->m_rcutsq )
                continue;

            // Access neighbor velocity; depends on fluid or fictitious solid particle
            Scalar3 vj  = make_scalar3(0.0, 0.0, 0.0);
            Scalar mj   = h_velocity.data[k].w;
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

            // Read particle k pressure
            Scalar Pj = h_pressure.data[k];

            // Compute velocity difference
            Scalar3 dv;
            dv.x = vi.x - vj.x;
            dv.y = vi.y - vj.y;
            dv.z = vi.z - vj.z;

            // Calculate absolute and normalized distance
            Scalar r = sqrt(rsq);

            // Mean smoothing length and denominator modifier
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar eps    = Scalar(0.1)*meanh;
            Scalar epssqr = eps*eps;

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = (r > Scalar(1e-8)*meanh) ? dwdr/r : Scalar(0);

            // ── Inter-particle pressure force ────────────────────────────────────
            // Symmetric volume formulation (Adami et al. 2013):
            //   $F_i^p = -\sum_j (V_i^2 + V_j^2) \cdot \bar{p}_{ij} \cdot (\partial W/\partial r / r) \cdot \mathbf{r}_{ij}$
            //
            // DENSITYSUMMATION -- density-weighted average pressure (Colagrossi 2003):
            //   $\bar{p}_{ij} = (\rho_j p_i + \rho_i p_j) / (\rho_i + \rho_j)$
            //
            //   With consistent interface pressure (CIP, Hu & Adams 2009), cross-phase
            //   pairs use rest-density weighting + hydrostatic correction:
            //     $\bar{p}_{ij} = (\rho_{0j} p_i + \rho_{0i} p_j + \rho_{0i} \rho_{0j} (\mathbf{g} \cdot \mathbf{r}_{ij})) / (\rho_{0i} + \rho_{0j})$
            //   The g.r_ij term makes the SPH pressure-gradient force reproduce gravity
            //   exactly for a hydrostatic column, eliminating parasitic interfacial
            //   currents that are amplified at large density ratios (e.g. water/air).
            //
            // DENSITYCONTINUITY -- mass-flux-consistent form:
            //   prefactor = m_i m_j;   $\bar{p}_{ij} = (p_i + p_j) / (\rho_i \rho_j)$
            Scalar prefactor = 0.0;
            if ( this->m_density_method == DENSITYSUMMATION )
            {
                // Consistent interface pressure (Hu & Adams 2009): rest-density weighting
                // + hydrostatic correction for cross-phase pairs; standard formula otherwise
                if ( m_consistent_interface_pressure && !j_issolid && (i_isfluid1 != j_isfluid1) )
                    temp0 = (rho0j*Pi + rho0i*Pj + rho0i*rho0j*dot(gvec, dx)) / (rho0i + rho0j);
                else
                    temp0 = (rhoj*Pi+rhoi*Pj)/(rhoi+rhoj);
                prefactor = Vi*Vi + Vj*Vj;
            }
            else if ( this->m_density_method == DENSITYCONTINUITY )
            {
                temp0 = (Pi+Pj)/(rhoi*rhoj);
                prefactor = mi * mj;
            }

            // ── Momentum dissipation (fluid–fluid pairs only) ────────────────────
            // Exactly one branch is active at a time (else-if).
            //
            // [A] Monaghan artificial viscosity (Monaghan 1992):
            //     $\Pi_{ij} = (-\alpha c_\mathrm{max} \mu_{ij} + \beta \mu_{ij}^2) / \bar{\rho}_{ij}$
            //     $\mu_{ij} = \bar{h} (\mathbf{v}_{ij} \cdot \mathbf{r}_{ij}) / (r_{ij}^2 + \eta^2)$   [has units of velocity]
            //   Activated via activateArtificialViscosity(alpha, beta).
            //
            // [B] Riemann-based dissipation (Zhang, Hu & Adams 2017):
            //     $Z^*_{ij} = Z_i Z_j / (Z_i + Z_j)$,  $Z = \rho c$   [harmonic mean impedance]
            //     $u_{ij}  = (\mathbf{v}_{ij} \cdot \mathbf{r}_{ij}) / (|\mathbf{r}_{ij}| + \eta)$       [signed radial velocity]
            //     $\mathrm{avc} = -\beta_R \cdot Z^*_{ij} \cdot u_{ij}^- / \bar{\rho}_{ij}$        (only if $\mathbf{v}_{ij} \cdot \mathbf{r}_{ij} < 0$)
            //   Impedance mismatch at the interface is handled automatically:
            //   $Z^* \to Z_\mathrm{lighter} / 2$ when $Z_\mathrm{heavy} \gg Z_\mathrm{lighter}$ (e.g. water/air).
            //   Activated via activateRiemannDissipation(beta).
            // Dissipation term "diss" carries its own scaling so that its units
            // match a force regardless of the density-method branch:
            //   [A] Monaghan AV:  Pi_ij ~ pressure/rho^2  ->  F -= m_i m_j Pi_ij grad W
            //   [B] Riemann:      p_d = -beta Z* u        ->  F -= (V_i^2+V_j^2) p_d grad W
            Scalar diss = 0.0;
            // [A] Monaghan AV — Monaghan (1992) Annu. Rev. Astron. Astrophys. 30, 543–574
            if ( this->m_artificial_viscosity && !j_issolid )
                {
                Scalar dotdvdx = dot(dv,dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar muij    = meanh*dotdvdx/(rsq+epssqr);
                    Scalar meanrho = Scalar(0.5)*(rhoi+rhoj);
                    diss = mi*mj*(-this->m_avalpha*this->m_cmax*muij+this->m_avbeta*muij*muij)/meanrho;
                    }
                }
            // [B] Riemann dissipation — Zhang, Hu & Adams (2017) J. Comput. Phys. 340, 439–455
            // The dissipative pair pressure is p_d = -beta Z* u^- (Z u is already a
            // pressure; no division by the mean density).
            else if ( m_riemann_dissipation && !j_issolid )
                {
                Scalar dotdvdx = dot(dv, dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar uij   = dotdvdx / (r + eps);
                    Scalar Zi    = rhoi * ci;
                    Scalar Zj    = rhoj * cj;
                    Scalar Zstar = (Zi * Zj) / (Zi + Zj);
                    Scalar pd    = -m_riemann_beta * Zstar * uij;
                    diss = (Vi*Vi + Vj*Vj) * pd;
                    }
                }

            // Add pressure + dissipation force contribution to fluid particle
            h_force.data[i].x -= ( prefactor*temp0 + diss ) * dwdr_r * dx.x;
            h_force.data[i].y -= ( prefactor*temp0 + diss ) * dwdr_r * dx.y;
            h_force.data[i].z -= ( prefactor*temp0 + diss ) * dwdr_r * dx.z;

            // Evaluate viscous interaction forces. Non-Newtonian phases use
            // the per-particle frame-invariant shear rate stored in the
            // energy array by compute_strain_rate() (zero when no NN model is
            // active, in which case computeNNViscosity ignores it anyway).
            {
            bool nn_active = (m_nn_model1 != NEWTONIAN || m_nn_model2 != NEWTONIAN);
            Scalar gdot_i = nn_active ? h_gdot.data[i] : Scalar(0);
            Scalar gdot_j = (nn_active && !j_issolid) ? h_gdot.data[k] : gdot_i;
            NonNewtonianModel nn_model_i = i_isfluid1 ? m_nn_model1 : m_nn_model2;
            Scalar mu_eff_i = computeNNViscosity(mui, gdot_i, nn_model_i,
                i_isfluid1 ? m_nn_K1 : m_nn_K2,
                i_isfluid1 ? m_nn_n1 : m_nn_n2,
                i_isfluid1 ? m_nn_mu0_1 : m_nn_mu0_2,
                i_isfluid1 ? m_nn_muinf_1 : m_nn_muinf_2,
                i_isfluid1 ? m_nn_lambda1 : m_nn_lambda2,
                i_isfluid1 ? m_nn_tauy1 : m_nn_tauy2,
                i_isfluid1 ? m_nn_m1 : m_nn_m2,
                i_isfluid1 ? m_nn_mu_min1 : m_nn_mu_min2);
            Scalar mu_eff_j;
            if (j_issolid)
                mu_eff_j = mu_eff_i;
            else
                {
                NonNewtonianModel nn_model_j = j_isfluid1 ? m_nn_model1 : m_nn_model2;
                mu_eff_j = computeNNViscosity(muj, gdot_j, nn_model_j,
                    j_isfluid1 ? m_nn_K1 : m_nn_K2,
                    j_isfluid1 ? m_nn_n1 : m_nn_n2,
                    j_isfluid1 ? m_nn_mu0_1 : m_nn_mu0_2,
                    j_isfluid1 ? m_nn_muinf_1 : m_nn_muinf_2,
                    j_isfluid1 ? m_nn_lambda1 : m_nn_lambda2,
                    j_isfluid1 ? m_nn_tauy1 : m_nn_tauy2,
                    j_isfluid1 ? m_nn_m1 : m_nn_m2,
                    j_isfluid1 ? m_nn_mu_min1 : m_nn_mu_min2);
                }
            Scalar mu_harm = Scalar(2) * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
            temp0 = mu_harm * (Vi*Vi+Vj*Vj) * dwdr_r;
            }
            h_force.data[i].x  += temp0 * dv.x;
            h_force.data[i].y  += temp0 * dv.y;
            h_force.data[i].z  += temp0 * dv.z;

            // Evaluate rate of change of density if CONTINUITY approach is used
            if ( this->m_density_method == DENSITYCONTINUITY )
                {
                if ( j_issolid )
                    {
                    // Use physical advection velocity rather than fictitious velocity here
                    vj.x = h_velocity.data[k].x;
                    vj.y = h_velocity.data[k].y;
                    vj.z = h_velocity.data[k].z;

                    // Recompute velocity difference
                    dv.x = vi.x - vj.x;
                    dv.y = vi.y - vj.y;
                    dv.z = vi.z - vj.z;

                    //Vj = mj / m_rho0;
                    }

                // Compute density rate of change
                h_ratedpe.data[i].x += rhoi*Vj*dot(dv,dwdr_r*dx);
                //h_ratedpe.data[i].x += mj*dot(dv,dwdr_r*dx);

                // Molteni–Colagrossi density diffusion (fluid–fluid pairs only).
                // Ref: Molteni & Colagrossi (2009) Comput. Phys. Commun. 180, 861–872.
                //
                // Drive term is $(\rho_i/\rho_{0i} - \rho_j/\rho_{0j})$ -- rest-density normalised.
                // The original term $(\rho_i/\rho_j - 1)$ is non-zero at equilibrium when
                // $\rho_{01} \neq \rho_{02}$ (different-phase rest densities), generating unphysical
                // density drift across the interface in stratified-flow setups.
                // The normalised form equals zero at equilibrium for both phases.
                // Corrected sign (Laplacian smoothing of the normalized density):
                //   drho_i/dt += 2 delta h c V_j rho0_i (rho_i/rho0_i - rho_j/rho0_j) (dW/dr)/r
                if ( !j_issolid && this->m_density_diffusion )
                    h_ratedpe.data[i].x += Scalar(2)*m_ddiff*meanh*m_cmax*(mj/rhoj)*rho0i*(rhoi/rho0i-rhoj/rho0j)*dwdr_r;
                }

            } // Closing Neighbor Loop

        // NOTE: pressure is re-evaluated from the per-phase EOS every step in
        // computeForces() (DENSITYCONTINUITY), so no dp/dt is integrated.

        // Add surface force
        h_force.data[i].x  += h_sf.data[i].x;
        h_force.data[i].y  += h_sf.data[i].y;
        h_force.data[i].z  += h_sf.data[i].z;

        } // Closing Fluid Particle Loop

    this->m_max_vel = Scalar(max_vel);
    } // End GPU Array Scope

    // Add volumetric force (gravity)
    this->applyBodyForce(timestep, this->m_fluidgroup);
    // if ( m_compute_solid_forces )
    //     this->applyBodyForce(timestep, m_solidgroup);

    } // end forcecomputation 



/*! Compute forces definition
*/
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::computeForces(uint64_t timestep)
{

    // start by updating the neighborlist
    this->m_nlist->compute(timestep);

    // This is executed once to initialize protected/private variables
    if (!m_params_set)
        {
        this->m_exec_conf->msg->error() << "sph.models.TwoPhaseFlow requires parameters to be set before run()" << std::endl;
        throw std::runtime_error("Error computing TwoPhaseFlow forces");
        }

    // m_solid_removed flag is set to False initially, so this 
    // only executes at timestep 0
    if (!m_solid_removed)
        {
        this->m_nlist->forceUpdate();
        this->m_nlist->compute(timestep);
        mark_solid_particles_toremove(timestep);
        this->m_solid_removed = true;
        }

    // Apply Shepard density renormalization if requested.
    // Density is reset from scratch, so pressure must also be reinitialized from EOS.
    if ( m_shepard_renormalization && timestep % m_shepardfreq == 0 )
        {
        renormalize_density(timestep);
        compute_pressure(timestep);
        m_pressure_initialized = true;
#ifdef ENABLE_MPI
         // Update ghost particle densities and pressures.
        update_ghost_density_pressure(timestep);
#endif
        }

    // Apply density reinitialization from summation if requested (DENSITYCONTINUITY only).
    // Density is reset from scratch, so pressure must also be reinitialized from EOS.
    if ( m_density_reinitialization && timestep % m_densityreinitfreq == 0 )
        {
        compute_ndensity(timestep);
        compute_pressure(timestep);
        m_pressure_initialized = true;
#ifdef ENABLE_MPI
        update_ghost_density_pressure(timestep);
#endif
        }

    if ( m_fickian_shifting )
    {
        compute_particle_concentration_gradient(timestep);
#ifdef ENABLE_MPI
        update_ghost_density_pressure_energy(timestep);
#endif
    }

    if (m_density_method == DENSITYSUMMATION)
    {
        // Density re-evaluated from kernel summation every step; derive pressure from EOS.
        compute_ndensity(timestep);
        compute_pressure(timestep);
    }
    else // DENSITYCONTINUITY
    {
        // Density is time-integrated via the continuity equation. Pressure is
        // re-evaluated from the per-phase EOS every step so it stays exactly
        // consistent with the integrated density.
        compute_pressure(timestep);
    }

#ifdef ENABLE_MPI
    // Update ghost particle densities and pressures.
    update_ghost_density_pressure(timestep);
#endif

    // Compute particle pressures
    compute_noslip(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    update_ghost_density_pressure(timestep);
#endif

    // Compute particle interfacial color gradient
    compute_colorgradients(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    update_ghost_aux123(timestep);
#endif

    // $\delta^+$-SPH particle shifting (Sun et al. 2017).
    // Interface normals in aux3 must be up-to-date before calling.
    // No forced neighbor-list rebuild: shifts are far below the nlist buffer
    // skin, and the standard displacement check picks them up next step. A
    // forced full rebuild every step costs ~7x in throughput at 1e5 particles.
    if ( m_particle_shifting )
        {
        compute_particle_shift(timestep);
        }

    compute_surfaceforce(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    update_ghost_aux4(timestep);
#endif


    // Non-Newtonian rheology: compute the per-particle shear rate. Uses the
    // energy array, which Fickian shifting also claims for |grad C|^2 — the
    // two features are therefore mutually exclusive.
    if ( m_nn_model1 != NEWTONIAN || m_nn_model2 != NEWTONIAN )
        {
        if ( m_fickian_shifting )
            throw std::runtime_error(
                "TwoPhaseFlow: non-Newtonian rheology and Fickian shifting are "
                "mutually exclusive (both store per-particle data in the energy "
                "array). Disable one of the two.");
        compute_strain_rate(timestep);
#ifdef ENABLE_MPI
        // energy piggybacks on the density/pressure ghost exchange
        update_ghost_density_pressure_energy(timestep);
#endif
        }

    // Execute the force computation
    // This includes the computation of the density if
    // DENSITYCONTINUITY method is used
    forcecomputation(timestep);

    if ( m_compute_solid_forces )
        {
        compute_solid_forces(timestep);
        }

} // End computeForces



template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::activateDensityReinitialization(unsigned int densityreinitfreq)
    {
    if (densityreinitfreq <= 0)
        {
        this->m_exec_conf->msg->error() << "sph.models.TwoPhaseFlow: Density reinitialization period has to be a positive real number" << std::endl;
        throw std::runtime_error("Error initializing TwoPhaseFlow.");
        }
    m_density_reinitialization = true;
    m_densityreinitfreq = densityreinitfreq;
    }


/*! Activate \f$\delta^+\f$-SPH particle shifting.
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::activateParticleShifting(Scalar A, Scalar R,
                                                                unsigned int n,
                                                                bool interface_condition)
    {
    if (A <= Scalar(0))
        {
        this->m_exec_conf->msg->error() << "sph.models.TwoPhaseFlow: shift amplitude A must be > 0" << std::endl;
        throw std::runtime_error("Error initializing TwoPhaseFlow particle shifting.");
        }
    m_particle_shifting         = true;
    m_shift_A                   = A;
    m_shift_R                   = R;
    m_shift_n                   = n;
    m_shift_interface_condition = interface_condition;
    }


/*! Compute and apply \f$\delta^+\f$-SPH particle position shifts (Sun et al. 2017).
 *
 * Three-pass algorithm:
 *   Pass 1 -- compute shift vector \f$\delta r_i\f$ for every fluid particle using
 *             the Sun et al. kernel-gradient formula with enhancement factor,
 *             then project out the interface-normal component if requested.
 *   Pass 2 -- (DENSITYCONTINUITY only) apply ALE density remapping correction:
 *             \f$\Delta\rho_i = \rho_i \sum_j V_j (\delta r_i - \delta r_j) \cdot \nabla W_{ij}\f$
 *   Pass 3 — update particle positions and wrap periodic boundaries.
 *
 * ArrayHandle scopes are separated so the same array is never opened twice.
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_particle_shift(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::ParticleShift." << endl;

    const BoxDim& box          = this->m_pdata->getGlobalBox();
    const unsigned int N_local  = this->m_pdata->getN();
    const unsigned int N_total  = N_local + this->m_pdata->getNGhosts();
    const unsigned int fluid_size = this->m_fluidgroup->getNumMembers();
    const Scalar eps      = Scalar(1e-10);
    const Scalar eps_norm = Scalar(1e-6);

    // Shift vectors for all slots; ghost slots remain zero (used as $\delta r_k = 0$ approx).
    std::vector<Scalar3> shift_vec(N_total,
                                   make_scalar3(Scalar(0), Scalar(0), Scalar(0)));

    { // ── scope: read positions / kernel data; readwrite density (ALE) ─────────
    ArrayHandle<Scalar4> h_pos     (this->m_pdata->getPositions(),     access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(),    access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_density (this->m_pdata->getDensities(),     access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar>  h_h       (this->m_pdata->getSlengths(),      access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_fn      (this->m_pdata->getAuxiliaries3(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(m_type_property_map,           access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_n_neigh  (this->m_nlist->getNNeighArray(),        access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist_arr(this->m_nlist->getNListArray(),         access_location::host, access_mode::read);
    ArrayHandle<size_t>       h_head_list(this->m_nlist->getHeadList(),           access_location::host, access_mode::read);

    // ── PASS 1: shift vectors ──────────────────────────────────────────────────
    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp14(this->m_fluidgroup->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for
    for (unsigned int group_idx = 0; group_idx < fluid_size; group_idx++)
        {
        unsigned int i = h_members_omp14.data[group_idx];

        Scalar hi    = m_const_slength ? m_ch : h_h.data[i];
        // W_ref: kernel at approx. initial inter-particle spacing $\Delta p \approx 0.5 h$
        Scalar w_ref = this->m_skernel->wij(hi, Scalar(0.5)*hi);
        if (w_ref < eps) w_ref = eps;

        Scalar3 grad_sum = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
        unsigned int n_neigh = h_n_neigh.data[i];
        size_t       head    = h_head_list.data[i];

        for (unsigned int neigh_idx = 0; neigh_idx < n_neigh; neigh_idx++)
            {
            unsigned int k = h_nlist_arr.data[head + neigh_idx];

            // Solid dummy particles PARTICIPATE in the shift sum: with their
            // Adami-extrapolated densities they close the kernel support at
            // walls, so the shift sees no spurious "free surface" that would
            // otherwise pump wall-adjacent fluid into the solid every step.
            // Guard against marked-removed solids carrying zero density.
            Scalar mk   = h_velocity.data[k].w;
            Scalar rhok = h_density.data[k];
            if (rhok < Scalar(1e-12)) continue;
            Scalar hk   = m_const_slength ? m_ch : h_h.data[k];
            Scalar Vk   = mk / rhok;

            Scalar3 dx;
            dx.x = h_pos.data[i].x - h_pos.data[k].x;
            dx.y = h_pos.data[i].y - h_pos.data[k].y;
            dx.z = h_pos.data[i].z - h_pos.data[k].z;
            dx = box.minImage(dx);

            Scalar rsq = dx.x*dx.x + dx.y*dx.y + dx.z*dx.z;
            if (rsq > m_rcutsq) continue;
            Scalar r = sqrt(rsq);

            Scalar meanh  = Scalar(0.5)*(hi + hk);
            Scalar dwdr   = this->m_skernel->dwijdr(meanh, r);
            Scalar wij_   = this->m_skernel->wij(meanh, r);
            Scalar dwdr_r = (r > Scalar(1e-8)*meanh) ? dwdr/r : Scalar(0);

            // Sun et al. 2017 enhancement factor: [1 + R*(W_ij/W_ref)^n]
            Scalar ratio = wij_ / w_ref;
            Scalar Rpow  = Scalar(1);
            for (unsigned int p = 0; p < m_shift_n; p++) Rpow *= ratio;
            Scalar enhance = Scalar(1) + m_shift_R * Rpow;

            grad_sum.x += enhance * Vk * dwdr_r * dx.x;
            grad_sum.y += enhance * Vk * dwdr_r * dx.y;
            grad_sum.z += enhance * Vk * dwdr_r * dx.z;
            }

        // $\delta r_i = -A\,\mathrm{Ma}_i (2h_i)^2 \sum_j [1 + R(W/W_\mathrm{ref})^n] V_j \nabla W_{ij}$
        // Sun et al. 2017 delta^+ scaling: the per-particle Mach number
        // Ma_i = |v_i|/c makes the shift vanish for quiescent fluid, so the
        // one-sided truncation next to excluded solid neighbors cannot pump
        // static particles into the walls.
        Scalar ci = checkfluid1(h_type_property_map.data, h_pos.data[i].w) ? m_c1 : m_c2;
        Scalar vmagi = sqrt(h_velocity.data[i].x*h_velocity.data[i].x +
                            h_velocity.data[i].y*h_velocity.data[i].y +
                            h_velocity.data[i].z*h_velocity.data[i].z);
        Scalar Mai   = vmagi / ci;
        Scalar coeff = -m_shift_A * Mai * Scalar(4.0) * hi * hi;
        Scalar3 dr;
        dr.x = coeff * grad_sum.x;
        dr.y = coeff * grad_sum.y;
        dr.z = coeff * grad_sum.z;

        // Interface condition: project out normal component at fluid–fluid interface
        // so particles cannot cross between phases (Mokos 2017, Lyu 2021).
        if (m_shift_interface_condition)
            {
            Scalar3 fn_i   = h_fn.data[i];
            Scalar  fn_mag = sqrt(fn_i.x*fn_i.x + fn_i.y*fn_i.y + fn_i.z*fn_i.z);
            if (fn_mag > eps_norm)
                {
                Scalar  inv_mag = Scalar(1) / fn_mag;
                Scalar3 n_hat   = make_scalar3(fn_i.x*inv_mag, fn_i.y*inv_mag, fn_i.z*inv_mag);
                Scalar  dr_n    = dr.x*n_hat.x + dr.y*n_hat.y + dr.z*n_hat.z;
                dr.x -= dr_n * n_hat.x;
                dr.y -= dr_n * n_hat.y;
                dr.z -= dr_n * n_hat.z;
                }
            }

        // NaN-safe magnitude cap at 0.1 h_i: degenerate near-pairs make
        // dW/dr / r blow up, and box.wrap() only unwraps a single periodic
        // image, so an unbounded shift would strand the particle outside the
        // box and corrupt the cell binning.
        Scalar drmag2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
        const Scalar drcap = Scalar(0.1) * hi;
        if (!std::isfinite(drmag2))
            {
            dr = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
            }
        else if (drmag2 > drcap * drcap)
            {
            Scalar rescale = drcap / sqrt(drmag2);
            dr.x *= rescale;
            dr.y *= rescale;
            dr.z *= rescale;
            }

        shift_vec[i] = dr;
        } // end PASS 1

    // ── PASS 2: ALE density correction (DENSITYCONTINUITY only) ───────────────
    // $\Delta\rho_i = \rho_i \sum_j V_j (\delta r_i - \delta r_j) \cdot \nabla W_{ij}$
    // Ghost neighbor j gets $\delta r_j = 0$ (conservative approximation).
    if (m_density_method == DENSITYCONTINUITY)
        {
        // NOTE: deliberately NOT OpenMP-parallelized. This loop writes
        // h_density[i] while also reading neighbor densities h_density[k]
        // for the volume V_k — a cross-iteration read/write overlap.
        for (unsigned int group_idx = 0; group_idx < fluid_size; group_idx++)
            {
            unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

            Scalar hi   = m_const_slength ? m_ch : h_h.data[i];
            Scalar rhoi = h_density.data[i];
            unsigned int n_neigh = h_n_neigh.data[i];
            size_t       head    = h_head_list.data[i];
            Scalar delta_rho = Scalar(0);

            for (unsigned int neigh_idx = 0; neigh_idx < n_neigh; neigh_idx++)
                {
                unsigned int k = h_nlist_arr.data[head + neigh_idx];

                // Solid neighbors participate with delta r_k = 0 (stationary);
                // consistent with their inclusion in the PASS 1 support closure.
                Scalar mk   = h_velocity.data[k].w;
                Scalar rhok = h_density.data[k];
                if (rhok < Scalar(1e-12)) continue;
                Scalar hk   = m_const_slength ? m_ch : h_h.data[k];
                Scalar Vk   = mk / rhok;

                Scalar3 dx;
                dx.x = h_pos.data[i].x - h_pos.data[k].x;
                dx.y = h_pos.data[i].y - h_pos.data[k].y;
                dx.z = h_pos.data[i].z - h_pos.data[k].z;
                dx = box.minImage(dx);

                Scalar rsq = dx.x*dx.x + dx.y*dx.y + dx.z*dx.z;
                if (rsq > m_rcutsq) continue;
                Scalar r = sqrt(rsq);

                Scalar meanh  = Scalar(0.5)*(hi + hk);
                Scalar dwdr   = this->m_skernel->dwijdr(meanh, r);
                Scalar dwdr_r = (r > Scalar(1e-8)*meanh) ? dwdr/r : Scalar(0);

                // $\delta r_i - \delta r_k$; ghost particles (k >= N_local) get $\delta r_k = 0$
                Scalar3 ddr;
                ddr.x = shift_vec[i].x - (k < N_local ? shift_vec[k].x : Scalar(0));
                ddr.y = shift_vec[i].y - (k < N_local ? shift_vec[k].y : Scalar(0));
                ddr.z = shift_vec[i].z - (k < N_local ? shift_vec[k].z : Scalar(0));

                delta_rho += rhoi * Vk * (ddr.x*dwdr_r*dx.x +
                                          ddr.y*dwdr_r*dx.y +
                                          ddr.z*dwdr_r*dx.z);
                }
            h_density.data[i] += delta_rho;
            }
        } // end PASS 2

    } // ── end scope: read handles released ──────────────────────────────────────

    // ── PASS 3: apply position updates and wrap periodic boundaries ───────────
    {
    ArrayHandle<Scalar4> h_pos_rw(this->m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<int3>    h_image (this->m_pdata->getImages(),    access_location::host, access_mode::readwrite);

    // Acquire the group index array once: getMemberIndex() acquires an
    // ArrayHandle per call, which is not thread-safe inside the parallel loop
    ArrayHandle<unsigned int> h_members_omp15(this->m_fluidgroup->getIndexArray(), access_location::host, access_mode::read);
    #pragma omp parallel for
    for (unsigned int group_idx = 0; group_idx < fluid_size; group_idx++)
        {
        unsigned int i = h_members_omp15.data[group_idx];
        h_pos_rw.data[i].x += shift_vec[i].x;
        h_pos_rw.data[i].y += shift_vec[i].y;
        h_pos_rw.data[i].z += shift_vec[i].z;
        box.wrap(h_pos_rw.data[i], h_image.data[i]);
        }
    }
    } // end compute_particle_shift


/*! Returns provided timestep quantities for use in adaptive timestep controllers
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
std::vector<double> TwoPhaseFlow<KT_, SET1_, SET2_>::getProvidedTimestepQuantities(uint64_t timestep)
{
    m_timestep_list[0] = m_rho01;
    m_timestep_list[1] = m_rho02;
    m_timestep_list[2] = m_c1;
    m_timestep_list[3] = m_c2;
    m_timestep_list[4] = m_ch;

    Scalar3 acc = this->getAcceleration(timestep);
    Scalar acc_total = sqrt((acc.x * acc.x) + (acc.y * acc.y) + (acc.z * acc.z));
    m_timestep_list[5] = acc_total;

    m_timestep_list[6] = m_mu1;
    m_timestep_list[7] = m_mu2;

    // Maximum fluid speed of the last force computation (MPI-reduced).
    m_timestep_list[8] = this->getMaxVelocity();

    return m_timestep_list;
}


/*! Compute fluid-induced forces on solid particles
 */
template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_solid_forces(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Compute Solid Forces." << endl;

    const BoxDim& box = this->m_pdata->getGlobalBox();
    const unsigned int group_size = m_solidgroup->getNumMembers();

        { // GPU Array Scope
        ArrayHandle<Scalar4> h_force(this->m_force, access_location::host, access_mode::readwrite);

        ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_density(this->m_pdata->getDensities(), access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_pressure(this->m_pdata->getPressures(), access_location::host, access_mode::read);
        ArrayHandle<Scalar3> h_vf(this->m_pdata->getAuxiliaries1(), access_location::host, access_mode::read);
        ArrayHandle<Scalar>  h_gdot(this->m_pdata->getEnergies(), access_location::host, access_mode::read); // per-particle shear rate
        ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

        ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
        ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

        assert(h_pos.data != NULL);

        unsigned int size;
        size_t myHead;
        Scalar temp0 = 0;

        // For each solid particle, accumulate the exact reaction of the pair
        // forces applied to its fluid neighbors in forcecomputation(). The pair
        // expressions are identical to the fluid loop (symmetric in i<->j with
        // dx and dv flipping sign), so Newton's third law holds without any
        // extra sign flip or mass-ratio scaling.
        // Acquire the group index array once: getMemberIndex() acquires an
        // ArrayHandle per call, which is not thread-safe inside the parallel loop
        ArrayHandle<unsigned int> h_members_omp16(m_solidgroup->getIndexArray(), access_location::host, access_mode::read);
        #pragma omp parallel for private(size, myHead) firstprivate(temp0)
        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            unsigned int i = h_members_omp16.data[group_idx];

            Scalar3 pi;
            pi.x = h_pos.data[i].x;
            pi.y = h_pos.data[i].y;
            pi.z = h_pos.data[i].z;

            // Fictitious (Adami) velocity: the fluid loop computed its viscous
            // pair force with dv = v_f - v~_s, so the reaction uses the same.
            Scalar3 vi;
            vi.x = h_vf.data[i].x;
            vi.y = h_vf.data[i].y;
            vi.z = h_vf.data[i].z;
            Scalar mi = h_velocity.data[i].w;

            Scalar Pi   = h_pressure.data[i];
            Scalar rhoi = h_density.data[i];
            Scalar Vi   = mi / rhoi;

            myHead = h_head_list.data[i];
            size = (unsigned int)h_n_neigh.data[i];

            for (unsigned int j = 0; j < size; j++)
                {
                unsigned int k = h_nlist.data[myHead + j];
                assert(k < this->m_pdata->getN() + this->m_pdata->getNGhosts());

                bool issolid = checksolid(h_type_property_map.data, h_pos.data[k].w);
                if ( issolid ) { continue; }

                bool j_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[k].w);

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

                if ( this->m_const_slength && rsq > this->m_rcutsq )
                    continue;

                Scalar3 vj;
                vj.x = h_velocity.data[k].x;
                vj.y = h_velocity.data[k].y;
                vj.z = h_velocity.data[k].z;
                Scalar mj   = h_velocity.data[k].w;
                Scalar rhoj = h_density.data[k];
                Scalar Vj   = mj / rhoj;
                Scalar Pj   = h_pressure.data[k];

                Scalar3 dv;
                dv.x = vi.x - vj.x;
                dv.y = vi.y - vj.y;
                dv.z = vi.z - vj.z;

                Scalar r = sqrt(rsq);

                Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);

                Scalar dwdr   = this->m_skernel->dwijdr(meanh, r);
                Scalar dwdr_r = (r > Scalar(1e-8)*meanh) ? dwdr/r : Scalar(0);

                // Same expression as the fluid loop: symmetric in i<->j, so this
                // yields the exact reaction force through the flipped dx.
                if ( m_density_method == DENSITYSUMMATION )
                    temp0 = (Vi*Vi+Vj*Vj)*((rhoj*Pi+rhoi*Pj)/(rhoi+rhoj));
                else
                    temp0 = mi*mj*(Pi+Pj)/(rhoi*rhoj);

                h_force.data[i].x -= temp0 * dwdr_r * dx.x;
                h_force.data[i].y -= temp0 * dwdr_r * dx.y;
                h_force.data[i].z -= temp0 * dwdr_r * dx.z;

                // Use viscosity of the fluid neighbor (with NN rheology,
                // evaluated at the fluid particle's per-particle shear rate)
                Scalar muj_base = j_isfluid1 ? this->m_mu1 : this->m_mu2;
                {
                bool nn_active = (m_nn_model1 != NEWTONIAN || m_nn_model2 != NEWTONIAN);
                Scalar gdot_j = nn_active ? h_gdot.data[k] : Scalar(0);
                NonNewtonianModel nn_model_j = j_isfluid1 ? m_nn_model1 : m_nn_model2;
                Scalar mu_eff_j = computeNNViscosity(muj_base, gdot_j, nn_model_j,
                    j_isfluid1 ? m_nn_K1 : m_nn_K2,
                    j_isfluid1 ? m_nn_n1 : m_nn_n2,
                    j_isfluid1 ? m_nn_mu0_1 : m_nn_mu0_2,
                    j_isfluid1 ? m_nn_muinf_1 : m_nn_muinf_2,
                    j_isfluid1 ? m_nn_lambda1 : m_nn_lambda2,
                    j_isfluid1 ? m_nn_tauy1 : m_nn_tauy2,
                    j_isfluid1 ? m_nn_m1 : m_nn_m2,
                    j_isfluid1 ? m_nn_mu_min1 : m_nn_mu_min2);
                temp0 = mu_eff_j * (Vi*Vi+Vj*Vj) * dwdr_r;
                }
                // Viscous reaction (dv = v~_s - v_f, the negative of the fluid
                // loop's dv, so this is -F_fluid).
                h_force.data[i].x += temp0 * dv.x;
                h_force.data[i].y += temp0 * dv.y;
                h_force.data[i].z += temp0 * dv.z;

                } // End neighbor loop

            } // End solid particle loop
        } // End GPU Array Scope
    }


namespace detail
{

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void export_TwoPhaseFlow(pybind11::module& m, std::string name)
{
    pybind11::class_<TwoPhaseFlow<KT_, SET1_, SET2_>, SPHBaseClass<KT_, SET1_>, std::shared_ptr<TwoPhaseFlow<KT_, SET1_, SET2_>>>(m, name.c_str())
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
        .def("setParams", &TwoPhaseFlow<KT_, SET1_, SET2_>::setParams)
        .def("setHysteresis", &TwoPhaseFlow<KT_, SET1_, SET2_>::setHysteresis)
        .def("getDensityMethod", &TwoPhaseFlow<KT_, SET1_, SET2_>::getDensityMethod)
        .def("setDensityMethod", &TwoPhaseFlow<KT_, SET1_, SET2_>::setDensityMethod)
        .def("getViscosityMethod", &TwoPhaseFlow<KT_, SET1_, SET2_>::getViscosityMethod)
        .def("setViscosityMethod", &TwoPhaseFlow<KT_, SET1_, SET2_>::setViscosityMethod)
        .def("getColorGradientMethod", &TwoPhaseFlow<KT_, SET1_, SET2_>::getColorGradientMethod)
        .def("setColorGradientMethod", &TwoPhaseFlow<KT_, SET1_, SET2_>::setColorGradientMethod)
        .def("setConstSmoothingLength", &TwoPhaseFlow<KT_, SET1_, SET2_>::setConstSmoothingLength)
        .def("computeSolidForces", &TwoPhaseFlow<KT_, SET1_, SET2_>::computeSolidForces)
        .def("activateArtificialViscosity", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateArtificialViscosity)
        .def("activateConsistentInterfacePressure", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateConsistentInterfacePressure)
        .def("deactivateConsistentInterfacePressure", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateConsistentInterfacePressure)
        .def("activateRiemannDissipation", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateRiemannDissipation,
             pybind11::arg("beta") = Scalar(1.0))
        .def("deactivateRiemannDissipation", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateRiemannDissipation)
        .def("activateDensityDiffusion", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateShepardRenormalization)
        .def("activateDensityReinitialization", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateDensityReinitialization)
        .def("deactivateDensityReinitialization", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateDensityReinitialization)
        .def("activateFickianShifting", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateFickianShifting)
        .def("deactivateFickianShifting", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateFickianShifting)
        .def("activateParticleShifting", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateParticleShifting,
             pybind11::arg("A")                   = Scalar(0.2),
             pybind11::arg("R")                   = Scalar(0.2),
             pybind11::arg("n")                   = 4,
             pybind11::arg("interface_condition") = true)
        .def("deactivateParticleShifting", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateParticleShifting)
        .def("getProvidedTimestepQuantities", &TwoPhaseFlow<KT_, SET1_, SET2_>::getProvidedTimestepQuantities)
        .def("activatePowerLaw1", &TwoPhaseFlow<KT_, SET1_, SET2_>::activatePowerLaw1,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("mu_min") = Scalar(0))
        .def("activateCarreau1", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateCarreau1)
        .def("activateBingham1", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateBingham1,
             pybind11::arg("mu_p"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("activateHerschelBulkley1", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateHerschelBulkley1,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("deactivateNonNewtonian1", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateNonNewtonian1)
        .def("activatePowerLaw2", &TwoPhaseFlow<KT_, SET1_, SET2_>::activatePowerLaw2,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("mu_min") = Scalar(0))
        .def("activateCarreau2", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateCarreau2)
        .def("activateBingham2", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateBingham2,
             pybind11::arg("mu_p"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("activateHerschelBulkley2", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateHerschelBulkley2,
             pybind11::arg("K"), pybind11::arg("n"), pybind11::arg("tauy"), pybind11::arg("m_reg"),
             pybind11::arg("mu_min") = Scalar(0))
        .def("deactivateNonNewtonian2", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateNonNewtonian2)
        .def("setAcceleration", &SPHBaseClass<KT_, SET1_>::setAcceleration)
        .def("setRCut", &TwoPhaseFlow<KT_, SET1_, SET2_>::setRCutPython)
        ;
}

} // end namespace detail

//! Explicit template instantiations
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc2, linear, linear>;
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc2, linear, tait>;
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc2, tait, linear>;
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc2, tait, tait>;
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc4, linear, linear>;
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc4, linear, tait>;
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc4, tait, linear>;
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc4, tait, tait>;
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc6, linear, linear>;
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc6, linear, tait>;
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc6, tait, linear>;
template class PYBIND11_EXPORT TwoPhaseFlow<wendlandc6, tait, tait>;
template class PYBIND11_EXPORT TwoPhaseFlow<quintic, linear, linear>;
template class PYBIND11_EXPORT TwoPhaseFlow<quintic, linear, tait>;
template class PYBIND11_EXPORT TwoPhaseFlow<quintic, tait, linear>;
template class PYBIND11_EXPORT TwoPhaseFlow<quintic, tait, tait>;
template class PYBIND11_EXPORT TwoPhaseFlow<cubicspline, linear, linear>;
template class PYBIND11_EXPORT TwoPhaseFlow<cubicspline, linear, tait>;
template class PYBIND11_EXPORT TwoPhaseFlow<cubicspline, tait, linear>;
template class PYBIND11_EXPORT TwoPhaseFlow<cubicspline, tait, tait>;


namespace detail
{

    template void export_TwoPhaseFlow<wendlandc2, linear, linear>(pybind11::module& m, std::string name = "TwoPF_WC2_LL");
    template void export_TwoPhaseFlow<wendlandc2, linear, tait>(pybind11::module& m, std::string name = "TwoPF_WC2_LT");
    template void export_TwoPhaseFlow<wendlandc2, tait, linear>(pybind11::module& m, std::string name = "TwoPF_WC2_TL");
    template void export_TwoPhaseFlow<wendlandc2, tait, tait>(pybind11::module& m, std::string name = "TwoPF_WC2_TT");
    
    template void export_TwoPhaseFlow<wendlandc4, linear, linear>(pybind11::module& m, std::string name = "TwoPF_WC4_LL");
    template void export_TwoPhaseFlow<wendlandc4, linear, tait>(pybind11::module& m, std::string name = "TwoPF_WC4_LT");
    template void export_TwoPhaseFlow<wendlandc4, tait, linear>(pybind11::module& m, std::string name = "TwoPF_WC4_TL");
    template void export_TwoPhaseFlow<wendlandc4, tait, tait>(pybind11::module& m, std::string name = "TwoPF_WC4_TT");
    
    template void export_TwoPhaseFlow<wendlandc6, linear, linear>(pybind11::module& m, std::string name = "TwoPF_WC6_LL");
    template void export_TwoPhaseFlow<wendlandc6, linear, tait>(pybind11::module& m, std::string name = "TwoPF_WC6_LT");
    template void export_TwoPhaseFlow<wendlandc6, tait, linear>(pybind11::module& m, std::string name = "TwoPF_WC6_TL");
    template void export_TwoPhaseFlow<wendlandc6, tait, tait>(pybind11::module& m, std::string name = "TwoPF_WC6_TT");
    
    template void export_TwoPhaseFlow<quintic, linear, linear>(pybind11::module& m, std::string name = "TwoPF_Q_LL");
    template void export_TwoPhaseFlow<quintic, linear, tait>(pybind11::module& m, std::string name = "TwoPF_Q_LT");
    template void export_TwoPhaseFlow<quintic, tait, linear>(pybind11::module& m, std::string name = "TwoPF_Q_TL");
    template void export_TwoPhaseFlow<quintic, tait, tait>(pybind11::module& m, std::string name = "TwoPF_Q_TT");
    
    template void export_TwoPhaseFlow<cubicspline, linear, linear>(pybind11::module& m, std::string name = "TwoPF_CS_LL");
    template void export_TwoPhaseFlow<cubicspline, linear, tait>(pybind11::module& m, std::string name = "TwoPF_CS_LT");
    template void export_TwoPhaseFlow<cubicspline, tait, linear>(pybind11::module& m, std::string name = "TwoPF_CS_TL");
    template void export_TwoPhaseFlow<cubicspline, tait, tait>(pybind11::module& m, std::string name = "TwoPF_CS_TT");

}  // end namespace detail
} // end namespace sph
} // end namespace hoomd
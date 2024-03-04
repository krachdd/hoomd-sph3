/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "TwoPhaseFlow.h"

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
        m_density_diffusion = false;
        m_shepard_renormalization = false;
        m_fickian_shifting = false;
        m_ch = Scalar(0.0);
        m_rcut = Scalar(0.0);
        m_rcutsq = Scalar(0.0);
        m_avalpha = Scalar(0.0);
        m_avbeta = Scalar(0.0);
        m_ddiff = Scalar(0.0);
        m_shepardfreq = 1;

        m_solid_removed = false;

        // Sanity checks
        assert(this->m_pdata);
        assert(this->m_nlist);
        assert(this->m_skernel);
        assert(this->m_eos1);
        assert(this->m_eos2);

        // If fluid phase speed of sounds are different, backpressures will be different
        // Apply larger backbressure to entire system
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
                h_type_property_map.data[m_fluidtypes[i]] |= SolidFluidTypeBit::FLUID1;
            }
            for (unsigned int i = 0; i < m_fluidtypes2.size(); i++) {
                h_type_property_map.data[m_fluidtypes[i]] |= SolidFluidTypeBit::FLUID2;
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

        m_r_cut_nlist = std::make_shared<GlobalArray<Scalar>>(m_typpair_idx.getNumElements(), this->m_exec_conf);
        this->m_nlist->addRCutMatrix(m_r_cut_nlist);

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

    // Compute solid-fluid interfacial tension using Young's equation
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
        // store r_cut**2 for use internally
        // ArrayHandle<Scalar> h_rcutsq(m_rcutsq, access_location::host, access_mode::readwrite);
        // h_rcutsq.data[m_typpair_idx(typ1, typ2)] = rcut * rcut;
        // h_rcutsq.data[m_typpair_idx(typ2, typ1)] = rcut * rcut;

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

// Grab handles for particle and neighbor data
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    unsigned int size;
    size_t myHead;

    // For all solid particles
    unsigned int group_size = this->m_solidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_solidgroup->getMemberIndex(group_idx);

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
    memset((void*)h_pressure.data,0,sizeof(Scalar)*this->m_force.getNumElements());

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    unsigned int size;
    size_t myHead;

    // Precompute self-density for homogeneous smoothing lengths
    // Scalar w0 = this->m_skernel->w0(m_ch);

    // Particle loop to compute the particle concentration
    // For each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);
        
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
            // Write \sum_j mj/rhoj * W_ij to h_pressure
            h_pressure.data[i] += (mj/rhoj)*(this->m_const_slength ? this->m_skernel->wij(m_ch,r) : this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r));

        } // End neighbour loop

    } // End fluid group loop

    Scalar3 gradCi;
    Scalar  temp0; 
    // Second loop to compute the actual gradient, stored in h_energy
    group_size = this->m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);
        
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

            // Mean smoothing length and denominator modifier
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar eps    = Scalar(0.1)*meanh;

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);
            
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

    // Precompute self-density for homogeneous smoothing lengths
    Scalar w0 = this->m_skernel->w0(m_ch);

    // Particle loop
    // For each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);
        
        // set temp variable to zero 
        ni = w0;

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
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int i = this->m_fluidgroup1->getMemberIndex(group_idx);
        // Evaluate pressure
        h_pressure.data[i] = this->m_eos1->Pressure(h_density.data[i]);
    
    } // End fluid group 1 loop

    // For each fluid particle of fluidgroup2 
    group_size = this->m_fluidgroup2->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
    {
        // Read particle index
        unsigned int i = this->m_fluidgroup2->getMemberIndex(group_idx);
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
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_solidgroup->getMemberIndex(group_idx);

        // Access the particle's position, velocity, mass and type
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;

        // Read acceleration of solid particle i if content is not NaN
        Scalar3 accel_i = make_scalar3(0,0,0);
        if ( h_accel.data[i].x != h_accel.data[i].x ||
             h_accel.data[i].y != h_accel.data[i].y ||
             h_accel.data[i].z != h_accel.data[i].z )
            {
            }
        else
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
        // if (fluidneighbors > 0 && h_density.data[i] > 0 )
        if (fluidneighbors > 0 && wij_c0 > 0 )
            {
            //  Compute zeroth order normalization constant
            // Scalar norm_constant = 1./h_density.data[i];
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
    // ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
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
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

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
    // ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Zero data before calculation
    memset((void*)h_sn.data,0,sizeof(Scalar3)*this->m_pdata->getAuxiliaries2().getNumElements());
    memset((void*)h_fn.data,0,sizeof(Scalar3)*this->m_pdata->getAuxiliaries3().getNumElements());

    // Particle loop
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

            // Mean smoothing length and denominator modifier
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar eps    = Scalar(0.1)*meanh;

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);


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
    ArrayHandle<Scalar>  h_energy(this->m_pdata->getDensities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    // ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);


    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Zero data before calculation
    memset((void*)h_sf.data,0,sizeof(Scalar3)*this->m_pdata->getAuxiliaries4().getNumElements());

    // for each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

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
            bool j_issolid  = checksolid(h_type_property_map.data, h_pos.data[k].w);
            bool j_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[k].w);
            if ( !(j_issolid) && i_isfluid1 != j_isfluid1 )
                {
                    nearfluidinterface = true;
                    break;
                }
            }
        if ( !nearfluidinterface )
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

        // Fluid phase 1 - Solid interface
        if ( i_isfluid1 && this->m_sigma01 > 0.0 && normsni > 0.0 )
        {
            temp0 = this->m_sigma01/normsni;
            istress[0] += temp0 * ( temp1 - sni.x * sni.x); // xx
            istress[1] += temp0 * ( temp1 - sni.y * sni.y); // yy
            istress[2] += temp0 * ( temp1 - sni.z * sni.z); // zz
            istress[3] -= temp0 * ( sni.x * sni.y);         // xy yx
            istress[4] -= temp0 * ( sni.x * sni.z);         // xz zx
            istress[5] -= temp0 * ( sni.y * sni.z);         // yz zy
        }

        // Fluid phase 2 - Solid interface
        if ( i_isfluid2 && this->m_sigma02 > 0.0 && normsni > 0.0 )
        {
            temp0 = this->m_sigma02/normsni;
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
            bool j_isfluid2 = checkfluid1(h_type_property_map.data, h_pos.data[k].w);

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

            // Mean smoothing length and denominator modifier
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar eps    = Scalar(0.1)*meanh;

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);

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
                istress[0] += temp0 * ( temp1 - fnj.x * fnj.x); // xx
                istress[1] += temp0 * ( temp1 - fnj.y * fnj.y); // yy
                istress[2] += temp0 * ( temp1 - fnj.z * fnj.z); // zz
                istress[3] -= temp0 * ( fnj.x * fnj.y);         // xy yx
                istress[4] -= temp0 * ( fnj.x * fnj.z);         // xz zx
                istress[5] -= temp0 * ( fnj.y * fnj.z);         // yz zy
            }

            if ( !m_fickian_shifting )
            {
                temp1 = 1./3. * normsnj * normsnj;
            }

            // Fluid phase 1 - Solid interface
            if ( j_isfluid1 && this->m_sigma01 > 0.0 && normsnj > 0.0 )
            {
                temp0 = this->m_sigma01/normsnj;
                istress[0] += temp0 * ( temp1 - snj.x * snj.x); // xx
                istress[1] += temp0 * ( temp1 - snj.y * snj.y); // yy
                istress[2] += temp0 * ( temp1 - snj.z * snj.z); // zz
                istress[3] -= temp0 * ( snj.x * snj.y);         // xy yx
                istress[4] -= temp0 * ( snj.x * snj.z);         // xz zx
                istress[5] -= temp0 * ( snj.y * snj.z);         // yz zy
            }

            // Fluid phase 2 - Solid interface
            if ( j_isfluid2 && this->m_sigma02 > 0.0 && normsnj > 0.0 )
            {
                temp0 = this->m_sigma02/normsnj;
                istress[0] += temp0 * ( temp1 - snj.x * snj.x); // xx
                istress[1] += temp0 * ( temp1 - snj.y * snj.y); // yy
                istress[2] += temp0 * ( temp1 - snj.z * snj.z); // zz
                istress[3] -= temp0 * ( snj.x * snj.y);         // xy yx
                istress[4] -= temp0 * ( snj.x * snj.z);         // xz zx
                istress[5] -= temp0 * ( snj.y * snj.z);         // yz zy
            }

            bool oldm = true;


            if ( oldm )
            {
                // Add contribution to surface force
                h_sf.data[i].x += dwdr_r*dx.x*(Vi*Vi*istress[0]+Vj*Vj*jstress[0])+
                                  dwdr_r*dx.y*(Vi*Vi*istress[3]+Vj*Vj*jstress[3])+
                                  dwdr_r*dx.z*(Vi*Vi*istress[4]+Vj*Vj*jstress[4]);
                h_sf.data[i].y += dwdr_r*dx.x*(Vi*Vi*istress[3]+Vj*Vj*jstress[3])+
                                  dwdr_r*dx.y*(Vi*Vi*istress[1]+Vj*Vj*jstress[1])+
                                  dwdr_r*dx.z*(Vi*Vi*istress[5]+Vj*Vj*jstress[5]);
                h_sf.data[i].z += dwdr_r*dx.x*(Vi*Vi*istress[4]+Vj*Vj*jstress[4])+
                                  dwdr_r*dx.y*(Vi*Vi*istress[5]+Vj*Vj*jstress[5])+
                                  dwdr_r*dx.z*(Vi*Vi*istress[2]+Vj*Vj*jstress[2]);
            }

            else
            {
                // Add contribution to surface force
                h_sf.data[i].x += dwdr_r * dx.x * ( istress[0] + jstress[0] )/(rhoi * rhoj)+
                                  dwdr_r * dx.y * ( istress[3] + jstress[3] )/(rhoi * rhoj)+
                                  dwdr_r * dx.z * ( istress[4] + jstress[4] )/(rhoi * rhoj);

                h_sf.data[i].y += dwdr_r * dx.x * ( istress[3] + jstress[3] )/(rhoi * rhoj)+
                                  dwdr_r * dx.y * ( istress[1] + jstress[1] )/(rhoi * rhoj)+
                                  dwdr_r * dx.z * ( istress[5] + jstress[5] )/(rhoi * rhoj);
                
                h_sf.data[i].z += dwdr_r * dx.x * ( istress[4] + jstress[4] )/(rhoi * rhoj)+
                                  dwdr_r * dx.y * ( istress[5] + jstress[5] )/(rhoi * rhoj)+
                                  dwdr_r * dx.z * ( istress[2] + jstress[2] )/(rhoi * rhoj);
            }




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
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_sf(this->m_pdata->getAuxiliaries4(), access_location::host,access_mode::read);

    // access the neighbor list
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    // ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
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

    // maximum velocity variable for adaptive timestep
    double max_vel = 0.0;

    // for each fluid particle
    unsigned int group_size = m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = m_fluidgroup->getMemberIndex(group_idx);

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
        if (i == 0) { max_vel = vi_total; }
        else if (vi_total > max_vel) { max_vel = vi_total; }

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
            Scalar dwdr_r = dwdr/(r+eps);

            // Evaluate inter-particle pressure force
            // Transport formulation proposed by Adami 2013
            Scalar prefactor = 0.0; 
            if ( this->m_densitymethod == DENSITYSUMMATION ) 
            {
                temp0 = (rhoj*Pi+rhoi*Pj)/(rhoi+rhoj);
                prefactor = Vi*Vi + Vj*Vj;
            }
            else if ( this->m_densitymethod == DENSITYCONTINUITY )
            {
                temp0 = (Pi+Pj)/(rhoi*rhoj);
                prefactor = mi * mj;
            }

            Scalar avc = 0.0;
            // Optionally add artificial viscosity
            // Monaghan (1983) J. Comput. Phys. 52 (2) 374–389
            if ( this->m_artificial_viscosity && !j_issolid )
                {
                Scalar dotdvdx = dot(dv,dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar muij    = meanh*dotdvdx/(rsq+epssqr);
                    Scalar meanrho = Scalar(0.5)*(rhoi+rhoj);
                    avc = (-this->m_avalpha*this->m_cmax*muij+this->m_avbeta*muij*muij)/meanrho;
                    }
                }

            // Add contribution to fluid particle
            h_force.data[i].x -= prefactor * ( temp0 + avc )* dwdr_r * dx.x;
            h_force.data[i].y -= prefactor * ( temp0 + avc )* dwdr_r * dx.y;
            h_force.data[i].z -= prefactor * ( temp0 + avc )* dwdr_r * dx.z;

            // Evaluate viscous interaction forces
            temp0 = ((Scalar(2)*(mui*muj))/(mui+muj)) * (Vi*Vi+Vj*Vj) * dwdr_r;
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

                // Add density diffusion if requested
                // Molteni and Colagrossi, Computer Physics Communications 180 (2009) 861–872
                if ( !j_issolid && this->m_density_diffusion )
                    h_ratedpe.data[i].x -= (Scalar(2)*m_ddiff*meanh*m_cmax*mj*(rhoi/rhoj-Scalar(1))*dot(dx,dwdr_r*dx))/(rsq+epssqr);
                }

            } // Closing Neighbor Loop

        // Add surface force
        h_force.data[i].x  += h_sf.data[i].x;
        h_force.data[i].y  += h_sf.data[i].y;
        h_force.data[i].z  += h_sf.data[i].z;

        } // Closing Fluid Particle Loop

    m_timestep_list[5] = max_vel;
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

    // Apply density renormalization if requested
    if ( m_shepard_renormalization && timestep % m_shepardfreq == 0 )
        {
        renormalize_density(timestep);
#ifdef ENABLE_MPI
         // Update ghost particle densities and pressures.
        update_ghost_density(timestep);
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
        compute_ndensity(timestep);
        // compute_particlenumberdensity(timestep);
    }

    // Compute fluid pressure based on m_eos;
    // Only working on the fluidgroup
    compute_pressure(timestep);

#ifdef ENABLE_MPI
    // Update ghost particle densities and pressures.
    update_ghost_density_pressure_energy(timestep);
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

    compute_surfaceforce(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    update_ghost_aux4(timestep);
#endif

    // Execute the force computation
    // This includes the computation of the density if 
    // DENSITYCONTINUITY method is used
    forcecomputation(timestep);

} // End computeForces



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
        .def("activateDensityDiffusion", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateShepardRenormalization)
        .def("activateFickianShifting", &TwoPhaseFlow<KT_, SET1_, SET2_>::activateFickianShifting)
        .def("deactivateFickianShifting", &TwoPhaseFlow<KT_, SET1_, SET2_>::deactivateFickianShifting)
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
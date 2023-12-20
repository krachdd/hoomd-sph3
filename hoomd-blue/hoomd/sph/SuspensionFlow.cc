/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SuspensionFlow.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

using namespace std;

#include <cmath>

namespace hoomd 
{
namespace sph
{
/*! Constructor
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
SuspensionFlow<KT_, SET_>::SuspensionFlow(std::shared_ptr<SystemDefinition> sysdef,
                                 std::shared_ptr<SmoothingKernel<KT_> > skernel,
                                 std::shared_ptr<StateEquation<SET_> > equationofstate,
                                 std::shared_ptr<nsearch::NeighborList> nlist,
                                 std::shared_ptr<ParticleGroup> fluidgroup,
                                 std::shared_ptr<ParticleGroup> solidgroup,
                                 std::shared_ptr<ParticleGroup> suspendedgroup,                                 
                                 DensityMethod mdensitymethod,
                                 ViscosityMethod mviscositymethod)
    : SPHBaseClass<KT_, SET_>(sysdef,skernel,equationofstate,nlist), m_fluidgroup(fluidgroup), m_solidgroup(solidgroup), m_suspendedgroup(suspendedgroup), m_typpair_idx(this->m_pdata->getNTypes())
      {
        this->m_exec_conf->msg->notice(5) << "Constructing SuspensionFlow" << std::endl;

        // Set private attributes to default values
        m_const_slength = false;
        m_params_set = false;
        m_compute_solid_forces = false;
        m_artificial_viscosity = false;
        m_density_diffusion = false;
        m_shepard_renormalization = false;
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
        assert(this->m_eos);

        // Contruct type vectors
        this->constructTypeVectors(fluidgroup,&m_fluidtypes);
        this->constructTypeVectors(solidgroup,&m_solidtypes);
        this->constructTypeVectors(suspendedgroup,&m_suspendedtypes);


        // all particle groups are based on the same particle data
        unsigned int num_types = this->m_sysdef->getParticleData()->getNTypes();

        m_type_property_map = GPUArray<unsigned int>(num_types, this->m_exec_conf);
        {
            ArrayHandle<unsigned int> h_type_property_map(m_type_property_map, access_location::host, access_mode::overwrite);
            fill_n(h_type_property_map.data, num_types, SolidFluidTypeBit::NONE);
            // no need to parallelize this as there should only be a few particle types
            for (unsigned int i = 0; i < m_fluidtypes.size(); i++) {
                h_type_property_map.data[m_fluidtypes[i]] |= SolidFluidTypeBit::FLUID;
            }
            for (unsigned int i = 0; i < m_solidtypes.size(); i++) {
                h_type_property_map.data[m_solidtypes[i]] |= SolidFluidTypeBit::SOLID;
            }
            for (unsigned int i = 0; i < m_suspendedtypes.size(); i++) {
                h_type_property_map.data[m_suspendedtypes[i]] |= SolidFluidTypeBit::SUSPENDED;
            }
        }

        this->m_exec_conf->msg->notice(2) << "Computing SuspensionFlow::number of solids in m_suspendedtypes " << m_suspendedtypes.size() << std::endl;

        // Init vectors for suspension model
        Scalar3 zeros3 = make_scalar3(0.0,0.0,0.0);
        Scalar4 zeros4 = make_scalar4(0.0,0.0,0.0,0.0);
        for (unsigned int i = 0; i < m_suspendedtypes.size(); ++i)
            {
            m_centerofmasses.push_back(zeros4);
            m_repulsiveforces.push_back(zeros3);
            m_velocities.push_back(zeros3);
            m_radii.push_back(0.0);
            }

        // Set simulations methods
        m_density_method = mdensitymethod;
        m_viscosity_method = mviscositymethod;

        // Get necessary variables from kernel and EOS classes
        m_rho0  = equationofstate->getRestDensity();
        m_c     = equationofstate->getSpeedOfSound();
        m_kappa = skernel->getKernelKappa();

        m_r_cut_nlist = std::make_shared<GlobalArray<Scalar>>(m_typpair_idx.getNumElements(), this->m_exec_conf);
        this->m_nlist->addRCutMatrix(m_r_cut_nlist);

      }

/*! Destructor
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
SuspensionFlow<KT_, SET_>::~SuspensionFlow()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying SuspensionFlow" << std::endl;
    }


/*! Returns provided timestep Quantities to Compute
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
std::vector<double> SuspensionFlow<KT_, SET_>::getProvidedTimestepQuantities(uint64_t timestep)
{
    m_timestep_list[0] = m_rho0;
    m_timestep_list[1] = m_c;
    m_timestep_list[2] = m_ch;

    Scalar3 acc =  this->getAcceleration(timestep);

    Scalar acc_total = sqrt((acc.x * acc.x) + (acc.y * acc.y) + (acc.z * acc.z));

    m_timestep_list[3] = acc_total;

    m_timestep_list[4] = m_mu;


    return m_timestep_list;
}

template<SmoothingKernelType KT_,StateEquationType SET_>
void  SuspensionFlow<KT_, SET_>::activateShepardRenormalization(unsigned int shepardfreq){
            if (shepardfreq <= 0)
                {
                this->m_exec_conf->msg->error() << "sph.models.SuspensionFlow: Shepard density reinitialization period has to be a positive real number" << std::endl;
                throw std::runtime_error("Error initializing SuspensionFlow.");
            }
            m_shepard_renormalization = true;
            m_shepardfreq = shepardfreq;
            }


/*! \post Set model parameters
 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::setParams(Scalar mu, Scalar rho0_S, Scalar f0)
    {
    m_mu   = mu;
    m_rhoS = rho0_S;
    m_f0 = f0;
    if (m_mu <= 0)
         {
         this->m_exec_conf->msg->error() << "sph.models.SuspensionFlow: Dynamic viscosity has to be a positive real number" << std::endl;
         throw std::runtime_error("Error initializing SuspensionFlow.");
         }
    //m_rhoS = rho0_S;
    if (m_rhoS <= 0)
        {
        this->m_exec_conf->msg->error() << "sph.models.SuspensionFlow: Solid density has to be a positive real number" << std::endl;
        throw std::runtime_error("Error initializing SuspensionFlow.");
        }
    //m_f0 = f0;
    //if (m_f0 <= 0)
    if (m_f0 < 0)
        {
        this->m_exec_conf->msg->error() << "sph.models.SuspensionFlow: Contact force magnitude has to be a positive real number" << std::endl;
        throw std::runtime_error("Error initializing SuspensionFlow.");
        }
    m_params_set = true;
    }

/*! Communicate (update) ghost particle fields
 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::update_ghost_density_pressure(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Update Ghost density, pressure" << std::endl;

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

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::update_ghost_density(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Update Ghost density" << std::endl;

#ifdef ENABLE_MPI
    if (this->m_comm)
        {
        CommFlags flags(0);
        flags[comm_flag::tag] = 0;
        flags[comm_flag::position] = 0;
        flags[comm_flag::velocity] = 0;
        // flags[comm_flag::dpe] = 1;
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
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::update_ghost_aux1(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Update Ghost aux1" << std::endl;

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
        flags[comm_flag::auxiliary1] = 1; // ficticios velocity
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
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::update_ghost_aux2(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Update Ghost aux2" << std::endl;

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
        flags[comm_flag::auxiliary1] = 0; // ficticios velocity
        flags[comm_flag::auxiliary2] = 1;
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
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::update_ghost_aux34(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Update Ghost aux34" << std::endl;

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
        flags[comm_flag::auxiliary1] = 0; // ficticios velocity
        flags[comm_flag::auxiliary2] = 0;
        flags[comm_flag::auxiliary3] = 1;
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

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::validateTypes(unsigned int typ1,
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
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::setRcut(unsigned int typ1, unsigned int typ2, Scalar rcut)
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

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::setRCutPython(pybind11::tuple types, Scalar r_cut)
    {
    auto typ1 = this->m_pdata->getTypeByName(types[0].cast<std::string>());
    auto typ2 = this->m_pdata->getTypeByName(types[1].cast<std::string>());
    setRcut(typ1, typ2, r_cut);
    }


/*! Mark solid particles to remove
    set mass of a particle to -999.0

 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::mark_solid_particles_toremove(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Mark solid Particles to remove at timestep " << timestep << std::endl;

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

/*! Perform number density computation
 * This method computes and stores
     - the density based on a real mass density ( rho_i = m_i * n_neighbours / volume_sphere_of_influnece ) for fluid particles
   in the density Array.
   It overestimates density, but is superfast in comparison to compute_ndensity.
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::compute_particlenumberdensity(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Particle Number Density" << std::endl;

    // Grab handles for particle data
    ArrayHandle<Scalar> h_density(this->m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);

    unsigned int size;

    // Volume of sphere with radius m_rcut
    double sphere_vol = Scalar(4.0)/Scalar(3.0) * Scalar(3.14159265359) * pow(m_rcut, 3);

    // Particle loop
    // For each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);
        
        // All of the neighbors of this particle
        size = (unsigned int)h_n_neigh.data[i];
        // +1 because particle itself also contributes to density
        h_density.data[i] = (size + 1) * h_velocity.data[i].w / sphere_vol;

        } // End of particle loop

    } // End particle number density


/*! Perform number density computation
 * This method computes and stores
   the number density based mass density ( rho_i = m_i * \sum w_ij ) for fluid particles 
   if the SUMMATION approach is being used in the density Array.
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::compute_ndensity(uint64_t timestep)
{
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Number Density" << std::endl;

    // Grab handles for particle data
    ArrayHandle<Scalar>  h_density(this->m_pdata->getDensities(), access_location::host, access_mode::readwrite);
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
    Scalar ni;

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

        //Scalar rhoi =  h_density.data[i];

        //Scalar w_sum = 0;
        //Scalar neighbors = 0;

        for (unsigned int j = 0; j < size; j++)
        {
            // Index of neighbor
            unsigned int k = h_nlist.data[myHead + j];

            // if (checkfluid1(h_type_property_map.data, h_pos.data[k].w))
            //     continue;

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

            // Calculate distance
            Scalar r = sqrt(rsq);

            // If particle distance is too large, continue with next neighbor in loop
            if ( this->m_const_slength && rsq > m_rcutsq )
                continue;

            ni += this->m_const_slength ? this->m_skernel->wij(m_ch,r) : this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);

        } // End neighbour loop

        // Compute mass density from number density if particle i is a fluid particle
        // rho_i = m_i * \sum_j wij
        h_density.data[i] = ni * h_velocity.data[i].w;

    } // End fluid group loop

} // End compute number density


/*! Perform pressure computation
 * Simply done by a equation of state
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::compute_pressure(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Pressure" << std::endl;

    // Define ArrayHandles
    ArrayHandle<Scalar> h_density(this->m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_pressure(this->m_pdata->getPressures(), access_location::host, access_mode::readwrite);

    // For each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);
        // Evaluate pressure
        h_pressure.data[i] = this->m_eos->Pressure(h_density.data[i]);
        
        } // End fluid group loop

    } // End compute pressure

/*! Compute fictitious solid particle properties
    This method updates fictitious solid particle pressures and velocities to account for
    no-slip boundary conditions. Method follows Adami et. al. (2012).
 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::compute_noslipsolid(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::NoSlip NoPenetration" << std::endl;

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
        // // Read particle index
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
            if ( checkfluid(h_type_property_map.data, h_pos.data[k].w) )
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
            h_pressure.data[i] = this->m_eos->getBackgroundPressure();
            // Density to rest density
            h_density.data[i] = this->m_rho0;

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

            // if ( this->m_body_acceleration )
            //     {

            //     // ehemals rakulan

            //     ph_c0.x += h_density.data[k] * dx.x * wij;
            //     ph_c0.y += h_density.data[k] * dx.y * wij;
            //     ph_c0.z += h_density.data[k] * dx.z * wij;

            //     }
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
            h_density.data[i] = this->m_eos->Density(h_pressure.data[i]);
            }
        else
            {
            // Set fictitious solid velocity to zero
            h_vf.data[i].x = 0.0;
            h_vf.data[i].y = 0.0;
            h_vf.data[i].z = 0.0;

            // If no fluid neighbors are present,
            // Set pressure to background pressure
            h_pressure.data[i] = this->m_eos->getBackgroundPressure();
            // Density to rest density
            // Here: overwrite the normalisation constant
            h_density.data[i] = this->m_rho0;
            }

        } // End solid particle loop

    } // End compute noslip computation

/*! Compute fictitious solid particle properties
    This method updates fictitious solid particle pressures and velocities to account for
    no-slip boundary conditions. Method follows Adami et. al. (2012).
 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::compute_noslipsuspended(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::NoSlip NoPenetration" << std::endl;

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

    // For all suspended particles

    unsigned int group_size_suspended = this->m_suspendedgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size_suspended; group_idx++)
        {
        // // Read particle index
        unsigned int i = this->m_suspendedgroup->getMemberIndex(group_idx);

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
            if ( checkfluid(h_type_property_map.data, h_pos.data[k].w) )
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
            h_pressure.data[i] = this->m_eos->getBackgroundPressure();
            // Density to rest density
            h_density.data[i] = this->m_rho0;
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

            // if ( this->m_body_acceleration )
            //     {

            //     // ehemals rakulan

            //     ph_c0.x += h_density.data[k] * dx.x * wij;
            //     ph_c0.y += h_density.data[k] * dx.y * wij;
            //     ph_c0.z += h_density.data[k] * dx.z * wij;

            //     }
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
            h_density.data[i] = this->m_eos->Density(h_pressure.data[i]);
            }
        else
            {
            // Set fictitious solid velocity to zero
            h_vf.data[i].x = 0.0;
            h_vf.data[i].y = 0.0;
            h_vf.data[i].z = 0.0;

            // If no fluid neighbors are present,
            // Set pressure to background pressure
            h_pressure.data[i] = this->m_eos->getBackgroundPressure();
            // Density to rest density
            // Here: overwrite the normalisation constant
            h_density.data[i] = this->m_rho0;
            }

        } // End solid particle loop



    } // End compute noslip computation


/*! Perform Shepard density renormalization
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::renormalize_density(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Density renormalization" << std::endl;

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
    long unsigned int myHead;


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

/*! Compute Center of Mass of Solids in the system
 * This method computes and stores
     - the center of mass of each type in the solidgroup
   in the local variable m_centerofmasses.
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::compute_Centerofmasses(uint64_t timestep, bool print)
    {
    this->m_exec_conf->msg->notice(5) << "Suspended Object Compute Center of Mass" << std::endl;

    // Handle data on multiple cores
    // int world_size;
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Count solid types
    m_maxSuspendedID = *max_element(m_suspendedtypes.begin(), m_suspendedtypes.end());

    // Read-Handle to velocity mass and position array
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);

    // grab the box dimensions
    BoxDim globalbox = this->m_pdata->getGlobalBox();

    // Limit values
    Scalar3 lo = globalbox.getLo();
    Scalar3 Ld = globalbox.getL();

    // Find max solid type number
    unsigned int particle_typeID;

    // Loop over all solid bodies --> if there are walls, they have to have the typeID 0 --> counter starts at 1
    unsigned int start = 0;
    if (m_walls)
        start = 1;

    for (unsigned int i = start; i < m_suspendedtypes.size(); i++)
        {
        // Average position on unit circle
        Scalar angles[6] = {0,0,0,0,0,0};

        // Loop over all solid particles ---> detecte those who belong to the subgroup
        unsigned int group_size = m_suspendedgroup->getNumMembers();
        unsigned int suspendedtype_size = 0;
            for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
                {
                // Read particle index
                unsigned int j  = m_suspendedgroup->getMemberIndex(group_idx);
                particle_typeID = __scalar_as_int(h_pos.data[j].w); 
                if( particle_typeID == m_suspendedtypes[i] )
                    {
                    suspendedtype_size = suspendedtype_size+1;
                    Scalar3 pos = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);

                    // Compute mapped angles
                    Scalar theta_x = ( ( pos.x - lo.x ) / Ld.x ) * Scalar(M_TWOPI);
                    Scalar theta_y = ( ( pos.y - lo.y ) / Ld.y ) * Scalar(M_TWOPI);
                    Scalar theta_z = ( ( pos.z - lo.z ) / Ld.z ) * Scalar(M_TWOPI);

                    // Add contribution to mass averaged points on unit circle
                    angles[0] += ( Ld.x / Scalar(M_TWOPI) ) * cos(theta_x);
                    angles[1] += ( Ld.x / Scalar(M_TWOPI) ) * sin(theta_x);
                    angles[2] += ( Ld.y / Scalar(M_TWOPI) ) * cos(theta_y);
                    angles[3] += ( Ld.y / Scalar(M_TWOPI) ) * sin(theta_y);
                    angles[4] += ( Ld.z / Scalar(M_TWOPI) ) * cos(theta_z);
                    angles[5] += ( Ld.z / Scalar(M_TWOPI) ) * sin(theta_z);
                    }
                }

        // Count total number of particles in group
        unsigned int totalgroupN = suspendedtype_size;

    #ifdef ENABLE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &angles[0], 6, MPI_HOOMD_SCALAR, MPI_SUM, this->m_exec_conf->getMPICommunicator());
        MPI_Allreduce(MPI_IN_PLACE, &totalgroupN, 1, MPI_UNSIGNED, MPI_SUM, this->m_exec_conf->getMPICommunicator());
    #endif

        // Averaging
        angles[0] /= totalgroupN;
        angles[1] /= totalgroupN;
        angles[2] /= totalgroupN;
        angles[3] /= totalgroupN;
        angles[4] /= totalgroupN;
        angles[5] /= totalgroupN;

        // Set attribute
        Scalar4 centerofmass;
        // bereits initialisiert
        if ( angles[1] != 0 && angles[0] != 0 )
            centerofmass.x = lo.x + ( atan2(-angles[1], -angles[0]) + Scalar(M_PI) ) * ( Ld.x / Scalar(M_TWOPI) );
        if ( angles[3] != 0 && angles[2] != 0 )
            centerofmass.y = lo.y + ( atan2(-angles[3], -angles[2]) + Scalar(M_PI) ) * ( Ld.y / Scalar(M_TWOPI) );
        if ( angles[5] != 0 && angles[4] != 0 )
            centerofmass.z = lo.z + ( atan2(-angles[5], -angles[4]) + Scalar(M_PI) ) * ( Ld.z / Scalar(M_TWOPI) );
        
        // Set typeID
        centerofmass.w = m_suspendedtypes[i];

        if ( print and timestep==1)
            {
            this->m_exec_conf->msg->notice(2) << "SuspendedObject " << m_suspendedtypes[i] << " consists of " << totalgroupN << " slave particles" << endl;
            this->m_exec_conf->msg->notice(2) << "SuspendedObject " << m_suspendedtypes[i] << " center of Mass x: " << centerofmass.x << endl;
            this->m_exec_conf->msg->notice(2) << "SuspendedObject " << m_suspendedtypes[i] << " center of mass y: " << centerofmass.y << endl;
            this->m_exec_conf->msg->notice(2) << "SuspendedObject " << m_suspendedtypes[i] << " center of mass z: " << centerofmass.z << endl;
            this->m_exec_conf->msg->notice(2) << "SuspendedObject " << m_suspendedtypes[i] << " type ID: " << centerofmass.w << endl;
            }

    #ifdef ENABLE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &centerofmass.x, 1, MPI_HOOMD_SCALAR, MPI_MAX, this->m_exec_conf->getMPICommunicator());
        MPI_Allreduce(MPI_IN_PLACE, &centerofmass.y, 1, MPI_HOOMD_SCALAR, MPI_MAX, this->m_exec_conf->getMPICommunicator());
        MPI_Allreduce(MPI_IN_PLACE, &centerofmass.z, 1, MPI_HOOMD_SCALAR, MPI_MAX, this->m_exec_conf->getMPICommunicator());
        MPI_Allreduce(MPI_IN_PLACE, &centerofmass.w, 1, MPI_HOOMD_SCALAR, MPI_MAX, this->m_exec_conf->getMPICommunicator());
    #endif

        // Handle data on multiple cores
        // centerofmass.x = centerofmass.x/world_size;
        // centerofmass.y = centerofmass.y/world_size;
        // centerofmass.z = centerofmass.z/world_size;
        // centerofmass.w = centerofmass.w/world_size;

        // Save data in global variable     
        m_centerofmasses[i]=centerofmass;

        }

    for (unsigned int i = start; i < m_suspendedtypes.size(); i++)
        {
        if ( print and timestep==1)
            {
            this->m_exec_conf->msg->notice(2) << "SuspendedObject " << m_suspendedtypes[i] << " center of Mass x: " << m_centerofmasses[i].x << endl;
            this->m_exec_conf->msg->notice(2) << "SuspendedObject " << m_suspendedtypes[i] << " center of mass y: " << m_centerofmasses[i].y << endl;
            this->m_exec_conf->msg->notice(2) << "SuspendedObject " << m_suspendedtypes[i] << " center of mass z: " << m_centerofmasses[i].z << endl;
            this->m_exec_conf->msg->notice(2) << "SuspendedObject " << m_suspendedtypes[i] << " type ID: " << m_centerofmasses[i].w << endl;
            }
        }

    }

/*! Compute equivalent radi of solids in the system
 * This method computes and stores
     - an equivalent radi of each type in the solidgroup
   in the local variable m_radi.
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::compute_equivalentRadii(uint64_t timestep, bool print)
    {
    this->m_exec_conf->msg->notice(5) << "Suspended Object Compute equivalent radi" << std::endl;

    // Initialize fields
    // ArrayHandle<Scalar> h_density(this->m_pdata->getDensities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // access the neighbor list
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);
     
    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // unsigned int myHead, size;
    unsigned int size;
    size_t myHead;


    // Local variable to store things
    Scalar3 temp3;

    // COMPUTE EQUIVALENT RADI FOR ALL BODIES
    unsigned int start = 0;
    if (m_walls)
        start = 1;

    if ( print and timestep==1)
        {
        this->m_exec_conf->msg->notice(2) << "maximum number of solid bodies: " << m_maxSuspendedID << endl;
        this->m_exec_conf->msg->notice(2) << "Size of vector of equivalent radii: " << m_radii.size() << endl;
        }

    // Loop over all solid bodies --> if there are walls, they have to have the typeID 0 --> counter starts at 1
    for (unsigned int i = start; i < m_suspendedtypes.size(); i++)
        {
        // Access COM position
        Scalar3  pi;
        pi.x = m_centerofmasses[i].x;
        pi.y = m_centerofmasses[i].y;
        pi.z = m_centerofmasses[i].z;
        Scalar typeID_COM = m_centerofmasses[i].w;
        
        // For mean distance to COM
        // Scalar meanradi = 0.0;
        // Scalar minradi = 100.0;
        Scalar maxradi = 0.000001;
        // unsigned int containing_particles = 0;
        unsigned int typeID_j;

        // Loop over all solid particles ---> detecte those who belong to the subgroup
        unsigned int group_size = m_suspendedgroup->getNumMembers();
        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            // Read particle index
            unsigned int j  = m_suspendedgroup->getMemberIndex(group_idx);
            typeID_j = __scalar_as_int(h_pos.data[j].w);

            if (typeID_j == (start-1) or (typeID_j != typeID_COM))
                continue;

            // Only check radi for particles with fluid neighbors
            // Loop over neighbors
            bool solid_w_fluid_neigh = false;
            myHead = h_head_list.data[j];
            size = (unsigned int)h_n_neigh.data[j];
            for (unsigned int k = 0; k < size; k++)
                {
                // Index of neighbor (MEM TRANSFER: 1 scalar)
                unsigned int a = h_nlist.data[myHead + k];

                // Check if neighbor is fluid
                if ( checkfluid(h_type_property_map.data, h_pos.data[a].w) )
                    {
                    solid_w_fluid_neigh = true;
                    //break;
                    }
                }

            // If fluid neighbors exist add part to mean radi
            if (solid_w_fluid_neigh)
                {
                // Access the particle's position, velocity, mass and type
                // Scalar3 pj = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);
                Scalar3 pj;
                pj.x = h_pos.data[j].x;
                pj.y = h_pos.data[j].y;
                pj.z = h_pos.data[j].z;

                // Compute distance vector
                Scalar3 dx = pi - pj;

                // Apply periodic boundary conditions
                dx = box.minImage(dx);

                // Calculate squared distance
                Scalar rsq = dot(dx, dx);

                // Calculate absolute and normalized distance
                Scalar r = sqrt(rsq);
                //minradi = r; // ==========>>>>>>>>>>>>> TEST DELETE LATER

                // Add part to sum of radii
                // meanradi = meanradi + r;
                // containing_particles = containing_particles + 1;
                // if (r < minradi)
                //     minradi = r;
                if (r > maxradi)
                    maxradi = r;


                // If velocity for this group was not saved do it here =====> Necessary for other contact model (based on distance of surface particles) TODO NK
                if (m_velocities[i].x == 0.0)
                    {
                    temp3.x = h_velocity.data[j].x;
                    temp3.y = h_velocity.data[j].y;
                    temp3.z = h_velocity.data[j].z;
                    m_velocities[i] = temp3;
                    }

                }

            }

    // Reduce on all processors
    #ifdef ENABLE_MPI
        // MPI_Allreduce(MPI_IN_PLACE, &meanradi, 1, MPI_HOOMD_SCALAR, MPI_SUM, this->m_exec_conf->getMPICommunicator());
        // MPI_Allreduce(MPI_IN_PLACE, &minradi, 1, MPI_HOOMD_SCALAR, MPI_MIN, this->m_exec_conf->getMPICommunicator());
        MPI_Allreduce(MPI_IN_PLACE, &maxradi, 1, MPI_HOOMD_SCALAR, MPI_MAX, this->m_exec_conf->getMPICommunicator());
        // MPI_Allreduce(MPI_IN_PLACE, &containing_particles, 1, MPI_UNSIGNED, MPI_SUM, this->m_exec_conf->getMPICommunicator());
    #endif

        // Calculate form factor out of minradi and maxradi
        // Scalar formfactor = (minradi+maxradi)/2.0;
        // Calculate meanradi
        // m_radii[i] = ((meanradi/containing_particles)+maxradi)/2.0;
        m_radii[i] = maxradi;

        if ( print and timestep==1)
            {
            this->m_exec_conf->msg->notice(2) << "Radi of solid body: " << typeID_COM << ": " << maxradi << endl;
            }

        }
        
    }


/*! Compute repulsive force due to contact of two solids based on Bian, Ellero, CompPhysComm (2014)
 * This method computes and stores
     - the repulsive force of each type in the solidgroup
   in the global variable aux2.
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::compute_repulsiveForce(uint64_t timestep, bool print)
    {

    this->m_exec_conf->msg->notice(5) << "Suspended Object Compute Repulsive Forces" << std::endl;


    // Check input data, can be omitted if need be
    //assert(h_conforce.data);

    // Read-Handle to velocity mass and position array
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_conforce(this->m_pdata->getAuxiliaries2(), access_location::host,access_mode::readwrite);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // access the neighbor list
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Zero data before calculation
    memset((void*)h_conforce.data,0,sizeof(Scalar3)*this->m_pdata->getAuxiliaries2().getNumElements());

    // Handle data on multiple cores
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Set start point if there are walls or not
    unsigned int start = 0;
    if (m_walls)
        start = 1;

    // Local variable to store things
    Scalar temp0;
    Scalar temp1;
    Scalar r_ref;
    Scalar F_total[3] = {0,0,0};
    // Scalar dist = 100.0;
    Scalar F_0 = m_f0;

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Loop over all solid bodies --> if there are walls, they have to have the typeID 0 --> counter starts at 1
    for (unsigned int i = start; i < m_suspendedtypes.size(); i++)
        {
        // Access COM position of solid i
        Scalar3  pi;
        pi.x = m_centerofmasses[i].x;
        pi.y = m_centerofmasses[i].y;
        pi.z = m_centerofmasses[i].z;
        Scalar typeID_i = m_centerofmasses[i].w;
        Scalar radi_i = m_radii[i];

        Scalar F_rep[3] = {0,0,0};

        for (unsigned int j = start; j < m_suspendedtypes.size(); j++)
            {   
            if (i == j)
                continue;
            
            // Access COM position of solid k
            Scalar3  pj;
            pj.x = m_centerofmasses[j].x;
            pj.y = m_centerofmasses[j].y;
            pj.z = m_centerofmasses[j].z;
            Scalar radi_j = m_radii[j];
            r_ref = (radi_i + radi_j)/Scalar(2.0);
            // F_0   = r_ref; // F_0 should be in the range of the effective radius

            // Compute distance vector
            Scalar3 dx = pi - pj;

            // Apply periodic boundary conditions
            dx = box.minImage(dx);

            // Calculate squared distance
            Scalar rsq = dot(dx, dx);

            // Calculate absolute and normalized distance
            Scalar r = sqrt(rsq);

            // s from the theoretical model
            Scalar dist = r - radi_i - radi_j;

            // If particle distance is too large, skip this loop
            // if (dist > 0.2*radi_i)
            //     continue;

            // tau from the theoretical model
            temp0 = 1/(0.01*r_ref);
            temp1 = (F_0 * temp0 * exp(-temp0*dist)) / (1 - exp(-temp0*dist));

            // Only add contact force if distance is smaller than 20%  /// 10 %  of radi and F_rep > e-5
            // Add F_rep in opposite direction of the vector between the centers of mass of two bodies
            if (dist < 0.2*radi_j && sqrt(temp1*temp1) > 0.00001)
                {
                F_rep[0] =+ temp1*dx.x;
                F_rep[1] =+ temp1*dx.y;
                F_rep[2] =+ temp1*dx.z;
                if (print and (sqrt(F_rep[0]*F_rep[0]+F_rep[1]*F_rep[1]+F_rep[2]*F_rep[2]))>0.0000001 and timestep%10)
                    {
                    this->m_exec_conf->msg->notice(2) << "Found contact between solid " << i << " and solid " << j <<"!" << endl;
                    this->m_exec_conf->msg->notice(2) << "Contact force on solid particle "<< i << ": Fx = " << F_rep[0] << " Fy = " << F_rep[1] << ": Fz = " << F_rep[2] << endl;
                    //this->m_exec_conf->msg->notice(2) << "Total contact force on solid particle "<< i << ": F_rep = " << (sqrt(F_rep[0]*F_rep[0]+F_rep[1]*F_rep[1]+F_rep[2]*F_rep[2])) << endl; 
                    }
                }

            } // end over solid bodies that might be in contact

    #ifdef ENABLE_MPI
        MPI_Allreduce(&F_rep, &F_total, 3, MPI_HOOMD_SCALAR, MPI_SUM, this->m_exec_conf->getMPICommunicator());
    #endif

        // Handle data on multiple cores
        F_total[0] = F_total[0]/world_size;
        F_total[1] = F_total[1]/world_size;
        F_total[2] = F_total[2]/world_size;

        // if (print and (sqrt(F_total[0]*F_total[0]+F_total[1]*F_total[1]+F_total[2]*F_total[2]))>0.0 and timestep%10 )
        //     {
        //     this->m_exec_conf->msg->notice(2) << "Total contact force on solid particle "<< i << ": Fx = " << F_total[0] << " Fy = " << F_total[1] << ": Fz = " << F_total[2] << endl; 
        //     // this->m_exec_conf->msg->notice(2) << "Total contact force on solid particle "<< i << ": F_total = " << (sqrt(F_total[0]*F_total[0]+F_total[1]*F_total[1]+F_total[2]*F_total[2])) << endl; 
        //     }

        // Loop over solid particles
        unsigned int group_size = m_suspendedgroup->getNumMembers();
        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            unsigned int k = m_suspendedgroup->getMemberIndex(group_idx);
            Scalar typeID_k = __scalar_as_int(h_pos.data[k].w);

            // Solid-solid contact force
            // If type id of solid particle matches id of solid body, force is applied. Counterpart is applied when i gets to this id
            if (typeID_k == typeID_i)
                {
                // Repulsive force
                h_conforce.data[k].x =+ F_total[0];
                h_conforce.data[k].y =+ F_total[1];
                h_conforce.data[k].z =+ F_total[2];

                // Wall contact force
                if (m_walls)
                    {
                    bool foundcontact = false;

                    // Loop over all of the neighbors of this particle
                    // const unsigned int myHead = h_head_list.data[i];
                    // const unsigned int size = (unsigned int)h_n_neigh.data[k];
                    const size_t myHead = h_head_list.data[i];
                    const unsigned int size = (unsigned int)h_n_neigh.data[k];
                    for (unsigned int f = 0; f < size; f++)
                        {

                        // Index of neighbor 
                        unsigned int g = h_nlist.data[myHead + f];

                        // Sanity check
                        assert(g < this->m_pdata->getN() + this->m_pdata->getNGhosts());

                        // Determine neighbor type
                        Scalar typeID_g = __scalar_as_int(h_pos.data[g].w);

                        if(foundcontact)
                            continue;
                        
                        // Only compute wall force if neighbor is a wall particle (= typeID 0)
                        if (typeID_g == 0)
                            {
                            // Access neighbor position
                            Scalar3 pk;
                            pk.x = h_pos.data[k].x;
                            pk.y = h_pos.data[k].y;
                            pk.z = h_pos.data[k].z;

                            Scalar3 pg;
                            pg.x = h_pos.data[g].x;
                            pg.y = h_pos.data[g].y;
                            pg.z = h_pos.data[g].z;

                            // Compute distance vector
                            Scalar3 walldist = pg - pk;

                            // Apply periodic boundary conditions
                            walldist = box.minImage(walldist);

                            // Calculate squared distance
                            Scalar wdistsq = dot(walldist, walldist);

                            // Calculate absolute and normalized distance
                            Scalar wd = sqrt(wdistsq);

                            // If particle distance is too large, skip this loop
                            if (wd > 0.2*radi_i)
                                continue;

                            if (print and timestep%10 and (wd<0.0001))
                                {
                                this->m_exec_conf->msg->notice(2) << "Found wall contact!" << endl; 
                                }

                            // tau from the theoretical model
                            temp0 = 1/(0.01*radi_i);
                            temp1 = (0.1*F_0 * temp0 * exp(-temp0*wd)) / (1 - exp(-temp0*wd));

                            // Add repulsive force
                            h_conforce.data[k].x =+ temp1*walldist.x;
                            h_conforce.data[k].y =+ temp1*walldist.y;
                            h_conforce.data[k].z =+ temp1*walldist.z;

                            foundcontact = true;

                            }

                        } // End neighbor loop

                    } // End m_walls

                } // End type_i == type_k    

            } // end loop over solid particles

        } // end loop over solid bodies

    }

/*! Perform force computation
 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::forcecomputation(uint64_t timestep)
    {

    if ( m_density_method == DENSITYSUMMATION )
        this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Forces using SUMMATION approach " << m_density_method << endl;
    else if ( m_density_method == DENSITYCONTINUITY )
        this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlow::Forces using CONTINUITY approach " << m_density_method << endl;

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
    ArrayHandle<Scalar3> h_conforce(this->m_pdata->getAuxiliaries2(), access_location::host,access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

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

        // // Total velocity of particle
        Scalar vi_total = sqrt((vi.x * vi.x) + (vi.y * vi.y) + (vi.z * vi.z));

        // Properties needed for adaptive timestep
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
            bool issolid = checksolid(h_type_property_map.data, h_pos.data[k].w);
            bool issuspended = checksuspended(h_type_property_map.data, h_pos.data[k].w);
            // If particle is another suspended particle continue (TODO: check if right)
            // if ( issuspended )
            //     continue;

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
            if ( m_const_slength && rsq > m_rcutsq )
                continue;

            // Access neighbor velocity; depends on fluid or fictitious solid particle
            Scalar3 vj  = make_scalar3(0.0, 0.0, 0.0);
            Scalar mj   = h_velocity.data[k].w;
            if ( issolid || issuspended )
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
            Scalar meanh  = m_const_slength ? m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar eps    = Scalar(0.1)*meanh;
            Scalar epssqr = eps*eps;

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);

            // Evaluate inter-particle pressure forces
            //temp0 = -((mi*mj)/(rhoj*rhoi))*(Pi+Pj);
            //temp0 = -Vi*Vj*( Pi + Pj );
            //temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);
            //temp0 = -mi*mj*( Pi/(rhoi*rhoj) + Pj/(rhoj*rhoj) );
            if ( m_density_method == DENSITYSUMMATION )
            {
                // Transport formulation proposed by Adami 2013
                temp0 = -(Vi*Vi+Vj*Vj)*((rhoj*Pi+rhoi*Pj)/(rhoi+rhoj)); 
            }
            else if ( m_density_method == DENSITYCONTINUITY) 
            { 
                temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);
            }

            // Optionally add artificial viscosity
            // Monaghan (1983) J. Comput. Phys. 52 (2) 374–389
            // TODO: Hier Klammer?
            if ( m_artificial_viscosity && (!issolid || !issuspended))
                {
                Scalar dotdvdx = dot(dv,dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar muij    = meanh*dotdvdx/(rsq+epssqr);
                    Scalar meanrho = Scalar(0.5)*(rhoi+rhoj);
                    temp0 += mi*mj*(m_avalpha*m_c*muij+m_avbeta*muij*muij)/meanrho;
                    }
                }

            // Add contribution to suspended particle
            h_force.data[i].x += temp0*dwdr_r*dx.x;
            h_force.data[i].y += temp0*dwdr_r*dx.y;
            h_force.data[i].z += temp0*dwdr_r*dx.z;

            // // Add contribution to solid particle
            // if ( issuspended && m_compute_solid_forces )
            //     {
            //     h_force.data[k].x -= (mj/mi)*temp0*dwdr_r*dx.x;
            //     h_force.data[k].y -= (mj/mi)*temp0*dwdr_r*dx.y;
            //     h_force.data[k].z -= (mj/mi)*temp0*dwdr_r*dx.z;
            //     }

            // Evaluate viscous interaction forces
            temp0 = m_mu * (Vi*Vi+Vj*Vj) * dwdr_r;
            h_force.data[i].x  += temp0*dv.x;
            h_force.data[i].y  += temp0*dv.y;
            h_force.data[i].z  += temp0*dv.z;

            // // Add contribution to solid particle
            // if ( issuspended && m_compute_solid_forces )
            //     {
            //     h_force.data[k].x -= (mj/mi)*temp0*dv.x;
            //     h_force.data[k].y -= (mj/mi)*temp0*dv.y;
            //     h_force.data[k].z -= (mj/mi)*temp0*dv.z;
            //     }

            // Evaluate rate of change of density if CONTINUITY approach is used
            if ( m_density_method == DENSITYCONTINUITY )
                {
                if ( issolid || issuspended )
                    {
                    // Use physical advection velocity rather than fictitious velocity here
                    vj.x = h_velocity.data[k].x;
                    vj.y = h_velocity.data[k].y;
                    vj.z = h_velocity.data[k].z;

                    // Recompute velocity difference
                    // dv = vi - vj;
                    dv.x = vi.x - vj.x;
                    dv.y = vi.y - vj.y;
                    dv.z = vi.z - vj.z;
                    //Vj = mj / m_rho0;
                    }

                // Compute density rate of change
                // std::cout << "Compute density rate of change: rhoi " << rhoi << " Vj " << Vj << " dot(dv,dwdr_r*dx) " << dot(dv,dwdr_r*dx) << std::endl;
                h_ratedpe.data[i].x += rhoi*Vj*dot(dv,dwdr_r*dx);
                // std::cout << "Compute density rate of change: h_ratedpe.data[i].x " << h_ratedpe.data[i].x << std::endl;

                //h_ratedpe.data[i].x += mj*dot(dv,dwdr_r*dx);

                // Add density diffusion if requested
                // Molteni and Colagrossi, Computer Physics Communications 180 (2009) 861–872
                if ( (!issolid || !issuspended) && m_density_diffusion )
                    h_ratedpe.data[i].x -= (Scalar(2)*m_ddiff*meanh*m_c*mj*(rhoi/rhoj-Scalar(1))*dot(dx,dwdr_r*dx))/(rsq+epssqr);
                }

            } // Closing Neighbor Loop

        } // Closing Fluid Particle Loop

    // for each suspended particle
    unsigned int group_size_suspended = m_suspendedgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size_suspended; group_idx++)
        {
        // Read particle index
        unsigned int i = m_suspendedgroup->getMemberIndex(group_idx);

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

        // // Total velocity of particle
        Scalar vi_total = sqrt((vi.x * vi.x) + (vi.y * vi.y) + (vi.z * vi.z));

        // Properties needed for adaptive timestep
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

            // TODO: Check if solid is valid but not suspended
            // Determine neighbor type
            bool issolid = checksolid(h_type_property_map.data, h_pos.data[k].w);
            bool issuspended = checksuspended(h_type_property_map.data, h_pos.data[k].w);
            if (issuspended)
                continue;

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
            if ( m_const_slength && rsq > m_rcutsq )
                continue;

            // Access neighbor velocity; depends on fluid or fictitious solid particle
            Scalar3 vj  = make_scalar3(0.0, 0.0, 0.0);
            Scalar mj   = h_velocity.data[k].w;
            if ( issolid )
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
            Scalar meanh  = m_const_slength ? m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar eps    = Scalar(0.1)*meanh;
            Scalar epssqr = eps*eps;

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);

            // Evaluate inter-particle pressure forces
            //temp0 = -((mi*mj)/(rhoj*rhoi))*(Pi+Pj);
            //temp0 = -Vi*Vj*( Pi + Pj );
            //temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);
            //temp0 = -mi*mj*( Pi/(rhoi*rhoj) + Pj/(rhoj*rhoj) );
            if ( m_density_method == DENSITYSUMMATION )
            {
                // Transport formulation proposed by Adami 2013
                temp0 = -(Vi*Vi+Vj*Vj)*((rhoj*Pi+rhoi*Pj)/(rhoi+rhoj)); 
            }
            else if ( m_density_method == DENSITYCONTINUITY) 
            { 
                temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);
            }


            // // TODO: hier überlegen
            // // Optionally add artificial viscosity
            // // Monaghan (1983) J. Comput. Phys. 52 (2) 374–389
            // if ( m_artificial_viscosity && !issolid )
            //     {
            //     Scalar dotdvdx = dot(dv,dx);
            //     if ( dotdvdx < Scalar(0) )
            //         {
            //         Scalar muij    = meanh*dotdvdx/(rsq+epssqr);
            //         Scalar meanrho = Scalar(0.5)*(rhoi+rhoj);
            //         temp0 += mi*mj*(m_avalpha*m_c*muij+m_avbeta*muij*muij)/meanrho;
            //         }
            //     }

            // Add contribution to fluid particle
            h_force.data[i].x += temp0*dwdr_r*dx.x;
            h_force.data[i].y += temp0*dwdr_r*dx.y;
            h_force.data[i].z += temp0*dwdr_r*dx.z;

            // // Add contribution to solid particle
            // if ( issolid && m_compute_solid_forces )
            //     {
            //     h_force.data[k].x -= (mj/mi)*temp0*dwdr_r*dx.x;
            //     h_force.data[k].y -= (mj/mi)*temp0*dwdr_r*dx.y;
            //     h_force.data[k].z -= (mj/mi)*temp0*dwdr_r*dx.z;
            //     }

            // Evaluate viscous interaction forces
            temp0 = m_mu * (Vi*Vi+Vj*Vj) * dwdr_r;
            h_force.data[i].x  += temp0*dv.x;
            h_force.data[i].y  += temp0*dv.y;
            h_force.data[i].z  += temp0*dv.z;

            // // Add contribution to solid particle
            // if ( issolid && m_compute_solid_forces )
            //     {
            //     h_force.data[k].x -= (mj/mi)*temp0*dv.x;
            //     h_force.data[k].y -= (mj/mi)*temp0*dv.y;
            //     h_force.data[k].z -= (mj/mi)*temp0*dv.z;
            //     }

            // // Evaluate rate of change of density if CONTINUITY approach is used
            // if ( m_density_method == DENSITYCONTINUITY )
            //     {
            //     if ( issolid )
            //         {
            //         // Use physical advection velocity rather than fictitious velocity here
            //         vj.x = h_velocity.data[k].x;
            //         vj.y = h_velocity.data[k].y;
            //         vj.z = h_velocity.data[k].z;

            //         // Recompute velocity difference
            //         // dv = vi - vj;
            //         dv.x = vi.x - vj.x;
            //         dv.y = vi.y - vj.y;
            //         dv.z = vi.z - vj.z;
            //         //Vj = mj / m_rho0;
            //         }

            //     // Compute density rate of change
            //     // std::cout << "Compute density rate of change: rhoi " << rhoi << " Vj " << Vj << " dot(dv,dwdr_r*dx) " << dot(dv,dwdr_r*dx) << std::endl;
            //     h_ratedpe.data[i].x += rhoi*Vj*dot(dv,dwdr_r*dx);
            //     // std::cout << "Compute density rate of change: h_ratedpe.data[i].x " << h_ratedpe.data[i].x << std::endl;

            //     //h_ratedpe.data[i].x += mj*dot(dv,dwdr_r*dx);

            //     // Add density diffusion if requested
            //     // Molteni and Colagrossi, Computer Physics Communications 180 (2009) 861–872
            //     if ( !issolid && m_density_diffusion )
            //         h_ratedpe.data[i].x -= (Scalar(2)*m_ddiff*meanh*m_c*mj*(rhoi/rhoj-Scalar(1))*dot(dx,dwdr_r*dx))/(rsq+epssqr);
            //     }

            } // Closing Neighbor Loop

        } // Closing Suspended Particle Loop

    m_timestep_list[5] = max_vel;
    // Add volumetric force (gravity)
    this->applyBodyForce(timestep, m_fluidgroup);
    if ( m_compute_solid_forces )
        {
        this->applyBodyForce(timestep, m_suspendedgroup);

        // for each particle in solid group
        unsigned int group_size = m_suspendedgroup->getNumMembers();
        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            // Read particle index
            unsigned int j = m_suspendedgroup->getMemberIndex(group_idx);

            // Add contribution of contact force to total force
            h_force.data[j].x += h_conforce.data[j].x;
            h_force.data[j].y += h_conforce.data[j].y;
            h_force.data[j].z += h_conforce.data[j].z;
            }
        }


    }



/*! Compute forces definition
*/

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::computeForces(uint64_t timestep)
    {

    // start by updating the neighborlist
    this->m_nlist->compute(timestep);

    // This is executed once to initialize protected/private variables
    if (!m_params_set)
        {
        this->m_exec_conf->msg->error() << "sph.models.SuspensionFlow requires parameters to be set before run()"
            << std::endl;
        throw std::runtime_error("Error computing SuspensionFlow forces");
        }

    // start by updating the neighborlist
    this->m_nlist->compute(timestep);

    // m_solid_removed flag is set to False initially, so this 
    // only executes at timestep 0
    if (!m_solid_removed)
        {
        this->m_nlist->forceUpdate();
        this->m_nlist->compute(timestep);
        mark_solid_particles_toremove(timestep);
        m_solid_removed = true;
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
    update_ghost_density_pressure(timestep);
#endif

    // Compute particle pressures
    // Includes the computation of the density of solid particles
    // based on ficticios pressure p_i^\ast
    compute_noslipsolid(timestep);

    // Compute particle pressures
    // Includes the computation of the density of solid particles
    // based on ficticios pressure p_i^\ast
    compute_noslipsuspended(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    update_ghost_aux1(timestep);
#endif

    // Compute center of mass of solid bodies
    compute_Centerofmasses(timestep, true);

    // Compute equivalent radi of solid bodies ============= >>> evtl als Vektor speichern (x,y,z)
    compute_equivalentRadii(timestep, true);

    // Compute repulsive forces
    compute_repulsiveForce(timestep, true);

    // communicate repulsive forces
#ifdef ENABLE_MPI
    // Update ghost particles
    update_ghost_aux2(timestep);
#endif

    // Execute the force computation
    // This includes the computation of the density if 
    // DENSITYCONTINUITY method is used
    forcecomputation(timestep);

    }

namespace detail 
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SuspensionFlow(pybind11::module& m, std::string name)
{
    pybind11::class_<SuspensionFlow<KT_, SET_>, SPHBaseClass<KT_, SET_> , std::shared_ptr<SuspensionFlow<KT_, SET_>>>(m, name.c_str()) 
        .def(pybind11::init< std::shared_ptr<SystemDefinition>,
                             std::shared_ptr<SmoothingKernel<KT_> >,
                             std::shared_ptr<StateEquation<SET_> >,
                             std::shared_ptr<nsearch::NeighborList>,
                             std::shared_ptr<ParticleGroup>,
                             std::shared_ptr<ParticleGroup>,
                             std::shared_ptr<ParticleGroup>,
                             DensityMethod,
                             ViscosityMethod >())
        .def("setParams", &SuspensionFlow<KT_, SET_>::setParams)
        .def("getDensityMethod", &SuspensionFlow<KT_, SET_>::getDensityMethod)
        .def("setDensityMethod", &SuspensionFlow<KT_, SET_>::setDensityMethod)
        .def("getViscosityMethod", &SuspensionFlow<KT_, SET_>::getViscosityMethod)
        .def("setViscosityMethod", &SuspensionFlow<KT_, SET_>::setViscosityMethod)
        .def("setConstSmoothingLength", &SuspensionFlow<KT_, SET_>::setConstSmoothingLength)
        .def("computeSolidForces", &SuspensionFlow<KT_, SET_>::computeSolidForces)
        .def("computeWallForces", &SuspensionFlow<KT_, SET_>::computeWallForces)
        .def("activateArtificialViscosity", &SuspensionFlow<KT_, SET_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity", &SuspensionFlow<KT_, SET_>::deactivateArtificialViscosity)
        .def("activateDensityDiffusion", &SuspensionFlow<KT_, SET_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion", &SuspensionFlow<KT_, SET_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &SuspensionFlow<KT_, SET_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization", &SuspensionFlow<KT_, SET_>::deactivateShepardRenormalization)
        .def("setAcceleration", &SPHBaseClass<KT_, SET_>::setAcceleration)
        .def("setRCut", &SuspensionFlow<KT_, SET_>::setRCutPython)
        ;

    }

} // end namespace detail

//! Explicit template instantiations
template class PYBIND11_EXPORT SuspensionFlow<wendlandc2, linear>;
template class PYBIND11_EXPORT SuspensionFlow<wendlandc2, tait>;
template class PYBIND11_EXPORT SuspensionFlow<wendlandc4, linear>;
template class PYBIND11_EXPORT SuspensionFlow<wendlandc4, tait>;
template class PYBIND11_EXPORT SuspensionFlow<wendlandc6, linear>;
template class PYBIND11_EXPORT SuspensionFlow<wendlandc6, tait>;
template class PYBIND11_EXPORT SuspensionFlow<quintic, linear>;
template class PYBIND11_EXPORT SuspensionFlow<quintic, tait>;
template class PYBIND11_EXPORT SuspensionFlow<cubicspline, linear>;
template class PYBIND11_EXPORT SuspensionFlow<cubicspline, tait>;


namespace detail
{

    template void export_SuspensionFlow<wendlandc2, linear>(pybind11::module& m, std::string name = "SuspensionPF_WC2_L");
    template void export_SuspensionFlow<wendlandc2, tait>(pybind11::module& m, std::string name = "SuspensionF_WC2_T");
    template void export_SuspensionFlow<wendlandc4, linear>(pybind11::module& m, std::string name = "SuspensionF_WC4_L");
    template void export_SuspensionFlow<wendlandc4, tait>(pybind11::module& m, std::string name = "SuspensionF_WC4_T");
    template void export_SuspensionFlow<wendlandc6, linear>(pybind11::module& m, std::string name = "SuspensionF_WC6_L");
    template void export_SuspensionFlow<wendlandc6, tait>(pybind11::module& m, std::string name = "SuspensionF_WC6_T");
    template void export_SuspensionFlow<quintic, linear>(pybind11::module& m, std::string name = "SuspensionF_Q_L");
    template void export_SuspensionFlow<quintic, tait>(pybind11::module& m, std::string name = "SuspensionF_Q_T");
    template void export_SuspensionFlow<cubicspline, linear>(pybind11::module& m, std::string name = "SuspensionF_CS_L");
    template void export_SuspensionFlow<cubicspline, tait>(pybind11::module& m, std::string name = "SuspensionF_CS_T");

} // end namespace detail
} // end namespace sph
} // end namespace hoomd

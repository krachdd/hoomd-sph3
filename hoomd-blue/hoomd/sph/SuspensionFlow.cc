/* ---------------------------------------------------------
maintainer: drostan, daniel.rostan@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SuspensionFlow.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include "hoomd/VectorMath.h"

#include <map>
#include <sstream>
#include <string.h>

using namespace std;

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
                                 std::shared_ptr<ParticleGroup> aggregategroup,
                                 std::shared_ptr<ParticleGroup> comgroup,
                                 DensityMethod mdensitymethod,
                                 ViscosityMethod mviscositymethod)
    : SPHBaseClassConstraint<KT_, SET_>(sysdef,skernel,equationofstate,nlist), m_fluidgroup(fluidgroup),
      m_solidgroup(solidgroup), m_aggregategroup(aggregategroup), m_comgroup(comgroup), m_typpair_idx(this->m_pdata->getNTypes()),
      m_bodies_changed(false), m_particles_added_removed(false)
      {
        // this->m_pdata->getGlobalParticleNumberChangeSignal()
        // .connect<SuspensionFlow<KT_, SET_>, &SuspensionFlow<KT_, SET_>::slotPtlsAddedRemoved>(this);

        // this->m_exec_conf->msg->notice(5) << "Constructing SuspensionFlow" << std::endl;
        this->m_exec_conf->msg->notice(5) << "Constructing SuspensionFlow" << std::endl;

        //std::cout << "Types:" << this->m_pdata->getNTypes() << std::endl;

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
        m_rhoS = Scalar(1000.0);

        m_solid_removed = false;

        // Sanity checks
        assert(this->m_pdata);
        assert(this->m_nlist);
        assert(this->m_skernel);
        assert(this->m_eos);

        // Contruct type vectors
        this->constructTypeVectors(fluidgroup,&m_fluidtypes);
        this->constructTypeVectors(solidgroup,&m_solidtypes);
        this->constructTypeVectors(aggregategroup,&m_aggregatetypes);
        this->constructTypeVectors(comgroup,&m_comtypes);


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
            for (unsigned int i = 0; i < m_aggregatetypes.size(); i++) {
                h_type_property_map.data[m_aggregatetypes[i]] |= SolidFluidTypeBit::FLUID1;
            }
            for (unsigned int i = 0; i < m_comtypes.size(); i++) {
                h_type_property_map.data[m_comtypes[i]] |= SolidFluidTypeBit::FLUID2;
            }
        }

        this->m_exec_conf->msg->notice(2) << "Computing SuspensionFlow::number of aggregates in m_aggregatetypes " << m_aggregatetypes.size() << std::endl;

        // // Init vectors for suspension model
        // Scalar3 zeros3 = make_scalar3(0.0,0.0,0.0);
        // Scalar4 zeros4 = make_scalar4(0.0,0.0,0.0,0.0);
        // for (unsigned int i = 0; i < m_solidtypes.size(); ++i)
        //     {
        //     m_centerofmasses.push_back(zeros4);
        //     m_repulsiveforces.push_back(zeros3);
        //     m_velocities.push_back(zeros3);
        //     m_radii.push_back(0.0);
        //     }

        //'Set simulations methods
        m_density_method = mdensitymethod;
        m_viscosity_method = mviscositymethod;

        // Get necessary variables from kernel and EOS classes
        m_rho0  = equationofstate->getRestDensity();
        m_c     = equationofstate->getSpeedOfSound();
        m_kappa = skernel->getKernelKappa();

        m_r_cut_nlist = std::make_shared<GlobalArray<Scalar>>(m_typpair_idx.getNumElements(), this->m_exec_conf);
        this->m_nlist->addRCutMatrix(m_r_cut_nlist);


        GlobalArray<unsigned int> body_types(this->m_pdata->getNTypes(), 1, this->m_exec_conf);
        m_body_types.swap(body_types);
        TAG_ALLOCATION(m_body_types);

        GlobalArray<Scalar3> body_pos(this->m_pdata->getNTypes(), 1, this->m_exec_conf);
        m_body_pos.swap(body_pos);
        TAG_ALLOCATION(m_body_pos);

        // GlobalArray<Scalar4> body_orientation(m_pdata->getNTypes(), 1, m_exec_conf);
        // m_body_orientation.swap(body_orientation);
        // TAG_ALLOCATION(m_body_orientation);

        GlobalArray<unsigned int> body_len(this->m_pdata->getNTypes(), this->m_exec_conf);
        m_body_len.swap(body_len);
        TAG_ALLOCATION(m_body_len);

        // After the constructor is called
        // Count the number of elements in m_body_types array
        unsigned int numBodyTypes = m_body_types.getNumElements();

        // Count the number of elements in m_body_len array
        unsigned int numBodyLen = m_body_len.getNumElements();

        // reset elements to zero
        ArrayHandle<unsigned int> h_body_len(m_body_len, access_location::host, access_mode::readwrite);
        for (unsigned int i = 0; i < this->m_pdata->getNTypes(); ++i)
            {
            h_body_len.data[i] = 0;
            }

        // m_body_charge.resize(m_pdata->getNTypes());
        // m_body_diameter.resize(m_pdata->getNTypes());

        m_d_max.resize(this->m_pdata->getNTypes(), Scalar(0.0));
        m_d_max_changed.resize(this->m_pdata->getNTypes(), false);

    #ifdef ENABLE_MPI
        if (this->m_sysdef->isDomainDecomposed())
            {
            auto comm_weak = this->m_sysdef->getCommunicator();
            assert(comm_weak.lock());
            m_comm = comm_weak.lock();

            // register this class with the communicator
            m_comm->getBodyGhostLayerWidthRequestSignal()
                .connect<SuspensionFlow, &SuspensionFlow::requestBodyGhostLayerWidth>(this);
            }
    #endif

      }

/*! Destructor
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
SuspensionFlow<KT_, SET_>::~SuspensionFlow()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying SuspensionFlow" << std::endl;

    // // disconnect from signal in ParticleData;
    // this->m_pdata->getGlobalParticleNumberChangeSignal()
    //     .disconnect<SuspensionFlow, &SuspensionFlow::slotPtlsAddedRemoved>(this);
    // #ifdef ENABLE_MPI
    //     if (this->m_sysdef->isDomainDecomposed())
    //         {
    //         m_comm->getBodyGhostLayerWidthRequestSignal()
    //             .disconnect<SuspensionFlow, &SuspensionFlow::requestBodyGhostLayerWidth>(this);
    //         }
    // #endif
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

// not sure if rho0S needed
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::setParams(Scalar mu, Scalar rhoS0)
    {
    m_mu   = mu;
    m_rhoS = rhoS0;

    if (m_mu <= 0)
         {
         this->m_exec_conf->msg->error() << "sph.models.SuspensionFlow: Dynamic viscosity has to be a positive real number" << std::endl;
         throw std::runtime_error("Error initializing SuspensionFlow.");
         }

    if (m_rhoS <= 0)
         {
         this->m_exec_conf->msg->error() << "sph.models.SuspensionFlow: Solid density has to be a positive real number" << std::endl;
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

        Scalar rhoi =  h_density.data[i];

        Scalar w_sum = 0;
        Scalar neighbors = 0;

        for (unsigned int j = 0; j < size; j++)
        {
            // Index of neighbor
            unsigned int k = h_nlist.data[myHead + j];

            if (checkfluid1(h_type_property_map.data, h_pos.data[k].w))
                continue;

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
// TODO: Muss angepasst werden zwecks Aggregate - Wall Interaktion

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::compute_noslip(uint64_t timestep)
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

            //     ph_c0.x += h_dpe.data[k].x * dx.x * wij;
            //     ph_c0.y += h_dpe.data[k].x * dx.y * wij;
            //     ph_c0.z += h_dpe.data[k].x * dx.z * wij;

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

// TODO: Compute equivalentRadii

// TODO: Compute repulsiveForce
template<SmoothingKernelType KT_, StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::setRigidParams(unsigned int body_typeid,
                              std::vector<unsigned int>& type,
                              std::vector<Scalar3>& pos)
                              // std::vector<Scalar4>& orientation,
                              // std::vector<Scalar>& charge,
                              // std::vector<Scalar>& diameter)
    {
    assert(m_body_types.getPitch() >= this->m_pdata->getNTypes());
    assert(m_body_pos.getPitch() >= this->m_pdata->getNTypes());
    // assert(m_body_orientation.getPitch() >= this->m_pdata->getNTypes());
    // assert(m_body_charge.size() >= this->m_pdata->getNTypes());
    // assert(m_body_diameter.size() >= this->m_pdata->getNTypes());

    if (body_typeid >= this->m_pdata->getNTypes())
        {
        throw std::runtime_error("Error initializing SuspensionFlow: Invalid rigid body type.");
        }

    //if (type.size() != pos.size() || orientation.size() != pos.size())
    if (type.size() != pos.size())
        {
        std::ostringstream error_msg;
        error_msg << "Error initializing SuspensionFlow: Constituent particle lists"
                  << " (position, orientation, type) are of unequal length.";
        throw std::runtime_error(error_msg.str());
        }
    // if (charge.size() && charge.size() != pos.size())
    //     {
    //     std::ostringstream error_msg;
    //     error_msg << "Error initializing SuspensionFlow: Charges are non-empty but of different "
    //               << "length than the positions.";
    //     throw std::runtime_error(error_msg.str());
    //     }
    // if (diameter.size() && diameter.size() != pos.size())
    //     {
    //     std::ostringstream error_msg;
    //     error_msg << "Error initializing SuspensionFlow: Diameters are non-empty but of different "
    //               << "length than the positions.";
    //     throw std::runtime_error(error_msg.str());
    //     }

    bool body_updated = false;

    bool body_len_changed = false;

        // detect if bodies have changed

        {
        ArrayHandle<unsigned int> h_body_type(m_body_types,
                                              access_location::host,
                                              access_mode::read);
        ArrayHandle<Scalar3> h_body_pos(m_body_pos, access_location::host, access_mode::read);
        // ArrayHandle<Scalar4> h_body_orientation(m_body_orientation,
        //                                         access_location::host,
        //                                         access_mode::read);
        ArrayHandle<unsigned int> h_body_len(m_body_len,
                                             access_location::host,
                                             access_mode::readwrite);

        assert(body_typeid < m_body_len.getNumElements());
        if (type.size() != h_body_len.data[body_typeid])
            {
            body_updated = true;

            h_body_len.data[body_typeid] = (unsigned int)type.size();
            body_len_changed = true;
            }
        else
            {
            for (unsigned int i = 0; i < type.size(); ++i)
                {
                // auto body_index = m_body_idx(body_typeid, i);
                // if (type[i] != h_body_type.data[body_index] || pos[i] != h_body_pos.data[body_index]
                //     || orientation[i] != h_body_orientation.data[body_index])
                //     {
                //     body_updated = true;
                //     }
                auto body_index = m_body_idx(body_typeid, i);
                if (type[i] != h_body_type.data[body_index] || pos[i] != h_body_pos.data[body_index])
                    {
                    body_updated = true;
                    }
                }
            }
        }

    if (body_len_changed)
        {
        if (type.size() > m_body_types.getHeight())
            {
            // resize per-type arrays
            m_body_types.resize(this->m_pdata->getNTypes(), type.size());
            m_body_pos.resize(this->m_pdata->getNTypes(), type.size());
            // m_body_orientation.resize(this->m_pdata->getNTypes(), type.size());

            m_body_idx = Index2D((unsigned int)m_body_types.getPitch(),
                                 (unsigned int)m_body_types.getHeight());
            }
        }

    if (body_updated)
        {
            {
            ArrayHandle<unsigned int> h_body_type(m_body_types,
                                                  access_location::host,
                                                  access_mode::readwrite);
            ArrayHandle<Scalar3> h_body_pos(m_body_pos,
                                            access_location::host,
                                            access_mode::readwrite);
            // ArrayHandle<Scalar4> h_body_orientation(m_body_orientation,
            //                                         access_location::host,
            //                                         access_mode::readwrite);

            // m_body_charge[body_typeid].resize(type.size());
            // m_body_diameter[body_typeid].resize(type.size());

            // store body data in GlobalArray
            for (unsigned int i = 0; i < type.size(); ++i)
                {
                h_body_type.data[m_body_idx(body_typeid, i)] = type[i];
                h_body_pos.data[m_body_idx(body_typeid, i)] = pos[i];
                // h_body_orientation.data[m_body_idx(body_typeid, i)] = orientation[i];

                // m_body_charge[body_typeid][i] = charge[i];
                // m_body_diameter[body_typeid][i] = diameter[i];
                }
            }
        m_bodies_changed = true;
        assert(m_d_max_changed.size() > body_typeid);

        // make sure central particle will be communicated
        m_d_max_changed[body_typeid] = true;

        // also update diameter on constituent particles
        for (unsigned int i = 0; i < type.size(); ++i)
            {
            m_d_max_changed[type[i]] = true;
            }
        }
    }

/** Compute the needed body ghost layer width.

    For central particles, the body ghost layer width is the maximum of [d_i + r_ghost_i] where d_i
    is distance of particle i from the center of the body and r_ghost_i is the ghost width for
    particle i determined by cutoff and other interactions.

    The body ghost layer width for constituent particles *should be* the maximum diameter d_i among
    all rigid body types that have this particle as a consituent. However, this must be larger due
    to limitations in the way that individual rigid body particles are indexed relative to the
    aggregates in MolecularForceCompute. In the worst case, for a ghost particle within the
    interaction ghost width r_ghost_i of a boundary, *ALL* other particles in that body must be
    included. The ghost layer width needed to satisfy this condition is the maximum of [2*d_i +
    r_ghost_i], allowing for enough distance to communicate another particle placed at -r_i.
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
Scalar SuspensionFlow<KT_, SET_>::requestBodyGhostLayerWidth(unsigned int type, Scalar* h_r_ghost)
    {
    assert(m_body_len.getNumElements() > type);
    ArrayHandle<unsigned int> h_body_len(m_body_len, access_location::host, access_mode::read);

    if (m_d_max_changed[type])
        {
        m_d_max[type] = Scalar(0.0);
        ArrayHandle<Scalar3> h_body_pos(m_body_pos, access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_body_type(m_body_types,
                                              access_location::host,
                                              access_mode::read);

        if (h_body_len.data[type] != 0)
            {
            // central particles
            for (unsigned int body_particle = 0; body_particle < h_body_len.data[type];
                 body_particle++)
                {
                unsigned int constituent_typeid = h_body_type.data[m_body_idx(type, body_particle)];
                vec3<Scalar> constituent_position(h_body_pos.data[m_body_idx(type, body_particle)]);
                Scalar d = slow::sqrt(dot(constituent_position, constituent_position));
                m_d_max[type] = std::max(m_d_max[type], d + h_r_ghost[constituent_typeid]);
                }
            }
        else
            {
            // constituent particles
            for (unsigned int body_type = 0; body_type < this->m_pdata->getNTypes(); body_type++)
                {
                if (h_body_len.data[body_type] == 0)
                    {
                    continue;
                    }

                for (unsigned int body_particle = 0; body_particle < h_body_len.data[body_type];
                     body_particle++)
                    {
                    unsigned int constituent_typeid
                        = h_body_type.data[m_body_idx(body_type, body_particle)];

                    if (constituent_typeid == type)
                        {
                        vec3<Scalar> constituent_position(
                            h_body_pos.data[m_body_idx(body_type, body_particle)]);
                        Scalar d = slow::sqrt(dot(constituent_position, constituent_position));
                        m_d_max[type] = std::max(m_d_max[type],
                                                 d * Scalar(2.0) + h_r_ghost[constituent_typeid]);
                        }
                    }
                }
            }
        }

    m_d_max_changed[type] = false;
    this->m_exec_conf->msg->notice(7) << "SuspensionFlow: requesting ghost layer for type "
                                << this->m_pdata->getNameByType(type) << ": " << m_d_max[type]
                                << std::endl;

    return m_d_max[type];
    }

template<SmoothingKernelType KT_, StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::validateRigidBodies()
    {
    if (!(m_bodies_changed || m_particles_added_removed))
        {
        return;
        }

    // check validity of rigid body types: no nested rigid bodies
    unsigned int ntypes = this->m_pdata->getNTypes();
    assert(m_body_types.getPitch() >= ntypes);
        {
        ArrayHandle<unsigned int> h_body_type(m_body_types,
                                              access_location::host,
                                              access_mode::read);
        ArrayHandle<unsigned int> h_body_len(m_body_len, access_location::host, access_mode::read);
        for (unsigned int itype = 0; itype < ntypes; ++itype)
            {
            for (unsigned int j = 0; j < h_body_len.data[itype]; ++j)
                {
                assert(h_body_type.data[m_body_idx(itype, j)] <= ntypes);
                if (h_body_len.data[h_body_type.data[m_body_idx(itype, j)]] != 0)
                    {
                    throw std::runtime_error(
                        "Error initializing SuspensionFlow: A rigid body type "
                        "may not contain constituent particles that are also rigid bodies.");
                    }
                }
            }
        }

    SnapshotParticleData<Scalar> snap;

    // take a snapshot on rank 0
    this->m_pdata->takeSnapshot(snap);

    std::vector<unsigned int> aggregate_tag;

    // number of bodies in system
    unsigned int nbodies = 0;

    // number of free particles in the system
    m_n_free_particles_global = 0;

    // Validate the body tags in the system and assign aggregates into aggregate tag
    if (this->m_exec_conf->getRank() == 0)
        {
        // access body data
        ArrayHandle<unsigned int> h_body_len(m_body_len, access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_body_type(m_body_types,
                                              access_location::host,
                                              access_mode::read);

        typedef std::map<unsigned int, unsigned int> map_t;
        // This will count the length of all aggregates (realized rigid bodies) in the system.
        map_t body_particle_count;

        aggregate_tag.resize(snap.size, NO_AGGREGATE);

        // count number of constituent particles to add
        for (unsigned i = 0; i < snap.size; ++i)
            {
            assert(snap.type[i] < ntypes);

            // If a particle is of a type with a non-zero body definition it should be a central
            // particle and the body value should equal its particle tag.
            if (h_body_len.data[snap.type[i]] != 0)
                {
                if (snap.body[i] != i)
                    {
                    throw std::runtime_error(
                        "Error validating rigid bodies: Particles of types defining rigid bodies "
                        "must have a body tag identical to their particle tag to be considered a "
                        "central particle.");
                    }

                // Create a new aggregate count for central particle i.
                body_particle_count.insert(std::make_pair(i, 0));
                aggregate_tag[i] = nbodies++;
                }
            // validate constituent particles. MIN_FLOPPY defines the maximum tag for a particle
            // that is in a rigid body. Tags higher than this can be in a floppy body.
            else if (snap.body[i] < MIN_FLOPPY)
                {
                // check if particle body tag correctly points to the central particle and less than
                // the number of particles in the system. This first check is to ensure that no
                // unallocated memory is attempted to be accessed in the second check.
                if (snap.body[i] >= snap.size || snap.body[snap.body[i]] != snap.body[i])
                    {
                    throw std::runtime_error(
                        "Error validating rigid bodies: Constituent particle body tags must "
                        "be the tag of their central particles.");
                    }

                unsigned int central_particle_index = snap.body[i];
                map_t::iterator it = body_particle_count.find(central_particle_index);
                // If find returns the end of the map, then the central particle for this
                // particular composite particle has not been seen yet. Since we are iterating
                // over the snapshot, this means that the tag for the central particle is higher
                // than the composite particles which is not allowed.
                if (it == body_particle_count.end())
                    {
                    throw std::runtime_error(
                        "Error validating rigid bodies: Central particle must have a lower "
                        "tag than all constituent particles.");
                    }

                unsigned int current_aggregate_size = it->second;
                unsigned int body_type = snap.type[central_particle_index];
                if (current_aggregate_size == h_body_len.data[body_type])
                    {
                    throw std::runtime_error(
                        "Error validating rigid bodies: Too many constituent particles for "
                        "rigid body.");
                    }

                if (h_body_type.data[m_body_idx(body_type, current_aggregate_size)] != snap.type[i])
                    {
                    throw std::runtime_error(
                        "Error validating rigid bodies: Constituent particle types must be "
                        "consistent with the rigid body definition.  Rigid body definition "
                        "is in order of particle tag.");
                    }
                // increase aggregate size by one as particle is validated
                it->second++;
                // Mark consistent particle in aggregate as belonging to its central particle.
                aggregate_tag[i] = aggregate_tag[snap.body[i]];
                }
            else
                {
                m_n_free_particles_global++;
                }
            }
        for (auto it = body_particle_count.begin(); it != body_particle_count.end(); ++it)
            {
            const auto central_particle_tag = it->first;
            const auto aggregate_size = it->second;
            unsigned int central_particle_type = snap.type[central_particle_tag];
            if (aggregate_size != h_body_len.data[central_particle_type])
                {
                std::ostringstream error_msg;
                error_msg << "Error validating rigid bodies: Incomplete rigid body with only "
                          << aggregate_size << " constituent particles "
                          << "instead of " << h_body_len.data[central_particle_type] << " for body "
                          << central_particle_tag;
                throw std::runtime_error(error_msg.str());
                }
            }
        }

#ifdef ENABLE_MPI
    if (this->m_pdata->getDomainDecomposition())
        {
        bcast(aggregate_tag, 0, this->m_exec_conf->getMPICommunicator());
        bcast(nbodies, 0, this->m_exec_conf->getMPICommunicator());
        bcast(m_n_free_particles_global, 0, this->m_exec_conf->getMPICommunicator());
        }
#endif

    // resize Molecular tag member array
    this->m_aggregate_tag.resize(aggregate_tag.size());
        {
        ArrayHandle<unsigned int> h_aggregate_tag(this->m_aggregate_tag,
                                                 access_location::host,
                                                 access_mode::overwrite);
        std::copy(aggregate_tag.begin(), aggregate_tag.end(), h_aggregate_tag.data);
        }

    // store number of aggregates in all ranks
    this->m_n_aggregates_global = nbodies;

    // reset flags
    m_bodies_changed = false;
    m_particles_added_removed = false;
    }

template<SmoothingKernelType KT_, StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::createRigidBodies()
    {
    SnapshotParticleData<Scalar> snap;
    const BoxDim& global_box = this->m_pdata->getGlobalBox();

    // take a snapshot on rank 0
    this->m_pdata->takeSnapshot(snap);
    bool remove_existing_constituents = false;
    unsigned int n_constituent_particles_to_add = 0;
    unsigned int n_free_particles = 0;

    if (this->m_exec_conf->getRank() == 0)
        {
        ArrayHandle<unsigned int> h_body_len(m_body_len, access_location::host, access_mode::read);
        for (unsigned int particle_tag = 0; particle_tag < snap.size; ++particle_tag)
            {
            // Determine whether each particle is rigid or free based on the particle type and the
            // rigid body definition.
            // TODO: Hier eventuell auch noch 1 für Wall Boundarys einführen

            // // Print the value of h_body_len.data for each particle_tag
            // std::cout << "h_body_len.data[" << snap.type[particle_tag] << "] = "
            //           << h_body_len.data[snap.type[particle_tag]] << std::endl;


            if (h_body_len.data[snap.type[particle_tag]] == 0)
                {
                n_free_particles++;
                }
            else
                {
                // Determine whether we need to remove existing constituent particles. These will
                // be recreated below.
                if (snap.body[particle_tag] != particle_tag)
                    {
                    remove_existing_constituents = true;
                    }

                // Increase the number of particles we need to add by the number of constituent
                // particles this rigid body center has based on its type.
                n_constituent_particles_to_add += h_body_len.data[snap.type[particle_tag]];
                }
            }
        }

#ifdef ENABLE_MPI
    if (this->m_pdata->getDomainDecomposition())
        {
        bcast(remove_existing_constituents, 0, this->m_exec_conf->getMPICommunicator());
        bcast(n_free_particles, 0, this->m_exec_conf->getMPICommunicator());
        }
#endif

    if (remove_existing_constituents)
        {
        this->m_exec_conf->msg->notice(7)
            << "SuspensionFlow reinitialize particle data without rigid bodies" << std::endl;
        this->m_pdata->initializeFromSnapshot(snap, true);
        this->m_pdata->takeSnapshot(snap);
        }

    std::vector<unsigned int> aggregate_tag;

    unsigned int n_central_particles = snap.size - n_free_particles;

    if (this->m_exec_conf->getRank() == 0)
        {
        unsigned int initial_snapshot_size = snap.size;
        snap.insert(snap.size, n_constituent_particles_to_add);

        ArrayHandle<unsigned int> h_body_len(m_body_len, access_location::host, access_mode::read);
        ArrayHandle<Scalar3> h_body_pos(m_body_pos, access_location::host, access_mode::read);
        // ArrayHandle<Scalar4> h_body_orientation(m_body_orientation,
        //                                         access_location::host,
        //                                         access_mode::read);
        ArrayHandle<unsigned int> h_body_type(m_body_types,
                                              access_location::host,
                                              access_mode::read);
        aggregate_tag.resize(snap.size, NO_AGGREGATE);

        unsigned int constituent_particle_tag = initial_snapshot_size;
        for (unsigned int particle_tag = 0; particle_tag < initial_snapshot_size; ++particle_tag)
            {
            assert(snap.type[particle_tag] < this->m_pdata->getNTypes());

            // If the length of the body definition is zero it must be a free particle because all
            // constituent particles have been removed.
            if (h_body_len.data[snap.type[particle_tag]] == 0)
                {
                snap.body[particle_tag] = NO_BODY;
                continue;
                }
            snap.body[particle_tag] = particle_tag;
            aggregate_tag[particle_tag] = particle_tag;

            unsigned int body_type = snap.type[particle_tag];
            unsigned int n_body_particles = h_body_len.data[body_type];

            std::cout << "body_type" << body_type << std::endl;
            std::cout << "n_body_particles" << n_body_particles << std::endl;

            for (unsigned int current_body_index = 0; current_body_index < n_body_particles;
                 ++current_body_index)
                {
                size_t body_idx = m_body_idx(body_type, current_body_index);

                // Update constituent particle snapshot properties from default.
                snap.type[constituent_particle_tag] = h_body_type.data[body_idx];
                snap.body[constituent_particle_tag] = particle_tag;
                // snap.charge[constituent_particle_tag]
                //     = m_body_charge[body_type][current_body_index];
                // snap.diameter[constituent_particle_tag]
                //     = m_body_diameter[body_type][current_body_index];

                // Set position and orientation of constituents
                vec3<Scalar> body_position(snap.pos[particle_tag]);
                //quat<Scalar> body_orientation(snap.orientation[particle_tag]);
                vec3<Scalar> local_position(h_body_pos.data[body_idx]);
                //quat<Scalar> local_orientation(h_body_orientation.data[body_idx]);

                std::cout << "local_position.x: " << local_position.x << "local_position.y: " << local_position.y << "local_position.z: " << local_position.z << std::endl;


                // vec3<Scalar> constituent_position = body_position + rotate(body_orientation, local_position);
                vec3<Scalar> constituent_position = body_position + local_position;

                //quat<Scalar> constituent_orientation = body_orientation * local_orientation;

                //vec3<Scalar> constituent_position(snap.pos[particle_tag]);

                std::cout << "constituent_position.x: " << constituent_position.x << "constituent_position.y: " << constituent_position.y << "constituent_position.z: " << constituent_position.z << std::endl;

                snap.pos[constituent_particle_tag] = constituent_position;
                snap.image[constituent_particle_tag] = snap.image[particle_tag];
                //snap.orientation[constituent_particle_tag] = constituent_orientation;

                // wrap back into the box
                global_box.wrap(snap.pos[constituent_particle_tag],
                                snap.image[constituent_particle_tag]);

                // Since the central particle tags here will be [0, n_central_particles), we know
                // that the aggregate number will be the same as the central particle tag.
                aggregate_tag[constituent_particle_tag] = particle_tag;

                ++constituent_particle_tag;
                }
            }
        }

    // Keep rigid bodies this time when initializing.
    this->m_pdata->initializeFromSnapshot(snap, false);

#ifdef ENABLE_MPI
    if (this->m_pdata->getDomainDecomposition())
        {
        bcast(aggregate_tag, 0, this->m_exec_conf->getMPICommunicator());
        }
#endif

    this->m_aggregate_tag.resize(aggregate_tag.size());
        {
        // store global aggregate information in GlobalArray
        ArrayHandle<unsigned int> h_aggregate_tag(this->m_aggregate_tag,
                                                 access_location::host,
                                                 access_mode::overwrite);
        std::copy(aggregate_tag.begin(), aggregate_tag.end(), h_aggregate_tag.data);
        }
    this->m_n_aggregates_global = n_central_particles;
    m_n_free_particles_global = n_free_particles;

    std::cout << "m_n_aggregates_global" << this->m_n_aggregates_global << std::endl;
    std::cout << std::endl;
    std::cout << "m_n_free_particles_global" << m_n_free_particles_global << std::endl;
    std::cout << std::endl;

    m_bodies_changed = false;
    m_particles_added_removed = false;
    }

// #ifdef ENABLE_MPI
// /*! \param timestep Current time step
//  */
// template<SmoothingKernelType KT_, StateEquationType SET_>
// CommFlags SuspensionFlow<KT_, SET_>::getRequestedCommFlags(uint64_t timestep)
//     {
//     CommFlags flags = CommFlags(0);

//     // // request orientations
//     // flags[comm_flag::orientation] = 1;

//     // // only communicate net virial if needed
//     // PDataFlags pdata_flags = this->m_pdata->getFlags();
//     // if (pdata_flags[pdata_flag::pressure_tensor])
//     //     {
//     //     flags[comm_flag::net_virial] = 1;
//     //     }

//     // request body ids
//     flags[comm_flag::body] = 1;

//     // we need central particle images
//     flags[comm_flag::image] = 1;

//     flags |= SPHBaseClassConstraint::getRequestedCommFlags(timestep);

//     return flags;
//     }
// #endif

/* Set position and velocity of constituent particles in rigid bodies in the 1st or second half of
 * integration on the CPU based on the body center of mass and particle relative position in each
 * body frame.
 */

template<SmoothingKernelType KT_, StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::updateCompositeParticles(uint64_t timestep)
    {
    // If no rigid bodies exist return early. This also prevents accessing arrays assuming that this
    // is non-zero.
    if (this->m_n_aggregates_global == 0)
        {
        return;
        }

    // access aggregate order (this needs to be on top because of ArrayHandle scope) and its
    // pervasive use across this function.
    ArrayHandle<unsigned int> h_aggregate_order(this->getAggregateOrder(),
                                               access_location::host,
                                               access_mode::read);
    ArrayHandle<unsigned int> h_aggregate_len(this->getAggregateLengths(),
                                             access_location::host,
                                             access_mode::read);
    ArrayHandle<unsigned int> h_aggregate_idx(this->getAggregateIndex(),
                                             access_location::host,
                                             access_mode::read);

    // access the particle data arrays
    ArrayHandle<Scalar4> h_postype(this->m_pdata->getPositions(),
                                   access_location::host,
                                   access_mode::readwrite);
    // ArrayHandle<Scalar4> h_orientation(this->m_pdata->getOrientationArray(),
    //                                    access_location::host,
    //                                    access_mode::readwrite);
    ArrayHandle<int3> h_image(this->m_pdata->getImages(), access_location::host, access_mode::readwrite);

    ArrayHandle<unsigned int> h_body(this->m_pdata->getBodies(),
                                     access_location::host,
                                     access_mode::read);
    ArrayHandle<unsigned int> h_rtag(this->m_pdata->getRTags(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_tag(this->m_pdata->getTags(), access_location::host, access_mode::read);

    // access body positions and orientations
    ArrayHandle<Scalar3> h_body_pos(m_body_pos, access_location::host, access_mode::read);
    // ArrayHandle<Scalar4> h_body_orientation(m_body_orientation,
    //                                         access_location::host,
    //                                         access_mode::read);
    ArrayHandle<unsigned int> h_body_len(m_body_len, access_location::host, access_mode::read);

    const BoxDim& box = this->m_pdata->getBox();
    const BoxDim& global_box = this->m_pdata->getGlobalBox();

    // we need to update both local and ghost particles
    unsigned int n_particles_local = this->m_pdata->getN() + this->m_pdata->getNGhosts();
    for (unsigned int particle_index = 0; particle_index < n_particles_local; particle_index++)
        {
        unsigned int central_tag = h_body.data[particle_index];

        // Do nothing with floppy bodies, since we don't need to update their positions or
        // orientations here.
        if (central_tag >= MIN_FLOPPY)
            {
            continue;
            }

        // body tag equals tag for central particle
        assert(central_tag <= this->m_pdata->getMaximumTag());
        unsigned int central_idx = h_rtag.data[central_tag];

        // If this is a rigid body center continue, since we do not need to update its position or
        // orientation (the integrator methods do this).
        if (particle_index == central_idx)
            {
            continue;
            }

        // If the central particle is not local, then we cannot update the position and orientation
        // of this particle. Ideally, this would perform an error check. However, that is not
        // feasible as SuspensionFlow does not have knowledge of which ghost particles are within
        // the interaction ghost width (and need therefore need to be updated) vs those that are
        // communicated to make bodies whole.
        if (central_idx == NOT_LOCAL)
            {
            continue;
            }

        // central particle position and orientation
        assert(central_idx <= this->m_pdata->getN() + this->m_pdata->getNGhosts());

        Scalar4 postype = h_postype.data[central_idx];
        vec3<Scalar> pos(postype);
        //quat<Scalar> orientation(h_orientation.data[central_idx]);

        // body type
        unsigned int type = __scalar_as_int(postype.w);

        unsigned int body_len = h_body_len.data[type];
        unsigned int agg_idx = h_aggregate_idx.data[particle_index];
        // Checks if the number of local particle in a aggregate denoted by
        // h_aggregate_len.data[particle_index] is equal to the number of particles in the rigid body
        // definition `body_len`. As above, this error check *should* be performed for all local and
        // ghost particles within the interaction ghost width. However, that check is not feasible
        // here. At least catch this error for particles local to this rank.
        if (body_len != h_aggregate_len.data[agg_idx] - 1)
            {
            if (particle_index < this->m_pdata->getN())
                {
                // if the aggregate is incomplete and has local members, this is an error
                std::ostringstream error_msg;
                error_msg << "Error while updating constituent particles:"
                          << "Composite particle with body tag " << central_tag << " incomplete: "
                          << "body_len=" << body_len
                          << ", aggregate_len=" << h_aggregate_len.data[agg_idx] - 1;
                throw std::runtime_error(error_msg.str());
                }

            // otherwise we must ignore it
            continue;
            }

        int3 img = h_image.data[central_idx];

        // fetch relative index in body from aggregate list
        assert(h_aggregate_order.data[particle_index] > 0);
        unsigned int idx_in_body = h_aggregate_order.data[particle_index] - 1;

        vec3<Scalar> local_pos(h_body_pos.data[m_body_idx(type, idx_in_body)]);
        //vec3<Scalar> dr_space = rotate(orientation, local_pos);

        // TODO: Eventuell muss Orientation in Aux2 Feld geschrieben werden 
        // update position and orientation
        vec3<Scalar> updated_pos(pos);
        //quat<Scalar> local_orientation(h_body_orientation.data[m_body_idx(type, idx_in_body)]);

        //updated_pos += dr_space;
        //quat<Scalar> updated_orientation = orientation * local_orientation;

        // this runs before the ForceComputes,
        // wrap into box, allowing rigid bodies to span multiple images
        int3 imgi = box.getImage(vec_to_scalar3(updated_pos));
        int3 negimgi = make_int3(-imgi.x, -imgi.y, -imgi.z);
        updated_pos = global_box.shift(updated_pos, negimgi);

        h_postype.data[particle_index] = make_scalar4(updated_pos.x,
                                                      updated_pos.y,
                                                      updated_pos.z,
                                                      h_postype.data[particle_index].w);
        // h_orientation.data[particle_index] = quat_to_scalar4(updated_orientation);
        h_image.data[particle_index] = img + imgi;
        }
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
            bool isaggregate = checkfluid1(h_type_property_map.data, h_pos.data[k].w);

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
            // TODO: Hier checken ob Berechnung mit Soliddichte oder Noslip-Dichte für Aggregates
            Scalar rhoj = h_density.data[k];
            Scalar Vj   = mj / rhoj;

            // Read particle k pressure
            // TODO: Hier checken ob Berechnung mit NoslipDruck --> Wahrscheinlich ja
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
                temp0 = -(Pi*Vi*Vi+Pj*Vj*Vj); // TODO: Aus Twophaseflow alt, keine Dichte enthalten
                //temp0 = -(Vi*Vi+Vj*Vj)*((rhoj*Pi+rhoi*Pj)/(rhoi+rhoj)); 
            }
            else if ( m_density_method == DENSITYCONTINUITY) 
            { 
                temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);
            }


            // Optionally add artificial viscosity
            // Monaghan (1983) J. Comput. Phys. 52 (2) 374–389
            if ( m_artificial_viscosity && !issolid )
                {
                Scalar dotdvdx = dot(dv,dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar muij    = meanh*dotdvdx/(rsq+epssqr);
                    Scalar meanrho = Scalar(0.5)*(rhoi+rhoj);
                    temp0 += mi*mj*(m_avalpha*m_c*muij+m_avbeta*muij*muij)/meanrho;
                    }
                }

            // Add contribution to fluid particle
            h_force.data[i].x += temp0*dwdr_r*dx.x;
            h_force.data[i].y += temp0*dwdr_r*dx.y;
            h_force.data[i].z += temp0*dwdr_r*dx.z;

            // Add contribution to solid particle
            if ( isaggregate && m_compute_solid_forces )
                {
                h_force.data[k].x -= (mj/mi)*temp0*dwdr_r*dx.x;
                h_force.data[k].y -= (mj/mi)*temp0*dwdr_r*dx.y;
                h_force.data[k].z -= (mj/mi)*temp0*dwdr_r*dx.z;
                }

            // Evaluate viscous interaction forces
            temp0 = m_mu * (Vi*Vi+Vj*Vj) * dwdr_r;
            h_force.data[i].x  += temp0*dv.x;
            h_force.data[i].y  += temp0*dv.y;
            h_force.data[i].z  += temp0*dv.z;

            // Add contribution to solid particle
            if ( isaggregate && m_compute_solid_forces )
                {
                h_force.data[k].x -= (mj/mi)*temp0*dv.x;
                h_force.data[k].y -= (mj/mi)*temp0*dv.y;
                h_force.data[k].z -= (mj/mi)*temp0*dv.z;
                }

            // Evaluate rate of change of density if CONTINUITY approach is used
            if ( m_density_method == DENSITYCONTINUITY )
                {
                if ( issolid )
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
                if ( !issolid && !isaggregate && m_density_diffusion )
                    h_ratedpe.data[i].x -= (Scalar(2)*m_ddiff*meanh*m_c*mj*(rhoi/rhoj-Scalar(1))*dot(dx,dwdr_r*dx))/(rsq+epssqr);
                }

            } // Closing Neighbor Loop

        } // Closing Fluid Particle Loop

    m_timestep_list[5] = max_vel;
    // Add volumetric force (gravity)
    this->applyBodyForce(timestep, m_fluidgroup);
    if ( m_compute_solid_forces )
        this->applyBodyForce(timestep, m_aggregategroup);

    }

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlow<KT_, SET_>::rigidforcecomputation(uint64_t timestep)
    {
    // If no rigid bodies exist return early. This also prevents accessing arrays assuming that this
    // is non-zero.
    if (this->m_n_aggregates_global == 0)
        {
        return;
        }

    std::cout << "Rigid Force Computation" << std::endl;

    // access local aggregate data
    // need to move this on top because of scoping issues
    Index2D aggregate_indexer = this->getAggregateIndexer();
    unsigned int nmol = aggregate_indexer.getH();

    ArrayHandle<unsigned int> h_aggregate_length(this->getAggregateLengths(),
                                                access_location::host,
                                                access_mode::read);
    ArrayHandle<unsigned int> h_aggregate_list(this->getAggregateList(),
                                              access_location::host,
                                              access_mode::read);

    // access particle data
    ArrayHandle<unsigned int> h_body(this->m_pdata->getBodies(),
                                     access_location::host,
                                     access_mode::read);
    ArrayHandle<unsigned int> h_rtag(this->m_pdata->getRTags(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_tag(this->m_pdata->getTags(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_postype(this->m_pdata->getPositions(),
                                   access_location::host,
                                   access_mode::read);
    // ArrayHandle<Scalar4> h_orientation(this->m_pdata->getOrientationArray(),
    //                                    access_location::host,
    //                                    access_mode::read);

    // access net force and torque acting on constituent particles
    ArrayHandle<Scalar4> h_net_force(this->m_pdata->getNetForce(),
                                     access_location::host,
                                     access_mode::readwrite);
    // TODO: Torque in Aux
    // ArrayHandle<Scalar4> h_net_torque(this->m_pdata->getNetTorqueArray(),
    //                                   access_location::host,
    //                                   access_mode::readwrite);


    // access the force and torque array for the central particle
    ArrayHandle<Scalar4> h_force(this->m_force, access_location::host, access_mode::overwrite);
    // ArrayHandle<Scalar4> h_torque(this->m_torque, access_location::host, access_mode::overwrite);

    // access rigid body definition
    ArrayHandle<Scalar3> h_body_pos(this->m_body_pos, access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_body_len(this->m_body_len, access_location::host, access_mode::read);

    // reset constraint forces and torques
    memset(h_force.data, 0, sizeof(Scalar4) * this->m_pdata->getN());
    //memset(h_torque.data, 0, sizeof(Scalar4) * this->m_pdata->getN());

    unsigned int n_particles_local = this->m_pdata->getN() + this->m_pdata->getNGhosts();

    PDataFlags flags = this->m_pdata->getFlags();

    // loop over all aggregates, also incomplete ones
    for (unsigned int ibody = 0; ibody < nmol; ibody++)
        {
        // get central particle tag from first particle in aggregate
        assert(h_aggregate_length.data[ibody] > 0);
        unsigned int first_idx = h_aggregate_list.data[aggregate_indexer(0, ibody)];

        assert(first_idx < this->m_pdata->getN() + this->m_pdata->getNGhosts());
        unsigned int central_tag = h_body.data[first_idx];

        assert(central_tag <= this->m_pdata->getMaximumTag());
        unsigned int central_idx = h_rtag.data[central_tag];

        if (central_idx >= n_particles_local)
            continue;

        // the central particle must be present
        assert(central_tag == h_tag.data[first_idx]);

        // central particle position and orientation
        Scalar4 postype = h_postype.data[central_idx];
        // quat<Scalar> orientation(U.data[central_idx]);

        // body type
        unsigned int type = __scalar_as_int(postype.w);

        // sum up forces and torques from constituent particles
        for (unsigned int constituent_index = 0; constituent_index < h_aggregate_length.data[ibody];
             ++constituent_index)
            {
            unsigned int idxj = h_aggregate_list.data[aggregate_indexer(constituent_index, ibody)];
            assert(idxj < this->m_pdata->getN() + this->m_pdata->getNGhosts());

            assert(idxj == central_idx || constituent_index > 0);
            if (idxj == central_idx)
                continue;

            // force and torque on particle
            Scalar4 net_force = h_net_force.data[idxj];
            // Scalar4 net_torque = h_net_torque.data[idxj];
            vec3<Scalar> f(net_force);

            // zero net energy on constituent particles to avoid double counting
            // also zero net force and torque for consistency
            h_net_force.data[idxj] = make_scalar4(0.0, 0.0, 0.0, 0.0);
            //h_net_torque.data[idxj] = make_scalar4(0.0, 0.0, 0.0, 0.0);

            // only add forces for local central particles
            if (central_idx < this->m_pdata->getN())
                {
                // if the central particle is local, the aggregate should be complete
                if (h_aggregate_length.data[ibody] != h_body_len.data[type] + 1)
                    {
                    std::ostringstream error_msg;
                    error_msg << "Composite particle with body tag " << central_tag
                              << " is incomplete.";
                    throw std::runtime_error(error_msg.str());
                    }

                // sum up center of mass force
                h_force.data[central_idx].x += f.x;
                h_force.data[central_idx].y += f.y;
                h_force.data[central_idx].z += f.z;

                // sum up energy
                h_force.data[central_idx].w += net_force.w;

                // fetch relative position from rigid body definition
                vec3<Scalar> dr(h_body_pos.data[m_body_idx(type, constituent_index - 1)]);

                // rotate into space frame
                //vec3<Scalar> dr_space = rotate(orientation, dr);

                // // torque = r x f
                // vec3<Scalar> delta_torque(cross(dr_space, f));
                // h_torque.data[central_idx].x += delta_torque.x;
                // h_torque.data[central_idx].y += delta_torque.y;
                // h_torque.data[central_idx].z += delta_torque.z;

                // /* from previous rigid body implementation: Access Torque elements from a single
                //    particle. Right now I will am assuming that the particle and rigid body reference
                //    frames are the same. Probably have to rotate first.
                //  */
                // h_torque.data[central_idx].x += net_torque.x;
                // h_torque.data[central_idx].y += net_torque.y;
                // h_torque.data[central_idx].z += net_torque.z;

                }

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
    compute_noslip(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    update_ghost_aux1(timestep);
#endif

#ifdef ENABLE_MPI
    // Update ghost particles
    update_ghost_aux2(timestep);
#endif

    rigidforcecomputation(timestep);


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
    pybind11::class_<SuspensionFlow<KT_, SET_>, SPHBaseClassConstraint<KT_, SET_> , std::shared_ptr<SuspensionFlow<KT_, SET_>>>(m, name.c_str()) 
        .def(pybind11::init< std::shared_ptr<SystemDefinition>,
                             std::shared_ptr<SmoothingKernel<KT_> >,
                             std::shared_ptr<StateEquation<SET_> >,
                             std::shared_ptr<nsearch::NeighborList>,
                             std::shared_ptr<ParticleGroup>,
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
        .def("activateArtificialViscosity", &SuspensionFlow<KT_, SET_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity", &SuspensionFlow<KT_, SET_>::deactivateArtificialViscosity)
        .def("activateDensityDiffusion", &SuspensionFlow<KT_, SET_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion", &SuspensionFlow<KT_, SET_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &SuspensionFlow<KT_, SET_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization", &SuspensionFlow<KT_, SET_>::deactivateShepardRenormalization)
        .def("setAcceleration", &SPHBaseClassConstraint<KT_, SET_>::setAcceleration)
        .def("setRCut", &SuspensionFlow<KT_, SET_>::setRCutPython)
        .def("setBody", &SuspensionFlow<KT_, SET_>::setBody)
        .def("getBody", &SuspensionFlow<KT_, SET_>::getBody)
        .def("validateRigidBodies", &SuspensionFlow<KT_, SET_>::validateRigidBodies)
        .def("createRigidBodies", &SuspensionFlow<KT_, SET_>::createRigidBodies)
        .def("updateCompositeParticles", &SuspensionFlow<KT_, SET_>::updateCompositeParticles);
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

    template void export_SuspensionFlow<wendlandc2, linear>(pybind11::module& m, std::string name = "SuspensionF_WC2_L");
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

/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SinglePhaseFlow.h"

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
template<SmoothingKernelType KT_,StateEquationType SET_>
SinglePhaseFlow<KT_, SET_>::SinglePhaseFlow(std::shared_ptr<SystemDefinition> sysdef,
                                 std::shared_ptr<SmoothingKernel<KT_> > skernel,
                                 std::shared_ptr<StateEquation<SET_> > equationofstate,
                                 std::shared_ptr<nsearch::NeighborList> nlist,
                                 std::shared_ptr<ParticleGroup> fluidgroup,
                                 std::shared_ptr<ParticleGroup> solidgroup,
                                 DensityMethod mdensitymethod,
                                 ViscosityMethod mviscositymethod)
    : SPHBaseClass<KT_, SET_>(sysdef,skernel,equationofstate,nlist), m_fluidgroup(fluidgroup), m_solidgroup(solidgroup), m_typpair_idx(this->m_pdata->getNTypes())
      {
        this->m_exec_conf->msg->notice(5) << "Constructing SinglePhaseFlow" << std::endl;

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
        m_log_computed_last_timestep = -1;

        // Sanity checks
        assert(this->m_pdata);
        assert(this->m_nlist);
        assert(this->m_skernel);
        assert(this->m_eos);

        // Contruct type vectors
        this->constructTypeVectors(fluidgroup,&m_fluidtypes);
        this->constructTypeVectors(solidgroup,&m_solidtypes);

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
        }

        // Set simulations methods
        m_density_method = mdensitymethod;
        m_viscosity_method = mviscositymethod;

        // Get necessary variables from kernel and EOS classes
        m_rho0  = equationofstate->getRestDensity();
        m_c     = equationofstate->getSpeedOfSound();
        m_kappa = skernel->getKernelKappa();

        // Allocate space for log quantities
        // GPUArray< Scalar > properties(singlephaseflow_logger_index::num_quantities, this->m_exec_conf);
        // m_properties.swap(properties);

        m_r_cut_nlist = std::make_shared<GlobalArray<Scalar>>(m_typpair_idx.getNumElements(), this->m_exec_conf);
        this->m_nlist->addRCutMatrix(m_r_cut_nlist);

//         // Define log quantities
//         m_logname_list.push_back(std::string("sum_fluid_velocity_x"));
//         m_logname_list.push_back(std::string("sum_fluid_velocity_y"));
//         m_logname_list.push_back(std::string("sum_fluid_velocity_z"));
//         m_logname_list.push_back(std::string("total_fluid_particles"));
//         m_logname_list.push_back(std::string("kinetic_energy"));
//         // m_logname_list.push_back(std::string("dt_adapt"));

// #ifdef ENABLE_MPI
//         m_properties_reduced = false;
// #endif
      }

/*! Destructor
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
SinglePhaseFlow<KT_, SET_>::~SinglePhaseFlow()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying SinglePhaseFlow" << std::endl;
    }

/*! Returns provided Log Quantities to Logger
*/
// template<SmoothingKernelType KT_,StateEquationType SET_>
// std::vector< std::string > SinglePhaseFlow<KT_, SET_>::getProvidedLogQuantities()
//     {
//     return m_logname_list;
// }

/*! Returns provided timestep Quantities to Compute
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
std::vector<double> SinglePhaseFlow<KT_, SET_>::getProvidedTimestepQuantities(uint64_t timestep)
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
void  SinglePhaseFlow<KT_, SET_>::activateShepardRenormalization(unsigned int shepardfreq){
            if (shepardfreq <= 0)
                {
                this->m_exec_conf->msg->error() << "sph.models.SinglePhaseFlow: Shepard density reinitialization period has to be a positive real number" << std::endl;
                throw std::runtime_error("Error initializing SinglePhaseFlow.");
            }
            m_shepard_renormalization = true;
            m_shepardfreq = shepardfreq;
            }


/*! Return requested Log quantity
*/
// template<SmoothingKernelType KT_,StateEquationType SET_>
// Scalar SinglePhaseFlow<KT_, SET_>::getLogValue(const std::string& quantity, uint64_t timestep)
//     {
//     if ( m_log_computed_last_timestep != timestep )
//         {
//         // Compute log quantities
//         computeProperties();
//         m_log_computed_last_timestep = timestep;
//         }

//     if (quantity == m_logname_list[0])
//         {
//         return GetSumFluidXVelocity();
//         }
//     else if (quantity == m_logname_list[1])
//         {
//         return GetSumFluidYVelocity();
//         }
//     else if (quantity == m_logname_list[2])
//         {
//         return GetSumFluidZVelocity();
//         }
//     else if (quantity == m_logname_list[3])
//         {
//         return GetFluidParticleNum();
//         }
//     else if (quantity == m_logname_list[4])
//         {
//         return GetKineticEnergy();
//         }
//     // else if (quantity == m_logname_list[4])
//     //     {
//     //     return GetAdaptTimestep();
//     //     }
//     else
//         {
//         this->m_exec_conf->msg->error() << "sph.SinglePhaseFlow: " << quantity << " is not a valid log quantity" << std::endl;
//         throw std::runtime_error("Error getting log value");
//         }
//     }

/*! Computes all log quantities
*/
// template<SmoothingKernelType KT_,StateEquationType SET_>
// void SinglePhaseFlow<KT_, SET_>::computeProperties()
//     {
//     // if (this->m_prof)
//     //     this->m_prof->push("SinglePhaseFlowlogger");

//     ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
//     ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);

//     double fluid_vel_x_sum  = 0.0;
//     double fluid_vel_y_sum  = 0.0;
//     double fluid_vel_z_sum  = 0.0;
//     double fluid_prtl = 0;
//     double kinetic_energy = 0.0;
//     // double adaptive_tstep = 0.0;
    

//     // for each fluid particle
//     unsigned int group_size = m_fluidgroup->getNumMembers();
//     for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
//         {
//         // Read particle index
//         unsigned int i = m_fluidgroup->getMemberIndex(group_idx);
//         // Count fluid particles
//         fluid_prtl += Scalar(1);
//         // Sum velocities
//         fluid_vel_x_sum += h_velocity.data[i].x;
//         fluid_vel_y_sum += h_velocity.data[i].y;
//         fluid_vel_z_sum += h_velocity.data[i].z;
//         kinetic_energy  += abs(0.5*h_velocity.data[i].w*sqrt(pow(h_velocity.data[i].x,2)+pow(h_velocity.data[i].y,2)+pow(h_velocity.data[i].z,2)));
//     }

//     ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::overwrite);
//     h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_x]  = Scalar(fluid_vel_x_sum);
//     h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_y]  = Scalar(fluid_vel_y_sum);
//     h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_z]  = Scalar(fluid_vel_z_sum);
//     h_properties.data[singlephaseflow_logger_index::total_fluid_particles] = Scalar(fluid_prtl);
//     h_properties.data[singlephaseflow_logger_index::kinetic_energy]        = Scalar(kinetic_energy);
//     // h_properties.data[singlephaseflow_logger_index::dt_adapt] = Scalar(adaptive_tstep);

// #ifdef ENABLE_MPI
//     this->m_properties_reduced = !this->m_pdata->getDomainDecomposition();
// #endif

//     // if (this->m_prof)
//     //     this->m_prof->pop();
//     }

// #ifdef ENABLE_MPI

// template<SmoothingKernelType KT_,StateEquationType SET_>
// void SinglePhaseFlow<KT_, SET_>::reduceProperties()
//     {
//     if (m_properties_reduced) return;

//     // reduce properties
//     ArrayHandle<Scalar> h_properties(m_properties, access_location::host, access_mode::readwrite);
//     MPI_Allreduce(MPI_IN_PLACE, h_properties.data, singlephaseflow_logger_index::num_quantities, MPI_HOOMD_SCALAR,
//             MPI_SUM, this->m_exec_conf->getMPICommunicator());

//     m_properties_reduced = true;
//     }
// #endif

/*! \post Set model parameters
 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlow<KT_, SET_>::setParams(Scalar mu)
    {
    m_mu   = mu;
    if (m_mu <= 0)
         {
         this->m_exec_conf->msg->error() << "sph.models.SinglePhaseFlow: Dynamic viscosity has to be a positive real number" << std::endl;
         throw std::runtime_error("Error initializing SinglePhaseFlow.");
         }

    m_params_set = true;
    }

/*! Communicate (update) ghost particle fields
 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlow<KT_, SET_>::update_ghost_dpe(uint64_t timestep)
    {
#ifdef ENABLE_MPI
    if (this->m_comm)
        {
        CommFlags flags(0);
        flags[comm_flag::tag] = 0;
        flags[comm_flag::position] = 0;
        flags[comm_flag::velocity] = 0;
        flags[comm_flag::dpe] = 1;
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
void SinglePhaseFlow<KT_, SET_>::update_ghost_aux1(uint64_t timestep)
    {
#ifdef ENABLE_MPI
    if (this->m_comm)
        {
        CommFlags flags(0);
        flags[comm_flag::tag] = 0;
        flags[comm_flag::position] = 0;
        flags[comm_flag::velocity] = 0;
        flags[comm_flag::dpe] = 1;
        flags[comm_flag::auxiliary1] = 1;
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
void SinglePhaseFlow<KT_, SET_>::validateTypes(unsigned int typ1,
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
void SinglePhaseFlow<KT_, SET_>::setRcut(unsigned int typ1, unsigned int typ2, Scalar rcut)
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
void SinglePhaseFlow<KT_, SET_>::setRCutPython(pybind11::tuple types, Scalar r_cut)
    {
    auto typ1 = this->m_pdata->getTypeByName(types[0].cast<std::string>());
    auto typ2 = this->m_pdata->getTypeByName(types[1].cast<std::string>());
    setRcut(typ1, typ2, r_cut);
    }

/*! Perform number density computation
 * This method computes and stores
     - the density based on a real mass density ( rho_i = m_i * n_neighbours / volume_sphere_of_influnece ) for fluid particles
   in the x-component of the h_dpe Array.
   It overestimates density, but is superfast in comparison to compute_ndensity.
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlow<KT_, SET_>::compute_particlenumberdensity(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Particle Number Density" << std::endl;

    // if (this->m_prof)
    //     this->m_prof->push("SinglePhaseFlowNDensity");

    // Grab handles for particle data
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    // ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    // ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    unsigned int size;
    long unsigned int myHead;

    // Volume of sphere with radius m_rcut
    double sphere_vol = Scalar(4.0)/Scalar(3.0) * Scalar(3.14159265359) * pow(m_rcut, 3);

    // Particle loop
    for (unsigned int i = 0; i < this->m_pdata->getN(); i++)
        {

        // Determine particle i type
        bool i_issolid = checksolid(h_type_property_map.data, h_pos.data[i].w);

        // Skip neighbor loop if this solid solid particle does not have fluid neighbors
        bool solid_w_fluid_neigh = false;
        if ( i_issolid )
            {
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
            }

        if ( i_issolid && !(solid_w_fluid_neigh) )
            {
            h_dpe.data[i].x = this->m_rho0;
            continue;
            }

        // All of the neighbors of this particle
        unsigned int size = (unsigned int)h_n_neigh.data[i];

        h_dpe.data[i].x = size * h_velocity.data[i].w / sphere_vol;
        std::cout << "Mass: " << h_velocity.data[i].w << " size " << size << " sphere_vol " << sphere_vol << " density " << h_dpe.data[i].x << std::endl;


        } // End of particle loop
    }




/*! Perform number density computation
 * This method computes and stores
     - the number density based mass density ( rho_i = m_i * \sum w_ij ) for fluid particles
       if the SUMMATION approach is being used.
     - the zeroth order normalization constant for solid particles
   in the x-component of the h_dpe Array.
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlow<KT_, SET_>::compute_ndensity(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Number Density" << std::endl;

    // if (this->m_prof)
    //     this->m_prof->push("SinglePhaseFlowNDensity");

    // Grab handles for particle data
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    // ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Precompute self-density for homogeneous smoothing lengths
    Scalar w0 = this->m_skernel->w0(m_ch);

    unsigned int size;
    long unsigned int myHead;

    // Particle loop
    for (unsigned int i = 0; i < this->m_pdata->getN(); i++)
        {
        // Access the particle's position
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;

        // Determine particle i type
        bool i_issolid = checksolid(h_type_property_map.data, h_pos.data[i].w);
        bool i_isfluid = checkfluid(h_type_property_map.data, h_pos.data[i].w);

        // Do not compute number density based mass density for fluid particle
        // if anything other that DensityMethod SUMMATION is used.
        if ( this->m_density_method != DENSITYSUMMATION && !(i_issolid) )
            {
            continue;
            }

        // Initialize number density with self density of kernel
        if ( i_issolid )
            h_dpe.data[i].x = 0;
        else
            h_dpe.data[i].x = this->m_const_slength ? w0 : this->m_skernel->w0(h_h.data[i]);

        // Skip neighbor loop if this solid solid particle does not have fluid neighbors
        bool solid_w_fluid_neigh = false;
        if ( i_issolid )
            {
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
            }

        if ( i_issolid && !(solid_w_fluid_neigh) )
            {
            h_dpe.data[i].x = this->m_rho0;
            continue;
            }

        // Loop over all of the neighbors of this particle
        size_t myHead = h_head_list.data[i];
        unsigned int size = (unsigned int)h_n_neigh.data[i];

        for (unsigned int j = 0; j < size; j++)
            {   
                // this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Number Density : Start neighbor loop" << std::endl;

                // Index of neighbor
                unsigned int k = h_nlist.data[myHead + j];

                // Read particle j type
                bool j_issolid = checksolid(h_type_property_map.data, h_pos.data[k].w);

                // If both particles are solid, continue with next neighbor in loop
                if ( i_issolid && j_issolid )
                    continue;

                // Access neighbor position
                Scalar3 pj;
                pj.x = h_pos.data[k].x;
                pj.y = h_pos.data[k].y;
                pj.z = h_pos.data[k].z;

                // Compute distance vector
                Scalar3 dx = pj - pi;

                // Apply periodic boundary conditions
                dx = box.minImage(dx);

                // Calculate squared distance
                Scalar rsq = dot(dx, dx);

                // If particle distance is too large, continue with next neighbor in loop
                if ( this->m_const_slength && rsq > m_rcutsq )
                    continue;

                // Calculate distance
                Scalar r = sqrt(rsq);

                // If i_issolid and j_issolid - This part is not computed
                // If i_issolid and j_isfluid - Add contribution to normalization constant
                // If j_isfluid and j_issolid - Add contribution to particle number density
                // If j_isfluid and j_isfluid - Add contribution to particle number density

                h_dpe.data[i].x += this->m_const_slength ? this->m_skernel->wij(m_ch,r) : this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);

            } // End of neighbor loop

        // Compute mass density from number density if particle i is a fluid particle
        // rho_i = m_i * \sum_j wij
        if ( this->m_density_method == DENSITYSUMMATION && i_isfluid )
            {
            h_dpe.data[i].x= h_dpe.data[i].x * h_velocity.data[i].w;

            }

        } // End of particle loop

    // if (this->m_prof)
    //     this->m_prof->pop();
    }


/*! Perform pressure computation
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlow<KT_, SET_>::compute_pressure(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Pressure" << std::endl;

    // if (this->m_prof)
    //     this->m_prof->push("SinglePhaseFlowPressure");

    // Define ArrayHandles
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::readwrite);

    // For each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);
        // Evaluate pressure
        h_dpe.data[i].y = this->m_eos->Pressure(h_dpe.data[i].x);
        }

    // if (this->m_prof)
    //     this->m_prof->pop();
    }

/*! Compute fictitious solid particle properties
    This method updates fictitious solid particle pressures and velocities to account for
    no-slip boundary conditions. Method follows Adami et. al. (2012).
 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlow<KT_, SET_>::compute_noslip(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::NoSlip" << std::endl;

    // if (this->m_prof)
    //     this->m_prof->push("SinglePhaseFlowNoSlip");

    // Grab handles for particle and neighbor data
    ArrayHandle<Scalar3> h_vf(this->m_pdata->getAuxiliaries1(), access_location::host,access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_accel(this->m_pdata->getAccelerations(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    // ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    unsigned int size;
    long unsigned int myHead;


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
            // Solid particles which do not have fluid neighbors are marked
            // using energy=1 so that they can be deleted during simulation
            h_dpe.data[i].z = Scalar(1);
            // Set fictitious solid velocity to zero
            h_vf.data[i].x = 0;
            h_vf.data[i].y = 0;
            h_vf.data[i].z = 0;
            // If no fluid neighbors are present,
            // Set pressure to background pressure
            h_dpe.data[i].y = this->m_eos->getBackgroundPressure();
            // Density to rest density
            h_dpe.data[i].x = this->m_rho0;

            continue;
            }

        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];
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
            Scalar3 dx = pi - pj;

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
            Scalar Pj = h_dpe.data[k].y;

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
            if ( this->m_body_acceleration )
                {
                ph_c0.x += h_dpe.data[k].x * dx.x * wij;
                ph_c0.y += h_dpe.data[k].x * dx.y * wij;
                ph_c0.z += h_dpe.data[k].x * dx.z * wij;
                }
            } // End neighbor loop

        // Store fictitious solid particle velocity
        if (fluidneighbors > 0 && h_dpe.data[i].x > 0 )
            {
            //  Compute zeroth order normalization constant
            Scalar norm_constant = 1./h_dpe.data[i].x;
            // Set fictitious velocity
            h_vf.data[i].x = 2 * h_velocity.data[i].x - norm_constant * uf_c0.x;
            h_vf.data[i].y = 2 * h_velocity.data[i].y - norm_constant * uf_c0.y;
            h_vf.data[i].z = 2 * h_velocity.data[i].z - norm_constant * uf_c0.z;
            h_dpe.data[i].y = norm_constant * pf_c0 + dot( this->getAcceleration(timestep) - accel_i , norm_constant * ph_c0);
            // Compute solid densities by inverting equation of state
            h_dpe.data[i].x = this->m_eos->Density(h_dpe.data[i].y);
            }
        else
            {
            // Set fictitious solid velocity to zero
            h_vf.data[i].x = 0;
            h_vf.data[i].y = 0;
            h_vf.data[i].z = 0;

            // If no fluid neighbors are present,
            // Set pressure to background pressure
            h_dpe.data[i].y = this->m_eos->getBackgroundPressure();
            // Density to rest density
            h_dpe.data[i].x = this->m_rho0;
            }

        } // End solid particle loop

    // if (this->m_prof)
    //     this->m_prof->pop();
    }


/*! Perform Shepard density renormalization
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlow<KT_, SET_>::renormalize_density(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Density renormalization" << std::endl;

    // if (this->m_prof)
    //     this->m_prof->push("SinglePhaseFlowDensityRenormalization");

    // Grab handles for particle data
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    // ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);

    auto tmp_pde = this->m_pdata->getDPEs();
    ArrayHandle<Scalar3> h_dpe_old(tmp_pde, access_location::host, access_mode::read);

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
        Scalar rhoi = h_dpe.data[i].x;

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
                Scalar3 dx = pj - pi;

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
                Scalar Vj =  h_velocity.data[k].w / h_dpe_old.data[k].x ;
                normalization += this->m_const_slength ? Vj*this->m_skernel->wij(m_ch,r) : Vj*this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);

            } // End of neighbor loop

        normalization = Scalar(1.0)/normalization;

        // Initialize density with normalized kernel self density
        h_dpe.data[i].x = this->m_const_slength ? w0*(mi*normalization): this->m_skernel->w0(h_h.data[i])*(mi*normalization);

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
            Scalar3 dx = pj - pi;

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
            h_dpe.data[i].x += this->m_const_slength ? factor*this->m_skernel->wij(m_ch,r) : factor*this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);
            }

        } // End of particle loop

    // if (this->m_prof)
    //     this->m_prof->pop();
    }

/*! Perform force computation
 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlow<KT_, SET_>::forcecomputation(uint64_t timestep)
    {

    if ( m_density_method == DENSITYSUMMATION )
        this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Forces using SUMMATION approach " << m_density_method << endl;
    else if ( m_density_method == DENSITYCONTINUITY )
        this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Forces using CONTINUITY approach " << m_density_method << endl;

    // Start the profile for this compute
    // if (this->m_prof)
    //     this->m_prof->push("SinglePhaseFlowForces");

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
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
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

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Local variable to store things
    Scalar temp0;

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
        Scalar Pi = h_dpe.data[i].y;

        // Read particle i density and volume
        Scalar rhoi = h_dpe.data[i].x;
        Scalar Vi   = mi / rhoi;

        // // Total velocity of particle
        Scalar vi_total = sqrt((vi.x * vi.x) + (vi.y * vi.y) + (vi.z * vi.z));

        if (i == 0) {
            max_vel = vi_total;
        } else if (vi_total > max_vel) {
            max_vel = vi_total;
        }

        // Loop over all of the neighbors of this particle
        const long unsigned int myHead = h_head_list.data[i];
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
            bool issolid = checksolid(h_type_property_map.data, h_pos.data[k].w);

            // Compute distance vector (FLOPS: 3)
            Scalar3 dx = pi - pj;

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
            Scalar rhoj = h_dpe.data[k].x;
            Scalar Vj   = mj / rhoj;

            // Read particle k pressure
            Scalar Pj = h_dpe.data[k].y;

            // Compute velocity difference
            Scalar3 dv = vi - vj;

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
                temp0 = -(Vi*Vi+Vj*Vj)*((rhoj*Pi+rhoi*Pj)/(rhoi+rhoj));
            else
                temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);


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
            if ( issolid && m_compute_solid_forces )
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
            if ( issolid && m_compute_solid_forces )
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
                    dv = vi - vj;

                    //Vj = mj / m_rho0;
                    }

                // Compute density rate of change
                h_ratedpe.data[i].x += rhoi*Vj*dot(dv,dwdr_r*dx);
                //h_ratedpe.data[i].x += mj*dot(dv,dwdr_r*dx);

                // Add density diffusion if requested
                // Molteni and Colagrossi, Computer Physics Communications 180 (2009) 861–872
                if ( !issolid && m_density_diffusion )
                    h_ratedpe.data[i].x -= (Scalar(2)*m_ddiff*meanh*m_c*mj*(rhoi/rhoj-Scalar(1))*dot(dx,dwdr_r*dx))/(rsq+epssqr);
                }

            } // Closing Neighbor Loop

        } // Closing Fluid Particle Loop

    m_timestep_list[5] = max_vel;
    // Add volumetric force (gravity)
    this->applyBodyForce(timestep, m_fluidgroup);
    if ( m_compute_solid_forces )
        this->applyBodyForce(timestep, m_solidgroup);

    // if (this->m_prof)
    //     this->m_prof->pop();
    }

/*! Compute forces definition
*/

template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlow<KT_, SET_>::computeForces(uint64_t timestep)
    {

    // start by updating the neighborlist
    this->m_nlist->compute(timestep);

    // This is executed once to initialize protected/private variables
    if (!m_params_set)
        {
        this->m_exec_conf->msg->error() << "sph.models.SinglePhaseFlow requires parameters to be set before run()"
            << std::endl;
        throw std::runtime_error("Error computing SinglePhaseFlow forces");
        }

    // Apply density renormalization if requested
    if ( m_shepard_renormalization && timestep % m_shepardfreq == 0 )
        {
        renormalize_density(timestep);
#ifdef ENABLE_MPI
         // Update ghost particle densities and pressures.
        update_ghost_dpe(timestep);
#endif
        }

    // Compute solid renormalization constant and, provided that SUMMATION approach is used,
    // particle mass densities.
    compute_ndensity(timestep);
    // compute_particlenumberdensity(timestep);

    // Compute fluid pressure based on m_eos;
    compute_pressure(timestep);

#ifdef ENABLE_MPI
    // Update ghost particle densities and pressures.
    update_ghost_dpe(timestep);
#endif

    // Compute particle pressures
    compute_noslip(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    update_ghost_aux1(timestep);
#endif

    // Execute the force computation
    forcecomputation(timestep);

    // if (this->m_prof)
    //     this->m_prof->pop();
    }

namespace detail 
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlow(pybind11::module& m, std::string name)
{
    pybind11::class_<SinglePhaseFlow<KT_, SET_>, SPHBaseClass<KT_, SET_> , std::shared_ptr<SinglePhaseFlow<KT_, SET_>>>(m, name.c_str()) 
        .def(pybind11::init< std::shared_ptr<SystemDefinition>,
                             std::shared_ptr<SmoothingKernel<KT_> >,
                             std::shared_ptr<StateEquation<SET_> >,
                             std::shared_ptr<nsearch::NeighborList>,
                             std::shared_ptr<ParticleGroup>,
                             std::shared_ptr<ParticleGroup>,
                             DensityMethod,
                             ViscosityMethod >())
        .def("setParams", &SinglePhaseFlow<KT_, SET_>::setParams)
        .def("getDensityMethod", &SinglePhaseFlow<KT_, SET_>::getDensityMethod)
        .def("setDensityMethod", &SinglePhaseFlow<KT_, SET_>::setDensityMethod)
        .def("getViscosityMethod", &SinglePhaseFlow<KT_, SET_>::getViscosityMethod)
        .def("setViscosityMethod", &SinglePhaseFlow<KT_, SET_>::setViscosityMethod)
        .def("setConstSmoothingLength", &SinglePhaseFlow<KT_, SET_>::setConstSmoothingLength)
        .def("computeSolidForces", &SinglePhaseFlow<KT_, SET_>::computeSolidForces)
        .def("activateArtificialViscosity", &SinglePhaseFlow<KT_, SET_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity", &SinglePhaseFlow<KT_, SET_>::deactivateArtificialViscosity)
        .def("activateDensityDiffusion", &SinglePhaseFlow<KT_, SET_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion", &SinglePhaseFlow<KT_, SET_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &SinglePhaseFlow<KT_, SET_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization", &SinglePhaseFlow<KT_, SET_>::deactivateShepardRenormalization)
        .def("setAcceleration", &SPHBaseClass<KT_, SET_>::setAcceleration)
        .def("setRCut", &SinglePhaseFlow<KT_, SET_>::setRCutPython)
        ;

    }

// void export_SinglePhaseFlow()
// {
//     #define EXPORT(r,seq) export_SinglePhaseFlow_templ<BOOST_PP_SEQ_ENUM(seq)>();
//     BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPORT, (KERNELTYPES)(SEQTYPES))
//     #undef EXPORT
// }


// template<SmoothingKernelType KT_, StateEquationType SET_>
// void export_SinglePhaseFlow_templ()
// {
//     std::string name="SinglePF_" + get_SE_name<SET_>() + "_" + get_kernel_name<KT_>();
//     pybind11::class_<SinglePhaseFlow<KT_, SET_> , std::shared_ptr<SinglePhaseFlow<KT_, SET_>>>(m, name.c_str()) 
//         .def(pybind11::init< std::shared_ptr<SystemDefinition>,
//                              std::shared_ptr<SmoothingKernel<KT_> >,
//                              std::shared_ptr<StateEquation<SET_> >,
//                              std::shared_ptr<nsearch::NeighborList>,
//                              std::shared_ptr<ParticleGroup>,
//                              std::shared_ptr<ParticleGroup>,
//                              DensityMethod,
//                              ViscosityMethod >())
//         .def("setParams", &SinglePhaseFlow<KT_, SET_>::setParams)
//         .def("getDensityMethod", &SinglePhaseFlow<KT_, SET_>::getDensityMethod)
//         .def("setDensityMethod", &SinglePhaseFlow<KT_, SET_>::setDensityMethod)
//         .def("getViscosityMethod", &SinglePhaseFlow<KT_, SET_>::getViscosityMethod)
//         .def("setViscosityMethod", &SinglePhaseFlow<KT_, SET_>::setViscosityMethod)
//         .def("setConstSmoothingLength", &SinglePhaseFlow<KT_, SET_>::setConstSmoothingLength)
//         .def("computeSolidForces", &SinglePhaseFlow<KT_, SET_>::computeSolidForces)
//         .def("activateArtificialViscosity", &SinglePhaseFlow<KT_, SET_>::activateArtificialViscosity)
//         .def("deactivateArtificialViscosity", &SinglePhaseFlow<KT_, SET_>::deactivateArtificialViscosity)
//         .def("activateDensityDiffusion", &SinglePhaseFlow<KT_, SET_>::activateDensityDiffusion)
//         .def("deactivateDensityDiffusion", &SinglePhaseFlow<KT_, SET_>::deactivateDensityDiffusion)
//         .def("activateShepardRenormalization", &SinglePhaseFlow<KT_, SET_>::activateShepardRenormalization)
//         .def("deactivateShepardRenormalization", &SinglePhaseFlow<KT_, SET_>::deactivateShepardRenormalization)
//         .def("setAcceleration", &SPHBaseClass<KT_, SET_>::setAcceleration)
//         ;
//     }

// void export_SinglePhaseFlow()
// {
//     #define EXPORT(r,seq) export_SinglePhaseFlow_templ<BOOST_PP_SEQ_ENUM(seq)>();
//     BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPORT, (KERNELTYPES)(SEQTYPES))
//     #undef EXPORT
// }


} // end namespace detail

//! Explicit template instantiations
template class PYBIND11_EXPORT SinglePhaseFlow<wendlandc2, linear>;
template class PYBIND11_EXPORT SinglePhaseFlow<wendlandc2, tait>;
template class PYBIND11_EXPORT SinglePhaseFlow<wendlandc4, linear>;
template class PYBIND11_EXPORT SinglePhaseFlow<wendlandc4, tait>;
template class PYBIND11_EXPORT SinglePhaseFlow<wendlandc6, linear>;
template class PYBIND11_EXPORT SinglePhaseFlow<wendlandc6, tait>;
template class PYBIND11_EXPORT SinglePhaseFlow<quintic, linear>;
template class PYBIND11_EXPORT SinglePhaseFlow<quintic, tait>;
template class PYBIND11_EXPORT SinglePhaseFlow<cubicspline, linear>;
template class PYBIND11_EXPORT SinglePhaseFlow<cubicspline, tait>;


namespace detail
{

    template void export_SinglePhaseFlow<wendlandc2, linear>(pybind11::module& m, std::string name = "SinglePF_WC2_L");
    template void export_SinglePhaseFlow<wendlandc2, tait>(pybind11::module& m, std::string name = "SinglePF_WC2_T");
    template void export_SinglePhaseFlow<wendlandc4, linear>(pybind11::module& m, std::string name = "SinglePF_WC4_L");
    template void export_SinglePhaseFlow<wendlandc4, tait>(pybind11::module& m, std::string name = "SinglePF_WC4_T");
    template void export_SinglePhaseFlow<wendlandc6, linear>(pybind11::module& m, std::string name = "SinglePF_WC6_L");
    template void export_SinglePhaseFlow<wendlandc6, tait>(pybind11::module& m, std::string name = "SinglePF_WC6_T");
    template void export_SinglePhaseFlow<quintic, linear>(pybind11::module& m, std::string name = "SinglePF_Q_L");
    template void export_SinglePhaseFlow<quintic, tait>(pybind11::module& m, std::string name = "SinglePF_Q_T");
    template void export_SinglePhaseFlow<cubicspline, linear>(pybind11::module& m, std::string name = "SinglePF_CS_L");
    template void export_SinglePhaseFlow<cubicspline, tait>(pybind11::module& m, std::string name = "SinglePF_CS_T");

}  // end namespace detail
} // end namespace sph
} // end namespace hoomd
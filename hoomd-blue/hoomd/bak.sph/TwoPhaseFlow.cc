// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

#include "TwoPhaseFlow.h"

// #include <boost/python.hpp>
// using namespace boost::python;

// #include <boost/shared_ptr.hpp>
// #include <boost/bind.hpp>
// using namespace boost;

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

using namespace std;

namespace hoomd
{
namespace sph
{
/*! Constructor
*/
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
TwoPhaseFlow<KT_, SET1_, SET2_>::TwoPhaseFlow(std::shared_ptr<SystemDefinition> sysdef,
                           std::shared_ptr<SmoothingKernel<KT_> > skernel,
                           std::shared_ptr<StateEquation<SET1_> > equationofstate1,
                           std::shared_ptr<StateEquation<SET2_> > equationofstate2,
                           std::shared_ptr<nsearch::NeighborList> nlist,
                           std::shared_ptr<ParticleGroup> fluidgroup1,
                           std::shared_ptr<ParticleGroup> fluidgroup2,
                           std::shared_ptr<ParticleGroup> solidgroup,
                           DensityMethod mdensitymethod,
                           ViscosityMethod mviscositymethod)
    : SPHBaseClass<KT_, SET1_>(sysdef,skernel,equationofstate1,nlist), m_fluidgroup1(fluidgroup1), m_fluidgroup2(fluidgroup2),
      m_solidgroup(solidgroup), m_eos1(equationofstate1), m_eos2(equationofstate2)
      {
        this->m_exec_conf->msg->notice(5) << "Constructing TwoPhaseFlow" << std::endl;

        // Set private attributes to default values
        this->m_const_slength = false;
        this->m_params_set = false;
        this->m_compute_solid_forces = false;
        this->m_artificial_viscosity = false;
        this->m_density_diffusion = false;
        this->m_shepard_renormalization = false;
        this->m_ch = Scalar(0.0);
        this->m_rcut = Scalar(0.0);
        this->m_rcutsq = Scalar(0.0);

        this->m_avalpha = Scalar(0.0);
        this->m_avbeta = Scalar(0.0);
        this->m_ddiff = Scalar(0.0);
        this->m_shepardfreq = 1;

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
        this->constructTypeVectors(solidgroup,&m_solidtypes);
        this->constructTypeVectors(this->m_fluidgroup,&m_fluidtypes);

        // all particle groups are based on the same particle data
        unsigned int num_types = this->m_sysdef->getParticleData()->getNTypes();
        m_type_property_map = GPUArray<unsigned int>(num_types, this->m_exec_conf);
        {
            ArrayHandle<unsigned int> h_type_property_map(m_type_property_map, access_location::host, access_mode::overwrite);
            fill_n(h_type_property_map.data, num_types, SolidFluidTypeBit::NONE);
            // no need to parallelize this as there should only be a few particle types
            for (unsigned int i = 0; i < m_fluidtypes1.size(); i++) {
                h_type_property_map.data[m_fluidtypes[i]] |= SolidFluidTypeBit::FLUID | SolidFluidTypeBit::FLUID1;
            }
            for (unsigned int i = 0; i < m_fluidtypes2.size(); i++) {
                h_type_property_map.data[m_fluidtypes[i]] |= SolidFluidTypeBit::FLUID | SolidFluidTypeBit::FLUID2;
            }
            for (unsigned int i = 0; i < m_solidtypes.size(); i++) {
                h_type_property_map.data[m_solidtypes[i]] |= SolidFluidTypeBit::SOLID;
            }
        }

        // Set simulations methods
        this->m_density_method = mdensitymethod;
        this->m_viscosity_method = mviscositymethod;

        // Get necessary variables from kernel and EOS classes
        this->m_rho01 = equationofstate1->getRestDensity();
        this->m_rho02 = equationofstate2->getRestDensity();
        this->m_c1    = equationofstate1->getSpeedOfSound();
        this->m_c2    = equationofstate1->getSpeedOfSound();
        this->m_cmax  = this->m_c1 > this->m_c2 ? this->m_c1 : this->m_c2;
        this->m_kappa = skernel->getKernelKappa();
      }

/*! Destructor
*/

template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
TwoPhaseFlow<KT_, SET1_, SET2_>::~TwoPhaseFlow()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying TwoPhaseFlow" << std::endl;
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

/*! Communicate (update) ghost particle fields
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::update_ghost_dpe(unsigned int timestep)
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
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::update_ghost_aux123(unsigned int timestep)
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
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::update_ghost_aux4(unsigned int timestep)
    {
#ifdef ENABLE_MPI
    if (this->m_comm)
        {
        CommFlags flags(0);
        flags[comm_flag::tag] = 0;
        flags[comm_flag::position] = 0;
        flags[comm_flag::velocity] = 0;
        flags[comm_flag::dpe] = 0;
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

/*! Perform number density computation
 * This method computes and stores
     - the number density based mass density ( rho_i = this->m_i * \sum w_ij ) for fluid particles
       if the SUMMATION approach is being used.
     - the zeroth order normalization constant for solid particles
   in the x-component of the h_dpe Array.
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_ndensity(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Number Density" << std::endl;

    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowNDensity");

    // Grab handles for particle data
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Precompute self-density for homogeneous smoothing lengths
    Scalar w0 = this->m_skernel->w0(this->m_ch);

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
            continue;

        // Initialize number density with self density of kernel
        if ( i_issolid )
            h_dpe.data[i].x = 0;
        else
            h_dpe.data[i].x = this->m_const_slength ? w0 : this->m_skernel->w0(h_h.data[i]);

        // Loop over all of the neighbors of this particle
        const unsigned int myHead = h_head_list.data[i];
        const unsigned int size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {
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
                if ( this->m_const_slength && rsq > this->m_rcutsq )
                    continue;

                // Calculate distance
                Scalar r = sqrt(rsq);

                // If i_issolid and j_issolid - This part is not computed
                // If i_issolid and j_isfluid - Add contribution to normalization constant
                // If j_isfluid and j_issolid - Add contribution to particle number density
                // If j_isfluid and j_isfluid - Add contribution to particle number density
                h_dpe.data[i].x += this->m_const_slength ? this->m_skernel->wij(this->m_ch,r) : this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);

            } // End of neighbor loop

        // Compute mass density from number density if particle i is a fluid particle
        // rho_i = this->m_i * \sum_j wij
        if ( this->m_density_method == DENSITYSUMMATION && i_isfluid )
            h_dpe.data[i].x= h_dpe.data[i].x * h_velocity.data[i].w;

        } // End of particle loop

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Perform pressure computation
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_pressure(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Pressure" << std::endl;

    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowPressure");

    // Define ArrayHandles
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::readwrite);

    // For each fluid particle in phase 1
    unsigned int group_size = this->m_fluidgroup1->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_fluidgroup1->getMemberIndex(group_idx);
        // Evaluate pressure
        h_dpe.data[i].y = this->m_eos1->Pressure(h_dpe.data[i].x);
        }

    // For each fluid particle in phase 2
    group_size = this->m_fluidgroup2->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_fluidgroup2->getMemberIndex(group_idx);
        // Evaluate pressure
        h_dpe.data[i].y = this->m_eos1->Pressure(h_dpe.data[i].x);
        }

    if (this->m_prof)
        this->m_prof->pop();
    }


/*! Compute fictitious solid particle properties
    This method updates fictitious solid particle pressures and velocities to account for
    no-slip boundary conditions. Method follows Adami et. al. (2012).
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_noslip(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::NoSlip" << std::endl;

    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowNoSlip");

    // Grab handles for particle and neighbor data
    ArrayHandle<Scalar3> h_vf(this->m_pdata->getAuxiliaries1(), access_location::host,access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_accel(this->m_pdata->getAccelerations(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Count fluid neighbors before setting solid particle properties
    unsigned int fluidneighbors;

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

        // Read normalization constant of solid particle i
        Scalar norm_constant = 1./h_dpe.data[i].x;

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
        fluidneighbors = 0;
        const unsigned int myHead = h_head_list.data[i];
        const unsigned int size = (unsigned int)h_n_neigh.data[i];
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
            Scalar wij = this->m_const_slength ? this->m_skernel->wij(this->m_ch,r) : this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);

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
        if (fluidneighbors > 0)
            {
            h_vf.data[i].x = 2 * h_velocity.data[i].x - norm_constant * uf_c0.x;
            h_vf.data[i].y = 2 * h_velocity.data[i].y - norm_constant * uf_c0.y;
            h_vf.data[i].z = 2 * h_velocity.data[i].z - norm_constant * uf_c0.z;
            }
        else
            {
            h_vf.data[i].x = 0;
            h_vf.data[i].y = 0;
            h_vf.data[i].z = 0;
            }

        // Store fictitious solid particle pressure
        if (fluidneighbors > 0)
            h_dpe.data[i].y = norm_constant * pf_c0 + dot( this->getAcceleration(timestep) - accel_i , norm_constant * ph_c0);
        else
            h_dpe.data[i].y = this->m_eos1->getBackgroundPressure();

        // Compute solid densities by inverting equation of state
        h_dpe.data[i].x = this->m_eos1->Density(h_dpe.data[i].y);

        } // End solid particle loop

    if (this->m_prof)
        this->m_prof->pop();
    }


/*! Perform Shepard density renormalization
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::renormalize_density(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Density renormalization" << std::endl;

    if (this->m_prof)
        this->m_prof->push("SinglePhaseFlowDensityRenormalization");

    // Grab handles for particle data
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Precompute self-density for homogeneous smoothing lengths
    Scalar w0 = this->m_skernel->w0(this->m_ch);

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

        // Check if fluid is of phase 1
        bool isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[i].w);

        // Loop over all of the neighbors of this particle
        unsigned int myHead = h_head_list.data[i];
        unsigned int size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {
                // Index of neighbor
                unsigned int k = h_nlist.data[myHead + j];

                // Only interpolate over fluid domain of same phase
                if ( checkfluid1(h_type_property_map.data, h_pos.data[k].w) != isfluid1 )
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
                if ( this->m_const_slength && rsq > this->m_rcutsq )
                    continue;

                // Calculate distance
                Scalar r = sqrt(rsq);

                // Add contribution to renormalization
                Scalar Vj =  h_velocity.data[k].w / h_dpe.data[k].x ;
                normalization += this->m_const_slength ? Vj*this->m_skernel->wij(this->m_ch,r) : Vj*this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);

            } // End of neighbor loop

        normalization = Scalar(1.0)/normalization;

        // Initialize density with normalized kernel self density
        h_dpe.data[i].x = this->m_const_slength ? w0*(mi*normalization): this->m_skernel->w0(h_h.data[i])*(mi*normalization);

        // Loop over all of the neighbors of this particle
        myHead = h_head_list.data[i];
        size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int j = 0; j < size; j++)
            {
                // Index of neighbor
                unsigned int k = h_nlist.data[myHead + j];

                // Only interpolate over fluid domain of same phase
                if ( checkfluid1(h_type_property_map.data, h_pos.data[k].w) != isfluid1 )
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
                if ( this->m_const_slength && rsq > this->m_rcutsq )
                    continue;

                // Calculate distance
                Scalar r = sqrt(rsq);

                // Add contribution to normalized density interpolation
                Scalar factor =  h_velocity.data[k].w * normalization ;
                h_dpe.data[i].x += this->m_const_slength ? factor*this->m_skernel->wij(this->m_ch,r) : factor*this->m_skernel->wij(Scalar(0.5)*(h_h.data[i]+h_h.data[k]),r);
            }

        } // End of particle loop

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Compute interfacial color gradients
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_colorgradients(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Normals" << std::endl;

    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowNormals");

    // Grab handles for particle and neighbor data
    ArrayHandle<Scalar3> h_sn(this->m_pdata->getAuxiliaries2(), access_location::host,access_mode::readwrite);
    ArrayHandle<Scalar3> h_fn(this->m_pdata->getAuxiliaries3(), access_location::host,access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
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
        Scalar rhoi = h_dpe.data[i].x;
        Scalar Vi   = mi / rhoi;

        // Detect particle i type
        bool i_issolid = checksolid(h_type_property_map.data, h_pos.data[i].w);
        bool i_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[i].w);

        // Loop over all of the neighbors of this particle
        const unsigned int myHead = h_head_list.data[i];
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

            // Skip color gradient computation if both particles belong to same phase
            if ( (i_issolid==true && j_issolid==true) || i_isfluid1 == j_isfluid1 )
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
            Scalar rhoj = h_dpe.data[k].x;
            Scalar Vj   = mj / rhoj;

            // Calculate absolute and normalized distance
            Scalar r = sqrt(rsq);

            // Mean smoothing length and denominator modifier
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar eps    = Scalar(0.1)*meanh;

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);

            // Evaluate interpolation contribution
            //Scalar temp0 = (Vi*Vi+Vj*Vj)*(rhoi/(rhoi+rhoj))*(mi/rhoj);
            Scalar temp0 = (Vj*Vj/Vi);

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

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Compute surface force vectors
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::compute_surfaceforce(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::SurfaceForce" << std::endl;

    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowSurfaceForce");

    // Grab handles for particle and neighbor data
    ArrayHandle<Scalar3> h_sf(this->m_pdata->getAuxiliaries4(), access_location::host,access_mode::readwrite);
    ArrayHandle<Scalar3> h_sn(this->m_pdata->getAuxiliaries2(), access_location::host,access_mode::read);
    ArrayHandle<Scalar3> h_fn(this->m_pdata->getAuxiliaries3(), access_location::host,access_mode::read);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_dpe(this->m_pdata->getDPEs(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
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

        // Check if there is any fluid particle near the current particle, if not continue
        // This makes sure that only particle near a fluid interface experience an interfacial force.
        // In other words, fluid particles only near solid interfaces are omitted.
        bool nearfluidinterface = false;
        // Loop over all of the neighbors of this particle
        unsigned int myHead = h_head_list.data[i];
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
        Scalar rhoi = h_dpe.data[i].x;
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
        // Fluid phase 1 - Fluid phase 2 interface
        if ( this->m_sigma12 > 0 && normfni > 0 )
            {
            Scalar temp0 = (this->m_sigma12/normfni);
            istress[0] +=  temp0*(normfni*normfni-fni.x*fni.x); // xx
            istress[1] +=  temp0*(normfni*normfni-fni.y*fni.y); // yy
            istress[2] +=  temp0*(normfni*normfni-fni.z*fni.z); // zz
            istress[3] += -temp0*(fni.x*fni.y);                 // xy yx
            istress[4] += -temp0*(fni.x*fni.z);                 // xz zx
            istress[5] += -temp0*(fni.y*fni.z);                 // yz zy
            }
        // Fluid phase 1 - Solid interface
        if ( i_isfluid1 && this->m_sigma01 > 0 && normsni > 0 )
            {
            Scalar temp0 = (this->m_sigma01/normsni);
            istress[0] +=  temp0*(normsni*normsni-sni.x*sni.x); // xx
            istress[1] +=  temp0*(normsni*normsni-sni.y*sni.y); // yy
            istress[2] +=  temp0*(normsni*normsni-sni.z*sni.z); // zz
            istress[3] += -temp0*(sni.x*sni.y);                 // xy yx
            istress[4] += -temp0*(sni.x*sni.z);                 // xz zx
            istress[5] += -temp0*(sni.y*sni.z);                 // yz zy
            }
        // Fluid phase 2 - Solid interface
        if ( !i_isfluid1 && this->m_sigma02 > 0 && normsni > 0 )
            {
            Scalar temp0 = (this->m_sigma02/normsni);
            istress[0] +=  temp0*(normsni*normsni-sni.x*sni.x); // xx
            istress[1] +=  temp0*(normsni*normsni-sni.y*sni.y); // yy
            istress[2] +=  temp0*(normsni*normsni-sni.z*sni.z); // zz
            istress[3] += -temp0*(sni.x*sni.y);                 // xy yx
            istress[4] += -temp0*(sni.x*sni.z);                 // xz zx
            istress[5] += -temp0*(sni.y*sni.z);                 // yz zy
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
            Scalar rhoj = h_dpe.data[k].x;
            Scalar Vj   = mj / rhoj;

            // Mean smoothing length and denominator modifier
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar eps    = Scalar(0.1)*meanh;

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);

            // Evaluate particle i interfacial stress tensor
            Scalar jstress[6] = {0};
            // Fluid phase 1 - Fluid phase 2 interface
            if ( !(j_issolid) && this->m_sigma12 > 0 && normfnj > 0 )
                {
                Scalar temp0 = (this->m_sigma12/normfnj);
                jstress[0] +=  temp0*(normfnj*normfnj-fnj.x*fnj.x); // xx
                jstress[1] +=  temp0*(normfnj*normfnj-fnj.y*fnj.y); // yy
                jstress[2] +=  temp0*(normfnj*normfnj-fnj.z*fnj.z); // zz
                jstress[3] += -temp0*(fnj.x*fnj.y);                 // xy yx
                jstress[4] += -temp0*(fnj.x*fnj.z);                 // xz zx
                jstress[5] += -temp0*(fnj.y*fnj.z);                 // yz zy
                }
            // Fluid phase 1 - Solid interface
            if ( j_isfluid1 && this->m_sigma01 > 0 && normsnj > 0 )
                {
                Scalar temp0 = (this->m_sigma01/normsnj);
                jstress[0] +=  temp0*(normsnj*normsnj-snj.x*snj.x); // xx
                jstress[1] +=  temp0*(normsnj*normsnj-snj.y*snj.y); // yy
                jstress[2] +=  temp0*(normsnj*normsnj-snj.z*snj.z); // zz
                jstress[3] += -temp0*(snj.x*snj.y);                 // xy yx
                jstress[4] += -temp0*(snj.x*snj.z);                 // xz zx
                jstress[5] += -temp0*(snj.y*snj.z);                 // yz zy
                }
            // Fluid phase 2 - Solid interface
            if ( !j_isfluid1 && this->m_sigma02 > 0 && normsnj > 0 )
                {
                Scalar temp0 = (this->m_sigma02/normsnj);
                jstress[0] +=  temp0*(normsnj*normsnj-snj.x*snj.x); // xx
                jstress[1] +=  temp0*(normsnj*normsnj-snj.y*snj.y); // yy
                jstress[2] +=  temp0*(normsnj*normsnj-snj.z*snj.z); // zz
                jstress[3] += -temp0*(snj.x*snj.y);                 // xy yx
                jstress[4] += -temp0*(snj.x*snj.z);                 // xz zx
                jstress[5] += -temp0*(snj.y*snj.z);                 // yz zy
                }

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

            } // End of neighbor loop

        // Set component normal to solid surface at solid interface to zero

        } // Closing Fluid Particle Loop

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Perform force computation
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::forcecomputation(unsigned int timestep)
    {

    if ( this->m_density_method == DENSITYSUMMATION )
        this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Forces using SUMMATION approach" << endl;
    else if ( this->m_density_method == DENSITYCONTINUITY )
        this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Forces using CONTINUITY approach" << endl;

    // Start the profile for this compute
    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowForces");

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
    ArrayHandle<Scalar3> h_sf(this->m_pdata->getAuxiliaries4(), access_location::host,access_mode::read);

    // access the neighbor list
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Check input data
    assert(h_pos.data != NULL);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Local variable to store things
    Scalar temp0;

    // for each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

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

        // Read particle i type, viscosity, speed of sound and rest density
        bool i_isfluid1 = checkfluid1(h_type_property_map.data, h_pos.data[i].w);
        Scalar mui   = i_isfluid1 ? this->m_mu1 : this->m_mu2;
        Scalar rho0i = i_isfluid1 ? this->m_rho01 : this->m_rho02;
        Scalar ci    = i_isfluid1 ? this->m_c1 : this->m_c2;

        // Loop over all of the neighbors of this particle
        const unsigned int myHead = h_head_list.data[i];
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

            // Read particle j viscosity, speed of sound and rest density
            Scalar muj   = j_isfluid1 ? this->m_mu1 : this->m_mu2;
            Scalar rho0j = j_isfluid1 ? this->m_rho01 : this->m_rho02;
            Scalar cj    = j_isfluid1 ? this->m_c1 : this->m_c2;
            // If particle j is solid, set parameters equal to those of particle i
            muj   = j_issolid ? mui : muj;
            rho0j = j_issolid ? rho0i : rho0j;
            cj    = j_issolid ? ci : cj;

            // Compute distance vector (FLOPS: 3)
            Scalar3 dx = pi - pj;

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
            Scalar rhoj = h_dpe.data[k].x;
            Scalar Vj   = mj / rhoj;

            // Read particle k pressure
            Scalar Pj = h_dpe.data[k].y;

            // Compute velocity difference
            Scalar3 dv = vi - vj;

            // Calculate absolute and normalized distance
            Scalar r = sqrt(rsq);

            // Mean smoothing length and denominator modifier
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar eps    = Scalar(0.1)*meanh;
            Scalar epssqr = eps*eps;

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);

            // Evaluate inter-particle pressure forces
            if ( this->m_density_method == DENSITYSUMMATION )
                temp0 = -(Vi*Vi+Vj*Vj)*((rhoj*Pi+rhoi*Pj)/(rhoi+rhoj));
            else
                temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);

            // Optionally add artificial viscosity
            // Monaghan (1983) J. Comput. Phys. 52 (2) 374–389
            if ( this->m_artificial_viscosity && !j_issolid )
                {
                Scalar dotdvdx = dot(dv,dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar muij    = meanh*dotdvdx/(rsq+epssqr);
                    Scalar meanrho = Scalar(0.5)*(rhoi+rhoj);
                    temp0 += mi*mj*(this->m_avalpha*m_cmax*muij+m_avbeta*muij*muij)/meanrho;
                    }
                }

            // Add contribution to fluid particle
            h_force.data[i].x += temp0*dwdr_r*dx.x;
            h_force.data[i].y += temp0*dwdr_r*dx.y;
            h_force.data[i].z += temp0*dwdr_r*dx.z;

            // Add contribution to solid particle
            if ( j_issolid && this->m_compute_solid_forces )
                {
                h_force.data[k].x -= (mj/mi)*temp0*dwdr_r*dx.x;
                h_force.data[k].y -= (mj/mi)*temp0*dwdr_r*dx.y;
                h_force.data[k].z -= (mj/mi)*temp0*dwdr_r*dx.z;
                }

            // Evaluate viscous interaction forces
            temp0 = ((Scalar(2)*(mui*muj))/(mui+muj)) * (Vi*Vi+Vj*Vj) * dwdr_r;
            h_force.data[i].x  += temp0*dv.x;
            h_force.data[i].y  += temp0*dv.y;
            h_force.data[i].z  += temp0*dv.z;

            // Add contribution to solid particle
            if ( j_issolid && this->m_compute_solid_forces )
                {
                h_force.data[k].x -= (mj/mi)*temp0*dv.x;
                h_force.data[k].y -= (mj/mi)*temp0*dv.y;
                h_force.data[k].z -= (mj/mi)*temp0*dv.z;
                }

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
                    dv = vi - vj;
                    }

                // Compute density rate of change
                h_ratedpe.data[i].x += rhoi*Vj*dot(dv,dwdr_r*dx);

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

    // Add volumetric force (gravity)
    this->applyBodyForce(timestep, this->m_fluidgroup);

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Compute forces definition
*/
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlow<KT_, SET1_, SET2_>::computeForces(unsigned int timestep)
    {
    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlow");

    // This is executed once to initialize protected/private variables
    if (!m_params_set)
        {
        this->m_exec_conf->msg->error() << "sph.models.TwoPhaseFlow requires parameters to be set before run()"
            << std::endl;
        throw std::runtime_error("Error computing TwoPhaseFlow forces");
        }

    // Make sure neighbor list is up-to-date
    this->m_nlist->compute(timestep);

    // Apply density renormalization if requested
    if ( this->m_shepard_renormalization && timestep % this->m_shepardfreq == 0 )
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

    // Compute fluid pressure based on this->m_eos;
    compute_pressure(timestep);

#ifdef ENABLE_MPI
    // Update ghost particle densities and pressures.
    update_ghost_dpe(timestep);
#endif

    // Compute particle pressures
    compute_noslip(timestep);

#ifdef ENABLE_MPI
    // Update ghost particle densities and pressures.
    update_ghost_dpe(timestep);
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
    forcecomputation(timestep);

    if (this->m_prof)
        this->m_prof->pop();
    }

template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void export_TwoPhaseFlow(pybind11::module& m)
    {
    std::string name="TwoPF_" + get_SE_name<SET1_>()+ get_SE_name<SET2_>() + "_" + get_kernel_name<KT_>();
    pybind11::class_<TwoPhaseFlow<KT_,SET1_,SET2_>, std::shared_ptr<TwoPhaseFlow<KT_,SET1_,SET2_>>(m, name.c_str())
        .def(pybind::init< std::shared_ptr<SystemDefinition>,
                           std::shared_ptr<SmoothingKernel<KT_> >,
                           std::shared_ptr<StateEquation<SET1_> >,
                           std::shared_ptr<StateEquation<SET2_> >,
                           std::shared_ptr<nsearch::NeighborList>,
                           std::shared_ptr<ParticleGroup>,
                           std::shared_ptr<ParticleGroup>,
                           std::shared_ptr<ParticleGroup>,
                           DensityMethod,
                           ViscosityMethod >())
        .def("setParams", &TwoPhaseFlow<KT_,SET1_,SET2_>::setParams)
        .def("getDensityMethod", &TwoPhaseFlow<KT_,SET1_,SET2_>::getDensityMethod)
        .def("setDensityMethod", &TwoPhaseFlow<KT_,SET1_,SET2_>::setDensityMethod)
        .def("getViscosityMethod", &TwoPhaseFlow<KT_,SET1_,SET2_>::getViscosityMethod)
        .def("setViscosityMethod", &TwoPhaseFlow<KT_,SET1_,SET2_>::setViscosityMethod)
        .def("setConstSmoothingLength", &TwoPhaseFlow<KT_,SET1_,SET2_>::setConstSmoothingLength)
        .def("computeSolidForces", &TwoPhaseFlow<KT_,SET1_,SET2_>::computeSolidForces)
        .def("activateArtificialViscosity", &TwoPhaseFlow<KT_,SET1_,SET2_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity", &TwoPhaseFlow<KT_,SET1_,SET2_>::deactivateArtificialViscosity)
        .def("activateDensityDiffusion", &TwoPhaseFlow<KT_,SET1_,SET2_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion", &TwoPhaseFlow<KT_,SET1_,SET2_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &TwoPhaseFlow<KT_,SET1_,SET2_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization", &TwoPhaseFlow<KT_,SET1_,SET2_>::deactivateShepardRenormalization)
        .def("setAcceleration", &SPHBaseClass<KT_,SET2_>::setAcceleration)
        ;
    }


// void export_TwoPhaseFlow(pybind11::module& m)
// {
//     #define EXPORT(r,seq) export_TwoPhaseFlow_templ<BOOST_PP_SEQ_ENUM(seq)>();
//     BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPORT, (KERNELTYPES)(SEQTYPES)(SEQTYPES))
//     #undef EXPORT
// }

} // end namespace sph
} // end namespace hoomd

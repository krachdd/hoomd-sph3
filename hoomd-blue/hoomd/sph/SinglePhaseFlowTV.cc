/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SinglePhaseFlowTV.h"

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
SinglePhaseFlowTV<KT_, SET_>::SinglePhaseFlowTV(std::shared_ptr<SystemDefinition> sysdef,
                                 std::shared_ptr<SmoothingKernel<KT_> > skernel,
                                 std::shared_ptr<StateEquation<SET_> > equationofstate,
                                 std::shared_ptr<nsearch::NeighborList> nlist,
                                 std::shared_ptr<ParticleGroup> fluidgroup,
                                 std::shared_ptr<ParticleGroup> solidgroup,
                                 DensityMethod mdensitymethod,
                                 ViscosityMethod mviscositymethod)
    : SinglePhaseFlow<KT_, SET_>(sysdef,skernel,equationofstate,nlist,fluidgroup,solidgroup,mdensitymethod,mviscositymethod)
      {
        this->m_exec_conf->msg->notice(5) << "Constructing SinglePhaseFlowTV" << std::endl;

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

        this->m_solid_removed = false;

        // Sanity checks
        assert(this->m_pdata);
        assert(this->m_nlist);
        assert(this->m_skernel);
        assert(this->m_eos);

        // Contruct type vectors
        this->constructTypeVectors(fluidgroup,&this->m_fluidtypes);
        this->constructTypeVectors(solidgroup,&this->m_solidtypes);

        // all particle groups are based on the same particle data
        unsigned int num_types = this->m_sysdef->getParticleData()->getNTypes();

        this->m_type_property_map = GPUArray<unsigned int>(num_types, this->m_exec_conf);
        {
            ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::overwrite);
            fill_n(h_type_property_map.data, num_types, SolidFluidTypeBit::NONE);
            // no need to parallelize this as there should only be a few particle types
            for (unsigned int i = 0; i < this->m_fluidtypes.size(); i++) {
                h_type_property_map.data[this->m_fluidtypes[i]] |= SolidFluidTypeBit::FLUID;
            }
            for (unsigned int i = 0; i < this->m_solidtypes.size(); i++) {
                h_type_property_map.data[this->m_solidtypes[i]] |= SolidFluidTypeBit::SOLID;
            }
        }

        // Set simulations methods
        this->m_density_method = mdensitymethod;
        this->m_viscosity_method = mviscositymethod;

        // Get necessary variables from kernel and EOS classes
        this->m_rho0  = equationofstate->getRestDensity();
        this->m_c     = equationofstate->getSpeedOfSound();
        this->m_kappa = skernel->getKernelKappa();

        this->m_r_cut_nlist = std::make_shared<GlobalArray<Scalar>>(this->m_typpair_idx.getNumElements(), this->m_exec_conf);
        this->m_nlist->addRCutMatrix(this->m_r_cut_nlist);

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
template<SmoothingKernelType KT_,StateEquationType SET_>
SinglePhaseFlowTV<KT_, SET_>::~SinglePhaseFlowTV()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying SinglePhaseFlowTV" << std::endl;
    }



template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlowTV<KT_, SET_>::update_ghost_aux123(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlowTV::Update Ghost aux123" << std::endl;

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
        flags[comm_flag::auxiliary2] = 2;
        flags[comm_flag::auxiliary3] = 3;
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

/*! Perform force computation
 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlowTV<KT_, SET_>::forcecomputation(uint64_t timestep)
    {

    if ( this->m_density_method == DENSITYSUMMATION )
        this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlowTV::Forces using SUMMATION approach " << this->m_density_method << endl;
    else if ( this->m_density_method == DENSITYCONTINUITY )
        this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlowTV::Forces using CONTINUITY approach " << this->m_density_method << endl;

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
    ArrayHandle<Scalar3> h_bpc(this->m_pdata->getAuxiliaries2(), access_location::host,access_mode::readwrite); // background pressure contribution to tv
    ArrayHandle<Scalar3> h_tv(this->m_pdata->getAuxiliaries3(), access_location::host,access_mode::read); // transport velocity of the particle tv
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
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

        // set background pressure contribution for tv to zero 
        h_bpc.data[i].x = 0.0;
        h_bpc.data[i].y = 0.0;
        h_bpc.data[i].z = 0.0;

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

        Scalar3 tvi;
        tvi.x = h_tv.data[i].x;
        tvi.y = h_tv.data[i].y;
        tvi.z = h_tv.data[i].z;

        // Read particle i pressure
        Scalar Pi = h_pressure.data[i];

        // Read particle i density and volume
        Scalar rhoi = h_density.data[i];
        Scalar Vi   = mi / rhoi;

        // // Total velocity of particle
        Scalar vi_total = sqrt((vi.x * vi.x) + (vi.y * vi.y) + (vi.z * vi.z));

        // compute \mathbf{A}_i, artificial stress tensor of particle i 
        Scalar A11i = rhoi * vi.x * ( tvi.x - vi.x );
        Scalar A12i = rhoi * vi.x * ( tvi.y - vi.y );
        Scalar A13i = rhoi * vi.x * ( tvi.z - vi.z );
        Scalar A21i = rhoi * vi.y * ( tvi.x - vi.x );
        Scalar A22i = rhoi * vi.y * ( tvi.y - vi.y );
        Scalar A23i = rhoi * vi.y * ( tvi.z - vi.z );
        Scalar A31i = rhoi * vi.z * ( tvi.x - vi.x );
        Scalar A32i = rhoi * vi.z * ( tvi.y - vi.y );
        Scalar A33i = rhoi * vi.z * ( tvi.z - vi.z );

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

            Scalar3 tvj;
            tvj.x = h_tv.data[k].x;
            tvj.y = h_tv.data[k].y;
            tvj.z = h_tv.data[k].z;

            // compute \mathbf{A}_j, artificial stress tensor of particle j 
            Scalar A11j = rhoj * vj.x * ( tvj.x - vj.x );
            Scalar A12j = rhoj * vj.x * ( tvj.y - vj.y );
            Scalar A13j = rhoj * vj.x * ( tvj.z - vj.z );
            Scalar A21j = rhoj * vj.y * ( tvj.x - vj.x );
            Scalar A22j = rhoj * vj.y * ( tvj.y - vj.y );
            Scalar A23j = rhoj * vj.y * ( tvj.z - vj.z );
            Scalar A31j = rhoj * vj.z * ( tvj.x - vj.x );
            Scalar A32j = rhoj * vj.z * ( tvj.y - vj.y );
            Scalar A33j = rhoj * vj.z * ( tvj.z - vj.z );


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
            temp0 = (rhoj*Pi+rhoi*Pj)/(rhoi+rhoj); 

            Scalar avc = 0.0;
            // Optionally add artificial viscosity
            // Monaghan (1983) J. Comput. Phys. 52 (2) 374–389
            if ( this->m_artificial_viscosity && !issolid )
                {
                Scalar dotdvdx = dot(dv,dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar muij    = meanh*dotdvdx/(rsq+epssqr);
                    Scalar meanrho = Scalar(0.5)*(rhoi+rhoj);
                    avc = (-this->m_avalpha*this->m_c*muij+this->m_avbeta*muij*muij)/meanrho;
                    }
                }

            // Add contribution to fluid particle; pressure interaction force
            Scalar vijsqr = Vi*Vi+Vj*Vj;
            h_force.data[i].x -= vijsqr * ( temp0 + avc )* dwdr_r * dx.x;
            h_force.data[i].y -= vijsqr * ( temp0 + avc )* dwdr_r * dx.y;
            h_force.data[i].z -= vijsqr * ( temp0 + avc )* dwdr_r * dx.z;

            // Evaluate viscous interaction forces
            temp0 = this->m_mu * vijsqr * dwdr_r;
            h_force.data[i].x  += temp0 * dv.x;
            h_force.data[i].y  += temp0 * dv.y;
            h_force.data[i].z  += temp0 * dv.z;

            // Evaluate and add artificial stress part
            temp0 = 0.5 * vijsqr * dwdr_r;
            Scalar A1ij = ( A11i + A11j ) * dx.x + ( A12i + A12j ) * dx.y + ( A13i + A13j ) * dx.z;
            Scalar A2ij = ( A21i + A21j ) * dx.x + ( A22i + A22j ) * dx.y + ( A23i + A23j ) * dx.z;
            Scalar A3ij = ( A31i + A31j ) * dx.x + ( A32i + A32j ) * dx.y + ( A33i + A33j ) * dx.z;

            h_force.data[k].x += temp0 * A1ij; 
            h_force.data[k].y += temp0 * A2ij; 
            h_force.data[k].z += temp0 * A3ij; 

            // Evaluate background pressure contribution in aux2
            h_bpc.data[i].x -= vijsqr * this->m_eos->getTransportVelocityPressure()/mi * dwdr_r * dx.x;
            h_bpc.data[i].y -= vijsqr * this->m_eos->getTransportVelocityPressure()/mi * dwdr_r * dx.y;
            h_bpc.data[i].z -= vijsqr * this->m_eos->getTransportVelocityPressure()/mi * dwdr_r * dx.z;

            // Evaluate rate of change of density if CONTINUITY approach is used
            if ( this->m_density_method == DENSITYCONTINUITY )
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
                if ( !issolid && this->m_density_diffusion )
                    h_ratedpe.data[i].x -= (Scalar(2)*this->m_ddiff*meanh*this->m_c*mj*(rhoi/rhoj-Scalar(1))*dot(dx,dwdr_r*dx))/(rsq+epssqr);
                }

            } // Closing Neighbor Loop

        } // Closing Fluid Particle Loop

    this->m_timestep_list[5] = max_vel;
    // Add volumetric force (gravity)
    this->applyBodyForce(timestep, this->m_fluidgroup);

    }

/*! Compute forces definition
*/

template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlowTV<KT_, SET_>::computeForces(uint64_t timestep)
    {

    // start by updating the neighborlist
    this->m_nlist->compute(timestep);

    // This is executed once to initialize protected/private variables
    if (!this->m_params_set)
        {
        this->m_exec_conf->msg->error() << "sph.models.SinglePhaseFlowTV requires parameters to be set before run()"
            << std::endl;
        throw std::runtime_error("Error computing SinglePhaseFlowTV forces");
        }

    // m_solid_removed flag is set to False initially, so this 
    // only executes at timestep 0
    if (!this->m_solid_removed)
        {
        this->m_nlist->forceUpdate();
        this->m_nlist->compute(timestep);
        this->mark_solid_particles_toremove(timestep);
        this->m_solid_removed = true;
        }

    // Apply density renormalization if requested
    if ( this->m_shepard_renormalization && timestep % this->m_shepardfreq == 0 )
        {
        this->renormalize_density(timestep);
#ifdef ENABLE_MPI
         // Update ghost particle densities and pressures.
        this->update_ghost_density(timestep);
#endif
        }

    if (this->m_density_method == DENSITYSUMMATION)
    {
        this->compute_ndensity(timestep);
        // compute_particlenumberdensity(timestep);
    }

    // Compute fluid pressure based on m_eos;
    // Only working on the fluidgroup
    this->compute_pressure(timestep);

#ifdef ENABLE_MPI
    // Update ghost particle densities and pressures.
    this->update_ghost_density_pressure(timestep);
#endif

    // Compute particle pressures
    // Includes the computation of the density of solid particles
    // based on ficticios pressure p_i^\ast
    this->compute_noslip(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    this->update_ghost_aux1(timestep);
#endif

    // Execute the force computation
    // This includes the computation of the density if 
    // DENSITYCONTINUITY method is used
    forcecomputation(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    update_ghost_aux123(timestep);
#endif

    }

namespace detail 
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlowTV(pybind11::module& m, std::string name)
{
    pybind11::class_<SinglePhaseFlowTV<KT_, SET_>, SPHBaseClass<KT_, SET_> , std::shared_ptr<SinglePhaseFlowTV<KT_, SET_>>>(m, name.c_str()) 
        .def(pybind11::init< std::shared_ptr<SystemDefinition>,
                             std::shared_ptr<SmoothingKernel<KT_> >,
                             std::shared_ptr<StateEquation<SET_> >,
                             std::shared_ptr<nsearch::NeighborList>,
                             std::shared_ptr<ParticleGroup>,
                             std::shared_ptr<ParticleGroup>,
                             DensityMethod,
                             ViscosityMethod >())
        .def("setParams", &SinglePhaseFlowTV<KT_, SET_>::setParams)
        .def("getDensityMethod", &SinglePhaseFlowTV<KT_, SET_>::getDensityMethod)
        .def("setDensityMethod", &SinglePhaseFlowTV<KT_, SET_>::setDensityMethod)
        .def("getViscosityMethod", &SinglePhaseFlowTV<KT_, SET_>::getViscosityMethod)
        .def("setViscosityMethod", &SinglePhaseFlowTV<KT_, SET_>::setViscosityMethod)
        .def("setConstSmoothingLength", &SinglePhaseFlowTV<KT_, SET_>::setConstSmoothingLength)
        .def("computeSolidForces", &SinglePhaseFlowTV<KT_, SET_>::computeSolidForces)
        .def("activateArtificialViscosity", &SinglePhaseFlowTV<KT_, SET_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity", &SinglePhaseFlowTV<KT_, SET_>::deactivateArtificialViscosity)
        .def("activateDensityDiffusion", &SinglePhaseFlowTV<KT_, SET_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion", &SinglePhaseFlowTV<KT_, SET_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &SinglePhaseFlowTV<KT_, SET_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization", &SinglePhaseFlowTV<KT_, SET_>::deactivateShepardRenormalization)
        .def("setAcceleration", &SPHBaseClass<KT_, SET_>::setAcceleration)
        .def("setRCut", &SinglePhaseFlowTV<KT_, SET_>::setRCutPython)
        ;

    }

} // end namespace detail

//! Explicit template instantiations
template class PYBIND11_EXPORT SinglePhaseFlowTV<wendlandc2, linear>;
template class PYBIND11_EXPORT SinglePhaseFlowTV<wendlandc2, tait>;
template class PYBIND11_EXPORT SinglePhaseFlowTV<wendlandc4, linear>;
template class PYBIND11_EXPORT SinglePhaseFlowTV<wendlandc4, tait>;
template class PYBIND11_EXPORT SinglePhaseFlowTV<wendlandc6, linear>;
template class PYBIND11_EXPORT SinglePhaseFlowTV<wendlandc6, tait>;
template class PYBIND11_EXPORT SinglePhaseFlowTV<quintic, linear>;
template class PYBIND11_EXPORT SinglePhaseFlowTV<quintic, tait>;
template class PYBIND11_EXPORT SinglePhaseFlowTV<cubicspline, linear>;
template class PYBIND11_EXPORT SinglePhaseFlowTV<cubicspline, tait>;


namespace detail
{

    template void export_SinglePhaseFlowTV<wendlandc2, linear>(pybind11::module& m, std::string name = "SinglePFTV_WC2_L");
    template void export_SinglePhaseFlowTV<wendlandc2, tait>(pybind11::module& m, std::string name = "SinglePFTV_WC2_T");
    template void export_SinglePhaseFlowTV<wendlandc4, linear>(pybind11::module& m, std::string name = "SinglePFTV_WC4_L");
    template void export_SinglePhaseFlowTV<wendlandc4, tait>(pybind11::module& m, std::string name = "SinglePFTV_WC4_T");
    template void export_SinglePhaseFlowTV<wendlandc6, linear>(pybind11::module& m, std::string name = "SinglePFTV_WC6_L");
    template void export_SinglePhaseFlowTV<wendlandc6, tait>(pybind11::module& m, std::string name = "SinglePFTV_WC6_T");
    template void export_SinglePhaseFlowTV<quintic, linear>(pybind11::module& m, std::string name = "SinglePFTV_Q_L");
    template void export_SinglePhaseFlowTV<quintic, tait>(pybind11::module& m, std::string name = "SinglePFTV_Q_T");
    template void export_SinglePhaseFlowTV<cubicspline, linear>(pybind11::module& m, std::string name = "SinglePFTV_CS_L");
    template void export_SinglePhaseFlowTV<cubicspline, tait>(pybind11::module& m, std::string name = "SinglePFTV_CS_T");

} // end namespace detail
} // end namespace sph
} // end namespace hoomd

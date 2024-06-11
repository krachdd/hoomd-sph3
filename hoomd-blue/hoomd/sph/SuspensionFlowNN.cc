/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SuspensionFlowNN.h"

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
SuspensionFlowNN<KT_, SET_>::SuspensionFlowNN(std::shared_ptr<SystemDefinition> sysdef,
                                 std::shared_ptr<SmoothingKernel<KT_> > skernel,
                                 std::shared_ptr<StateEquation<SET_> > equationofstate,
                                 std::shared_ptr<nsearch::NeighborList> nlist,
                                 std::shared_ptr<ParticleGroup> fluidgroup,
                                 std::shared_ptr<ParticleGroup> solidgroup,
                                 std::shared_ptr<ParticleGroup> suspendedgroup,
                                 DensityMethod mdensitymethod,
                                 ViscosityMethod mviscositymethod)
    : SuspensionFlow<KT_, SET_>(sysdef, skernel, equationofstate, nlist, fluidgroup, solidgroup, suspendedgroup, mdensitymethod, mviscositymethod)
      {
        this->m_exec_conf->msg->notice(5) << "Constructing SuspensionFlowNN" << std::endl;

        // Set private attributes to default values
        this->m_const_slength = false;
        this->m_params_set = false;
        this->m_compute_solid_forces = false;
        this->m_compute_fluid_acceleration = true;
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
        this->constructTypeVectors(suspendedgroup,&this->m_suspendedtypes);

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
            for (unsigned int i = 0; i < this->m_suspendedtypes.size(); i++) {
                h_type_property_map.data[this->m_suspendedtypes[i]] |= SolidFluidTypeBit::SUSPENDED;
            }
        }

        // Init vectors for suspension model
        Scalar3 zeros3 = make_scalar3(0.0,0.0,0.0);
        Scalar4 zeros4 = make_scalar4(0.0,0.0,0.0,0.0);
        for (unsigned int i = 0; i < this->m_suspendedtypes.size(); ++i)
            {
            this->m_centerofmasses.push_back(zeros4);
            this->m_repulsiveforces.push_back(zeros3);
            this->m_velocities.push_back(zeros3);
            this->m_radii.push_back(0.0);
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

      }

/*! Destructor
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
SuspensionFlowNN<KT_, SET_>::~SuspensionFlowNN()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying SuspensionFlowNN" << std::endl;
    }


/*! \post Set model parameters
 */

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlowNN<KT_, SET_>::setParams(Scalar k, Scalar tau, Scalar m, Scalar rhoS, Scalar f0)
    {
    m_k = k;
    m_tau = tau;
    m_m = m;
    this->m_rhoS = rhoS;
    this->m_f0 = f0;
    if (m_k <= 0)
         {
         this->m_exec_conf->msg->error() << "sph.models.SinglePhaseFlowNN: Dynamic viscosity has to be a positive real number" << std::endl;
         throw std::runtime_error("Error initializing SinglePhaseFlowNN.");
         }
    if (m_tau <= 0)
         {
         this->m_exec_conf->msg->error() << "sph.models.SinglePhaseFlowNN: Tau has to be a positive real number" << std::endl;
         throw std::runtime_error("Error initializing SinglePhaseFlowNN.");
         }
    if (m_m <= 0)
         {
         this->m_exec_conf->msg->error() << "sph.models.SinglePhaseFlowNN: m value has to be a positive real number" << std::endl;
         throw std::runtime_error("Error initializing SinglePhaseFlowNN.");
         }
    if (this->m_rhoS <= 0)
         {
         this->m_exec_conf->msg->error() << "sph.models.SuspensionFlowNN: Solid density has to be a positive real number" << std::endl;
         throw std::runtime_error("Error initializing SuspensionFlowNN.");
         }
    if (this->m_f0 < 0)
        {
        this->m_exec_conf->msg->error() << "sph.models.SuspensionFlowNN: Spring stiffness has to be a positive real number" << std::endl;
        throw std::runtime_error("Error initializing SuspensionFlowNN.");
        }
    this->m_params_set = true;
    }

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlowNN<KT_, SET_>::update_ghost_aux4(uint64_t timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlowNN::Update Ghost aux1" << std::endl;

#ifdef ENABLE_MPI
    if (this->m_comm)
        {
        CommFlags flags(0);
        flags[comm_flag::tag] = 0;
        flags[comm_flag::position] = 0;
        flags[comm_flag::velocity] = 0;
        // flags[comm_flag::dpe] = 1;
        flags[comm_flag::density] = 0;
        flags[comm_flag::pressure] = 0;
        flags[comm_flag::energy] = 0;
        flags[comm_flag::auxiliary1] = 0; // ficticios velocity
        flags[comm_flag::auxiliary2] = 0; // repulsive force
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

/*! Perform viscosity computation
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlowNN<KT_, SET_>::compute_viscosity(uint64_t timestep)
    {
    // access the particle data
    ArrayHandle<Scalar>  h_density(this->m_pdata->getDensities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);

    ArrayHandle<Scalar3> h_vf(this->m_pdata->getAuxiliaries1(), access_location::host,access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_nn(this->m_pdata->getAuxiliaries4(), access_location::host,access_mode::readwrite);

    // access the neighbor list
    ArrayHandle<unsigned int> h_n_neigh(this->m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(this->m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<size_t> h_head_list(this->m_nlist->getHeadList(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_type_property_map(this->m_type_property_map, access_location::host, access_mode::read);

    // Check input data
    //assert(h_pos.data != NULL);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    unsigned int size;
    size_t myHead;

    // Particle loop
    // For each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_fluidgroup->getMemberIndex(group_idx);

        // TODO: Aus Schleife raus initialisieren hier nur null setzen
        // Initialize velocity gradient
        Scalar L[9] = {0};
        // Initialize shear strain rate
        Scalar D[9] = {0};
        // Initialize L2-norm of shear rate
        Scalar norm_shear_rate = 0.0;
        // Initialize variable to compute current viscosity
        Scalar mu_i = 0.0;
        // Initialize variable to compute current shear stress
        Scalar tau_i = 0.0;

        // Access the particle's position, velocity, mass and type
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;

        Scalar3 vi;
        vi.x = h_velocity.data[i].x;
        vi.y = h_velocity.data[i].y;
        vi.z = h_velocity.data[i].z;
        
        // Loop over all of the neighbors of this particle
        // const long unsigned int myHead = h_head_list.data[i];
        // const unsigned int size = (unsigned int)h_n_neigh.data[i];
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

            // Compute distance vector (FLOPS: 3)
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

            // Determine neighbor type
            bool issolid = checksolid(h_type_property_map.data, h_pos.data[k].w);
            bool issuspended = checksuspended(h_type_property_map.data, h_pos.data[k].w);

            // Access neighbor velocity; depends on fluid or fictitious solid particle
            Scalar3 vj  = make_scalar3(0.0, 0.0, 0.0);
            Scalar mj   = h_velocity.data[k].w;


            if (issolid || issuspended)
                continue;

            vj.x = h_velocity.data[k].x;
            vj.y = h_velocity.data[k].y;
            vj.z = h_velocity.data[k].z;

            Scalar rhoj = h_density.data[k];
            Scalar Vj   = mj / rhoj;    // =============>>> USE NUMBER DENSITY ???
            // HINT:
            // mass_density = m_i * \sum_j wij ==> h_dpe.data[i].x = h_dpe.data[i].x * h_velocity.data[i].w;

            // Compute velocity difference
            Scalar3 dv;
            dv.x = vj.x - vi.x;
            dv.y = vj.y - vi.y;
            dv.z = vj.z - vi.z;

            // Calculate absolute and normalized distance
            Scalar r = sqrt(rsq);

            // Mean smoothing length and denominator modifier
            Scalar meanh  = this->m_const_slength ? this->m_ch : Scalar(0.5)*(h_h.data[i]+h_h.data[k]);
            Scalar eps    = Scalar(0.1)*meanh;

            // Kernel function derivative evaluation
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);

            // Evaluate grad u (zeroth-order consistency --> doesn't need local momenmtum conservation)
            L[0] += Vj * dv.x * dx.x * dwdr_r;
            L[1] += Vj * dv.x * dx.y * dwdr_r;
            L[2] += Vj * dv.x * dx.z * dwdr_r;

            L[3] += Vj * dv.y * dx.x * dwdr_r;
            L[4] += Vj * dv.y * dx.y * dwdr_r;
            L[5] += Vj * dv.y * dx.z * dwdr_r;

            L[6] += Vj * dv.z * dx.x * dwdr_r;
            L[7] += Vj * dv.z * dx.y * dwdr_r;
            L[8] += Vj * dv.z * dx.z * dwdr_r;

            // close neighbor loop
            }

        // COMPUTATION VARIANT 1: Following Rütten (Verallgemeinerte newtonsche Fluide, Springer 2019)
        // Evaluate shear strain rate tensor D = 1/2*(grad u + grad^T u) ( Sym. part of velocity gradient )
        D[0] = Scalar(0.5) * ( L[0] + L[0] );
        D[1] = Scalar(0.5) * ( L[1] + L[3] );
        D[2] = Scalar(0.5) * ( L[2] + L[6] );

        D[3] = Scalar(0.5) * ( L[3] + L[1] );
        D[4] = Scalar(0.5) * ( L[4] + L[4] );
        D[5] = Scalar(0.5) * ( L[5] + L[7] );

        D[6] = Scalar(0.5) * ( L[6] + L[2] );
        D[7] = Scalar(0.5) * ( L[7] + L[5] );
        D[8] = Scalar(0.5) * ( L[8] + L[8] );

        // Evaluate trace of D*D  (tr(D^2))
        norm_shear_rate = D[0] * D[0] + D[4] * D[4] + D[8] * D[8] + 2.0 * D[1] * D[1] + 2.0 * D[2] * D[2] + 2.0 * D[5] * D[5];
        norm_shear_rate = sqrt(Scalar(2.0) * norm_shear_rate);

        // Save norm of shear rate in aux3.x
        h_nn.data[i].x = norm_shear_rate;

        // // TODO: Eventuell threshold groesser als 0 für shear rate und mu0, das ist quatsch hier -- jetzt geaendert
        // if ( norm_shear_rate == Scalar(0.0) ) 
        //     {
        //     mu_i  = this->m_k+m_tau*m_m;
        //     tau_i = 0.0;
        //     } 
        // else 
        //     {
        //     mu_i  = this->m_k + ( this->m_tau / norm_shear_rate ) * ( Scalar(1.0) - exp( -(norm_shear_rate * this->m_m) ) );
        //     tau_i = this->m_k * norm_shear_rate + this->m_tau * ( Scalar(1.0) - exp( -(norm_shear_rate * this->m_m) ) ) ;
        //     }


        // TODO: Eventuell threshold groesser als 0 für shear rate und mu0, das ist quatsch hier -- jetzt geaendert
        if ( norm_shear_rate == Scalar(0.0) ) 
            {
            mu_i  = this->m_k+m_tau*m_m;
            tau_i = 0.0;
            } 
        else 
            {
            mu_i  = this->m_k + ( this->m_tau / norm_shear_rate ) * ( Scalar(1.0) - exp( -(norm_shear_rate * this->m_m) ) );
            tau_i = this->m_k * norm_shear_rate + this->m_tau * ( Scalar(1.0) - exp( -(norm_shear_rate * this->m_m) ) ) ;
            }

        if ( tau_i <= m_tau)
            {
            mu_i   = this->m_k+m_tau*m_m;
            }

        // Save norm of shear rate in aux3.x
        h_nn.data[i].x = norm_shear_rate;
        // Save viscosity in aux3.y
        h_nn.data[i].y = mu_i;
        // Save viscosity in aux3.z
        h_nn.data[i].z = tau_i;

        }
    }


template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlowNN<KT_, SET_>::forcecomputation(uint64_t timestep)
    {

    if ( this->m_density_method == DENSITYSUMMATION )
        this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlowNN::Forces using SUMMATION approach " << this->m_density_method << endl;
    else if ( this->m_density_method == DENSITYCONTINUITY )
        this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlowNN::Forces using CONTINUITY approach " << this->m_density_method << endl;

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
    ArrayHandle<Scalar3> h_nn(this->m_pdata->getAuxiliaries4(), access_location::host,access_mode::read);
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

    // Local variable to store mean viscosity
    Scalar meanmu = 0;

    // maximum velocity variable for adaptive timestep
    double max_vel = 0.0;

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
        Scalar Pi = h_pressure.data[i];

        // Read particle i density and volume
        Scalar rhoi = h_density.data[i];
        Scalar Vi   = mi / rhoi;

        // // Total velocity of particle
        Scalar vi_total = sqrt((vi.x * vi.x) + (vi.y * vi.y) + (vi.z * vi.z));

        // Read particle viscosity
        Scalar mui  = h_nn.data[i].y;

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
            if (issolid || issuspended)
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

            // Scalar Vj = Scalar(0.0);

            Scalar Vj = mj/rhoj;

            // Read particle k viscosity
            Scalar muj  = h_nn.data[k].y;


            // // Issolid wohl nicht noetig wenn nur BC - Masse m bleibt immer gleich deshalb mi
            // if (issolid || issuspended)
            //     {
            //     Vj = mi / rhoj;
            //     mj = mi;
            //     }
            // else
            //     {
            //     Vj = mj / rhoj; 
            //     }

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

            // Evaluate inter-particle pressure forces
            //temp0 = -((mi*mj)/(rhoj*rhoi))*(Pi+Pj);
            //temp0 = -Vi*Vj*( Pi + Pj );
            //temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);
            //temp0 = -mi*mj*( Pi/(rhoi*rhoj) + Pj/(rhoj*rhoj) );
            // if ( m_density_method == DENSITYSUMMATION )
            // {
            //     temp0 = -(Vi*Vi+Vj*Vj)*((rhoj*Pi+rhoi*Pj)/(rhoi+rhoj)); 
            // }
            // // Hier noch Anpassung noetig wegen Masse
            // else if ( m_density_method == DENSITYCONTINUITY) 
            // { 
            //     temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);
            // }

            temp0 = -(Vi*Vi+Vj*Vj)*((rhoj*Pi+rhoi*Pj)/(rhoi+rhoj)); 


            // Optionally add artificial viscosity
            // Monaghan (1983) J. Comput. Phys. 52 (2) 374–389
            if ( this->m_artificial_viscosity && (!issolid || !issuspended))
                {
                Scalar dotdvdx = dot(dv,dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar muij    = meanh*dotdvdx/(rsq+epssqr);
                    Scalar meanrho = Scalar(0.5)*(rhoi+rhoj);
                    temp0 += mi*mj*(this->m_avalpha*this->m_c*muij+this->m_avbeta*muij*muij)/meanrho;
                    }
                }

            // Add contribution to fluid particle
            h_force.data[i].x += temp0*dwdr_r*dx.x;
            h_force.data[i].y += temp0*dwdr_r*dx.y;
            h_force.data[i].z += temp0*dwdr_r*dx.z;

            // Evaluate viscous interaction forces
            // meanmu = (Scalar(2.0) * mui * muj) / (mui + muj);
            // Evaluate viscous interaction forces
            if (!issolid && !issuspended)
            {
                meanmu = (Scalar(2.0) * mui * muj) / (mui + muj);
            }
            else
            {
                meanmu = mui;
            }

            temp0 = meanmu * (Vi*Vi+Vj*Vj) * dwdr_r;
            h_force.data[i].x  += temp0*dv.x;
            h_force.data[i].y  += temp0*dv.y;
            h_force.data[i].z  += temp0*dv.z;

            // Evaluate rate of change of density if CONTINUITY approach is used
            if ( this->m_density_method == DENSITYCONTINUITY )
                {
                if ( issolid || issuspended)
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
                if ( (!issolid || !issuspended) && this->m_density_diffusion )
                    h_ratedpe.data[i].x -= (Scalar(2)*this->m_ddiff*meanh*this->m_c*mj*(rhoi/rhoj-Scalar(1))*dot(dx,dwdr_r*dx))/(rsq+epssqr);
                }

            } // Closing Neighbor Loop

        } // Closing Fluid Particle Loop

    this->m_timestep_list[5] = max_vel;
    // Add volumetric force (gravity)
    // Add volumetric force (gravity)
    if ( this->m_compute_fluid_acceleration)
        {
        this->applyBodyForce(timestep, this->m_fluidgroup);
        }
    }


template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlowNN<KT_, SET_>::forcecomputationsolids(uint64_t timestep)
    {

    if ( this->m_density_method == DENSITYSUMMATION )
        this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlowNN::Forces using SUMMATION approach " << this->m_density_method << endl;
    else if ( this->m_density_method == DENSITYCONTINUITY )
        this->m_exec_conf->msg->notice(7) << "Computing SuspensionFlowNN::Forces using CONTINUITY approach " << this->m_density_method << endl;

    // Grab handles for particle data
    // Access mode overwrite implies that data does not need to be read in
    ArrayHandle<Scalar4> h_force(this->m_force,access_location::host, access_mode::readwrite);

    // Check input data, can be omitted if need be
    assert(h_force.data);

    // // Zero data before force calculation
    // memset((void*)h_force.data,0,sizeof(Scalar4)*this->m_force.getNumElements());

    // access the particle data
    ArrayHandle<Scalar4> h_pos(this->m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_velocity(this->m_pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar>  h_density(this->m_pdata->getDensities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar>  h_pressure(this->m_pdata->getPressures(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_vf(this->m_pdata->getAuxiliaries1(), access_location::host,access_mode::read);
    ArrayHandle<Scalar3> h_nn(this->m_pdata->getAuxiliaries4(), access_location::host,access_mode::read);
    ArrayHandle<Scalar>  h_h(this->m_pdata->getSlengths(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_conforce(this->m_pdata->getAuxiliaries2(), access_location::host,access_mode::read);

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

    // Local variable to store mean viscosity
    Scalar meanmu = 0;

    // maximum velocity variable for adaptive timestep
    double max_vel = 0.0;

    // for each suspended particle
    unsigned int group_size_suspended = this->m_suspendedgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size_suspended; group_idx++)
        {
        // Read particle index
        unsigned int i = this->m_suspendedgroup->getMemberIndex(group_idx);

        h_force.data[i].x = Scalar(0.0);
        h_force.data[i].y = Scalar(0.0);
        h_force.data[i].z = Scalar(0.0);

        //unsigned int typeid_i = __scalar_as_int(h_pos.data[i].w);

        // Access the particle's position, velocity, mass and type
        Scalar3 pi;
        pi.x = h_pos.data[i].x;
        pi.y = h_pos.data[i].y;
        pi.z = h_pos.data[i].z;

        // Use ficitious velocity here for momentum conservation
        Scalar3 vi;
        // vi.x = h_velocity.data[i].x;
        // vi.y = h_velocity.data[i].y;
        // vi.z = h_velocity.data[i].z;
        vi.x = h_vf.data[i].x;
        vi.y = h_vf.data[i].y;
        vi.z = h_vf.data[i].z;

        Scalar mi = h_velocity.data[i].w;

        // Read particle i pressure
        Scalar Pi = h_pressure.data[i];

        // Read particle i density and volume
        Scalar rhoi = h_density.data[i];
        // Commented since using fluid mass for force calculation
        Scalar Vi   = mi / rhoi;

        // // Total velocity of particle
        Scalar vi_total = sqrt((vi.x * vi.x) + (vi.y * vi.y) + (vi.z * vi.z));

        // Properties needed for adaptive timestep
        if (i == 0) { max_vel = vi_total; }
        else if (vi_total > max_vel) { max_vel = vi_total; }

        // Read particle viscosity
        //Scalar mui  = h_nn.data[i].y;

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

            //unsigned int typeid_j = __scalar_as_int(h_pos.data[k].w);

            // TODO: Check if solid is valid but not suspended
            // Determine neighbor type
            bool issolid = checksolid(h_type_property_map.data, h_pos.data[k].w);
            bool issuspended = checksuspended(h_type_property_map.data, h_pos.data[k].w);
            // if (issuspended && typeid_i == typeid_j)
            //     continue;
            if (issuspended || issuspended)
                continue;

            // Compute distance vector (FLOPS: 3)
            // Scalar3 dx = pi - pj;
            Scalar3 dx;
            dx.x = pi.x - pj.x;
            dx.y = pi.y - pj.y;
            dx.z = pi.z - pj.z;
            // // Quatsch Richtung aendern, andere VZ bleiben gleich
            // dx.x = pj.x - pi.x;
            // dx.y = pj.y - pi.y;
            // dx.z = pj.z - pi.z;

            // Apply periodic boundary conditions (FLOPS: 9)
            dx = box.minImage(dx);

            // Calculate squared distance (FLOPS: 5)
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, skip this loop
            if ( this->m_const_slength && rsq > this->m_rcutsq )
                continue;

            // Access neighbor velocity; depends on fluid or fictitious solid particle
            Scalar3 vj  = make_scalar3(0.0, 0.0, 0.0);
            // // TODO: Hier wohl auch suspended?!
            // if ( issolid || issuspended)
            //     {
            //     vj.x = h_vf.data[k].x;
            //     vj.y = h_vf.data[k].y;
            //     vj.z = h_vf.data[k].z;
            //     }
            // else
            //     {
            //     vj.x = h_velocity.data[k].x;
            //     vj.y = h_velocity.data[k].y;
            //     vj.z = h_velocity.data[k].z;
            //     }

            vj.x = h_velocity.data[k].x;
            vj.y = h_velocity.data[k].y;
            vj.z = h_velocity.data[k].z;

            // // Vi hier damit identisch mit vorher fluiddensity nutzen
            // Scalar Vi   = mj / rhoi;
            // Scalar rhoj = h_density.data[k];
            // Scalar Vj   = mj / rhoj;
            // Scalar mi   = mj;

            Scalar rhoj = h_density.data[k];
            Scalar mj   = h_velocity.data[k].w;
            Scalar Vj = mj/rhoj;

            // Read particle k viscosity
            Scalar muj  = h_nn.data[k].y;

            // Read particle k pressure
            Scalar Pj = h_pressure.data[k];

            // Compute velocity difference
            // Richtung aendern, andere VZ bleiben gleich
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
            // Radius bleibt gleich -  da Skalar und immer positiv
            Scalar dwdr   = this->m_skernel->dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);

            // Evaluate inter-particle pressure forces
            //temp0 = -((mi*mj)/(rhoj*rhoi))*(Pi+Pj);
            //temp0 = -Vi*Vj*( Pi + Pj );
            //temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);
            //temp0 = -mi*mj*( Pi/(rhoi*rhoj) + Pj/(rhoj*rhoj) );
            // if ( m_density_method == DENSITYSUMMATION )
            // {
            //     temp0 = -(Vi*Vi+Vj*Vj)*((rhoj*Pi+rhoi*Pj)/(rhoi+rhoj)); 
            // }
            // else if ( m_density_method == DENSITYCONTINUITY) 
            // { 
            //     temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);
            // }

            temp0 = -(Vi*Vi+Vj*Vj)*((rhoj*Pi+rhoi*Pj)/(rhoi+rhoj)); 


            // Add contribution to suspended particle
            h_force.data[i].x += temp0*dwdr_r*dx.x;
            h_force.data[i].y += temp0*dwdr_r*dx.y;
            h_force.data[i].z += temp0*dwdr_r*dx.z;

            // // // Add contribution to solid particle
            // // if ( issolid && m_compute_solid_forces )
            // //     {
            // //     h_force.data[k].x -= (mj/mi)*temp0*dwdr_r*dx.x;
            // //     h_force.data[k].y -= (mj/mi)*temp0*dwdr_r*dx.y;
            // //     h_force.data[k].z -= (mj/mi)*temp0*dwdr_r*dx.z;
            // //     }

            // // Evaluate viscous interaction forces
            // // meanmu = (Scalar(2.0) * mui * muj) / (mui + muj);
            // // Evaluate viscous interaction forces
            // if (!issolid && !issuspended)
            // {
            //     meanmu = (Scalar(2.0) * mui * muj) / (mui + muj);
            // }
            // else
            // {
            //     meanmu = mui;
            // }
            meanmu = muj;

            temp0 = meanmu * (Vi*Vi+Vj*Vj) * dwdr_r;
            h_force.data[i].x  += temp0*dv.x;
            h_force.data[i].y  += temp0*dv.y;
            h_force.data[i].z  += temp0*dv.z;

            // // // Add contribution to solid particle
            // // if ( issolid && m_compute_solid_forces )
            // //     {
            // //     h_force.data[k].x -= (mj/mi)*temp0*dv.x;
            // //     h_force.data[k].y -= (mj/mi)*temp0*dv.y;
            // //     h_force.data[k].z -= (mj/mi)*temp0*dv.z;
            // //     }

            // // Add contribution to suspended particle
            // h_force.data[i].x += (mj/mi)*temp0*dwdr_r*dx.x;
            // h_force.data[i].y += (mj/mi)*temp0*dwdr_r*dx.y;
            // h_force.data[i].z += (mj/mi)*temp0*dwdr_r*dx.z;

            // // Evaluate viscous interaction forces
            // temp0 = m_mu * (Vi*Vi+Vj*Vj) * dwdr_r;
            // h_force.data[i].x  += (mj/mi)*temp0*dv.x;
            // h_force.data[i].y  += (mj/mi)*temp0*dv.y;
            // h_force.data[i].z  += (mj/mi)*temp0*dv.z;


            } // Closing Neighbor Loop

        // Add contribution of contact force to total force
        h_force.data[i].x += h_conforce.data[i].x;
        h_force.data[i].y += h_conforce.data[i].y;
        h_force.data[i].z += h_conforce.data[i].z;

        } // Closing Suspended Particle Loop

    this->m_timestep_list[5] = max_vel;

    // Add volumetric force (gravity)
    if ( this->m_compute_solid_forces )
        this->applyBodyForce(timestep, this->m_suspendedgroup);
     }

/*! Compute forces definition
*/

template<SmoothingKernelType KT_,StateEquationType SET_>
void SuspensionFlowNN<KT_, SET_>::computeForces(uint64_t timestep)
    {

    // start by updating the neighborlist
    this->m_nlist->compute(timestep);

    // This is executed once to initialize protected/private variables
    if (!this->m_params_set)
        {
        this->m_exec_conf->msg->error() << "sph.models.SuspensionFlowNN requires parameters to be set before run()"
            << std::endl;
        throw std::runtime_error("Error computing SuspensionFlowNN forces");
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

    // Compute particle pressures
    // Includes the computation of the density of suspended  particles
    // based on ficticios pressure p_i^\ast
    this->compute_noslipsuspended(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    this->update_ghost_aux1(timestep);
#endif

    // // Compute center of mass of solid bodies
    // compute_Centerofmasses(timestep, true);
    // Compute center of mass of solid bodies
    this->compute_Centerofmasses(timestep);

    // Compute equivalent radi of solid bodies ============= >>> evtl als Vektor speichern (x,y,z)
    this->compute_equivalentRadii(timestep);

    // Compute repulsive forces
    this->compute_repulsiveForce(timestep);

    // communicate repulsive forces
#ifdef ENABLE_MPI
    // Update ghost particles
    this->update_ghost_aux2(timestep);
#endif

    // Compute fluid viscosity
    this->compute_viscosity(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    update_ghost_aux4(timestep);
#endif

    // Execute the force computation
    // This includes the computation of the density if 
    // DENSITYCONTINUITY method is used
    forcecomputation(timestep);

    // Execute the force computation
    // This includes the computation of the density if 
    // DENSITYCONTINUITY method is used
    forcecomputationsolids(timestep);

    }

namespace detail 
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SuspensionFlowNN(pybind11::module& m, std::string name)
{
    pybind11::class_<SuspensionFlowNN<KT_, SET_>, SPHBaseClass<KT_, SET_> , std::shared_ptr<SuspensionFlowNN<KT_, SET_>>>(m, name.c_str()) 
        .def(pybind11::init< std::shared_ptr<SystemDefinition>,
                             std::shared_ptr<SmoothingKernel<KT_> >,
                             std::shared_ptr<StateEquation<SET_> >,
                             std::shared_ptr<nsearch::NeighborList>,
                             std::shared_ptr<ParticleGroup>,
                             std::shared_ptr<ParticleGroup>,
                             std::shared_ptr<ParticleGroup>,
                             DensityMethod,
                             ViscosityMethod >())
        .def("setParams", &SuspensionFlowNN<KT_, SET_>::setParams)
        .def("getDensityMethod", &SuspensionFlowNN<KT_, SET_>::getDensityMethod)
        .def("setDensityMethod", &SuspensionFlowNN<KT_, SET_>::setDensityMethod)
        .def("getViscosityMethod", &SuspensionFlowNN<KT_, SET_>::getViscosityMethod)
        .def("setViscosityMethod", &SuspensionFlowNN<KT_, SET_>::setViscosityMethod)
        .def("setConstSmoothingLength", &SuspensionFlowNN<KT_, SET_>::setConstSmoothingLength)
        .def("computeSolidForces", &SuspensionFlowNN<KT_, SET_>::computeSolidForces)
        .def("deactivateFluidAcceleration", &SuspensionFlowNN<KT_, SET_>::deactivateFluidAcceleration)
        .def("activateArtificialViscosity", &SuspensionFlowNN<KT_, SET_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity", &SuspensionFlowNN<KT_, SET_>::deactivateArtificialViscosity)
        .def("activateDensityDiffusion", &SuspensionFlowNN<KT_, SET_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion", &SuspensionFlowNN<KT_, SET_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &SuspensionFlowNN<KT_, SET_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization", &SuspensionFlowNN<KT_, SET_>::deactivateShepardRenormalization)
        .def("setAcceleration", &SPHBaseClass<KT_, SET_>::setAcceleration)
        .def("setRCut", &SuspensionFlowNN<KT_, SET_>::setRCutPython)
        ;

    }

} // end namespace detail

//! Explicit template instantiations
template class PYBIND11_EXPORT SuspensionFlowNN<wendlandc2, linear>;
template class PYBIND11_EXPORT SuspensionFlowNN<wendlandc2, tait>;
template class PYBIND11_EXPORT SuspensionFlowNN<wendlandc4, linear>;
template class PYBIND11_EXPORT SuspensionFlowNN<wendlandc4, tait>;
template class PYBIND11_EXPORT SuspensionFlowNN<wendlandc6, linear>;
template class PYBIND11_EXPORT SuspensionFlowNN<wendlandc6, tait>;
template class PYBIND11_EXPORT SuspensionFlowNN<quintic, linear>;
template class PYBIND11_EXPORT SuspensionFlowNN<quintic, tait>;
template class PYBIND11_EXPORT SuspensionFlowNN<cubicspline, linear>;
template class PYBIND11_EXPORT SuspensionFlowNN<cubicspline, tait>;


namespace detail
{

    template void export_SuspensionFlowNN<wendlandc2, linear>(pybind11::module& m, std::string name = "SuspensionFNN_WC2_L");
    template void export_SuspensionFlowNN<wendlandc2, tait>(pybind11::module& m, std::string name = "SuspensionFNN_WC2_T");
    template void export_SuspensionFlowNN<wendlandc4, linear>(pybind11::module& m, std::string name = "SuspensionFNN_WC4_L");
    template void export_SuspensionFlowNN<wendlandc4, tait>(pybind11::module& m, std::string name = "SuspensionFNN_WC4_T");
    template void export_SuspensionFlowNN<wendlandc6, linear>(pybind11::module& m, std::string name = "SuspensionFNN_WC6_L");
    template void export_SuspensionFlowNN<wendlandc6, tait>(pybind11::module& m, std::string name = "SuspensionFNN_WC6_T");
    template void export_SuspensionFlowNN<quintic, linear>(pybind11::module& m, std::string name = "SuspensionFNN_Q_L");
    template void export_SuspensionFlowNN<quintic, tait>(pybind11::module& m, std::string name = "SuspensionFNN_Q_T");
    template void export_SuspensionFlowNN<cubicspline, linear>(pybind11::module& m, std::string name = "SuspensionFNN_CS_L");
    template void export_SuspensionFlowNN<cubicspline, tait>(pybind11::module& m, std::string name = "SuspensionFNN_CS_T");

} // end namespace detail
} // end namespace sph
} // end namespace hoomd

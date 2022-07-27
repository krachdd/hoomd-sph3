// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "SinglePhaseFlowGPU.h"

#include <boost/python.hpp>
using namespace boost::python;

#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
using namespace boost;

using namespace std;


/*! Constructor
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
SinglePhaseFlowGPU<KT_, SET_>::SinglePhaseFlowGPU(boost::shared_ptr<SystemDefinition> sysdef,
                                       boost::shared_ptr<SmoothingKernel<KT_> > skernel,
                                       boost::shared_ptr<StateEquation<SET_> > equationofstate,
                                       boost::shared_ptr<nsearch::NeighborList> nlist,
                                       boost::shared_ptr<ParticleGroup> fluidgroup,
                                       boost::shared_ptr<ParticleGroup> solidgroup,
                                       DensityMethod mdensitymethod,
                                       ViscosityMethod mviscositymethod)
    : SinglePhaseFlow<KT_, SET_>(sysdef,skernel,equationofstate,nlist, fluidgroup, solidgroup, mdensitymethod, mviscositymethod)
      {
        this->m_exec_conf->msg->notice(5) << "Constructing SinglePhaseFlow GPU" << std::endl;
      }

/*! Destructor
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
SinglePhaseFlowGPU<KT_, SET_>::~SinglePhaseFlowGPU()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying SinglePhaseFlow GPU" << std::endl;
    }


/*! Computes all log quantities
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlowGPU<KT_, SET_>::computeProperties()
    {
    if (this->m_prof)
        this->m_prof->push("SinglePhaseFlowlogger");

    ArrayHandle<Scalar4> d_velocity(this->m_pdata->getVelocities(), access_location::device, access_mode::read);
    
    //double fluid_vel_x_sum  = 0.0;
    //double fluid_vel_y_sum  = 0.0;
    //double fluid_vel_z_sum  = 0.0;
    double fluid_prtl = 0;
    Scalar3 *sum;
    Scalar3 h_sum;
    cudaMalloc((void**)&sum, 3*sizeof(Scalar));
    cudaMemset(sum, 0, 3*sizeof(Scalar));
    
    // for each fluid particle
    
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    fluid_prtl = group_size;
    ArrayHandle< unsigned int > d_index_array(this->m_fluidgroup->getIndexArray(), access_location::device, access_mode::read);

    gpu_spf_compute_properties(d_velocity.data, d_index_array.data, group_size, sum);
    cudaMemcpy(&h_sum, sum, 3*sizeof(Scalar), cudaMemcpyDeviceToHost);
    
    ArrayHandle<Scalar> h_properties(this->m_properties, access_location::host, access_mode::overwrite);
    h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_x] = Scalar(h_sum.x);
    h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_y] = Scalar(h_sum.y);
    h_properties.data[singlephaseflow_logger_index::sum_fluid_velocity_z] = Scalar(h_sum.z);
    h_properties.data[singlephaseflow_logger_index::total_fluid_particles]      = Scalar(fluid_prtl);

#ifdef ENABLE_MPI
    this->m_properties_reduced = !this->m_pdata->getDomainDecomposition();
#endif

    if (this->m_prof)
        this->m_prof->pop();
    }


#ifdef ENABLE_MPI_CUDA

template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlow<KT_, SET_>::reduceProperties()
    {
    if (m_properties_reduced) return;

    // reduce properties
    ArrayHandle<Scalar> d_properties(m_properties, access_location::device, access_mode::readwrite);
    MPI_Allreduce(MPI_IN_PLACE, d_properties.data, singlephaseflow_logger_index::num_quantities, MPI_HOOMD_SCALAR,
            MPI_SUM, this->m_exec_conf->getMPICommunicator());

    m_properties_reduced = true;
    }
#endif

/*! Perform number density computation
 * This method computes and stores
     - the number density based mass density ( rho_i = m_i * \sum w_ij ) for fluid particles
       if the SUMMATION approach is being used.
     - the zeroth order normalization constant for solid particles
   in the x-component of the d_dpe Array.
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlowGPU<KT_, SET_>::compute_ndensity(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Number Density" << std::endl;

    if (this->m_prof)
        this->m_prof->push("SinglePhaseFlowNDensity");

    // Grab handles for particle data
    ArrayHandle<Scalar3> d_dpe(this->m_pdata->getDPEs(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_pos(this->m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_velocity(this->m_pdata->getVelocities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h(this->m_pdata->getSlengths(), access_location::device, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(this->m_nlist->getNListArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_head_list(this->m_nlist->getHeadList(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_type_property_map(this->m_type_property_map, access_location::device, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Particle loop
    gpu_spf_compute_ndensity(d_pos.data,
                         d_dpe.data,
                         d_velocity.data,
                         d_h.data,
                         d_n_neigh.data,
                         d_nlist.data,
                         d_head_list.data,
                         this->m_pdata->getN(),
                         box,
                         this->m_ch,
                         this->m_rcutsq,
                         this->m_rho0,
                         this->m_const_slength,
                         this->m_skernel.get(),
                         d_type_property_map.data,
                         this->m_density_method == DENSITYSUMMATION
    );

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Perform pressure computation
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlowGPU<KT_, SET_>::compute_pressure(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Pressure" << std::endl;

    if (this->m_prof)
        this->m_prof->push("SinglePhaseFlowPressure");

    // Define ArrayHandles
    ArrayHandle<Scalar3> d_dpe(this->m_pdata->getDPEs(), access_location::device, access_mode::readwrite);

    // For each fluid particle
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    ArrayHandle< unsigned int > d_index_array(this->m_fluidgroup->getIndexArray(), access_location::device, access_mode::read);

    gpu_spf_compute_pressure(d_dpe.data,
                             d_index_array.data,
                             group_size,
                             this->m_eos.get());

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Compute fictitious solid particle properties
    This method updates fictitious solid particle pressures and velocities to account for
    no-slip boundary conditions. Method follows Adami et. al. (2012).
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlowGPU<KT_, SET_>::compute_noslip(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::NoSlip" << std::endl;

    if (this->m_prof)
        this->m_prof->push("SinglePhaseFlowNoSlip");

    // Grab handles for particle and neighbor data
    ArrayHandle<Scalar3> d_vf(this->m_pdata->getAuxiliaries1(), access_location::device,access_mode::readwrite);
    ArrayHandle<Scalar4> d_pos(this->m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_dpe(this->m_pdata->getDPEs(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_velocity(this->m_pdata->getVelocities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_accel(this->m_pdata->getAccelerations(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h(this->m_pdata->getSlengths(), access_location::device, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(this->m_nlist->getNListArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_head_list(this->m_nlist->getHeadList(), access_location::device, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // For all solid particles
    ArrayHandle< unsigned int > d_index_array(this->m_solidgroup->getIndexArray(), access_location::device, access_mode::read);
    unsigned int group_size = this->m_solidgroup->getNumMembers();
    ArrayHandle<unsigned int> d_type_property_map(this->m_type_property_map, access_location::device, access_mode::read);

    gpu_spf_compute_noslip(d_vf.data,
                       d_pos.data,
                       d_dpe.data,
                       d_velocity.data,
                       d_accel.data,
                       d_h.data,
                       d_n_neigh.data,
                       d_nlist.data,
                       d_head_list.data,
                       d_index_array.data,
                       group_size,
                       box,
                       this->m_rcutsq,
                       this->m_ch,
                       this->m_rho0,
                       this->m_const_slength,
                       this->getAcceleration(timestep),
                       this->m_eos.get(),
                       this->m_skernel.get(),
                       d_type_property_map.data,
                       this->m_body_acceleration
    );

    if (this->m_prof)
        this->m_prof->pop();
    }


/*! Perform Shepard density renormalization
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlowGPU<KT_, SET_>::renormalize_density(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Density renormalization" << std::endl;

    if (this->m_prof)
        this->m_prof->push("SinglePhaseFlowDensityRenormalization");

    // Grab handles for particle data
    ArrayHandle<Scalar3> d_dpe(this->m_pdata->getDPEs(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_pos(this->m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_velocity(this->m_pdata->getVelocities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h(this->m_pdata->getSlengths(), access_location::device, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(this->m_nlist->getNListArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_head_list(this->m_nlist->getHeadList(), access_location::device, access_mode::read);

    auto tmp_pde = this->m_pdata->getDPEs();
    ArrayHandle<Scalar3> d_dpe_old(tmp_pde, access_location::device, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Particle loop
    // For each fluid particle
    ArrayHandle< unsigned int > d_index_array(this->m_fluidgroup->getIndexArray(), access_location::device, access_mode::read);
    unsigned int group_size = this->m_fluidgroup->getNumMembers();

    gpu_spf_compute_ndensityrenormalization(d_pos.data,
                                        d_dpe.data,
                                        d_dpe_old.data,
                                        d_velocity.data,
                                        d_h.data,
                                        d_n_neigh.data,
                                        d_nlist.data,
                                        d_head_list.data,
                                        d_index_array.data,
                                        group_size,
                                        box,
                                        this->m_rcutsq,
                                        this->m_ch,
                                        this->m_const_slength,
                                        this->m_skernel.get()
    );

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Perform force computation
 */
template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlowGPU<KT_, SET_>::forcecomputation(unsigned int timestep)
    {

    if ( this->m_density_method == DENSITYSUMMATION )
        this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Forces using SUMMATION approach" << endl;
    else if ( this->m_density_method == DENSITYCONTINUITY )
        this->m_exec_conf->msg->notice(7) << "Computing SinglePhaseFlow::Forces using CONTINUITY approach" << endl;

    // Start the profile for this compute
    if (this->m_prof)
        this->m_prof->push("SinglePhaseFlowForces");

    // Grab handles for particle data
    // Access mode overwrite implies that data does not need to be read in
    ArrayHandle<Scalar4> d_force(this->m_force,access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_ratedpe(this->m_ratedpe,access_location::device, access_mode::readwrite);

    // Check input data, can be omitted if need be
    assert(d_force.data);
    assert(d_ratedpe.data);

    // Zero data before force calculation
    cudaMemset((void*)d_force.data,0,sizeof(Scalar4)*this->m_force.getNumElements());
    cudaMemset((void*)d_ratedpe.data,0,sizeof(Scalar4)*this->m_ratedpe.getNumElements());

    // access the particle data
    ArrayHandle<Scalar4> d_pos(this->m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_velocity(this->m_pdata->getVelocities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_dpe(this->m_pdata->getDPEs(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3> d_vf(this->m_pdata->getAuxiliaries1(), access_location::device,access_mode::read);
    ArrayHandle<Scalar>  d_h(this->m_pdata->getSlengths(), access_location::device, access_mode::read);

    // access the neighbor list
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(this->m_nlist->getNListArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_head_list(this->m_nlist->getHeadList(), access_location::device, access_mode::read);

    // Check input data
    assert(d_pos.data != NULL);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // for each fluid particle
    ArrayHandle< unsigned int > d_index_array(this->m_fluidgroup->getIndexArray(), access_location::device, access_mode::read);
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    ArrayHandle<unsigned int> d_type_property_map(this->m_type_property_map, access_location::device, access_mode::read);

    gpu_spf_forcecomputation(d_pos.data,
                                d_velocity.data,
                                d_dpe.data,
                                d_vf.data,
                                d_h.data,
                                d_force.data,
                                d_ratedpe.data,
                                d_n_neigh.data,
                                d_nlist.data,
                                d_head_list.data,
                                d_index_array.data,
                                group_size,
                                box,
                                this->m_rcutsq,
                                this->m_ch,
                                this->m_mu,
                                this->m_const_slength,
                                this->m_artificial_viscosity,
                                this->m_avalpha,
                                this->m_c,
                                this->m_avbeta,
                                this->m_ddiff,
                                this->m_skernel.get(),
                                d_type_property_map.data,
                                this->m_density_method == DENSITYSUMMATION,
                                this->m_density_method == DENSITYCONTINUITY,
                                this->m_compute_solid_forces,
                                this->m_density_diffusion
                            );

    // Add volumetric force (gravity)
    this->applyBodyForceGPU(timestep, this->m_fluidgroup);
    if ( this->m_compute_solid_forces )
        this->applyBodyForceGPU(timestep, this->m_solidgroup);

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Compute forces definition
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
void SinglePhaseFlowGPU<KT_, SET_>::computeForces(unsigned int timestep)
    {
    if (this->m_prof)
        this->m_prof->push("SinglePhaseFlow");

    // This is executed once to initialize protected/private variables
    if (!this->m_params_set)
        {
        this->m_exec_conf->msg->error() << "sph.models.SinglePhaseFlow requires parameters to be set before run()"
            << std::endl;
        throw std::runtime_error("Error computing SinglePhaseFlow forces");
        }

    // Make sure neighbor list is up-to-date
    this->m_nlist->compute(timestep);

    // Apply density renormalization if requested
    if (this-> m_shepard_renormalization && timestep % this->m_shepardfreq == 0 )
        {
        this->renormalize_density(timestep);
#ifdef ENABLE_MPI
         // Update ghost particle densities and pressures.
        this->update_ghost_dpe(timestep);
#endif
        }

    // Compute solid renormalization constant and, provided that SUMMATION approach is used,
    // particle mass densities.
    this->compute_ndensity(timestep);

    // Compute fluid pressure based on m_eos;
    this->compute_pressure(timestep);

#ifdef ENABLE_MPI
    // Update ghost particle densities and pressures.
    this->update_ghost_dpe(timestep);
#endif

    // Compute particle pressures
    this->compute_noslip(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    this->update_ghost_aux1(timestep);
#endif

    // Execute the force computation
    this->forcecomputation(timestep);

    if (this->m_prof)
        this->m_prof->pop();
    }

template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlowGPU_templ()
    {
    std::string name="SinglePFGPU_" + get_SE_name<SET_>() + "_" + get_kernel_name<KT_>();
    class_<SinglePhaseFlowGPU<KT_, SET_>, boost::shared_ptr<SinglePhaseFlowGPU<KT_, SET_> >, bases<ForceCompute>, boost::noncopyable>
        (name.c_str(), init< boost::shared_ptr<SystemDefinition>,
                                  boost::shared_ptr<SmoothingKernel<KT_> >,
                                  boost::shared_ptr<StateEquation<SET_> >,
                                  boost::shared_ptr<nsearch::NeighborList>,
                                  boost::shared_ptr<ParticleGroup>,
                                  boost::shared_ptr<ParticleGroup>,
                                  DensityMethod,
                                  ViscosityMethod >())
        .def("setParams", &SinglePhaseFlowGPU<KT_, SET_>::setParams)
        .def("getDensityMethod", &SinglePhaseFlowGPU<KT_, SET_>::getDensityMethod)
        .def("setDensityMethod", &SinglePhaseFlowGPU<KT_, SET_>::setDensityMethod)
        .def("getViscosityMethod", &SinglePhaseFlowGPU<KT_, SET_>::getViscosityMethod)
        .def("setViscosityMethod", &SinglePhaseFlowGPU<KT_, SET_>::setViscosityMethod)
        .def("setConstSmoothingLength", &SinglePhaseFlowGPU<KT_, SET_>::setConstSmoothingLength)
        .def("computeSolidForces", &SinglePhaseFlowGPU<KT_, SET_>::computeSolidForces)
        .def("activateArtificialViscosity", &SinglePhaseFlowGPU<KT_, SET_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity", &SinglePhaseFlowGPU<KT_, SET_>::deactivateArtificialViscosity)
        .def("activateDensityDiffusion", &SinglePhaseFlowGPU<KT_, SET_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion", &SinglePhaseFlowGPU<KT_, SET_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &SinglePhaseFlowGPU<KT_, SET_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization", &SinglePhaseFlowGPU<KT_, SET_>::deactivateShepardRenormalization)
        .def("setAcceleration", &SinglePhaseFlowGPU<KT_, SET_>::setAcceleration)
        ;

    }


void export_SinglePhaseFlowGPU()
{
    #define EXPORT(r,seq) export_SinglePhaseFlowGPU_templ<BOOST_PP_SEQ_ENUM(seq)>();
    BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPORT, (KERNELTYPES)(SEQTYPES))
    #undef EXPORT
}
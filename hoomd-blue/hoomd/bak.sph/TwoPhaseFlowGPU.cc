// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "TwoPhaseFlowGPU.h"

#include <boost/python.hpp>
using namespace boost::python;

#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
using namespace boost;

using namespace std;


/*! Constructor
*/
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
TwoPhaseFlowGPU<KT_, SET1_, SET2_>::TwoPhaseFlowGPU(boost::shared_ptr<SystemDefinition> sysdef,
                                            boost::shared_ptr<SmoothingKernel<KT_> > skernel,
                                            boost::shared_ptr<StateEquation<SET1_> > equationofstate1,
                                            boost::shared_ptr<StateEquation<SET2_> > equationofstate2,
                                            boost::shared_ptr<nsearch::NeighborList> nlist,
                                            boost::shared_ptr<ParticleGroup> fluidgroup1,
                                            boost::shared_ptr<ParticleGroup> fluidgroup2,
                                            boost::shared_ptr<ParticleGroup> solidgroup,
                                            DensityMethod mdensitymethod,
                                            ViscosityMethod mviscositymethod)
    : TwoPhaseFlow<KT_, SET1_, SET2_>(sysdef,skernel,equationofstate1,equationofstate2,nlist, fluidgroup1,fluidgroup2, solidgroup, mdensitymethod, mviscositymethod)
      {
        this->m_exec_conf->msg->notice(5) << "Constructing TwoPhaseFlow GPU" << std::endl;
      }

/*! Destructor
*/
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
TwoPhaseFlowGPU<KT_, SET1_, SET2_>::~TwoPhaseFlowGPU()
    {
    this->m_exec_conf->msg->notice(5) << "Destroying TwoPhaseFlow GPU" << std::endl;
    }


/*! Perform number density computation
 * This method computes and stores
     - the number density based mass density ( rho_i = m_i * \sum w_ij ) for fluid particles
       if the SUMMATION approach is being used.
     - the zeroth order normalization constant for solid particles
   in the x-component of the d_dpe Array.
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::compute_ndensity(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Number Density" << std::endl;

    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowNDensity");

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
    gpu_tpf_compute_ndensity(d_pos.data,
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
                         this->m_rho01,
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
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::compute_pressure(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Pressure" << std::endl;

    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowPressure");

    // Define ArrayHandles
    ArrayHandle<Scalar3> d_dpe(this->m_pdata->getDPEs(), access_location::device, access_mode::readwrite);

    // For each fluid particle
    unsigned int group_size1 = this->m_fluidgroup1->getNumMembers();
    unsigned int group_size2 = this->m_fluidgroup2->getNumMembers();
    ArrayHandle< unsigned int > d_index_array1(this->m_fluidgroup1->getIndexArray(), access_location::device, access_mode::read);
    ArrayHandle< unsigned int > d_index_array2(this->m_fluidgroup2->getIndexArray(), access_location::device, access_mode::read);

    gpu_tpf_compute_pressure(d_dpe.data,
                             d_index_array1.data,
                             group_size1,
                             this->m_eos1.get());

    gpu_tpf_compute_pressure(d_dpe.data,
                             d_index_array2.data,
                             group_size2,
                             this->m_eos1.get());

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Compute fictitious solid particle properties
    This method updates fictitious solid particle pressures and velocities to account for
    no-slip boundary conditions. Method follows Adami et. al. (2012).
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::compute_noslip(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::NoSlip" << std::endl;

    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowNoSlip");

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

    gpu_tpf_compute_noslip(d_vf.data,
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
                       this->m_rho01,
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
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::renormalize_density(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Density renormalization" << std::endl;

    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowDensityRenormalization");

    // Grab handles for particle data
    ArrayHandle<Scalar3> d_dpe(this->m_pdata->getDPEs(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_pos(this->m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_velocity(this->m_pdata->getVelocities(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h(this->m_pdata->getSlengths(), access_location::device, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(this->m_nlist->getNListArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_head_list(this->m_nlist->getHeadList(), access_location::device, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Particle loop
    // For each fluid particle
    ArrayHandle< unsigned int > d_index_array(this->m_fluidgroup->getIndexArray(), access_location::device, access_mode::read);
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    ArrayHandle<unsigned int> d_type_property_map(this->m_type_property_map, access_location::device, access_mode::read);

    gpu_tpf_compute_ndensityrenormalization(d_pos.data,
                                        d_dpe.data,
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
                                        this->m_skernel.get(),
                                        d_type_property_map.data
    );

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Compute interfacial color gradients
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::compute_colorgradients(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::Normals" << std::endl;

    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowNormals");

    // Grab handles for particle and neighbor data
    ArrayHandle<Scalar3> d_sn(this->m_pdata->getAuxiliaries2(), access_location::device,access_mode::readwrite);
    ArrayHandle<Scalar3> d_fn(this->m_pdata->getAuxiliaries3(), access_location::device,access_mode::readwrite);
    ArrayHandle<Scalar4> d_pos(this->m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_dpe(this->m_pdata->getDPEs(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h(this->m_pdata->getSlengths(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_velocity(this->m_pdata->getVelocities(), access_location::device, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(this->m_nlist->getNListArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_head_list(this->m_nlist->getHeadList(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_type_property_map(this->m_type_property_map, access_location::device, access_mode::read);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Zero data before calculation
    cudaMemset((void*)d_sn.data,0,sizeof(Scalar3)*this->m_pdata->getAuxiliaries2().getNumElements());
    cudaMemset((void*)d_fn.data,0,sizeof(Scalar3)*this->m_pdata->getAuxiliaries3().getNumElements());

    gpu_tpf_colorgradients(d_pos.data,
                           d_velocity.data,
                           d_dpe.data,
                           d_h.data,
                           d_sn.data,
                           d_fn.data,
                           d_n_neigh.data,
                           d_nlist.data,
                           d_head_list.data,
                           this->m_pdata->getN(),
                           box,
                           this->m_rcutsq,
                           this->m_const_slength,
                           this->m_artificial_viscosity,
                           this->m_avalpha,
                           this->m_ch,
                           this->m_avbeta,
                           this->m_ddiff,
                           this->m_skernel.get(),
                           d_type_property_map.data,
                           this->m_compute_solid_forces,
                           this->m_density_diffusion);

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Compute surface force vectors
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::compute_surfaceforce(unsigned int timestep)
    {
    this->m_exec_conf->msg->notice(7) << "Computing TwoPhaseFlow::SurfaceForce" << std::endl;

    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlowSurfaceForce");

    // Grab handles for particle and neighbor data
    ArrayHandle<Scalar3> d_sf(this->m_pdata->getAuxiliaries4(), access_location::device,access_mode::readwrite);
    ArrayHandle<Scalar3> d_sn(this->m_pdata->getAuxiliaries2(), access_location::device,access_mode::read);
    ArrayHandle<Scalar3> d_fn(this->m_pdata->getAuxiliaries3(), access_location::device,access_mode::read);
    ArrayHandle<Scalar4> d_pos(this->m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<Scalar3> d_dpe(this->m_pdata->getDPEs(), access_location::device, access_mode::read);
    ArrayHandle<Scalar>  d_h(this->m_pdata->getSlengths(), access_location::device, access_mode::read);
    ArrayHandle<Scalar4> d_velocity(this->m_pdata->getVelocities(), access_location::device, access_mode::read);

    // Grab handles for neighbor data
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(this->m_nlist->getNListArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_head_list(this->m_nlist->getHeadList(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_type_property_map(this->m_type_property_map, access_location::device, access_mode::read);

    ArrayHandle< unsigned int > d_index_array(this->m_fluidgroup->getIndexArray(), access_location::device, access_mode::read);
    unsigned int group_size = this->m_fluidgroup->getNumMembers();
    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // Zero data before calculation
    cudaMemset((void*)d_sf.data,0,sizeof(Scalar3)*this->m_pdata->getAuxiliaries4().getNumElements());

    gpu_tpf_compute_surfaceforce(d_pos.data,
                                 d_velocity.data,
                                 d_dpe.data,
                                 d_sf.data,
                                 d_h.data,
                                 d_sn.data,
                                 d_fn.data,
                                 d_n_neigh.data,
                                 d_nlist.data,
                                 d_head_list.data,
                                 d_index_array.data,
                                 group_size,
                                 box,
                                 this->m_rcutsq,
                                 this->m_const_slength,
                                 this->m_sigma01,
                                 this->m_sigma02,
                                 this->m_sigma12,
                                 this->m_ch,
                                 this->m_skernel.get(),
                                 d_type_property_map.data,
                                 this->m_compute_solid_forces,
                                 this->m_density_diffusion);

    if (this->m_prof)
        this->m_prof->pop();
    }


/*! Perform force computation
 */
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::forcecomputation(unsigned int timestep)
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
    ArrayHandle<Scalar3> d_sf(this->m_pdata->getAuxiliaries4(), access_location::device,access_mode::read);


    // access the neighbor list
    ArrayHandle<unsigned int> d_n_neigh(this->m_nlist->getNNeighArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_nlist(this->m_nlist->getNListArray(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_head_list(this->m_nlist->getHeadList(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_type_property_map(this->m_type_property_map, access_location::device, access_mode::read);

    // Check input data
    assert(d_pos.data != NULL);

    // Local copy of the simulation box
    const BoxDim& box = this->m_pdata->getGlobalBox();

    // for each fluid particle
    ArrayHandle< unsigned int > d_index_array(this->m_fluidgroup->getIndexArray(), access_location::device, access_mode::read);
    unsigned int group_size = this->m_fluidgroup->getNumMembers();

    gpu_tpf_forcecomputation(d_pos.data,
                                d_velocity.data,
                                d_dpe.data,
                                d_vf.data,
                                d_h.data,
                                d_force.data,
                                d_ratedpe.data,
                                d_sf.data,
                                d_n_neigh.data,
                                d_nlist.data,
                                d_head_list.data,
                                d_index_array.data,
                                group_size,
                                box,
                                this->m_rcutsq,
                                this->m_ch,
                                this->m_mu1,
                                this->m_mu2,
                                this->m_rho01,
                                this->m_rho02,
                                this->m_c1,
                                this->m_c2,
                                this->m_const_slength,
                                this->m_artificial_viscosity,
                                this->m_avalpha,
                                this->m_cmax,
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
    this->applyBodyForce(timestep, this->m_fluidgroup);

    if (this->m_prof)
        this->m_prof->pop();
    }

/*! Compute forces definition
*/
template<SmoothingKernelType KT_,StateEquationType SET1_,StateEquationType SET2_>
void TwoPhaseFlowGPU<KT_, SET1_, SET2_>::computeForces(unsigned int timestep)
    {
    if (this->m_prof)
        this->m_prof->push("TwoPhaseFlow");

    // This is executed once to initialize protected/private variables
    if (!this->m_params_set)
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
        this->renormalize_density(timestep);
#ifdef ENABLE_MPI
         // Update ghost particle densities and pressures.
        this->update_ghost_dpe(timestep);
#endif
        }

    // Compute solid renormalization constant and, provided that SUMMATION approach is used,
    // particle mass densities.
    this->compute_ndensity(timestep);

    // Compute fluid pressure based on this->m_eos;
    this->compute_pressure(timestep);

#ifdef ENABLE_MPI
    // Update ghost particle densities and pressures.
    this->update_ghost_dpe(timestep);
#endif

    // Compute particle pressures
    this->compute_noslip(timestep);

#ifdef ENABLE_MPI
    // Update ghost particle densities and pressures.
    this->update_ghost_dpe(timestep);
#endif

    // Compute particle interfacial color gradient
    this->compute_colorgradients(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    this->update_ghost_aux123(timestep);
#endif

    this->compute_surfaceforce(timestep);

#ifdef ENABLE_MPI
    // Update ghost particles
    this->update_ghost_aux4(timestep);
#endif

    // Execute the force computation
    this->forcecomputation(timestep);

    if (this->m_prof)
        this->m_prof->pop();
    }

template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
void export_TwoPhaseFlowGPU_templ()
    {
    std::string name="TwoPFGPU_" + get_SE_name<SET1_>()+ get_SE_name<SET2_>() + "_" + get_kernel_name<KT_>();
    class_<TwoPhaseFlowGPU<KT_,SET1_,SET2_>, boost::shared_ptr<TwoPhaseFlowGPU<KT_,SET1_,SET2_> >, bases<ForceCompute>, boost::noncopyable>
        (name.c_str(), init< boost::shared_ptr<SystemDefinition>,
                               boost::shared_ptr<SmoothingKernel<KT_> >,
                               boost::shared_ptr<StateEquation<SET1_> >,
                               boost::shared_ptr<StateEquation<SET2_> >,
                               boost::shared_ptr<nsearch::NeighborList>,
                               boost::shared_ptr<ParticleGroup>,
                               boost::shared_ptr<ParticleGroup>,
                               boost::shared_ptr<ParticleGroup>,
                               DensityMethod,
                               ViscosityMethod >())
        .def("setParams", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::setParams)
        .def("getDensityMethod", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::getDensityMethod)
        .def("setDensityMethod", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::setDensityMethod)
        .def("getViscosityMethod", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::getViscosityMethod)
        .def("setViscosityMethod", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::setViscosityMethod)
        .def("setConstSmoothingLength", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::setConstSmoothingLength)
        .def("computeSolidForces", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::computeSolidForces)
        .def("activateArtificialViscosity", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::activateArtificialViscosity)
        .def("deactivateArtificialViscosity", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::deactivateArtificialViscosity)
        .def("activateDensityDiffusion", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::activateDensityDiffusion)
        .def("deactivateDensityDiffusion", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::deactivateDensityDiffusion)
        .def("activateShepardRenormalization", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::activateShepardRenormalization)
        .def("deactivateShepardRenormalization", &TwoPhaseFlowGPU<KT_,SET1_,SET2_>::deactivateShepardRenormalization)
        .def("setAcceleration", &SPHBaseClass<KT_,SET2_>::setAcceleration)
        ;

    }
    
void export_TwoPhaseFlowGPU()
{
    #define EXPORT(r,seq) export_TwoPhaseFlowGPU_templ<BOOST_PP_SEQ_ENUM(seq)>();
    BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPORT, (KERNELTYPES)(SEQTYPES)(SEQTYPES))
    #undef EXPORT
}

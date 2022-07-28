// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "hoomd/BoxDim.h"
#include <vector>

#include "SmoothingKernel.h"
#include "StateEquations.h"

//! Helper function to compute particle pressures
template<StateEquationType SET_>
cudaError_t gpu_spf_compute_pressure(
                                            Scalar3 *d_dpe,
                                            unsigned int *d_index_array,
                                            unsigned int group_size,
                                            StateEquation<SET_> *m_eos);

/*! Helper function to compute particle number density
 * \post Particle number densities are stores in charge array
 */
template<SmoothingKernelType KT_>
cudaError_t gpu_spf_compute_ndensity(Scalar4 *d_pos,
                                 Scalar3 *d_dpe,
                                 Scalar4 *d_vel,
                                 Scalar *d_h,
                                 unsigned int *d_n_neigh,
                                 unsigned int *d_nlist,
                                 unsigned int *d_head_list,
                                 unsigned int N,
                                 const BoxDim& box,
                                 Scalar m_ch,
                                 Scalar m_rcutsq,
                                 Scalar m_rho0,
                                 bool m_const_slength,
                                 SmoothingKernel<KT_> *m_skernel,
                                 unsigned int * d_type_property_map,
                                 bool density_sum
                                    );


/*! Helper function to compute renormalized number density
 * \post Fluid particle number density field is renormalized
 */
template<SmoothingKernelType KT_>
cudaError_t gpu_spf_compute_ndensityrenormalization(Scalar4 *d_pos,
                                                Scalar3 *d_dpe,
                                                Scalar3 *d_dpe_old,
                                                Scalar4 *d_vel,
                                                Scalar *d_h,
                                                unsigned int *d_n_neigh,
                                                unsigned int *d_nlist,
                                                unsigned int *d_head_list,
                                                unsigned int *d_index_array,
                                                unsigned int group_size,
                                                const BoxDim& box,
                                                Scalar m_rcutsq,
                                                Scalar m_ch,
                                                bool m_const_slength,
                                                SmoothingKernel<KT_> *m_skernel);

/*! Helper function to compute fictitious solid particle properties (pressures and velocities)
 * \pre Ghost particle number densities (i.e. charge array) must be up-to-date
 * \pre Normalization constant \sum_fluids Wij must be computed and stored in orientation array
 * \post Fictitious particle properties are computed and stored in orentiation array
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
cudaError_t gpu_spf_compute_noslip(Scalar3 *d_vf,
                               Scalar4 *d_pos,
                               Scalar3 *d_dpe,
                               Scalar4 *d_velocity,
                               Scalar3 *d_accel,
                               Scalar *d_h,
                               unsigned int *d_n_neigh,
                               unsigned int *d_nlist,
                               unsigned int *d_head_list,
                               unsigned int *d_index_array,
                               unsigned int group_size,
                               const BoxDim& box,
                               Scalar m_rcutsq,
                               Scalar m_ch,
                               Scalar m_rho0,
                               bool m_const_slength,
                               Scalar3 g_accel,
                               StateEquation<SET_> *m_eos,
                               SmoothingKernel<KT_> *m_skernel,
                               unsigned int *d_type_property_map,
                               bool m_body_acceleration
                                  );

/*! Helper function to set communication flags and update ghosts charges
 * \param timestep The time step
 * \post Ghost particle charges are up-to-date
 */
cudaError_t gpu_spf_compute_properties(Scalar4 *d_vel,
                                   unsigned int *d_group_members,
                                   unsigned int group_size,
                                   Scalar3 *sum);

//! Kernel driver for the first part of the NVE update called by RigidBodyIntegratorGPU
template<SmoothingKernelType KT_>
cudaError_t gpu_spf_forcecomputation(Scalar4 *d_pos,
                                     Scalar4 *d_velocity,
                                     Scalar3 *d_dpe,
                                     Scalar3 *d_vf,
                                     Scalar *d_h,
                                        Scalar4 *d_force,
                                        Scalar4 *d_ratedpe,   
                                        unsigned int *d_n_neigh,
                                        unsigned int *d_nlist,
                                        unsigned int *d_head_list,
                                        unsigned int *d_group_members,
                                        unsigned int group_size,
                                        const BoxDim& box,
                                        Scalar m_rcutsq,
                                        Scalar m_h,
                                        Scalar m_mu,
                                        bool m_const_slength,
					Scalar m_artificial_viscosity,
				   	Scalar m_avalpha,
				 	Scalar m_c,
					Scalar m_avbeta,
					Scalar m_ddiff,
                                        SmoothingKernel<KT_> *m_skernel,
                                        unsigned int *d_type_property_map,
                                        bool density_sum,
                                        bool density_cont,
					bool m_compute_solid_forces,
					bool m_density_diffusion
                                        );

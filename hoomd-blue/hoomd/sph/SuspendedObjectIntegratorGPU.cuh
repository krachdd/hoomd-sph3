// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

//#ifndef __SUSPENDED_OBJECT_INTEGRATOR_H__
//#define __SUSPENDED_OBJECT_INTEGRATOR_H__

#include <cuda_runtime.h>

#include "hoomd/ParticleData.cuh"
#include "hoomd/HOOMDMath.h"

/*! \file SuspendedObjectIntegratorGPU.cuh
 *    \brief Declares GPU kernel code 
 */

//! Kernel driver for the the first step of the computation
cudaError_t gpu_SuspendedOIntegrator_step_one(Scalar4 *d_pos,
                                              Scalar4 *d_vel,
                                              unsigned int *d_group_members,
                                              unsigned int group_size,
                                              const BoxDim& box,
                                              Scalar3 centerofmass,
                                              vec3<Scalar> angularvel,
                                              Scalar3 translationvel, Scalar deltaT);

//! Kernel driver for the the second step of the computation called by NPTUpdaterGPU
cudaError_t gpu_SuspendedOIntegrator_step_two_a(Scalar4 *d_pos,
                                                Scalar4 *d_net_force,
                                                unsigned int *d_group_members,
                                                unsigned int group_size,
                                                const BoxDim& box,
                                                Scalar3 centerofmass,
                                                Scalar* totalforcetorque);

//! Kernel driver for the the second step of the computation called by NPTUpdaterGPU
cudaError_t gpu_SuspendedOIntegrator_step_two_b(Scalar4 *d_vel,
                                                Scalar4 *d_pos,
                                                unsigned int *d_group_members,
                                                unsigned int group_size,
                                                const BoxDim& box,
                                                Scalar3 centerofmass,
                                                Scalar3 translationvel,
                                                vec3<Scalar> angularvel);

//! Kernel driver
cudaError_t gpu_SuspendedOIntegrator_compute_center_of_mass(const Scalar4* pos,
                                                            const unsigned int* indices, unsigned int count,
                                                            Scalar* dst_angles,
                                                            Scalar3 lo, Scalar3 Ld);

//! Kernel driver
cudaError_t gpu_SuspendedOIntegrator_compute_translation_velocity(const Scalar4* vel,
                                                                  Scalar totalmass,
                                                                  const unsigned int* indices, unsigned int count,
                                                                  Scalar* translationvel);

//! Kernel driver
cudaError_t gpu_SuspendedOIntegrator_compute_moment_of_intertia(const Scalar4* pos, const Scalar4* vel,
                                                                const unsigned int* indices, unsigned int count,
                                                                Scalar3 centerofmass, BoxDim box,
                                                                Scalar* J);

//! Kernel driver
cudaError_t gpu_suspendedOIntegrator_compute_angular_velocity(Scalar4* pos, Scalar4* vel,
                                                              Scalar3 centerofmass, Scalar3 translationvel, BoxDim box,
                                                              const unsigned int* indices, unsigned int count,
                                                              Scalar* angmomentum);

//#endif

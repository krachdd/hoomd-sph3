// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

/*! \file RigidBodyIntegratorGPU.cuh
 \ brief Declares GPU kernel code for* NVE integration on the GPU. Used by RigidBodyIntegratorGPU.
 */

#include "hoomd/ParticleData.cuh"

//#ifndef __RIGID_BODY_INTEGRATOR_H__
//#define __RIGID_BODY_INTEGRATOR_H__

//! Kernel driver for the first part of the NVE update called by RigidBodyIntegratorGPU
cudaError_t gpu_RigidBody_step_one(Scalar4 *d_pos,
                                   Scalar4 *d_vel,
                                   Scalar3 *d_accel,
                                   unsigned int *d_group_members,
                                   unsigned int group_size,
                                   const BoxDim& box,
                                   Scalar deltaT,
                                   Scalar3 pivotpnt,
                                   Scalar3 rotaxis,
                                   Scalar transvel_x,
                                   Scalar transvel_y,
                                   Scalar transvel_z,
                                   Scalar rotatvel,
                                   Scalar transvel_x_dt,
                                   Scalar transvel_y_dt,
                                   Scalar transvel_z_dt,
                                   Scalar rotatvel_dt,
                                   unsigned int block_size);

//! Kernel driver for the second part of the NVE update called by RigidBodyIntegratorGPU
cudaError_t gpu_RigidBody_step_two(Scalar4 *d_pos,
                                   Scalar4 *d_vel,
                                   unsigned int *d_group_members,
                                   unsigned int group_size,
                                   Scalar deltaT,
                                   unsigned int block_size);

//#endif // __RIGID_BODY_INTEGRATOR_H__
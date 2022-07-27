// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

/*! \file VelocityVerletGPU.cuh
 \ brief Declares GPU kernel c*ode
 */

#include "hoomd/ParticleData.cuh"
#include "hoomd/HOOMDMath.h"

//#ifndef __VELOCITY_VERLET_H__
//#define __VELOCITY_VERLET_H__

//! Kernel driver for the first part of the NVE update called by TwoStepNVEGPU
cudaError_t gpu_VelocityVerlet_step_one(Scalar4 *d_pos,
                                        Scalar4 *d_vel,
                                        const Scalar3 *d_accel,
                                        int3 *d_image,
                                        Scalar3 *d_dpe,
                                        Scalar3 *d_dpedt,
                                        unsigned int *d_group_members,
                                        unsigned int group_size,
                                        const BoxDim& box,
                                        Scalar deltaT,
                                        unsigned int block_size);

//! Kernel driver for the second part of the NVE update called by TwoStepNVEGPU
cudaError_t gpu_VelocityVerlet_step_two(Scalar4 *d_pos,
                                        Scalar4 *d_vel,
                                        Scalar3 *d_accel,
                                        Scalar3 *d_dpe,
                                        Scalar3 *d_dpedt,
                                        unsigned int *d_group_members,
                                        unsigned int group_size,
                                        Scalar4 *d_net_force,
                                        Scalar4 *d_net_ratedpe,
                                        Scalar deltaT,
                                        unsigned int block_size
);

//#endif //__VELOCITY_VERLET_H__
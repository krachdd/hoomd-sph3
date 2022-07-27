// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "RigidBodyIntegratorGPU.cuh"

/*! \file RigidBodyIntegratorGPU.cu
 *    \brief Defines GPU kernel code for NVE integration on the GPU. Used by RigidBodyIntegratorGPU.
 */

//! Takes the first half-step forward in the velocity-verlet NVE integration on a group of particles
/*! \param d_pos array of particle positions
 *    \param d_vel array of particle velocities
 *    \param d_accel array of particle accelerations
 *    \param d_image array of particle images
 *    \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 *    \param group_size Number of members in the group
 *    \param box Box dimensions for periodic boundary condition handling
 *    \param deltaT timestep
 *    \param limit If \a limit is true, then the dynamics will be limited so that particles do not move
 *        a distance further than \a limit_val in one step.
 *    \param limit_val Length to limit particle distance movement to
 *    \param zero_force Set to true to always assign an acceleration of 0 to all particles in the group
 * 
 *    This kernel must be executed with a 1D grid of any block size such that the number of threads is greater than or
 *    equal to the number of members in the group. The kernel's implementation simply reads one particle in each thread
 *    and updates that particle.
 * 
 *    <b>Performance notes:</b>
 *    Particle properties are read via the texture cache to optimize the bandwidth obtained with sparse groups. The writes
 *    in sparse groups will not be coalesced. However, because ParticleGroup sorts the index list the writes will be as
 *    contiguous as possible leading to fewer memory transactions on compute 1.3 hardware and more cache hits on Fermi.
 */
extern "C" __global__
void gpu_RigidBody_step_one_kernel(Scalar4 *d_pos,
                                   Scalar4 *d_vel,
                                   Scalar3 *d_accel,
                                   unsigned int *d_group_members,
                                   unsigned int group_size,
                                   const BoxDim box,
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
                                   Scalar rotatvel_dt)
{
    // determine which particle this thread works on (MEM TRANSFER: 4 bytes)
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (group_idx < group_size)
    {
        vec3<Scalar> angularvel(rotatvel*rotaxis.x, rotatvel*rotaxis.y, rotatvel*rotaxis.z);
        vec3<Scalar> angularaccel(rotatvel_dt*rotaxis.x, rotatvel_dt*rotaxis.y, rotatvel_dt*rotaxis.z);
        
        unsigned int idx = d_group_members[group_idx];
        
        Scalar3 pos = make_scalar3(d_pos[idx].x, d_pos[idx].y, d_pos[idx].z);
        
        // Relative position to rotational pivot point
        Scalar3 rdiff = pos-pivotpnt;
        
        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);
        
        // Relative position to center of mass
        vec3<Scalar> rdiffv(rdiff.x, rdiff.y, rdiff.z);
        
        // Angular velocity
        vec3<Scalar> rot = cross(angularvel,rdiffv);
        
        // Angular acceleration
        vec3<Scalar> rotaccel = cross(angularaccel,rdiffv);
        
        // Centrifugal acceleration
        vec3<Scalar> centrifugal = cross(angularvel,rot);
        
        // Update slave particle velocities
        d_vel[idx].x = transvel_x + rot.x;
        d_vel[idx].y = transvel_y + rot.y;
        d_vel[idx].z = transvel_z + rot.z;
        
        // Update acceleration of slave particles
        d_accel[idx].x = transvel_x_dt + rotaccel.x + centrifugal.x;
        d_accel[idx].y = transvel_y_dt + rotaccel.y + centrifugal.y;
        d_accel[idx].z = transvel_z_dt + rotaccel.z + centrifugal.z;
        
        // Update slave particle positions
        d_pos[idx].x += Scalar(1.0/2.0)*d_vel[idx].x*deltaT;
        d_pos[idx].y += Scalar(1.0/2.0)*d_vel[idx].y*deltaT;
        d_pos[idx].z += Scalar(1.0/2.0)*d_vel[idx].z*deltaT;
    }
}

/*! \param d_pos array of particle positions
 *    \param d_vel array of particle velocities
 *    \param d_accel array of particle accelerations
 *    \param d_image array of particle images
 *    \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 *    \param group_size Number of members in the group
 *    \param box Box dimensions for periodic boundary condition handling
 *    \param deltaT timestep
 *    \param limit If \a limit is true, then the dynamics will be limited so that particles do not move
 *        a distance further than \a limit_val in one step.
 *    \param limit_val Length to limit particle distance movement to
 *    \param zero_force Set to true to always assign an acceleration of 0 to all particles in the group
 * 
 *    See gpu_RigidBody_step_one_kernel() for full documentation, this function is just a driver.
 */
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
                                   unsigned int block_size)
{
    static unsigned int max_block_size = UINT_MAX;
    if (max_block_size == UINT_MAX)
    {
        cudaFuncAttributes attr;
        cudaFuncGetAttributes(&attr, (const void*)gpu_RigidBody_step_one_kernel);
        max_block_size = attr.maxThreadsPerBlock;
    }
    
    unsigned int run_block_size = min(block_size, max_block_size);
    
    // setup the grid to run the kernel
    dim3 grid( (group_size/run_block_size) + 1, 1, 1);
    dim3 threads(run_block_size, 1, 1);
    
    // run the kernel
    gpu_RigidBody_step_one_kernel<<< grid, threads >>>(d_pos,
                                                       d_vel, d_accel, d_group_members,
                                                       group_size, box, deltaT, pivotpnt, rotaxis, transvel_x, transvel_y, transvel_z, rotatvel,
                                                       transvel_x_dt, transvel_y_dt, transvel_z_dt, rotatvel_dt);
    
    return cudaSuccess;
}


//! Takes the second half-step forward in the velocity-verlet NVE integration on a group of particles
/*! \param d_vel array of particle velocities
 *    \param d_accel array of particle accelerations
 *    \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 *    \param group_size Number of members in the group
 *    \param d_net_force Net force on each particle
 *    \param deltaT Amount of real time to step forward in one time step
 *    \param limit If \a limit is true, then the dynamics will be limited so that particles do not move
 *        a distance further than \a limit_val in one step.
 *    \param limit_val Length to limit particle distance movement to
 *    \param zero_force Set to true to always assign an acceleration of 0 to all particles in the group
 * 
 *    This kernel is implemented in a very similar manner to gpu_RigidBody_step_one_kernel(), see it for design details.
 */
extern "C" __global__
void gpu_RigidBody_step_two_kernel(
    Scalar4 *d_pos,
    Scalar4 *d_vel,
    unsigned int *d_group_members,
    unsigned int group_size,
    Scalar deltaT)
{
    // determine which particle this thread works on (MEM TRANSFER: 4 bytes)
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (group_idx < group_size)
    {
        unsigned int idx = d_group_members[group_idx];
        
        d_pos[idx].x += Scalar(1.0/2.0)*d_vel[idx].x*deltaT;
        d_pos[idx].y += Scalar(1.0/2.0)*d_vel[idx].y*deltaT;
        d_pos[idx].z += Scalar(1.0/2.0)*d_vel[idx].z*deltaT;
    }
}

/*! \param d_vel array of particle velocities
 *    \param d_accel array of particle accelerations
 *    \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 *    \param group_size Number of members in the group
 *    \param d_net_force Net force on each particle
 *    \param deltaT Amount of real time to step forward in one time step
 *    \param limit If \a limit is true, then the dynamics will be limited so that particles do not move
 *        a distance further than \a limit_val in one step.
 *    \param limit_val Length to limit particle distance movement to
 *    \param zero_force Set to true to always assign an acceleration of 0 to all particles in the group
 * 
 *    This is just a driver for gpu_RigidBody_step_two_kernel(), see it for details.
 */
cudaError_t gpu_RigidBody_step_two(Scalar4 *d_pos,
                                   Scalar4 *d_vel,
                                   unsigned int *d_group_members,
                                   unsigned int group_size,
                                   Scalar deltaT,
                                   unsigned int block_size)
{
    static unsigned int max_block_size = UINT_MAX;
    if (max_block_size == UINT_MAX)
    {
        cudaFuncAttributes attr;
        cudaFuncGetAttributes(&attr, (const void *)gpu_RigidBody_step_two_kernel);
        max_block_size = attr.maxThreadsPerBlock;
    }
    
    unsigned int run_block_size = min(block_size, max_block_size);
    
    // setup the grid to run the kernel
    dim3 grid( (group_size/run_block_size) + 1, 1, 1);
    dim3 threads(run_block_size, 1, 1);
    
    // run the kernel
    gpu_RigidBody_step_two_kernel<<< grid, threads >>>(d_pos,
                                                       d_vel,
                                                       d_group_members,
                                                       group_size,
                                                       deltaT);
    
    return cudaSuccess;
}

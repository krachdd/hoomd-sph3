// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "VelocityVerletGPU.cuh"
#include "hoomd/VectorMath.h"


/*! \file VelocityVerletGPU.cu
 \ *brief Defines GPU kernel code for NVE integration on the GPU. Used by VelocityVerletGPU.
 */

//! Takes the first half-step forward in the velocity-verlet NVE integration on a group of particles
/*! \param d_pos array of particle positions
 \ *param d_vel array of particle velocities
 \param d_accel array of particle accelerations
 \param d_image array of particle images
 \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 \param group_size Number of members in the group
 \param box Box dimensions for periodic boundary condition handling
 \param deltaT timestep
 
 This kernel must be executed with a 1D grid of any block size such that the number of threads is greater than or
 equal to the number of members in the group. The kernel's implementation simply reads one particle in each thread
 and updates that particle.
 
 <b>Performance notes:</b>
 Particle properties are read via the texture cache to optimize the bandwidth obtained with sparse groups. The writes
 in sparse groups will not be coalesced. However, because ParticleGroup sorts the index list the writes will be as
 contiguous as possible leading to fewer memory transactions on compute 1.3 hardware and more cache hits on Fermi.
 */
extern "C" __global__
void gpu_VelocityVerlet_step_one_kernel(Scalar4 *d_pos,
                                        Scalar4 *d_vel,
                                        const Scalar3 *d_accel,
                                        int3 *d_image,
                                        Scalar3 *d_dpe,
                                        Scalar3 *d_dpedt,
                                        unsigned int *d_group_members,
                                        unsigned int group_size,
                                        BoxDim box,
                                        Scalar deltaT)
{
    // determine which particle this thread works on (MEM TRANSFER: 4 bytes)
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (group_idx < group_size)
    {
        unsigned int idx = d_group_members[group_idx];
        
        // dpe(t+deltaT/2) = dpe(t) + (1/2)*dpedt(t)*deltaT
        d_dpe[idx].x += Scalar(1.0/2.0)*d_dpedt[idx].x*deltaT;
        d_dpe[idx].y += Scalar(1.0/2.0)*d_dpedt[idx].y*deltaT;
        d_dpe[idx].z += Scalar(1.0/2.0)*d_dpedt[idx].z*deltaT;

        // v(t+deltaT/2) = v(t) + (1/2)*a(t)*deltaT
        d_vel[idx].x += Scalar(1.0/2.0)*d_accel[idx].x*deltaT;
        d_vel[idx].y += Scalar(1.0/2.0)*d_accel[idx].y*deltaT;
        d_vel[idx].z += Scalar(1.0/2.0)*d_accel[idx].z*deltaT;

        // r(t+deltaT/2) = r(t) + v(t+deltaT/2)*deltaT/2
        d_pos[idx].x += Scalar(1.0/2.0)*d_vel[idx].x*deltaT;
        d_pos[idx].y += Scalar(1.0/2.0)*d_vel[idx].y*deltaT;
        d_pos[idx].z += Scalar(1.0/2.0)*d_vel[idx].z*deltaT;
        
        // read in the particle's image (MEM TRANSFER: 16 bytes)
        int3 image = d_image[idx];
        
        // fix the periodic boundary conditions (FLOPS: 15)
        box.wrap(d_pos[idx], image);
        
        d_image[idx] = image;
    }
}

/*! \param d_pos array of particle positions
 \ *param d_vel array of particle velocities
 \param d_accel array of particle accelerations
 \param d_image array of particle images
 \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 \param group_size Number of members in the group
 \param box Box dimensions for periodic boundary condition handling
 \param deltaT timestep
 
 See gpu_VelocityVerlet_step_one_kernel() for full documentation, this function is just a driver.
 */
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
                                        unsigned int block_size)
{
    static unsigned int max_block_size = UINT_MAX;
    if (max_block_size == UINT_MAX)
    {
        cudaFuncAttributes attr;
        cudaFuncGetAttributes(&attr, (const void*)gpu_VelocityVerlet_step_one_kernel);
        max_block_size = attr.maxThreadsPerBlock;
    }
    
    unsigned int run_block_size = min(block_size, max_block_size);
    
    // setup the grid to run the kernel
    dim3 grid( (group_size/run_block_size) + 1, 1, 1);
    dim3 threads(run_block_size, 1, 1);
    
    // run the kernel
    gpu_VelocityVerlet_step_one_kernel<<< grid, threads >>>(d_pos, d_vel, d_accel, d_image, d_dpe, d_dpedt, d_group_members, group_size, box, deltaT);
    
    return cudaSuccess;
}

//! Takes the second half-step forward in the velocity-verlet NVE integration on a group of particles
/*! \param d_vel array of particle velocities
 \ *param d_accel array of particle accelerations
 \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 \param group_size Number of members in the group
 \param d_net_force Net force on each particle
 \param deltaT Amount of real time to step forward in one time step
 
 This kernel is implemented in a very similar manner to gpu_VelocityVerlet_step_one_kernel(), see it for design details.
 */
extern "C" __global__
void gpu_VelocityVerlet_step_two_kernel(
    Scalar4 *d_pos,
    Scalar4 *d_vel,
    Scalar3 *d_accel,
    Scalar3 *d_dpe,
    Scalar3 *d_dpedt,
    unsigned int *d_group_members,
    unsigned int group_size,
    Scalar4 *d_net_force,
    Scalar4 *d_net_ratedpe,
    Scalar deltaT
)
{
    // determine which particle this thread works on (MEM TRANSFER: 4 bytes)
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (group_idx < group_size)
    {
        unsigned int idx = d_group_members[group_idx];
        
        // first, calculate acceleration from the net force
        Scalar minv = Scalar(1.0) / d_vel[idx].w;
        d_accel[idx].x = d_net_force[idx].x*minv;
        d_accel[idx].y = d_net_force[idx].y*minv;
        d_accel[idx].z = d_net_force[idx].z*minv;

        d_dpedt[idx].x = d_net_ratedpe[idx].x;
        d_dpedt[idx].y = d_net_ratedpe[idx].y;
        d_dpedt[idx].z = d_net_ratedpe[idx].z;

        // dpe(t+deltaT) = dpe(t+deltaT/2) + 1/2 * dpedt(t+deltaT)*deltaT
        d_dpe[idx].x += Scalar(1.0/2.0)*d_dpedt[idx].x*deltaT;
        d_dpe[idx].y += Scalar(1.0/2.0)*d_dpedt[idx].y*deltaT;
        d_dpe[idx].z += Scalar(1.0/2.0)*d_dpedt[idx].z*deltaT;

        // r(t+deltaT) = r(t+deltaT/2) + v(t+deltaT/2)*deltaT/2
        d_pos[idx].x += Scalar(1.0/2.0)*d_vel[idx].x*deltaT;
        d_pos[idx].y += Scalar(1.0/2.0)*d_vel[idx].y*deltaT;
        d_pos[idx].z += Scalar(1.0/2.0)*d_vel[idx].z*deltaT;

        // v(t+deltaT) = v(t+deltaT/2) + 1/2 * a(t+deltaT)*deltaT
        d_vel[idx].x += Scalar(1.0/2.0)*d_accel[idx].x*deltaT;
        d_vel[idx].y += Scalar(1.0/2.0)*d_accel[idx].y*deltaT;
        d_vel[idx].z += Scalar(1.0/2.0)*d_accel[idx].z*deltaT;
        
    }
}

/*! \param d_vel array of particle velocities
 \ *param d_accel array of particle accelerations
 \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 \param group_size Number of members in the group
 \param d_net_force Net force on each particle
 \param deltaT Amount of real time to step forward in one time step
 
 This is just a driver for gpu_VelocityVerlet_step_two_kernel(), see it for details.
 */
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
)
{
    static unsigned int max_block_size = UINT_MAX;
    if (max_block_size == UINT_MAX)
    {
        cudaFuncAttributes attr;
        cudaFuncGetAttributes(&attr, (const void *)gpu_VelocityVerlet_step_two_kernel);
        max_block_size = attr.maxThreadsPerBlock;
    }
    
    unsigned int run_block_size = min(block_size, max_block_size);
    
    // setup the grid to run the kernel
    dim3 grid( (group_size/run_block_size) + 1, 1, 1);
    dim3 threads(run_block_size, 1, 1);
    
    // run the kernel
    gpu_VelocityVerlet_step_two_kernel<<< grid, threads >>>(d_pos,
                                                            d_vel,
                                                            d_accel,
                                                            d_dpe,
                                                            d_dpedt,
                                                            d_group_members,
                                                            group_size,
                                                            d_net_force,
                                                            d_net_ratedpe,
                                                            deltaT
    );
    
    return cudaSuccess;
}
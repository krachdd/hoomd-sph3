// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "SPHBaseClass.cuh"

//! Helper function to compute particle pressures
__global__ void gpu_sphbase_apply_body_force_kernel(
                                 Scalar4 *d_force,
                                 Scalar4 *d_velocity,
                                 Scalar3 bforce,
                                 unsigned int *d_index_array,
                                 unsigned int group_size)
{
    // determine which particle this thread works on
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (group_idx < group_size)
    {
        // Read particle index
        unsigned int i = d_index_array[group_idx];
        
        // Evaluate pressure
        Scalar mi = d_velocity[i].w;

        // Add contribution to force
        d_force[i].x += bforce.x*mi;
        d_force[i].y += bforce.y*mi;
        d_force[i].z += bforce.z*mi;
    }
}
cudaError_t gpu_sphbase_apply_body_force(
                                 Scalar4 *d_force,
                                 Scalar4 *d_velocity,
                                 Scalar3 bforce,
                                 unsigned int *d_index_array,
                                 unsigned int group_size)
{
    // setup the grid to run the kernel
    unsigned int block_size=256;
    dim3 grid( (group_size / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    gpu_sphbase_apply_body_force_kernel<<< grid, threads >>>(
        d_force,
        d_velocity,
        bforce,
        d_index_array,
        group_size);
    
    return cudaSuccess;
}

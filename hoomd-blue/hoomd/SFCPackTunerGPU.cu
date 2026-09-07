// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*! \file SFCPackTunerGPU.cu
    \brief Defines GPU kernel code for generating the space-filling curve sorted order on the GPU.
   Used by SFCPackTunerGPU.
*/

#include <hip/hip_runtime.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#pragma GCC diagnostic pop

#include "SFCPackTunerGPU.cuh"

namespace hoomd
    {
namespace kernel
    {
//! Kernel to bin particles
template<bool twod>
__global__ void gpu_sfc_bin_particles_kernel(unsigned int N,
                                             const Scalar4* d_pos,
                                             unsigned int* d_particle_bins,
                                             const unsigned int* d_traversal_order,
                                             unsigned int n_grid,
                                             unsigned int* d_sorted_order,
                                             const BoxDim box)
    {
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;

    if (idx >= N)
        return;

    // fetch particle position
    Scalar4 postype = d_pos[idx];
    Scalar3 p = make_scalar3(postype.x, postype.y, postype.z);

    Scalar3 f = box.makeFraction(p);
    int ib = (unsigned int)(f.x * n_grid) % n_grid;
    int jb = (unsigned int)(f.y * n_grid) % n_grid;
    int kb = (unsigned int)(f.z * n_grid) % n_grid;

    // if the particle is slightly outside, move back into grid
    if (ib < 0)
        ib = 0;
    if (ib >= n_grid)
        ib = n_grid - 1;

    if (jb < 0)
        jb = 0;
    if (jb >= n_grid)
        jb = n_grid - 1;

    if (kb < 0)
        kb = 0;
    if (kb >= n_grid)
        kb = n_grid - 1;

    // record its bin
    unsigned int bin;
    if (twod)
        {
        // do not use Hilbert curve in 2D
        bin = ib * n_grid + jb;
        d_particle_bins[idx] = bin;
        }
    else
        {
        bin = ib * (n_grid * n_grid) + jb * n_grid + kb;
        d_particle_bins[idx] = d_traversal_order[bin];
        }

    // store index of ptl
    d_sorted_order[idx] = idx;
    }

/*! \param N number of local particles
    \param d_pos Device array of positions
    \param d_particle_bins Device array of particle bins
    \param d_traversal_order Device array of Hilbert-curve bins
    \param n_grid Number of grid elements along one edge
    \param d_sorted_order Sorted order of particles
    \param box Box dimensions
    \param twod If true, bin particles in two dimensions
    */
void gpu_generate_sorted_order(unsigned int N,
                               const Scalar4* d_pos,
                               unsigned int* d_particle_bins,
                               unsigned int* d_traversal_order,
                               unsigned int n_grid,
                               unsigned int* d_sorted_order,
                               const BoxDim& box,
                               bool twod,
                               CachedAllocator& alloc)
    {
    // maybe need to autotune, but SFCPackTuner is called infrequently
    unsigned int block_size = 256;
    unsigned int n_blocks = N / block_size + 1;

    if (twod)
        hipLaunchKernelGGL(HIP_KERNEL_NAME(gpu_sfc_bin_particles_kernel<true>),
                           dim3(n_blocks),
                           dim3(block_size),
                           0,
                           0,
                           N,
                           d_pos,
                           d_particle_bins,
                           d_traversal_order,
                           n_grid,
                           d_sorted_order,
                           box);
    else
        hipLaunchKernelGGL(HIP_KERNEL_NAME(gpu_sfc_bin_particles_kernel<false>),
                           dim3(n_blocks),
                           dim3(block_size),
                           0,
                           0,
                           N,
                           d_pos,
                           d_particle_bins,
                           d_traversal_order,
                           n_grid,
                           d_sorted_order,
                           box);

    // Sort particles
    if (N)
        {
        // Clear any sticky CUDA error from the prior binning kernel before calling Thrust.
        hipGetLastError();
        thrust::device_ptr<unsigned int> particle_bins(d_particle_bins);
        thrust::device_ptr<unsigned int> sorted_order(d_sorted_order);
#ifdef __HIP_PLATFORM_HCC__
        thrust::sort_by_key(thrust::hip::par(alloc),
#else
        thrust::sort_by_key(thrust::cuda::par(alloc),
#endif
                            particle_bins,
                            particle_bins + N,
                            sorted_order);
        }
    }

//! Kernel to apply sorted order (SPH field set)
__global__ void gpu_apply_sorted_order_kernel(unsigned int N,
                                              unsigned int n_ghost,
                                              const unsigned int* d_sorted_order,
                                              const Scalar4* d_pos,
                                              Scalar4* d_pos_alt,
                                              const Scalar4* d_vel,
                                              Scalar4* d_vel_alt,
                                              const Scalar* d_density,
                                              Scalar* d_density_alt,
                                              const Scalar* d_pressure,
                                              Scalar* d_pressure_alt,
                                              const Scalar* d_energy,
                                              Scalar* d_energy_alt,
                                              const Scalar3* d_aux1,
                                              Scalar3* d_aux1_alt,
                                              const Scalar3* d_aux2,
                                              Scalar3* d_aux2_alt,
                                              const Scalar3* d_aux3,
                                              Scalar3* d_aux3_alt,
                                              const Scalar3* d_aux4,
                                              Scalar3* d_aux4_alt,
                                              const Scalar* d_slength,
                                              Scalar* d_slength_alt,
                                              const Scalar3* d_accel,
                                              Scalar3* d_accel_alt,
                                              const Scalar3* d_dpedt,
                                              Scalar3* d_dpedt_alt,
                                              const int3* d_image,
                                              int3* d_image_alt,
                                              const unsigned int* d_body,
                                              unsigned int* d_body_alt,
                                              const unsigned int* d_tag,
                                              unsigned int* d_tag_alt,
                                              const Scalar4* d_net_force,
                                              Scalar4* d_net_force_alt,
                                              const Scalar4* d_net_ratedpe,
                                              Scalar4* d_net_ratedpe_alt,
                                              unsigned int* d_rtag)
    {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= N + n_ghost)
        return;

    // apply sorted order only for local ptls
    unsigned int old_idx = (idx < N ? d_sorted_order[idx] : idx);

    // permute and copy over particle data
    d_pos_alt[idx]         = d_pos[old_idx];
    d_vel_alt[idx]         = d_vel[old_idx];
    d_density_alt[idx]     = d_density[old_idx];
    d_pressure_alt[idx]    = d_pressure[old_idx];
    d_energy_alt[idx]      = d_energy[old_idx];
    d_aux1_alt[idx]        = d_aux1[old_idx];
    d_aux2_alt[idx]        = d_aux2[old_idx];
    d_aux3_alt[idx]        = d_aux3[old_idx];
    d_aux4_alt[idx]        = d_aux4[old_idx];
    d_slength_alt[idx]     = d_slength[old_idx];
    d_accel_alt[idx]       = d_accel[old_idx];
    d_dpedt_alt[idx]       = d_dpedt[old_idx];
    d_image_alt[idx]       = d_image[old_idx];
    d_body_alt[idx]        = d_body[old_idx];
    unsigned int tag       = d_tag[old_idx];
    d_tag_alt[idx]         = tag;
    d_net_force_alt[idx]   = d_net_force[old_idx];
    d_net_ratedpe_alt[idx] = d_net_ratedpe[old_idx];

    if (idx < N)
        {
        // update rtag to point to particle position in new arrays
        d_rtag[tag] = idx;
        }
    }

void gpu_apply_sorted_order(unsigned int N,
                            unsigned int n_ghost,
                            const unsigned int* d_sorted_order,
                            const Scalar4* d_pos,
                            Scalar4* d_pos_alt,
                            const Scalar4* d_vel,
                            Scalar4* d_vel_alt,
                            const Scalar* d_density,
                            Scalar* d_density_alt,
                            const Scalar* d_pressure,
                            Scalar* d_pressure_alt,
                            const Scalar* d_energy,
                            Scalar* d_energy_alt,
                            const Scalar3* d_aux1,
                            Scalar3* d_aux1_alt,
                            const Scalar3* d_aux2,
                            Scalar3* d_aux2_alt,
                            const Scalar3* d_aux3,
                            Scalar3* d_aux3_alt,
                            const Scalar3* d_aux4,
                            Scalar3* d_aux4_alt,
                            const Scalar* d_slength,
                            Scalar* d_slength_alt,
                            const Scalar3* d_accel,
                            Scalar3* d_accel_alt,
                            const Scalar3* d_dpedt,
                            Scalar3* d_dpedt_alt,
                            const int3* d_image,
                            int3* d_image_alt,
                            const unsigned int* d_body,
                            unsigned int* d_body_alt,
                            const unsigned int* d_tag,
                            unsigned int* d_tag_alt,
                            const Scalar4* d_net_force,
                            Scalar4* d_net_force_alt,
                            const Scalar4* d_net_ratedpe,
                            Scalar4* d_net_ratedpe_alt,
                            unsigned int* d_rtag)
    {
    unsigned int block_size = 256;
    unsigned int n_blocks = (N + n_ghost) / block_size + 1;

    hipLaunchKernelGGL(gpu_apply_sorted_order_kernel,
                       dim3(n_blocks),
                       dim3(block_size),
                       0,
                       0,
                       N,
                       n_ghost,
                       d_sorted_order,
                       d_pos,
                       d_pos_alt,
                       d_vel,
                       d_vel_alt,
                       d_density,
                       d_density_alt,
                       d_pressure,
                       d_pressure_alt,
                       d_energy,
                       d_energy_alt,
                       d_aux1,
                       d_aux1_alt,
                       d_aux2,
                       d_aux2_alt,
                       d_aux3,
                       d_aux3_alt,
                       d_aux4,
                       d_aux4_alt,
                       d_slength,
                       d_slength_alt,
                       d_accel,
                       d_accel_alt,
                       d_dpedt,
                       d_dpedt_alt,
                       d_image,
                       d_image_alt,
                       d_body,
                       d_body_alt,
                       d_tag,
                       d_tag_alt,
                       d_net_force,
                       d_net_force_alt,
                       d_net_ratedpe,
                       d_net_ratedpe_alt,
                       d_rtag);
    }

    } // end namespace kernel

    } // end namespace hoomd

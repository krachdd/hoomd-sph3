// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "ParticleData.cuh"

/*! \file ParticleData.cu
    \brief Implements GPU kernel code and data structure functions used by ParticleData
*/

#ifdef ENABLE_MPI

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <hipcub/hipcub.hpp>

#include <thrust/device_ptr.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/scatter.h>
#pragma GCC diagnostic pop

namespace hoomd
    {
namespace kernel
    {
//! Kernel to partition particle data (SPH version)
__global__ void gpu_scatter_particle_data_kernel(const unsigned int nwork,
                                                 const Scalar4* d_pos,
                                                 const Scalar4* d_vel,
                                                 const Scalar* d_density,
                                                 const Scalar* d_pressure,
                                                 const Scalar* d_energy,
                                                 const Scalar3* d_aux1,
                                                 const Scalar3* d_aux2,
                                                 const Scalar3* d_aux3,
                                                 const Scalar3* d_aux4,
                                                 const Scalar* d_slength,
                                                 const Scalar3* d_accel,
                                                 const Scalar3* d_dpedt,
                                                 const int3* d_image,
                                                 const unsigned int* d_body,
                                                 const Scalar4* d_net_force,
                                                 const Scalar4* d_net_ratedpe,
                                                 const unsigned int* d_tag,
                                                 unsigned int* d_rtag,
                                                 Scalar4* d_pos_alt,
                                                 Scalar4* d_vel_alt,
                                                 Scalar* d_density_alt,
                                                 Scalar* d_pressure_alt,
                                                 Scalar* d_energy_alt,
                                                 Scalar3* d_aux1_alt,
                                                 Scalar3* d_aux2_alt,
                                                 Scalar3* d_aux3_alt,
                                                 Scalar3* d_aux4_alt,
                                                 Scalar* d_slength_alt,
                                                 Scalar3* d_accel_alt,
                                                 Scalar3* d_dpedt_alt,
                                                 int3* d_image_alt,
                                                 unsigned int* d_body_alt,
                                                 Scalar4* d_net_force_alt,
                                                 Scalar4* d_net_ratedpe_alt,
                                                 unsigned int* d_tag_alt,
                                                 detail::pdata_element* d_out,
                                                 unsigned int* d_comm_flags,
                                                 unsigned int* d_comm_flags_out,
                                                 const unsigned int* d_scan)
    {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= nwork)
        return;
    bool remove = d_comm_flags[idx];

    unsigned int scan_remove = d_scan[idx];
    unsigned int scan_keep = idx - scan_remove;

    if (remove)
        {
        detail::pdata_element p;
        p.pos         = d_pos[idx];
        p.vel         = d_vel[idx];
        p.density     = d_density[idx];
        p.pressure    = d_pressure[idx];
        p.energy      = d_energy[idx];
        p.aux1        = d_aux1[idx];
        p.aux2        = d_aux2[idx];
        p.aux3        = d_aux3[idx];
        p.aux4        = d_aux4[idx];
        p.slength     = d_slength[idx];
        p.accel       = d_accel[idx];
        p.dpedt       = d_dpedt[idx];
        p.image       = d_image[idx];
        p.body        = d_body[idx];
        p.net_force   = d_net_force[idx];
        p.net_ratedpe = d_net_ratedpe[idx];
        p.tag         = d_tag[idx];
        d_out[scan_remove] = p;
        d_comm_flags_out[scan_remove] = d_comm_flags[idx];

        // reset communication flags
        d_comm_flags[idx] = 0;

        // reset rtag
        d_rtag[p.tag] = NOT_LOCAL;
        }
    else
        {
        d_pos_alt[scan_keep]         = d_pos[idx];
        d_vel_alt[scan_keep]         = d_vel[idx];
        d_density_alt[scan_keep]     = d_density[idx];
        d_pressure_alt[scan_keep]    = d_pressure[idx];
        d_energy_alt[scan_keep]      = d_energy[idx];
        d_aux1_alt[scan_keep]        = d_aux1[idx];
        d_aux2_alt[scan_keep]        = d_aux2[idx];
        d_aux3_alt[scan_keep]        = d_aux3[idx];
        d_aux4_alt[scan_keep]        = d_aux4[idx];
        d_slength_alt[scan_keep]     = d_slength[idx];
        d_accel_alt[scan_keep]       = d_accel[idx];
        d_dpedt_alt[scan_keep]       = d_dpedt[idx];
        d_image_alt[scan_keep]       = d_image[idx];
        d_body_alt[scan_keep]        = d_body[idx];
        d_net_force_alt[scan_keep]   = d_net_force[idx];
        d_net_ratedpe_alt[scan_keep] = d_net_ratedpe[idx];
        unsigned int tag = d_tag[idx];
        d_tag_alt[scan_keep] = tag;

        // update rtag
        d_rtag[tag] = scan_keep;
        }
    }

__global__ void
gpu_select_sent_particles(unsigned int N, unsigned int* d_comm_flags, unsigned int* d_tmp)
    {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= N)
        return;
    d_tmp[idx] = d_comm_flags[idx] ? 1 : 0;
    }

unsigned int gpu_pdata_remove(const unsigned int N,
                              const Scalar4* d_pos,
                              const Scalar4* d_vel,
                              const Scalar* d_density,
                              const Scalar* d_pressure,
                              const Scalar* d_energy,
                              const Scalar3* d_aux1,
                              const Scalar3* d_aux2,
                              const Scalar3* d_aux3,
                              const Scalar3* d_aux4,
                              const Scalar* d_slength,
                              const Scalar3* d_accel,
                              const Scalar3* d_dpedt,
                              const int3* d_image,
                              const unsigned int* d_body,
                              const Scalar4* d_net_force,
                              const Scalar4* d_net_ratedpe,
                              const unsigned int* d_tag,
                              unsigned int* d_rtag,
                              Scalar4* d_pos_alt,
                              Scalar4* d_vel_alt,
                              Scalar* d_density_alt,
                              Scalar* d_pressure_alt,
                              Scalar* d_energy_alt,
                              Scalar3* d_aux1_alt,
                              Scalar3* d_aux2_alt,
                              Scalar3* d_aux3_alt,
                              Scalar3* d_aux4_alt,
                              Scalar* d_slength_alt,
                              Scalar3* d_accel_alt,
                              Scalar3* d_dpedt_alt,
                              int3* d_image_alt,
                              unsigned int* d_body_alt,
                              Scalar4* d_net_force_alt,
                              Scalar4* d_net_ratedpe_alt,
                              unsigned int* d_tag_alt,
                              detail::pdata_element* d_out,
                              unsigned int* d_comm_flags,
                              unsigned int* d_comm_flags_out,
                              unsigned int max_n_out,
                              unsigned int* d_tmp,
                              CachedAllocator& alloc)
    {
    if (!N)
        return 0;

    assert(d_pos);       assert(d_vel);
    assert(d_density);   assert(d_pressure);  assert(d_energy);
    assert(d_aux1);      assert(d_aux2);       assert(d_aux3);      assert(d_aux4);
    assert(d_slength);   assert(d_accel);      assert(d_dpedt);
    assert(d_image);     assert(d_body);
    assert(d_net_force); assert(d_net_ratedpe);
    assert(d_tag);       assert(d_rtag);
    assert(d_pos_alt);   assert(d_vel_alt);
    assert(d_density_alt); assert(d_pressure_alt); assert(d_energy_alt);
    assert(d_aux1_alt);  assert(d_aux2_alt);   assert(d_aux3_alt);  assert(d_aux4_alt);
    assert(d_slength_alt); assert(d_accel_alt);  assert(d_dpedt_alt);
    assert(d_image_alt); assert(d_body_alt);
    assert(d_net_force_alt); assert(d_net_ratedpe_alt);
    assert(d_tag_alt);
    assert(d_out);
    assert(d_comm_flags); assert(d_comm_flags_out); assert(d_tmp);

    unsigned int n_out;

    unsigned int block_size = 256;
    unsigned int n_blocks = N / block_size + 1;

    // select nonzero communication flags
    hipLaunchKernelGGL(gpu_select_sent_particles,
                       dim3(n_blocks),
                       dim3(block_size),
                       0,
                       0,
                       N,
                       d_comm_flags,
                       d_tmp);

    // perform an exclusive prefix sum
    void* d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;

    unsigned int* d_scan = alloc.getTemporaryBuffer<unsigned int>(N);
    assert(d_scan);

    hipcub::DeviceScan::ExclusiveSum(d_temp_storage, temp_storage_bytes, d_tmp, d_scan, N);
    d_temp_storage = alloc.getTemporaryBuffer<char>(temp_storage_bytes);
    hipcub::DeviceScan::ExclusiveSum(d_temp_storage, temp_storage_bytes, d_tmp, d_scan, N);
    alloc.deallocate((char*)d_temp_storage);

    // total number of removed particles
    d_temp_storage = NULL;
    temp_storage_bytes = 0;
    unsigned int* d_n_out = (unsigned int*)alloc.getTemporaryBuffer<unsigned int>(1);
    assert(d_n_out);
    hipcub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_tmp, d_n_out, N);
    d_temp_storage = alloc.allocate(temp_storage_bytes);
    hipcub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_tmp, d_n_out, N);
    alloc.deallocate((char*)d_temp_storage);
    hipMemcpy(&n_out, d_n_out, sizeof(unsigned int), hipMemcpyDeviceToHost);
    alloc.deallocate((char*)d_n_out);

    if (n_out <= max_n_out)
        {
        unsigned int nwork = N;
        n_blocks = nwork / block_size + 1;

        hipLaunchKernelGGL(gpu_scatter_particle_data_kernel,
                           dim3(n_blocks),
                           dim3(block_size),
                           0,
                           0,
                           nwork,
                           d_pos, d_vel,
                           d_density, d_pressure, d_energy,
                           d_aux1, d_aux2, d_aux3, d_aux4,
                           d_slength, d_accel, d_dpedt,
                           d_image, d_body,
                           d_net_force, d_net_ratedpe,
                           d_tag, d_rtag,
                           d_pos_alt, d_vel_alt,
                           d_density_alt, d_pressure_alt, d_energy_alt,
                           d_aux1_alt, d_aux2_alt, d_aux3_alt, d_aux4_alt,
                           d_slength_alt, d_accel_alt, d_dpedt_alt,
                           d_image_alt, d_body_alt,
                           d_net_force_alt, d_net_ratedpe_alt,
                           d_tag_alt,
                           d_out,
                           d_comm_flags, d_comm_flags_out,
                           d_scan);
        }

    alloc.deallocate((char*)d_scan);
    return n_out;
    }

__global__ void gpu_pdata_add_particles_kernel(unsigned int old_nparticles,
                                               unsigned int num_add_ptls,
                                               Scalar4* d_pos,
                                               Scalar4* d_vel,
                                               Scalar* d_density,
                                               Scalar* d_pressure,
                                               Scalar* d_energy,
                                               Scalar3* d_aux1,
                                               Scalar3* d_aux2,
                                               Scalar3* d_aux3,
                                               Scalar3* d_aux4,
                                               Scalar* d_slength,
                                               Scalar3* d_accel,
                                               Scalar3* d_dpedt,
                                               int3* d_image,
                                               unsigned int* d_body,
                                               Scalar4* d_net_force,
                                               Scalar4* d_net_ratedpe,
                                               unsigned int* d_tag,
                                               unsigned int* d_rtag,
                                               const detail::pdata_element* d_in,
                                               unsigned int* d_comm_flags)
    {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= num_add_ptls)
        return;

    detail::pdata_element p = d_in[idx];
    unsigned int add_idx = old_nparticles + idx;

    d_pos[add_idx]         = p.pos;
    d_vel[add_idx]         = p.vel;
    d_density[add_idx]     = p.density;
    d_pressure[add_idx]    = p.pressure;
    d_energy[add_idx]      = p.energy;
    d_aux1[add_idx]        = p.aux1;
    d_aux2[add_idx]        = p.aux2;
    d_aux3[add_idx]        = p.aux3;
    d_aux4[add_idx]        = p.aux4;
    d_slength[add_idx]     = p.slength;
    d_accel[add_idx]       = p.accel;
    d_dpedt[add_idx]       = p.dpedt;
    d_image[add_idx]       = p.image;
    d_body[add_idx]        = p.body;
    d_net_force[add_idx]   = p.net_force;
    d_net_ratedpe[add_idx] = p.net_ratedpe;
    d_tag[add_idx]         = p.tag;
    d_rtag[p.tag]          = add_idx;
    d_comm_flags[add_idx]  = 0;
    }

void gpu_pdata_add_particles(const unsigned int old_nparticles,
                             const unsigned int num_add_ptls,
                             Scalar4* d_pos,
                             Scalar4* d_vel,
                             Scalar* d_density,
                             Scalar* d_pressure,
                             Scalar* d_energy,
                             Scalar3* d_aux1,
                             Scalar3* d_aux2,
                             Scalar3* d_aux3,
                             Scalar3* d_aux4,
                             Scalar* d_slength,
                             Scalar3* d_accel,
                             Scalar3* d_dpedt,
                             int3* d_image,
                             unsigned int* d_body,
                             Scalar4* d_net_force,
                             Scalar4* d_net_ratedpe,
                             unsigned int* d_tag,
                             unsigned int* d_rtag,
                             const detail::pdata_element* d_in,
                             unsigned int* d_comm_flags)
    {
    assert(d_pos); assert(d_vel);
    assert(d_density); assert(d_pressure); assert(d_energy);
    assert(d_aux1); assert(d_aux2); assert(d_aux3); assert(d_aux4);
    assert(d_slength); assert(d_accel); assert(d_dpedt);
    assert(d_image); assert(d_body);
    assert(d_net_force); assert(d_net_ratedpe);
    assert(d_tag); assert(d_rtag); assert(d_in); assert(d_comm_flags);

    if (!num_add_ptls)
        return;

    const int block_size = 256;
    hipLaunchKernelGGL(gpu_pdata_add_particles_kernel,
                       dim3(num_add_ptls / block_size + 1),
                       dim3(block_size),
                       0,
                       0,
                       old_nparticles,
                       num_add_ptls,
                       d_pos, d_vel,
                       d_density, d_pressure, d_energy,
                       d_aux1, d_aux2, d_aux3, d_aux4,
                       d_slength, d_accel, d_dpedt,
                       d_image, d_body,
                       d_net_force, d_net_ratedpe,
                       d_tag, d_rtag,
                       d_in,
                       d_comm_flags);
    }

    } // end namespace kernel

    } // end namespace hoomd

#endif // ENABLE_MPI

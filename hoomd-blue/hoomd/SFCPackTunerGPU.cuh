// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __SFC_PACK_UPDATER_GPU_CUH__
#define __SFC_PACK_UPDATER_GPU_CUH__

#include "BoxDim.h"
#include "CachedAllocator.h"
#include "HOOMDMath.h"

/*! \file SFCPackTunerGPU.cuh
    \brief Defines GPU functions for generating the space-filling curve sorted order on the GPU.
   Used by SFCPackTunerGPU.
*/

namespace hoomd
    {
namespace kernel
    {
//! Generate sorted order on GPU
void gpu_generate_sorted_order(unsigned int N,
                               const Scalar4* d_pos,
                               unsigned int* d_particle_bins,
                               unsigned int* d_traversal_order,
                               unsigned int n_grid,
                               unsigned int* d_sorted_order,
                               const BoxDim& box,
                               bool twod,
                               CachedAllocator& alloc);

//! Reorder particle data (GPU driver function) — SPH field set
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
                            unsigned int* d_rtag);

    } // namespace kernel

    } // end namespace hoomd

#endif // __SFC_PACK_UPDATER_GPU_CUH__

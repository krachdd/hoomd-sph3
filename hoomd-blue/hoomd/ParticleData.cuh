// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#pragma once

#ifdef ENABLE_HIP
#include "BoxDim.h"

#include "hoomd/CachedAllocator.h"

/*! \file ParticleData.cuh
    \brief Declares GPU kernel code and data structure functions used by ParticleData
*/

#ifdef __HIPCC__
//! Sentinel value in \a body to signify that this particle does not belong to a body
const unsigned int NO_BODY = 0xffffffff;

//! Unsigned value equivalent to a sign flip in a signed int.
const unsigned int MIN_FLOPPY = 0x80000000;

//! Sentinel value in \a r_tag to signify that this particle is not currently present on the local
//! processor
const unsigned int NOT_LOCAL = 0xffffffff;
#endif

namespace hoomd
    {
namespace detail
    {
#ifdef __HIPCC__
//! Compact particle data storage (SPH version)
struct pdata_element
    {
    Scalar4 pos;          //!< Position
    Scalar4 vel;          //!< Velocity
    Scalar density;       //!< Density
    Scalar pressure;      //!< Pressure
    Scalar energy;        //!< Energy
    Scalar3 aux1;         //!< Auxiliary vector field 1
    Scalar3 aux2;         //!< Auxiliary vector field 2
    Scalar3 aux3;         //!< Auxiliary vector field 3
    Scalar3 aux4;         //!< Auxiliary vector field 4
    Scalar slength;       //!< Smoothing length
    Scalar3 accel;        //!< Acceleration
    Scalar3 dpedt;        //!< Density, pressure and energy rate of change
    int3 image;           //!< Image
    unsigned int body;    //!< Body id
    unsigned int tag;     //!< global tag
    Scalar4 net_force;    //!< net force
    Scalar4 net_ratedpe;  //!< net rate of DPE
    };
#else
//! Forward declaration
struct pdata_element;
#endif

    } // end namespace detail

namespace kernel
    {
//! Pack particle data into output buffer and remove marked particles
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
                              CachedAllocator& alloc);

//! Update particle data with new particles
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
                             unsigned int* d_comm_flags);
    } // end namespace kernel

    } // end namespace hoomd

#endif // ENABLE_HIP

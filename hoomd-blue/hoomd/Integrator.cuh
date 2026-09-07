// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*! \file Integrator.cuh
    \brief Declares methods and data structures used by the Integrator class on the GPU
*/

#ifndef __INTEGRATOR_CUH__
#define __INTEGRATOR_CUH__

#include "HOOMDMath.h"
#include "ParticleData.cuh"

namespace hoomd
    {
namespace kernel
    {
//! struct to pack up several force and ratedpe arrays for addition
/*! To keep the argument count down to gpu_integrator_sum_net_force, up to 6 force/ratedpe array
   pairs are packed up in this struct for addition to the net force/ratedpe in a single kernel
   call. If there is not a multiple of 6 forces to sum, set some of the pointers to NULL and they
   will be ignored.
*/
struct gpu_force_list
    {
    //! Initializes to NULL
    gpu_force_list()
        : f0(NULL), f1(NULL), f2(NULL), f3(NULL), f4(NULL), f5(NULL),
          r0(NULL), r1(NULL), r2(NULL), r3(NULL), r4(NULL), r5(NULL)
        {
        }

    Scalar4* f0; //!< Pointer to force array 0
    Scalar4* f1; //!< Pointer to force array 1
    Scalar4* f2; //!< Pointer to force array 2
    Scalar4* f3; //!< Pointer to force array 3
    Scalar4* f4; //!< Pointer to force array 4
    Scalar4* f5; //!< Pointer to force array 5

    Scalar4* r0; //!< Pointer to ratedpe array 0
    Scalar4* r1; //!< Pointer to ratedpe array 1
    Scalar4* r2; //!< Pointer to ratedpe array 2
    Scalar4* r3; //!< Pointer to ratedpe array 3
    Scalar4* r4; //!< Pointer to ratedpe array 4
    Scalar4* r5; //!< Pointer to ratedpe array 5
    };

//! Driver for gpu_integrator_sum_net_force_kernel()
hipError_t gpu_integrator_sum_net_force(Scalar4* d_net_force,
                                        Scalar4* d_net_ratedpe,
                                        const gpu_force_list& force_list,
                                        unsigned int nparticles,
                                        bool clear);

    } // end namespace kernel

    } // end namespace hoomd

#endif

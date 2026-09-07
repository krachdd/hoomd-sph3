// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "Integrator.cuh"
#include <hip/hip_runtime.h>

#include <assert.h>

/*! \file Integrator.cu
    \brief Defines methods and data structures used by the Integrator class on the GPU
*/

namespace hoomd
    {
namespace kernel
    {
//! Helper to add a force and ratedpe from one ForceCompute to the running totals
__device__ void add_force_total(Scalar4& net_force,
                                Scalar4& net_ratedpe,
                                Scalar4* d_f,
                                Scalar4* d_r,
                                int idx)
    {
    if (d_f != NULL)
        {
        Scalar4 f = d_f[idx];
        net_force.x += f.x;
        net_force.y += f.y;
        net_force.z += f.z;
        net_force.w += f.w;
        }
    if (d_r != NULL)
        {
        Scalar4 r = d_r[idx];
        net_ratedpe.x += r.x;
        net_ratedpe.y += r.y;
        net_ratedpe.z += r.z;
        net_ratedpe.w += r.w;
        }
    }

//! Kernel for summing forces on the GPU
__global__ void gpu_integrator_sum_net_force_kernel(Scalar4* d_net_force,
                                                    Scalar4* d_net_ratedpe,
                                                    const gpu_force_list force_list,
                                                    unsigned int nwork,
                                                    bool clear)
    {
    int idx = blockDim.x * blockIdx.x + threadIdx.x;

    if (idx < nwork)
        {
        Scalar4 net_force;
        Scalar4 net_ratedpe;
        if (clear)
            {
            net_force   = make_scalar4(Scalar(0.0), Scalar(0.0), Scalar(0.0), Scalar(0.0));
            net_ratedpe = make_scalar4(Scalar(0.0), Scalar(0.0), Scalar(0.0), Scalar(0.0));
            }
        else
            {
            net_force   = d_net_force[idx];
            net_ratedpe = d_net_ratedpe[idx];
            }

        add_force_total(net_force, net_ratedpe, force_list.f0, force_list.r0, idx);
        add_force_total(net_force, net_ratedpe, force_list.f1, force_list.r1, idx);
        add_force_total(net_force, net_ratedpe, force_list.f2, force_list.r2, idx);
        add_force_total(net_force, net_ratedpe, force_list.f3, force_list.r3, idx);
        add_force_total(net_force, net_ratedpe, force_list.f4, force_list.r4, idx);
        add_force_total(net_force, net_ratedpe, force_list.f5, force_list.r5, idx);

        d_net_force[idx]   = net_force;
        d_net_ratedpe[idx] = net_ratedpe;
        }
    }

hipError_t gpu_integrator_sum_net_force(Scalar4* d_net_force,
                                        Scalar4* d_net_ratedpe,
                                        const gpu_force_list& force_list,
                                        unsigned int nparticles,
                                        bool clear)
    {
    assert(d_net_force);
    assert(d_net_ratedpe);

    const int block_size = 256;
    unsigned int nwork = nparticles;

    hipLaunchKernelGGL(gpu_integrator_sum_net_force_kernel,
                       dim3(nwork / block_size + 1),
                       dim3(block_size),
                       0,
                       0,
                       d_net_force,
                       d_net_ratedpe,
                       force_list,
                       nwork,
                       clear);

    return hipSuccess;
    }

    } // end namespace kernel

    } // end namespace hoomd

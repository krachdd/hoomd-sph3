// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

//GeneralGPUDeviceFunctions.cu

/*! \file GeneralGPUDeviceFunctions.cu
    \brief __device__ funciton need to be in the same compilation unit as the __global__ kernel
*/

#include "GeneralGPUDeviceFunctions.cuh"

#ifndef __GENERALGPUDEVICEFUNCTIONS_CU__
#define __GENERALGPUDEVICEFUNCTIONS_CU__


// atomicAdd(double*,double) is not defined for cuda verions prior to 6.x
#ifdef __CUDA_ARCH__
#if __CUDA_ARCH__ < 600
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                        __longlong_as_double(assumed)));

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif
#endif

__device__ bool gpu_find(unsigned int *array, unsigned int len, unsigned int val)
{
    bool r(false);
        for(int k(0); k <len; k++){
            if(array[k]==val)
            {
                r=true;
                break;
            }
        }
    return r;
}


// shared should have size of >= blockdim.x
// blockDim should have a power of two
//
// val will be ignored if tid >= count
template<typename T>
__device__ void _gpu_sum(T* shared, unsigned int count, T val, T* dst) {

    unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

    //// Example
    //// val(i) = i
    //// count = 2 * blockdim = 8

    // put all in shared memory TODO first reduction step
    if (tid < count) {
        shared[tid] = val;
    } else {
        shared[tid] = 0;
    }
    __syncthreads();

    //// now we have
    //// [0, 1, 2, 3, 4, 5, 6, 7]

    // Reduction
    int stride;
    for (stride = blockDim.x / 2; stride >= 1; stride /= 2) {
        if (tid < stride) {
            shared[tid] += shared[tid + stride];
        }
        __syncthreads();
    }
    //// this leads to
    //// [0+4, 1+5, 2+6, 3+7] = [4, 6, 8, 10]
    //// [4+8, 6+10, x, x] = [12, 16, x, x]
    //// [12+16, x, x, x] = [28, x, x, x]

    if (tid == 0) {
        atomicAdd(dst, shared[0]);
    }
}
template __device__ void _gpu_sum<double>(double* shared, unsigned int count, double val, double* dst);
#endif
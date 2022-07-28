// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

//GeneralGPUDeviceFunctions.cu

/*! \file GeneralGPUDeviceFunctions.cu
    \brief __device__ funciton need to be in the same compilation unit as the __global__ kernel
*/

#ifndef __GENERALGPUDEVICEFUNCTIONS_CUH__
#define __GENERALGPUDEVICEFUNCTIONS_CUH__

// atomicAdd(double*,double) is not defined for cuda verions prior to 6.x
__device__ double atomicAdd(double* address, double val);


__device__ bool gpu_find(unsigned int *array, unsigned int len, unsigned int val);

// shared should have size of >= blockdim.x
// blockDim should have a power of two
//
// val will be ignored if tid >= count
template<typename T>
__device__ void _gpu_sum(T* shared, unsigned int count, T val, T* dst);

#endif
// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

//GeneralGPUFunctions.cuh

#ifndef __GENERALGPUFUNCTIONS_CUH__
#define __GENERALGPUFUNCTIONS_CUH__

#include <stdio.h>
#include "hoomd/HOOMDMath.h"

#define REDUCTION_BLOCK_SIZE 512


__global__ void gpu_SimpleFind_kernel(unsigned int *pszData, unsigned int pszTarget, unsigned int* pFound, unsigned int Len);


template<typename Sin>
cudaError_t gpu_sum_xyzw(Sin* input, unsigned int len, Scalar4* sum, unsigned int* indices = 0);

template<typename Sin>
cudaError_t gpu_sum_xyz(Sin* input, unsigned int len, Scalar3* sum, unsigned int* indices = 0);

template<typename Sin>
cudaError_t gpu_sum_xy(Sin* input, unsigned int len, Scalar2* sum, unsigned int* indices = 0);

template<typename Sin>
cudaError_t gpu_sum_x(Sin* input, unsigned int len, Scalar* sum, unsigned int* indices = 0);

#endif
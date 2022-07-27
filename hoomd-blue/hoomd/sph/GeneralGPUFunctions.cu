// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

//GeneralGPUFunctions.cu

#include <stdio.h>
#include "GeneralGPUFunctions.cuh"
#include "GeneralGPUDeviceFunctions.cuh"

//unsigned int
__global__ void gpu_SimpleFind_kernel(unsigned int *pszData, unsigned int pszTarget, unsigned int* pFound, unsigned int Len)
{
    __shared__ unsigned int test;
    unsigned int startIndex = blockDim.x*blockIdx.x + threadIdx.x;

    if(threadIdx.x==0)
        test = *pFound;
    __syncthreads();

    if(startIndex < Len)
    {

        if(*pFound > startIndex)
        {
            // only continue if an earlier instance hasn't already been found
            if (pszData[startIndex] == pszTarget){
                atomicMin(&test, startIndex);
                atomicMin(pFound, startIndex);
            }

        }
    }
    __syncthreads();

}
/*
__device__ int weh;
__global__ void test(int *pszData, int Len, int* w)
{
        *w = Len;
        int *was;
        (cudaMalloc((void**)&was, sizeof(int)));
        *was=Len;
        weh=Len;
        printf("KERNEL KERNEL LAUNCH: POS: %i\t %p\n", *was, was);
        int numBlocks = ceil((float)Len/(REDUCTION_BLOCK_SIZE));
        dim3 dimGrid(numBlocks, 1, 1);
        dim3 dimBlock(REDUCTION_BLOCK_SIZE, 1, 1);
        SimpleFindKind<<<dimGrid, dimBlock>>>(pszData, 1909, was, Len);
        printf("");
        printf("KERNEL KERNEL LAUNCH: POS: %i\t%p\n", *was, was);
        printf("KERNEL KERNEL LAUNCH: POS: %i\t%p\n", *was, was);
        *w=*was;
}
*/

struct IndexerLinear {

    IndexerLinear() {
    }

    __host__ __device__ inline
    unsigned int ix(unsigned int i) {
        return i;
    }
};

struct IndexerIndirect {
    unsigned int* indirect;

    IndexerIndirect(unsigned int* indices) : indirect(indices) {
    }

    __host__ __device__ inline
    unsigned int ix(unsigned int i) {
        return indirect[i];
    }
};


// reduce(+) operation  on `input`. `len` is the number of ScalarN in input
// the result is _added_ to `sum`.
//
// Algorithm:
// chunks of 2*REDUCTION_BLOCK_SIZE elements are loaded into shared memory.
// then a scan operation is performed in shared memory
// and the last element (== sum of all input values in the current chunk) is added to the total result
template<typename Sin, typename Ix>
__global__ void gpu_sum_xyzw_kernel(Sin * input, Ix indexer, unsigned int len, Scalar4 *sum) {
    //assert(blockDim.x == REDUCTION_BLOCK_SIZE);

    // Load a segment of the input vector into shared memory
    __shared__ Scalar scan_array_x[REDUCTION_BLOCK_SIZE];
    __shared__ Scalar scan_array_y[REDUCTION_BLOCK_SIZE];
    __shared__ Scalar scan_array_z[REDUCTION_BLOCK_SIZE];
    __shared__ Scalar scan_array_w[REDUCTION_BLOCK_SIZE];

    // start is the offset<ScalarN> for the current chunk
    unsigned int t = threadIdx.x, start = 2 * blockIdx.x * REDUCTION_BLOCK_SIZE;

    if (start + t < len) {
        scan_array_x[t] = input[indexer.ix(start + t)].x;
        scan_array_y[t] = input[indexer.ix(start + t)].y;
        scan_array_z[t] = input[indexer.ix(start + t)].z;
        scan_array_w[t] = input[indexer.ix(start + t)].w;
    } else {
        scan_array_x[t] = 0;
        scan_array_y[t] = 0;
        scan_array_z[t] = 0;
        scan_array_w[t] = 0;
    }
    if (start + REDUCTION_BLOCK_SIZE + t < len) {
        scan_array_x[t] += input[indexer.ix(start + REDUCTION_BLOCK_SIZE + t)].x;
        scan_array_y[t] += input[indexer.ix(start + REDUCTION_BLOCK_SIZE + t)].y;
        scan_array_z[t] += input[indexer.ix(start + REDUCTION_BLOCK_SIZE + t)].z;
        scan_array_w[t] += input[indexer.ix(start + REDUCTION_BLOCK_SIZE + t)].w;
    }
    __syncthreads();

    // Reduction
    int stride;
    for (stride = REDUCTION_BLOCK_SIZE / 2; stride >= 1; stride /= 2) {
        if (t < stride) {
            scan_array_x[t] += scan_array_x[t + stride];
            scan_array_y[t] += scan_array_y[t + stride];
            scan_array_z[t] += scan_array_z[t + stride];
            scan_array_w[t] += scan_array_w[t + stride];
        }
        __syncthreads();
    }

    if (t == 0) {
        atomicAdd(&(sum->x), scan_array_x[0]);
        atomicAdd(&(sum->y), scan_array_y[0]);
        atomicAdd(&(sum->z), scan_array_z[0]);
        atomicAdd(&(sum->w), scan_array_w[0]);
    }
}

template<typename Sin, typename Ix>
__global__ void gpu_sum_xyz_kernel(Sin * input, Ix indexer, unsigned int len, Scalar3 *sum) {
    //assert(blockDim.x == REDUCTION_BLOCK_SIZE);

    // Load a segment of the input vector into shared memory
    __shared__ Scalar scan_array_x[REDUCTION_BLOCK_SIZE];
    __shared__ Scalar scan_array_y[REDUCTION_BLOCK_SIZE];
    __shared__ Scalar scan_array_z[REDUCTION_BLOCK_SIZE];

    // start is the offset<ScalarN> for the current chunk
    unsigned int t = threadIdx.x;
    unsigned int start = 2 * blockIdx.x * REDUCTION_BLOCK_SIZE;

    // first reduction step
    if (start + t < len) {
        scan_array_x[t] = input[indexer.ix(start + t)].x;
        scan_array_y[t] = input[indexer.ix(start + t)].y;
        scan_array_z[t] = input[indexer.ix(start + t)].z;
    } else {
        scan_array_x[t] = 0;
        scan_array_y[t] = 0;
        scan_array_z[t] = 0;
    }
    if (start + REDUCTION_BLOCK_SIZE + t < len) {
        scan_array_x[t] += input[indexer.ix(start + REDUCTION_BLOCK_SIZE + t)].x;
        scan_array_y[t] += input[indexer.ix(start + REDUCTION_BLOCK_SIZE + t)].y;
        scan_array_z[t] += input[indexer.ix(start + REDUCTION_BLOCK_SIZE + t)].z;
    }
    __syncthreads();

    // Reduction
    int stride;
    for (stride = REDUCTION_BLOCK_SIZE / 2; stride >= 1; stride /= 2) {
        if (t < stride) {
            scan_array_x[t] += scan_array_x[t + stride];
            scan_array_y[t] += scan_array_y[t + stride];
            scan_array_z[t] += scan_array_z[t + stride];
        }
        __syncthreads();
    }

    if (t == 0) {
        atomicAdd(&(sum->x), scan_array_x[0]);
        atomicAdd(&(sum->y), scan_array_y[0]);
        atomicAdd(&(sum->z), scan_array_z[0]);
    }
}

template<typename Sin, typename Ix>
__global__ void gpu_sum_xy_kernel(Sin * input, Ix indexer, unsigned int len, Scalar2 *sum) {
    //assert(blockDim.x == REDUCTION_BLOCK_SIZE);

    // Load a segment of the input vector into shared memory
    __shared__ Scalar scan_array_x[REDUCTION_BLOCK_SIZE];
    __shared__ Scalar scan_array_y[REDUCTION_BLOCK_SIZE];

    // start is the offset<ScalarN> for the current chunk
    unsigned int t = threadIdx.x;
    unsigned int start = 2 * blockIdx.x * REDUCTION_BLOCK_SIZE;

    // first reduction step
    if (start + t < len) {
        scan_array_x[t] = input[indexer.ix(start + t)].x;
        scan_array_y[t] = input[indexer.ix(start + t)].y;
    } else {
        scan_array_x[t] = 0;
        scan_array_y[t] = 0;
    }
    if (start + REDUCTION_BLOCK_SIZE + t < len) {
        scan_array_x[t] += input[indexer.ix(start + REDUCTION_BLOCK_SIZE + t)].x;
        scan_array_y[t] += input[indexer.ix(start + REDUCTION_BLOCK_SIZE + t)].y;
    }
    __syncthreads();

    // Reduction
    int stride;
    for (stride = REDUCTION_BLOCK_SIZE / 2; stride >= 1; stride /= 2) {
        if (t < stride) {
            scan_array_x[t] += scan_array_x[t + stride];
            scan_array_y[t] += scan_array_y[t + stride];
        }
        __syncthreads();
    }

    if (t == 0) {
        atomicAdd(&(sum->x), scan_array_x[0]);
        atomicAdd(&(sum->y), scan_array_y[0]);
    }
}


// this is needed to use the first element in both Scalar and ScalarN in templated code
template<typename S> __device__ inline
Scalar x(S s) {
    return s.x;
}

template<> __device__ inline
Scalar x<Scalar>(Scalar s) {
    return s;
}

template<typename Sin, typename Ix>
__global__ void gpu_sum_x_kernel(Sin * input, Ix indexer, unsigned int len, Scalar *sum) {
    //assert(blockDim.x == REDUCTION_BLOCK_SIZE);

    // Load a segment of the input vector into shared memory
    __shared__ Scalar scan_array_x[REDUCTION_BLOCK_SIZE];

    // start is the offset<ScalarN> for the current chunk
    unsigned int t = threadIdx.x;
    unsigned int start = 2 * blockIdx.x * REDUCTION_BLOCK_SIZE;

    // first reduction step
    if (start + t < len) {
        scan_array_x[t] = x(input[indexer.ix(start + t)]);
    } else {
        scan_array_x[t] = 0;
    }
    if (start + REDUCTION_BLOCK_SIZE + t < len) {
        // no need to sync as index t is only accessed by thread t
        scan_array_x[t] += x(input[indexer.ix(start + REDUCTION_BLOCK_SIZE + t)]);
    }
    __syncthreads();

    // Reduction
    int stride;
    for (stride = REDUCTION_BLOCK_SIZE / 2; stride >= 1; stride /= 2) {
        if (t < stride) {
            scan_array_x[t] += scan_array_x[t + stride];
        }
        __syncthreads();
    }

    if (t == 0) {
        atomicAdd(sum, scan_array_x[0]);
    }
}


template<typename Sin>
cudaError_t gpu_sum_xyzw(Sin* input, unsigned int len, Scalar4* sum, unsigned int* indices) {
    unsigned int block_size=REDUCTION_BLOCK_SIZE;
    dim3 grid( (len / block_size / 2) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);

    if (!indices) {
        gpu_sum_xyzw_kernel<<<grid, threads>>>(input, IndexerLinear(), len, sum);
    } else {
        // TODO better way to do indirect indices
        gpu_sum_xyzw_kernel<<<grid, threads>>>(input, IndexerIndirect(indices), len, sum);
    }

    return cudaSuccess;
}

template<typename Sin>
cudaError_t gpu_sum_xyz(Sin* input, unsigned int len, Scalar3* sum, unsigned int* indices) {
    unsigned int block_size=REDUCTION_BLOCK_SIZE;
    dim3 grid( (len / block_size / 2) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);

    if (!indices) {
        gpu_sum_xyz_kernel<<<grid, threads>>>(input, IndexerLinear(), len, sum);
    } else {
        // TODO better way to do indirect indices
        gpu_sum_xyz_kernel<<<grid, threads>>>(input, IndexerIndirect(indices), len, sum);
    }

    return cudaSuccess;
}

template<typename Sin>
cudaError_t gpu_sum_xy(Sin* input, unsigned int len, Scalar2* sum, unsigned int* indices) {
    unsigned int block_size=REDUCTION_BLOCK_SIZE;
    dim3 grid( (len / block_size / 2) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);

    if (!indices) {
        gpu_sum_xy_kernel<<<grid, threads>>>(input, IndexerLinear(), len, sum);
    } else {
        // TODO better way to do indirect indices
        gpu_sum_xy_kernel<<<grid, threads>>>(input, IndexerIndirect(indices), len, sum);
    }

    return cudaSuccess;
}

template<typename Sin>
cudaError_t gpu_sum_x(Sin* input, unsigned int len, Scalar* sum, unsigned int* indices) {
    unsigned int block_size=REDUCTION_BLOCK_SIZE;
    dim3 grid( (len / block_size / 2) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);

    if (!indices) {
        gpu_sum_x_kernel<<<grid, threads>>>(input, IndexerLinear(), len, sum);
    } else {
        // TODO better way to do indirect indices
        gpu_sum_x_kernel<<<grid, threads>>>(input, IndexerIndirect(indices), len, sum);
    }

    return cudaSuccess;
}

template cudaError_t gpu_sum_xyz<Scalar4>(Scalar4* input, unsigned int len, Scalar3* sum, unsigned int* indices);
template cudaError_t gpu_sum_x<Scalar4>(Scalar4* input, unsigned int len, Scalar* sum, unsigned int* indices);

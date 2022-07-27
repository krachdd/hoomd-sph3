// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

//GeneralGPUFunctions_test.cu

// run with: nvcc -O3 --std c++11 --compiler-options -O3,-march=native  -I../.. GeneralGpuFunctions_test.cu --run

#include <iostream>
#include "GeneralGPUFunctions.cu"
#include "GeneralGPUDeviceFunctions.cu"

void checkForErrors(const cudaError_t status, const char *label, const int line, const char *file)
{
    if (status != cudaSuccess) {
        std::cerr << "CUDA ERROR (" << label << ") ";
        std::cerr << "at " << file << ":" << line << std::endl;
        std::cerr << cudaGetErrorString(status) << ". Exiting..." << std::endl;
        exit(1);
    }
}

#define cuchck(func) checkForErrors(func, #func, __LINE__, __FILE__)

int tests_failed = 0;
void assert_small(const char* name, double val, double bound = 1e-10) {
    if (abs(val) > bound) {
        std::cout << "! Test FAILED: " << name << "\twith value=" << val << std::endl;
        tests_failed++;
    } else {
        std::cout << "  Test passed: " << name << std::endl;
    }
}

int main() {
    size_t N = 6000000;
    double* data = new double[N];

    double expected_sum_x = N*(N-1)/2.0; // sum [0,N)
    double expected_sum_x0 = 2*((N/2)*((N/2)-1)/2.0); // 2 * sum [0,N/2)
    double expected_sum_x00 = 3*((N/3)*((N/3)-1)/2.0); // 3 * sum [0,N/3)
    double expected_sum_x000 = 4*((N/4)*((N/4)-1)/2.0); // 4 * sum [0,N/4)
    for (size_t i = 0; i < N; i++) {
        data[i] = i;
    }

    double* d_data = 0;
    cuchck(cudaMalloc(&d_data, sizeof(double) * N));
    cuchck(cudaMemcpy(d_data, data, sizeof(double) * N, cudaMemcpyHostToDevice));
    Scalar* d_sum_x = 0;
    Scalar2* d_sum_xy = 0;
    Scalar3* d_sum_xyz = 0;
    Scalar4* d_sum_xyzw = 0;
    cuchck(cudaMalloc(&d_sum_x, sizeof(Scalar)));
    cuchck(cudaMalloc(&d_sum_xy, sizeof(Scalar2)));
    cuchck(cudaMalloc(&d_sum_xyz, sizeof(Scalar3)));
    cuchck(cudaMalloc(&d_sum_xyzw, sizeof(Scalar4)));
    Scalar sum_x = 0;
    Scalar2 sum_xy = make_scalar2(0,0);
    Scalar3 sum_xyz = make_scalar3(0,0,0);
    Scalar4 sum_xyzw = make_scalar4(0,0,0,0);
    const Scalar4 zero4 = make_scalar4(0,0,0,0);

    // ZERO accumulators
    cuchck(cudaMemcpy(d_sum_x, &zero4, sizeof(Scalar), cudaMemcpyHostToDevice));
    cuchck(cudaMemcpy(d_sum_xy, &zero4, sizeof(Scalar2), cudaMemcpyHostToDevice));
    cuchck(cudaMemcpy(d_sum_xyz, &zero4, sizeof(Scalar3), cudaMemcpyHostToDevice));
    cuchck(cudaMemcpy(d_sum_xyzw, &zero4, sizeof(Scalar4), cudaMemcpyHostToDevice));

    // T -> T reduction
    std::cout << "Test Reduction: 1,2,3,4-dim" << std::endl;
    gpu_sum_x((Scalar*)d_data, N, d_sum_x, 0);
    gpu_sum_xy((Scalar2*)d_data, N/2, d_sum_xy, 0);
    gpu_sum_xyz((Scalar3*)d_data, N/3, d_sum_xyz, 0);
    gpu_sum_xyzw((Scalar4*)d_data, N/4, d_sum_xyzw, 0);

    cuchck(cudaMemcpy(&sum_x, d_sum_x, sizeof(Scalar), cudaMemcpyDeviceToHost));
    cuchck(cudaMemcpy(&sum_xy, d_sum_xy, sizeof(Scalar2), cudaMemcpyDeviceToHost));
    cuchck(cudaMemcpy(&sum_xyz, d_sum_xyz, sizeof(Scalar3), cudaMemcpyDeviceToHost));
    cuchck(cudaMemcpy(&sum_xyzw, d_sum_xyzw, sizeof(Scalar4), cudaMemcpyDeviceToHost));

    assert_small("sum_x", sum_x - expected_sum_x);
    assert_small("sum_xy x", sum_xy.x - expected_sum_x0);
    assert_small("sum_xy y", sum_xy.y - expected_sum_x0 - N/2);
    assert_small("sum_xyz x", sum_xyz.x - expected_sum_x00);
    assert_small("sum_xyz y", sum_xyz.y - expected_sum_x00 - N/3);
    assert_small("sum_xyz z", sum_xyz.z - expected_sum_x00 - 2*N/3);
    assert_small("sum_xyzw x", sum_xyzw.x - expected_sum_x000);
    assert_small("sum_xyzw x", sum_xyzw.y - expected_sum_x000 - N/4);
    assert_small("sum_xyzw x", sum_xyzw.z - expected_sum_x000 - 2*N/4);
    assert_small("sum_xyzw x", sum_xyzw.w - expected_sum_x000 - 3*N/4);


    // ZERO accumulators
    cuchck(cudaMemcpy(d_sum_x, &zero4, sizeof(Scalar), cudaMemcpyHostToDevice));
    cuchck(cudaMemcpy(d_sum_xy, &zero4, sizeof(Scalar2), cudaMemcpyHostToDevice));
    cuchck(cudaMemcpy(d_sum_xyz, &zero4, sizeof(Scalar3), cudaMemcpyHostToDevice));
    cuchck(cudaMemcpy(d_sum_xyzw, &zero4, sizeof(Scalar4), cudaMemcpyHostToDevice));

    // subset reduction
    std::cout << "Test Reduction on Subsets: 2->1, 3->2, 4->3 dim" << std::endl;
    gpu_sum_x((Scalar2*)d_data, N/2, d_sum_x, 0);
    gpu_sum_xy((Scalar3*)d_data, N/3, d_sum_xy, 0);
    gpu_sum_xyz((Scalar4*)d_data, N/4, d_sum_xyz, 0);

    cuchck(cudaMemcpy(&sum_x, d_sum_x, sizeof(Scalar), cudaMemcpyDeviceToHost));
    cuchck(cudaMemcpy(&sum_xy, d_sum_xy, sizeof(Scalar2), cudaMemcpyDeviceToHost));
    cuchck(cudaMemcpy(&sum_xyz, d_sum_xyz, sizeof(Scalar3), cudaMemcpyDeviceToHost));
    cuchck(cudaMemcpy(&sum_xyzw, d_sum_xyzw, sizeof(Scalar4), cudaMemcpyDeviceToHost));

    assert_small("sum_x", sum_x - expected_sum_x0);
    assert_small("sum_xy x", sum_xy.x - expected_sum_x00);
    assert_small("sum_xy y", sum_xy.y - expected_sum_x00 - N/3);
    assert_small("sum_xyz x", sum_xyz.x - expected_sum_x000);
    assert_small("sum_xyz y", sum_xyz.y - expected_sum_x000 - N/4);
    assert_small("sum_xyz z", sum_xyz.z - expected_sum_x000 - 2*N/4);


    cuchck(cudaFree(d_sum_xyzw));
    cuchck(cudaFree(d_sum_xyz));
    cuchck(cudaFree(d_sum_xy));
    cuchck(cudaFree(d_sum_x));
    cuchck(cudaFree(d_data));
    delete[] data;

    return tests_failed;
}
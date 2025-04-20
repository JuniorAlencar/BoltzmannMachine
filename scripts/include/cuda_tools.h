#ifndef CUDA_TOOLS_H
#define CUDA_TOOLS_H

#include <cuda_runtime.h>
#include <iostream>

inline void check_cuda_error(const char* msg = "") {
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: " << cudaGetErrorString(err) << " " << msg << std::endl;
        exit(EXIT_FAILURE);
    }
}

#endif

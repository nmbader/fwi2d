#pragma once

#include <stdexcept>
#include <cuda.h>
#include <cuda_runtime_api.h>

#define cudaCheckError(ans) { _cudaCheckError((ans), __FILE__, __LINE__); }
#define cudaKernelError() { _cudaKernelError(__FILE__, __LINE__); }

inline void _cudaCheckError(cudaError_t code, const char *file, int line)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"\n========CUDA API ERROR in file %s, line %d========\n",file,line);
        throw std::logic_error(cudaGetErrorString(code));
    }
}
inline void _cudaKernelError(const char *file, int line)
{
    cudaError err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr,"\n========CUDA KERNEL ERROR in file %s, line %d========\n",file,line);
        throw std::logic_error(cudaGetErrorString(err));
    }
}

inline void printGpuInfo(int iGpu){
    cudaDeviceProp  prop;
    int count;
    cudaGetDeviceCount( &count );
    fprintf(stderr,"Number of GPUs available: %d\n",count);
    cudaGetDeviceProperties( &prop, iGpu );
    fprintf(stderr, "   --- General Information for device %d ---\n", iGpu );
    fprintf(stderr, "Name:  %s\n", prop.name );
    fprintf(stderr, "Compute capability:  %d.%d\n", prop.major, prop.minor );
    fprintf(stderr, "Clock rate:  %d\n", prop.clockRate );
    fprintf(stderr, "Device copy overlap:  " );
    if (prop.deviceOverlap) fprintf(stderr, "Enabled\n" );
    else fprintf(stderr, "Disabled\n" );
    fprintf(stderr, "Kernel execition timeout :  " );
    if (prop.kernelExecTimeoutEnabled) fprintf(stderr, "Enabled\n" );
    else fprintf(stderr, "Disabled\n" );
    fprintf(stderr, "   --- Memory Information for device %d ---\n", iGpu );
    fprintf(stderr, "Total global mem:  %ld\n", prop.totalGlobalMem );
    fprintf(stderr, "Total constant Mem:  %ld\n", prop.totalConstMem );
    fprintf(stderr, "Max mem pitch:  %ld\n", prop.memPitch );
    fprintf(stderr, "Texture Alignment:  %ld\n", prop.textureAlignment );
    fprintf(stderr, "   --- MP Information for device %d ---\n", iGpu );
    fprintf(stderr, "Multiprocessor count:  %d\n",prop.multiProcessorCount );
    fprintf(stderr, "Shared mem per mp:  %ld\n", prop.sharedMemPerBlock );
    fprintf(stderr, "Registers per mp:  %d\n", prop.regsPerBlock );
    fprintf(stderr, "Threads in warp:  %d\n", prop.warpSize );
    fprintf(stderr, "Max threads per block:  %d\n",prop.maxThreadsPerBlock );
    fprintf(stderr, "Max thread dimensions:  (%d, %d, %d)\n",prop.maxThreadsDim[0], prop.maxThreadsDim[1],prop.maxThreadsDim[2] );
    fprintf(stderr, "Max grid dimensions:  (%d, %d, %d)\n",prop.maxGridSize[0], prop.maxGridSize[1],prop.maxGridSize[2] );
    fprintf(stderr, "\n" );
}
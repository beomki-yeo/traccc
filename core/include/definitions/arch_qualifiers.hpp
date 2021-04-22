/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

// CUDA macros
#ifdef __CUDACC__      
#include <cuda.h>
#include <cuda_runtime.h>
#define __CUDA_QUALIFIER__ inline __device__ __host__
#else
#define __CUDA_QUALIFIER__
#endif   

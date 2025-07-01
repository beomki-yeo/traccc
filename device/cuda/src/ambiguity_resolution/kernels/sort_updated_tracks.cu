/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/barrier.hpp"
#include "../../utils/global_index.hpp"
#include "sort_updated_tracks.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void sort_updated_tracks(
    device::sort_updated_tracks_payload payload) {

    if (*(payload.terminate) == 1 || *(payload.n_updated_tracks) == 0) {
        return;
    }

    __shared__ unsigned int N;
    extern __shared__ unsigned int shared_mem_tracks[];

    vecmem::device_vector<const traccc::scalar> rel_shared(
        payload.rel_shared_view);
    vecmem::device_vector<const traccc::scalar> pvals(payload.pvals_view);
    vecmem::device_vector<unsigned int> updated_tracks(
        payload.updated_tracks_view);

    const unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid == 0) {
        N = (*(payload.n_updated_tracks) == 0)
                ? 1
                : 1 << (32 - __clz(*(payload.n_updated_tracks) - 1));
    }

    __syncthreads();

    // Load to shared memory
    if (tid < *(payload.n_updated_tracks)) {
        shared_mem_tracks[tid] = updated_tracks[tid];
    } else {
        shared_mem_tracks[tid] = UINT_MAX;  // dummy index
    }

    __syncthreads();

    if (tid == 0) {
        printf("raw shared_mem_tracks: ");
        for (std::size_t i = 0; i < *(payload.n_updated_tracks); ++i) {
            printf("%u ", shared_mem_tracks[i]);
        }
        printf("\n");
    }

    if (tid == 0) {
        printf("Before sorting \n");
        for (std::size_t i = 0; i < *(payload.n_updated_tracks); i++) {
            unsigned int trk_i = shared_mem_tracks[i];
            traccc::scalar rel_i = rel_shared[trk_i];
            traccc::scalar pv_i = pvals[trk_i];
            printf("(%f %f) ", rel_i, pv_i);
        }
        printf("\n");
    }

    __syncthreads();

    // Bitonic sort on shared_mem_tracks using rel_shared & pvals
    for (unsigned int k = 2; k <= N; k <<= 1) {
        for (unsigned int j = k >> 1; j > 0; j >>= 1) {
            unsigned int ixj = tid ^ j;

            if (ixj > tid && ixj < N && tid < N) {
                // printf("%d %d %d \n", ixj, tid, N);

                unsigned int trk_i = shared_mem_tracks[tid];
                unsigned int trk_j = shared_mem_tracks[ixj];

                bool should_swap = false;
                if (trk_i == UINT_MAX && trk_j != UINT_MAX) {
                    should_swap = true;
                } else if (trk_i != UINT_MAX && trk_j == UINT_MAX) {
                    should_swap = false;
                } else if (trk_i == UINT_MAX && trk_j == UINT_MAX) {
                    should_swap = false;
                } else {

                    traccc::scalar rel_i = rel_shared[trk_i];
                    traccc::scalar rel_j = rel_shared[trk_j];
                    traccc::scalar pv_i = pvals[trk_i];
                    traccc::scalar pv_j = pvals[trk_j];

                    bool ascending = ((tid & k) == 0);
                    should_swap = (rel_i > rel_j || (rel_i == rel_j &&
                                                     pv_i < pv_j)) == ascending;
                }
                
                if (should_swap) {

                    shared_mem_tracks[tid] = trk_j;
                    shared_mem_tracks[ixj] = trk_i;
                }
            }

            __syncthreads();
        }
    }

    __syncthreads();

    if (tid == 0) {
        printf("After sorting \n");
        for (std::size_t i = 0; i < *(payload.n_updated_tracks); i++) {
            unsigned int trk_i = shared_mem_tracks[i];
            // printf("%d \n", trk_i);
            traccc::scalar rel_i = rel_shared[trk_i];
            traccc::scalar pv_i = pvals[trk_i];
            printf("(%f %f) ", rel_i, pv_i);
        }
        printf("\n");
    }

    if (tid < *(payload.n_updated_tracks)) {
        updated_tracks[tid] = shared_mem_tracks[tid];
    }
}

}  // namespace traccc::cuda::kernels

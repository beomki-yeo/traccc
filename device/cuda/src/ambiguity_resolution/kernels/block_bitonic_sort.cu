/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "block_bitonic_sort.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void block_bitonic_sort(device::block_bitonic_sort_payload payload) {

    if (*(payload.terminate) == 1 || *(payload.n_meas_to_remove) <= 1) {
        return;
    }

    vecmem::device_vector<measurement_id_type> meas_to_remove(
        payload.meas_to_remove_view);
    vecmem::device_vector<unsigned int> threads(payload.threads_view);

    __shared__ measurement_id_type sh_meas_ids[512];
    __shared__ unsigned int sh_threads[512];
    __shared__ unsigned int N;

    auto threadIndex = threadIdx.x;

    sh_meas_ids[threadIndex] = std::numeric_limits<measurement_id_type>::max();
    sh_threads[threadIndex] = std::numeric_limits<unsigned int>::max();

    if (threadIndex <= *(payload.n_meas_to_remove)){
        sh_meas_ids[threadIndex] = meas_to_remove[threadIndex];
        sh_threads[threadIndex] = threads[threadIndex];
    }

    if (threadIndex == 0) {

        // Padding N to the power of 2
        N = 1 << (32 - __clz(*(payload.n_meas_to_remove) - 1));
    }
    __syncthreads();

    // Bitonic sort on meas_to_thread w.r.t. measurement id
    const auto tid = threadIdx.x;
    for (int k = 2; k <= N; k <<= 1) {
        for (int j = k >> 1; j > 0; j >>= 1) {
            int ixj = tid ^ j;

            if (ixj > tid && ixj < N && tid < N) {
                auto meas_i = sh_meas_ids[tid];
                auto meas_j = sh_meas_ids[ixj];
                auto thread_i = sh_threads[tid];
                auto thread_j = sh_threads[ixj];

                bool ascending = ((tid & k) == 0);
                bool should_swap =
                    (meas_i > meas_j ||
                     (meas_i == meas_j && thread_i > thread_j)) == ascending;

                if (should_swap) {
                    sh_meas_ids[tid] = meas_j;
                    sh_meas_ids[ixj] = meas_i;
                    sh_threads[tid] = thread_j;
                    sh_threads[ixj] = thread_i;
                }
            }
            __syncthreads();
        }
    }
}

}  // namespace traccc::cuda::kernels

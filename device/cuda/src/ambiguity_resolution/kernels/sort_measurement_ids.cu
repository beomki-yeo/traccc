/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "sort_measurement_ids.cuh"

// VecMem include(s).
#include <vecmem/containers/jagged_device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void sort_measurement_ids(
    device::sort_measurement_ids_payload payload) {

    __shared__ unsigned int sh_meas_ids[1024];

    vecmem::jagged_device_vector<measurement_id_type> all_meas_ids(
        payload.meas_ids_view);

    auto meas_ids = all_meas_ids.at(blockIdx.x);
    const unsigned int tid = threadIdx.x;
    const unsigned int n_meas = meas_ids.size();

    sh_meas_ids[tid] = std::numeric_limits<measurement_id_type>::max();

    if (tid < n_meas) {
        sh_meas_ids[tid] = meas_ids[tid];
    }

    // Bitonic sort
    const unsigned int N = 1 << (32 - __clz(n_meas - 1));
    for (int k = 2; k <= N; k <<= 1) {

        bool ascending = ((tid & k) == 0);

        for (int j = k >> 1; j > 0; j >>= 1) {
            int ixj = tid ^ j;

            if (ixj > tid && ixj < N && tid < N) {
                auto meas_i = sh_meas_ids[tid];
                auto meas_j = sh_meas_ids[ixj];

                bool should_swap = (meas_i > meas_j) == ascending;

                if (should_swap) {
                    sh_meas_ids[tid] = meas_j;
                    sh_meas_ids[ixj] = meas_i;
                }
            }
            __syncthreads();
        }
    }

    if (tid < n_meas) {
        meas_ids[tid] = sh_meas_ids[tid];
    }
}
}  // namespace traccc::cuda::kernels

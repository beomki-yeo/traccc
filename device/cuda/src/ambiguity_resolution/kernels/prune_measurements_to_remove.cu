/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "prune_measurements_to_remove.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void prune_measurements_to_remove(
    device::prune_measurements_to_remove_payload payload) {

    if (*(payload.terminate) == 1) {
        return;
    }

    auto threadIndex = threadIdx.x;

    __shared__ measurement_id_type sh_meas_ids[512];
    __shared__ unsigned int sh_threads[512];
    __shared__ int prefix[512];
    __shared__ unsigned int min_thread;

    if (threadIndex == 0) {
        min_thread = std::numeric_limits<unsigned int>::max();
    }

    vecmem::device_vector<const unsigned int> n_accepted_tracks_per_measurement(
        payload.n_accepted_tracks_per_measurement_view);
    vecmem::device_vector<const unsigned int> meas_id_to_unique_id(
        payload.meas_id_to_unique_id_view);
    vecmem::device_vector<measurement_id_type> meas_to_remove(
        payload.meas_to_remove_view);
    vecmem::device_vector<unsigned int> threads(payload.threads_view);

    sh_meas_ids[threadIndex] = std::numeric_limits<measurement_id_type>::max();
    sh_threads[threadIndex] = std::numeric_limits<unsigned int>::max();

    if (threadIndex < *(payload.n_meas_to_remove)) {
        sh_meas_ids[threadIndex] = meas_to_remove[threadIndex];
        sh_threads[threadIndex] = threads[threadIndex];
    }

    __syncthreads();

    // Find starting point
    if (threadIndex < *(payload.n_meas_to_remove)) {
        auto mid = sh_meas_ids[threadIndex];
        bool is_start =
            (threadIndex == 0) || (sh_meas_ids[threadIndex - 1] != mid);
        const auto unique_meas_idx = meas_id_to_unique_id.at(mid);
        const auto its_accepted_tracks =
            n_accepted_tracks_per_measurement.at(unique_meas_idx);

        if (is_start) {

            int i = threadIndex + 1;
            int n_sharing_tracks = 1;

            while (i < *(payload.n_meas_to_remove) && sh_meas_ids[i] == mid) {
                if (sh_threads[i] != sh_threads[i - 1]) {
                    n_sharing_tracks++;

                    if (n_sharing_tracks == its_accepted_tracks) {
                        atomicMin(&min_thread, sh_threads[i - 1]);
                        break;
                    }
                }
                i++;
            }
        }
    }

    __syncthreads();

    if (threadIndex == 0) {
        if (min_thread == 0) {
            *(payload.n_removable_tracks) = 1;
        } else if (min_thread == std::numeric_limits<unsigned int>::max()) {
            *(payload.n_removable_tracks) = *(payload.n_tracks_to_iterate);
        } else {
            *(payload.n_removable_tracks) = min_thread;
        }
    }

    __syncthreads();

    meas_to_remove[threadIndex] = sh_meas_ids[threadIndex];
    threads[threadIndex] = sh_threads[threadIndex];

    __syncthreads();

    int is_valid =
        (threads[threadIndex] < *(payload.n_removable_tracks)) ? 1 : 0;

    // TODO: Use better reduction algorithm
    if (is_valid) {
        atomicAdd(payload.n_valid_threads, 1);
    }

    __syncthreads();

    // Exclusive scan (Hillis-Steele)
    prefix[threadIndex] = is_valid;  // copy input
    __syncthreads();

    for (int offset = 1; offset < *(payload.n_meas_to_remove); offset <<= 1) {
        int val = 0;
        if (threadIndex >= offset) {
            val = prefix[threadIndex - offset];
        }
        __syncthreads();
        prefix[threadIndex] += val;
        __syncthreads();
    }

    if (is_valid) {
        prefix[threadIndex] -= 1;
        sh_meas_ids[prefix[threadIndex]] = meas_to_remove[threadIndex];
        sh_threads[prefix[threadIndex]] = threads[threadIndex];
    }

    __syncthreads();

    meas_to_remove[threadIndex] = sh_meas_ids[threadIndex];
    threads[threadIndex] = sh_threads[threadIndex];
}

}  // namespace traccc::cuda::kernels
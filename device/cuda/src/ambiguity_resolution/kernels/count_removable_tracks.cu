/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "count_removable_tracks.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

// Thrust include(s).
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>

namespace traccc::cuda::kernels {

__device__ void count_tracks(int tid, int* sh_n_meas, int n_tracks,
                             unsigned int& bound, unsigned int& count,
                             bool& stop) {

    unsigned int add = 0;
    unsigned int offset = 0;
    for (unsigned int stride = 1; stride < (n_tracks - count); stride *= 2) {
        if ((count + tid + stride) < n_tracks) {
            sh_n_meas[count + tid] += sh_n_meas[count + tid + stride];
        }
        __syncthreads();

        if (sh_n_meas[count] < bound) {
            if (tid == 0) {
                offset = sh_n_meas[count];
                add = stride * 2;
            }
        }

        __syncthreads();
    }

    if (tid == 0) {
        bound -= offset;
        count += add;

        if (add == 0) {
            stop = true;
        }
    }

    __syncthreads();
}

__launch_bounds__(512) __global__ void count_removable_tracks(
    device::count_removable_tracks_payload payload) {

    if (threadIdx.x == 0) {
        if (*(payload.max_shared) == 0) {
            *(payload.terminate) = 1;
        }
    }

    __syncthreads();

    if (*(payload.terminate) == 1) {
        return;
    }

    __shared__ int shared_n_meas[512];
    __shared__ unsigned int n_meas_total;
    __shared__ unsigned int bound;
    __shared__ unsigned int n_tracks_to_iterate;
    __shared__ bool stop;

    vecmem::device_vector<const unsigned int> sorted_ids(
        payload.sorted_ids_view);
    vecmem::jagged_device_vector<const measurement_id_type> meas_ids(
        payload.meas_ids_view);
    vecmem::device_vector<const unsigned int> n_meas(payload.n_meas_view);
    vecmem::device_vector<measurement_id_type> meas_to_remove(
        payload.meas_to_remove_view);
    vecmem::device_vector<unsigned int> threads(payload.threads_view);
    vecmem::device_vector<const unsigned int> n_accepted_tracks_per_measurement(
        payload.n_accepted_tracks_per_measurement_view);
    vecmem::device_vector<const unsigned int> meas_id_to_unique_id(
        payload.meas_id_to_unique_id_view);

    auto threadIndex = threadIdx.x;

    int gid = static_cast<int>(*payload.n_accepted) - 1 - threadIndex;
    shared_n_meas[threadIndex] = 0;

    if (threadIndex == 0) {
        *(payload.n_removable_tracks) = 0;
        *(payload.n_meas_to_remove) = 0;
        *(payload.n_valid_threads) = 0;
        *(payload.n_tracks_to_iterate) = 0;
        n_meas_total = 0;
        bound = 512;
        n_tracks_to_iterate = 0;
        stop = false;
    }

    __syncthreads();

    unsigned int trk_id = 0;
    unsigned int n_m = 0;
    if (gid >= 0) {
        trk_id = sorted_ids[gid];
        n_m = n_meas[trk_id];
        shared_n_meas[threadIndex] = n_m;
    }

    __syncthreads();

    auto n_tracks_total = min(bound, *payload.n_accepted);

    // @TODO: Improve the logic
    count_tracks(threadIdx.x, shared_n_meas, n_tracks_total, bound,
                 n_tracks_to_iterate, stop);
    /*
    for (int i = 0; i < 100; i++) {
        count_tracks(threadIdx.x, shared_n_meas, n_tracks_total, bound,
                     n_tracks_to_iterate, stop);
        __syncthreads();
        if (stop)
            break;

        if (gid >= 0 && static_cast<unsigned int>(gid) < sorted_ids.size()) {
            const auto trk_id = sorted_ids[gid];
            if (trk_id < n_meas.size()) {
                shared_n_meas[threadIndex] = n_meas[trk_id];
            }
        }
        __syncthreads();
    }
    */

    if (threadIndex == 0 && n_tracks_to_iterate == 0) {
        n_tracks_to_iterate = 1;
    }

    // @TODO: Improve the logic
    if (threadIndex < n_tracks_to_iterate && gid >= 0) {
        const unsigned int pos = atomicAdd(&n_meas_total, n_m);

        const auto& mids = meas_ids[trk_id];
        for (int i = 0; i < n_m; i++) {
            meas_to_remove[pos + i] = mids[i];
            threads[pos + i] = threadIndex;
        }
    }

    __syncthreads();

    if (threadIndex == 0) {
        *(payload.n_meas_to_remove) = n_meas_total;
        *(payload.n_tracks_to_iterate) = n_tracks_to_iterate;
    }
}

}  // namespace traccc::cuda::kernels

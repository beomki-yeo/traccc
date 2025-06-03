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
            __syncthreads();
        }
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

__device__ void bitonic_sort_shared(
    const int tid, traccc::pair<std::size_t, unsigned int>* shared_data,
    const int count, const int N) {

    if (tid >= count && tid < N) {
        shared_data[tid] = {std::numeric_limits<std::size_t>::max(),
                            std::numeric_limits<unsigned int>::max()};
    }

    __syncthreads();

    for (int k = 2; k <= N; k <<= 1) {
        for (int j = k >> 1; j > 0; j >>= 1) {
            int ixj = tid ^ j;

            if (ixj > tid && ixj < N && tid < N) {
                auto elem_i = shared_data[tid];
                auto elem_j = shared_data[ixj];

                bool ascending = ((tid & k) == 0);
                bool should_swap =
                    (elem_i.first > elem_j.first ||
                     (elem_i.first == elem_j.first &&
                      elem_i.second > elem_j.second)) == ascending;

                if (should_swap) {
                    shared_data[tid] = elem_j;
                    shared_data[ixj] = elem_i;
                }
            }
            __syncthreads();
        }
    }
}

__global__ void count_removable_tracks(
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

    __shared__ int shared_n_meas[1024];
    __shared__ traccc::pair<std::size_t, unsigned int> meas_to_thread[1024];
    __shared__ unsigned int n_meas_total;
    __shared__ unsigned int bound;
    __shared__ unsigned int n_tracks_to_iterate;
    __shared__ unsigned int min_thread;
    __shared__ unsigned int N;
    __shared__ bool stop;

    vecmem::device_vector<const unsigned int> sorted_ids(
        payload.sorted_ids_view);
    vecmem::jagged_device_vector<const std::size_t> meas_ids(
        payload.meas_ids_view);
    vecmem::device_vector<const std::size_t> n_meas(payload.n_meas_view);
    vecmem::device_vector<traccc::pair<std::size_t, unsigned int>>
        meas_to_remove(payload.meas_to_remove_view);
    vecmem::device_vector<const std::size_t> unique_meas(
        payload.unique_meas_view);
    vecmem::device_vector<const unsigned int> n_accepted_tracks_per_measurement(
        payload.n_accepted_tracks_per_measurement_view);

    auto threadIndex = threadIdx.x;

    int gid = static_cast<int>(*payload.n_accepted) - 1 - threadIndex;
    shared_n_meas[threadIndex] = 0;
    meas_to_thread[threadIndex] = {std::numeric_limits<std::size_t>::max(),
                                   std::numeric_limits<unsigned int>::max()};

    if (threadIndex == 0) {
        *(payload.n_removable_tracks) = 0;
        *(payload.n_meas_to_remove) = 0;
        n_meas_total = 0;
        bound = 1024;
        N = 1;
        n_tracks_to_iterate = 0;
        min_thread = std::numeric_limits<unsigned int>::max();
        stop = false;
    }

    __syncthreads();

    if (gid >= 0) {
        shared_n_meas[threadIndex] = n_meas[sorted_ids[gid]];
    }

    __syncthreads();

    auto n_tracks_total = min(bound, *payload.n_accepted);

    // @TODO: Improve the logic
    count_tracks(threadIdx.x, shared_n_meas, n_tracks_total, bound,
                 n_tracks_to_iterate, stop);

    if (threadIndex == 0 && n_tracks_to_iterate == 0) {
        n_tracks_to_iterate = 1;
    }

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

    // @TODO: Improve the logic
    if (threadIndex < n_tracks_to_iterate && gid >= 0) {
        const auto& mids = meas_ids[sorted_ids[gid]];
        for (const auto& id : mids) {
            const unsigned int pos = atomicAdd(&n_meas_total, 1);
            meas_to_thread[pos] = {id, threadIndex};
        }
    }

    __syncthreads();

    // Bitonic sort on meas_to_thread w.r.t. measurement id
    if (threadIndex == 0) {
        N = (n_meas_total == 0) ? 1 : 1 << (32 - __clz(n_meas_total - 1));
    }
    __syncthreads();

    bitonic_sort_shared(threadIndex, meas_to_thread, n_meas_total, N);
    /*
    if (threadIndex == 0) {
        for (const auto& e: unique_meas){
            printf("%lu ", e);
        }
        printf("\n");
    }
    */
    // Find starting point
    if (threadIndex < n_meas_total) {
        auto curr = meas_to_thread[threadIndex];
        bool is_start = (threadIndex == 0) ||
                        (meas_to_thread[threadIndex - 1].first != curr.first);

        if (is_start) {

            int i = threadIndex + 1;
            int n_sharing_tracks = 1;

            const std::size_t unique_meas_idx =
                thrust::lower_bound(thrust::seq, unique_meas.begin(),
                                    unique_meas.end(), curr.first) -
                unique_meas.begin();

            while (i < n_meas_total && meas_to_thread[i].first == curr.first) {
                if (meas_to_thread[i].second != meas_to_thread[i - 1].second) {
                    n_sharing_tracks++;

                    /*
                    printf("%d %d %d \n", threadIndex, n_sharing_tracks,
                           n_accepted_tracks_per_measurement.at(
                               static_cast<unsigned int>(unique_meas_idx)));
                    */

                    /*
                    printf(
                        "thread index %d n sharing %d unique meas idx %lu curr "
                        "first %lu\n",
                        threadIndex, n_sharing_tracks, unique_meas_idx,
                        curr.first);
                    */

                    // atomicMin(&min_thread, meas_to_thread[i].second);

                    if (n_sharing_tracks ==
                        n_accepted_tracks_per_measurement.at(unique_meas_idx)) {
                        atomicMin(&min_thread, meas_to_thread[i - 1].second);
                        break;
                    }
                }
                /*
                if (meas_to_thread[i].second != curr.second) {
                    atomicMin(&min_thread, meas_to_thread[i].second);
                }
                */
                i++;
            }

            /*
            if (n_sharing_tracks >= 2 &&
                (n_sharing_tracks ==
                 n_accepted_tracks_per_measurement.at(unique_meas_idx))) {
                atomicMin(&min_thread, curr.second);
            }
            */
        }
    }

    __syncthreads();

    if (threadIndex == 0) {
        n_meas_total = 0;

        if (min_thread == 0) {
            *(payload.n_removable_tracks) = 1;
        } else if (min_thread == std::numeric_limits<unsigned int>::max()) {
            *(payload.n_removable_tracks) = n_tracks_to_iterate;
        } else {
            *(payload.n_removable_tracks) = min_thread;
        }
    }

    __syncthreads();

    // Count n_meas_total again
    if (meas_to_thread[threadIndex].second < *(payload.n_removable_tracks)) {
        const unsigned int pos = atomicAdd(&n_meas_total, 1);
        meas_to_remove[pos] = meas_to_thread[threadIndex];
    }

    __syncthreads();

    if (threadIndex == 0) {
        *(payload.n_meas_to_remove) = n_meas_total;
    }

    if (threadIndex == 0) {
        printf(
            "min thread %d removable tracks %d max shared %d n meas to remove "
            "%d\n",
            min_thread, *(payload.n_removable_tracks), *(payload.max_shared),
            *(payload.n_meas_to_remove));

        for (int i = 0; i < *(payload.n_meas_to_remove); i++) {
            printf("(%lu %d) ", meas_to_remove[i].first,
                   meas_to_remove[i].second);
        }
        printf("\n");

        printf("n accepted track per meas \n");
        for (int i = 0; i < unique_meas.size(); i++) {
            printf("(%lu %d) ", unique_meas.at(i),
                   n_accepted_tracks_per_measurement.at(i));
        }
        printf("\n");
    }
}

}  // namespace traccc::cuda::kernels

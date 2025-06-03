/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/utils/pair.hpp"

// VecMem include(s).
#include <vecmem/containers/data/jagged_vector_view.hpp>
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

inline __device__ void bitonic_sort_shared(
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

/// (Event Data) Payload for the @c
/// traccc::device::count_removable_tracks function
struct count_removable_tracks_payload {

    /**
     * @brief Whether to terminate the calculation
     */
    int* terminate;

    /**
     * @brief The number of max shared
     */
    unsigned int* max_shared;

    /**
     * @brief View object to the sorted track
     */
    vecmem::data::vector_view<const unsigned int> sorted_ids_view;

    /**
     * @brief The number of accepted tracks
     */
    unsigned int* n_accepted;

    /**
     * @brief View object to the vector of measured ids per track
     */
    vecmem::data::jagged_vector_view<const std::size_t> meas_ids_view;

    /**
     * @brief View object to the vector of number of measurements
     */
    vecmem::data::vector_view<const std::size_t> n_meas_view;

    /**
     * @brief View object to the unique measurement ids
     */
    vecmem::data::vector_view<const std::size_t> unique_meas_view;

    /**
     * @brief View object to the number of accepted tracks per measurement
     */
    vecmem::data::vector_view<const unsigned int>
        n_accepted_tracks_per_measurement_view;

    /**
     * @brief The number of worst tracks removable
     */
    unsigned int* n_removable_tracks;

    /**
     * @brief The number of measurements to remove
     */
    unsigned int* n_meas_to_remove;

    /**
     * @brief View object to measurements to remove
     */
    vecmem::data::vector_view<traccc::pair<std::size_t, unsigned int>>
        meas_to_remove_view;
};

}  // namespace traccc::device

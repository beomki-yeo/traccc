/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Event Data) Payload for the @c
/// traccc::device::prune_measurements_to_remove function
struct prune_measurements_to_remove_payload {

    /**
     * @brief Whether to terminate the calculation
     */
    int* terminate;

    /**
     * @brief The number of worst tracks removable
     */
    unsigned int* n_removable_tracks;

    /**
     * @brief The number of measurements to remove
     */
    unsigned int* n_meas_to_remove;

    /**
     * @brief The number of threads that can remove its corresponding track
     */
    unsigned int* n_valid_threads;

    /**
     * @brief The number of tracks to iterate
     */
    unsigned int* n_tracks_to_iterate;

    /**
     * @brief View object to the meas id to unique id map
     */
    vecmem::data::vector_view<const unsigned int> meas_id_to_unique_id_view;

    /**
     * @brief View object to the number of accepted tracks per measurement
     */
    vecmem::data::vector_view<const unsigned int>
        n_accepted_tracks_per_measurement_view;

    /**
     * @brief View object to measurements to remove
     */
    vecmem::data::vector_view<measurement_id_type> meas_to_remove_view;

    /**
     * @brief View object to thread id of measurements to remove
     */
    vecmem::data::vector_view<unsigned int> threads_view;
};

}  // namespace traccc::device

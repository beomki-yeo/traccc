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
/// traccc::device::block_bitonic_sort function
struct block_bitonic_sort_payload {

    /**
     * @brief Whether to terminate the calculation
     */
    int* terminate;

    /**
     * @brief The number of measurements to remove
     */
    unsigned int* n_meas_to_remove;

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

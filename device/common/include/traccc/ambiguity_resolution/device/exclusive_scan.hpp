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
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Event Data) Payload for the @c
/// traccc::device::exclusive_scan function
struct exclusive_scan_payload {

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
    vecmem::data::vector_view<traccc::pair<std::size_t, unsigned int>>
        meas_to_remove_view;
};

}  // namespace traccc::device

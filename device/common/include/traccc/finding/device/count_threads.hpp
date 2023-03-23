/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

namespace traccc::device {

template <typename config_t>
TRACCC_DEVICE inline void count_threads(
    std::size_t globalIndex, const config_t cfg,
    vecmem::data::vector_view<unsigned int> n_measurements_view,
    vecmem::data::vector_view<unsigned int> n_threads_view,
    const unsigned int& n_total_measurements,
    unsigned int& n_measurements_per_thread, unsigned int& n_total_threads);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/count_threads.ipp"
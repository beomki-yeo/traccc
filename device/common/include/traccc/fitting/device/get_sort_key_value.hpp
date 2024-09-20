/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/device/sort_key.hpp"
#include "traccc/edm/track_candidate.hpp"

// System include(s).
#include <cstddef>

namespace traccc::device {

/// Function used for fitting a track for a given track candidates
///
/// @param[in] globalIndex   The index of the current thread
/// @param[out] keys_view    The key values
/// @param[out] ids_view     The param ids
///
TRACCC_HOST_DEVICE inline void get_sort_key_value(
    std::size_t globalIndex,
    track_candidate_container_types::const_view track_candidates_view,
    vecmem::data::vector_view<device::sort_key> keys_view,
    vecmem::data::vector_view<unsigned int> ids_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/fitting/device/impl/get_sort_key_value.ipp"

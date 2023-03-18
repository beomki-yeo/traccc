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

TRACCC_DEVICE inline void build_tracks(
    std::size_t globalIndex,
    measurement_container_types::const_view measurements_view,
    bound_track_parameters_collection_types::view seeds_view,
    vecmem::data::jagged_vector_view<candidate_link> link_view,
    vecmem::data::vector_view<thrust::pair<unsigned int, unsigned int>>
        tips_view,
    track_candidate_container_types::view track_candidates_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/build_tracks.ipp"
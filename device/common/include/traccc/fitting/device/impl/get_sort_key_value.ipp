/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::device {

TRACCC_HOST_DEVICE inline void get_sort_key_value(
    std::size_t globalIndex,
    track_candidate_container_types::const_view track_candidates_view,
    vecmem::data::vector_view<device::sort_key> keys_view,
    vecmem::data::vector_view<unsigned int> ids_view) {

    track_candidate_container_types::const_device track_candidates(
        track_candidates_view);

    if (globalIndex >= track_candidates.size()) {
        return;
    }

    // Keys
    vecmem::device_vector<traccc::scalar> keys_device(keys_view);

    // Param id
    vecmem::device_vector<unsigned int> ids_device(ids_view);

    keys_device.at(globalIndex) =
        device::get_sort_key(track_candidates.at(globalIndex).header);
    ids_device.at(globalIndex) = globalIndex;
}

}  // namespace traccc::device
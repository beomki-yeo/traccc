/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::device {

TRACCC_HOST_DEVICE inline void fill_sort_keys(
    std::size_t globalIndex,
    const track_candidate_container_types::const_view& track_candidates_view,
    vecmem::data::vector_view<device::sort_key> keys_view,
    vecmem::data::vector_view<unsigned int> ids_view) {

    track_candidate_container_types::const_device track_candidates(
        track_candidates_view);

    // Keys
    vecmem::device_vector<device::sort_key> keys_device(keys_view);

    // Param id
    vecmem::device_vector<unsigned int> ids_device(ids_view);

    if (globalIndex >= keys_device.size()) {
        return;
    }

    // Key = The number of measurements
    keys_device.at(globalIndex) =
        device::get_sort_key(track_candidates.at(globalIndex).header);
    ids_device.at(globalIndex) = globalIndex;
}

}  // namespace traccc::device
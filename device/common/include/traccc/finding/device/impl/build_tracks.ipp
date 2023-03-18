/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::device {

TRACCC_DEVICE inline void build_tracks(
    std::size_t globalIndex,
    measurement_container_types::const_view measurements_view,
    bound_track_parameters_collection_types::view seeds_view,
    vecmem::data::jagged_vector_view<candidate_link> link_view,
    vecmem::data::vector_view<thrust::pair<unsigned int, unsigned int>>
        tips_view,
    track_candidate_container_types::view track_candidates_view) {

    measurement_container_types::const_device measurements(measurements_view);

    bound_track_parameters_collection_types::device seeds(seeds_view);

    vecmem::jagged_device_vector<candidate_link> links(link_view);

    vecmem::device_vector<thrust::pair<unsigned int, unsigned int>> tips(
        tips_view);

    track_candidate_container_types::device track_candidates(
        track_candidates_view);

    if (globalIndex >= tips.size()) {
        return;
    }

    const auto tip = tips.at(globalIndex);
    auto& seed = track_candidates[globalIndex].header;
    auto cands_per_track = track_candidates[globalIndex].items;

    auto L = links[tip.first][tip.second];

    for (auto it = cands_per_track.rbegin(); it != cands_per_track.rend();
         it++) {

        auto& cand = *it;
        cand = {L.surface_link, measurements.at(L.meas_link)};

        if (it == cands_per_track.rend() - 1) {
            seed = seeds.at(L.last_link.second);
            break;
        }

        L = links[L.last_link.first][L.last_link.second];
    }
}

}  // namespace traccc::device
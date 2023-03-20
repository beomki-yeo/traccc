/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/track_parameters.hpp"

namespace traccc::device {

template <typename propagator_t, typename config_t>
TRACCC_DEVICE inline void find_tracks(
    std::size_t globalIndex, const config_t cfg,
    typename propagator_t::detector_type::detector_view_type det_data,
    vecmem::data::jagged_vector_view<typename propagator_t::intersection_type>
        nav_candidates_buffer,
    measurement_container_types::const_view measurements_view,
    vecmem::data::vector_view<thrust::pair<unsigned int, unsigned int>>
        module_map_view,
    bound_track_parameters_collection_types::view in_params_view,
    bound_track_parameters_collection_types::view out_params_view,
    vecmem::data::vector_view<candidate_link> links_view,
    vecmem::data::vector_view<unsigned int> param_to_link_view,
    vecmem::data::vector_view<thrust::pair<unsigned int, unsigned int>>
        tips_view,
    vecmem::data::vector_view<unsigned int> n_threads_view,
    const unsigned int& iteration,
    const unsigned int& n_measurements_per_thread,
    const unsigned int& n_total_threads, unsigned int& n_measurements_per_step,
    unsigned int& n_out_params);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/find_tracks.ipp"
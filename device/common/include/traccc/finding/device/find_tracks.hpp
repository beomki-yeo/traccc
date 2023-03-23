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

/// Function for combinatorial finding:
/// If the chi2 of the measurement < chi2_max, its measurement index and the
/// index of the link from the previous step are added to the link container.
/// Then the tracks are spawned for the measurements satisfying the condition
/// and propagated to the next surface. If a track finds a surface that contains
/// measurements, its bound track parameter on the surface will be used for the
/// next step. Otherwise, the link is added into the tip link container so that
/// we can know which links in the link container are the final measurements of
/// full tracks
///
/// @param[in] globalIndex       The index of the current thread
/// @param[in] cfg               Track finding config object
/// @param[in] det_data          Detector view object
/// @param[in] nav_candidates_buffer Navgation buffer
/// @param[in] measurements_view Measurements container view
/// @param[in] module_map_view   Module map view
/// @param[in] in_params_view    Input parameters
/// @param[in] n_threads_view    The number of threads per tracks
/// @param[in] step              Step index
/// @param[in] n_measurements_per_thread  Number of measurements per thread
/// @param[in] n_total_threads   Number of total threads
/// @param[out] out_params_view  Output parameters
/// @param[out] links_view       link container for the current step
/// @param[out] param_to_link_view  Container for param index -> link index
/// @param[out] tips_view        Tip link container for the current step
///
template <typename propagator_t, typename config_t>
TRACCC_DEVICE inline void find_tracks(
    std::size_t globalIndex, const config_t cfg,
    typename propagator_t::detector_type::detector_view_type det_data,
    vecmem::data::jagged_vector_view<typename propagator_t::intersection_type>
        nav_candidates_buffer,
    measurement_container_types::const_view measurements_view,
    vecmem::data::vector_view<const thrust::pair<unsigned int, unsigned int>>
        module_map_view,
    bound_track_parameters_collection_types::const_view in_params_view,
    vecmem::data::vector_view<const unsigned int> n_threads_view,
    const unsigned int step, const unsigned int& n_measurements_per_thread,
    const unsigned int& n_total_threads,
    bound_track_parameters_collection_types::view out_params_view,
    vecmem::data::vector_view<candidate_link> links_view,
    vecmem::data::vector_view<unsigned int> param_to_link_view,
    vecmem::data::vector_view<thrust::pair<unsigned int, unsigned int>>
        tips_view,
    unsigned int& n_measurements_per_step, unsigned int& n_out_params);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/find_tracks.ipp"
/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/fitting/kalman_filter/gain_matrix_updater.hpp"

// System include(s).
#include <limits>

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
    vecmem::data::vector_view<thrust::pair<unsigned int, unsigned int>>
        tips_view,
    vecmem::data::vector_view<unsigned int> n_threads_view,
    const unsigned int& iteration,
    const unsigned int& n_measurements_per_thread,
    const unsigned int& n_total_threads, unsigned int& n_candidates,
    unsigned int& n_out_params) {

    if (globalIndex >= n_total_threads) {
        return;
    }

    // Detector
    typename propagator_t::detector_type det(det_data);

    // Navigation candidate buffer
    vecmem::jagged_device_vector<typename propagator_t::intersection_type>
        nav_candidates(nav_candidates_buffer);

    // Measurement
    measurement_container_types::const_device measurements(measurements_view);

    // module map
    vecmem::device_vector<thrust::pair<unsigned int, unsigned int>> module_map(
        module_map_view);

    // Input parameters
    bound_track_parameters_collection_types::device in_params(in_params_view);

    // Output parameters
    bound_track_parameters_collection_types::device out_params(out_params_view);

    // links
    vecmem::device_vector<candidate_link> links(links_view);

    // tips
    vecmem::device_vector<thrust::pair<unsigned int, unsigned int>> tips(
        tips_view);

    // n threads
    vecmem::device_vector<unsigned int> n_threads(n_threads_view);

    // Search for out_param index
    const auto lo1 = thrust::lower_bound(thrust::seq, n_threads.begin(),
                                         n_threads.end(), globalIndex + 1);
    const auto in_param_id = std::distance(n_threads.begin(), lo1);

    // Get module id
    const auto module_id = in_params.at(in_param_id).surface_link();

    // Search for measuremetns header ID
    const auto lo2 = thrust::lower_bound(
        thrust::seq, module_map.begin(), module_map.end(), module_id,
        compare_pair_int<thrust::pair, unsigned int>());
    const auto header_id = (*lo2).second;

    // Get measurements on surface
    const auto measurements_on_surface = measurements.at(header_id).items;

    unsigned int ref;
    if (lo1 == n_threads.begin()) {
        ref = 0;
    } else {
        ref = *(lo1 - 1);
    }

    const unsigned int offset = globalIndex - ref;
    const unsigned int stride = offset * n_measurements_per_thread;
    const unsigned int n_meas_on_surface = measurements_on_surface.size();

    // Iterate over the measurements
    const auto& mask_store = det.mask_store();
    const auto& surface = det.surfaces(module_id);

    // Create propagator
    propagator_t propagator({}, {});

    // Last iteration ID
    const unsigned int last_iteration =
        (iteration == 0) ? std::numeric_limits<unsigned int>::max()
                         : iteration - 1;

    for (unsigned int i = 0; i < n_measurements_per_thread; i++) {
        if (i + stride >= n_meas_on_surface) {
            break;
        }

        bound_track_parameters in_par = in_params.at(in_param_id);
        const auto meas = measurements_on_surface.at(i + stride);
        track_state<typename propagator_t::transform3_type> trk_state(
            {module_id, meas});

        mask_store.template visit<
            gain_matrix_updater<typename propagator_t::transform3_type>>(
            surface.mask(), trk_state, in_par);

        const auto chi2 = trk_state.filtered_chi2();

        if (chi2 < cfg.chi2_max) {

            // Create propagator state
            typename propagator_t::state propagation(
                in_par, det.get_bfield(), det,
                std::move(nav_candidates.at(globalIndex)));

            // Actor state
            // @TODO: simplify the syntax here
            using actor_list_type =
                typename propagator_t::actor_chain_type::actor_list_type;
            typename detray::detail::tuple_element<
                0, actor_list_type>::type::state s0{};
            typename detray::detail::tuple_element<
                1, actor_list_type>::type::state s1{};
            typename detray::detail::tuple_element<
                2, actor_list_type>::type::state s2{
                cfg.min_step_length_for_surface_aborter};

            // Run propagation
            // @TODO test with propagate_sync()
            propagator.propagate(propagation, std::tie(s0, s1, s2));

            vecmem::device_atomic_ref<unsigned int> num_measurements_per_layer(
                n_candidates);
            const unsigned int l_pos = num_measurements_per_layer.fetch_add(1);

            // @TODO; Consider max_num_branches_per_surface

            links[l_pos] = {{last_iteration, in_param_id},
                            {header_id, i + stride},
                            module_id};

            if (s2.success) {
                vecmem::device_atomic_ref<unsigned int> num_out_params(
                    n_out_params);
                const unsigned int p_pos = num_out_params.fetch_add(1);

                out_params[p_pos] = propagation._stepping._bound_params;
            }
            // Unless the track found a surface, it is considered a tip
            else if (!s2.success &&
                     iteration >= cfg.min_track_candidates_per_track) {
                tips.push_back({iteration, l_pos});
            }
        }
    }
}

}  // namespace traccc::device
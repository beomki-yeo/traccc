/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::device {

template <typename propagator_t>
TRACCC_DEVICE inline void apply_interaction(
    std::size_t globalIndex,
    typename propagator_t::detector_type::detector_view_type det_data,
    vecmem::data::jagged_vector_view<typename propagator_t::intersection_type>
        nav_candidates_buffer,
    const int n_params,
    bound_track_parameters_collection_types::view params_view) {

    // Detector
    typename propagator_t::detector_type det(det_data);

    // Navigation candidate buffer
    vecmem::jagged_device_vector<typename propagator_t::intersection_type>
        nav_candidates(nav_candidates_buffer);

    // in param
    bound_track_parameters_collection_types::device params(params_view);

    if (globalIndex >= n_params) {
        return;
    }

    // Create propagator
    propagator_t propagator({}, {});

    // Create propagator state
    typename propagator_t::state propagation(
        params.at(globalIndex), det.get_bfield(), det,
        std::move(nav_candidates.at(globalIndex)));

    // Actor state
    // @TODO: simplify the syntax here
    using actor_list_type =
        typename propagator_t::actor_chain_type::actor_list_type;
    typename detray::detail::tuple_element<0, actor_list_type>::type::state
        s0{};

    // Run propagation
    // @TODO test with propagate_sync()
    propagator.propagate(propagation, std::tie(s0));

    // Replace the parameter
    params.at(globalIndex) = propagation._stepping._bound_params;
}

}  // namespace traccc::device

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc {

/// Track Finding algorithm for a set of tracks
template <typename stepper_t, typename navigator_t>
class finding_algorithm
    : public algorithm<track_candidate_container_types::buffer(
          const typename navigator_t::detector_type::detector_view_type&,
          const typename measurement_container_types::const_view&,
          bound_track_parameters_collection_types::buffer&&)> {

    // Transform3 type
    using transform3_type = typename stepper_t::transform3_type;

    // Detector type
    using detector_type = typename navigator_t::detector_type;

    // Actor types
    using transporter = detray::parameter_transporter<transform3_type>;
    using interactor = detray::pointwise_material_interactor<transform3_type>;

    // scalar type
    using scalar_type = typename transform3_type::scalar_type;

    // Actor chain for material interactor and its propagator type
    using actor_for_interaction = detray::actor_chain<std::tuple, interactor>;

    using propagator_for_interaction =
        detray::propagator<stepper_t, navigator_t, actor_for_interaction>;

    // Actor chain for propagate to the next surface and its propagator type
    using actor_for_propagation =
        detray::actor_chain<std::tuple, detray::pathlimit_aborter, transporter,
                            detray::next_surface_aborter>;

    using propagator_for_propagation =
        detray::propagator<stepper_t, navigator_t, actor_for_propagation>;

    public:
    /// Constructor for the finding algorithm
    ///
    /// @param mr The memory resource to use
    finding_algorithm(const traccc::memory_resource& mr);
};

}  // namespace traccc
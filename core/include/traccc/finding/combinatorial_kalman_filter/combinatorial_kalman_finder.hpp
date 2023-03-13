/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/finding/combinatorial_kalman_filter/ckf_actor.hpp"

// detray include(s).
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/propagator.hpp"

namespace traccc {

/// Combinatorial Kalman finder algorithm to find tracks from measurement
/// container and seed parameters
template <typename stepper_t, typename navigator_t>
class combinatorial_kalman_finder {

    public:
    // scalar type
    using scalar_type = typename stepper_t::scalar_type;

    // vector type
    template <typename T>
    using vector_type = typename navigator_t::template vector_type<T>;

    // Combinatorial Kalman finder configuration
    struct config {
        scalar_type pathlimit = std::numeric_limits<scalar>::max();
        scalar_type overstep_tolerance = -10 * detray::unit<scalar>::um;
        scalar_type step_constraint = 5. * detray::unit<scalar>::mm;

        // maximum number of branches per seed
        // std::size_t max_num_branches_per_seed = 128;

        // max chi2 for selection
        scalar_type max_chi2 = 15.;

        // maximum number of branches per module
        std::size_t max_num_branches_per_module = 10;
    };

    // navigator candidate type
    using intersection_type = typename navigator_t::intersection_type;

    using source_link_type =
        typename measurement_container_types::host::link_type;

    // transform3 type
    using transform3_type = typename stepper_t::transform3_type;

    // Detector type
    using detector_type = typename navigator_t::detector_type;

    // Actor types
    using aborter = detray::pathlimit_aborter;
    using transporter = detray::parameter_transporter<transform3_type>;
    using interactor = detray::pointwise_material_interactor<transform3_type>;
    using find_actor = traccc::ckf_actor<transform3_type, vector_type>;
    using resetter = detray::parameter_resetter<transform3_type>;

    using actor_chain_type =
        detray::actor_chain<std::tuple, aborter, transporter, interactor,
                            fit_actor, resetter>;

    // Propagator type
    using propagator_type =
        detray::propagator<stepper_t, navigator_t, actor_chain_type>;

    /// Constructor with a detector
    ///
    /// @param det the detector object
    TRACCC_HOST_DEVICE
    combinatorial_kalman_finder(const detector_type& det) : m_detector(det) {}

    /// Combinatorial Kalman finder state
    struct state {

        /// State constructor
        state() = default;

        /// @return the actor chain state
        TRACCC_HOST_DEVICE
        typename actor_chain_type::state operator()() {
            return std::tie(m_aborter_state, m_transporter_state,
                            m_interactor_state, m_fit_actor_state,
                            m_resetter_state);
        }

        /// Individual actor states
        typename aborter::state m_aborter_state{};
        typename transporter::state m_transporter_state{};
        typename interactor::state m_interactor_state{};
        typename find_actor::state m_find_actor_state{};
        typename resetter::state m_resetter_state{};
        vector_type<std::pair<source_link_type, source_link_type>> m_links{};
    };

    TRACCC_HOST_DEVICE void find() {}

    TRACCC_HOST_DEVICE void find_per_module(
        const bound_track_parameters& seed_params, state& finder_state,
        vector_type<intersection_type>&& nav_candidates = {}) {

        // Create propagator
        propagator_type propagator({}, {});

        // Set path limit
        finder_state.m_aborter_state.set_path_limit(m_cfg.pathlimit);

        // Create propagator state
        typename propagator_type::state propagation(
            seed_params, m_detector.get_bfield(), m_detector,
            std::move(nav_candidates));

        // Set overstep tolerance and stepper constraint
        propagation._stepping().set_overstep_tolerance(
            m_cfg.overstep_tolerance);
        propagation._stepping
            .template set_constraint<detray::step::constraint::e_accuracy>(
                m_cfg.step_constraint);

        // Run forward filtering
        propagator.propagate(propagation, finder_state());

        return;
    }

    private:
    // Detector object
    const detector_type& m_detector;
    // Configuration object
    config m_cfg;
    // Source links
    vector_type<source_link_type> m_links;
};

}  // namespace traccc 
/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/track_candidate.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/edm/track_state.hpp"
#include "traccc/fitting/kalman_filter/kalman_actor.hpp"

// detray include(s).
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/propagator.hpp"

namespace traccc {

/// kalman fitter for CPU
template <typename stepper_t, typename navigator_t>
class kalman_fitter {

    public:
    struct config {
        std::size_t n_iterations = 1;
    };

    // transform3 type
    using transform3_type = typename stepper_t::transform3_type;

    // Detector type
    using detector_type = typename navigator_t::detector_type;

    // Actor types
    using transporter = detray::parameter_transporter<transform3_type>;
    using interactor = detray::pointwise_material_interactor<transform3_type>;
    using fit_actor = traccc::kalman_actor<transform3_type, vecmem::vector>;
    using resetter = detray::parameter_resetter<transform3_type>;

    using actor_chain_type =
        detray::actor_chain<std::tuple, transporter, interactor, fit_actor,
                            resetter>;

    // Propagator type
    using propagator_type =
        detray::propagator<stepper_t, navigator_t, actor_chain_type>;

    kalman_fitter(const detector_type& det)
        : m_detector(std::make_unique<detector_type>(det)) {}

    template <typename seed_parameters_t>
    void run(const seed_parameters_t& seed_params,
             vecmem::vector<track_state<transform3_type>>&& track_states) {

        // Kalman actor state that takes track candidates
        typename fit_actor::state fit_actor_state(std::move(track_states));

        for (std::size_t i = 0; i < m_cfg.n_iterations; i++) {
            fit_actor_state.reset();

            if (i == 0) {
                fit(seed_params, fit_actor_state);
            } else {
                const auto& new_seed_params =
                    fit_actor_state.m_track_states[0].smoothed();

                fit(new_seed_params, fit_actor_state);
            }
        }

        m_track_states = std::move(fit_actor_state.m_track_states);
    }

    template <typename seed_parameters_t>
    void fit(const seed_parameters_t& seed_params,
             typename fit_actor::state& fit_actor_state) {

        // Create propagator
        propagator_type propagator({}, {});

        // Create actor chain states
        typename actor_chain_type::state actor_states =
            std::tie(m_transporter_state, m_interactor_state, fit_actor_state,
                     m_resetter_state);

        // Create propagator state
        typename propagator_type::state propagation(
            seed_params, m_detector->get_bfield(), *m_detector, actor_states);

        // Run forward filtering
        propagator.propagate(propagation);

        // Run smoothing
        smooth(fit_actor_state.m_track_states);

        // Backward Propagation (Currently Doesn't work)
        /*
        // Flip the stepper direction and momentum for backward filtering
        propagation._stepping().flip();

        // Flip the navigation direction for backward filtering
        propagation._navigation.set_direction(
            detray::navigation::direction::e_backward);

        // Run backward filtering
        propagator.propagate(propagation);
        */

        // TODO: Write track info
    }

    // For Eq. (3.36 - 3.37) of "Pattern Recognition, Tracking and Vertex
    // Reconstruction in Particle Detectors" from R.Frühwirth and A. Strandlie
    template <typename track_state_collection_t>
    void smooth(track_state_collection_t& track_states) {

        // Last track state's smoothed parameter is eqaul to filtered one
        auto& last = track_states.back();
        last.smoothed().set_vector(last.filtered().vector());
        last.smoothed().set_covariance(last.filtered().covariance());

        const auto& mask_store = m_detector.get()->mask_store();

        for (typename track_state_collection_t::reverse_iterator it =
                 track_states.rbegin() + 1;
             it != track_states.rend(); ++it) {

            // Surface
            const auto& surface =
                m_detector.get()->surface_by_index(it->surface_link());

            // Run kalman smoother
            mask_store.template call<gain_matrix_smoother<transform3_type>>(
                surface.mask(), surface, *it, *(it - 1));
        }
    }

    /*
    // For Eq. (3.38 - 3.40) of "Pattern Recognition, Tracking and Vertex
    // Reconstruction in Particle Detectors" from R.Frühwirth and A. Strandlie
    template <typename track_state_collection_t>
    void smooth2(track_state_collection_t& track_states) {

        const auto& mask_store = m_detector.get()->mask_store();

        for (auto it = track_states.begin(); it != track_states.end(); it++) {
            auto& trk_state = *it;

            // Surface
            const auto& surface =
                m_detector.get()->surface_by_index(trk_state.surface_link());

            // Run kalman smoother
            mask_store.template call<gain_matrix_smoother<transform3_type>>(
                surface.mask(), surface, trk_state);
        }
    }
    */

    fitter_info<transform3_type> get_fitter_info() const { return m_fit_info; }
    vecmem::vector<track_state<transform3_type>> get_track_states() const {
        return m_track_states;
    }

    private:
    std::unique_ptr<detector_type> m_detector;
    fitter_info<transform3_type> m_fit_info;
    vecmem::vector<track_state<transform3_type>> m_track_states;

    typename transporter::state m_transporter_state{};
    typename interactor::state m_interactor_state{};
    typename resetter::state m_resetter_state{};

    config m_cfg;
};

}  // namespace traccc
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
#include "traccc/edm/track_candidate.hpp"
#include "traccc/edm/track_state.hpp"
#include "traccc/fitting/kalman_filter/gain_matrix_updater.hpp"

// detray include(s).
#include "detray/propagator/base_actor.hpp"

// Thrust include(s).
#include <thrust/pair.h>

// System include(s).
#include <limits>

namespace traccc {

/// Detray actor for Combinatorial Kalman Filtering counter
template <typename algebra_t, template <typename...> class vector_t>
struct ckf_actor : detray::actor {

    using scalar_type = typename algebra_t::scalar_type;

    // Actor state
    struct state {
        /*
        unsigned int max_num_measurements_per_module = 10;
        const measurement_container_types::host& measurements;
        vector_t<bound_track_parameter> seeds;
        vector_t<track_candidate> candidates{};
        vector_t<thrust::pair<scalar_type, scalar_type>> chi2_index_buffer;
        */
    };

    template <typename propagator_state_t>
    TRACCC_HOST_DEVICE void operator()(state& actor_state,
                                       propagator_state_t& propagation) const {

        const auto& navigation = propagation._navigation;

        // triggered only for sensitive surfaces
        if (navigation.is_on_sensitive()) {

            const auto& stepping = propagation._stepping;
            const auto& measurements = actor_state.measurements;

            // Geometry id
            const auto geo_id = navigation.current_object();
            // Container index
            const auto idx = prefix[geo_id];
            // Get measurements per module
            const auto measurements_per_module = measurements.at(idx).items;

            scalar_type min_chi2 = std::numeric_limits<scalar_type>::max();
            std::size_t min_idx = 0;

            const auto n_measurements = measurements_per_module.size();

            for (const auto& ms : measurements_per_module) {
            }

            /*
            for (std::size_t i_ms = 0; i_ms < n_measurements; i_ms++) {
                const auto& ms = measurements_per_module[i_ms];
                track_candidate cand{geo_id, ms};
                track_state trk_st{cand};
                // Set full jacobian
                trk_state.jacobian() = stepping._full_jacobian;
                auto det = navigation.detector();
                const auto& mask_store = det->mask_store();
                // Surface
                const auto& surface =
                    det->surface_by_index(trk_state.surface_link());
                // Run kalman updater
                mask_store.template call<gain_matrix_updater<algebra_t>>(
                    surface.mask(), trk_state, propagation);
                const auto chi2 = trk_state.filtered_chi2();
                if (chi2 < min_chi2) {
                    min_chi2 = chi2;
                    min_idx = i_ms;
                }
            }
            */
        }
    }
};

}  // namespace traccc 
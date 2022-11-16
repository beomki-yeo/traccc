/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/track_candidate.hpp"
#include "traccc/edm/track_state.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"
#include "traccc/utils/algorithm.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// traccc library include(s).
#include "traccc/utils/memory_resource.hpp"

namespace traccc::cuda {

/// Fitting algorithm for a set of tracks
template <typename fitter_t, typename seed_parameter_t>
class fitting_algorithm : public algorithm<track_state_container_types::buffer(
                              const typename track_candidates_container_types<
                                  seed_parameter_t>::const_view&)> {

    public:
    using detector_type = typename fitter_t::detector_type;
    using transform3_type = typename fitter_t::transform3_type;

    /// Constructor with a detector
    fitting_algorithm(const detector_type& det)
        : m_detector(std::make_unique<detector_type>(det)) {}

    /// Operator executing the algorithm
    ///
    /// @param track_candidates_view is a view of candidate measurements from
    /// track finding
    /// @return the buffer of the fitted track parameters
    track_state_container_types::buffer operator()(
        const typename track_candidates_container_types<seed_parameter_t>::
            const_view& track_candidates_view) const override {

        fitter_t fitter(*m_detector.get());

        track_state_container_types::host trk_states;

        /*
        track_state_container_types::host trk_states;

        // The number of tracks
        std::size_t n_tracks = track_candidates.size();

        // Iterate over tracks
        for (std::size_t i = 0; i < n_tracks; i++) {

            // Seed parameter
            const auto& seed_param = track_candidates[i].header;

            // Make a vector of track state
            auto& cands = track_candidates[i].items;
            vecmem::vector<track_state<transform3_type>> track_states;
            for (auto& cand : cands) {
                track_states.emplace_back(cand);
            }

            // Run fitter
            fitter.run(seed_param, std::move(track_states));

            trk_states.push_back(std::move(fitter.get_fitter_info()),
                                 std::move(fitter.get_track_states()));
        }

        return trk_states;
        */
    }

    private:
    traccc::memory_resource m_mr;
    std::unique_ptr<detector_type> m_detector;
};

}  // namespace traccc::cuda
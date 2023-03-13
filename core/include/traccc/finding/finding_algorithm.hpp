/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/track_candidate.hpp"
#include "traccc/finding/combinatorial_kalman_filter/combinatorial_kalman_finder.hpp"
#include "traccc/utils/algorithm.hpp"

namespace traccc {

/// Track Finding algorithm for a set of tracks
template <typename finder_t>
class finding_algorithm
    : public algorithm<track_candidate_container_types::host(
          const typename finder_t::detector_type&,
          const typename measurement_container_types::host&,
          const typename bound_track_parameters_collection_types::host&)> {

    public:
    using transform3_type = typename finder_t::transform3_type;

    // Measurement link_type
    using link_type = typename measurement_container_types::host::link_type;

    struct config {
        std::size_t max_num_branches_per_seed = 128;
    };

    // Run the algorithm
    track_candidate_container_types::host operator()(
        const typename finder_t::detector_type& det,
        const typename measurement_container_types::host& measurements,
        const typename bound_track_parameters_collection_types::host&
            seed_params),
        const override {

        finder_t finder(det);

        track_candidate_container_types::host output_candidates;

        // The number of seeds
        const std::size_t n_seeds = seed_params.size();

        // Iterate over tracks
        for (std::size_t i = 0; i < n_seeds; i++) {

            // Seed parameter
            const auto& seed_param = seed_params[i];

            // Make a finder state
            typename finder_t::state finder_state{};

            finder.find(seed_param, measurements, output_candidates,
                        finder_state);
        }

        return output_candidates;
    }

    private:
    // Configuration object
    config m_cfg;
};

}  // namespace traccc 
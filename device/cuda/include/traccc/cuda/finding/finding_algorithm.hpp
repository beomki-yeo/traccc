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
#include "traccc/utils/algorithm.hpp"

// detray include(s).
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/propagator.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// traccc library include(s).
#include "traccc/utils/memory_resource.hpp"

// Thrust Library
#include <thrust/pair.h>

namespace traccc::cuda {

/// Track Finding algorithm for a set of tracks
template <typename stepper_t, typename navigator_t>
class finding_algorithm
    : public algorithm<track_candidate_container_types::buffer(
          const typename navigator_t::detector_type::detector_view_type&,
          const vecmem::data::jagged_vector_view<
              typename navigator_t::intersection_type>&,
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

    struct config {
        unsigned int max_num_branches_per_seed = 10;
        unsigned int max_num_branches_per_surface = 10;
        unsigned int max_track_candidates_per_track = 10;
        unsigned int min_track_candidates_per_track = 6;
        unsigned int n_avg_threads_per_track = 4;
        scalar_type min_step_length_for_surface_aborter =
            0.1f * detray::unit<scalar_type>::mm;
        scalar_type chi2_max = 15.;
    };

    /// Get config object
    config& get_config() { return m_cfg; }

    const config& get_config() const { return m_cfg; }

    /// Run the algorithm
    track_candidate_container_types::buffer operator()(
        const typename detector_type::detector_view_type& det_view,
        const vecmem::data::jagged_vector_view<
            typename navigator_t::intersection_type>& navigation_buffer,
        const typename measurement_container_types::const_view& measurements,
        bound_track_parameters_collection_types::buffer&& seeds) const override;

    private:
    /// Memory resource used by the algorithm
    traccc::memory_resource m_mr;
    /// Copy object used by the algorithm
    std::unique_ptr<vecmem::copy> m_copy;
    /// Config object
    config m_cfg;
};

}  // namespace traccc::cuda
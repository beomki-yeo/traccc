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
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// traccc library include(s).
#include "traccc/utils/memory_resource.hpp"

namespace traccc::cuda {

template <typename metadata>
void fit(detray::detector<metadata>& det_data,
         track_candidate_container_types::const_view track_candidates_view,
         track_state_container_types::view track_states_view);

/// Fitting algorithm for a set of tracks
template <typename fitter_t, typename host_detector_t>
class fitting_algorithm
    : public algorithm<track_state_container_types::buffer(
          host_detector_t&&,
          const typename track_candidate_container_types::const_view&)> {

    public:
    using device_detector_type = typename fitter_t::detector_type;
    using transform3_type = typename fitter_t::transform3_type;

    /// Constructor for the fitting algorithm
    ///
    /// @param mr The memory resource to use
    fitting_algorithm(const traccc::memory_resource& mr) : m_mr(mr) {

        // Initialize m_copy ptr based on memory resources that were given
        if (mr.host) {
            m_copy = std::make_unique<vecmem::cuda::copy>();
        } else {
            m_copy = std::make_unique<vecmem::copy>();
        }
    };

    /// Run the algorithm
    track_state_container_types::buffer operator()(
        host_detector_t&& det,
        const typename track_candidate_container_types::const_view&
            track_candidates_view) const override {

        // Number of tracks
        const track_candidate_container_types::const_device::header_vector::
            size_type n_tracks =
                m_copy->get_size(track_candidates_view.headers);

        // Get the sizes of the track candidates in each track
        const std::vector<track_candidate_container_types::const_device::
                              item_vector::value_type::size_type>
            candidate_sizes = m_copy->get_sizes(track_candidates_view.items);

        track_state_container_types::buffer track_states_buffer{
            {n_tracks, m_mr.main},
            {std::vector<std::size_t>(n_tracks),
             std::vector<std::size_t>(candidate_sizes.begin(),
                                      candidate_sizes.end()),
             m_mr.main, m_mr.host}};

        m_copy->setup(track_states_buffer.headers);
        m_copy->setup(track_states_buffer.items);

        fit(det, track_candidates_view, track_states_buffer);

        return track_states_buffer;
    };

    private:
    /// Memory resource used by the algorithm
    traccc::memory_resource m_mr;
    /// Copy object used by the algorithm
    std::unique_ptr<vecmem::copy> m_copy;
};

}  // namespace traccc::cuda
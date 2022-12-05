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

// traccc library include(s).
#include "traccc/utils/memory_resource.hpp"

namespace traccc::cuda {

/// Fitting algorithm for a set of tracks
template <typename fitter_t>
class fitting_algorithm
    : public algorithm<track_state_container_types::buffer(
          const typename fitter_t::detector_type&,
          const typename track_candidate_container_types::const_view&)> {

    public:
    using transform3_type = typename fitter_t::transform3_type;

    /// Constructor for the fitting algorithm
    ///
    /// @param mr The memory resource to use
    ///
    fitting_algorithm(const traccc::memory_resource& mr) : m_mr(mr) {}

    /// Run the algorithm
    track_state_container_types::buffer operator()(
        const typename fitter_t::detector_type& det,
        const typename track_candidate_container_types::const_view&
            track_candidates_view) const override {

        fitter_t fitter(det);

        // Number of tracks
        const track_candidate_container_types::const_device::header_vector::
            size_type n_tracks =
                m_copy->get_size(track_candidates_view.headers);

        // Get the sizes of the track candidates in each track
        const std::vector<track_candidate_containter_types::const_device::
                              item_vector::value_type::size_type>
            candidate_sizes = m_copy->get_sizes(track_candidates_view.items);

        track_state_container_types::buffer trk_states_buffer{
            {n_tracks, m_mr.main},
            {std::vector<fitter_info>(n_tracks),
             std::vector<track_state>(candidate_sizes.begin(),
                                      candidate_sizes.end()),
             m_mr.main, m_mr.host}};

        m_copy->setup(trk_states_buffer.headers);
        m_copy->setup(trk_states_buffer.items);
    }

    private:
    /// Memory resource used by the algorithm
    traccc::memory_resource m_mr;
    /// Copy object used by the algorithm
    std::unique_ptr<vecmem::copy> m_copy;
};

}  // namespace traccc::cuda
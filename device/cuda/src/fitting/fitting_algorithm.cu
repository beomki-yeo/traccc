/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/fitting/fitting_algorithm.hpp"

// System include(s).
#include <vector>

namespace traccc::cuda {

template <typename fitter_t>
fitting_algorithm<fitter_t>::fitting_algorithm(
    const traccc::memory_resource& mr)
    : m_mr(mr) {}

template <typename fitter_t>
track_state_container_types::buffer fitting_algorithm<fitter_t>::operator()(
    const typename fitter_t::detector_type& det,
    const typename track_candidate_container_types::const_view&
        track_candidates_view) const {

    fitter_t fitter(det);

    // Number of tracks
    const track_candidate_container_types::const_device::header_vector::
        size_type n_tracks = m_copy->get_size(track_candidates_view.headers);

    // Get the sizes of the track candidates in each track
    const std::vector<track_candidate_container_types::const_device::
                          item_vector::value_type::size_type>
        candidate_sizes = m_copy->get_sizes(track_candidates_view.items);

    track_state_container_types::buffer trk_states_buffer{
        {n_tracks, m_mr.main},
        {std::vector<fitter_info<transform3_type>>(n_tracks),
         std::vector<track_state<transform3_type>>(candidate_sizes.begin(),
                                                   candidate_sizes.end()),
         m_mr.main, m_mr.host}};

    m_copy->setup(trk_states_buffer.headers);
    m_copy->setup(trk_states_buffer.items);

    // Calculate the number of threads and thread blocks to run the track
    // fitting
    /*
    const unsigned int nThreads = WARP_SIZE * 2;
    const unsigned int nBlocks =
        (track_states_buffer.size() + nThreads - 1) / nThreads;

    // Run the track fitting

    kernels::fit<<<nBlocks, nThreads>>>(track_states_buffer,
                                        track_candidates_view);
    */

    return trk_states_buffer;
}

}  // namespace traccc::cuda
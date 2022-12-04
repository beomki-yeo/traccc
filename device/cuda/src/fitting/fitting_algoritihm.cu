/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/fitting/fitting_algorithm.hpp"

namespace traccc::cuda {

template <typename fitter_t>
fitting_algorithm<fitter_t>::fitting_algorithm(const detector_type& det)
    : m_detector(std::make_unique<detector_type>(det)) {}

template <typename fitter_t>
track_state_container_types::buffer fitting_algorithm<fitter_t>::operator()(
    const typename track_candidates_container_types::const_view&
        track_candidates_view) const {

    // Get the sizes from the track candidates view
    auto track_sizes = m_copy->get_sizes(track_candidates_view);

    // Set up the track states buffer
    track_state_container_types::buffer track_states_buffer(
        track_sizes, *m_copy, m_mr.main, m_mr.host);

    // Calculate the number of threads and thread blocks to run the track
    // fitting
    const unsigned int nThreads = WARP_SIZE * 2;
    const unsigned int nBlocks =
        (track_states_buffer.size() + nThreads - 1) / nThreads;

    // Run the track fitting
    /*
    kernels::fit<<<nBlocks, nThreads>>>(track_states_buffer,
                                        track_candidates_view);
    */

    return track_states_buffer;
}

}  // namespace traccc::cuda
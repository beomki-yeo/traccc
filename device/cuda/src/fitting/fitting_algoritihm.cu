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

    // Get the size of the track candidates view
    const std::size_t track_size = m_copy->get_size(track_candidates_view);

    // Create device buffer for the track states
    track_state_container_types::buffer track_states_buffer(track_size,
                                                            m_mr.main);

    return track_states_buffer;
}

}  // namespace traccc::cuda
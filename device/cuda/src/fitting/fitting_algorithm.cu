/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/fitting/fitting_algorithm.hpp"
#include "traccc/cuda/utils/definitions.hpp"

// detray include(s).
#include "detray/detectors/detector_metadata.hpp"
#include "detray/propagator/rk_stepper.hpp"

// System include(s).
#include <vector>

namespace traccc::cuda {

namespace kernels {

template <typename metadata>
__global__ void fit(
    detray::detector_view<detray::detector<metadata>> det_data,
    track_candidate_container_types::const_view track_candidates_view,
    track_state_container_types::view track_states_view) {

    using device_detector_type =
        detray::detector<metadata, covfie::field_view,
                         detray::device_container_types>;
}

}  // namespace kernels

template <typename metadata>
void fit(detray::detector<metadata>& det,
         track_candidate_container_types::const_view track_candidates_view,
         track_state_container_types::view track_states_view) {

    // Calculate the number of threads and thread blocks to run the track
    // fitting
    const unsigned int nThreads = WARP_SIZE * 2;
    const unsigned int nBlocks =
        (track_candidates_view.headers.size() + nThreads - 1) / nThreads;

    auto det_data = detray::get_data(det);

    // Run the track fitting
    kernels::fit<<<nBlocks, nThreads>>>(det_data, track_candidates_view,
                                        track_states_view);
}

// Explicit instantiation of template function
template void fit<detray::detector_registry::telescope_detector>(
    detray::detector<detray::detector_registry::telescope_detector>& det,
    track_candidate_container_types::const_view track_candidates_view,
    track_state_container_types::view track_states_view);

}  // namespace traccc::cuda
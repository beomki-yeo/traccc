/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/fitting/fitting_algorithm.hpp"
#include "traccc/cuda/utils/definitions.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"

// detray include(s).
#include "detray/detectors/detector_metadata.hpp"
#include "detray/propagator/rk_stepper.hpp"

// System include(s).
#include <vector>

namespace traccc::cuda {

namespace kernels {

template <typename fitter_t, typename detector_view_t>
__global__ void fit(
    detector_view_t det_data,
    vecmem::data::jagged_vector_view<typename fitter_t::intersection_type>
        nav_candidates_buffer,
    track_candidate_container_types::const_view track_candidates_view,
    track_state_container_types::view track_states_view) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    typename fitter_t::detector_type det(det_data);
    vecmem::jagged_device_vector<typename fitter_t::intersection_type>
        nav_candidates(nav_candidates_buffer);

    track_candidate_container_types::const_device track_candidates(
        track_candidates_view);

    track_state_container_types::device track_states(track_states_view);

    fitter_t fitter(det);

    if (gid >= track_states.size()) {
        return;
    }

    // Track candidates per track
    const auto& track_candidates_per_track = track_candidates[gid].items;

    // Seed parameter
    const auto& seed_param = track_candidates[gid].header;

    // Track states per track
    auto track_states_per_track = track_states[gid].items;

    for (auto& cand : track_candidates_per_track) {
        track_states_per_track.emplace_back(cand);
    }

    typename fitter_t::state fitter_state(track_states_per_track);

    fitter.fit(seed_param, fitter_state, nav_candidates.at(gid));
}

}  // namespace kernels

template <typename fitter_t, typename host_detector_t>
fitting_algorithm<fitter_t, host_detector_t>::fitting_algorithm(
    const traccc::memory_resource& mr)
    : m_mr(mr) {

    // Initialize m_copy ptr based on memory resources that were given
    if (mr.host) {
        m_copy = std::make_unique<vecmem::cuda::copy>();
    } else {
        m_copy = std::make_unique<vecmem::copy>();
    }
};

template <typename fitter_t, typename host_detector_t>
track_state_container_types::buffer
fitting_algorithm<fitter_t, host_detector_t>::operator()(
    host_detector_t&& det,
    const typename track_candidate_container_types::const_view&
        track_candidates_view) const {

    // Number of tracks
    const track_candidate_container_types::const_device::header_vector::
        size_type n_tracks = m_copy->get_size(track_candidates_view.headers);

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

    // Calculate the number of threads and thread blocks to run the track
    // fitting
    const unsigned int nThreads = WARP_SIZE * 2;
    const unsigned int nBlocks = (n_tracks + nThreads - 1) / nThreads;

    auto det_data = detray::get_data(det);

    // Create navigator candidates buffer
    auto nav_candidates_buffer =
        detray::create_candidates_buffer(det, n_tracks, m_mr.main, m_mr.host);
    m_copy->setup(nav_candidates_buffer);

    // Run the track fitting
    kernels::fit<fitter_t>
        <<<nBlocks, nThreads>>>(det_data, nav_candidates_buffer,
                                track_candidates_view, track_states_buffer);
    return track_states_buffer;
}

// Explicit template instantiation
using host_detector_type =
    detray::detector<detray::detector_registry::telescope_detector,
                     covfie::field>;
using device_detector_type =
    detray::detector<detray::detector_registry::telescope_detector,
                     covfie::field_view, detray::device_container_types>;
using b_field_t = typename host_detector_type::bfield_type;
using rk_stepper_type = detray::rk_stepper<b_field_t::view_t, transform3,
                                           detray::constrained_step<>>;
using device_navigator_type = detray::navigator<const device_detector_type>;
using device_fitter_type =
    kalman_fitter<rk_stepper_type, device_navigator_type>;
template class fitting_algorithm<device_fitter_type, host_detector_type>;

/*
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
*/

}  // namespace traccc::cuda
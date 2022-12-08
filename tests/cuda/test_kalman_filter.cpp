/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "tests/seed_generator.hpp"
#include "traccc/cuda/fitting/fitting_algorithm.hpp"
#include "traccc/device/container_h2d_copy_alg.hpp"
#include "traccc/edm/track_state.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/resolution/fitting_performance_writer.hpp"
#include "traccc/utils/memory_resource.hpp"

// Test include(s).
#include "tests/kalman_fitting_test.hpp"

// detray include(s).
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/simulation/simulator.hpp"
#include "detray/simulation/track_generators.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <climits>

using namespace traccc;

// This defines the local frame test suite
TEST_P(KalmanFittingTests, Run) {

    // Test Parameters
    const scalar p0 = std::get<0>(GetParam());
    const scalar phi0 = std::get<1>(GetParam());

    // File path
    std::string file_path =
        std::to_string(p0) + "_GeV_" + std::to_string(phi0) + "_phi/";
    std::string full_path =
        "detray_simulation/telescope/kf_validation/" + file_path;

    /*****************************
     * Build a telescope geometry
     *****************************/

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::device_memory_resource device_mr;
    traccc::memory_resource mr{device_mr, &host_mr};
    vecmem::cuda::managed_memory_resource mng_mr;

    host_detector_type host_det = create_telescope_detector(
        mng_mr,
        b_field_t(b_field_t::backend_t::configuration_t{B[0], B[1], B[2]}),
        plane_positions, traj, std::numeric_limits<scalar>::infinity(),
        std::numeric_limits<scalar>::infinity(), mat, thickness);

    /***************
     * Run fitting
     ***************/

    vecmem::cuda::copy copy;

    traccc::device::container_h2d_copy_alg<
        traccc::track_candidate_container_types>
        track_candidate_h2d{mr, copy};

    // Seed generator
    seed_generator<rk_stepper_type, host_navigator_type> sg(host_det, stddevs);

    // Fitting algorithm object
    traccc::fitting_algorithm<host_fitter_type> host_fitting;

    traccc::cuda::fitting_algorithm<device_fitter_type, host_detector_type>
        device_fitting(mr);

    std::size_t n_events = 100;

    // Iterate over events
    for (std::size_t i_evt = 0; i_evt < n_events; i_evt++) {
        // Event map
        traccc::event_map2 evt_map(i_evt, full_path, full_path, full_path);

        // Truth Track Candidates
        traccc::track_candidate_container_types::host track_candidates =
            evt_map.generate_truth_candidates(sg, mng_mr);

        // n_trakcs = 100
        ASSERT_EQ(track_candidates.size(), 100);

        // Run fitting (host)
        auto host_track_states = host_fitting(host_det, track_candidates);

        // Run fitting (device)
        const traccc::track_candidate_container_types::buffer
            track_candidates_cuda_buffer =
                track_candidate_h2d(traccc::get_data(track_candidates));

        auto device_track_states =
            device_fitting(std::move(host_det), track_candidates_cuda_buffer);
    }
}

INSTANTIATE_TEST_SUITE_P(
    KalmanFitValidation, KalmanFittingTests,
    ::testing::Values(std::make_tuple(1 * detray::unit<scalar>::GeV, 0),
                      std::make_tuple(10 * detray::unit<scalar>::GeV, 0),
                      std::make_tuple(100 * detray::unit<scalar>::GeV, 0)));

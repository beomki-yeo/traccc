/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "tests/seed_generator.hpp"
#include "traccc/edm/track_state.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/resolution/fitting_performance_writer.hpp"

// Test include(s).
#include "tests/kalman_fitting_test.hpp"

// detray include(s).
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/simulation/simulator.hpp"
#include "detray/simulation/track_generators.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

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

    // Memory resource(s)
    vecmem::cuda::managed_memory_resource mng_mr;

    const host_detector_type host_det = create_telescope_detector(
        mng_mr,
        b_field_t(b_field_t::backend_t::configuration_t{B[0], B[1], B[2]}),
        plane_positions, traj, std::numeric_limits<scalar>::infinity(),
        std::numeric_limits<scalar>::infinity(), mat, thickness);

    /***************
     * Run fitting
     ***************/

    // Seed generator
    seed_generator<rk_stepper_type, host_navigator_type> sg(host_det, stddevs);

    // Fitting algorithm object
    fitting_algorithm<host_fitter_type> host_fitting;
    // fitting_algorithm<device_fitter_type> device_fitting(mng_mr);

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
        // auto device_track_states = device_fitting();
    }
}

INSTANTIATE_TEST_SUITE_P(
    KalmanFitValidation, KalmanFittingTests,
    ::testing::Values(std::make_tuple(1 * detray::unit<scalar>::GeV, 0),
                      std::make_tuple(10 * detray::unit<scalar>::GeV, 0),
                      std::make_tuple(100 * detray::unit<scalar>::GeV, 0)));

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/finding/finding_algorithm.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/efficiency/finding_performance_writer.hpp"
#include "traccc/simulation/simulator.hpp"
#include "traccc/utils/ranges.hpp"

// Test include(s).
#include "tests/ckf_toy_detector_test.hpp"
#include "traccc/utils/seed_generator.hpp"

// detray include(s).
#include "detray/io/common/detector_reader.hpp"
#include "detray/io/common/detector_writer.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <filesystem>
#include <string>

using namespace traccc;

TEST_P(CkfEfficiencyToyDetectorTests, Run) {

    // Get the parameters
    const std::string name = std::get<0>(GetParam());
    const std::array<scalar, 3u> origin = std::get<1>(GetParam());
    const std::array<scalar, 3u> origin_stddev = std::get<2>(GetParam());
    const std::array<scalar, 2u> mom_range = std::get<3>(GetParam());
    const std::array<scalar, 2u> eta_range = std::get<4>(GetParam());
    const std::array<scalar, 2u> theta_range = eta_to_theta_range(eta_range);
    const std::array<scalar, 2u> phi_range = std::get<5>(GetParam());
    const unsigned int n_truth_tracks = std::get<6>(GetParam());
    const unsigned int n_events = std::get<7>(GetParam());

    // Performance writer
    traccc::find_performance_writer::config find_writer_cfg;
    find_writer_cfg.file_path = "performance_track_finding_" + name + ".root";
    traccc::fitting_performance_writer fit_performance_writer(fit_writer_cfg);

    // Remove the data
    std::filesystem::remove_all(full_path);
}

/*
INSTANTIATE_TEST_SUITE_P(
    CkfEfficiencyToyDetectorValidation0, CkfEfficiencyToyDetectorTests,
    ::testing::Values(std::make_tuple(
        "toy_detector", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 200.f, 200.f},
        std::array<scalar, 2u>{1.f, 1.f}, std::array<scalar, 2u>{0.f, 0.f},
        std::array<scalar, 2u>{0.f, 0.f}, 1, 5000)));
*/
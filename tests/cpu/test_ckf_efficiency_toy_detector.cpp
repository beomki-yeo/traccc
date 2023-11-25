/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/efficiency/finding_performance_writer.hpp"
#include "traccc/finding/finding_algorithm.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/utils.hpp"
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
    const scalar charge = std::get<6>(GetParam());
    const unsigned int n_truth_tracks = std::get<7>(GetParam());
    const unsigned int n_events = std::get<8>(GetParam());

    // Performance writer
    traccc::finding_performance_writer::config find_writer_cfg;
    find_writer_cfg.file_path = "performance_track_finding_" + name + ".root";
    traccc::finding_performance_writer find_performance_writer(find_writer_cfg);

    /*****************************
     * Build a toy detector
     *****************************/

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;

    // Read back detector file
    const std::string path = name + "/";
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg.add_file(path + "toy_detector_geometry.json")
        .add_file(path + "toy_detector_homogeneous_material.json")
        .add_file(path + "toy_detector_surface_grids.json");

    const auto [host_det, names] =
        detray::io::read_detector<host_detector_type>(host_mr, reader_cfg);

    auto field = detray::bfield::create_const_field(B);

    /***************************
     * Generate simulation data
     ***************************/

    // Track generator
    using generator_type =
        detray::random_track_generator<traccc::free_track_parameters,
                                       uniform_gen_t>;
    generator_type::configuration gen_cfg{};
    gen_cfg.n_tracks(n_truth_tracks);
    gen_cfg.origin(origin);
    gen_cfg.origin_stddev(origin_stddev);
    gen_cfg.phi_range(phi_range[0], phi_range[1]);
    gen_cfg.theta_range(theta_range[0], theta_range[1]);
    gen_cfg.mom_range(mom_range[0], mom_range[1]);
    gen_cfg.charge(charge);
    generator_type generator(gen_cfg);

    // Smearing value for measurements
    traccc::measurement_smearer<transform3> meas_smearer(smearing[0],
                                                         smearing[1]);

    using writer_type =
        traccc::smearing_writer<traccc::measurement_smearer<transform3>>;

    typename writer_type::config smearer_writer_cfg{meas_smearer};

    // Run simulator
    const std::string full_path = io::data_directory() + path;
    std::filesystem::create_directories(full_path);
    auto sim = traccc::simulator<host_detector_type, b_field_t, generator_type,
                                 writer_type>(
        n_events, host_det, field, std::move(generator),
        std::move(smearer_writer_cfg), full_path);

    // Set constrained step size to 2 mm
    sim.get_config().step_constraint = step_constraint;
    sim.get_config().overstep_tolerance = overstep_tolerance;

    sim.run();

    /*****************************
     * Do the reconstruction
     *****************************/

    // Seed generator
    seed_generator<host_detector_type> sg(host_det, stddevs);

    // Finding algorithm configuration
    typename traccc::finding_algorithm<rk_stepper_type,
                                       host_navigator_type>::config_type cfg;

    // Finding algorithm object
    traccc::finding_algorithm<rk_stepper_type, host_navigator_type>
        host_finding(cfg);

    // Iterate over events
    for (std::size_t i_evt = 0; i_evt < n_events; i_evt++) {

        // Truth Track Candidates
        traccc::event_map2 evt_map(i_evt, path, path, path);

        traccc::track_candidate_container_types::host truth_track_candidates =
            evt_map.generate_truth_candidates(sg, host_mr);

        ASSERT_EQ(truth_track_candidates.size(), n_truth_tracks);

        // Prepare truth seeds
        traccc::bound_track_parameters_collection_types::host seeds(&host_mr);
        for (unsigned int i_trk = 0; i_trk < n_truth_tracks; i_trk++) {
            seeds.push_back(truth_track_candidates.at(i_trk).header);
        }
        ASSERT_EQ(seeds.size(), n_truth_tracks);

        // Read measurements
        traccc::io::measurement_reader_output readOut(&host_mr);
        traccc::io::read_measurements(readOut, i_evt, path,
                                      traccc::data_format::csv);
        traccc::measurement_collection_types::host& measurements_per_event =
            readOut.measurements;

        // Run finding
        auto track_candidates =
            host_finding(host_det, field, measurements_per_event, seeds);

        find_performance_writer.write(traccc::get_data(track_candidates),
                                      evt_map);
    }

    find_performance_writer.finalize();

    // Check the efficiency
    ckf_efficiency_tests(find_writer_cfg.file_path, eta_range[0], eta_range[1], 0.95f);
}

INSTANTIATE_TEST_SUITE_P(
    CkfEfficiencyToyDetectorValidation0, CkfEfficiencyToyDetectorTests,
    ::testing::Values(std::make_tuple(
        "toy_1_GeV", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f}, std::array<scalar, 2u>{1.f, 1.f},
        std::array<scalar, 2u>{-3.f, 3.f},
        std::array<scalar, 2u>{0.f, 2.0f * detray::constant<scalar>::pi}, -1.f,
        4000, 1)));

INSTANTIATE_TEST_SUITE_P(
    CkfEfficiencyToyDetectorValidation1, CkfEfficiencyToyDetectorTests,
    ::testing::Values(std::make_tuple(
        "toy_10_GeV", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 2u>{10.f, 10.f}, std::array<scalar, 2u>{-3.f, 3.f},
        std::array<scalar, 2u>{0.f, 2.0f * detray::constant<scalar>::pi}, -1.f,
        4000, 1)));

INSTANTIATE_TEST_SUITE_P(
    CkfEfficiencyToyDetectorValidation2, CkfEfficiencyToyDetectorTests,
    ::testing::Values(std::make_tuple(
        "toy_100_GeV", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 2u>{100.f, 100.f}, std::array<scalar, 2u>{-3.f, 3.f},
        std::array<scalar, 2u>{0.f, 2.0f * detray::constant<scalar>::pi}, -1.f,
        4000, 1)));

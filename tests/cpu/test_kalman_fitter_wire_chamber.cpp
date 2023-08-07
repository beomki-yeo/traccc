/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/edm/track_state.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/resolution/fitting_performance_writer.hpp"
#include "traccc/utils/ranges.hpp"
#include "traccc/utils/seed_generator.hpp"

// Test include(s).
#include "tests/kalman_fitting_wire_chamber_test.hpp"

// detray include(s).
#include "detray/detectors/create_wire_chamber.hpp"
#include "detray/io/common/detector_reader.hpp"
#include "detray/io/common/detector_writer.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/simulation/simulator.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <filesystem>
#include <string>

using namespace traccc;

TEST_P(KalmanFittingWireChamberTests, Run) {

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
    traccc::fitting_performance_writer::config fit_writer_cfg;
    fit_writer_cfg.file_path = "performance_track_fitting_" + name + ".root";
    traccc::fitting_performance_writer fit_performance_writer(fit_writer_cfg);

    /*****************************
     * Build a Wire Chamber
     *****************************/

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;

    detray::wire_chamber_config wire_chamber_cfg;
    wire_chamber_cfg.n_layers(n_wire_layers);
    wire_chamber_cfg.bfield_vec(B);

    // Create telescope detector
    auto [det, name_map] = create_wire_chamber(host_mr, wire_chamber_cfg);

    // Write detector file
    auto writer_cfg = detray::io::detector_writer_config{}
                          .format(detray::io::format::json)
                          .replace_files(true);
    detray::io::write_detector(det, name_map, writer_cfg);

    // Read back detector file
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg.add_file("wire_chamber_geometry.json")
        .add_file("wire_chamber_homogeneous_material.json")
        .bfield_vec(B[0], B[1], B[2]);

    const auto [host_det, names] =
        detray::io::read_detector<host_detector_type>(host_mr, reader_cfg);

    /***************************
     * Generate simulation data
     ***************************/

    auto generator =
        detray::random_track_generator<traccc::free_track_parameters,
                                       uniform_gen_t>(n_truth_tracks, origin,
                                                      origin_stddev, mom_range,
                                                      theta_range, phi_range);

    // Smearing value for measurements
    detray::measurement_smearer<transform3> meas_smearer(smearing[0],
                                                         smearing[1]);

    // Run simulator
    const std::string path = name + "/";
    const std::string full_path = io::data_directory() + path;
    std::filesystem::create_directories(full_path);
    auto sim = detray::simulator(n_events, host_det, std::move(generator),
                                 meas_smearer, full_path);
    sim.run();

    /***************
     * Run fitting
     ***************/

    // Seed generator
    seed_generator<host_detector_type> sg(host_det, stddevs);

    // Fitting algorithm object
    typename traccc::fitting_algorithm<host_fitter_type>::config_type fit_cfg;
    fitting_algorithm<host_fitter_type> fitting(fit_cfg);

    // Iterate over events
    for (std::size_t i_evt = 0; i_evt < n_events; i_evt++) {

        std::cout << "Event ID: " << i_evt << std::endl;

        // Event map
        traccc::event_map2 evt_map(i_evt, path, path, path);

        // Truth Track Candidates
        traccc::track_candidate_container_types::host track_candidates =
            evt_map.generate_truth_candidates(sg, host_mr);

        // n_trakcs = 100
        ASSERT_EQ(track_candidates.size(), n_truth_tracks);

        // Run fitting
        auto track_states = fitting(host_det, track_candidates);

        // Iterator over tracks
        const std::size_t n_tracks = track_states.size();

        // n_trakcs = 100
        ASSERT_EQ(n_tracks, n_truth_tracks);
        for (std::size_t i_trk = 0; i_trk < n_tracks; i_trk++) {

            std::cout << i_trk << std::endl;

            const auto& track_states_per_track = track_states[i_trk].items;
            ASSERT_EQ(track_states_per_track.size(), n_wire_layers);
            const auto& fit_info = track_states[i_trk].header;
            ASSERT_FLOAT_EQ(fit_info.ndf, n_wire_layers - 5.f);

            fit_performance_writer.write(track_states_per_track, fit_info,
                                         host_det, evt_map);
        }
    }

    fit_performance_writer.finalize();

    /********************
     * Pull value test
     ********************/

    static const std::vector<std::string> pull_names{
        "pull_d0", "pull_z0", "pull_phi", "pull_theta", "pull_qop"};
    pull_value_tests(fit_writer_cfg.file_path, pull_names);

    // Remove the data
    std::filesystem::remove_all(full_path);
}

/*
INSTANTIATE_TEST_SUITE_P(
    KalmanFitWireChamberValidation0, KalmanFittingWireChamberTests,
    ::testing::Values(std::make_tuple(
        "1_GeV_0_phi", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 2u>{1.f, 1.f}, std::array<scalar, 2u>{0.f, 0.f},
        std::array<scalar, 2u>{0.f, 0.f}, 100, 100)));

INSTANTIATE_TEST_SUITE_P(
    KalmanFitWireChamberValidation1, KalmanFittingWireChamberTests,
    ::testing::Values(std::make_tuple(
        "10_GeV_0_phi", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 2u>{10.f, 10.f}, std::array<scalar, 2u>{0.f, 0.f},
        std::array<scalar, 2u>{0.f, 0.f}, 100, 100)));
*/
INSTANTIATE_TEST_SUITE_P(
    KalmanFitWireChamberValidation2, KalmanFittingWireChamberTests,
    ::testing::Values(std::make_tuple("100_GeV_0_phi",
                                      std::array<scalar, 3u>{0.f, 0.f, 0.f},
                                      std::array<scalar, 3u>{0.f, 0.f, 0.f},
                                      std::array<scalar, 2u>{100.f, 100.f},
                                      std::array<scalar, 2u>{-0.5f, 0.5f},
                                      std::array<scalar, 2u>{0.f, 0.f}, 100,
                                      100)));

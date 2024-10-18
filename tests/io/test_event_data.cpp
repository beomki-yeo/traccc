/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/clusterization/measurement_creation_algorithm.hpp"
#include "traccc/clusterization/sparse_ccl_algorithm.hpp"
#include "traccc/io/data_format.hpp"
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/utils/event_data.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

TEST(event_data, acts_odd) {

    /// Type declarations
    using host_detector_type = detray::detector<detray::default_metadata,
                                                detray::host_container_types>;

    vecmem::host_memory_resource resource;

    const std::string path = "odd/geant4_1muon_100GeV";
    const std::string det_file =
        "geometries/odd/odd-detray_geometry_detray.json";

    // Read detector file
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg.add_file(traccc::io::data_directory() + det_file);

    auto [host_det, names] =
        detray::io::read_detector<host_detector_type>(resource, reader_cfg);

    {
        // without cell
        traccc::event_data evt_data(path, 0u, resource, true, &host_det,
                                    traccc::data_format::csv, false);
        EXPECT_EQ(evt_data.m_particle_map.size(), 4515u);
        EXPECT_EQ(evt_data.m_meas_to_ptc_map.size(), 58u);
        EXPECT_EQ(evt_data.m_meas_to_param_map.size(), 58u);
    }
    {
        // with cell
        traccc::event_data evt_data(path, 0u, resource, true, &host_det,
                                    traccc::data_format::csv, true);
        EXPECT_EQ(evt_data.m_particle_map.size(), 4515u);
        EXPECT_EQ(evt_data.m_meas_to_ptc_map.size(), 58u);
        EXPECT_EQ(evt_data.m_meas_to_param_map.size(), 58u);
    }
}

TEST(event_data, mock_data) {

    /***
     * Mock data test
     *
     * Mock data consists of three particles each of which has one hit
     *
     * first particle: one hit, three cells
     * second particle: one hit, four cells
     * thrid particle: one hit, three cells
     *
     * [ ] [1] [ ] [ ] [ ] [ ] [ ] [ ]
     * [1][1,2][2] [ ] [ ] [ ] [ ] [ ]
     * [ ] [2] [2] [ ] [ ] [ ] [ ] [ ]
     * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
     * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
     * [ ] [ ] [ ] [ ] [ ] [3] [3] [3]
     * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
     * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
     *
     * Current traccc's CCA algorithm will make two clusters
     *
     * first cluster is a set of cells generated by 1st and 2nd particles
     * second cluster is a set of cells generated by 3rd particle
     *
     */

    /// Type declarations
    using host_detector_type = detray::detector<detray::default_metadata,
                                                detray::host_container_types>;

    vecmem::host_memory_resource resource;

    const std::string path = TRACCC_TEST_IO_MOCK_DATA_DIR;

    // Dummy detector and digitization file
    const std::string det_file =
        "geometries/odd/odd-detray_geometry_detray.json";
    const std::string digi_file =
        "geometries/odd/odd-digi-geometric-config.json";

    // Read detector file
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg.add_file(traccc::io::data_directory() + det_file);

    auto [host_det, names] =
        detray::io::read_detector<host_detector_type>(resource, reader_cfg);

    traccc::event_data evt_data(path, 0u, resource, true, &host_det,
                                traccc::data_format::csv, true);

    // There are three measurements
    EXPECT_EQ(evt_data.m_meas_to_ptc_map.size(), 3u);
    EXPECT_EQ(evt_data.m_meas_to_param_map.size(), 3u);
    // There are three particles
    EXPECT_EQ(evt_data.m_particle_map.size(), 3u);
    EXPECT_EQ(evt_data.m_ptc_to_meas_map.size(), 3u);

    for (auto const& [meas, ptcs] : evt_data.m_meas_to_ptc_map) {
        EXPECT_EQ(ptcs.size(), 1);

        for (auto const& [ptc, count] : ptcs) {
            if (ptc.particle_id == 4503599644147712) {
                // number of cells from 1st particle
                EXPECT_EQ(count, 3);
            } else if (ptc.particle_id == 4503599660924928) {
                // number of cells from 2nd particle
                EXPECT_EQ(count, 4);
            } else if (ptc.particle_id == 4503599744811008) {
                // number of cells from 3rd particle
                EXPECT_EQ(count, 3);
            }
        }
    }

    /// Test with CCA

    // Construct the detector description object.
    traccc::silicon_detector_description::host det_descr{resource};
    traccc::io::read_detector_description(det_descr, det_file, digi_file,
                                          traccc::data_format::json);
    traccc::silicon_detector_description::data det_descr_data{
        vecmem::get_data(det_descr)};

    // Algorithms
    traccc::host::sparse_ccl_algorithm cc(resource);
    traccc::host::measurement_creation_algorithm mc(resource);

    // Read cells
    traccc::edm::silicon_cell_collection::host cells{resource};
    traccc::io::read_cells(cells, 0u, path, &det_descr,
                           traccc::data_format::csv);
    const auto cells_view = vecmem::get_data(cells);

    auto clusters = cc(cells_view);
    auto measurements =
        mc(cells_view, vecmem::get_data(clusters), det_descr_data);

    evt_data.fill_cca_result(cells, clusters, measurements, det_descr);

    EXPECT_EQ(evt_data.m_found_meas_to_ptc_map.size(), 2u);
    EXPECT_EQ(evt_data.m_found_meas_to_param_map.size(), 2u);

    bool has_first_cluster = false;
    bool has_second_cluster = false;

    bool has_first_particle = false;
    bool has_second_particle = false;
    bool has_third_particle = false;

    for (auto const& [meas, ptcs] : evt_data.m_found_meas_to_ptc_map) {

        // first measurement (or cluster) is contributed by 1st and 2nd
        // particles
        if (ptcs.size() == 2) {
            for (auto const& [ptc, count] : ptcs) {
                if (ptc.particle_id == 4503599644147712) {
                    // number of cells from 1st particle
                    EXPECT_EQ(count, 3);
                    has_first_particle = true;
                } else if (ptc.particle_id == 4503599660924928) {
                    // number of cells from 2nd particle
                    EXPECT_EQ(count, 4);
                    has_second_particle = true;
                }
            }

            has_first_cluster = true;
        }

        // second measurement (or cluster) is contributed by 3rd particle
        else if (ptcs.size() == 1) {
            for (auto const& [ptc, count] : ptcs) {
                if (ptc.particle_id == 4503599744811008) {
                    // number of cells from 3rd particle
                    EXPECT_EQ(count, 3);
                    has_third_particle = true;
                }
            }

            has_second_cluster = true;
        }
    }

    EXPECT_EQ(has_first_cluster, true);
    EXPECT_EQ(has_second_cluster, true);
    EXPECT_EQ(has_first_particle, true);
    EXPECT_EQ(has_second_particle, true);
    EXPECT_EQ(has_third_particle, true);

    bool has_first_param = false;
    bool has_second_param = false;

    for (auto const& [meas, param] : evt_data.m_found_meas_to_param_map) {
        if (meas.measurement_id == 0u) {
            has_first_param = true;
        } else if (meas.measurement_id == 1u) {
            has_second_param = true;
        }
    }

    EXPECT_EQ(has_first_param, true);
    EXPECT_EQ(has_second_param, true);
}

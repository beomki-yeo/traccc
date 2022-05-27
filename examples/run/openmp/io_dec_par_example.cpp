/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/clusterization/component_connection.hpp"
#include "traccc/clusterization/measurement_creation.hpp"
#include "traccc/clusterization/spacepoint_formation.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/spacepoint.hpp"
#include "traccc/io/demonstrator_edm.hpp"
#include "traccc/io/reader.hpp"
#include "traccc/io/writer.hpp"

// Boost
#include <boost/program_options.hpp>

#ifdef _OPENMP
#include "omp.h"
#endif

// System include(s).
#include <chrono>
#include <exception>
#include <iostream>

namespace po = boost::program_options;

traccc::demonstrator_result run(traccc::demonstrator_input input_data,
                                vecmem::host_memory_resource resource) {

    // Algorithms
    traccc::component_connection cc(resource);
    traccc::measurement_creation mc(resource);
    traccc::spacepoint_formation sf(resource);

    // Output stats
    uint64_t n_modules = 0;
    uint64_t n_cells = 0;
    uint64_t n_measurements = 0;
    uint64_t n_spacepoints = 0;

    auto startAlgorithms = std::chrono::system_clock::now();

    traccc::demonstrator_result aggregated_results;

#pragma omp parallel for reduction (+:n_modules, n_cells, n_measurements, n_spacepoints)
    for (size_t event = 0; event < input_data.size(); ++event) {
        traccc::cell_container_types::const_view cells_view =
            input_data.operator[](event);

        /*-------------------
            CCL
          -------------------*/

        auto clusters_view = cc(cells_view);

        /*------------------------
            Measurement Creation
          ------------------------*/

        auto measurements_view = mc(clusters_view);

        /*------------------------
            Spacepoint formation
          ------------------------*/

        auto spacepoints_view = sf(measurements_view);

        /*----------------------------
          Statistics
          ----------------------------*/

        const traccc::cell_container_types::const_device cells_per_event(
            cells_view);
        const traccc::measurement_container_types::const_device
            measurements_per_event(measurements_view);
        const traccc::spacepoint_container_types::const_device
            spacepoints_per_event(spacepoints_view);

        n_modules += cells_per_event.size();
        n_cells += cells_per_event.total_size();
        n_measurements += measurements_per_event.total_size();
        n_spacepoints += spacepoints_per_event.total_size();

#pragma omp critical
        aggregated_results.push_back(
            traccc::result({measurements_view, spacepoints_view}));
    }

    auto endAlgorithms = std::chrono::system_clock::now();
    std::chrono::duration<double> diffAlgo = endAlgorithms - startAlgorithms;
    std::cout << "Algorithms time: " << diffAlgo.count() << " sec."
              << std::endl;

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read    " << n_cells << " cells from " << n_modules
              << " modules" << std::endl;
    std::cout << "- created " << n_measurements << " measurements. "
              << std::endl;
    std::cout << "- created " << n_spacepoints << " spacepoints. " << std::endl;

    return aggregated_results;
}

// The main routine
int main(int argc, char *argv[]) {

    // Set up the program options.
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Give some help with the program's options");
    desc.add_options()("detector_file", po::value<std::string>()->required(),
                       "specify detector file");
    desc.add_options()("cell_directory", po::value<std::string>()->required(),
                       "specify the directory of cell files");
    desc.add_options()("events", po::value<int>()->required(),
                       "number of events");

    // Interpret the program options.
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // Print a help message if the user asked for it.
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    // Handle any and all errors.
    try {
        po::notify(vm);
    } catch (const std::exception &ex) {
        std::cerr << "Couldn't interpret command line options because of:\n\n"
                  << ex.what() << "\n\n"
                  << desc << std::endl;
        return 1;
    }

    auto detector_file = vm["detector_file"].as<std::string>();
    auto cell_directory = vm["cell_directory"].as<std::string>();
    auto data_format = vm["data_format"].as<std::string>();
    auto events = vm["events"].as<int>();

    std::cout << "Running " << argv[0] << " " << detector_file << " "
              << cell_directory << " " << events << std::endl;

    auto start = std::chrono::system_clock::now();
    vecmem::host_memory_resource resource;
    set_default_resource(&resource);

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Total execution time: " << diff.count() << " sec."
              << std::endl;
}

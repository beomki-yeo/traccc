/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"

// io
#include "traccc/io/read_geometry.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/utils.hpp"

// algorithms
#include "traccc/seeding/experimental/spacepoint_formation.hpp"

// options
#include "traccc/options/common_options.hpp"
#include "traccc/options/detector_input_options.hpp"
#include "traccc/options/handle_argument_errors.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/io/common/detector_reader.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <iostream>

using namespace traccc;
namespace po = boost::program_options;

// The main routine
//
int main(int argc, char* argv[]) {
    // Set up the program options
    po::options_description desc("Allowed options");

    // Add options
    desc.add_options()("help,h", "Give some help with the program's options");
    traccc::common_options common_opts(desc);
    traccc::detector_input_options det_opts(desc);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // Check errors
    traccc::handle_argument_errors(vm, desc);

    // Read options
    common_opts.read(vm);
    det_opts.read(vm);

    std::cout << "Running " << argv[0] << " " << common_opts.input_directory
              << " " << common_opts.events << std::endl;

    // Output stats
    uint64_t n_spacepoints = 0;
    uint64_t n_measurements = 0;

    /// Type declarations
    using host_detector_type = detray::detector<detray::default_metadata,
                                                detray::host_container_types>;

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;

    // Read the detector
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg.add_file(traccc::io::data_directory() + det_opts.detector_file);
    if (!det_opts.material_file.empty()) {
        reader_cfg.add_file(traccc::io::data_directory() +
                            det_opts.material_file);
    }
    if (!det_opts.grid_file.empty()) {
        reader_cfg.add_file(traccc::io::data_directory() + det_opts.grid_file);
    }
    const auto [host_det, names] =
        detray::io::read_detector<host_detector_type>(host_mr, reader_cfg);

    // Loop over events
    for (unsigned int event = common_opts.skip;
         event < common_opts.events + common_opts.skip; ++event) {

        // Read measurements
        traccc::io::measurement_reader_output meas_read_out(&host_mr);
        traccc::io::read_measurements(meas_read_out, event,
                                      common_opts.input_directory,
                                      traccc::data_format::csv);
        traccc::measurement_collection_types::host& measurements_per_event =
            meas_read_out.measurements;
        n_measurements += measurements_per_event.size();

        // Run spacepoint formation
        experimental::spacepoint_formation<decltype(host_det)> sp_formation(
            host_mr);

        auto spacepoints = sp_formation(host_det, measurements_per_event);
        n_spacepoints += spacepoints.size();
    }

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read    " << n_measurements << " measurements" << std::endl;
    std::cout << "- created " << n_spacepoints << " spacepoints" << std::endl;
}

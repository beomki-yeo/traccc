/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/options/common_options.hpp"
#include "traccc/options/handle_argument_errors.hpp"
#include "traccc/options/options.hpp"
#include "traccc/options/particle_gen_options.hpp"
#include "traccc/options/propagation_options.hpp"
#include "traccc/options/simulation_options.hpp"
#include "traccc/utils/digitization_writer.hpp"

// detray include(s).
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/detector_writer.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/simulation/simulator.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// Boost include(s).
#include <boost/filesystem.hpp>

using namespace traccc;
namespace po = boost::program_options;

int simulate(std::string output_directory, unsigned int events,
             const traccc::particle_gen_options<scalar>& pg_opts,
             const traccc::propagation_options<scalar>& propagation_opts,
             const traccc::simulation_options& simulation_opts) {

    // Use deterministic random number generator for testing
    using uniform_gen_t =
        detray::random_numbers<scalar, std::uniform_real_distribution<scalar>,
                               std::seed_seq>;

    // Memory resource
    vecmem::host_memory_resource host_mr;

    /*****************************
     * Build a toy geometry
     *****************************/

    // Detector type
    using detector_type = detray::detector<detray::toy_metadata<>>;

    // B field value and its type
    // @TODO: Set B field as argument
    const vector3 B{0, 0, 2 * detray::unit<scalar>::T};
    using field_type = typename detector_type::bfield_type;

    // Create the toy geometry
    const auto [det, name_map] =
        detray::create_toy_geometry<detray::host_container_types>(
            host_mr,
            field_type(
                field_type::backend_t::configuration_t{B[0], B[1], B[2]}),
            4u, 7u);

    /***************************
     * Generate simulation data
     ***************************/

    // Origin of particles
    auto generator =
        detray::random_track_generator<traccc::free_track_parameters,
                                       uniform_gen_t>(
            pg_opts.gen_nparticles, pg_opts.vertex, pg_opts.vertex_stddev,
            pg_opts.mom_range, pg_opts.theta_range, pg_opts.phi_range);
    using generator_type = decltype(generator);

    // Smearing value for measurements
    detray::measurement_smearer<transform3> meas_smearer(
        50 * detray::unit<scalar>::um, 50 * detray::unit<scalar>::um);

    // Create directory
    const std::string full_path = io::data_directory() + output_directory;
    boost::filesystem::create_directories(full_path);

    if (!simulation_opts.run_digitization) {
        // Writer config
        using writer_type =
            detray::smearing_writer<detray::measurement_smearer<transform3>>;
        typename writer_type::config smearer_writer_cfg{meas_smearer};

        // Run simulator
        auto sim =
            detray::simulator<detector_type, generator_type, writer_type>(
                events, det, std::move(generator),
                std::move(smearer_writer_cfg), full_path);
        sim.get_config().step_constraint = propagation_opts.step_constraint;
        sim.get_config().overstep_tolerance =
            propagation_opts.overstep_tolerance;
        sim.run();
    } else {
        // Writer config
        using writer_type = traccc::digitization_writer;

        // Read the digitization configuration file
        auto digi_map = traccc::io::experimental::read_digitization_config(
            simulation_opts.digitization_config_file);

        typename writer_type::config digitization_writer_cfg{digi_map};

        // Run simulator
        auto sim =
            detray::simulator<detector_type, generator_type, writer_type>(
                events, det, std::move(generator),
                std::move(digitization_writer_cfg), full_path);
        sim.get_config().step_constraint = propagation_opts.step_constraint;
        sim.get_config().overstep_tolerance =
            propagation_opts.overstep_tolerance;
        sim.run();
    }

    // Create detector file
    auto writer_cfg = detray::io::detector_writer_config{}
                          .format(detray::io::format::json)
                          .replace_files(true);
    detray::io::write_detector(det, name_map, writer_cfg);

    return 1;
}

// The main routine
//
int main(int argc, char* argv[]) {
    // Set up the program options
    po::options_description desc("Allowed options");

    // Add options
    desc.add_options()("help,h", "Give some help with the program's options");
    desc.add_options()("output_directory", po::value<std::string>()->required(),
                       "specify the directory of output data");
    desc.add_options()("events", po::value<unsigned int>()->required(),
                       "number of events");
    traccc::particle_gen_options<scalar> pg_opts(desc);
    traccc::propagation_options<scalar> propagation_opts(desc);
    traccc::simulation_options simulation_opts(desc);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // Check errors
    traccc::handle_argument_errors(vm, desc);

    // Read options
    auto output_directory = vm["output_directory"].as<std::string>();
    auto events = vm["events"].as<unsigned int>();
    pg_opts.read(vm);
    propagation_opts.read(vm);
    simulation_opts.read(vm);

    std::cout << "Running " << argv[0] << " " << output_directory << " "
              << events << std::endl;

    return simulate(output_directory, events, pg_opts, propagation_opts,
                    simulation_opts);
}

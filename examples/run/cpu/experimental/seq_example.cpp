/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// io
#include "traccc/io/experimental/event_map.hpp"
#include "traccc/io/experimental/read_cells.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/read_geometry.hpp"
#include "traccc/io/utils.hpp"

// algorithms
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/finding/finding_algorithm.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/seeding/experimental/spacepoint_formation.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

// performance
#include "traccc/efficiency/finding_performance_writer.hpp"
#include "traccc/efficiency/seeding_performance_writer.hpp"
#include "traccc/resolution/fitting_performance_writer.hpp"

// options
#include "traccc/options/common_options.hpp"
#include "traccc/options/finding_input_options.hpp"
#include "traccc/options/full_tracking_input_options.hpp"
#include "traccc/options/handle_argument_errors.hpp"
#include "traccc/options/propagation_options.hpp"
#include "traccc/options/seeding_input_options.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/io/common/detector_reader.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <exception>
#include <iostream>

namespace po = boost::program_options;

// The main routine
//
int main(int argc, char* argv[]) {
    // Set up the program options
    po::options_description desc("Allowed options");

    // Add options
    desc.add_options()("help,h", "Give some help with the program's options");
    traccc::common_options common_opts(desc);
    traccc::full_tracking_input_config full_tracking_input_cfg(desc);
    traccc::seeding_input_config seeding_input_cfg(desc);
    traccc::finding_input_config finding_input_cfg(desc);
    traccc::propagation_options<traccc::scalar> propagation_opts(desc);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // Check errors
    traccc::handle_argument_errors(vm, desc);

    // Read options
    common_opts.read(vm);
    seeding_input_cfg.read(vm);
    finding_input_cfg.read(vm);
    full_tracking_input_cfg.read(vm);
    propagation_opts.read(vm);

    // Output stats
    uint64_t n_cells = 0;
    uint64_t n_measurements = 0;
    uint64_t n_spacepoints = 0;
    uint64_t n_seeds = 0;
    uint64_t n_found_tracks = 0;
    uint64_t n_fitted_tracks = 0;

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    /// Type declarations
    using host_detector_type = detray::detector<detray::toy_metadata<>>;

    using b_field_t = typename host_detector_type::bfield_type;
    using rk_stepper_type =
        detray::rk_stepper<b_field_t::view_t,
                           typename host_detector_type::transform3,
                           detray::constrained_step<>>;
    using host_navigator_type = detray::navigator<const host_detector_type>;
    using host_fitter_type =
        traccc::kalman_fitter<rk_stepper_type, host_navigator_type>;

    // B field value and its type
    // @TODO: Set B field as argument
    const traccc::vector3 B{0, 0, 2 * detray::unit<traccc::scalar>::T};

    // Read the detector
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg
        .add_file(traccc::io::data_directory() + common_opts.detector_file)
        .add_file(traccc::io::data_directory() + common_opts.material_file)
        .bfield_vec(B[0], B[1], B[2]);

    auto [host_det, names] =
        detray::io::read_detector<host_detector_type>(host_mr, reader_cfg);

    // Read the digitization configuration file
    auto digi_cfg = traccc::io::experimental::read_digitization_config(
        full_tracking_input_cfg.digitization_config_file);

    // Performance writer
    traccc::seeding_performance_writer sd_performance_writer(
        traccc::seeding_performance_writer::config{});
    traccc::finding_performance_writer find_performance_writer(
        traccc::finding_performance_writer::config{});
    traccc::fitting_performance_writer fit_performance_writer(
        traccc::fitting_performance_writer::config{});

    // Algorithms
    traccc::clusterization_algorithm ca(host_mr);

    traccc::experimental::spacepoint_formation<host_detector_type> sf(host_mr);

    traccc::seedfinder_config finder_config;
    finder_config.sigmaScattering = 10.f;
    finder_config.maxPtScattering = 1.f * traccc::unit<traccc::scalar>::GeV;
    traccc::spacepoint_grid_config grid_config(finder_config);
    traccc::seedfilter_config filter_config;
    traccc::seeding_algorithm sa(finder_config, grid_config, filter_config,
                                 host_mr);

    traccc::track_params_estimation tp(host_mr);

    typename traccc::finding_algorithm<rk_stepper_type,
                                       host_navigator_type>::config_type cfg;
    cfg.min_track_candidates_per_track =
        finding_input_cfg.track_candidates_range[0];
    cfg.max_track_candidates_per_track =
        finding_input_cfg.track_candidates_range[1];
    cfg.constrained_step_size = propagation_opts.step_constraint;
    traccc::finding_algorithm<rk_stepper_type, host_navigator_type>
        host_finding(cfg);

    typename traccc::fitting_algorithm<host_fitter_type>::config_type fit_cfg;
    fit_cfg.step_constraint = propagation_opts.step_constraint;
    traccc::fitting_algorithm<host_fitter_type> host_fitting(fit_cfg);

    // Loop over events
    for (unsigned int event = common_opts.skip;
         event < common_opts.events + common_opts.skip; ++event) {

        traccc::io::cell_reader_output readOut(&host_mr);

        // Read the cells from the relevant event file
        traccc::io::experimental::read_cells(
            readOut, event, common_opts.input_directory,
            common_opts.input_data_format, &digi_cfg);
        traccc::cell_collection_types::host& cells_per_event = readOut.cells;
        traccc::cell_module_collection_types::host& modules_per_event =
            readOut.modules;

        /*-------------------
            Clusterization
          -------------------*/

        auto measurements_per_event = ca(cells_per_event, modules_per_event);

        /*------------------------
            Spacepoint formation
          ------------------------*/

        auto spacepoints_per_event = sf(host_det, measurements_per_event);

        /*----------------
             Seeding
          ---------------*/

        auto seeds = sa(spacepoints_per_event);

        /*----------------------------
           Track Parameter Estimation
          ----------------------------*/

        auto params = tp(spacepoints_per_event, seeds,
                         {0.f, 0.f, finder_config.bFieldInZ});

        /*------------------------
           Track Finding with CKF
          ------------------------*/

        auto track_candidates =
            host_finding(host_det, std::move(measurements_per_event), params);

        /*------------------------
           Track Fitting with KF
          ------------------------*/

        auto track_states = host_fitting(host_det, track_candidates);

        /*----------------------------
          Statistics
          ----------------------------*/

        n_cells += cells_per_event.size();
        n_measurements += measurements_per_event.size();
        n_spacepoints += spacepoints_per_event.size();
        n_seeds += seeds.size();
        n_found_tracks += track_candidates.size();
        n_fitted_tracks += track_states.size();

        /*------------
          Writer
          ------------*/

        if (common_opts.check_performance) {

            traccc::io::experimental::event_map evt_map(
                event, host_mr,
                full_tracking_input_cfg.digitization_config_file,
                common_opts.input_directory, common_opts.input_directory,
                common_opts.input_directory);
            sd_performance_writer.write(vecmem::get_data(seeds),
                                        vecmem::get_data(spacepoints_per_event),
                                        evt_map);

            /*
            find_performance_writer.write(traccc::get_data(track_candidates),
                                          evt_map);

            for (unsigned int i = 0; i < track_states.size(); i++) {
                const auto& trk_states_per_track = track_states.at(i).items;

                const auto& fit_info = track_states[i].header;

                fit_performance_writer.write(trk_states_per_track, fit_info,
                                             host_det, evt_map);
            }
            */
        }
    }

    if (common_opts.check_performance) {
        sd_performance_writer.finalize();
        find_performance_writer.finalize();
        fit_performance_writer.finalize();
    }

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read    " << n_cells << " cells" << std::endl;
    std::cout << "- created (cpu)  " << n_measurements << " measurements"
              << std::endl;
    std::cout << "- created (cpu)  " << n_spacepoints << " spacepoints"
              << std::endl;
    std::cout << "- created (cpu)  " << n_seeds << " seeds" << std::endl;
    std::cout << "- created (cpu)  " << n_found_tracks << " found tracks"
              << std::endl;
    std::cout << "- created (cpu)  " << n_fitted_tracks << " fitted tracks"
              << std::endl;

    return 0;
}

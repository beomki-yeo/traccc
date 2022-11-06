/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "tests/parameter_smearer.hpp"
#include "traccc/edm/track_state.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/io/reader.hpp"
#include "traccc/resolution/fitting_performance_writer.hpp"

// detray include(s).
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/propagator.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"
#include "tests/common/tools/simulator.hpp"
#include "tests/common/tools/track_generators.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace traccc;
using matrix_operator = typename transform3::matrix_actor;

// This defines the local frame test suite
TEST(kalman_filter, telescope_truth_tracking) {

    // Memory resource
    vecmem::host_memory_resource host_mr;

    /*****************************
     *  Setup values
     *****************************/

    // Numer of events
    std::size_t n_events = 10000;

    // Standard deviations for seed track parameters
    std::array<scalar, e_bound_size> stddevs = {
        0.03 * detray::unit_constants::mm,
        0.03 * detray::unit_constants::mm,
        0.017,
        0.017,
        0.001 / detray::unit_constants::GeV,
        1 * detray::unit_constants::ns};

    // Material and thickness for detector planes
    const auto mat = detray::silicon_tml<scalar>();
    const scalar thickness = 0.5 * detray::unit_constants::mm;

    /*****************************
     * Build a telescope geometry
     *****************************/

    // Plane alignment direction
    detray::detail::ray<transform3> traj{{0, 0, 0}, 0, {1, 0, 0}, -1};
    // Position of planes (in mm unit)
    std::vector<scalar> plane_positions = {-10., 20., 40., 60.,  80., 100.,
                                           120., 140, 160, 180., 200.};

    // Detector type
    using detector_type =
        detray::detector<detray::detector_registry::telescope_detector,
                         covfie::field>;

    // B field value and its type
    const vector3 B{2, 0, 0 * detray::unit_constants::T};
    using b_field_t = typename detector_type::bfield_type;

    // Create the detector
    const detector_type det = create_telescope_detector(
        host_mr,
        b_field_t(b_field_t::backend_t::configuration_t{B[0], B[1], B[2]}),
        plane_positions, traj, 100000. * detray::unit_constants::mm,
        100000. * detray::unit_constants::mm, mat, thickness);

    /***************************
     * Generate simulation data
     ***************************/

    // Basic setups
    // 1. one track per event
    // 2. vertex: (0,0,0)
    // 3. momentum: 1 GeV/c
    constexpr unsigned int theta_steps{1};
    constexpr unsigned int phi_steps{1};
    const vector3 x_0{0, 0, 0};
    const scalar mom_0 = 1 * detray::unit_constants::GeV;

    // Track direction: {theta = M_PI/2, phi = 0}
    auto generator =
        detray::uniform_track_generator<traccc::free_track_parameters>(
            theta_steps, phi_steps, x_0, mom_0, {M_PI / 2., M_PI / 2.},
            {0., 0.});

    // Smearing value for measurements
    detray::measurement_smearer<scalar> meas_smearer(
        50 * detray::unit_constants::um, 50 * detray::unit_constants::um);

    auto sim = detray::simulator(n_events, data_directory(), det, generator,
                                 meas_smearer);
    sim.run();

    /***************************
     * Prepare track candidates
     ***************************/

    // Navigator, stepper, and fitter
    using navigator_type = detray::navigator<decltype(det)>;
    using rk_stepper_type = detray::rk_stepper<b_field_t::view_t, transform3>;
    using fitter_type = kalman_fitter<rk_stepper_type, navigator_type>;
    // using seed_parameter_type = free_track_parameters;
    using seed_parameter_type = bound_track_parameters;

    seed_generator<rk_stepper_type, navigator_type> sg(det);

    // Track candidates for multiple tracks
    traccc::track_candidates_container_types<seed_parameter_type>::host
        track_candidates(&host_mr);

    for (std::size_t event = 0; event < n_events; event++) {

        // Read the measurements from the relevant event file
        traccc::measurement_container_types::host measurements_per_event =
            traccc::read_measurements_from_event(event, "", data_format::csv,
                                                 host_mr);

        // Read the vertex information
        const auto io_particles_file =
            data_directory() + get_event_filename(0, "-particles.csv");
        particle_reader preader(io_particles_file);

        csv_particle io_particle;
        preader.read(io_particle);

        point3 pos{io_particle.vx, io_particle.vy, io_particle.vz};
        vector3 mom{io_particle.px, io_particle.py, io_particle.pz};

        // Make a seed parameter
        // @TODO: Smear the seed
        free_track_parameters vertex(pos, io_particle.vt, mom, io_particle.q);
        auto seed_params = sg(vertex, stddevs);
        // auto seed_params = parameter_smearer()(vertex, stddevs);

        // Make a track candidate vector
        vecmem::vector<track_candidate> candidates(&host_mr);

        // The number of measurement should be eqaul to the number of physical
        // planes (add -2 because the first and last surface are portals)
        const std::size_t n_headers = measurements_per_event.size();
        ASSERT_EQ(n_headers, plane_positions.size() - 2);

        // The number of measurements per plane should be one
        for (std::size_t i = 0; i < n_headers; i++) {
            ASSERT_EQ(measurements_per_event[i].items.size(), 1);

            candidates.push_back({measurements_per_event[i].header.module,
                                  measurements_per_event[i].items[0]});
        }

        track_candidates.push_back(std::move(seed_params),
                                   std::move(candidates));
    }

    /********************
     * Run Track fitting
     ********************/

    fitting_algorithm<fitter_type, seed_parameter_type> fitting(det);

    // Run fitting
    auto track_states = fitting(track_candidates);

    ASSERT_EQ(track_states.size(), n_events);

    // performance writer
    traccc::fitting_performance_writer fit_performance_writer(
        traccc::fitting_performance_writer::config{});

    fit_performance_writer.add_cache("CPU");

    for (std::size_t i = 0; i < n_events; i++) {
        const auto& track_states_per_track = track_states[i].items;

        ASSERT_EQ(track_states_per_track.size(), plane_positions.size() - 2);

        traccc::event_map2<detector_type> evt_map(det, i, "", "", "");
        fit_performance_writer.write("CPU", track_states_per_track, evt_map);
    }

    fit_performance_writer.finalize();
}
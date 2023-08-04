/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"

// detray include(s).
#include "detray/core/detector.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/masks/masks.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>
#include <string_view>
#include <vector>

namespace traccc {

/// Kalman Fitting Test with Telescope Geometry
///
/// Test parameters:
/// (1) name
/// (2) origin
/// (3) origin stddev
/// (4) momentum range
/// (5) eta range
/// (6) phi range
/// (7) number of tracks per event
/// (8) number of events
class KalmanFittingTests
    : public ::testing::TestWithParam<std::tuple<
          std::string, std::array<scalar, 3u>, std::array<scalar, 3u>,
          std::array<scalar, 2u>, std::array<scalar, 2u>,
          std::array<scalar, 2u>, unsigned int, unsigned int>> {

    public:
    /// Type declarations
    using host_detector_type =
        detray::detector<detray::default_metadata, covfie::field,
                         detray::host_container_types>;
    using device_detector_type =
        detray::detector<detray::default_metadata, covfie::field_view,
                         detray::device_container_types>;

    using b_field_t = typename host_detector_type::bfield_type;
    using rk_stepper_type = detray::rk_stepper<b_field_t::view_t, transform3,
                                               detray::constrained_step<>>;
    using host_navigator_type = detray::navigator<const host_detector_type>;
    using host_fitter_type =
        kalman_fitter<rk_stepper_type, host_navigator_type>;
    using device_navigator_type = detray::navigator<const device_detector_type>;
    using device_fitter_type =
        kalman_fitter<rk_stepper_type, device_navigator_type>;

    // Use deterministic random number generator for testing
    using uniform_gen_t =
        detray::random_numbers<scalar, std::uniform_real_distribution<scalar>,
                               std::seed_seq>;

    /// Verify that pull distribtions follow the normal distribution
    ///
    /// @param file_name The name of the file holding the distributions
    /// @param hist_names The names of the histograms to process
    ///
    void pull_value_tests(std::string_view file_name,
                          const std::vector<std::string>& hist_names) const;
};

}  // namespace traccc
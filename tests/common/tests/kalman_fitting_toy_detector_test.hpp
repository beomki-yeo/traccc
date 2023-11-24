/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "kalman_fitting_test.hpp"

// Detray include(s).
#include "detray/detectors/bfield.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/detector_reader.hpp"
#include "detray/io/common/detector_writer.hpp"

namespace traccc {

/// Kalman Finding Test for toy geometry
class KalmanFittingToyDetectorTests : public KalmanFittingTests {

    public:
    /// Number of barrel layers
    static const inline unsigned int n_barrels{4u};

    /// Number of endcap layers
    static const inline unsigned int n_endcaps{7u};

    /// B field value and its type
    static constexpr vector3 B{0, 0, 2 * detray::unit<scalar>::T};

    /*
    /// Step constraint (only for simulation)
    static const inline scalar step_constraint = 2 * detray::unit<scalar>::mm;

    /// Overstep tolerance
    static const inline scalar overstep_tolerance =
        -100.f * detray::unit<scalar>::um;

    // Set mask tolerance to a large value not to miss the surface during KF
    static const inline scalar mask_tolerance = 50.f * detray::unit<scalar>::um;
    */
    /// Measurement smearing parameters
    static constexpr std::array<scalar, 2u> smearing{
        50.f * detray::unit<scalar>::um, 50.f * detray::unit<scalar>::um};

    /// Standard deviations for seed track parameters
    static constexpr std::array<scalar, e_bound_size> stddevs = {
        0.01 * detray::unit<scalar>::mm,
        0.01 * detray::unit<scalar>::mm,
        0.001,
        0.001,
        0.01 / detray::unit<scalar>::GeV,
        0.01 * detray::unit<scalar>::ns};

    protected:
    virtual void SetUp() override {
        vecmem::host_memory_resource host_mr;

        // Create the toy geometry
        auto [det, name_map] =
            detray::create_toy_geometry(host_mr, {n_barrels, n_endcaps});

        // Write detector file
        auto writer_cfg = detray::io::detector_writer_config{}
                              .format(detray::io::format::json)
                              .replace_files(true)
                              .path(std::get<0>(GetParam()));
        detray::io::write_detector(det, name_map, writer_cfg);
    }
};

}  // namespace traccc
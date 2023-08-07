/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "kalman_fitting_test.hpp"

namespace traccc {

/// Combinatorial Kalman Finding Test with Sparse tracks
class KalmanFittingWireChamberTests : public KalmanFittingTests {

    public:
    /// Number of layers
    static const inline unsigned int n_wire_layers{20u};

    /// B field value and its type
    static constexpr vector3 B{0, 0, 2 * detray::unit<scalar>::T};

    /// Measurement smearing parameters
    static constexpr std::array<scalar, 2u> smearing{
        50 * detray::unit<scalar>::um, 50 * detray::unit<scalar>::um};

    /// Standard deviations for seed track parameters
    static constexpr std::array<scalar, e_bound_size> stddevs = {
        0.03 * detray::unit<scalar>::mm,
        0.03 * detray::unit<scalar>::mm,
        0.017,
        0.017,
        0.001 / detray::unit<scalar>::GeV,
        1 * detray::unit<scalar>::ns};
};

}  // namespace traccc
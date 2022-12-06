/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray include(s).
#include "detray/detectors/create_telescope_detector.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// ROOT include(s).
#include <TF1.h>

// Kalman Fitting Test with Telescope Geometry
//
// Tuple parameter made of (1) initial particle momentum and (2) initial phi
class KalmanFittingTests
    : public ::testing::TestWithParam<std::tuple<scalar, scalar>> {

    using detector_type =
        detray::detector<detray::detecgtor_registry::telescope_detector,
                         covfie::field>;

    void pull_value_test(TFile& f) {

        std::array<std::string, 5> pull_names{"pull_d0", "pull_z0", "pull_phi",
                                              "pull_theta", "pull_qop"};
        TF1 gaus{"gaus", "gaus", -5, 5};
        double fit_par[3];

        for (auto name : pull_names) {

            auto pull_d0 = (TH1F*)f->Get(name.c_str());

            // Set the mean seed to 0
            gaus.SetParameters(1, 0.);
            gaus.SetParLimits(1, -1., 1.);
            // Set the standard deviation seed to 1
            gaus.SetParameters(2, 1.0);
            gaus.SetParLimits(2, 0.5, 2.);

            auto res = pull_d0->Fit("gaus", "Q0S");

            gaus.GetParameters(&fit_par[0]);

            // Mean check
            EXPECT_NEAR(fit_par[1], 0, 0.05) << name << " mean value error";

            // Sigma check
            EXPECT_NEAR(fit_par[2], 1, 0.1) << name << " sigma value error";
        }
    }
};

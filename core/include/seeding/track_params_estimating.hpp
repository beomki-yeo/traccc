/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <edm/internal_spacepoint.hpp>
#include <edm/seed.hpp>
#include <edm/track_parameters.hpp>
#include <iostream>
#include <seeding/detail/seeding_config.hpp>
#include <seeding/detail/statistics.hpp>
#include <seeding/doublet_finding.hpp>
#include <seeding/seed_filtering.hpp>
#include <seeding/triplet_finding.hpp>

#include "Utilities/Definitions.hpp"

namespace traccc {

/// track parameters estimation
/// Algorithms from Acts/Seeding/EstimateTrackParamsFromSeed.hpp

struct track_params_estimating {
    track_params_estimating() {}

    host_track_parameters_container operator()(const host_measurement_container& measurements,
					       const host_seed_container& seed_container,
					       const Acts::GeometryContext& gctx){
        host_track_parameters_container track_parameters_container(
            {host_track_parameters_container::header_vector(1, 0),
             host_track_parameters_container::item_vector(1)});
	this->operator()(measurements, seed_container, gctx, track_parameters_container);

	return track_parameters_container;
    }
    
    /// Callable operator for track parameter estimation
    /// naively assumed that bfield is 2 T...    
    void operator()(const host_measurement_container& measurements,
                    const host_seed_container& seed_container,
		    const Acts::GeometryContext& gctx,
		    host_track_parameters_container& track_parameters_container,
                    vector3 bfield = {0, 0, 2},
                    scalar mass = 139.57018 * Acts::UnitConstants::MeV) {
	
        auto& seeds = seed_container.items[0];

        for (auto seed : seeds) {
            array<vector3, 3> sp_global_positions;

            const auto& spB = seed.spB;
            const auto& spM = seed.spM;
            const auto& spT = seed.spT;

            sp_global_positions[0] = spB.global_position();
            sp_global_positions[1] = spM.global_position();
            sp_global_positions[2] = spT.global_position();

            // Define a new coordinate frame with its origin at the bottom space
            // point, z axis long the magnetic field direction and y axis
            // perpendicular to vector from the bottom to middle space point.
            // Hence, the projection of the middle space point on the tranverse
            // plane will be located at the x axis of the new frame.
            vector3 relVec = sp_global_positions[1] - sp_global_positions[0];
            vector3 newZAxis = vector::normalize(bfield);
            vector3 newYAxis =
                vector::normalize(vector::cross(newZAxis, relVec));
            vector3 newXAxis = vector::cross(newYAxis, newZAxis);

            // The center of the new frame is at the bottom space point
            vector3 translation = sp_global_positions[0];

            transform3 trans(translation, newZAxis, newXAxis);

	    
            // The coordinate of the middle and top space point in the new frame
            auto local1 = trans.point_to_local(sp_global_positions[1]);
            auto local2 = trans.point_to_local(sp_global_positions[2]);

            // std::cout << local1[0] << "  " << local1[1] << "  "  << local2[0]
            // << "  " << local2[1] << std::endl;

            // Lambda to transform the coordinates to the (u, v) space
            auto uvTransform = [](const vector3& local) -> vector2 {
                vector2 uv;
                scalar denominator = local[0] * local[0] + local[1] * local[1];
                uv[0] = local[0] / denominator;
                uv[1] = local[1] / denominator;
                return uv;
            };

            // The uv1.y() should be zero
            vector2 uv1 = uvTransform(local1);
            vector2 uv2 = uvTransform(local2);

            // A,B are slope and intercept of the straight line in the u,v plane
            // connecting the three points
            scalar A = (uv2[1] - uv1[1]) / (uv2[0] - uv1[0]);
            scalar B = uv2[1] - A * uv2[0];

            // Curvature (with a sign) estimate
            scalar rho = -2.0 * B / std::hypot(1., A);
            // The projection of the top space point on the transverse plane of
            // the new frame
            scalar rn = local2[0] * local2[0] + local2[1] * local2[1];
            // The (1/tanTheta) of momentum in the new frame,
            scalar invTanTheta =
                local2[2] * std::sqrt(1. / rn) / (1. + rho * rho * rn);

            // The momentum direction in the new frame (the center of the circle
            // has the coordinate (-1.*A/(2*B), 1./(2*B)))
            vector3 transDirection(1., A, std::hypot(1, A) * invTanTheta);
            // Transform it back to the original frame
            vector3 direction = transform3::rotate(
                trans._data, vector::normalize(transDirection));

	    Acts::BoundVector params;

            // The estimated phi and theta
            params[Acts::eBoundPhi] = getter::phi(direction);
            params[Acts::eBoundTheta] = getter::theta(direction);

            // The measured loc0 and loc1
            const auto& meas_for_spB =
                measurements.items[spB.m_idx[0]][spB.m_idx[1]];
            params[Acts::eBoundLoc0] = meas_for_spB.local[0];
            params[Acts::eBoundLoc1] = meas_for_spB.local[1];

            // The estimated q/pt in [GeV/c]^-1 (note that the pt is the
            // projection of momentum on the transverse plane of the new frame)
            scalar qOverPt =
                rho * (Acts::UnitConstants::m) / (0.3 * getter::norm(bfield));
            // The estimated q/p in [GeV/c]^-1
            params[Acts::eBoundQOverP] = qOverPt / std::hypot(1., invTanTheta);

            // The estimated momentum, and its projection along the magnetic
            // field diretion
            scalar pInGeV = std::abs(1.0 / params[Acts::eBoundQOverP]);
            scalar pzInGeV = 1.0 / std::abs(qOverPt) * invTanTheta;
            scalar massInGeV = mass / Acts::UnitConstants::GeV;

            // The estimated velocity, and its projection along the magnetic
            // field diretion
            scalar v = pInGeV / std::hypot(pInGeV, massInGeV);
            scalar vz = pzInGeV / std::hypot(pInGeV, massInGeV);
            // The z coordinate of the bottom space point along the magnetic
            // field direction
            scalar pathz = vector::dot(sp_global_positions[0], bfield) /
                           getter::norm(bfield);

            // The estimated time (use path length along magnetic field only if
            // it's not zero)
            if (pathz != 0) {
                params[Acts::eBoundTime] = pathz / vz;
            } else {
                params[Acts::eBoundTime] =
                    getter::norm(sp_global_positions[0]) / v;
            }
	    
	    Acts::BoundMatrix cov = Acts::BoundSymMatrix::Zero();
	    cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = meas_for_spB.variance[0];
	    cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = meas_for_spB.variance[1];

	    // Note: cov params for theta,phi and p should be investigated in detail
	    cov(Acts::eBoundPhi, Acts::eBoundPhi) = 1*Acts::UnitConstants::degree;
	    cov(Acts::eBoundTheta, Acts::eBoundTheta) = 1*Acts::UnitConstants::degree;
	    ActsScalar p = 1/params[Acts::eBoundQOverP];
	    ActsScalar sigmaP = 0.001 * p;
	    ActsScalar sigmaQoP = sigmaP / (p*p);
	    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = sigmaQoP * sigmaQoP;

	    bound_parameters bd_params(gctx, cov, params, nullptr);
        }
    }
};

}  // namespace traccc

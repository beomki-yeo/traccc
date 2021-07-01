/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <edm/internal_spacepoint.hpp>
#include <edm/seed.hpp>
#include <iostream>
#include <seeding/detail/seeding_config.hpp>
#include <seeding/detail/statistics.hpp>
#include <seeding/doublet_finding.hpp>
#include <seeding/seed_filtering.hpp>
#include <seeding/triplet_finding.hpp>

namespace traccc{

/// track parameters estimation
/// Originated from Acts/Seeding/EstimateTrackParamsFromSeed.hpp

    
struct track_params_estimating {
    
    /// Callable operator for track parameter estimation
    /// naively assumed that bfield is 2 T...
    void operator()(host_seed_container& seed_container, vector3 bfield = {0,0,2}) {
	auto& seeds = seed_container.items[0];

	for (auto seed: seeds){
	    array<vector3, 3> sp_global_positions;
	    sp_global_positions[0] = seed.spB.global_position();
	    sp_global_positions[1] = seed.spM.global_position();
	    sp_global_positions[2] = seed.spT.global_position();
	    

	    // Define a new coordinate frame with its origin at the bottom space point, z
	    // axis long the magnetic field direction and y axis perpendicular to vector
	    // from the bottom to middle space point. Hence, the projection of the middle
	    // space point on the tranverse plane will be located at the x axis of the new
	    // frame.
	    vector3 relVec = sp_global_positions[1] - sp_global_positions[0];
	    vector3 newZAxis = vector::normalize(bfield);
	    vector3 newYAxis = vector::normalize(vector::cross(newZAxis,relVec));
	    vector3 newXAxis = vector::cross(newYAxis, newZAxis);
	    
	    Acts::RotationMatrix3 rotation;
	    rotation.col(0) = std::array<scalar, 3>(newXAxis);
	    //rotation.col(1) = newYAxis;
	    //rotation.col(2) = newZAxis;
	    
	}
    }


    

};

} // traccc

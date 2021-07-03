/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "edm/measurement.hpp"
#include "edm/spacepoint.hpp"

namespace traccc {

/// Connected component labeling.
struct spacepoint_formation {

    /// Callable operator for the space point formation
    ///
    /// @param measurements are the input measurements for an event
    ///
    /// C++20 piping interface
    ///
    /// @return a spacepoint container
    host_spacepoint_container operator()(const host_measurement_container& measurements) const {
	host_spacepoint_container spacepoints;
        this->operator()(measurements, spacepoints);
        return spacepoints;
	
    }

    /// Callable operator for the space point formation, based on one single
    /// module
    ///
    /// @param measurements are the input measurements for an event
    ///
    /// void interface
    ///
    /// @return a spacepoint container
    void operator()(const host_measurement_container& measurements,
                    host_spacepoint_container& spacepoints) const {
        // Run the algorithm
        spacepoints.headers.reserve(measurements.headers.size());
	spacepoints.items.reserve(measurements.items.size());
	
	for (size_t i_h = 0; i_h < measurements.headers.size(); i_h++){

	    const auto& module = measurements.headers[i_h];
	    const auto& measurements_per_module = measurements.items[i_h];

	    host_spacepoint_collection spacepoints_per_module;
	    spacepoints_per_module.reserve(measurements_per_module.size());
	    
	    for (size_t i_m = 0; i_m < measurements_per_module.size(); ++i_m) {
		const auto& m = measurements_per_module[i_m];
		spacepoint s;
		point3 local_3d = {m.local[0], m.local[1], 0.};
		s.global = module.placement.point_to_global(local_3d);
		s.m_idx = spacepoint::measurement_index({i_h,i_m});
		// @todo add variance estimation
		spacepoints_per_module.push_back(std::move(s));
	    }

	    spacepoints.headers.push_back(module.module);
	    spacepoints.items.push_back(spacepoints_per_module);
	}       
    }
};

}  // namespace traccc

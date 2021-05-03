/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "edm/spacepoint.hpp"
#include "definitions/primitives.hpp"

// spacepoint grid configuration
struct spacepoint_grid_config {
    // magnetic field in kTesla
    float bFieldInZ;
    // minimum pT to be found by seedfinder in MeV
    float minPt;
    // maximum extension of sensitive detector layer relevant for seeding as
    // distance from x=y=0 (i.e. in r) in mm
    float rMax;
    // maximum extension of sensitive detector layer relevant for seeding in
    // positive direction in z in mm
    float zMax;
    // maximum extension of sensitive detector layer relevant for seeding in
    // negative direction in z in mm
    float zMin;
    // maximum distance in r from middle space point to bottom or top spacepoint
    // in mm                                                
    float deltaRMax;
    // maximum forward direction expressed as cot(theta)
    float cotThetaMax;
};

struct axis {

}

struct spacepoint_grid{
    spacepoint_grid(const spacepoint_grid_config& config){

	// calculate circle intersections of helix and max detector radius
	float minHelixRadius = config.minPt / (300. * config.bFieldInZ);  // in mm
	float maxR2 = config.rMax * config.rMax;
	float xOuter = maxR2 / (2 * minHelixRadius);
	float yOuter = std::sqrt(maxR2 - xOuter * xOuter);
	float outerAngle = std::atan(xOuter / yOuter);
	// intersection of helix and max detector radius minus maximum R distance from
	// middle SP to top SP
	float innerAngle = 0;
	if (config.rMax > config.deltaRMax) {
	    float innerCircleR2 =
		(config.rMax - config.deltaRMax) * (config.rMax - config.deltaRMax);
	    float xInner = innerCircleR2 / (2 * minHelixRadius);
	    float yInner = std::sqrt(innerCircleR2 - xInner * xInner);
	    innerAngle = std::atan(xInner / yInner);
	}

	// FIXME: phibin size must include max impact parameters
	// divide 2pi by angle delta to get number of phi-bins
	// size is always 2pi even for regions of interest
	int phiBins = std::floor(2 * M_PI / (outerAngle - innerAngle));
	Acts::detail::Axis<detail::AxisType::Equidistant,
			   detail::AxisBoundaryType::Closed>
	    phiAxis(-M_PI, M_PI, phiBins);
	
	// TODO: can probably be optimized using smaller z bins
	// and returning (multiple) neighbors only in one z-direction for forward
	// seeds
	// FIXME: zBinSize must include scattering
	
	float zBinSize = config.cotThetaMax * config.deltaRMax;
	int zBins = std::floor((config.zMax - config.zMin) / zBinSize);
	detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound>
	    zAxis(config.zMin, config.zMax, zBins);
       	
    }
    
    vecmem::jagged_vector< spacepoint > binned_sp;
}

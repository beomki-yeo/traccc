/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"


namespace traccc{

class track_finding_algorithm {

public:

    track_finding_algorithm(host_measurement_container& measurements,
			    host_seed_container& seeds):
	m_measurements(measurements),
	m_seeds(seeds)
    {}

    void operator()(){

    }
    
    // acts_cpu
    using TrackFinderOptions =
	Acts::CombinatorialKalmanFilterOptions<ActsExamples::MeasurementCalibrator,
					       Acts::MeasurementSelector>;

    using TrackFinderResult = std::vector<
	Acts::Result<Acts::CombinatorialKalmanFilterResult<ActsExamples::IndexSourceLink>>>;
    using TrackFinderFunction = std::function<TrackFinderResult(
	const ActsExamples::IndexSourceLinkContainer&,
	const ActsExamples::TrackParametersContainer&,
	const TrackFinderOptions&)>;

private:
    host_measurement_container& m_measurements;
    host_seed_container& m_seeds;
};

}

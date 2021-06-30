/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "edm/trajectory.hpp"

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

    track_finding_algorithm(){}
    
    host_trajectory_container operator()(const host_measurement_container& measurement_container){
        host_trajectory_container trajectory_container(
            {host_trajectory_container::header_vector(1, 0),
             host_trajectory_container::item_vector(1)});
        this->operator()(measurement_container, trajectory_container);
	
    }

    void operator()(const host_measurement_container& measurement_container,
		    host_trajectory_container& trajectory_container){


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
};

}

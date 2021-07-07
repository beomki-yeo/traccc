/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// gpuKalmanFilter (Actscore)
#include "Utilities/Definitions.hpp"
#include "Surfaces/Surface.hpp"
#include "Surfaces/PlaneSurface.hpp"
#include "EventData/TrackParameters.hpp"


namespace traccc {

using bound_parameters = Acts::BoundParameters<Acts::PlaneSurface<Acts::InfiniteBounds>>;
    
/// Header: number of tracks
/// Item: bound track parameters per track
    
/// Container of seeds for an event
template <template <typename> class vector_t>
using track_parameters_collection = vector_t<bound_parameters>;

/// Convenience declaration for the track_parameters collection type to use
/// in host code
using host_track_parameters_collection = track_parameters_collection<vecmem::vector>;

/// Convenience declaration for the track_parameters collection type to use
/// in device code
using device_track_parameters_collection = track_parameters_collection<vecmem::device_vector>;

/// Convenience declaration for the track_parameters container type to use in
/// host code
using host_track_parameters_container = host_container<unsigned int, bound_parameters>;

/// Convenience declaration for the track_parameters container type to use in
/// device code
using device_track_parameters_container = device_container<unsigned int, bound_parameters>;

/// Convenience declaration for the track_parameters container data type to
/// use in host code
using track_parameters_container_data = container_data<unsigned int, bound_parameters>;

/// Convenience declaration for the track_parameters container buffer type to
/// use in host code
using track_parameters_container_buffer = container_buffer<unsigned int, bound_parameters>;

/// Convenience declaration for the track_parameters container view type to
/// use in host code
using track_parameters_container_view = container_view<unsigned int, bound_parameters>;

};  // namespace traccc
    

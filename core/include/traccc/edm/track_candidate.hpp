/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/track_parameters.hpp"

namespace traccc {

struct track_candidate {
    geometry_id surface_link;
    measurement meas;
};

/// Declare a track candidates container type
template <typename parameter_t>
using track_candidates_container_types =
    container_types<parameter_t, track_candidate>;

}  // namespace traccc
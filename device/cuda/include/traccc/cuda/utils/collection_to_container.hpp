/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/measurement.hpp"

namespace traccc::cuda {

/// Measurement collection to container
measurement_container_types::buffer measurement_collection_to_container(
    measurement_collection_types::buffer measurement_buffer);

}  // namespace traccc::cuda
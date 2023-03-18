/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/measurement.hpp"

namespace traccc::device {

struct candidate_link {
    thrust::pair<unsigned int, unsigned int> last_link;
    typename measurement_container_types::host::link_type meas_link;
    std::size_t surface_link;
};

}  // namespace traccc::device
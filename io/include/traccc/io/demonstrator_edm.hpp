/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/edm/cell.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/spacepoint.hpp"

namespace traccc {
struct result {
    measurement_container_types::const_view measurements;
    spacepoint_container_types::const_view spacepoints;
};

using demonstrator_input = vecmem::vector<cell_container_types::const_view>;
using demonstrator_result = vecmem::vector<result>;
}  // namespace traccc

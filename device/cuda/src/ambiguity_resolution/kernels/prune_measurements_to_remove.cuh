/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/ambiguity_resolution/device/prune_measurements_to_remove.hpp"

namespace traccc::cuda::kernels {

__global__ void prune_measurements_to_remove(
    device::prune_measurements_to_remove_payload payload);
}

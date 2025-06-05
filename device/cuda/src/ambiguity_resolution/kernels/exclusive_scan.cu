/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "exclusive_scan.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void exclusive_scan(device::exclusive_scan_payload payload) {

    if (*(payload.terminate) == 1) {
        return;
    }

}

}  // namespace traccc::cuda::kernels
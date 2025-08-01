/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "sort_tracks_per_measurement.cuh"

// VecMem include(s).
#include <vecmem/containers/jagged_device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void sort_tracks_per_measurement(
    device::sort_tracks_per_measurement_payload payload) {

    vecmem::jagged_device_vector<unsigned int> tracks_per_measurement(
        payload.tracks_per_measurement_view);
}
}  // namespace traccc::cuda::kernels

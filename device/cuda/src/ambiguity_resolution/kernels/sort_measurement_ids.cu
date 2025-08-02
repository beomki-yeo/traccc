/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "sort_measurement_ids.cuh"

// VecMem include(s).
#include <vecmem/containers/jagged_device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void sort_measurement_ids(
    device::sort_measurement_ids_payload payload) {

    //__shared__ measurement_id_type sh_meas_ids[1024];

    vecmem::jagged_device_vector<measurement_id_type> meas_ids(payload.meas_ids_view);
}
}  // namespace traccc::cuda::kernels

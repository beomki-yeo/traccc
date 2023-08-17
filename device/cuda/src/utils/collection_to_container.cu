/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/utils/collection_to_container.hpp"

// Vecmem include(s).
#include <vecmem/memory/memory_resource.hpp>

// Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#include <thrust/unique.h>

namespace traccc::cuda {

/// Measurement collection to container
measurement_container_types::buffer measurement_collection_to_container(
    measurement_collection_types::buffer measurement_buffer,
    vecmem::memory_resource& mr) {

    // Create the measurement collection from buffer
    measurement_collection_types::device measurements(buffer);

    // Sort the measurement collection
    thrust::sort(thrust::device, measurements.begin(), measurements.end());

    // Get Copy of uniques
    measurement_collection_types::buffer uniques_buffer{measurements.size(),
                                                        mr};
    measurement_collection_types::device uniques(uniques_buffer);

    thrust::unique_copy(thrust::device, measurements.begin(),
                        measurements.end(), uniques.begin());

    // Get lower bound from unique element
}

}  // namespace traccc::cuda

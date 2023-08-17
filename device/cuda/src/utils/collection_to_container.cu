/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/utils/collection_to_container.hpp"

// Vecmem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>

// Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#include <thrust/unique.h>

namespace traccc::cuda {

/// Measurement collection to container
// measurement_container_types::buffer measurement_collection_to_container(
void measurement_collection_to_container(
    measurement_collection_types::buffer measurement_buffer,
    vecmem::host_memory_resourcce& host_mr;
    vecmem::cuda::device_memory_resource & dev_mr, vecmem::cuda::copy& copy) {

    // Create the measurement collection from buffer
    measurement_collection_types::device measurements(buffer);

    // Sort the measurement collection
    thrust::sort(thrust::device, measurements.begin(), measurements.end());

    // Get Copy of uniques
    measurement_collection_types::buffer uniques_buffer{measurements.size(),
                                                        dev_mr};
    measurement_collection_types::device uniques(uniques_buffer);
    thrust::unique_copy(thrust::device, measurements.begin(),
                        measurements.end(), uniques.begin());

    // Get lower bounds from unique element (Used as an unique element counter)
    vecmem::vector_buffer<unsigned int> lower_bounds_buffer{uniques.size(),
                                                            dev_mr};
    vecmem::device_vector<unsigned int> lower_bounds(lower_bounds_buffer);
    thrust::lower_bound(measurements.begin(), measurements.end(),
                        uniques.begin(), uniques.end(), lower_bounds.begin());

    vecmem::vector<unsigned int> lower_bounds_host{uniques.size(), host_mr};
    auto lower_bounds_host_buffer = vecmem::get_data(lower_bounds_host);

    // Copy device lower bounds to host lower bounds
    copy(lower_bounds_buffer, lower_bounds_host_buffer,
         vecmem::copy::type::device_to_host);

    /*
    thrust::lower_bound(input.begin(), input.end(), values.begin(),
                        values.end(), output.begin(), thrust::less<int>());
    */
}

}  // namespace traccc::cuda

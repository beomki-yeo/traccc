/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// VecMem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// CUB include(s).
#include <cub/cub.cuh>

// GTest include(s).
#include <gtest/gtest.h>

// This defines the local frame test suite

namespace {
vecmem::cuda::copy copy;
vecmem::host_memory_resource host_resource;
vecmem::cuda::device_memory_resource device_resource;

}  // namespace

TEST(CUDAcub, sort) {

    vecmem::vector<unsigned int> host_vector{{3, 2, 1, 8, 4}, &host_resource};

    auto host_buffer = vecmem::get_data(host_vector);
    auto device_buffer = copy.to(vecmem::get_data(host_vector), device_resource,
                                 vecmem::copy::type::host_to_device);

    const auto n_size = static_cast<unsigned int>(host_vector.size());

    vecmem::data::vector_buffer<unsigned int> temp_device_buffer{
        n_size, device_resource};

    void* temp_storage = nullptr;
    std::size_t temp_bytes = 0;

    // Calculate temporary meomory size for sorting
    cub::DeviceRadixSort::SortKeys(nullptr, temp_bytes, device_buffer.ptr(),
                                   temp_device_buffer.ptr(), n_size);

    // Sort
    cudaMalloc(&temp_storage, temp_bytes);
    cub::DeviceRadixSort::SortKeys(&temp_storage, temp_bytes,
                                   device_buffer.ptr(),
                                   temp_device_buffer.ptr(), n_size);

    copy(temp_device_buffer, host_buffer, vecmem::copy::type::device_to_host)
        ->wait();

    ASSERT_EQ(host_vector[0], 1);
    ASSERT_EQ(host_vector[1], 2);
    ASSERT_EQ(host_vector[2], 3);
    ASSERT_EQ(host_vector[3], 4);
    ASSERT_EQ(host_vector[4], 8);
}

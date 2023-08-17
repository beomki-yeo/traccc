/**
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Projection include(s).
#include "traccc/cuda/utils/collection_to_container.hpp"
#include "traccc/edm/measurement.hpp"

// Detray include(s).
#include "detray/geometry/barcode.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// Google test include(s).
#include <gtest/gtest.h>

using namespace traccc;

// This defines the local frame test suite

namespace {
vecmem::cuda::copy copy;
vecmem::host_memory_resource host_resource;
vecmem::cuda::device_memory_resource device_resource;

}  // namespace

TEST(CollectionToContainer, Measurement) {

    traccc::measurement_collection_types::host measurements_host{
        &host_resource};

    measurements_host.push_back(
        {{1.0, 2.0}, {1.5, 2.5}, detray::geometry::barcode{0u}});
    measurements_host.push_back(
        {{9.5, 1.7}, {1.5, 2.5}, detray::geometry::barcode{3u}});
    measurements_host.push_back(
        {{5.2, 7.9}, {1.5, 2.5}, detray::geometry::barcode{1u}});
    measurements_host.push_back(
        {{8.8, 8.2}, {1.5, 2.5}, detray::geometry::barcode{2u}});
    measurements_host.push_back(
        {{1.3, 3.3}, {1.5, 2.5}, detray::geometry::barcode{2u}});
    measurements_host.push_back(
        {{1.9, 3.2}, {1.5, 2.5}, detray::geometry::barcode{1u}});
    measurements_host.push_back(
        {{7.2, 2.5}, {1.5, 2.5}, detray::geometry::barcode{3u}});
    measurements_host.push_back(
        {{6.5, 9.0}, {1.5, 2.5}, detray::geometry::barcode{2u}});
    measurements_host.push_back(
        {{6.6, 1.0}, {1.5, 2.5}, detray::geometry::barcode{2u}});

    ASSERT_EQ(measurements_host.size(), 9u);

}
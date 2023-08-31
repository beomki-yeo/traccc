/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/utils/digitization_algorithm.hpp"

// detray include(s).
#include "detray/detectors/create_toy_geometry.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace traccc;

// This defines the local frame test suite
TEST(digitization, planar_surface) {

    // Memory resource
    vecmem::host_memory_resource host_mr;

    // Detector type
    using detector_type = detray::detector<detray::toy_metadata<>>;

    // B field value and its type
    // @TODO: Set B field as argument
    const vector3 B{0, 0, 2 * detray::unit<scalar>::T};
    using field_type = typename detector_type::bfield_type;

    // Create the toy geometry
    const auto [det, name_map] =
        detray::create_toy_geometry<detray::host_container_types>(
            host_mr,
            field_type(
                field_type::backend_t::configuration_t{B[0], B[1], B[2]}),
            4u, 7u);

    // Read the digitization configuration file
    auto digi_map = traccc::io::experimental::read_digitization_config(
        "../tests/utils/default-input-config-toy.json");

    digitization_algorithm<detector_type> digi(det, digi_map);

    // Create bound track parameter
    traccc::bound_track_parameters bound_param;
    // Volume 1 surface
    bound_param.set_surface_link(detray::geometry::barcode{4503599627379199});
    EXPECT_EQ(bound_param.surface_link().volume(), 1u);

    // Run digitization
    const auto cells = digi(bound_param);

    EXPECT_EQ(cells.size(), 5u);
    EXPECT_EQ(cells[0][0], 480);
    EXPECT_EQ(cells[0][1], 108);
    EXPECT_EQ(cells[1][0], 479);
    EXPECT_EQ(cells[1][1], 108);
    EXPECT_EQ(cells[2][0], 480);
    EXPECT_EQ(cells[2][1], 107);
    EXPECT_EQ(cells[3][0], 481);
    EXPECT_EQ(cells[3][1], 108);
    EXPECT_EQ(cells[4][0], 480);
    EXPECT_EQ(cells[4][1], 109);
}
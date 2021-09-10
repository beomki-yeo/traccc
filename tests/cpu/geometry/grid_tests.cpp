/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// traccc include
#include "edm/spacepoint.hpp"
#include "edm/internal_spacepoint.hpp"
#include "seeding/detail/seeding_config.hpp"

// detray include
#include "tests/common/test_defs.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"
#include "utils/indexing.hpp"

// vecmem include
#include "vecmem/containers/vector.hpp"

#include <gtest/gtest.h>

#include <climits>

// This defines the local frame test suite
TEST(algorithms, cuda_grid) {

    using scalar = detray::scalar;
    using spacepoint_t = traccc::internal_spacepoint<traccc::spacepoint>;
    
    traccc::spacepoint_grid_config grid_config;

    /*---------------
      grid setup
      ---------------*/
    
    // calculate circle intersections of helix and max detector radius
    scalar minHelixRadius =
        grid_config.minPt / (300. * grid_config.bFieldInZ);  // in mm
    scalar maxR2 = grid_config.rMax * grid_config.rMax;
    scalar xOuter = maxR2 / (2 * minHelixRadius);
    scalar yOuter = std::sqrt(maxR2 - xOuter * xOuter);
    scalar outerAngle = std::atan(xOuter / yOuter);
    // intersection of helix and max detector radius minus maximum R distance
    // from middle SP to top SP
    scalar innerAngle = 0;
    if (grid_config.rMax > grid_config.deltaRMax) {
        scalar innerCircleR2 = (grid_config.rMax - grid_config.deltaRMax) *
                               (grid_config.rMax - grid_config.deltaRMax);
        scalar xInner = innerCircleR2 / (2 * minHelixRadius);
        scalar yInner = std::sqrt(innerCircleR2 - xInner * xInner);
        innerAngle = std::atan(xInner / yInner);
    }

    // FIXME: phibin size must include max impact parameters
    // divide 2pi by angle delta to get number of phi-bins
    // size is always 2pi even for regions of interest
    int phiBins = std::floor(2 * M_PI / (outerAngle - innerAngle));


    detray::attach_populator< false, spacepoint_t> replacer;
    detray::serializer2 serializer;
    
    detray::axis::regular<> xaxis{10, -5., 5.};
    detray::axis::regular<> yaxis{10, -5., 5.};
    
    /*
    detray::replace_populator<> replacer;
    detray::serializer2 serializer;

    detray::axis::regular<> xaxis{10, -5., 5.};
    detray::axis::regular<> yaxis{10, -5., 5.};
    using grid2r = detray::grid2<decltype(replacer), decltype(xaxis), decltype(yaxis), decltype(serializer)>;
    
    grid2r g2(std::move(xaxis), std::move(yaxis));

    // Test the initialization
    detray::test::point2 p = {-4.5, -4.5};
    for (unsigned int ib0 = 0; ib0 < 10; ++ib0)
    {
        for (unsigned int ib1 = 0; ib1 < 10; ++ib1)
        {
            p = {static_cast<detray::scalar>(-4.5 + ib0), static_cast<detray::scalar>(-4.5 + ib1)};
            EXPECT_EQ(g2.bin(p), std::numeric_limits<detray::dindex>::max());
        }
    }
    
    p = {-4.5, -4.5};
    // Fill and read
    g2.populate(p, 3u);
    EXPECT_EQ(g2.bin(p), 3u);
    
    // Fill and read two times, fill first 0-99, then 100-199
    for (unsigned int il = 0; il < 2; ++il)
    {
        unsigned int counter = il * 100;
        for (unsigned int ib0 = 0; ib0 < 10; ++ib0)
        {
            for (unsigned int ib1 = 0; ib1 < 10; ++ib1)
            {
                p = {static_cast<detray::scalar>(-4.5 + ib0), static_cast<detray::scalar>(-4.5 + ib1)};
                g2.populate(p, counter);
                EXPECT_EQ(g2.bin(p), counter++);
            }
        }
    }
    */
    /*
    // A zone test w/o neighbour hood
    p = {-4.5, -4.5};
    auto test = g2.zone(p);
    dvector<dindex> expect = {100u};
    EXPECT_EQ(test, expect);

    // A zone test with neighbour hood
    p = {0.5, 0.5};

    darray<dindex, 2> zone11 = {1u, 1u};
    darray<dindex, 2> zone22 = {2u, 2u};

    test = g2.zone(p, {zone11, zone22}, true);
    expect = {143u, 144u, 145u, 146u, 147u, 153u, 154u, 155u, 156u, 157u, 163u, 164u, 165u, 166u, 167u};
    EXPECT_EQ(test, expect);

    axis::circular<> circular{4, -2., 2.};
    axis::regular<> closed{5, 0., 5.};
    using grid2cc = grid2<decltype(replacer), decltype(circular), decltype(closed), decltype(serializer)>;

    grid2cc g2cc(std::move(circular), std::move(closed));
    unsigned int counter = 0;
    for (unsigned icl = 0; icl < 5; ++icl)
    {
        for (unsigned ici = 0; ici < 4; ++ici)
        {
            p = {static_cast<scalar>(-1.5 + ici), static_cast<scalar>(0.5 + icl)};
            g2cc.populate(p, counter++);
        }
    }

    // A zone test for circular testing
    p = {1.5, 2.5};
    test = g2cc.zone(p, {zone11, zone11}, true);
    expect = {4u, 6u, 7u, 8u, 10u, 11u, 12u, 14u, 15u};
    EXPECT_EQ(test, expect);
    */
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

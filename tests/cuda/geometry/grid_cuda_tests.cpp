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
#include "clusterization/clusterization_algorithm.hpp"

// detray include
#include "tests/common/test_defs.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"
#include "utils/indexing.hpp"

// vecmem include
#include "vecmem/containers/vector.hpp"
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// io
#include "io/csv.hpp"
#include "io/reader.hpp"
#include "io/utils.hpp"
#include "io/writer.hpp"

#include <gtest/gtest.h>

#include <climits>

// This defines the local frame test suite
TEST(algorithms, grid_cuda) {

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;
    
    // Read the surface transforms
    auto surface_transforms = traccc::read_geometry("tml_detector/trackml-detector.csv");

    // Read the cells from the relevant event file
    std::string io_cells_file =
	traccc::data_directory() + "tml_pixels/" + "/" +
	traccc::get_event_filename(0, "-cells.csv");
    
    traccc::cell_reader creader(io_cells_file, {"geometry_id", "hit_id", "cannel0", "channel1","activation", "time"});
    traccc::host_cell_container cells_per_event =
	traccc::read_cells(creader, mng_mr, &surface_transforms);
    
    /*-------------------
      Run clusterization
      -------------------*/

    traccc::clusterization_algorithm ca;
    auto ca_result = ca(cells_per_event);
    auto& measurements_per_event = ca_result.first;
    auto& spacepoints_per_event = ca_result.second;
    
    /*---------------
      Axis setup
      ---------------*/

    using scalar = detray::scalar;
    using spacepoint_t = traccc::internal_spacepoint<traccc::spacepoint>;
    
    traccc::spacepoint_grid_config grid_config;
    
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
    int phi_bins = std::floor(2 * M_PI / (outerAngle - innerAngle));
    scalar half_module = 2 * M_PI_2 / 2*phi_bins;
    scalar phi_min = -M_PI + half_module;
    scalar phi_max = M_PI - half_module;
    detray::axis::circular<> phi_axis{phi_bins, phi_min, phi_max};    

    // TODO: can probably be optimized using smaller z bins
    // and returning (multiple) neighbors only in one z-direction for forward
    // seeds
    // FIXME: zBinSize must include scattering

    scalar zBinSize = grid_config.cotThetaMax * grid_config.deltaRMax;
    int zBins = std::floor((grid_config.zMax - grid_config.zMin) / zBinSize);
    detray::axis::regular<> z_axis{zBins, grid_config.zMin, grid_config.zMax};
    
    /*---------------
      Grid setup
      ---------------*/
    
    detray::attach_populator< false, spacepoint_t > replacer;
    detray::serializer2 serializer;
    
    using grid2r = detray::grid2<decltype(replacer), decltype(phi_axis), decltype(z_axis), decltype(serializer)>;

    /*------------------------
      Run spacepoint grouping
      ------------------------*/

    
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

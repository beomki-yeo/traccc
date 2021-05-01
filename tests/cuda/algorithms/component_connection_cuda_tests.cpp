/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "edm/cell.hpp"
#include "edm/cluster.hpp"
#include "algorithms/component_connection.hpp"
#include "cuda/algorithms/component_connection_kernels.cuh"
#include <gtest/gtest.h>
#include <iostream>

// vecmem include
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>


// This defines the local frame test suite
TEST(algorithms, component_connection){

    // Host memory resource used in the test.
    vecmem::host_memory_resource host_mr;

    // cuda managed memory resource used in the test.
    vecmem::cuda::managed_memory_resource mng_mr;

    // mock data
    traccc::host_cell_collection host_cells =
      { { {1, 0, 1., 0. },
          {8, 4, 2., 0.},
          {10, 4, 3., 0.},
          {9, 5, 4., 0.},
          {10, 5, 5., 0},
          {12, 12, 6, 0},
          {3, 13, 7, 0},
          {11, 13, 8, 0},
          {4, 14, 9, 0 } },
        &host_mr };
    
    traccc::cell_module module;
    module.module = 0;

    // cpu
    traccc::component_connection ccl;
    auto clusters = ccl(host_cells, module);

    ASSERT_EQ(clusters.items.size(), 4u);

    // cuda - prepare data
    traccc::host_cell_container mng_cells = {
	traccc::host_cell_container::header_vector(&mng_mr),
        traccc::host_cell_container::item_vector(&mng_mr) };
    
    mng_cells.items.push_back(host_cells);
    mng_cells.headers.push_back(module);    

    traccc::cuda::detail::host_label_container mng_labels={
	vecmem::vector<unsigned int>(1,&mng_mr),
	vecmem::jagged_vector<unsigned int>(1,&mng_mr)
    };
    mng_labels.items[0] = vecmem::vector< unsigned int >(mng_cells.items[0].size(),0);

    // cuda - component connection algorithm
    traccc::cuda::component_connection(mng_cells, mng_labels, &mng_mr);

    ASSERT_EQ(mng_labels.headers[0], 4u); 
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

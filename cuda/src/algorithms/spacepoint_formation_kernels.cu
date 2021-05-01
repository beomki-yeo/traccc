/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "cuda/algorithms/spacepoint_formation_kernels.cuh"
#include "cuda/utils/definitions.hpp"
#include "definitions/algebra.hpp"

namespace traccc {
namespace cuda{
    
    // Declare kernel function for spacepoint formation per module
    __global__
    void spacepoint_formation_kernel(measurement_container_view ms_view,
				     spacepoint_container_view sp_view);

    
    void spacepoint_formation(host_measurement_container& ms_container,
			      host_spacepoint_container& sp_container,
			      vecmem::memory_resource* resource){

	auto ms_data = get_data(ms_container, resource);
	auto sp_data = get_data(sp_container, resource);

	measurement_container_view ms_view(ms_data);
	spacepoint_container_view sp_view(sp_data);
	
	unsigned int num_threads = WARP_SIZE*2; 
	unsigned int num_blocks = ms_data.headers.m_size/num_threads + 1;
	
	spacepoint_formation_kernel<<< num_blocks, num_threads >>>(ms_view, sp_view);
	
	CUDA_ERROR_CHECK(cudaGetLastError());
	CUDA_ERROR_CHECK(cudaDeviceSynchronize());		
      
    }

    __global__
    void spacepoint_formation_kernel(measurement_container_view ms_view,
				     spacepoint_container_view sp_view){

	int gid = blockDim.x * blockIdx.x + threadIdx.x;
	if (gid>=ms_view.items.m_size) return;
	
	device_measurement_container ms_data({ms_view.headers, ms_view.items});
	device_spacepoint_container sp_data({sp_view.headers, sp_view.items});

	auto ms_per_module = ms_data.items.at(gid);
	auto sp_per_module = sp_data.items.at(gid);
       
	for (int i=0; i<ms_per_module.size(); ++i){
	    auto& ms = ms_per_module[i];
	    auto& sp = sp_per_module[i];
	    point3 local_3d = {std::get<0>(ms.local), std::get<1>(ms.local), 0.};

	    // transform the local position to global position
	    sp.global= ms_data.headers[gid].placement.point_to_global(local_3d);

	    // TODO: local error transformation to global one
	}

	sp_data.headers.at(gid) = ms_data.headers.at(gid).module;
    }
}
}

/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "cuda/algorithms/measurement_creation_kernels.cuh"
#include "cuda/utils/definitions.hpp"

namespace traccc {
namespace cuda{

    __global__
    void count_measurements_kernel(cell_container_view cell_view,
				   detail::label_container_view cc_label_view,
				   detail::label_container_view ms_label_view);
    
    void count_measurements(host_cell_container& cells,
			    detail::host_label_container& cc_labels,
			    detail::host_label_container& ms_labels,
			    vecmem::memory_resource* resource){
	auto cell_data = get_data(cells, resource);
	auto cc_label_data = detail::get_data(cc_labels, resource);
	auto ms_label_data = detail::get_data(ms_labels, resource);

	cell_container_view cell_view(cell_data);
	detail::label_container_view cc_label_view(cc_label_data);
	detail::label_container_view ms_label_view(ms_label_data);
	
	unsigned int num_threads = WARP_SIZE*2; 
	unsigned int num_blocks = cell_data.modules.m_size/num_threads + 1;
	
	count_measurements_kernel<<< num_blocks, num_threads >>>(cell_view,
								 cc_label_view,
								 ms_label_view);
	
	CUDA_ERROR_CHECK(cudaGetLastError());
	CUDA_ERROR_CHECK(cudaDeviceSynchronize());	
	
    }

    __global__
    void count_measurements_kernel(cell_container_view cell_view,
				   detail::label_container_view cc_label_view,
				   detail::label_container_view ms_label_view){

	int gid = blockDim.x * blockIdx.x + threadIdx.x;
	if (gid>=cell_view.cells.m_size) return;
	
	device_cell_container cells_data({cell_view.modules, cell_view.cells});
	detail::device_label_container cc_label_data({cc_label_view.counts, cc_label_view.labels});
	detail::device_label_container ms_label_data({ms_label_view.counts, ms_label_view.labels});
	
	auto cells_per_module = cells_data.cells.at(gid);
	auto cc_counts_per_module = cc_label_data.counts.at(gid);
	auto cc_labels_per_module = cc_label_data.labels.at(gid);	
	auto ms_counts_per_module = ms_label_data.counts.at(gid);	
	auto ms_labels_per_module = ms_label_data.labels.at(gid);
	
	for(int i=1; i<=cc_counts_per_module; ++i){
	    scalar weight = 0;
	    for(int j=0; j<cc_labels_per_module.size(); ++j){
		auto& label = cc_labels_per_module[j];
		auto& cell  = cells_per_module[j];
		if ( i == label ){
		    weight += cell.activation;
		}
	    }
	    if (weight > 0){ // need to pass threshold later
		ms_labels_per_module[ms_counts_per_module] = i;
		ms_counts_per_module++;
	    }	    
	}
    }    
}
}

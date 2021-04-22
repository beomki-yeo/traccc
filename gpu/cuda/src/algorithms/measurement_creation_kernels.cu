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

    __global__
    void measurement_creation_kernel(cell_container_view cell_view,
				     detail::label_container_view cc_label_view,
				     detail::label_container_view ms_label_view,
				     measurement_container_view ms_view);
    
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

    void measurement_creation(host_cell_container& ce_container,
			      detail::host_label_container& cc_labels,
			      detail::host_label_container& ms_labels,
			      host_measurement_container& ms_container,
			      vecmem::memory_resource* resource){

	auto cell_data = get_data(ce_container, resource);
	auto cc_label_data = detail::get_data(cc_labels, resource);
	auto ms_label_data = detail::get_data(ms_labels, resource);
	auto ms_data = get_data(ms_container, resource);

	cell_container_view cell_view(cell_data);
	detail::label_container_view cc_label_view(cc_label_data);
	detail::label_container_view ms_label_view(ms_label_data);
	measurement_container_view ms_view(ms_data);
	
	unsigned int num_threads = WARP_SIZE*2; 
	unsigned int num_blocks = cell_data.modules.m_size/num_threads + 1;

	measurement_creation_kernel<<< num_blocks, num_threads >>>(cell_view,
								   cc_label_view,
								   ms_label_view,
								   ms_view);
	
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
	auto cc_labels_per_module = cc_label_data.labels.at(gid);	
	auto ms_labels_per_module = ms_label_data.labels.at(gid);
	
	for(int i=1; i<=cc_label_data.counts.at(gid); ++i){
	    scalar weight = 0;
	    for(int j=0; j<cc_labels_per_module.size(); ++j){
		auto& label = cc_labels_per_module[j];
		auto& cell  = cells_per_module[j];
		if ( i == label ){
		    weight += cell.activation;
		}
	    }
	    if (weight > 0){ // need to pass threshold later
		ms_labels_per_module[ms_label_data.counts.at(gid)] = i;
		ms_label_data.counts.at(gid)++;
	    }	    
	}
    }

    __global__
    void measurement_creation_kernel(cell_container_view cell_view,
				     detail::label_container_view cc_label_view,
				     detail::label_container_view ms_label_view,
				     measurement_container_view ms_view){

	int gid = blockDim.x * blockIdx.x + threadIdx.x;
	if (gid>=cell_view.cells.m_size) return;
	
	device_cell_container cells_data({cell_view.modules, cell_view.cells});
	detail::device_label_container cc_label_data({cc_label_view.counts, cc_label_view.labels});
	detail::device_label_container ms_label_data({ms_label_view.counts, ms_label_view.labels});
	device_measurement_container ms_data({ms_view.modules, ms_view.measurements});
	
	auto cells_per_module = cells_data.cells.at(gid);
	
	auto cc_counts = cc_label_data.counts.at(gid);
	auto cc_labels_per_module = cc_label_data.labels.at(gid);
	
	auto ms_counts = ms_label_data.counts.at(gid);	
	auto ms_labels_per_module = ms_label_data.labels.at(gid);

	auto ms_per_module = ms_data.measurements.at(gid);

	auto pix = traccc::pixel_segmentation{-8.425, -36.025, 0.05, 0.05};
	
	for(int i=0; i<ms_counts; ++i){	    
	    int clabel = ms_labels_per_module[i];
	    scalar total_weight = 0;

	    for (int j=0; j<cc_labels_per_module.size(); ++j){
		if ( clabel == cc_labels_per_module[j] ){
		    auto& cell = cells_per_module[j];
		    scalar weight = cell.activation;
		    total_weight+=weight;
		    auto cell_position = pix(cell.channel0, cell.channel1);
		    ms_per_module[clabel-1].local = ms_per_module[clabel-1].local + weight * cell_position;
		}
	    }
	    ms_per_module[clabel-1].local = 1./total_weight * ms_per_module[clabel-1].local;
	}

	ms_data.modules[gid] = cells_data.modules[gid];
    }
}
}

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

    // Declare kernel function for couting measurement per module
    __global__
    void count_measurements_kernel(cell_container_view cell_view,
				   detail::label_container_view cc_label_view,
				   detail::label_container_view ms_label_view);
    
    // Declare kernel function for measurement creation per module
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
	auto cc_label_data = get_data(cc_labels, resource);
	auto ms_label_data = get_data(ms_labels, resource);

	cell_container_view cell_view(cell_data);
	detail::label_container_view cc_label_view(cc_label_data);
	detail::label_container_view ms_label_view(ms_label_data);
	
	unsigned int num_threads = WARP_SIZE*2; 
	unsigned int num_blocks = cell_data.headers.m_size/num_threads + 1;
	
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
	auto cc_label_data = get_data(cc_labels, resource);
	auto ms_label_data = get_data(ms_labels, resource);
	auto ms_data = get_data(ms_container, resource);

	cell_container_view cell_view(cell_data);
	detail::label_container_view cc_label_view(cc_label_data);
	detail::label_container_view ms_label_view(ms_label_data);
	measurement_container_view ms_view(ms_data);
	
	unsigned int num_threads = WARP_SIZE*2; 
	unsigned int num_blocks = cell_data.headers.m_size/num_threads + 1;

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
	if (gid>=cell_view.items.m_size) return;
	
	device_cell_container cells_data({cell_view.headers, cell_view.items});
	detail::device_label_container cc_label_data({cc_label_view.headers, cc_label_view.items});
	detail::device_label_container ms_label_data({ms_label_view.headers, ms_label_view.items});
	
	auto cells_per_module = cells_data.items.at(gid);
	auto cc_labels_per_module = cc_label_data.items.at(gid);	
	auto ms_labels_per_module = ms_label_data.items.at(gid);

	// Loop over unique labels
	for(int i=1; i<=cc_label_data.headers.at(gid); ++i){
	    scalar weight = 0;
	    // Loop over the labels of cells
	    for(int j=0; j<cc_labels_per_module.size(); ++j){
		auto& label = cc_labels_per_module[j];
		auto& cell  = cells_per_module[j];
		// if the label of a cell (label) matches with unique label (i),
		// sum the weight
		if ( i == label ){
		    weight += cell.activation;
		}
	    }
	    if (weight > 0){ // TODO: pass the threshold from argument
		// if the total weight is larger than threshold,
		// record the unique label (i) as valid label
		ms_labels_per_module[ms_label_data.headers.at(gid)] = i;
		ms_label_data.headers.at(gid)++;
	    }	    
	}
    }

    __global__
    void measurement_creation_kernel(cell_container_view cell_view,
				     detail::label_container_view cc_label_view,
				     detail::label_container_view ms_label_view,
				     measurement_container_view ms_view){

	int gid = blockDim.x * blockIdx.x + threadIdx.x;
	if (gid>=cell_view.items.m_size) return;
	
	device_cell_container cells_data({cell_view.headers, cell_view.items});
	detail::device_label_container cc_label_data({cc_label_view.headers, cc_label_view.items});
	detail::device_label_container ms_label_data({ms_label_view.headers, ms_label_view.items});
	device_measurement_container ms_data({ms_view.headers, ms_view.items});
	
	auto cells_per_module = cells_data.items.at(gid);
	auto module = cells_data.headers.at(gid);
	auto cc_counts = cc_label_data.headers.at(gid);
	auto cc_labels_per_module = cc_label_data.items.at(gid);	
	auto ms_counts = ms_label_data.headers.at(gid);	
	auto ms_labels_per_module = ms_label_data.items.at(gid);
	auto ms_per_module = ms_data.items.at(gid);

	// Width of a pixel
	auto pitch = module.pixel.get_pitch();

	// Loop over the number of measurements per module
	for(int i=0; i<ms_counts; ++i){
	    // clabel: valid label which can form a measurement out of cluster
	    int clabel = ms_labels_per_module[i];
	    scalar total_weight = 0;

	    // Loop over the labels of cells
	    for (int j=0; j<cc_labels_per_module.size(); ++j){
		// Find the cells whose label is same with the valid label (clabel)
		// and compute the weighted average of local position and error
		if ( clabel == cc_labels_per_module[j] ){
		    auto& cell = cells_per_module[j];
		    scalar weight = cell.activation;
		    total_weight+=weight;
		    auto cell_position = module.pixel(cell.channel0, cell.channel1);
		    point2 square_pos = {powf(std::get<0>(cell_position),2),
					 powf(std::get<1>(cell_position),2)};
		    // weighted average of local position
		    ms_per_module[clabel-1].local = ms_per_module[clabel-1].local + weight * cell_position;
		    // weighted average of variance
		    ms_per_module[clabel-1].variance = ms_per_module[clabel-1].variance + weight * square_pos;
		}
	    }
	    // normalize the cell position
	    ms_per_module[clabel-1].local = 1./total_weight * ms_per_module[clabel-1].local;
	    // normalize the variance
	    ms_per_module[clabel-1].variance = 1./total_weight * ms_per_module[clabel-1].variance;
	    // plus pitch^2 / 12
	    ms_per_module[clabel-1].variance = ms_per_module[clabel-1].variance +
		point2{powf(std::get<0>(pitch),2)/12, powf(std::get<1>(pitch),2)/12};
	    // minus <x>^2
	    ms_per_module[clabel-1].variance = ms_per_module[clabel-1].variance -
		point2{powf(std::get<0>(ms_per_module[clabel-1].local),2),
		       powf(std::get<1>(ms_per_module[clabel-1].local),2)};

	}

	// pass the cell_module from cell container
	ms_data.headers[gid] = cells_data.headers[gid];
    }
}
}

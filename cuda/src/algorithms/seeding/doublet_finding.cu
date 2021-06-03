/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include <cuda/algorithms/seeding/doublet_finding.cuh>
#include <cuda/utils/definitions.hpp>

namespace traccc{    
namespace cuda{

__global__
void doublet_finding_kernel(const seedfinder_config config,
			    internal_spacepoint_container_view internal_sp_data,
			    doublet_counter_container_view doublet_counter_view,
			    doublet_container_view mid_bot_doublet_view,
			    doublet_container_view mid_top_doublet_view);    
    
void doublet_finding(const seedfinder_config& config,
		     host_internal_spacepoint_container& internal_sp_container,
		     host_doublet_counter_container& doublet_counter_container,	 
		     host_doublet_container& mid_bot_doublet_container,
		     host_doublet_container& mid_top_doublet_container,
		     vecmem::memory_resource* resource){
    auto internal_sp_data = get_data(internal_sp_container, resource);
    auto doublet_counter_view = get_data(doublet_counter_container, resource);    
    auto mid_bot_doublet_view = get_data(mid_bot_doublet_container, resource);
    auto mid_top_doublet_view = get_data(mid_top_doublet_container, resource);
    
    unsigned int num_threads = WARP_SIZE*4; 
    unsigned int num_blocks = internal_sp_data.headers.m_size;
    
    doublet_finding_kernel<<< num_blocks, num_threads >>>(config,
							  internal_sp_data,
							  doublet_counter_view,
							  mid_bot_doublet_view,
							  mid_top_doublet_view);   
    CUDA_ERROR_CHECK(cudaGetLastError());
    CUDA_ERROR_CHECK(cudaDeviceSynchronize());	        
}
    
__global__
void doublet_finding_kernel(const seedfinder_config config,
			    internal_spacepoint_container_view internal_sp_view,
			    doublet_counter_container_view doublet_counter_view,
			    doublet_container_view mid_bot_doublet_view,
			    doublet_container_view mid_top_doublet_view){

    device_internal_spacepoint_container internal_sp_device({internal_sp_view.headers, internal_sp_view.items});
    device_doublet_counter_container doublet_counter_device({doublet_counter_view.headers, doublet_counter_view.items});
    
    device_doublet_container mid_bot_doublet_device({mid_bot_doublet_view.headers, mid_bot_doublet_view.items});
    device_doublet_container mid_top_doublet_device({mid_top_doublet_view.headers, mid_top_doublet_view.items});
    
    
    auto bin_info = internal_sp_device.headers.at(blockIdx.x);
    auto internal_sp_per_bin = internal_sp_device.items.at(blockIdx.x);

    auto& num_compat_spM_per_bin = doublet_counter_device.headers.at(blockIdx.x);
    auto doublet_counter_per_bin = doublet_counter_device.items.at(blockIdx.x);
    
    auto& num_mid_bot_doublets_per_bin = mid_bot_doublet_device.headers.at(blockIdx.x);
    auto mid_bot_doublets_per_bin = mid_bot_doublet_device.items.at(blockIdx.x);

    auto& num_mid_top_doublets_per_bin = mid_top_doublet_device.headers.at(blockIdx.x);
    auto mid_top_doublets_per_bin = mid_top_doublet_device.items.at(blockIdx.x);
    
    size_t n_iter = internal_sp_per_bin.size()/blockDim.x + 1;
    
    for (size_t i_it = 0; i_it < n_iter; ++i_it){

	auto sp_idx = i_it*blockDim.x + threadIdx.x;
	
	if (sp_idx >= internal_sp_per_bin.size()) {
	    continue;
	}

	auto spM_loc = sp_location({blockIdx.x, sp_idx});
	auto isp = internal_sp_per_bin[sp_idx];

	if (doublet_counter_per_bin[sp_idx].n_mid_bot == 0 ||
	    doublet_counter_per_bin[sp_idx].n_mid_top == 0){
	    continue;
	}	

	size_t n_mid_bot_per_spM = 0;
	size_t n_mid_top_per_spM = 0;

	size_t mid_bot_start_idx = 0;
	size_t mid_top_start_idx = 0;
	
	for (size_t i=0; i<sp_idx; i++){
	    if (doublet_counter_per_bin[i].n_mid_bot == 0 ||
		doublet_counter_per_bin[i].n_mid_top == 0){
		continue;
	    }
	    mid_bot_start_idx += doublet_counter_per_bin[i].n_mid_bot;
	    mid_top_start_idx += doublet_counter_per_bin[i].n_mid_top;	    
	}	
	
		
	for(size_t i_n=0; i_n<bin_info.bottom_idx.counts; ++i_n){		
	    auto neigh_bin = bin_info.bottom_idx.vector_indices[i_n];	    
	    auto neigh_internal_sp_per_bin = internal_sp_device.items.at(neigh_bin);
	    
	    for (size_t spB_idx=0; spB_idx<neigh_internal_sp_per_bin.size(); ++spB_idx){	       		
		auto neigh_isp = neigh_internal_sp_per_bin[spB_idx];		
		
		if (doublet_finding_helper::isCompatible(isp, neigh_isp, config, true)){
		    
		    auto spB_loc = sp_location({neigh_bin, spB_idx});
		    auto lin = doublet_finding_helper::transform_coordinates(isp, neigh_isp, true);
		    
		    if (n_mid_bot_per_spM < doublet_counter_per_bin[sp_idx].n_mid_bot &&
			num_mid_bot_doublets_per_bin < mid_bot_doublets_per_bin.size()){			
			size_t pos = mid_bot_start_idx + n_mid_bot_per_spM;	  
			if (pos>=mid_bot_doublets_per_bin.size()) {
			    continue;
			}
			
			mid_bot_doublets_per_bin[pos] = doublet({spM_loc,
								 spB_loc,
								 lin});
			
			atomicAdd(&num_mid_bot_doublets_per_bin,1);
			n_mid_bot_per_spM++;
			
		    }
		    
		}

		if (doublet_finding_helper::isCompatible(isp, neigh_isp, config, false)){
		    
		    auto spT_loc = sp_location({neigh_bin, spB_idx});
		    auto lin = doublet_finding_helper::transform_coordinates(isp, neigh_isp, false);

		    if (n_mid_top_per_spM < doublet_counter_per_bin[sp_idx].n_mid_top &&
			num_mid_top_doublets_per_bin < mid_top_doublets_per_bin.size()){
			
			size_t pos = mid_top_start_idx + n_mid_top_per_spM;
			if (pos>=mid_top_doublets_per_bin.size()) {
			    continue;
			}
			
			mid_top_doublets_per_bin[pos] = doublet({spM_loc,
								 spT_loc,
								 lin});
			
			atomicAdd(&num_mid_top_doublets_per_bin,1);	      
			n_mid_top_per_spM++;

		    }		    
		}
	    }				    	    
	}

	
	/*
	bool hasCompatBottom = false;
	bool hasCompatTop = false;
     
	for(size_t i_n=0; i_n<bin_info.bottom_idx.counts; ++i_n){	   

	    if (hasCompatBottom && hasCompatTop){
		break;
	    }
	    
	    auto neigh_bin = bin_info.bottom_idx.vector_indices[i_n];
	    
	    auto neigh_internal_sp_per_bin = internal_sp_device.items.at(neigh_bin);		
	    for (size_t spB_idx=0; spB_idx<neigh_internal_sp_per_bin.size(); ++spB_idx){		
		auto neigh_isp = neigh_internal_sp_per_bin[spB_idx];
		if (!hasCompatBottom){
		    hasCompatBottom = doublet_finding_helper::isCompatible(isp, neigh_isp, config, true);
		}
		if (!hasCompatTop){
		    hasCompatTop = doublet_finding_helper::isCompatible(isp, neigh_isp, config, false);
		}
	    }
	}

	// Skip if there is not compatible hits on any of bottom and top side
	if (!hasCompatBottom || !hasCompatTop) {
	    continue;
	}
	
	for(size_t i_n=0; i_n<bin_info.bottom_idx.counts; ++i_n){		
	    auto neigh_bin = bin_info.bottom_idx.vector_indices[i_n];
	    
	    auto neigh_internal_sp_per_bin = internal_sp_device.items.at(neigh_bin);		
	    for (size_t spB_idx=0; spB_idx<neigh_internal_sp_per_bin.size(); ++spB_idx){	       		
		auto neigh_isp = neigh_internal_sp_per_bin[spB_idx];		

		if (!doublet_finding_helper::isCompatible(isp, neigh_isp, config, true)) continue;
		
		auto spB_loc = sp_location({neigh_bin, spB_idx});
		auto lin = doublet_finding_helper::transform_coordinates(isp, neigh_isp, true);
		auto pos = atomicAdd(&num_mid_bot_doublets_per_bin,1);
	       
		if (pos>=mid_bot_doublets_per_bin.size()){
		    num_mid_bot_doublets_per_bin = mid_bot_doublets_per_bin.size();
		    continue;
		}
		
		mid_bot_doublets_per_bin[pos] = doublet({spM_loc, spB_loc, lin});
		
		
	    }				
	}

	for(size_t i_n=0; i_n<bin_info.top_idx.counts; ++i_n){		
	    auto neigh_bin = bin_info.top_idx.vector_indices[i_n];
	    
	    auto neigh_internal_sp_per_bin = internal_sp_device.items.at(neigh_bin);		
	    for (size_t spT_idx=0; spT_idx<neigh_internal_sp_per_bin.size(); ++spT_idx){		
		auto neigh_isp = neigh_internal_sp_per_bin[spT_idx];		
	    
		if (!doublet_finding_helper::isCompatible(isp, neigh_isp, config, false)) continue;		
		auto lin = doublet_finding_helper::transform_coordinates(isp, neigh_isp, false);
		auto spT_loc = sp_location({neigh_bin, spT_idx});

		auto pos = atomicAdd(&num_mid_top_doublets_per_bin,1);
		
		if (pos>=mid_top_doublets_per_bin.size()){
		    num_mid_top_doublets_per_bin = mid_top_doublets_per_bin.size();
		    continue;
		}
		
		mid_top_doublets_per_bin[pos] = doublet({spM_loc, spT_loc, lin});   
	    }				
	}	   
	*/
    }    
}
    
}// namespace cuda
}// namespace traccc

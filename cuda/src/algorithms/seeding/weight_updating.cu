/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include <cuda/algorithms/seeding/weight_updating.cuh>
#include <cuda/utils/definitions.hpp>

namespace traccc{    
namespace cuda{

__global__
void weight_updating_kernel(const seedfilter_config filter_config,
			    internal_spacepoint_container_view internal_sp_view,
			    triplet_container_view triplet_view,
			    vecmem::data::jagged_vector_view< size_t > n_triplets_per_mb_view,
			    vecmem::data::jagged_vector_view< float > compatseed_view);
    
void weight_updating(const seedfilter_config& filter_config,
		     host_internal_spacepoint_container& internal_sp_container,
		     host_triplet_container& triplet_container,
		     vecmem::jagged_vector< size_t >& n_triplets_per_mb,
		     vecmem::jagged_vector< float >& compatseed_container,
		     vecmem::memory_resource* resource
		     ){

    auto internal_sp_data = get_data(internal_sp_container, resource);
    auto triplet_view = get_data(triplet_container, resource);
    auto n_triplets_per_mb_view = vecmem::get_data(n_triplets_per_mb, resource);
    auto compatseed_view = vecmem::get_data(compatseed_container, resource);
    
    unsigned int num_threads = WARP_SIZE*2; 
    unsigned int num_blocks = internal_sp_data.headers.m_size;

    
    weight_updating_kernel<<< num_blocks, num_threads >>>(filter_config,
							  internal_sp_data,
							  triplet_view,
							  n_triplets_per_mb_view,
							  compatseed_view);   
    
    CUDA_ERROR_CHECK(cudaGetLastError());
    CUDA_ERROR_CHECK(cudaDeviceSynchronize());	        
    
}

__global__
void weight_updating_kernel(const seedfilter_config filter_config,
			    internal_spacepoint_container_view internal_sp_view,
			    triplet_container_view triplet_view,
			    vecmem::data::jagged_vector_view< size_t > n_triplets_per_mb_view,
			    vecmem::data::jagged_vector_view< float > compatseed_view){

    device_internal_spacepoint_container internal_sp_device({internal_sp_view.headers, internal_sp_view.items});
    
    device_triplet_container triplet_device({triplet_view.headers, triplet_view.items});
    vecmem::jagged_device_vector< size_t > n_triplets_per_mb_device(n_triplets_per_mb_view);
    vecmem::jagged_device_vector< float > compatseed_device(compatseed_view);
   
    auto bin_info = internal_sp_device.headers.at(blockIdx.x);
    auto internal_sp_per_bin = internal_sp_device.items.at(blockIdx.x);

    auto& num_triplets_per_bin = triplet_device.headers.at(blockIdx.x);
    auto triplets_per_bin = triplet_device.items.at(blockIdx.x);

    auto n_triplets_per_mb = n_triplets_per_mb_device.at(blockIdx.x);
    auto compatseed = compatseed_device.at(blockIdx.x);
    
    size_t n_iter = num_triplets_per_bin/blockDim.x + 1;

    for (size_t i_it = 0; i_it < n_iter; ++i_it){
	auto tr_idx = i_it*blockDim.x + threadIdx.x;
	auto& triplet = triplets_per_bin[tr_idx];

	if (tr_idx >= num_triplets_per_bin){
	    continue;
	}

	if (n_triplets_per_mb.size() == 0){
	    continue;
	}

	size_t start_idx = 0;
	size_t end_idx = 0;	

	for (auto n_triplets : n_triplets_per_mb){
	    
	    end_idx += n_triplets;	    
	    if (end_idx > tr_idx){
		break;
	    }
	    start_idx += n_triplets;	    
	}

	auto& spT_idx = triplet.sp3;
	auto& current_spT = internal_sp_per_bin[spT_idx.sp_idx];
	float currentTop_r = current_spT.radius();
	
	// if two compatible seeds with high distance in r are found, compatible
	// seeds span 5 layers
	// -> very good seed		
	float lowerLimitCurv = triplet.curvature - filter_config.deltaInvHelixDiameter;
	float upperLimitCurv = triplet.curvature + filter_config.deltaInvHelixDiameter;	
	int n_compatseed = 0;
	
	// iterate over triplets
	for (auto tr_it = triplets_per_bin.begin()+start_idx;
	     tr_it!= triplets_per_bin.begin()+end_idx ;
	     tr_it++){

	    if (triplet == *tr_it){
		continue;
	    }

	    auto& other_triplet = *tr_it;
	    auto other_spT = internal_sp_per_bin[(*tr_it).sp3.sp_idx];
	    
	    // compared top SP should have at least deltaRMin distance
	    float otherTop_r = other_spT.radius();
	    float deltaR = currentTop_r - otherTop_r;
	    if (std::abs(deltaR) < filter_config.deltaRMin) {
		continue;
	    }
	    
	    // curvature difference within limits?
	    // TODO: how much slower than sorting all vectors by curvature
	    // and breaking out of loop? i.e. is vector size large (e.g. in jets?)
	    if (other_triplet.curvature < lowerLimitCurv) {
		continue;
	    }
	    if (other_triplet.curvature > upperLimitCurv) {
		continue;
	    }

	    bool newCompSeed = true;
	    
	    //for (float previousDiameter : compatibleSeedR) {
	    for (size_t i_s = start_idx; i_s < start_idx+n_compatseed; ++i_s){
		float previousDiameter = compatseed[i_s];
		
		// original ATLAS code uses higher min distance for 2nd found compatible
		// seed (20mm instead of 5mm)
		// add new compatible seed only if distance larger than rmin to all
		// other compatible seeds
		if (std::abs(previousDiameter - otherTop_r) < filter_config.deltaRMin) {
		    newCompSeed = false;
		    break;
		}
	    }
	    	    
	    if (newCompSeed) {
		n_compatseed++;
		compatseed[start_idx+n_compatseed] = otherTop_r;
		triplet.weight += filter_config.compatSeedWeight;
	    }
	    
	    if (n_compatseed >= filter_config.compatSeedLimit) {
		break;
	    }	   	    
	}
    }

    
    
}
    
}// namespace cuda
}// namespace traccc    

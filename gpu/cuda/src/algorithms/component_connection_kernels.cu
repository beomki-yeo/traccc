/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "cuda/algorithms/component_connection_kernels.cuh"
#include "cuda/utils/definitions.hpp"

__device__
unsigned int find_root(const vecmem::device_vector< unsigned int >& L, unsigned int e){
  unsigned int r = e;
  while (L[r] != r){
    r = L[r];
  }
  return r;
}

__device__
unsigned int make_union(vecmem::device_vector< unsigned int >& L, unsigned int e1, unsigned int e2){
  int e;
  if (e1 < e2){
    e = e1;
    L[e2] = e;
  } else {
    e = e2;
    L[e1] = e;
  }
  return e;
}

__device__
bool is_adjacent(traccc::cell a, traccc::cell b){
  return (a.channel0 - b.channel0)*(a.channel0 - b.channel0) <= 1 
    and (a.channel1 - b.channel1)*(a.channel1 - b.channel1) <= 1; 
}

__device__
bool is_far_enough(traccc::cell a, traccc::cell b){
  return (a.channel1 - b.channel1) > 1; 
}

__device__
unsigned int sparse_ccl(const vecmem::device_vector< traccc::cell > cells,
			vecmem::device_vector< unsigned int >& L){    
    unsigned int start_j = 0;
    for (unsigned int i=0; i<cells.size(); ++i){
	L[i] = i;
	int ai = i;
	if (i > 0){
	    for (unsigned int j = start_j; j < i; ++j){
		if (is_adjacent(cells[i], cells[j])){
		    ai = make_union(L, ai, find_root(L, j));
		} else if (is_far_enough(cells[i], cells[j])){
		    ++start_j;
		}
	    }
	}    
    }
    
    // second scan: transitive closure
    unsigned int labels = 0;
    
    for (unsigned int i = 0; i < cells.size(); ++i){
	unsigned int l = 0;
	if (L[i] == i){
	    ++labels;
	    l = labels; 
	} else {
      l = L[L[i]];
	}
	L[i] = l;
    }
    
    return labels;
}


namespace traccc {
namespace cuda{

    __global__
    void cc_kernel(cell_container_view cell_data,
		   detail::label_container_view label_data);
    
    void component_connection(host_cell_container& cells_per_event,
			      detail::host_label_container& labels_per_event,
			      vecmem::memory_resource* resource){
	auto cell_data = get_data(cells_per_event, resource);
	auto label_data = detail::get_data(labels_per_event, resource);
	cell_container_view cell_view(cell_data);
	detail::label_container_view label_view(label_data);

	unsigned int num_threads = WARP_SIZE * 2; 
	unsigned int num_blocks = cell_data.modules.m_size/num_threads + 1;

	cc_kernel<<< num_blocks, num_threads >>>(cell_view, label_view);
	
	CUDA_ERROR_CHECK(cudaGetLastError());
	CUDA_ERROR_CHECK(cudaDeviceSynchronize());	
    }

    __global__
    void cc_kernel(cell_container_view cell_data,
		   detail::label_container_view label_data){
	int gid = blockDim.x * blockIdx.x + threadIdx.x;
	if (gid>=cell_data.cells.m_size) return;

	device_cell_container cells_device({cell_data.modules, cell_data.cells});
	detail::device_label_container labels_device({label_data.num_label, label_data.labels});
	auto num_label = labels_device.num_label;
	auto cells_per_module = cells_device.cells.at(gid);
	auto labels_per_module = labels_device.labels.at(gid);
		
	num_label[gid] = sparse_ccl(cells_per_module, labels_per_module);
	return;	
    }
}
}

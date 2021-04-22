/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "edm/cell.hpp"
#include "edm/cluster.hpp"
#include "edm/measurement.hpp"
#include "edm/spacepoint.hpp"
#include "geometry/pixel_segmentation.hpp"
#include "algorithms/component_connection.hpp"
#include "algorithms/measurement_creation.hpp"
#include "algorithms/spacepoint_formation.hpp"
#include "csv/csv_io.hpp"

// cuda include
#include "cuda/algorithms/component_connection_kernels.cuh"
#include "cuda/algorithms/measurement_creation_kernels.cuh"
#include "cuda/algorithms/spacepoint_formation_kernels.cuh"

// vecmem include
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

// std 
#include <iostream>
#include <chrono>

int seq_run(const std::string& detector_file, const std::string& cells_dir, unsigned int events)
{
    auto env_d_d = std::getenv("TRACCC_TEST_DATA_DIR");
    if (env_d_d == nullptr)
    {
        throw std::ios_base::failure("Test data directory not found. Please set TRACCC_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d) + std::string("/");

    // Read the surface transforms
    std::string io_detector_file = data_directory + detector_file;
    traccc::surface_reader sreader(io_detector_file, {"geometry_id", "cx", "cy", "cz", "rot_xu", "rot_xv", "rot_xw", "rot_zu", "rot_zv", "rot_zw"});
    auto surface_transforms = traccc::read_surfaces(sreader);

    // Algorithms
    traccc::component_connection cc;
    traccc::measurement_creation mt;
    traccc::spacepoint_formation sp;

    // Output stats
    uint64_t n_cells = 0;
    uint64_t m_modules = 0;
    uint64_t n_clusters = 0;
    uint64_t n_measurements = 0;
    uint64_t n_space_points = 0;

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::managed_memory_resource mng_mr;

    // Elapsed time
    float read_cpu(0), ccl_cpu(0), ms_cpu(0), sp_cpu(0), total_cpu(0);
    float read_cuda(0), ccl_cuda(0), ms_cuda(0), sp_cuda(0), total_cuda(0);

    
    // Loop over events
    for (unsigned int event = 0; event < events; ++event){

	/*time*/ auto start_total_cpu = std::chrono::system_clock::now();
	
        // Read the cells from the relevant event file
        std::string event_string = "000000000";
        std::string event_number = std::to_string(event);
        event_string.replace(event_string.size()-event_number.size(), event_number.size(), event_number);

        std::string io_cells_file = data_directory+cells_dir+std::string("/event")+event_string+std::string("-cells.csv");

	/*time*/ auto start_read_cpu = std::chrono::system_clock::now();
        traccc::cell_reader creader(io_cells_file, {"geometry_id", "hit_id", "cannel0", "channel1", "activation", "time"});	
        traccc::host_cell_container cells_per_event = traccc::read_cells(creader, host_mr, surface_transforms);

	/*time*/ auto end_read_cpu = std::chrono::system_clock::now();
	/*time*/ std::chrono::duration<double> time_read_cpu = end_read_cpu - start_read_cpu; 
	/*time*/ read_cpu += time_read_cpu.count();
	
        m_modules += cells_per_event.modules.size();

        // Output containers
        traccc::host_measurement_container measurements_per_event;
        traccc::host_spacepoint_container spacepoints_per_event;
        measurements_per_event.modules.reserve(cells_per_event.modules.size());
	measurements_per_event.measurements.reserve(cells_per_event.modules.size());
        spacepoints_per_event.modules.reserve(cells_per_event.modules.size());
	spacepoints_per_event.spacepoints.reserve(cells_per_event.modules.size()); 
	
        for (std::size_t i = 0; i < cells_per_event.cells.size(); ++i )
        {
            // The algorithmic code part: start
	    /*time*/ auto start_ccl_cpu = std::chrono::system_clock::now();
	    
	    auto& module = cells_per_event.modules[i];
            traccc::cluster_collection clusters_per_module = cc(cells_per_event.cells[i], cells_per_event.modules[i]);

	    /*time*/ auto end_ccl_cpu = std::chrono::system_clock::now();
	    /*time*/ std::chrono::duration<double> time_ccl_cpu = end_ccl_cpu - start_ccl_cpu; 
	    /*time*/ ccl_cpu += time_ccl_cpu.count();
	    
            clusters_per_module.position_from_cell = traccc::pixel_segmentation{-8.425, -36.025, 0.05, 0.05};


	    /*time*/ auto start_ms_cpu = std::chrono::system_clock::now();
            traccc::host_measurement_collection measurements_per_module = mt(clusters_per_module);

	    /*time*/ auto end_ms_cpu = std::chrono::system_clock::now();
	    /*time*/ std::chrono::duration<double> time_ms_cpu = end_ms_cpu - start_ms_cpu; 
	    /*time*/ ms_cpu += time_ms_cpu.count();

	    /*time*/ auto start_sp_cpu = std::chrono::system_clock::now();
	    
            traccc::host_spacepoint_collection spacepoints_per_module = sp(module, measurements_per_module);
            // The algorithmnic code part: end

	    /*time*/ auto end_sp_cpu = std::chrono::system_clock::now();
	    /*time*/ std::chrono::duration<double> time_sp_cpu = end_sp_cpu - start_sp_cpu; 
	    /*time*/ sp_cpu += time_sp_cpu.count();
	    
            n_cells += cells_per_event.cells[i].size();
            n_clusters += clusters_per_module.items.size();
            n_measurements += measurements_per_module.size();
            n_space_points += spacepoints_per_module.size();

            measurements_per_event.measurements.push_back(std::move(measurements_per_module));
	    measurements_per_event.modules.push_back(module);
	    
            spacepoints_per_event.spacepoints.push_back(std::move(spacepoints_per_module));
	    spacepoints_per_event.modules.push_back(module.module);
        }

	/*time*/ auto end_total_cpu = std::chrono::system_clock::now();
	/*time*/ std::chrono::duration<double> time_total_cpu = end_total_cpu - start_total_cpu; 
	/*time*/ total_cpu += time_total_cpu.count();

	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	/*time*/ auto start_total_cuda = std::chrono::system_clock::now();	
  
	//// cuda - read cell data with managed memory resource
	/*time*/ auto start_read_cuda = std::chrono::system_clock::now();
	
	traccc::cell_reader creader_for_cuda(io_cells_file, {"geometry_id", "hit_id", "cannel0", "channel1", "activation", "time"});
        traccc::host_cell_container ce_container = traccc::read_cells(creader_for_cuda, mng_mr, surface_transforms);

	// prepare ccl label vector
	int n_modules = ce_container.cells.size();
	traccc::cuda::detail::host_label_container cc_labels = {
	    vecmem::vector< unsigned int >(n_modules, 0, &mng_mr),
	    vecmem::jagged_vector< unsigned int >(n_modules, &mng_mr)
	};
	for (int i=0; i<n_modules; ++i){
	    cc_labels.labels[i]=vecmem::vector<unsigned int>(ce_container.cells[i].size(),0);
	}
	
	/*time*/ auto end_read_cuda = std::chrono::system_clock::now();
	/*time*/ std::chrono::duration<double> time_read_cuda = end_read_cuda - start_read_cuda; 
	/*time*/ read_cuda += time_read_cuda.count();

	//// cuda - component connection algorithm
	
	/*time*/ auto start_ccl_cuda = std::chrono::system_clock::now();

	traccc::cuda::component_connection(ce_container,
					   cc_labels,
					   &mng_mr);

	/*time*/ auto end_ccl_cuda = std::chrono::system_clock::now();
	/*time*/ std::chrono::duration<double> time_ccl_cuda = end_ccl_cuda - start_ccl_cuda; 
	/*time*/ ccl_cuda += time_ccl_cuda.count();

	//// cuda - count measurements

	/*time*/ auto start_ms_cuda = std::chrono::system_clock::now();
	traccc::cuda::detail::host_label_container ms_labels = {
	    vecmem::vector< unsigned int >(n_modules, 0, &mng_mr),
	    vecmem::jagged_vector< unsigned int >(n_modules, &mng_mr)
	};
	for (int i=0; i<n_modules; ++i){
	    ms_labels.labels[i] = vecmem::vector< unsigned int >(cc_labels.labels[i].size(),0);
	}

	traccc::cuda::count_measurements(ce_container, cc_labels,
					 ms_labels, &mng_mr);
	
	//// cuda - mesaurements creation
	
	traccc::host_measurement_container ms_container = {
	     traccc::host_measurement_container::cell_module_vector( n_modules, &mng_mr ),	     
	     traccc::host_measurement_container::measurement_vector( n_modules, &mng_mr ) 
	};	
	for ( int i=0; i< n_modules; ++i ){
	    ms_container.measurements[i] =
		vecmem::vector< traccc::measurement >(ms_labels.counts[i]);
	}	
	
	traccc::cuda::measurement_creation(ce_container,
					   cc_labels,
					   ms_labels,
					   ms_container,
					   &mng_mr);

	/*time*/ auto end_ms_cuda = std::chrono::system_clock::now();
	/*time*/ std::chrono::duration<double> time_ms_cuda = end_ms_cuda - start_ms_cuda; 
	/*time*/ ms_cuda += time_ms_cuda.count();
	
	//// cuda - spacepoint formation

	/*time*/ auto start_sp_cuda = std::chrono::system_clock::now();
	traccc::host_spacepoint_container sp_container = {
	     traccc::host_spacepoint_container::module_vector( n_modules, &mng_mr ),	     
	     traccc::host_spacepoint_container::spacepoint_vector( n_modules, &mng_mr ) 
	};	
	for ( int i=0; i< n_modules; ++i ){
	    sp_container.spacepoints[i] =
		vecmem::vector< traccc::spacepoint >(ms_labels.counts[i]);
	}	

	traccc::cuda::spacepoint_formation(ms_container,
					   sp_container,
					   &mng_mr);

	/*time*/ auto end_sp_cuda = std::chrono::system_clock::now();
	/*time*/ std::chrono::duration<double> time_sp_cuda = end_sp_cuda - start_sp_cuda; 
	/*time*/ sp_cuda += time_sp_cuda.count();

	/*time*/ auto end_total_cuda = std::chrono::system_clock::now();
	/*time*/ std::chrono::duration<double> time_total_cuda = end_total_cuda - start_total_cuda; 
	/*time*/ total_cuda += time_total_cuda.count();

	
	// assure we have same results for spacepoints	
	for(int i=0; i<sp_container.spacepoints.size(); i++){
	    for(int j=0; j<sp_container.spacepoints[i].size(); j++){
		auto sp_cuda = sp_container.spacepoints[i][j];
		auto sp_cpu = spacepoints_per_event.spacepoints[i][j];
		if (sp_cuda.global[0] != sp_cpu.global[0] ||
		    sp_cuda.global[1] != sp_cpu.global[1] ||
		    sp_cuda.global[2] != sp_cpu.global[2]){
		    std::cout << "WARNING: spacepoint mismatch detected" << std::endl;
		}
	    }
	}
	
        traccc::measurement_writer mwriter{std::string("event")+event_number+"-measurements.csv"};
	for (int i=0; i<measurements_per_event.measurements.size(); ++i){
	    auto measurements_per_module = measurements_per_event.measurements[i];
            auto module = measurements_per_event.modules[i];
            for (const auto& measurement : measurements_per_module){
                const auto& local = measurement.local;
                mwriter.append({ module.module, local[0], local[1], 0., 0.});
            }
        }	
	
        traccc::spacepoint_writer spwriter{std::string("event")+event_number+"-spacepoints.csv"};
	for (int i=0; i<spacepoints_per_event.spacepoints.size(); ++i){
	    auto spacepoints_per_module = spacepoints_per_event.spacepoints[i];
            auto module = spacepoints_per_event.modules[i];
	    
            for (const auto& spacepoint : spacepoints_per_module){
                const auto& pos = spacepoint.global;
                spwriter.append({ module, pos[0], pos[1], pos[2], 0., 0., 0.});
            }
        }
    }

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read    " << n_cells << " cells from " << m_modules << " modules" << std::endl;
    std::cout << "- created " << n_clusters << " clusters. " << std::endl;
    std::cout << "- created " << n_measurements << " measurements. " << std::endl;
    std::cout << "- created " << n_space_points << " space points. " << std::endl;

    std::cout << "==> Elapsed time ... " << std::endl;
    std::cout << "- reading time | " << " cpu: " << read_cpu << " |  cuda: " << read_cuda << std::endl;
    std::cout << "- ccl time     | " << " cpu: " << ccl_cpu << " |  cuda: " << ccl_cuda << std::endl;
    std::cout << "- ms time      | " << " cpu: " << ms_cpu << " |  cuda: " << ms_cuda << std::endl;
    std::cout << "- sp time      | " << " cpu: " << sp_cpu << " |  cuda: " << sp_cuda << std::endl;
    std::cout << "- total time      | " << " cpu: " << total_cpu << " |  cuda: " << total_cuda << std::endl;
    
    return 0;
}

// The main routine
//
int main(int argc, char *argv[])
{
    if (argc < 4){
        std::cout << "Not enough arguments, minimum requirement: " << std::endl;
        std::cout << "./seq_example <detector_file> <cell_directory> <events>" << std::endl;
        return -1;
    }

    auto detector_file = std::string(argv[1]);
    auto cell_directory = std::string(argv[2]);
    auto events = std::atoi(argv[3]);

    std::cout << "Running ./seq_exammple " << detector_file << " " << cell_directory << " " << events << std::endl;
    return seq_run(detector_file, cell_directory, events);
}

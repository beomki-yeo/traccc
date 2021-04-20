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
    float read_cpu(0), ccl_cpu(0), total_cpu(0);
    float read_cuda(0), ccl_cuda(0), total_cuda(0);
    
    // Loop over events
    for (unsigned int event = 0; event < events; ++event){

        // Read the cells from the relevant event file
        std::string event_string = "000000000";
        std::string event_number = std::to_string(event);
        event_string.replace(event_string.size()-event_number.size(), event_number.size(), event_number);

        std::string io_cells_file = data_directory+cells_dir+std::string("/event")+event_string+std::string("-cells.csv");

	/**/ auto start_read_cpu = std::chrono::system_clock::now();
        traccc::cell_reader creader(io_cells_file, {"geometry_id", "hit_id", "cannel0", "channel1", "activation", "time"});	
        traccc::host_cell_container cells_per_event = traccc::read_cells(creader, host_mr, surface_transforms);

	/**/ auto end_read_cpu = std::chrono::system_clock::now();
	/**/ std::chrono::duration<double> time_read_cpu = end_read_cpu - start_read_cpu; 
	/**/ read_cpu += time_read_cpu.count();
	
        m_modules += cells_per_event.modules.size();

        // Output containers
        traccc::measurement_container measurements_per_event;
        traccc::spacepoint_container spacepoints_per_event;
        measurements_per_event.reserve(cells_per_event.modules.size());
        spacepoints_per_event.reserve(cells_per_event.modules.size());

        for (std::size_t i = 0; i < cells_per_event.cells.size(); ++i )
        {
            // The algorithmic code part: start
	    /**/ auto start_ccl_cpu = std::chrono::system_clock::now();
	    
            traccc::cluster_collection clusters_per_module = cc(cells_per_event.cells[i], cells_per_event.modules[i]);

	    /**/ auto end_ccl_cpu = std::chrono::system_clock::now();
	    /**/ std::chrono::duration<double> time_ccl_cpu = end_ccl_cpu - start_ccl_cpu; 
	    /**/ ccl_cpu += time_ccl_cpu.count();
	    
            clusters_per_module.position_from_cell = traccc::pixel_segmentation{-8.425, -36.025, 0.05, 0.05};
            traccc::measurement_collection measurements_per_module = mt(clusters_per_module);
            traccc::spacepoint_collection spacepoints_per_module = sp(measurements_per_module);
            // The algorithmnic code part: end

            n_cells += cells_per_event.cells[i].size();
            n_clusters += clusters_per_module.items.size();
            n_measurements += measurements_per_module.items.size();
            n_space_points += spacepoints_per_module.items.size();

            measurements_per_event.push_back(std::move(measurements_per_module));
            spacepoints_per_event.push_back(std::move(spacepoints_per_module));
        }

	// cuda - read cell data with managed memory resource
	/**/ auto start_read_cuda = std::chrono::system_clock::now();
	
	traccc::cell_reader creader_for_cuda(io_cells_file, {"geometry_id", "hit_id", "cannel0", "channel1", "activation", "time"});
        traccc::host_cell_container mng_cells = traccc::read_cells(creader_for_cuda, mng_mr, surface_transforms);
	auto mng_labels = traccc::cuda::detail::get_label_from_cell(mng_cells, &mng_mr);
	/**/ auto end_read_cuda = std::chrono::system_clock::now();
	/**/ std::chrono::duration<double> time_read_cuda = end_read_cuda - start_read_cuda; 
	/**/ read_cuda += time_read_cuda.count();

	// cuda - component connection algorithm
	
	/**/ auto start_ccl_cuda = std::chrono::system_clock::now();

	traccc::cuda::component_connection(mng_cells, mng_labels, &mng_mr);

	/**/ auto end_ccl_cuda = std::chrono::system_clock::now();
	/**/ std::chrono::duration<double> time_ccl_cuda = end_ccl_cuda - start_ccl_cuda; 
	/**/ ccl_cuda += time_ccl_cuda.count();
	
        traccc::measurement_writer mwriter{std::string("event")+event_number+"-measurements.csv"};
        for (const auto& measurements_per_module : measurements_per_event){
            auto module = measurements_per_module.module;
            for (const auto& measurement : measurements_per_module.items){
                const auto& local = measurement.local;
                mwriter.append({ module, local[0], local[1], 0., 0.});
            }
        }

        traccc::spacepoint_writer spwriter{std::string("event")+event_number+"-spacepoints.csv"};
        for (const auto& spacepoint_per_module : spacepoints_per_event){
            auto module = spacepoint_per_module.module;
            for (const auto& spacepoint : spacepoint_per_module.items){
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

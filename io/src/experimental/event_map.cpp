/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/experimental/event_map.hpp"

#include "../csv/make_cell_reader.hpp"
#include "../csv/make_hit_reader.hpp"
#include "../csv/make_measurement_hit_id_reader.hpp"
#include "../csv/make_measurement_reader.hpp"
#include "../csv/make_particle_reader.hpp"
#include "traccc/io/experimental/read_cells.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/utils.hpp"

// Project include(s).
#include "traccc/clusterization/component_connection.hpp"
#include "traccc/clusterization/measurement_creation.hpp"

namespace traccc::io::experimental {

event_map::event_map(std::size_t event, vecmem::memory_resource& resource,
                     const std::string& digi_config_file,
                     const std::string& particle_dir,
                     const std::string& hit_dir, const std::string& cell_dir) {

    /**************************
     * Read the input files
     **************************/

    // Containers
    std::vector<traccc::io::csv::particle> particles;
    std::vector<traccc::io::csv::hit> hits;
    std::vector<traccc::io::csv::measurement_hit_id> meas_hit_ids;
    std::vector<traccc::io::csv::cell> cells;

    // Read Particles
    std::string io_particle_file =
        io::data_directory() + particle_dir +
        io::get_event_filename(event, "-particles.csv");

    auto preader = io::csv::make_particle_reader(io_particle_file);

    traccc::io::csv::particle io_particle;
    while (preader.read(io_particle)) {
        point3 pos{io_particle.vx, io_particle.vy, io_particle.vz};
        vector3 mom{io_particle.px, io_particle.py, io_particle.pz};

        particles.push_back(io_particle);
        ptc_map[io_particle.particle_id] =
            particle{io_particle.particle_id, io_particle.particle_type,
                     io_particle.process,     pos,
                     io_particle.vt,          mom,
                     io_particle.m,           io_particle.q};
    }

    // Read Meas-to-Hit IDs
    std::string io_meas_hit_id_file =
        io::data_directory() + hit_dir +
        io::get_event_filename(event, "-measurement-simhit-map.csv");

    auto mhid_reader =
        io::csv::make_measurement_hit_id_reader(io_meas_hit_id_file);

    traccc::io::csv::measurement_hit_id io_mh_id;
    while (mhid_reader.read(io_mh_id)) {
        meas_hit_ids.push_back(io_mh_id);
    }

    // Read Hits
    std::string io_hit_file = io::data_directory() + hit_dir +
                              io::get_event_filename(event, "-hits.csv");

    auto hreader = io::csv::make_hit_reader(io_hit_file);

    traccc::io::csv::hit io_hit;

    hit_id hid = 0u;
    while (hreader.read(io_hit)) {
        hits.push_back(io_hit);

        spacepoint sp;
        sp.global = {io_hit.tx, io_hit.ty, io_hit.tz};
        h_map[meas_hit_ids[hid].hit_id] = sp;
        hid++;
    }

    // Read Cells
    std::string io_cell_file = io::data_directory() + cell_dir +
                               io::get_event_filename(event, "-cells.csv");

    auto creader = io::csv::make_cell_reader(io_cell_file);

    traccc::io::csv::cell io_cell;
    while (creader.read(io_cell)) {
        cells.push_back(io_cell);
    }

    /***************************************
     * Generate measurement to cluster map
     ***************************************/

    // CCA algorithms
    component_connection cc(resource);
    measurement_creation mc(resource);

    // Read the digitization configuration file
    auto digi_cfg =
        io::experimental::read_digitization_config(digi_config_file);

    // Read the cells from the relevant event file
    traccc::io::cell_reader_output readOut(&resource);
    io::experimental::read_cells(readOut, event, cell_dir,
                                 traccc::data_format::csv, &digi_cfg);
    cell_collection_types::host& cells_per_event = readOut.cells;
    cell_module_collection_types::host& modules_per_event = readOut.modules;

    auto clusters_per_event = cc(cells_per_event);
    auto measurements_per_event = mc(clusters_per_event, modules_per_event);

    assert(measurements_per_event.size() == clusters_per_event.size());

    // Fill the measurement-to-cluster map
    for (unsigned int i = 0; i < measurements_per_event.size(); ++i) {
        const auto& clus = clusters_per_event.get_items()[i];

        meas_clus_map[measurements_per_event[i]] = clus;
    }

    // Fill the barcode link map
    for (unsigned int i = 0; i < modules_per_event.size(); ++i) {
        bcd_link_map[modules_per_event[i].surface_link] = i;
    }

    /***************************************
     * Generate hit to particle map
     ***************************************/

    for (const auto& hit : hits) {
        spacepoint sp;
        sp.global = {hit.tx, hit.ty, hit.tz};

        /*
        unsigned int link = 0;
        auto it = link_map.find(hit.geometry_id);
        if (it != link_map.end()) {
            link = (*it).second;
        }
        */

        hit_ptc_map[sp] = ptc_map[hit.particle_id];
    }

    /***************************************
     * Generate hit to cluster map
     ***************************************/

    for (const auto& iocell : cells) {
        unsigned int link = 0;
        auto it =
            bcd_link_map.find(detray::geometry::barcode{iocell.geometry_id});
        if (it != bcd_link_map.end()) {
            link = (*it).second;
        }

        hit_clus_map[h_map[iocell.hit_id]].push_back(
            cell{iocell.channel0, iocell.channel1, iocell.value,
                 iocell.timestamp, link});
    }

    /***************************************
     * Generate cell to particle map
     ***************************************/

    for (auto const& [hit, ptc] : hit_ptc_map) {
        auto& clus = hit_clus_map[hit];

        for (auto& c : clus) {
            cell_ptc_map[c] = ptc;
        }
    }

    /***************************************
     * Generate measurement to particle map
     ***************************************/

    for (auto const& [meas, clus] : meas_clus_map) {
        for (const auto& c : clus) {
            meas_ptc_map[meas][cell_ptc_map[c]]++;
        }
    }
}

event_map::event_map(std::size_t event, const std::string& particle_dir,
                     const std::string& hit_dir,
                     const std::string& measurement_dir) {

    // Containers
    std::vector<traccc::io::csv::particle> particles;
    std::vector<traccc::io::csv::hit> hits;
    std::vector<traccc::io::csv::measurement> measurements;
    std::vector<traccc::io::csv::measurement_hit_id> meas_hit_ids;

    // Read Particles
    std::string io_particle_file =
        io::data_directory() + particle_dir +
        io::get_event_filename(event, "-particles.csv");

    auto preader = io::csv::make_particle_reader(io_particle_file);

    traccc::io::csv::particle io_particle;
    while (preader.read(io_particle)) {
        point3 pos{io_particle.vx, io_particle.vy, io_particle.vz};
        vector3 mom{io_particle.px, io_particle.py, io_particle.pz};

        particles.push_back(io_particle);
        ptc_map[io_particle.particle_id] =
            particle{io_particle.particle_id, io_particle.particle_type,
                     io_particle.process,     pos,
                     io_particle.vt,          mom,
                     io_particle.m,           io_particle.q};
    }

    // Read Hits
    std::string io_hit_file = io::data_directory() + hit_dir +
                              io::get_event_filename(event, "-hits.csv");

    auto hreader = io::csv::make_hit_reader(io_hit_file);

    traccc::io::csv::hit io_hit;
    while (hreader.read(io_hit)) {
        hits.push_back(io_hit);
    }

    // Read Measurements
    std::string io_measurement_file =
        io::data_directory() + measurement_dir +
        io::get_event_filename(event, "-measurements.csv");

    auto mreader = io::csv::make_measurement_reader(io_measurement_file);

    traccc::io::csv::measurement io_measurement;
    while (mreader.read(io_measurement)) {
        measurements.push_back(io_measurement);
    }

    // Read Meas-to-Hit IDs
    std::string io_meas_hit_id_file =
        io::data_directory() + hit_dir +
        io::get_event_filename(event, "-measurement-simhit-map.csv");

    auto mhid_reader =
        io::csv::make_measurement_hit_id_reader(io_meas_hit_id_file);

    traccc::io::csv::measurement_hit_id io_mh_id;
    while (mhid_reader.read(io_mh_id)) {
        meas_hit_ids.push_back(io_mh_id);
    }

    // Iterate over Measurements
    for (const auto& csv_meas : measurements) {

        // Hit index
        const auto h_id = meas_hit_ids[csv_meas.measurement_id].hit_id;

        // Make spacepoint
        const auto csv_hit = hits[h_id];
        point3 global_pos{csv_hit.tx, csv_hit.ty, csv_hit.tz};
        point3 global_mom{csv_hit.tpx, csv_hit.tpy, csv_hit.tpz};

        // Make particle
        const auto csv_ptc = particles[csv_hit.particle_id];
        point3 pos{csv_ptc.vx, csv_ptc.vy, csv_ptc.vz};
        vector3 mom{csv_ptc.px, csv_ptc.py, csv_ptc.pz};
        particle ptc{csv_ptc.particle_id, csv_ptc.particle_type,
                     csv_ptc.process,     pos,
                     csv_ptc.vt,          mom,
                     csv_ptc.m,           csv_ptc.q};

        // Make measurement
        point2 local{csv_meas.local0, csv_meas.local1};
        variance2 var{csv_meas.var_local0, csv_meas.var_local1};
        measurement meas{local, var,
                         detray::geometry::barcode{csv_meas.geometry_id}};

        // Fill measurement to truth global position and momentum map
        meas_param_map[meas] = std::make_pair(global_pos, global_mom);

        // Fill particle to measurement map
        ptc_meas_map[ptc].push_back(meas);

        // Fill measurement to particle map
        auto& contributing_particles = meas_ptc_map[meas];
        contributing_particles[ptc]++;
    }
}

}  // namespace traccc::io::experimental

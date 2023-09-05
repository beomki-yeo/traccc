/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/cell.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/particle.hpp"
#include "traccc/edm/spacepoint.hpp"
#include "traccc/edm/track_candidate.hpp"

// Detray
#include "detray/geometry/barcode.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

namespace traccc::io::experimental {

class event_map {

    public:
    /// Constructor with cells
    event_map(std::size_t event, vecmem::memory_resource& resource,
              const std::string& digi_config_file,
              const std::string& particle_dir, const std::string& hit_dir,
              const std::string& cell_dir);

    /// Constructor with measurements
    event_map(std::size_t event, const std::string& particle_dir,
              const std::string& hit_dir, const std::string& measurement_dir);

    /// Map for measurement to truth global position and momentum
    using measurement_parameter_map =
        std::map<measurement, std::pair<point3, point3>>;
    /// Map for measurement to the contributing particles
    using measurement_particle_map =
        std::map<measurement, std::map<particle, uint64_t>>;
    /// Map for particle to the vector of (geometry_id, measurement)
    using particle_measurement_map =
        std::map<particle, std::vector<measurement>>;
    /// Map for measurement to cluster (cells)
    using measurement_cluster_map = std::map<measurement, vecmem::vector<cell>>;
    /// Map for cell to particle map
    using cell_particle_map = std::map<cell, particle>;
    /// Particle ID type
    using particle_id = uint64_t;
    /// Map for particles and their IDs
    using particle_map = std::map<particle_id, particle>;
    /// Hit ID type
    using hit_id = uint64_t;
    /// Map for hits and their IDs
    using hit_map = std::map<hit_id, spacepoint>;
    /// Map for Hit to particle
    using hit_particle_map = std::map<spacepoint, particle>;
    /// Map for hit to cluster
    using hit_cluster_map = std::map<spacepoint, std::vector<cell>>;
    /// Map for barcode to link
    using barcode_link_map = std::map<detray::geometry::barcode, unsigned int>;

    particle_map ptc_map;
    hit_map h_map;
    measurement_parameter_map meas_param_map;
    measurement_particle_map meas_ptc_map;
    measurement_cluster_map meas_clus_map;
    particle_measurement_map ptc_meas_map;
    cell_particle_map cell_ptc_map;
    hit_particle_map hit_ptc_map;
    hit_cluster_map hit_clus_map;
    barcode_link_map bcd_link_map;
};

}  // namespace traccc::io::experimental

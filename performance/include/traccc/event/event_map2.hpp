/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/edm/measurement.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/edm/track_state.hpp"
#include "traccc/io/csv.hpp"
#include "traccc/io/utils.hpp"

namespace traccc {

inline bool operator==(const csv_measurement& lhs,
                       const track_state<transform3>& rhs) {
    if (lhs.geometry_id == rhs.surface_link() &&
        lhs.local0 == rhs.m_measurement.local[0] &&
        lhs.local1 == rhs.m_measurement.local[1]) {
        return true;
    }
    return false;
}

template <typename detector_type>
struct event_map2 {

    using point3 = typename detector_type::point3;
    using vector3 = typename detector_type::vector3;

    event_map2(const detector_type& det, size_t event,
               const std::string& hit_dir, const std::string& measurement_dir,
               const std::string& particle_dir)
        : m_detector(std::make_unique<detector_type>(det)) {

        std::string io_meas_hit_id_file =
            data_directory() + hit_dir +
            get_event_filename(event, "-measurement-simhit-map.csv");

        std::string io_particle_file =
            data_directory() + particle_dir +
            get_event_filename(event, "-particles.csv");

        std::string io_hit_file =
            data_directory() + hit_dir + get_event_filename(event, "-hits.csv");

        std::string io_measurement_file =
            data_directory() + measurement_dir +
            get_event_filename(event, "-measurements.csv");

        meas_hit_id_reader mhid_reader(io_meas_hit_id_file,
                                       {"measurement_id", "hit_id"});

        particle_reader preader(
            io_particle_file, {"particle_id", "particle_type", "process", "vx",
                               "vy", "vz", "vt", "px", "py", "pz", "m", "q"});

        fatras_hit_reader hreader(
            io_hit_file,
            {"particle_id", "geometry_id", "tx", "ty", "tz", "tt", "tpx", "tpy",
             "tpz", "te", "deltapx", "deltapy", "deltapz", "deltae", "index"});

        traccc::measurement_reader mreader(
            io_measurement_file,
            {"geometry_id", "local_key", "local0", "local1", "phi", "theta",
             "time", "var_local0", "var_local1", "var_phi", "var_theta",
             "var_time"});

        csv_meas_hit_id io_mh_id;
        while (mhid_reader.read(io_mh_id)) {
            m_meas_hit_ids.push_back(io_mh_id);
        }

        csv_particle io_particle;
        while (preader.read(io_particle)) {
            m_particles.push_back(io_particle);
        }

        csv_fatras_hit io_hit;
        while (hreader.read(io_hit)) {
            m_hits.push_back(io_hit);
        }

        csv_measurement io_measurement;
        while (mreader.read(io_measurement)) {
            m_measurements.push_back(io_measurement);
        }

        // Check if the size of measurements is equal to measurement-simhit-map
        assert(m_measurements.size() == m_meas_hit_ids.size());
    }

    bound_track_parameters find_truth_param(
        const track_state<transform3>& trk_state) const {

        // Find the corresponding measurement
        auto it =
            std::find(m_measurements.begin(), m_measurements.end(), trk_state);
        assert(it != m_measurements.end());

        // Measurement index
        const auto m_id = std::distance(m_measurements.begin(), it);

        // Hit index
        const auto h_id = m_meas_hit_ids[m_id].hit_id;

        // Hit objects
        const auto hit = m_hits[h_id];

        // Particle index
        const auto p_id = hit.particle_id;

        // Particle obejcts
        const auto ptc = m_particles[p_id];

        // Get truth local position
        const point3 global_pos{hit.tx, hit.ty, hit.tz};
        const vector3 global_mom{hit.tpx, hit.tpy, hit.tpz};
        const auto truth_local = m_detector.get()->global_to_local(
            hit.geometry_id, global_pos, vector::normalize(global_mom));

        // Return value
        bound_track_parameters ret;
        auto& ret_vec = ret.vector();
        getter::element(ret_vec, e_bound_loc0, 0) = truth_local[0];
        getter::element(ret_vec, e_bound_loc1, 0) = truth_local[1];
        getter::element(ret_vec, e_bound_phi, 0) = it->phi;
        getter::element(ret_vec, e_bound_theta, 0) = it->theta;
        getter::element(ret_vec, e_bound_time, 0) = it->time;
        getter::element(ret_vec, e_bound_qoverp, 0) =
            ptc.q / std::sqrt(hit.tpx * hit.tpx + hit.tpy * hit.tpy +
                              hit.tpz * hit.tpz);
        return ret;
    }

    private:
    std::unique_ptr<detector_type> m_detector;
    std::vector<csv_meas_hit_id> m_meas_hit_ids;
    std::vector<csv_particle> m_particles;
    std::vector<csv_fatras_hit> m_hits;
    std::vector<csv_measurement> m_measurements;
};

}  // namespace traccc
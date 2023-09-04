/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/utils/cell.hpp"
#include "traccc/utils/digitization_algorithm.hpp"

// Detray include(s).
#include "detray/io/csv/csv_io_types.hpp"
#include "detray/propagator/base_actor.hpp"

// DFE include(s).
#include <dfe/dfe_io_dsv.hpp>
#include <dfe/dfe_namedtuple.hpp>

namespace traccc {

struct digitization_writer : detray::actor {

    struct config {
        digitization_map digi_map;
    };

    struct state {
        state(std::size_t event_id, config&& writer_cfg,
              const std::string directory)
            : m_particle_writer(directory + io::get_event_filename(
                                                event_id, "-particles.csv")),
              m_hit_writer(directory +
                           io::get_event_filename(event_id, "-hits.csv")),
              m_cell_writer(directory +
                            io::get_event_filename(event_id, "-cells.csv")),
              m_digi_map(writer_cfg.digi_map) {}

        uint64_t particle_id = 0u;
        detray::particle_writer m_particle_writer;
        detray::hit_writer m_hit_writer;
        traccc::cell_writer m_cell_writer;
        uint64_t m_hit_count = 0u;
        digitization_map m_digi_map;

        void set_seed(const uint_fast64_t /*sd*/) {}

        void write_particle(const free_track_parameters& track) {
            detray::csv_particle particle;
            const auto pos = track.pos();
            const auto mom = track.mom();

            particle.particle_id = particle_id;
            particle.vx = pos[0];
            particle.vy = pos[1];
            particle.vz = pos[2];
            particle.vt = track.time();
            particle.px = mom[0];
            particle.py = mom[1];
            particle.pz = mom[2];
            particle.q = track.charge();

            m_particle_writer.append(particle);
        }
    };

    template <typename propagator_state_t>
    void operator()(state& writer_state,
                    propagator_state_t& propagation) const {

        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // triggered only for sensitive surfaces
        if (navigation.is_on_sensitive()) {

            // Write hits
            detray::csv_hit hit;

            const auto track = stepping();
            const auto pos = track.pos();
            const auto mom = track.mom();

            const auto sf = navigation.get_surface();

            hit.particle_id = writer_state.particle_id;
            hit.geometry_id = sf.barcode().value();
            hit.tx = pos[0];
            hit.ty = pos[1];
            hit.tz = pos[2];
            hit.tt = track.time();
            hit.tpx = mom[0];
            hit.tpy = mom[1];
            hit.tpz = mom[2];

            writer_state.m_hit_writer.append(hit);

            /*
            // Write measurements
            // NOTE: Works for 2D measurement only
            csv_measurement meas;

            const auto bound_params = stepping._bound_params;

            const auto local = sf.template visit_mask<measurement_kernel>(
                bound_params, writer_state.m_meas_smearer);

            meas.measurement_id = writer_state.m_hit_count;
            meas.geometry_id = hit.geometry_id;
            meas.local_key = "unknown";
            meas.local0 = local[0];
            meas.local1 = local[1];
            auto stddev_0 = writer_state.m_meas_smearer.stddev[0];
            auto stddev_1 = writer_state.m_meas_smearer.stddev[1];
            meas.var_local0 = stddev_0 * stddev_0;
            meas.var_local1 = stddev_1 * stddev_1;
            meas.phi = bound_params.phi();
            meas.theta = bound_params.theta();
            meas.time = bound_params.time();

            writer_state.m_meas_writer.append(meas);
            */
           
            // Write cells 
            // NOTE: Works for 2D measurement only
            digitization_algorithm digi_alg(navigation.detector(),
                                            writer_state.m_digi_map);

            const auto cells = digi_alg(stepping._bound_params);

            for (const auto& cell : cells) {
                traccc::csv_cell ce;
                ce.geometry_id = hit.geometry_id;
                ce.hit_id = writer_state.m_hit_count;
                ce.channel0 = cell[0];
                ce.channel1 = cell[1];
                ce.value = 1.f;
                writer_state.m_cell_writer.append(ce);
            }

            writer_state.m_hit_count++;
        }
    }
};

}  // namespace traccc

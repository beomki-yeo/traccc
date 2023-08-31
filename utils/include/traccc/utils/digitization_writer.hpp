/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/utils/digitization_algorithm.hpp"

// Detray include(s).
#include "detray/propagator/base_actor.hpp"

namespace traccc {

struct digitization_writer : detray::actor {

    struct state {
        std::unique_ptr<digitization_map> m_digi_map;
    };

    template <typename propagator_state_t>
    void operator()(state& writer_state,
                    propagator_state_t& propagation) const {

        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // triggered only for sensitive surfaces
        if (navigation.is_on_sensitive()) {

            // Write hits
            csv_hit hit;

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

            // Write cells
            digitization_algorithm digi_alg(navigation.detector(),
                                            writer_state.m_digi_map);

            // csv_cell ce;
            /*
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

            // Write hit measurement map
            csv_meas_hit_id meas_hit_id;
            meas_hit_id.hit_id = writer_state.m_hit_count;
            meas_hit_id.measurement_id = writer_state.m_hit_count;
            writer_state.m_meas_hit_id_writer.append(meas_hit_id);
            writer_state.m_hit_count++;
            */
        }
    }
};

}  // namespace traccc

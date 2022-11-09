/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/edm/track_parameters.hpp"

// detray include(s).
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/propagator.hpp"

// System include(s).
#include <random>

namespace traccc {

using matrix_operator = typename transform3::matrix_actor;

/// Track parameter smearer
struct parameter_smearer {

    std::random_device rd{};
    std::mt19937 generator{rd()};

    bound_track_parameters operator()(
        bound_track_parameters& param,
        const std::array<scalar, e_bound_size>& stddevs) {
        // New vector
        auto new_vec = matrix_operator().template zero<e_bound_size, 1>();

        // New covariance
        auto new_cov =
            matrix_operator().template zero<e_bound_size, e_bound_size>();

        for (std::size_t i = 0; i < e_bound_size; i++) {

            matrix_operator().element(new_vec, i, 0) =
                std::normal_distribution<scalar>(
                    matrix_operator().element(param.vector(), i, 0),
                    stddevs[i])(generator);

            matrix_operator().element(new_cov, i, i) = stddevs[i] * stddevs[i];
        }

        // Set vector and covariance
        param.set_vector(new_vec);
        param.set_covariance(new_cov);

        return param;
    }

    free_track_parameters operator()(
        free_track_parameters& param,
        const std::array<scalar, e_free_size>& stddevs) {
        // New vector
        auto new_vec = matrix_operator().template zero<e_free_size, 1>();

        // New covariance
        auto new_cov =
            matrix_operator().template zero<e_free_size, e_free_size>();

        for (std::size_t i = 0; i < e_free_size; i++) {

            matrix_operator().element(new_vec, i, 0) =
                std::normal_distribution<scalar>(
                    matrix_operator().element(param.vector(), i, 0),
                    stddevs[i])(generator);

            matrix_operator().element(new_cov, i, i) = stddevs[i] * stddevs[i];
        }

        // Set vector and covariance
        param.set_vector(new_vec);
        param.set_covariance(new_cov);

        // Normalize dir
        const auto dir = param.dir();
        param.set_dir(vector::normalize(dir));

        return param;
    }
};

template <typename stepper_t, typename navigator_t>
struct seed_generator {

    struct aborter : detray::actor {
        struct state {};

        template <typename propagator_state_t>
        void operator()(state& /*abrt_state*/,
                        propagator_state_t& propagation) const {

            auto& navigation = propagation._navigation;

            //// DELETE ME ////////
            /*
            auto& stepping = propagation._stepping;
            printf("Propagating... ");
            for (std::size_t i = 0; i < e_free_size; i++) {
                printf("%f ", getter::element(stepping().vector(), i, 0));
            }
            printf("\n");
            */
            //// DELETE ME ////////

            if (navigation.is_on_module()) {
                propagation._heartbeat &= navigation.abort();
            }
        }
    };

    using transform3_type = typename stepper_t::transform3_type;
    using detector_type = typename navigator_t::detector_type;

    using transporter = detray::parameter_transporter<transform3_type>;
    using resetter = detray::parameter_resetter<transform3_type>;
    using actor_chain_type =
        detray::actor_chain<std::tuple, transporter, resetter, aborter>;
    using propagator_type =
        detray::propagator<stepper_t, navigator_t, actor_chain_type>;

    seed_generator(const detector_type& det)
        : m_detector(std::make_unique<detector_type>(det)) {}

    bound_track_parameters operator()(
        const free_track_parameters& vertex,
        const std::array<scalar, e_bound_size>& stddevs) {
        propagator_type propagator({}, {});
        typename actor_chain_type::state actor_states =
            std::tie(m_transporter_state, m_resetter_state, m_aborter_state);
        typename propagator_type::state propagation(
            vertex, m_detector->get_bfield(), *m_detector, actor_states);

        auto& stepping = propagation._stepping;

        //// DELETE ME ////////
        /*
        printf("Initial \n");
        for (std::size_t i = 0; i < e_free_size; i++) {
            printf("%f ", getter::element(stepping().vector(), i, 0));
        }
        printf("\n");
        */
        propagator.propagate(propagation);
        /*
        printf("Seed Parameter \n");

        for (std::size_t i = 0; i < e_free_size; i++) {
            printf("%f ", getter::element(stepping().vector(), i, 0));
        }
        printf("\n");

        printf("%lu ", stepping._bound_params.surface_link());
        for (std::size_t i = 0; i < e_bound_size; i++) {
            printf("%f ",
                   getter::element(stepping._bound_params.vector(), i, 0));
        }
        printf("\n");

        for (std::size_t i = 0; i < e_bound_size; i++) {
            for (std::size_t j = 0; j < e_bound_size; j++) {
                printf("%f ", getter::element(
                                  stepping._bound_params.covariance(), i, j));
            }
            printf("\n");
        }
        printf("\n");
        */
        //// DELETE ME ////////

        return parameter_smearer()(stepping._bound_params, stddevs);
    }

    std::random_device rd{};
    std::mt19937 generator{rd()};

    std::unique_ptr<detector_type> m_detector;
    typename transporter::state m_transporter_state{};
    typename resetter::state m_resetter_state{};
    typename aborter::state m_aborter_state{};
};

}  // namespace traccc
/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/track_parametrization.hpp"

namespace traccc {

template <typename algebra_t>
struct gain_matrix_updater {

    using output_type = bool;
    using matrix_operator = typename algebra_t::matrix_actor;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;

    template <typename mask_group_t, typename index_t, typename surface_t,
              typename track_state_t, typename propagator_state_t>
    TRACCC_HOST_DEVICE inline output_type operator()(
        const mask_group_t& mask_group, const index_t& /*index*/,
        const surface_t& surface, track_state_t& trk_state,
        propagator_state_t& propagation) const {

        // Reference: Application of Kalman filtering to track and vertex
        // fitting, R.Frühwirth, NIM A
        auto& stepping = propagation._stepping;

        // Mask associated with the track state
        const auto& mask = mask_group[surface.mask().index()];

        // Some identity matrices
        const auto I66 =
            matrix_operator().template identity<e_bound_size, e_bound_size>();
        const auto I22 = matrix_operator().template identity<2, 2>();

        // projection matrix
        // Dimension: (2 X 6)
        const auto H = mask.template projection_matrix<e_bound_size>();

        // Measurement data on surface
        // Dimension: (2 X 1)
        const auto& meas_local = trk_state.measurement_local();

        // Predicted vector of bound track parameters
        // Dimension: (6 X 1)
        const auto& predicted_vec = stepping._bound_params.vector();

        // Predicted covaraince of bound track parameters
        // Dimension: (6 X 6)
        const auto& predicted_cov = stepping._bound_params.covariance();

        // Set track state parameters
        trk_state.predicted().set_vector(predicted_vec);
        trk_state.predicted().set_covariance(predicted_cov);

        // Spatial resolution (Measurement covariance)
        // Dimension: (2 X 2)
        const auto V = trk_state.measurement_covariance();

        // Dimension: (2 X 6) * (6 X 6) * (6 X 2)
        const matrix_type<2, 2> M =
            H * predicted_cov * matrix_operator().transpose(H) + V;

        // Kalman gain matrix
        // Dimension: (6 X 6) * (6 X 2) * (2 X 2)
        const auto K = predicted_cov * matrix_operator().transpose(H) *
                       matrix_operator().inverse(M);

        // Dimension: (6 X 1) + (6 X 2) * [ (2 X 1) - (2 X 6) * (6 X 1) ]
        const auto filtered_vec =
            predicted_vec + K * (meas_local - H * predicted_vec);

        // Dimension: [ (6 X 6) - (6 X 2) * (2 X 6) ] * (6 X 6)
        const auto filtered_cov = (I66 - K * H) * predicted_cov;

        // Residual between measurement and (projected) filtered vector
        // Dimension: (2 X 1) - (2 X 6) * (6 X 1);
        const matrix_type<2, 1> residual = meas_local - H * filtered_vec;

        // Calculate the chi square
        // Dimension: [ (2 X 2) - (2 X 6) * (6 X 2) ] * (2 X 2)
        const matrix_type<2, 2> R = (I22 - H * K) * V;
        // Dimension: (1 X 2) * ( 2 X 2 ) * (2 X 1)
        const matrix_type<1, 1> chi2 = matrix_operator().transpose(residual) *
                                       matrix_operator().inverse(R) * residual;

        // Set Stepper parameter
        stepping._bound_params.set_vector(filtered_vec);
        stepping._bound_params.set_covariance(filtered_cov);
        // Set track state parameters
        trk_state.filtered().set_vector(filtered_vec);
        trk_state.filtered().set_covariance(filtered_cov);
        trk_state.filtered_chi2() = matrix_operator().element(chi2, 0, 0);

        //////////////////// Delete ME /////////////////////////
        /*
        printf("Surface: %lu \n", trk_state.surface_link());

        printf("Predicted: ");
        for (auto i = 0; i < e_bound_size; i++){
            printf("%f ", getter::element(trk_state.predicted().vector(),i,0));
        }
        printf("\n");

        printf("Filtered: ");
        for (auto i = 0; i < e_bound_size; i++){
            printf("%f ", getter::element(trk_state.filtered().vector(),i,0));
        }
        printf("Filtered Chi2: %f", trk_state.filtered_chi2());

        printf("\n");
        */

        return true;
    }
};

}  // namespace traccc

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/track_parameters.hpp"

// detray include(s).
#include "detray/propagator/navigator.hpp"

namespace traccc {

/// Type unrolling functor to smooth the track parameters after the Kalman
/// filtering
template <typename algebra_t>
struct gain_matrix_smoother {

    // Type declarations
    using output_type = bool;
    using matrix_operator = typename algebra_t::matrix_actor;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;

    /// Gain matrix smoother operation
    ///
    /// @brief Based on "Application of Kalman filtering to track and vertex
    /// fitting", R.Frühwirth, NIM A
    ///
    /// @param mask_group mask group that contains the mask of the current
    /// surface
    /// @param index mask index of the current surface
    /// @param cur_state track state of the current surface
    /// @param next_state track state of the next surface
    ///
    /// @return true if the update succeeds
    template <typename mask_group_t, typename index_t>
    TRACCC_HOST_DEVICE inline output_type operator()(
        const mask_group_t& mask_group, const index_t& index,
        track_state<algebra_t>& cur_state,
        const track_state<algebra_t>& next_state) {

        const auto& next_smoothed = next_state.smoothed();
        const auto& next_predicted = next_state.predicted();
        const auto& cur_filtered = cur_state.filtered();

        // Next track state parameters
        const auto& next_jacobian = next_state.jacobian();
        const auto& next_smoothed_vec = next_smoothed.vector();
        const auto& next_smoothed_cov = next_smoothed.covariance();
        const auto& next_predicted_vec = next_predicted.vector();
        const auto& next_predicted_cov = next_predicted.covariance();

        // Current track state parameters
        const auto& cur_filtered_vec = cur_filtered.vector();
        const auto& cur_filtered_cov = cur_filtered.covariance();

        // Regularization matrix for numerical stability
        static constexpr double epsilon = 1e-13;
        const matrix_type<6, 6> regularization =
            matrix_operator().template identity<e_bound_size, e_bound_size>() *
            epsilon;
        const matrix_type<6, 6> regularized_predicted_cov =
            next_predicted_cov + regularization;

        // Calculate smoothed parameter for current state
        // Dimension: (6 X 6) * (6 X 6) * (6 X 6)
        const matrix_type<6, 6> A =
            cur_filtered_cov * matrix_operator().transpose(next_jacobian) *
            matrix_operator().inverse(regularized_predicted_cov);
        const auto smt_vec =
            cur_filtered_vec + A * (next_smoothed_vec - next_predicted_vec);
        const auto smt_cov =
            cur_filtered_cov + A * (next_smoothed_cov - next_predicted_cov) *
                                   matrix_operator().transpose(A);

        cur_state.smoothed().set_vector(smt_vec);
        cur_state.smoothed().set_covariance(smt_cov);

        //////////////////// Delete ME /////////////////////////
        /*
        printf("Surface: %lu \n", cur_state.surface_link());

        printf("Predicted: ");
        for (auto i = 0; i < e_bound_size; i++){
            printf("%f ", getter::element(cur_state.predicted().vector(),i,0));
        }
        printf("\n");

        printf("Filtered: ");
        for (auto i = 0; i < e_bound_size; i++){
            printf("%f ", getter::element(cur_filtered_vec,i,0));
        }
        printf("\n");

        printf("Smoothed: ");
        for (auto i = 0; i < e_bound_size; i++){
            printf("%f ", getter::element(smt_vec,i,0));
        }
        printf("\n");
        */

        // projection matrix
        // Dimension: (2 X 6)
        const auto H =
            mask_group[index].template projection_matrix<e_bound_size>();

        // Calculate smoothed chi square
        const auto& meas_local = cur_state.measurement_local();
        const auto& V = cur_state.measurement_covariance();

        // Dimension: (2 X 1) - (2 X 6) * (6 X 1);
        const matrix_type<2, 1> residual = meas_local - H * smt_vec;
        // Dimension: (2 X 2) - (2 X 6) * (6 X 6) * (6 X 2)
        const matrix_type<2, 2> R =
            V - H * smt_cov * matrix_operator().transpose(H);

        // Dimension: (1 X 2) * ( 2 X 2 ) * (2 X 1)
        const matrix_type<1, 1> chi2 = matrix_operator().transpose(residual) *
                                       matrix_operator().inverse(R) * residual;

        cur_state.smoothed_chi2() = matrix_operator().element(chi2, 0, 0);

        return true;
    }

    // Reference: Eq. (3.38 - 3.40) of "Pattern Recognition, Tracking and Vertex
    // Reconstruction in Particle Detectors" from R.Frühwirth and A. Strandlie
    // NOTE: Kalman fitter of ACTS just uses the backward filtering parameters
    // as a final result. This smoothering might be very expensive, therefore,
    // we need to decide if we are going to use the smoothering or not.
    /*
    template <typename mask_group_t, typename index_t, typename surface_t>
    TRACCC_HOST_DEVICE inline output_type operator()(
        const mask_group_t& mask_group, const index_t& index,
        const surface_t surface, track_state<algebra_t>& trk_state) {

        const auto param_1 = trk_state.filtered();

        const auto param_2 = trk_state.backward_predicted();

        // Todo: flip param_2

        // Mask associated with the track state
        const auto& mask = mask_group[surface.mask().index()];

        // projection matrix
        // Dimension: (2 X 6)
        const auto H = mask.template projection_matrix<e_bound_size>();

        // Get vectors and covarainces (and their inverses)
        const auto& param_vec_1 = param_1.vector();
        const auto& param_cov_1 = param_1.covariance();
        const auto& param_vec_2 = param_2.vector();
        const auto& param_cov_2 = param_2.covariance();

        const auto param_cov_inv_1 = matrix_operator().inverse(param_cov_1);

        const auto param_cov_inv_2 = matrix_operator().inverse(param_cov_2);

        // Calculate Eq. 3.38
        const auto smt_cov_inv = param_cov_inv_1 + param_cov_inv_2;
        const auto smt_cov = matrix_operator().inverse(smt_cov_inv);

        // Calculate Eq. 3.38 for smoothed vector
        const auto smt_vec = smt_cov * (param_cov_inv_1 * param_vec_1 +
                                        param_cov_inv_2 * param_vec_2);

        // Spatial resolution (Measurement covariance)
        // Dimension: (2 X 2)
        const auto V = trk_state.measurement_covariance();

        // Calculate Eq. 3.39 for smoothed covariance
        // Dimension: (2 X 2) - (2 X 6) * (6 X 6) * (6 * 2)
        const auto R = V - H * smt_cov * matrix_operator().transpose(H);

        // Measurement data on surface
        // Dimension: (2 X 1)
        const auto& meas_local = trk_state.measurement_local();

        // Calculate Eq. 3.39 for residual
        // Dimension: (2 X 1) - (2 X 6) * (6 X 1);
        const auto residual = meas_local - H * smt_vec;

        // Calculate Eq. 3.40 for chi square
        // Dimension: (1 X 2) * ( 2 X 2 ) * (2 X 1)
        const auto chi2 = matrix_operator().transpose(residual) *
                          matrix_operator().inverse(R) * residual;

        trk_state.smoothed() =
            bound_track_parameters(trk_state.surface_link(), smt_vec, smt_cov);
        trk_state.smoothed_chi2() = matrix_operator().element(chi2, 0, 0);

        return true;
    }
    */
};

}  // namespace traccc
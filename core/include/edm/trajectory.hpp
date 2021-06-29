/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// traccc include(s).
#include "definitions/algebra.hpp"
#include "definitions/primitives.hpp"
#include "edm/container.hpp"

// VecMem include(s).
#include <vecmem/containers/data/jagged_vector_buffer.hpp>
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>
#include <vecmem/containers/jagged_vector.hpp>
#include <vecmem/containers/vector.hpp>

// standard
#include <bitset>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <vector>

// Acts
#include "Acts/EventData/Measurement.hpp"

// Eigen
#include <Eigen/Core>

namespace traccc {

/// Either type T or const T depending on the boolean.
template <typename T, bool select>
using ConstIf = std::conditional_t<select, const T, T>;

/// Type construction helper for coefficients and associated covariances.
template <size_t Size>
struct Types {
    enum {
        Flags = Eigen::ColMajor | Eigen::AutoAlign,
        SizeIncrement = 8,
    };
    using Scalar = double;
    // single items
    using Coefficients = Eigen::Matrix<Scalar, Size, 1, Flags>;
    using Covariance = Eigen::Matrix<Scalar, Size, Size, Flags>;
};

struct trajectory {
    enum {
        MeasurementSizeMax = Acts::eBoundSize,
        // 25 for the maximum number of measurements in trackML detector
        NumMeasurementsMax = 25,
    };
    unsigned int n_measurements;
    unsigned int seed_idx;
    array<Types<Acts::eBoundSize>::Coefficients, NumMeasurementsMax> m_params;
    array<Types<Acts::eBoundSize>::Covariance, NumMeasurementsMax> m_cov;
    array<Types<MeasurementSizeMax>::Coefficients, NumMeasurementsMax> m_meas;
    array<Types<MeasurementSizeMax>::Covariance, NumMeasurementsMax> m_meas_cov;
};

/// Container of trajectorys belonging to one detector module
template <template <typename> class vector_t>
using trajectory_collection = vector_t<trajectory>;

/// Convenience declaration for the trajectory collection type to use in host
/// code
using host_trajectory_collection = trajectory_collection<vecmem::vector>;
/// Convenience declaration for the trajectory collection type to use in device
/// code
using device_trajectory_collection =
    trajectory_collection<vecmem::device_vector>;

/// Convenience declaration for the trajectory container type to use in host
/// code
using host_trajectory_container = host_container<unsigned int, trajectory>;

/// Convenience declaration for the trajectory container type to use in device
/// code
using device_trajectory_container = device_container<unsigned int, trajectory>;

/// Convenience declaration for the trajectory container data type to use in
/// host code
using trajectory_container_data = container_data<unsigned int, trajectory>;

/// Convenience declaration for the trajectory container buffer type to use in
/// host code
using trajectory_container_buffer = container_buffer<unsigned int, trajectory>;

/// Convenience declaration for the trajectory container view type to use in
/// host code
using trajectory_container_view = container_view<unsigned int, trajectory>;

}  // namespace traccc

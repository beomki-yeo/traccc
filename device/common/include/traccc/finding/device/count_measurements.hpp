/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/finding/device/make_module_map.hpp"

namespace traccc::device {

// @TODO: usu detector_type as a template parameter
template <typename detector_t>
TRACCC_DEVICE inline void count_measurements(
    std::size_t globalIndex, typename detector_t::detector_view_type det_data,
    measurement_container_types::const_view measurements_view,
    vecmem::data::vector_view<thrust::pair<unsigned int, unsigned int>>
        module_map_view,
    const int n_params,
    bound_track_parameters_collection_types::view params_view,
    vecmem::data::vector_view<unsigned int> n_measurements_view,
    unsigned int& n_total_measurements);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/count_measurements.ipp"
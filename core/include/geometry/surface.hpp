/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "definitions/primitives.hpp"
#include "edm/collection.hpp"

// used surface class of gpuKalmanFilter

namespace traccc {

template < typename surface_t >    
struct surface_link{
    geometry_id geometry;
    surface_t surface;
};
    
/// Convenience declaration for the surface collection type to use in host
/// code
template <typename surface_t>    
using host_surface_collection = host_collection<surface_link<surface_t>>;

/// Convenience declaration for the surface collection type to use in device
/// code
template <typename surface_t>      
using device_surface_collection = device_collection<surface_link<surface_t>>;
    
} // namespace traccc

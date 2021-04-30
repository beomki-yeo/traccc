/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "edm/cell.hpp"
#include "edm/measurement.hpp"
#include "edm/spacepoint.hpp"
#include "detail/label.hpp"

namespace traccc {
namespace cuda{

    /// Wrapper function for counting the number of measurements per event
    ///
    /// @param ms_container: The input measurements of an event
    ///
    /// @param sp_container: The output spacepoints of an event
    ///
    /// @param resource: memery resource for cuda
    ///
    /// The global positions and errors of measurements are obtained
    void spacepoint_formation(host_measurement_container& ms_container,
			      host_spacepoint_container& sp_container,
			      vecmem::memory_resource* resource);
    
}
}

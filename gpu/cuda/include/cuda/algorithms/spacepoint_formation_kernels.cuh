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
    
    void spacepoint_formation(host_measurement_container& ms_container,
			      host_spacepoint_container& sp_container,
			      vecmem::memory_resource* resource);
    
}
}

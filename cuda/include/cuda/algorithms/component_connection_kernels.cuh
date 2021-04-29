/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "edm/cell.hpp"
#include "detail/label.hpp"

namespace traccc {
namespace cuda{
    
    void component_connection(host_cell_container& cells_per_event,
			      detail::host_label_container& labels_per_event,
			      vecmem::memory_resource* resource);

}
}

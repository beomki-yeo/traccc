/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "edm/cell.hpp"
#include "detail/label.hpp"
#include "definitions/algebra.hpp"
#include "definitions/primitives.hpp"

namespace traccc {
namespace cuda{
    void check_valid_label(host_cell_container& cells_per_event,
			   detail::host_label_container& labels_per_event,
			   vecmem::jagged_vector< bool >& valid_labels,
			   vecmem::memory_resource* resource);
}
}

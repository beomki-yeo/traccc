/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "edm/cell.hpp"
#include "edm/measurement.hpp"
#include "detail/label.hpp"

namespace traccc {
namespace cuda{
    
    void count_measurements(host_cell_container& cells,
			    detail::host_label_container& cc_labels,
			    detail::host_label_container& ms_labels,
			    vecmem::memory_resource* resource);

    void measurement_creation(host_cell_container& ce_container,
			      detail::host_label_container& cc_labels,
			      detail::host_label_container& ms_labels,
			      host_measurement_container& ms_container,
			      vecmem::memory_resource* resource);
}
}

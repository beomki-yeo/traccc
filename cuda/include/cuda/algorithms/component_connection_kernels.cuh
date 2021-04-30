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

    /// Wrapper function for the connected component kernel function per event
    ///
    /// @param cells_per_event: The input cells of an event
    /// @param labels_per_event: The labeling result from connected component.
    ///                          count the number of unique labels and record labels
    ///                          for all cells
    /// @param resource: memery resource for cuda

    /// labeling example result for a module
    ///        < c1, c2, c3, c4, c5, c6, c7, c8 >   : cells_per_module  --- input
    ///  { 4 , <  2,  1,  3,  1,  1,  1,  2,  4 > } : labels_per_module --- output
    ///   --- -----------------------------------
    ///  # of         labels for each cell
    ///  unique
    ///  labels
    void component_connection(host_cell_container& cells_per_event,
			      detail::host_label_container& labels_per_event,
			      vecmem::memory_resource* resource);

}
}

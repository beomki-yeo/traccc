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
#include "geometry/pixel_segmentation.hpp"

namespace traccc {
namespace cuda{
    
    /// Wrapper function for counting the number of measurements per event
    ///
    /// @param cells: The input cells of an event
    /// @param cc_labels: The labeling result from connected component.
    ///                   Count the number of unique labels and record labels
    ///                   for all cells
    /// @param ms_labels: The measurement couting result.
    ///                   Count the number of measurements and record valid labels
    /// @param resource: memery resource for cuda

    /// counting result for a module
    ///        < c1, c2, c3, c4, c5, c6, c7, c8 >   : cells     --- input
    ///  { 4 , <  2,  1,  3,  1,  1,  1,  2,  4 > } : cc_labels --- input
    ///   --- -----------------------------------
    ///  # of         labels for each cell
    ///  unique
    ///  labels
    ///
    ///  { 4 , <  1,  2,  3,  4,  0,  0,  0,  0 > } : ms_labels --- output
    ///   ---  -----------------  -------------
    ///  # of    valid labels        not used
    ///  meas-
    ///  ments
    ///
    ///
    ///  # of measurements can be smaller than # of uniqe labels
    ///  if sum of weights of clusters is less than threshold.
    
    void count_measurements(host_cell_container& cells,
			    detail::host_label_container& cc_labels,
			    detail::host_label_container& ms_labels,
			    vecmem::memory_resource* resource);

    /// Wrapper function for measurement creation per event
    ///
    /// @param cells: The input cells of an event
    /// @param cc_labels: The labeling result from connected component.
    ///                   Count the number of unique labels and record labels
    ///                   for all cells
    /// @param ms_labels: The measurement couting result.
    ///                   Count the number of measurements and record valid labels
    /// @param ms_container: The output measurements of an event
    /// @param resource: memery resource for cuda

    /// The size of ms_continaer is based on number of_measurements.
    /// In the kernel function, the weighted averages of local position and erros
    /// are calculated.
    void measurement_creation(host_cell_container& ce_container,
			      detail::host_label_container& cc_labels,
			      detail::host_label_container& ms_labels,
			      host_measurement_container& ms_container,
			      vecmem::memory_resource* resource);
}
}

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/io/data_format.hpp"
#include "traccc/io/reader_edm.hpp"

// Project include(s).
#include "traccc/edm/cell.hpp"
#include "traccc/io/digitization_config.hpp"

// System include(s).
#include <cstddef>
#include <string_view>

namespace traccc::io::experimental {

/// Read cell data into memory
///
/// The file to read is selected according the naming conventions used in
/// our data.
///
/// @param out A cell & a cell_module (host) collections
/// @param event The event ID to read in the cells for
/// @param directory The directory holding the cell data files
/// @param format The format of the cell data files (to read)
/// @param dconfig The detector's digitization configuration
///
void read_cells(cell_reader_output &out, std::size_t event,
                std::string_view directory,
                data_format format = data_format::csv,
                const digitization_map *dconfig = nullptr);

/// Read cell data into memory
///
/// The file name is selected explicitly by the user.
///
/// @param out A cell & a cell_module (host) collections
/// @param filename The file to read the cell data from
/// @param format The format of the cell data files (to read)
/// @param dconfig The detector's digitization configuration
///
void read_cells(cell_reader_output &out, std::string_view filename,
                data_format format = data_format::csv,
                const digitization_map *dconfig = nullptr);

}  // namespace traccc::io::experimental

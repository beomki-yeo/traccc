/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/io/experimental/csv/hit.hpp"

// DFE include(s).
#include <dfe/dfe_io_dsv.hpp>

// System include(s).
#include <string_view>

namespace traccc::io::experimental::csv {

/// Set up an object for reading a CSV file containing hit information
///
/// @param filename The name of the file to read
/// @return An object that can read the specified CSV file
///
dfe::NamedTupleCsvReader<hit> make_hit_reader(std::string_view filename);

}  // namespace traccc::io::experimental::csv

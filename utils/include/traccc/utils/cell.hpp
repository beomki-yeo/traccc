/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// DFE include(s).
#include <dfe/dfe_io_dsv.hpp>
#include <dfe/dfe_namedtuple.hpp>

// System include(s).
#include <cstdint>

namespace traccc {

/// Type used in reading CSV data into memory
struct csv_cell {

    uint64_t geometry_id = 0;
    uint64_t hit_id = 0;
    uint32_t channel0 = 0;
    uint32_t channel1 = 0;
    float timestamp = 0.;
    float value = 0.;

    // geometry_id,hit_id,channel0,channel1,timestamp,value
    DFE_NAMEDTUPLE(csv_cell, geometry_id, hit_id, channel0, channel1, timestamp,
                   value);
};

using cell_writer = dfe::NamedTupleCsvWriter<csv_cell>;

}  // namespace traccc

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/cuda/seeding/seeding_algorithm.hpp"

// Project include(s).
#include "traccc/seeding/detail/seeding_config.hpp"

// System include(s).
#include <cmath>

namespace traccc::cuda {

seeding_algorithm::seeding_algorithm(const seedfinder_config& finder_config,
                                     const spacepoint_grid_config& grid_config,
                                     const seedfilter_config& filter_config,
                                     const traccc::memory_resource& mr,
                                     vecmem::copy& copy, stream& str)
    : m_spacepoint_binning(finder_config, grid_config, mr, copy, str),
      m_seed_finding(finder_config, filter_config, mr, copy, str) {}

seeding_algorithm::output_type seeding_algorithm::operator()(
    const spacepoint_collection_types::const_view& spacepoints_view) const {

    /// Sub-algorithm performing the spacepoint binning
    spacepoint_binning binning_alg(m_finder_config, m_grid_config, m_mr, m_copy,
                                   m_stream);
    /// Sub-algorithm performing the seed finding
    seed_finding finding_alg(m_finder_config, m_filter_config, m_mr, m_copy,
                             m_stream);

    return finding_alg(spacepoints_view, binning_alg(spacepoints_view));
}

}  // namespace traccc::cuda

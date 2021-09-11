/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <algorithm>
#include <edm/internal_spacepoint.hpp>
#include <edm/spacepoint.hpp>
#include <seeding/detail/seeding_config.hpp>
#include <seeding/detail/spacepoint_grid.hpp>
#include "utils/algorithm.hpp"

namespace traccc {
    
struct spacepoint_binning
    : public algorithm<const host_spacepoint_container&,
		       phi_z_grid>{

    spacepoint_binning(const seedfinder_config& config,
		       const spacepoint_grid_config& grid_config,
		       vecmem::memory_resource& mr)
        : m_config(config), m_grid_config(grid_config), m_mr(mr){


    }

private:
    seedfinder_config m_config;
    spacepoint_grid_config m_grid_config;
    vecmem::memory_resource& m_mr;
};
    
} // namespace traccc

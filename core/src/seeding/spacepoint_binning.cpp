/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/seeding/spacepoint_binning.hpp"

#include "traccc/definitions/primitives.hpp"
#include "traccc/seeding/spacepoint_binning_helper.hpp"

namespace traccc {

spacepoint_binning::spacepoint_binning(
    const seedfinder_config& config, const spacepoint_grid_config& grid_config,
    vecmem::memory_resource& mr)
    : m_config(config),
      m_grid_config(grid_config),
      m_axes(get_axes(grid_config, mr)),
      m_mr(mr) {}

spacepoint_binning::output_type spacepoint_binning::operator()(
    const spacepoint_container_types::const_view& spacepoints_view) const {

    output_type g2(m_axes.first, m_axes.second, m_mr.get());

    const spacepoint_container_types::const_device spacepoints(
        spacepoints_view);

    djagged_vector<sp_location> rbins(m_config.get_num_rbins());

    for (unsigned int i = 0; i < spacepoints.size(); i++) {
        for (unsigned int j = 0; j < spacepoints.get_items()[i].size(); j++) {
            sp_location sp_loc{i, j};
            fill_radius_bins<spacepoint_container_types::const_device,
                             djagged_vector>(m_config, spacepoints, sp_loc,
                                             rbins);
        }
    }

    // fill rbins into grid such that each grid bin is sorted in r
    // space points with delta r < rbin size can be out of order
    for (auto& rbin : rbins) {
        for (auto& sp_loc : rbin) {

            auto isp = internal_spacepoint<spacepoint>(
                spacepoints, {sp_loc.bin_idx, sp_loc.sp_idx}, m_config.beamPos);

            point2 sp_position = {isp.phi(), isp.z()};
            g2.populate(sp_position, std::move(isp));
        }
    }

    return g2;
}

}  // namespace traccc

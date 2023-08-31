/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/cell.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/io/digitization_config.hpp"

// System include(s).
#include <array>

namespace traccc {

template <typename detector_t>
struct digitization_algorithm {

    digitization_algorithm(const detector_t det,
                           const digitization_map& digi_cfg)
        : m_detector(std::make_unique<detector_t>(det)),
          m_digi_map(std::make_unique<digitization_map>(digi_cfg)){};

    std::vector<std::array<channel_id, 2>> operator()(
        const traccc::bound_track_parameters& bound_param) const {

        std::vector<std::array<channel_id, 2>> ret;

        const auto vid = bound_param.surface_link().volume();
        const auto segmentation = (*m_digi_map)[vid];

        // Binning boundaries
        std::size_t n_bins0 = segmentation.binningData()[0].bins();
        std::size_t n_bins1 = segmentation.binningData()[1].bins();

        const auto local_tmp = bound_param.bound_local();
        Acts::Vector2 local{local_tmp[0], local_tmp[1]};

        const auto bin0 = static_cast<channel_id>(segmentation.bin(local, 0));
        const auto bin1 = static_cast<channel_id>(segmentation.bin(local, 1));

        // Generate five cells (center, left, right, up, and down)
        ret.push_back({bin0, bin1});
        if (bin0 > 0u) {
            ret.push_back({bin0 - 1u, bin1});
        }
        if (bin1 > 0u) {
            ret.push_back({bin0, bin1 - 1u});
        }
        if (bin0 < n_bins0 - 1u) {
            ret.push_back({bin0 + 1u, bin1});
        }
        if (bin1 < n_bins1 - 1u) {
            ret.push_back({bin0, bin1 + 1u});
        }

        return ret;
    };

    std::unique_ptr<detector_t> m_detector;
    std::unique_ptr<digitization_map> m_digi_map;
};

}  // namespace traccc
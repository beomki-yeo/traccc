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

namespace traccc {

template <typename detector_t>
struct digitization_algorithm {

    digitization_algorithm(const detector_t det,
                           const digitization_map& digi_cfg)
        : m_detector(std::make_unique<detector_t>(det)),
          m_digi_cfg(std::make_unique<digitization_map>(digi_cfg)){};

    std::vector<cell> operator()(
        const traccc::bound_track_parameters& bound_param) const {

        std::vector<cell> ret;

        const auto vid = bound_param.surface_link().volume();
        const auto segmentation = (*m_digi_cfg)[vid];

        return ret;
    };

    std::unique_ptr<detector_t> m_detector;
    std::unique_ptr<digitization_map> m_digi_cfg;
};

}  // namespace traccc
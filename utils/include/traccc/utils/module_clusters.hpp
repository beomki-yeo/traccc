/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/track_parametrization.hpp"

namespace traccc {

struct module_value {
    std::vector<bound_indices> param_indices = {};
    std::vector<scalar> param_values = {};
    std::vector<scalar> param_variances = {};
};

struct module_clusters {

    module_clusters(Acts::BinUtility segmentation,
                    std::vector<bound_indices> geoIndices, bool merge,
                    double nsigma, bool commonCorner)
        : m_segmentation(std::move(segmentation)),
          m_geoIndices(std::move(geoIndices)),
          m_merge(merge),
          m_nsigma(nsigma),
          m_commonCorner(commonCorner) {}

    private:
    Acts::BinUtility m_segmentation;
    std::vector<bound_indices> m_geoIndices;
    std::vector<module_value> m_moduleValues;
    bool m_merge;
    double m_nsigma;
    bool m_commonCorner;
};

}  // namespace traccc
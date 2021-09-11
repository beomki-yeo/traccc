/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "edm/spacepoint.hpp"
#include "edm/internal_spacepoint.hpp"

// ActsCore
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

// detray include
#include "tests/common/test_defs.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"

namespace traccc {

// define spacepoint_grid
using spacepoint_grid = Acts::detail::Grid<
    int,
    Acts::detail::Axis<Acts::detail::AxisType::Equidistant,
                       Acts::detail::AxisBoundaryType::Closed>,
    Acts::detail::Axis<Acts::detail::AxisType::Equidistant,
                       Acts::detail::AxisBoundaryType::Bound> >;

    using phi_z_grid = detray::grid2<detray::attach_populator<false, internal_spacepoint<spacepoint>>, detray::axis::circular<>, detray::axis::regular<>, detray::serializer2>;
    
}  // namespace traccc

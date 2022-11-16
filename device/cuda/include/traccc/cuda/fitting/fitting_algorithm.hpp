/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/track_candidate.hpp"
#include "traccc/edm/track_state.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

// traccc library include(s).
#include "traccc/utils/memory_resource.hpp"

namespace traccc::cuda {

/// Fitting algorithm for a set of tracks
template <typename fitter_t>
class fitting_algorithm
    : public algorithm<track_state_container_types::buffer(
          const typename track_candidates_container_types::const_view&)> {

    public:

    /// Type delcarations
    using detector_type = typename fitter_t::detector_type;
    using transform3_type = typename fitter_t::transform3_type;

    /// Constructor with a detector
    fitting_algorithm(const detector_type& det);

    /// Operator executing the algorithm
    ///
    /// @param track_candidates_view is a view of candidate measurements from
    /// track finding
    /// @return the buffer of the fitted track parameters
    track_state_container_types::buffer operator()(
        const typename track_candidates_container_types::const_view&
            track_candidates_view) const override;

    private:
    /// Detector object
    std::unique_ptr<detector_type> m_detector;
    /// Memory resource used by the algorithm
    traccc::memory_resource m_mr;
    /// Copy object used by the algorithm
    std::unique_ptr<vecmem::copy> m_copy;
};

}  // namespace traccc::cuda
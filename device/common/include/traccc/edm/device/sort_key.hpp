/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/track_parameters.hpp"

namespace traccc::device {

//using sort_key = traccc::scalar;

struct sort_key{
    traccc::scalar key;
};


TRACCC_HOST_DEVICE
inline sort_key get_sort_key(const bound_track_parameters& params) {
    // key = |theta - pi/2|
    return sort_key{math::abs(params.theta() - constant<traccc::scalar>::pi_2)};
}


/// Comparator based on detray barcode value
struct sort_key_comp {
    TRACCC_HOST_DEVICE
    bool operator()(const sort_key& lhs, const sort_key& rhs) {
        return lhs.key < rhs.key;
    }
};

}  // namespace traccc::device

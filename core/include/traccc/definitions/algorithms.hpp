/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if defined(__CUDACC__)
#include <thrust/binary_search.h>
#else
#include <algorithm>
#endif

namespace traccc::detail {
/*
#if defined(__CUDACC__)
namespace algorithm = thrust;
#else
namespace algorithm = std;
#endif
*/
}  // namespace traccc::detail

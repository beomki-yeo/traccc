/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/cell.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

// System include(s).
#include <cassert>

namespace traccc {

/// Implemementation of SparseCCL, following
/// [DOI: 10.1109/DASIP48288.2019.9049184]
///
/// Requires cells to be sorted in column major
namespace detail {

/// Find root of the tree for entry @param e
///
/// @param L an equivalance table
///
/// @return the root of @param e
template <template <typename> class vector_t>
TRACCC_HOST_DEVICE inline unsigned int find_root(
    const vector_t<unsigned int>& L, unsigned int e) {

    unsigned int r = e;
    assert(r < L.size());
    while (L[r] != r) {
        r = L[r];
        assert(r < L.size());
    }
    return r;
}

/// Create a union of two entries @param e1 and @param e2
///
/// @param L an equivalance table
///
/// @return the rleast common ancestor of the entries
template <template <typename> class vector_t>
TRACCC_HOST_DEVICE inline unsigned int make_union(vector_t<unsigned int>& L,
                                                  unsigned int e1,
                                                  unsigned int e2) {

    int e;
    if (e1 < e2) {
        e = e1;
        assert(e2 < L.size());
        L[e2] = e;
    } else {
        e = e2;
        assert(e1 < L.size());
        L[e1] = e;
    }
    return e;
}

/// Helper method to find adjacent cells
///
/// @param a the first cell
/// @param b the second cell
///
/// @return boolan to indicate 8-cell connectivity
TRACCC_HOST_DEVICE inline bool is_adjacent(traccc::cell a, traccc::cell b) {
    return (a.channel0 - b.channel0) * (a.channel0 - b.channel0) <= 1 and
           (a.channel1 - b.channel1) * (a.channel1 - b.channel1) <= 1;
}

/// Helper method to find define distance,
/// does not need abs, as channels are sorted in
/// column major
///
/// @param a the first cell
/// @param b the second cell
///
/// @return boolan to indicate !8-cell connectivity
TRACCC_HOST_DEVICE inline bool is_far_enough(traccc::cell a, traccc::cell b) {
    return (a.channel1 - b.channel1) > 1;
}

/// Sparce CCL algorithm
///
/// @param cells is the cell collection
/// @param L is the vector of the output indices (to which cluster a cell
/// belongs to)
/// @param labels is the number of clusters found
template <template <typename> class vector_t>
TRACCC_HOST_DEVICE inline unsigned int sparse_ccl(
    const cell_collection_types::const_device& cells,
    vector_t<unsigned int>& L) {

    unsigned int labels = 0;

    // The number of cells.
    const unsigned int n_cells = cells.size();

    // first scan: pixel association
    unsigned int start_j = 0;
    for (unsigned int i = 0; i < n_cells; ++i) {
        L[i] = i;
        int ai = i;
        if (i > 0) {
            for (unsigned int j = start_j; j < i; ++j) {
                if (is_adjacent(cells[i], cells[j])) {
                    ai = make_union(L, ai, find_root(L, j));
                } else if (is_far_enough(cells[i], cells[j])) {
                    ++start_j;
                }
            }
        }
    }

    // second scan: transitive closure
    for (unsigned int i = 0; i < n_cells; ++i) {
        unsigned int l = 0;
        if (L[i] == i) {
            ++labels;
            l = labels;
        } else {
            l = L[L[i]];
        }
        L[i] = l;
    }

    return labels;
}
}  // namespace detail

}  // namespace traccc

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/clusterization/component_connection.hpp"

#include "traccc/clusterization/detail/sparse_ccl.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>

// System include(s)
#include <iterator>

namespace traccc {

component_connection::output_type component_connection::operator()(
    const cell_container_types::const_view& cells_view) const {

    // Create the result container.
    cluster_container_types::host result(&(m_mr.get()));
    auto& result_headers = result.get_headers();
    auto& result_items = result.get_items();

    unsigned int num_clusters = 0;

    const cell_container_types::const_device cells_device(cells_view);

    for (std::size_t i = 0; i < cells_device.size(); i++) {

        const auto module = cells_device.get_headers()[i];
        const cell_collection_types::const_device cells_per_module =
            cells_device.get_items()[i];

        // Set up the index vector.
        std::vector<unsigned int> CCL_indices(cells_per_module.size());

        // Run the algorithm
        auto num_clusters_per_module =
            detail::sparse_ccl(cells_per_module, CCL_indices);
        num_clusters += num_clusters_per_module;
        result.reserve(num_clusters);

        cluster_container_types::host clusters_per_module(
            num_clusters_per_module, &(m_mr.get()));

        auto& cluster_modules = clusters_per_module.get_headers();
        for (auto& cl_id : cluster_modules) {
            cl_id.module = module.module;
            cl_id.placement = module.placement;
        }

        auto& cluster_cells = clusters_per_module.get_items();
        unsigned int icell = 0;
        for (auto cell_label : CCL_indices) {
            auto cindex = static_cast<unsigned int>(cell_label - 1);
            if (cindex < cluster_cells.size()) {
                cluster_cells[cindex].push_back(cells_per_module[icell++]);
            }
        }

        result_headers.insert(result_headers.end(),
                              std::make_move_iterator(cluster_modules.begin()),
                              std::make_move_iterator(cluster_modules.end()));

        result_items.insert(result_items.end(),
                            std::make_move_iterator(cluster_cells.begin()),
                            std::make_move_iterator(cluster_cells.end()));
    }

    return get_data(result);

    /*
    // Create the result container.
    output_type result(&(m_mr.get()));

    // Set up a device collection on top of the host collection.
    const cell_collection_types::const_view cells_view =
        vecmem::get_data(cells);
    const cell_collection_types::const_device cells_device(cells_view);

    // Set up the index vector.
    vecmem::vector<unsigned int> connected_cells(cells.size(), &(m_mr.get()));
    vecmem::device_vector<unsigned int> connected_cells_device(
        vecmem::get_data(connected_cells));

    // Run the algorithm
    unsigned int num_clusters = 0;
    detail::sparse_ccl(cells_device, connected_cells_device, num_clusters);

    result.resize(num_clusters);
    for (auto& cl_id : result.get_headers()) {
        cl_id.module = module.module;
        cl_id.placement = module.placement;
    }

    auto& cluster_items = result.get_items();
    unsigned int icell = 0;
    for (auto cell_label : connected_cells) {
        auto cindex = static_cast<unsigned int>(cell_label - 1);
        if (cindex < cluster_items.size()) {
            cluster_items[cindex].push_back(cells[icell++]);
        }
    }

    // Return the cluster container.
    return result;
    */
}

}  // namespace traccc

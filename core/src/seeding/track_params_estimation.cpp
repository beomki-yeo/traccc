/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/seeding/track_params_estimation.hpp"

#include "traccc/seeding/track_params_estimation_helper.hpp"

namespace traccc {

track_params_estimation::track_params_estimation(vecmem::memory_resource& mr)
    : m_mr(mr) {}

track_params_estimation::output_type track_params_estimation::operator()(
    const spacepoint_container_types::const_view& spacepoints_view,
    const seed_collection_types::const_view& seeds_view) const {

    bound_track_parameters_collection_types::host result(&m_mr.get());
    const spacepoint_container_types::const_device spacepoints(
        spacepoints_view);
    const traccc::seed_collection_types::const_device seeds(seeds_view);

    // convenient assumption on bfield and mass
    // TODO: Make use of bfield extenstion in the future
    vector3 bfield = {0, 0, 2};

    for (const auto& seed : seeds) {
        bound_track_parameters track_params;
        track_params.vector() =
            seed_to_bound_vector(spacepoints, seed, bfield, PION_MASS_MEV);

        result.push_back(track_params);
    }

    return vecmem::get_data(result);
}

}  // namespace traccc

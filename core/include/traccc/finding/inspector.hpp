/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/tuple_helpers.hpp"

// System include(s)
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

namespace traccc {

namespace propagation {

/// Print inspector that runs as actor in the propagation
struct print_inspector : detray::actor {

    struct state {
        std::stringstream stream{};

        std::string to_string() const { return stream.str(); }
    };

    template <typename propagation_state_t>
    DETRAY_HOST_DEVICE void operator()(
        state &printer, const propagation_state_t &prop_state) const {
        const auto &navigation = prop_state._navigation;
        const auto &stepping = prop_state._stepping;

        printer.stream << std::left << std::setw(30);
        switch (navigation.status()) {
            case detray::navigation::status::e_abort:
                printer.stream << "status: abort";
                break;
            case detray::navigation::status::e_on_target:
                printer.stream << "status: e_on_target";
                break;
            case detray::navigation::status::e_unknown:
                printer.stream << "status: unknowm";
                break;
            case detray::navigation::status::e_towards_object:
                printer.stream << "status: towards_surface";
                break;
            case detray::navigation::status::e_on_module:
                printer.stream << "status: on_module";
                break;
            case detray::navigation::status::e_on_portal:
                printer.stream << "status: on_portal";
                break;
        };

        if (navigation.volume() == detray::dindex_invalid) {
            printer.stream << "volume: " << std::setw(10) << "invalid";
        } else {
            printer.stream << "volume: " << std::setw(10)
                           << navigation.volume();
        }

        if (navigation.current_object().is_invalid()) {
            printer.stream << "surface: " << std::setw(14) << "invalid";
        } else {
            printer.stream << "surface: " << std::setw(14)
                           << navigation.current_object();
        }

        printer.stream << "step_size: " << std::setw(10) << stepping._step_size
                       << std::endl;

        printer.stream
            << std::setw(10)
            << detray::detail::ray<
                   typename propagation_state_t::detector_type::transform3>(
                   stepping())
            << std::endl;
    }
};

}  // namespace propagation

}  // namespace traccc
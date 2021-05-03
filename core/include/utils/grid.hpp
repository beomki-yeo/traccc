/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "utils/axis.hpp"

template < clase... Axes >
class grid final{
public:
    /// number of dimensions of the grid
    static constexpr size_t DIM = sizeof...(Axes);

    /// type for points in d-dimensional grid space
    using point_t = std::array<scalar, DIM>;
    /// index type using local bin indices along each axis
    using index_t = std::array<size_t, DIM>;
    
    /// @brief default constructor
    ///
    /// @param [in] axes actual axis objects spanning the grid
    grid(std::tuple<Axes...> axes) : m_axes(std::move(axes)) {
	m_values.resize(size());
    }

};

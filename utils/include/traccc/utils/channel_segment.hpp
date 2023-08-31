/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc {

/// Nested struct for representing channel steps.
struct channel_segment {
    /// Shorthand for a 2D segment
    using Segment2D = std::array<point2, 2>;
    /// Shorthand for a 2D bin
    using Bin2D = std::array<unsigned int, 2>;

    /// The bin of this segment
    Bin2D bin = {0, 0};
    /// The segment start, end points
    Segment2D path2D;
    /// The (clipped) value (uncorrected: path length)
    double activation = 0.;

    /// Constructor with arguments
    ///
    /// @param bin_ The bin corresponding to this step
    /// @param path2D_ The start/end 2D position of the segement
    /// @param activation_ The segment activation (clean: length) for this bin
    channel_segment(Bin2D bin_, Segment2D path2D_, double activation_)
        : bin(bin_), path2D(std::move(path2D_)), activation(activation_) {}
};

}  // namespace traccc
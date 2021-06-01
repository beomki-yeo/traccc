/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc{
    
    struct sp_location{
	/// index of the bin of the spacepoint grid
	size_t bin_idx;
	/// index of the spacepoint in the bin
	size_t sp_idx; 
    };       	
    
    struct lin_circle {
	float Zo;
	float cotTheta;
	float iDeltaR;
	float Er;
	float U;
	float V;
    };

    // Middle - Bottom or Middle - Top
    struct doublet{	
	sp_location sp1;
	sp_location sp2;
	lin_circle lin;
    };

    /// Container of doublet belonging to one detector module
    template< template< typename > class vector_t >
    using doublet_collection = vector_t< doublet >;

    /// Convenience declaration for the doublet collection type to use in host code
    using host_doublet_collection
    = doublet_collection< vecmem::vector >;

    /// Convenience declaration for the doublet collection type to use in device code
    using device_doublet_collection
    = doublet_collection< vecmem::device_vector >;

    /// Convenience declaration for the doublet container type to use in host code
    using host_doublet_container
    = host_container< int, doublet >;

    /// Convenience declaration for the doublet container type to use in device code
    using device_doublet_container
    = device_container< int, doublet>;

    /// Convenience declaration for the doublet container data type to use in host code
    using doublet_container_data
    = container_data< int, doublet >;

    /// Convenience declaration for the doublet container buffer type to use in host code
    using doublet_container_buffer
    = container_buffer< int, doublet >;

    /// Convenience declaration for the doublet container view type to use in host code
    using doublet_container_view
    = container_view< int, doublet >;    
    
} // namespace traccc

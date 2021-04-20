/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <vector>

namespace traccc {
    
    /// A cell definition: maximum two channel identifiers
    /// and one activiation value;
    struct spacepoint {
        point3 global = { 0., 0., 0.};
        variance3 variance = { 0., 0., 0.};
    };

    /// Container of spacepoints belonging to one detector module
    template< template< typename > class vector_t >
    using spacepoint_collection = vector_t< spacepoint >;

    /// Convenience declaration for the spacepoint collection type to use in host code
    using host_spacepoint_collection = spacepoint_collection< vecmem::vector >;
    /// Convenience declaration for the spacepoint collection type to use in device code
    using device_spacepoint_collection = spacepoint_collection< vecmem::device_vector >;

    /// Container describing all of the spacepoints in a given event
    ///
    /// This is the "main" spacepoint container of the code, holding all relevant
    /// information about all of the spacepoints in a given event.
    ///
    /// It can be instantiated with different vector types, to be able to use
    /// the same container type in both host and device code.
    ///
    template< template< typename > class vector_t,
              template< typename > class jagged_vector_t >
    class spacepoint_container {

    public:
        /// @name Type definitions
        /// @{

        /// Vector type used by the spacepoint container
        template< typename T >
        using vector_type = vector_t< T >;
        /// Jagged vector type used by the spacepoint container
        template< typename T >
        using jagged_vector_type = jagged_vector_t< T >;

        /// The spacepoint module vector type
        using module_vector = vector_type< geometry_id >;
        /// The spacepoint vector type
        using spacepoint_vector = jagged_vector_type< spacepoint >;

        /// @}

        /// Headers for all of the modules (holding spacepoints) in the event
        module_vector modules;
        /// All of the spacepoints in the event
        spacepoint_vector spacepoints;

	/// Reserve
	void reserve(size_t val){
	    modules.reserve(val);
	    spacepoints.reserve(val);
	}
    }; // class spacepoint_container

    /// Convenience declaration for the spacepoint container type to use in host code
    using host_spacepoint_container =
        spacepoint_container< vecmem::vector, vecmem::jagged_vector >;
    /// Convenience declaration for the spacepoint container type to use in device code
    using device_spacepoint_container =
        spacepoint_container< vecmem::device_vector, vecmem::jagged_device_vector >;

    /// @}

    /// @name Types used to send data back and forth between host and device code
    /// @{

    /// Structure holding (some of the) data about the spacepoints in host code
    struct spacepoint_container_data {
        vecmem::data::vector_view< geometry_id > modules;
        vecmem::data::jagged_vector_data< spacepoint > spacepoints;
    }; // struct spacepoint_container_data

    /// Structure holding (all of the) data about the spacepoints in host code
    struct spacepoint_container_buffer {
        vecmem::data::vector_buffer< geometry_id > modules;
        vecmem::data::jagged_vector_buffer< spacepoint > spacepoints;
    }; // struct spacepoint_container_data

    /// Structure used to send the data about the spacepoints to device code
    ///
    /// This is the type that can be passed to device code as-is. But since in
    /// host code one needs to manage the data describing a
    /// @c traccc::spacepoint_container either using @c traccc::spacepoint_container_data or
    /// @c traccc::spacepoint_container_buffer, it needs to have constructors from
    /// both of those types.
    ///
    /// In fact it needs to be created from one of those types, as such an
    /// object can only function if an instance of one of those types exists
    /// alongside it as well.
    ///
    struct spacepoint_container_view {

        /// Constructor from a @c spacepoint_container_data object
        spacepoint_container_view( const spacepoint_container_data& data )
        : modules( data.modules ), spacepoints( data.spacepoints ) {}

        /// Constructor from a @c spacepoint_container_buffer object
        spacepoint_container_view( const spacepoint_container_buffer& buffer )
        : modules( buffer.modules ), spacepoints( buffer.spacepoints ) {}

        /// View of the data describing the headers of the spacepoint holding modules
        vecmem::data::vector_view< geometry_id > modules;
        /// View of the data describing all of the spacepoints
        vecmem::data::jagged_vector_view< spacepoint > spacepoints;

    }; // struct spacepoint_container_view

    /// Helper function for making a "simple" object out of the spacepoint container
    spacepoint_container_data get_data( host_spacepoint_container& cc, vecmem::memory_resource* resource = nullptr ) {
        return { { vecmem::get_data( cc.modules ) },
                 { vecmem::get_data( cc.spacepoints, resource ) } };
    }

    /// @}
}

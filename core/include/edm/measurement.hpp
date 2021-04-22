/** TRACCC library, part of the ACTS projWect (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "definitions/primitives.hpp"
#include "definitions/algebra.hpp"
#include <vector>

namespace traccc {
    
    /// A measurement definition, fix to two-dimensional here
    struct measurement {
        point2 local = {0., 0.};
        variance2 variance = { 0., 0.};
    };

    /// Container of measurements belonging to one detector module
    template< template< typename > class vector_t >
    using measurement_collection = vector_t< measurement >;

    /// Convenience declaration for the measurement collection type to use in host code
    using host_measurement_collection = measurement_collection< vecmem::vector >;
    /// Convenience declaration for the measurement collection type to use in device code
    using device_measurement_collection = measurement_collection< vecmem::device_vector >;

    /// Container describing all of the measurements in a given event
    ///
    /// This is the "main" measurement container of the code, holding all relevant
    /// information about all of the measurements in a given event.
    ///
    /// It can be instantiated with different vector types, to be able to use
    /// the same container type in both host and device code.
    ///
    template< template< typename > class vector_t,
              template< typename > class jagged_vector_t >
    class measurement_container {

    public:
        /// @name Type definitions
        /// @{

        /// Vector type used by the measurement container
        template< typename T >
        using vector_type = vector_t< T >;
        /// Jagged vector type used by the measurement container
        template< typename T >
        using jagged_vector_type = jagged_vector_t< T >;

        /// The measurement module vector type
        using cell_module_vector = vector_type< cell_module >;
        /// The measurement vector type
        using measurement_vector = jagged_vector_type< measurement >;

        /// @}

        /// Headers for all of the modules (holding measurements) in the event
        cell_module_vector modules;
        /// All of the measurements in the event
        measurement_vector measurements;

    }; // class measurement_container

    /// Convenience declaration for the measurement container type to use in host code
    using host_measurement_container =
        measurement_container< vecmem::vector, vecmem::jagged_vector >;
    /// Convenience declaration for the measurement container type to use in device code
    using device_measurement_container =
        measurement_container< vecmem::device_vector, vecmem::jagged_device_vector >;

    /// @}

    /// @name Types used to send data back and forth between host and device code
    /// @{

    /// Structure holding (some of the) data about the measurements in host code
    struct measurement_container_data {
        vecmem::data::vector_view< cell_module > modules;
        vecmem::data::jagged_vector_data< measurement > measurements;
    }; // struct measurement_container_data

    /// Structure holding (all of the) data about the measurements in host code
    struct measurement_container_buffer {
        vecmem::data::vector_buffer< cell_module > modules;
        vecmem::data::jagged_vector_buffer< measurement > measurements;
    }; // struct measurement_container_data

    /// Structure used to send the data about the measurements to device code
    ///
    /// This is the type that can be passed to device code as-is. But since in
    /// host code one needs to manage the data describing a
    /// @c traccc::measurement_container either using @c traccc::measurement_container_data or
    /// @c traccc::measurement_container_buffer, it needs to have constructors from
    /// both of those types.
    ///
    /// In fact it needs to be created from one of those types, as such an
    /// object can only function if an instance of one of those types exists
    /// alongside it as well.
    ///
    struct measurement_container_view {

        /// Constructor from a @c measurement_container_data object
        measurement_container_view( const measurement_container_data& data )
        : modules( data.modules ), measurements( data.measurements ) {}

        /// Constructor from a @c measurement_container_buffer object
        measurement_container_view( const measurement_container_buffer& buffer )
        : modules( buffer.modules ), measurements( buffer.measurements ) {}

        /// View of the data describing the headers of the measurement holding modules
        vecmem::data::vector_view< cell_module > modules;
        /// View of the data describing all of the measurements
        vecmem::data::jagged_vector_view< measurement > measurements;

    }; // struct measurement_container_view

    /// Helper function for making a "simple" object out of the measurement container
    __inline__
    measurement_container_data get_data( host_measurement_container& cc, vecmem::memory_resource* resource = nullptr ) {
        return { { vecmem::get_data( cc.modules ) },
                 { vecmem::get_data( cc.measurements, resource ) } };
    }
    
    /// @}

    
}

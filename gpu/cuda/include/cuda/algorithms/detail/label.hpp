/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "vecmem/memory/cuda/managed_memory_resource.hpp"
#include "edm/cell.hpp"

namespace traccc {
namespace cuda{
namespace detail{

// Label container that stores the sparse ccl information
    template< template< typename > class vector_t,
              template< typename > class jagged_vector_t >
    class label_container {

    public:
        /// @name Type definitions
        /// @{

        /// Vector type used by the label container
        template< typename T >
        using vector_type = vector_t< T >;
        /// Jagged vector type used by the label container
        template< typename T >
        using jagged_vector_type = jagged_vector_t< T >;

        /// The label module vector type
        using ccl_label = vector_type< unsigned int >;
        /// The label vector type
        using ccl_label_vector = jagged_vector_type< unsigned int >;
	
        /// @}

        /// number of labels per module
	ccl_label counts;
        /// All of the labels in the event
        ccl_label_vector labels;

    };

    using host_label_container =
        label_container< vecmem::vector, vecmem::jagged_vector >;
    
    using device_label_container =
        label_container< vecmem::device_vector, vecmem::jagged_device_vector >;

    /// Structure holding (some of the) data about the labels in host code
    struct label_container_data {
        vecmem::data::vector_view< unsigned int > counts;
        vecmem::data::jagged_vector_data< unsigned int > labels;
    }; 

    struct label_container_view {

        /// Constructor from a @c label_container_data object
        label_container_view( const label_container_data& data )
        : counts( data.counts ), labels( data.labels ) {}

	/// View of the data describing the headers of the label holding modules
        vecmem::data::vector_view< unsigned int > counts;
        /// View of the data describing all of the labels
        vecmem::data::jagged_vector_view< unsigned int > labels;

    };
    
    /// Helper function for making a "simple" object out of the label container
    label_container_data get_data( host_label_container& lc, vecmem::memory_resource* resource = nullptr  ) {
        return { { vecmem::get_data( lc.counts ) },
                 { vecmem::get_data( lc.labels, resource ) } };
    }

}
}
}

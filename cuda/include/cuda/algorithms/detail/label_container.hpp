/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "vecmem/memory/cuda/managed_memory_resource.hpp"

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
	ccl_label num_label;
        /// All of the labels in the event
        ccl_label_vector labels;
	
    };

    using host_label_container =
        label_container< vecmem::vector, vecmem::jagged_vector >;

    using managed_label_container =
        label_container< vecmem::vector, vecmem::jagged_vector >;
    
    using device_label_container =
        label_container< vecmem::device_vector, vecmem::jagged_device_vector >;

    /// Structure holding (some of the) data about the labels in host code
    struct label_container_data {
        vecmem::data::vector_view< unsigned int > num_label;
        vecmem::data::jagged_vector_data< unsigned int > labels;
    }; 
    
    /// Helper function for making a "simple" object out of the label container
    label_container_data get_data( host_label_container& lc, vecmem::memory_resource* resource = nullptr  ) {
        return { { vecmem::get_data( lc.num_label ) },
                 { vecmem::get_data( lc.labels, resource ) } };
    }
    /*
    label_container_data get_data( managed_label_container& lc, vecmem::cuda::managed_memory_resource& resource ) {
        return { { vecmem::get_data( lc.num_label ) },
                 { vecmem::get_data( lc.labels, resource ) } };
    } 
    */       
}
}
}

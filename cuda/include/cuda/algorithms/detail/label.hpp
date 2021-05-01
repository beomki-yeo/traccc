/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "vecmem/memory/cuda/managed_memory_resource.hpp"

namespace traccc {
namespace cuda{
namespace detail{

    // Define label container
    using label = unsigned int;
    using host_label_container = host_container< label, label >;   
    using device_label_container = device_container< label, label >;
    using label_container_data = container_data< label, label >;
    using label_container_buffer = container_buffer< label, label >;
    using label_container_view = container_view< label, label >;
    
}
}
}


# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Guard against multiple includes.
include_guard( GLOBAL )

# Include the helper function(s).
include( traccc-functions )

set( CMAKE_CUDA_STANDARD 17 CACHE STRING "CUDA C++ standard to use" )

foreach( mode RELEASE RELWITHDEBINFO MINSIZEREL DEBUG )
  traccc_add_flag( CMAKE_CUDA_FLAGS_${mode} "--expt-relaxed-constexpr" )
  traccc_add_flag( CMAKE_CUDA_FLAGS_${mode} "-fmad=false" )
endforeach()
  

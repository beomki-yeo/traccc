# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Guard against multiple includes.
include_guard( GLOBAL )

# Tell the user what's happening.
message( STATUS "Building detray as part of the traccc project" )

# Declare where to get VecMem from.
FetchContent_Declare( Detray
  GIT_REPOSITORY "https://github.com/acts-project/detray.git"
  GIT_TAG        552d5622a0dc4aeb99a13576b4ebc12ca469eadd)

# Prevent Detray from building its tests and benchmarks
# builds/uses GoogleTest.
set( DETRAY_UNIT_TESTS OFF )
set( DETRAY_BENCHMARKS OFF )

# Get it into the current directory.
FetchContent_Populate( Detray )
add_subdirectory( "${detray_SOURCE_DIR}" "${detray_BINARY_DIR}"
  EXCLUDE_FROM_ALL )

# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the "build" of the traccc::option library.
traccc_add_library( traccc_utils utils TYPE SHARED
  # header files
  "include/traccc/utils/digitization_algorithm.hpp"
  "include/traccc/utils/digitization_writer.hpp"
  )
target_link_libraries( traccc_utils 
  PUBLIC traccc::io detray::core detray::utils )
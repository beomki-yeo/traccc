# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Silence FetchContent warnings with CMake >=3.24.
if( POLICY CMP0135 )
   cmake_policy( SET CMP0135 NEW )
endif()

# Tell the user what's happening.
message( STATUS "Building oneDPL as part of the TRACCC project" )

# Declare where to get DPL from.
set( TRACCC_DPL_SOURCE
   "URL;https://github.com/oneapi-src/oneDPL/archive/refs/tags/oneDPL-2022.0.0-release.tar.gz;URL_MD5;46b60bd87325042df716d1df4686eb3b"
   CACHE STRING "Source for DPL, when built as part of this project" )
mark_as_advanced( TRACCC_DPL_SOURCE )
FetchContent_Declare( DPL ${TRACCC_DPL_SOURCE} )

# Exclude TBB
set( ONEDPL_BACKEND dpcpp_only )

# Get it into the current directory.
FetchContent_MakeAvailable( DPL )

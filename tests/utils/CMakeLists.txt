# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Declare the io library test(s).
traccc_add_test( utils
   "test_digitization.cpp"
   LINK_LIBRARIES GTest::gtest_main traccc_tests_common
                  traccc::utils detray::utils )

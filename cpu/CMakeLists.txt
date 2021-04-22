add_library(traccc_cpu INTERFACE)

target_include_directories(traccc_cpu
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(traccc_cpu
  INTERFACE traccc::core)

install(TARGETS traccc_cpu
        PUBLIC_HEADER
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cpu
)

install(
  DIRECTORY include/algorithms
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(traccc::cpu ALIAS traccc_cpu)

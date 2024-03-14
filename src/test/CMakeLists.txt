enable_testing()

add_executable(
    google_test
    google_test.cpp
)
target_link_libraries(
    google_test
    GTest::gtest_main
)

set(vmc_src
  "../"
  )
include_directories(${vmc_src})

include(GoogleTest)
gtest_discover_tests(google_test)

target_link_libraries(google_test VMC_TYPE LATTICE)
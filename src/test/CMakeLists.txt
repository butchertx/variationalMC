

include(Dart)
set(CTEST_OUTPUT_ON_FAILURE 1)

# To allow type definitions for VSCode debugging
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")

# Set data path env variables
add_definitions(-DEXAMPLES_PATH="${PROJECT_SOURCE_DIR}/examples/")

# UNIT TESTS

# Set executable
set(UNIT_TEST_SOURCE UnitTest.cpp)
set(UNIT_TEST_TARGET unit_test)

set(EXAMPLE_TEST_SOURCE ExampleTest.cpp)
set(EXAMPLE_TEST_TARGET example_test)

set(zheev_TEST_SOURCE test_zheev.cpp)
set(zheev_TEST_TARGET test_zheev)

# Set compile and link flags
set(COMPILE_FLAGS "-Wall --std=c++17 -DMKL_ILP64")
set(LINK_FLAGS "-v")

set(vmc_src
  "../"
  )
include_directories(${vmc_src})

# Set executables

add_executable(${UNIT_TEST_TARGET} ${UNIT_TEST_SOURCE})
set_target_properties(${UNIT_TEST_TARGET} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
set_target_properties(${UNIT_TEST_TARGET} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
target_link_libraries(
    ${UNIT_TEST_TARGET}
    GTest::gtest_main
    VMC_IO
    VMC_TYPE
    MATRIX
    LATTICE
    MEAN_FIELD_ANSATZ
    # PROJECTED_STATE
    # WAVEFUNCTION
    $<LINK_ONLY:MKL::MKL>
)

add_executable(${EXAMPLE_TEST_TARGET} ${EXAMPLE_TEST_SOURCE})
set_target_properties(${EXAMPLE_TEST_TARGET} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
set_target_properties(${EXAMPLE_TEST_TARGET} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
target_link_libraries(
    ${EXAMPLE_TEST_TARGET}
    GTest::gtest_main
    VMC_DRIVER
    VMC_IO
    VMC_TYPE
    LATTICE
    MEAN_FIELD_ANSATZ
    PROJECTED_STATE
    WAVEFUNCTION
    SPIN_MODEL
    VMC_RUNNER
    VMC_NUMERICAL
    VMC_RESULTS
    $<LINK_ONLY:MKL::MKL>
)

add_executable(${zheev_TEST_TARGET} ${zheev_TEST_SOURCE})
set_target_properties(${zheev_TEST_TARGET} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
set_target_properties(${zheev_TEST_TARGET} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
target_link_libraries(
    ${zheev_TEST_TARGET}
    MATRIX
    $<LINK_ONLY:MKL::MKL>
)

# Discover tests

include(GoogleTest)
gtest_discover_tests(unit_test)
gtest_discover_tests(example_test)
enable_testing()

# To allow type definitions for VSCode debugging
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")

# Set data path env variables
add_definitions(-DEXAMPLES_PATH="${PROJECT_SOURCE_DIR}/examples/")

# UNIT TESTS

# Set executable
set(SOURCE_FILE unit_test.cpp)
set(TARGET_NAME unit_test)

# Set compile and link flags
set(COMPILE_FLAGS "-Wall --std=c++17" )
set(LINK_FLAGS "-v")

set(vmc_src
  "../"
  )
include_directories(${vmc_src})

add_executable(${TARGET_NAME} ${SOURCE_FILE})
set_target_properties(${TARGET_NAME} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
set_target_properties(${TARGET_NAME} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
target_link_libraries(
    unit_test
    GTest::gtest_main
    VMC_IO
    VMC_TYPE
    LATTICE
    MEAN_FIELD_ANSATZ
    PROJECTED_STATE
    WAVEFUNCTION
    $<LINK_ONLY:MKL::MKL>
)


# BENCHMARKS


# # Set executable
# set(SOURCE_FILE benchmark.cpp)
# set(TARGET_NAME benchmark)

# # Set compile and link flags
# set(COMPILE_FLAGS "-Wall")
# set(LINK_FLAGS "-v")

# add_executable(${TARGET_NAME} ${SOURCE_FILE})
# set_target_properties(${TARGET_NAME} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
# set_target_properties(${TARGET_NAME} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# target_link_libraries(${TARGET_NAME} GTest::gtest_main VMC_TYPE LATTICE)

# target_link_libraries(
#     benchmark
#     GTest::gtest_main
#     VMC_TYPE
#     LATTICE
# )

include(GoogleTest)
gtest_discover_tests(unit_test)
# gtest_discover_tests(benchmark)
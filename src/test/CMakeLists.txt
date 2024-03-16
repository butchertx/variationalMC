enable_testing()

# To allow type definitions for VSCode debugging
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")

# Set executable
set(SOURCE_FILE google_test.cpp)
set(TARGET_NAME google_test)

# Set compile and link flags
set(COMPILE_FLAGS "-Wall")
set(LINK_FLAGS "-v")

set(vmc_src
  "../"
  )
include_directories(${vmc_src})

# To compile in a single command:
#    icpx <file>.cpp -o <file>
# CMake executes:
#    [compile] icpx -o <file>.cpp.o -c <file>.cpp
#    [link]    icpx <file>.cpp.o -o <file>
add_executable(${TARGET_NAME} ${SOURCE_FILE})
set_target_properties(${TARGET_NAME} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
set_target_properties(${TARGET_NAME} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
target_link_libraries(${TARGET_NAME} GTest::gtest_main VMC_TYPE LATTICE)

target_link_libraries(
    google_test
    GTest::gtest_main
    VMC_TYPE
    LATTICE
)

include(GoogleTest)
gtest_discover_tests(google_test)
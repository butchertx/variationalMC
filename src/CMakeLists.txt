# To allow type definitions for VSCode debugging
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")

# Set executable
set(SOURCE_FILE example.cpp)
set(TARGET_NAME run_example)

# Set libraries
add_library(VMC_IO
    vmc_io.cpp
)
add_library(VMC_TYPE
    vmctype.cpp
)
add_library(LATTICE
    Lattice.cpp
)
add_library(MEAN_FIELD_ANSATZ
    MeanFieldAnsatz.cpp
)

# Set compile and link flags
set(COMPILE_FLAGS "-Wall")
set(LINK_FLAGS "-v")

# To compile in a single command:
#    icpx <file>.cpp -o <file>
# CMake executes:
#    [compile] icpx -o <file>.cpp.o -c <file>.cpp
#    [link]    icpx <file>.cpp.o -o <file>
add_executable(${TARGET_NAME} ${SOURCE_FILE})
set_target_properties(${TARGET_NAME} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
set_target_properties(${TARGET_NAME} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
add_custom_target(cpu-gpu DEPENDS ${TARGET_NAME})


target_link_libraries(${TARGET_NAME} VMC_IO VMC_TYPE LATTICE MEAN_FIELD_ANSATZ)
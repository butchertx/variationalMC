

# Set executable
set(SOURCE_FILE latticex.cpp)
set(TARGET_NAME latticex)

# Set compile and link flags
set(COMPILE_FLAGS "-Wall")
# set(LINK_FLAGS "-fsycl")

# To compile in a single command:
#    icpx <file>.cpp -o <file>
# CMake executes:
#    [compile] icpx -o <file>.cpp.o -c <file>.cpp
#    [link]    icpx <file>.cpp.o -o <file>
add_executable(${TARGET_NAME} ${SOURCE_FILE})
set_target_properties(${TARGET_NAME} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
set_target_properties(${TARGET_NAME} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
add_custom_target(cpu-gpu DEPENDS ${TARGET_NAME})
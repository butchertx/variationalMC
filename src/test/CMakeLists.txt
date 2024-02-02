enable_testing()

add_executable(
    google_test
    google_test.cc
)
target_link_libraries(
    google_test
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(google_test)
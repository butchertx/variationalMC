#include <gtest/gtest.h>
#include <vmctype.h>

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions2) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello1", "world2");
  // Expect equality.
  EXPECT_EQ(7 * 8, 56);
}

TEST(VMCHello, VMCHello1) {
    vmctype::lattice_options lat;
}
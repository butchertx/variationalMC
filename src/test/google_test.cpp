#include <gtest/gtest.h>
#include <vmctype.h>
#include <Lattice.h>

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

TEST(LatticeTest, ValidateLattices) {
  Lattice chainLattice(Lattice_type_t::CHAIN, vec3<int>(2, 1, 1), vec3<int>(0, 0, 0));
  auto n = chainLattice.get_N();
  EXPECT_EQ(n, 2);
}
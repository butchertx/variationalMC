#include <gtest/gtest.h>
#include <filesystem>
#include <vmctype.h>
#include <vmc_io.h>
#include <Lattice.h>
#include <MeanFieldAnsatz.h>
#include <ProjectedState.h>
#include "mkl_types.h"

// Global Variables

#ifdef EXAMPLES_PATH
    const std::string EXAMPLES_DIR = EXAMPLES_PATH;
#else
    const std::string EXAMPLES_DIR = std::filesystem::absolute("../../examples/");
#endif

// Hello World Tests

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

// Lattice Tests

TEST(LatticeTest, CreateBasicLattices) {
    Lattice chainLattice(Lattice_type_t::CHAIN, vec3<int>(2, 1, 1), vec3<int>(0, 0, 0));
    auto n = chainLattice.get_N();
    EXPECT_EQ(n, 2);
}

// MFAnsatz Tests

class MFAnsatzTest : public ::testing::Test {

protected:
    void SetUp() override {
        chainLatticeHalf = Lattice(Lattice_type_t::CHAIN, vec3<int>(2, 1, 1), vec3<int>(0, 0, 0));
        chainLatticeOne = Lattice(Lattice_type_t::CHAIN, vec3<int>(3, 1, 1), vec3<int>(1, 0, 0));
    }

    Lattice chainLatticeHalf, chainLatticeOne;
};

TEST_F(MFAnsatzTest, CheckFixture) {
    EXPECT_EQ(chainLatticeHalf.get_N(), 2);
    EXPECT_EQ(chainLatticeOne.get_N(), 3);
}

// TEST_F(MFAnsatzTest, SpinOneChain) {
//     mf_ansatz = MeanFieldAnsatz_ONE()
// }

// ProjectedState Tests

class ProjectedStateTest : public ::testing::Test {

protected:
    void SetUp() override {
        chainLatticeHalf = Lattice(Lattice_type_t::CHAIN, vec3<int>(2, 1, 1), vec3<int>(0, 0, 0));
        chainLatticeOne = Lattice(Lattice_type_t::CHAIN, vec3<int>(3, 1, 1), vec3<int>(1, 0, 0));
        
        chain_half_options = read_json_wavefunction_from_dir(EXAMPLES_DIR + "spin_half/1d/trivial");
        chain_one_options = read_json_wavefunction_from_dir(EXAMPLES_DIR + "spin_one/1d/trivial");

        mf_ansatz_half = std::shared_ptr<MeanFieldAnsatz>(new MeanFieldAnsatz_HALF(chain_half_options, chainLatticeHalf));
        mf_ansatz_one = std::shared_ptr<MeanFieldAnsatz>(new MeanFieldAnsatz_ONE(chain_one_options, chainLatticeOne));
        
    }

    Lattice chainLatticeHalf, chainLatticeOne;

    WavefunctionOptions chain_half_options, chain_one_options;

    std::shared_ptr<MeanFieldAnsatz> mf_ansatz_half, mf_ansatz_one;
};

TEST_F(ProjectedStateTest, CheckFixture) {
    EXPECT_EQ(chainLatticeHalf.get_N(), 2);
    EXPECT_EQ(chainLatticeOne.get_N(), 3);

    RandomEngine r(0, chainLatticeHalf.get_N(), chainLatticeHalf.get_neighbor_counts()[0]);
    ProjectedState wf_half(*mf_ansatz_half, r);
    std::vector<int> flips({0, 1});
    wf_half.update(flips);
    wf_half.print_matrix("Slater");
    wf_half.print_matrix("LU");
    wf_half.print_matrix("Winv");
    wf_half.print_matrix("UP1");
    wf_half.print_matrix("UP2");
    wf_half.print_matrix("UP3");
    wf_half.print_matrix("ipiv");

    RandomEngine r2(0, chainLatticeOne.get_N(), chainLatticeOne.get_neighbor_counts()[0]);
    ProjectedState wf_one(*mf_ansatz_one, r2);
    wf_one.update(flips);
    wf_one.print_matrix("Slater");
    wf_one.print_matrix("LU");
    wf_one.print_matrix("Winv");
    wf_one.print_matrix("UP1");
    wf_one.print_matrix("UP2");
    wf_one.print_matrix("UP3");
    wf_one.print_matrix("ipiv");
    EXPECT_TRUE(true);
}
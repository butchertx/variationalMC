#include <gtest/gtest.h>
#include <filesystem>
#include <vmctype.h>
#include <vmc_io.h>
#include <Matrix.h>
#include <Lattice.h>
#include <MeanFieldAnsatz.h>
// #include <ProjectedState.h>
#include "mkl_types.h"

// Global Variables

#ifdef EXAMPLES_PATH
    const std::string EXAMPLES_DIR = EXAMPLES_PATH;
#else
    const std::string EXAMPLES_DIR = std::filesystem::absolute("../../examples/");
#endif

// Custom Assertions

inline void AssertMKLComplexEqual__(const MKL_Complex16& base, const MKL_Complex16& other)
{
    ASSERT_DOUBLE_EQ(base.real, other.real);
    ASSERT_DOUBLE_EQ(base.imag, other.imag);
}

#define ASSERT_COMPLEX_EQUAL(base__, other__)  \
    SCOPED_TRACE("Complex values don't match"); \
    AssertMKLComplexEqual__(base__, other__)

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

// Matrix Tests

TEST(MatrixTest, CreateMatrix) {
    Matrix<double> m(2, 2);
    m(0, 0) = 1.0;
    m(1, 1) = 2.0;
    EXPECT_EQ(m(0, 0), 1.0);
    EXPECT_EQ(m(1, 1), 2.0);
}

TEST(MatrixTest, MatrixMultiplication) {
    ComplexDoubleMatrix<MKL_Complex16> m1(2, 2);
    ComplexDoubleMatrix<MKL_Complex16> m2(2, 2);
    m1(0, 0) = {1.0, 0.0};
    m1(0, 1) = {2.0, 0.0};
    m1(1, 0) = {3.0, 0.0};
    m1(1, 1) = {4.0, 0.0};

    m2(0, 0) = {5.0, 0.0};
    m2(0, 1) = {6.0, 0.0};
    m2(1, 0) = {7.0, 0.0};
    m2(1, 1) = {8.0, 0.0};

    ComplexDoubleMatrix<MKL_Complex16> result(2, 2);
    result(0, 0) = {19.0, 0.0};
    result(0, 1) = {22.0, 0.0};
    result(1, 0) = {43.0, 0.0};
    result(1, 1) = {50.0, 0.0};

    EXPECT_EQ(m1 * m2, result);
}

TEST(MatrixTest, MatrixDiagonalize){
    ComplexDoubleMatrix<MKL_Complex16> m(2, 2);
    m(0, 0) = {0.0, 0.0};
    m(0, 1) = {1.0, 0.0};
    m(1, 0) = {1.0, 0.0};
    m(1, 1) = {0.0, 0.0};

    auto eigensystem = m.hermitian_diagonalize();
    Matrix<double> result(2, 1);
    result(0, 0) = -1.0;
    result(1, 0) = 1.0;

    EXPECT_EQ(eigensystem.second, result);
}

TEST(MatrixTest, MatrixNorm){
    ComplexDoubleMatrix<MKL_Complex16> m(3, 3);
    m(0, 0) = {1.0, 0.0}, m(0, 1) = {2.0, 0.0}, m(0, 2) = {3.0, 0.0};
    m(1, 0) = {4.0, 0.0}, m(1, 1) = {5.0, 0.0}, m(1, 2) = {6.0, 0.0};
    m(2, 0) = {7.0, 0.0}, m(2, 1) = {8.0, 0.0}, m(2, 2) = {9.0, 0.0};

    Matrix<double> result = m.vector_norm();
    EXPECT_EQ(result(0, 0), sqrt(66.0));
    EXPECT_EQ(result(0, 1), sqrt(93.0));
    EXPECT_EQ(result(0, 2), sqrt(126.0));

    Matrix<double> result2 = m.vector_norm(0, 2);
    EXPECT_EQ(result2(0, 0), sqrt(17.0));
    EXPECT_EQ(result2(0, 1), sqrt(29.0));
    EXPECT_EQ(result2(0, 2), sqrt(45.0));
}

TEST(MatrixTest, MatrixCopyRow) {
    ComplexDoubleMatrix<MKL_Complex16> m1(2, 2);
    ComplexDoubleMatrix<MKL_Complex16> m2(2, 2);
    m1(0, 0) = {1.0, 0.0};
    m1(0, 1) = {2.0, 0.0};
    m1(1, 0) = {3.0, 0.0};
    m1(1, 1) = {4.0, 0.0};

    // copy row 0 of m1 to row 1 of m2
    m2.copy_row(m2, m1, 1, 0);

    EXPECT_EQ(m2(1, 0), m1(0, 0));
    EXPECT_EQ(m2(1, 1), m1(0, 1));
}

TEST(MatrixTest, MatrixCopyRow2) {
    ComplexDoubleMatrix<MKL_Complex16> m1(4, 4);
    ComplexDoubleMatrix<MKL_Complex16> m2(2, 2);
    m1(0, 0) = {1.0, 0.0};
    m1(0, 1) = {2.0, 0.0};
    m1(0, 2) = {3.0, 0.0};
    m1(0, 3) = {4.0, 0.0};
    m1(1, 0) = {3.0, 0.0};
    m1(1, 1) = {4.0, 0.0};
    m1(1, 2) = {5.0, 0.0};
    m1(1, 3) = {6.0, 0.0};

    // copy row 0 of m1 to row 1 of m2
    // only first two columns
    m2.copy_row(m2, m1, 1, 0, 2);

    EXPECT_EQ(m2(1, 0), m1(0, 0));
    EXPECT_EQ(m2(1, 1), m1(0, 1));
}

TEST(MatrixTest, MatrixInverse) {
    ComplexDoubleMatrix<MKL_Complex16> m(2, 2);
    m(0, 0) = {1.0, 0.0};
    m(0, 1) = {2.0, 0.0};
    m(1, 0) = {3.0, 0.0};
    m(1, 1) = {4.0, 0.0};

    ComplexDoubleMatrix<MKL_Complex16> inv = m.compute_inverse();
    std::cout << "Inverse matrix:\n";   
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << inv(i, j).real << " ";
        }
        std::cout << "\n";
    }
    ASSERT_COMPLEX_EQUAL(inv(0, 0), MKL_Complex16({-2.0, 0.0}));
    ASSERT_COMPLEX_EQUAL(inv(0, 1), MKL_Complex16({1.0, 0.0}));
    ASSERT_COMPLEX_EQUAL(inv(1, 0), MKL_Complex16({1.5, 0.0}));
    ASSERT_COMPLEX_EQUAL(inv(1, 1), MKL_Complex16({-0.5, 0.0}));
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

// // ProjectedState Tests

// class ProjectedStateTest : public ::testing::Test {

// protected:
//     void SetUp() override {
//         chainLatticeHalf = Lattice(Lattice_type_t::CHAIN, vec3<int>(2, 1, 1), vec3<int>(0, 0, 0));
//         chainLatticeOne = Lattice(Lattice_type_t::CHAIN, vec3<int>(3, 1, 1), vec3<int>(1, 0, 0));
        
//         chain_half_options = read_json_wavefunction_from_dir(EXAMPLES_DIR + "spin_half/1d/trivial");
//         chain_one_options = read_json_wavefunction_from_dir(EXAMPLES_DIR + "spin_one/1d/trivial");

//         mf_ansatz_half = std::shared_ptr<MeanFieldAnsatz>(new MeanFieldAnsatz_HALF(chain_half_options, chainLatticeHalf));
//         mf_ansatz_one = std::shared_ptr<MeanFieldAnsatz>(new MeanFieldAnsatz_ONE(chain_one_options, chainLatticeOne));
        
//     }

//     Lattice chainLatticeHalf, chainLatticeOne;

//     WavefunctionOptions chain_half_options, chain_one_options;

//     std::shared_ptr<MeanFieldAnsatz> mf_ansatz_half, mf_ansatz_one;
// };

// TEST_F(ProjectedStateTest, CheckFixture) {
//     EXPECT_EQ(chainLatticeHalf.get_N(), 2);
//     EXPECT_EQ(chainLatticeOne.get_N(), 3);

//     RandomEngine r(0, chainLatticeHalf.get_N(), chainLatticeHalf.get_neighbor_counts()[0]);
//     ProjectedState wf_half(*mf_ansatz_half, r);
//     std::vector<int> flips({0, 1});
// 	std::cout << "Starting with psi = " << wf_half.get_det() << "\n";
//     std::cout << "configuration = " << vec2str(wf_half.get_configuration()) << "\n";
//     std::cout << "ratio = " << wf_half.psi_over_psi(flips) << "\n";
//     wf_half.update(flips);
// 	std::cout << "After update psi = " << wf_half.get_det() << "\n";
//     std::cout << "configuration = " << vec2str(wf_half.get_configuration()) << "\n";
//     wf_half.print_matrix("Slater");
//     wf_half.print_matrix("LU");
//     wf_half.print_matrix("Winv");
//     wf_half.print_matrix("UP1");
//     wf_half.print_matrix("UP2");
//     wf_half.print_matrix("UP3");
//     wf_half.print_matrix("ipiv");

//     RandomEngine r2(0, chainLatticeOne.get_N(), chainLatticeOne.get_neighbor_counts()[0]);
//     ProjectedState wf_one(*mf_ansatz_one, r2);
//     std::vector<int> flips1({0, 1});
// 	std::cout << "Starting with psi = " << wf_one.get_det() << "\n";
//     std::cout << "ratio = " << wf_one.psi_over_psi(flips1) << "\n";
//     wf_one.update(flips1);
// 	std::cout << "After update psi = " << wf_one.get_det() << "\n";
//     std::cout << "configuration = " << vec2str(wf_one.get_configuration()) << "\n";
//     wf_one.print_matrix("Slater");
//     wf_one.print_matrix("LU");
//     wf_one.print_matrix("Winv");
//     wf_one.print_matrix("UP1");
//     wf_one.print_matrix("UP2");
//     wf_one.print_matrix("UP3");
//     wf_one.print_matrix("ipiv");
//     EXPECT_TRUE(true);
// }
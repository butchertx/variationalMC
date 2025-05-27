#pragma once
#include <complex>
#include <vector>
#include <assert.h>
#include <numeric>
#include "RandomEngine.h"
#include "MeanFieldAnsatz.h"
#include "Wavefunction.h"
#include "Lattice.h"
#include "vmc_io.h"
#include "mkl_types.h"

class ProjectedState : public Wavefunction {

	// This uses MKL_Complex16 for complex numbers, which is compatible with the Intel MKL library.
	// The interface converts them back to std::complex<double> when needed for external access.

	MeanFieldAnsatz& ansatz;
	RandomEngine& rand;
	JastrowTable jastrow;
	int N, DIM; // number of sites/particles, and state space dimension
	std::vector<int> parton_labels;
	ComplexDoubleMatrix<MKL_Complex16> Slater, Winv;
	static const int CONFIG_ATTEMPTS = 50;

	// helpers
	int Spin_t_to_row(int spin_idx);

	// initialization
	void clear_matrices();
	void initialize_configuration();
	bool try_configuration();
	void set_configuration(std::vector<int> conf);

	// updates
	void update_(std::vector<int>& flips, std::vector<int>& new_sz);
	void update(int site1, int site2);
	void updateMatrixInverse(int, int, int, int);

	// matrix elements
	MKL_Complex16 psi_over_psi2(int site1, int site2, int new_sz1, int new_sz2); // swap 2 sites with specified sz values, with jastrow
	MKL_Complex16 psi_over_psi_swap(int site1, int site2, int site3); // 3-site ring exchange, with jastrow

public:

	ProjectedState(MeanFieldAnsatz& M, RandomEngine& rand_in);

	ProjectedState(MeanFieldAnsatz& M, RandomEngine& rand_in, JastrowTable jastrow_in);

	~ProjectedState() {};

	// printing
	void print_matrix(std::string name);

	void print_timers() {
		if (jastrow.exist()) {
			jastrow.print_timers();
		}
	}

	// Override Parent Virtual Functions
	
	void f() override {};
	std::complex<double> basis_element(const std::vector<int>&) override { return { 0.0, 0.0 }; }

	// can swap spins at 2 or 3 sites given in ring_swap
	std::complex<double> psi_over_psi(std::vector<int>& ring_swap) override;
	// chooses a ring swap or a 2-site swap, potentially with an additional spin flip for the 2-site swap
	std::complex<double> psi_over_psi(std::vector<int>& flips, std::vector<int>& new_sz) override;

	void update(std::vector<int>& ring_swap) override;
	void update(std::vector<int>& flips, std::vector<int>& new_sz) override;	

	std::vector<double> log_derivative() override { return jastrow.log_derivative(); }
	std::vector<double> greedy_log_derivative() override { return jastrow.greedy_log_derivative(configuration);	}
	void update_parameters(std::vector<double> new_params) override { jastrow.set_params(new_params); }
	std::vector<double> get_parameters() override { return jastrow.get_params(); }
	
	// TODO: this should be moved to VariationalMonteCarlo and should only write the Markov chain updates
	void write_configuration(std::ofstream* f) override {
		*f << configuration[0];
		for (int i = 1; i < configuration.size(); ++i) {
			*f << "," << configuration[i];
		}
	}

	//Additional Functions

	std::complex<double> basis_element(int site, int sz) {
		// to implement this for real, first we retrieve the current determinant
		// then if conf[site] == sz, we return the current determinant
		// if conf[site] != sz, compute the ratio of determinants for the swap and
		// return the ratio times the jastrow factor times the current determinant
		return { 1.0, 0.0 };
	}

	//Tests
	// TODO: move these to an appropriate test module

	bool test_2_spin_swap_pop(bool);

	bool test_2_spin_flip_pop(std::vector<int>& flips, std::vector<int>& new_sz);

	bool test_3_spin_swap_pop(bool);
	
};

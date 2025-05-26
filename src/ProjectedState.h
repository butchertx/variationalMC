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

	MeanFieldAnsatz& ansatz;
	RandomEngine& rand;
	JastrowTable jastrow;
	int N, DIM; // number of sites/particles, and state space dimension
	std::vector<int> parton_labels;
	ComplexDoubleMatrix<MKL_Complex16> Slater, Winv;
	lapack_complex_double *UP1, *UP2, *UP3;
	static const int CONFIG_ATTEMPTS = 50;

	// helpers
	int Spin_t_to_row(int spin_idx);

	// initialization
	void malloc_matrices() {
		UP1 = (lapack_complex_double*)mkl_malloc(DIM * 2 * sizeof(lapack_complex_double), 64);
		UP2 = (lapack_complex_double*)mkl_malloc(N * 2 * sizeof(lapack_complex_double), 64);
		UP3 = (lapack_complex_double*)mkl_malloc(N * 2 * sizeof(lapack_complex_double), 64);
	}
	void clear_matrices();
	void initialize_configuration();
	bool try_configuration();
	void set_configuration(std::vector<int> conf);
	MKL_Complex16 calc_det();

	// updates
	void update(std::vector<int>& flips, std::vector<int>& new_sz, MKL_Complex16 pop);
	void update(int site1, int site2, MKL_Complex16 pop);
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
	MKL_Complex16 basis_element(const std::vector<int>&) override { return { 0.0, 0.0 }; }

	// can swap spins at 2 or 3 sites given in ring_swap
	MKL_Complex16 psi_over_psi(std::vector<int>& ring_swap) override;
	// chooses a ring swap or a 2-site swap, potentially with an additional spin flip for the 2-site swap
	MKL_Complex16 psi_over_psi(std::vector<int>& flips, std::vector<int>& new_sz) override;

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

	MKL_Complex16 basis_element(int site, int sz) {
		return { 1.0, 0.0 };
	}	

	MKL_Complex16 get_det() {
		return det;
	}

	//Tests
	// TODO: move these to an appropriate test module

	bool test_2_spin_swap_pop(bool);

	bool test_2_spin_flip_pop(std::vector<int>& flips, std::vector<int>& new_sz);

	bool test_3_spin_swap_pop(bool);
	
};

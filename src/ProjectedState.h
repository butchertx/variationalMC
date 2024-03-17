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
	lapack_complex_double *Slater, *LU, *Winv, *UP1, *UP2, *UP3;
	int *ipiv;
	std::complex<double> det;
	static const int CONFIG_ATTEMPTS = 50;

	// helpers
	int Spin_t_to_row(int spin_idx);

	// initialization
	void malloc_matrices() {
		Slater = (lapack_complex_double*)mkl_malloc(N * N * sizeof(lapack_complex_double), 64);
		LU = (lapack_complex_double*)mkl_malloc(N * N * sizeof(lapack_complex_double), 64);
		Winv = (lapack_complex_double*)mkl_malloc(DIM * N * sizeof(lapack_complex_double), 64);
		UP1 = (lapack_complex_double*)mkl_malloc(DIM * 2 * sizeof(lapack_complex_double), 64);
		UP2 = (lapack_complex_double*)mkl_malloc(N * 2 * sizeof(lapack_complex_double), 64);
		UP3 = (lapack_complex_double*)mkl_malloc(N * 2 * sizeof(lapack_complex_double), 64);
		ipiv = (int *)mkl_malloc(N * N * sizeof(int), 64);
	}
	void clear_matrices();
	void initialize_configuration();
	bool try_configuration();
	void set_configuration(std::vector<int> conf);
	std::complex<double> calc_det();

	// updates
	void update(std::vector<int>& flips, std::vector<int>& new_sz, std::complex<double> pop);
	void update(int site1, int site2, std::complex<double> pop);
	void upinvhop2(int, int, int, int);

	// matrix elements
	std::complex<double> psi_over_psi2(int site1, int site2, int new_sz1, int new_sz2); // swap 2 sites with specified sz values, with jastrow
	std::complex<double> psi_over_psi_swap(int site1, int site2, int site3); // 3-site ring exchange, with jastrow

public:

	ProjectedState(MeanFieldAnsatz& M, RandomEngine& rand_in);

	ProjectedState(MeanFieldAnsatz& M, RandomEngine& rand_in, JastrowTable jastrow_in);

	~ProjectedState() {
		mkl_free(Slater);
		mkl_free(LU);
		mkl_free(Winv);
		mkl_free(UP1);
		mkl_free(UP2);
		mkl_free(UP3);
		mkl_free(ipiv);
	}

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

	std::vector<int>& state_ref() {
		return configuration;
	}

	std::complex<double> basis_element(int site, int sz) {
		return { 1.0, 0.0 };
	}	

	std::complex<double> get_det() {
		return det;
	}

	//Tests
	// TODO: move these to an appropriate test module

	bool test_2_spin_swap_pop(bool);

	bool test_2_spin_flip_pop(std::vector<int>& flips, std::vector<int>& new_sz);

	bool test_3_spin_swap_pop(bool);
	
};

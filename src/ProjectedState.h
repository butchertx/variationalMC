#pragma once
#include <complex>
#include <vector>
#include <assert.h>
#include <numeric>
#include "RandomEngine.h"
#include "MeanFieldAnsatz.h"
#include "Wavefunction.h"
#include "Lattice.h"

class ProjectedState : public Wavefunction {

	MeanFieldAnsatz& ansatz;
	RandomEngine& rand;
	JastrowTable jastrow;
	MKL_INT N, DIM; // number of sites/particles, and state space dimension
	std::vector<int> parton_labels;
	lapack_complex_double *Slater, *LU, *Winv, *UP1, *UP2, *UP3;
	lapack_int *ipiv;
	std::complex<double> det;
	static const int CONFIG_ATTEMPTS = 50;

	// helpers
	int Spin_t_to_row(int spin_idx);

	// initialization
	void malloc_matrices();
	void clear_matrices();
	void initialize_configuration();
	bool try_configuration();
	void set_configuration(std::vector<int> conf);
	std::complex<double> calc_det();

	// updates
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

	void print_timers() {
		if (jastrow.exist()) {
			jastrow.print_timers();
		}
	}

	// Override Parent Virtual Functions
	
	void f() override {};

	std::complex<double> basis_element(const std::vector<int>&) override { return { 0.0, 0.0 }; }

	std::complex<double> psi_over_psi(std::vector<int>& ring_swap) override {
		std::vector<int> sz(ring_swap.size());
		for (int i = 0; i < ring_swap.size(); ++i) {
			if (i == ring_swap.size() - 1) {
				sz[i] = configuration[ring_swap[0]];
			}
			else {
				sz[i] = configuration[ring_swap[i + 1]];
			}
		}
		return psi_over_psi(ring_swap, sz);
	}

	std::complex<double> psi_over_psi(std::vector<int>& flips, std::vector<int>& new_sz) override;

	void update(std::vector<int>& flips, std::vector<int>& new_sz) override;
	void update(std::vector<int>& flips, std::vector<int>& new_sz, std::complex<double> pop);

	void update(std::vector<int>& ring_swap) override {
		std::vector<int> sz(ring_swap.size());
		for (int i = 0; i < ring_swap.size(); ++i) {
			if (i == ring_swap.size() - 1) {
				sz[i] = configuration[ring_swap[0]];
			}
			else {
				sz[i] = configuration[ring_swap[i + 1]];
			}
		}
		update(ring_swap, sz);
	}

	std::vector<double> log_derivative() override { 
		return jastrow.log_derivative(); 
	}

	std::vector<double> greedy_log_derivative() override {
		return jastrow.greedy_log_derivative(configuration);
	}

	void update_parameters(std::vector<double> new_params) override {
		jastrow.set_params(new_params);
	};

	std::vector<double> get_parameters() override {
		return jastrow.get_params();
	};
	
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

	void update(int, int, std::complex<double>);

	//Tests

	bool test_2_spin_swap_pop(bool);

	bool test_2_spin_flip_pop(std::vector<int>& flips, std::vector<int>& new_sz);

	bool test_3_spin_swap_pop(bool);
	
};

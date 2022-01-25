#pragma once
#include <complex>
#include <vector>
#include <assert.h>
#include <numeric>
#include "RandomEngine.h"
#include "MeanFieldAnsatz.h"
#include "Wavefunction.h"
//#include "su3_irreps.h"
#include "Lattice.h"


//class SU3RingJastrow : public Jastrow<double> {
//
//	RingList rings;
//	double logpsi_value;
//
//public:
//
//	SU3RingJastrow(RingList rings_in) : rings(rings_in) , logpsi_value(0.0) {
//		num_params = 6;
//		for (int i = 0; i < 3; ++i) {
//			param_list.push_back({});
//			param_list[i].push_back(0.0);//uptri
//			param_list[i].push_back(0.0);//downtri
//		}
//	}
//
//	double logpsi_over_psi(std::vector<int> flips, std::vector<int> new_sz, const std::vector<int>& conf) {
//		return logpsi_over_psi(flips, conf);
//	}
//
//	double logpsi_over_psi(std::vector<int> flips, const std::vector<int>& conf) {
//		//assert(flips.size() == 3);
//		double oldlp = 0.0, newlp = 0.0;
//		std::vector<int> tri_list;
//		std::vector<int> tri_labels(2);
//		bool found = false;
//		for (int site : flips) {
//			tri_labels = rings.get_ring_labels(site);
//			for (int label : tri_labels) {
//				tri_list.push_back(label);
//			}
//		}
//		std::sort(tri_list.begin(), tri_list.end());
//		auto last = std::unique(tri_list.begin(), tri_list.end());
//		tri_list.erase(last, tri_list.end());
//		//assert(tri_list.size() == 4);
//		std::vector<int> sites;
//		std::vector<int> newconf(3);
//		found = false;
//		for (int label : tri_list) {
//			sites = rings.get_ring(label);
//			for (int i = 0; i < sites.size(); ++i) {
//				found = false;
//				for (int j = 0; j < flips.size()-1; ++j) {
//					if (sites[i] == flips[j]) {
//						newconf[i] = conf[flips[j+1]] + 1;
//						found = true;
//					}
//				}
//				if (sites[i] == flips[flips.size()-1]) {
//					newconf[i] = conf[flips[0]] + 1;
//					found = true;
//				}
//				if(!found) {
//					newconf[i] = conf[sites[i]] + 1;
//				}
//			}
//			newlp += logpsi_spins(newconf[0], newconf[1], newconf[2], rings.get_ring_name(label));
//			oldlp += logpsi_tri(label, conf);
//		}
//		return newlp - oldlp;
//	}
//
//	double logpsi_spins(int s1, int s2, int s3, std::string updown) {
//		double param = 0.0;
//		double exp_factor = 0.0;
//		std::complex<double> overlap;
//		if (updown.compare("up") == 0) {
//			param = param_list[0][0]; //singlet
//			overlap = overlap_table_singlet[s1][s2][s3];
//			exp_factor += param * std::abs(overlap) * std::abs(overlap);
//
//			param = param_list[1][0]; //8
//			overlap = overlap_table_left8[s1][s2][s3];
//			exp_factor += 2.0 * param * std::abs(overlap) * std::abs(overlap);
//
//			param = param_list[2][0]; //10
//			overlap = overlap_table10[s1][s2][s3];
//			exp_factor += param * std::abs(overlap) * std::abs(overlap);
//		}
//		else if (updown.compare("down") == 0) {
//			param = param_list[0][1]; //singlet
//			overlap = overlap_table_singlet[s1][s2][s3];
//			exp_factor += param * std::abs(overlap) * std::abs(overlap);
//
//			param = param_list[1][1]; //8
//			overlap = overlap_table_left8[s1][s2][s3];
//			exp_factor += 2.0 * param * std::abs(overlap) * std::abs(overlap);
//
//			param = param_list[2][1]; //10
//			overlap = overlap_table10[s1][s2][s3];
//			exp_factor += param * std::abs(overlap) * std::abs(overlap);
//		}
//		else {
//			std::cout << "Error: no valid ring name for ring label \n";
//			exit(0);
//		}
//		return exp_factor;
//	}
//
//	double logpsi_tri(int tri_label, const std::vector<int>& conf) {
//		int s1, s2, s3;
//		std::vector<int> ring;
//		ring = rings.get_ring(tri_label);
//		s1 = conf[ring[0]] + 1, s2 = conf[ring[1]] + 1, s3 = conf[ring[2]] + 1;
//		
//		return logpsi_spins(s1, s2, s3, rings.get_ring_name(tri_label));
//	}
//
//	double logpsi(const std::vector<int>& conf) {
//		double exp_factor = 0.0;
//		for (int r = 0; r < rings.get_size(); ++r) {
//			exp_factor += logpsi_tri(r, conf);
//		}
//		return exp_factor;
//	}
//
//	void set_logpsi(const std::vector<int>& conf) {
//		logpsi_value = logpsi(conf);
//	}
//
//	std::vector<double> dlogpsi(const std::vector<int>& conf) { return {}; }
//
//};


class ProjectedState : public Wavefunction {

	MeanFieldAnsatz& ansatz;
	RandomEngine& rand;
	JastrowTable jastrow;
	int N;
	std::vector<int> parton_labels;
	//std::vector<int> configuration;
	lapack_complex_double *Slater, *LU, *Winv, *UP1, *UP2, *UP3;
	lapack_int *ipiv;
	std::complex<double> det;
	std::vector<int> conf_copy;

	void initialize_matrices();
	std::complex<double> calc_det();
	void upinvhop2(int, int, int, int);
	void upinvhop2_flip(int, int, int, int);
	void set_configuration(std::vector<int> conf);
	bool try_configuration();
	std::complex<double> psi_over_psi2(int site1, int site2, int new_sz1, int new_sz2);
	std::complex<double> psi_over_psi_swap(int site1, int site2, int site3);
	std::complex<double> psi_over_psi(int, int, int);

public:

	void print_timers() {
		if (jastrow.exist()) {
			jastrow.print_timers();
		}
	}

	~ProjectedState() {
		mkl_free(Slater);
		mkl_free(LU);
		mkl_free(Winv);
		mkl_free(UP1);
		mkl_free(UP2);
		mkl_free(UP3);
		mkl_free(ipiv);
	}

	ProjectedState(MeanFieldAnsatz& M, RandomEngine& rand_in);

	ProjectedState(MeanFieldAnsatz& M, RandomEngine& rand_in, JastrowTable jastrow_in);

	void reset_configuration() {
		for (int i = 0; i < N * N; ++i) {
			Slater[i] = { 0,0 };
			LU[i] = { 0,0 };
			ipiv[i] = 0;
			for (int j = 0; j < 3; ++j) {
				Winv[i + j * N * N] = { 0,0 };
			}
			if (i < 3 * N * 2) {
				UP1[i] = { 0.0, 0.0 };
				if (i < N * 2) {
					UP2[i] = { 0.0, 0.0 };
					UP3[i] = { 0.0, 0.0 };
				}
			}
		}
		int config_attempt = 0;
		while (!try_configuration() && config_attempt < 50) {
			det = { 0, 0 };
			++config_attempt;
		}
		assert(config_attempt < 50);
		if (jastrow.exist()) {
			jastrow.initialize_tables(configuration);
		}
	}

	//Overload Parent Virtual Functions
	
	void f() {};

	std::complex<double> basis_element(const std::vector<int>&) { return { 0.0, 0.0 }; }

	std::complex<double> psi_over_psi(std::vector<int>& ring_swap) {
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

	std::complex<double> psi_over_psi(std::vector<int>& flips, std::vector<int>& new_sz);

	void update(std::vector<int>& flips, std::vector<int>& new_sz);
	void update(std::vector<int>& flips, std::vector<int>& new_sz, std::complex<double> pop);

	void update(std::vector<int>& ring_swap) {
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

	std::vector<double> log_derivative() { 
		return jastrow.log_derivative(); 
	}

	void update_parameters(std::vector<double> new_params) {
		jastrow.set_params(new_params);
	};

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

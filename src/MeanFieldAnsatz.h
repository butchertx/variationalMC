#pragma once
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !_USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <complex>
#ifndef MKL_Complex16
#define MKL_Complex16 std::complex<double>
#endif // !MKL_Complex16
#include "mkl.h"
#include <assert.h>
#include "vmctype.h"

using namespace vmctype;

class Lattice;

class RandomEngine;

class TightBindingSitePair {

protected:
	int x1, x2;//site indices (ordered)

public:
	TightBindingSitePair(int x1_in, int x2_in)
		: x1(x1_in), x2(x2_in) {};

	void set_sites(int first, int second) {
		x1 = first;
		x2 = second;
	}

	void get_sites(int& first, int& second) {
		first = std::abs(x1);
		second = std::abs(x2);
	}

	virtual void get_couplings(double&, double&) = 0;

	virtual void get_couplings(double&, double&, double&, double&) = 0;

	virtual void conjugate() = 0;

	virtual std::string to_string() = 0;

};

class TightBindingSitePair_HALF : public TightBindingSitePair {

	double t;//coupling
	double t_phase;//phase (in units of 2pi, ordered)

public:
	TightBindingSitePair_HALF(int x1_in, int x2_in)
		: TightBindingSitePair(x1_in, x2_in), t(0.0), t_phase(0.0) {}

	TightBindingSitePair_HALF(int x1_in, int x2_in, double t_in, double t_phase_in)
		: TightBindingSitePair(x1_in, x2_in), t(t_in), t_phase(t_phase_in) {}

	void get_couplings(double& t, double& t_phase) override;

	void get_couplings(double&, double&, double&, double&) override{
		throw std::invalid_argument("Spin half hoppings have 2 parameters");
	};

	void conjugate() override;

	std::string to_string() override;
};

class TightBindingSitePair_ONE : public TightBindingSitePair {

	double tz, txy;//couplings
	double tz_phase, txy_phase;//phases (in units of 2pi, ordered)

public:
	TightBindingSitePair_ONE(int x1_in, int x2_in)
		: TightBindingSitePair(x1_in, x2_in), tz(0.0), txy(0.0), tz_phase(0.0), txy_phase(0.0) {}

	TightBindingSitePair_ONE(int x1_in, int x2_in, double tz_in, double txy_in, double tz_phase_in, double txy_phase_in)
		: TightBindingSitePair(x1_in, x2_in), tz(tz_in), txy(txy_in), tz_phase(tz_phase_in), txy_phase(txy_phase_in) {}

	void get_couplings(double&, double&) override {
		throw std::invalid_argument("Spin one hoppings have 4 parameters");
	};

	void get_couplings(double& tzr, double& tzi, double& txyr, double& txyi) override;

	void conjugate() override;

	std::string to_string() override;
};

class TightBindingUnitCell {
public:
	std::vector<std::unique_ptr<TightBindingSitePair>> site_pairs;
};

class Orbital {

	std::vector<double> occupation_numbers;
	vec3<double> k;
	double energy;
	int index;

public:

	Orbital() {};
	Orbital(std::vector<double> occ_in, double e_in, int ind_in)
		: occupation_numbers(occ_in), energy(e_in), index(ind_in) {};

	int get_index() {
		return index;
	};

	void update_index(int new_index) {
		index = new_index;
	}

	double get_overlap(int sz) {
		return occupation_numbers[1 - sz];
	}
};

class FermiSurface {

	std::vector<double> inner_shell;//number of each flavor in the inner shell
	std::vector<Orbital> orbitals;
	double energy;

public:

	FermiSurface() {};
	FermiSurface(double e_in, std::vector<double> shell_in) 
	: inner_shell(shell_in), energy(e_in) {};

	void update_index(int fs_index, int new_orb_index) {
		orbitals[fs_index].update_index(new_orb_index);
	}

	void add_orbital(Orbital o) { 
		orbitals.push_back(o);
	}

	Orbital get_orbital(int fs_index) {
		return orbitals[fs_index];
	}

	std::vector<Orbital> get_orbitals() {
		return orbitals;
	}

	int get_size() {
		return orbitals.size();
	}

	double get_inner_shell_count(int sz) {
		return inner_shell[1 - sz];
	}
};

class MeanFieldAnsatz {

protected:

	int N, DIM = 0, info, fermi_surface_start, fermi_surface_end;
	vmctype::Spin_t SPIN_TYPE;
	double field;
	lapack_complex_double *HMF, *Phi; // , * Pair_Eig, * PhiR;
	std::vector<lapack_complex_double*> del_H; //each element corresponds to dH for a given variational param
	double *Energy;
	std::vector<std::vector<std::shared_ptr<TightBindingSitePair>>> site_pair_list;//each hopping vmc param has its own vector of site pairs
	std::vector<std::vector<std::complex<double>>> mean_field_hamiltonian;
	FermiSurface fermi;

	virtual void set_hamiltonian() = 0;
	void diagonalize_hamiltonian();

public:

	// Implemented functions

	MeanFieldAnsatz(int N_in, double field_in) : N(N_in), field(field_in) {};

	~MeanFieldAnsatz() {
		mkl_free(HMF);
		mkl_free(Phi);
		mkl_free(Energy);
		for (auto p : del_H) {
			mkl_free(p);
		}
	}

	lapack_complex_double* get_H() { return HMF; }

	lapack_complex_double* get_Phi() { return Phi; }

	double* get_Energy() { return Energy; }

	int get_N() { return N; }
	
	int get_dim() { return DIM; }

	vmctype::Spin_t get_spin_type() { return SPIN_TYPE; }

	int get_num_hop_classes() { return site_pair_list.size(); }

	std::vector<std::pair<int,int>> get_tb_pairs(int hop_class);

	std::string get_tb_string();

	// pure functions

	virtual void print_levels(bool print_all = false) = 0;

	virtual int get_N0F() = 0;

};

class MeanFieldAnsatz_HALF : public MeanFieldAnsatz {

	vmctype::SpecificWFOptions opts;

	virtual void set_hamiltonian() override;

public:

	MeanFieldAnsatz_HALF(WavefunctionOptions& mf_in, Lattice& lat_in);

	virtual void print_levels(bool) override;

	virtual int get_N0F() override { return 0; }

};

class MeanFieldAnsatz_ONE : public MeanFieldAnsatz {

	int n0_F;
	vmctype::SpecificWFOptions opts;
	std::vector<vec3<std::complex<double>>> directors;

	virtual void set_hamiltonian() override;
	std::complex<double> get_director_element(vec3<std::complex<double>>, int m1, int m2);
	double get_su3_element(std::string, int, int, int);

public:

	MeanFieldAnsatz_ONE(WavefunctionOptions& mf_in, Lattice& lat_in);

	virtual void print_levels(bool) override;
	void print_fermi_level();
	void write_levels(std::ofstream* f);
	void write_directors(std::ofstream* f);
	void set_fermi_surface();

	void shuffle_FS(int n0, int n1, RandomEngine* rand);

	bool check_orbital_overlap(int fs_index, int sz) {
		return std::abs(fermi.get_orbital(fs_index).get_overlap(sz)) > EPSILON;
	}

	int get_orbital_index(int fs_index) {
		return fermi.get_orbital(fs_index).get_index();
	}

	int get_N0F() override {
		if (opts.su3_symmetry) {
			assert(3 * (N / 3) == N);
			return N/3;
		}
		else {
			return n0_F;
		}
	}

	int get_fermi_surface_start() {
		return fermi_surface_start;
	}

	int get_fermi_surface_end() {
		return fermi_surface_end;
	}

};


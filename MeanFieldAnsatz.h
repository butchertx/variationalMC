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
#include "mkl.h" // "C:\Program Files (x86)\Intel\oneAPI\mkl\2021.1.1\include\mkl.h"
#include <assert.h>
#include "vmctype.h"

using namespace vmctype;

class Lattice;

class RandomEngine;

class TightBindingSitePair {
	int x1, x2;//site indices (ordered)
	double tz, txy;//couplings
	double tz_phase, txy_phase;//phases (in units of 2pi, ordered)

public:
	TightBindingSitePair(int x1_in, int x2_in)
		: x1(x1_in), x2(x2_in) {
		tz = 0.0;
		txy = 0.0;
		tz_phase = 0.0;
		txy_phase = 0.0;
	}

	TightBindingSitePair(int x1_in, int x2_in, double tz_in, double txy_in, double tz_phase_in, double txy_phase_in)
		: x1(x1_in), x2(x2_in), tz(tz_in), txy(txy_in), tz_phase(tz_phase_in), txy_phase(txy_phase_in) {}

	void set_sites(int first, int second) {
		x1 = first;
		x2 = second;
	}

	void get_sites(int * first, int * second) {
		*first = std::abs(x1);
		*second = std::abs(x2);
	}

	void get_couplings(double * tzr, double * tzi, double * txyr, double * txyi);

	void conjugate() {
		int temp = x2;
		x2 = temp < 0 ? -x1 : x1;
		x1 = std::abs(temp);
		tz_phase *= -1;
		txy_phase *= -1;
	}

	std::string to_string() {
		std::stringstream ss;
		ss << "Sites: (" << x1 << "," << x2 << "); Tz, Txy (polar, units of 2pi) = {(" << tz << "," << tz_phase << "), (" << txy << "," << txy_phase << ")}";
		return ss.str();
	}
};

class TightBindingUnitCell {
public:
	std::vector<TightBindingSitePair> site_pairs;
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
	: energy(e_in), inner_shell(shell_in) {};

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

	int N, info, n0_F, fermi_surface_start, fermi_surface_end;
	double field;
	lapack_complex_double *HMF, *Phi; // , * Pair_Eig, * PhiR;
	std::vector<lapack_complex_double*> del_H; //each element corresponds to dH for a given variational param
	double *Energy;
	std::vector<std::vector<TightBindingSitePair>> site_pair_list;//each hopping vmc param has its own vector of site pairs
	//std::vector<std::vector<SingletPairingSitePair>> singlet_pair_list;
	std::vector<vec3<std::complex<double>>> directors;
	std::vector<std::vector<std::complex<double>>> mean_field_hamiltonian;
	FermiSurface fermi;

	void set_hamiltonian();
	void diagonalize_hamiltonian();
	std::complex<double> get_director_element(vec3<std::complex<double>>, int m1, int m2);
	double get_su3_element(std::string, int, int, int);


public:

	~MeanFieldAnsatz() {
		mkl_free(HMF);
		mkl_free(Phi);
		mkl_free(Energy);
		for (auto p : del_H) {
			mkl_free(p);
		}
	}

	//MeanFieldAnsatz(int Lx_in, int Ly_in, TightBindingUnitCell uc_in);

	MeanFieldAnsatz(mean_field_options& mf_in, Lattice& lat_in, bool unit_cell_construction);

	lapack_complex_double* get_H() {
		return HMF;
	}

	lapack_complex_double* get_Phi() {
		return Phi;
	}

	double* get_Energy() {
		return Energy;
	}

	int get_dim() {
		return 3 * N;
	}

	int get_num_hop_classes() {
		return site_pair_list.size();
	}

	std::string get_tb_string() {
		std::stringstream ss;
		ss << "Listing Hopping pairs:\n";
		for (auto tb_list : site_pair_list) {
			for (auto tb_pair : tb_list) {
				ss << tb_pair.to_string() << "\n";
			}
		}
		return ss.str();
	}

	std::vector<std::pair<int,int>> get_tb_pairs(int hop_class) {
		std::vector<std::pair<int,int>> result;
		int first, second;
		assert(hop_class < site_pair_list.size());
		for (auto tb_pair : site_pair_list[hop_class]) {
			tb_pair.get_sites(&first, &second);
			result.push_back(std::pair<int,int>(std::abs(first),std::abs(second)));
		}
		return result;
	}

	void print_levels();
	void print_fermi_level();
	void set_fermi_surface();

	void shuffle_FS(int n0, int n1, RandomEngine* rand);

	bool check_orbital_overlap(int fs_index, int sz) {
		return std::abs(fermi.get_orbital(fs_index).get_overlap(sz)) > EPSILON;
	}

	int get_orbital_index(int fs_index) {
		return fermi.get_orbital(fs_index).get_index();
	}

	int get_N0F() {
		return n0_F;
	}

	int get_fermi_surface_start() {
		return fermi_surface_start;
	}

	int get_fermi_surface_end() {
		return fermi_surface_end;
	}

};


#pragma once
#include <complex>
#include <vector>
#include <set>
#include <unordered_set>
#include <assert.h>
#include <numeric>
#include <algorithm>
#include "MemTimeTester.h"


class JastrowDensity {

	double strength;
	double sz0_density;

public:

	JastrowDensity() {};

	JastrowDensity(double strength_);

	void initialize_Sz0(std::vector<int>&);

	double greedy_eval(std::vector<int>&);

	double greedy_log_derivative(std::vector<int>&);

	double conf_sum_difference(std::vector<int>&, std::vector<int>&, std::vector<int>&);

	double lazy_eval(std::vector<int>&, std::vector<int>&, std::vector<int>&);

	void update_Sz0(std::vector<int>&, std::vector<int>&, std::vector<int>&);

	double log_derivative() {
		return sz0_density;
	}

	double get_param() {
		return strength;
	}

	void set_param(double strength_) {
		strength = strength_;
	}
};

class JastrowFactor {
	bool sz2 = false; // true if the jastrow factor couples to (S^z)^2
	double strength = 0.0;
	std::vector<std::vector<int>> neighbor_table;
	std::vector<std::vector<bool>> neighbor_bool;
	std::vector<std::set<int>> neighbor_set;
	std::vector<double> exp_table = {};
	double conf_sum = 0.0;
	MemTimeTester timer;

public:

	void print_timers() {
		timer.print_timers();
	}

	std::vector<std::vector<int>> get_neighbor_table() {
		return neighbor_table;
	}

	bool check_neighbors(int site1, int site2) {
		return neighbor_bool[site1][site2];
	}

	JastrowFactor() {};

	JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_);
	
	JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_, std::vector<int>& configuration_);

	JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_, bool sz2_flag);

	JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_, std::vector<int>& configuration_, bool sz2_flag);

	void initialize_table(std::vector<int>&);

	double greedy_eval(std::vector<int>&);

	double greedy_log_derivative(std::vector<int>&);

	double conf_sum_difference(std::vector<int>&, std::vector<int>&, std::vector<int>&, bool);

	double lazy_eval(std::vector<int>&, std::vector<int>&, std::vector<int>&);

	//double lazy_eval(std::vector<int>&, std::vector<int>&, std::vector<int>&, bool);

	void update_tables(std::vector<int>&, std::vector<int>&, std::vector<int>&);

	double log_derivative() {
		return conf_sum;
	}

	double get_param() {
		return strength;
	}

	void set_param(double strength_) {
		strength = strength_;
	}
};

class JastrowTable {

	std::vector<JastrowFactor> jastrows;
	JastrowDensity jdens;
	bool jdens_flag = false;
	MemTimeTester timer;

public:

	void print_timers() {
		std::cout << "Jastrow timers:\n";
		timer.print_timers();
	}

	JastrowTable() {};
	
	JastrowTable(std::vector<JastrowFactor> jastrows_) : jastrows(jastrows_) {
		jdens_flag = false;
	};

	JastrowTable(std::vector<JastrowFactor> jastrows_, JastrowDensity jdens_) : jastrows(jastrows_), jdens(jdens_) {
		jdens_flag = true;
	};

	void initialize_tables(std::vector<int>& configuration) {

		std::vector<std::vector<int>> neighbor_table;

		for (int i = 0; i < jastrows.size(); ++i) {

			jastrows[i].initialize_table(configuration);

		}

		if (jdens_flag) {
			jdens.initialize_Sz0(configuration);
		}
	}

	bool exist() {
		return jastrows.size() > 0 || jdens_flag;
	}

	double greedy_eval(std::vector<int>& configuration_) {
		double result = 1.0;
		for (int j = 0; j < jastrows.size(); ++j) {
			result *= jastrows[j].greedy_eval(configuration_);
		}
		if (jdens_flag) {
			result *= jdens.greedy_eval(configuration_);
		}
		return result;
	}

	double lazy_eval(std::vector<int>& a, std::vector<int>& b, std::vector<int>& c) {
		timer.flag_start_time("lazy eval");
		double result = 1.0;
		for (int j = 0; j < jastrows.size(); ++j) {
			result *= jastrows[j].lazy_eval(a, b, c);
		}
		if (jdens_flag) {
			result *= jdens.lazy_eval(a, b, c);
		}
		timer.flag_end_time("lazy eval");
		return result;
	}

	void update_tables(std::vector<int>& a, std::vector<int>& b, std::vector<int>& c) {
		timer.flag_start_time("update tables");
		for (int j = 0; j < jastrows.size(); ++j) {
			jastrows[j].update_tables(a, b, c);
		}
		if (jdens_flag) {
			jdens.update_Sz0(a, b, c);
		}
		timer.flag_end_time("update tables");
	}

	std::vector<double> log_derivative() {
		std::vector<double> result;
		if (jdens_flag) {
			result.push_back(jdens.log_derivative());
		}
		for (int j = 0; j < jastrows.size(); ++j) {
			result.push_back(jastrows[j].log_derivative());
		}
		return result;
	}

	std::vector<double> greedy_log_derivative(std::vector<int>& conf) {
		std::vector<double> result;
		if (jdens_flag) {
			result.push_back(jdens.greedy_log_derivative(conf));
		}
		for (int j = 0; j < jastrows.size(); ++j) {
			result.push_back(jastrows[j].greedy_log_derivative(conf));
		}
		return result;
	}

	std::vector<double> get_params() {
		std::vector<double> result;
		if (jdens_flag) {
			result.push_back(jdens.get_param());
		}
		for (int j = 0; j < jastrows.size(); ++j) {
			result.push_back(jastrows[j].get_param());
		}
		return result;
	}

	void set_params(std::vector<double> params_) {
		assert(params_.size() == jastrows.size() + jdens_flag);
		if (jdens_flag) {
			jdens.set_param(params_[0]);
		}
		for (int i = jdens_flag; i < params_.size(); ++i) {
			jastrows[i-jdens_flag].set_param(params_[i]);
		}
	}

};

//template <class T>
//class GreedyJastrow : public Jastrow<T> {
//
//public:
//
//	GreedyJastrow() {}
//
//	void f() {}
//
//};



class Wavefunction {
	//A Wavefunction's main role is to define how to calculate matrix elements and correlations
	//Defined on a Lattice, has a set of basis states (usually Sz eigenstates?), calculates overlaps and matrix elements for basis states
protected:
	std::vector<int> configuration;

public:

	virtual void f() = 0;

	virtual std::complex<double> basis_element(const std::vector<int>&) = 0;

	virtual std::complex<double> psi_over_psi(std::vector<int>&, std::vector<int>&) = 0;
	virtual std::complex<double> psi_over_psi(std::vector<int>&) = 0;

	virtual void update(std::vector<int>&, std::vector<int>&) = 0;
	virtual void update(std::vector<int>&) = 0;

	const std::vector<int>& conf_ref() { return configuration; }

	virtual std::vector<double> greedy_log_derivative() { return { 0.0 }; }
	virtual std::vector<double> log_derivative() { return { 0.0 }; }

	virtual void update_parameters(std::vector<double>) {};
	virtual std::vector<double> get_parameters() { return {}; };

	

	//virtual void write_amplitudes(std::ofstream*) = 0;
	//virtual void write_configuration(std::ofstream*) = 0;

};




//
//class DirectorProductState : public Wavefunction {
//
//	Lattice& lat;
//	std::vector<vec3<double>> u;
//	std::vector<vec3<double>> v;
//
//
//
//public:
//
//	std::complex<double> basis_element(int site, int sz) {
//		assert(site >= 0 && site < u.size());
//		//assert(sz == 0 || sz == -1 || sz == 1);
//
//		switch (sz) {
//		case 1:
//			return 1 / sqrt(2) * (u[site].y - v[site].x, u[site].x + v[site].y);
//		case -1:
//			return 1 / sqrt(2) * (u[site].y + v[site].x, -u[site].x + v[site].y);
//		case 0:
//			return (v[site].z, -u[site].z);
//		default:
//			return (0.0, 0.0);//could throw an error here
//		}
//
//	}
//
//	void f() {}
//
//	DirectorProductState(Lattice& lat_in, vec3<vec3<double>> Q_in, std::vector<double> basis_phases_in)
//		: lat(lat_in), u(std::vector<vec3<double>>(lat_in.get_N())), v(std::vector<vec3<double>>(lat.get_N())) {
//		//Momentum initialization
//		assert(basis_phases_in.size() == lat.get_num_base());
//	}
//
//	DirectorProductState(Lattice& lat_in, std::vector<vec3<double>> u_in, std::vector<vec3<double>> v_in)
//		: lat(lat_in) {
//		//sublattice initialization
//		assert(u_in.size() == lat.get_num_base() && v_in.size() == lat.get_num_base());
//
//		for (int cell = 0; cell < lat.get_N() / u_in.size(); ++cell) {
//			for (int base = 0; base < u_in.size(); ++base) {
//				u.push_back(u_in[base]);
//				v.push_back(v_in[base]);
//			}
//		}
//	}
//
//	DirectorProductState(Lattice& lat_in, std::vector<vec3<double>> u_in, std::vector<vec3<double>> v_in, bool greedy)
//		: lat(lat_in) {
//		//greedy initialization
//		assert(u_in.size() == lat.get_N() && v_in.size() == lat.get_N());
//
//		for (int base = 0; base < lat.get_N(); ++base) {
//			u.push_back(u_in[base]);
//			v.push_back(v_in[base]);
//		}
//
//	}
//	
//	//compute <x'|psi>/<x|psi>
//	std::complex<double> psi_over_psi(std::vector<int>& state, std::vector<int>& flips, std::vector<int>& new_sz);
//
//	std::vector<double>  expectation_Nz();
//
//	double expectation_Nz_total();
//
//};


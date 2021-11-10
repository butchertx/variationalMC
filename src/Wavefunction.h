#pragma once
#include <complex>
#include <vector>
#include <assert.h>
#include <numeric>


template <class T>//param class: double or complex
class Jastrow {

protected:

	std::vector<std::vector<T>> param_list; //each param type has its own vector (ex. sz.sz is first list, sz^2.sz^2 is second list, etc.)
	int num_params;

public:

	virtual T logpsi(const std::vector<int>&) = 0;
	virtual T logpsi_over_psi(std::vector<int>, std::vector<int>, const std::vector<int>&) = 0;
	virtual std::vector<T> dlogpsi(const std::vector<int>&) = 0;

	std::vector<T> unroll_params() {
		std::vector<T> result;
		for (auto ptype : param_list) {
			for (auto param : ptype) {
				result.push_back(param);
			}
		}
		return result;
	}

	void set_params(std::vector<T> p_in) {
		assert(p_in.size() == num_params);
		int p_idx = 0;
		std::vector<std::vector<T>> new_p_in;
		for (int i = 0; i < param_list.size(); ++i) {
			new_p_in.push_back({});
			for (int j = 0; j < param_list[i].size(); ++j) {
				new_p_in[i].push_back(p_in[p_idx]);
				++p_idx;
			}
		}
		set_params(new_p_in);
	}

	void set_params(std::vector<std::vector<T>> p_in) {
		assert(p_in.size() == param_list.size());
		for (int i = 0; i < p_in.size(); ++i) {
			assert(p_in[i].size() == param_list[i].size());
		}
		param_list = p_in;
	}

};

template <class T>
class GreedyJastrow : public Jastrow<T> {

public:

	GreedyJastrow() {}

	void f() {}

};



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

	//virtual void update(int, int, std::complex<double>) = 0;
	

	//virtual std::complex<double> psi_over_psi(int, int, int, int) = 0;

	//virtual std::complex<double> psi_over_psi_full_calculation(int, int, int, int) = 0;
	
	/*virtual void update(int, int, std::complex<double>) = 0;
	virtual void update(int, bool, std::complex<double>) = 0;*/

	

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


#pragma once
#include <math.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <assert.h>
#include <complex>
#include "Lattice.h"

class FlipList {
public:
	std::vector<std::complex<double>> multipliers;
	std::vector<std::vector<int>> flips;
	std::vector<std::vector<int>> new_sz;

	void empty() {
		multipliers.clear();
		flips.clear();
		new_sz.clear();
	}

};

static class FlipListDiag : public FlipList {

	//dummy FlipList object that doesn't have off-diagonal terms

};

class Interaction {
	//has an action on a basis state with a multiplier, diagonal, and off-diagonal elements (multiplied from the left) (Use Ising basis)
protected:

	double coefficient; // coefficient for a hamiltonian or group of operators in an order parameter

	FlipList flip_buffer;

public:
	Interaction(double coefficient_) : coefficient(coefficient_) {};

	virtual std::complex<double> diag(const std::vector<int>& state) = 0;

	virtual FlipList& off_diag(const std::vector<int>& state) = 0;

	virtual void print_info(const std::vector<int>& state) = 0;


};

class SwapExchange : public Interaction {

	int i, j;

public:

	SwapExchange(int i_in, int j_in, double coefficient_)
		: Interaction(coefficient), i(i_in), j(j_in) {
		flip_buffer.multipliers.push_back(coefficient_);
		flip_buffer.flips.push_back({ i,j });
		flip_buffer.new_sz.push_back({ 0, 0 });
	};

	std::complex<double> diag(const std::vector<int>& state) { 
		return { 0.0, 0.0 }; // coupling;
	}

	FlipList& off_diag(const std::vector<int>& state) {		
		flip_buffer.new_sz[0][0] = state[j];
		flip_buffer.new_sz[0][1] = state[i];
		return flip_buffer;
	}

	void print_info(const std::vector<int>& state) {
		std::cout << "Sites (i,j) = (" << i << "," << j << ")\n";
		std::cout << "Spins (szi, szj) = (" << state[i] << "," << state[j] << ")\n";
		std::cout << "Flips:\n";
		for (int i = 0; i < flip_buffer.flips.size(); ++i) {
			std::cout << "Mult = " << flip_buffer.multipliers[i] << "\n";
			std::cout << "Flip sites = (" << flip_buffer.flips[i][0] << "," << flip_buffer.flips[i][1] << ")\n";
		}
	}
};

class HeisenbergExchange : public Interaction {
	int i, j;
	double S;

public:

	HeisenbergExchange(int i_in, int j_in, double S_in, double coefficient_)
		: Interaction(coefficient_), i(i_in), j(j_in), S(S_in) {
		flip_buffer.multipliers.push_back(coefficient_);
		flip_buffer.multipliers.push_back(coefficient_);
		flip_buffer.flips.push_back({ i, j });
		flip_buffer.flips.push_back({ i, j });
		flip_buffer.new_sz.push_back({ 0, 0 });
		flip_buffer.new_sz.push_back({ 0, 0 });
	};

	std::complex<double> diag(const std::vector<int>& state) {
		return coefficient * S * S * state[i] * state[j];
	}

	FlipList& off_diag(const std::vector<int>& state) {
		flip_buffer.multipliers[0] = std::complex<double>{ 0.0, 0.0 };
		flip_buffer.multipliers[1] = std::complex<double>{ 0.0, 0.0 };
		if (state[i] != -S && state[j] != S) {
			flip_buffer.multipliers[0] = coefficient;
			flip_buffer.new_sz[0][0] = state[i] - 1;
			flip_buffer.new_sz[0][1] = state[j] + 1;
		}
		if (state[i] != S && state[j] != -S) {
			flip_buffer.multipliers[1] = coefficient;
			flip_buffer.new_sz[1][0] = state[i] + 1;
			flip_buffer.new_sz[1][1] = state[j] - 1;
		}

		return flip_buffer;
	}

	void print_info(const std::vector<int>& state) {
		std::cout << "Sites (i,j) = (" << i << "," << j << ")\n";
		std::cout << "Spins (szi, szj) = (" << state[i] << "," << state[j] << ")\n";
		std::cout << "Flips:\n";
		for (int i = 0; i < flip_buffer.flips.size(); ++i) {
			std::cout << "Mult = " << flip_buffer.multipliers[i] << "\n";
			std::cout << "Flip sites = (" << flip_buffer.flips[i][0] << "," << flip_buffer.flips[i][1] << ")\n";
		}
	}
};

class IsingExchange : public Interaction {

	int i, j;
	double S;

	FlipList flip_off_diag = FlipListDiag();

public:

	IsingExchange(int i_in, int j_in, double S_in)
		: Interaction(1.0), i(i_in), j(j_in), S(S_in) {};

	std::complex<double> diag(const std::vector<int>& state) {
		return S * S * state[i] * state[j];
	}

	FlipList& off_diag(const std::vector<int>& state) {
		return flip_off_diag;
	}

	void print_info(const std::vector<int>& state) {
		std::cout << "Sites (i,j) = (" << i << "," << j << ")\n";
		std::cout << "Spins (szi, szj) = (" << state[i] << "," << state[j] << ")\n";
	}
};

class LadderExchange : public Interaction {
	int i, j;
	double S;

public:

	LadderExchange(int i_in, int j_in, double S_in)
		//only implemented for spin=1
		: Interaction(1.0), i(i_in), j(j_in), S(S_in) {
		flip_buffer.multipliers.push_back(1.0);
		flip_buffer.multipliers.push_back(1.0);
		flip_buffer.flips.push_back({ i, j });
		flip_buffer.flips.push_back({ i, j });
		flip_buffer.new_sz.push_back({ 0, 0 });
		flip_buffer.new_sz.push_back({ 0, 0 });
	};

	std::complex<double> diag(const std::vector<int>& state) {
		return 0.0;
	}

	FlipList& off_diag(const std::vector<int>& state) {
		flip_buffer.multipliers[0] = std::complex<double>{ 0.0, 0.0 };
		flip_buffer.multipliers[1] = std::complex<double>{ 0.0, 0.0 };
		if (state[i] != -S && state[j] != S) {
			flip_buffer.multipliers[0] = coefficient;
			flip_buffer.new_sz[0][0] = state[i] - 1;
			flip_buffer.new_sz[0][1] = state[j] + 1;
		}
		if (state[i] != S && state[j] != -S) {
			flip_buffer.multipliers[1] = coefficient;
			flip_buffer.new_sz[1][0] = state[i] + 1;
			flip_buffer.new_sz[1][1] = state[j] - 1;
		}

		return flip_buffer;
	}

	void print_info(const std::vector<int>& state) {
		std::cout << "Sites (i,j) = (" << i << "," << j << ")\n";
		std::cout << "Spins (szi, szj) = (" << state[i] << "," << state[j] << ")\n";
		std::cout << "Flips:\n";
		for (int i = 0; i < flip_buffer.flips.size(); ++i) {
			std::cout << "Mult = " << flip_buffer.multipliers[i] << "\n";
			std::cout << "Flip sites = (" << flip_buffer.flips[i][0] << "," << flip_buffer.flips[i][1] << ")\n";
		}
	}
};

class RingExchange : public Interaction {

	int i, j, k;

	int tri_label;

	bool CC;

public:

	RingExchange( RingList& r, int tri_label_in, bool CC_in, double coefficient_)
		: Interaction(coefficient_), tri_label(tri_label_in), CC(CC_in){
		std::vector<int> sites = r.get_ring(tri_label);
		if (!CC) {
			i = sites[0];
			j = sites[1];
			k = sites[2];
		}
		else {
			i = sites[0];
			k = sites[1];
			j = sites[2];
		}
		flip_buffer.multipliers.push_back(coefficient_);
		flip_buffer.flips.push_back({ i,j,k });
		flip_buffer.new_sz.push_back({ 0, 0, 0 });
	}

	std::complex<double> diag(const std::vector<int>& state) { return { 0.0, 0.0 }; }

	FlipList& off_diag(const std::vector<int>& state) {
		//flip_buffer.empty();
		//flip_buffer.multipliers.push_back(coupling);
		//flip_buffer.flips.push_back({ i,j,k });
		flip_buffer.new_sz[0][0] = state[j];
		flip_buffer.new_sz[0][1] = state[k];
		flip_buffer.new_sz[0][2] = state[i];
		return flip_buffer;
	}

	void print_info(const std::vector<int>& state) {
		std::cout << "Sites (i,j,k) = (" << i << "," << j << "," << k << ")\n";
		std::cout << "Spins (szi, szj) = (" << state[i] << "," << state[j] << "," << state[k] << ")\n";
		std::cout << "Flips:\n";
		for (int i = 0; i < flip_buffer.flips.size(); ++i) {
			std::cout << "Mult = " << flip_buffer.multipliers[i] << "\n";
			std::cout << "Flip sites = (" << flip_buffer.flips[i][0] << "," << flip_buffer.flips[i][1] << "," << flip_buffer.flips[i][2] << ")\n";
		}
	}

	bool get_CC() {
		return CC;
	}

	int get_tri_label() {
		return tri_label;
	}

	


};

class SingleIonAnisotropy : public Interaction {

	int i;
	double S;

	FlipList flip_off_diag = FlipListDiag();

public:

	SingleIonAnisotropy(int i_in, double S_in)
		: Interaction(1.0), i(i_in), S(S_in) {};

	std::complex<double> diag(const std::vector<int>& state) {
		return S * S * state[i] * state[i];
	}

	FlipList& off_diag(const std::vector<int>& state) {
		return flip_off_diag;
	}

	void print_info(const std::vector<int>& state) {
		std::cout << "Site (i) = (" << i << ")\n";
		std::cout << "Spin (szi) = (" << state[i]  << ")\n";
	}
};

//class BilinearSpinInteraction : public Interaction {
//protected:
//	double S; //spin
//	int i, j; //two sites
//
//public:
//	BilinearSpinInteraction(double coupling_in, int i_in, int j_in, double S_in) : Interaction(coupling_in), i(i_in), j(j_in), S(S_in) {
//		assert(S == 1 / 2 || S == 1);//limit the cases for now
//	}
//
//	double diag(std::vector<int>& state) { 
//		return coupling * S * S * state[i]*state[j]; 
//	}
//
//};

//class HeisenbergInteraction : public Interaction {
//
//	double S; //spin
//	double diag_factor;
//	int i, j; //two sites
//
//public:
//	HeisenbergInteraction(double coupling_in, int i_in, int j_in, double S_in) : Interaction(coupling_in), i(i_in), j(j_in), S(S_in) {
//		diag_factor = 0.5 * coupling * S * S;
//	}
//
//	double diag(std::vector<int>& state) {
//		return diag_factor * state[i] * state[j];
//	}
//
//	FlipList off_diag(std::vector<int>& state);
//
//	void print_info(std::vector<int>& state);
//
//};
//
//class BiquadraticInteraction : public Interaction {
//
//	double S;
//	double diag_factor;
//	int i, j;
//
//public:
//	BiquadraticInteraction(double coupling_in, int i_in, int j_in, double S_in) : Interaction(coupling_in), i(i_in), j(j_in), S(S_in) {
//		diag_factor = 0.5 * coupling * S * S * S * S;
//	}
//
//	double diag(std::vector<int>& state) {
//		return (state[i] == -1 || state[j] == 1) ? 
//					diag_factor * state[i] * state[i] * state[j] * state[j] :
//					diag_factor * state[i] * state[i] * state[j] * state[j] + coupling;
//	}
//
//	FlipList off_diag(std::vector<int>& state);
//
//	void print_info(std::vector<int>& state);
//};

class Observable {

	std::string name;
	//List of interactions
	std::vector<Interaction*> interactions;

public:
	Observable(std::string name_) : name(name_) {};

	Observable() {};

	void add_interaction(Interaction * interaction_in) {
		interactions.push_back(interaction_in);
	}

	int total() {
		return interactions.size();
	}

	std::vector<Interaction*>& get() {
		return interactions;
	}
};

class SpinModel {

	// map of names to observables that add up to get the energy
	std::map<std::string, Observable> terms;
	std::map<std::string, std::complex<double>> couplings;
	std::complex<double> E0;

public:
	SpinModel() {};

	void add_term(std::string name_in, Observable term_in, std::complex<double> coupling_in) {
		terms.insert_or_assign(name_in, term_in);
		couplings.insert_or_assign(name_in, coupling_in);
	}

	void add_constant(std::complex<double> E0) {
		this->E0 = E0;
	}

	std::complex<double> get_E0() {
		return E0;
	}

	std::vector<std::string> get_terms() {
		std::vector<std::string> term_names;
		for (std::map<std::string, Observable>::iterator it = terms.begin(); it != terms.end(); ++it) {
			term_names.push_back(it->first);
		}
		return term_names;
	}

	std::vector<Interaction*>& get_interactions(std::string name_in) {
		return terms[name_in].get();
	}

	std::complex<double> get_coupling(std::string name_in) {
		return couplings[name_in];
	}
};
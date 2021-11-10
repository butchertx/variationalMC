#pragma once
#include <math.h>
#include <vector>
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

class Interaction {
	//has a coupling and an action on a basis state (multiplied from the left) (Use Ising basis for now)
protected:

	std::complex<double> coupling;

	FlipList flip_buffer;

public:
	Interaction(std::complex<double> coupling_in) : coupling(coupling_in) {}

	virtual std::complex<double> diag(const std::vector<int>& state) = 0;

	virtual FlipList off_diag(const std::vector<int>& state) = 0;

	virtual void print_info(const std::vector<int>& state) = 0;

	std::complex<double> get_coupling() {
		return coupling;
	}


};

class SwapExchange : public Interaction {

	int i, j;

public:

	SwapExchange(std::complex<double> coupling_in, int i_in, int j_in)
		: Interaction(coupling_in), i(i_in), j(j_in) {
		flip_buffer.multipliers.push_back(coupling);
		flip_buffer.flips.push_back({ i,j });
		flip_buffer.new_sz.push_back({ 0, 0 });
	};

	std::complex<double> diag(const std::vector<int>& state) { 
		return { 0.0, 0.0 }; // coupling;
	}

	FlipList off_diag(const std::vector<int>& state) {		
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
	//have to add both sites, only implemented Si+Sj-
	int i, j;
	double S;

public:

	HeisenbergExchange(std::complex<double> coupling_in, int i_in, int j_in, double S_in)
		: Interaction(coupling_in), i(i_in), j(j_in), S(S_in) {
		flip_buffer.multipliers.push_back(coupling / 2.0);
		flip_buffer.multipliers.push_back(coupling / 2.0);
		flip_buffer.flips.push_back({ i, j });
		flip_buffer.flips.push_back({ i, j });
		flip_buffer.new_sz.push_back({ 0, 0 });
		flip_buffer.new_sz.push_back({ 0, 0 });
	};

	std::complex<double> diag(const std::vector<int>& state) {
		return 0.5 * coupling.real() * S * S * state[i] * state[j];
	}

	FlipList off_diag(const std::vector<int>& state) {
		flip_buffer.multipliers[0] = std::complex<double>{ 0.0, 0.0 };
		flip_buffer.multipliers[1] = std::complex<double>{ 0.0, 0.0 };
		if (state[i] != -S && state[j] != S) {
			flip_buffer.multipliers[0] = coupling / 2.0;
			flip_buffer.new_sz[0][0] = state[i] - 1;
			flip_buffer.new_sz[0][1] = state[j] + 1;
		}
		if (state[i] != S && state[j] != -S) {
			flip_buffer.multipliers[0] = coupling / 2.0;
			flip_buffer.new_sz[0][0] = state[i] + 1;
			flip_buffer.new_sz[0][1] = state[j] - 1;
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

	RingExchange(std::complex<double> coupling_in, RingList& r, int tri_label_in, bool CC_in) 
		: Interaction(coupling_in), tri_label(tri_label_in), CC(CC_in){
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
		flip_buffer.multipliers.push_back(coupling);
		flip_buffer.flips.push_back({ i,j,k });
		flip_buffer.new_sz.push_back({ 0, 0, 0 });
	}

	std::complex<double> diag(const std::vector<int>& state) { return { 0.0, 0.0 }; }

	FlipList off_diag(const std::vector<int>& state) {
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

class SpinModel {
	//List of interactions
	std::vector<Interaction*> interactions;

public:
	SpinModel() {};

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
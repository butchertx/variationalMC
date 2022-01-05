#include "Wavefunction.h"
#include "vmc_io.h"

JastrowFactor::JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_)
	: strength(strength_), neighbor_table(neighbor_table_) {}

JastrowFactor::JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_, std::vector<int>& configuration_)
	: strength(strength_), neighbor_table(neighbor_table_) {
	assert(neighbor_table_.size() == configuration_.size());
	initialize_table(configuration_);
}

void JastrowFactor::initialize_table(std::vector<int>& configuration_) {
	exp_table.clear();
	conf_sum = 0.0;
	std::vector<int> site_neighbors;
	for (int site = 0; site < configuration_.size(); ++site) {
		exp_table.push_back(0.0);
		site_neighbors = neighbor_table[site];
		for (auto neigh : site_neighbors) {
			exp_table[site] += configuration_[neigh];
		}
		conf_sum += 0.5 * exp_table[site] * configuration_[site];
	}
}

double JastrowFactor::greedy_eval(std::vector<int>& configuration_) {
	assert(neighbor_table.size() == configuration_.size());
	initialize_table(configuration_);
	double result = 0.0;
	for (int site = 0; site < configuration_.size(); ++site) {
		result += 0.5 * strength * configuration_[site] * exp_table[site];
	}
	return result;
}

double JastrowFactor::lazy_eval(std::vector<int>& flip_sites, std::vector<int>& new_val, std::vector<int>& configuration_) {
	//Note: this assumes only two sites have been changed - see M.W.Butcher thesis to generalize to any # of sites
	assert(flip_sites.size() == 2);
	double flip_sum = 0.0;
	double neighbor_sum = 0.0;
	int del_s = 0;
	//check if site2 is a neighbor of site1 for this Jastrow factor and add the associated term to the calculation
	//this can be sped up dramatically by pre-sorting the neighbor sites and using binary search instead
	//better way to do this might be a hard-coded matrix where the factors are just looked up
	if (std::find(neighbor_table[flip_sites[0]].begin(), neighbor_table[flip_sites[0]].end(), flip_sites[1]) != neighbor_table[flip_sites[0]].end()) {
		flip_sum = strength * (new_val[0] * new_val[1] - configuration_[flip_sites[0]] * configuration_[flip_sites[1]]
			- 0.5 * new_val[0] * configuration_[flip_sites[1]] - 0.5 * new_val[1] * configuration_[flip_sites[0]]);
	}
	for (int flip_idx = 0; flip_idx < flip_sites.size(); ++flip_idx) {
		del_s = (new_val[flip_idx] - configuration_[flip_sites[flip_idx]]);
		neighbor_sum += 0.5 * strength * exp_table[flip_sites[flip_idx]] * del_s;
	}
	return exp(flip_sum + neighbor_sum);
}

void JastrowFactor::update_tables(std::vector<int>& flip_sites, std::vector<int>& new_val, std::vector<int>& old_configuration_) {
	for (int flip_idx = 0; flip_idx < flip_sites.size(); ++flip_idx) {
		for (auto site : neighbor_table[flip_sites[flip_idx]]) {
			exp_table[site] += new_val[flip_idx] - old_configuration_[flip_sites[flip_idx]]; // update the spin in exp_table
			conf_sum += 0.5 * exp_table[site] * (new_val[flip_idx] - old_configuration_[flip_sites[flip_idx]]); // update conf_sum with the new exp_table value
		}
	}
}

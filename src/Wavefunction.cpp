#include "Wavefunction.h"
#include "vmc_io.h"

JastrowFactor::JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_)
	: strength(strength_), neighbor_table(neighbor_table_) {

	for (auto neighlist : neighbor_table_) {
		neighbor_set.push_back(std::set<int>(neighlist.begin(), neighlist.end()));
	}
}

JastrowFactor::JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_, std::vector<int>& configuration_)
	: strength(strength_), neighbor_table(neighbor_table_) {
	for (auto neighlist : neighbor_table_) {
		neighbor_set.push_back(std::set<int>(neighlist.begin(), neighlist.end()));
	}
	assert(neighbor_table.size() == configuration_.size());
	assert(neighbor_set.size() == configuration_.size());
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

double JastrowFactor::lazy_eval(std::vector<int>& flip_sites, std::vector<int>& new_val, std::vector<int>& configuration_, bool are_neighbors) {
	//Note: this assumes only two sites have been changed - see M.W.Butcher thesis to generalize to any # of sites
	assert(flip_sites.size() == 2);
	double flip_sum = 0.0;
	double neighbor_sum = 0.0;
	int del_s = 0;

	timer.flag_start_time("are neighbors");
	if (are_neighbors) {
		flip_sum = strength * (new_val[0] * new_val[1] - configuration_[flip_sites[0]] * configuration_[flip_sites[1]]
			- 0.5 * new_val[0] * configuration_[flip_sites[1]] - 0.5 * new_val[1] * configuration_[flip_sites[0]]);
	}
	timer.flag_end_time("are neighbors");
	timer.flag_start_time("neighbor sum");
	for (int flip_idx = 0; flip_idx < flip_sites.size(); ++flip_idx) {
		del_s = (new_val[flip_idx] - configuration_[flip_sites[flip_idx]]);
		neighbor_sum += 0.5 * strength * exp_table[flip_sites[flip_idx]] * del_s;
	}
	timer.flag_end_time("neighbor sum");
	//timer.flag_start_time("exp calc");
	return exp(flip_sum + neighbor_sum);
}

double JastrowFactor::lazy_eval(std::vector<int>& flip_sites, std::vector<int>& new_val, std::vector<int>& configuration_) {
	bool are_neighbors = neighbor_set[flip_sites[0]].find(flip_sites[1]) != neighbor_set[flip_sites[0]].end();
	double result = lazy_eval(flip_sites, new_val, configuration_, are_neighbors);
	//timer.flag_end_time("exp calc");
	return result;
}

void JastrowFactor::update_tables(std::vector<int>& flip_sites, std::vector<int>& new_val, std::vector<int>& old_configuration_) {
	for (int flip_idx = 0; flip_idx < flip_sites.size(); ++flip_idx) {
		for (auto site : neighbor_table[flip_sites[flip_idx]]) {
			exp_table[site] += new_val[flip_idx] - old_configuration_[flip_sites[flip_idx]]; // update the spin in exp_table
			conf_sum += 0.5 * exp_table[site] * (new_val[flip_idx] - old_configuration_[flip_sites[flip_idx]]); // update conf_sum with the new exp_table value
		}
	}
}

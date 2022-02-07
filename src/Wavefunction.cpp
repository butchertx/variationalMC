#include "Wavefunction.h"
#include "vmc_io.h"

JastrowDensity::JastrowDensity(double strength_) : strength(strength_){}

void JastrowDensity::initialize_Sz0(std::vector<int>& configuration_) {
	sz0_density = 0.0;
	for (int i = 0; i < configuration_.size(); ++i) {
		sz0_density += (1 - configuration_[i] * configuration_[i]);
	}
	//sz0_density /= configuration_.size();
}

double JastrowDensity::greedy_eval(std::vector<int>& configuration_) {
	int num_Sz0_greedy = 0;
	for (int i = 0; i < configuration_.size(); ++i) {
		num_Sz0_greedy += (1 - configuration_[i] * configuration_[i]);
	}
	return exp(strength * num_Sz0_greedy);// / configuration_.size());
}

double JastrowDensity::greedy_log_derivative(std::vector<int>& configuration_) {
	int num_Sz0_greedy = 0;
	for (int i = 0; i < configuration_.size(); ++i) {
		num_Sz0_greedy += (1 - configuration_[i] * configuration_[i]);
	}
	return num_Sz0_greedy;
}

double JastrowDensity::conf_sum_difference(std::vector<int>& flips, std::vector<int>& new_sz, std::vector<int>& configuration_) {
	
	double new_sum = 0.0;
	for (int s = 0; s < flips.size(); ++s) {
		new_sum += configuration_[flips[s]] * configuration_[flips[s]] - new_sz[s] * new_sz[s];
	}
	
	return new_sum; // / configuration_.size();
}

double JastrowDensity::lazy_eval(std::vector<int>& flips, std::vector<int>& new_sz, std::vector<int>& configuration_) {
	return exp(strength * (sz0_density + conf_sum_difference(flips, new_sz, configuration_)));
}

void JastrowDensity::update_Sz0(std::vector<int>& flips, std::vector<int>& new_sz, std::vector<int>& configuration_) {
	sz0_density += conf_sum_difference(flips, new_sz, configuration_);
}

JastrowFactor::JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_)
	: strength(strength_), neighbor_table(neighbor_table_) {

	neighbor_bool = std::vector<std::vector<bool>>(neighbor_table_.size(), std::vector<bool>(neighbor_table_.size(), false));
	for (int i = 0; i < neighbor_table_.size(); ++i) {
		for (int jdx = 0; jdx < neighbor_table_[i].size(); ++jdx) {
			neighbor_bool[i][neighbor_table_[i][jdx]] = true;
		}
	}

	for (auto neighlist : neighbor_table_) {
		neighbor_set.push_back(std::set<int>(neighlist.begin(), neighlist.end()));
	}
}

JastrowFactor::JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_, bool sz2_flag)
	: strength(strength_), neighbor_table(neighbor_table_), sz2(sz2_flag) {
	neighbor_bool = std::vector<std::vector<bool>>(neighbor_table_.size(), std::vector<bool>(neighbor_table_.size(), false));
	for (int i = 0; i < neighbor_table_.size(); ++i) {
		for (int jdx = 0; jdx < neighbor_table_[i].size(); ++jdx) {
			neighbor_bool[i][neighbor_table_[i][jdx]] = true;
		}
	}
	for (auto neighlist : neighbor_table_) {
		neighbor_set.push_back(std::set<int>(neighlist.begin(), neighlist.end()));
	}
}

JastrowFactor::JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_, std::vector<int>& configuration_)
	: strength(strength_), neighbor_table(neighbor_table_) {
	neighbor_bool = std::vector<std::vector<bool>>(neighbor_table_.size(), std::vector<bool>(neighbor_table_.size(), false));
	for (int i = 0; i < neighbor_table_.size(); ++i) {
		for (int jdx = 0; jdx < neighbor_table_[i].size(); ++jdx) {
			neighbor_bool[i][neighbor_table_[i][jdx]] = true;
		}
	}
	for (auto neighlist : neighbor_table_) {
		neighbor_set.push_back(std::set<int>(neighlist.begin(), neighlist.end()));
	}
	assert(neighbor_table.size() == configuration_.size());
	assert(neighbor_set.size() == configuration_.size());
	initialize_table(configuration_);
}

JastrowFactor::JastrowFactor(double strength_, std::vector<std::vector<int>> neighbor_table_, std::vector<int>& configuration_, bool sz2_flag)
	: strength(strength_), neighbor_table(neighbor_table_), sz2(sz2_flag) {
	neighbor_bool = std::vector<std::vector<bool>>(neighbor_table_.size(), std::vector<bool>(neighbor_table_.size(), false));
	for (int i = 0; i < neighbor_table_.size(); ++i) {
		for (int jdx = 0; jdx < neighbor_table_[i].size(); ++jdx) {
			neighbor_bool[i][neighbor_table_[i][jdx]] = true;
		}
	}
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
			exp_table[site] += sz2 ? configuration_[neigh]*configuration_[neigh] : configuration_[neigh];
		}
	}
	for (int site = 0; site < configuration_.size(); ++site) {
		conf_sum += 0.5 * exp_table[site] * (sz2 ? configuration_[site] * configuration_[site] : configuration_[site]);
	}
}

double JastrowFactor::greedy_eval(std::vector<int>& configuration_) {
	assert(neighbor_table.size() == configuration_.size());
	initialize_table(configuration_);
	double result = 0.0;
	for (int site = 0; site < configuration_.size(); ++site) {
		result += 0.5 * strength * (sz2 ? configuration_[site] * configuration_[site] : configuration_[site]) * exp_table[site];
	}
	return exp(result);
}

double JastrowFactor::greedy_log_derivative(std::vector<int>& configuration_) {
	assert(neighbor_table.size() == configuration_.size());
	std::vector<int> site_neighbors;
	double conf_sum_greedy = 0.0;
	std::vector<double> exp_table_greedy;
	for (int site = 0; site < configuration_.size(); ++site) {
		exp_table_greedy.push_back(0.0);
		site_neighbors = neighbor_table[site];
		for (auto neigh : site_neighbors) {
			exp_table_greedy[site] += sz2 ? configuration_[neigh] * configuration_[neigh] : configuration_[neigh];
		}
	}
	for (int site = 0; site < configuration_.size(); ++site) {
		conf_sum_greedy += 0.5 * exp_table_greedy[site] * (sz2 ? configuration_[site] * configuration_[site] : configuration_[site]);
	}
	return conf_sum_greedy;
}

double JastrowFactor::conf_sum_difference(std::vector<int>& flip_sites, std::vector<int>& new_val, std::vector<int>& configuration_, bool are_neighbors) {
	//See M.W.Butcher thesis to generalize to any # of sites
	assert(flip_sites.size() == new_val.size());
	double flip_sum = 0.0, neighbor_sum = 0.0;

	//compute Sz or Sz^2
	std::vector<int> new_confs(new_val.size()), old_confs(flip_sites.size()), del_s(flip_sites.size());
	for (int s = 0; s < flip_sites.size(); ++s) {
		if (sz2) {
			new_confs[s] = new_val[s] * new_val[s];
			old_confs[s] = configuration_[flip_sites[s]] * configuration_[flip_sites[s]];
			del_s[s] = new_confs[s] - old_confs[s];
		}
		else {
			new_confs[s] = new_val[s];
			old_confs[s] = configuration_[flip_sites[s]];
			del_s[s] = new_confs[s] - old_confs[s];
		}
	}
	if (flip_sites.size() == 2) {
		if (are_neighbors) {
			flip_sum = (new_confs[0] * new_confs[1] + old_confs[0] * old_confs[1]
				- new_confs[0] * old_confs[1] - new_confs[1] * old_confs[0]);
		}
		for (int flip_idx = 0; flip_idx < flip_sites.size(); ++flip_idx) {
			neighbor_sum += exp_table[flip_sites[flip_idx]] * del_s[flip_idx];
		}
	}
	else if (flip_sites.size() == 3) {
		//assume v12 = v13 = v23
		//also assume if are_neighbors is true (for sites 1 and 2) it is true for every pair of the three
		if (are_neighbors) {
			flip_sum = (old_confs[0]*old_confs[1] + old_confs[1]*old_confs[2] + old_confs[2]*old_confs[0]
				- old_confs[0]*old_confs[0] - old_confs[1]*old_confs[1] - old_confs[2]*old_confs[2]);
		}
		for (int flip_idx = 0; flip_idx < flip_sites.size(); ++flip_idx) {
			neighbor_sum += exp_table[flip_sites[flip_idx]] * del_s[flip_idx];
		}
	}
	else {
		throw std::runtime_error("Jastrow only implemented for 2- or 3- site updates");
	}
	return flip_sum + neighbor_sum;
}

double JastrowFactor::lazy_eval(std::vector<int>& flip_sites, std::vector<int>& new_val, std::vector<int>& configuration_) {
	bool are_neighbors = neighbor_bool[flip_sites[0]][flip_sites[1]]; 
	double result = conf_sum_difference(flip_sites, new_val, configuration_, are_neighbors);
	return exp(strength * result);
}

void JastrowFactor::update_tables(std::vector<int>& flip_sites, std::vector<int>& new_val, std::vector<int>& old_configuration_) {
	//compute Sz or Sz^2
	std::vector<int> new_confs(new_val.size()), old_confs(flip_sites.size());
	for (int s = 0; s < flip_sites.size(); ++s) {
		if (sz2) {
			new_confs[s] = new_val[s] * new_val[s];
			old_confs[s] = old_configuration_[flip_sites[s]] * old_configuration_[flip_sites[s]];
		}
		else {
			new_confs[s] = new_val[s];
			old_confs[s] = old_configuration_[flip_sites[s]];
		}
	}

	conf_sum += conf_sum_difference(flip_sites, new_val, old_configuration_, neighbor_bool[flip_sites[0]][flip_sites[1]]);

	for (int flip_idx = 0; flip_idx < flip_sites.size(); ++flip_idx) {
		for (auto site : neighbor_table[flip_sites[flip_idx]]) {
			exp_table[site] += new_confs[flip_idx] - old_confs[flip_idx]; // update the spin in exp_table
		}
	}

	
}

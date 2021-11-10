#pragma once
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>

class RandomEngine {

	static unsigned const TEST_SEED = 8271617;

	std::default_random_engine generator_default;

	std::mt19937 generator_twister;

	std::uniform_real_distribution<> distribution_prob;//Roll double between 0 and 1

	std::uniform_int_distribution<> distribution_spin;//Roll Sz value in [-1, 1]

	std::uniform_int_distribution<> distribution_lat;//Roll lattice site in [0, N]

	std::uniform_int_distribution<> distribution_neighbor; //Roll neighbor site in [0, N]

public:

	RandomEngine(int seed_in, int N_in, int swap_range_in) {
		unsigned seed;
		if (seed_in < 0) {
			seed = std::chrono::system_clock::now().time_since_epoch().count();
		}
		else if (seed_in == 0) {
			seed = TEST_SEED;
		}
		else {
			seed = seed_in;
		}

		generator_twister.seed(seed);
		generator_twister.discard(10000);//"Equilibrate" the generator
		generator_default.seed(seed);
		generator_default.discard(10000);//"Equilibrate" the generator

		distribution_prob = std::uniform_real_distribution<double>(0.0, 1.0);
		distribution_spin = std::uniform_int_distribution<int>(-1, 1);
		distribution_lat = std::uniform_int_distribution<int>(0, N_in - 1);
		distribution_neighbor = std::uniform_int_distribution<int>(0, swap_range_in - 1);
	}

	void set_lattice_distribution(int N) {
		distribution_lat = std::uniform_int_distribution<int>(0, N - 1);
	}

	double get_rand_prob() {
		return distribution_prob(generator_twister);
	}

	int get_rand_spin() {
		return distribution_spin(generator_twister);
	}

	int get_rand_site() {
		return distribution_lat(generator_twister);
	}

	int get_rand_neighbor() {
		return distribution_neighbor(generator_twister);
	}

	bool coin_flip() {
		return distribution_prob(generator_twister) > 0.5;
	}

	int get_rand_in_range(int max) {
		//return random int from 0 to max-1
		return std::floor(distribution_prob(generator_twister) * max);
	}

	std::vector<double> get_rand_vec_prob(int length) {
		std::vector<double> result(length);
		for (int i = 0; i < length; ++i) {
			result[i] = get_rand_prob();
		}
		return result;
	}

	std::vector<int> get_rand_vec_spin(int length) {
		std::vector<int> result(length);
		for (int i = 0; i < length; ++i) {
			result[i] = get_rand_spin();
		}
		return result;
	}

	std::vector<int> get_rand_vec_site(int length) {
		std::vector<int> result(length);
		for (int i = 0; i < length; ++i) {
			result[i] = get_rand_site();
		}
		return result;
	}

	std::vector<int> get_rand_vec_neighbor(int length) {
		std::vector<int> result(length);
		for (int i = 0; i < length; ++i) {
			result[i] = get_rand_neighbor();
		}
		return result;
	}

	std::vector<int> get_rand_spin_state(std::vector<int> spindist, int Nsite) {
		//spindist = {#(-s), #(-s+1), ... #(s-1), #(s)}
		std::vector<int> result(Nsite), sites;
		for (int i = 0; i < Nsite; ++i) {
			sites.push_back(i);
		}
		std::random_shuffle(sites.begin(), sites.end());

		int smax = spindist.size() % 2 == 0 ? spindist.size() / 2 : (spindist.size() - 1) / 2;
		int spinval = -smax;
		std::vector<int> spinvec;
		for (int sidx = 0; sidx < spindist.size(); ++sidx) {
			for (int i = 0; i < spindist[sidx]; ++i) {
				spinvec.push_back(spinval);
			}
			++spinval;
		}

		for (int i = 0; i < Nsite; ++i) {
			result[sites[i]] = spinvec[i];
		}

		return result;
	}
};
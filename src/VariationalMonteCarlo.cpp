#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !_USE_MATH_DEFINES
#include "mkl.h"
#include "VariationalMonteCarlo.h"
#include <exception>
#include "vmc_io.h"
#include "vmc_numerical.h"

std::vector<double> vec_r(std::vector<std::complex<double>> cvec) {
	std::vector<double> result;
	for (auto c : cvec) {
		result.push_back(c.real());
	}
	return result;
}

std::vector<double> vec_i(std::vector<std::complex<double>> cvec) {
	std::vector<double> result;
	for (auto c : cvec) {
		result.push_back(c.imag());
	}
	return result;
}


std::pair<std::vector<int>, std::vector<int>> MonteCarloEngine::step_su2() {

	//Propose move
	std::vector<int> sites(2, 0), spins(2, 0);
	auto wf_conf = WF.conf_ref();
	int tempspin = 0;

	//swap spins at two sites or change Nz-sector according to Bieri(2012) app. C
	while (sites[0] == sites[1]) {
		sites[0] = rand.get_rand_site();
		sites[1] = rand.get_rand_site();
	}
	spins[0] = wf_conf[sites[0]];
	spins[1] = wf_conf[sites[1]];

	if (spins[0] != spins[1]) {
		if (spins[0] != 0 && spins[1] != 0 && rand.get_rand_prob() < 0.5) {
			// |1>|-1> and coin flip
			//	flip to |0>|0>
			spins[0] = 0;
			spins[1] = 0;
		}
		else {
			// |0>|1>, |1>|0>, or (|1>|-1> and coin flip)
			// swap spins
			tempspin = spins[0];
			spins[0] = spins[1];
			spins[1] = tempspin;
		}
	}
	else if (spins[0] == 0){
		// |0>|0>
		// switch to |1>|-1>
		spins[0] = 1;
		spins[1] = -1;
	}
	else {
		// |1>|1> or |-1>|-1>
		//pick two different-spin sites and swap
		sites[1] = sites[0];
		spins[1] = spins[0];
		while (spins[0] == spins[1]) {
			sites[0] = rand.get_rand_site();
			spins[1] = wf_conf[sites[0]];
			sites[1] = rand.get_rand_site();
			spins[0] = wf_conf[sites[1]];
		}
	}

	double sqrt_p = std::abs(WF.psi_over_psi(sites, spins));
	if (rand.get_rand_prob() < sqrt_p * sqrt_p) {
		WF.update(sites, spins);
		return { sites, spins };
	}
	else {
		return { {}, {} };
	}



	
}

std::vector<int> MonteCarloEngine::step_su3(int num_site) {

	std::vector<int> swaplist(num_site);
	if (num_site == 3) {
		//Propose ring-exchange move
		int ring_label = rand.get_rand_in_range(lat.ring_ref().get_size());
		bool cc = rand.coin_flip();
		swaplist = lat.ring_ref().get_ring(ring_label);
		int tempsite;
		if (cc) {
			tempsite = swaplist[0];
			swaplist[0] = swaplist[1];
			swaplist[1] = tempsite;
		}
	}
	else if (num_site == 2) {
		//Propose swap move
		int site = rand.get_rand_site();
		int neigh = lat.get_neighbors(site, 0)[rand.get_rand_neighbor()];
		swaplist = { site, neigh };
	}

	std::complex<double> sqrt_p = WF.psi_over_psi(swaplist);
	if (rand.get_rand_prob() < std::abs(sqrt_p) * std::abs(sqrt_p)) {
		WF.update(swaplist);
		return swaplist;
	}
	else {
		return {};
	}
}

void MonteCarloEngine::measure_energy() {

	std::complex<double> energy = (0.0, 0.0), temp_O_val = (0.0, 0.0),  temp_E_val = (0.0, 0.0);

	FlipList f;

	//Run over all interactions in model
	for (std::string term_name : H.get_terms()) {
		temp_O_val = std::complex<double>(0.0, 0.0);
		auto interaction_list = H.get_interactions(term_name);
		for (auto interaction = interaction_list.begin(); interaction != interaction_list.end(); ++interaction) {

			timer.flag_start_time("Energy diag");
			temp_E_val = (*interaction)->diag(WF.conf_ref());
			temp_O_val += temp_E_val;
			energy += H.get_coupling(term_name) * temp_E_val;
			timer.flag_end_time("Energy diag");

			//get flips
			f = (*interaction)->off_diag(WF.conf_ref());

			timer.flag_start_time("Energy fliplist iteration");
			//iterate through flips
			for (int i = 0; i < f.multipliers.size(); ++i) {
				temp_E_val = f.multipliers[i] * WF.psi_over_psi(f.flips[i], f.new_sz[i]);
				temp_O_val += temp_E_val;
				energy += H.get_coupling(term_name) * temp_E_val;
			}
			timer.flag_end_time("Energy fliplist iteration");
		}
		observable_measures[term_name].push_back(temp_O_val);
	}
	//if (std::abs(energy.imag()) > 1e-8) {
	//	std::cout << "Imaginary Energy: " << std::abs(energy.imag()) << "\n";
	//	assert(std::abs(energy.imag()) < 1e-8);
	//}

	//std::cout << energy.real() / (2 * lat.get_N()) << "\n";
	E.push_back((energy.real() + H.get_E0().real()));
	E_imag.push_back((energy.imag() + H.get_E0().imag()));

	//if (!params.su3) {
	//	double nz_val = 0.0;
	//	auto conf = WF.conf_ref();
	//	for (int i = 0; i < conf.size(); ++i) {
	//		nz_val += (1 - conf[i]*conf[i]);
	//	}
	//	Nz.push_back(nz_val);
	//}
}

void MonteCarloEngine::measure_obs_functions() {

	std::complex<double> temp_O_val = (0.0, 0.0);
	std::vector<Observable> temp_function;
	std::vector<std::complex<double>> temp_function_eval;

	FlipList f;

	//Run over all observable functions
	timer.flag_start_time("Obs function calculation");
	for (std::string term_name : observable_function_names) {

		//Each observable function is a vector of Observables
		temp_function = observable_functions[term_name];
		temp_function_eval = std::vector<std::complex<double>>();

		for (int i = 0; i < temp_function.size(); ++i) {

			//Each Observable is itself a sum of Interactions
			temp_O_val = std::complex<double>(0.0, 0.0);
			auto interaction_list = temp_function[i].get();

			for (auto interaction = interaction_list.begin(); interaction != interaction_list.end(); ++interaction) {

				temp_O_val += (*interaction)->diag(WF.conf_ref());

				//get flips
				f = (*interaction)->off_diag(WF.conf_ref());

				//iterate through flips
				for (int i = 0; i < f.multipliers.size(); ++i) {
					temp_O_val += f.multipliers[i] * WF.psi_over_psi(f.flips[i], f.new_sz[i]);
				}

			}

			temp_function_eval.push_back(temp_O_val);
			
		}

		observable_function_measures[term_name].push_back(temp_function_eval);

	}
	timer.flag_end_time("Obs function calculation");

}

//find the error via the bin technique using a specified number of bins
double MonteCarloEngine::error(std::vector<double>::iterator vals_begin, std::vector<double>::iterator vals_end, double mean, int bins) {
	if (bins < 2) {
		return 0.0;
	}
	else {
		int NMC = bins;//Nb is bin size, NMC is number of bins
		int Nb = (vals_end - vals_begin) / NMC;
		std::vector<double> avgs(NMC);
		for (int i = 0; i < NMC; ++i) {
			double avg = 0;
			for (int j = 0; j < Nb; ++j) {
				avg += *(vals_begin + (Nb*i + j));
			}
			avg = avg / Nb;
			avgs[i] = avg;
		}
		double std_dev = 0;
		for (int i = 0; i < avgs.size(); ++i) {
			std_dev += (avgs[i] - mean)*(avgs[i] - mean);
		}
		std_dev = sqrt(std_dev / NMC / (NMC - 1));
		return std_dev;
	}
}

void MonteCarloEngine::run() {
	if (!params.optimization) {
		run_bin(false);
	}
	else {
		// run SR iterations
		// this would be much clearer with some simple linear algebra implementation (MKL?)
		std::vector<std::vector<double>> ok_measures;
		std::vector<double> E_bin;

		std::vector<std::vector<double>> sjk; // vparam covariation
		double e_mean;
		std::vector<double> ok_mean, E_ok_mean;
		std::vector<double> fk; // cross correlation with H

		std::vector<double> alpha;

		std::ofstream optfile;
		optfile.open("results/optimization.csv");
		optfile << "SR Bin, Total Energy, V params\n";

		for (int SR_bin = 0; SR_bin < params.sr.bins; ++SR_bin) {

			ok_measures.clear();
			E_bin.clear();

			// run bin and collect O_k measurements
			ok_measures = run_bin(true);
			for (int m = 0; m < ok_measures.size(); ++m) {
				ok_measures[m].insert(ok_measures[m].begin(), 1.0);
			}
			std::copy(E.end() - ok_measures.size(), E.end(), std::back_inserter(E_bin));
			ok_mean = std::vector<double>(ok_measures[0].size(), 0.0);
			E_ok_mean = std::vector<double>(ok_measures[0].size(), 0.0);
			e_mean = 0.0;
			for (int v = 0; v < ok_measures[0].size(); ++v) {
				E_ok_mean[v] = 0.0;
				for (int m = 0; m < ok_measures.size(); ++m) {
					if (v == 0) {
						e_mean += E_bin[m] / E_bin.size();
					}
					ok_mean[v] += ok_measures[m][v] / E_bin.size();
					E_ok_mean[v] += E_bin[m] * ok_measures[m][v] / E_bin.size();
				}
			}


			// calculate forces and cross-correlations
			// this is a naive implementation of the cross-correlation.  It is hopefully not too slow unless there are many variational parameters
			sjk = std::vector<std::vector<double>>(ok_measures[0].size(), std::vector<double>(ok_measures[0].size()));
			fk = std::vector<double>(ok_measures[0].size(), 0.0);
			for (int j = 0; j < sjk.size(); ++j) {
				for (int k = 0; k < sjk.size(); ++k) {
					sjk[j][k] = 0.0;
					for (int m = 0; m < ok_measures.size(); ++m) {
						sjk[j][k] += ok_measures[m][j] * ok_measures[m][k];
					}
					sjk[j][k] /= ok_measures.size();
				}
			}
			for (int j = 0; j < sjk.size(); ++j) {
				for (int k = 0; k < sjk.size(); ++k) {
					sjk[j][k] -= sjk[j][0] * sjk[k][0];
				}
			}

			for (int k = 0; k < fk.size(); ++k) {
				fk[k] = ok_mean[k] * e_mean - E_ok_mean[k];
			}




			// invert and solve for d\alpha
			std::vector<double> d_alpha = solve_linear_system(sjk, fk);

			// save values of old params and energy bins
			alpha = WF.get_parameters();
			v_params_record.push_back(alpha);
			E_bins_record.push_back(e_mean);

			// send new params to WF
			for (int k = 0; k < alpha.size(); ++k) {
				alpha[k] += params.sr.timestep * d_alpha[k+1];
			}
			WF.update_parameters(alpha);

			optfile << SR_bin+1 << "," << e_mean << "," << vec2str(alpha) << "\n";
			std::cout << "SR iteration " << SR_bin + 1 << " of " << params.sr.bins << "\n";
			std::cout << "Bin energy = " << e_mean << "; New params: " << vec2str(alpha) << "\n";
		}
		optfile.close();
	}
}

//void  FSSMonteCarlo::step() {
//	//Propose new triangle configuration
//	int ring_label = rand.get_rand_in_range(lat.ring_ref().get_size());
//	Color rand_perm1[6] = { R, R, G, G, B, B };
//	Color rand_perm2[6] = { G, B, R, B, G, R };
//	Color rand_perm3[6] = { B, G, B, R, R, G };
//	int roll = 0;
//	std::vector<int> ring = lat.ring_ref().get_ring(ring_label);
//	ThreeBar s1 = conf.read(ring[0]), s2 = conf.read(ring[1]), s3 = conf.read(ring[2]),
//		s1prop = s1, s2prop = s2, s3prop = s3;
//	if (ring_label % 2 == 0) {
//		while (s1prop.x == s1.x && s2prop.x == s2.x && s3prop.x == s3.x) {
//			//roll new configuration
//			roll = rand.get_rand_in_range(6);
//			s1prop.x = rand_perm1[roll];
//			s2prop.x = rand_perm2[roll];
//			s3prop.x = rand_perm3[roll];
//		}
//	}
//	else {
//		while (s1prop.y == s1.y && s2prop.y == s2.y && s3prop.y == s3.y) {
//			//roll new configuration
//			roll = rand.get_rand_in_range(6);
//			s1prop.y = rand_perm1[roll];
//			s2prop.y = rand_perm2[roll];
//			s3prop.y = rand_perm3[roll];
//		}
//	}
//
//	//Check probability and accept/reject move
//	bool accept = WF.check_basis_element(s1prop, s2prop, s3prop);
//	if (accept) {
//		WF.update({ s1prop, s2prop, s3prop }, ring_label);
//	}
//}
//
//void FSSMonteCarlo::measure_energy() {
//	std::complex<double> energy = (0.0, 0.0), econt = (0.0, 0.0);
//	double K;
//
//	FlipList f;
//
//	//Run over all interactions in model
//	for (auto interaction = H.get().begin(); interaction != H.get().end(); ++interaction) {
//		econt = WF.psi_over_psi((*interaction)->get_tri_label(), (*interaction)->get_CC());
//		K = (*interaction)->get_coupling();
//		energy += (*interaction)->get_CC() ?
//			std::complex<double>(K * econt.real() - h * econt.imag(), K * econt.imag() + h*econt.real()) :
//			std::complex<double>(K * econt.real() + h * econt.imag(), K * econt.imag() - h*econt.real());
//	}
//
//	//std::cout << energy.real() / (2 * lat.get_N()) << "\n";
//	energy = std::complex<double>{ energy.real() / lat.get_N(), energy.imag() / lat.get_N() };
//	E.push_back(energy);
//}
//
//double FSSMonteCarlo::error(std::vector<double>::iterator vals_begin, std::vector<double>::iterator vals_end, double mean, int bins) {
//	if (bins < 2) {
//		return 0.0;
//	}
//	else {
//		int NMC = bins;//Nb is bin size, NMC is number of bins
//		int Nb = (vals_end - vals_begin) / NMC;
//		std::vector<double> avgs(NMC);
//		for (int i = 0; i < NMC; ++i) {
//			double avg = 0;
//			for (int j = 0; j < Nb; ++j) {
//				avg += *(vals_begin + (Nb*i + j));
//			}
//			avg = avg / Nb;
//			avgs[i] = avg;
//		}
//		double std_dev = 0;
//		for (int i = 0; i < avgs.size(); ++i) {
//			std_dev += (avgs[i] - mean)*(avgs[i] - mean);
//		}
//		std_dev = sqrt(std_dev / NMC / (NMC - 1));
//		return std_dev;
//	}
//}
//
//std::complex<double> FSSMonteCarlo::error(std::vector<std::complex<double>>::iterator start, std::vector<std::complex<double>>::iterator end, std::complex<double> mean, int bins) {
//	std::vector<double> re, imag;
//	double mean_r = mean.real(), mean_i = mean.imag();
//	for (auto it = start; it != end; ++it) {
//		re.push_back(it->real());
//		imag.push_back(it->imag());
//	}
//	return std::complex<double>{error(re.begin(), re.end(), mean_r, bins), error(imag.begin(), imag.end(), mean_i, bins) };
//}

//double VariationalMonteCarlo::bootstrap(std::vector<double> x, int n_boot, std::string func) {
//	//compute the error in <x^2> - <x>^2 using the bootstrap method
//	//use n_boot as the number of bootstrap bins and a random seed
//	std::vector<double> boots(n_boot);
//	rand_init_(&n_boot);
//	double xB, xB_sq, x_temp;
//	int N = x.size();
//	for (int boot = 0; boot < n_boot; ++boot) {
//		//find a bootstrap value for x^B_boot and x^2 ^B_boot
//		xB = 0;
//		xB_sq = 0;
//		for (int i = 0; i < N; ++i) {
//			x_temp = x[(int)(drand1_()*N)];
//			xB += x_temp;
//			xB_sq += x_temp * x_temp;
//		}
//		xB /= N;
//		xB_sq /= N;
//		//find the values for f^B_boot
//		if (func.compare("binder") == 0) {
//			//for this case, input x should be m^2 measurements
//			boots[boot] = xB_sq / (xB*xB);
//		}
//		else if (func.compare("susceptibility") == 0) {
//			boots[boot] = xB_sq - xB * xB;
//		}
//		else {
//			std::cout << "Error: invalid bootstrap function option \"" << func << "\". Returning 0 for error estimate.\n";
//			return 0;
//		}
//	}
//	//average the values for f^B_boot to get susceptibility
//	double av = mean(boots);
//	//std::cout << "Bootstrap mean: " << av << "\n";
//
//	//std dev for f^B is std dev for f
//	double boots_sq = 0;
//	for (int i = 0; i < n_boot; ++i) {
//		boots_sq += boots[i] * boots[i];
//	}
//	boots_sq /= n_boot;
//	return sqrt(boots_sq - av * av);
//}
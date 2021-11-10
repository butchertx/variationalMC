#include "VariationalMonteCarlo.h"
#include <exception>

//void MonteCarloEngine::step_su2(std::string step_type) {
//	
//	//Propose move
//	std::vector<int> sites(2, 0), spins(2, 0);
//	int tempspin = 0;
//	if (strcmp(step_type.c_str(), "NZ_ADAPT") == 0) {
//		//swap spins at two sites or change Nz-sector according to Bieri(2012) app. C
//		while (sites[0] == sites[1]) {
//			sites[0] = rand.get_rand_site();
//			sites[1] = rand.get_rand_site();
//		}
//		std::cout << "Remember to change so that spins[i] holds the NEW spin at site [i]\n";
//		spins[0] = conf.read_spin(sites[0]);
//		spins[1] = conf.read_spin(sites[1]);
//		std::cout << "Spins are : (" << spins[0] << "," << spins[1] << ")\n";
//
//		if (spins[0] != spins[1] || spins[0] == 0) {
//			if (spins[0] == spins[1]) {
//				spins[0] = 1;
//				spins[1] = -1;
//				std::cout << "Flip up and down\n";
//			}
//			else if (spins[0] != 0 && spins[1] != 0 && rand.get_rand_prob() < 0.5) {
//				spins[0] = 0;
//				spins[1] = 0;
//				std::cout << "Flip to 0\n";
//			}
//			else {
//				tempspin = spins[0];
//				spins[0] = spins[1];
//				spins[1] = tempspin;
//				std::cout << "Swap\n";
//			}
//		}
//		else {
//			//pick two different-spin sites and swap
//			sites[1] = sites[0];
//			spins[1] = spins[0];
//			while (sites[0] == sites[1] || spins[0] == spins[1]) {
//				sites[0] = rand.get_rand_site();
//				spins[1] = conf.read_spin(sites[0]);
//				sites[1] = rand.get_rand_site();
//				spins[0] = conf.read_spin(sites[1]);
//			}
//			std::cout << "Swapping two different sites\n";
//		}
//
//		double sqrt_p = std::abs(WF.psi_over_psi(conf.ref(), sites));//, spins));
//		//std::cout << "sqrt_p = " << sqrt_p << "\n";
//		if (rand.get_rand_prob() < sqrt_p * sqrt_p) {
//			conf.set(sites[0], sites[1], spins[0], spins[1]);
//			//std::cout << "Move Accepted, Nz = " << conf.get_Nz() << "\n";
//		}
//		
//	}
//	else if (strcmp(step_type.c_str(), "NZ_CONST") == 0) {
//		while (sites[0] == sites[1] || spins[0] == spins[1]) {
//			sites[0] = rand.get_rand_site();
//			sites[1] = lat.get_neighbor(sites[0], rand.get_rand_neighbor());
//			spins[0] = conf.read_spin(sites[1]);
//			spins[1] = conf.read_spin(sites[0]);
//		}
//		
//		//std::cout << "Spins are : (" << spins[0] << "," << spins[1] << ")\n";
//		//Check probability and accept/reject move
//		std::complex<double> sqrt_p = WF.psi_over_psi(sites, spins);
//		//if ((double)std::abs(std::abs(sqrt_p) - std::abs(WF.psi_over_psi_full_calculation(sites[0], sites[1], spins[0], spins[1]))) > 1e-10) {
//			//std::cout << "Sqrt_p and Full Det: " << sqrt_p << "; " << WF.psi_over_psi_full_calculation(sites[0], sites[1], spins[0], spins[1]) << "\n";
//		//}
//		//std::cout << "sqrt_p = " << sqrt_p << "\n";
//		//std::cout << "psi/psi full calculation = " << WF.psi_over_psi_full_calculation(sites[0], sites[1], spins[0], spins[1]) << "\n";
//		if (rand.get_rand_prob() < std::abs(sqrt_p) * std::abs(sqrt_p)) {
//			WF.update(sites[0], sites[1], sqrt_p);
//			//std::cout << "Move Accepted, Nz = " << conf.get_Nz() << "\n";
//		}
//	}
//
//
//
//	
//}

//void MonteCarloEngine::propose_ring_exchange(int tri_label, bool CC, std::vector<int>& affected_sites, std::vector<int>& affected_rings) {
//
//	std::vector<int> sites = lat.ring_ref().get_ring(tri_label);
//	//also need the sites and labels of neighboring rings
//	std::vector<int> ring_labels1 = lat.ring_ref().get_ring_labels(sites[0]),
//		ring_labels2 = lat.ring_ref().get_ring_labels(sites[1]),
//		ring_labels3 = lat.ring_ref().get_ring_labels(sites[2]);
//
//
//	affected_rings[0] = tri_label;
//	affected_sites[0] = sites[0];
//	if (CC) {
//		affected_sites[1] = sites[1];
//		affected_sites[2] = sites[2];
//	}
//	else {
//		affected_sites[2] = sites[1];
//		affected_sites[1] = sites[2];
//	}
//
//	if (ring_labels1[0] == tri_label) {
//		affected_rings[1] = ring_labels1[1];
//	}
//	else {
//		affected_rings[1] = ring_labels1[0];
//	}
//	sites = lat.ring_ref().get_ring(affected_rings[1]);
//
//	for (int i = 0; i < sites.size(); ++i) {
//		//if site is in original ring, replace with new site after ring exchange
//		if (sites[i] == affected_sites[0]) {
//			sites[i] = affected_sites[2];
//		}
//	}
//
//	affected_sites[3] = sites[0];
//	affected_sites[4] = sites[1];
//	affected_sites[5] = sites[2];
//
//
//	if (ring_labels2[0] == tri_label) {
//		affected_rings[2] = ring_labels2[1];
//	}
//	else {
//		affected_rings[2] = ring_labels2[0];
//	}
//	sites = lat.ring_ref().get_ring(affected_rings[2]);
//	for (int i = 0; i < sites.size(); ++i) {
//		//if site is in original ring, replace with new site after ring exchange
//		if (sites[i] == affected_sites[1]) {
//			sites[i] = affected_sites[0];
//		}
//	}
//
//	affected_sites[6] = sites[0];
//	affected_sites[7] = sites[1];
//	affected_sites[8] = sites[2];
//
//	if (ring_labels3[0] == tri_label) {
//		affected_rings[3] = ring_labels3[1];
//	}
//	else {
//		affected_rings[3] = ring_labels3[0];
//	}
//	sites = lat.ring_ref().get_ring(affected_rings[3]);
//	for (int i = 0; i < sites.size(); ++i) {
//		//if site is in original ring, replace with new site after ring exchange
//		if (sites[i] == affected_sites[2]) {
//			sites[i] = affected_sites[1];
//		}
//	}
//
//	affected_sites[9] = sites[0];
//	affected_sites[10] = sites[1];
//	affected_sites[11] = sites[2];
//
//}

void MonteCarloEngine::step_su3(int num_site) {

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
	}
}

//void MonteCarloEngine::measure_energy_su2() {
//
//	std::complex<double> energy = (0.0, 0.0);
//
//	FlipList f;
//
//	//Run over all interactions in model
//	for (auto interaction = H.get().begin(); interaction != H.get().end(); ++interaction) {
//
//		//get diagonal contribution and add to energy
//		energy.real(energy.real() + (*interaction)->diag(conf.ref()));
//
//		//get flips
//		f = (*interaction)->off_diag(conf.ref());
//
//		//iterate through flips
//		for (int i = 0; i < f.multipliers.size(); ++i){
//			energy += f.multipliers[i] * WF.psi_over_psi(f.flips[i], f.new_sz[i]);
//		}
//		//std::cout << "Spot check: Interaction output: \n";
//		//(*interaction)->print_info(conf.ref());
//		//std::cout << "Energy result of first flip = " << f.multipliers[0] * WF.psi_over_psi(conf.ref(), f.flips[0], f.new_sz[0]) << "\n";
//
//	}
//	if (std::abs(energy.imag()) > 1e-8) {
//		std::cout << "Imaginary Energy: " << std::abs(energy.imag()) << "\n";
//		assert(std::abs(energy.imag()) < 1e-8);
//	}
//	
//	//std::cout << energy.real() / (2 * lat.get_N()) << "\n";
//	E.push_back(energy.real() / (lat.get_N()));
//}

//void MonteCarloEngine::measure_energy_su3() {
//
//	std::complex<double> energy = (0.0, 0.0);
//
//	FlipList f;
//	std::vector<int> affected_rings(4);
//	std::vector<int> affected_sites(12);
//
//	//Run over all interactions in model
//	for (auto interaction = H.get().begin(); interaction != H.get().end(); ++interaction) {
//
//		//set ring labels
//		//propose_ring_exchange((*interaction)->get_tri_label(), (*interaction)->get_CC(), affected_sites, affected_rings);
//
//		//get diagonal contribution and add to energy
//		energy.real(energy.real() + (*interaction)->diag(conf.ref()));
//
//		//get flips
//		f = (*interaction)->off_diag(conf.ref());
//
//		//iterate through flips
//		for (int i = 0; i < f.multipliers.size(); ++i) {
//			energy += f.multipliers[i] * WF.psi_over_psi((*interaction)->get_tri_label(), (*interaction)->get_CC());
//		}
//		//std::cout << "Spot check: Interaction output: \n";
//		//(*interaction)->print_info(conf.ref());
//		//std::cout << "Energy result of first flip = " << f.multipliers[0] * WF.psi_over_psi(conf.ref(), f.flips[0], f.new_sz[0]) << "\n";
//
//	}
//	if (std::abs(energy.imag()) > 1e-8) {
//		std::cout << "Imaginary Energy: " << std::abs(energy.imag()) << "\n";
//		assert(std::abs(energy.imag()) < 1e-8);
//	}
//
//	//std::cout << energy.real() / (2 * lat.get_N()) << "\n";
//	E.push_back(energy.real() / (lat.get_N()));
//	E_imag.push_back(energy.imag() / lat.get_N());
//}

void MonteCarloEngine::measure_energy() {

	std::complex<double> energy = (0.0, 0.0);

	FlipList f;
	std::vector<int> affected_rings(4);
	std::vector<int> affected_sites(12);

	//Run over all interactions in model
	for (auto interaction = H.get().begin(); interaction != H.get().end(); ++interaction) {

		energy += (*interaction)->diag(WF.conf_ref());

		//get flips
		f = (*interaction)->off_diag(WF.conf_ref());

		//iterate through flips
		for (int i = 0; i < f.multipliers.size(); ++i) {
			energy += f.multipliers[i] * WF.psi_over_psi(f.flips[i], f.new_sz[i]);
		}
		//std::cout << "Spot check: Interaction output: \n";
		//(*interaction)->print_info(conf.ref());
		//std::cout << "Energy result of first flip = " << f.multipliers[0] * WF.psi_over_psi(conf.ref(), f.flips[0], f.new_sz[0]) << "\n";

	}
	//if (std::abs(energy.imag()) > 1e-8) {
	//	std::cout << "Imaginary Energy: " << std::abs(energy.imag()) << "\n";
	//	assert(std::abs(energy.imag()) < 1e-8);
	//}

	//std::cout << energy.real() / (2 * lat.get_N()) << "\n";
	E.push_back(energy.real() / (lat.get_N()));
	E_imag.push_back(energy.imag() / lat.get_N());
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
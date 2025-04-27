#include "ProjectedState.h"
#include "mkl_types.h"

// printing

void ProjectedState::print_matrix(std::string name){
	if (std::strcmp(name.c_str(), "Slater") == 0){
		vmc_io::print_matrix("Slater", N, N, Slater, N);
	}
	else if (std::strcmp(name.c_str(), "LU") == 0){
		vmc_io::print_matrix("LU", N, N, LU, N);
	}
	else if (std::strcmp(name.c_str(), "Winv") == 0){
		vmc_io::print_matrix("Winv", DIM, N, Winv, N);
	}
	else if (std::strcmp(name.c_str(), "UP1") == 0){
		vmc_io::print_matrix("UP1", DIM, 2, UP1, 2);
	}
	else if (std::strcmp(name.c_str(), "UP2") == 0){
		vmc_io::print_matrix("UP2", 2, N, UP2, N);
	}
	else if (std::strcmp(name.c_str(), "UP3") == 0){
		vmc_io::print_matrix("UP3", 2, N, UP3, N);
	}
	else if (std::strcmp(name.c_str(), "ipiv") == 0){
		vmc_io::print_matrix("ipiv", N, N, ipiv, N);
	}
	else if (std::strcmp(name.c_str(), "Phi") == 0){
		vmc_io::print_matrix("Phi", DIM, DIM, ansatz.get_Phi(), DIM);
	}
	else {
		std::cerr << "Matrix name " << name << " not recognized.\n";
	}
}

// constructors

ProjectedState::ProjectedState(MeanFieldAnsatz& M_, RandomEngine& rand_)
	: ansatz(M_), rand(rand_), N(ansatz.get_N()), DIM(ansatz.get_dim()){
	// initialize the Slater-Jastrow state

	conserve_sz2 = M_.get_conserve_sz2();
	malloc_matrices();
	clear_matrices();
	initialize_configuration();
}

ProjectedState::ProjectedState(MeanFieldAnsatz& M_, RandomEngine& rand_, JastrowTable jastrow_)
	: ProjectedState(M_, rand_){
	conserve_sz2 = M_.get_conserve_sz2();
	jastrow = jastrow_;
	jastrow.initialize_tables(configuration);
}

void ProjectedState::clear_matrices(){
	for (int i = 0; i < N * N; ++i) {
		Slater[i] = { 0,0 };
		LU[i] = { 0,0 };
		ipiv[i] = 0;
	}
	for (int i = 0; i < DIM * N; ++i) {
		Winv[i] = { 0,0 };
	}
	for (int i = 0; i < 2 * DIM; ++i) {
		UP1[i] = { 0.0, 0.0 };
		if (i < 2 * N) {
			UP2[i] = { 0.0, 0.0 };
			UP3[i] = { 0.0, 0.0 };
		}
	}
}

void ProjectedState::initialize_configuration(){
	int config_attempt = 0;
	while (!try_configuration() && config_attempt < CONFIG_ATTEMPTS) {
		det = { 0, 0 };
		++config_attempt;
	}
}

bool ProjectedState::try_configuration() {
	// TODO: do this in a more robust manner, maybe type checking on initialization
	int N0 = ansatz.get_N0F();
	if (2 * ((N - N0) / 2) != N - N0) {
		N0 += 1;
	}
	set_configuration(rand.get_rand_spin_state(std::vector<int>{ (N - N0) / 2, N0, (N - N0) / 2 }, N));
	
	CBLAS_INDEX low = 0;
	low = cblas_izamin(N, LU, N+1);
	return (cblas_dcabs1(&(LU[low*(N+1)])) > 10e-10);
}

// TODO: this implementation can be sped up with a lookup table
int ProjectedState::Spin_t_to_row(int spin_idx){
	if (ansatz.get_spin_type() == vmctype::Spin_t::ONE){
		// spin_idx = -1, 0, or 1
		return (-spin_idx + 1) * N;
	}
	else if (ansatz.get_spin_type() == vmctype::Spin_t::HALF){
		// spin_idx = -1, or 1
		return ((-spin_idx + 1) / 2) * N;
	}
	throw vmctype::NotImplemented("Spins other than 1/2 and 1 not implemented.");
}

void ProjectedState::set_configuration(std::vector<int> conf) {
	configuration = conf;
	int row = 0;
	lapack_complex_double* phi = ansatz.get_Phi();
	parton_labels.clear();
	int info;
	MKL_Complex16 alpha = { 1.0, 0.0 }, beta = { 0.0, 0.0 };
	MKL_INT64 N_64 = N, DIM_64 = DIM; // needed for use with intel ilp64 interface / libraries

	// move the relevant rows of Phi into Slater
	for (int i = 0; i < N; ++i) {
		parton_labels.push_back(i);
		row = Spin_t_to_row(configuration[i]) + i;
		std::memcpy(&(Slater[parton_labels[i] * N]), &(phi[row * DIM]), N * sizeof(lapack_complex_double));
	}
	info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, Slater, N, ipiv);
	std::memcpy(LU, Slater, N * N * sizeof(lapack_complex_double));
	if (info == 0) {
		info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, N, Slater, N, ipiv);
		// zgemm3m("N", "N", &DIM, &N, &N, &alpha, phi, &N, Slater, &N, &beta, Winv, &N);
		cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, DIM, N, N, &alpha, phi, DIM, Slater, N, &beta, Winv, N);
		// cblas_zgemm3m_64(CblasRowMajor, CblasNoTrans, CblasNoTrans, DIM_64, N_64, N_64, &alpha, phi, DIM_64, Slater, N_64, &beta, Winv, N_64);
	}
	det = calc_det();
}

std::complex<double> ProjectedState::calc_det() {
	std::complex<double> result = { 1.0, 0.0 };
	for (int i = 0; i < N; ++i) {
		result *= LU[i*N + i];
	}
	return result;
}

void ProjectedState::upinvhop2(int rowk, int colk, int rowl, int coll) {
	//perform the update of Winv according to the Woodbury Matrix identity
	//Winv' = Winv - Winv * U (I_k + V A^-1 U)^-1 V A^-1
	//where U and V are defined so A' = A + UV (A is the Slater matrix)
	//UP1 is Winv * U; the i and j columns of Winv where i and j are the sites to update
	//UP2 is V A^-1 (which can be computed easily from rows of Winv)
	//UP3 is (I_k + V A^-1 U)^-1 (2x2) times UP2 (2xN)

	std::complex<double>
		c11 = Winv[rowl * N + coll],
		c22 = Winv[rowk * N + colk],
		c12 = -Winv[rowk * N + coll],
		c21 = -Winv[rowl * N + colk];

	std::complex<double> g = c11 * c22 - c12 * c21, beta = { 1.0, 0.0 };

	cblas_zcopy(DIM, &(Winv[colk]), N, UP1, 2);
	cblas_zcopy(DIM, &(Winv[coll]), N, &(UP1[1]), 2);
	cblas_zcopy(N, &(Winv[rowk * N]), 1, UP2, 1);
	UP2[colk] -= std::complex<double>(1.0, 0.0);
	cblas_zcopy(N, &(Winv[rowl * N]), 1, &(UP2[N]), 1);
	UP2[N + coll] -= std::complex<double>(1.0, 0.0);

	g = std::complex<double>(-1.0, 0.0) / g;

	for (int i = 0; i < N; ++i) {
		UP3[i] = c11 * UP2[i] + c12 * UP2[N + i];
		UP3[N + i] = c21 * UP2[i] + c22 * UP2[N + i];
	}

	MKL_INT64 N_64 = N, DIM_64 = DIM; // needed for use with intel ilp64 interface / libraries
	cblas_zgemm3m_64(CblasRowMajor, CblasNoTrans, CblasNoTrans, DIM_64, N_64, 2, &g, UP1, 2, UP3, N_64, &beta, Winv, N_64);
	// cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, DIM, N, 2, &g, UP1, 2, UP3, N, &beta, Winv, N);
}

/// PRIVATE MATRIX ELEMENTS

// swap 2 sites with specified sz values, with jastrow
std::complex<double> ProjectedState::psi_over_psi2(int site1, int site2, int new_sz1, int new_sz2) {
	MKL_Complex16 result;
	//new_sz1 = configuration[site2], new_sz2 = configuration[site1];

	int spin_row1 = Spin_t_to_row(new_sz1), spin_row2 = Spin_t_to_row(new_sz2);

	result = Winv[(site1 + spin_row1) * N + parton_labels[site1]] * Winv[(site2 + spin_row2) * N + parton_labels[site2]]
		- Winv[(site1 + spin_row1) * N + parton_labels[site2]] * Winv[(site2 + spin_row2) * N + parton_labels[site1]];

	std::vector<int> flip_sites = { site1, site2 }, new_sz = { new_sz1, new_sz2 };
	return result * jastrow.lazy_eval(flip_sites, new_sz, configuration);
}

// 3-site ring exchange, with jastrow
std::complex<double> ProjectedState::psi_over_psi_swap(int site1, int site2, int site3) {
	MKL_Complex16 result;
	int new_sz1 = configuration[site2], new_sz2 = configuration[site3], new_sz3 = configuration[site1];
	int spin_row1 = Spin_t_to_row(new_sz1), spin_row2 = Spin_t_to_row(new_sz2), spin_row3 = Spin_t_to_row(new_sz3);

	std::complex<double> 
		row1_1 = Winv[(site1 + spin_row1) * N + parton_labels[site1]],
		row1_2 = Winv[(site1 + spin_row1) * N + parton_labels[site2]],
		row1_3 = Winv[(site1 + spin_row1) * N + parton_labels[site3]];

	result = row1_1 * (Winv[(site2 + spin_row2) * N + parton_labels[site2]] * Winv[(site3 + spin_row3) * N + parton_labels[site3]]
					 - Winv[(site3 + spin_row3) * N + parton_labels[site2]] * Winv[(site2 + spin_row2) * N + parton_labels[site3]])
		- row1_2 * (Winv[(site2 + spin_row2) * N + parton_labels[site1]] * Winv[(site3 + spin_row3) * N + parton_labels[site3]]
				- Winv[(site3 + spin_row3) * N + parton_labels[site1]] * Winv[(site2 + spin_row2) * N + parton_labels[site3]])
		+ row1_3 * (Winv[(site2 + spin_row2) * N + parton_labels[site1]] * Winv[(site3 + spin_row3) * N + parton_labels[site2]]
				- Winv[(site3 + spin_row3) * N + parton_labels[site1]] * Winv[(site2 + spin_row2) * N + parton_labels[site2]]);

	std::vector<int> flip_sites = { site1, site2, site3 }, new_sz = { new_sz1, new_sz2, new_sz3 };
	return result * jastrow.lazy_eval(flip_sites, new_sz, configuration);
}

/// OVERRIDE MATRIX ELEMENTS

// can swap spins at 2 or 3 sites given in ring_swap
std::complex<double> ProjectedState::psi_over_psi(std::vector<int>& ring_swap) {
		std::vector<int> sz(ring_swap.size());
		for (int i = 0; i < ring_swap.size(); ++i) {
			if (i == ring_swap.size() - 1) {
				sz[i] = configuration[ring_swap[0]];
			}
			else {
				sz[i] = configuration[ring_swap[i + 1]];
			}
		}
		return psi_over_psi(ring_swap, sz);
	}

// chooses a ring swap or a 2-site swap, potentially with an additional spin flip for the 2-site swap
std::complex<double> ProjectedState::psi_over_psi(std::vector<int>& flips, std::vector<int>& new_sz) {

	std::complex<double> result(1.0, 0.0);
	if (flips.size() == 2) {
		result *= psi_over_psi2(flips[0], flips[1], new_sz[0], new_sz[1]);
	}
	else if (flips.size() == 3) {
		result *= psi_over_psi_swap(flips[0], flips[1], flips[2]);
	}
	else {
		assert(flips.size() == 2);
	}
	return result;

}

/// PRIVATE UPDATES - ACTUAL CALLERS

void ProjectedState::update(int site1, int site2, std::complex<double> psioverpsi) {
	//only for swapping spins at two different sites

	if (configuration[site1] < configuration[site2]) {
		upinvhop2((site1 + Spin_t_to_row(configuration[site2])), parton_labels[site2], (site2 + Spin_t_to_row(configuration[site1])), parton_labels[site1]);
	}
	else {
		upinvhop2((site2 + Spin_t_to_row(configuration[site1])), parton_labels[site1], (site1 + Spin_t_to_row(configuration[site2])), parton_labels[site2]);
	}
	//std::cout << "labels: " << vec2str(parton_labels) << "\n";
	int templabel = parton_labels[site1];
	parton_labels[site1] = parton_labels[site2];
	parton_labels[site2] = templabel;
	templabel = configuration[site1];
	configuration[site1] = configuration[site2];
	configuration[site2] = templabel;
	det *= psioverpsi;
}

void ProjectedState::update(std::vector<int>& sites, std::vector<int>& new_sz, std::complex<double> psioverpsi) {

	assert(sites.size() == 2);
	assert(new_sz.size() == 2);

	if (configuration[sites[0]] < configuration[sites[1]]) {
		upinvhop2((sites[0] + Spin_t_to_row(new_sz[0])), parton_labels[sites[0]], (sites[1] + Spin_t_to_row(new_sz[1])), parton_labels[sites[1]]);
	}
	else {
		upinvhop2((sites[1] + Spin_t_to_row(new_sz[1])), parton_labels[sites[1]], (sites[0] + Spin_t_to_row(new_sz[0])), parton_labels[sites[0]]);
	}

	configuration[sites[0]] = new_sz[0];
	configuration[sites[1]] = new_sz[1];
	det *= psioverpsi;
}

/// OVERRIDE UPDATES - PUBLIC INTERFACE

// pass-thru, analogous to psi_over_psi(ring_swap)
void ProjectedState::update(std::vector<int>& ring_swap) {
	std::vector<int> sz(ring_swap.size());
	for (int i = 0; i < ring_swap.size(); ++i) {
		if (i == ring_swap.size() - 1) {
			sz[i] = configuration[ring_swap[0]];
		}
		else {
			sz[i] = configuration[ring_swap[i + 1]];
		}
	}
	update(ring_swap, sz);
}

void ProjectedState::update(std::vector<int>& flips, std::vector<int>& new_sz) {
	std::complex<double> pop;

	if (flips.size() == 2) {
		jastrow.update_tables(flips, new_sz, configuration);
		pop = psi_over_psi(flips, new_sz);
		if (configuration[flips[0]] == new_sz[1] && configuration[flips[1]] == new_sz[0]) {
			update(flips[0], flips[1], pop);
		}
		else {
			update(flips, new_sz, pop);
		}
	}
	else if (flips.size() == 3) {
		assert(!jastrow.exist()); //jastrow not implemented for ring exchanges
		std::vector<int> fliplist(2), spinlist(2);
		fliplist = { flips[0], flips[1] };
		spinlist = { new_sz[0], new_sz[2] };
		pop = psi_over_psi(fliplist, spinlist);
		update(flips[0], flips[1], pop);

		fliplist = { flips[1], flips[2] };
		spinlist = { new_sz[1], new_sz[2] };
		pop = psi_over_psi(fliplist, spinlist);
		update(flips[1], flips[2], pop);
	}
	else {
		std::stringstream ss;
		ss << "Spin exchanges with <2 or >3 sites not implemented.\nFlip list:\n";
		for (int i = 0; i < flips.size(); ++i) {
			ss << flips[i] << ",";
		}
		ss << "\n";
		throw vmctype::NotImplemented(ss.str());
	}
}

//Tests

bool ProjectedState::test_2_spin_swap_pop(bool output) {
	assert(!jastrow.exist()); //test jastrow separately
	bool success = false;

	//Choose two random sites
	std::vector<int> flips = rand.get_rand_vec_site(2);
	std::vector<int> new_sz = { configuration[flips[1]], configuration[flips[0]] };

	//Calculate psi fast and slow
	std::complex<double> pfast, pslow, oldpsi = calc_det();
	//fast
	pfast = psi_over_psi2(flips[0], flips[1], configuration[flips[1]], configuration[flips[0]]);
	//slow
	if (std::abs(pfast) > 1e-10) {
		update(flips[0], flips[1], pfast);
		set_configuration(configuration);
		pslow = calc_det() / oldpsi;

		if (output && std::abs(std::abs(pslow) - std::abs(pfast)) > 1e-8 && std::abs(pfast) > 1e-16) {
			std::cout << "Test 2 spin swap psi over psi\n";
			std::cout << "Swap sites " << flips[0] << ", " << flips[1] << "\n";
			std::cout << "With sz = " << configuration[flips[0]] << ", " << configuration[flips[1]] << "\n";
			std::cout << "And new_sz = " << new_sz[0] << ", " << new_sz[1] << "\n";
			std::cout << "psi fast = " << pfast.real() << " + " << pfast.imag() << "i\n";
			std::cout << "psi slow = " << pslow.real() << " + " << pslow.imag() << "i\n";
		}
	}

	

	return success;
}

bool ProjectedState::test_2_spin_flip_pop(std::vector<int>& flips, std::vector<int>& new_sz) {
	//std::cout << "Warning: set all Jastrow factors = 0 for accurate results\n";
	bool success = false;
	std::vector<int> old_sz = { configuration[flips[0]], configuration[flips[1]] };

	//Calculate psi fast and slow
	std::complex<double> pfast, pslow, oldpsi = calc_det();
	//fast
	pfast = psi_over_psi(flips, new_sz);
	//slow
	if (std::abs(pfast) > 1e-10) {
		update(flips, new_sz, pfast);
		set_configuration(configuration);
		pslow = calc_det() / oldpsi;

		if (std::abs(std::abs(pslow) - std::abs(pfast)) > 1e-8 && std::abs(pfast) > 1e-16) {
			std::cout << "Test 2 spin flip psi over psi\n";
			std::cout << "Swap sites " << flips[0] << ", " << flips[1] << "\n";
			std::cout << "With sz = " << old_sz[0] << ", " << old_sz[1] << "\n";
			std::cout << "And new_sz = " << new_sz[0] << ", " << new_sz[1] << "\n";
			std::cout << "psi fast = " << pfast.real() << " + " << pfast.imag() << "i\n";
			std::cout << "psi slow = " << pslow.real() << " + " << pslow.imag() << "i\n";
		}
	}



	return success;
}

bool ProjectedState::test_3_spin_swap_pop(bool output) {
	assert(!jastrow.exist()); //test jastrow separately
	int config_attempt = 0;
	while (!try_configuration() && config_attempt < 50) {
		det = { 0, 0 };
		++config_attempt;
	}
	assert(config_attempt < 50);
	bool success = false;

	//Choose three random sites
	int site = rand.get_rand_site();
	std::vector<int> flips = { site, (site + 1) % N, (site + 2) % N };
	std::vector<int> new_sz = { configuration[flips[1]], configuration[flips[2]], configuration[flips[0]] };

	//Calculate psi fast and slow
	std::complex<double> pfast, pslow, oldpsi = calc_det();
	//fast
	pfast = psi_over_psi_swap(flips[0], flips[1], flips[2]);

	//slow
	if (std::abs(pfast) > 1e-10) {
		update(flips);
		set_configuration(configuration);
		pslow = calc_det() / oldpsi;

		if (output && std::abs(std::abs(pslow) - std::abs(pfast)) > 1e-8 && std::abs(pfast) > 1e-16) {
			std::cout << "Test 3 spin swap psi over psi\n";
			std::cout << "Starting with psi = " << det << "\n";
			std::cout << "Swap sites " << flips[0] << ", " << flips[1] << ", " << flips[2] << "\n";
			std::cout << "With sz = " << configuration[flips[0]] << ", " << configuration[flips[1]] << ", " << configuration[flips[2]] << "\n";
			std::cout << "And new_sz = " << new_sz[0] << ", " << new_sz[1] << ", " << new_sz[2] << "\n";
			std::cout << "psi fast = " << pfast.real() << " + " << pfast.imag() << "i\n";
			std::cout << "psi slow = " << pslow.real() << " + " << pslow.imag() << "i\n";
		}
	}

	return success;
}

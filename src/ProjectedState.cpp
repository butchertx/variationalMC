#include "ProjectedState.h"

void print_matrix(const char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda) {
	std::cout << desc << ":\n";
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n-1; ++j) {
			std::cout << a[i * lda + j] << ",";
		}
		std::cout << a[i * lda + n - 1] << "\n";
	}
}

ProjectedState::ProjectedState(MeanFieldAnsatz& M_, RandomEngine& rand_)
	: ansatz(M_), rand(rand_) {
	initialize_matrices();
	jastrow = {};
	assert(!jastrow.exist());
}

ProjectedState::ProjectedState(MeanFieldAnsatz& M, RandomEngine& rand_in, JastrowTable jastrow_)
	: ansatz(M), rand(rand_in), jastrow(jastrow_) {
	initialize_matrices();
	jastrow.initialize_tables(configuration);
}

void ProjectedState::initialize_matrices(){
	//construct Slater Matrix and set up matrix element calculation
	N = ansatz.get_dim() / 3;
	Slater = (lapack_complex_double*)mkl_malloc(N * N * sizeof(lapack_complex_double), 64);
	LU = (lapack_complex_double*)mkl_malloc(N * N * sizeof(lapack_complex_double), 64);
	Winv = (lapack_complex_double*)mkl_malloc(3 * N * N * sizeof(lapack_complex_double), 64);
	UP1 = (lapack_complex_double*)mkl_malloc(3 * N * 2 * sizeof(lapack_complex_double), 64);
	UP2 = (lapack_complex_double*)mkl_malloc(N * 2 * sizeof(lapack_complex_double), 64);
	UP3 = (lapack_complex_double*)mkl_malloc(N * 2 * sizeof(lapack_complex_double), 64);
	ipiv = (lapack_int*)mkl_malloc(N * N * sizeof(lapack_int), 64);
	for (int i = 0; i < N * N; ++i) {
		Slater[i] = { 0,0 };
		LU[i] = { 0,0 };
		ipiv[i] = 0;
		for (int j = 0; j < 3; ++j) {
			Winv[i + j * N * N] = { 0,0 };
		}
		if (i < 3 * N * 2) {
			UP1[i] = { 0.0, 0.0 };
			if (i < N * 2) {
				UP2[i] = { 0.0, 0.0 };
				UP3[i] = { 0.0, 0.0 };
			}
		}
	}
	int config_attempt = 0;
	while (!try_configuration() && config_attempt < 50) {
		det = { 0, 0 };
		++config_attempt;
	}
	//assert(config_attempt < 50);
	std::cout << "Starting with psi = " << det << "\n";
}

std::complex<double> ProjectedState::psi_over_psi2(int site1, int site2, int new_sz1, int new_sz2) {
	MKL_Complex16 result;
	//new_sz1 = configuration[site2], new_sz2 = configuration[site1];

	result = Winv[(site1 + (1 - new_sz1) * N) * N + parton_labels[site1]] * Winv[(site2 + (1 - new_sz2) * N) * N + parton_labels[site2]]
		- Winv[(site1 + (1 - new_sz1) * N) * N + parton_labels[site2]] * Winv[(site2 + (1 - new_sz2) * N) * N + parton_labels[site1]];

	std::vector<int> flip_sites = { site1, site2 }, new_sz = { new_sz1, new_sz2 };
	return result * jastrow.lazy_eval(flip_sites, new_sz, configuration);
}

std::complex<double> ProjectedState::psi_over_psi_swap(int site1, int site2, int site3) {
	MKL_Complex16 result;
	int new_sz1 = configuration[site2], new_sz2 = configuration[site3], new_sz3 = configuration[site1];

	std::complex<double> row1_1 = Winv[(site1 + (1 - new_sz1) * N) * N + parton_labels[site1]],
		row1_2 = Winv[(site1 + (1 - new_sz1) * N) * N + parton_labels[site2]],
		row1_3 = Winv[(site1 + (1 - new_sz1) * N) * N + parton_labels[site3]];

	result = row1_1 * (Winv[(site2 + (1 - new_sz2) * N) * N + parton_labels[site2]] * Winv[(site3 + (1 - new_sz3) * N) * N + parton_labels[site3]]
		- Winv[(site3 + (1 - new_sz3) * N) * N + parton_labels[site2]] * Winv[(site2 + (1 - new_sz2) * N) * N + parton_labels[site3]])
		- row1_2 * (Winv[(site2 + (1 - new_sz2) * N) * N + parton_labels[site1]] * Winv[(site3 + (1 - new_sz3) * N) * N + parton_labels[site3]]
			- Winv[(site3 + (1 - new_sz3) * N) * N + parton_labels[site1]] * Winv[(site2 + (1 - new_sz2) * N) * N + parton_labels[site3]])
		+ row1_3 * (Winv[(site2 + (1 - new_sz2) * N) * N + parton_labels[site1]] * Winv[(site3 + (1 - new_sz3) * N) * N + parton_labels[site2]]
			- Winv[(site3 + (1 - new_sz3) * N) * N + parton_labels[site1]] * Winv[(site2 + (1 - new_sz2) * N) * N + parton_labels[site2]]);

	std::vector<int> flip_sites = { site1, site2, site3 }, new_sz = { new_sz1, new_sz2, new_sz3 };
	return result * jastrow.lazy_eval(flip_sites, new_sz, configuration);
}

std::complex<double> ProjectedState::psi_over_psi(int site1, int site2, int site3) {
	assert(!jastrow.exist()); //jastrow not implemented for ring exchanges

	//swap sites 1 and 2, with labels 1 and 2, newsz 1 and 3
	std::complex<double> detswap1 = psi_over_psi2(site1, site2, configuration[site2], configuration[site1]);

	if (std::abs(detswap1) > 10e-16) {
		update(site1, site2, detswap1);
		std::complex<double> detswap2 = psi_over_psi2(site2, site3, configuration[site3], configuration[site2]);
		update(site2, site1, 1.0 / detswap1);
		return detswap1 * detswap2;
	}
	else { 
		detswap1 = psi_over_psi2(site2, site3, configuration[site3], configuration[site2]);
		if (std::abs(detswap1) > 10e-16) {
			update(site2, site3, detswap1);
			//std::cout << "successful update 1, detswap = " << detswap1 << ", 1/detswap = " << 1.0 / detswap1 << "\n";
			//swap sites 2 and 3, with labels 1 and 3, newsz 2 and 3
			std::complex<double> detswap2 = psi_over_psi2(site3, site1, configuration[site1], configuration[site3]);
			update(site3, site2, 1.0 / detswap1);
			//std::cout << "successful update 2, detswap = " << detswap2 << "\n";
			return detswap1 * detswap2;
		}
		else {
			return { 0.0, 0.0 };
		}
	}
}

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
	//conf_copy = configuration;
	//for (int i = 0; i < flips.size(); ++i) {
	//	conf_copy[flips[i]] = new_sz[i];
	//}
	//double newj, oldj;
	//oldj = jastrow.logpsi(configuration);
	//newj = jastrow.logpsi(conf_copy);
	return result;// *exp(jastrow.logpsi_over_psi(flips, configuration));

}


void ProjectedState::update(std::vector<int>& flips, std::vector<int>& new_sz) {
	std::complex<double> pop;

	if (flips.size() == 2) {
		/*std::cout << "Lazy and Greedy jastrows: \n" << jastrow.lazy_eval(flips, new_sz, configuration); 
		double greedy_old = jastrow.greedy_eval(configuration);*/
		jastrow.update_tables(flips, new_sz, configuration);
		pop = psi_over_psi(flips, new_sz);
		if (configuration[flips[0]] == new_sz[1] && configuration[flips[1]] == new_sz[0]) {
			update(flips[0], flips[1], pop);
		}
		else {
			//test_2_spin_flip_pop(flips, new_sz);
			update(flips, new_sz, pop);
		}
		/*double greedy_new = jastrow.greedy_eval(configuration);
		std::cout << ", " << greedy_new / greedy_old << "\n";*/
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
		std::cout << "Error: did not provide 2 or 3 elements to update Projected State. Flip list:\n";
		for (int i = 0; i < flips.size(); ++i) {
			std::cout << flips[i] << ",";
		}
		std::cout << "\n";
	}
}

void ProjectedState::update(int site1, int site2, std::complex<double> psioverpsi) {
	//only for swapping spins at two different sites

	if (configuration[site1] < configuration[site2]) {
		upinvhop2((site1 + (1 - configuration[site2]) * N), parton_labels[site2], (site2 + (1 - configuration[site1]) * N), parton_labels[site1]);
	}
	else {
		upinvhop2((site2 + (1 - configuration[site1]) * N), parton_labels[site1], (site1 + (1 - configuration[site2]) * N), parton_labels[site2]);
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
		upinvhop2_flip((sites[0] + (1 - new_sz[0]) * N), parton_labels[sites[0]], (sites[1] + (1 - new_sz[1]) * N), parton_labels[sites[1]]);
	}
	else {
		upinvhop2_flip((sites[1] + (1 - new_sz[1]) * N), parton_labels[sites[1]], (sites[0] + (1 - new_sz[0]) * N), parton_labels[sites[0]]);
	}

	configuration[sites[0]] = new_sz[0];
	configuration[sites[1]] = new_sz[1];
	det *= psioverpsi;
}

bool ProjectedState::try_configuration() {
	int N0 = ansatz.get_N0F();
	if (2 * ((N - N0) / 2) != N - N0) {
		N0 += 1;
	}
	set_configuration(rand.get_rand_spin_state(std::vector<int>{ (N - N0) / 2, N0, (N - N0) / 2 }, N));
	
	CBLAS_INDEX low = 0;
	low = cblas_izamin(N, LU, N+1);
	return (cblas_dcabs1(&(LU[low*(N+1)])) > 10e-10);
}

void ProjectedState::set_configuration(std::vector<int> conf) {
	configuration = conf;
	int row = 0, N3 = 3 * N, mi = 0, fs_start = ansatz.get_fermi_surface_start(), fs_size = ansatz.get_fermi_surface_end()-fs_start, rand_orbital = 0;
	std::vector<int> fs_orbitals;
	for (int i = 0; i < fs_size; ++i) {
		fs_orbitals.push_back(i);
	}
	bool selected = false;
	lapack_complex_double* phi = ansatz.get_Phi();
	parton_labels.clear();
	//print_matrix("Orbitals 3Nx10", 3*N, 10, phi, 3*N);
	lapack_int info;
	MKL_Complex16 alpha = { 1.0, 0.0 }, beta = { 0.0, 0.0 };
	for (int i = 0; i < N; ++i) {
		parton_labels.push_back(i);
		mi = configuration[i];
		row = (-mi + 1) * N + i;
		//std::cout << "Site " << i << ", Sz = " << conf[i] << ", row " << row << "\n";
		//print_matrix("Row: ", 1, 10, &(phi[row*3*N]), 3*N);
		std::memcpy(&(Slater[parton_labels[i] * N]), &(phi[row * 3 * N]), N * sizeof(lapack_complex_double));
	}
	info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, Slater, N, ipiv);
	std::memcpy(LU, Slater, N * N * sizeof(lapack_complex_double));
	if (info == 0) {
		info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, N, Slater, N, ipiv);
		//zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3 * N, N, N, &alpha, phi, 3 * N, Slater, N, &beta, Winv, 3 * N);
		zgemm3m("N", "N", &N, &N3, &N, &alpha, Slater, &N, phi, &N3, &beta, Winv, &N);
		//print_matrix("Winv: ", 10,10, Winv, N);
	}
	det = calc_det();
}

std::complex<double> ProjectedState::calc_det() {
	std::complex<double> result = { 1.0, 0.0 };
	for (int i = 0; i < N; ++i) {
		//std::cout << "(" << LU[i*N + i].real << "," << LU[i*N + i].imag << ")\n";
		result *= LU[i*N + i];
	}
	return result;
}

void ProjectedState::upinvhop2(int rowk, int colk, int rowl, int coll) {

	int N3 = 3 * N;

	std::complex<double> 
		c11 = Winv[rowl*N + coll],
		c22 = Winv[rowk*N + colk],
		c12 = -Winv[rowk*N + coll],
		c21 = -Winv[rowl*N + colk];

	std::complex<double> g = c11 * c22 - c12 * c21, beta = { 1.0, 0.0 };

	//if (std::abs(g) < 10e-20) {
	//	std::cout << "Invalid Update Requested, g = " << g << "; Check Acceptance Probability\n";
	//	std::cout << "Program Finished, Press Enter to exit\n";
	//	std::string trash;
	//	std::getline(std::cin, trash);
	//	exit(0);
	//}

	cblas_zcopy(N3, &(Winv[colk]), N, UP1, 2);
	cblas_zcopy(N3, &(Winv[coll]), N, &(UP1[1]), 2);
	//print_matrix("UP1", 3 * N, 2, UP1, 2);
	cblas_zcopy(N, &(Winv[rowk*N]), 1, UP2, 1);
	UP2[colk] -= std::complex<double>(1.0, 0.0);
	cblas_zcopy(N, &(Winv[rowl*N]), 1, &(UP2[N]), 1);
	UP2[N + coll] -= std::complex<double>(1.0, 0.0);

	g = std::complex<double>({ -1.0, 0.0 }) / g;

	for (int i = 0; i < N; ++i) {
		UP3[i] = c11 * UP2[i] + c12 * UP2[N + i];
		UP3[N + i] = c21 * UP2[i] + c22 * UP2[N + i];
	}
	//print_matrix("UP1", 3 * N, 2, UP1, 2);
	//print_matrix("UP2", 2, N, UP2, N);
	cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, N3, N, 2, &g, UP1, 2, UP3, N, &beta, Winv, N);
}

void ProjectedState::upinvhop2_flip(int rowk, int colk, int rowl, int coll) {
	//perform the update of Winv according to the Woodbury Matrix identity
	//Winv' = Winv - Winv * U (I_k + V A^-1 U)^-1 V A^-1
	//where U and V are defined so A' = A + UV (A is the Slater matrix)
	//UP1 is Winv * U; the i and j columns of Winv where i and j are the sites to update
	//UP2 is V A^-1 (which can be computed easily from rows of Winv)
	//UP3 is (I_k + V A^-1 U)^-1 (2x2) times UP2 (2xN)

	int N3 = 3 * N;

	std::complex<double>
		c11 = Winv[rowl * N + coll],
		c22 = Winv[rowk * N + colk],
		c12 = -Winv[rowk * N + coll],
		c21 = -Winv[rowl * N + colk];

	std::complex<double> g = c11 * c22 - c12 * c21, beta = { 1.0, 0.0 };

	cblas_zcopy(N3, &(Winv[colk]), N, UP1, 2);
	cblas_zcopy(N3, &(Winv[coll]), N, &(UP1[1]), 2);
	cblas_zcopy(N, &(Winv[rowk * N]), 1, UP2, 1);
	UP2[colk] -= std::complex<double>(1.0, 0.0);
	cblas_zcopy(N, &(Winv[rowl * N]), 1, &(UP2[N]), 1);
	UP2[N + coll] -= std::complex<double>(1.0, 0.0);

	g = std::complex<double>({ -1.0, 0.0 }) / g;

	for (int i = 0; i < N; ++i) {
		UP3[i] = c11 * UP2[i] + c12 * UP2[N + i];
		UP3[N + i] = c21 * UP2[i] + c22 * UP2[N + i];
	}

	cblas_zgemm3m(CblasRowMajor, CblasNoTrans, CblasNoTrans, N3, N, 2, &g, UP1, 2, UP3, N, &beta, Winv, N);
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

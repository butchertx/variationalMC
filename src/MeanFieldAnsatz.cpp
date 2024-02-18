#include "MeanFieldAnsatz.h"
#include "Lattice.h"
#include "RandomEngine.h"
#include <cstring>
#include <cstdlib>

//void print_matrix(const char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda);

class origin_direction_pair {
public:
	int origin;
	double phase;
	vec3<double> direction;

	origin_direction_pair(int o_, double x_, double y_, double z_, double phase_) : origin(o_), direction(vec3<double>(x_, y_, z_)), phase(phase_) {};
	origin_direction_pair(int o_, vec3<double> vec_, double phase_) : origin(o_), direction(vec_), phase(phase_) {};

};

void TightBindingSitePair::get_couplings(double * tzr, double * tzi, double * txyr, double * txyi) {
	int pbc_sign = x2 >= 0 ? 1 : -1;
	*tzr = pbc_sign * tz * cos(2 * M_PI * tz_phase);
	*tzi = pbc_sign * tz * sin(2 * M_PI * tz_phase);
	*txyr = pbc_sign * txy * cos(2 * M_PI * txy_phase);
	*txyi = pbc_sign * txy * sin(2 * M_PI * txy_phase);
}

//
//MeanFieldAnsatz::MeanFieldAnsatz(int Lx_in, int Ly_in, TightBindingUnitCell uc_in){
//	//for Lx,Ly, for site pairs in uc_in:
//	//shift indices, create new site pair with the given params
//	int N = Lx_in * Ly_in;
//	int shift;
//	int site1, site2;
//	site_pair_list.push_back({});
//	for (int x = 0; x < Lx_in; ++x) {
//		for (int y = 0; y < Ly_in; ++y) {
//			shift = x*Ly_in + y;
//			for (auto pair : uc_in.site_pairs) {
//				pair.get_sites(&site1, &site2);
//				pair.set_sites((site1 + shift)%N, (site2 + shift)%N);
//				site_pair_list[0].push_back(pair);
//			}
//		}
//	}
//}

MeanFieldAnsatz::MeanFieldAnsatz(WavefunctionOptions& mf_in, Lattice& lat_in, bool unit_cell_construction) 
	: N(lat_in.get_N()), su3_symmetry(mf_in.su3_symmetry), mu_z(mf_in.mu_z){

	//compatibility conditions:
	//1.  Lattice types are the same
	//2.  all hopping elements have valid connections
	assert(mf_in.lattice_type.compare(Lattice_type_to_string(lat_in.get_lattice_type())) == 0);

	//std::vector<std::vector<origin_direction_pair>> hopping_list(1);
	//std::vector<vec3<std::complex<double>>> directors;

	//TightBindingUnitCell tb_cell;
	double tz, tzp, txy, txyp;
	int origin, neighbor;
	TightBindingSitePair tbp(0,0);
	std::vector<std::vector<int>> basis_partition = lat_in.basis_partition(mf_in.basis);

	for (auto hopterm : mf_in.hopping_list) {//int hop = 0; hop < mf_in.hopping_list.size(); ++hop) {
		site_pair_list.push_back({});
		tz = hopterm.spin_row == 0 ? hopterm.strength : 0.0;
		tzp = 0.0;
		txy = hopterm.spin_row == 1 ? hopterm.strength : 0.0;
		txyp = 0.0;
		for (int uc = 0; uc < basis_partition.size(); ++uc) {
			for (int termind = 0; termind < hopterm.origins.size(); ++termind) {
				origin = basis_partition[uc][hopterm.origins[termind]];
				neighbor = lat_in.get_neighbor_with_pbc(origin, hopterm.distance, hopterm.neighbor_index[termind]);
				neighbor = mf_in.match_lattice_pbc ? neighbor : abs(neighbor);
				tbp = TightBindingSitePair(origin, neighbor, tz, txy, tzp + hopterm.phases[termind] / 360.0, txyp + hopterm.phases[termind] / 360.0);
				site_pair_list[site_pair_list.size() - 1].push_back(tbp);
				tbp.conjugate();
				site_pair_list[site_pair_list.size() - 1].push_back(tbp);
			}
		}

		std::cout << "TB pair list size: " << site_pair_list[site_pair_list.size() - 1].size() << "\n";

	}


	if (mf_in.directors.unit_cell_u_polar.size() > 0) {
		field = mf_in.field;
		for (auto uc : basis_partition) {
			for (int site = 0; site < uc.size(); ++site) {
				directors.push_back(mf_in.directors.eval_at(site, lat_in.get_coordinate(uc[0])));
			}
		}
	}

	HMF = (lapack_complex_double *)mkl_malloc(9 * N * N * sizeof(lapack_complex_double), 64);
	Phi = (lapack_complex_double *)mkl_malloc(9 * N * N * sizeof(lapack_complex_double), 64);
	Energy = (double *)mkl_malloc(3 * N * sizeof(double), 64);
	set_hamiltonian();
	diagonalize_hamiltonian();
	set_fermi_surface();
}

void MeanFieldAnsatz::set_hamiltonian() {
	//Use row major form (n = row*num_cols + col)
	//initialize to zero
	for (int i = 0; i < 9 * N * N; ++i) {
		HMF[i] = { 0.0, 0.0 };
	}

	//iterate TB pairs and add upper triangle values
	int row, col, dim = 3*N;
	double tzr, tzi, txyr, txyi;
	for (int hop_class = 0; hop_class < site_pair_list.size(); ++hop_class) {
		for (auto tb_pair : site_pair_list[hop_class]) {
			tb_pair.get_sites(&row, &col);
			if (row <= col) {
				tb_pair.get_couplings(&tzr, &tzi, &txyr, &txyi);
				HMF[row*dim + col] -= std::complex<double>(txyr, txyi);
				HMF[(row + N)*dim + col + N] -= std::complex<double>(tzr, tzi);
				HMF[(row + 2 * N)*dim + col + 2 * N] -= std::complex<double>(txyr,txyi);
			}
		}
	}

	//SU(3) irrep terms (kagome only)
	//std::string irrep = "octet";
	//for (int i = 0; i < N; ++i) {
	//	for (int m1 = 0; m1 < 3; ++m1) {
	//		for (int m2 = 0; m2 < 3; ++m2) {
	//			row = i + m1 * N;
	//			col = i + m2 * N;
	//			if (row <= col) {
	//				HMF[row * dim + col] -= 0.0;// -0.05 * get_su3_element(irrep, m1, m2, i);
	//			}
	//		}
	//	}
	//}

	//add director terms
	if (directors.size() > 0) {
		for (int i = 0; i < N; ++i) {
			for (int m1 = 0; m1 < 3; ++m1) {
				for (int m2 = 0; m2 < 3; ++m2) {
					row = i + m1 * N;
					col = i + m2 * N;
					if (row <= col) {
						HMF[row*dim + col] -= field * get_director_element(directors[i], m1, m2);
						if (m1 == 1 && m2 == 1) {
							HMF[row * dim + col] -= mu_z;
						}
					}
				}
			}
		}
	}

}

double MeanFieldAnsatz::get_su3_element(std::string irrep, int m1, int m2, int i) {
	double s1[3][3] = { {8.0, 4.0, -2.0}, {4.0, 8.0, -2.0}, {-2.0, -2.0, 8.0} };
	double s2[3][3] = { {8.0, -2.0, 4.0}, {-2.0, 8.0, -2.0}, {4.0, -2.0, 8.0} };
	double s3[3][3] = { {8.0, -2.0, -2.0}, {-2.0, 8.0, 4.0}, {-2.0, 4.0, 8.0} };
	if (i % 3 == 0) {
		return s1[m1][m2];
	}
	else if (i % 3 == 1) {
		return s2[m1][m2];
	}
	else if (i % 3 == 2) {
		return s3[m1][m2];
	}
	else {
		return 0.0;
	}
}

std::complex<double> MeanFieldAnsatz::get_director_element(vec3<std::complex<double>> d, int m1, int m2) {
	assert(m1 >= 0 && m1 <= 2);
	assert(m2 >= 0 && m2 <= 2);
	vec3<double> u(std::real(d.x), std::real(d.y), std::real(d.z)), v(std::imag(d.x), std::imag(d.y), std::imag(d.z));
	switch (m1) {
	case 0:
		switch (m2) {
		case 0:
			return { 0.5 * ((u.y - v.x)*(u.y - v.x) + (u.x + v.y)*(u.x + v.y)), 0.0 };
		case 1:
			return { -1 / sqrt(2)*((u.x + v.y)*u.z + v.z*(v.x - u.y)), -1 / sqrt(2)*(-v.z*(u.x + v.y) + u.z*(v.x - u.y)) };
		case 2:
			return { 0.5 * (-u.x*u.x - v.x*v.x + u.y*u.y + v.y*v.y), u.y*u.x + v.y*v.x };
		}
	case 1:
		switch (m2) {
		case 0:
			return { -1 / sqrt(2)*((u.x + v.y)*u.z + v.z*(v.x - u.y)), 1 / sqrt(2)*(-v.z*(u.x + v.y) + u.z*(v.x - u.y)) };
		case 1:
			return { u.z*u.z + v.z*v.z, 0.0 };
		case 2:
			return { 1 / sqrt(2)*((u.x - v.y)*u.z + v.z*(v.x + u.y)), 1 / sqrt(2)*(v.z*(u.x - v.y) - u.z*(v.x + u.y)) };
		}
	case 2:
		switch (m2) {
		case 0:
			return { 0.5 * (-u.x*u.x - v.x*v.x + u.y*u.y + v.y*v.y), -u.y*u.x + -v.y*v.x };
		case 1:
			return { 1 / sqrt(2)*((u.x - v.y)*u.z + v.z*(v.x + u.y)), -1 / sqrt(2)*(v.z*(u.x - v.y) - u.z*(v.x + u.y)) };
		case 2:
			return { 0.5 * ((u.y + v.x)*(u.y + v.x) + (u.x - v.y)*(u.x - v.y)), 0.0 };
		}
	default:
		return { 0.0, 0.0 };
	}
}

void MeanFieldAnsatz::diagonalize_hamiltonian() {
	std::memcpy(Phi, HMF, 9 * N * N * sizeof(lapack_complex_double));
	info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', 3*N, Phi, 3*N, Energy);
}

void MeanFieldAnsatz::print_levels() {
	MKL_INT i;
	double E, n1 = 0, n0 = 0, n_1 = 0,
		n1_tot = 0, n0_tot = 0, n_1_tot = 0;
	printf("\n %s\n", "Energies and Occupation numbers of Single-Particle Orbitals");
	for (i = 0; i < N+1; i++) {
		E = Energy[i];
		n1 = cblas_dznrm2(N, &(Phi[i]), 3 * N);
		n0 = cblas_dznrm2(N, &(Phi[i+3*N*N]), 3 * N);
		n_1 = cblas_dznrm2(N, &(Phi[i+6*N*N]), 3 * N);
		n1_tot += n1*n1;
		n0_tot += n0*n0;
		n_1_tot += n_1*n_1;
		if (i == N - 1) {
			n0_F = std::round(n0_tot);
		}
		printf(" (%6.2f,%6.2f,%6.2f,%6.2f)", E, n1_tot, n0_tot, n_1_tot);
		printf("\n");
	}
}

void MeanFieldAnsatz::print_fermi_level() {
	MKL_INT i;
	double E, n1 = 0, n0 = 0, n_1 = 0,
		n1_tot = 0, n0_tot = 0, n_1_tot = 0, Ef = Energy[N-1];
	bool subsurface = true, endsurface = true;
	printf("\n %s\n", "Energies and Occupation numbers of Single-Particle Orbitals At the Fermi Level");
	for (i = 0; i < 2*N; i++) {
		E = Energy[i];
		n1 = cblas_dznrm2(N, &(Phi[i]), 3 * N);
		n0 = cblas_dznrm2(N, &(Phi[i + 3 * N*N]), 3 * N);
		n_1 = cblas_dznrm2(N, &(Phi[i + 6 * N*N]), 3 * N);
		n1_tot += n1*n1;
		n0_tot += n0*n0;
		n_1_tot += n_1*n_1;
		if (std::abs(E - Ef) < EPSILON || i == N) {
			if (subsurface) {
				printf(" (%6d, %6.2f,%6.2f,%6.2f,%6.2f)", i, Energy[i-1], n1_tot - n1 * n1, n0_tot - n0 * n0, n_1_tot - n_1 * n_1);
				printf("\n");
				std::cout << "------------- FS starts here\n";
				subsurface = false;
			}

			printf(" (%6d, %6.2f,%6.2f,%6.2f,%6.2f)", i+1, E, n1*n1, n0*n0, n_1*n_1);
			printf("\n");
			if (i == N - 1) {
				std::cout << "------------- 1 Particle Per Site\n";
			}
		}
		if ( E - Ef > EPSILON && endsurface) {
			std::cout << "------------- FS ends here\n";
			printf(" (%6d, %6.2f,%6.2f,%6.2f,%6.2f)", i + 1, E, n1 * n1, n0 * n0, n_1 * n_1);
			printf("\n");
			endsurface = false;
		}
	}
}

void MeanFieldAnsatz::set_fermi_surface() {
	MKL_INT i;
	double E, n1 = 0, n0 = 0, n_1 = 0,
		n1_tot = 0, n0_tot = 0, n_1_tot = 0, Ef = Energy[N - 1];
	bool subsurface = true, endsurface = true;
	for (i = 0; i < 2 * N; i++) {
		E = Energy[i];
		n1 = cblas_dznrm2(N, &(Phi[i]), 3 * N);
		n0 = cblas_dznrm2(N, &(Phi[i + 3 * N * N]), 3 * N);
		n_1 = cblas_dznrm2(N, &(Phi[i + 6 * N * N]), 3 * N);
		n1_tot += n1 * n1;
		n0_tot += n0 * n0;
		n_1_tot += n_1 * n_1;
		if (std::abs(E - Ef) < EPSILON || i == N) {
			if (subsurface) {
				subsurface = false;
				fermi_surface_start = i;
				fermi = FermiSurface(Ef, { n1_tot, n0_tot, n_1_tot });
			}

			fermi.add_orbital(Orbital({ n1 * n1, n0 * n0, n_1 * n_1 }, E, i));
			if (i == N - 1) {
				n0_F = std::round(n0_tot);
			}
		}
		if (E - Ef > EPSILON&& endsurface) {
			fermi_surface_end = i;//index of first gapped orbital
			endsurface = false;
		}
	}
}

void MeanFieldAnsatz::shuffle_FS(int n0, int n1, RandomEngine* rand) {
	if (fermi_surface_end != n0 + 2 * n1) {
		lapack_complex_double* fs_temp_matrix = (lapack_complex_double*)mkl_malloc(3 * N * (fermi.get_size()) * sizeof(lapack_complex_double), 64);
		double n0fermi = n0 - fermi.get_inner_shell_count(0), n1fermi = n1 - fermi.get_inner_shell_count(1), n_1fermi = n1fermi;
		int fermisize = fermi.get_size();

		std::vector<int> fs_indices(fermi.get_size());
		for (int i = 0; i < fs_indices.size(); ++i) {
			fs_indices[i] = i;
		}

		std::vector<int> occ_indices;
		int temp_orb_index = 0;
		Orbital temp_orb;

		//select random orbitals from FS until all number requirements are met
		while (n0fermi + n1fermi + n_1fermi > EPSILON) {
			temp_orb_index = rand->get_rand_in_range(fs_indices.size());
			temp_orb = fermi.get_orbital(temp_orb_index);

			if (n1fermi - temp_orb.get_overlap(1) > -EPSILON
				&& n0fermi - temp_orb.get_overlap(0) > -EPSILON
				&& n_1fermi - temp_orb.get_overlap(-1) > -EPSILON) {

				occ_indices.push_back(temp_orb.get_index());
				fs_indices.erase(fs_indices.begin() + temp_orb_index);
				n1fermi -= temp_orb.get_overlap(1);
				n0fermi -= temp_orb.get_overlap(0);
				n_1fermi -= temp_orb.get_overlap(-1);
			}
		}

		assert(fermisize == fs_indices.size() + occ_indices.size());
		//copy selected columns (occ_indices) of phi into first columns of fs_temp_matrix and update orbitals
		for (int fs_index = 0; fs_index < occ_indices.size(); ++fs_index) {
			for (int row = 0; row < 3 * N; ++row) {
				fs_temp_matrix[row * fermisize + fs_index] = Phi[row * 3 * N + occ_indices[fs_index]];
			}
			fermi.update_index(fs_index, fermi_surface_start + fs_index);
		}

		//copy remaining fs_indices into remaining columns of fs_temp_matrix and update orbitals
		for (int fs_index = occ_indices.size(); fs_index < fs_indices.size() + occ_indices.size(); ++fs_index) {
			for (int row = 0; row < 3 * N; ++row) {
				fs_temp_matrix[row * fermisize + fs_index] = Phi[row * 3 * N + fs_indices[fs_index - occ_indices.size()]];
			}
			fermi.update_index(fs_index, fermi_surface_start + fs_index);
		}

		//copy columns of fs_temp_matrix back into corresponding columns of phi
		for (int row = 0; row < 3 * N; ++row) {
			std::memcpy(&(Phi[row * 3 * N + fermi_surface_start]), &(fs_temp_matrix[row * fermisize]), fermisize * sizeof(lapack_complex_double));
		}



		mkl_free(fs_temp_matrix);
	}
}

void MeanFieldAnsatz::write_levels(std::ofstream *f) {
	MKL_INT i;
	double E, n1 = 0, n0 = 0, n_1 = 0,
		n1_tot = 0, n0_tot = 0, n_1_tot = 0;
	*f << "Band Energy, Total N, cumul. N_-1, cumul. N_0, cumul. N_1, N_-1, N_0, N_1\n";
	char buf[1024];
	for (i = 0; i < N + 1; i++) {
		E = Energy[i];
		n1 = cblas_dznrm2(N, &(Phi[i]), 3 * N);
		n0 = cblas_dznrm2(N, &(Phi[i + 3 * N * N]), 3 * N);
		n_1 = cblas_dznrm2(N, &(Phi[i + 6 * N * N]), 3 * N);
		n1_tot += n1 * n1;
		n0_tot += n0 * n0;
		n_1_tot += n_1 * n_1;
		if (i == N - 1) {
			n0_F = std::round(n0_tot);
		}
		snprintf(buf, 200, "%6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f", E, n1_tot + n0_tot + n_1_tot, n1_tot, n0_tot, n_1_tot, n1 * n1, n0 * n0, n_1 * n_1);
		*f << std::string(buf) << "\n";
	}
}
void MeanFieldAnsatz::write_directors(std::ofstream* f) {
	*f << "site, ux, uy, uz, vx, vy, vz\n";
	for (auto site = 0; site < directors.size(); ++site) {
		*f << site << ", " << directors[site].x.real() << ", " << directors[site].y.real() << ", " << directors[site].z.real() << ", ";
		*f << directors[site].x.imag() << ", " << directors[site].y.imag() << ", " << directors[site].z.imag() << "\n";
	}
}

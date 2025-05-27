#include "ProjectedState.h"
#include "mkl_types.h"

// constructors

ProjectedState::ProjectedState(MeanFieldAnsatz& M_, RandomEngine& rand_)
	: ansatz(M_), rand(rand_), N(ansatz.get_N()), DIM(ansatz.get_dim()){
	// initialize the Slater-Jastrow state

	conserve_sz2 = M_.get_conserve_sz2();
	Slater = ComplexDoubleMatrix<MKL_Complex16>(N, N);
	Winv = ComplexDoubleMatrix<MKL_Complex16>(ansatz.get_dim(), N);
	clear_matrices();
	initialize_configuration();
}

ProjectedState::ProjectedState(MeanFieldAnsatz& M_, RandomEngine& rand_, JastrowTable jastrow_)
	: ProjectedState(M_, rand_){
	conserve_sz2 = M_.get_conserve_sz2();
	Slater = ComplexDoubleMatrix<MKL_Complex16>(N, N);
	Winv = ComplexDoubleMatrix<MKL_Complex16>(ansatz.get_dim(), N);
	jastrow = jastrow_;
	jastrow.initialize_tables(configuration);
}

void ProjectedState::clear_matrices(){
	Slater.clear_matrix();
	Winv.clear_matrix();
}

void ProjectedState::initialize_configuration(){
	int config_attempt = 0;
	while (!try_configuration() && config_attempt < CONFIG_ATTEMPTS) {
		++config_attempt;
	}
}

bool ProjectedState::try_configuration() {
	// TODO: do this in a more robust manner, maybe type checking on initialization
	// return True if we have a valid configuration, False otherwise
	int N0 = ansatz.get_N0F();
	if (2 * ((N - N0) / 2) != N - N0) {
		N0 += 1;
	}
	set_configuration(rand.get_rand_spin_state(std::vector<int>{ (N - N0) / 2, N0, (N - N0) / 2 }, N));

	return !(Winv.is_determinant_zero());
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
	auto phi = ansatz.get_Phi();
	parton_labels.clear();

	// move the relevant rows of Phi into Slater
	for (int i = 0; i < N; ++i) {
		parton_labels.push_back(i);
		row = Spin_t_to_row(configuration[i]) + i;
		Slater.copy_row(Slater, phi, parton_labels[i], row, N);
	}

	// invert Slater matrix
	auto Slater_inverse = Slater.compute_inverse();

	// Multiply phi into Slater^{-1}
	Winv = phi.get_slice(0, phi.rows(), 0, N) * Slater_inverse;

}

void ProjectedState::updateMatrixInverse(int rowk, int colk, int rowl, int coll) {
	//perform the update of Winv according to the Woodbury Matrix identity
	//Winv' = Phi * I_{dNxN} * A^{-1} * U (I_k + V A^-1 U)^-1 V A^-1
	//      = (Winv * U) * (I_k + V A^-1 U)^-1 * RatioDiff
	//where U and V are defined so A' = A + UV (A is the Slater matrix)
	//WinvU is Winv * U; the i and j columns of Winv where i and j are the sites to update
	//RatioDiff is V A^-1 (which can be computed easily from rows of Winv)
	//WoodburyDiff is (I_k + V A^-1 U)^-1 (2x2) times RatioDiff (2xN)

	ComplexDoubleMatrix<MKL_Complex16> U(2, N);
	U(colk, 0) = MKL_Complex16({1.0, 0.0});
	U(coll, 1) = MKL_Complex16({1.0, 0.0});
	ComplexDoubleMatrix<MKL_Complex16> WinvU = Winv * U;

	ComplexDoubleMatrix<MKL_Complex16> RatioDiff(2, N);
	RatioDiff.copy_row(RatioDiff, WinvU, rowk, 0);
	RatioDiff.copy_row(RatioDiff, WinvU, rowl, 0);
	RatioDiff(0, colk) = RatioDiff(0, colk) - MKL_Complex16({1.0, 0.0});
	RatioDiff(1, coll) = RatioDiff(1, coll) - MKL_Complex16({1.0, 0.0});

	ComplexDoubleMatrix<MKL_Complex16> Woodbury = ComplexDoubleMatrix<MKL_Complex16>::identity(2) - RatioDiff * U;
	Woodbury = Woodbury.compute_inverse();

	Winv = WinvU * Woodbury * RatioDiff;
}

/// PRIVATE MATRIX ELEMENTS

// swap 2 sites with specified sz values, with jastrow
MKL_Complex16 ProjectedState::psi_over_psi2(int site1, int site2, int new_sz1, int new_sz2) {
	MKL_Complex16 result;
	//new_sz1 = configuration[site2], new_sz2 = configuration[site1];

	int spin_row1 = Spin_t_to_row(new_sz1), spin_row2 = Spin_t_to_row(new_sz2);

	result = Winv((site1 + spin_row1), parton_labels[site1]) * Winv((site2 + spin_row2),  parton_labels[site2])
				- Winv((site1 + spin_row1), parton_labels[site2]) * Winv((site2 + spin_row2), parton_labels[site1]);

	std::vector<int> flip_sites = { site1, site2 }, new_sz = { new_sz1, new_sz2 };
	return result * jastrow.lazy_eval(flip_sites, new_sz, configuration);
}

// 3-site ring exchange, with jastrow
MKL_Complex16 ProjectedState::psi_over_psi_swap(int site1, int site2, int site3) {
	MKL_Complex16 result;
	int new_sz1 = configuration[site2], new_sz2 = configuration[site3], new_sz3 = configuration[site1];
	int spin_row1 = Spin_t_to_row(new_sz1), spin_row2 = Spin_t_to_row(new_sz2), spin_row3 = Spin_t_to_row(new_sz3);

	MKL_Complex16 
		row1_1 = Winv((site1 + spin_row1), parton_labels[site1]),
		row1_2 = Winv((site1 + spin_row1), parton_labels[site2]),
		row1_3 = Winv((site1 + spin_row1), parton_labels[site3]);

	result = row1_1 * (Winv((site2 + spin_row2), parton_labels[site2]) * Winv((site3 + spin_row3), parton_labels[site3])
					 - Winv((site3 + spin_row3), parton_labels[site2]) * Winv((site2 + spin_row2), parton_labels[site3]))
		- row1_2 * (Winv((site2 + spin_row2), parton_labels[site1]) * Winv((site3 + spin_row3), parton_labels[site3])
				- Winv((site3 + spin_row3),  parton_labels[site1]) * Winv((site2 + spin_row2), parton_labels[site3]))
		+ row1_3 * (Winv((site2 + spin_row2), parton_labels[site1]) * Winv((site3 + spin_row3), parton_labels[site2])
				- Winv((site3 + spin_row3), parton_labels[site1]) * Winv((site2 + spin_row2), parton_labels[site2]));

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

	MKL_Complex16 result({1.0, 0.0});
	if (flips.size() == 2) {
		result *= psi_over_psi2(flips[0], flips[1], new_sz[0], new_sz[1]);
	}
	else if (flips.size() == 3) {
		result *= psi_over_psi_swap(flips[0], flips[1], flips[2]);
	}
	else {
		assert(flips.size() == 2);
	}
	return to_std_complex(result);

}

/// PRIVATE UPDATES - ACTUAL CALLERS

void ProjectedState::update(int site1, int site2) {
	//only for swapping spins at two different sites

	if (configuration[site1] < configuration[site2]) {
		updateMatrixInverse((site1 + Spin_t_to_row(configuration[site2])), parton_labels[site2], (site2 + Spin_t_to_row(configuration[site1])), parton_labels[site1]);
	}
	else {
		updateMatrixInverse((site2 + Spin_t_to_row(configuration[site1])), parton_labels[site1], (site1 + Spin_t_to_row(configuration[site2])), parton_labels[site2]);
	}
	//std::cout << "labels: " << vec2str(parton_labels) << "\n";
	int templabel = parton_labels[site1];
	parton_labels[site1] = parton_labels[site2];
	parton_labels[site2] = templabel;
	templabel = configuration[site1];
	configuration[site1] = configuration[site2];
	configuration[site2] = templabel;
}

void ProjectedState::update_(std::vector<int>& sites, std::vector<int>& new_sz) {

	assert(sites.size() == 2);
	assert(new_sz.size() == 2);

	if (configuration[sites[0]] < configuration[sites[1]]) {
		updateMatrixInverse((sites[0] + Spin_t_to_row(new_sz[0])), parton_labels[sites[0]], (sites[1] + Spin_t_to_row(new_sz[1])), parton_labels[sites[1]]);
	}
	else {
		updateMatrixInverse((sites[1] + Spin_t_to_row(new_sz[1])), parton_labels[sites[1]], (sites[0] + Spin_t_to_row(new_sz[0])), parton_labels[sites[0]]);
	}

	configuration[sites[0]] = new_sz[0];
	configuration[sites[1]] = new_sz[1];
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

	if (flips.size() == 2) {
		jastrow.update_tables(flips, new_sz, configuration);
		if (configuration[flips[0]] == new_sz[1] && configuration[flips[1]] == new_sz[0]) {
			update(flips[0], flips[1]);
		}
		else {
			update_(flips, new_sz);
		}
	}
	else if (flips.size() == 3) {
		assert(!jastrow.exist()); //jastrow not implemented for ring exchanges
		std::vector<int> fliplist(2), spinlist(2);
		fliplist = { flips[0], flips[1] };
		spinlist = { new_sz[0], new_sz[2] };
		update(flips[0], flips[1]);

		fliplist = { flips[1], flips[2] };
		spinlist = { new_sz[1], new_sz[2] };
		update(flips[1], flips[2]);
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

#include "Wavefunction.h"
#include "vmc_io.h"







//std::complex<double> DirectorProductState::psi_over_psi(std::vector<int>& state, std::vector<int>& flips, std::vector<int>& new_sz) {
//
//	assert(flips.size() == 2);//only flip two spins in this iteration of the code
//
//	//more than 2 spin flips should still work in this part
//	std::complex<double> result(1.0, 0.0);
//
//	for (int i = 0; i < flips.size(); ++i) {
//		result *= basis_element(flips[i], new_sz[i]) / basis_element(flips[i], state[flips[i]]);
//	}
//
//	return result;
//
//}
//
//std::vector<double> DirectorProductState::expectation_Nz() {
//	std::vector<double> result(lat.get_N());
//
//	for (int i = 0; i < result.size(); ++i) {
//		result[i] = std::norm(basis_element(i, 0));
//	}
//
//	return result;
//}
//
//double DirectorProductState::expectation_Nz_total() {
//
//	std::vector<double> per_site = expectation_Nz();
//
//	return std::accumulate(per_site.begin(), per_site.end(), 0.0);
//}
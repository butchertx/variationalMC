#include "Lattice.h"
#include "vmc_io.h"
#include "SpinModel.h"

SpinModel create_su2_Hamiltonian(Lattice l, std::vector<double> J, vmctype::Spin_t spin_val) {
    SpinModel ham;
    std::cout << "Creating Heisenberg Hamiltonian with J = [" << vec2str(J) << "]\n";

    //2-site exchanges
    std::vector<int> nei;
    Observable *s12;
    std::string obs_name;
    std::complex<double> E0 = { 0.0, 0.0 };
    for (int range = 0; range < J.size(); ++range) {
        //Bilinear terms
        obs_name = "<S_i . S_j>_" + std::to_string(range);
        s12 = new Observable(obs_name);
        for (int i = 0; i < l.get_N(); ++i) {
            nei = l.get_neighbors(i, range);
            for (int n = 0; n < nei.size(); ++n) {
                HeisenbergExchange* I0;
                if (spin_val == vmctype::Spin_t::HALF){
                    I0 = new HeisenbergExchange(i, nei[n], 0.5, 0.5);
                }
                else if (spin_val == vmctype::Spin_t::ONE){
                    I0 = new HeisenbergExchange(i, nei[n], 1.0, 0.5);
                }
                s12->add_interaction(I0);
            }
        }
        ham.add_term(obs_name, *s12, { J[range], 0.0 });
    }
    ham.add_constant(E0);

    return ham;
}
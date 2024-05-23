/**
 * Example of a driver program, which can be used to run the pre-made example problems.
*/

#include <complex>
#include "vmctype.h"
#include "vmc_io.h"
#include "MemTimeTester.h"
#include "Lattice.h"
#include "MeanFieldAnsatz.h"
#include "RandomEngine.h"
#include "Wavefunction.h"
#include "ProjectedState.h"
#include "SpinModel.h"
#include "VariationalMonteCarlo.h"
// #include "model_and_calculation_helper.h"

MemTimeTester timer;
SpinModel create_su2_Hamiltonian(Lattice l, std::vector<double> J, vmctype::Spin_t spin_val);

struct results_struct {
    std::complex<double> E;
    std::complex<double> E_err;
    std::map<std::string, std::complex<double>> observables;
    std::map<std::string, std::complex<double>> observables_err;
};

int main(int argc, char* argv[]) {
    timer.flag_start_time("Total Program Time");
    std::string example_path, outfile_name;
    if (argc != 2) {
        std::cout << "Mandatory command line argument: <example directory>\n";
        std::cout << "Exiting...";
        return 1;
    }

    example_path = argv[1];
    std::cout << example_path << "\n";

    // Read in options
    LatticeOptions lat_options = read_json_lattice_from_dir(example_path);
    WavefunctionOptions wf_options = read_json_wavefunction_from_dir(example_path);
    ModelOptions mdl_options = read_json_model_from_dir(example_path);
    VMCOptions mc_options = read_json_vmc_from_dir(example_path);

    // Create objects
    Lattice lattice(Lattice_type_from_string(lat_options.type), vec3<int>(lat_options.L), vec3<int>(lat_options.pbc));
    std::shared_ptr<MeanFieldAnsatz> mf_ansatz;
    if (wf_options.other_options.spin == vmctype::Spin_t::HALF){
        mf_ansatz = std::shared_ptr<MeanFieldAnsatz>(new MeanFieldAnsatz_HALF(wf_options, lattice));
    }
    else if (wf_options.other_options.spin == vmctype::Spin_t::ONE){
        mf_ansatz = std::shared_ptr<MeanFieldAnsatz>(new MeanFieldAnsatz_ONE(wf_options, lattice));
    }
    else {
        std::stringstream ss;
        ss << "Spin value " << wf_options.other_options.spin << " not implemented\n";
        throw vmctype::NotImplemented(ss.str());
    }
    mf_ansatz->print_levels(true);

    RandomEngine r(-1, lattice.get_N(), lattice.get_neighbor_counts()[0]);
    ProjectedState wf(*mf_ansatz, r);
    std::vector<int> flips({0, 1});
    wf.update(flips);
    wf.print_matrix("Winv");

    SpinModel Ham = create_su2_Hamiltonian(lattice, mdl_options.get_su2_terms("J"), wf_options.other_options.spin);
    MonteCarloEngine sampler(Ham, wf, lattice, r, mc_options);
    sampler.run();

    results_struct results;
    results.E = sampler.get_energy() * (1.0/lattice.get_N());
    results.E_err = sampler.get_energy_err() * (1.0 / lattice.get_N());

    std::cout << "Results: E = " << results.E << " +- " << results.E_err << "\n";
    std::cout << "\n\n";
    timer.flag_end_time("Total Program Time");
    timer.print_timers();

    return 0;

}

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
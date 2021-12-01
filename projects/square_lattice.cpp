#include <complex>
#include "../src/vmctype.h"
#include "../src/vmc_io.h"
#include "../src/MemTimeTester.h"
#include "../src/Lattice.h"
#include "../src/MeanFieldAnsatz.h"
#include "../src/RandomEngine.h"
#include "../src/Wavefunction.h"
#include "../src/ProjectedState.h"
#include "../src/SpinModel.h"
#include "../src/VariationalMonteCarlo.h"

using namespace vmctype;

SpinModel create_Hamiltonian(Lattice, double J, double K);

struct results_struct {
    std::complex<double> E;
    std::complex<double> E_err;
    std::complex<double> J;
    std::complex<double> J_err;
    std::complex<double> K;
    std::complex<double> K_err;
};

results_struct run_mc(Lattice, mean_field_options, double, double);

MemTimeTester timer;

int main(int argc, char* argv[]) {
    timer.flag_start_time("Total Program Time");
    std::string infile_name, outfile_name;
    if (argc == 1) {
        std::cout << "Mandatory command line argument: <input_filename>\n";
        std::cout << "Optional command line argument: <results label> to get output file results/<results label>.csv\n";
        std::cout << "Exiting...";
        return 0;
    }
    else if (argc == 2) {
        infile_name = argv[1];
        outfile_name = "results/observables.csv";
        std::cout << "Using input file: " << infile_name << "\n";
        std::cout << "Using default output file: results/observables.csv\n";
    }
    else {
        infile_name = argv[1];
        std::string outfile_label = argv[2];
        outfile_name = "results/" + outfile_label + ".csv";
        std::cout << "Using input file: " << infile_name << "\n";
        std::cout << "Using output file: " << outfile_name << "\n";
    }

    lattice_options lat_options;
    mean_field_options wf_options;
    model_options mdl_options;
    vmc_options mc_options;

    read_json_full_input(&lat_options, &wf_options, &mdl_options, &mc_options, infile_name);

    //lattice_options lat_options = read_json_lattice(infile_name);
    Lattice lattice(Lattice_type_from_string(lat_options.type), vec3<int>(lat_options.L), vec3<int>(lat_options.pbc));

    makePath("./data");
    std::ofstream neighborfile;
    neighborfile.open("data/neighbors.txt");
    lattice.print_neighbors(&neighborfile);
    neighborfile.close();

    std::ofstream ringfile;
    ringfile.open("data/rings.txt");
    lattice.print_rings(&ringfile);
    ringfile.close();
    lattice.print_timers();

    //mean_field_options wf_options = read_json_wavefunction(infile_name);

    results_struct results;
    results = run_mc(lattice, wf_options, mdl_options.bilinear_terms[0].coupling, mdl_options.bilinear_terms[1].coupling);
    makePath("./results");
    std::ofstream results_file;
    results_file.open(outfile_name);
    results_file << ", E, E_err, <S_i * S_j>, <S_i * S_j>_err, <(S_i * S_j)^2>, <(S_i * S_j)^2>_err\n";
    results_file << "real, " << results.E.real() << ", " << results.E_err.real() << ", " << results.J.real() << ", " << results.J_err.real()
        << ", " << results.K.real() << ", " << results.K_err.real() << "\n";
    results_file << "imag, " << results.E.imag() << ", " << results.E_err.imag() << ", " << results.J.imag() << ", " << results.J_err.imag()
        << ", " << results.K.imag() << ", " << results.K_err.imag() << "\n";

    timer.flag_end_time("Total Program Time");
    timer.print_timers();
    return 1;

}

results_struct run_mc(Lattice lattice, mean_field_options mf_options, double J, double K) {
    MeanFieldAnsatz mf(mf_options, lattice, true);
    mf.print_levels();
    mf.print_fermi_level();
    std::ofstream director_file;
    director_file.open("data/directors.csv");
    mf.print_directors(&director_file);
    director_file.close();

    RandomEngine r(-1, lattice.get_N(), lattice.get_neighbor_counts()[0]);
    ProjectedState wf(mf, r);
    results_struct results;

    if (std::abs(wf.get_det()) == 0) {
        std::cout << "No valid initializations\n";
    }
    else {

        SpinModel Ham = create_Hamiltonian(lattice, J, K);
        MonteCarloEngine sampler(Ham, wf, lattice, r, VMCParams(10, 1000, 0));
        sampler.run();

        results.E = sampler.get_energy();
        results.E_err = sampler.get_energy_err();
        results.J = sampler.get_observable("Bilinear");
        results.J_err = sampler.get_observable_err("Bilinear");
        results.K = 1.0 + sampler.get_observable("SU(3) Exchange") - sampler.get_observable("Bilinear");
        results.K_err = sampler.get_observable_err("SU(3) Exchange") + sampler.get_observable_err("Bilinear");
        sampler.print_timers();
        std::cout << "Results: E = " << results.E << " +- " << results.E_err << "\n";
    }
    return results;
}


SpinModel create_Hamiltonian(Lattice l, double J, double K) {
    SpinModel ham;
    std::cout << "Creating BLBQ Hamiltonian with J = " << J << ", K = " << K << "\n";

    //2-site exchanges
    std::vector<int> nei;
    Observable s12("Bilinear");
    for (int i = 0; i < l.get_N(); ++i) {
        nei = l.get_neighbors(i, 0);
        for (int n = 0; n < nei.size(); ++n) {
            auto* I0 = new HeisenbergExchange(i, nei[n], 1.0, 0.5);
            s12.add_interaction(I0);
        }
    }
    ham.add_term("Bilinear", s12, { J - K, 0.0 });

    //2-site exchanges
    Observable p12("SU(3) Exchange");
    std::complex<double> E0 = { 0.0, 0.0 };
    for (int i = 0; i < l.get_N(); ++i) {
        nei = l.get_neighbors(i, 0);
        for (int n = 0; n < nei.size(); ++n) {
            auto* I0 = new SwapExchange(i, nei[n], 0.5);
            p12.add_interaction(I0);
            E0 += { 0.5*K, 0.0 };
        }
    }
    ham.add_term("SU(3) Exchange", p12, { K, 0.0 });
    ham.add_constant(E0);

    return ham;
}
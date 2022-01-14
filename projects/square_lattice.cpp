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

std::vector<Observable> create_Correlation_SZ(Lattice);

std::vector<Observable> create_Correlation_ladder(Lattice);

JastrowTable create_Jastrow(Lattice, JastrowTableOptions);

struct results_struct {
    std::complex<double> E;
    std::complex<double> E_err;
    std::complex<double> J;
    std::complex<double> J_err;
    std::complex<double> K;
    std::complex<double> K_err;
};

results_struct run_mc(lattice_options, mean_field_options, model_options, vmc_options);

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

    results_struct results = run_mc(lat_options, wf_options, mdl_options, mc_options);
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

results_struct run_mc(lattice_options lat_options, mean_field_options mf_options, model_options mdl_options, vmc_options v_options) {

    Lattice lattice(Lattice_type_from_string(lat_options.type), vec3<int>(lat_options.L), vec3<int>(lat_options.pbc));

    makePath("./data");
    std::ofstream neighborfile;
    neighborfile.open("data/neighbors.txt");
    lattice.print_neighbors(&neighborfile);
    neighborfile.close();

    MeanFieldAnsatz mf(mf_options, lattice, true);
    mf.print_levels();
    mf.print_fermi_level();

    std::ofstream director_file;
    director_file.open("data/directors.csv");
    mf.print_directors(&director_file);
    director_file.close();

    RandomEngine r(-1, lattice.get_N(), lattice.get_neighbor_counts()[0]);
    ProjectedState wf(mf, r, create_Jastrow(lattice, mf_options.jastrow));
    results_struct results;

    if (std::abs(wf.get_det()) == 0) {
        std::cout << "No valid initializations\n";
    }
    else {

        SpinModel Ham = create_Hamiltonian(lattice, mdl_options.bilinear_terms[0].coupling, mdl_options.bilinear_terms[1].coupling);
        MonteCarloEngine sampler(Ham, wf, lattice, r, v_options);
        sampler.add_observable_function(create_Correlation_SZ(lattice), "SzSz_correlation");
        sampler.add_observable_function(create_Correlation_ladder(lattice), "Spm_correlation");
        sampler.run();

        std::ofstream obs_file_r, obs_file_i;
        obs_file_r.open("results/SzSz_correlation_real.csv");
        obs_file_i.open("results/SzSz_correlation_imag.csv");
        sampler.write_observable_functions(&obs_file_r, &obs_file_i, "SzSz_correlation");
        obs_file_r.close();
        obs_file_i.close();

        obs_file_r.open("results/Spm_correlation_real.csv");
        obs_file_i.open("results/Spm_correlation_imag.csv");
        sampler.write_observable_functions(&obs_file_r, &obs_file_i, "Spm_correlation");
        obs_file_r.close();
        obs_file_i.close();
        

        results.E = sampler.get_energy();
        results.E_err = sampler.get_energy_err();
        results.J = sampler.get_observable("Bilinear");
        results.J_err = sampler.get_observable_err("Bilinear");
        results.K = 1.0 + sampler.get_observable("SU(3) Exchange") - sampler.get_observable("Bilinear");
        results.K_err = sampler.get_observable_err("SU(3) Exchange") + sampler.get_observable_err("Bilinear");
        std::cout << "Results: E = " << results.E << " +- " << results.E_err << "\n";
        std::cout << "\n\n";
        sampler.print_timers();
        
    }
    wf.print_timers();
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

std::vector<Observable> create_Correlation_SZ(Lattice l) {
    std::vector<Observable> correlation;

    for (int i = 0; i < l.get_N(); ++i) {
        auto* s12 = new Observable("Szi Szj");
        auto* I0 = new IsingExchange(0, i, 1.0);
        s12->add_interaction(I0);
        correlation.push_back(*s12);
    }
    
    return correlation;
}

std::vector<Observable> create_Correlation_ladder(Lattice l) {
    std::vector<Observable> correlation;

    for (int i = 0; i < l.get_N(); ++i) {
        auto* s12 = new Observable("Si+- Sj-+");
        auto* I0 = new LadderExchange(0, i, 1.0);
        s12->add_interaction(I0);
        correlation.push_back(*s12);
    }

    return correlation;
}

JastrowTable create_Jastrow(Lattice lattice, JastrowTableOptions jopt) {
    JastrowFactorOptions j(jopt.sz);

    //assume isotropic=true
    
    std::vector<double> params;
    if (j.distance_max == j.values.size()) {
        params = j.values;
    } 
    else{
        std::cout << "Initializing Jastrow with all zeros\n";
        params = std::vector<double>(j.distance_max, 0.0);
    }

    std::vector<JastrowFactor> jlist;
    std::vector<std::vector<int>> neighbors;
    for (int i = 0; i < params.size(); ++i) {
        neighbors = lattice.get_neighbors(i);
        jlist.push_back(JastrowFactor(params[i], neighbors));
    }

    return JastrowTable(jlist);
}
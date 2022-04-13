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
#include "../src/model_and_calculation_helper.h"

// Contains default behavior for accepting input files and writing/choosing output files

// lattice and model selection are included in input parameters and offloaded to "model_and_calculation_helper.h"

// Default output files:
//// data/neighbors.csv: List of sites and neighbors in the lattice
//// data/rings.csv: List of ring-exchange groupings in the lattice (if ring exchange is implemented for the chosen lattice/model)
// data/wavefunction.json: Wavefunction parameters that can be read back in for subsequent calculations
//// data/directors.csv: directors for ordered wavefunctions
//// data/mean_field_energies.csv: Info about the non-interacting orbitals and Fermi surface
//// data/conf_init.csv: initial configuration
//// data/conf_final.csv: final configuration
// data/markov.csv: list of markov chain updates leading to the final configuration
//// results/observables.csv: expectation values and errors for observables calculated in the Hamiltonian
//// results/<output_file>: alternate name for observables.csv that can be specified by command line argument
//// results/optimization.csv: Energy bins and variational parameter values in an optimization calculation
//// results/<obs>_correlation_<imag/real>.csv: Correlation function connecting <obs> at sites i and j, real and imaginary part

MemTimeTester timer;

void output_results(std::string outfile_name, results_struct res);

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
        outfile_name = "results/observables";
        std::cout << "Using input file: " << infile_name << "\n";
        std::cout << "Using default output file: results/observables.csv\n";
    }
    else {
        infile_name = argv[1];
        std::string outfile_label = argv[2];
        outfile_name = "results/" + outfile_label;
        std::cout << "Using input file: " << infile_name << "\n";
        std::cout << "Using output file: " << outfile_name << "\n";
    }

    lattice_options lat_options;
    mean_field_options wf_options;
    model_options mdl_options;
    vmc_options mc_options;

    read_json_full_input(&lat_options, &wf_options, &mdl_options, &mc_options, infile_name);
    makePath("./results");

    if (mc_options.search) {
        std::vector<results_struct> results_list;

        //unroll all psearch vectors
        //assume mu and g are the only valid params
        std::vector<double> psearch_mu, psearch_g;
        assert(mc_options.psearch.num_parallel == mc_options.psearch.mu_z.size() * mc_options.psearch.g.size());
        assert(mc_options.psearch.txy.size() == 0);
        assert(mc_options.psearch.tz.size() == 0);
        for (int mdx = 0; mdx < mc_options.psearch.mu_z.size(); ++mdx) {
            for (int gdx = 0; gdx < mc_options.psearch.g.size(); ++gdx) {
                psearch_mu.push_back(mc_options.psearch.mu_z[mdx]);
                psearch_g.push_back(mc_options.psearch.g[gdx]);
            }
        }

        for (int i = 0; i < mc_options.psearch.num_parallel; ++i) {
            wf_options.mu_z = psearch_mu[i];
            wf_options.field = psearch_g[i];
            try {
                results_list.push_back(run_mc(lat_options, wf_options, mdl_options, mc_options));
            }
            catch (...) {
                results_struct dummy;
                dummy.E = std::complex<double>(1000.0, 0.0);
                results_list.push_back(dummy);
                // move on to the next simulation
            }
        }

        std::vector<double> energy_vals, n0_vals;
        int opt_res_idx = 0;
        double e_min = results_list[0].E.real();
        std::ofstream param_file;
        param_file.open(outfile_name + "_params.csv");
        param_file << "idx, E, E_err, mu_z, g" << "\n";

        for (int i = 0; i < results_list.size(); ++i) {
            energy_vals.push_back(results_list[i].E.real());
            n0_vals.push_back(1.0 - results_list[i].observables["<(S_i^z)^2>"].real());

            output_results(outfile_name + std::to_string(i) + ".csv", results_list[i]);
            if (results_list[i].E.real() < e_min) {
                opt_res_idx = i;
                e_min = results_list[i].E.real();
            }

            param_file << i << ", " << results_list[i].E.real() << ", "  << results_list[i].E_err.real()  << ", " << psearch_mu[i] << ", " << psearch_g[i] << "\n";

        }
        param_file.close();

        output_results(outfile_name + "_opt.csv", results_list[opt_res_idx]);
        std::ofstream results_file;
        results_file.open(outfile_name + "_opt.csv", std::ios_base::app);
        results_file << "mu_z_opt," << psearch_mu[opt_res_idx] << "\n";
        results_file << "g_opt," << psearch_g[opt_res_idx];
        results_file.close();

    }
    else {
        results_struct results = run_mc(lat_options, wf_options, mdl_options, mc_options);    
        output_results(outfile_name + ".csv", results);        
    }

    timer.flag_end_time("Total Program Time");
    timer.print_timers();

    return 1;

}

void output_results(std::string outfile_name, results_struct results) {
    std::ofstream results_file;
    results_file.open(outfile_name);
    std::vector<std::string> obs_names;
    std::vector<double> obs_real, obs_real_err, obs_imag, obs_imag_err;
    for (auto item : results.observables) {
        obs_names.push_back(item.first);
        obs_real.push_back(item.second.real());
        obs_real_err.push_back(results.observables_err[item.first].real());
        obs_imag.push_back(item.second.imag());
        obs_imag_err.push_back(results.observables_err[item.first].imag());
    }
    results_file << ", E, " << vec2str(obs_names) << "\n";
    results_file << "real, " << results.E.real() << ", " << vec2str(obs_real) << "\n";
    results_file << "real_err, " << results.E_err.real() << ", " << vec2str(obs_real_err) << "\n";
    results_file << "imag, " << results.E.imag() << ", " << vec2str(obs_imag) << "\n";
    results_file << "imag_err, " << results.E_err.imag() << ", " << vec2str(obs_imag_err) << "\n";
    results_file.close();
}
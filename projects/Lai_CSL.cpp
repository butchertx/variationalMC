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

SpinModel create_Hamiltonian(Lattice, double J, double K, double h);

struct results_struct {
    std::complex<double> E;
    std::complex<double> E_err;
    std::complex<double> J;
    std::complex<double> J_err;
    std::complex<double> K;
    std::complex<double> K_err;
    std::complex<double> h;
    std::complex<double> h_err;
};

results_struct run_mc(Lattice, mean_field_options, double, double, double);

MemTimeTester timer;

//Inputs:
//read vmc options
//read ansatz params
//read model params

int main(int argc, char* argv[]) {
    timer.flag_start_time("Total Program Time");
    std::string infile_name, outfile_name;
    if(argc == 1){
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
    int steps = 10, measures = 10000, throwaway = 100;

    lattice_options lat_options = read_json_lattice(infile_name);
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

    mean_field_options wf_options = read_json_wavefunction(infile_name);

    results_struct results;
    results = run_mc(lattice, wf_options, 1.0, 0.0, 0.0);
    makePath("./results");
    std::ofstream results_file;
    results_file.open(outfile_name);
    results_file << ", E, E_err, <P_12>, <P_12>_err, <P_sym>, <P_sym>_err, <chi>, <chi>_err\n";
    results_file << "real, " << results.E.real() << ", " << results.E_err.real() << ", " << results.J.real() << ", " << results.J_err.real()
        << ", " << results.K.real() << ", " << results.K_err.real() << ", " << results.h.real() << ", " << results.h_err.real() << "\n";
    results_file << "imag, " << results.E.imag() << ", " << results.E_err.imag() << ", " << results.J.imag() << ", " << results.J_err.imag()
        << ", " << results.K.imag() << ", " << results.K_err.imag() << ", " << results.h.imag() << ", " << results.h_err.imag() << "\n";
    //std::cout << "E = " << results.E << " +- " << results.err << "\n";

    //for (int sim = 0; sim < num_t2; ++sim) {
    //    t2 = t2min + sim * t2change;
    //    wf_options.hopping_list[2].strength = t2;
    //    wf_options.hopping_list[3].strength = t2;
    //    results = run_mc(lattice, wf_options, p, 1.0, 0.0);
    //    t2_list.push_back(t2);
    //    k_list.push_back(1.0);
    //    h_list.push_back(0.0);
    //    e_list.push_back(results.E);
    //    err_list.push_back(results.err);

    //    results = run_mc(lattice, wf_options, p, 0.0, 1.0);
    //    t2_list.push_back(t2);
    //    k_list.push_back(0.0);
    //    h_list.push_back(1.0);
    //    e_list.push_back(results.E);
    //    err_list.push_back(results.err);
    //}

    //std::ofstream file;
    //file << "K, h, t2, E_real, E_imag, err_real, err_imag\n";
    //for (int i = 0; i < 2*num_t2; ++i) {
    //    file << k_list[i] << ", " << h_list[i] << ", " << t2_list[i] << ", " << e_list[i].real() << ", " << e_list[i].imag() << ", " << err_list[i].real() << ", " << err_list[i].imag() << "\n";
    //}
    //file.close();
    timer.flag_end_time("Total Program Time");
    timer.print_timers();
    return 1;

}

results_struct run_mc(Lattice lattice, mean_field_options mf_options, double J, double K, double h) {
    MeanFieldAnsatz mf(mf_options, lattice, true);
    mf.print_levels();
    mf.print_fermi_level();

    RandomEngine r(-1, lattice.get_N(), lattice.get_neighbor_counts()[0]);
    ProjectedState wf(mf, r);
    results_struct results;

    if (std::abs(wf.get_det()) == 0) {
        std::cout << "No valid initializations\n";
    }
    else {

        SpinModel Ham = create_Hamiltonian(lattice, J, K, h);
        MonteCarloEngine sampler(Ham, wf, lattice, r, VMCParams(10, 1000, 0));
        sampler.run();

        results.E = sampler.get_energy();
        results.E_err = sampler.get_energy_err();
        results.J = sampler.get_observable("J");
        results.J_err = sampler.get_observable_err("J");
        results.K = sampler.get_observable("K"), 
        results.K_err = sampler.get_observable_err("K"), 
        results.h = sampler.get_observable("h"),
        results.h_err = sampler.get_observable_err("h");
        sampler.print_timers();
        std::cout << "Results: E = " << results.E << " +- " << results.E_err << "\n";
    }
    return results;
}


SpinModel create_Hamiltonian(Lattice l, double J, double K, double h) {
    SpinModel ham;
    std::cout << "Creating SU(3) Hamiltonian with J = " << J << ", K = " << K << ", h = " << h << "\n";

    //2-site exchanges
    std::vector<int> nei;
    Observable p12("J");
    for (int i = 0; i < l.get_N(); ++i) {
        nei = l.get_neighbors(i, 0);
        for (int n = 0; n < nei.size(); ++n) {
            auto* I0 = new SwapExchange(i, nei[n], 0.5);
            p12.add_interaction(I0);
        }
    }
    ham.add_term("J", p12, { J, 0.0 });

    //all ring exchanges
    std::vector<int> ring_sites;
    RingList& r = l.ring_ref();
    Observable p123_left("K"), p123_right("h");
    for (int tri = 0; tri < r.get_size(); tri += 1) {
        auto* I1 = new RingExchange(r, tri, true, 1.0);
        p123_left.add_interaction(I1);
        I1 = new RingExchange(r, tri, false, 1.0);
        p123_left.add_interaction(I1);

        auto* I2 = new RingExchange(r, tri, true, 1.0);
        p123_right.add_interaction(I2);
        I2 = new RingExchange(r, tri, false, -1.0);
        p123_right.add_interaction(I2);
    }
    ham.add_term("K", p123_left, { K, 0.0 });
    ham.add_term("h", p123_right, { 0.0, h });

    return ham;
}

#include <assert.h>

#include "vmc_io.h"
#include "Lattice.h"
#include "MeanFieldAnsatz.h"
#include "Wavefunction.h"
#include "ProjectedState.h"
#include "SpinModel.h"
#include "VariationalMonteCarlo.h"
#include "model_and_calculation_helper.h"

results_struct run_mc(LatticeOptions lat_options, MeanFieldOptions mf_options, ModelOptions mdl_options, VMCOptions v_options) {

    Lattice lattice(Lattice_type_from_string(lat_options.type), vec3<int>(lat_options.L), vec3<int>(lat_options.pbc));

    makePath("./data");
    std::ofstream latticefile;
    latticefile.open("data/neighbors.csv");
    lattice.write_neighbors(&latticefile);
    latticefile.close();

    if (lattice.has_rings()) {
        latticefile.open("data/rings.csv");
        lattice.write_rings(&latticefile);
        latticefile.close();
    }

    MeanFieldAnsatz mf(mf_options, lattice, true);

    std::ofstream wavefunctionfile;
    wavefunctionfile.open("data/mean_field_energies.csv");
    mf.write_levels(&wavefunctionfile);
    wavefunctionfile.close();
    // mf.print_fermi_level();

    wavefunctionfile.open("data/directors.csv");
    mf.write_directors(&wavefunctionfile);
    wavefunctionfile.close();

    RandomEngine r(-1, lattice.get_N(), lattice.get_neighbor_counts()[0]);
    ProjectedState wf(mf, r, create_Jastrow(lattice, mf_options.jastrow));

    std::ofstream conf_file;
    conf_file.open("data/conf_init.csv");
    wf.write_configuration(&conf_file);
    conf_file.close();

    results_struct results;

    if (std::abs(wf.get_det()) == 0) {
        std::cout << "No valid initializations\n";
    }
    else {

        makePath("./results");
        SpinModel Ham;
        if (mdl_options.model_type.compare("blbq") == 0) {
            Ham = create_su2_Hamiltonian(lattice, mdl_options.get_su2_terms("J"), mdl_options.get_su2_terms("K"), mdl_options.single_ion);
        }
        else if (mdl_options.model_type.compare("su3") == 0) {
            Ham = create_su3_Hamiltonian(lattice, mdl_options.bilinear_terms[0].coupling, mdl_options.ring3_terms[0].coupling_real, mdl_options.ring3_terms[0].coupling_imag, mdl_options.single_ion);
        }
        MonteCarloEngine sampler(Ham, wf, lattice, r, v_options);
        sampler.add_observable_function(create_Correlation_SZ(lattice), "SzSz_correlation");
        sampler.add_observable_function(create_Correlation_ladder(lattice), "Spm_correlation");
        sampler.add_observable_function(create_Correlation_Swap(lattice), "Swap_correlation");
        sampler.run();
        conf_file.open("data/conf_final.csv");
        wf.write_configuration(&conf_file);
        conf_file.close();

        std::ofstream obs_file_r, obs_file_i;
        //obs_file_r.open("results/SzSz_correlation_real.csv");
        //obs_file_i.open("results/SzSz_correlation_imag.csv");
        //sampler.write_observable_functions(&obs_file_r, &obs_file_i, "SzSz_correlation");
        //obs_file_r.close();
        //obs_file_i.close();

        //obs_file_r.open("results/Spm_correlation_real.csv");
        //obs_file_i.open("results/Spm_correlation_imag.csv");
        //sampler.write_observable_functions(&obs_file_r, &obs_file_i, "Spm_correlation");
        //obs_file_r.close();
        //obs_file_i.close();

        //obs_file_r.open("results/Swap_correlation_real.csv");
        //obs_file_i.open("results/Swap_correlation_imag.csv");
        //sampler.write_observable_functions(&obs_file_r, &obs_file_i, "Swap_correlation");
        //obs_file_r.close();
        //obs_file_i.close();


        results.E = sampler.get_energy() * (1.0/lattice.get_N());
        results.E_err = sampler.get_energy_err() * (1.0 / lattice.get_N());
        results.observables = sampler.get_all_observables();
        results.observables_err = sampler.get_all_observables_err();
        // normalize results by lattice size
        for (auto it = results.observables.begin(); it != results.observables.end(); ++it) {
            it->second *= (1.0 / lattice.get_N());
        }
        for (auto it = results.observables_err.begin(); it != results.observables_err.end(); ++it) {
            it->second *= (1.0 / lattice.get_N());
        }

        std::cout << "Results: E = " << results.E << " +- " << results.E_err << "\n";
        std::cout << "\n\n";

        sampler.print_timers();

    }
    wf.print_timers();
    return results;
}


SpinModel create_su2_Hamiltonian(Lattice l, std::vector<double> J, std::vector<double> K, double D) {
    SpinModel ham;
    std::cout << "Creating BLBQ Hamiltonian with J = [" << vec2str(J) << "], K = [" << vec2str(K) << "], D = " << D << "\n";
    if (J.size() != K.size()) {
        //need to broadcast so they are the same size
        if (J.size() < K.size()) {
            for (int jsize = J.size(); jsize < K.size(); ++jsize) {
                J.push_back(0.0);
            }
        }
        else {
            for (int ksize = K.size(); ksize < J.size(); ++ksize) {
                K.push_back(0.0);
            }
        }
    }
    assert(J.size() == K.size());

    //2-site exchanges
    std::vector<int> nei;
    Observable *s12, *p12;
    std::string obs_name;
    std::complex<double> E0 = { 0.0, 0.0 };
    for (int range = 0; range < J.size(); ++range) {
        //Bilinear terms
        obs_name = "<S_i . S_j>_" + std::to_string(range);
        s12 = new Observable(obs_name);
        for (int i = 0; i < l.get_N(); ++i) {
            nei = l.get_neighbors(i, range);
            for (int n = 0; n < nei.size(); ++n) {
                auto* I0 = new HeisenbergExchange(i, nei[n], 1.0, 0.5);
                s12->add_interaction(I0);
            }
        }
        ham.add_term(obs_name, *s12, { J[range] - K[range], 0.0 });

        //Biquadratic terms
        obs_name = "<S_i <-> S_j>_" + std::to_string(range);
        p12 = new Observable(obs_name);
        
        for (int i = 0; i < l.get_N(); ++i) {
            nei = l.get_neighbors(i, range);
            for (int n = 0; n < nei.size(); ++n) {
                auto* I0 = new SwapExchange(i, nei[n], 0.5);
                p12->add_interaction(I0);
                E0 += std::complex<double>(0.5 * K[range], 0.0);
            }
        }
        ham.add_term(obs_name, *p12, { K[range], 0.0 });
    }
    ham.add_constant(E0);

    

    Observable ds2("<(S_i^z)^2>");
    for (int i = 0; i < l.get_N(); ++i) {
        auto* I0 = new SingleIonAnisotropy(i, 1.0);
        ds2.add_interaction(I0);
    }
    ham.add_term("<(S_i^z)^2>", ds2, { D, 0.0 });

    return ham;
}

SpinModel create_su3_Hamiltonian(Lattice l, double J, double K, double h, double D) {
    SpinModel ham;
    std::cout << "Creating SU(3) Hamiltonian with J = " << J << ", K = " << K << ", h = " << h << "\n";

    //2-site exchanges
    std::vector<int> nei;
    Observable p12("<S_i <-> S_j>");
    for (int i = 0; i < l.get_N(); ++i) {
        nei = l.get_neighbors(i, 0);
        for (int n = 0; n < nei.size(); ++n) {
            auto* I0 = new SwapExchange(i, nei[n], 0.5);
            p12.add_interaction(I0);
        }
    }
    ham.add_term("<S_i <-> S_j>", p12, { J, 0.0 });

    //all ring exchanges
    assert(l.get_lattice_type() == TRIANGLE); // only implemented for triangle lattice
    std::vector<int> ring_sites;
    RingList& r = l.ring_ref();
    Observable p123_left("<2Re(C_3)>"), p123_right("<2Im(C_3)>");
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
    ham.add_term("<2Re(C_3)>", p123_left, { K, 0.0 });
    ham.add_term("<2Im(C_3)>", p123_right, { 0.0, h });

    if (D != 0.0) {
        Observable ds2("<(S_i^z)^2>");
        for (int i = 0; i < l.get_N(); ++i) {
            auto* I0 = new SingleIonAnisotropy(i, 1.0);
            ds2.add_interaction(I0);
        }
        ham.add_term("<(S_i^z)^2>", ds2, { D, 0.0 });
    }

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

std::vector<Observable> create_Correlation_Swap(Lattice l) {
    std::vector<Observable> correlation;

    for (int i = 0; i < l.get_N(); ++i) {
        auto* s12 = new Observable("Swap");
        auto* I0 = new SwapExchange(0, i, 1.0);
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
    //
    //  Sz options
    //
    JastrowFactorOptions j(jopt.sz);

    //assume isotropic=true

    std::vector<double> params;
    if (j.distance_max == j.values.size()) {
        params = j.values;
    }
    else {
        std::cout << "Initializing Jastrow Sz with all zeros\n";
        params = std::vector<double>(j.distance_max, 0.0);
    }

    std::vector<JastrowFactor> jlist;
    std::vector<std::vector<int>> neighbors;
    for (int i = 0; i < params.size(); ++i) {
        neighbors = lattice.get_neighbors(i);
        jlist.push_back(JastrowFactor(params[i], neighbors));
    }

    //
    //  Sz2 options
    //
    JastrowFactorOptions j2(jopt.sz2);

    //assume isotropic=true
    if (j2.distance_max > -1) {

        if (j2.distance_max == j2.values.size()) {
            params = j2.values;
        }
        else {
            std::cout << "Initializing Jastrow Sz2 with all zeros\n";
            params = std::vector<double>(j2.distance_max, 0.0);
        }

        for (int i = 0; i < params.size(); ++i) {
            neighbors = lattice.get_neighbors(i);
            jlist.push_back(JastrowFactor(params[i], neighbors, true));
        }
    }

    if (jopt.density_flag) {
        JastrowDensity jdens(jopt.density_coupling);
        std::cout << "Warning: Jastrow Density not implemented!\n";
        //return JastrowTable(jlist, jdens);
        return JastrowTable(jlist);
    }
    else {
        return JastrowTable(jlist);
    }

}

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
#include "VMCDriver.h"
#include "VMCResults.h"

MemTimeTester timer;

VMCResults VMCDriver::run() {
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

    RandomEngine r(-1, lattice.get_N(), lattice.get_neighbor_counts()[0]);
    ProjectedState wf(*mf_ansatz, r);
    wf.print_matrix("Phi");
    wf.print_matrix("Slater");
    wf.print_matrix("LU");
    wf.print_matrix("Winv");
    wf.print_matrix("UP1");
    wf.print_matrix("UP2");
    wf.print_matrix("UP3");
    wf.print_matrix("ipiv");

    SpinModel Ham = create_su2_Hamiltonian(lattice, mdl_options.get_su2_terms("J"), wf_options.other_options.spin);
    MonteCarloEngine sampler(Ham, wf, lattice, r, mc_options);
    sampler.run();


    VMCResults results(sampler.get_energy() * (1.0/lattice.get_N()),
                        sampler.get_energy_err() * (1.0 / lattice.get_N()),
                        sampler.get_all_observables(),
                        sampler.get_all_observables_err());

    return results;
}


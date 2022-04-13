#pragma once
#include <complex>
#include <vector>

#include "vmctype.h"

using namespace vmctype;

class SpinModel;
class Observable;
class Lattice;
class JastrowTable;

using namespace vmctype;

SpinModel create_su2_Hamiltonian(Lattice, std::vector<double> J, std::vector<double> K, double D);
SpinModel create_su3_Hamiltonian(Lattice, double J, double K, double h, double D = 0.0);

std::vector<Observable> create_Correlation_SZ(Lattice);

std::vector<Observable> create_Correlation_Swap(Lattice);

std::vector<Observable> create_Correlation_ladder(Lattice);

JastrowTable create_Jastrow(Lattice, JastrowTableOptions);

struct results_struct {
    std::complex<double> E;
    std::complex<double> E_err;
    std::map<std::string, std::complex<double>> observables;
    std::map<std::string, std::complex<double>> observables_err;
};

results_struct run_mc(lattice_options, mean_field_options, model_options, vmc_options);
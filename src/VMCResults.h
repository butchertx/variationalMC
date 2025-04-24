/**
 * @file VMCResults.h
 * @brief Header file for the VMCResults class.
 * @details This file contains the definition of the VMCResults class, which holds the data resulting from a VMC calculation
 */
#pragma once
#include <complex>
#include <map>
#include <string>
#include <iostream>
#include <vector>

class VMCResults {

    std::complex<double> E;
    std::complex<double> E_err;

    std::map<std::string, std::complex<double>> observables;
    std::map<std::string, std::complex<double>> observables_err;

public:

    VMCResults() : E(0.0), E_err(0.0) {};

    VMCResults(std::complex<double> E_in, std::complex<double> E_err_in,
                   std::map<std::string, std::complex<double>> observables_in,
                   std::map<std::string, std::complex<double>> observables_err_in) :
        E(E_in), E_err(E_err_in), observables(observables_in), observables_err(observables_err_in) {}

    // Getters
    std::complex<double> get_energy() const { return E; }
    std::complex<double> get_energy_err() const { return E_err; }
    std::map<std::string, std::complex<double>> get_observables() const { return observables; }
    std::map<std::string, std::complex<double>> get_observables_err() const { return observables_err; }
    
    // Setters
    void set_energy(std::complex<double> E_in) { E = E_in; }
    void set_energy_err(std::complex<double> E_err_in) { E_err = E_err_in; }
    void set_observables(std::map<std::string, std::complex<double>> observables_in) { observables = observables_in; }
    void set_observables_err(std::map<std::string, std::complex<double>> observables_err_in) { observables_err = observables_err_in; }

    void print_results() { 
        std::cout << "Results: E = " << E << " +- " << E_err << "\n";
        for (const auto& obs : observables) {
            std::cout << obs.first << " = " << obs.second << " +- " << observables_err[obs.first] << "\n";
        }
    }
};
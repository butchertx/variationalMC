#pragma once
#include <random>
#include <string>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <map>
#include "SpinModel.h"
#include "Wavefunction.h"
#include "RandomEngine.h"
#include <typeinfo>
#include <fstream>
#include <iostream>
#include "vmctype.h"
#include "vmc_io.h"

std::vector<double> vec_r(std::vector<std::complex<double>> cvec);

std::vector<double> vec_i(std::vector<std::complex<double>> cvec); 

class Lattice;

class MonteCarloEngine {

	MemTimeTester timer;

	SpinModel& H;
	std::map<std::string, Observable> observables; // separate from Hamiltonian observables
	std::vector<std::string> observable_names; // separate from Hamiltonian observables
	std::map<std::string, std::vector<Observable>> observable_functions;
	std::vector<std::string> observable_function_names;

	Wavefunction& WF;

	Lattice& lat;

	std::vector<double> E;
	std::vector<double> E_imag;

	std::vector<double> Nz;

	RandomEngine& rand;

	vmc_options params;

	std::map<std::string, std::vector<std::complex<double>>> observable_measures; // Hamiltonian observables are automatically added here
	std::map<std::string, std::vector<std::vector<std::complex<double>>>> observable_function_measures;

	std::vector<std::vector<double>> v_params_record;
	std::vector<double> E_bins_record;

	std::pair<std::vector<int>, std::vector<int>> step_su2();
	std::vector<int> step_su3(int num_swap);
	void measure_energy();
	void measure_obs_functions();

public:
	
	MonteCarloEngine(SpinModel& H_in, Wavefunction& WF_in, Lattice& lat_in, RandomEngine& rand_in, vmc_options params_in)
		: H(H_in), WF(WF_in), lat(lat_in), rand(rand_in), params(params_in) {

		for (std::string term_name : H.get_terms()) {
			observable_measures.insert(std::pair<std::string, std::vector<std::complex<double>>> (term_name, std::vector<std::complex<double>>({})));
		}
	}

	void add_observable(Observable obs_) {

	}

	void add_observable_function(std::vector<Observable> obs_, std::string name_) {
		observable_function_measures.insert(std::pair<std::string, std::vector<std::vector<std::complex<double>>>>(name_, std::vector<std::vector<std::complex<double>>>()));
		observable_functions.insert(std::pair<std::string, std::vector<Observable>>(name_, obs_));
		observable_function_names.push_back(name_);
	}

	std::vector<std::vector<double>> run_bin(bool opt) {

		std::vector<std::vector<double>> v_param_derivatives;

		for (int m = 0; m < params.num_measures; ++m) {
			timer.flag_start_time("Steps");
			for (int s = 0; s < params.steps_per_measure; ++s) {
				if (params.su3) {
					step_su3(2);
				}
				else {
					step_su2();
				}
			}
			timer.flag_end_time("Steps");
			timer.flag_start_time("Energy Calculation");
			measure_energy();
			timer.flag_end_time("Energy Calculation");
			timer.flag_start_time("Obs Calculation");
			measure_obs_functions();
			timer.flag_end_time("Obs Calculation");

			if (opt) {
				v_param_derivatives.push_back(WF.log_derivative());
			}
			else {
				if (m % 1000 == 0) {
					std::cout << "measurement " << m << " of " << params.num_measures << "\n";
				}
			}
		}

		return v_param_derivatives;
	}

	void run();

	void output_energy() {
		double e_avg;
		e_avg = std::accumulate(E.begin() + params.throwaway_measures, E.end(), 0.0) / (E.size() - params.throwaway_measures);
		std::cout << "Average Energy: \n" << "E = " << e_avg / lat.get_N() << " +- " << error(E.begin() + params.throwaway_measures, E.end(), e_avg, 10) / lat.get_N() << "\n";
	}

	void output_Nz() {
		std::cout << "Average Nz: \n" << "<Nz> = " << std::accumulate(Nz.begin() + params.throwaway_measures, Nz.end(), 0.0) / (Nz.size() - params.throwaway_measures) 
			<< "\n<Nz>/N = " << std::accumulate(Nz.begin() + params.throwaway_measures, Nz.end(), 0.0) / (Nz.size() - params.throwaway_measures) / lat.get_N() << "\n";
	}

	std::complex<double> get_energy() {
		std::complex<double> e_avg;
		if (params.optimization) {
			int num_skip = params.sr.throwaway_bins * params.num_measures;
			e_avg = { std::accumulate(E.begin() + num_skip, E.end(), 0.0) / (E.size() - num_skip),
				std::accumulate(E_imag.begin() + num_skip, E_imag.end(), 0.0) / (E.size() - num_skip) };
		}
		else {
			e_avg = { std::accumulate(E.begin() + params.throwaway_measures, E.end(), 0.0) / (E.size() - params.throwaway_measures),
				std::accumulate(E_imag.begin() + params.throwaway_measures, E_imag.end(), 0.0) / (E.size() - params.throwaway_measures) };
		}
		return e_avg;
	}



	std::complex<double> get_observable(std::string term) {
		std::vector<std::complex<double>> vals = observable_measures[term];
		std::complex<double> O_avg(0.0, 0.0);
		int num_skip = params.throwaway_measures;
		if (params.optimization) {
			num_skip = params.sr.throwaway_bins * params.num_measures;
		}
		for (auto v = vals.begin() + num_skip; v != vals.end(); ++v) {
			O_avg += *v;
		}
		return std::complex<double>(O_avg.real() / (vals.size() - num_skip), O_avg.imag() / (vals.size() - num_skip));
	}

	std::map<std::string, std::complex<double>> get_all_observables() {
		std::map<std::string, std::complex<double>> result;
		for (auto item : observable_measures) {
			result.insert({ item.first, get_observable(item.first) });
		}
		return result;
	}

	void write_observable_functions(std::ofstream* f_r, std::ofstream* f_i, std::string function_name) {
		for (int i = 0; i < observable_function_measures[function_name].size(); ++i) {
			for (int j = 0; j < observable_function_measures[function_name][i].size()-1; ++j) {
				*f_r << observable_function_measures[function_name][i][j].real() << ",";
				*f_i << observable_function_measures[function_name][i][j].imag() << ",";
			}
			*f_r << observable_function_measures[function_name][i][observable_function_measures[function_name][i].size() - 1].real() << "\n";
			*f_i << observable_function_measures[function_name][i][observable_function_measures[function_name][i].size() - 1].imag() << "\n";
		}
	}

	std::complex<double> get_energy_err() {
		double var_r, var_i;
		std::complex<double> e_avg;
		int num_skip = params.throwaway_measures;
		if (params.optimization) {
			num_skip = params.sr.throwaway_bins * params.num_measures;
		}
		e_avg = get_energy();
		var_r = error(E.begin() + num_skip, E.end(), get_energy().real(), 10);
		var_i = error(E_imag.begin() + num_skip, E_imag.end(), get_energy().imag(), 10);
		return { var_r, var_i };
	}

	std::complex<double> get_observable_err(std::string term) {
		double var_r, var_i;
		std::vector<std::complex<double>> vals = observable_measures[term];
		std::vector<double> vals_r = vec_r(vals), vals_i = vec_i(vals);
		std::complex<double> O_avg = get_observable(term);
		int num_skip = params.throwaway_measures;
		if (params.optimization) {
			num_skip = params.sr.throwaway_bins * params.num_measures;
		}
		var_r = error(vals_r.begin() + num_skip, vals_r.end(), O_avg.real(), 10);
		var_i = error(vals_i.begin() + num_skip, vals_i.end(), O_avg.imag(), 10);
		return { var_r, var_i };
	}

	std::map<std::string, std::complex<double>> get_all_observables_err() {
		std::map<std::string, std::complex<double>> result;
		for (auto item : observable_measures) {
			result.insert({ item.first, get_observable_err(item.first) });
		}
		return result;
	}

	void print_timers() {
		std::cout << "VMC Simulation timers:\n";
		timer.print_timers();
	}

	//find the error via the bin technique using a specified number of bins
	double error(std::vector<double>::iterator, std::vector<double>::iterator, double mean, int bins);
	//double error(std::vector<std::complex<double>>::iterator, std::vector<std::complex<double>>::iterator, std::complex<double> mean, int bins);

	//double bootstrap(std::vector<double>, int, std::string);
};
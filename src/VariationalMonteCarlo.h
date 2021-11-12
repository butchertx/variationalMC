#pragma once
#include <random>
#include <string>
#include <chrono>
#include <numeric>
#include <map>
#include "SpinModel.h"
#include "Wavefunction.h"
#include "RandomEngine.h"
#include <typeinfo>
#include <fstream>
#include <iostream>
#include "vmctype.h"

std::vector<double> vec_r(std::vector<std::complex<double>> cvec);

std::vector<double> vec_i(std::vector<std::complex<double>> cvec); 

class Lattice;

class MonteCarloEngine {

	MemTimeTester timer;

	SpinModel& H;

	Wavefunction& WF;

	Lattice& lat;

	std::vector<double> E;
	std::vector<double> E_imag;

	std::vector<double> Nz;

	RandomEngine& rand;

	VMCParams params;

	std::map<std::string, std::vector<std::complex<double>>> observable_measures;

	//void step_su2(std::string step_type);
	void step_su3(int num_swap);

	//void propose_ring_exchange(int, bool, std::vector<int>& sites, std::vector<int>& tris);

	//void measure_energy_su2();
	//void measure_energy_su3();
	void measure_energy();

public:
	
	MonteCarloEngine(SpinModel& H_in, Wavefunction& WF_in, Lattice& lat_in, RandomEngine& rand_in, VMCParams params_in)
		: H(H_in), WF(WF_in), lat(lat_in), rand(rand_in), params(params_in) {

		for (std::string term_name : H.get_terms()) {
			observable_measures.insert(std::pair<std::string, std::vector<std::complex<double>>> (term_name, std::vector<std::complex<double>>({})));
		}
	}

	void run() {
		for (int m = 0; m < params.total_measures; ++m) {
			timer.flag_start_time("Steps");
			for (int s = 0; s < params.steps_per_measure; ++s) {
				//if (step_outputs && stepnum < 100) {
				//	//write configuration
				//	stepnum = params.steps_per_measure * m + s;
				//	conffile = outdir + "conf" + std::to_string(stepnum) + ".csv";
				//	outfile.open(conffile);
				//	WF.write_configuration(&outfile);
				//	outfile.close();
				//	//write amplitudes on trimers
				//	ampfile = outdir + "amplitude" + std::to_string(stepnum) + ".csv";
				//	outfile.open(ampfile);
				//	WF.write_amplitudes(&outfile);
				//	outfile.close();
				//}
				step_su3(2);
			}
			timer.flag_end_time("Steps");
			timer.flag_start_time("Energy Calculation");
			measure_energy();
			timer.flag_end_time("Energy Calculation");
			if (m % 1000 == 0) {
				std::cout << "measurement " << m << " of " << params.total_measures << "\n";
			}
		}
	}

	//void run(bool step_outputs = false) {
	//	std::ofstream outfile;
	//	std::string outdir = "../../markov_chain/";
	//	if (step_outputs) {
	//		if (!isDirExist(outdir)) {
	//			makePath(outdir);
	//		}
	//		std::string latfile = outdir + "lattice.csv";
	//		outfile.open(latfile);
	//		lat.write_coordinates(&outfile);
	//		outfile.close();
	//		std::string ringfile = outdir + "rings.csv";
	//		outfile.open(ringfile);
	//		lat.write_rings(&outfile);
	//		outfile.close();
	//	}
	//	if (typeid(WF) == typeid(ProjectedState)) {
	//		for (int m = 0; m < params.total_measures; ++m) {
	//			timer.flag_start_time("Steps");
	//			for (int s = 0; s < params.steps_per_measure; ++s) {
	//				step_su2("NZ_CONST");
	//			}
	//			timer.flag_end_time("Steps");
	//			timer.flag_start_time("Energy Calculation");
	//			measure_energy_su2();
	//			timer.flag_end_time("Energy Calculation");
	//			Nz.push_back(conf.get_Nz());
	//		}
	//		//std::cout << "Run finished.  Continue?\n";
	//		//std::string trash;
	//		//std::getline(std::cin, trash);
	//	}
	//	else if (typeid(WF) == typeid(SU3ProductState)) {
	//		int stepnum = 0;
	//		std::string conffile;
	//		std::string ampfile;
	//		for (int m = 0; m < params.total_measures; ++m) {
	//			timer.flag_start_time("Steps");
	//			for (int s = 0; s < params.steps_per_measure; ++s) {
	//				if (step_outputs && stepnum < 100) {
	//					//write configuration
	//					stepnum = params.steps_per_measure * m + s;
	//					conffile = outdir + "conf" + std::to_string(stepnum) + ".csv";
	//					outfile.open(conffile);
	//					WF.write_configuration(&outfile);
	//					outfile.close();
	//					//write amplitudes on trimers
	//					ampfile = outdir + "amplitude" + std::to_string(stepnum) + ".csv";
	//					outfile.open(ampfile);
	//					WF.write_amplitudes(&outfile);
	//					outfile.close();
	//				}
	//				step_su3();
	//			}
	//			timer.flag_end_time("Steps");
	//			timer.flag_start_time("Energy Calculation");
	//			measure_energy_su3();
	//			timer.flag_end_time("Energy Calculation");
	//		}
	//	}
	//}

	void output_energy() {
		double e_avg;
		e_avg = std::accumulate(E.begin() + params.throwaway_measures, E.end(), 0.0) / (E.size() - params.throwaway_measures);
		std::cout << "Average Energy: \n" << "E = " << e_avg << " +- " << error(E.begin() + params.throwaway_measures, E.end(), e_avg, 10) << "\n";
	}

	void output_Nz() {
		std::cout << "Average Nz: \n" << "<Nz> = " << std::accumulate(Nz.begin() + params.throwaway_measures, Nz.end(), 0.0) / (Nz.size() - params.throwaway_measures) 
			<< "\n<Nz>/N = " << std::accumulate(Nz.begin() + params.throwaway_measures, Nz.end(), 0.0) / (Nz.size() - params.throwaway_measures) / lat.get_N() << "\n";
	}

	std::complex<double> get_energy() {
		std::complex<double> e_avg;
		e_avg = { std::accumulate(E.begin() + params.throwaway_measures, E.end(), 0.0) / (E.size() - params.throwaway_measures),
			std::accumulate(E_imag.begin() + params.throwaway_measures, E_imag.end(), 0.0) / (E.size() - params.throwaway_measures) };
		return e_avg;
	}

	std::complex<double> get_observable(std::string term) {
		std::vector<std::complex<double>> vals = observable_measures[term];
		std::complex<double> O_avg(0.0, 0.0);
		for (auto v = vals.begin() + params.throwaway_measures; v != vals.end(); ++v) {
			O_avg += *v;
		}
		return std::complex<double>(O_avg.real() / (vals.size() - params.throwaway_measures), O_avg.imag() / (vals.size() - params.throwaway_measures));
	}

	std::complex<double> get_energy_err() {
		double var_r, var_i;
		std::complex<double> e_avg;
		e_avg = get_energy();
		var_r = error(E.begin() + params.throwaway_measures, E.end(), get_energy().real(), 10);
		var_i = error(E_imag.begin() + params.throwaway_measures, E_imag.end(), get_energy().imag(), 10);
		return { var_r, var_i };
	}

	std::complex<double> get_observable_err(std::string term) {
		double var_r, var_i;
		std::vector<std::complex<double>> vals = observable_measures[term];
		std::vector<double> vals_r = vec_r(vals), vals_i = vec_i(vals);
		std::complex<double> O_avg = get_observable(term);
		var_r = error(vals_r.begin() + params.throwaway_measures, vals_r.end(), O_avg.real(), 10);
		var_i = error(vals_i.begin() + params.throwaway_measures, vals_i.end(), O_avg.imag(), 10);
		return { var_r, var_i };
	}

	void print_timers() {
		timer.print_timers();
	}

	//find the error via the bin technique using a specified number of bins
	double error(std::vector<double>::iterator, std::vector<double>::iterator, double mean, int bins);
	//double error(std::vector<std::complex<double>>::iterator, std::vector<std::complex<double>>::iterator, std::complex<double> mean, int bins);

	//double bootstrap(std::vector<double>, int, std::string);
};
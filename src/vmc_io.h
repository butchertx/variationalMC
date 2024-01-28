#pragma once
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include "json.hpp"
#include "vmctype.h"
#include <sys/stat.h> // stat
#include <errno.h>    // errno, ENOENT, EEXIST
#if defined(_WIN32)
    #include <direct.h> //mkdir
#endif

// for convenience
using json = nlohmann::json;
using namespace vmctype;

//------------------------------------------
//	Forward Declarations
//------------------------------------------

//------------------------------------------
//	General Helper Functions
//------------------------------------------

bool doesDirExist(const std::string& path);

bool makePath(const std::string& path);


//------------------------------------------
//	I/O Functions
//------------------------------------------

template <class T>
std::string vec2str(std::vector<T> vec) {
	std::stringstream ss;
	for (int i = 0; i < vec.size() - 1; ++i) {
		ss << vec[i] << ", ";
	}
	ss << vec.back();
	return ss.str();
}

lattice_options read_json_lattice(json j);
lattice_options read_json_lattice(std::string infile_name);
lattice_options read_json_lattice_from_dir(const std::string& dir_name);

mean_field_options read_json_wavefunction(json j);
mean_field_options read_json_wavefunction(std::string infile_name);
mean_field_options read_json_wavefunction_from_dir(const std::string& dir_name);

model_options read_json_model(json j);
model_options read_json_model(std::string infile_name);
model_options read_json_model_from_dir(const std::string& dir_name);

vmc_options read_json_vmc(json j);
vmc_options read_json_vmc(std::string infile_name);
vmc_options read_json_vmc_from_dir(const std::string& dir_name);

void read_json_full_input(lattice_options* lat, mean_field_options* wf, model_options* H, vmc_options* vmc, std::string infile_name);


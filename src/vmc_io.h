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

LatticeOptions read_json_lattice(json j);
LatticeOptions read_json_lattice(std::string infile_name);
LatticeOptions read_json_lattice_from_dir(const std::string& dir_name);

MeanFieldOptions read_json_wavefunction(json j);
MeanFieldOptions read_json_wavefunction(std::string infile_name);
MeanFieldOptions read_json_wavefunction_from_dir(const std::string& dir_name);

ModelOptions read_json_model(json j);
ModelOptions read_json_model(std::string infile_name);
ModelOptions read_json_model_from_dir(const std::string& dir_name);

VMCOptions read_json_vmc(json j);
VMCOptions read_json_vmc(std::string infile_name);
VMCOptions read_json_vmc_from_dir(const std::string& dir_name);

void read_json_full_input(LatticeOptions* lat, MeanFieldOptions* wf, ModelOptions* H, VMCOptions* vmc, std::string infile_name);


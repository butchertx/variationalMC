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
#include "mkl.h"

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
//	I/O Print Functions
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

namespace vmc_io {

	template <class T>
	void print_matrix(const char* desc, int m, int n, T* a, int lda){
		std::cout << desc << ":\n";
		// this is gross but I don't know how else to format print a complex number
		int WIDTH = 20;
		int token_width = 0;
		int num_spaces = 0;
		std::stringstream ss;
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n-1; ++j) {
				ss.str("");
				ss << a[i * lda + j];
				token_width = ss.str().length();
				num_spaces = WIDTH - token_width;
				if (num_spaces < 1) {
					num_spaces = 1;
				}
				std::cout << std::string(num_spaces, ' ') << a[i * lda + j] << ",";
			}
			std::cout << a[i * lda + n - 1] << "\n";
		}
	}
	
}

//------------------------------------------
//	I/O File Functions
//------------------------------------------

LatticeOptions read_json_lattice(json j);
LatticeOptions read_json_lattice(std::string infile_name);
LatticeOptions read_json_lattice_from_dir(const std::string& dir_name);

WavefunctionOptions read_json_wavefunction(json j);
WavefunctionOptions read_json_wavefunction(std::string infile_name);
WavefunctionOptions read_json_wavefunction_from_dir(const std::string& dir_name);

ModelOptions read_json_model(json j);
ModelOptions read_json_model(std::string infile_name);
ModelOptions read_json_model_from_dir(const std::string& dir_name);

VMCOptions read_json_vmc(json j);
VMCOptions read_json_vmc(std::string infile_name);
VMCOptions read_json_vmc_from_dir(const std::string& dir_name);

void read_json_full_input(LatticeOptions* lat, WavefunctionOptions* wf, ModelOptions* H, VMCOptions* vmc, std::string infile_name);


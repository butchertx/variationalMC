#include "vmc_io.h"

//------------------------------------------
//	Forward Declarations
//------------------------------------------


//------------------------------------------
//	General Helper Functions
//------------------------------------------



bool isDirExist(const std::string& path)
{
#if defined(_WIN32)
	struct _stat info;
	if (_stat(path.c_str(), &info) != 0)
	{
		return false;
	}
	return (info.st_mode & _S_IFDIR) != 0;
#else 
	struct stat info;
	if (stat(path.c_str(), &info) != 0)
	{
		return false;
	}
	return (info.st_mode & S_IFDIR) != 0;
#endif
}

bool makePath(const std::string& path)
{
#if defined(_WIN32)
	int ret = _mkdir(path.c_str());
#else
	mode_t mode = 0755;
	int ret = mkdir(path.c_str(), mode);
#endif
	if (ret == 0)
		return true;

	switch (errno)
	{
	case ENOENT:
		// parent didn't exist, try to create it
	{
		int pos = path.find_last_of('/');
		if (pos == std::string::npos)
#if defined(_WIN32)
			pos = path.find_last_of('\\');
		if (pos == std::string::npos)
#endif
			return false;
		if (!makePath(path.substr(0, pos)))
			return false;
	}
	// now, try to create again
#if defined(_WIN32)
	return 0 == _mkdir(path.c_str());
#else 
	return 0 == mkdir(path.c_str(), mode);
#endif

	case EEXIST:
		// done!
		return isDirExist(path);

	default:
		return false;
	}
}



lattice_options read_json_lattice(json j) {
	lattice_options lat_opt;
	lat_opt.type = j["lattice"]["type"].get<std::string>();
	lat_opt.dimension = j["lattice"]["dimension"];
	lat_opt.L = j["lattice"]["L"].get<std::vector<int>>();
	lat_opt.pbc = j["lattice"]["pbc"].get<std::vector<int>>();

	lat_opt.Lx = lat_opt.L[0];
	lat_opt.Ly = lat_opt.L[1];
	lat_opt.Lz = lat_opt.L[2];
	return lat_opt;
}

lattice_options read_json_lattice(std::string infile_name) {
	std::ifstream i(infile_name);
	json j;
	i >> j;
	return read_json_lattice(j);
}

mean_field_options read_json_wavefunction(json j) {
	mean_field_options wf_opt;
	wf_opt.lattice_type = j["wavefunction"]["lattice type"].get<std::string>();
	wf_opt.wf_type = j["wavefunction"]["wavefunction type"].get<std::string>();

	if (wf_opt.wf_type == "meanfield") {
		if (j["wavefunction"].contains("basis")) {
			wf_opt.basis = j["wavefunction"]["basis"].get<std::vector<vec3<int>>>();
			wf_opt.inequivalent_sites = j["wavefunction"]["inequivalent sites"];
		}
		wf_opt.spin = j["wavefunction"]["spin"];
		wf_opt.field = j["wavefunction"]["field"];
		wf_opt.num_spin_orbit = j["wavefunction"]["num spin-orbit"];
		wf_opt.match_lattice_pbc = j["wavefunction"]["match lattice pbc"];
		wf_opt.su3_symmetry = j["wavefunction"]["su3_symmetry"];
		wf_opt.hopping_list = j["wavefunction"]["hopping terms"].get<std::vector<HoppingTerm>>();
		if (j["wavefunction"].contains("directors")) {
			wf_opt.directors = j["wavefunction"]["directors"].get<QuadrupoleOrder>();
		}
		if (j["wavefunction"].contains("jastrow")) {
			wf_opt.jastrow_flag = true;
			if (j["wavefunction"]["jastrow"].contains("sz")) {
				wf_opt.jastrow.sz = j["wavefunction"]["jastrow"]["sz"].get<JastrowFactorOptions>();
			}
		}
	}
	else {
		throw;
	}
	return wf_opt;
}

mean_field_options read_json_wavefunction(std::string infile_name) {
	std::ifstream i(infile_name);
	json j;
	i >> j;
	return read_json_wavefunction(j);
}

model_options read_json_model(json j) {
	model_options model_opt;
	if (j["model"].contains("bilinear")) {
		model_opt.bilinear_terms = j["model"]["bilinear"].get<std::vector<BilinearOptions>>();
	}
	if (j["model"].contains("trilinear")) {
		model_opt.ring3_terms = j["model"]["trilinear"].get<std::vector<TrilinearOptions>>();
	}
	return model_opt;
}

model_options read_json_model(std::string infile_name) {
	std::ifstream i(infile_name);
	json j;
	i >> j;
	return read_json_model(j);
}

vmc_options read_json_vmc(json j) {
	vmc_options vmc_opt;
	if (j["vmc"].contains("optimization")) {
		vmc_opt.optimization = j["vmc"]["optimization"].get<bool>();
	}
	else {
		vmc_opt.optimization = false;
	}
	if (j["vmc"].contains("su3")) {
		vmc_opt.su3 = j["vmc"]["su3"].get<bool>();
	}
	else {
		vmc_opt.su3 = true;
	}
	
	if (!vmc_opt.optimization) {
		vmc_opt.num_measures = j["vmc"]["measures"].get<int>();
		vmc_opt.steps_per_measure = j["vmc"]["steps"].get<int>();
		vmc_opt.throwaway_measures = j["vmc"]["throwaway"].get<int>();
	}
	return vmc_opt;
}

vmc_options read_json_vmc(std::string infile_name) {
	std::ifstream i(infile_name);
	json j;
	i >> j;
	return read_json_vmc(j);
}

void read_json_full_input(lattice_options* lat, mean_field_options* wf, model_options* H, vmc_options* vmc, std::string infile_name) {
	std::ifstream i(infile_name);
	json j;
	i >> j;
	*lat = read_json_lattice(j);
	*wf = read_json_wavefunction(j);
	*H = read_json_model(j);
	*vmc = read_json_vmc(j);
}


std::string lattice_options::to_string() {
		std::stringstream ss;
		ss << "Lattice=" << type << "\n"
			<< "dimension=" << dimension << "\n"
			<< "Lx=" << Lx << "\n" << "Ly=" << Ly << "\n" << "Lz=" << Lz << "\n";
		return ss.str();
}


std::string mean_field_options::to_string() {
		std::stringstream ss;
		ss << "Lattice=" << lattice_type << "\n"
			<< "Inequivalent sites=" << inequivalent_sites << "\n";
			//<< "tz magnitudes: " << vec2str(tz_list) << "\n"
			//<< "tz phases (degrees): " << vec2str(tz_phase_list) << "\n"
			//<< "txy magnitudes: " << vec2str(txy_list) << "\n"
			//<< "tz phases (degrees): " << vec2str(txy_phase_list) << "\n";
		return ss.str();
}

//void read_lattice_input(std::ifstream* file_p, lattice_options* params) {
//	std::string line;
//	while (std::getline(*file_p, line)) {
//		std::stringstream str(line);
//		std::string variable_name;
//		std::string string_in;
//		//double double_in;
//		int int_in;
//		if ((std::getline(str, variable_name, '='))) {
//
//			if (variable_name.compare("dimension") == 0 || variable_name.compare("Dimension") == 0 || variable_name.compare("dim") == 0) {
//				if (str >> int_in) {
//					params->dimension = int_in;
//				}
//				else {
//					std::cout << "Error: Invalid value for parameter " << variable_name << "\n";
//				}
//			}
//			else if (variable_name.compare("lattice") == 0 || variable_name.compare("Lattice") == 0 || variable_name.compare("lat") == 0) {
//				if (str >> string_in) {
//					if (string_in.compare("chain") == 0) {
//						params->type = CHAIN;
//					}
//					else if (string_in.compare("kagome") == 0) {
//						params->type = KAGOME;
//					}
//					else if (string_in.compare("square") == 0) {
//						params->type = SQUARE;
//					}
//					else if (string_in.compare("diamond") == 0) {
//						params->type = DIAMOND;
//					}
//					else if (string_in.compare("diamond111") == 0) {
//						params->type = DIAMOND111;
//					}
//					else if (string_in.compare("diamondcc") == 0) {
//						params->type = DIAMONDCC;
//					}
//					else {
//						std::cout << "Error: Invalid lattice type entered, " << string_in << " is not a valid lattice type\n";
//					}
//				}
//				else {
//					std::cout << "Error: Invalid value for parameter " << variable_name << "\n";
//				}
//			}
//			else if (variable_name.compare("Lx") == 0 || variable_name.compare("lx") == 0) {
//				if (str >> int_in) {
//					params->Lx = int_in;
//				}
//				else {
//					std::cout << "Error: Invalid value for parameter " << variable_name << "\n";
//				}
//			}
//			else if (variable_name.compare("Ly") == 0 || variable_name.compare("ly") == 0) {
//				if (str >> int_in) {
//					params->Ly = int_in;
//				}
//				else {
//					std::cout << "Error: Invalid value for parameter " << variable_name << "\n";
//				}
//			}
//			else if (variable_name.compare("Lz") == 0 || variable_name.compare("lz") == 0) {
//				if (str >> int_in) {
//					params->Lz = int_in;
//				}
//				else {
//					std::cout << "Error: Invalid value for parameter " << variable_name << "\n";
//				}
//			}
//			else {
//				std::cout << "Error: " << variable_name << " is not a valid input variable.  Make sure there are no spaces.\n";
//			}
//		}
//		else {
//			std::cout << "Error: must have input data in the form 'variable name = value'\n";
//		}
//	}
//
//	assert(params->is_valid());
//}
//
//void read_mean_field_input(std::ifstream* file_p, mean_field_options* params) {
//
//	std::vector<std::string> line_list;
//
//	//lattice type
//	line_list = getNextLineAndSplitIntoTokens(*file_p);
//	if (line_list[1].compare("diamond") == 0) {
//		params->lattice_type = DIAMOND;
//	}
//	else {
//		std::cout << "Error: Invalid value " << line_list[1] << " for parameter " << line_list[0] << "\n";
//	}
//	
//	//Inequivalent sites
//	line_list = getNextLineAndSplitIntoTokens(*file_p);
//	if (line_list[0].compare("Inequivalent sites") == 0) {
//		params->inequivalent_sites = std::stod(line_list[1]);
//	}
//	else {
//		std::cout << "Error: Invalid value " << line_list[1] << " for parameter " << line_list[0] << "\n";
//	}
//	
//	while (getNextLineAndSplitIntoTokens(*file_p)[0].compare("Hopping") == 0) {
//		std::cout << "Made it past \"Hopping\" line\n";
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		if (line_list.size() == 4) {
//			params->tz_list.push_back(std::stod(line_list[0]));
//			params->tz_phase_list.push_back(std::stod(line_list[1]));
//			params->txy_list.push_back(std::stod(line_list[2]));
//			params->txy_phase_list.push_back(std::stod(line_list[3]));
//			params->hopping_list.push_back({});
//		}
//		else {
//			std::cout << "Need 4 values for hopping params\n";
//		}
//
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		while (line_list.size() == 4) {
//			params->hopping_list[params->hopping_list.size() - 1]
//				.push_back(origin_direction_pair(std::stoi(line_list[0]), std::stod(line_list[1]), std::stod(line_list[2]), std::stod(line_list[3]), std::stod(line_list[4])));
//			line_list = getNextLineAndSplitIntoTokens(*file_p);
//		}
//	}
//}
//
//void read_mean_field_input(std::ifstream* file_p, mean_field_options* params, Lattice* lat) {
//
//	std::vector<std::string> line_list;
//
//	//lattice type
//	line_list = getNextLineAndSplitIntoTokens(*file_p);
//	if (line_list[1].compare(Lattice_type_to_string(lat->get_lattice_type())) == 0) {
//		params->lattice_type = lat->get_lattice_type();
//	}
//	else {
//		std::cout << "Error: Invalid value " << line_list[1] << " for parameter " << line_list[0] << "\n";
//	}
//
//	//Inequivalent sites
//	line_list = getNextLineAndSplitIntoTokens(*file_p);
//	if (line_list[0].compare("Inequivalent sites") == 0) {
//		params->inequivalent_sites = std::stod(line_list[1]);
//	}
//	else {
//		std::cout << "Error: Invalid value " << line_list[1] << " for parameter " << line_list[0] << "\n";
//	}
//
//	line_list = getNextLineAndSplitIntoTokens(*file_p);
//	while (line_list[0].compare("Hopping") == 0) {
//		std::cout << "Made it past \"Hopping\" line\n";
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		if (line_list.size() == 4) {
//			params->tz_list.push_back(std::stod(line_list[0]));
//			params->tz_phase_list.push_back(std::stod(line_list[1]));
//			params->txy_list.push_back(std::stod(line_list[2]));
//			params->txy_phase_list.push_back(std::stod(line_list[3]));
//			params->hopping_list.push_back({});
//		}
//		else {
//			std::cout << "Need 4 values for hopping params\n";
//		}
//
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		while (line_list.size() == 3) {
//			//std::cout << std::stoi(line_list[0]) << "," << std::stoi(line_list[1]) << vec3str(lat->nearest_periodic_translation(std::stoi(line_list[0]), std::stoi(line_list[1]))) << "\n";
//			params->hopping_list[params->hopping_list.size() - 1]
//				.push_back(origin_direction_pair(std::stoi(line_list[0]), lat->nearest_periodic_translation(std::stoi(line_list[0]), std::stoi(line_list[1])),std::stod(line_list[2])));
//			line_list = getNextLineAndSplitIntoTokens(*file_p);
//			//std::cout << line_list[0] << "\n";
//		}
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//	}
//
//	
//	while (line_list[0].compare("Singlet Pairing") == 0) {
//		std::cout << "Made it past \"Singlet Pairing\" line\n";
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		if (line_list.size() == 2) {
//			params->delta_s_list.push_back(std::stod(line_list[0]));
//			params->delta_s_phase_list.push_back(std::stod(line_list[1]));
//			params->singlet_pair_list.push_back({});
//		}
//		else {
//			std::cout << "Need 2 values for singlet pairing params\n";
//		}
//
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		line_list = getNextLineAndSplitIntoTokens(*file_p);
//		while (line_list.size() == 3) {
//			//std::cout << std::stoi(line_list[0]) << "," << std::stoi(line_list[1]) << vec3str(lat->nearest_periodic_translation(std::stoi(line_list[0]), std::stoi(line_list[1]))) << "\n";
//			params->singlet_pair_list[params->singlet_pair_list.size() - 1]
//				.push_back(origin_direction_pair(std::stoi(line_list[0]), lat->nearest_periodic_translation(std::stoi(line_list[0]), std::stoi(line_list[1])), std::stod(line_list[2])));
//			line_list = getNextLineAndSplitIntoTokens(*file_p);
//		}
//	}
//}
//
//void read_vmc_input(std::ifstream* file_p, vmc_options* params) {
//
//
//}
//
//std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str) {
//	std::vector<std::string>   result;
//	std::string                line;
//	std::getline(str, line);
//
//	std::stringstream          lineStream(line);
//	std::string                cell;
//
//	while (std::getline(lineStream, cell, ','))
//	{
//		result.push_back(cell);
//	}
//	// This checks for a trailing comma with no data after it.
//	if (!lineStream && cell.empty())
//	{
//		// If there was a trailing comma then add an empty element.
//		result.push_back("");
//	}
//	return result;
//}

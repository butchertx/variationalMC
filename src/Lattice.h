#pragma once
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <algorithm>
#include <functional>
#include "MemTimeTester.h"
#include "vmctype.h"

using namespace vmctype;

/*! \brief Enum for specifying a lattice type
*/
enum Lattice_type_t { CHAIN, TRIANGLE, KAGOME, DICE, HONEYCOMB, SQUARE, CUBIC, DIAMOND, DIAMOND111, DIAMONDCC };//require SQUARE is the last 2-D lattice in list

/*! \brief Get string for Lattice_type_t
 */
std::string Lattice_type_to_string(Lattice_type_t lat);

/*! \brief Get Lattice_type_t for string
 */
Lattice_type_t Lattice_type_from_string(std::string str);

/*! \brief List of sites for Ring Exchanges
*/
class RingList {
	//list groups of indices for ring-exchange type terms
	//always listed in counter-clockwise direction
	int num_sites;//number of sites in each ring
	int total_sites;
	int num_rings;//number of rings
	std::vector<std::vector<int>> rings;
	std::vector<std::vector<int>> ring_labels;//ring labels for each site
	std::vector<std::string> ring_names;//label the rings up or down

public:
	RingList() {}

	RingList(int num_sites_in, int num_rings_in, int total_sites_in) 
		: num_sites(num_sites_in), num_rings(num_rings_in), total_sites(total_sites_in){
		for (int r = 0; r < num_rings; ++r) {
			if (r % 2 == 0) { ring_names.push_back("up"); }
			else { ring_names.push_back("down"); }//TODO: make this more robust (currently is fine for kagome)
			rings.push_back({});
			for (int s = 0; s < num_sites; ++s) {
				rings[r].push_back(0);
			}
			rings[r].shrink_to_fit();
		}
		rings.shrink_to_fit();
		
		for (int i = 0; i < total_sites; ++i) {
			ring_labels.push_back({});
		}
		ring_labels.shrink_to_fit();
	}

	void set_ring(int ring_label, int site1, int site2, int site3) {
		rings[ring_label] = { site1, site2, site3 };
		ring_labels[site1].push_back(ring_label);
		ring_labels[site2].push_back(ring_label);
		ring_labels[site3].push_back(ring_label);
	}

	std::vector<int> get_ring_labels(int site) {
		return ring_labels[site];
	}

	std::vector<int> get_ring(int index) {
		return rings[index];
	}

	std::string get_ring_name(int index) {
		return ring_names[index];
	}

	int get_size() {
		return rings.size();
	}

	int get_total_sites() {
		return total_sites;
	}

};

/*! \brief Lattice Class
*
*  Contains sites, neighbors, and various utilities for defining and using a Lattice
*/
class Lattice {
private:
	MemTimeTester timer;

	Lattice_type_t lat_type;
	vec3<vec3<double>> a;//lattice vectors
	double a_const;//lattice constant
	std::vector<vec3<double>> basis;
	vec3<vec3<double>> b;//reciprocal lattice vectors
	vec3<int> L;//dimensions
	int N;//total number of sites
    vec3<int> pbc_signs;//periodic/antiperiodic boundary conditions in each direction
	std::vector<vec3<double>> coordinates;
	std::vector<double> distances;
	std::vector<std::vector<std::vector<int>>> neighbors;
    std::vector<std::vector<std::vector<int>>> pbc;
	RingList rings1;

	double nearest_periodic_distance(int, int, vec3<int>&);
	double nearest_periodic_distance(int, int);
	void set_neighbors(int);
	void set_neighbors_faster();
	void set_coordinates();
	void set_reciprocal_vectors();
	void set_distances();
	void set_rings();
	std::vector<int> get_ring_sites(int base_vertex, std::vector<vec3<double>> vertices, int neighbor_distance);

public:

	vec3<double> nearest_periodic_translation(int, int);

	Lattice();
	Lattice(Lattice_type_t, vec3<int>, vec3<int>);
	~Lattice();

	RingList& ring_ref() {
		return rings1;
	}



	Lattice_type_t get_lattice_type() {
		return lat_type;
	}

    int get_N(){
        return N;
    }

    int get_Lx(){
        return L.x;
    }

    int get_Ly(){
        return L.y;
    }

    int get_Lz(){
        return L.z;
    }

    int get_num_base(){
        return basis.size();
    }

    int get_num_distance(){
        return neighbors[0].size();
    }

    int get_max_neighbor_count(){
        int max = 1;
        for(int i = 0; i < neighbors.size(); ++i){
            for(int j = 0; j < neighbors[i].size(); ++j){
                max = neighbors[i][j].size() > max ? neighbors[i][j].size() : max;
            }
        }
        return max;
    }

    std::vector<int> get_neighbor_counts(){
        std::vector<int> imulti(neighbors[0].size());
        for(int i = 0; i < imulti.size(); ++i){
            imulti[i] = neighbors[0][i].size();
        }
        return imulti;
    }

	int get_neighbor(int base_site, int index) {
		if (index < neighbors[base_site][0].size()) {
			return neighbors[base_site][0][index];
		}
		else if (index < neighbors[base_site][1].size()) {
			return neighbors[base_site][1][index];
		}
		else {
			std::cout << "Invalid neighbor request; Exiting Code.\n";
			exit(0);
		}
	}

    std::vector<int> get_neighbors(int base_site, int distance){
		assert(base_site < N && base_site >= 0);
        assert(distance >= 0 && distance < neighbors[0].size());
        std::vector<int> neighbor_list(neighbors[base_site][distance].size());
        for(int j = 0; j < neighbor_list.size(); ++j){
            neighbor_list[j] = neighbors[base_site][distance][j];
        }
        return neighbor_list;
    }

	std::vector<std::vector<int>> get_neighbors(int distance) {
		assert(distance >= 0 && distance < neighbors[0].size());
		std::vector<std::vector<int>> neighbor_list;
		for (int j = 0; j < neighbors.size(); ++j) {
			neighbor_list.push_back(neighbors[j][distance]);
		}
		return neighbor_list;
	}

	std::vector<int> get_neighbor_pbcs(int base_site, int distance) {
        assert(base_site < N && base_site >= 0);
        assert(distance >= 0 && distance < neighbors[0].size());
        std::vector<int> neighbor_list(neighbors[base_site][distance].size());
        for(int j = 0; j < neighbor_list.size(); ++j){
            neighbor_list[j] = pbc[base_site][distance][j];
        }
        return neighbor_list;
    }

	int get_neighbor_with_pbc(int base_site, int distance, int neighbor_index) {
		int neighbor, pbc;
		neighbor = get_neighbors(base_site, distance)[neighbor_index];
		pbc = get_neighbor_pbcs(base_site, distance)[neighbor_index];
		return neighbor * pbc;
	}

    vec3<double> get_coordinate(int site){
        assert(site >=0 && site < N);
        return coordinates[site];
    }

    vec3<double> get_momentum(int kx, int ky, int kz){
        assert(kx >= 0 && kx <= L.x);
        assert(ky >= 0 && ky <= L.y);
        assert(kz >= 0 && kz <= L.z);
        return b.x*(((double)kx)/L.x) + b.y*(((double)ky)/L.y) + b.z*(((double)kz)/L.z);
    }

	void print_neighbors() {
		std::cout << "Table of Neighbors:\n";
		std::stringstream ss;
		for (int neighbor = 0; neighbor < neighbors[0].size(); ++neighbor) {
			std::cout << "At distance k = " << neighbor+1 << "\n\n";
			for (int i = 0; i < N; ++i) {
				ss.str("");
				ss << "i = " << i << ";  ";
				for (int j = 0; j < neighbors[i][neighbor].size(); ++j) {
					ss << neighbors[i][neighbor][j]*pbc[i][neighbor][j] << "  ";
				}
				ss << "\n";
				std::cout << ss.str();
			}
		}
	}

	void print_neighbors(std::ofstream* f) {
		*f << "Table of Neighbors:\n";
		std::stringstream ss;
		for (int neighbor = 0; neighbor < neighbors[0].size(); ++neighbor) {
			*f << "At distance k = " << neighbor + 1 << "\n\n";
			for (int i = 0; i < N; ++i) {
				ss.str("");
				ss << "i = " << i << ";  ";
				for (int j = 0; j < neighbors[i][neighbor].size(); ++j) {
					ss << neighbors[i][neighbor][j] * pbc[i][neighbor][j] << "  ";
				}
				ss << "\n";
				*f << ss.str();
			}
		}
	}

	void print_neighbors(int max_distance) {
		std::cout << "Table of Neighbors:\n";
		std::stringstream ss;
		for (int neighbor = 0; neighbor < max_distance; ++neighbor) {
			std::cout << "At distance k = " << neighbor + 1 << "\n\n";
			for (int i = 0; i < N; ++i) {
				ss.str("");
				ss << "i = " << i << ";  ";
				for (int j = 0; j < neighbors[i][neighbor].size(); ++j) {
					ss << neighbors[i][neighbor][j] * pbc[i][neighbor][j] << "  ";
				}
				ss << "\n";
				std::cout << ss.str();
			}
		}
	}

	void print_neighbors(int origin, int distance) {
		std::cout << "Table of Neighbors for site " << origin << " to distance " << distance << ":\n";
		std::stringstream ss;
		int neighbor = distance - 1;
		int i = origin;
		ss.str("");
		ss << "i = " << i << ":\n";
		for (int j = 0; j < neighbors[i][neighbor].size(); ++j) {
			ss << neighbors[i][neighbor][j] * pbc[i][neighbor][j] << "    " << vec3str(nearest_periodic_translation(i, neighbors[i][neighbor][j])) << "\n";
		}
		ss << "\n";
		std::cout << ss.str();
	}

	void print_rings(std::ofstream* f) {
		*f << "Table of Ring Exchanges:\n";
		std::stringstream ss;
		std::vector<int> ring_vals;
		for (int ring = 0; ring < rings1.get_size(); ++ring) {
			ring_vals = rings1.get_ring(ring);
			ss.str("");
			ss << "ring index = " << ring << ";  ";
			for (int i = 0; i < ring_vals.size(); ++i) {
				ss << ring_vals[i] << ",  ";				
			}
			ss << "\n";
			*f << ss.str();
		}
	}

	void print_coordinates() {
		std::cout << "Lattice Coordinates:\n";
		for (int i = 0; i < N; ++i) {
			std::cout << "i=" << i << "    <" << coordinates[i].x << "," << coordinates[i].y << "," << coordinates[i].z << ">\n";
		}
		std::cout << "\n";
	}

	void write_coordinates(std::ofstream* file_out) {
		*file_out << "site, x, y, z\n";
		for (int i = 0; i < N; ++i) {
			*file_out << i << "," << coordinates[i].x << "," << coordinates[i].y << "," << coordinates[i].z << "\n";
		}
	}

	void write_rings(std::ofstream* file_out) {
		*file_out << "ring_idx, site1, site2, site3\n";
		for (int i = 0; i < rings1.get_size(); ++i) {
			std::vector<int> ring = rings1.get_ring(i);
			*file_out << i << "," << ring[0] << "," << ring[1] << "," << ring[2] << "\n";
		}
	}

	void print_lattice_vectors(bool real_space) {
		vec3<vec3<double>> out;
		if (real_space) {
			out = a;
		}
		else {
			out = b;
		}
		std::cout << (real_space ? "Lattice Vectors:\n" : "Reciprocal Lattice Vectors:\n")
			<< "<" << out.x.x << "," << out.x.y << "," << out.x.z << ">\n"
			<< "<" << out.y.x << "," << out.y.y << "," << out.y.z << ">\n"
			<< "<" << out.z.x << "," << out.z.y << "," << out.z.z << ">\n\n";
	}

	void print_distances() {
		std::cout << "Distances:\n";
		for (int i = 0; i < distances.size(); ++i) {
			std::cout << distances[i] << "\n";
		}
		std::cout << "\n";
	}

	int get_neighbor_label_with_pbc(int origin, vec3<double> shift) {
		int result = origin;
		//shift = a.x * shift.x + a.y * shift.y + a.z * shift.z; USE THIS FOR OLD VERSION OF MFUC
		double dist = sqrt(shift * shift) - EPSILON;
		int dist_ind = std::lower_bound(distances.begin(), distances.end(), dist) - distances.begin();
		//std::cout << "(dist_ind, Dist) = " << dist_ind << "," << dist << "\n";
		std::vector<int> neighbor_list = get_neighbors(origin, dist_ind);
		std::vector<int> pbc_list = get_neighbor_pbcs(origin, dist_ind);
		//std::cout << "origin,shift: " << origin << "," << vec3str(shift) << "\n";
		for (auto label = neighbor_list.begin(); label != neighbor_list.end(); ++label) {
			//std::cout << "neighbor: " << *label << "\n";
			if (check_shift(origin, *label, shift)) {
				result = *label * pbc_list[label - neighbor_list.begin()];
			}
		}
		//std::cout << "result, origin " << result << "," << origin << "\n";
		assert(std::abs(result) != origin);
		return result;
	}

	bool check_shift(int origin, int neighbor, vec3<double> shift) {
		//check periodic images to see if origin and neighbor are separated by vector shift
		bool result = false;
		vec3<double> lattice_shift;
		for (int i = -1; i <= 1; ++i) {
			for (int j = -1; j <= 1; ++j) {
				for (int k = -1; k <= 1; ++k) {
					lattice_shift = a.x * (i*L.x) + a.y * (j*L.y) + a.z * (k*L.z);
					result = result || (get_coordinate(origin) + lattice_shift + shift == get_coordinate(neighbor));
				}
			}
		}
		return result;
	}

	int get_label(int nx, int ny, int nz, int base) {
		assert(nx >= 0 && nx < L.x);
		assert(ny >= 0 && ny < L.y);
		assert(nz >= 0 && nz < L.z);
		assert(base >= 0 && base < basis.size());
		return basis.size() * (nx*L.y*L.z + ny * L.z + nz) + base;
	}

	int get_nx(int label) {
		assert(label >= 0 && label < N);
		return (label / basis.size()) / (L.y*L.x);
	}

	int get_ny(int label) {
		assert(label >= 0 && label < N);
		return ((label / basis.size()) % (L.y*L.z)) / L.z;
	}

	int get_nz(int label) {
		assert(label >= 0 && label < N);
		return (label / basis.size()) % L.z;
	}

	int get_base(int label) {
		assert(label >= 0 && label < N);
		return label % basis.size();
	}

	void print_timers() {
		timer.print_timers();
	}

    bool check_in_plane(int,int);

	std::vector<std::vector<int>> basis_partition(std::vector<vec3<int>> basis);

};
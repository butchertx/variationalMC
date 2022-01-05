#include "Lattice.h"
#include <iostream>

/*! \brief Get string for Lattice_type_t
 */
std::string Lattice_type_to_string(Lattice_type_t lat) {
	switch (lat) {
	case CHAIN:
		return "chain";
	case TRIANGLE:
		return "triangle";
	case KAGOME:
		return "kagome";
	case DICE:
		return "dice";
	case HONEYCOMB:
		return "honeycomb";
	case SQUARE:
		return "square";
	case CUBIC:
		return "cubic";
	case DIAMOND:
		return "diamond";
	case DIAMOND111:
		return "diamond 111";
	case DIAMONDCC:
		return "diamond conventional cell";
	default:
		return "INVALID";
	}

}

/*! \brief Get Lattice_type_t for string
 */
Lattice_type_t Lattice_type_from_string(std::string str) {
	if (str.compare("chain") == 0) {
		return CHAIN;
	}
	else if (str.compare("triangle") == 0) {
		return TRIANGLE;
	}
	else if (str.compare("kagome") == 0) {
		return KAGOME;
	}
	else if (str.compare("square") == 0) {
		return SQUARE;
	}
	else if (str.compare("diamond conventional cell") == 0) {
		return DIAMONDCC;
	}
	else {
		std::cout << "Invalid lattice type from string conversion:\n" << "string: " << str << "\n";
		exit(0);
	}
}

bool vmctype::lattice_options::is_valid() {
	bool allowed = true;
	std::cout << dimension << " " << type << " " << Lx << " " << Ly << " " << Lz << "\n";
	allowed = (dimension == 1 && Lattice_type_from_string(type) == CHAIN) || (dimension == 2 && Lattice_type_from_string(type) <= SQUARE) || (dimension == 3 && Lattice_type_from_string(type) > SQUARE);
	allowed = allowed && Lattice_type_from_string(type) >= 0 && Lattice_type_from_string(type) <= DIAMONDCC;
	allowed = allowed && Lx > 0 && Ly > 0 && Lz > 0;
	return allowed;
}

///////////////////////////////////
///////////////////////////////////
// LATTICE 3D
///////////////////////////////////
///////////////////////////////////

Lattice::Lattice() {
}

Lattice::~Lattice() {}

Lattice::Lattice(Lattice_type_t lat_type_in, vec3<int> L_in, vec3<int> pbc_in)
	: lat_type(lat_type_in), L(L_in), pbc_signs(pbc_in) {
	assert(L.x > 0);
	assert(pbc_signs.x == 1 || pbc_signs.x == -1);
	assert(L.y > 0);
	assert(pbc_signs.y == 1 || pbc_signs.y == -1);
	assert(L.z > 0);
	assert(pbc_signs.z == 1 || pbc_signs.z == -1);

	vec3<double> delta_1;
	vec3<double> primitive_axy;
	vec3<double> primitive_ayz;
	vec3<double> primitive_axz;

	

	switch (lat_type) {
	case CHAIN:
		basis.push_back(vec3<double>(0.0, 0.0, 0.0));
		a_const = 1.0;
		a = vec3<vec3<double>>(vec3<double>(1.0, 0.0, 0.0), vec3<double>(0.0, 1.0, 0.0), vec3<double>(0.0, 0.0, 1.0));
		break;
	case TRIANGLE:
		a_const = 1.0;
		basis.push_back(vec3<double>(0.0, 0.0, 0.0));
		a = vec3<vec3<double>>(vec3<double>(a_const, 0.0, 0.0), vec3<double>(a_const / 2, a_const * sqrt(3) / 2, 0.0), vec3<double>(0.0, 0.0, a_const));
		break;
	case KAGOME:
		a_const = 1.0;
		basis.push_back(vec3<double>(0.0, 0.0, 0.0));
		basis.push_back(vec3<double>(a_const/2, 0.0, 0.0));
		basis.push_back(vec3<double>(a_const/4, a_const*sqrt(3)/4, 0.0));
		a = vec3<vec3<double>>(vec3<double>(a_const, 0.0, 0.0), vec3<double>(a_const/2, a_const*sqrt(3)/2, 0.0), vec3<double>(0.0, 0.0, a_const));
		break;
	case SQUARE:
		basis.push_back(vec3<double>(0.0, 0.0, 0.0));
		a_const = 1.0;
		a = vec3<vec3<double>>(vec3<double>(a_const, 0.0, 0.0), vec3<double>(0.0, a_const, 0.0), vec3<double>(0.0, 0.0, a_const));
		break;
	case CUBIC:
		basis.push_back(vec3<double>(0.0, 0.0, 0.0));
		a_const = 1.0;
		a = vec3<vec3<double>>(vec3<double>(1.0, 0.0, 0.0), vec3<double>(0.0,1.0,0.0), vec3<double>(0.0,0.0,1.0));
		break;
	case DIAMOND:
		basis.push_back(vec3<double>(0.0, 0.0, 0.0));
		basis.push_back(vec3<double>(1 / sqrt(3.0), 1 / sqrt(3.0), 1 / sqrt(3.0)));
		a_const = sqrt(16.0 / 3.0);//lattice constant: a=sqrt(16/3) maintains first neighbor distance = 1
		a = vec3<vec3<double>>(vec3<double>(0.5, 0.5, 0.0)*a_const, vec3<double>(0.0, 0.5, 0.5)*a_const, vec3<double>(0.5, 0.0, 0.5)*a_const);
		break;
    case DIAMOND111:
        a_const = sqrt(16.0/3.0);
        delta_1 = vec3<double>(1.0, 1.0, 1.0)*0.25*a_const;
        basis.push_back(vec3<double>(0.0, 0.0, 0.0));
        basis.push_back(delta_1);
        basis.push_back(vec3<double>(0.0, 1.0, 1.0)*0.5*a_const);
        basis.push_back(vec3<double>(0.0, 1.0, 1.0)*0.5*a_const + delta_1);
        basis.push_back(vec3<double>(0.0, 1.0, 1.0)*a_const);
        basis.push_back(vec3<double>(0.0, 1.0, 1.0)*a_const + delta_1);
        a = vec3<vec3<double>>(vec3<double>(-0.5, 0.5, 0.0)*a_const, vec3<double>(-0.5, 0.0, 0.5)*a_const, vec3<double>(1.0, 1.0, 1.0)*a_const);
        break;
	case DIAMONDCC:
		//Conventional unit cell.  A and B sublattices are even and odd labels, respectively
		a_const = sqrt(16.0 / 3.0);
		primitive_axy = vec3<double>(0.5, 0.5, 0.0)*a_const;
		primitive_ayz = vec3<double>(0.0, 0.5, 0.5)*a_const;
		primitive_axz = vec3<double>(0.5, 0.0, 0.5)*a_const;
		delta_1 = vec3<double>(1.0, 1.0, 1.0)*0.25*a_const;
		basis.push_back(vec3<double>(0.0, 0.0, 0.0));
		basis.push_back(delta_1);
		basis.push_back(primitive_axy);
		basis.push_back(primitive_axy + delta_1);
		basis.push_back(primitive_ayz);
		basis.push_back(primitive_ayz + delta_1);
		basis.push_back(primitive_axz);
		basis.push_back(primitive_axz + delta_1);
		a = vec3<vec3<double>>(vec3<double>(1,0,0)*a_const, vec3<double>(0,1,0)*a_const, vec3<double>(0,0,1)*a_const);
		break;
	default:
		std::cout << "Error, no correct lattice input\n";
		break;
	}

	
	N = basis.size() * L.x * L.y * L.z;

	timer.flag_start_time("Reciprocal Vectors");
	set_reciprocal_vectors();
	timer.flag_end_time("Reciprocal Vectors");

	timer.flag_start_time("Set Coordinates");
	set_coordinates();
	timer.flag_end_time("Set Coordinates");

	timer.flag_start_time("Set Distances");
	set_distances();
	timer.flag_end_time("Set Distances");

	timer.flag_start_time("Set Neighbors");
	set_neighbors_faster();
	timer.flag_end_time("Set Neighbors");

	set_rings();

}

std::vector<int> Lattice::get_ring_sites(int base_vertex, std::vector<vec3<double>> vertices, int neighbor_distance = 0) {
	std::vector<int> nei1 = get_neighbors(base_vertex, neighbor_distance);
	std::vector<int> ring;
	ring.push_back(base_vertex);
	
	for (int k = 0; k < vertices.size(); ++k) {
		for (int j = 0; j < nei1.size(); ++j) {
			if (check_shift(base_vertex, nei1[j], vertices[k])) {
				ring.push_back(nei1[j]);
			}
		}
	}

	assert(ring.size() == vertices.size() + 1);
	return ring;
}

void Lattice::set_rings() {
	std::vector<int> nei1(3);
	int ring_index = 0, s1, s2, s3;
	vec3<double> down1, down2;
	switch (lat_type){
		case KAGOME:
			rings1 = RingList(3, 2 * N / 3, N);
			//even rings are the three basis sites (up triangles)
			//odd rings are site 2(mod 3) + {basis[0], basis[2], basis[2] - basis[1]}
			down1 = basis[2];
			down2 = basis[2] - basis[1];
			for (int i = 0; i < N; i = i + 3) {
				rings1.set_ring(ring_index, i, i + 1, i + 2);
				s1 = i + 2;
				//cycle neighbors.  At each one check for down1 then down2
				nei1 = get_neighbors(s1, 0);
				for (int j = 0; j < nei1.size(); ++j) {
					if (check_shift(s1, nei1[j], down1)) {
						s2 = nei1[j];
					}
					else if (check_shift(s1, nei1[j], down2)) {
						s3 = nei1[j];
					}
				}
				rings1.set_ring(ring_index + 1, s1, s2, s3);
				ring_index += 2;
			}
			std::cout << "Ring list not fully implemented yet for " << Lattice_type_to_string(lat_type) << " lattice.\n";
			break;
		case TRIANGLE:
			rings1 = RingList(3, 2 * N, N);
			//even rings are up triangles (with neighbors at 0, a1, and a2)
			//odd rings are down triangles composed of base site, a2, and a2-a1
			for (int i = 0; i < N; ++i) {
				//up triangles
				down1 = a.x;
				down2 = a.y;
				nei1 = get_ring_sites(i, { down1, down2 });
				rings1.set_ring(ring_index, nei1[0], nei1[1], nei1[2]);

				//down triangles
				down1 = a.y;
				down2 = a.y - a.x;
				nei1 = get_ring_sites(i, { down1, down2 });
				rings1.set_ring(ring_index+1, nei1[0], nei1[1], nei1[2]);

				ring_index += 2;
			}
			break;
		default:
			std::cout << "Ring list not implemented yet for " << Lattice_type_to_string(lat_type) << " lattice.\n";
			break;
	}
}

void Lattice::set_distances() {
	double distance;
	bool found = false;
	for (int i = 1; i < N; ++i) {
		distance = nearest_periodic_distance(0, i);
		if (distances.size() == 0) {
			distances.push_back(distance);
		}
		else {
			for (int j = 0; j < distances.size(); ++j) {
				found = found || (std::abs(distances[j] - distance) < EPSILON);
			}
			if (!found) {
				distances.push_back(distance);
			}
			else {
				found = false;
			}
		}
	}
	std::sort(distances.begin(), distances.end());
}

double Lattice::nearest_periodic_distance(int reference_site, int other_site) {
	//cycle through all periodic images of "other_site" and return the distance to the closest image
	vec3<double> reference_position = coordinates[reference_site], other_position = coordinates[other_site];
	double distance = reference_position ^ other_position;
    double new_distance;
	for (int i = -1; i <= 1; ++i) {
		for (int j = -1; j <= 1; ++j) {
			for (int k = -1; k <= 1; ++k) {
                new_distance = reference_position ^ (other_position + (a.x*i*L.x) + (a.y*j*L.y) + (a.z*k*L.z));
                if(new_distance < distance){
                    distance = new_distance;
                }
			}
		}
	}

	return distance;
}

vec3<double> Lattice::nearest_periodic_translation(int reference_site, int other_site) {
	//cycle through all periodic images of "other_site" and return the distance to the closest image
	vec3<double> reference_position = coordinates[reference_site], other_position = coordinates[other_site], result;
	double distance = reference_position ^ other_position;
	double new_distance;
	result = other_position - reference_position;
	for (int i = -1; i <= 1; ++i) {
		for (int j = -1; j <= 1; ++j) {
			for (int k = -1; k <= 1; ++k) {
				new_distance = reference_position ^ (other_position + (a.x*i*L.x) + (a.y*j*L.y) + (a.z*k*L.z));
				if (new_distance < distance) {
					distance = new_distance;
					result = (other_position + (a.x*i*L.x) + (a.y*j*L.y) + (a.z*k*L.z)) - reference_position;
				}
			}
		}
	}

	return result;
}

double Lattice::nearest_periodic_distance(int reference_site, int other_site, vec3<int>& pbcs) {
	//cycle through all periodic images of "other_site" and return the distance to the closest image
	vec3<double> reference_position = coordinates[reference_site], other_position = coordinates[other_site];
	double distance = reference_position ^ other_position;
    double new_distance;
	for (int i = -1; i <= 1; ++i) {
		for (int j = -1; j <= 1; ++j) {
			for (int k = -1; k <= 1; ++k) {
                new_distance = reference_position ^ (other_position + (a.x*i*L.x) + (a.y*j*L.y) + (a.z*k*L.z));
                if(new_distance <= distance){
                    distance = new_distance;
                    pbcs.x = i;
                    pbcs.y = j;
                    pbcs.z = k;
                }
			}
		}
	}

	return distance;
}

void Lattice::set_neighbors(int max_distance) {
	double neighbor_distance;
    vec3<int> lt;//lattice translations
    int flipx = 1, flipy = 1, flipz = 1;
	for (int i = 0; i < N; ++i) {
		neighbors.push_back({});
        pbc.push_back({});
		for (int neighbor_number = 0; neighbor_number < fmin(distances.size(), max_distance); ++neighbor_number) {
			neighbors[i].push_back({});
            pbc[i].push_back({});
			neighbor_distance = distances[neighbor_number];
			for (int other_site = 0; other_site < N; ++other_site) {
				if (std::abs(nearest_periodic_distance(i, other_site, lt) - neighbor_distance ) < EPSILON){
					neighbors[i][neighbor_number].push_back(other_site);
                    flipx = abs(lt.x) == 0 ? 1 : pbc_signs.x;
					flipy = abs(lt.y) == 0 ? 1 : pbc_signs.y;
					flipz = abs(lt.z) == 0 ? 1 : pbc_signs.z;
					pbc[i][neighbor_number].push_back(flipx*flipy*flipz);
					assert(flipx*flipy*flipz == 1 || flipx*flipy*flipz == -1);
				}
			}
		}
	}
}

void Lattice::set_neighbors_faster() {
	double neighbor_distance;
    vec3<int> lt;//lattice translations
    int flipx = 1, flipy = 1, flipz = 1;
    auto distance_it = distances.begin();
	
	//for ordering
	struct neighbor_entry {
		int site, pbc;
		double angle;
		double z_coord;
	};
	neighbor_entry tempneigh_entry;
	std::vector<int> tempneighbors, temppbcs;
	vec3<double> direction;
	std::vector<neighbor_entry> tempneigh_list;


	//set basis site neighbors first
	for (int i = 0; i < N; ++i) {
		neighbors.push_back({});
        pbc.push_back({});
		for (int neighbor_number = 0; neighbor_number < distances.size(); ++neighbor_number) {
			neighbors[i].push_back({});
            pbc[i].push_back({});
		}
        for (int other_site = 0; other_site < N; ++other_site) {
            if(other_site != i){
                neighbor_distance = nearest_periodic_distance(i, other_site, lt);
                distance_it = std::lower_bound(distances.begin(), distances.end(), neighbor_distance - EPSILON);//should point to the element giving neighbor_distance
                neighbors[i][distance_it-distances.begin()].push_back(other_site);
                flipx = std::abs(lt.x) == 0 ? 1 : pbc_signs.x;
                flipy = std::abs(lt.y) == 0 ? 1 : pbc_signs.y;
                flipz = std::abs(lt.z) == 0 ? 1 : pbc_signs.z;
                pbc[i][distance_it-distances.begin()].push_back(flipx*flipy*flipz);
                assert(flipx*flipy*flipz == 1 || flipx*flipy*flipz == -1);
            }
        }
		//ensure proper ordering of neighbors
		for (int neighbor_number = 0; neighbor_number < distances.size(); ++neighbor_number) {
			tempneighbors = neighbors[i][neighbor_number];
			temppbcs = pbc[i][neighbor_number];
			tempneigh_list.clear();
			for (int n = 0; n < tempneighbors.size(); ++n) {
				direction = nearest_periodic_translation(i, tempneighbors[n]);
				tempneigh_entry.site = tempneighbors[n];
				tempneigh_entry.pbc = temppbcs[n];
				tempneigh_entry.angle = direction.angle_xy();
				tempneigh_entry.z_coord = direction.z;
				tempneigh_list.push_back(tempneigh_entry);
			}
			std::sort(tempneigh_list.begin(), tempneigh_list.end(), [](auto const& a, auto const& b) { return std::abs(a.angle - b.angle) < EPSILON ? a.z_coord > b.z_coord : a.angle < b.angle; });
			for (int n = 0; n < neighbors[i][neighbor_number].size(); ++n) {
				neighbors[i][neighbor_number][n] = tempneigh_list[n].site;
				pbc[i][neighbor_number][n] = tempneigh_list[n].pbc;
			}
		}
	}

	////cycle through the rest of the sites adding based on corresponding basis site
	//int basis_corr, shift;
	//for (int i = basis.size(); i < N; ++i) {
	//	basis_corr = i % basis.size();
	//	shift = i - basis_corr;
	//	neighbors.push_back({});
	//	pbc.push_back({});
	//	for (int neighbor_number = 0; neighbor_number < distances.size(); ++neighbor_number) {
	//		neighbors[i].push_back({});
	//		pbc[i].push_back({});
	//	}
	//}
}

void Lattice::set_coordinates() {
	for (int x = 0; x < L.x; ++x) {
		for (int y = 0; y < L.y; ++y) {
			for (int z = 0; z < L.z; ++z) {
				for (int base = 0; base < basis.size(); ++base) {
					coordinates.push_back(a.x*x + a.y*y + a.z*z + basis[base]);
				}
			}
		}
	}
}

void Lattice::set_reciprocal_vectors() {
	double volume = a.x * (a.y % a.z);
	b.x =  (a.y % a.z) * 2 * M_PI * (1 / volume);
	b.y =  (a.z % a.x) * 2 * M_PI * (1 / volume);
	b.z =  (a.x % a.y) * 2 * M_PI * (1 / volume);

	assert(abs(2 * M_PI - a.x*b.x) < EPSILON);
	assert(abs(2 * M_PI - a.y*b.y) < EPSILON);
	assert(abs(2 * M_PI - a.z*b.z) < EPSILON);
}

bool Lattice::check_in_plane(int i1, int i2){
    vec3<double> shift = get_coordinate(i2) - get_coordinate(i1), new_shift(0,0,0);
    bool result = false;
	for (int i = -1; i <= 1; ++i) {
		for (int j = -1; j <= 1; ++j) {
			for (int k = -1; k <= 1; ++k) {
				new_shift = shift + (a.x*i*L.x) + (a.y*j*L.y) + (a.z*k*L.z);
				result = result || (std::abs(new_shift.z) < EPSILON);
			}
		}
	}
    return result;
}

std::vector<std::vector<int>> Lattice::basis_partition(std::vector<vec3<int>> basis) {
	std::vector<std::vector<int>> partition;
	std::vector<int> unit_cell;

	if (basis.size() == 0 || basis.size() == 1) {
		// return the partition of this lattice's basis
		int num_lattice_basis = this->basis.size();
		int site = 0;
		for (int uc = 0; uc < (this->N) / num_lattice_basis; ++uc) {
			partition.push_back({});
			for (int b = 0; b < num_lattice_basis; ++b) {
				partition[uc].push_back(site);
				++site;
			}
		}
	}
	else {
		//create list of all possible enlarged unit cells
		for (int i = 0; i < N; ++i) {
			unit_cell.clear();
			for (int b = 0; b < basis.size(); ++b) {
				for (int b0 = 0; b0 < this->basis.size(); ++b0) {
					if (b == 0 && b0 == 0) {
						unit_cell.push_back(i);
					}
					else if ( this->basis.size() * (i / this->basis.size() ) == i){
						unit_cell.push_back(get_neighbor_label_with_pbc(i, basis[b] * a + this->basis[b0]));
					}
				}
			}
			partition.push_back(unit_cell);
		}

		if (basis.size() > 1) {
			//remove unit cells that duplicate existing cells
			for (int i = 0; i < partition.size(); ++i) {
				if (partition[i][0] != N) {
					for (int b = 1; b < partition[i].size(); ++b) {
						partition[partition[i][b]][0] = N;//mark the unit cell that is already accounted in another unit cell
					}
				}
			}
			for (int i = partition.size() - 1; i >= 0; --i) {
				if (partition[i][0] == N) {
					partition.erase(partition.begin() + i);
				}
			}
		}
	}

	return partition;
}
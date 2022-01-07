#pragma once
#include <vector>
#include <complex>
#include <iostream>
#include "json.hpp"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// for convenience
using json = nlohmann::json;

namespace vmctype {


	/*! \brief Small value for double comparisons
		*/
	const double EPSILON = 1e-12;

	class HoppingTerm {
	public:
		int spin_row; //|m|: in spin-1 this is 0 for tz and 1 for txy
		double strength;
		int distance; //neighbor distance
		std::vector<int> origins, neighbor_index;//within unit cell
		std::vector<double> phases;

	};

	/*! \brief 2-D Vector
*/
	template <class T>
	class vec2 {
	public:
		T x, y;

		vec2() {}

		vec2(T x_in, T y_in) : x(x_in), y(y_in) {}
	};



	/*! \brief 3-D Vector with operations
	*/
	template <class T>
	class vec3 {
	public:
		T x, y, z;

		vec3() {}

		vec3(std::vector<T> v_in) {
			assert(v_in.size() == 3);
			x = v_in[0];
			y = v_in[1];
			z = v_in[2];
		}

		vec3(T x_in, T y_in, T z_in) : x(x_in), y(y_in), z(z_in) {}

		double angle_xy() {
			//get the angle in rad in the xy plane (from the x-axis)
			if (y == 0) {
				return x > 0 ? 0.0 : M_PI;
			}
			else {
				return atan2(y, x) < 0 ? atan2(y, x) + 2.0*M_PI : atan2(y, x);//[0,2pi]
			}
		}

		vec3 operator+ (const vec3<T>& v) const {
			return vec3(this->x + v.x, this->y + v.y, this->z + v.z);
		}

		vec3 operator* (const int& a) const {
			return vec3(this->x * a, this->y * a, this->z * a);
		}

		vec3 operator* (const double& a) const {
			return vec3(this->x * a, this->y * a, this->z * a);
		}

		template <class R>
		R operator* (const vec3<R>& v) const {
			return (v.x * this->x + v.y * this->y + v.z * this->z);
		}

		vec3 operator% (const vec3<T>& v) const {
			//cross product
			return vec3(this->y * v.z - this->z * v.y, this->z * v.x - this->x * v.z, this->x * v.y - this->y * v.x);
		}

		vec3 operator- (const vec3<T>& v) const {
			return vec3(this->x - v.x, this->y - v.y, this->z - v.z);
		}

		double operator^ (const vec3<T>& v) const {
			//distance between the two vectors
			vec3 diff(this->x - v.x, this->y - v.y, this->z - v.z);
			return sqrt(diff * diff);
		}

		bool operator== (const vec3<T>& v) const {
			//compare vectors
			return (std::abs(this->x - v.x) < EPSILON
				&& std::abs(this->y - v.y) < EPSILON
				&& std::abs(this->z - v.z) < EPSILON);
		}

	};

	template <class T>
	std::string vec3str(vec3<T> vec_in) {
		std::stringstream ss;
		ss << "<" << vec_in.x << "," << vec_in.y << "," << vec_in.z << ">";
		return ss.str();
	}

	template <class T>
	void to_json(json& j, const vec3<T>& p) {
		j = json{ vec3str(p) };
	}

	template <class T>
	void from_json(const json& j, vec3<T>& p) {
		std::vector<T> vec(3);
		j.get_to<std::vector<T>>(vec);
		p.x = vec[0];
		p.y = vec[1];
		p.z = vec[2];
	}

	class QuadrupoleOrder {
		//As defined here, polar coordinates are in degrees and Q momenta are in full units (including factors of 2pi)
	public:
		std::vector<vec3<double>> unit_cell_u_polar, unit_cell_v_polar;
		vec3<double> Q_t_u, Q_p_u, Q_t_v, Q_p_v;

		vec3<std::complex<double>> eval_at(int uc_idx, vec3<double> position) {
			double r_u, r_v, theta_u, theta_v, phi_u, phi_v;
			r_u = unit_cell_u_polar[uc_idx].x;
			theta_u = M_PI / 180.0 * (unit_cell_u_polar[uc_idx].y) + Q_t_u * position;
			phi_u = M_PI / 180.0 * (unit_cell_u_polar[uc_idx].z) + Q_p_u * position;

			r_v = unit_cell_v_polar[uc_idx].x;
			theta_v = M_PI / 180.0 * (unit_cell_v_polar[uc_idx].y) + Q_t_v * position;
			phi_v = M_PI / 180.0 * (unit_cell_v_polar[uc_idx].z) + Q_p_v * position;

			double norm = r_u * r_u + r_v * r_v;
			r_u = r_u / sqrt(norm);
			r_v = r_v / sqrt(norm);
			return vec3<std::complex<double>>({ r_u * cos(phi_u) * sin(theta_u), r_v * cos(phi_v) * sin(theta_v) },
				{ r_u * sin(phi_u) * sin(theta_u), r_v * sin(phi_v) * sin(theta_v) },
				{ r_u * cos(theta_u), r_v * cos(theta_v) });
		}
	};

	struct JastrowFactorOptions {

		bool isotropic = true;
		int distance_max = -1;
		std::vector<double> values = {};

	};

	struct JastrowTableOptions {

		JastrowFactorOptions sz;

	};

	struct BilinearOptions {

		std::string interaction_type = "su2 exchange";
		double coupling = 0.0;
		int neighbor_index = 0;

	};

	struct TrilinearOptions {

		std::string interaction_type = "su3 ring exchange";
		bool hermitian = true;
		double coupling_real;
		double coupling_imag;

	};

	struct lattice_options {
		std::string type;
		int dimension;
		std::vector<int> L;
		std::vector<int> pbc;
		int Lx;
		int Ly;
		int Lz;

		std::string to_string();
		bool is_valid();
	};

	struct mean_field_options {
		std::string lattice_type;
		std::string wf_type;
		int inequivalent_sites = 1, spin, num_spin_orbit;
		bool match_lattice_pbc, su3_symmetry;
		double field = 0.0;
		std::vector<vec3<int>> basis;
		std::vector<vmctype::HoppingTerm> hopping_list;
		QuadrupoleOrder directors;
		bool jastrow_flag = false;
		JastrowTableOptions jastrow;

		std::string to_string();
	};

	struct model_options {
		std::vector<BilinearOptions> bilinear_terms;
		std::vector<TrilinearOptions> ring3_terms;
	};

	struct vmc_options {
		int steps_per_measure, num_measures, throwaway_measures;
		bool optimization = false;
	};

	void to_json(json& j, const HoppingTerm& p);

	void from_json(const json& j, HoppingTerm& p);

	void to_json(json& j, const QuadrupoleOrder& p);

	void from_json(const json& j, QuadrupoleOrder& p);

	void to_json(json& j, const JastrowFactorOptions& p);

	void from_json(const json& j, JastrowFactorOptions& p);

	void to_json(json& j, const BilinearOptions& p);

	void from_json(const json& j, BilinearOptions& p);

	void to_json(json& j, const TrilinearOptions& p);

	void from_json(const json& j, TrilinearOptions& p);

}
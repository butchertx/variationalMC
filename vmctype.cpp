#include "vmctype.h"

void vmctype::to_json(json& j, const vmctype::HoppingTerm& p) {
	j = json{
			{"spin row", p.spin_row},
			{"strength", p.strength},
			{"distance", p.distance},
			{"origins", p.origins},
			{"neighbor index", p.neighbor_index},
			{"phases", p.phases},
	};
}

void vmctype::from_json(const json& j, vmctype::HoppingTerm& p) {
	j.at("spin row").get_to(p.spin_row);
	j.at("strength").get_to(p.strength);
	j.at("distance").get_to(p.distance);
	j.at("origins").get_to<std::vector<int>>(p.origins);
	j.at("neighbor index").get_to<std::vector<int>>(p.neighbor_index);
	j.at("phases").get_to<std::vector<double>>(p.phases);
}

void vmctype::to_json(json& j, const vmctype::QuadrupoleOrder& p) {
	j = json{
			{"u_r_theta_phi", p.unit_cell_u_polar},
			{"v_r_theta_phi", p.unit_cell_v_polar},
			{"Qthetau", p.Q_t_u},
			{"Qphiu", p.Q_p_u},
			{"Qthetav", p.Q_t_v},
			{"Qphiv", p.Q_p_v},
	};
}

void vmctype::from_json(const json& j, vmctype::QuadrupoleOrder& p) {
	j.at("u_r_theta_phi").get_to<std::vector<vec3<double>>>(p.unit_cell_u_polar);
	j.at("v_r_theta_phi").get_to<std::vector<vec3<double>>>(p.unit_cell_v_polar);
	j.at("Qthetau").get_to<vec3<double>>(p.Q_t_u);
	j.at("Qphiu").get_to<vec3<double>>(p.Q_p_u);
	j.at("Qthetav").get_to<vec3<double>>(p.Q_t_v);
	j.at("Qphiv").get_to<vec3<double>>(p.Q_p_v);
}
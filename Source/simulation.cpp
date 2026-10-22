#include "sapling/simulation.hpp"
#include "Simulations/simulations.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace sapling {

bool Simulation::set_parameter(const std::string& id, double value) {
    for (Parameter& p : m_params) {
        if (p.id != id) continue;
        if (!std::isfinite(value)) return false;
        value = std::min(p.max, std::max(p.min, value));
        if (p.value == value) return true;
        p.value = value;
        parameter_changed(id);
        return true;
    }
    return false;
}

double Simulation::parameter(const std::string& id) const {
    for (const Parameter& p : m_params) {
        if (p.id == id) return p.value;
    }
    return 0.0;
}

void Simulation::add_parameter(const std::string& id, const std::string& label, double value, double min, double max, double step, const std::string& unit) {
    m_params.push_back(Parameter{ id, label, value, min, max, step, unit });
}

void Simulation::parameter_changed(const std::string&) {
    reset();
}

const std::vector<Section>& sections() {
    static const std::vector<Section> all = {
        { "mechanics", "Mechanics", "Pendulums, springs, rigid bodies and everything that moves under forces" },
        { "gravity", "Gravity", "Orbits in Newton's gravity and around black holes in general relativity" },
        { "quantum", "Quantum", "Wave functions, tunneling, bound states, atoms and molecules" },
        { "thermo", "Thermodynamics", "Heat flow, temperature, engines and statistical mechanics" },
        { "fluids", "Fluids and aero", "Flow around objects, lift, vortices and waves" },
        { "electro", "Electromagnetism", "Charges, fields, currents, magnets and light" }
    };
    return all;
}

const std::vector<SimulationInfo>& simulations() {
    static const std::vector<SimulationInfo> all = {
        { "double_pendulum", "mechanics", "Double pendulum", "Two rods hinged end to end: a tiny change in the start gives a completely different swing." },
        { "kepler", "gravity", "Kepler orbit", "A planet around a star in Newton's gravity, stepped with a symplectic integrator so the orbit never drifts." },
        { "black_hole", "gravity", "Orbit around a black hole", "A geodesic of the Schwarzschild metric: the orbit's closest point turns a little every lap, and close in it plunges." },
        { "tunneling", "quantum", "Quantum tunneling", "A wave packet hits a barrier taller than its energy, and part of it gets through." },
        { "heat_rod", "thermo", "Heat in a rod", "A rod with its ends held at two temperatures, settling toward a straight line." },
        { "cylinder_flow", "fluids", "Flow past a cylinder", "Ideal flow around a cylinder with circulation: the streamlines and the lift it makes." },
        { "field_lines", "electro", "Electric field lines", "The field of a few point charges, traced from each positive charge." }
    };
    return all;
}

const SimulationInfo* find_simulation(const std::string& id) {
    for (const SimulationInfo& s : simulations()) {
        if (s.id == id) return &s;
    }
    return nullptr;
}

std::unique_ptr<Simulation> make_simulation(const std::string& id) {
    std::unique_ptr<Simulation> s;
    if (id == "double_pendulum") s = make_double_pendulum();
    else if (id == "kepler") s = make_kepler();
    else if (id == "black_hole") s = make_black_hole();
    else if (id == "tunneling") s = make_tunneling();
    else if (id == "heat_rod") s = make_heat_rod();
    else if (id == "cylinder_flow") s = make_cylinder_flow();
    else if (id == "field_lines") s = make_field_lines();
    if (s) s->reset();
    return s;
}

std::string format_number(double v, int digits) {
    if (!std::isfinite(v)) return v != v ? "undefined" : (v > 0 ? "\xE2\x88\x9E" : "-\xE2\x88\x9E");
    if (std::fabs(v) < 1e-300) v = 0.0;
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*g", digits, v);
    return buf;
}

} // namespace sapling

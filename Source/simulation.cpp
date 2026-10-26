#include "sapling/simulation.hpp"

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

std::string format_number(double v, int digits) {
    if (!std::isfinite(v)) return v != v ? "undefined" : (v > 0 ? "\xE2\x88\x9E" : "-\xE2\x88\x9E");
    if (std::fabs(v) < 1e-300) v = 0.0;
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*g", digits, v);
    return buf;
}

} // namespace sapling

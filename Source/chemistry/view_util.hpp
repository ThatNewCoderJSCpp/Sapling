#ifndef SAPLING_CHEMISTRY_VIEW_UTIL_HPP
#define SAPLING_CHEMISTRY_VIEW_UTIL_HPP

#include "sapling/chemistry.hpp"

#include <fizmo/../chemistry.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string>
#include <vector>

namespace sapling {
namespace chem_util {

namespace fc = fizmo::chemistry;

inline std::string arg(const std::vector<std::string>& a, std::size_t i) {
    return i < a.size() ? a[i] : std::string();
}

inline double num(const std::string& s, double fallback) {
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    return end == s.c_str() || !std::isfinite(v) ? fallback : v;
}

inline bool finite(double v) {
    return v == v && std::isfinite(v);
}

inline std::string fixed(double v, int digits) {
    if (!finite(v)) return "n/a";
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*f", digits, v);
    return buf;
}

inline std::string trim(std::string s) {
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r' || s.back() == '\n')) s.pop_back();
    std::size_t i = 0;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
    return s.substr(i);
}

inline std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == sep) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    out.push_back(cur);
    return out;
}

inline const fc::Element* element_of(const std::string& text) {
    std::string t = trim(text);
    const std::size_t sp = t.find(' ');
    if (sp != std::string::npos) t = t.substr(0, sp);
    return t.empty() ? nullptr : fc::find_element(t);
}

inline double to_si_or(double value, const std::string& unit, fc::Dimension d, double fallback) {
    try {
        return fc::to_si(value, unit, d);
    } catch (...) {
        return fallback;
    }
}

inline double from_si_or(double value, const std::string& unit, fc::Dimension d) {
    try {
        return fc::from_si(value, unit, d);
    } catch (...) {
        return value;
    }
}

Color heat_color(double t, const Theme& theme);
Color signed_color(double v, const Theme& theme);
std::optional<fc::Molecule> molecule_from_text(const std::string& text, int* charge, std::string* name, std::string* smiles);
fc::PhaseProperties phase_from_saved(const std::string& kind, const std::string& json, std::string* error);
std::optional<fc::IonDescriptor> ion_from_text(const std::string& text, bool cation);
std::string dielectric_readouts(const fc::DielectricProperties& d);
std::string mechanics_readouts(const fc::MechanicalProperties& m);
std::string electron_gas_readouts(const fc::ElectronGasProperties& g, const UnitPrefs& u);
std::string susceptibility_readouts(const fc::SusceptibilityInfo& s);

} // namespace chem_util
} // namespace sapling

#endif // SAPLING_CHEMISTRY_VIEW_UTIL_HPP

#include "view_util.hpp"

#include <algorithm>

namespace sapling {
namespace chem_util {

Color heat_color(double t, const Theme& theme) {
    t = std::max(0.0, std::min(1.0, t));
    const Color stops[] = { theme.background, rgb(49, 50, 120), theme.accent, rgb(166, 227, 255), rgb(255, 255, 255) };
    const double x = t * 4.0;
    const int i = std::min(3, static_cast<int>(x));
    return mix(stops[i], stops[i + 1], x - i);
}

Color signed_color(double v, const Theme& theme) {
    const double a = std::min(1.0, std::fabs(v));
    const Color hi = v >= 0.0 ? rgb(243, 139, 168) : rgb(137, 180, 250);
    return mix(theme.background, hi, std::pow(a, 0.6));
}

std::optional<fc::IonDescriptor> ion_from_text(const std::string& text, bool cation) {
    const std::string t = trim(text);
    if (t.empty()) return std::nullopt;
    if (const fc::Element* e = fc::find_element(t)) {
        if (cation) {
            const std::vector<int> charges = fc::likely_cation_charges(*e);
            if (charges.empty()) return std::nullopt;
            return fc::monatomic_ion(*e, charges.front());
        }
        const int q = fc::likely_anion_charge(*e);
        if (q >= 0) return std::nullopt;
        return fc::monatomic_ion(*e, q);
    }
    try {
        std::optional<fc::IonDescriptor> d = fc::parse_ion(t);
        if (d && d->valid() && (d->charge > 0) == cation) return d;
    } catch (...) {
    }
    return std::nullopt;
}

std::optional<fc::Molecule> molecule_from_text(const std::string& text, int* charge, std::string* name, std::string* smiles) {
    const std::string t = trim(text);
    if (t.empty()) return std::nullopt;
    std::string s;
    int q = 0;
    std::string n;
    if (const fc::PolyatomicIonEntry* ion = fc::find_polyatomic_ion(t)) {
        s = ion->smiles;
        q = ion->charge;
        n = ion->name;
    }
    if (s.empty()) {
        const auto found = fc::find_common_compounds(t);
        for (const fc::CommonCompoundEntry* c : found) {
            if (c->smiles && *c->smiles) {
                s = c->smiles;
                n = c->name;
                break;
            }
        }
    }
    if (s.empty()) s = t;
    try {
        fc::Molecule m = fc::Molecule::from_smiles(s);
        if (m.atom_count() == 0) return std::nullopt;
        if (q == 0) q = m.charge();
        if (charge) *charge = q;
        if (name) *name = n;
        if (smiles) *smiles = s;
        return m;
    } catch (...) {
        return std::nullopt;
    }
}

fc::PhaseProperties phase_from_saved(const std::string& kind, const std::string& json, std::string* error) {
    fc::PhaseProperties props;
    JsonValue v;
    if (!parse_json(json, v, error)) return props;
    try {
        if (kind == "covalent") {
            fc::Molecule m;
            const JsonValue* atoms = v.find("atoms");
            const JsonValue* bonds = v.find("bonds");
            if (atoms && atoms->type == JsonValue::Array) {
                for (const JsonValue& a : atoms->array) {
                    const fc::Element* e = fc::find_element(a.string_or("element", "C"));
                    if (!e) continue;
                    const std::size_t i = m.add_atom(*e);
                    m.set_position(i, a.number_or("x", 0.0), a.number_or("y", 0.0), a.number_or("z", 0.0));
                }
            }
            if (bonds && bonds->type == JsonValue::Array) {
                for (const JsonValue& b : bonds->array) {
                    const int x = static_cast<int>(b.number_or("a", -1)), y = static_cast<int>(b.number_or("b", -1));
                    const int o = std::max(1, std::min(3, static_cast<int>(b.number_or("order", 1))));
                    if (x < 0 || y < 0 || x == y || x >= static_cast<int>(m.atom_count()) || y >= static_cast<int>(m.atom_count())) continue;
                    if (m.find_bond(static_cast<std::size_t>(x), static_cast<std::size_t>(y))) continue;
                    m.add_bond(static_cast<std::size_t>(x), static_cast<std::size_t>(y), static_cast<fc::BondOrder>(o));
                }
            }
            if (m.atom_count() == 0) {
                if (error) *error = "the saved molecule is empty";
                return props;
            }
            const std::string formula = m.formula().hill();
            props = fc::find_phase_properties(formula);
            if (!props.usable()) props = fc::phase_properties(m, v.string_or("name", formula));
        } else if (kind == "ionic") {
            const auto c = ion_from_text(v.string_or("cation", ""), true);
            const auto a = ion_from_text(v.string_or("anion", ""), false);
            if (!c || !a) {
                if (error) *error = "the saved ionic compound has ions I can't read";
                return props;
            }
            props = fc::phase_properties(fc::ionic_crystal(*c, *a));
        } else if (kind == "metallic") {
            std::vector<fc::AlloyComponent> comps;
            const JsonValue* list = v.find("components");
            if (list && list->type == JsonValue::Array) {
                for (const JsonValue& c : list->array) {
                    const fc::Element* e = fc::find_element(c.string_or("element", ""));
                    if (e) comps.push_back(fc::AlloyComponent{ e, c.number_or("fraction", 1.0) });
                }
            }
            if (comps.empty()) {
                if (error) *error = "the saved metal has no components";
                return props;
            }
            props = fc::phase_properties(fc::metallic_crystal(comps));
        } else {
            if (error) *error = "unknown kind " + kind;
            return props;
        }
        if (!props.usable()) fc::complete_phase_properties(props);
    } catch (const std::exception& e) {
        if (error) *error = e.what();
    }
    return props;
}

} // namespace chem_util
} // namespace sapling

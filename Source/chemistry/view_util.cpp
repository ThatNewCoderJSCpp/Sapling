#include "view_util.hpp"

#include <algorithm>
#include <sstream>
#include <string>

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

namespace {

std::string approx(bool estimated) {
    return estimated ? "\xE2\x89\x88 " : "";
}

} // namespace

std::string dielectric_readouts(const fc::DielectricProperties& d) {
    std::ostringstream o;
    o << "#Permittivity and light\n";
    const std::string a = approx(d.permittivity_estimated);
    if (d.conductor) {
        o << "Relative permittivity \xCE\xB5\xE1\xB5\xA3\teffectively infinite (a conductor)\n";
    } else if (finite(d.relative_permittivity)) {
        o << "Relative permittivity \xCE\xB5\xE1\xB5\xA3\t" << a << format_number(d.relative_permittivity, 5) << "\n";
        if (finite(d.absolute_permittivity)) o << "Permittivity \xCE\xB5 = \xCE\xB5\xE1\xB5\xA3\xCE\xB5\xE2\x82\x80\t" << a << format_number(d.absolute_permittivity, 4) << " F/m\n";
        if (finite(d.susceptibility)) o << "Electric susceptibility \xCF\x87\xE2\x82\x91\t" << a << format_number(d.susceptibility, 4) << "\n";
        o << "Charges inside it attract with 1/" << format_number(d.relative_permittivity, 3) << " of the force they would in vacuum\n";
    }
    if (finite(d.high_frequency_permittivity) && !d.conductor) o << "Optical permittivity \xCE\xB5\xE2\x88\x9E\t" << approx(d.refractive_index_estimated || d.permittivity_estimated) << format_number(d.high_frequency_permittivity, 4) << "\n";
    if (finite(d.optical_permittivity)) o << "Permittivity at 550 nm\t" << format_number(d.optical_permittivity, 4) << " (negative: light can't propagate inside)\n";
    if (finite(d.refractive_index)) o << "Refractive index n\t" << approx(d.refractive_index_estimated) << format_number(d.refractive_index, 4) << "\n";
    if (finite(d.speed_of_light)) o << "Speed of light inside\t" << fixed(d.speed_of_light / 1000.0, 0) << " km/s (" << format_number(d.speed_of_light / 2.99792458e8 * 100.0, 4) << "% of c)\n";
    if (finite(d.reflectance)) o << "Reflects at normal incidence\t" << format_number(d.reflectance * 100.0, 3) << "%\n";
    if (finite(d.brewster_angle)) o << "Brewster angle\t" << format_number(d.brewster_angle, 4) << "\xC2\xB0\n";
    if (finite(d.critical_angle)) o << "Critical angle (total internal reflection)\t" << format_number(d.critical_angle, 4) << "\xC2\xB0\n";
    if (finite(d.optical_penetration_depth)) o << "Light gets in about\t" << format_number(d.optical_penetration_depth * 1e9, 3) << " nm\n";
    if (finite(d.polarizability_volume)) o << "Polarizability volume\t" << format_number(d.polarizability_volume, 4) << " \xC3\x85\xC2\xB3\n";
    if (finite(d.polarizability)) o << "Polarizability \xCE\xB1\t" << format_number(d.polarizability, 4) << " C\xC2\xB7m\xC2\xB2/V\n";
    if (finite(d.molar_refraction)) o << "Molar refraction\t" << format_number(d.molar_refraction * 1e6, 4) << " cm\xC2\xB3/mol\n";
    if (finite(d.dielectric_strength)) o << "Dielectric strength\t" << format_number(d.dielectric_strength / 1e6, 3) << " MV/m (it breaks down and sparks above this)\n";
    if (!d.method.empty()) o << "Method\t" << d.method << "\n";
    for (const std::string& n : d.notes) o << n << "\n";
    return o.str();
}

std::string mechanics_readouts(const fc::MechanicalProperties& m) {
    std::ostringstream o;
    if (!finite(m.youngs_modulus) && !finite(m.bulk_modulus) && !finite(m.linear_expansion) && !finite(m.mohs_hardness)) return std::string();
    o << "#Stiffness, sound and expansion\n";
    const std::string a = approx(m.estimated);
    if (finite(m.youngs_modulus)) o << "Young's modulus\t" << a << format_number(m.youngs_modulus, 4) << " GPa\n";
    if (finite(m.bulk_modulus)) o << "Bulk modulus\t" << a << format_number(m.bulk_modulus, 4) << " GPa\n";
    if (finite(m.shear_modulus)) o << "Shear modulus\t" << a << format_number(m.shear_modulus, 4) << " GPa\n";
    if (finite(m.poisson_ratio)) o << "Poisson's ratio\t" << a << format_number(m.poisson_ratio, 3) << "\n";
    if (finite(m.compressibility)) o << "Compressibility\t" << a << format_number(m.compressibility * 1e9, 4) << " per GPa\n";
    if (finite(m.longitudinal_sound_speed)) o << "Speed of sound (pressure wave)\t" << a << format_number(m.longitudinal_sound_speed, 4) << " m/s\n";
    if (finite(m.transverse_sound_speed)) o << "Speed of sound (shear wave)\t" << a << format_number(m.transverse_sound_speed, 4) << " m/s\n";
    if (finite(m.bar_sound_speed)) o << "Speed of sound in a thin rod\t" << a << format_number(m.bar_sound_speed, 4) << " m/s\n";
    if (finite(m.acoustic_impedance)) o << "Acoustic impedance\t" << a << format_number(m.acoustic_impedance / 1e6, 4) << " MRayl\n";
    if (finite(m.linear_expansion)) o << "Thermal expansion (linear)\t" << a << format_number(m.linear_expansion * 1e6, 4) << " \xC2\xB5m/(m\xC2\xB7K)\n";
    if (finite(m.volume_expansion)) o << "Thermal expansion (volume)\t" << a << format_number(m.volume_expansion, 4) << " per K\n";
    if (finite(m.gruneisen)) o << "Gr\xC3\xBCneisen parameter\t" << a << format_number(m.gruneisen, 3) << "\n";
    if (finite(m.mohs_hardness)) o << "Mohs hardness\t" << a << format_number(m.mohs_hardness, 2) << "\n";
    for (const std::string& n : m.notes) o << n << "\n";
    return o.str();
}

std::string electron_gas_readouts(const fc::ElectronGasProperties& g, const UnitPrefs& u) {
    std::ostringstream o;
    if (!finite(g.fermi_energy) && !finite(g.work_function)) return std::string();
    o << "#Free electrons\n";
    if (finite(g.fermi_energy)) o << "Fermi energy\t" << format_number(g.fermi_energy, 4) << " eV\n";
    if (finite(g.fermi_temperature)) o << "Fermi temperature\t" << u.format(g.fermi_temperature, "temperature", 4) << "\n";
    if (finite(g.fermi_velocity)) o << "Fermi velocity\t" << format_number(g.fermi_velocity / 1000.0, 4) << " km/s\n";
    if (finite(g.fermi_wavevector)) o << "Fermi wavevector\t" << format_number(g.fermi_wavevector * 1e-9, 4) << " nm\xE2\x81\xBB\xC2\xB9\n";
    if (finite(g.density_of_states)) o << "Density of states at E_F\t" << format_number(g.density_of_states * 1.602176634e-19 * 1e-27, 4) << " per eV per nm\xC2\xB3\n";
    if (finite(g.work_function)) o << "Work function\t" << format_number(g.work_function, 3) << " eV\n";
    if (finite(g.photoelectric_threshold)) o << "Photoelectric threshold\t" << format_number(g.photoelectric_threshold * 1e9, 4) << " nm (shorter light knocks electrons out)\n";
    if (finite(g.richardson_current)) {
        if (g.richardson_current > 1e-12) o << "Thermionic emission\t" << format_number(g.richardson_current, 3) << " A/m\xC2\xB2 (Richardson)\n";
        else o << "Thermionic emission\tnone to speak of at this temperature (it takes well over 1000 K)\n";
    }
    if (finite(g.seebeck_coefficient)) o << "Seebeck coefficient\t" << format_number(g.seebeck_coefficient * 1e6, 3) << " \xC2\xB5V/K (free-electron estimate)\n";
    if (finite(g.thomas_fermi_length)) o << "Screening length\t" << format_number(g.thomas_fermi_length * 1e12, 3) << " pm (fields die out this fast inside)\n";
    if (finite(g.pauli_susceptibility)) o << "Pauli susceptibility\t" << format_number(g.pauli_susceptibility, 3) << "\n";
    return o.str();
}

std::string susceptibility_readouts(const fc::SusceptibilityInfo& s) {
    std::ostringstream o;
    if (!finite(s.volume_susceptibility) && !finite(s.molar_susceptibility_cgs)) return std::string();
    o << "#Magnetic susceptibility\n";
    const std::string a = approx(s.estimated);
    if (finite(s.volume_susceptibility)) {
        o << "Volume susceptibility \xCF\x87\t" << a << format_number(s.volume_susceptibility, 4) << (s.volume_susceptibility < 0 ? " (diamagnetic: pushed out of fields)" : " (paramagnetic: pulled into fields)") << "\n";
        o << "Relative permeability \xCE\xBC\xE1\xB5\xA3\t" << a << format_number(1.0 + s.volume_susceptibility, 8) << "\n";
    }
    if (finite(s.molar_susceptibility_cgs)) o << "Molar susceptibility\t" << a << format_number(s.molar_susceptibility_cgs * 1e6, 4) << "\xC3\x97" "10\xE2\x81\xBB\xE2\x81\xB6 cm\xC2\xB3/mol (cgs)\n";
    if (finite(s.mass_susceptibility)) o << "Mass susceptibility\t" << a << format_number(s.mass_susceptibility, 4) << " m\xC2\xB3/kg\n";
    if (s.unpaired_electrons > 0) o << "Unpaired electrons\t" << s.unpaired_electrons << "\n";
    if (finite(s.effective_moment)) o << "Effective moment\t" << format_number(s.effective_moment, 3) << " \xCE\xBC" "B (spin only)\n";
    if (finite(s.curie_constant)) o << "Curie constant\t" << format_number(s.curie_constant, 4) << " m\xC2\xB3\xC2\xB7K/mol\n";
    if (!s.method.empty()) o << "Method\t" << s.method << "\n";
    return o.str();
}

} // namespace chem_util
} // namespace sapling

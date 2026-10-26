#include "view_util.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <sstream>

namespace sapling {

using namespace chem_util;

namespace {

std::string clean(std::string s) {
    for (char& c : s) {
        if (c == '\t' || c == '\n' || c == '\r') c = ' ';
    }
    return s;
}

std::string join(const std::vector<std::string>& v, const std::string& sep, std::size_t max = 1000) {
    std::string out;
    for (std::size_t i = 0; i < v.size() && i < max; ++i) out += (i ? sep : std::string()) + v[i];
    return out;
}

std::vector<std::string> split_list(const char* s) {
    std::vector<std::string> out;
    if (!s || !*s) return out;
    for (const std::string& p : split(s, ';')) {
        const std::string t = trim(p);
        if (!t.empty()) out.push_back(t);
    }
    return out;
}

bool looks_like_name(const std::string& t) {
    bool lower = false, digit = false;
    for (unsigned char c : t) {
        if (std::isdigit(c)) digit = true;
        if (std::islower(c)) lower = true;
    }
    if (t.find(' ') != std::string::npos && !digit) return true;
    if (digit) return false;
    if (!t.empty() && std::islower(static_cast<unsigned char>(t[0]))) return true;
    return lower && t.size() > 3 && !fc::find_element(t);
}

double si_from(const std::string& value, const std::string& unit, fc::Dimension d) {
    const double v = num(value, fc::NO_VALUE);
    if (!finite(v)) return fc::NO_VALUE;
    if (d == fc::Dimension::Dimensionless && (unit.empty() || !fc::find_unit(unit, d))) return v;
    return to_si_or(v, unit, d, fc::NO_VALUE);
}

std::string name_tool(const std::string& text) {
    std::ostringstream o;
    const std::string t = trim(text);
    if (t.empty()) return "#Name or formula\nType a formula (Fe2(SO4)3, CuSO4\xC2\xB7" "5H2O, HNO3, C6H12O6) or a name (iron(III) sulfate, nitric acid).\n";
    fc::CompoundName r;
    bool from_name = false;
    try {
        if (looks_like_name(t)) {
            r = fc::formula_from_name(t);
            from_name = r.ok;
        }
        if (!r.ok) r = fc::name_compound(t);
        if (!r.ok && !from_name) {
            fc::CompoundName alt = fc::formula_from_name(t);
            if (alt.ok) {
                r = alt;
                from_name = true;
            }
        }
    } catch (const std::exception& e) {
        return std::string("#Name or formula\n!") + e.what() + "\n";
    }
    if (!r.ok) return "#Name or formula\n!" + (r.explanation.empty() ? std::string("I don't know how to name that") : r.explanation) + "\n";
    o << "#" << (from_name ? "Formula" : "Name") << "\n";
    o << "Name\t" << r.name << "\n";
    o << "Formula\t" << pretty_formula(r.formula) << "\n";
    o << "Kind\t" << fc::to_string(r.kind) << "\n";
    if (!r.other_names.empty()) o << "Also called\t" << join(r.other_names, ", ", 5) << "\n";
    if (r.hydrate_waters > 0) o << "Water of hydration\t" << r.hydrate_waters << " H\xE2\x82\x82O per formula unit\n";
    if (!r.explanation.empty()) o << r.explanation << "\n";
    if (!r.ions.empty()) {
        o << "#Ions\n";
        for (const fc::IonDescriptor& ion : r.ions) o << ion.text(true) << "\t" << ion.name << (ion.polyatomic ? " (polyatomic)" : "") << "\n";
    }
    try {
        const fc::Formula f = fc::Formula::parse(r.formula);
        o << "#Numbers\n";
        o << "Molar mass\t" << format_number(f.molar_mass(), 6) << " g/mol\n";
        const fc::SolubilityResult s = fc::solubility_of(r.formula);
        if (s.solubility != fc::Solubility::Unknown) o << "In water\t" << fc::to_string(s.solubility) << (s.rule.empty() ? "" : " (" + s.rule + ")") << "\n";
        if (fc::is_strong_acid(r.formula)) o << "Acid strength\tstrong acid\n";
        else if (fc::is_weak_acid(r.formula)) o << "Acid strength\tweak acid\n";
        if (fc::is_strong_base(r.formula)) o << "Base strength\tstrong base\n";
        else if (fc::is_weak_base(r.formula)) o << "Base strength\tweak base\n";
        if (fc::is_strong_electrolyte(r.formula)) o << "Electrolyte\tstrong\n";
    } catch (...) {
    }
    return o.str();
}

std::string ions_tool(const std::string& category, const std::string& search) {
    std::ostringstream o;
    std::vector<const fc::PolyatomicIonEntry*> list;
    if (!trim(search).empty()) {
        list = fc::search_polyatomic_ions(trim(search));
    } else {
        for (const fc::PolyatomicIonEntry& e : fc::polyatomic_ions()) list.push_back(&e);
    }
    for (const fc::PolyatomicIonEntry* e : list) {
        if (!category.empty() && category != "all" && category != e->category) continue;
        o << clean(e->name) << "\t" << fc::ion_formula_text(*e, true) << "\t" << e->charge << "\t" << clean(e->category) << "\n";
    }
    return o.str();
}

std::string ion_info(const std::string& text) {
    std::ostringstream o;
    const fc::PolyatomicIonEntry* e = fc::find_polyatomic_ion(trim(text));
    if (!e) {
        std::optional<fc::IonDescriptor> d;
        try {
            d = fc::parse_ion(trim(text));
        } catch (...) {
        }
        if (!d || !d->valid()) return "#Ion\n!no ion called '" + text + "'\n";
        o << "#" << d->name << "\n";
        o << "Formula\t" << d->text(true) << "\n";
        o << "Charge\t" << (d->charge > 0 ? "+" : "") << d->charge << "\n";
        if (!d->classical_name.empty()) o << "Older name\t" << d->classical_name << "\n";
        return o.str();
    }
    o << "#" << e->name << "\n";
    o << "Formula\t" << fc::ion_formula_text(*e, true) << "\n";
    o << "Charge\t" << (e->charge > 0 ? "+" : "") << e->charge << "\n";
    o << "Group\t" << e->category << "\n";
    const std::vector<std::string> others = split_list(e->other_names);
    if (!others.empty()) o << "Also called\t" << join(others, ", ") << "\n";
    if (e->acid && *e->acid) o << "Acid\t" << e->acid << (e->acid_formula && *e->acid_formula ? std::string(" (") + pretty_formula(e->acid_formula) + ")" : std::string()) << "\n";
    try {
        o << "Mass\t" << format_number(fc::ion_formula(*e).molar_mass(), 6) << " g/mol\n";
    } catch (...) {
    }
    const std::string n = e->name;
    if (n.size() > 3 && n.compare(n.size() - 3, 3, "ate") == 0) o << "-ate ions have one more oxygen than the matching -ite ion; their acids end in -ic.\n";
    else if (n.size() > 3 && n.compare(n.size() - 3, 3, "ite") == 0) o << "-ite ions have one fewer oxygen than the -ate ion; their acids end in -ous.\n";
    if (n.rfind("hydrogen ", 0) == 0) o << "Adding H\xE2\x81\xBA to an anion lowers its charge by one and puts \"hydrogen\" in front of the name.\n";
    return o.str();
}

std::string combine_tool(const std::string& cation, const std::string& anion) {
    std::ostringstream o;
    const auto c = ion_from_text(cation, true);
    const auto a = ion_from_text(anion, false);
    if (!c) return "#Combine ions\n!'" + cation + "' isn't a cation I know\n";
    if (!a) return "#Combine ions\n!'" + anion + "' isn't an anion I know\n";
    try {
        const fc::IonicCombination r = fc::combine_ions(*c, *a);
        o << "#" << r.name << "\n";
        o << "Formula\t" << pretty_formula(r.formula) << "\n";
        o << "Ions\t" << r.cation_count << " " << r.cation.text(true) << " + " << r.anion_count << " " << r.anion.text(true) << "\n";
        o << "Charges cancel\t" << r.cation_count << " \xC3\x97 (" << (r.cation.charge > 0 ? "+" : "") << r.cation.charge << ") + " << r.anion_count << " \xC3\x97 (" << r.anion.charge << ") = 0\n";
        o << "Molar mass\t" << format_number(r.composition().molar_mass(), 6) << " g/mol\n";
        const fc::SolubilityResult s = fc::solubility(r.cation, r.anion);
        o << "In water\t" << fc::to_string(s.solubility) << (s.rule.empty() ? "" : " (" + s.rule + ")") << "\n";
    } catch (const std::exception& e) {
        return std::string("#Combine ions\n!") + e.what() + "\n";
    }
    return o.str();
}

void reaction_details(std::ostringstream& o, const fc::Reaction& r, double t, const UnitPrefs& u) {
    try {
        o << "Type\t" << fc::to_string(r.classify()) << "\n";
    } catch (...) {
    }
    try {
        const double dn = r.delta_n_gas();
        if (finite(dn) && std::fabs(dn) > 1e-9) o << "\xCE\x94n (gas)\t" << format_number(dn, 4) << "\n";
    } catch (...) {
    }
    try {
        if (r.is_redox()) {
            o << "#Redox\n";
            for (const fc::OxidationChange& c : r.oxidation_changes()) {
                if (!c.oxidized() && !c.reduced()) continue;
                o << c.element->symbol << "\t" << format_number(c.before.value(), 3) << " \xE2\x86\x92 " << format_number(c.after.value(), 3) << " (" << (c.oxidized() ? "oxidized" : "reduced") << ")\n";
            }
            const double e = r.electrons_transferred();
            if (finite(e)) o << "Electrons transferred\t" << format_number(e, 4) << "\n";
        }
    } catch (...) {
    }
    try {
        if (r.has_thermo_data()) {
            o << "#Thermodynamics at " << u.format(t, "temperature", 5) << "\n";
            const double dh = r.delta_h(t), ds = r.delta_s(t), dg = r.delta_g(t);
            o << "\xCE\x94H\t" << u.format(dh, "molar energy", 5) << (dh < 0 ? " (exothermic)" : " (endothermic)") << "\n";
            o << "\xCE\x94S\t" << u.format(ds, "molar entropy", 5) << "\n";
            o << "\xCE\x94G\t" << u.format(dg, "molar energy", 5) << (dg < 0 ? " (spontaneous)" : " (not spontaneous)") << "\n";
            const double k = r.equilibrium_constant(t);
            if (finite(k)) o << "K\t" << format_number(k, 4) << "\n";
            const double ts = r.spontaneity_temperature();
            if (finite(ts) && ts > 0.0) o << "Switches spontaneity at\t" << u.format(ts, "temperature", 4) << "\n";
        }
    } catch (...) {
    }
    try {
        const fc::IonicEquation ie = fc::ionic_equations(r);
        if (!ie.complete.empty() && ie.complete != ie.net) {
            o << "#Ionic equations\n";
            o << "Complete ionic\t" << ie.complete << "\n";
            o << "Net ionic\t" << (ie.no_net_reaction ? std::string("no net reaction") : ie.net) << "\n";
            if (!ie.spectators.empty()) o << "Spectator ions\t" << join(ie.spectators, ", ") << "\n";
        }
    } catch (...) {
    }
}

std::string predict_tool(ChemShared& shared, const std::string& left, double t) {
    std::ostringstream o;
    if (trim(left).empty()) return "#Predict products\nType the reactants, like CH4 + O2, Zn + HCl, AgNO3 + NaCl, NaOH + H2SO4 or CaCO3.\n";
    fc::PredictedReaction p;
    try {
        p = fc::predict_products(trim(left));
    } catch (const std::exception& e) {
        return std::string("#Predict products\n!") + e.what() + "\n";
    }
    if (!p.ok) return "#Predict products\n!" + (p.explanation.empty() ? std::string("no rule I know fits these reactants") : p.explanation) + "\n";
    o << "#Prediction\n";
    if (p.no_reaction) {
        o << "Result\tno reaction\n";
        shared.last_prediction.clear();
        shared.last_prediction_sides.clear();
    } else {
        o << "Equation\t" << p.equation(true) << "\n";
        shared.last_prediction = p.equation(false);
        std::string lhs, rhs;
        for (const auto& term : p.reaction.reactants()) lhs += (lhs.empty() ? "" : " + ") + term.formula.display_string(false);
        for (const auto& term : p.reaction.products()) rhs += (rhs.empty() ? "" : " + ") + term.formula.display_string(false);
        shared.last_prediction_sides = lhs + "\t" + (p.reaction.reversible() ? "<=>" : "->") + "\t" + rhs;
    }
    if (!p.kind.empty()) o << "Rule used\t" << p.kind << "\n";
    o << "Reaction type\t" << fc::to_string(p.type) << "\n";
    if (!p.explanation.empty()) o << p.explanation << "\n";
    for (const std::string& n : p.notes) o << n << "\n";
    if (!p.no_reaction) {
        if (!p.complete_ionic.empty()) {
            o << "#Ionic equations\n";
            o << "Complete ionic\t" << p.complete_ionic << "\n";
            o << "Net ionic\t" << p.net_ionic << "\n";
            if (!p.spectators.empty()) o << "Spectator ions\t" << join(p.spectators, ", ") << "\n";
        }
        std::ostringstream more;
        reaction_details(more, p.reaction, t, shared.units);
        std::string m = more.str();
        const std::size_t ion = m.find("#Ionic equations");
        if (ion != std::string::npos && !p.complete_ionic.empty()) m = m.substr(0, ion);
        o << m;
    }
    return o.str();
}

std::string balance_tool(ChemShared& shared, const std::string& left, const std::string& arrow, const std::string& right, double t) {
    std::ostringstream o;
    if (trim(left).empty() || trim(right).empty()) return "#Balance\nType both sides, like Fe + O2 and Fe2O3, pick the arrow, and press Balance.\n";
    const bool reversible = arrow == "<=>" || arrow == "<->" || arrow == "\xE2\x87\x8C";
    const std::string text = trim(left) + (reversible ? " <=> " : " -> ") + trim(right);
    fc::Reaction r;
    try {
        r = fc::Reaction::parse(text);
        r.balance();
        r.set_reversible(reversible);
    } catch (const std::exception& e) {
        return std::string("#Balance\n!") + e.what() + "\n";
    }
    shared.last_balanced = r.to_string();
    o << "#Balanced\n";
    o << "Equation\t" << r.to_unicode_string() << "\n";
    try {
        o << "Check\t" << (r.is_balanced() ? "every element and the charge balance" : "not balanced") << "\n";
        if (!r.balance_is_unique()) o << "!More than one set of coefficients works (two independent reactions are mixed together); this is the smallest one.\n";
        std::string coefs;
        for (const auto& term : r.reactants()) coefs += (coefs.empty() ? "" : ", ") + format_number(term.coefficient, 6);
        coefs += " \xE2\x86\x92 ";
        bool first = true;
        for (const auto& term : r.products()) {
            coefs += (first ? "" : ", ") + format_number(term.coefficient, 6);
            first = false;
        }
        o << "Coefficients\t" << coefs << "\n";
        o << "Mass check\t" << format_number(r.mass_balance(), 4) << " g/mol difference\n";
    } catch (...) {
    }
    reaction_details(o, r, t, shared.units);
    return o.str();
}

std::string stoich_species(const std::string& equation) {
    std::ostringstream o;
    try {
        fc::Reaction r = fc::Reaction::parse(trim(equation));
        if (!r.is_balanced()) r.balance();
        std::size_t i = 0;
        for (const auto& term : r.reactants()) o << i++ << "\t" << term.formula.display_string(true) << "\treactant\t" << format_number(term.coefficient, 6) << "\t" << format_number(term.formula.molar_mass(), 6) << "\n";
        for (const auto& term : r.products()) o << i++ << "\t" << term.formula.display_string(true) << "\tproduct\t" << format_number(term.coefficient, 6) << "\t" << format_number(term.formula.molar_mass(), 6) << "\n";
    } catch (const std::exception& e) {
        return std::string("error\t") + e.what();
    }
    return o.str();
}

std::string stoich_tool(ChemShared& shared, const std::vector<std::string>& a) {
    std::ostringstream o;
    const UnitPrefs& u = shared.units;
    fc::Reaction r;
    try {
        r = fc::Reaction::parse(trim(arg(a, 1)));
        if (!r.is_balanced()) r.balance();
    } catch (const std::exception& e) {
        return std::string("#Stoichiometry\n!") + e.what() + "\n";
    }
    const double t = to_si_or(num(arg(a, 3), 298.15), arg(a, 4).empty() ? "K" : arg(a, 4), fc::Dimension::Temperature, 298.15);
    const double p = to_si_or(num(arg(a, 5), 1.0), arg(a, 6).empty() ? "atm" : arg(a, 6), fc::Dimension::Pressure, 101325.0);
    const std::string out_unit = arg(a, 7).empty() ? u.get("mass") : arg(a, 7);
    std::vector<fc::StoichiometryKnown> knowns;
    for (const std::string& entry : split(arg(a, 2), ';')) {
        const auto f = split(entry, '|');
        if (f.size() < 3 || trim(f[1]).empty()) continue;
        fc::StoichiometryKnown k;
        k.species = static_cast<std::size_t>(num(f[0], -1));
        const std::string unit = trim(f[2]);
        const double value = num(f[1], fc::NO_VALUE);
        if (!finite(value)) continue;
        if (fc::find_unit(unit, fc::Dimension::Mass)) {
            k.kind = fc::AmountKind::Mass;
            k.value = to_si_or(value, unit, fc::Dimension::Mass, value);
        } else if (unit == "particles") {
            k.kind = fc::AmountKind::Particles;
            k.value = value;
        } else if (fc::find_unit(unit, fc::Dimension::Amount)) {
            k.kind = fc::AmountKind::Moles;
            k.value = to_si_or(value, unit, fc::Dimension::Amount, value);
        } else if (fc::find_unit(unit, fc::Dimension::Volume)) {
            k.value = to_si_or(value, unit, fc::Dimension::Volume, value);
            if (f.size() >= 5 && finite(num(f[3], fc::NO_VALUE))) {
                k.kind = fc::AmountKind::SolutionVolume;
                k.molarity = to_si_or(num(f[3], 0.0), trim(f[4]).empty() ? "M" : trim(f[4]), fc::Dimension::Concentration, 0.0);
            } else {
                k.kind = fc::AmountKind::GasVolume;
            }
        } else {
            continue;
        }
        knowns.push_back(k);
    }
    if (knowns.empty()) return "#Stoichiometry\nEnter an amount for at least one species (a mass, moles, a gas volume, or a solution volume with its molarity).\n";
    fc::StoichiometryResult res;
    try {
        res = fc::solve_stoichiometry(r, knowns, t, p);
    } catch (const std::exception& e) {
        return std::string("#Stoichiometry\n!") + e.what() + "\n";
    }
    if (!res.ok) return "#Stoichiometry\n!" + res.message + "\n";
    o << "#Result\n";
    o << "Equation\t" << res.reaction.to_unicode_string() << "\n";
    if (res.limiting >= 0 && static_cast<std::size_t>(res.limiting) < res.species.size()) o << "Limiting reactant\t" << pretty_formula(res.species[static_cast<std::size_t>(res.limiting)].formula) << "\n";
    o << "Extent of reaction\t" << format_number(res.extent, 5) << " mol\n";
    if (!res.message.empty()) o << res.message << "\n";
    o << "#Amounts (" << out_unit << ")\n";
    for (const fc::StoichiometrySpecies& s : res.species) {
        const double grams = s.final_moles * s.molar_mass;
        std::string line;
        if (s.reactant) {
            line = "used " + u.format(s.reacted_moles * s.molar_mass / 1000.0, "mass", 5) + ", left " + fc::format_quantity(grams / 1000.0, fc::Dimension::Mass, out_unit, 5) + " (" + format_number(s.final_moles, 5) + " mol)";
        } else {
            line = "made " + fc::format_quantity(grams / 1000.0, fc::Dimension::Mass, out_unit, 5) + " (" + format_number(s.final_moles, 5) + " mol)";
            if (s.gas) line += ", " + u.format(s.final_moles * 8.314462618 * t / p, "volume", 4) + " as gas";
        }
        o << pretty_formula(s.formula) << (s.limiting ? " (limiting)" : "") << "\t" << line << "\n";
    }
    if (!res.steps.empty()) {
        o << "#Steps\n";
        for (const std::string& s : res.steps) o << s << "\n";
    }
    return o.str();
}

std::string equation_detail(const std::string& id) {
    const fc::ChemistryEquation* eq = fc::find_equation(id);
    if (!eq) return "error\tunknown equation";
    std::ostringstream o;
    o << "formula\t" << eq->formula << "\n";
    o << "name\t" << eq->name << "\n";
    if (!eq->notes.empty()) o << "notes\t" << clean(eq->notes) << "\n";
    for (const fc::SolverVariable& v : eq->variables) {
        std::string units;
        if (v.dimension != fc::Dimension::Dimensionless) {
            units = join(fc::unit_symbols(v.dimension), ",");
        } else if (!v.unit_label.empty()) {
            units = v.unit_label;
        }
        o << "var\t" << v.id << "\t" << clean(v.label) << "\t" << fc::to_string(v.dimension) << "\t" << (v.unit.empty() ? v.unit_label : v.unit) << "\t" << units << "\n";
    }
    return o.str();
}

std::string solve_tool(const std::vector<std::string>& a) {
    const fc::ChemistryEquation* eq = fc::find_equation(arg(a, 1));
    if (!eq) return "error\tpick an equation";
    std::vector<double> si(eq->variables.size(), fc::NO_VALUE);
    for (const std::string& entry : split(arg(a, 2), ';')) {
        const std::size_t eqpos = entry.find('=');
        if (eqpos == std::string::npos) continue;
        const std::string id = trim(entry.substr(0, eqpos));
        const auto parts = split(trim(entry.substr(eqpos + 1)), '|');
        const int idx = eq->index_of(id);
        if (idx < 0 || parts.empty() || trim(parts[0]).empty()) continue;
        const fc::SolverVariable& v = eq->variables[static_cast<std::size_t>(idx)];
        si[static_cast<std::size_t>(idx)] = si_from(trim(parts[0]), parts.size() > 1 ? trim(parts[1]) : v.unit, v.dimension);
    }
    int unknown = eq->index_of(arg(a, 3));
    if (unknown < 0) {
        int missing = 0;
        for (std::size_t i = 0; i < si.size(); ++i) {
            if (!finite(si[i])) {
                unknown = static_cast<int>(i);
                ++missing;
            }
        }
        if (missing != 1) return "error\tleave exactly one value blank (or pick what to solve for)";
    }
    si[static_cast<std::size_t>(unknown)] = fc::NO_VALUE;
    for (std::size_t i = 0; i < si.size(); ++i) {
        if (static_cast<int>(i) != unknown && !finite(si[i])) return "error\tfill in " + eq->variables[i].label;
    }
    fc::EquationSolution s;
    try {
        s = fc::solve_equation(*eq, si, unknown);
    } catch (const std::exception& e) {
        return std::string("error\t") + e.what();
    }
    if (!s.ok) return "error\t" + (s.message.empty() ? std::string("no solution found") : s.message);
    const fc::SolverVariable& v = eq->variables[static_cast<std::size_t>(unknown)];
    std::string unit = arg(a, 4).empty() ? v.unit : arg(a, 4);
    std::string shown;
    if (v.dimension == fc::Dimension::Dimensionless) shown = format_number(s.value, 6) + (v.unit_label.empty() ? "" : " " + v.unit_label);
    else shown = fc::format_quantity(s.value, v.dimension, unit, 6);
    std::string raw = v.dimension == fc::Dimension::Dimensionless ? format_number(s.value, 8) : format_number(from_si_or(s.value, unit, v.dimension), 8);
    return "ok\t" + v.id + "\t" + raw + "\t" + clean(v.label) + " = " + shown + (s.message.empty() ? "" : "\t" + clean(s.message));
}

std::string percent_tool(const std::string& formula) {
    std::ostringstream o;
    try {
        const fc::Formula f = fc::Formula::parse(trim(formula));
        o << "#" << pretty_formula(trim(formula)) << "\n";
        o << "Molar mass\t" << format_number(f.molar_mass(), 7) << " g/mol\n";
        for (const fc::PercentComposition& p : fc::percent_composition(f)) {
            o << p.element->name << " (" << p.element->symbol << ")\t" << format_number(p.percent, 5) << "% (" << format_number(p.count, 4) << " \xC3\x97 " << format_number(p.element->atomic_weight, 6) << " = " << format_number(p.mass, 6) << " g/mol)\n";
        }
    } catch (const std::exception& e) {
        return std::string("#Percent composition\n!") + e.what() + "\n";
    }
    return o.str();
}

std::string empirical_tool(const std::string& spec, const std::string& molar_mass) {
    std::ostringstream o;
    std::vector<std::pair<const fc::Element*, double>> parts;
    for (const std::string& entry : split(spec, ';')) {
        const std::size_t eq = entry.find_first_of("=:");
        if (eq == std::string::npos) continue;
        const fc::Element* e = fc::find_element(trim(entry.substr(0, eq)));
        const double v = num(entry.substr(eq + 1), fc::NO_VALUE);
        if (e && finite(v) && v > 0.0) parts.emplace_back(e, v);
    }
    if (parts.empty()) return "#Empirical formula\nEnter element amounts as C=40.0; H=6.7; O=53.3 (grams or percent).\n";
    try {
        const fc::Formula emp = fc::empirical_formula_from_masses(parts);
        o << "#Empirical formula\n";
        o << "Empirical formula\t" << pretty_formula(emp.to_string()) << "\n";
        o << "Empirical mass\t" << format_number(emp.molar_mass(), 6) << " g/mol\n";
        for (const auto& p : parts) o << p.first->symbol << "\t" << format_number(p.second, 5) << " \xC3\xB7 " << format_number(p.first->atomic_weight, 5) << " = " << format_number(p.second / p.first->atomic_weight, 5) << " mol\n";
        const double mm = num(molar_mass, fc::NO_VALUE);
        if (finite(mm) && mm > 0.0) {
            const fc::Formula mol = fc::molecular_formula(emp, mm);
            o << "Molecular formula\t" << pretty_formula(mol.to_string()) << " (\xC3\x97" << format_number(mm / emp.molar_mass(), 3) << ")\n";
            try {
                const fc::CompoundName n = fc::name_of(mol);
                if (n.ok) o << "Could be\t" << n.name << "\n";
            } catch (...) {
            }
        }
    } catch (const std::exception& e) {
        return std::string("#Empirical formula\n!") + e.what() + "\n";
    }
    return o.str();
}

std::string convert_tool(const std::vector<std::string>& a) {
    const double v = num(arg(a, 1), fc::NO_VALUE);
    if (!finite(v)) return "error\tnot a number";
    try {
        const fc::Dimension d = fc::dimension_named(arg(a, 4));
        const double r = arg(a, 4).empty() ? fc::convert(v, arg(a, 2), arg(a, 3)) : fc::convert(v, arg(a, 2), arg(a, 3), d);
        return "ok\t" + format_number(r, 8) + "\t" + format_number(v, 8) + " " + arg(a, 2) + " = " + format_number(r, 8) + " " + arg(a, 3);
    } catch (const std::exception& e) {
        return std::string("error\t") + e.what();
    }
}

} // namespace

std::string chemistry_tool(ChemShared& shared, const std::vector<std::string>& a) {
    const std::string c = arg(a, 0);
    if (c == "units") {
        shared.units.set(arg(a, 1), arg(a, 2));
        return shared.units.get(arg(a, 1));
    }
    if (c == "unit_list") {
        const fc::Dimension d = fc::dimension_named(arg(a, 1));
        return join(fc::unit_symbols(d), "\n");
    }
    if (c == "unit_prefs") {
        std::string out;
        for (const auto& kv : shared.units.units) out += kv.first + "\t" + kv.second + "\n";
        return out;
    }
    if (c == "dimensions") {
        std::string out;
        for (int i = 1; i < static_cast<int>(fc::Dimension::Count); ++i) {
            const fc::Dimension d = static_cast<fc::Dimension>(i);
            if (fc::unit_symbols(d).size() > 1) out += std::string(fc::to_string(d)) + "\n";
        }
        return out;
    }
    if (c == "convert") return convert_tool(a);
    if (c == "name") return name_tool(arg(a, 1));
    if (c == "ions") return ions_tool(arg(a, 1), arg(a, 2));
    if (c == "ion_categories") return join(fc::polyatomic_ion_categories(), "\n");
    if (c == "ion_info") return ion_info(arg(a, 1));
    if (c == "combine") return combine_tool(arg(a, 1), arg(a, 2));
    if (c == "compounds") {
        std::string out;
        const std::string s = trim(arg(a, 1));
        for (const fc::CommonCompoundEntry& e : fc::common_compounds()) {
            if (!s.empty()) {
                std::string hay = std::string(e.name) + " " + e.formula + " " + e.other_names + " " + e.category;
                std::string lo = hay, needle = s;
                std::transform(lo.begin(), lo.end(), lo.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
                std::transform(needle.begin(), needle.end(), needle.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
                if (lo.find(needle) == std::string::npos) continue;
            }
            out += clean(e.name) + "\t" + pretty_formula(e.formula) + "\t" + clean(e.category) + "\t" + (e.smiles && *e.smiles ? "1" : "0") + "\n";
        }
        return out;
    }
    if (c == "predict") return predict_tool(shared, arg(a, 1), to_si_or(num(arg(a, 2), 298.15), arg(a, 3).empty() ? "K" : arg(a, 3), fc::Dimension::Temperature, 298.15));
    if (c == "last_prediction") return shared.last_prediction;
    if (c == "last_prediction_sides") return shared.last_prediction_sides;
    if (c == "last_balanced") return shared.last_balanced;
    if (c == "balance") return balance_tool(shared, arg(a, 1), arg(a, 2), arg(a, 3), to_si_or(num(arg(a, 4), 298.15), arg(a, 5).empty() ? "K" : arg(a, 5), fc::Dimension::Temperature, 298.15));
    if (c == "stoich_species") return stoich_species(arg(a, 1));
    if (c == "stoich") return stoich_tool(shared, a);
    if (c == "equation_categories") return join(fc::equation_categories(), "\n");
    if (c == "equations") {
        std::string out;
        for (const fc::ChemistryEquation& e : fc::chemistry_equations()) {
            if (!arg(a, 1).empty() && arg(a, 1) != "all" && e.category != arg(a, 1)) continue;
            out += e.id + "\t" + clean(e.name) + "\t" + clean(e.category) + "\n";
        }
        return out;
    }
    if (c == "equation") return equation_detail(arg(a, 1));
    if (c == "solve") return solve_tool(a);
    if (c == "percent") return percent_tool(arg(a, 1));
    if (c == "empirical") return empirical_tool(arg(a, 1), arg(a, 2));
    if (c == "molar_mass") {
        try {
            return format_number(fc::Formula::parse(trim(arg(a, 1))).molar_mass(), 7);
        } catch (const std::exception& e) {
            return std::string("error\t") + e.what();
        }
    }
    if (c == "pretty") return pretty_formula(arg(a, 1));
    if (c == "save_dir") return chemistry_save_root().string();
    if (c == "saved") {
        std::string out;
        for (const std::string& kind : saved_chemistry_kinds()) {
            if (!arg(a, 1).empty() && arg(a, 1) != kind) continue;
            for (const std::string& n : list_saved_chemistry(kind)) out += kind + "\t" + n + "\n";
        }
        return out;
    }
    return "error\tunknown tool " + c;
}

} // namespace sapling

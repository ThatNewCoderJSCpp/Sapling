#include "bonding_view.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <sstream>

namespace sapling {

namespace {

constexpr double PI = 3.14159265358979323846;

const fc::Element& el(int z) {
    return fc::element(z);
}

std::string atom_label(const CovalentState& s, int i) {
    if (i < 0 || i >= static_cast<int>(s.atoms.size())) return "?";
    return std::string(el(s.atoms[static_cast<std::size_t>(i)].z).symbol) + std::to_string(i + 1);
}

std::vector<std::string> chem_split(const char* text) {
    std::vector<std::string> out;
    std::string cur;
    for (const char* p = text; p && *p; ++p) {
        if (*p == ';') {
            if (!cur.empty()) out.push_back(cur);
            cur.clear();
        } else if (!(cur.empty() && *p == ' ')) {
            cur += *p;
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

std::string fixed(double v, int digits) {
    if (!std::isfinite(v)) return "n/a";
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*f", digits, v);
    return buf;
}

fc::BondOrder to_order(int o) {
    return static_cast<fc::BondOrder>(std::max(1, std::min(3, o)));
}

int from_order(fc::BondOrder o) {
    return o == fc::BondOrder::Aromatic ? 1 : static_cast<int>(o);
}

double interpolated_energy(const fc::Element& a, const fc::Element& b, double order) {
    const int lo = std::max(1, std::min(3, static_cast<int>(std::floor(order + 1e-9))));
    const int hi = std::min(3, lo + 1);
    double el_lo = fc::NO_VALUE, el_hi = fc::NO_VALUE;
    try {
        el_lo = fc::bond_energy(a, b, to_order(lo));
    } catch (...) {
    }
    if (hi == lo || order - lo < 1e-6) return el_lo;
    try {
        el_hi = fc::bond_energy(a, b, to_order(hi));
    } catch (...) {
    }
    if (!(el_hi == el_hi)) return el_lo;
    if (!(el_lo == el_lo)) return el_hi;
    return el_lo + (el_hi - el_lo) * (order - lo);
}

double table_length(const fc::Element& a, const fc::Element& b, double order) {
    try {
        const int whole = static_cast<int>(std::lround(order));
        if (std::fabs(order - whole) < 1e-6) return fc::bond_length(a, b, to_order(whole));
        const double single = fc::bond_length(a, b, fc::BondOrder::Single);
        return fc::length_from_bond_order(single, order);
    } catch (...) {
        return fc::NO_VALUE;
    }
}

std::vector<Vec3> sphere_directions(int n) {
    std::vector<Vec3> out;
    const double golden = PI * (3.0 - std::sqrt(5.0));
    for (int i = 0; i < n; ++i) {
        const double y = 1.0 - 2.0 * (i + 0.5) / n;
        const double r = std::sqrt(std::max(0.0, 1.0 - y * y));
        const double t = golden * i;
        out.push_back(Vec3{ r * std::cos(t), y, r * std::sin(t) });
    }
    return out;
}

std::string composition_key(const fc::Formula& f) {
    std::vector<std::pair<int, long long>> parts;
    for (const auto& c : f.composition()) parts.emplace_back(c.first->atomic_number, std::llround(c.second * 1000.0));
    std::sort(parts.begin(), parts.end());
    std::string key;
    for (const auto& p : parts) key += std::to_string(p.first) + ":" + std::to_string(p.second) + ";";
    return key;
}

struct KnownName {
    std::string name;
    std::vector<std::string> others;
};

std::vector<KnownName> known_names(const fc::Formula& f, int charge) {
    static std::map<std::string, std::vector<KnownName>> table;
    static bool built = false;
    if (!built) {
        built = true;
        for (const fc::PolyatomicIonEntry& e : fc::polyatomic_ions()) {
            try {
                KnownName k;
                k.name = e.name;
                for (const std::string& o : chem_split(e.other_names)) k.others.push_back(o);
                table[composition_key(fc::ion_formula(e)) + "q" + std::to_string(e.charge)].push_back(k);
            } catch (...) {
            }
        }
        for (const fc::CommonCompoundEntry& e : fc::common_compounds()) {
            try {
                KnownName k;
                k.name = e.name;
                for (const std::string& o : chem_split(e.other_names)) k.others.push_back(o);
                table[composition_key(fc::Formula::parse(e.formula)) + "q0"].push_back(k);
            } catch (...) {
            }
        }
    }
    const auto it = table.find(composition_key(f) + "q" + std::to_string(charge));
    return it == table.end() ? std::vector<KnownName>{} : it->second;
}

Vec3 dipole_vector(const fc::Molecule& m, const std::vector<double>* lengths = nullptr) {
    Vec3 v;
    for (std::size_t bi = 0; bi < m.bonds().size(); ++bi) {
        const fc::MoleculeBond& b = m.bonds()[bi];
        const fc::MoleculeAtom& p = m.atoms()[b.a];
        const fc::MoleculeAtom& q = m.atoms()[b.b];
        const double xp = p.element->electronegativity_pauling == p.element->electronegativity_pauling ? p.element->electronegativity_pauling : 0.0;
        const double xq = q.element->electronegativity_pauling == q.element->electronegativity_pauling ? q.element->electronegativity_pauling : 0.0;
        const Vec3 d{ q.x - p.x, q.y - p.y, q.z - p.z };
        const double l = length(d);
        if (l <= 0.0) continue;
        const double ionic = fc::percent_ionic_character(std::fabs(xp - xq)) / 100.0;
        const double sign = xq > xp ? 1.0 : -1.0;
        const double len = lengths && bi < lengths->size() && (*lengths)[bi] == (*lengths)[bi] ? (*lengths)[bi] : l;
        v = v + d * (sign * fc::bond_dipole(ionic, len) / l);
    }
    return v;
}

} // namespace

void BondingView::ensure_analysis() {
    if (!m_analysis_dirty) return;
    m_analysis_dirty = false;
    CovalentAnalysis an;
    if (m_cov.atoms.empty()) {
        m_analysis = an;
        return;
    }
    try {
        fc::Molecule m;
        for (const SandboxAtom& a : m_cov.atoms) {
            const std::size_t i = m.add_atom(el(a.z));
            m.set_position(i, a.p.x, a.p.y, a.p.z);
        }
        for (const SandboxBond& b : m_cov.bonds) {
            if (b.a == b.b || m.find_bond(static_cast<std::size_t>(b.a), static_cast<std::size_t>(b.b))) continue;
            m.add_bond(static_cast<std::size_t>(b.a), static_cast<std::size_t>(b.b), to_order(b.order));
        }
        an.molecule = m;
        const std::size_t n = m.atom_count();
        const auto pieces = m.fragments();
        if (pieces.size() <= 1) {
            an.lewis = m.lewis_structure(m_cov.charge, true);
            if (m.bond_count() > 0 && n <= 40) {
                an.best = m.best_lewis_structure(m_cov.charge);
                if (an.best.ok && (!an.lewis.ok || an.best.score < an.lewis.score - 1e-6)) {
                    bool differs = !an.lewis.ok;
                    for (std::size_t i = 0; i < an.best.bond_orders.size() && i < m_cov.bonds.size(); ++i) {
                        differs = differs || from_order(an.best.bond_orders[i]) != m_cov.bonds[i].order;
                    }
                    an.better_exists = differs;
                }
            }
        } else {
            std::size_t largest = 0;
            for (std::size_t f = 1; f < pieces.size(); ++f) {
                if (pieces[f].size() > pieces[largest].size()) largest = f;
            }
            an.lewis = fc::LewisStructure{};
            an.lewis.ok = true;
            an.lewis.charge = m_cov.charge;
            an.lewis.atoms.assign(n, fc::LewisAtom{});
            an.best = fc::LewisStructure{};
            an.best.ok = true;
            for (const SandboxBond& b : m_cov.bonds) {
                an.lewis.bond_orders.push_back(to_order(b.order));
                an.best.bond_orders.push_back(to_order(b.order));
            }
            for (std::size_t f = 0; f < pieces.size(); ++f) {
                std::vector<int> local(n, -1);
                fc::Molecule sub;
                for (std::size_t atom : pieces[f]) {
                    local[atom] = static_cast<int>(sub.add_atom(*m.atoms()[atom].element));
                    sub.set_position(static_cast<std::size_t>(local[atom]), m.atoms()[atom].x, m.atoms()[atom].y, m.atoms()[atom].z);
                }
                std::vector<std::size_t> global_bonds;
                for (std::size_t bi = 0; bi < m_cov.bonds.size(); ++bi) {
                    const SandboxBond& b = m_cov.bonds[bi];
                    const int la = local[static_cast<std::size_t>(b.a)], lb = local[static_cast<std::size_t>(b.b)];
                    if (la < 0 || lb < 0 || sub.find_bond(static_cast<std::size_t>(la), static_cast<std::size_t>(lb))) continue;
                    sub.add_bond(static_cast<std::size_t>(la), static_cast<std::size_t>(lb), to_order(b.order));
                    global_bonds.push_back(bi);
                }
                const int q = f == largest ? m_cov.charge : 0;
                const fc::LewisStructure part = sub.lewis_structure(q, true);
                an.lewis.ok = an.lewis.ok && part.ok;
                an.lewis.valence_electrons += part.valence_electrons;
                an.lewis.bonding_electrons += part.bonding_electrons;
                an.lewis.lone_pair_electrons += part.lone_pair_electrons;
                an.lewis.unpaired_electrons += part.unpaired_electrons;
                an.lewis.score += part.score;
                for (const std::string& note : part.notes) an.lewis.notes.push_back(note);
                for (std::size_t atom : pieces[f]) {
                    if (static_cast<std::size_t>(local[atom]) < part.atoms.size()) an.lewis.atoms[atom] = part.atoms[static_cast<std::size_t>(local[atom])];
                }
                if (sub.bond_count() > 0 && sub.atom_count() <= 40) {
                    const fc::LewisStructure best = sub.best_lewis_structure(q);
                    if (best.ok && (!part.ok || best.score < part.score - 1e-6)) {
                        for (std::size_t k = 0; k < global_bonds.size() && k < best.bond_orders.size(); ++k) {
                            an.best.bond_orders[global_bonds[k]] = best.bond_orders[k];
                            if (from_order(best.bond_orders[k]) != m_cov.bonds[global_bonds[k]].order) an.better_exists = true;
                        }
                        if (!part.ok) an.better_exists = true;
                    }
                }
            }
        }
        if (pieces.size() <= 1 && m.bond_count() > 0 && n <= 28) an.resonance = m.resonance_structures(m_cov.charge, 12);
        an.avg_orders.assign(m_cov.bonds.size(), 1.0);
        for (std::size_t i = 0; i < m_cov.bonds.size(); ++i) an.avg_orders[i] = m_cov.bonds[i].order;
        bool user_matches_resonance = false;
        if (an.resonance.has_resonance() && an.resonance.average_bond_orders.size() == m_cov.bonds.size()) {
            for (const fc::LewisStructure& s : an.resonance.structures) {
                bool same = s.bond_orders.size() == m_cov.bonds.size();
                for (std::size_t i = 0; same && i < s.bond_orders.size(); ++i) same = from_order(s.bond_orders[i]) == m_cov.bonds[i].order;
                user_matches_resonance = user_matches_resonance || same;
            }
            if (user_matches_resonance) an.avg_orders = an.resonance.average_bond_orders;
        }
        an.lone_pairs.assign(n, 0);
        an.unpaired.assign(n, 0);
        an.formal.assign(n, 0);
        an.problem.assign(n, false);
        if (an.lewis.atoms.size() == n) {
            for (std::size_t i = 0; i < n; ++i) {
                const fc::LewisAtom& la = an.lewis.atoms[i];
                an.lone_pairs[i] = la.lone_pairs;
                an.unpaired[i] = la.unpaired;
                an.formal[i] = la.formal_charge;
                const bool h = m.atoms()[i].element->atomic_number == 1;
                an.problem[i] = (la.incomplete_octet && !h) || (la.expanded_octet && m.atoms()[i].element->period < 3) || m.is_overbonded(i);
            }
        }
        an.partial = m.partial_charges();
        const Vec3 drawn_dipole = dipole_vector(m);
        an.mu_drawn = length(drawn_dipole);
        if (n >= 2 && m.bond_count() > 0) {
            fc::Molecule relaxed = m;
            relaxed.layout_3d(400);
            std::vector<Vec3> from, to;
            for (std::size_t i = 0; i < n; ++i) {
                from.push_back(Vec3{ relaxed.atoms()[i].x, relaxed.atoms()[i].y, relaxed.atoms()[i].z });
                to.push_back(m_cov.atoms[i].p);
            }
            an.relaxed = from;
            std::vector<double> lengths;
            for (std::size_t i = 0; i < m_cov.bonds.size(); ++i) {
                const SandboxBond& b = m_cov.bonds[i];
                lengths.push_back(table_length(el(m_cov.atoms[static_cast<std::size_t>(b.a)].z), el(m_cov.atoms[static_cast<std::size_t>(b.b)].z), an.avg_orders[i]));
            }
            Vec3 v = dipole_vector(relaxed, &lengths);
            if (length(v) < 0.12) v = Vec3{};
            an.mu_shape = length(v);
            an.dipole_vector = aligned_vector(from, to, v);
        }
        an.polar = an.mu_shape > 0.12;
        an.bonds.resize(m_cov.bonds.size());
        an.intact_fraction = 1.0;
        for (std::size_t i = 0; i < m_cov.bonds.size(); ++i) {
            const SandboxBond& b = m_cov.bonds[i];
            const fc::Element& ea = el(m_cov.atoms[static_cast<std::size_t>(b.a)].z);
            const fc::Element& eb = el(m_cov.atoms[static_cast<std::size_t>(b.b)].z);
            BondReport& r = an.bonds[i];
            r.avg_order = an.avg_orders[i];
            r.table_length = table_length(ea, eb, r.avg_order);
            r.drawn_length = length(m_cov.atoms[static_cast<std::size_t>(b.a)].p - m_cov.atoms[static_cast<std::size_t>(b.b)].p);
            r.energy = interpolated_energy(ea, eb, r.avg_order);
            try {
                r.thermo = fc::bond_dissociation(ea, eb, to_order(static_cast<int>(std::lround(r.avg_order))), m_temperature, m_pressure, r.table_length, r.energy);
            } catch (...) {
            }
            const double xa = ea.electronegativity_pauling == ea.electronegativity_pauling ? ea.electronegativity_pauling : 0.0;
            const double xb = eb.electronegativity_pauling == eb.electronegativity_pauling ? eb.electronegativity_pauling : 0.0;
            r.en_difference = std::fabs(xa - xb);
            r.ionic_percent = fc::percent_ionic_character(r.en_difference);
            r.dipole = bond_dipole_debye(ea.atomic_number, eb.atomic_number, r.table_length == r.table_length ? r.table_length : r.drawn_length);
            if (r.thermo.fraction_broken == r.thermo.fraction_broken) an.intact_fraction *= (1.0 - r.thermo.fraction_broken);
        }

        const auto& frags = pieces;
        int radicals = 0;
        for (int u : an.unpaired) radicals += u;
        bool octet_issue = false;
        bool impossible = !an.lewis.ok;
        for (std::size_t i = 0; i < n; ++i) {
            if (m.is_overbonded(i) && an.lewis.atoms.size() == n && an.lewis.atoms[i].expanded_octet && m.atoms()[i].element->period < 3) impossible = true;
            if (an.problem[i]) octet_issue = true;
        }
        for (const std::string& note : an.lewis.notes) an.issues.push_back(note);
        if (frags.size() > 1) an.issues.push_back(std::to_string(frags.size()) + " separate pieces that aren't bonded to each other");
        int fc_sum = 0;
        for (int f : an.formal) fc_sum += std::abs(f);
        if (an.better_exists) an.issues.push_back("a structure with fewer or smaller formal charges exists (Fix Lewis applies it)");
        if (an.intact_fraction < 0.5) an.issues.push_back("at " + format_t() + " most of these molecules fall apart (" + fixed(100.0 * (1.0 - an.intact_fraction), 1) + "% have a broken bond)");
        if (impossible) an.verdict = "Impossible as drawn";
        else if (radicals > 0) an.verdict = "Radical: odd electron count, very reactive";
        else if (an.intact_fraction < 0.5) an.verdict = "Falls apart at this temperature";
        else if (octet_issue) an.verdict = "Unusual Lewis structure (octet rule broken)";
        else if (an.better_exists) an.verdict = "Plausible, but a better Lewis structure exists";
        else if (frags.size() > 1) an.verdict = "Separate pieces";
        else if (fc_sum > 2) an.verdict = "Stable but charge-separated";
        else an.verdict = "Stable Lewis structure";
        an.ok = true;
    } catch (const std::exception& e) {
        an.ok = false;
        an.error = e.what();
        an.verdict = std::string("can't analyze this: ") + e.what();
    }
    m_analysis = std::move(an);
}

bool BondingView::load_smiles(const std::string& smiles, int charge, const std::string& name) {
    fc::Molecule m;
    try {
        m = fc::Molecule::from_smiles(smiles);
    } catch (...) {
        return false;
    }
    if (m.atom_count() == 0) return false;
    try {
        if (m_3d) m.layout_3d(500);
        else m.layout_2d(140.0);
    } catch (...) {
    }
    CovalentState s;
    s.name = name;
    s.charge = charge != 0 ? charge : m.charge();
    for (const fc::MoleculeAtom& a : m.atoms()) s.atoms.push_back(SandboxAtom{ a.element->atomic_number, Vec3{ a.x, a.y, m_3d ? a.z : 0.0 } });
    bool aromatic = false;
    for (const fc::MoleculeBond& b : m.bonds()) {
        aromatic = aromatic || b.order == fc::BondOrder::Aromatic;
        s.bonds.push_back(SandboxBond{ static_cast<int>(b.a), static_cast<int>(b.b), from_order(b.order) });
    }
    m_cov = std::move(s);
    m_mode = BondingMode::Covalent;
    m_sel_atom = m_sel_bond = -1;
    m_analysis_dirty = true;
    if (aromatic) apply_best_lewis();
    if (!m_3d) m_cam.face_front();
    frame_view();
    changed(true);
    return true;
}

void BondingView::tidy() {
    if (m_cov.atoms.empty()) return;
    m_analysis_dirty = true;
    ensure_analysis();
    if (!m_analysis.ok) return;
    fc::Molecule m = m_analysis.molecule;
    try {
        if (m_3d) m.layout_3d(600);
        else m.layout_2d(140.0);
    } catch (...) {
        return;
    }
    for (std::size_t i = 0; i < m_cov.atoms.size() && i < m.atom_count(); ++i) {
        m_cov.atoms[i].p = Vec3{ m.atoms()[i].x, m.atoms()[i].y, m_3d ? m.atoms()[i].z : 0.0 };
    }
    m_analysis_dirty = true;
}

void BondingView::fill_hydrogens() {
    m_analysis_dirty = true;
    ensure_analysis();
    if (!m_analysis.ok) return;
    const fc::Molecule& m = m_analysis.molecule;
    const std::size_t heavy = m_cov.atoms.size();
    const std::vector<Vec3> dirs3 = sphere_directions(96);
    for (std::size_t i = 0; i < heavy; ++i) {
        if (m_cov.atoms[i].z == 1) continue;
        int missing = m.missing_bonds(i);
        if (missing <= 0) continue;
        std::vector<Vec3> used;
        for (const SandboxBond& b : m_cov.bonds) {
            if (b.a != static_cast<int>(i) && b.b != static_cast<int>(i)) continue;
            const int o = b.a == static_cast<int>(i) ? b.b : b.a;
            used.push_back(normalized(m_cov.atoms[static_cast<std::size_t>(o)].p - m_cov.atoms[i].p));
        }
        const int lone = i < m_analysis.lone_pairs.size() ? std::max(0, std::min(3, m_analysis.lone_pairs[i] - missing)) : 0;
        for (int k = 0; k < missing; ++k) {
            Vec3 best{ 1.0, 0.0, 0.0 };
            double best_score = -2.0;
            if (m_3d) {
                for (const Vec3& d : dirs3) {
                    double worst = 2.0;
                    for (const Vec3& u : used) worst = std::min(worst, 1.0 - dot(d, u));
                    if (used.empty()) worst = d.y;
                    if (worst > best_score) {
                        best_score = worst;
                        best = d;
                    }
                }
            } else {
                std::vector<double> angles;
                for (const Vec3& u : used) angles.push_back(std::atan2(-u.y, u.x));
                const std::vector<double> pick = lone_pair_angles(angles, 1 + (used.empty() ? 0 : 0));
                const double a = pick.empty() ? 0.0 : pick.front();
                best = Vec3{ std::cos(a), -std::sin(a), 0.0 };
            }
            (void)lone;
            const Vec3 p = m_cov.atoms[i].p + best * ideal_length(m_cov.atoms[i].z, 1, 1);
            const int h = add_atom(1, p);
            m_cov.bonds.push_back(SandboxBond{ static_cast<int>(i), h, 1 });
            used.push_back(best);
        }
    }
    m_analysis_dirty = true;
}

void BondingView::apply_best_lewis() {
    m_analysis_dirty = true;
    ensure_analysis();
    if (!m_analysis.ok || m_cov.bonds.empty()) return;
    fc::LewisStructure best = m_analysis.best;
    if (!best.ok) {
        try {
            best = m_analysis.molecule.best_lewis_structure(m_cov.charge);
        } catch (...) {
            return;
        }
    }
    if (!best.ok) return;
    for (std::size_t i = 0; i < m_cov.bonds.size() && i < best.bond_orders.size(); ++i) m_cov.bonds[i].order = from_order(best.bond_orders[i]);
    m_analysis_dirty = true;
}

void BondingView::apply_resonance(int index) {
    ensure_analysis();
    if (index < 0 || index >= static_cast<int>(m_analysis.resonance.structures.size())) return;
    const fc::LewisStructure& s = m_analysis.resonance.structures[static_cast<std::size_t>(index)];
    for (std::size_t i = 0; i < m_cov.bonds.size() && i < s.bond_orders.size(); ++i) m_cov.bonds[i].order = from_order(s.bond_orders[i]);
    m_analysis_dirty = true;
}

std::string BondingView::covalent_readouts() {
    ensure_analysis();
    std::ostringstream o;
    const UnitPrefs& u = m_shared.units;
    if (m_cov.atoms.empty()) {
        o << "#Covalent sandbox\n";
        o << "Pick an element from the table, then click the canvas to place atoms. Drag from one atom to another to bond them; click a bond with the Bond tool to make it double or triple.\n";
        o << "Fill H adds the hydrogens each atom needs, Tidy straightens the drawing (VSEPR shape in 3D), and Fix Lewis picks the bond orders with the smallest formal charges.\n";
        return o.str();
    }
    if (!m_analysis.ok) {
        o << "#Problem\n!" << m_analysis.verdict << "\n";
        return o.str();
    }
    const fc::Molecule& m = m_analysis.molecule;
    const fc::Formula f = m.formula();
    std::string formula = f.hill();
    o << "#Molecule\n";
    o << "Formula\t" << pretty_formula(formula) << superscript_charge(m_cov.charge) << "\n";
    const std::vector<KnownName> known = known_names(f, m_cov.charge);
    std::string systematic;
    std::string kind;
    try {
        fc::Formula named = fc::Formula::parse(formula);
        named.set_charge(m_cov.charge);
        const fc::CompoundName cn = fc::name_of(named);
        if (cn.ok) {
            systematic = cn.name;
            kind = fc::to_string(cn.kind);
        }
    } catch (...) {
    }
    if (!known.empty()) {
        o << "Name\t" << known.front().name << (m_cov.charge != 0 ? " ion" : "") << "\n";
        std::vector<std::string> others;
        for (std::size_t i = 1; i < known.size(); ++i) others.push_back(known[i].name + " (isomer)");
        for (const std::string& x : known.front().others) others.push_back(x);
        if (!systematic.empty() && systematic != known.front().name && m_cov.charge == 0) others.push_back(systematic);
        std::string list;
        for (std::size_t i = 0; i < others.size() && i < 4; ++i) list += (i ? ", " : "") + others[i];
        if (!list.empty()) o << "Also called\t" << list << "\n";
        o << "Kind\t" << (m_cov.charge != 0 ? std::string("polyatomic ion") : kind.empty() ? std::string("molecule") : kind) << "\n";
    } else if (!systematic.empty()) {
        o << "Name\t" << systematic << "\n";
        if (!kind.empty()) o << "Kind\t" << kind << "\n";
    }
    if (!m_cov.name.empty()) o << "Label\t" << m_cov.name << "\n";
    o << "Molar mass\t" << u.format(m.molar_mass(), "molar mass", 6) << "\n";
    o << "Atoms, bonds\t" << m.atom_count() << ", " << m.bond_count() << "\n";
    o << "Valence electrons\t" << m_analysis.lewis.valence_electrons << " (" << m_analysis.lewis.bonding_electrons << " bonding, " << m_analysis.lewis.lone_pair_electrons << " in lone pairs)\n";

    o << "#Stability\n";
    o << "Verdict\t" << m_analysis.verdict << "\n";
    for (const std::string& issue : m_analysis.issues) o << "!" << issue << "\n";

    o << "#Lewis structure\n";
    int lp_total = 0, fc_total = 0;
    for (int x : m_analysis.lone_pairs) lp_total += x;
    for (int x : m_analysis.formal) fc_total += std::abs(x);
    o << "Lone pairs\t" << lp_total << "\n";
    o << "Formal charges\t" << (fc_total == 0 ? std::string("all zero") : "total |q| = " + std::to_string(fc_total)) << "\n";
    for (std::size_t i = 0; i < m.atom_count() && i < m_analysis.lewis.atoms.size(); ++i) {
        const fc::LewisAtom& la = m_analysis.lewis.atoms[i];
        if (m_cov.atoms[i].z == 1 && la.formal_charge == 0) continue;
        std::string what = std::to_string(la.lone_pairs) + " LP, " + std::to_string(la.electrons) + " e\xE2\x81\xBB";
        if (la.formal_charge != 0) what += ", formal " + std::string(la.formal_charge > 0 ? "+" : "") + std::to_string(la.formal_charge);
        if (la.unpaired) what += ", " + std::to_string(la.unpaired) + " unpaired";
        if (la.expanded_octet) what += ", expanded octet";
        if (la.incomplete_octet && m_cov.atoms[i].z != 1) what += ", incomplete octet";
        o << atom_label(m_cov, static_cast<int>(i)) << "\t" << what << "\n";
    }
    if (m_analysis.better_exists) o << "!Fix Lewis gives a better set of bond orders\n";

    o << "#Resonance\n";
    if (m_analysis.resonance.has_resonance()) {
        o << "Structures\t" << m_analysis.resonance.count() << "\n";
        std::map<std::string, std::pair<double, int>> groups;
        for (std::size_t i = 0; i < m_cov.bonds.size() && i < m_analysis.resonance.average_bond_orders.size(); ++i) {
            const double avg = m_analysis.resonance.average_bond_orders[i];
            if (std::fabs(avg - std::round(avg)) < 0.02) continue;
            const SandboxBond& b = m_cov.bonds[i];
            std::string a = el(m_cov.atoms[static_cast<std::size_t>(b.a)].z).symbol, c = el(m_cov.atoms[static_cast<std::size_t>(b.b)].z).symbol;
            if (a > c) std::swap(a, c);
            std::string key = a + "\xE2\x80\x93" + c + " " + fixed(avg, 2);
            groups[key].first = avg;
            groups[key].second += 1;
        }
        for (const auto& g : groups) o << "Average order\t" << g.first << " (\xC3\x97" << g.second.second << ")\n";
        o << "The real molecule is a blend of these structures, so those bonds have non-integer orders and equal lengths.\n";
    } else {
        o << "Structures\t1 (no resonance)\n";
    }

    o << "#Shape and polarity\n";
    for (std::size_t i = 0; i < m.atom_count(); ++i) {
        if (m.degree(i) < 2) continue;
        const int lp = i < m_analysis.lone_pairs.size() ? m_analysis.lone_pairs[i] : 0;
        const fc::VseprResult v = fc::vsepr(m.degree(i), lp);
        std::string text = std::string(fc::to_string(v.molecular_geometry)) + " (" + v.axe_notation + ", " + fc::to_string(v.hybridization) + ")";
        if (!v.ideal_angles.empty()) text += ", " + fixed(v.ideal_angles.front(), 1) + "\xC2\xB0";
        o << atom_label(m_cov, static_cast<int>(i)) << "\t" << text << "\n";
    }
    o << "Polarity\t" << (m_analysis.polar ? "polar" : "nonpolar") << "\n";
    o << "Dipole (VSEPR shape)\t" << fixed(m_analysis.mu_shape, 2) << " D\n";
    if (m_3d) o << "Dipole (as drawn)\t" << fixed(m_analysis.mu_drawn, 2) << " D\n";

    o << "#Energy at " << format_t() << "\n";
    double total = 0.0;
    for (const BondReport& r : m_analysis.bonds) total += r.energy == r.energy ? r.energy : 0.0;
    o << "Bond energy total\t" << u.format_kj(total, 5) << "\n";
    try {
        const double hf = m.estimated_enthalpy_of_formation();
        if (hf == hf) o << "\xCE\x94H\xC2\xB0" "f (formation, est.)\t" << u.format_kj(hf, 4) << "\n";
    } catch (...) {
    }
    double dh_t = 0.0, dg_t = 0.0;
    bool have_t = false;
    int weakest = -1;
    double weakest_g = 1e300;
    for (std::size_t i = 0; i < m_analysis.bonds.size(); ++i) {
        const BondReport& r = m_analysis.bonds[i];
        if (r.thermo.enthalpy == r.thermo.enthalpy) {
            dh_t += r.thermo.enthalpy;
            dg_t += r.thermo.gibbs;
            have_t = true;
            if (r.thermo.gibbs < weakest_g) {
                weakest_g = r.thermo.gibbs;
                weakest = static_cast<int>(i);
            }
        }
    }
    if (have_t) {
        o << "\xCE\x94H to atomize\t" << u.format_kj(dh_t, 5) << "\n";
        o << "\xCE\x94G to atomize\t" << u.format_kj(dg_t, 5) << "\n";
        o << "Intact molecules\t" << fixed(100.0 * m_analysis.intact_fraction, m_analysis.intact_fraction > 0.999 ? 3 : 1) << "%\n";
        if (weakest >= 0) {
            const SandboxBond& b = m_cov.bonds[static_cast<std::size_t>(weakest)];
            o << "Weakest bond\t" << atom_label(m_cov, b.a) << "\xE2\x80\x93" << atom_label(m_cov, b.b) << " (\xCE\x94G " << u.format_kj(weakest_g, 4) << ")\n";
        }
    }
    return o.str();
}

std::string BondingView::selection_readouts() {
    ensure_analysis();
    if (m_sel_bond >= 0 && m_sel_bond < static_cast<int>(m_cov.bonds.size())) return bond_readout(m_sel_bond);
    if (m_sel_atom >= 0 && m_sel_atom < static_cast<int>(m_cov.atoms.size())) return atom_readout(m_sel_atom);
    return "#Selection\nClick a bond or an atom to see its details.\n";
}

std::string BondingView::bond_readout(int index) {
    std::ostringstream o;
    const UnitPrefs& u = m_shared.units;
    const SandboxBond& b = m_cov.bonds[static_cast<std::size_t>(index)];
    const fc::Element& ea = el(m_cov.atoms[static_cast<std::size_t>(b.a)].z);
    const fc::Element& eb = el(m_cov.atoms[static_cast<std::size_t>(b.b)].z);
    o << "#Bond " << atom_label(m_cov, b.a) << "\xE2\x80\x93" << atom_label(m_cov, b.b) << "\n";
    if (!m_analysis.ok || index >= static_cast<int>(m_analysis.bonds.size())) {
        o << "!" << m_analysis.verdict << "\n";
        return o.str();
    }
    const BondReport& r = m_analysis.bonds[static_cast<std::size_t>(index)];
    static const char* names[] = { "", "single", "double", "triple" };
    o << "Drawn as\t" << names[std::max(1, std::min(3, b.order))] << "\n";
    const bool fractional = std::fabs(r.avg_order - std::round(r.avg_order)) > 0.02;
    o << "Bond order\t" << (fractional ? fixed(r.avg_order, 3) + " (non-integer, from resonance)" : std::to_string(static_cast<int>(std::lround(r.avg_order))) + " (integer)") << "\n";
    o << "Character\t" << fc::to_string(r.en_difference < 0.4 ? fc::BondCharacter::NonpolarCovalent : r.en_difference < 1.8 ? fc::BondCharacter::PolarCovalent : fc::BondCharacter::Ionic) << "\n";
    o << "Average length\t" << u.format(r.table_length * 1e-12, "length", 4) << "\n";
    o << "Length as drawn\t" << u.format(r.drawn_length * 1e-12, "length", 4) << "\n";
    o << "Average bond energy\t" << u.format_kj(r.energy, 4) << "\n";
    o << "#At " << format_t() << "\n";
    if (r.thermo.enthalpy == r.thermo.enthalpy) {
        o << "\xCE\x94H to break\t" << u.format_kj(r.thermo.enthalpy, 4) << "\n";
        o << "\xCE\x94S to break\t" << u.format(r.thermo.entropy, "molar entropy", 4) << "\n";
        o << "\xCE\x94G to break\t" << u.format_kj(r.thermo.gibbs, 4) << "\n";
        if (r.thermo.equilibrium_constant == r.thermo.equilibrium_constant) o << "K (dissociation)\t" << format_number(r.thermo.equilibrium_constant, 4) << "\n";
        o << "Fraction broken\t" << format_number(100.0 * r.thermo.fraction_broken, 4) << "%\n";
        if (r.thermo.half_broken_temperature == r.thermo.half_broken_temperature) o << "Half broken at\t" << u.format(r.thermo.half_broken_temperature, "temperature", 4) << "\n";
        o << "Raising T makes \xCE\x94G = \xCE\x94H \xE2\x88\x92 T\xCE\x94S smaller, and lowering P favors the two separate atoms, so both break bonds.\n";
    } else {
        o << "No bond energy data for this pair.\n";
    }
    o << "#Polarity\n";
    o << "\xCE\x94\xCF\x87 (Pauling)\t" << fixed(r.en_difference, 2) << "\n";
    o << "Ionic character\t" << fixed(r.ionic_percent, 1) << "%\n";
    o << "Bond dipole\t" << fixed(r.dipole, 2) << " D\n";
    const double xa = ea.electronegativity_pauling, xb = eb.electronegativity_pauling;
    if (r.en_difference > 0.05) o << "\xCE\xB4\xE2\x88\x92 end\t" << (xa > xb ? ea.symbol : eb.symbol) << "\n";
    try {
        const double k = fc::gordy_force_constant(ea, eb, r.table_length, r.avg_order);
        const double mu = fc::reduced_mass(ea.atomic_weight, eb.atomic_weight);
        const double freq = fc::vibrational_frequency(k, mu);
        if (k == k && freq == freq && freq > 0.0) {
            o << "#Vibration\n";
            o << "Force constant\t" << format_number(k, 4) << " N/m\n";
            o << "Stretch\t" << format_number(freq / 2.99792458e10, 4) << " cm\xE2\x81\xBB\xC2\xB9\n";
        }
    } catch (...) {
    }
    return o.str();
}

std::string BondingView::atom_readout(int index) {
    std::ostringstream o;
    const SandboxAtom& a = m_cov.atoms[static_cast<std::size_t>(index)];
    const fc::Element& e = el(a.z);
    o << "#Atom " << atom_label(m_cov, index) << " \xC2\xB7 " << e.name << "\n";
    o << "Atomic number\t" << static_cast<int>(e.atomic_number) << "\n";
    o << "Electronegativity\t" << (e.electronegativity_pauling == e.electronegativity_pauling ? fixed(e.electronegativity_pauling, 2) : std::string("n/a")) << "\n";
    o << "Valence electrons\t" << e.valence_electrons() << "\n";
    if (m_analysis.ok && static_cast<std::size_t>(index) < m_analysis.lone_pairs.size()) {
        const fc::Molecule& m = m_analysis.molecule;
        o << "Lone pairs\t" << m_analysis.lone_pairs[static_cast<std::size_t>(index)] << "\n";
        o << "Formal charge\t" << m_analysis.formal[static_cast<std::size_t>(index)] << "\n";
        if (static_cast<std::size_t>(index) < m_analysis.partial.size()) o << "Partial charge\t" << fixed(m_analysis.partial[static_cast<std::size_t>(index)], 3) << " e\n";
        o << "Bonds (order sum)\t" << m.degree(static_cast<std::size_t>(index)) << " (" << fixed(m.bond_order_sum(static_cast<std::size_t>(index)), 1) << ")\n";
        o << "Usual valence\t" << m.standard_valence(static_cast<std::size_t>(index)) << "\n";
        if (m.degree(static_cast<std::size_t>(index)) >= 1) {
            const fc::VseprResult v = fc::vsepr(m.degree(static_cast<std::size_t>(index)), m_analysis.lone_pairs[static_cast<std::size_t>(index)]);
            o << "Electron geometry\t" << fc::to_string(v.electron_geometry) << "\n";
            o << "Shape\t" << fc::to_string(v.molecular_geometry) << " (" << v.axe_notation << ")\n";
            o << "Hybridization\t" << fc::to_string(v.hybridization) << "\n";
            const std::vector<std::size_t> nb = m.neighbors(static_cast<std::size_t>(index));
            for (std::size_t i = 0; i < nb.size() && i < 6; ++i) {
                for (std::size_t j = i + 1; j < nb.size() && j < 6; ++j) {
                    const Vec3 p = m_cov.atoms[static_cast<std::size_t>(index)].p;
                    const Vec3 u = normalized(m_cov.atoms[nb[i]].p - p), w = normalized(m_cov.atoms[nb[j]].p - p);
                    const double ang = std::acos(std::max(-1.0, std::min(1.0, dot(u, w)))) * 180.0 / PI;
                    o << "Angle " << atom_label(m_cov, static_cast<int>(nb[i])) << "\xE2\x80\x93" << atom_label(m_cov, index) << "\xE2\x80\x93" << atom_label(m_cov, static_cast<int>(nb[j])) << "\t" << fixed(ang, 1) << "\xC2\xB0 drawn";
                    if (!v.ideal_angles.empty()) o << ", " << fixed(v.ideal_angles.front(), 1) << "\xC2\xB0 ideal";
                    o << "\n";
                }
            }
        }
    }
    return o.str();
}

} // namespace sapling

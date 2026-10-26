#include "molecule_draw.hpp"
#include "view_util.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sapling {

using namespace chem_util;

namespace {

constexpr double PI = 3.14159265358979323846;

class StructureView : public ChemView {
public:
    explicit StructureView(ChemShared& shared) : m_shared(shared) {
        m_cam.perspective = 0.0008;
        show("carbonate");
    }

    void draw(Canvas& canvas, const Theme& theme) override {
        m_cam.width = canvas.width();
        m_cam.height = canvas.height();
        canvas.clear(theme.background);
        if (m_atoms.empty()) {
            draw_hud(canvas, theme, { "Structure", m_error.empty() ? std::string("pick an ion or a compound") : m_error }, false);
            return;
        }
        if (!m_framed) {
            std::vector<Vec3> pts;
            for (const DrawAtom& a : m_atoms) pts.push_back(a.p);
            m_cam.frame(pts, 70.0);
            m_cam.zoom = std::min(m_cam.zoom, std::min(canvas.width(), canvas.height()) / 380.0);
            m_framed = true;
        }
        DrawOptions opt;
        opt.flat = !m_3d;
        opt.lone_pairs = true;
        opt.formal_charges = true;
        opt.labels = true;
        draw_molecule(canvas, m_cam, theme, m_atoms, m_bonds, opt);
        std::vector<std::string> hud{ m_title };
        if (m_resonance > 1) hud.push_back(std::to_string(m_resonance) + " resonance structures (one shown)");
        draw_hud(canvas, theme, hud, false);
    }

    bool mouse_down(Vec2 p, int, unsigned) override {
        m_down = true;
        m_last = p;
        return false;
    }

    bool mouse_move(Vec2 p, bool) override {
        if (!m_down) return false;
        const Vec2 d = p - m_last;
        m_last = p;
        if (m_3d) m_cam.orbit(d.x, d.y);
        else m_cam.pan = m_cam.pan + d;
        return true;
    }

    bool mouse_up(Vec2, int) override {
        m_down = false;
        return false;
    }

    bool wheel(Vec2, double dy) override {
        m_cam.zoom_by(dy > 0 ? 1.12 : 1.0 / 1.12);
        return true;
    }

    bool double_click(Vec2) override {
        m_framed = false;
        return true;
    }

    bool step(double dt) override {
        if (m_spin && m_3d && !m_down && !m_atoms.empty()) {
            m_cam.yaw += dt * 0.6;
            return true;
        }
        return false;
    }

    std::string command(const std::vector<std::string>& a) override {
        const std::string c = arg(a, 0);
        if (c == "show") return show(arg(a, 1)) ? "ok\t" + m_title : "error\t" + m_error;
        if (c == "view") {
            m_3d = arg(a, 1) != "2d";
            relayout();
            return m_3d ? "3d" : "2d";
        }
        if (c == "spin") {
            m_spin = arg(a, 1) == "1";
            return arg(a, 1);
        }
        if (c == "smiles") return m_smiles + "\t" + std::to_string(m_charge) + "\t" + m_name;
        if (c == "readouts") return m_readouts;
        if (c == "version") return std::to_string(m_version);
        return "error\tunknown command " + c;
    }

private:
    bool show(const std::string& text) {
        m_error.clear();
        int q = 0;
        std::string name, smiles;
        std::optional<fc::Molecule> m = molecule_from_text(text, &q, &name, &smiles);
        if (!m) {
            m_error = "no structure for '" + text + "'";
            m_atoms.clear();
            m_bonds.clear();
            m_readouts = "#Structure\n!" + m_error + "\n";
            ++m_version;
            return false;
        }
        m_mol = *m;
        m_charge = q;
        m_name = name;
        m_smiles = smiles;
        try {
            const fc::LewisStructure best = m_mol.best_lewis_structure(q);
            if (best.ok) m_mol.apply_lewis(best);
            m_lewis = m_mol.lewis_structure(q, true);
            m_resonance = static_cast<int>(m_mol.atom_count() <= 28 ? m_mol.resonance_structures(q, 12).count() : 1);
        } catch (...) {
            m_resonance = 1;
        }
        const fc::Formula f = m_mol.formula();
        m_title = (name.empty() ? f.hill() : name) + "  \xC2\xB7  " + pretty_formula(f.hill()) + (f.charge() == 0 ? superscript_charge(q) : std::string());
        relayout();
        build_readouts();
        return true;
    }

    void relayout() {
        try {
            if (m_3d) m_mol.layout_3d(500);
            else m_mol.layout_2d(140.0);
        } catch (...) {
        }
        m_atoms.clear();
        m_bonds.clear();
        for (std::size_t i = 0; i < m_mol.atom_count(); ++i) {
            const fc::MoleculeAtom& at = m_mol.atoms()[i];
            DrawAtom d;
            d.p = Vec3{ at.x, at.y, m_3d ? at.z : 0.0 };
            d.z = at.element->atomic_number;
            if (m_lewis.atoms.size() == m_mol.atom_count()) {
                d.lone_pairs = m_lewis.atoms[i].lone_pairs;
                d.unpaired = m_lewis.atoms[i].unpaired;
                d.formal_charge = m_lewis.atoms[i].formal_charge;
            }
            m_atoms.push_back(d);
        }
        for (const fc::MoleculeBond& b : m_mol.bonds()) {
            DrawBond d;
            d.a = static_cast<int>(b.a);
            d.b = static_cast<int>(b.b);
            d.order = fc::numeric_bond_order(b.order);
            m_bonds.push_back(d);
        }
        if (!m_3d) m_cam.face_front();
        else {
            m_cam.yaw = 0.5;
            m_cam.pitch = 0.3;
        }
        m_framed = false;
        ++m_version;
    }

    void build_readouts() {
        std::ostringstream o;
        const fc::Formula f = m_mol.formula();
        o << "#Structure\n";
        o << "Formula\t" << pretty_formula(f.hill()) << (f.charge() == 0 ? superscript_charge(m_charge) : std::string()) << "\n";
        if (!m_name.empty()) o << "Name\t" << m_name << "\n";
        o << "Charge\t" << (m_charge > 0 ? "+" : "") << m_charge << "\n";
        o << "Molar mass\t" << m_shared.units.format(m_mol.molar_mass(), "molar mass", 6) << "\n";
        o << "Valence electrons\t" << m_lewis.valence_electrons << "\n";
        o << "Resonance structures\t" << m_resonance << "\n";
        int center = -1, best_degree = 0;
        for (std::size_t i = 0; i < m_mol.atom_count(); ++i) {
            if (m_mol.degree(i) > best_degree) {
                best_degree = m_mol.degree(i);
                center = static_cast<int>(i);
            }
        }
        if (center >= 0 && best_degree >= 2) {
            const int lp = m_lewis.atoms.size() > static_cast<std::size_t>(center) ? m_lewis.atoms[static_cast<std::size_t>(center)].lone_pairs : 0;
            const fc::VseprResult v = fc::vsepr(best_degree, lp);
            o << "Shape around " << m_mol.atoms()[static_cast<std::size_t>(center)].element->symbol << "\t" << fc::to_string(v.molecular_geometry) << " (" << v.axe_notation << ", " << fc::to_string(v.hybridization) << ")\n";
        }
        for (const std::string& n : m_lewis.notes) o << "!" << n << "\n";
        if (!m_smiles.empty()) o << "SMILES\t" << m_smiles << "\n";
        m_readouts = o.str();
    }

    ChemShared& m_shared;
    Camera3D m_cam;
    bool m_3d = true;
    bool m_spin = true;
    bool m_framed = false;
    bool m_down = false;
    Vec2 m_last;
    fc::Molecule m_mol;
    fc::LewisStructure m_lewis;
    std::vector<DrawAtom> m_atoms;
    std::vector<DrawBond> m_bonds;
    int m_charge = 0;
    int m_resonance = 1;
    std::string m_name;
    std::string m_smiles;
    std::string m_title;
    std::string m_error;
    std::string m_readouts;
    std::uint64_t m_version = 1;
};

class ElementView : public ChemView {
public:
    explicit ElementView(ChemShared& shared) : m_shared(shared) {}

    void draw(Canvas& canvas, const Theme& theme) override {
        canvas.clear(theme.background);
        const fc::Element& e = fc::element(m_z);
        const double w = canvas.width(), h = canvas.height();
        const bool both = m_style == "both";
        if (m_style == "bohr" || both) draw_bohr(canvas, theme, e, both ? Vec2{ w * 0.27, h * 0.5 } : Vec2{ w * 0.5, h * 0.5 }, both ? std::min(w * 0.25, h * 0.45) : std::min(w, h) * 0.45);
        if (m_style == "boxes" || both) draw_boxes(canvas, theme, e, both ? w * 0.54 : 30.0, 40.0, both ? w * 0.44 : w - 60.0);
        draw_hud(canvas, theme, { std::string(e.name) + " (" + e.symbol + ")  \xC2\xB7  Z = " + std::to_string(e.atomic_number), config_text() }, false);
    }

    bool step(double dt) override {
        if (m_style == "boxes") return false;
        m_time += dt;
        return true;
    }

    std::string command(const std::vector<std::string>& a) override {
        const std::string c = arg(a, 0);
        if (c == "element") {
            if (const fc::Element* e = element_of(arg(a, 1))) m_z = e->atomic_number;
            m_charge = static_cast<int>(num(arg(a, 2), 0.0));
            return fc::element(m_z).symbol;
        }
        if (c == "charge") {
            m_charge = static_cast<int>(num(arg(a, 1), 0.0));
            return std::to_string(m_charge);
        }
        if (c == "style") {
            m_style = arg(a, 1).empty() ? "both" : arg(a, 1);
            return m_style;
        }
        if (c == "readouts") return readouts();
        return "error\tunknown command " + c;
    }

private:
    fc::ElectronConfiguration config() const {
        const fc::Element& e = fc::element(m_z);
        try {
            if (m_charge != 0 && m_charge < m_z) return fc::ElectronConfiguration::of_ion(e, m_charge);
        } catch (...) {
        }
        return fc::ElectronConfiguration::ground_state(e);
    }

    std::string config_text() const {
        try {
            std::string s = config().to_unicode_string(true);
            if (m_charge != 0) s = fc::element(m_z).symbol + superscript_charge(m_charge) + ": " + s;
            return s;
        } catch (...) {
            return fc::element(m_z).electron_configuration;
        }
    }

    void draw_bohr(Canvas& canvas, const Theme& theme, const fc::Element& e, Vec2 c, double radius) {
        const fc::ElectronConfiguration cfg = config();
        std::vector<int> shells = cfg.shells();
        while (!shells.empty() && shells.back() == 0) shells.pop_back();
        const int n = static_cast<int>(shells.size());
        const double r0 = radius * 0.2;
        canvas.sphere_px(c, r0, rgb(230, 120, 110), 1.0);
        canvas.text(Vec2{ c.x, c.y - 1.0 }, std::to_string(e.protons()) + "p\xE2\x81\xBA", std::max(9.0, r0 * 0.45), rgb(30, 30, 40), TextAlign::Center);
        canvas.text(Vec2{ c.x, c.y + r0 * 0.5 }, std::to_string(e.neutrons()) + "n\xE2\x81\xB0", std::max(9.0, r0 * 0.45), rgb(30, 30, 40), TextAlign::Center);
        for (int k = 0; k < n; ++k) {
            const double r = r0 + (radius - r0) * (k + 1) / std::max(1, n);
            canvas.circle_px(c, r, theme.grid, 1.3);
            canvas.text(Vec2{ c.x + r * 0.71 + 4.0, c.y - r * 0.71 - 4.0 }, std::string(1, fc::shell_letter(k + 1)) + " " + std::to_string(shells[static_cast<std::size_t>(k)]), 11.0, theme.muted, TextAlign::Left);
            const int count = shells[static_cast<std::size_t>(k)];
            const double speed = 1.2 / (k + 1);
            for (int i = 0; i < count; ++i) {
                const double ang = m_time * speed + 2.0 * PI * i / std::max(1, count) + k * 0.4;
                const Vec2 p{ c.x + r * std::cos(ang), c.y + r * std::sin(ang) };
                canvas.disc_px(p, std::max(2.5, radius * 0.02), k == n - 1 ? theme.yellow : theme.accent);
            }
        }
    }

    void draw_boxes(Canvas& canvas, const Theme& theme, const fc::Element&, double x0, double y0, double width) {
        const fc::ElectronConfiguration cfg = config();
        std::vector<fc::Subshell> subs = cfg.subshells();
        std::stable_sort(subs.begin(), subs.end(), [](const fc::Subshell& a, const fc::Subshell& b) {
            if (a.madelung() != b.madelung()) return a.madelung() < b.madelung();
            return a.n < b.n;
        });
        const double box = 22.0;
        double x = x0, y = y0 + 30.0;
        canvas.text(Vec2{ x0, y0 + 10.0 }, "Orbital diagram (lowest energy first)", 12.0, theme.muted);
        for (const fc::Subshell& s : subs) {
            const int orbitals = s.orbitals();
            const double need = 40.0 + orbitals * (box + 2.0) + 14.0;
            if (x + need > x0 + width && x > x0) {
                x = x0;
                y += box + 30.0;
            }
            canvas.text(Vec2{ x, y + box * 0.7 }, fc::subshell_name(s.n, s.l), 13.0, theme.text);
            const double bx = x + 36.0;
            for (int o = 0; o < orbitals; ++o) {
                const double ox = bx + o * (box + 2.0);
                canvas.stroke_rect_px(ox, y, box, box, theme.axes, 1.2);
                const int up = s.electrons > o ? 1 : 0;
                const int down = s.electrons > orbitals + o ? 1 : 0;
                if (up) canvas.arrow_px(Vec2{ ox + box * 0.36, y + box - 3.0 }, Vec2{ ox + box * 0.36, y + 3.0 }, theme.accent, 1.6, 5.0);
                if (down) canvas.arrow_px(Vec2{ ox + box * 0.66, y + 3.0 }, Vec2{ ox + box * 0.66, y + box - 3.0 }, theme.warm, 1.6, 5.0);
            }
            x += need;
        }
        y += box + 40.0;
        const int unpaired = cfg.unpaired_electrons();
        canvas.text(Vec2{ x0, y }, std::to_string(unpaired) + " unpaired electron" + (unpaired == 1 ? "" : "s") + (unpaired ? " (paramagnetic)" : " (diamagnetic)"), 12.0, theme.muted);
    }

    std::string readouts() const {
        std::ostringstream o;
        const fc::Element& e = fc::element(m_z);
        const UnitPrefs& u = m_shared.units;
        o << "#" << e.name << "\n";
        o << "Symbol\t" << e.symbol << "\n";
        o << "Atomic number\t" << static_cast<int>(e.atomic_number) << "\n";
        o << "Molar mass\t" << u.format(e.atomic_weight, "molar mass", 7) << "\n";
        o << "Category\t" << e.category_name() << " (" << e.classification_name() << ")\n";
        o << "Group, period, block\t" << static_cast<int>(e.group) << ", " << static_cast<int>(e.period) << ", " << e.block_name() << "\n";
        o << "Protons, neutrons\t" << e.protons() << ", " << e.neutrons() << " (mass number " << e.mass_number() << ")\n";
        o << "#Electrons\n";
        o << "Configuration\t" << config_text() << "\n";
        try {
            o << "Shells\t" << config().shell_string() << "\n";
            o << "Unpaired\t" << config().unpaired_electrons() << "\n";
        } catch (...) {
        }
        o << "Valence electrons\t" << e.valence_electrons() << "\n";
        try {
            const fc::ElectronConfiguration cfg = config();
            const fc::QuantumNumbers q = cfg.last_electron();
            o << "Last electron\tn=" << q.n << ", l=" << q.l << ", m\xE2\x82\x97=" << q.ml << ", m\xE2\x82\x9B=" << (q.ms > 0 ? "+\xC2\xBD" : "\xE2\x88\x92\xC2\xBD") << "\n";
            o << "Zeff (outer)\t" << fixed(cfg.effective_nuclear_charge(m_z), 2) << "\n";
        } catch (...) {
        }
        o << "#Energies\n";
        if (finite(e.electronegativity_pauling)) o << "Electronegativity\t" << fixed(e.electronegativity_pauling, 2) << "\n";
        for (int i = 1; i <= 3; ++i) {
            const double ie = e.ionization_energy_kj_per_mol(i);
            if (finite(ie)) o << "Ionization energy " << i << "\t" << u.format_kj(ie, 5) << "\n";
        }
        if (e.has_electron_affinity()) o << "Electron affinity\t" << u.format_kj(e.electron_affinity_kj_per_mol(), 4) << "\n";
        o << "#Size\n";
        if (finite(e.atomic_radius)) o << "Atomic radius\t" << u.format(e.atomic_radius * 1e-12, "length", 4) << "\n";
        if (finite(e.covalent_radius)) o << "Covalent radius\t" << u.format(e.covalent_radius * 1e-12, "length", 4) << "\n";
        if (finite(e.van_der_waals_radius)) o << "van der Waals radius\t" << u.format(e.van_der_waals_radius * 1e-12, "length", 4) << "\n";
        if (finite(e.metallic_radius)) o << "Metallic radius\t" << u.format(e.metallic_radius * 1e-12, "length", 4) << "\n";
        o << "#Bulk\n";
        o << "Standard state\t" << fc::to_string(e.standard_state) << "\n";
        if (finite(e.melting_point)) o << "Melting point\t" << u.format(e.melting_point, "temperature", 5) << "\n";
        if (finite(e.boiling_point)) o << "Boiling point\t" << u.format(e.boiling_point, "temperature", 5) << "\n";
        if (finite(e.density)) o << "Density\t" << u.format(e.density * 1000.0, "density", 4) << "\n";
        if (e.crystal_structure && *e.crystal_structure) o << "Crystal structure\t" << e.crystal_structure << "\n";
        std::string ox;
        for (const auto& s : e.oxidation_states()) {
            if (!s.common) continue;
            ox += (ox.empty() ? "" : ", ") + std::string(s.state > 0 ? "+" : "") + std::to_string(s.state);
        }
        if (!ox.empty()) o << "Common oxidation states\t" << ox << "\n";
        if (e.discovery_year > 0) o << "Discovered\t" << e.discovery_year << "\n";
        if (!e.is_stable) o << "!radioactive (no stable isotopes)\n";
        return o.str();
    }

    ChemShared& m_shared;
    int m_z = 8;
    int m_charge = 0;
    std::string m_style = "both";
    double m_time = 0.0;
};

} // namespace

std::unique_ptr<ChemView> make_structure_view(ChemShared& shared) {
    return std::unique_ptr<ChemView>(new StructureView(shared));
}

std::unique_ptr<ChemView> make_element_view(ChemShared& shared) {
    return std::unique_ptr<ChemView>(new ElementView(shared));
}

} // namespace sapling

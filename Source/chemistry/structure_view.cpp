#include "molecule_draw.hpp"
#include "view_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

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
        const double top = 62.0;
        if (m_style == "nucleus") {
            draw_nucleus(canvas, theme, Vec2{ w * 0.5, top + (h - top - 50.0) * 0.5 }, std::max(20.0, std::min(w, h - top - 50.0) * 0.42));
        } else {
            if (m_style == "bohr" || both) {
                const double area_w = both ? w * 0.46 : w;
                const double radius = std::max(20.0, std::min(area_w * 0.46, (h - top - 34.0) * 0.48));
                const Vec2 center{ both ? area_w * 0.5 : w * 0.5, top + (h - top - 34.0) * 0.5 };
                draw_bohr(canvas, theme, e, center, radius);
                draw_counts(canvas, theme, Vec2{ center.x, std::min(h - 10.0, center.y + radius + 24.0) });
            }
            if (m_style == "boxes" || both) draw_boxes(canvas, theme, both ? w * 0.48 : 30.0, top, both ? w * 0.5 : w - 60.0);
        }
        std::string head = std::string(e.name) + " (" + e.symbol + ")  \xC2\xB7  Z = " + std::to_string(e.atomic_number);
        if (m_extra != 0 || m_charge != 0) head += "  \xC2\xB7  " + isotope_label();
        draw_hud(canvas, theme, { head, config_text(true) }, false);
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
        m_yaw += d.x * 0.01;
        m_pitch = std::max(-1.4, std::min(1.4, m_pitch + d.y * 0.01));
        return m_nucleus == "3d" && m_style != "boxes";
    }

    bool mouse_up(Vec2, int) override {
        m_down = false;
        return false;
    }

    bool step(double dt) override {
        if (m_style == "boxes") return false;
        m_time += dt;
        if (m_nucleus == "3d" && !m_down) m_yaw += dt * 0.4;
        return true;
    }

    std::string command(const std::vector<std::string>& a) override {
        const std::string c = arg(a, 0);
        if (c == "element") {
            if (const fc::Element* e = element_of(arg(a, 1))) m_z = e->atomic_number;
            m_charge = static_cast<int>(num(arg(a, 2), 0.0));
            if (a.size() > 3) m_extra = static_cast<int>(num(arg(a, 3), 0.0));
            return fc::element(m_z).symbol;
        }
        if (c == "charge") {
            m_charge = static_cast<int>(num(arg(a, 1), 0.0));
            return std::to_string(m_charge);
        }
        if (c == "neutrons") {
            m_extra = std::max(-base_neutrons(), std::min(200, static_cast<int>(num(arg(a, 1), 0.0))));
            return std::to_string(m_extra);
        }
        if (c == "style") {
            m_style = arg(a, 1).empty() ? "both" : arg(a, 1);
            return m_style;
        }
        if (c == "nucleus") {
            const std::string v = arg(a, 1);
            m_nucleus = v == "2d" || v == "3d" ? v : "circle";
            return m_nucleus;
        }
        if (c == "temperature") {
            m_t = std::max(1.0, num(arg(a, 1), 298.15));
            return "ok";
        }
        if (c == "readouts") return readouts();
        return "error\tunknown command " + c;
    }

private:
    int base_neutrons() const {
        const fc::Element& e = fc::element(m_z);
        const int a = e.most_abundant_mass_number ? e.most_abundant_mass_number : e.reference_mass_number;
        return std::max(0, a - m_z);
    }

    int neutrons() const { return std::max(0, base_neutrons() + m_extra); }
    int mass_number() const { return m_z + neutrons(); }

    std::string isotope_label() const {
        return superscript_digits(std::to_string(mass_number())) + fc::element(m_z).symbol + superscript_charge(m_charge);
    }

    fc::ElectronConfiguration config() const {
        const fc::Element& e = fc::element(m_z);
        try {
            if (m_charge != 0 && m_charge < m_z) return fc::ElectronConfiguration::of_ion(e, m_charge);
        } catch (...) {
        }
        return fc::ElectronConfiguration::ground_state(e);
    }

    std::string config_text(bool noble) const {
        try {
            std::string s = config().to_unicode_string(noble);
            if (m_charge != 0) s = fc::element(m_z).symbol + superscript_charge(m_charge) + ": " + s;
            return s;
        } catch (...) {
            return fc::element(m_z).electron_configuration;
        }
    }

    void outlined(Canvas& canvas, Vec2 p, const std::string& s, double size) {
        const Color dark = rgb(20, 20, 30);
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                if (dx || dy) canvas.text(Vec2{ p.x + dx * 1.2, p.y + dy * 1.2 }, s, size, dark, TextAlign::Center, 1.2);
            }
        }
        canvas.text(p, s, size, rgb(250, 250, 255), TextAlign::Center, 1.2);
    }

    static std::string grouped(double v) {
        if (!(v > 0) || !std::isfinite(v)) return "?";
        const double mag = std::pow(10.0, std::floor(std::log10(v)) - 1.0);
        long long x = static_cast<long long>(std::llround(v / mag) * mag);
        std::string d = std::to_string(x), out;
        for (std::size_t i = 0; i < d.size(); ++i) {
            if (i > 0 && (d.size() - i) % 3 == 0) out += ',';
            out += d[i];
        }
        return out;
    }

    void nucleon_sites(bool flat, std::vector<Vec3>& sites, std::vector<bool>& is_proton) const {
        const int z = m_z, n = neutrons(), total = z + n;
        sites.clear();
        is_proton.clear();
        if (total <= 0) return;
        const int span = static_cast<int>(std::ceil(flat ? std::sqrt(static_cast<double>(total)) : std::cbrt(static_cast<double>(total)))) + 2;
        if (flat) {
            for (int j = -span; j <= span; ++j) {
                for (int i = -span; i <= span; ++i) {
                    sites.push_back(Vec3{ i + 0.5 * (j & 1) + 0.0013 * j, j * 0.8660254, 0.0 });
                }
            }
        } else {
            for (int i = -span; i <= span; ++i) {
                for (int j = -span; j <= span; ++j) {
                    for (int k = -span; k <= span; ++k) {
                        if (((i + j + k) & 1) != 0) continue;
                        sites.push_back(Vec3{ (i + 0.013 * j) * 0.70710678, (j + 0.007 * k) * 0.70710678, (k + 0.011 * i) * 0.70710678 });
                    }
                }
            }
        }
        std::sort(sites.begin(), sites.end(), [](const Vec3& p, const Vec3& q) { return p.x * p.x + p.y * p.y + p.z * p.z < q.x * q.x + q.y * q.y + q.z * q.z; });
        if (static_cast<int>(sites.size()) > total) sites.resize(static_cast<std::size_t>(total));
        std::vector<int> order(sites.size());
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
        std::uint32_t seed = 2463534242u + static_cast<std::uint32_t>(z * 131 + n);
        for (std::size_t i = order.size(); i > 1; --i) {
            seed ^= seed << 13;
            seed ^= seed >> 17;
            seed ^= seed << 5;
            std::swap(order[i - 1], order[seed % i]);
        }
        is_proton.assign(sites.size(), false);
        for (int i = 0; i < z && i < static_cast<int>(order.size()); ++i) is_proton[static_cast<std::size_t>(order[static_cast<std::size_t>(i)])] = true;
    }

    void draw_nucleons(Canvas& canvas, Vec2 c, double radius) const {
        const Color proton = rgb(235, 100, 95), neutron = rgb(170, 175, 190);
        const bool flat = m_nucleus == "2d";
        std::vector<Vec3> sites;
        std::vector<bool> is_proton;
        nucleon_sites(flat, sites, is_proton);
        if (sites.empty()) return;
        double reach = 0.0;
        for (const Vec3& q : sites) reach = std::max(reach, std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z));
        const double unit = radius / (reach + 0.5);
        if (flat) {
            for (std::size_t i = 0; i < sites.size(); ++i) canvas.sphere_px(Vec2{ c.x + sites[i].x * unit, c.y + sites[i].y * unit }, unit * 0.49, is_proton[i] ? proton : neutron, 1.0);
            return;
        }
        Camera3D cam;
        cam.width = canvas.width();
        cam.height = canvas.height();
        cam.yaw = m_yaw;
        cam.pitch = m_pitch;
        cam.zoom = 1.0;
        cam.pan = Vec2{ c.x - 0.5 * canvas.width(), c.y - 0.5 * canvas.height() };
        Scene3D scene(cam);
        for (std::size_t i = 0; i < sites.size(); ++i) scene.sphere(Vec3{ sites[i].x * unit, sites[i].y * unit, sites[i].z * unit }, unit * 0.5, is_proton[i] ? proton : neutron);
        scene.draw(canvas);
    }

    void draw_counts(Canvas& canvas, const Theme& theme, Vec2 at) const {
        const int electrons = std::max(0, m_z - m_charge);
        struct Part {
            std::string text;
            Color color;
        };
        const Part parts[] = { { std::to_string(m_z) + " protons", rgb(235, 100, 95) }, { std::to_string(neutrons()) + " neutrons", rgb(170, 175, 190) }, { std::to_string(electrons) + " electrons", theme.accent } };
        const double size = 12.0;
        double total = 0.0;
        for (const Part& p : parts) total += 14.0 + canvas.text_width(p.text, size) + 18.0;
        double x = std::max(8.0, std::min(canvas.width() - 8.0 - (total - 18.0), at.x - 0.5 * (total - 18.0)));
        for (const Part& p : parts) {
            canvas.disc_px(Vec2{ x + 5.0, at.y - 4.0 }, 5.0, p.color);
            canvas.text(Vec2{ x + 14.0, at.y }, p.text, size, theme.text, TextAlign::Left);
            x += 14.0 + canvas.text_width(p.text, size) + 18.0;
        }
    }

    void draw_nucleus(Canvas& canvas, const Theme& theme, Vec2 c, double radius) {
        const int z = m_z, n = neutrons(), total = z + n;
        if (total <= 0) return;
        if (m_nucleus == "circle") {
            canvas.sphere_px(c, radius * 0.8, rgb(200, 90, 85), 1.0);
            const double fs = std::max(14.0, std::min(30.0, radius * 0.2));
            outlined(canvas, Vec2{ c.x, c.y - 4.0 }, std::to_string(z) + "p\xE2\x81\xBA", fs);
            outlined(canvas, Vec2{ c.x, c.y + fs + 4.0 }, std::to_string(n) + "n\xE2\x81\xB0", fs);
        } else {
            draw_nucleons(canvas, c, radius * 0.9);
        }
        draw_counts(canvas, theme, Vec2{ c.x, c.y + radius + 24.0 });
        const double r_fm = 1.2 * std::cbrt(static_cast<double>(total));
        std::string line = "nuclear radius \xE2\x89\x88 1.2 fm \xC3\x97 A^(1/3) = " + fixed(r_fm, 2) + " fm" + (fc::element(m_z).atomic_radius > 0 ? ", about 1/" + grouped(fc::element(m_z).atomic_radius * 1000.0 / r_fm) + " of the atom's radius" : std::string());
        if (m_nucleus == "3d") line += "  \xC2\xB7  drag to turn";
        canvas.text(Vec2{ c.x, c.y + radius + 42.0 }, line, 11.0, theme.muted, TextAlign::Center);
    }

    void draw_bohr(Canvas& canvas, const Theme& theme, const fc::Element&, Vec2 c, double radius) {
        const fc::ElectronConfiguration cfg = config();
        std::vector<int> shells = cfg.shells();
        while (!shells.empty() && shells.back() == 0) shells.pop_back();
        const int n = static_cast<int>(shells.size());
        const double r0 = radius * (n > 4 ? 0.17 : 0.22) * (m_nucleus == "circle" ? 1.0 : 1.2);
        if (m_nucleus == "circle") {
            canvas.sphere_px(c, r0, rgb(200, 90, 85), 1.0);
            const double fs = std::max(10.0, std::min(18.0, r0 * 0.42));
            outlined(canvas, Vec2{ c.x, c.y - 2.0 }, std::to_string(m_z) + "p\xE2\x81\xBA", fs);
            outlined(canvas, Vec2{ c.x, c.y + fs + 2.0 }, std::to_string(neutrons()) + "n\xE2\x81\xB0", fs);
        } else {
            draw_nucleons(canvas, c, r0);
        }
        for (int k = 0; k < n; ++k) {
            const double r = r0 + 10.0 + (radius - r0 - 10.0) * (k + 1) / std::max(1, n);
            canvas.circle_px(c, r, theme.grid, 1.3);
            canvas.text(Vec2{ c.x + r * 0.71 + 4.0, c.y - r * 0.71 - 4.0 }, std::string(1, fc::shell_letter(k + 1)) + " " + std::to_string(shells[static_cast<std::size_t>(k)]), 11.0, theme.muted, TextAlign::Left);
            const int count = shells[static_cast<std::size_t>(k)];
            const double speed = 1.2 / (k + 1);
            for (int i = 0; i < count; ++i) {
                const double ang = m_time * speed + 2.0 * 3.14159265358979 * i / std::max(1, count) + k * 0.4;
                const Vec2 p{ c.x + r * std::cos(ang), c.y + r * std::sin(ang) };
                canvas.disc_px(p, std::max(2.5, radius * 0.018), k == n - 1 ? theme.yellow : theme.accent);
            }
        }
    }

    void draw_boxes(Canvas& canvas, const Theme& theme, double x0, double y0, double width) {
        const fc::ElectronConfiguration cfg = config();
        std::vector<fc::Subshell> subs = cfg.subshells();
        std::stable_sort(subs.begin(), subs.end(), [](const fc::Subshell& a, const fc::Subshell& b) {
            if (a.madelung() != b.madelung()) return a.madelung() < b.madelung();
            return a.n < b.n;
        });
        const double box = 20.0;
        double x = x0, y = y0 + 22.0;
        canvas.text(Vec2{ x0, y0 + 6.0 }, "Orbital diagram (lowest energy first)", 12.0, theme.muted);
        for (const fc::Subshell& s : subs) {
            const int orbitals = s.orbitals();
            const double need = 36.0 + orbitals * (box + 2.0) + 12.0;
            if (x + need > x0 + width && x > x0) {
                x = x0;
                y += box + 24.0;
            }
            canvas.text(Vec2{ x, y + box * 0.7 }, fc::subshell_name(s.n, s.l), 13.0, theme.text);
            const double bx = x + 32.0;
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
        y += box + 34.0;
        canvas.text(Vec2{ x0, y }, "Full:", 12.0, theme.muted);
        y += 18.0;
        std::string full = config_text(false);
        while (!full.empty()) {
            std::size_t cut = full.size();
            while (cut > 0 && canvas.text_width(full.substr(0, cut), 13.0) > width) {
                const std::size_t sp = full.rfind(' ', cut - 1);
                if (sp == std::string::npos || sp == 0) break;
                cut = sp;
            }
            canvas.text(Vec2{ x0, y }, full.substr(0, cut), 13.0, theme.text);
            y += 18.0;
            full = cut < full.size() ? full.substr(cut + 1) : std::string();
        }
        canvas.text(Vec2{ x0, y + 4.0 }, "Noble-gas form:", 12.0, theme.muted);
        canvas.text(Vec2{ x0, y + 24.0 }, config_text(true), 13.0, theme.text);
        const int unpaired = cfg.unpaired_electrons();
        canvas.text(Vec2{ x0, y + 48.0 }, std::to_string(unpaired) + " unpaired electron" + (unpaired == 1 ? "" : "s") + (unpaired ? " (paramagnetic atom)" : " (diamagnetic atom)"), 12.0, theme.muted);
    }

    std::string readouts() const {
        std::ostringstream o;
        const fc::Element& e = fc::element(m_z);
        const UnitPrefs& u = m_shared.units;
        o << "#" << e.name << "\n";
        o << "Symbol\t" << e.symbol << "\n";
        o << "Atomic number\t" << static_cast<int>(e.atomic_number) << "\n";
        o << "Molar mass (natural)\t" << u.format(e.atomic_weight, "molar mass", 7) << "\n";
        o << "Category\t" << e.category_name() << " (" << e.classification_name() << ")\n";
        o << "Group, period, block\t" << static_cast<int>(e.group) << ", " << static_cast<int>(e.period) << ", " << e.block_name() << "\n";
        o << "#Isotope " << isotope_label() << "\n";
        o << "Protons, neutrons\t" << m_z << ", " << neutrons() << (m_extra ? " (" + std::string(m_extra > 0 ? "+" : "") + std::to_string(m_extra) + " from the usual)" : std::string()) << "\n";
        o << "Mass number A\t" << mass_number() << "\n";
        const fc::Isotope* iso = e.isotope(mass_number());
        if (iso) {
            if (iso->mass == iso->mass) o << "Isotope mass\t" << format_number(iso->mass, 9) << " u\n";
            if (iso->is_natural()) o << "Natural abundance\t" << format_number(iso->abundance * (iso->abundance <= 1.0 ? 100.0 : 1.0), 5) << "%\n";
            if (iso->is_stable) o << "Stability\tstable\n";
            else if (iso->half_life_known()) {
                const double hl = iso->half_life;
                std::string shown = hl < 1e-3 ? format_number(hl * 1e6, 4) + " \xC2\xB5s" : hl < 60 ? format_number(hl, 4) + " s" : hl < 3600 ? format_number(hl / 60.0, 4) + " min" : hl < 86400 ? format_number(hl / 3600.0, 4) + " h" : hl < 3.15576e7 * 2 ? format_number(hl / 86400.0, 4) + " days" : format_number(hl / 3.15576e7, 4) + " years";
                o << "Half-life\t" << shown << "\n";
            }
            std::string decays;
            for (const auto& d : iso->decays()) {
                if (d.sub_branch) continue;
                std::string dn = d.name();
                if (dn == "beta-" || dn == "B-") dn = "\xCE\xB2\xE2\x81\xBB";
                else if (dn == "beta+" || dn == "B+") dn = "\xCE\xB2\xE2\x81\xBA";
                else if (dn == "alpha" || dn == "A") dn = "\xCE\xB1";
                decays += (decays.empty() ? "" : ", ") + dn + (d.percent_known() ? " " + format_number(d.percent, 3) + "%" : std::string());
            }
            if (!decays.empty()) o << "Decays by\t" << decays << "\n";
            if (iso->has_spin()) o << "Nuclear spin\t" << format_number(iso->spin(), 2) << (iso->is_nmr_active() ? " (NMR active)" : "") << "\n";
            const double be = iso->binding_energy_per_nucleon();
            if (finite(be)) o << "Binding energy\t" << format_number(be, 5) << " MeV per nucleon\n";
            const double mev = 931.49410242;
            o << "Mass defect\t" << format_number(iso->binding_energy() / mev, 5) << " u\n";
        } else {
            o << "!no known nucleus has " << m_z << " protons and " << neutrons() << " neutrons\n";
            const double be = fc::semi_empirical_binding_energy(m_z, mass_number());
            if (finite(be) && be > 0) o << "Binding energy (formula)\t" << format_number(be / mass_number(), 4) << " MeV per nucleon (semi-empirical)\n";
            else o << "The strong force can't hold this many " << (neutrons() > m_z ? "neutrons" : "protons") << " together; it would fall apart at once.\n";
        }
        o << "Neutron-to-proton ratio\t" << fixed(static_cast<double>(neutrons()) / std::max(1, m_z), 3) << "\n";
        o << "#Electrons\n";
        o << "Configuration\t" << config_text(false) << "\n";
        o << "Noble-gas form\t" << config_text(true) << "\n";
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
        try {
            const fc::TransportProperties t = fc::element_transport(e, m_t, 60.0);
            o << "#Heat at " << u.format(m_t, "temperature", 4) << "\n";
            if (finite(t.molar_heat_capacity)) o << "Molar heat capacity\t" << u.format(t.molar_heat_capacity, "molar entropy", 4) << "\n";
            if (finite(t.specific_heat)) o << "Specific heat\t" << format_number(t.specific_heat / 1000.0, 4) << " J/(g\xC2\xB7K)\n";
            if (finite(t.debye_temperature)) o << "Debye temperature\t" << u.format(t.debye_temperature, "temperature", 3) << " (estimate)\n";
            if (finite(t.thermal_conductivity)) o << "Thermal conductivity\t" << format_number(t.thermal_conductivity, 4) << " W/(m\xC2\xB7K)\n";
            if (finite(t.thermal_diffusivity)) o << "Thermal diffusivity\t" << format_number(t.thermal_diffusivity * 1e6, 4) << " mm\xC2\xB2/s\n";
            o << "#Electricity and magnetism\n";
            o << "Conductor type\t" << fc::to_string(t.kind) << "\n";
            if (finite(t.band_gap)) o << "Band gap\t" << format_number(t.band_gap, 3) << " eV\n";
            if (finite(t.resistivity)) o << "Resistivity\t" << format_number(t.resistivity, 4) << " \xCE\xA9\xC2\xB7m\n";
            if (finite(t.conductivity)) o << "Conductivity\t" << format_number(t.conductivity, 4) << " S/m\n";
            if (finite(t.temperature_coefficient) && t.kind == fc::ConductorKind::Metal) o << "Temperature coefficient\t" << format_number(t.temperature_coefficient, 3) << " per K\n";
            if (finite(t.mean_free_path)) o << "Electron mean free path\t" << format_number(t.mean_free_path * 1e9, 3) + " nm" << "\n";
            if (finite(t.skin_depth)) o << "Skin depth at 60 Hz\t" << format_number(t.skin_depth * 1000.0, 3) + " mm" << "\n";
            o << "Magnetism\t" << fc::to_string(t.magnetic.order);
            if (finite(t.magnetic.ordering_temperature)) o << " (orders below " << u.format(t.magnetic.ordering_temperature, "temperature", 4) << ")";
            o << "\n";
            if (!t.magnetic.note.empty()) o << t.magnetic.note << "\n";
            for (const std::string& n : t.notes) o << n << "\n";
        } catch (...) {
        }
        try {
            o << susceptibility_readouts(fc::element_susceptibility(e, m_t));
            o << electron_gas_readouts(fc::element_electron_gas(e, m_t), u);
            o << dielectric_readouts(fc::element_dielectric(e, m_t, 101325.0));
            o << mechanics_readouts(fc::element_mechanics(e, m_t));
        } catch (...) {
        }
        if (e.discovery_year > 0) o << "Discovered\t" << e.discovery_year << "\n";
        if (!e.is_stable) o << "!radioactive (no stable isotopes)\n";
        return o.str();
    }

    ChemShared& m_shared;
    int m_z = 8;
    int m_charge = 0;
    int m_extra = 0;
    double m_t = 298.15;
    std::string m_style = "both";
    std::string m_nucleus = "circle";
    double m_time = 0.0;
    double m_yaw = 0.4;
    double m_pitch = 0.3;
    bool m_down = false;
    Vec2 m_last;
};

} // namespace

std::unique_ptr<ChemView> make_structure_view(ChemShared& shared) {
    return std::unique_ptr<ChemView>(new StructureView(shared));
}

std::unique_ptr<ChemView> make_element_view(ChemShared& shared) {
    return std::unique_ptr<ChemView>(new ElementView(shared));
}

} // namespace sapling

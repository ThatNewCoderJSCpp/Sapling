#include "molecule_draw.hpp"
#include "view_util.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <set>
#include <sstream>

namespace sapling {

using namespace chem_util;

namespace {

namespace fc = fizmo::chemistry;

constexpr double NETWORK_PI = 3.14159265358979323846;

struct NetAtom {
    int z = 6;
    Vec3 p;
};

struct NetBond {
    int a = 0;
    int b = 0;
    int order = 1;
    bool weak = false;
};

struct NetworkInfo {
    const char* key;
    const char* name;
    const char* formula;
    const char* kind;
    const char* bonding;
    const char* hybridization;
    const char* bond_length;
    const char* bond_angle;
    const char* coordination;
    const char* density;
    const char* hardness;
    const char* electrical;
    const char* thermal;
    const char* melting;
    const char* uses;
    const char* why;
};

const NetworkInfo NETWORKS[] = {
    { "diamond", "Diamond", "C", "carbon allotrope (covalent network)", "every carbon makes four single bonds to four other carbons, so the whole crystal is one giant molecule", "sp\xC2\xB3", "154.4 pm", "109.5\xC2\xB0 (tetrahedral)", "4", "3.51 g/cm\xC2\xB3", "10 on the Mohs scale (the hardest natural material)", "insulator: every valence electron is locked in a bond (band gap 5.5 eV)", "about 2200 W/(m\xC2\xB7K), the best of any bulk solid, carried by lattice vibrations", "turns to graphite above about 1900 K without air; burns in air above about 1000 K", "cutting and drilling tools, abrasives, heat spreaders, gems", "to scratch or melt diamond you have to break strong C\xE2\x80\x93" "C bonds in every direction at once" },
    { "graphite", "Graphite", "C", "carbon allotrope (layered covalent network)", "each carbon makes three bonds in a flat hexagonal sheet; the fourth electron joins a delocalized \xCF\x80 system spread over the whole sheet. The sheets are held together only by weak London forces", "sp\xC2\xB2", "142 pm in a sheet, 335 pm between sheets", "120\xC2\xB0 (trigonal planar)", "3 in the sheet", "2.26 g/cm\xC2\xB3", "1\xE2\x80\x93" "2 on the Mohs scale: the sheets slide past each other", "conducts along the sheets (delocalized \xCF\x80 electrons), barely across them", "very high along the sheets, low across them", "sublimes at about 3900 K at 1 atm (the most stable form of carbon at normal pressure)", "pencils, lubricants, battery anodes, electrodes, crucibles", "strong bonds inside the sheets but weak forces between them: hard to melt, easy to smear" },
    { "graphene", "Graphene", "C", "single sheet of graphite", "one atom-thick honeycomb sheet; every carbon has three \xCF\x83 bonds and shares its fourth electron in a \xCF\x80 system over the whole sheet", "sp\xC2\xB2", "142 pm", "120\xC2\xB0", "3", "about 0.77 mg/m\xC2\xB2 (one layer)", "stronger than steel for its weight (tensile strength about 130 GPa)", "excellent conductor; its electrons behave as if massless", "about 5000 W/(m\xC2\xB7K) along the sheet", "the sheet is stable to well above 3000 K", "research, sensors, transparent conductors, composites", "a single layer of graphite: all the strength of the sheet, none of the weak interlayer forces" },
    { "nanotube", "Carbon nanotube (6,6)", "C", "graphene rolled into a tube", "a graphene sheet rolled into a seamless cylinder; this (6,6) armchair tube is about 0.8 nm across", "sp\xC2\xB2 (slightly bent)", "about 142 pm", "close to 120\xC2\xB0", "3", "about 1.3 g/cm\xC2\xB3 for bundles", "very stiff and strong along the tube", "armchair tubes like this one conduct like metals; other twists are semiconductors", "very high along the tube", "stable to about 2800 K without air", "composites, wires, field emitters, research", "rolling the sheet keeps all the strong sp\xC2\xB2 bonds, so the tube is strong along its length" },
    { "fullerene", "Buckminsterfullerene C\xE2\x82\x86\xE2\x82\x80", "C60", "carbon allotrope (molecular)", "60 carbons in 12 pentagons and 20 hexagons, like a soccer ball; unlike the others it's a real molecule, and solid C\xE2\x82\x86\xE2\x82\x80 is held together by London forces", "sp\xC2\xB2 (curved)", "140 pm (between hexagons) and 145 pm (hexagon to pentagon)", "108\xC2\xB0 in pentagons, 120\xC2\xB0 in hexagons", "3", "1.65 g/cm\xC2\xB3", "soft (a molecular solid)", "semiconductor; superconducts when doped with potassium", "low", "sublimes near 800 K", "research, organic solar cells, antioxidants", "a closed cage: every carbon is bonded, so there are no dangling edges" },
    { "lignite", "Coal: lignite (brown coal)", "C\xE2\x82\x86\xE2\x82\x80H\xE2\x82\x86\xE2\x82\x80O\xE2\x82\x82\xE2\x82\x80 (roughly; it varies)", "coal (a disordered organic solid)", "small aromatic ring clusters (sp\xC2\xB2 carbon, like tiny pieces of graphite) joined to each other and decorated with hydrogen, OH and other oxygen groups; the clusters are jumbled, not stacked", "mostly sp\xC2\xB2 in the rings, sp\xC2\xB3 in the links", "140\xE2\x80\x93" "154 pm", "120\xC2\xB0 in rings, 109.5\xC2\xB0 in links", "3 or 4", "about 1.2 g/cm\xC2\xB3", "soft and crumbly", "insulator", "low", "doesn't melt: heating drives off water and gases, then it burns", "power stations (lowest energy coal)", "about 65\xE2\x80\x93" "70% carbon: the youngest coal, with small ring clusters and lots of oxygen and water" },
    { "bituminous", "Coal: bituminous", "about 80\xE2\x80\x93" "85% C, with H, O, N and S", "coal (a disordered organic solid)", "larger aromatic ring clusters, partly stacked, with hydrogen and some oxygen groups at their edges; this is a simplified picture, real coal is far less regular", "mostly sp\xC2\xB2, some sp\xC2\xB3", "140\xE2\x80\x93" "154 pm", "120\xC2\xB0 and 109.5\xC2\xB0", "3 or 4", "about 1.3 g/cm\xC2\xB3", "medium", "insulator", "low", "softens around 600\xE2\x80\x93" "800 K and gives off tar and gas (that's how coke is made)", "power, coke for making steel", "older and more heated than lignite: oxygen has been driven off and the ring clusters have grown" },
    { "anthracite", "Coal: anthracite", "over 90% C", "coal (a disordered organic solid)", "big aromatic sheets, mostly stacked like disordered graphite, with little hydrogen or oxygen left", "almost all sp\xC2\xB2", "about 142 pm", "120\xC2\xB0", "3", "about 1.4\xE2\x80\x93" "1.5 g/cm\xC2\xB3", "hard for a coal, shiny", "slightly conductive (the sheets are big enough to start acting like graphite)", "low", "burns slowly with little smoke", "home heating, water filters", "the oldest, most carbon-rich coal: on its way to becoming graphite" },
    { "silicon", "Silicon", "Si", "covalent network (same structure as diamond)", "every silicon makes four bonds in the diamond structure; the bonds are longer and weaker than carbon's", "sp\xC2\xB3", "235 pm", "109.5\xC2\xB0", "4", "2.33 g/cm\xC2\xB3", "6.5 on the Mohs scale", "semiconductor (band gap 1.12 eV): the basis of computer chips", "about 150 W/(m\xC2\xB7K)", "melts at 1687 K", "chips, solar cells, alloys", "same network as diamond, but the bigger atoms bond less strongly, so the band gap is small enough for heat and light to free some electrons" },
    { "sic", "Silicon carbide (3C, carborundum)", "SiC", "covalent network (zinc blende: diamond structure with alternating Si and C)", "every Si bonds to four C and every C to four Si, so it's diamond-like but with two kinds of atoms", "sp\xC2\xB3", "189 pm", "109.5\xC2\xB0", "4", "3.21 g/cm\xC2\xB3", "9\xE2\x80\x93" "9.5 on the Mohs scale", "wide-gap semiconductor (2.4 eV), used for high-voltage electronics", "about 360 W/(m\xC2\xB7K)", "decomposes at about 2800 K", "abrasives, brake discs, power electronics, armor", "strong, slightly polar Si\xE2\x80\x93" "C bonds in a 3D network" },
    { "cbn", "Cubic boron nitride", "BN", "covalent network (zinc blende)", "boron and nitrogen alternate in the diamond structure; B\xE2\x80\x93" "N has the same number of electrons as C\xE2\x80\x93" "C", "sp\xC2\xB3", "157 pm", "109.5\xC2\xB0", "4", "3.45 g/cm\xC2\xB3", "about 9.5 on the Mohs scale (second only to diamond)", "insulator (band gap about 6 eV)", "about 740 W/(m\xC2\xB7K)", "stable in air to about 1600 K, better than diamond", "cutting steel, where diamond would react with the iron", "a diamond look-alike: B and N together have the same electrons as two carbons" },
    { "hbn", "Hexagonal boron nitride (white graphite)", "BN", "layered covalent network", "flat honeycomb sheets of alternating B and N held together by weak forces, like graphite; B\xE2\x80\x93" "N bonds are polar, so the \xCF\x80 electrons stay on the nitrogens", "sp\xC2\xB2", "145 pm in a sheet, 333 pm between sheets", "120\xC2\xB0", "3 in the sheet", "2.1 g/cm\xC2\xB3", "soft and slippery", "insulator, unlike graphite: the polar bonds keep the \xCF\x80 electrons from moving across the sheet", "high along the sheets", "stable to about 3200 K", "lubricants, cosmetics, crucibles, a base layer for graphene devices", "graphite's structure with B and N, but the electrons can't roam" }
};

const NetworkInfo* network_info(const std::string& key) {
    for (const NetworkInfo& n : NETWORKS) {
        if (key == n.key) return &n;
    }
    return nullptr;
}

std::string fixed_text(double v, int digits) {
    std::ostringstream o;
    o.setf(std::ios::fixed);
    o.precision(digits);
    o << v;
    return o.str();
}

class NetworkView : public ChemView {
public:
    explicit NetworkView(ChemShared& shared) : m_shared(shared) {
        m_cam.yaw = 0.55;
        m_cam.pitch = 0.35;
        build();
    }

    void draw(Canvas& canvas, const Theme& theme) override {
        canvas.clear(theme.background);
        m_cam.width = canvas.width();
        m_cam.height = canvas.height();
        if (!m_framed && !m_atoms.empty()) {
            std::vector<Vec3> pts;
            for (const NetAtom& a : m_atoms) pts.push_back(a.p);
            m_cam.target = Vec3{};
            m_cam.pan = Vec2{};
            m_cam.frame(pts, 40.0);
            m_framed = true;
        }
        Scene3D scene(m_cam);
        const double scale = m_atoms.size() > 900 ? 0.6 : 1.0;
        for (const NetBond& b : m_bonds) {
            const Vec3 pa = m_atoms[static_cast<std::size_t>(b.a)].p, pb = m_atoms[static_cast<std::size_t>(b.b)].p;
            if (b.weak) {
                if (m_show_weak) scene.dashed(pa, pb, 1.0, rgb(160, 160, 190, 110));
                continue;
            }
            scene.cylinder(pa, pb, 16.0 * scale, element_color(m_atoms[static_cast<std::size_t>(b.a)].z), element_color(m_atoms[static_cast<std::size_t>(b.b)].z));
        }
        for (std::size_t i = 0; i < m_atoms.size(); ++i) {
            const NetAtom& a = m_atoms[i];
            const double r = (a.z == 1 ? 18.0 : a.z == 6 ? 30.0 : a.z == 8 ? 28.0 : 34.0) * scale;
            Color c = element_color(a.z);
            scene.sphere(a.p, r, c);
            if (static_cast<int>(i) == m_selected) scene.ring(a.p, 10.0, 2.0, theme.accent, 1.0);
        }
        scene.draw(canvas);
        const NetworkInfo* info = network_info(m_key);
        draw_hud(canvas, theme, { info ? info->name : m_key, std::to_string(m_atoms.size()) + " atoms, " + std::to_string(strong_bonds()) + " bonds; drag to turn, click an atom" }, false);
    }

    bool mouse_down(Vec2 p, int, unsigned) override {
        m_down = true;
        m_moved = false;
        m_last = p;
        return false;
    }

    bool mouse_move(Vec2 p, bool) override {
        if (!m_down) return false;
        const Vec2 d = p - m_last;
        m_last = p;
        if (std::fabs(d.x) + std::fabs(d.y) > 0.5) m_moved = true;
        m_cam.orbit(d.x, d.y);
        return true;
    }

    bool mouse_up(Vec2 p, int) override {
        const bool click = m_down && !m_moved;
        m_down = false;
        if (!click) return false;
        int best = -1;
        double best_d = 12.0;
        for (std::size_t i = 0; i < m_atoms.size(); ++i) {
            const Vec2 s = m_cam.project(m_atoms[i].p);
            const double d = length(s - p);
            if (d < best_d) {
                best_d = d;
                best = static_cast<int>(i);
            }
        }
        m_selected = best;
        ++m_version;
        return true;
    }

    bool wheel(Vec2, double delta) override {
        m_cam.zoom_by(delta > 0 ? 1.1 : 1.0 / 1.1);
        return true;
    }

    std::string command(const std::vector<std::string>& a) override {
        const std::string c = arg(a, 0);
        if (c == "structure") {
            if (network_info(arg(a, 1))) m_key = arg(a, 1);
            build();
            return "ok\t" + m_key;
        }
        if (c == "size") {
            m_size = std::max(1, std::min(6, static_cast<int>(num(arg(a, 1), 2.0))));
            build();
            return "ok";
        }
        if (c == "weak") {
            m_show_weak = arg(a, 1) == "1";
            return "ok";
        }
        if (c == "tp") {
            m_temperature = std::max(1.0, fc::to_si(num(arg(a, 1), 298.15), arg(a, 2).empty() ? "K" : arg(a, 2), fc::Dimension::Temperature));
            m_pressure = std::max(1e-6, fc::to_si(num(arg(a, 3), 1.0), arg(a, 4).empty() ? "atm" : arg(a, 4), fc::Dimension::Pressure));
            return "ok";
        }
        if (c == "list") {
            std::string out;
            for (const NetworkInfo& n : NETWORKS) out += std::string(n.key) + "\t" + n.name + "\n";
            return out;
        }
        if (c == "version") return std::to_string(m_version);
        if (c == "readouts") return readouts();
        if (c == "selection") return selection();
        return "error\tunknown command " + c;
    }

private:
    ChemShared& m_shared;
    std::string m_key = "diamond";
    int m_size = 2;
    std::vector<NetAtom> m_atoms;
    std::vector<NetBond> m_bonds;
    Camera3D m_cam;
    bool m_framed = false;
    bool m_down = false;
    bool m_moved = false;
    Vec2 m_last;
    int m_selected = -1;
    bool m_show_weak = true;
    double m_temperature = 298.15;
    double m_pressure = 101325.0;
    int m_version = 0;

    int strong_bonds() const {
        int n = 0;
        for (const NetBond& b : m_bonds) n += b.weak ? 0 : 1;
        return n;
    }

    void center() {
        if (m_atoms.empty()) return;
        Vec3 c;
        for (const NetAtom& a : m_atoms) c = c + a.p;
        c = c * (1.0 / static_cast<double>(m_atoms.size()));
        for (NetAtom& a : m_atoms) a.p = a.p - c;
    }

    void connect(double lo, double hi, bool weak) {
        const std::size_t n = m_atoms.size();
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                const double d = length(m_atoms[i].p - m_atoms[j].p);
                if (d >= lo && d <= hi) m_bonds.push_back(NetBond{ static_cast<int>(i), static_cast<int>(j), 1, weak });
            }
        }
    }

    void diamond_type(int za, int zb, double a) {
        const double f[4][3] = { { 0, 0, 0 }, { 0, 0.5, 0.5 }, { 0.5, 0, 0.5 }, { 0.5, 0.5, 0 } };
        const int n = m_size + 1;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    for (const auto& s : f) {
                        const Vec3 base{ (i + s[0]) * a, (j + s[1]) * a, (k + s[2]) * a };
                        if (base.x <= m_size * a + 1e-6 && base.y <= m_size * a + 1e-6 && base.z <= m_size * a + 1e-6) m_atoms.push_back(NetAtom{ za, base });
                        const Vec3 second = base + Vec3{ 0.25 * a, 0.25 * a, 0.25 * a };
                        if (second.x <= m_size * a && second.y <= m_size * a && second.z <= m_size * a) m_atoms.push_back(NetAtom{ zb, second });
                    }
                }
            }
        }
        const double bond = a * std::sqrt(3.0) / 4.0;
        connect(bond * 0.9, bond * 1.1, false);
        prune();
    }

    void prune() {
        std::vector<int> degree(m_atoms.size(), 0);
        for (const NetBond& b : m_bonds) {
            if (b.weak) continue;
            ++degree[static_cast<std::size_t>(b.a)];
            ++degree[static_cast<std::size_t>(b.b)];
        }
        std::vector<int> remap(m_atoms.size(), -1);
        std::vector<NetAtom> kept;
        for (std::size_t i = 0; i < m_atoms.size(); ++i) {
            if (degree[i] == 0) continue;
            remap[i] = static_cast<int>(kept.size());
            kept.push_back(m_atoms[i]);
        }
        std::vector<NetBond> bonds;
        for (NetBond b : m_bonds) {
            if (remap[static_cast<std::size_t>(b.a)] < 0 || remap[static_cast<std::size_t>(b.b)] < 0) continue;
            b.a = remap[static_cast<std::size_t>(b.a)];
            b.b = remap[static_cast<std::size_t>(b.b)];
            bonds.push_back(b);
        }
        m_atoms = kept;
        m_bonds = bonds;
    }

    void honeycomb_layer(std::vector<NetAtom>& out, int za, int zb, double d, int rings, double z, Vec3 shift) const {
        const Vec3 a1{ std::sqrt(3.0) * d, 0.0, 0.0 };
        const Vec3 a2{ std::sqrt(3.0) * d / 2.0, 1.5 * d, 0.0 };
        for (int i = -rings; i <= rings; ++i) {
            for (int j = -rings; j <= rings; ++j) {
                const Vec3 cell = a1 * static_cast<double>(i) + a2 * static_cast<double>(j) + shift;
                const double limit = (rings + 0.3) * std::sqrt(3.0) * d;
                for (int s = 0; s < 2; ++s) {
                    const Vec3 p = cell + Vec3{ 0.0, s == 0 ? 0.0 : d, z };
                    if (std::fabs(p.x) > limit || std::fabs(p.y) > limit * 0.9) continue;
                    out.push_back(NetAtom{ s == 0 ? za : zb, p });
                }
            }
        }
    }

    void layered(int za, int zb, double d, double spacing, int layers, bool swap_layers) {
        const int rings = 2 + m_size;
        for (int l = 0; l < layers; ++l) {
            std::vector<NetAtom> layer;
            const bool odd = l % 2 == 1;
            const Vec3 shift = (odd && !swap_layers) ? Vec3{ 0.0, d, 0.0 } : Vec3{};
            honeycomb_layer(layer, odd && swap_layers ? zb : za, odd && swap_layers ? za : zb, d, rings, (l - 0.5 * (layers - 1)) * spacing, shift);
            for (const NetAtom& a : layer) m_atoms.push_back(a);
        }
        connect(d * 0.9, d * 1.1, false);
        if (layers > 1) {
            const std::size_t n = m_atoms.size();
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = i + 1; j < n; ++j) {
                    const Vec3 dv = m_atoms[i].p - m_atoms[j].p;
                    if (std::fabs(std::fabs(dv.z) - spacing) < 1.0 && std::fabs(dv.x) < 1.0 && std::fabs(dv.y) < 1.0) m_bonds.push_back(NetBond{ static_cast<int>(i), static_cast<int>(j), 1, true });
                }
            }
        }
    }

    void nanotube() {
        const double d = 142.0;
        const int n = 6;
        const Vec3 a1{ std::sqrt(3.0) * d, 0.0, 0.0 };
        const Vec3 a2{ std::sqrt(3.0) * d / 2.0, 1.5 * d, 0.0 };
        const Vec3 chiral = (a1 + a2) * static_cast<double>(n);
        const double circumference = length(chiral);
        const Vec3 u = chiral * (1.0 / circumference);
        const Vec3 v{ -u.y, u.x, 0.0 };
        const double tube_length = (4 + 2 * m_size) * std::sqrt(3.0) * d;
        const double radius = circumference / (2.0 * NETWORK_PI);
        std::vector<Vec3> flat;
        for (int i = -40; i <= 40; ++i) {
            for (int j = -40; j <= 40; ++j) {
                const Vec3 cell = a1 * static_cast<double>(i) + a2 * static_cast<double>(j);
                for (int s = 0; s < 2; ++s) {
                    const Vec3 p = cell + Vec3{ 0.0, s == 0 ? 0.0 : d, 0.0 };
                    const double cu = dot(p, u), cv = dot(p, v);
                    if (cu < -1e-6 || cu >= circumference - 1e-6 || cv < 0.0 || cv > tube_length) continue;
                    flat.push_back(Vec3{ cu, cv, 0.0 });
                }
            }
        }
        for (const Vec3& f : flat) {
            const double angle = 2.0 * NETWORK_PI * f.x / circumference;
            m_atoms.push_back(NetAtom{ 6, Vec3{ radius * std::cos(angle), f.y, radius * std::sin(angle) } });
        }
        connect(d * 0.85, d * 1.12, false);
    }

    void fullerene() {
        const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
        const double base[3][3] = { { 0.0, 1.0, 3.0 * phi }, { 1.0, 2.0 + phi, 2.0 * phi }, { phi, 2.0, 2.0 * phi + 1.0 } };
        std::vector<Vec3> pts;
        for (const auto& b : base) {
            for (int sx = -1; sx <= 1; sx += 2) {
                for (int sy = -1; sy <= 1; sy += 2) {
                    for (int sz = -1; sz <= 1; sz += 2) {
                        const double v[3] = { b[0] * sx, b[1] * sy, b[2] * sz };
                        for (int r = 0; r < 3; ++r) {
                            const Vec3 p{ v[r % 3], v[(r + 1) % 3], v[(r + 2) % 3] };
                            bool dup = false;
                            for (const Vec3& q : pts) dup = dup || length(q - p) < 1e-6;
                            if (!dup) pts.push_back(p);
                        }
                    }
                }
            }
        }
        const double scale = 143.0 / 2.0;
        for (const Vec3& p : pts) m_atoms.push_back(NetAtom{ 6, p * scale });
        connect(130.0, 155.0, false);
    }

    void coal(int rank) {
        std::mt19937 rng(static_cast<unsigned>(17 + rank));
        std::uniform_real_distribution<double> U(0.0, 1.0);
        const double d = 140.0;
        const int fragments = rank == 0 ? 7 : rank == 1 ? 6 : 4;
        const int min_rings = rank == 0 ? 1 : rank == 1 ? 3 : 8;
        const int max_rings = rank == 0 ? 3 : rank == 1 ? 6 : 14;
        const double oxygen = rank == 0 ? 0.35 : rank == 1 ? 0.12 : 0.03;
        const double hydrogen = rank == 0 ? 0.9 : rank == 1 ? 0.8 : 0.5;
        for (int f = 0; f < fragments; ++f) {
            const int rings = min_rings + static_cast<int>(U(rng) * (max_rings - min_rings + 1));
            std::set<std::pair<int, int>> centers{ { 0, 0 } };
            while (static_cast<int>(centers.size()) < rings) {
                auto it = centers.begin();
                std::advance(it, static_cast<long>(U(rng) * centers.size()) % static_cast<long>(centers.size()));
                static const int steps[6][2] = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 }, { 1, -1 }, { -1, 1 } };
                const int s = static_cast<int>(U(rng) * 6) % 6;
                centers.insert({ it->first + steps[s][0], it->second + steps[s][1] });
            }
            std::vector<Vec3> local;
            const Vec3 a1{ std::sqrt(3.0) * d, 0.0, 0.0 };
            const Vec3 a2{ std::sqrt(3.0) * d / 2.0, 1.5 * d, 0.0 };
            for (const auto& c : centers) {
                const Vec3 cc = a1 * static_cast<double>(c.first) + a2 * static_cast<double>(c.second);
                for (int k = 0; k < 6; ++k) {
                    const double ang = NETWORK_PI / 6.0 + k * NETWORK_PI / 3.0;
                    const Vec3 p = cc + Vec3{ d * std::cos(ang), d * std::sin(ang), 0.0 };
                    bool dup = false;
                    for (const Vec3& q : local) dup = dup || length(q - p) < 1.0;
                    if (!dup) local.push_back(p);
                }
            }
            const double yaw = U(rng) * 2.0 * NETWORK_PI;
            const double tilt = rank == 2 ? (U(rng) - 0.5) * 0.25 : (U(rng) - 0.5) * 1.6;
            const double spread = rank == 2 ? 600.0 : 1100.0;
            const Vec3 offset = rank == 2 ? Vec3{ (U(rng) - 0.5) * 300.0, (U(rng) - 0.5) * 300.0, (f - 0.5 * (fragments - 1)) * 345.0 } : Vec3{ (U(rng) - 0.5) * spread, (U(rng) - 0.5) * spread, (U(rng) - 0.5) * spread };
            const std::size_t first = m_atoms.size();
            for (const Vec3& p : local) {
                const double cy = std::cos(yaw), sy = std::sin(yaw), ct = std::cos(tilt), st = std::sin(tilt);
                const Vec3 r1{ cy * p.x - sy * p.y, sy * p.x + cy * p.y, p.z };
                const Vec3 r2{ r1.x, ct * r1.y - st * r1.z, st * r1.y + ct * r1.z };
                m_atoms.push_back(NetAtom{ 6, r2 + offset });
            }
            const std::size_t last = m_atoms.size();
            for (std::size_t i = first; i < last; ++i) {
                for (std::size_t j = i + 1; j < last; ++j) {
                    const double dist = length(m_atoms[i].p - m_atoms[j].p);
                    if (dist > d * 0.9 && dist < d * 1.1) m_bonds.push_back(NetBond{ static_cast<int>(i), static_cast<int>(j), 1, false });
                }
            }
            Vec3 centroid;
            for (std::size_t i = first; i < last; ++i) centroid = centroid + m_atoms[i].p;
            centroid = centroid * (1.0 / static_cast<double>(last - first));
            for (std::size_t i = first; i < last; ++i) {
                int degree = 0;
                for (const NetBond& b : m_bonds) degree += (b.a == static_cast<int>(i) || b.b == static_cast<int>(i)) && !b.weak;
                if (degree != 2) continue;
                const Vec3 out = normalized(m_atoms[i].p - centroid);
                const double roll = U(rng);
                if (roll < oxygen) {
                    const int o = static_cast<int>(m_atoms.size());
                    m_atoms.push_back(NetAtom{ 8, m_atoms[i].p + out * 136.0 });
                    m_bonds.push_back(NetBond{ static_cast<int>(i), o, 1, false });
                    if (U(rng) < 0.7) {
                        const int h = static_cast<int>(m_atoms.size());
                        m_atoms.push_back(NetAtom{ 1, m_atoms[static_cast<std::size_t>(o)].p + out * 96.0 });
                        m_bonds.push_back(NetBond{ o, h, 1, false });
                    }
                } else if (roll < oxygen + hydrogen) {
                    const int h = static_cast<int>(m_atoms.size());
                    m_atoms.push_back(NetAtom{ 1, m_atoms[i].p + out * 108.0 });
                    m_bonds.push_back(NetBond{ static_cast<int>(i), h, 1, false });
                }
            }
        }
        if (rank < 2) {
            std::vector<int> carbons;
            for (std::size_t i = 0; i < m_atoms.size(); ++i) if (m_atoms[i].z == 6) carbons.push_back(static_cast<int>(i));
            for (int link = 0; link < fragments; ++link) {
                const int a = carbons[static_cast<std::size_t>(U(rng) * carbons.size()) % carbons.size()];
                int best = -1;
                double best_d = 1e300;
                for (int b : carbons) {
                    if (b == a) continue;
                    const double dist = length(m_atoms[static_cast<std::size_t>(a)].p - m_atoms[static_cast<std::size_t>(b)].p);
                    bool same = false;
                    for (const NetBond& nb : m_bonds) same = same || ((nb.a == a && nb.b == b) || (nb.a == b && nb.b == a));
                    if (dist > 250.0 && dist < best_d && !same) {
                        best_d = dist;
                        best = b;
                    }
                }
                if (best < 0 || best_d > 700.0) continue;
                const Vec3 mid = (m_atoms[static_cast<std::size_t>(a)].p + m_atoms[static_cast<std::size_t>(best)].p) * 0.5;
                const int link_atom = static_cast<int>(m_atoms.size());
                m_atoms.push_back(NetAtom{ U(rng) < 0.3 ? 8 : 6, mid });
                m_bonds.push_back(NetBond{ a, link_atom, 1, false });
                m_bonds.push_back(NetBond{ link_atom, best, 1, false });
            }
        }
    }

    void build() {
        ++m_version;
        m_atoms.clear();
        m_bonds.clear();
        m_selected = -1;
        m_framed = false;
        if (m_key == "diamond") diamond_type(6, 6, 356.7);
        else if (m_key == "silicon") diamond_type(14, 14, 543.1);
        else if (m_key == "sic") diamond_type(14, 6, 435.9);
        else if (m_key == "cbn") diamond_type(5, 7, 361.5);
        else if (m_key == "graphite") layered(6, 6, 142.0, 335.0, 3, false);
        else if (m_key == "graphene") layered(6, 6, 142.0, 335.0, 1, false);
        else if (m_key == "hbn") layered(5, 7, 144.6, 333.0, 3, true);
        else if (m_key == "nanotube") nanotube();
        else if (m_key == "fullerene") fullerene();
        else if (m_key == "lignite") coal(0);
        else if (m_key == "bituminous") coal(1);
        else if (m_key == "anthracite") coal(2);
        center();
    }

    std::vector<int> neighbors_of(int i) const {
        std::vector<int> out;
        for (const NetBond& b : m_bonds) {
            if (b.weak) continue;
            if (b.a == i) out.push_back(b.b);
            else if (b.b == i) out.push_back(b.a);
        }
        return out;
    }

    std::string readouts() const {
        std::ostringstream o;
        const UnitPrefs& u = m_shared.units;
        const NetworkInfo* n = network_info(m_key);
        if (!n) return "#Network solid\n!pick a structure\n";
        o << "#" << n->name << "\n";
        o << "Formula\t" << n->formula << "\n";
        o << "Kind\t" << n->kind << "\n";
        o << n->bonding << "\n";
        o << "#Bonds\n";
        o << "Hybridization\t" << n->hybridization << "\n";
        o << "Bond length\t" << n->bond_length << "\n";
        o << "Bond angle\t" << n->bond_angle << "\n";
        o << "Neighbors per atom\t" << n->coordination << "\n";
        double total = 0.0;
        int count = 0;
        std::map<std::string, int> elements;
        for (const NetAtom& a : m_atoms) ++elements[fc::element(a.z).symbol];
        for (const NetBond& b : m_bonds) {
            if (b.weak) continue;
            total += length(m_atoms[static_cast<std::size_t>(b.a)].p - m_atoms[static_cast<std::size_t>(b.b)].p);
            ++count;
        }
        if (count > 0) o << "Average bond shown\t" << fixed_text(total / count, 1) << " pm\n";
        std::string list;
        for (const auto& kv : elements) list += (list.empty() ? "" : ", ") + std::to_string(kv.second) + " " + kv.first;
        o << "Atoms shown\t" << list << "\n";
        o << "#Properties\n";
        o << "Density\t" << n->density << "\n";
        o << "Hardness\t" << n->hardness << "\n";
        o << "Electricity\t" << n->electrical << "\n";
        o << "Heat conduction\t" << n->thermal << "\n";
        o << "Heating\t" << n->melting << "\n";
        o << "Uses\t" << n->uses << "\n";
        o << "Why\t" << n->why << "\n";
        const std::string symbol_c = "C";
        if (elements.size() == 1 && elements.begin()->first == symbol_c && m_key != "fullerene") {
            const double t = m_temperature;
            const double p_gpa = m_pressure / 1e9;
            const double line = 1.6 + 0.0027 * t;
            o << "#Which carbon is stable at " << u.format(t, "temperature", 4) << " and " << u.format(m_pressure, "pressure", 4) << "\n";
            if (t >= 4500.0 && p_gpa < 0.01) o << "Stable form\tcarbon vapor (graphite sublimes near 3900 K at 1 atm)\n";
            else if (t >= 4800.0) o << "Stable form\tliquid carbon\n";
            else if (p_gpa > line) o << "Stable form\tdiamond (the pressure is above the graphite\xE2\x80\x93" "diamond line, about " << fixed_text(line, 1) << " GPa here)\n";
            else o << "Stable form\tgraphite (diamond needs more than about " << fixed_text(line, 1) << " GPa at this temperature)\n";
            if (m_key == "diamond" && p_gpa <= line) o << "Diamond here is metastable: it should turn into graphite, but the bonds are so strong that it takes thousands of degrees to start. Diamonds really are not forever, just very slow.\n";
            if (m_key == "graphite" && p_gpa > line) o << "Graphite wants to become diamond here, but without heat (or a metal catalyst) the change is too slow; synthetic diamonds are made at about 5 GPa and 1700 K.\n";
            if (m_key == "lignite" || m_key == "bituminous" || m_key == "anthracite") o << "Coal turns into graphite only after millions of years of heat and pressure; anthracite is the closest to it.\n";
        }
        o << "Click an atom to see its bonds.\n";
        return o.str();
    }

    std::string selection() const {
        if (m_selected < 0 || static_cast<std::size_t>(m_selected) >= m_atoms.size()) return "";
        std::ostringstream o;
        const NetAtom& a = m_atoms[static_cast<std::size_t>(m_selected)];
        const std::vector<int> nb = neighbors_of(m_selected);
        o << "#Atom " << (m_selected + 1) << ": " << fc::element(a.z).name << "\n";
        o << "Bonded neighbors\t" << nb.size() << "\n";
        std::string kind = nb.size() == 4 ? "sp\xC2\xB3 (tetrahedral)" : nb.size() == 3 ? "sp\xC2\xB2 (trigonal planar)" : nb.size() == 2 ? (a.z == 6 ? "an edge carbon (its other bond is cut off at the edge of the picture)" : "two bonds") : nb.size() == 1 ? "a single bond" : "none";
        if (a.z == 1) kind = "hydrogen on an edge";
        o << "Bonding\t" << kind << "\n";
        for (std::size_t i = 0; i < nb.size(); ++i) {
            const NetAtom& b = m_atoms[static_cast<std::size_t>(nb[i])];
            o << "To " << fc::element(b.z).symbol << " " << (nb[i] + 1) << "\t" << fixed_text(length(a.p - b.p), 1) << " pm\n";
        }
        if (nb.size() >= 2) {
            const Vec3 u0 = normalized(m_atoms[static_cast<std::size_t>(nb[0])].p - a.p);
            const Vec3 u1 = normalized(m_atoms[static_cast<std::size_t>(nb[1])].p - a.p);
            o << "Angle\t" << fixed_text(std::acos(std::max(-1.0, std::min(1.0, dot(u0, u1)))) * 180.0 / NETWORK_PI, 1) << "\xC2\xB0\n";
        }
        return o.str();
    }
};

}

std::unique_ptr<ChemView> make_network_view(ChemShared& shared) {
    return std::make_unique<NetworkView>(shared);
}

} // namespace sapling

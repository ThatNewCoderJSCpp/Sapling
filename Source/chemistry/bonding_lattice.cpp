#include "bonding_view.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <optional>
#include <sstream>

namespace sapling {

namespace {

constexpr double PI = 3.14159265358979323846;

bool finite(double v) {
    return v == v && std::isfinite(v);
}

std::string trim(std::string s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
    return s;
}

std::optional<fc::IonDescriptor> ion_from_text(const std::string& text, bool cation) {
    const std::string t = trim(text);
    if (t.empty()) return std::nullopt;
    const fc::Element* e = fc::find_element(t);
    if (e) {
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

double ease(double t) {
    t = std::max(0.0, std::min(1.0, t));
    return t * t * (3.0 - 2.0 * t);
}

std::string fixed(double v, int digits) {
    if (!finite(v)) return "n/a";
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*f", digits, v);
    return buf;
}

} // namespace

void BondingView::set_ionic(const std::string& cation, const std::string& anion) {
    m_lattice_error.clear();
    if (!cation.empty()) m_cation = trim(cation);
    if (!anion.empty()) m_anion = trim(anion);
    fc::IonicCrystal crystal;
    try {
        const fc::Element* ce = fc::find_element(m_cation);
        const fc::Element* ae = fc::find_element(m_anion);
        if (ce && ae) {
            crystal = fc::ionic_crystal(*ce, *ae);
        } else {
            const auto c = ion_from_text(m_cation, true);
            const auto a = ion_from_text(m_anion, false);
            if (!c) m_lattice_error = "'" + m_cation + "' isn't a cation I know (try Na, Fe3+, NH4+, Hg2 2+)";
            else if (!a) m_lattice_error = "'" + m_anion + "' isn't an anion I know (try Cl, O, SO4 2-, nitrate)";
            else crystal = fc::ionic_crystal(*c, *a);
        }
    } catch (const std::exception& e) {
        m_lattice_error = e.what();
    }
    if (!crystal.ok && m_lattice_error.empty()) {
        m_lattice_error = crystal.notes.empty() ? "these two don't make an ionic compound" : crystal.notes.front();
    }
    m_ionic = crystal;
    m_sel_lattice = -1;
    rebuild_lattice();
    restart_assembly();
    frame_view();
    changed(true);
}

void BondingView::set_metal(const std::string& spec) {
    m_lattice_error.clear();
    std::vector<std::pair<int, double>> comps;
    std::string cur;
    std::vector<std::string> parts;
    for (char ch : spec + ";") {
        if (ch == ';' || ch == ',' || ch == ' ') {
            if (!trim(cur).empty()) parts.push_back(trim(cur));
            cur.clear();
        } else {
            cur += ch;
        }
    }
    for (const std::string& p : parts) {
        const std::size_t eq = p.find_first_of("=:");
        const std::string sym = trim(eq == std::string::npos ? p : p.substr(0, eq));
        double frac = 1.0;
        if (eq != std::string::npos) {
            char* end = nullptr;
            const std::string val = p.substr(eq + 1);
            frac = std::strtod(val.c_str(), &end);
            if (end == val.c_str() || !(frac > 0.0)) frac = 0.0;
        }
        const fc::Element* e = fc::find_element(sym);
        if (!e) {
            m_lattice_error = "'" + sym + "' isn't an element";
            continue;
        }
        if (frac <= 0.0) continue;
        bool merged = false;
        for (auto& c : comps) {
            if (c.first == e->atomic_number) {
                c.second += frac;
                merged = true;
            }
        }
        if (!merged) comps.emplace_back(e->atomic_number, frac);
    }
    double total = 0.0;
    for (const auto& c : comps) total += c.second;
    if (comps.empty() || total <= 0.0) {
        if (m_lattice_error.empty()) m_lattice_error = "pick at least one metal";
        m_metal = fc::MetallicCrystal{};
        m_lattice.clear();
        changed(true);
        return;
    }
    for (auto& c : comps) c.second /= total;
    m_components = comps;
    std::vector<fc::AlloyComponent> ac;
    for (const auto& c : comps) ac.push_back(fc::AlloyComponent{ &fc::element(c.first), c.second });
    try {
        m_metal = fc::metallic_crystal(ac);
        if (!m_metal.ok) m_lattice_error = m_metal.notes.empty() ? "can't build that metal" : m_metal.notes.front();
    } catch (const std::exception& e) {
        m_metal = fc::MetallicCrystal{};
        m_lattice_error = e.what();
    }
    m_sel_lattice = -1;
    rebuild_lattice();
    restart_assembly();
    frame_view();
    changed(true);
}

void BondingView::rebuild_lattice() {
    m_lattice.clear();
    m_neighbors.clear();
    m_electrons.clear();
    std::vector<fc::CrystalAtom> atoms;
    try {
        if (m_mode == BondingMode::Ionic && m_ionic.ok && !m_ionic.cell.sites.empty()) {
            atoms = fc::build_crystal(m_ionic.cell, m_cells, m_cells, m_cells, true);
        } else if (m_mode == BondingMode::Metallic && m_metal.ok) {
            atoms = fc::build_metal(m_metal, m_cells + 1, m_cells + 1, m_cells + 1, 3u);
        }
    } catch (const std::exception& e) {
        m_lattice_error = e.what();
        return;
    }
    if (atoms.empty()) return;
    Vec3 c;
    for (const fc::CrystalAtom& a : atoms) c = c + Vec3{ a.x, a.y, a.z };
    c = c * (1.0 / static_cast<double>(atoms.size()));
    m_origin = Vec3{} - c;
    m_box_lo = Vec3{ 1e300, 1e300, 1e300 };
    m_box_hi = Vec3{ -1e300, -1e300, -1e300 };
    for (const fc::CrystalAtom& a : atoms) {
        LatticeAtom la;
        la.site = Vec3{ a.x, a.y, a.z } - c;
        la.pos = la.site;
        la.species = a.species;
        la.interstitial = a.interstitial;
        m_box_lo = Vec3{ std::min(m_box_lo.x, la.site.x), std::min(m_box_lo.y, la.site.y), std::min(m_box_lo.z, la.site.z) };
        m_box_hi = Vec3{ std::max(m_box_hi.x, la.site.x), std::max(m_box_hi.y, la.site.y), std::max(m_box_hi.z, la.site.z) };
        m_lattice.push_back(la);
    }
    if (atoms.size() <= 2500) {
        try {
            m_neighbors = fc::nearest_neighbor_pairs(atoms, 0.08);
        } catch (...) {
        }
    }
    m_nn_distance = 250.0;
    if (!m_neighbors.empty()) {
        const auto& p = m_neighbors.front();
        m_nn_distance = length(m_lattice[p.first].site - m_lattice[p.second].site);
    }
    if (m_mode == BondingMode::Metallic && m_metal.ok) {
        const double valence = finite(m_metal.valence_electrons_per_atom) ? m_metal.valence_electrons_per_atom : 1.0;
        const std::size_t count = std::min<std::size_t>(600, static_cast<std::size_t>(std::max(1.0, valence) * static_cast<double>(m_lattice.size())));
        std::uniform_real_distribution<double> U(0.0, 1.0);
        for (std::size_t i = 0; i < count; ++i) {
            LatticeAtom e;
            e.pos = Vec3{ m_box_lo.x + U(m_rng) * (m_box_hi.x - m_box_lo.x), m_box_lo.y + U(m_rng) * (m_box_hi.y - m_box_lo.y), m_box_lo.z + U(m_rng) * (m_box_hi.z - m_box_lo.z) };
            const double th = U(m_rng) * 2.0 * PI, ph = std::acos(2.0 * U(m_rng) - 1.0);
            e.vel = Vec3{ std::sin(ph) * std::cos(th), std::cos(ph), std::sin(ph) * std::sin(th) } * 400.0;
            m_electrons.push_back(e);
        }
    }
}

void BondingView::restart_assembly() {
    m_assembly = 0.0;
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    const double r = std::max(400.0, 1.6 * length(m_box_hi - m_box_lo));
    for (LatticeAtom& a : m_lattice) {
        Vec3 d;
        do {
            d = Vec3{ U(m_rng), U(m_rng), U(m_rng) };
        } while (dot(d, d) > 1.0);
        a.start = d * r;
        a.pos = a.start;
        a.phase = (U(m_rng) + 1.0) * PI;
        a.vel = Vec3{ U(m_rng), U(m_rng), U(m_rng) } * 300.0;
    }
}

std::string BondingView::lattice_state() const {
    double mp = fc::NO_VALUE, bp = fc::NO_VALUE;
    if (m_mode == BondingMode::Ionic) {
        mp = m_ionic.melting_point;
        bp = m_ionic.boiling_point;
    } else {
        mp = m_metal.melting_point;
        bp = m_metal.boiling_point;
    }
    if (finite(bp) && m_temperature >= bp) return "gas";
    if (finite(mp) && m_temperature >= mp) return "liquid";
    return "solid";
}

void BondingView::step_lattice(double dt) {
    m_time += dt;
    dt = std::min(dt, 0.05);
    if (m_lattice.empty()) return;
    const std::string state = lattice_state();
    const int code = state == "solid" ? 0 : state == "liquid" ? 1 : 2;
    if (code != m_lattice_state) {
        if (code == 0) {
            for (LatticeAtom& a : m_lattice) a.start = a.pos;
            m_assembly = 0.0;
        }
        m_lattice_state = code;
    }
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    const Vec3 half = (m_box_hi - m_box_lo) * 0.5;
    if (code == 0) {
        m_assembly = std::min(1.0, m_assembly + dt / 2.4);
        const double k = ease(m_assembly);
        double mp = m_mode == BondingMode::Ionic ? m_ionic.melting_point : m_metal.melting_point;
        if (!finite(mp)) mp = 1500.0;
        const double amp = m_options["thermal"] ? 0.07 * m_nn_distance * std::sqrt(std::max(0.0, m_temperature) / mp) : 0.0;
        for (LatticeAtom& a : m_lattice) {
            const Vec3 j{ std::sin(a.phase + m_time * 11.0), std::sin(a.phase * 1.37 + m_time * 13.0), std::sin(a.phase * 1.71 + m_time * 9.0) };
            a.pos = a.start + (a.site - a.start) * k + j * amp;
        }
    } else {
        const double spread = code == 1 ? 1.15 : 2.6;
        const double speed = code == 1 ? 140.0 : 520.0;
        const Vec3 lim = half * spread + Vec3{ m_nn_distance, m_nn_distance, m_nn_distance } * 0.5;
        for (LatticeAtom& a : m_lattice) {
            if (code == 1) a.vel = a.vel * 0.96 + Vec3{ U(m_rng), U(m_rng), U(m_rng) } * (speed * 0.35);
            const double vl = length(a.vel);
            if (vl > speed * 1.5) a.vel = a.vel * (speed * 1.5 / vl);
            if (code == 2 && vl < speed * 0.5) a.vel = normalized(a.vel + Vec3{ U(m_rng), U(m_rng), U(m_rng) }) * speed;
            a.pos = a.pos + a.vel * dt;
            double* p[3] = { &a.pos.x, &a.pos.y, &a.pos.z };
            double* v[3] = { &a.vel.x, &a.vel.y, &a.vel.z };
            const double l[3] = { lim.x, lim.y, lim.z };
            for (int i = 0; i < 3; ++i) {
                if (*p[i] > l[i]) {
                    *p[i] = l[i];
                    *v[i] = -std::fabs(*v[i]);
                }
                if (*p[i] < -l[i]) {
                    *p[i] = -l[i];
                    *v[i] = std::fabs(*v[i]);
                }
            }
        }
    }
    if (!m_electrons.empty()) {
        const bool field = m_options["field"];
        const Vec3 lo = m_box_lo - Vec3{ 40.0, 40.0, 40.0 }, hi = m_box_hi + Vec3{ 40.0, 40.0, 40.0 };
        for (LatticeAtom& e : m_electrons) {
            if (U(m_rng) > 1.0 - 4.0 * dt) {
                const double th = (U(m_rng) + 1.0) * PI, ph = std::acos(U(m_rng));
                e.vel = Vec3{ std::sin(ph) * std::cos(th), std::cos(ph), std::sin(ph) * std::sin(th) } * 420.0;
            }
            if (field) e.vel.x += 900.0 * dt;
            e.pos = e.pos + e.vel * dt;
            if (e.pos.x > hi.x) e.pos.x = lo.x;
            if (e.pos.x < lo.x) e.pos.x = hi.x;
            if (e.pos.y > hi.y || e.pos.y < lo.y) e.vel.y = -e.vel.y;
            if (e.pos.z > hi.z || e.pos.z < lo.z) e.vel.z = -e.vel.z;
            e.pos.y = std::max(lo.y, std::min(hi.y, e.pos.y));
            e.pos.z = std::max(lo.z, std::min(hi.z, e.pos.z));
        }
    }
}

Color BondingView::species_color(int species) const {
    if (m_mode == BondingMode::Ionic) {
        const fc::IonDescriptor& ion = species == 0 ? m_ionic.cation : m_ionic.anion;
        if (ion.element) return element_color(ion.element->atomic_number);
        return species == 0 ? rgb(180, 140, 230) : rgb(120, 210, 170);
    }
    if (species >= 0 && species < static_cast<int>(m_metal.components.size()) && m_metal.components[static_cast<std::size_t>(species)].element) return element_color(m_metal.components[static_cast<std::size_t>(species)].element->atomic_number);
    return rgb(180, 180, 190);
}

double BondingView::species_radius(int species) const {
    double r = 100.0;
    if (m_mode == BondingMode::Ionic) {
        r = species == 0 ? m_ionic.cation_radius : m_ionic.anion_radius;
    } else if (species >= 0 && species < static_cast<int>(m_metal.components.size()) && m_metal.components[static_cast<std::size_t>(species)].element) {
        const fc::Element& e = *m_metal.components[static_cast<std::size_t>(species)].element;
        r = finite(e.metallic_radius) ? e.metallic_radius : finite(m_metal.metallic_radius) ? m_metal.metallic_radius : 130.0;
    }
    if (!finite(r) || r <= 0.0) r = 100.0;
    return r;
}

std::string BondingView::species_label(int species) const {
    if (m_mode == BondingMode::Ionic) {
        const fc::IonDescriptor& ion = species == 0 ? m_ionic.cation : m_ionic.anion;
        return ion.text(true);
    }
    if (species >= 0 && species < static_cast<int>(m_metal.components.size()) && m_metal.components[static_cast<std::size_t>(species)].element) return m_metal.components[static_cast<std::size_t>(species)].element->symbol;
    return "?";
}

int BondingView::pick_lattice(Vec2 screen) const {
    int best = -1;
    double best_depth = -1e300;
    const double scale = m_options.at("space_fill") ? 1.0 : 0.45;
    const double zmid = 0.0;
    for (std::size_t i = 0; i < m_lattice.size(); ++i) {
        const LatticeAtom& a = m_lattice[i];
        if (!m_3d && std::fabs(a.site.z - zmid) > 0.3 * m_nn_distance && !m_lattice.empty()) continue;
        double depth = 0.0, s = 1.0;
        const Vec2 p = m_cam.project(a.pos, &depth, &s);
        const double r = std::max(6.0, species_radius(a.species) * scale * s);
        if (length(p - screen) <= r && depth > best_depth) {
            best = static_cast<int>(i);
            best_depth = depth;
        }
    }
    return best;
}

void BondingView::draw_lattice(Canvas& canvas, const Theme& theme) {
    const bool ionic = m_mode == BondingMode::Ionic;
    const bool ok = ionic ? m_ionic.ok : m_metal.ok;
    if (!ok || m_lattice.empty()) {
        std::vector<std::string> lines{ ionic ? "Ionic lattice" : "Metallic lattice" };
        lines.push_back(m_lattice_error.empty() ? std::string(ionic ? "pick a metal and a nonmetal" : "pick one or more metals") : m_lattice_error);
        draw_hud(canvas, theme, lines, false);
        return;
    }
    if (!m_framed) {
        std::vector<Vec3> pts;
        for (const LatticeAtom& a : m_lattice) pts.push_back(a.site);
        if (!m_3d) m_cam.face_front();
        m_cam.frame(pts, 1.2 * std::max(species_radius(0), species_radius(1)));
        m_framed = true;
    }
    const double scale = m_options["space_fill"] ? 1.0 : 0.45;
    const bool solid = m_lattice_state == 0 && m_assembly > 0.97;
    const bool sticks = ionic ? m_options["neighbors"] : m_options["metal_neighbors"];
    std::vector<int> coordination(m_lattice.size(), 0);
    for (const auto& p : m_neighbors) {
        ++coordination[p.first];
        ++coordination[p.second];
    }
    if (m_3d) {
        Scene3D scene(m_cam);
        if (sticks && solid && !m_options["space_fill"]) {
            for (const auto& p : m_neighbors) {
                const LatticeAtom& a = m_lattice[p.first];
                const LatticeAtom& b = m_lattice[p.second];
                scene.cylinder(a.pos, b.pos, 9.0, mix(species_color(a.species), rgb(150, 150, 160), 0.5), mix(species_color(b.species), rgb(150, 150, 160), 0.5));
            }
        }
        for (std::size_t i = 0; i < m_lattice.size(); ++i) {
            const LatticeAtom& a = m_lattice[i];
            Color c = species_color(a.species);
            if (a.interstitial) c = mix(c, rgb(255, 255, 255), 0.15);
            scene.sphere(a.pos, species_radius(a.species) * scale * (a.interstitial ? 0.8 : 1.0), c);
        }
        if (m_sel_lattice >= 0 && m_sel_lattice < static_cast<int>(m_lattice.size())) {
            const LatticeAtom& s = m_lattice[static_cast<std::size_t>(m_sel_lattice)];
            double depth = 0.0, sc = 1.0;
            m_cam.project(s.pos, &depth, &sc);
            scene.ring(s.pos, species_radius(s.species) * scale * sc + 3.0, 2.5, theme.accent, 1e6);
            for (const auto& p : m_neighbors) {
                if (static_cast<int>(p.first) != m_sel_lattice && static_cast<int>(p.second) != m_sel_lattice) continue;
                const LatticeAtom& o = m_lattice[static_cast<int>(p.first) == m_sel_lattice ? p.second : p.first];
                scene.line(s.pos, o.pos, 2.5, theme.yellow, 1e6);
                scene.ring(o.pos, species_radius(o.species) * scale * sc + 2.0, 1.8, theme.yellow, 1e6);
            }
        }
        if (m_mode == BondingMode::Metallic && m_options["electron_sea"]) {
            for (const LatticeAtom& e : m_electrons) scene.disc(e.pos, 2.2, rgb(theme.yellow.r, theme.yellow.g, theme.yellow.b, 210), 30.0);
        }
        if (m_options["unit_cell"] && solid) {
            const fc::UnitCell& cell = ionic ? m_ionic.cell : m_metal.cell;
            if (!cell.sites.empty()) {
                const Vec3 origin = m_origin;
                const auto corner = [&](double fx, double fy, double fz) {
                    const std::array<double, 3> v = cell.to_cartesian(fx, fy, fz);
                    return origin + Vec3{ v[0], v[1], v[2] };
                };
                const double f[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
                const int e[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
                for (const auto& edge : e) {
                    scene.dashed(corner(f[edge[0]][0], f[edge[0]][1], f[edge[0]][2]), corner(f[edge[1]][0], f[edge[1]][1], f[edge[1]][2]), 1.4, theme.accent, 5e5);
                }
            }
        }
        scene.draw(canvas);
    } else {
        const double slab = 0.3 * m_nn_distance;
        double zlayer = 0.0, best = 1e300;
        for (const LatticeAtom& a : m_lattice) {
            if (std::fabs(a.site.z) < best) {
                best = std::fabs(a.site.z);
                zlayer = a.site.z;
            }
        }
        const auto in_layer = [&](const LatticeAtom& a) { return std::fabs(a.site.z - zlayer) <= slab; };
        if (sticks && solid) {
            for (const auto& p : m_neighbors) {
                const LatticeAtom& a = m_lattice[p.first];
                const LatticeAtom& b = m_lattice[p.second];
                if (!in_layer(a) || !in_layer(b)) continue;
                canvas.line_px(m_cam.project(a.pos), m_cam.project(b.pos), theme.grid, 2.0);
            }
        }
        for (std::size_t i = 0; i < m_lattice.size(); ++i) {
            const LatticeAtom& a = m_lattice[i];
            if (!in_layer(a) && solid) continue;
            if (!solid && std::fabs(a.pos.z - zlayer) > 2.0 * m_nn_distance) continue;
            double s = 1.0;
            const Vec2 p = m_cam.project(a.pos, nullptr, &s);
            const double r = species_radius(a.species) * scale * s;
            const Color c = species_color(a.species);
            canvas.disc_px(p, r, mix(c, theme.background, 0.25));
            canvas.circle_px(p, r, c, 1.5);
            if (r > 9.0) canvas.text(Vec2{ p.x, p.y + std::min(r * 0.4, 7.0) }, species_label(a.species), std::min(18.0, r * 0.8), theme.text, TextAlign::Center);
            if (static_cast<int>(i) == m_sel_lattice) canvas.circle_px(p, r + 3.0, theme.accent, 2.5);
        }
        if (m_mode == BondingMode::Metallic && m_options["electron_sea"]) {
            for (const LatticeAtom& e : m_electrons) {
                if (std::fabs(e.pos.z - zlayer) > m_nn_distance) continue;
                canvas.disc_px(m_cam.project(e.pos), 2.2, rgb(theme.yellow.r, theme.yellow.g, theme.yellow.b, 210));
            }
        }
    }
    std::vector<std::string> hud;
    if (ionic) {
        hud.push_back(pretty_formula(m_ionic.formula) + "  \xC2\xB7  " + m_ionic.name + "  \xC2\xB7  " + fc::to_string(m_ionic.structure));
    } else {
        hud.push_back(m_metal.name + "  \xC2\xB7  " + fc::to_string(m_metal.type) + "  \xC2\xB7  " + fc::to_string(m_metal.structure));
    }
    std::string st = lattice_state() + " at " + m_shared.units.format(m_temperature, "temperature", 5);
    if (m_assembly < 1.0 && m_lattice_state == 0) st += "  \xC2\xB7  assembling " + fixed(100.0 * m_assembly, 0) + "%";
    if (!m_3d) st += "  \xC2\xB7  one layer (z = 0)";
    hud.push_back(st);
    draw_hud(canvas, theme, hud, false);
    std::vector<std::string> bottom;
    if (m_sel_lattice >= 0 && m_sel_lattice < static_cast<int>(m_lattice.size())) {
        const LatticeAtom& s = m_lattice[static_cast<std::size_t>(m_sel_lattice)];
        bottom.push_back(species_label(s.species) + "  \xC2\xB7  " + std::to_string(coordination[static_cast<std::size_t>(m_sel_lattice)]) + " nearest neighbors");
    } else {
        bottom.push_back(m_3d ? "drag to rotate, wheel to zoom, Shift+drag to pan, click an atom" : "drag to pan, wheel to zoom, click an atom");
    }
    draw_hud(canvas, theme, bottom, true);
}

std::string BondingView::ionic_readouts() {
    std::ostringstream o;
    const UnitPrefs& u = m_shared.units;
    if (!m_ionic.ok) {
        o << "#Ionic compound\n!" << (m_lattice_error.empty() ? "pick a metal and a nonmetal" : m_lattice_error) << "\n";
        return o.str();
    }
    const fc::IonicCrystal& c = m_ionic;
    o << "#Compound\n";
    o << "Formula\t" << pretty_formula(c.formula) << "\n";
    o << "Name\t" << c.name << "\n";
    o << "Ions\t" << c.cation_count << " " << c.cation.text(true) << " (" << c.cation.name << ") + " << c.anion_count << " " << c.anion.text(true) << " (" << c.anion.name << ")\n";
    o << "Structure\t" << fc::to_string(c.structure) << "\n";
    if (!c.structure_reason.empty()) o << c.structure_reason << "\n";
    o << "#Geometry\n";
    o << "Cation radius\t" << u.format(c.cation_radius * 1e-12, "length", 4) << "\n";
    o << "Anion radius\t" << u.format(c.anion_radius * 1e-12, "length", 4) << "\n";
    o << "Radius ratio r\xE2\x82\x8A/r\xE2\x82\x8B\t" << fixed(c.radius_ratio, 3) << "\n";
    o << "Coordination\t" << c.cation_coordination << " : " << c.anion_coordination;
    if (c.cell.coordination_shape[0] && *c.cell.coordination_shape[0]) o << " (" << c.cell.coordination_shape[0] << ")";
    o << "\n";
    o << "Nearest distance\t" << u.format(c.nearest_distance * 1e-12, "length", 4) << "\n";
    o << "Cell edge a\t" << u.format(c.cell.a * 1e-12, "length", 4) << "\n";
    if (c.cell.hexagonal() || std::fabs(c.cell.c - c.cell.a) > 1e-6) o << "Cell edge c\t" << u.format(c.cell.c * 1e-12, "length", 4) << "\n";
    o << "Formula units per cell\t" << c.cell.formula_units << "\n";
    o << "#Lattice energy\n";
    o << "Lattice energy U\t" << u.format_kj(c.lattice_energy(), 4) << "\n";
    if (finite(c.lattice_energy_born_haber)) o << "Born\xE2\x80\x93Haber cycle\t" << u.format_kj(c.lattice_energy_born_haber, 4) << "\n";
    if (finite(c.lattice_energy_born_lande)) o << "Born\xE2\x80\x93Land\xC3\xA9\t" << u.format_kj(c.lattice_energy_born_lande, 4) << "\n";
    if (finite(c.lattice_energy_born_mayer)) o << "Born\xE2\x80\x93Mayer\t" << u.format_kj(c.lattice_energy_born_mayer, 4) << "\n";
    if (finite(c.lattice_energy_kapustinskii)) o << "Kapustinskii\t" << u.format_kj(c.lattice_energy_kapustinskii, 4) << "\n";
    if (finite(c.madelung_constant)) o << "Madelung constant\t" << fixed(c.madelung_constant, 4) << "\n";
    if (finite(c.born_exponent)) o << "Born exponent n\t" << fixed(c.born_exponent, 1) << "\n";
    o << "U is the energy released when the gas-phase ions come together into the solid, the step the animation shows.\n";
    if (!c.steps.empty()) {
        o << "#Born\xE2\x80\x93Haber cycle\n";
        for (const std::string& s : c.steps) o << s << "\n";
    }
    o << "#Properties\n";
    o << "Molar mass\t" << u.format(c.molar_mass, "molar mass", 6) << "\n";
    if (finite(c.density)) o << "Density\t" << u.format(c.density * 1000.0, "density", 4) << "\n";
    if (finite(c.electronegativity_difference)) o << "\xCE\x94\xCF\x87\t" << fixed(c.electronegativity_difference, 2) << " (" << fixed(c.percent_ionic, 0) << "% ionic)\n";
    if (finite(c.enthalpy_of_formation)) o << "\xCE\x94H\xC2\xB0" "f\t" << u.format_kj(c.enthalpy_of_formation, 4) << "\n";
    if (finite(c.melting_point)) o << "Melting point\t" << u.format(c.melting_point, "temperature", 5) << "\n";
    if (finite(c.boiling_point)) o << "Boiling point\t" << u.format(c.boiling_point, "temperature", 5) << "\n";
    o << "Solubility in water\t" << fc::to_string(c.solubility.solubility) << "\n";
    if (!c.solubility.rule.empty()) o << c.solubility.rule << "\n";
    o << "#At " << format_t() << "\n";
    o << "State\t" << lattice_state() << "\n";
    if (!c.notes.empty()) {
        o << "#Notes\n";
        for (const std::string& n : c.notes) o << n << "\n";
    }
    return o.str();
}

std::string BondingView::metallic_readouts() {
    std::ostringstream o;
    const UnitPrefs& u = m_shared.units;
    if (!m_metal.ok) {
        o << "#Metal\n!" << (m_lattice_error.empty() ? "pick one or more metals" : m_lattice_error) << "\n";
        return o.str();
    }
    const fc::MetallicCrystal& m = m_metal;
    o << "#" << (m.components.size() > 1 ? "Alloy" : "Metal") << "\n";
    o << "Name\t" << m.name << "\n";
    if (!m.other_names.empty()) {
        std::string others;
        for (std::size_t i = 0; i < m.other_names.size() && i < 3; ++i) others += (i ? ", " : "") + m.other_names[i];
        o << "Also called\t" << others << "\n";
    }
    if (!m.formula.empty()) o << "Formula\t" << pretty_formula(m.formula) << "\n";
    std::string comp;
    for (const fc::AlloyComponent& c : m.components) comp += (comp.empty() ? "" : ", ") + std::string(c.element->symbol) + " " + fixed(100.0 * c.fraction, 1) + "%";
    o << "Composition (atoms)\t" << comp << "\n";
    o << "Type\t" << fc::to_string(m.type) << "\n";
    o << "Structure\t" << fc::to_string(m.structure) << "\n";
    if (!m.structure_note.empty()) o << m.structure_note << "\n";
    if (!m.checks.empty()) {
        o << "#Hume-Rothery rules\n";
        for (const fc::HumeRotheryCheck& c : m.checks) o << c.rule << "\t" << (c.passes ? "passes" : "fails") << (c.detail.empty() ? "" : ": " + c.detail) << "\n";
        if (finite(m.size_mismatch)) o << "Size mismatch\t" << fixed(100.0 * m.size_mismatch, 1) << "%\n";
        if (finite(m.electronegativity_spread)) o << "\xCF\x87 spread\t" << fixed(m.electronegativity_spread, 2) << "\n";
    }
    o << "#Lattice\n";
    o << "Lattice constant a\t" << u.format(m.lattice_constant * 1e-12, "length", 4) << "\n";
    if (finite(m.lattice_c)) o << "Lattice constant c\t" << u.format(m.lattice_c * 1e-12, "length", 4) << "\n";
    o << "Metallic radius\t" << u.format(m.metallic_radius * 1e-12, "length", 4) << "\n";
    o << "Coordination\t" << m.coordination << "\n";
    if (finite(m.packing_fraction)) o << "Packing fraction\t" << fixed(100.0 * m.packing_fraction, 1) << "%\n";
    o << "Molar mass\t" << u.format(m.molar_mass, "molar mass", 5) << "\n";
    if (finite(m.density)) o << "Density\t" << u.format(m.density * 1000.0, "density", 4) << "\n";
    o << "#Electron sea\n";
    o << "Valence e\xE2\x81\xBB per atom\t" << fixed(m.valence_electrons_per_atom, 2) << "\n";
    if (finite(m.electron_density)) o << "Electron density\t" << format_number(m.electron_density, 4) << " m\xE2\x81\xBB\xC2\xB3\n";
    if (finite(m.fermi_energy)) o << "Fermi energy\t" << fixed(m.fermi_energy, 2) << " eV\n";
    if (finite(m.fermi_temperature)) o << "Fermi temperature\t" << u.format(m.fermi_temperature, "temperature", 4) << "\n";
    if (finite(m.fermi_velocity)) o << "Fermi velocity\t" << format_number(m.fermi_velocity, 4) << " m/s\n";
    o << "Metal cations sit in a sea of shared valence electrons; turn on Electric field to see them drift (a current).\n";
    o << "#Thermal\n";
    if (finite(m.melting_point)) o << "Melting point\t" << u.format(m.melting_point, "temperature", 5) << "\n";
    if (finite(m.boiling_point)) o << "Boiling point\t" << u.format(m.boiling_point, "temperature", 5) << "\n";
    if (finite(m.cohesive_energy)) o << "Cohesive energy\t" << u.format_kj(m.cohesive_energy, 4) << "\n";
    o << "State at " << u.format(m_temperature, "temperature", 5) << "\t" << lattice_state() << "\n";
    if (!m.notes.empty()) {
        o << "#Notes\n";
        for (const std::string& n : m.notes) o << n << "\n";
    }
    return o.str();
}

std::string BondingView::lattice_selection() {
    std::ostringstream o;
    if (m_sel_lattice < 0 || m_sel_lattice >= static_cast<int>(m_lattice.size())) {
        o << "#Selection\nClick an atom or ion in the lattice to see its neighbors.\n";
        return o.str();
    }
    const LatticeAtom& a = m_lattice[static_cast<std::size_t>(m_sel_lattice)];
    int cn = 0;
    double d = 0.0;
    for (const auto& p : m_neighbors) {
        if (static_cast<int>(p.first) != m_sel_lattice && static_cast<int>(p.second) != m_sel_lattice) continue;
        ++cn;
        d = length(m_lattice[p.first].site - m_lattice[p.second].site);
    }
    o << "#" << species_label(a.species) << "\n";
    if (m_mode == BondingMode::Ionic) {
        const fc::IonDescriptor& ion = a.species == 0 ? m_ionic.cation : m_ionic.anion;
        o << "Ion\t" << ion.name << "\n";
        o << "Charge\t" << (ion.charge > 0 ? "+" : "") << ion.charge << "\n";
        const char* shape = m_ionic.cell.coordination_shape[a.species == 0 ? 0 : 1];
        if (shape && *shape) o << "Surrounded by\t" << shape << "\n";
    } else if (a.species >= 0 && a.species < static_cast<int>(m_metal.components.size()) && m_metal.components[static_cast<std::size_t>(a.species)].element) {
        const fc::Element& e = *m_metal.components[static_cast<std::size_t>(a.species)].element;
        o << "Element\t" << e.name << "\n";
        if (a.interstitial) o << "Site\tinterstitial (in a hole between the host atoms)\n";
    }
    o << "Radius\t" << m_shared.units.format(species_radius(a.species) * 1e-12, "length", 4) << "\n";
    o << "Neighbors here\t" << cn << (cn > 0 ? "" : " (on the edge of the block)") << "\n";
    if (d > 0.0) o << "Neighbor distance\t" << m_shared.units.format(d * 1e-12, "length", 4) << "\n";
    o << "Atoms on the outside of the block have fewer neighbors than ones inside.\n";
    return o.str();
}

} // namespace sapling

#include "bonding_view.hpp"


#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

namespace sapling {

namespace {

std::string arg(const std::vector<std::string>& a, std::size_t i) {
    return i < a.size() ? a[i] : std::string();
}

double num(const std::string& s, double fallback) {
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    return end == s.c_str() || !std::isfinite(v) ? fallback : v;
}

const char* tool_name(SandboxTool t) {
    switch (t) {
        case SandboxTool::Atom: return "atom";
        case SandboxTool::Bond: return "bond";
        case SandboxTool::Move: return "move";
        case SandboxTool::Erase: return "erase";
        default: return "inspect";
    }
}

const char* mode_name(BondingMode m) {
    switch (m) {
        case BondingMode::Covalent: return "covalent";
        case BondingMode::Ionic: return "ionic";
        default: return "metallic";
    }
}

} // namespace

BondingView::BondingView(ChemShared& shared) : m_shared(shared) {
    m_options["lone_pairs"] = true;
    m_options["formal_charges"] = true;
    m_options["partial_charges"] = false;
    m_options["bond_dipoles"] = false;
    m_options["molecular_dipole"] = false;
    m_options["labels"] = true;
    m_options["space_fill"] = false;
    m_options["unit_cell"] = true;
    m_options["electron_sea"] = true;
    m_options["neighbors"] = true;
    m_options["metal_neighbors"] = false;
    m_options["spin"] = false;
    m_options["thermal"] = true;
    m_options["field"] = false;
    m_options["snap"] = true;
    m_cam.face_front();
    load_smiles("O", 0, "water");
    m_undo.clear();
    m_components = { { 29, 0.7 }, { 30, 0.3 } };
}

double BondingView::current_scale() const {
    return m_options.count("space_fill") && m_options.at("space_fill") ? 2.2 : 1.0;
}

void BondingView::changed(bool structure) {
    ++m_version;
    if (structure) m_analysis_dirty = true;
}

void BondingView::push_undo() {
    m_undo.push_back(m_cov);
    if (m_undo.size() > 120) m_undo.erase(m_undo.begin());
    m_redo.clear();
}

void BondingView::frame_view() {
    m_framed = false;
}

std::vector<DrawAtom> BondingView::draw_atoms() const {
    std::vector<DrawAtom> out;
    out.reserve(m_cov.atoms.size());
    for (std::size_t i = 0; i < m_cov.atoms.size(); ++i) {
        DrawAtom a;
        a.p = m_cov.atoms[i].p;
        a.z = m_cov.atoms[i].z;
        if (m_analysis.ok && i < m_analysis.lone_pairs.size()) {
            a.lone_pairs = m_analysis.lone_pairs[i];
            a.unpaired = m_analysis.unpaired[i];
            a.formal_charge = m_analysis.formal[i];
            a.partial = i < m_analysis.partial.size() ? m_analysis.partial[i] : 0.0;
            a.problem = m_analysis.problem[i];
        }
        a.selected = static_cast<int>(i) == m_sel_atom;
        a.highlight = static_cast<int>(i) == m_hover_atom && m_hover_atom != m_sel_atom;
        out.push_back(a);
    }
    return out;
}

std::vector<DrawBond> BondingView::draw_bonds() const {
    std::vector<DrawBond> out;
    out.reserve(m_cov.bonds.size());
    for (std::size_t i = 0; i < m_cov.bonds.size(); ++i) {
        DrawBond b;
        b.a = m_cov.bonds[i].a;
        b.b = m_cov.bonds[i].b;
        b.order = m_cov.bonds[i].order;
        b.selected = static_cast<int>(i) == m_sel_bond;
        out.push_back(b);
    }
    return out;
}

void BondingView::draw(Canvas& canvas, const Theme& theme) {
    m_cam.width = canvas.width();
    m_cam.height = canvas.height();
    canvas.clear(theme.background);
    if (m_mode != BondingMode::Covalent) {
        draw_lattice(canvas, theme);
        return;
    }
    ensure_analysis();
    if (!m_framed) {
        std::vector<Vec3> pts;
        for (const SandboxAtom& a : m_cov.atoms) pts.push_back(a.p);
        if (pts.empty()) {
            m_cam.target = Vec3{};
            m_cam.zoom = std::min(canvas.width(), canvas.height()) / 900.0;
            m_cam.pan = Vec2{};
        } else {
            m_cam.frame(pts, 160.0);
            m_cam.zoom = std::min(m_cam.zoom, std::min(canvas.width(), canvas.height()) / 300.0);
        }
        m_framed = true;
    }
    if (flat()) {
        const double step = 50.0 * m_cam.zoom;
        if (step > 12.0) {
            const Vec2 o = m_cam.project(Vec3{});
            const Color dotc = mix(theme.background, theme.grid, 0.9);
            for (double x = std::fmod(o.x, step); x < canvas.width(); x += step) {
                for (double y = std::fmod(o.y, step); y < canvas.height(); y += step) canvas.put_pixel(static_cast<int>(x), static_cast<int>(y), dotc);
            }
        }
    }
    std::vector<DrawAtom> atoms = draw_atoms();
    std::vector<DrawBond> bonds = draw_bonds();
    if (m_down && m_dragged && m_drag_atom >= 0 && (m_tool == SandboxTool::Atom || m_tool == SandboxTool::Bond)) {
        const int over = pick_atom(m_cam, atoms, m_mouse, flat(), current_scale());
        DrawBond b;
        b.a = m_drag_atom;
        b.preview = true;
        if (over >= 0 && over != m_drag_atom) {
            b.b = over;
        } else {
            DrawAtom ghost;
            ghost.z = m_element;
            const Vec3 from = m_cov.atoms[static_cast<std::size_t>(m_drag_atom)].p;
            Vec3 to = place_point(m_mouse, m_drag_depth);
            if (m_options["snap"] && !(m_mods & k_shift)) {
                const Vec3 d = to - from;
                if (length(d) > 1e-6) to = from + normalized(d) * ideal_length(m_cov.atoms[static_cast<std::size_t>(m_drag_atom)].z, m_element, 1);
            }
            ghost.p = to;
            ghost.highlight = true;
            atoms.push_back(ghost);
            b.b = static_cast<int>(atoms.size()) - 1;
        }
        bonds.push_back(b);
    }
    DrawOptions opt;
    opt.flat = flat();
    opt.lone_pairs = m_options["lone_pairs"];
    opt.formal_charges = m_options["formal_charges"];
    opt.partial_charges = m_options["partial_charges"];
    opt.bond_dipoles = m_options["bond_dipoles"];
    opt.molecular_dipole = m_options["molecular_dipole"] && m_analysis.ok;
    opt.labels = m_options["labels"];
    opt.atom_scale = current_scale();
    opt.dipole = m_analysis.dipole_vector;
    opt.dipole_debye = m_analysis.mu_shape;
    draw_molecule(canvas, m_cam, theme, atoms, bonds, opt);

    std::vector<std::string> hud;
    const fc::Element* e = fc::find_element(m_element);
    std::string head = std::string(m_3d ? "3D" : "2D") + " covalent  \xC2\xB7  tool: " + tool_name(m_tool) + "  \xC2\xB7  element: " + (e ? e->symbol : "?");
    if (m_cov.charge != 0) head += "  \xC2\xB7  charge " + std::to_string(m_cov.charge);
    hud.push_back(head);
    if (m_cov.atoms.empty()) hud.push_back("click to place an atom; drag from an atom to grow a bond");
    else if (!m_analysis.verdict.empty()) hud.push_back(m_analysis.verdict);
    draw_hud(canvas, theme, hud, false);
    std::vector<std::string> bottom;
    if (m_hover_atom >= 0 && m_hover_atom < static_cast<int>(m_cov.atoms.size())) {
        const fc::Element* he = fc::find_element(m_cov.atoms[static_cast<std::size_t>(m_hover_atom)].z);
        std::string line = std::string(he ? he->name : "?") + " (" + (he ? he->symbol : "?") + ")";
        if (m_analysis.ok && static_cast<std::size_t>(m_hover_atom) < m_analysis.lone_pairs.size()) {
            line += "  \xC2\xB7  lone pairs " + std::to_string(m_analysis.lone_pairs[static_cast<std::size_t>(m_hover_atom)]);
            line += "  \xC2\xB7  formal charge " + std::to_string(m_analysis.formal[static_cast<std::size_t>(m_hover_atom)]);
        }
        bottom.push_back(line);
    } else if (m_hover_bond >= 0 && m_hover_bond < static_cast<int>(m_analysis.bonds.size())) {
        const SandboxBond& b = m_cov.bonds[static_cast<std::size_t>(m_hover_bond)];
        const BondReport& r = m_analysis.bonds[static_cast<std::size_t>(m_hover_bond)];
        const fc::Element* ea = fc::find_element(m_cov.atoms[static_cast<std::size_t>(b.a)].z);
        const fc::Element* eb = fc::find_element(m_cov.atoms[static_cast<std::size_t>(b.b)].z);
        char buf[64];
        std::snprintf(buf, sizeof buf, "%.2f", r.avg_order);
        bottom.push_back(std::string(ea ? ea->symbol : "?") + "\xE2\x80\x93" + (eb ? eb->symbol : "?") + "  order " + buf + "  \xC2\xB7  " + m_shared.units.format(r.table_length * 1e-12, "length", 4) + "  \xC2\xB7  " + m_shared.units.format_kj(r.energy, 4) + " (click for more)");
    } else if (!m_hint.empty()) {
        bottom.push_back(m_hint);
    }
    draw_hud(canvas, theme, bottom, true);
}

Vec3 BondingView::place_point(Vec2 screen, double depth) const {
    Vec3 p = m_cam.unproject(screen, depth);
    if (flat()) p.z = 0.0;
    return p;
}

bool BondingView::mouse_down(Vec2 p, int button, unsigned mods) {
    m_down = true;
    m_button = button;
    m_mods = mods;
    m_down_at = p;
    m_last = p;
    m_mouse = p;
    m_dragged = false;
    m_drag_atom = -1;
    m_drag_bond = -1;
    m_drag_pushed = false;
    if (m_mode != BondingMode::Covalent) return false;
    if (button != 1) return false;
    const std::vector<DrawAtom> atoms = draw_atoms();
    m_drag_atom = pick_atom(m_cam, atoms, p, flat(), current_scale());
    if (m_drag_atom < 0) m_drag_bond = pick_bond(m_cam, atoms, draw_bonds(), p);
    if (m_drag_atom >= 0) {
        double depth = 0.0;
        m_cam.project(m_cov.atoms[static_cast<std::size_t>(m_drag_atom)].p, &depth);
        m_drag_depth = depth;
    } else {
        m_drag_depth = 0.0;
    }
    return true;
}

bool BondingView::mouse_move(Vec2 p, bool) {
    const Vec2 delta = p - m_last;
    m_last = p;
    m_mouse = p;
    if (!m_down) {
        if (m_mode != BondingMode::Covalent) {
            return false;
        }
        const std::vector<DrawAtom> atoms = draw_atoms();
        const int ha = pick_atom(m_cam, atoms, p, flat(), current_scale());
        const int hb = ha < 0 ? pick_bond(m_cam, atoms, draw_bonds(), p) : -1;
        if (ha != m_hover_atom || hb != m_hover_bond) {
            m_hover_atom = ha;
            m_hover_bond = hb;
            return true;
        }
        return false;
    }
    if (!m_dragged && length(p - m_down_at) > 4.0) m_dragged = true;
    if (!m_dragged) return false;
    const bool atom_action = m_mode == BondingMode::Covalent && m_button == 1 && m_drag_atom >= 0 && (m_tool == SandboxTool::Atom || m_tool == SandboxTool::Bond || m_tool == SandboxTool::Move);
    if (!atom_action) {
        const bool pan = m_button == 2 || (m_mods & k_shift) || !m_3d;
        if (pan) m_cam.pan = m_cam.pan + delta;
        else m_cam.orbit(delta.x, delta.y);
        return true;
    }
    if (m_tool == SandboxTool::Move && m_drag_atom >= 0) {
        if (!m_drag_pushed) {
            push_undo();
            m_drag_pushed = true;
        }
        m_cov.atoms[static_cast<std::size_t>(m_drag_atom)].p = place_point(p, m_drag_depth);
        changed(true);
        return true;
    }
    return true;
}

bool BondingView::mouse_up(Vec2 p, int) {
    if (!m_down) return false;
    m_down = false;
    m_mouse = p;
    if (m_mode != BondingMode::Covalent) {
        if (!m_dragged && m_button == 1) {
            m_sel_lattice = pick_lattice(p);
            changed(false);
        }
        return true;
    }
    if (m_button != 1) return true;
    const std::vector<DrawAtom> atoms = draw_atoms();
    const bool drew = m_dragged && m_drag_atom >= 0;
    if (m_tool == SandboxTool::Atom || m_tool == SandboxTool::Bond) {
        if (drew) {
            const int over = pick_atom(m_cam, atoms, p, flat(), current_scale());
            push_undo();
            if (over >= 0 && over != m_drag_atom) {
                connect(m_drag_atom, over, true);
            } else if (over < 0) {
                const Vec3 from = m_cov.atoms[static_cast<std::size_t>(m_drag_atom)].p;
                Vec3 to = place_point(p, m_drag_depth);
                if (m_options["snap"] && !(m_mods & k_shift)) {
                    const Vec3 d = to - from;
                    if (length(d) > 1e-6) to = from + normalized(d) * ideal_length(m_cov.atoms[static_cast<std::size_t>(m_drag_atom)].z, m_element, 1);
                }
                const int n = add_atom(m_element, to);
                connect(m_drag_atom, n, false);
            }
            changed(true);
            return true;
        }
        if (m_dragged) return true;
        if (m_drag_atom >= 0) {
            if (m_tool == SandboxTool::Atom && m_cov.atoms[static_cast<std::size_t>(m_drag_atom)].z != m_element) {
                push_undo();
                m_cov.atoms[static_cast<std::size_t>(m_drag_atom)].z = m_element;
                changed(true);
            } else if (m_tool == SandboxTool::Bond) {
                if (m_sel_atom >= 0 && m_sel_atom != m_drag_atom) {
                    push_undo();
                    connect(m_sel_atom, m_drag_atom, true);
                    m_sel_atom = -1;
                    changed(true);
                } else {
                    m_sel_atom = m_drag_atom;
                    m_sel_bond = -1;
                    changed(false);
                }
            } else {
                m_sel_atom = m_drag_atom;
                m_sel_bond = -1;
                changed(false);
            }
            return true;
        }
        if (m_drag_bond >= 0) {
            if (m_tool == SandboxTool::Bond) {
                push_undo();
                SandboxBond& b = m_cov.bonds[static_cast<std::size_t>(m_drag_bond)];
                b.order = b.order >= 3 ? 1 : b.order + 1;
                changed(true);
            }
            m_sel_bond = m_drag_bond;
            m_sel_atom = -1;
            changed(false);
            return true;
        }
        if (m_tool == SandboxTool::Atom) {
            push_undo();
            m_sel_atom = add_atom(m_element, place_point(p, 0.0));
            m_sel_bond = -1;
            changed(true);
        } else {
            m_sel_atom = -1;
            m_sel_bond = -1;
            changed(false);
        }
        return true;
    }
    if (m_dragged) return true;
    if (m_tool == SandboxTool::Erase) {
        if (m_drag_atom >= 0) {
            push_undo();
            remove_atom(m_drag_atom);
            changed(true);
        } else if (m_drag_bond >= 0) {
            push_undo();
            m_cov.bonds.erase(m_cov.bonds.begin() + m_drag_bond);
            m_sel_bond = -1;
            changed(true);
        }
        m_hover_atom = -1;
        m_hover_bond = -1;
        return true;
    }
    m_sel_atom = m_drag_atom;
    m_sel_bond = m_drag_atom < 0 ? m_drag_bond : -1;
    changed(false);
    return true;
}

bool BondingView::wheel(Vec2 p, double dy) {
    if (dy == 0.0) return false;
    const double f = std::pow(1.15, dy > 0 ? 1.0 : -1.0);
    const Vec2 center{ 0.5 * m_cam.width, 0.5 * m_cam.height };
    const Vec2 c = center + m_cam.pan;
    const Vec2 next = p - (p - c) * f;
    m_cam.pan = next - center;
    m_cam.zoom_by(f);
    return true;
}

bool BondingView::double_click(Vec2) {
    frame_view();
    return true;
}

bool BondingView::key(const std::string& k, unsigned mods) {
    if (m_mode != BondingMode::Covalent) return false;
    if ((mods & k_control) && (k == "Z" || k == "z")) {
        if (mods & k_shift) command({ "redo" });
        else command({ "undo" });
        return true;
    }
    if ((mods & k_control) && (k == "Y" || k == "y")) {
        command({ "redo" });
        return true;
    }
    if (k == "Delete" || k == "Backspace") {
        if (m_sel_atom >= 0) {
            push_undo();
            remove_atom(m_sel_atom);
            changed(true);
            return true;
        }
        if (m_sel_bond >= 0 && m_sel_bond < static_cast<int>(m_cov.bonds.size())) {
            push_undo();
            m_cov.bonds.erase(m_cov.bonds.begin() + m_sel_bond);
            m_sel_bond = -1;
            changed(true);
            return true;
        }
        return false;
    }
    if ((k == "1" || k == "2" || k == "3") && m_sel_bond >= 0 && m_sel_bond < static_cast<int>(m_cov.bonds.size())) {
        push_undo();
        m_cov.bonds[static_cast<std::size_t>(m_sel_bond)].order = k[0] - '0';
        changed(true);
        return true;
    }
    if (k.size() == 1 && std::isalpha(static_cast<unsigned char>(k[0])) && !(mods & k_control)) {
        std::string sym(1, static_cast<char>(std::toupper(static_cast<unsigned char>(k[0]))));
        const fc::Element* e = fc::find_element(sym);
        if (e) {
            m_element = e->atomic_number;
            changed(false);
            return true;
        }
    }
    if (k == "Escape") {
        m_sel_atom = -1;
        m_sel_bond = -1;
        changed(false);
        return true;
    }
    return false;
}

bool BondingView::step(double dt) {
    bool redraw = false;
    if (m_options["spin"] && (m_3d || m_mode != BondingMode::Covalent) && !m_down) {
        m_cam.yaw += dt * 0.5;
        redraw = true;
    }
    if (m_mode != BondingMode::Covalent) {
        step_lattice(dt);
        redraw = true;
    }
    return redraw;
}

int BondingView::add_atom(int z, Vec3 p) {
    m_cov.atoms.push_back(SandboxAtom{ z, p });
    return static_cast<int>(m_cov.atoms.size()) - 1;
}

void BondingView::remove_atom(int index) {
    if (index < 0 || index >= static_cast<int>(m_cov.atoms.size())) return;
    std::vector<SandboxBond> keep;
    for (const SandboxBond& b : m_cov.bonds) {
        if (b.a == index || b.b == index) continue;
        SandboxBond c = b;
        if (c.a > index) --c.a;
        if (c.b > index) --c.b;
        keep.push_back(c);
    }
    m_cov.bonds = std::move(keep);
    m_cov.atoms.erase(m_cov.atoms.begin() + index);
    m_sel_atom = -1;
    m_sel_bond = -1;
    m_hover_atom = -1;
    m_hover_bond = -1;
}

int BondingView::find_bond(int a, int b) const {
    for (std::size_t i = 0; i < m_cov.bonds.size(); ++i) {
        const SandboxBond& x = m_cov.bonds[i];
        if ((x.a == a && x.b == b) || (x.a == b && x.b == a)) return static_cast<int>(i);
    }
    return -1;
}

void BondingView::connect(int a, int b, bool cycle) {
    if (a == b || a < 0 || b < 0) return;
    const int existing = find_bond(a, b);
    if (existing >= 0) {
        if (cycle) {
            SandboxBond& x = m_cov.bonds[static_cast<std::size_t>(existing)];
            x.order = x.order >= 3 ? 1 : x.order + 1;
        }
        m_sel_bond = existing;
        return;
    }
    m_cov.bonds.push_back(SandboxBond{ a, b, 1 });
    m_sel_bond = static_cast<int>(m_cov.bonds.size()) - 1;
    m_sel_atom = -1;
}

double BondingView::ideal_length(int za, int zb, int order) const {
    const fc::Element* a = fc::find_element(za);
    const fc::Element* b = fc::find_element(zb);
    if (!a || !b) return 150.0;
    try {
        const double l = fc::bond_length(*a, *b, static_cast<fc::BondOrder>(std::max(1, std::min(3, order))));
        if (l == l && l > 0.0) return l;
    } catch (...) {
    }
    return 150.0;
}

std::string BondingView::format_t() const {
    return m_shared.units.format(m_temperature, "temperature", 5) + ", " + m_shared.units.format(m_pressure, "pressure", 5);
}

std::string BondingView::kind_name() const {
    return mode_name(m_mode);
}

std::string BondingView::state_text() const {
    std::ostringstream o;
    const fc::Element* e = fc::find_element(m_element);
    o << "mode\t" << mode_name(m_mode) << "\n";
    o << "view\t" << (m_3d ? "3d" : "2d") << "\n";
    o << "tool\t" << tool_name(m_tool) << "\n";
    o << "element\t" << (e ? e->symbol : "C") << "\n";
    o << "charge\t" << m_cov.charge << "\n";
    o << "cation\t" << m_cation << "\n";
    o << "anion\t" << m_anion << "\n";
    o << "cells\t" << m_cells << "\n";
    std::string comp;
    for (const auto& c : m_components) {
        const fc::Element* ce = fc::find_element(c.first);
        if (!ce) continue;
        if (!comp.empty()) comp += ";";
        comp += std::string(ce->symbol) + "=" + format_number(c.second * 100.0, 4);
    }
    o << "metal\t" << comp << "\n";
    o << "name\t" << m_cov.name << "\n";
    for (const auto& kv : m_options) o << "option." << kv.first << "\t" << (kv.second ? 1 : 0) << "\n";
    return o.str();
}

std::string BondingView::to_json() const {
    std::ostringstream o;
    o << "{\n  \"type\": " << json_quote(kind_name()) << ",\n  \"sapling\": \"chemistry\",\n";
    if (m_mode == BondingMode::Covalent) {
        o << "  \"name\": " << json_quote(m_cov.name) << ",\n";
        std::string formula;
        if (m_analysis.ok) formula = m_analysis.molecule.formula().hill();
        o << "  \"formula\": " << json_quote(formula) << ",\n";
        o << "  \"charge\": " << m_cov.charge << ",\n  \"atoms\": [";
        for (std::size_t i = 0; i < m_cov.atoms.size(); ++i) {
            const SandboxAtom& a = m_cov.atoms[i];
            const fc::Element* e = fc::find_element(a.z);
            o << (i ? ",\n" : "\n") << "    {\"element\": " << json_quote(e ? e->symbol : "C") << ", \"x\": " << format_number(a.p.x, 8) << ", \"y\": " << format_number(a.p.y, 8) << ", \"z\": " << format_number(a.p.z, 8) << (a.neutrons ? ", \"neutrons\": " + std::to_string(a.neutrons) : std::string()) << "}";
        }
        o << "\n  ],\n  \"bonds\": [";
        for (std::size_t i = 0; i < m_cov.bonds.size(); ++i) {
            const SandboxBond& b = m_cov.bonds[i];
            o << (i ? ",\n" : "\n") << "    {\"a\": " << b.a << ", \"b\": " << b.b << ", \"order\": " << b.order << "}";
        }
        o << "\n  ]\n}";
    } else if (m_mode == BondingMode::Ionic) {
        o << "  \"cation\": " << json_quote(m_cation) << ",\n  \"anion\": " << json_quote(m_anion) << ",\n";
        o << "  \"formula\": " << json_quote(m_ionic.ok ? m_ionic.formula : std::string()) << ",\n";
        o << "  \"name\": " << json_quote(m_ionic.ok ? m_ionic.name : std::string()) << "\n}";
    } else {
        o << "  \"name\": " << json_quote(m_metal.ok ? m_metal.name : std::string()) << ",\n  \"components\": [";
        for (std::size_t i = 0; i < m_components.size(); ++i) {
            const fc::Element* e = fc::find_element(m_components[i].first);
            o << (i ? ", " : "") << "{\"element\": " << json_quote(e ? e->symbol : "Fe") << ", \"fraction\": " << format_number(m_components[i].second, 8) << "}";
        }
        o << "]\n}";
    }
    return o.str();
}

bool BondingView::from_json(const std::string& json, std::string* error) {
    JsonValue v;
    if (!parse_json(json, v, error)) return false;
    const std::string type = v.string_or("type", "covalent");
    if (type == "covalent") {
        CovalentState s;
        s.name = v.string_or("name", "");
        s.charge = static_cast<int>(v.number_or("charge", 0.0));
        const JsonValue* atoms = v.find("atoms");
        const JsonValue* bonds = v.find("bonds");
        if (atoms && atoms->type == JsonValue::Array) {
            for (const JsonValue& a : atoms->array) {
                const fc::Element* e = fc::find_element(a.string_or("element", "C"));
                if (!e) continue;
                s.atoms.push_back(SandboxAtom{ e->atomic_number, Vec3{ a.number_or("x", 0.0), a.number_or("y", 0.0), a.number_or("z", 0.0) }, static_cast<int>(a.number_or("neutrons", 0.0)) });
            }
        }
        if (bonds && bonds->type == JsonValue::Array) {
            for (const JsonValue& b : bonds->array) {
                SandboxBond sb{ static_cast<int>(b.number_or("a", -1)), static_cast<int>(b.number_or("b", -1)), static_cast<int>(b.number_or("order", 1)) };
                if (sb.a < 0 || sb.b < 0 || sb.a >= static_cast<int>(s.atoms.size()) || sb.b >= static_cast<int>(s.atoms.size()) || sb.a == sb.b) continue;
                sb.order = std::max(1, std::min(3, sb.order));
                s.bonds.push_back(sb);
            }
        }
        push_undo();
        m_cov = std::move(s);
        m_mode = BondingMode::Covalent;
        bool depth = false;
        for (const SandboxAtom& a : m_cov.atoms) depth = depth || std::fabs(a.p.z) > 1e-6;
        m_3d = depth;
        if (!m_3d) m_cam.face_front();
        m_sel_atom = m_sel_bond = -1;
        frame_view();
        changed(true);
        return true;
    }
    if (type == "ionic") {
        m_mode = BondingMode::Ionic;
        set_ionic(v.string_or("cation", "Na+"), v.string_or("anion", "Cl-"));
        return true;
    }
    if (type == "metallic") {
        std::string spec;
        const JsonValue* comps = v.find("components");
        if (comps && comps->type == JsonValue::Array) {
            for (const JsonValue& c : comps->array) {
                if (!spec.empty()) spec += ";";
                spec += c.string_or("element", "Fe") + "=" + format_number(c.number_or("fraction", 1.0), 8);
            }
        }
        m_mode = BondingMode::Metallic;
        set_metal(spec);
        return true;
    }
    if (error) *error = "unknown type '" + type + "'";
    return false;
}

namespace {

std::string central_atom_smiles(const std::string& text) {
    try {
        const fc::Formula f = fc::Formula::parse(text);
        const auto parts = f.composition();
        if (parts.size() < 2) return std::string();
        const fc::Element* central = nullptr;
        int ligands = 0;
        for (const auto& c : parts) {
            const int n = static_cast<int>(std::llround(c.second));
            if (std::fabs(c.second - n) > 1e-9 || n <= 0) return std::string();
            if (n == 1 && c.first->atomic_number != 1) {
                const double en = c.first->electronegativity();
                if (!central || (en == en && (central->electronegativity() != central->electronegativity() || en < central->electronegativity()))) central = c.first;
            }
        }
        if (!central) return std::string();
        for (const auto& c : parts) {
            if (c.first != central) ligands += static_cast<int>(std::llround(c.second));
        }
        if (ligands < 1 || ligands > 8) return std::string();
        std::vector<std::string> subs;
        for (const auto& c : parts) {
            if (c.first == central) continue;
            const bool chalcogen = c.first->atomic_number == 8 || c.first->atomic_number == 16;
            for (int i = 0; i < static_cast<int>(std::llround(c.second)); ++i) subs.push_back(std::string(chalcogen ? "=" : "") + "[" + c.first->symbol + "]");
        }
        std::string out = std::string("[") + central->symbol + "]";
        for (std::size_t i = 0; i + 1 < subs.size(); ++i) out += "(" + subs[i] + ")";
        out += subs.back();
        fc::Molecule m = fc::Molecule::from_smiles(out);
        return m.atom_count() > 0 ? out : std::string();
    } catch (...) {
        return std::string();
    }
}

} // namespace

std::string BondingView::command(const std::vector<std::string>& a) {
    const std::string c = arg(a, 0);
    if (c == "version") return std::to_string(m_version);
    if (c == "state") return state_text();
    if (c == "mode") {
        const std::string m = arg(a, 1);
        BondingMode next = m == "ionic" ? BondingMode::Ionic : m == "metallic" ? BondingMode::Metallic : BondingMode::Covalent;
        if (next != m_mode) {
            m_mode = next;
            m_sel_lattice = -1;
            if (m_mode != BondingMode::Covalent) {
                if (m_mode == BondingMode::Ionic && !m_ionic.ok) set_ionic(m_cation, m_anion);
                if (m_mode == BondingMode::Metallic && !m_metal.ok) {
                    std::string spec;
                    for (const auto& cp : m_components) {
                        const fc::Element* e = fc::find_element(cp.first);
                        if (!e) continue;
                        if (!spec.empty()) spec += ";";
                        spec += std::string(e->symbol) + "=" + format_number(cp.second, 6);
                    }
                    set_metal(spec);
                }
                rebuild_lattice();
                restart_assembly();
            }
            frame_view();
            if (!m_3d && m_mode == BondingMode::Covalent) m_cam.face_front();
            changed(true);
        }
        return mode_name(m_mode);
    }
    if (c == "view") {
        const bool want = arg(a, 1) == "3d";
        if (want != m_3d) {
            if (m_mode == BondingMode::Covalent) {
                if (!want) {
                    push_undo();
                    for (SandboxAtom& at : m_cov.atoms) {
                        const Vec3 v = m_cam.to_view(at.p);
                        at.p = Vec3{ v.x + m_cam.target.x, v.y + m_cam.target.y, 0.0 };
                    }
                    m_cam.face_front();
                    m_cam.perspective = 0.0;
                } else {
                    bool depth = false;
                    for (const SandboxAtom& at : m_cov.atoms) depth = depth || std::fabs(at.p.z) > 1e-6;
                    m_3d = true;
                    if (!depth && m_cov.atoms.size() > 2) tidy();
                    m_cam.yaw = 0.5;
                    m_cam.pitch = 0.3;
                    m_cam.perspective = 0.0008;
                }
            } else {
                if (!want) {
                    m_cam.face_front();
                    m_cam.perspective = 0.0;
                } else {
                    m_cam.yaw = 0.6;
                    m_cam.pitch = 0.45;
                    m_cam.perspective = 0.0005;
                }
            }
            m_3d = want;
            frame_view();
            changed(true);
        }
        return m_3d ? "3d" : "2d";
    }
    if (c == "tool") {
        const std::string t = arg(a, 1);
        m_tool = t == "bond" ? SandboxTool::Bond : t == "move" ? SandboxTool::Move : t == "erase" ? SandboxTool::Erase : t == "inspect" ? SandboxTool::Inspect : SandboxTool::Atom;
        static const char* hints[] = {
            "Atom: click empty space to place, click an atom to change its element, drag from an atom to bond a new one",
            "Bond: drag between atoms (or click two) to bond them; click a bond to cycle single, double, triple",
            "Move: drag atoms around; drag empty space to rotate (3D) or pan (2D); Shift+drag pans",
            "Erase: click an atom or a bond to remove it",
            "Inspect: click an atom or a bond to see its details"
        };
        m_hint = hints[static_cast<int>(m_tool)];
        m_sel_atom = -1;
        changed(false);
        return tool_name(m_tool);
    }
    if (c == "element") {
        std::string s = arg(a, 1);
        const std::size_t sp = s.find(' ');
        if (sp != std::string::npos) s = s.substr(0, sp);
        const fc::Element* e = fc::find_element(s);
        if (e) m_element = e->atomic_number;
        changed(false);
        return e ? e->symbol : "";
    }
    if (c == "charge") {
        push_undo();
        m_cov.charge = static_cast<int>(num(arg(a, 1), 0.0));
        changed(true);
        return std::to_string(m_cov.charge);
    }
    if (c == "neutrons") {
        if (m_sel_atom < 0 || m_sel_atom >= static_cast<int>(m_cov.atoms.size())) return "error\tselect an atom first";
        push_undo();
        m_cov.atoms[static_cast<std::size_t>(m_sel_atom)].neutrons = std::max(-200, std::min(200, static_cast<int>(num(arg(a, 1), 0.0))));
        changed(true);
        return std::to_string(m_cov.atoms[static_cast<std::size_t>(m_sel_atom)].neutrons);
    }
    if (c == "selected_neutrons") {
        if (m_sel_atom < 0 || m_sel_atom >= static_cast<int>(m_cov.atoms.size())) return "";
        return std::to_string(m_cov.atoms[static_cast<std::size_t>(m_sel_atom)].neutrons);
    }
    if (c == "frequency") {
        m_frequency = std::max(1e-3, num(arg(a, 1), 60.0));
        changed(false);
        return format_number(m_frequency, 6);
    }
    if (c == "order") {
        if (m_sel_bond >= 0 && m_sel_bond < static_cast<int>(m_cov.bonds.size())) {
            push_undo();
            m_cov.bonds[static_cast<std::size_t>(m_sel_bond)].order = std::max(1, std::min(3, static_cast<int>(num(arg(a, 1), 1.0))));
            changed(true);
            return "ok";
        }
        return "";
    }
    if (c == "option") {
        std::string key = arg(a, 1);
        if (key == "neighbors" && m_mode == BondingMode::Metallic) key = "metal_neighbors";
        m_options[key] = arg(a, 2) == "1" || arg(a, 2) == "true";
        changed(false);
        return arg(a, 2);
    }
    if (c == "tp") {
        try {
            m_temperature = fc::to_si(num(arg(a, 1), 298.15), arg(a, 2).empty() ? "K" : arg(a, 2), fc::Dimension::Temperature);
            m_pressure = fc::to_si(num(arg(a, 3), 1.0), arg(a, 4).empty() ? "atm" : arg(a, 4), fc::Dimension::Pressure);
        } catch (const std::exception& e) {
            return std::string("error\t") + e.what();
        }
        m_temperature = std::max(0.01, m_temperature);
        m_pressure = std::max(1e-6, m_pressure);
        changed(true);
        return format_t();
    }
    if (c == "clear") {
        if (m_mode == BondingMode::Covalent) {
            push_undo();
            m_cov = CovalentState{};
            m_sel_atom = m_sel_bond = -1;
            frame_view();
            changed(true);
        }
        return "ok";
    }
    if (c == "undo") {
        if (m_undo.empty()) return "";
        m_redo.push_back(m_cov);
        m_cov = m_undo.back();
        m_undo.pop_back();
        m_sel_atom = m_sel_bond = -1;
        changed(true);
        return "ok";
    }
    if (c == "redo") {
        if (m_redo.empty()) return "";
        m_undo.push_back(m_cov);
        m_cov = m_redo.back();
        m_redo.pop_back();
        m_sel_atom = m_sel_bond = -1;
        changed(true);
        return "ok";
    }
    if (c == "fill_h") {
        push_undo();
        fill_hydrogens();
        changed(true);
        return "ok";
    }
    if (c == "snap") {
        push_undo();
        const std::string what = snap_to_correct();
        changed(true);
        return "ok\t" + what;
    }
    if (c == "tidy") {
        push_undo();
        tidy();
        frame_view();
        changed(true);
        return "ok";
    }
    if (c == "fix_lewis") {
        push_undo();
        apply_best_lewis();
        changed(true);
        return "ok";
    }
    if (c == "resonance") {
        push_undo();
        apply_resonance(static_cast<int>(num(arg(a, 1), 0.0)));
        changed(true);
        return "ok";
    }
    if (c == "frame") {
        frame_view();
        changed(false);
        return "ok";
    }
    if (c == "load_smiles") {
        push_undo();
        const bool ok = load_smiles(arg(a, 1), static_cast<int>(num(arg(a, 2), 0.0)), arg(a, 3));
        return ok ? "ok" : "error\tcouldn't read that SMILES";
    }
    if (c == "load_ion" || c == "load_compound") {
        const std::string text = arg(a, 1);
        std::string smiles;
        int charge = 0;
        std::string name;
        if (const fc::PolyatomicIonEntry* ion = fc::find_polyatomic_ion(text)) {
            smiles = ion->smiles;
            charge = ion->charge;
            name = ion->name;
        } else {
            const auto found = fc::find_common_compounds(text);
            if (!found.empty()) {
                smiles = found.front()->smiles;
                name = found.front()->name;
            }
        }
        if (smiles.empty()) {
            try {
                fc::Molecule m = fc::Molecule::from_smiles(text);
                if (m.atom_count() > 0) smiles = text;
            } catch (...) {
            }
        }
        if (smiles.empty()) smiles = central_atom_smiles(text);
        if (smiles.empty()) return "error\tno structure for '" + text + "'";
        if (name.empty()) name = text;
        push_undo();
        m_mode = BondingMode::Covalent;
        return load_smiles(smiles, charge, name) ? "ok\t" + name : "error\tcouldn't build " + text;
    }
    if (c == "ionic") {
        set_ionic(arg(a, 1), arg(a, 2));
        return m_ionic.ok ? "ok\t" + m_ionic.formula : "error\t" + m_lattice_error;
    }
    if (c == "metal") {
        set_metal(arg(a, 1));
        return m_metal.ok ? "ok\t" + m_metal.name : "error\t" + m_lattice_error;
    }
    if (c == "cells") {
        m_cells = std::max(1, std::min(5, static_cast<int>(num(arg(a, 1), 2.0))));
        rebuild_lattice();
        restart_assembly();
        frame_view();
        changed(true);
        return std::to_string(m_cells);
    }
    if (c == "assemble") {
        restart_assembly();
        changed(false);
        return "ok";
    }
    if (c == "readouts") {
        if (m_mode == BondingMode::Ionic) return ionic_readouts();
        if (m_mode == BondingMode::Metallic) return metallic_readouts();
        return covalent_readouts();
    }
    if (c == "selection") {
        if (m_mode != BondingMode::Covalent) return lattice_selection();
        return selection_readouts();
    }
    if (c == "resonance_list") {
        ensure_analysis();
        std::string out;
        for (std::size_t i = 0; i < m_analysis.resonance.structures.size(); ++i) out += std::to_string(i) + "\tstructure " + std::to_string(i + 1) + "\n";
        return out;
    }
    if (c == "name") {
        m_cov.name = arg(a, 1);
        return m_cov.name;
    }
    if (c == "json") return to_json();
    if (c == "save") {
        std::string name = arg(a, 1);
        if (name.empty() && m_mode == BondingMode::Covalent) name = m_cov.name;
        if (name.empty() && m_mode == BondingMode::Ionic && m_ionic.ok) name = m_ionic.name;
        if (name.empty() && m_mode == BondingMode::Metallic && m_metal.ok) name = m_metal.name;
        if (m_mode == BondingMode::Covalent) {
            ensure_analysis();
            m_cov.name = name;
        }
        std::string error;
        if (!save_chemistry(kind_name(), name, to_json(), &error)) return "error\t" + error;
        return "ok\t" + (chemistry_save_root() / kind_name() / (safe_file_name(name) + ".json")).string();
    }
    if (c == "saved") {
        std::string out;
        const std::string kind = arg(a, 1).empty() ? kind_name() : arg(a, 1);
        for (const std::string& n : list_saved_chemistry(kind)) out += n + "\n";
        return out;
    }
    if (c == "load") {
        std::string json;
        if (!load_chemistry(arg(a, 1), arg(a, 2), json)) return "error\tcan't read " + arg(a, 2);
        std::string error;
        if (!from_json(json, &error)) return "error\t" + error;
        if (m_mode == BondingMode::Covalent && m_cov.name.empty()) m_cov.name = arg(a, 2);
        frame_view();
        changed(true);
        return "ok\t" + kind_name();
    }
    if (c == "delete") {
        return delete_saved_chemistry(arg(a, 1), arg(a, 2)) ? "ok" : "error\tnothing to delete";
    }
    if (c == "save_dir") return chemistry_save_root().string();
    if (c == "hover") return m_hint;
    return "error\tunknown command " + c;
}

std::unique_ptr<ChemView> make_bonding_view(ChemShared& shared) {
    return std::unique_ptr<ChemView>(new BondingView(shared));
}

} // namespace sapling

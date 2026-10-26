#include "molecule_draw.hpp"
#include "view_util.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

namespace sapling {

using namespace chem_util;

namespace {

struct PropField {
    const char* key;
    const char* label;
    const char* dimension;
};

const PropField PROP_FIELDS[] = {
    { "critical_temperature", "Critical temperature", "temperature" },
    { "critical_pressure", "Critical pressure", "pressure" },
    { "triple_temperature", "Triple-point temperature", "temperature" },
    { "triple_pressure", "Triple-point pressure", "pressure" },
    { "melting_point", "Normal melting point", "temperature" },
    { "boiling_point", "Normal boiling point", "temperature" },
    { "vaporization_enthalpy", "\xCE\x94H vaporization", "molar energy" },
    { "fusion_enthalpy", "\xCE\x94H fusion", "molar energy" },
    { "liquid_density", "Liquid density", "density" },
    { "solid_density", "Solid density", "density" },
    { "molar_mass", "Molar mass", "molar mass" },
    { "acentric_factor", "Acentric factor \xCF\x89", "dimensionless" },
};

double* prop_slot(fc::PhaseProperties& p, const std::string& key) {
    if (key == "critical_temperature") return &p.critical_temperature;
    if (key == "critical_pressure") return &p.critical_pressure;
    if (key == "triple_temperature") return &p.triple_temperature;
    if (key == "triple_pressure") return &p.triple_pressure;
    if (key == "melting_point") return &p.melting_point;
    if (key == "boiling_point") return &p.boiling_point;
    if (key == "vaporization_enthalpy") return &p.vaporization_enthalpy;
    if (key == "fusion_enthalpy") return &p.fusion_enthalpy;
    if (key == "liquid_density") return &p.liquid_density;
    if (key == "solid_density") return &p.solid_density;
    if (key == "molar_mass") return &p.molar_mass;
    if (key == "acentric_factor") return &p.acentric_factor;
    return nullptr;
}

const char* estimate_name(const std::string& key) {
    if (key == "critical_temperature") return "critical temperature";
    if (key == "critical_pressure") return "critical pressure";
    if (key == "triple_temperature") return "triple point temperature";
    if (key == "triple_pressure") return "triple point pressure";
    if (key == "melting_point") return "melting point";
    if (key == "boiling_point") return "boiling point";
    if (key == "vaporization_enthalpy") return "enthalpy of vaporization";
    if (key == "fusion_enthalpy") return "enthalpy of fusion";
    if (key == "liquid_density") return "liquid density";
    if (key == "solid_density") return "solid density";
    if (key == "acentric_factor") return "acentric factor";
    if (key == "molar_mass") return "molar mass";
    return "";
}

bool is_estimated(const fc::PhaseProperties& p, const std::string& key) {
    const std::string name = estimate_name(key);
    for (const std::string& e : p.estimated) {
        if (e == key || e == name) return true;
    }
    return false;
}

constexpr int PLASMA_CELL = 99;
constexpr double H_PLANCK = 6.62607015e-34;
constexpr double K_BOLTZ = 1.380649e-23;
constexpr double M_ELECTRON = 9.1093837139e-31;
constexpr double AVOGADRO = 6.02214076e23;

Color phase_color(fc::Phase ph) {
    switch (ph) {
        case fc::Phase::Solid: return rgb(116, 134, 170);
        case fc::Phase::Liquid: return rgb(76, 140, 220);
        case fc::Phase::Gas: return rgb(222, 196, 120);
        case fc::Phase::Supercritical: return rgb(186, 140, 230);
        default: return rgb(90, 90, 100);
    }
}

Color cell_color(int code) {
    if (code == PLASMA_CELL) return rgb(245, 140, 110);
    return phase_color(static_cast<fc::Phase>(code));
}

Color region_color(fc::PvtRegion r) {
    switch (r) {
        case fc::PvtRegion::Solid: return rgb(116, 134, 170);
        case fc::PvtRegion::Liquid: return rgb(76, 140, 220);
        case fc::PvtRegion::Gas: return rgb(222, 196, 120);
        case fc::PvtRegion::Supercritical: return rgb(186, 140, 230);
        case fc::PvtRegion::LiquidVapor: return rgb(120, 190, 170);
        case fc::PvtRegion::SolidLiquid: return rgb(100, 150, 200);
        default: return rgb(170, 160, 150);
    }
}

const char* cell_word(int code) {
    if (code == PLASMA_CELL) return "PLASMA";
    switch (static_cast<fc::Phase>(code)) {
        case fc::Phase::Solid: return "SOLID";
        case fc::Phase::Liquid: return "LIQUID";
        case fc::Phase::Gas: return "GAS";
        case fc::Phase::Supercritical: return "SUPERCRITICAL FLUID";
        default: return "";
    }
}

class PhaseView : public ChemView {
public:
    explicit PhaseView(ChemShared& shared) : m_shared(shared) {
        m_cam.yaw = 0.7;
        m_cam.pitch = 0.35;
        set_substance("water");
    }

    void draw(Canvas& canvas, const Theme& theme) override {
        canvas.clear(theme.background);
        m_cam.width = canvas.width();
        m_cam.height = canvas.height();
        if (!m_diag.ok) {
            draw_hud(canvas, theme, { "Phase diagram", m_error.empty() ? std::string("type a substance, pick an element, or use something you saved") : m_error }, false);
            return;
        }
        if (m_3d) draw_3d(canvas, theme);
        else draw_2d(canvas, theme);
    }

    bool mouse_down(Vec2 p, int button, unsigned mods) override {
        m_down = true;
        m_button = button;
        m_mods = mods;
        m_last = p;
        m_follow = -1;
        if (m_3d) {
            if (button == 1 && !(mods & k_shift)) {
                probe(p);
                return true;
            }
            return false;
        }
        if (button != 1 || !m_frame.inside(p)) return false;
        const std::string hit = marker_at(p);
        if (!hit.empty()) {
            m_selection = hit;
            m_down = false;
            ++m_version;
            return true;
        }
        const int curve = curve_at(p);
        if (curve >= 0) {
            m_follow = curve;
            m_selection = curve == 0 ? "sublimation" : curve == 1 ? "vaporization" : "fusion";
            follow_to(p);
            return true;
        }
        set_point_from_screen(p);
        return true;
    }

    bool mouse_move(Vec2 p, bool) override {
        const Vec2 d = p - m_last;
        m_last = p;
        m_mouse = p;
        if (m_3d) {
            if (!m_down) return false;
            if (m_button == 1 && !(m_mods & k_shift)) probe(p);
            else m_cam.orbit(d.x, d.y);
            return true;
        }
        if (m_down && m_button == 1) {
            if (m_follow >= 0) follow_to(p);
            else if (m_frame.inside(p)) set_point_from_screen(p);
        }
        return true;
    }

    bool mouse_up(Vec2, int) override {
        m_down = false;
        return false;
    }

    bool wheel(Vec2, double dy) override {
        if (!m_3d) return false;
        m_cam.zoom_by(dy > 0 ? 1.12 : 1.0 / 1.12);
        return true;
    }

    bool double_click(Vec2) override {
        m_cam.yaw = 0.7;
        m_cam.pitch = 0.35;
        m_cam.zoom = 0.0;
        return true;
    }

    std::string command(const std::vector<std::string>& a) override {
        const std::string c = arg(a, 0);
        if (c == "substance") return set_substance(arg(a, 1)) ? "ok\t" + m_props.name : "error\t" + m_error;
        if (c == "element") {
            const fc::Element* e = element_of(arg(a, 1));
            if (!e) return "error\tpick an element";
            m_overrides.clear();
            m_base = fc::phase_properties(*e);
            m_base_source = std::string("element ") + e->name;
            return build() ? "ok\t" + m_props.name : "error\t" + m_error;
        }
        if (c == "use_saved") {
            std::string json, error;
            if (!load_chemistry(arg(a, 1), arg(a, 2), json)) return "error\tcan't read " + arg(a, 2);
            fc::PhaseProperties p = phase_from_saved(arg(a, 1), json, &error);
            if (!error.empty() && !p.usable()) return "error\t" + error;
            m_overrides.clear();
            m_base = p;
            if (m_base.name.empty()) m_base.name = arg(a, 2);
            m_base_source = "saved " + arg(a, 1) + " \"" + arg(a, 2) + "\"";
            return build() ? "ok\t" + m_props.name : "error\t" + m_error;
        }
        if (c == "view") {
            m_3d = arg(a, 1) == "3d";
            if (m_3d) ensure_surface();
            return m_3d ? "3d" : "2d";
        }
        if (c == "point") {
            m_point_t = to_si_or(num(arg(a, 1), 298.15), arg(a, 2).empty() ? "K" : arg(a, 2), fc::Dimension::Temperature, m_point_t);
            m_point_p = to_si_or(num(arg(a, 3), 1.0), arg(a, 4).empty() ? "atm" : arg(a, 4), fc::Dimension::Pressure, m_point_p);
            m_point_p = std::max(1e-9, m_point_p);
            ++m_version;
            return describe(m_point_t, m_point_p);
        }
        if (c == "point_values") {
            const UnitPrefs& u = m_shared.units;
            const std::string tu = arg(a, 1).empty() ? u.get("temperature") : arg(a, 1);
            const std::string pu = arg(a, 2).empty() ? u.get("pressure") : arg(a, 2);
            const bool si = tu == "K" && pu == "Pa";
            return (si ? format_number(m_point_t, 9) : fixed(from_si_or(m_point_t, tu, fc::Dimension::Temperature), 2)) + "\t" + format_number(from_si_or(m_point_p, pu, fc::Dimension::Pressure), si ? 9 : 5);
        }
        if (c == "props") {
            std::string out;
            for (const PropField& f : PROP_FIELDS) {
                double* slot = prop_slot(m_props, f.key);
                const std::string dim = f.dimension;
                const std::string unit = dim == "dimensionless" ? std::string() : m_shared.units.get(dim);
                double v = slot ? *slot : fc::NO_VALUE;
                const double si = v;
                std::string shown = finite(si) ? (dim == "dimensionless" ? format_number(si, 5) : format_number(from_si_or(si, unit, fc::dimension_named(dim)), 6)) : std::string();
                const bool est = is_estimated(m_props, f.key);
                const bool over = m_overrides.count(f.key) > 0;
                out += std::string(f.key) + "\t" + f.label + "\t" + shown + "\t" + dim + "\t" + unit + "\t" + (over ? "yours" : est ? "estimated" : "data") + "\n";
            }
            return out;
        }
        if (c == "set_prop") {
            const std::string key = arg(a, 1);
            if (arg(a, 2).empty()) {
                m_overrides.erase(key);
            } else {
                std::string dim;
                for (const PropField& f : PROP_FIELDS) {
                    if (key == f.key) dim = f.dimension;
                }
                if (dim.empty()) return "error\tunknown property";
                double v = num(arg(a, 2), fc::NO_VALUE);
                if (!finite(v)) return "error\tnot a number";
                if (dim != "dimensionless") v = to_si_or(v, arg(a, 3), fc::dimension_named(dim), v);
                m_overrides[key] = v;
            }
            return build() ? "ok" : "error\t" + m_error;
        }
        if (c == "reset_props") {
            m_overrides.clear();
            return build() ? "ok" : "error\t" + m_error;
        }
        if (c == "readouts") return readouts();
        if (c == "selection") return selection_readouts();
        if (c == "axes") {
            m_tmin = arg(a, 1).empty() ? fc::NO_VALUE : std::max(0.0, num(arg(a, 1), fc::NO_VALUE));
            m_tmax = arg(a, 2).empty() ? fc::NO_VALUE : std::min(1e7, num(arg(a, 2), fc::NO_VALUE));
            m_pmin = arg(a, 3).empty() ? fc::NO_VALUE : std::max(0.0, num(arg(a, 3), fc::NO_VALUE));
            m_pmax = arg(a, 4).empty() ? fc::NO_VALUE : std::min(1e13, num(arg(a, 4), fc::NO_VALUE));
            m_region_key.clear();
            ++m_version;
            return "ok";
        }
        if (c == "logp" || c == "logt" || c == "extreme") {
            const bool on = arg(a, 1) == "1";
            if (c == "logp") m_logp = on;
            else if (c == "logt") m_logt = on;
            else m_extreme = on;
            m_region_key.clear();
            ++m_version;
            return arg(a, 1);
        }
        if (c == "axis_values") {
            double t0, t1, p0, p1;
            ranges(t0, t1, p0, p1);
            return format_number(t0, 6) + "\t" + format_number(t1, 6) + "\t" + format_number(p0, 6) + "\t" + format_number(p1, 6);
        }
        if (c == "use_json") {
            JsonValue v;
            std::string error;
            if (!parse_json(arg(a, 1), v, &error)) return "error\t" + error;
            const std::string kind = v.string_or("type", "covalent");
            fc::PhaseProperties p = phase_from_saved(kind, arg(a, 1), &error);
            if (!p.usable() && !error.empty()) return "error\t" + error;
            m_overrides.clear();
            m_base = p;
            if (m_base.name.empty()) m_base.name = v.string_or("name", v.string_or("formula", "your " + kind + " build"));
            m_base_source = "built here (" + kind + ")";
            return build() ? "ok\t" + m_props.name : "error\t" + m_error;
        }
        if (c == "version") return std::to_string(m_version);
        if (c == "saved") {
            std::string out;
            for (const std::string& kind : saved_chemistry_kinds()) {
                for (const std::string& n : list_saved_chemistry(kind)) out += kind + "\t" + n + "\n";
            }
            return out;
        }
        return "error\tunknown command " + c;
    }

private:
    bool set_substance(const std::string& text) {
        m_error.clear();
        m_overrides.clear();
        const std::string t = trim(text);
        if (t.empty()) {
            m_error = "type a substance";
            return false;
        }
        fc::PhaseProperties p;
        try {
            p = fc::find_phase_properties(t);
        } catch (...) {
        }
        if (!p.usable()) {
            const fc::Element* e = fc::find_element(t);
            if (e) p = fc::phase_properties(*e);
        }
        if (!p.usable()) {
            int q = 0;
            std::string name;
            std::optional<fc::Molecule> m = molecule_from_text(t, &q, &name, nullptr);
            if (m) {
                try {
                    p = fc::phase_properties(*m, name.empty() ? t : name);
                } catch (...) {
                }
            }
        }
        if (!p.usable()) {
            m_error = "no phase data for '" + t + "' (try a name like water, CO2, methane, a SMILES string, or an element)";
            m_diag = fc::PhaseDiagram{};
            ++m_version;
            return false;
        }
        m_base = p;
        m_base_source = p.source.empty() ? "data" : p.source;
        return build();
    }

    bool build() {
        m_error.clear();
        fc::PhaseProperties p = m_base;
        for (const auto& kv : m_overrides) {
            if (double* slot = prop_slot(p, kv.first)) *slot = kv.second;
            const std::string name = estimate_name(kv.first);
            p.estimated.erase(std::remove_if(p.estimated.begin(), p.estimated.end(), [&](const std::string& e) { return e == kv.first || e == name; }), p.estimated.end());
        }
        try {
            fc::complete_phase_properties(p);
            m_props = p;
            m_diag = fc::phase_diagram(p, 240);
        } catch (const std::exception& e) {
            m_error = e.what();
            m_diag = fc::PhaseDiagram{};
        }
        if (!m_diag.ok && m_error.empty()) m_error = m_diag.notes.empty() ? "couldn't draw a phase diagram from these properties" : m_diag.notes.front();
        m_surface = fc::PvtSurface{};
        m_region_key.clear();
        ++m_version_build;
        m_cam.zoom = 0.0;
        if (m_diag.ok) {
            m_point_t = std::max(m_diag.t_min, std::min(m_diag.t_max, 298.15));
            m_point_p = 101325.0;
            if (m_3d) ensure_surface();
        }
        ++m_version;
        return m_diag.ok;
    }

    void ensure_surface() {
        if (m_surface.ok || !m_diag.ok) return;
        try {
            m_surface = fc::pvt_surface(m_diag, 44, 60);
        } catch (...) {
        }
    }

    double disp_t(double si) const { return from_si_or(si, m_shared.units.get("temperature"), fc::Dimension::Temperature); }
    double disp_p(double si) const { return from_si_or(si, m_shared.units.get("pressure"), fc::Dimension::Pressure); }
    double si_t(double d) const { return to_si_or(d, m_shared.units.get("temperature"), fc::Dimension::Temperature, d); }
    double si_p(double d) const { return to_si_or(d, m_shared.units.get("pressure"), fc::Dimension::Pressure, d); }

    void ranges(double& t0, double& t1, double& p0, double& p1) const {
        t0 = finite(m_tmin) ? m_tmin : m_diag.t_min;
        t1 = finite(m_tmax) ? m_tmax : m_diag.t_max;
        p0 = finite(m_pmin) ? m_pmin : (m_logp ? m_diag.p_min : 0.0);
        p1 = finite(m_pmax) ? m_pmax : m_diag.p_max;
        if (m_logt) t0 = std::max(t0, 1e-9);
        if (m_logp) p0 = std::max(p0, 1e-12);
        if (!(t1 > t0)) t1 = t0 + 1.0;
        if (!(p1 > p0)) p1 = m_logp ? p0 * 10.0 : p0 + 1.0;
    }

    double ax_t(double si) const { return m_logt ? si : disp_t(si); }
    double ax_t_inv(double v) const { return m_logt ? v : si_t(v); }
    double ax_p(double si) const { return disp_p(si); }
    double ax_p_inv(double v) const { return si_p(v); }

    Vec2 screen_of(double t, double p) const { return m_frame.to_screen(ax_t(t), ax_p(p)); }

    void set_point_from_screen(Vec2 p) {
        m_point_t = std::max(0.0, ax_t_inv(m_frame.x_of(p.x)));
        m_point_p = std::max(0.0, ax_p_inv(m_frame.y_of(p.y)));
        m_selection.clear();
        ++m_version;
    }

    std::string marker_at(Vec2 p) const {
        const std::pair<const char*, fc::PhasePoint> marks[] = { { "triple", m_diag.triple }, { "critical", m_diag.critical }, { "normal_melting", m_diag.normal_melting }, { "normal_boiling", m_diag.normal_boiling } };
        for (const auto& m : marks) {
            if (!finite(m.second.t) || !finite(m.second.p) || m.second.p <= 0) continue;
            if (length(screen_of(m.second.t, m.second.p) - p) < 10.0) return m.first;
        }
        return std::string();
    }

    int curve_at(Vec2 p) const {
        const std::vector<fc::PhasePoint>* curves[] = { &m_diag.sublimation, &m_diag.vaporization, &m_diag.fusion };
        int best = -1;
        double bd = 7.0;
        for (int c = 0; c < 3; ++c) {
            const auto& pts = *curves[c];
            for (std::size_t i = 0; i + 1 < pts.size(); ++i) {
                if (!(pts[i].p > 0) || !(pts[i + 1].p > 0)) continue;
                const Vec2 a = screen_of(pts[i].t, pts[i].p), b = screen_of(pts[i + 1].t, pts[i + 1].p);
                const Vec2 d = b - a;
                const double dd = d.x * d.x + d.y * d.y;
                double t = dd > 0 ? ((p.x - a.x) * d.x + (p.y - a.y) * d.y) / dd : 0.0;
                t = std::max(0.0, std::min(1.0, t));
                const double dist = length(p - (a + d * t));
                if (dist < bd) {
                    bd = dist;
                    best = c;
                }
            }
        }
        return best;
    }

    void follow_to(Vec2 p) {
        const double tt = m_props.triple_temperature, tc = m_props.critical_temperature;
        double t = ax_t_inv(m_frame.x_of(p.x));
        try {
            if (m_follow == 0) {
                t = std::max(m_diag.t_min, std::min(tt, t));
                m_point_t = t;
                m_point_p = m_diag.sublimation_pressure(t);
            } else if (m_follow == 1) {
                t = std::max(tt, std::min(tc, t));
                m_point_t = t;
                m_point_p = m_diag.vapor_pressure(t);
            } else {
                const double pp = ax_p_inv(m_frame.y_of(p.y));
                m_point_p = std::max(m_props.triple_pressure, pp);
                m_point_t = m_diag.melting_temperature(m_point_p);
            }
        } catch (...) {
        }
        ++m_version;
    }

    void probe(Vec2 p) {
        int best = -1;
        double bd = 14.0, depth = -1e300;
        for (std::size_t i = 0; i < m_probe_vertices.size(); ++i) {
            const double d = length(m_probe_vertices[i].s - p);
            if (d < bd && (d < bd - 3.0 || m_probe_vertices[i].depth > depth)) {
                bd = d;
                depth = m_probe_vertices[i].depth;
                best = static_cast<int>(i);
            }
        }
        if (best < 0) return;
        m_probe_t = m_probe_vertices[static_cast<std::size_t>(best)].it;
        m_probe_v = m_probe_vertices[static_cast<std::size_t>(best)].iv;
        m_point_t = m_surface.temperatures[static_cast<std::size_t>(m_probe_t)];
        m_point_p = std::max(0.0, m_surface.pressure(m_probe_t, m_probe_v));
        m_selection = "probe";
        ++m_version;
    }

    double ionization_energy_ev() const {
        double best = 1e300;
        try {
            const fc::Formula f = fc::Formula::parse(m_props.formula);
            for (const auto& c : f.composition()) {
                const double ie = c.first->ionization_energy(1);
                if (ie == ie && ie > 0) best = std::min(best, ie);
            }
        } catch (...) {
        }
        return best < 1e300 ? best : 12.0;
    }

    double ionization_fraction(double t, double p) const {
        if (!(t > 0) || !(p > 0)) return 0.0;
        const double chi = ionization_energy_ev() * 1.602176634e-19;
        const double a = std::pow(2.0 * 3.14159265358979 * M_ELECTRON * K_BOLTZ * t / (H_PLANCK * H_PLANCK), 1.5);
        const double expo = -chi / (K_BOLTZ * t);
        if (expo < -700.0) return 0.0;
        const double k = K_BOLTZ * t / p * a * std::exp(expo);
        return std::sqrt(k / (1.0 + k));
    }

    bool boson() const {
        long long neutrons = 0;
        try {
            const fc::Formula f = fc::Formula::parse(m_props.formula);
            for (const auto& c : f.composition()) {
                const int a = c.first->most_abundant_mass_number ? c.first->most_abundant_mass_number : c.first->reference_mass_number;
                neutrons += std::llround(c.second) * (a - c.first->atomic_number);
            }
        } catch (...) {
        }
        return neutrons % 2 == 0;
    }

    double particle_mass() const {
        return (finite(m_props.molar_mass) && m_props.molar_mass > 0 ? m_props.molar_mass : 18.0) / 1000.0 / AVOGADRO;
    }

    double quantum_pressure(double t) const {
        const double m = particle_mass();
        if (boson()) return 2.612 * K_BOLTZ * t * std::pow(2.0 * 3.14159265358979 * m * K_BOLTZ * t / (H_PLANCK * H_PLANCK), 1.5);
        const double hbar = H_PLANCK / (2.0 * 3.14159265358979);
        const double n = std::pow(2.0 * m * K_BOLTZ * t / (hbar * hbar), 1.5) / (3.0 * 3.14159265358979 * 3.14159265358979);
        return n * K_BOLTZ * t;
    }

    std::string point_details(const std::string& key) const {
        const UnitPrefs& u = m_shared.units;
        std::ostringstream o;
        fc::PhasePoint pt;
        std::string title, meaning;
        if (key == "triple") {
            pt = m_diag.triple;
            title = "Triple point";
            meaning = "Solid, liquid and gas all coexist here. Below this pressure the liquid can't exist, so the solid sublimes straight to gas.";
        } else if (key == "critical") {
            pt = m_diag.critical;
            title = "Critical point";
            meaning = "The end of the liquid\xE2\x80\x93gas line: above this temperature and pressure the liquid and gas become one supercritical fluid, with no boiling.";
        } else if (key == "normal_melting") {
            pt = m_diag.normal_melting;
            title = "Normal melting point";
            meaning = "Where the melting line crosses 1 atm: the melting point you'd measure in the open air.";
        } else {
            pt = m_diag.normal_boiling;
            title = "Normal boiling point";
            meaning = "Where the vapor pressure equals 1 atm: the boiling point at sea level.";
        }
        o << "#" << title << "\n";
        o << "Temperature\t" << u.format(pt.t, "temperature", 5) << "\n";
        o << "Pressure\t" << u.format(pt.p, "pressure", 5) << "\n";
        o << meaning << "\n";
        try {
            if (key == "critical") {
                if (finite(m_props.critical_volume)) o << "Critical molar volume\t" << u.format(m_props.critical_volume, "molar volume", 4) << "\n";
                if (finite(m_props.critical_volume) && m_props.molar_mass > 0) o << "Critical density\t" << u.format(m_props.molar_mass / 1000.0 / m_props.critical_volume, "density", 4) << "\n";
                o << "Compressibility Z\xE1\xB6\x9C\t" << format_number(pt.p * m_props.critical_volume / (8.314462618 * pt.t), 4) << "\n";
            }
            const double h = 0.002 * pt.t;
            if (key == "triple" || key == "normal_boiling") {
                const double dpdt = (m_diag.vapor_pressure(pt.t + h) - m_diag.vapor_pressure(std::max(pt.t - h, 1e-6))) / (2.0 * h);
                o << "Slope of the boiling line\t" << u.format(dpdt, "pressure", 4) << " per K\n";
                if (finite(m_props.vaporization_enthalpy)) o << "\xCE\x94H vaporization\t" << u.format(m_props.vaporization_enthalpy, "molar energy", 4) << "\n";
            }
            if (key == "triple" || key == "normal_melting") {
                const double dv = m_props.fusion_volume_change();
                if (finite(dv) && finite(m_props.fusion_enthalpy) && std::fabs(dv) > 0) o << "Slope of the melting line\t" << u.format(m_props.fusion_enthalpy / (pt.t * dv), "pressure", 4) << " per K (Clapeyron)\n";
                if (finite(m_props.fusion_enthalpy)) o << "\xCE\x94H fusion\t" << u.format(m_props.fusion_enthalpy, "molar energy", 4) << "\n";
            }
            if (key == "triple") {
                const double sub = finite(m_props.vaporization_enthalpy) && finite(m_props.fusion_enthalpy) ? m_props.vaporization_enthalpy + m_props.fusion_enthalpy : fc::NO_VALUE;
                if (finite(sub)) o << "\xCE\x94H sublimation\t" << u.format(sub, "molar energy", 4) << " (fusion + vaporization)\n";
            }
            if (key == "normal_boiling" && finite(m_props.vaporization_enthalpy)) o << "Trouton's ratio \xCE\x94S\t" << format_number(m_props.vaporization_enthalpy / pt.t, 4) << " J/(mol\xC2\xB7K) (about 88 for most liquids)\n";
        } catch (...) {
        }
        return o.str();
    }

    std::string selection_readouts() const {
        if (!m_diag.ok) return std::string();
        const UnitPrefs& u = m_shared.units;
        std::ostringstream o;
        if (m_selection == "triple" || m_selection == "critical" || m_selection == "normal_melting" || m_selection == "normal_boiling") return point_details(m_selection);
        if (m_selection == "sublimation" || m_selection == "vaporization" || m_selection == "fusion") {
            static const char* what[] = { "solid \xE2\x87\x8C gas", "liquid \xE2\x87\x8C gas", "solid \xE2\x87\x8C liquid" };
            const int c = m_selection == "sublimation" ? 0 : m_selection == "vaporization" ? 1 : 2;
            o << "#On the " << m_selection << " line\n";
            o << "Coexisting\t" << what[c] << "\n";
            o << "Temperature\t" << u.format(m_point_t, "temperature", 5) << "\n";
            o << "Pressure\t" << u.format(m_point_p, "pressure", 5) << "\n";
            if (c == 1) o << "This liquid boils at " << u.format(m_point_t, "temperature", 4) << " when the pressure is " << u.format(m_point_p, "pressure", 4) << ".\n";
            if (c == 0) o << "Below the triple point the solid turns straight into gas (like dry ice or frost vanishing).\n";
            if (c == 2) o << "Drag up the melting line to see how pressure changes the melting point.\n";
            return o.str();
        }
        if (m_selection == "probe" && m_probe_t >= 0 && m_probe_v >= 0) {
            const double v = m_surface.specific_volumes[static_cast<std::size_t>(m_probe_v)];
            o << "#Point on the surface\n";
            o << "Temperature\t" << u.format(m_point_t, "temperature", 5) << "\n";
            o << "Specific volume\t" << u.format(v, "specific volume", 4) << "\n";
            o << "Pressure\t" << u.format(m_point_p, "pressure", 5) << "\n";
            o << "Region\t" << fc::to_string(m_surface.region(m_probe_t, m_probe_v)) << "\n";
            if (m_props.molar_mass > 0) o << "Z = Pv/RT\t" << format_number(m_point_p * v * m_props.molar_mass / 1000.0 / (8.314462618 * m_point_t), 4) << "\n";
            return o.str();
        }
        o << "#Selection\nClick the triple point, critical point or a normal point for details; drag along a line to follow it. In 3D, drag across the surface to read points (Shift+drag or right-drag rotates).\n";
        return o.str();
    }

    std::string describe(double t, double p) const {
        const UnitPrefs& u = m_shared.units;
        const fc::Phase ph = m_diag.phase_at(t, p);
        std::string out = fc::to_string(ph);
        try {
            if (ph == fc::Phase::Liquid && p < m_props.critical_pressure) out += "; boils at " + u.format(m_diag.boiling_temperature(p), "temperature", 4) + " at this pressure";
            if (ph == fc::Phase::Solid) {
                if (p >= m_props.triple_pressure) out += "; melts at " + u.format(m_diag.melting_temperature(p), "temperature", 4) + " here";
                else out += "; sublimes at " + u.format(m_diag.sublimation_temperature(p), "temperature", 4) + " here";
            }
            if (ph == fc::Phase::Gas && t < m_props.critical_temperature && t > m_props.triple_temperature) out += "; condenses above " + u.format(m_diag.vapor_pressure(t), "pressure", 4) + " at this temperature";
        } catch (...) {
        }
        return out;
    }

    void draw_2d(Canvas& canvas, const Theme& theme) {
        const UnitPrefs& u = m_shared.units;
        PlotFrame& f = m_frame;
        f.x0 = 78.0;
        f.y0 = 84.0;
        f.w = std::max(50.0, canvas.width() - f.x0 - 24.0);
        f.h = std::max(50.0, canvas.height() - f.y0 - 58.0);
        double t0, t1, p0, p1;
        ranges(t0, t1, p0, p1);
        f.xmin = ax_t(t0);
        f.xmax = ax_t(t1);
        f.ymin = ax_p(p0);
        f.ymax = ax_p(p1);
        if (m_logp) f.ymin = std::max(1e-300, f.ymin);
        f.log_y = m_logp;
        f.log_x = m_logt;
        draw_plot_axes(canvas, f, theme, m_logt ? "Temperature (K, log scale)" : "Temperature (" + u.get("temperature") + ")", "Pressure (" + u.get("pressure") + (m_logp ? ", log scale)" : ")"));
        const std::string key = std::to_string(canvas.width()) + "x" + std::to_string(canvas.height()) + ":" + std::to_string(m_version_build) + u.get("temperature") + u.get("pressure") + (m_extreme ? "x" : "") + format_number(f.xmin, 8) + format_number(f.xmax, 8) + format_number(f.ymin, 8) + format_number(f.ymax, 8) + (m_logp ? "L" : "") + (m_logt ? "T" : "");
        const int step = 4;
        if (key != m_region_key) {
            m_region_key = key;
            m_cells.clear();
            m_label_pos.clear();
            std::map<int, std::pair<Vec2, int>> sums;
            for (double y = f.y0; y < f.y0 + f.h; y += step) {
                for (double x = f.x0; x < f.x0 + f.w; x += step) {
                    const double t = ax_t_inv(f.x_of(x + 0.5 * step));
                    const double p = ax_p_inv(f.y_of(y + 0.5 * step));
                    int code = static_cast<int>(fc::Phase::Unknown);
                    try {
                        code = static_cast<int>(m_diag.phase_at(t, std::max(p, 1e-300)));
                    } catch (...) {
                    }
                    if (m_extreme && (code == static_cast<int>(fc::Phase::Gas) || code == static_cast<int>(fc::Phase::Supercritical)) && ionization_fraction(t, p) > 0.5) code = PLASMA_CELL;
                    m_cells.push_back(code);
                    auto& sm = sums[code];
                    sm.first = sm.first + Vec2{ x, y };
                    sm.second += 1;
                }
            }
            for (const auto& kv : sums) {
                if (kv.second.second < 30) continue;
                m_label_pos[kv.first] = kv.second.first * (1.0 / kv.second.second);
            }
        }
        std::size_t idx = 0;
        for (double y = f.y0; y < f.y0 + f.h; y += step) {
            for (double x = f.x0; x < f.x0 + f.w; x += step) {
                if (idx >= m_cells.size()) break;
                Color c = cell_color(m_cells[idx++]);
                c.a = 70;
                canvas.fill_rect_px(x, y, std::min<double>(step, f.x0 + f.w - x), std::min<double>(step, f.y0 + f.h - y), c);
            }
        }
        for (const auto& kv : m_label_pos) {
            const char* w = cell_word(kv.first);
            const double half = 0.5 * canvas.text_width(w, 14.0);
            const double lx = std::max(f.x0 + half + 4.0, std::min(f.x0 + f.w - half - 4.0, kv.second.x));
            double ly = kv.second.y;
            for (const fc::PhasePoint* mp : { &m_diag.triple, &m_diag.critical, &m_diag.normal_melting, &m_diag.normal_boiling }) {
                if (!finite(mp->t) || !finite(mp->p) || mp->p <= 0.0) continue;
                const Vec2 q = screen_of(mp->t, mp->p);
                if (std::fabs(q.y - ly) < 22.0 && std::fabs(q.x - lx) < half + 70.0) ly = q.y + (ly >= q.y ? 26.0 : -26.0);
            }
            if (*w) canvas.text(Vec2{ lx, ly + 5.0 }, w, 14.0, mix(cell_color(kv.first), theme.text, 0.5), TextAlign::Center, 1.2);
        }
        const auto clip = [&](Vec2 q) { return q.y >= f.y0 - 2.0 && q.y <= f.y0 + f.h + 2.0 && q.x >= f.x0 - 2.0 && q.x <= f.x0 + f.w + 2.0; };
        const auto curve = [&](const std::vector<fc::PhasePoint>& pts, Color c, double width) {
            Vec2 prev{};
            bool have = false;
            for (const fc::PhasePoint& pt : pts) {
                if (!finite(pt.t) || !finite(pt.p) || pt.p <= 0.0) {
                    have = false;
                    continue;
                }
                const Vec2 q = screen_of(pt.t, pt.p);
                if (!clip(q)) {
                    have = false;
                    continue;
                }
                if (have) canvas.line_px(prev, q, c, width);
                prev = q;
                have = true;
            }
        };
        const auto hot = [&](int c) { return m_follow == c || (m_follow < 0 && !m_down && curve_at(m_mouse) == c); };
        curve(m_diag.sublimation, hot(0) ? theme.yellow : theme.text, hot(0) ? 3.4 : 2.4);
        curve(m_diag.vaporization, hot(1) ? theme.yellow : theme.text, hot(1) ? 3.4 : 2.4);
        curve(m_diag.fusion, hot(2) ? theme.yellow : theme.text, hot(2) ? 3.4 : 2.4);
        if (m_extreme) {
            std::vector<fc::PhasePoint> half, quantum;
            for (int i = 0; i <= 160; ++i) {
                const double pp = m_logp ? std::pow(10.0, std::log10(std::max(p0, 1e-12)) + (std::log10(p1) - std::log10(std::max(p0, 1e-12))) * i / 160.0) : p0 + (p1 - p0) * i / 160.0;
                double lo = 10.0, hi = 1e7;
                if (ionization_fraction(hi, pp) < 0.5) continue;
                for (int it = 0; it < 80; ++it) {
                    const double mid = std::sqrt(lo * hi);
                    if (ionization_fraction(mid, pp) > 0.5) hi = mid;
                    else lo = mid;
                }
                half.push_back(fc::PhasePoint{ hi, pp });
            }
            for (int i = 0; i <= 200; ++i) {
                const double t = m_logt ? std::pow(10.0, std::log10(t0) + (std::log10(t1) - std::log10(t0)) * i / 200.0) : t0 + (t1 - t0) * i / 200.0;
                if (t > 0) quantum.push_back(fc::PhasePoint{ t, quantum_pressure(t) });
            }
            curve(half, rgb(245, 140, 110), 1.6);
            Vec2 prev{};
            bool have = false;
            for (const fc::PhasePoint& pt : quantum) {
                const Vec2 q = screen_of(pt.t, pt.p);
                if (clip(q) && have) canvas.dashed_line_px(prev, q, theme.green, 1.6, 6.0);
                prev = q;
                have = clip(q);
            }
            for (const fc::PhasePoint& pt : quantum) {
                const Vec2 q = screen_of(pt.t, pt.p);
                if (clip(q) && q.x > f.x0 + 20.0) {
                    canvas.text(Vec2{ q.x + 4.0, q.y - 4.0 }, boson() ? "BEC line (ideal Bose gas)" : "Fermi degeneracy line", 10.0, theme.green);
                    break;
                }
            }
        }
        const double atm = ax_p(101325.0);
        if (atm > f.ymin && atm < f.ymax) {
            const Vec2 l = f.to_screen(f.xmin, atm), r = f.to_screen(f.xmax, atm);
            canvas.dashed_line_px(l, r, theme.muted, 1.2, 6.0);
            canvas.text(Vec2{ r.x - 4.0, r.y - 5.0 }, "1 atm", 10.0, theme.muted, TextAlign::Right);
        }
        std::vector<std::array<double, 4>> placed;
        const auto marker = [&](const char* key, const fc::PhasePoint& pt, const std::string& label, Color c, bool left) {
            if (!finite(pt.t) || !finite(pt.p) || pt.p <= 0.0) return;
            const Vec2 q = screen_of(pt.t, pt.p);
            if (!f.inside(q)) return;
            const bool sel = m_selection == key;
            const bool hover = !m_down && length(q - m_mouse) < 10.0;
            canvas.disc_px(q, sel || hover ? 7.0 : 5.0, c);
            canvas.circle_px(q, sel || hover ? 9.5 : 7.0, sel ? theme.text : theme.background, 1.5);
            const double tw = canvas.text_width(label, 11.0);
            const bool flip = left || q.x + 9.0 + tw > f.x0 + f.w;
            const double lx0 = flip ? q.x - 10.0 - tw : q.x + 10.0;
            double ly = q.y - 8.0;
            for (int tries = 0; tries < 6; ++tries) {
                bool hit = false;
                for (const auto& r : placed) {
                    if (lx0 < r[2] && lx0 + tw > r[0] && ly - 10.0 < r[3] && ly + 2.0 > r[1]) {
                        hit = true;
                        break;
                    }
                }
                if (!hit) break;
                ly -= 13.0;
            }
            placed.push_back({ lx0, ly - 10.0, lx0 + tw, ly + 2.0 });
            canvas.text(Vec2{ flip ? q.x - 10.0 : q.x + 10.0, ly }, label, 11.0, c, flip ? TextAlign::Right : TextAlign::Left);
        };
        marker("triple", m_diag.triple, "triple point", theme.yellow, false);
        marker("critical", m_diag.critical, "critical point", theme.red, false);
        marker("normal_melting", m_diag.normal_melting, "normal mp", theme.green, true);
        marker("normal_boiling", m_diag.normal_boiling, "normal bp", theme.warm, false);
        const Vec2 sp = screen_of(m_point_t, m_point_p);
        if (f.inside(sp)) {
            canvas.line_px(Vec2{ sp.x, f.y0 }, Vec2{ sp.x, f.y0 + f.h }, rgb(theme.accent.r, theme.accent.g, theme.accent.b, 90), 1.0);
            canvas.line_px(Vec2{ f.x0, sp.y }, Vec2{ f.x0 + f.w, sp.y }, rgb(theme.accent.r, theme.accent.g, theme.accent.b, 90), 1.0);
            canvas.circle_px(sp, 8.0, theme.accent, 2.4);
            canvas.disc_px(sp, 2.5, theme.accent);
        }
        std::vector<std::string> hud{ m_props.name.empty() ? m_props.formula : m_props.name };
        std::string state = "state point: " + u.format(m_point_t, "temperature", 4) + ", " + u.format(m_point_p, "pressure", 4) + " \xE2\x86\x92 " + fc::to_string(m_diag.phase_at(m_point_t, std::max(m_point_p, 1e-300)));
        if (m_follow >= 0) state += "  (following the line)";
        hud.push_back(state);
        draw_hud(canvas, theme, hud, false);
        if (f.inside(m_mouse) && !m_down) {
            const double t = ax_t_inv(f.x_of(m_mouse.x)), p = ax_p_inv(f.y_of(m_mouse.y));
            std::string line = u.format(t, "temperature", 4) + ", " + u.format(p, "pressure", 4) + ": " + fc::to_string(m_diag.phase_at(t, std::max(p, 1e-300)));
            if (m_extreme) line += ", " + format_number(100.0 * ionization_fraction(t, p), 3) + "% ionized";
            if (!marker_at(m_mouse).empty()) line = "click for details about this point";
            else if (curve_at(m_mouse) >= 0) line = "drag along this line to follow it";
            draw_hud(canvas, theme, { line }, true);
        }
    }

    void draw_3d(Canvas& canvas, const Theme& theme) {
        ensure_surface();
        const UnitPrefs& u = m_shared.units;
        if (!m_surface.ok || m_surface.nt < 2 || m_surface.nv < 2) {
            draw_hud(canvas, theme, { "3D P\xE2\x80\x93v\xE2\x80\x93T surface", "not enough data for this substance" }, false);
            return;
        }
        const fc::PvtSurface& s = m_surface;
        const double tmin = s.temperatures.front(), tmax = s.temperatures.back();
        const double lvmin = std::log10(s.specific_volumes.front()), lvmax = std::log10(s.specific_volumes.back());
        double lpmin = 1e300, lpmax = -1e300;
        for (double p : s.pressures) {
            if (p > 0.0 && finite(p)) {
                lpmin = std::min(lpmin, std::log10(p));
                lpmax = std::max(lpmax, std::log10(p));
            }
        }
        lpmin = std::max(lpmin, std::log10(std::max(1e-3, m_diag.p_min)));
        lpmax = std::min(lpmax, std::log10(std::max(m_diag.p_max, m_props.critical_pressure * 3.0)));
        if (!(lpmax > lpmin)) lpmax = lpmin + 1.0;
        const auto map = [&](double t, double v, double p) {
            const double x = (t - tmin) / (tmax - tmin) * 300.0 - 150.0;
            const double z = (std::log10(v) - lvmin) / (lvmax - lvmin) * 300.0 - 150.0;
            const double lp = std::max(lpmin, std::min(lpmax, std::log10(std::max(1e-300, p))));
            const double y = (lp - lpmin) / (lpmax - lpmin) * 220.0 - 110.0;
            return Vec3{ x, y, z };
        };
        if (m_cam.zoom <= 0.0) {
            m_cam.target = Vec3{};
            m_cam.pan = Vec2{};
            m_cam.zoom = 0.42 * std::min(canvas.width(), canvas.height()) / 230.0;
            m_cam.perspective = 0.0012;
        }
        struct Quad {
            Vec2 p[4];
            double depth;
            Color c;
        };
        m_probe_vertices.clear();
        for (int it = 0; it < s.nt; ++it) {
            for (int iv = 0; iv < s.nv; ++iv) {
                double depth = 0.0;
                const Vec2 sp = m_cam.project(map(s.temperatures[static_cast<std::size_t>(it)], s.specific_volumes[static_cast<std::size_t>(iv)], s.pressure(it, iv)), &depth);
                m_probe_vertices.push_back(ProbeVertex{ sp, depth, it, iv });
            }
        }
        std::vector<Quad> quads;
        quads.reserve(static_cast<std::size_t>(s.nt * s.nv));
        const Vec3 light = normalized(Vec3{ -0.4, 0.8, 0.5 });
        for (int it = 0; it + 1 < s.nt; ++it) {
            for (int iv = 0; iv + 1 < s.nv; ++iv) {
                const Vec3 a = map(s.temperatures[static_cast<std::size_t>(it)], s.specific_volumes[static_cast<std::size_t>(iv)], s.pressure(it, iv));
                const Vec3 b = map(s.temperatures[static_cast<std::size_t>(it + 1)], s.specific_volumes[static_cast<std::size_t>(iv)], s.pressure(it + 1, iv));
                const Vec3 c = map(s.temperatures[static_cast<std::size_t>(it + 1)], s.specific_volumes[static_cast<std::size_t>(iv + 1)], s.pressure(it + 1, iv + 1));
                const Vec3 d = map(s.temperatures[static_cast<std::size_t>(it)], s.specific_volumes[static_cast<std::size_t>(iv + 1)], s.pressure(it, iv + 1));
                Vec3 n = normalized(cross(c - a, d - b));
                if (n.y < 0) n = n * -1.0;
                const double shade = 0.45 + 0.55 * std::fabs(dot(n, light));
                Color col = region_color(s.region(it, iv));
                col = rgb(static_cast<int>(col.r * shade), static_cast<int>(col.g * shade), static_cast<int>(col.b * shade));
                Quad q;
                double da = 0, db = 0, dc = 0, dd = 0;
                q.p[0] = m_cam.project(a, &da);
                q.p[1] = m_cam.project(b, &db);
                q.p[2] = m_cam.project(c, &dc);
                q.p[3] = m_cam.project(d, &dd);
                q.depth = 0.25 * (da + db + dc + dd);
                q.c = col;
                quads.push_back(q);
            }
        }
        std::sort(quads.begin(), quads.end(), [](const Quad& x, const Quad& y) { return x.depth < y.depth; });
        Scene3D axes(m_cam);
        const Vec3 o{ -150.0, -110.0, -150.0 };
        axes.line(o, Vec3{ 150.0, -110.0, -150.0 }, 1.5, theme.axes, -1e6);
        axes.line(o, Vec3{ -150.0, 110.0, -150.0 }, 1.5, theme.axes, -1e6);
        axes.line(o, Vec3{ -150.0, -110.0, 150.0 }, 1.5, theme.axes, -1e6);
        axes.draw(canvas);
        for (const Quad& q : quads) {
            canvas.fill_polygon(std::vector<Vec2>{ q.p[0], q.p[1], q.p[2], q.p[3] }, q.c);
        }
        Scene3D over(m_cam);
        for (std::size_t i = 0; i + 1 < s.dome_t.size(); ++i) {
            over.line(map(s.dome_t[i], s.dome_liquid_v[i], s.dome_p[i]), map(s.dome_t[i + 1], s.dome_liquid_v[i + 1], s.dome_p[i + 1]), 2.0, theme.text, 1e6);
            over.line(map(s.dome_t[i], s.dome_vapor_v[i], s.dome_p[i]), map(s.dome_t[i + 1], s.dome_vapor_v[i + 1], s.dome_p[i + 1]), 2.0, theme.text, 1e6);
        }
        if (finite(s.triple_t) && finite(s.triple_solid_v) && finite(s.triple_vapor_v)) {
            over.line(map(s.triple_t, s.triple_solid_v, s.triple_p), map(s.triple_t, s.triple_vapor_v, s.triple_p), 2.5, theme.yellow, 1e6);
            over.text(map(s.triple_t, s.triple_vapor_v, s.triple_p), "triple line", 11.0, theme.yellow, 1e6, Vec2{ 0.0, -8.0 });
        }
        if (finite(s.critical_t) && finite(s.critical_v)) {
            over.disc(map(s.critical_t, s.critical_v, s.critical_p), 5.0, theme.red, 1e6);
            over.text(map(s.critical_t, s.critical_v, s.critical_p), "critical point", 11.0, theme.red, 1e6, Vec2{ 0.0, -10.0 });
        }
        if (m_probe_t >= 0 && m_probe_v >= 0 && m_probe_t < s.nt && m_probe_v < s.nv) {
            const Vec3 pp = map(s.temperatures[static_cast<std::size_t>(m_probe_t)], s.specific_volumes[static_cast<std::size_t>(m_probe_v)], s.pressure(m_probe_t, m_probe_v));
            over.disc(pp, 6.0, theme.accent, 3e6);
            over.ring(pp, 9.0, 2.0, theme.text, 3e6);
            over.text(pp, u.format(s.temperatures[static_cast<std::size_t>(m_probe_t)], "temperature", 4) + ", " + u.format(s.pressure(m_probe_t, m_probe_v), "pressure", 4) + ", " + fc::to_string(s.region(m_probe_t, m_probe_v)), 11.0, theme.text, 3e6, Vec2{ 0.0, -14.0 });
        }
        over.text(Vec3{ 165.0, -110.0, -150.0 }, "T (" + u.get("temperature") + ")", 12.0, theme.text, 2e6);
        over.text(Vec3{ -150.0, 125.0, -150.0 }, "log P", 12.0, theme.text, 2e6);
        over.text(Vec3{ -150.0, -110.0, 170.0 }, "log v (specific volume)", 12.0, theme.text, 2e6);
        over.text(Vec3{ -150.0, -128.0, -150.0 }, fixed(disp_t(tmin), 0), 10.0, theme.muted, 2e6);
        over.text(Vec3{ 150.0, -128.0, -150.0 }, fixed(disp_t(tmax), 0), 10.0, theme.muted, 2e6);
        over.draw(canvas);
        std::vector<std::string> hud{ (m_props.name.empty() ? m_props.formula : m_props.name) + "  \xC2\xB7  P\xE2\x80\x93v\xE2\x80\x93T surface (Peng\xE2\x80\x93Robinson)" };
        hud.push_back("v from " + u.format(s.specific_volumes.front(), "specific volume", 3) + " to " + u.format(s.specific_volumes.back(), "specific volume", 3));
        draw_hud(canvas, theme, hud, false);
        draw_hud(canvas, theme, { "drag on the surface to read points  \xC2\xB7  Shift+drag or right-drag to rotate  \xC2\xB7  wheel to zoom" }, true);
    }

    std::string readouts() const {
        std::ostringstream o;
        const UnitPrefs& u = m_shared.units;
        if (!m_diag.ok) {
            o << "#Phase diagram\n!" << (m_error.empty() ? "nothing loaded" : m_error) << "\n";
            return o.str();
        }
        const fc::PhaseProperties& p = m_props;
        const auto est = [&](const char* key) {
            if (m_overrides.count(key)) return std::string(" (yours)");
            return is_estimated(p, key) ? std::string(" (estimated)") : std::string();
        };
        o << "#" << (p.name.empty() ? p.formula : p.name) << "\n";
        if (!p.formula.empty()) o << "Formula\t" << pretty_formula(p.formula) << "\n";
        o << "Source\t" << m_base_source << "\n";
        if (finite(p.molar_mass)) o << "Molar mass\t" << u.format(p.molar_mass, "molar mass", 5) << est("molar_mass") << "\n";
        o << "#Key points\n";
        o << "Triple point\t" << u.format(p.triple_temperature, "temperature", 5) << ", " << u.format(p.triple_pressure, "pressure", 4) << est("triple_temperature") << "\n";
        o << "Critical point\t" << u.format(p.critical_temperature, "temperature", 5) << ", " << u.format(p.critical_pressure, "pressure", 4) << est("critical_temperature") << "\n";
        if (finite(p.melting_point)) o << "Normal melting point\t" << u.format(p.melting_point, "temperature", 5) << est("melting_point") << "\n";
        if (finite(p.boiling_point)) o << "Normal boiling point\t" << u.format(p.boiling_point, "temperature", 5) << est("boiling_point") << "\n";
        else if (p.triple_pressure > 101325.0) o << "Sublimes at 1 atm\t" << u.format(m_diag.sublimation_temperature(101325.0), "temperature", 5) << "\n";
        if (finite(p.vaporization_enthalpy)) o << "\xCE\x94H vaporization\t" << u.format(p.vaporization_enthalpy, "molar energy", 4) << est("vaporization_enthalpy") << "\n";
        if (finite(p.fusion_enthalpy)) o << "\xCE\x94H fusion\t" << u.format(p.fusion_enthalpy, "molar energy", 4) << est("fusion_enthalpy") << "\n";
        if (finite(p.acentric_factor)) o << "Acentric factor\t" << fixed(p.acentric_factor, 3) << est("acentric_factor") << "\n";
        if (finite(p.liquid_density)) o << "Liquid density\t" << u.format(p.liquid_density, "density", 4) << est("liquid_density") << "\n";
        if (finite(p.solid_density)) o << "Solid density\t" << u.format(p.solid_density, "density", 4) << est("solid_density") << "\n";
        const double dv = p.fusion_volume_change();
        if (finite(dv)) o << "Melting line\t" << (dv < 0 ? "leans left (solid is less dense than liquid, like ice)" : "leans right (solid is denser)") << "\n";
        o << "#State point\n";
        o << "Temperature\t" << u.format(m_point_t, "temperature", 5) << "\n";
        o << "Pressure\t" << u.format(m_point_p, "pressure", 5) << "\n";
        o << "Phase\t" << describe(m_point_t, m_point_p) << "\n";
        try {
            if (m_point_t > p.triple_temperature && m_point_t < p.critical_temperature) o << "Vapor pressure here\t" << u.format(m_diag.vapor_pressure(m_point_t), "pressure", 4) << "\n";
            else if (m_point_t <= p.triple_temperature) o << "Sublimation pressure\t" << u.format(m_diag.sublimation_pressure(m_point_t), "pressure", 4) << "\n";
        } catch (...) {
        }
        try {
            const fc::Phase ph = m_diag.phase_at(m_point_t, std::max(m_point_p, 1e-300));
            if (p.substance) {
                const fc::Phase use = ph == fc::Phase::Supercritical ? fc::Phase::Gas : ph;
                const double cp = p.substance->heat_capacity(use, m_point_t);
                if (finite(cp) && cp > 0) {
                    o << "Heat capacity C\xE2\x82\x9A\t" << u.format(cp, "molar entropy", 4) << " (" << fc::to_string(use) << ")\n";
                    if (p.molar_mass > 0) o << "Specific heat\t" << format_number(cp / p.molar_mass, 4) << " J/(g\xC2\xB7K)\n";
                }
            }
        } catch (...) {
        }
        o << "#Extreme states\n";
        const double ie = ionization_energy_ev();
        o << "Lowest ionization energy\t" << format_number(ie, 4) << " eV\n";
        double lo = 10.0, hi = 1e7;
        for (int it = 0; it < 80; ++it) {
            const double mid = std::sqrt(lo * hi);
            if (ionization_fraction(mid, 101325.0) > 0.5) hi = mid;
            else lo = mid;
        }
        o << "Half ionized at 1 atm\t" << u.format(hi, "temperature", 3) << " (Saha equation)\n";
        o << "Particle statistics\t" << (boson() ? "boson (even number of neutrons)" : "fermion (odd number of neutrons)") << "\n";
        const double hbar = H_PLANCK / (2.0 * 3.14159265358979);
        const double n_cold = 1e20;
        const double tq = boson() ? 2.0 * 3.14159265358979 * hbar * hbar / (particle_mass() * K_BOLTZ) * std::pow(n_cold / 2.612, 2.0 / 3.0) : hbar * hbar / (2.0 * particle_mass() * K_BOLTZ) * std::pow(3.0 * 3.14159265358979 * 3.14159265358979 * n_cold, 2.0 / 3.0);
        o << (boson() ? "BEC temperature" : "Fermi temperature") << " at 10\xC2\xB2\xE2\x81\xB0 atoms/m\xC2\xB3\t" << format_number(tq * 1e9, 4) << " nK\n";
        o << "Turn on Show plasma and quantum gases to shade where the gas is mostly ionized and draw the line below which an ideal gas of these particles would be quantum degenerate. Real substances freeze long before; only very dilute, laser-cooled gases reach it.\n";
        if (!p.estimated.empty() || !m_diag.notes.empty()) {
            o << "#Notes\n";
            if (!p.estimated.empty()) o << "Values marked estimated come from group contributions, corresponding states or the element data, so treat them as rough.\n";
            for (const std::string& n : m_diag.notes) o << n << "\n";
        }
        return o.str();
    }

    ChemShared& m_shared;
    fc::PhaseProperties m_base;
    std::string m_base_source;
    std::map<std::string, double> m_overrides;
    fc::PhaseProperties m_props;
    fc::PhaseDiagram m_diag;
    fc::PvtSurface m_surface;
    std::string m_error;
    bool m_3d = false;
    double m_point_t = 298.15;
    double m_point_p = 101325.0;
    PlotFrame m_frame;
    Camera3D m_cam;
    bool m_down = false;
    int m_button = 0;
    Vec2 m_last;
    Vec2 m_mouse{ -1.0, -1.0 };
    std::string m_region_key;
    std::vector<int> m_cells;
    std::map<int, Vec2> m_label_pos;
    std::uint64_t m_version = 1;
    std::uint64_t m_version_build = 0;
    unsigned m_mods = 0;
    int m_follow = -1;
    std::string m_selection;
    bool m_logp = true;
    bool m_logt = false;
    bool m_extreme = false;
    double m_tmin = fc::NO_VALUE, m_tmax = fc::NO_VALUE, m_pmin = fc::NO_VALUE, m_pmax = fc::NO_VALUE;
    struct ProbeVertex {
        Vec2 s;
        double depth;
        int it, iv;
    };
    std::vector<ProbeVertex> m_probe_vertices;
    int m_probe_t = -1, m_probe_v = -1;
};

} // namespace

std::unique_ptr<ChemView> make_phase_view(ChemShared& shared) {
    return std::unique_ptr<ChemView>(new PhaseView(shared));
}

} // namespace sapling

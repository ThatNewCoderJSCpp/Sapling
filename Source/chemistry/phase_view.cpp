#include "molecule_draw.hpp"
#include "view_util.hpp"

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

Color phase_color(fc::Phase ph) {
    switch (ph) {
        case fc::Phase::Solid: return rgb(116, 134, 170);
        case fc::Phase::Liquid: return rgb(76, 140, 220);
        case fc::Phase::Gas: return rgb(222, 196, 120);
        case fc::Phase::Supercritical: return rgb(186, 140, 230);
        default: return rgb(90, 90, 100);
    }
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

const char* phase_word(fc::Phase ph) {
    switch (ph) {
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

    bool mouse_down(Vec2 p, int button, unsigned) override {
        m_down = true;
        m_button = button;
        m_last = p;
        if (!m_3d && button == 1 && m_frame.inside(p)) {
            set_point_from_screen(p);
            return true;
        }
        return false;
    }

    bool mouse_move(Vec2 p, bool) override {
        const Vec2 d = p - m_last;
        m_last = p;
        m_mouse = p;
        if (m_3d) {
            if (m_down) {
                m_cam.orbit(d.x, d.y);
                return true;
            }
            return false;
        }
        if (m_down && m_button == 1 && m_frame.inside(p)) set_point_from_screen(p);
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
            return fixed(from_si_or(m_point_t, tu, fc::Dimension::Temperature), 2) + "\t" + format_number(from_si_or(m_point_p, pu, fc::Dimension::Pressure), 5);
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
                const bool est = std::find(m_props.estimated.begin(), m_props.estimated.end(), f.key) != m_props.estimated.end();
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
            p.estimated.erase(std::remove(p.estimated.begin(), p.estimated.end(), kv.first), p.estimated.end());
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

    void set_point_from_screen(Vec2 p) {
        m_point_t = si_t(m_frame.x_of(p.x));
        m_point_p = si_p(m_frame.y_of(p.y));
        ++m_version;
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
        f.xmin = disp_t(m_diag.t_min);
        f.xmax = disp_t(m_diag.t_max);
        f.ymin = std::max(1e-12, disp_p(m_diag.p_min));
        f.ymax = std::max(f.ymin * 10.0, disp_p(m_diag.p_max));
        f.log_y = true;
        f.log_x = false;
        draw_plot_axes(canvas, f, theme, "Temperature (" + u.get("temperature") + ")", "Pressure (" + u.get("pressure") + ", log scale)");
        const std::string key = std::to_string(canvas.width()) + "x" + std::to_string(canvas.height()) + ":" + std::to_string(m_version_build) + u.get("temperature") + u.get("pressure");
        const int step = 4;
        if (key != m_region_key) {
            m_region_key = key;
            m_cells.clear();
            m_label_pos.clear();
            std::map<int, std::pair<Vec2, int>> sums;
            for (double y = f.y0; y < f.y0 + f.h; y += step) {
                for (double x = f.x0; x < f.x0 + f.w; x += step) {
                    const double t = si_t(f.x_of(x + 0.5 * step));
                    const double p = si_p(f.y_of(y + 0.5 * step));
                    fc::Phase ph = fc::Phase::Unknown;
                    try {
                        ph = m_diag.phase_at(t, p);
                    } catch (...) {
                    }
                    m_cells.push_back(static_cast<int>(ph));
                    auto& s = sums[static_cast<int>(ph)];
                    s.first = s.first + Vec2{ x, y };
                    s.second += 1;
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
                Color c = phase_color(static_cast<fc::Phase>(m_cells[idx++]));
                c.a = 70;
                canvas.fill_rect_px(x, y, std::min<double>(step, f.x0 + f.w - x), std::min<double>(step, f.y0 + f.h - y), c);
            }
        }
        for (const auto& kv : m_label_pos) {
            const char* w = phase_word(static_cast<fc::Phase>(kv.first));
            const double half = 0.5 * canvas.text_width(w, 14.0);
            const double lx = std::max(f.x0 + half + 4.0, std::min(f.x0 + f.w - half - 4.0, kv.second.x));
            if (*w) canvas.text(Vec2{ lx, kv.second.y + 5.0 }, w, 14.0, mix(phase_color(static_cast<fc::Phase>(kv.first)), theme.text, 0.5), TextAlign::Center, 1.2);
        }
        const auto curve = [&](const std::vector<fc::PhasePoint>& pts, Color c) {
            std::vector<Vec2> s;
            for (const fc::PhasePoint& p : pts) {
                if (!finite(p.t) || !finite(p.p) || p.p <= 0.0) continue;
                const Vec2 q = f.to_screen(disp_t(p.t), disp_p(p.p));
                if (q.y < f.y0 - 2.0 || q.y > f.y0 + f.h + 2.0 || q.x < f.x0 - 2.0 || q.x > f.x0 + f.w + 2.0) {
                    if (s.size() > 1) {
                        for (std::size_t i = 0; i + 1 < s.size(); ++i) canvas.line_px(s[i], s[i + 1], c, 2.4);
                    }
                    s.clear();
                    continue;
                }
                s.push_back(q);
            }
            for (std::size_t i = 0; i + 1 < s.size(); ++i) canvas.line_px(s[i], s[i + 1], c, 2.4);
        };
        curve(m_diag.sublimation, theme.text);
        curve(m_diag.vaporization, theme.text);
        curve(m_diag.fusion, theme.text);
        const double atm = disp_p(101325.0);
        if (atm > f.ymin && atm < f.ymax) {
            const Vec2 l = f.to_screen(f.xmin, atm), r = f.to_screen(f.xmax, atm);
            canvas.dashed_line_px(l, r, theme.muted, 1.2, 6.0);
            canvas.text(Vec2{ r.x - 4.0, r.y - 5.0 }, "1 atm", 10.0, theme.muted, TextAlign::Right);
        }
        const auto marker = [&](const fc::PhasePoint& p, const std::string& label, Color c, bool left = false) {
            if (!finite(p.t) || !finite(p.p) || p.p <= 0.0) return;
            const Vec2 q = f.to_screen(disp_t(p.t), disp_p(p.p));
            if (!f.inside(q)) return;
            canvas.disc_px(q, 5.0, c);
            canvas.circle_px(q, 7.0, theme.background, 1.5);
            const double tw = canvas.text_width(label, 11.0);
            const bool flip = left || q.x + 9.0 + tw > f.x0 + f.w;
            canvas.text(Vec2{ flip ? q.x - 9.0 : q.x + 9.0, q.y - 7.0 }, label, 11.0, c, flip ? TextAlign::Right : TextAlign::Left);
        };
        marker(m_diag.triple, "triple point", theme.yellow);
        marker(m_diag.critical, "critical point", theme.red);
        marker(m_diag.normal_melting, "normal mp", theme.green, true);
        marker(m_diag.normal_boiling, "normal bp", theme.warm);
        const Vec2 sp = f.to_screen(disp_t(m_point_t), disp_p(m_point_p));
        if (f.inside(sp)) {
            canvas.line_px(Vec2{ sp.x, f.y0 }, Vec2{ sp.x, f.y0 + f.h }, rgb(theme.accent.r, theme.accent.g, theme.accent.b, 90), 1.0);
            canvas.line_px(Vec2{ f.x0, sp.y }, Vec2{ f.x0 + f.w, sp.y }, rgb(theme.accent.r, theme.accent.g, theme.accent.b, 90), 1.0);
            canvas.circle_px(sp, 8.0, theme.accent, 2.4);
            canvas.disc_px(sp, 2.5, theme.accent);
        }
        std::vector<std::string> hud{ m_props.name.empty() ? m_props.formula : m_props.name };
        hud.push_back("state point: " + u.format(m_point_t, "temperature", 4) + ", " + u.format(m_point_p, "pressure", 4) + " \xE2\x86\x92 " + fc::to_string(m_diag.phase_at(m_point_t, m_point_p)));
        draw_hud(canvas, theme, hud, false);
        if (f.inside(m_mouse) && !m_down) {
            const double t = si_t(f.x_of(m_mouse.x)), p = si_p(f.y_of(m_mouse.y));
            draw_hud(canvas, theme, { u.format(t, "temperature", 4) + ", " + u.format(p, "pressure", 4) + ": " + fc::to_string(m_diag.phase_at(t, p)) + "  (click to move the state point)" }, true);
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
        over.text(Vec3{ 165.0, -110.0, -150.0 }, "T (" + u.get("temperature") + ")", 12.0, theme.text, 2e6);
        over.text(Vec3{ -150.0, 125.0, -150.0 }, "log P", 12.0, theme.text, 2e6);
        over.text(Vec3{ -150.0, -110.0, 170.0 }, "log v (specific volume)", 12.0, theme.text, 2e6);
        over.text(Vec3{ -150.0, -128.0, -150.0 }, fixed(disp_t(tmin), 0), 10.0, theme.muted, 2e6);
        over.text(Vec3{ 150.0, -128.0, -150.0 }, fixed(disp_t(tmax), 0), 10.0, theme.muted, 2e6);
        over.draw(canvas);
        std::vector<std::string> hud{ (m_props.name.empty() ? m_props.formula : m_props.name) + "  \xC2\xB7  P\xE2\x80\x93v\xE2\x80\x93T surface (Peng\xE2\x80\x93Robinson)" };
        hud.push_back("v from " + u.format(s.specific_volumes.front(), "specific volume", 3) + " to " + u.format(s.specific_volumes.back(), "specific volume", 3));
        draw_hud(canvas, theme, hud, false);
        draw_hud(canvas, theme, { "solid  \xC2\xB7  liquid  \xC2\xB7  gas  \xC2\xB7  supercritical  \xC2\xB7  liquid+vapor dome  \xC2\xB7  drag to rotate" }, true);
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
            return std::find(p.estimated.begin(), p.estimated.end(), key) != p.estimated.end() ? std::string(" (estimated)") : std::string();
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
};

} // namespace

std::unique_ptr<ChemView> make_phase_view(ChemShared& shared) {
    return std::unique_ptr<ChemView>(new PhaseView(shared));
}

} // namespace sapling

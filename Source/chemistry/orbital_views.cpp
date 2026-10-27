#include "molecule_draw.hpp"
#include "view_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <random>
#include <sstream>

namespace sapling {

using namespace chem_util;

namespace {

constexpr double BOHR_PM = 52.917721;

Color l_color(int l, int n, int max_n) {
    static const Color base[] = { rgb(137, 180, 250), rgb(166, 227, 161), rgb(250, 179, 135), rgb(203, 166, 247), rgb(243, 139, 168) };
    const Color c = base[std::max(0, std::min(4, l))];
    const double dim = max_n > 1 ? 0.45 * (max_n - n) / (max_n - 1.0) : 0.0;
    return mix(c, rgb(60, 60, 70), dim);
}

fc::OrbitalModel model_of(const std::string& s) {
    return s == "slater" ? fc::OrbitalModel::Slater : fc::OrbitalModel::Hydrogenic;
}

class RadialView : public ChemView {
public:
    explicit RadialView(ChemShared& shared) : m_shared(shared) { rebuild(); }

    void draw(Canvas& canvas, const Theme& theme) override {
        canvas.clear(theme.background);
        if (m_curves.empty()) {
            draw_hud(canvas, theme, { "Radial probability", m_error.empty() ? std::string("pick an element") : m_error }, false);
            return;
        }
        const UnitPrefs& u = m_shared.units;
        const std::string lu = u.get("length");
        PlotFrame& f = m_frame;
        f.x0 = 70.0;
        f.y0 = 64.0;
        f.w = std::max(50.0, canvas.width() - f.x0 - 170.0);
        f.h = std::max(50.0, canvas.height() - f.y0 - 56.0);
        f.xmin = 0.0;
        f.xmax = from_si_or(m_rmax * 1e-12, lu, fc::Dimension::Length);
        double ymax = 0.0;
        for (const fc::RadialCurve& c : visible()) {
            for (double v : (m_cumulative ? c.cumulative : c.probability)) ymax = std::max(ymax, v);
        }
        if (m_total && !m_cumulative && m_show == "all") {
            for (double v : m_total_curve.probability) ymax = std::max(ymax, v);
        }
        if (m_cumulative) ymax = 1.0;
        f.ymin = 0.0;
        f.ymax = ymax > 0.0 ? ymax * 1.08 : 1.0;
        draw_plot_axes(canvas, f, theme, "distance from the nucleus r (" + lu + ")", m_cumulative ? "probability inside r" : "P(r) = r\xC2\xB2R\xC2\xB2 per pm");
        int max_n = 1;
        for (const fc::RadialCurve& c : m_curves) max_n = std::max(max_n, c.n);
        const auto plot = [&](const fc::RadialCurve& c, Color col, double width) {
            const std::vector<double>& ys = m_cumulative ? c.cumulative : c.probability;
            Vec2 prev{};
            for (std::size_t i = 0; i < c.r_pm.size() && i < ys.size(); ++i) {
                const Vec2 p = f.to_screen(from_si_or(c.r_pm[i] * 1e-12, lu, fc::Dimension::Length), ys[i]);
                if (i > 0) canvas.line_px(prev, p, col, width);
                prev = p;
            }
        };
        if (m_total && m_show == "all") plot(m_total_curve, theme.text, 2.6);
        for (const fc::RadialCurve& c : visible()) plot(c, l_color(c.l, c.n, max_n), 2.0);
        const std::vector<fc::RadialCurve> vis = visible();
        if (vis.size() == 1) {
            const fc::RadialCurve& c = vis.front();
            const auto vline = [&](double rpm, Color col, const std::string& label) {
                const Vec2 a = f.to_screen(from_si_or(rpm * 1e-12, lu, fc::Dimension::Length), f.ymin);
                const Vec2 b = f.to_screen(from_si_or(rpm * 1e-12, lu, fc::Dimension::Length), f.ymax);
                canvas.dashed_line_px(a, b, col, 1.2, 5.0);
                canvas.text(Vec2{ b.x + 4.0, b.y + 14.0 }, label, 10.0, col, TextAlign::Left);
            };
            vline(c.most_probable_pm, theme.yellow, "most probable");
            for (double r : c.node_radii_pm) vline(r, theme.red, "node");
        }
        double ly = f.y0 + 6.0;
        const double lx = f.x0 + f.w + 16.0;
        if (m_total && m_show == "all") {
            canvas.line_px(Vec2{ lx, ly + 6.0 }, Vec2{ lx + 18.0, ly + 6.0 }, theme.text, 2.6);
            canvas.text(Vec2{ lx + 24.0, ly + 10.0 }, "total", 11.0, theme.text);
            ly += 18.0;
        }
        for (const fc::RadialCurve& c : m_curves) {
            const bool shown = m_show == "all" || m_show == c.label;
            const Color col = shown ? l_color(c.l, c.n, max_n) : theme.grid;
            canvas.line_px(Vec2{ lx, ly + 6.0 }, Vec2{ lx + 18.0, ly + 6.0 }, col, 2.0);
            canvas.text(Vec2{ lx + 24.0, ly + 10.0 }, c.label + subscript_digits("") + "  " + std::to_string(c.occupancy) + " e\xE2\x81\xBB", 11.0, shown ? theme.text : theme.muted);
            ly += 18.0;
            if (ly > f.y0 + f.h) break;
        }
        const fc::Element& e = fc::element(m_z);
        draw_hud(canvas, theme, { std::string(e.name) + (m_charge ? " " + std::string(e.symbol) + superscript_charge(m_charge) : std::string()) + "  \xC2\xB7  " + (m_model == "slater" ? "Slater orbitals" : "hydrogen-like orbitals with Slater Zeff") }, false);
        if (f.inside(m_mouse)) {
            canvas.line_px(Vec2{ m_mouse.x, f.y0 }, Vec2{ m_mouse.x, f.y0 + f.h }, rgb(theme.accent.r, theme.accent.g, theme.accent.b, 120), 1.0);
            const double r_disp = f.x_of(m_mouse.x);
            const double rpm = to_si_or(r_disp, lu, fc::Dimension::Length, r_disp * 1e-12) * 1e12;
            std::string line = "r = " + u.format(rpm * 1e-12, "length", 4) + ":";
            int shown = 0;
            for (const fc::RadialCurve& c : vis) {
                if (c.r_pm.size() < 2) continue;
                const double step = c.r_pm[1] - c.r_pm[0];
                const std::size_t i = static_cast<std::size_t>(std::max(0.0, std::min(static_cast<double>(c.r_pm.size() - 1), std::round(rpm / step))));
                line += "  " + c.label + " " + format_number((m_cumulative ? c.cumulative : c.probability)[i], 3);
                if (++shown >= 5) break;
            }
            draw_hud(canvas, theme, { line }, true);
        }
    }

    bool mouse_move(Vec2 p, bool) override {
        m_mouse = p;
        return true;
    }

    std::string command(const std::vector<std::string>& a) override {
        const std::string c = arg(a, 0);
        if (c == "element") {
            if (const fc::Element* e = element_of(arg(a, 1))) m_z = e->atomic_number;
            m_auto_r = true;
            m_show = "all";
            rebuild();
            return fc::element(m_z).symbol;
        }
        if (c == "charge") {
            m_charge = std::max(-3, std::min(m_z - 1, static_cast<int>(num(arg(a, 1), 0.0))));
            m_auto_r = true;
            rebuild();
            return std::to_string(m_charge);
        }
        if (c == "model") {
            m_model = arg(a, 1) == "slater" ? "slater" : "hydrogenic";
            rebuild();
            return m_model;
        }
        if (c == "show") {
            m_show = arg(a, 1).empty() ? "all" : arg(a, 1);
            if (m_auto_r) rebuild();
            return m_show;
        }
        if (c == "total") {
            m_total = arg(a, 1) == "1";
            return arg(a, 1);
        }
        if (c == "cumulative") {
            m_cumulative = arg(a, 1) == "1";
            return arg(a, 1);
        }
        if (c == "rmax") {
            if (arg(a, 1) == "auto" || arg(a, 1).empty()) {
                m_auto_r = true;
            } else {
                const double v = to_si_or(num(arg(a, 1), 0.0), arg(a, 2).empty() ? "pm" : arg(a, 2), fc::Dimension::Length, 0.0) * 1e12;
                if (v > 1.0) {
                    m_rmax = v;
                    m_auto_r = false;
                }
            }
            rebuild();
            return format_number(m_rmax, 5);
        }
        if (c == "subshells") {
            std::string out = "all\tall subshells\n";
            for (const fc::RadialCurve& r : m_curves) out += r.label + "\t" + r.label + " (" + std::to_string(r.occupancy) + " e\xE2\x81\xBB)\n";
            return out;
        }
        if (c == "readouts") return readouts();
        return "error\tunknown command " + c;
    }

private:
    std::vector<fc::RadialCurve> visible() const {
        std::vector<fc::RadialCurve> out;
        for (const fc::RadialCurve& c : m_curves) {
            if (m_show == "all" || m_show == c.label) out.push_back(c);
        }
        return out;
    }

    void rebuild() {
        m_error.clear();
        m_curves.clear();
        const fc::Element& e = fc::element(m_z);
        try {
            if (m_auto_r) m_rmax = fc::suggested_radial_range_pm(e, m_charge);
            m_curves = fc::atom_radial_distributions(e, m_rmax, 500, model_of(m_model), m_charge);
            if (m_auto_r && m_show != "all") {
                const auto find = [&]() -> const fc::RadialCurve* {
                    for (const fc::RadialCurve& c : m_curves) {
                        if (c.label == m_show && !c.r_pm.empty()) return &c;
                    }
                    return nullptr;
                };
                const fc::RadialCurve* c = find();
                if (c) {
                    double guess = std::max(c->mean_pm, c->most_probable_pm) * 4.0;
                    if (!(guess > 0.5)) guess = m_rmax;
                    m_rmax = guess;
                    m_curves = fc::atom_radial_distributions(e, m_rmax, 500, model_of(m_model), m_charge);
                    c = find();
                    if (c && !c->cumulative.empty() && c->cumulative.back() > 0.0) {
                        const double full = c->cumulative.back();
                        double reach = c->r_pm.back();
                        for (std::size_t i = 0; i < c->cumulative.size() && i < c->r_pm.size(); ++i) {
                            if (c->cumulative[i] >= 0.9995 * full) {
                                reach = c->r_pm[i];
                                break;
                            }
                        }
                        m_rmax = std::max(reach * 1.15, c->most_probable_pm * 2.5);
                        m_curves = fc::atom_radial_distributions(e, m_rmax, 500, model_of(m_model), m_charge);
                    }
                }
            }
            m_total_curve = fc::total_radial_density(e, m_rmax, 500, model_of(m_model), m_charge);
        } catch (const std::exception& ex) {
            m_error = ex.what();
        }
    }

    std::string readouts() const {
        std::ostringstream o;
        const UnitPrefs& u = m_shared.units;
        const fc::Element& e = fc::element(m_z);
        o << "#" << e.name << (m_charge ? " " + std::string(e.symbol) + superscript_charge(m_charge) : std::string()) << "\n";
        o << "Model\t" << (m_model == "slater" ? "Slater-type orbitals" : "hydrogen-like with Slater's rules") << "\n";
        o << "P(r) is the chance of finding the electron in a thin shell at distance r. The peaks are where each subshell's electrons usually are; nodes are radii where P(r) drops to zero.\n";
        for (const fc::RadialCurve& c : m_curves) {
            o << "#" << c.label << " (" << c.occupancy << " e\xE2\x81\xBB)\n";
            o << "Zeff\t" << fixed(c.z_eff, 2) << "\n";
            o << "Most probable r\t" << u.format(c.most_probable_pm * 1e-12, "length", 4) << "\n";
            o << "Average r\t" << u.format(c.mean_pm * 1e-12, "length", 4) << "\n";
            o << "90% inside\t" << u.format(c.radius_90_pm * 1e-12, "length", 4) << "\n";
            o << "Radial nodes\t" << c.node_radii_pm.size() << "\n";
            if (finite(c.energy_ev)) o << "Orbital energy\t" << fixed(c.energy_ev, 2) << " eV\n";
        }
        return o.str();
    }

    ChemShared& m_shared;
    int m_z = 11;
    int m_charge = 0;
    std::string m_model = "hydrogenic";
    std::string m_show = "all";
    bool m_total = true;
    bool m_cumulative = false;
    bool m_auto_r = true;
    double m_rmax = 400.0;
    std::vector<fc::RadialCurve> m_curves;
    fc::RadialCurve m_total_curve;
    std::string m_error;
    PlotFrame m_frame;
    Vec2 m_mouse{ -1.0, -1.0 };
};

class CloudView : public ChemView {
public:
    explicit CloudView(ChemShared& shared) : m_shared(shared) {
        m_cam.yaw = 0.6;
        m_cam.pitch = 0.4;
        m_cam.perspective = 0.0;
        rebuild();
    }

    void draw(Canvas& canvas, const Theme& theme) override {
        canvas.clear(theme.background);
        m_cam.width = canvas.width();
        m_cam.height = canvas.height();
        if (!m_error.empty()) {
            draw_hud(canvas, theme, { "Electron cloud", m_error }, false);
            return;
        }
        if (m_3d) draw_3d(canvas, theme);
        else draw_2d(canvas, theme);
        std::string caption;
        if (m_3d) {
            if (m_style3 == "balloon") caption = m_target == "orbital" ? "surface holding 90% of the probability; colors are the signs of \xCF\x88" : "one sphere per subshell, holding 90% of its electrons";
            else if (m_style3 == "balloon_dots") caption = "90% surface with " + std::to_string(m_points.size()) + " sampled positions inside";
            else caption = std::to_string(m_points.size()) + " sampled positions; denser = more likely";
        } else {
            static const std::map<std::string, std::string> words = { { "heat", "brightness = probability density" }, { "gray", "brightness = probability density" }, { "contour", "lines of equal probability density" }, { "dots", "dots placed with probability |\xCF\x88|\xC2\xB2" }, { "outline", "the region holding 90% of the probability in this plane" } };
            caption = std::string("slice through the nucleus (") + m_plane + " plane), " + words.at(m_style2);
        }
        draw_hud(canvas, theme, { title(), caption }, false);
    }

    bool mouse_down(Vec2 p, int, unsigned) override {
        m_down = true;
        m_last = p;
        return false;
    }

    bool mouse_move(Vec2 p, bool) override {
        const Vec2 d = p - m_last;
        m_last = p;
        if (m_down && m_3d) {
            m_cam.orbit(d.x, d.y);
            return true;
        }
        return false;
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

    bool step(double dt) override {
        if (m_spin && m_3d && !m_down) {
            m_cam.yaw += dt * 0.4;
            return true;
        }
        return false;
    }

    std::string command(const std::vector<std::string>& a) override {
        const std::string c = arg(a, 0);
        if (c == "element") {
            if (const fc::Element* e = element_of(arg(a, 1))) m_z = e->atomic_number;
            m_charge = 0;
            m_target = "atom";
            rebuild();
            return fc::element(m_z).symbol;
        }
        if (c == "charge") {
            m_charge = std::max(-3, std::min(m_z - 1, static_cast<int>(num(arg(a, 1), 0.0))));
            rebuild();
            return std::to_string(m_charge);
        }
        if (c == "orbital") {
            const std::string t = arg(a, 1);
            if (t == "atom" || t == "valence") {
                m_target = t;
            } else {
                const auto parts = split(t, ',');
                if (parts.size() >= 3) {
                    m_target = "orbital";
                    m_n = std::max(1, static_cast<int>(num(parts[0], 1.0)));
                    m_l = std::max(0, std::min(m_n - 1, static_cast<int>(num(parts[1], 0.0))));
                    m_m = std::max(-m_l, std::min(m_l, static_cast<int>(num(parts[2], 0.0))));
                }
            }
            rebuild();
            return title();
        }
        if (c == "orbitals") {
            std::string out = "atom\twhole atom (all electrons)\nvalence\tvalence electrons only\n";
            try {
                const fc::ElectronConfiguration cfg = m_charge ? fc::ElectronConfiguration::of_ion(fc::element(m_z), m_charge) : fc::ElectronConfiguration::ground_state(fc::element(m_z));
                for (const fc::Subshell& s : cfg.subshells()) {
                    for (int m = -static_cast<int>(s.l); m <= static_cast<int>(s.l); ++m) {
                        out += std::to_string(s.n) + "," + std::to_string(s.l) + "," + std::to_string(m) + "\t" + fc::orbital_name(s.n, s.l, m) + " (" + fc::orbital_shape_name(s.l, m) + ")\n";
                    }
                }
                const int nx = cfg.highest_shell() + 1;
                for (int l = 0; l < std::min(nx, 4); ++l) {
                    for (int m = -l; m <= l; ++m) out += std::to_string(nx) + "," + std::to_string(l) + "," + std::to_string(m) + "\t" + fc::orbital_name(nx, l, m) + " (empty, excited)\n";
                }
            } catch (...) {
            }
            return out;
        }
        if (c == "view") {
            m_3d = arg(a, 1) == "3d";
            m_framed = false;
            return m_3d ? "3d" : "2d";
        }
        if (c == "plane") {
            const std::string p = arg(a, 1);
            if (p == "xy" || p == "xz" || p == "yz") m_plane = p;
            rebuild();
            return m_plane;
        }
        if (c == "signed") {
            m_signed = arg(a, 1) == "1";
            rebuild();
            return arg(a, 1);
        }
        if (c == "points") {
            m_count = std::max(200, std::min(60000, static_cast<int>(num(arg(a, 1), 8000.0))));
            rebuild();
            return std::to_string(m_count);
        }
        if (c == "model") {
            m_model = arg(a, 1) == "slater" ? "slater" : "hydrogenic";
            rebuild();
            return m_model;
        }
        if (c == "spin") {
            m_spin = arg(a, 1) == "1";
            return arg(a, 1);
        }
        if (c == "style3") {
            const std::string v = arg(a, 1);
            if (v == "sign" || v == "mono" || v == "type" || v == "balloon" || v == "balloon_dots") m_style3 = v;
            m_cache_key.clear();
            return "ok";
        }
        if (c == "style2") {
            const std::string v = arg(a, 1);
            if (v == "heat" || v == "gray" || v == "contour" || v == "dots" || v == "outline") m_style2 = v;
            return "ok";
        }
        if (c == "readouts") return readouts();
        return "error\tunknown command " + c;
    }

private:
    std::string title() const {
        const fc::Element& e = fc::element(m_z);
        std::string base = std::string(e.symbol) + superscript_charge(m_charge);
        if (m_target == "atom") return base + " whole atom";
        if (m_target == "valence") return base + " valence electrons";
        return base + " " + fc::orbital_name(m_n, m_l, m_m) + " orbital";
    }

    double zeff() const {
        try {
            return fc::effective_nuclear_charge(fc::element(m_z), m_n, m_l, m_charge);
        } catch (...) {
            return 1.0;
        }
    }

    void rebuild() {
        m_error.clear();
        m_points.clear();
        m_grid_dirty = true;
        m_cache_key.clear();
        m_shells.clear();
        m_slice = fc::DensitySlice{};
        const fc::Element& e = fc::element(m_z);
        try {
            if (m_target == "orbital") {
                const double z = zeff();
                m_extent = std::max(40.0, fc::orbital_extent_bohr(m_n, m_l, z, 0.985, model_of(m_model)) * BOHR_PM * 1.05);
                m_points = fc::sample_orbital(m_n, m_l, m_m, z, m_count, 7u, model_of(m_model));
                m_slice = fc::orbital_slice(m_n, m_l, m_m, z, m_plane[0], m_plane[1], m_extent, 200, m_signed, model_of(m_model));
            } else {
                m_extent = fc::suggested_radial_range_pm(e, m_charge) * 0.7;
                m_points = fc::sample_atom(e, m_count, 7u, model_of(m_model), m_charge, m_target == "valence");
                m_slice = fc::atom_slice(e, m_extent, 200, model_of(m_model), m_charge);
            }
        } catch (const std::exception& ex) {
            m_error = ex.what();
        }
        m_framed = false;
    }

    void draw_2d(Canvas& canvas, const Theme& theme) {
        if (m_slice.width <= 0 || m_slice.values.empty()) return;
        const double side = std::max(60.0, std::min(canvas.width(), canvas.height()) - 80.0);
        const double x0 = 0.5 * (canvas.width() - side), y0 = 0.5 * (canvas.height() - side) + 10.0;
        const int W = m_slice.width, H = m_slice.height;
        const double maxabs = std::max(std::fabs(m_slice.max_value), std::fabs(m_slice.min_value));
        const bool signed_values = m_signed && m_target == "orbital";
        const auto at = [&](int a, int b) { return m_slice.values[static_cast<std::size_t>(b) * W + a]; };
        const auto density = [&](double v) { return signed_values ? v * v : std::max(0.0, v); };
        double max_density = 0.0;
        for (double v : m_slice.values) max_density = std::max(max_density, density(v));
        const int px = static_cast<int>(side);
        const bool field = m_style2 == "heat" || m_style2 == "gray" || m_style2 == "contour";
        if (field) {
            for (int j = 0; j < px; ++j) {
                const double fy = (j + 0.5) / side * H - 0.5;
                const int j0 = std::max(0, std::min(H - 1, static_cast<int>(std::floor(fy))));
                const int j1 = std::min(H - 1, j0 + 1);
                const double ty = std::max(0.0, std::min(1.0, fy - j0));
                for (int i = 0; i < px; ++i) {
                    const double fx = (i + 0.5) / side * W - 0.5;
                    const int i0 = std::max(0, std::min(W - 1, static_cast<int>(std::floor(fx))));
                    const int i1 = std::min(W - 1, i0 + 1);
                    const double tx = std::max(0.0, std::min(1.0, fx - i0));
                    const double v = (at(i0, j0) * (1 - tx) + at(i1, j0) * tx) * (1 - ty) + (at(i0, j1) * (1 - tx) + at(i1, j1) * tx) * ty;
                    Color c;
                    if (m_style2 == "gray") {
                        const double t = max_density > 0 ? std::pow(density(v) / max_density, 0.45) : 0.0;
                        const auto g = static_cast<std::uint8_t>(std::lround(255.0 * t));
                        c = rgb(g, g, g);
                    } else if (signed_values) {
                        c = signed_color(maxabs > 0 ? v / maxabs : 0.0, theme);
                    } else {
                        c = heat_color(max_density > 0 ? std::pow(density(v) / max_density, 0.45) : 0.0, theme);
                    }
                    if (m_style2 == "contour") c.a = 90;
                    canvas.put_pixel(static_cast<int>(x0) + i, static_cast<int>(y0) + j, c);
                }
            }
        }
        const auto to_screen = [&](double gx, double gy) { return Vec2{ x0 + (gx + 0.5) * side / W, y0 + (gy + 0.5) * side / H }; };
        const auto contour = [&](const std::function<double(int, int)>& f, double level, Color color, double width) {
            for (int j = 0; j + 1 < H; ++j) {
                for (int i = 0; i + 1 < W; ++i) {
                    const double v[4] = { f(i, j) - level, f(i + 1, j) - level, f(i + 1, j + 1) - level, f(i, j + 1) - level };
                    const double cx[4] = { static_cast<double>(i), static_cast<double>(i + 1), static_cast<double>(i + 1), static_cast<double>(i) };
                    const double cy[4] = { static_cast<double>(j), static_cast<double>(j), static_cast<double>(j + 1), static_cast<double>(j + 1) };
                    Vec2 hits[4];
                    int count = 0;
                    for (int e = 0; e < 4; ++e) {
                        const int f2 = (e + 1) % 4;
                        if ((v[e] < 0) == (v[f2] < 0)) continue;
                        const double t = v[e] / (v[e] - v[f2]);
                        hits[count++] = to_screen(cx[e] + (cx[f2] - cx[e]) * t, cy[e] + (cy[f2] - cy[e]) * t);
                    }
                    if (count >= 2) canvas.line_px(hits[0], hits[1], color, width);
                    if (count == 4) canvas.line_px(hits[2], hits[3], color, width);
                }
            }
        };
        std::vector<double> sorted;
        double level90 = 0.0;
        if (m_style2 == "outline" || m_style2 == "contour") {
            for (double v : m_slice.values) sorted.push_back(density(v));
            std::sort(sorted.begin(), sorted.end(), std::greater<double>());
            double total = 0.0;
            for (double d : sorted) total += d;
            double run = 0.0;
            for (double d : sorted) {
                run += d;
                level90 = d;
                if (run >= 0.9 * total) break;
            }
        }
        if (m_style2 == "contour" && max_density > 0) {
            for (double frac : { 0.7, 0.4, 0.2, 0.1, 0.05, 0.02 }) {
                const double level = frac * max_density;
                if (signed_values) {
                    const double amp = std::sqrt(level);
                    contour([&](int a, int b) { return at(a, b); }, amp, rgb(243, 139, 168), 1.2);
                    contour([&](int a, int b) { return -at(a, b); }, amp, rgb(137, 180, 250), 1.2);
                } else {
                    contour([&](int a, int b) { return density(at(a, b)); }, level, theme.text, 1.0);
                }
            }
        }
        if (m_style2 == "outline" && level90 > 0) {
            for (int j = 0; j < px; ++j) {
                const int gj = std::max(0, std::min(H - 1, static_cast<int>((j + 0.5) / side * H)));
                for (int i = 0; i < px; ++i) {
                    const int gi = std::max(0, std::min(W - 1, static_cast<int>((i + 0.5) / side * W)));
                    const double v = at(gi, gj);
                    if (density(v) < level90) continue;
                    Color c = signed_values ? (v >= 0 ? rgb(243, 139, 168, 120) : rgb(137, 180, 250, 120)) : rgb(166, 227, 161, 110);
                    canvas.put_pixel(static_cast<int>(x0) + i, static_cast<int>(y0) + j, c);
                }
            }
            if (signed_values) {
                contour([&](int a, int b) { return at(a, b); }, std::sqrt(level90), rgb(243, 139, 168), 2.0);
                contour([&](int a, int b) { return -at(a, b); }, std::sqrt(level90), rgb(137, 180, 250), 2.0);
            } else {
                contour([&](int a, int b) { return density(at(a, b)); }, level90, rgb(166, 227, 161), 2.0);
            }
        }
        if (m_style2 == "dots" && max_density > 0) {
            std::mt19937 rng(11u);
            std::uniform_real_distribution<double> u01(0.0, 1.0);
            const int wanted = std::max(500, m_count / 2);
            int placed = 0;
            for (int tries = 0; tries < wanted * 400 && placed < wanted; ++tries) {
                const double gx = u01(rng) * (W - 1), gy = u01(rng) * (H - 1);
                const int gi = static_cast<int>(gx), gj = static_cast<int>(gy);
                const double v = at(gi, gj);
                if (u01(rng) * max_density > density(v)) continue;
                ++placed;
                const Color c = signed_values ? (v >= 0 ? rgb(243, 139, 168, 200) : rgb(137, 180, 250, 200)) : theme.text;
                canvas.disc_px(to_screen(gx, gy), 1.1, c);
            }
        }
        canvas.stroke_rect_px(x0, y0, side, side, theme.axes, 1.0);
        const UnitPrefs& u = m_shared.units;
        canvas.text(Vec2{ x0 + side * 0.5, y0 + side + 18.0 }, std::string(1, m_plane[0]) + " axis, \xC2\xB1" + u.format(m_extent * 1e-12, "length", 3), 11.0, theme.muted, TextAlign::Center);
        canvas.text(Vec2{ x0 - 8.0, y0 + side * 0.5 }, std::string(1, m_plane[1]), 12.0, theme.muted, TextAlign::Right);
        canvas.disc_px(Vec2{ x0 + side * 0.5, y0 + side * 0.5 }, 2.0, theme.red);
    }

    void build_grid() {
        m_grid_dirty = false;
        m_grid.clear();
        m_grid_n = 0;
        m_iso = 0.0;
        m_cache_key.clear();
        if (m_target != "orbital") {
            try {
                m_shells = fc::atom_radial_distributions(fc::element(m_z), std::max(60.0, m_extent * 2.5), 400, model_of(m_model), m_charge);
            } catch (...) {
                m_shells.clear();
            }
            return;
        }
        const int G = 56;
        const double E = m_extent;
        const double z = zeff();
        m_grid.assign(static_cast<std::size_t>(G) * G * G, 0.0f);
        std::vector<double> squares;
        squares.reserve(m_grid.size());
        for (int k = 0; k < G; ++k) {
            for (int j = 0; j < G; ++j) {
                for (int i = 0; i < G; ++i) {
                    const double x = -E + 2.0 * E * i / (G - 1), y = -E + 2.0 * E * j / (G - 1), w = -E + 2.0 * E * k / (G - 1);
                    double psi = 0.0;
                    try {
                        psi = fc::wavefunction(m_n, m_l, m_m, z, x / BOHR_PM, y / BOHR_PM, w / BOHR_PM, model_of(m_model));
                    } catch (...) {
                    }
                    m_grid[(static_cast<std::size_t>(k) * G + j) * G + i] = static_cast<float>(psi);
                    squares.push_back(psi * psi);
                }
            }
        }
        m_grid_n = G;
        std::sort(squares.begin(), squares.end(), std::greater<double>());
        double total = 0.0;
        for (double q : squares) total += q;
        double run = 0.0;
        for (double q : squares) {
            run += q;
            m_iso = std::sqrt(q);
            if (run >= 0.9 * total) break;
        }
    }

    double grid_at(Vec3 p) const {
        const int G = m_grid_n;
        if (G < 2) return 0.0;
        const double E = m_extent;
        const double fx = (p.x + E) / (2.0 * E) * (G - 1), fy = (p.y + E) / (2.0 * E) * (G - 1), fz = (p.z + E) / (2.0 * E) * (G - 1);
        if (fx < 0 || fy < 0 || fz < 0 || fx > G - 1 || fy > G - 1 || fz > G - 1) return 0.0;
        const int i = std::min(G - 2, static_cast<int>(fx)), j = std::min(G - 2, static_cast<int>(fy)), k = std::min(G - 2, static_cast<int>(fz));
        const double tx = fx - i, ty = fy - j, tz = fz - k;
        const auto g = [&](int a, int b, int c) { return static_cast<double>(m_grid[(static_cast<std::size_t>(c) * G + b) * G + a]); };
        const double c00 = g(i, j, k) * (1 - tx) + g(i + 1, j, k) * tx;
        const double c10 = g(i, j + 1, k) * (1 - tx) + g(i + 1, j + 1, k) * tx;
        const double c01 = g(i, j, k + 1) * (1 - tx) + g(i + 1, j, k + 1) * tx;
        const double c11 = g(i, j + 1, k + 1) * (1 - tx) + g(i + 1, j + 1, k + 1) * tx;
        return (c00 * (1 - ty) + c10 * ty) * (1 - tz) + (c01 * (1 - ty) + c11 * ty) * tz;
    }

    void draw_balloon(Canvas& canvas, bool see_through) {
        if (m_grid_dirty) build_grid();
        if (m_grid_n < 2 || !(m_iso > 0.0)) return;
        const int W = canvas.width(), H = canvas.height();
        const int step = m_down ? 4 : 2;
        const int bw = (W + step - 1) / step, bh = (H + step - 1) / step;
        char key[256];
        std::snprintf(key, sizeof(key), "%.5f|%.5f|%.6f|%.2f|%.2f|%d|%d|%d|%d", m_cam.yaw, m_cam.pitch, m_cam.zoom, m_cam.pan.x, m_cam.pan.y, W, H, step, see_through ? 1 : 0);
        if (m_cache_key != key) {
            m_cache_key = key;
            m_cache.assign(static_cast<std::size_t>(bw) * bh, Color{ 0, 0, 0, 0 });
            const double E = m_extent;
            const double R = E * 1.8;
            const double spacing = 2.0 * E / (m_grid_n - 1);
            const double dz = spacing * 0.6;
            const Vec3 light = normalized(Vec3{ -0.45, 0.55, 0.7 });
            const Vec3 half = normalized(light + Vec3{ 0.0, 0.0, 1.0 });
            const Vec3 origin_view = m_cam.to_view(Vec3{});
            for (int by = 0; by < bh; ++by) {
                for (int bx = 0; bx < bw; ++bx) {
                    const Vec2 sp{ bx * step + 0.5 * step, by * step + 0.5 * step };
                    Vec3 a = m_cam.unproject(sp, R);
                    Vec3 b = m_cam.unproject(sp, -R);
                    double t0 = 0.0, t1 = 1.0;
                    const Vec3 d = b - a;
                    bool miss = false;
                    for (int axis = 0; axis < 3 && !miss; ++axis) {
                        const double o = axis == 0 ? a.x : axis == 1 ? a.y : a.z;
                        const double v = axis == 0 ? d.x : axis == 1 ? d.y : d.z;
                        if (std::fabs(v) < 1e-12) {
                            if (o < -E || o > E) miss = true;
                            continue;
                        }
                        double ta = (-E - o) / v, tb = (E - o) / v;
                        if (ta > tb) std::swap(ta, tb);
                        t0 = std::max(t0, ta);
                        t1 = std::min(t1, tb);
                        if (t0 > t1) miss = true;
                    }
                    if (miss) continue;
                    const double len = length(d) * (t1 - t0);
                    const int steps = std::max(2, static_cast<int>(len / dz));
                    Vec3 prev = a + d * t0;
                    double prev_v = grid_at(prev);
                    for (int s = 1; s <= steps; ++s) {
                        const Vec3 p = a + d * (t0 + (t1 - t0) * s / steps);
                        const double v = grid_at(p);
                        if (std::fabs(v) < m_iso) {
                            prev = p;
                            prev_v = v;
                            continue;
                        }
                        Vec3 lo = prev, hi = p;
                        for (int it = 0; it < 6; ++it) {
                            const Vec3 mid = (lo + hi) * 0.5;
                            if (std::fabs(grid_at(mid)) >= m_iso) hi = mid;
                            else lo = mid;
                        }
                        const Vec3 h = hi;
                        const double e = spacing * 0.5;
                        Vec3 grad{
                            std::fabs(grid_at(h + Vec3{ e, 0, 0 })) - std::fabs(grid_at(h - Vec3{ e, 0, 0 })),
                            std::fabs(grid_at(h + Vec3{ 0, e, 0 })) - std::fabs(grid_at(h - Vec3{ 0, e, 0 })),
                            std::fabs(grid_at(h + Vec3{ 0, 0, e })) - std::fabs(grid_at(h - Vec3{ 0, 0, e }))
                        };
                        Vec3 normal = normalized(grad * -1.0);
                        const Vec3 nv = normalized(m_cam.to_view(normal) - origin_view);
                        const double lambert = std::max(0.0, dot(nv, light));
                        const double spec = std::pow(std::max(0.0, dot(nv, half)), 32.0) * 0.45;
                        const bool positive = grid_at(h) >= 0.0;
                        const double br = positive ? 243.0 : 137.0, bg = positive ? 139.0 : 180.0, bb = positive ? 168.0 : 250.0;
                        const double shade = 0.28 + 0.72 * lambert;
                        const auto channel = [&](double base) { return static_cast<std::uint8_t>(std::max(0.0, std::min(255.0, base * shade + 255.0 * spec))); };
                        m_cache[static_cast<std::size_t>(by) * bw + bx] = Color{ channel(br), channel(bg), channel(bb), static_cast<std::uint8_t>(see_through ? 150 : 255) };
                        (void)prev_v;
                        break;
                    }
                }
            }
        }
        for (int by = 0; by < bh; ++by) {
            for (int bx = 0; bx < bw; ++bx) {
                const Color c = m_cache[static_cast<std::size_t>(by) * bw + bx];
                if (c.a == 0) continue;
                for (int yy = 0; yy < step; ++yy) {
                    for (int xx = 0; xx < step; ++xx) canvas.put_pixel(bx * step + xx, by * step + yy, c);
                }
            }
        }
    }

    void draw_3d(Canvas& canvas, const Theme& theme) {
        if (!m_framed) {
            m_cam.target = Vec3{};
            m_cam.pan = Vec2{};
            m_cam.zoom = 0.45 * std::min(canvas.width(), canvas.height()) / std::max(30.0, m_extent);
            m_framed = true;
        }
        Scene3D scene(m_cam);
        const double L = m_extent;
        scene.line(Vec3{ -L, 0, 0 }, Vec3{ L, 0, 0 }, 1.0, theme.grid);
        scene.line(Vec3{ 0, -L, 0 }, Vec3{ 0, L, 0 }, 1.0, theme.grid);
        scene.line(Vec3{ 0, 0, -L }, Vec3{ 0, 0, L }, 1.0, theme.grid);
        scene.text(Vec3{ L * 1.05, 0, 0 }, "x", 11.0, theme.muted);
        scene.text(Vec3{ 0, L * 1.05, 0 }, "y", 11.0, theme.muted);
        scene.text(Vec3{ 0, 0, L * 1.05 }, "z", 11.0, theme.muted);
        const bool balloon = m_style3 == "balloon" || m_style3 == "balloon_dots";
        const bool dots = !balloon || m_style3 == "balloon_dots";
        scene.sphere_px(Vec3{}, 4.0, theme.red);
        int max_n = 1;
        for (const fc::CloudPoint& p : m_points) max_n = std::max(max_n, p.n);
        const double dot = m_points.size() > 20000 ? 1.0 : 1.4;
        if (dots) {
            for (const fc::CloudPoint& p : m_points) {
                Color c;
                if (m_style3 == "mono") c = rgb(205, 214, 244, 120);
                else if (m_style3 == "type") c = l_color(p.l, p.n, max_n), c.a = 150;
                else if (m_target == "orbital" || m_signed) c = p.sign >= 0 ? rgb(243, 139, 168, 170) : rgb(137, 180, 250, 170);
                else {
                    c = l_color(p.l, p.n, max_n);
                    c.a = 150;
                }
                if (balloon) c.a = 110;
                scene.disc(Vec3{ p.x, p.y, p.z }, dot, c);
            }
        }
        if (balloon && m_target != "orbital") {
            if (m_grid_dirty) build_grid();
            int shell_max = 1;
            for (const fc::RadialCurve& r : m_shells) shell_max = std::max(shell_max, r.n);
            std::vector<const fc::RadialCurve*> order_list;
            for (const fc::RadialCurve& r : m_shells) {
                if (m_target == "valence" && r.n < shell_max && r.l < 2) continue;
                if (r.occupancy > 0 && r.radius_90_pm > 0) order_list.push_back(&r);
            }
            std::sort(order_list.begin(), order_list.end(), [](const fc::RadialCurve* x, const fc::RadialCurve* y) { return x->radius_90_pm > y->radius_90_pm; });
            for (std::size_t i = 0; i < order_list.size(); ++i) {
                Color c = l_color(order_list[i]->l, order_list[i]->n, shell_max);
                c.a = m_style3 == "balloon_dots" ? 55 : 90;
                scene.sphere(Vec3{}, order_list[i]->radius_90_pm, c, -static_cast<double>(i) * 1e-3);
            }
        }
        scene.draw(canvas);
        if (balloon && m_target == "orbital") {
            draw_balloon(canvas, m_style3 == "balloon_dots");
            canvas.disc_px(m_cam.project(Vec3{}), 3.0, theme.red);
        }
    }

    std::string readouts() const {
        std::ostringstream o;
        const UnitPrefs& u = m_shared.units;
        const fc::Element& e = fc::element(m_z);
        o << "#" << title() << "\n";
        if (m_target == "orbital") {
            const double z = zeff();
            o << "Quantum numbers\tn = " << m_n << ", l = " << m_l << ", m\xE2\x82\x97 = " << m_m << "\n";
            o << "Shape\t" << fc::orbital_shape_name(m_l, m_m) << "\n";
            o << "Zeff\t" << fixed(z, 2) << "\n";
            o << "Radial nodes\t" << fc::radial_node_count(m_n, m_l) << "\n";
            o << "Angular nodes\t" << fc::angular_node_count(m_l) << "\n";
            try {
                o << "Most probable r\t" << u.format(fc::most_probable_radius_bohr(m_n, m_l, z, model_of(m_model)) * BOHR_PM * 1e-12, "length", 4) << "\n";
                o << "Average r\t" << u.format(fc::mean_radius_bohr(m_n, m_l, z) * BOHR_PM * 1e-12, "length", 4) << "\n";
                o << "Energy (approx.)\t" << fixed(fc::orbital_energy_ev(m_n, z), 2) << " eV\n";
            } catch (...) {
            }
            o << "The two colors are the + and \xE2\x88\x92 lobes of the wavefunction \xCF\x88; the probability is \xCF\x88\xC2\xB2 either way.\n";
        } else {
            o << "Configuration\t" << e.electron_configuration << "\n";
            o << "Colors\ts blue, p green, d orange, f purple (darker = inner shells)\n";
            o << "Each dot is one possible position of an electron, picked with probability |\xCF\x88|\xC2\xB2. Where dots crowd together, the electron is most likely to be.\n";
        }
        o << "Picture radius\t" << u.format(m_extent * 1e-12, "length", 4) << "\n";
        return o.str();
    }

    ChemShared& m_shared;
    int m_z = 6;
    int m_charge = 0;
    std::string m_target = "orbital";
    int m_n = 2;
    int m_l = 1;
    int m_m = 0;
    bool m_3d = true;
    std::string m_plane = "xz";
    bool m_signed = true;
    int m_count = 9000;
    std::string m_model = "hydrogenic";
    bool m_spin = true;
    std::string m_style3 = "sign";
    std::string m_style2 = "heat";
    std::vector<float> m_grid;
    int m_grid_n = 0;
    double m_iso = 0.0;
    bool m_grid_dirty = true;
    std::vector<Color> m_cache;
    std::string m_cache_key;
    std::vector<fc::RadialCurve> m_shells;
    double m_extent = 300.0;
    std::vector<fc::CloudPoint> m_points;
    fc::DensitySlice m_slice;
    std::string m_error;
    Camera3D m_cam;
    bool m_framed = false;
    bool m_down = false;
    Vec2 m_last;
};

} // namespace

std::unique_ptr<ChemView> make_orbital_view(ChemShared& shared) {
    return std::unique_ptr<ChemView>(new RadialView(shared));
}

std::unique_ptr<ChemView> make_cloud_view(ChemShared& shared) {
    return std::unique_ptr<ChemView>(new CloudView(shared));
}

} // namespace sapling

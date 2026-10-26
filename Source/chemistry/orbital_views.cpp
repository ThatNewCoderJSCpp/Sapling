#include "molecule_draw.hpp"
#include "view_util.hpp"

#include <algorithm>
#include <cmath>
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
        if (m_total && !m_cumulative) {
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
        if (m_total) plot(m_total_curve, theme.text, 2.6);
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
        if (m_total) {
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
        draw_hud(canvas, theme, { title(), m_3d ? std::to_string(m_points.size()) + " sampled positions; denser = more likely" : std::string("slice through the nucleus (") + m_plane + " plane)" + (m_signed ? ", colored by the sign of \xCF\x88" : ", brightness = probability density") }, false);
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
        const int px = static_cast<int>(side);
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
                const auto at = [&](int a, int b) { return m_slice.values[static_cast<std::size_t>(b) * W + a]; };
                const double v = (at(i0, j0) * (1 - tx) + at(i1, j0) * tx) * (1 - ty) + (at(i0, j1) * (1 - tx) + at(i1, j1) * tx) * ty;
                Color c;
                if (m_signed) c = signed_color(maxabs > 0 ? v / maxabs : 0.0, theme);
                else c = heat_color(m_slice.max_value > 0 ? std::pow(std::max(0.0, v) / m_slice.max_value, 0.45) : 0.0, theme);
                canvas.put_pixel(static_cast<int>(x0) + i, static_cast<int>(y0) + j, c);
            }
        }
        canvas.stroke_rect_px(x0, y0, side, side, theme.axes, 1.0);
        const UnitPrefs& u = m_shared.units;
        canvas.text(Vec2{ x0 + side * 0.5, y0 + side + 18.0 }, std::string(1, m_plane[0]) + " axis, \xC2\xB1" + u.format(m_extent * 1e-12, "length", 3), 11.0, theme.muted, TextAlign::Center);
        canvas.text(Vec2{ x0 - 8.0, y0 + side * 0.5 }, std::string(1, m_plane[1]), 12.0, theme.muted, TextAlign::Right);
        canvas.disc_px(Vec2{ x0 + side * 0.5, y0 + side * 0.5 }, 2.0, theme.red);
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
        scene.sphere_px(Vec3{}, 4.0, theme.red);
        int max_n = 1;
        for (const fc::CloudPoint& p : m_points) max_n = std::max(max_n, p.n);
        const double dot = m_points.size() > 20000 ? 1.0 : 1.4;
        for (const fc::CloudPoint& p : m_points) {
            Color c;
            if (m_target == "orbital" || m_signed) c = p.sign >= 0 ? rgb(243, 139, 168, 170) : rgb(137, 180, 250, 170);
            else {
                c = l_color(p.l, p.n, max_n);
                c.a = 150;
            }
            scene.disc(Vec3{ p.x, p.y, p.z }, dot, c);
        }
        scene.draw(canvas);
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

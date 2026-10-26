#include "molecule_draw.hpp"
#include "view_util.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

namespace sapling {

using namespace chem_util;

namespace {

constexpr double R_GAS = 8.314462618;
constexpr double K_B = 1.380649e-23;
constexpr double N_A = 6.02214076e23;

struct Series {
    std::string name;
    Color color;
    std::vector<std::pair<double, double>> points;
    bool dashed = false;
    bool dots = false;
};

struct Marker {
    double x = 0.0;
    double y = 0.0;
    std::string label;
    Color color;
};

Color palette(int i) {
    static const Color c[] = { rgb(137, 180, 250), rgb(166, 227, 161), rgb(250, 179, 135), rgb(203, 166, 247), rgb(243, 139, 168), rgb(249, 226, 175), rgb(148, 226, 213), rgb(245, 194, 231) };
    return c[((i % 8) + 8) % 8];
}

void plot_series(Canvas& canvas, const Theme& theme, PlotFrame& f, const std::vector<Series>& series, const std::vector<Marker>& markers, const std::string& xl, const std::string& yl, bool auto_range, double top) {
    f.x0 = 74.0;
    f.y0 = top;
    f.w = std::max(60.0, canvas.width() - f.x0 - 24.0);
    f.h = std::max(60.0, canvas.height() - f.y0 - 56.0);
    if (auto_range) {
        double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
        for (const Series& s : series) {
            for (const auto& p : s.points) {
                if (!finite(p.first) || !finite(p.second)) continue;
                x0 = std::min(x0, p.first);
                x1 = std::max(x1, p.first);
                y0 = std::min(y0, p.second);
                y1 = std::max(y1, p.second);
            }
        }
        if (!(x1 > x0)) {
            x0 = 0.0;
            x1 = 1.0;
        }
        if (!(y1 > y0)) {
            y0 -= 0.5;
            y1 += 0.5;
        }
        const double pad = 0.06 * (y1 - y0);
        f.xmin = x0;
        f.xmax = x1;
        f.ymin = y0 >= 0.0 && y0 - pad < 0.0 ? 0.0 : y0 - pad;
        f.ymax = y1 + pad;
    }
    draw_plot_axes(canvas, f, theme, xl, yl);
    for (const Series& s : series) {
        Vec2 prev{};
        bool have = false;
        for (const auto& p : s.points) {
            if (!finite(p.first) || !finite(p.second)) {
                have = false;
                continue;
            }
            const Vec2 q = f.to_screen(p.first, std::max(f.ymin - (f.ymax - f.ymin), std::min(f.ymax + (f.ymax - f.ymin), p.second)));
            if (s.dots) {
                canvas.disc_px(q, 3.5, s.color);
                continue;
            }
            if (have) {
                if (s.dashed) canvas.dashed_line_px(prev, q, s.color, 1.6, 5.0);
                else canvas.line_px(prev, q, s.color, 2.2);
            }
            prev = q;
            have = true;
        }
    }
    for (const Marker& m : markers) {
        const Vec2 q = f.to_screen(m.x, m.y);
        if (!f.inside(q)) continue;
        canvas.disc_px(q, 4.5, m.color);
        canvas.text(Vec2{ q.x + 7.0, q.y - 6.0 }, m.label, 11.0, m.color);
    }
    double ly = f.y0 + 8.0;
    for (const Series& s : series) {
        if (s.name.empty()) continue;
        const double w = canvas.text_width(s.name, 11.0);
        const double lx = f.x0 + f.w - w - 34.0;
        canvas.fill_rect_px(lx - 4.0, ly - 2.0, w + 34.0, 16.0, rgb(theme.background.r, theme.background.g, theme.background.b, 180));
        canvas.line_px(Vec2{ lx, ly + 6.0 }, Vec2{ lx + 20.0, ly + 6.0 }, s.color, 2.2);
        canvas.text(Vec2{ lx + 26.0, ly + 10.0 }, s.name, 11.0, theme.text);
        ly += 17.0;
    }
}

std::vector<std::string> pieces(const std::string& s, const std::string& seps) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (seps.find(c) != std::string::npos) {
            if (!trim(cur).empty()) out.push_back(trim(cur));
            cur.clear();
        } else {
            cur += c;
        }
    }
    if (!trim(cur).empty()) out.push_back(trim(cur));
    return out;
}

class KineticsView : public ChemView {
public:
    explicit KineticsView(ChemShared& shared) : m_shared(shared) { rebuild(); }

    void draw(Canvas& canvas, const Theme& theme) override {
        canvas.clear(theme.background);
        if (!m_error.empty()) {
            draw_hud(canvas, theme, { "Kinetics", m_error }, false);
            return;
        }
        plot_series(canvas, theme, m_frame, m_series, m_markers, m_xlabel, m_ylabel, true, 58.0);
        draw_hud(canvas, theme, { m_title }, false);
        if (m_frame.inside(m_mouse)) {
            canvas.line_px(Vec2{ m_mouse.x, m_frame.y0 }, Vec2{ m_mouse.x, m_frame.y0 + m_frame.h }, rgb(theme.accent.r, theme.accent.g, theme.accent.b, 110), 1.0);
            const double x = m_frame.x_of(m_mouse.x);
            std::string line = m_xlabel.substr(0, m_xlabel.find(' ')) + " = " + format_number(x, 4) + ":";
            int shown = 0;
            for (const Series& s : m_series) {
                if (s.points.empty() || s.dots) continue;
                double best = 0.0, bd = 1e300;
                for (const auto& p : s.points) {
                    if (std::fabs(p.first - x) < bd) {
                        bd = std::fabs(p.first - x);
                        best = p.second;
                    }
                }
                line += "  " + (s.name.empty() ? std::string("y") : s.name) + " " + format_number(best, 4);
                if (++shown >= 4) break;
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
        if (c == "mode") {
            m_mode = arg(a, 1);
        } else if (c == "rate") {
            m_mode = "rate";
            m_order = num(arg(a, 1), 1.0);
            m_k = num(arg(a, 2), 0.1);
            m_a0 = num(arg(a, 3), 1.0);
            m_tmax = num(arg(a, 4), 0.0);
            m_plot = arg(a, 5).empty() ? "conc" : arg(a, 5);
        } else if (c == "mechanism") {
            m_mode = "mechanism";
            m_steps_text = arg(a, 1);
            m_init_text = arg(a, 2);
            m_tmax = num(arg(a, 3), 0.0);
        } else if (c == "arrhenius") {
            m_mode = "arrhenius";
            m_ea = num(arg(a, 1), 50000.0);
            m_prefactor = num(arg(a, 2), 1e10);
            m_t1 = num(arg(a, 3), 298.15);
            m_t2 = num(arg(a, 4), 348.15);
            m_plot = arg(a, 5).empty() ? "k" : arg(a, 5);
        } else if (c == "profile") {
            m_mode = "profile";
            m_dh = num(arg(a, 1), -50000.0);
            m_ea = num(arg(a, 2), 80000.0);
            m_ea_cat = num(arg(a, 3), fc::NO_VALUE);
        } else if (c == "boltzmann") {
            m_mode = "boltzmann";
            m_molar = num(arg(a, 1), fc::NO_VALUE);
            if (!finite(m_molar)) {
                try {
                    m_molar = fc::Formula::parse(trim(arg(a, 1))).molar_mass();
                } catch (...) {
                    return "error\tnot a formula or molar mass: " + arg(a, 1);
                }
            }
            if (!(m_molar > 0.0)) m_molar = 28.0;
            m_t1 = num(arg(a, 2), 298.15);
            m_t2 = num(arg(a, 3), 600.0);
            m_ea = num(arg(a, 4), 50000.0);
        } else if (c == "fit") {
            m_mode = "fit";
            m_data_text = arg(a, 1);
        } else if (c == "readouts") {
            return m_readouts;
        } else {
            return "error\tunknown command " + c;
        }
        rebuild();
        return m_error.empty() ? "ok" : "error\t" + m_error;
    }

private:
    void rebuild() {
        m_error.clear();
        m_series.clear();
        m_markers.clear();
        try {
            if (m_mode == "rate") build_rate();
            else if (m_mode == "mechanism") build_mechanism();
            else if (m_mode == "arrhenius") build_arrhenius();
            else if (m_mode == "profile") build_profile();
            else if (m_mode == "boltzmann") build_boltzmann();
            else build_fit();
        } catch (const std::exception& e) {
            m_error = e.what();
            m_readouts = std::string("#Kinetics\n!") + e.what() + "\n";
        }
    }

    double concentration(double t) const {
        if (std::fabs(m_order) < 1e-9) return fc::zero_order_concentration(m_a0, m_k, t);
        if (std::fabs(m_order - 1.0) < 1e-9) return fc::first_order_concentration(m_a0, m_k, t);
        if (std::fabs(m_order - 2.0) < 1e-9) return fc::second_order_concentration(m_a0, m_k, t);
        return fc::nth_order_concentration(m_a0, m_k, t, m_order);
    }

    double half_life_at(double c0) const {
        if (std::fabs(m_order) < 1e-9) return fc::zero_order_half_life(c0, m_k);
        if (std::fabs(m_order - 1.0) < 1e-9) return fc::first_order_half_life(m_k);
        if (std::fabs(m_order - 2.0) < 1e-9) return fc::second_order_half_life(c0, m_k);
        return fc::nth_order_half_life(c0, m_k, m_order);
    }

    void build_rate() {
        if (!(m_k > 0) || !(m_a0 > 0)) throw std::invalid_argument("k and [A]\xE2\x82\x80 have to be positive");
        const double t12 = half_life_at(m_a0);
        double tmax = m_tmax > 0 ? m_tmax : 5.0 * t12;
        if (std::fabs(m_order) < 1e-9) tmax = m_tmax > 0 ? m_tmax : 1.2 * m_a0 / m_k;
        Series s{ "[A]", palette(0), {}, false, false };
        Series p{ "[P]", palette(1), {}, false, false };
        for (int i = 0; i <= 400; ++i) {
            const double t = tmax * i / 400.0;
            const double c = std::max(0.0, concentration(t));
            double y = c;
            if (m_plot == "ln") y = c > 0 ? std::log(c) : fc::NO_VALUE;
            else if (m_plot == "inv") y = c > 0 ? 1.0 / c : fc::NO_VALUE;
            s.points.emplace_back(t, y);
            if (m_plot == "conc") p.points.emplace_back(t, m_a0 - c);
        }
        m_series.push_back(s);
        if (m_plot == "conc") m_series.push_back(p);
        double c0 = m_a0;
        double t = 0.0;
        for (int i = 1; i <= 4 && finite(t12); ++i) {
            const double h = half_life_at(c0);
            if (!finite(h) || t + h > tmax) break;
            t += h;
            c0 *= 0.5;
            const double y = m_plot == "ln" ? std::log(c0) : m_plot == "inv" ? 1.0 / c0 : c0;
            m_markers.push_back(Marker{ t, y, std::to_string(i) + (i == 1 ? " half-life" : " half-lives"), palette(5) });
        }
        static const char* names[] = { "zero", "first", "second", "third" };
        const int whole = static_cast<int>(std::lround(m_order));
        const std::string order_name = std::fabs(m_order - whole) < 1e-9 && whole >= 0 && whole <= 3 ? std::string(names[whole]) + "-order" : "order " + format_number(m_order, 3);
        m_title = order_name + " reaction: rate = k[A]" + (whole == 1 ? std::string() : superscript_digits(format_number(m_order, 3)));
        m_xlabel = "time (s)";
        m_ylabel = m_plot == "ln" ? "ln [A]" : m_plot == "inv" ? "1/[A] (1/M)" : "concentration (M)";
        std::ostringstream o;
        o << "#" << order_name << " kinetics\n";
        o << "Rate law\trate = k[A]^" << format_number(m_order, 3) << "\n";
        o << "k\t" << format_number(m_k, 5) << " " << (whole == 1 ? "1/s" : whole == 0 ? "M/s" : "M^" + format_number(1.0 - m_order, 3) + "/s") << "\n";
        o << "Initial rate\t" << format_number(m_k * std::pow(m_a0, m_order), 5) << " M/s\n";
        o << "First half-life\t" << format_number(t12, 5) << " s\n";
        if (whole == 1) o << "Half-life is constant: each one takes " << format_number(t12, 4) << " s no matter the concentration.\n";
        else if (whole == 0) o << "Half-lives get shorter as A runs out (t\xC2\xBD = [A]\xE2\x82\x80/2k).\n";
        else o << "Half-lives get longer as A runs out (t\xC2\xBD = 1/(k[A]\xE2\x82\x80) for second order).\n";
        for (double frac : { 0.5, 0.9, 0.99 }) {
            double lo = 0.0, hi = tmax * 100.0;
            for (int it = 0; it < 200; ++it) {
                const double mid = 0.5 * (lo + hi);
                if (concentration(mid) > m_a0 * (1.0 - frac)) lo = mid;
                else hi = mid;
            }
            o << "Time to " << format_number(frac * 100.0, 3) << "% done\t" << format_number(hi, 5) << " s\n";
        }
        o << "#Straight-line test\n";
        o << "Zero order\t[A] vs t is a line\n";
        o << "First order\tln[A] vs t is a line (slope \xE2\x88\x92k)\n";
        o << "Second order\t1/[A] vs t is a line (slope k)\n";
        m_readouts = o.str();
    }

    void build_mechanism() {
        std::map<std::string, std::size_t> index;
        std::vector<std::string> names;
        const auto species = [&](const std::string& n) {
            auto it = index.find(n);
            if (it != index.end()) return it->second;
            index[n] = names.size();
            names.push_back(n);
            return names.size() - 1;
        };
        const auto side = [&](const std::string& text) {
            std::vector<std::pair<std::size_t, double>> out;
            for (const std::string& term : pieces(text, "+")) {
                std::size_t i = 0;
                while (i < term.size() && (std::isdigit(static_cast<unsigned char>(term[i])) || term[i] == '.')) ++i;
                const double coef = i > 0 ? num(term.substr(0, i), 1.0) : 1.0;
                const std::string n = trim(term.substr(i));
                if (!n.empty()) out.emplace_back(species(n), coef);
            }
            return out;
        };
        fc::ReactionNetwork net;
        std::vector<std::string> labels;
        for (const std::string& step : pieces(m_steps_text, ";\n")) {
            const std::size_t colon = step.find(':');
            const std::string eq = trim(step.substr(0, colon));
            const auto ks = colon == std::string::npos ? std::vector<std::string>{} : pieces(step.substr(colon + 1), ",");
            bool reversible = false;
            std::size_t arrow = eq.find("<=>");
            std::size_t alen = 3;
            if (arrow == std::string::npos) {
                arrow = eq.find("<->");
            }
            if (arrow != std::string::npos) reversible = true;
            else {
                arrow = eq.find("->");
                alen = 2;
            }
            if (arrow == std::string::npos) throw std::invalid_argument("each step needs -> or <=> (like A -> B : 0.5)");
            const auto lhs = side(eq.substr(0, arrow));
            const auto rhs = side(eq.substr(arrow + alen));
            const double kf = ks.empty() ? 1.0 : num(ks[0], 1.0);
            const double kr = reversible ? (ks.size() > 1 ? num(ks[1], 0.1) : 0.1) : 0.0;
            net.add(lhs, rhs, kf, kr);
            labels.push_back(eq);
        }
        if (net.steps().empty()) throw std::invalid_argument("type at least one step, like A -> B : 0.5; B -> C : 0.1");
        net.set_species(names.size());
        std::vector<double> c(names.size(), 0.0);
        for (const std::string& entry : pieces(m_init_text, ",;")) {
            const std::size_t eq = entry.find('=');
            if (eq == std::string::npos) continue;
            const std::string n = trim(entry.substr(0, eq));
            if (index.count(n)) c[index[n]] = num(entry.substr(eq + 1), 0.0);
        }
        double kmin = 1e300;
        for (const fc::ElementaryStep& s : net.steps()) {
            if (s.forward_rate > 0) kmin = std::min(kmin, s.forward_rate);
        }
        const double tmax = m_tmax > 0 ? m_tmax : (kmin < 1e300 ? 6.0 / kmin : 10.0);
        std::vector<Series> series(names.size());
        for (std::size_t i = 0; i < names.size(); ++i) series[i] = Series{ names[i], palette(static_cast<int>(i)), {}, false, false };
        const int steps = 2000;
        const double dt = tmax / steps;
        std::vector<double> peak(names.size(), 0.0), peak_t(names.size(), 0.0);
        for (int n = 0; n <= steps; ++n) {
            if (n % 5 == 0) {
                for (std::size_t i = 0; i < names.size(); ++i) series[i].points.emplace_back(n * dt, c[i]);
            }
            for (std::size_t i = 0; i < names.size(); ++i) {
                if (c[i] > peak[i]) {
                    peak[i] = c[i];
                    peak_t[i] = n * dt;
                }
            }
            if (n == steps) break;
            std::vector<double> k1 = net.derivatives(c), tmp(c.size());
            for (std::size_t i = 0; i < c.size(); ++i) tmp[i] = c[i] + 0.5 * dt * k1[i];
            std::vector<double> k2 = net.derivatives(tmp);
            for (std::size_t i = 0; i < c.size(); ++i) tmp[i] = c[i] + 0.5 * dt * k2[i];
            std::vector<double> k3 = net.derivatives(tmp);
            for (std::size_t i = 0; i < c.size(); ++i) tmp[i] = c[i] + dt * k3[i];
            std::vector<double> k4 = net.derivatives(tmp);
            for (std::size_t i = 0; i < c.size(); ++i) c[i] = std::max(0.0, c[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
        }
        m_series = series;
        m_title = "mechanism with " + std::to_string(net.steps().size()) + " elementary step" + (net.steps().size() == 1 ? "" : "s");
        m_xlabel = "time (s)";
        m_ylabel = "concentration (M)";
        std::ostringstream o;
        o << "#Mechanism\n";
        for (std::size_t i = 0; i < labels.size(); ++i) {
            const fc::ElementaryStep& s = net.steps()[i];
            o << "Step " << i + 1 << "\t" << labels[i] << "  (k = " << format_number(s.forward_rate, 4) << (s.reverse_rate > 0 ? ", k\xE2\x82\x8B = " + format_number(s.reverse_rate, 4) : std::string()) << ")\n";
        }
        o << "Each step's rate law comes from its own reactants (elementary steps), and the curves come from integrating them all together (RK4).\n";
        o << "#At t = " << format_number(tmax, 4) << " s\n";
        for (std::size_t i = 0; i < names.size(); ++i) o << "[" << names[i] << "]\t" << format_number(c[i], 5) << " M\n";
        o << "#Peaks\n";
        for (std::size_t i = 0; i < names.size(); ++i) {
            if (peak_t[i] > 0.0 && peak_t[i] < tmax * 0.98) o << names[i] << " (intermediate)\tpeaks at " << format_number(peak[i], 4) << " M, t = " << format_number(peak_t[i], 4) << " s\n";
        }
        m_readouts = o.str();
    }

    void build_arrhenius() {
        if (!(m_prefactor > 0) || !(m_t1 > 0)) throw std::invalid_argument("A and the temperatures have to be positive");
        const double lo = std::min(m_t1, m_t2) * 0.7, hi = std::max(m_t1, m_t2) * 1.3;
        Series s{ m_plot == "ln" ? "ln k" : "k", palette(0), {}, false, false };
        for (int i = 0; i <= 300; ++i) {
            const double t = lo + (hi - lo) * i / 300.0;
            const double k = fc::arrhenius(m_prefactor, m_ea, t);
            if (m_plot == "ln") s.points.emplace_back(1000.0 / t, std::log(k));
            else s.points.emplace_back(t, k);
        }
        m_series.push_back(s);
        const double k1 = fc::arrhenius(m_prefactor, m_ea, m_t1), k2 = fc::arrhenius(m_prefactor, m_ea, m_t2);
        if (m_plot == "ln") {
            m_markers.push_back(Marker{ 1000.0 / m_t1, std::log(k1), "T\xE2\x82\x81", palette(1) });
            m_markers.push_back(Marker{ 1000.0 / m_t2, std::log(k2), "T\xE2\x82\x82", palette(2) });
        } else {
            m_markers.push_back(Marker{ m_t1, k1, "T\xE2\x82\x81", palette(1) });
            m_markers.push_back(Marker{ m_t2, k2, "T\xE2\x82\x82", palette(2) });
        }
        m_title = "Arrhenius: k = A e^(\xE2\x88\x92" "E\xE2\x82\x90/RT)";
        m_xlabel = m_plot == "ln" ? "1000/T (1/K)" : "T (K)";
        m_ylabel = m_plot == "ln" ? "ln k (slope = \xE2\x88\x92" "E\xE2\x82\x90/R)" : "k";
        const UnitPrefs& u = m_shared.units;
        std::ostringstream o;
        o << "#Arrhenius\n";
        o << "E\xE2\x82\x90\t" << u.format(m_ea, "molar energy", 4) << "\n";
        o << "A\t" << format_number(m_prefactor, 4) << "\n";
        o << "k at T\xE2\x82\x81 (" << u.format(m_t1, "temperature", 4) << ")\t" << format_number(k1, 5) << "\n";
        o << "k at T\xE2\x82\x82 (" << u.format(m_t2, "temperature", 4) << ")\t" << format_number(k2, 5) << "\n";
        o << "k\xE2\x82\x82 / k\xE2\x82\x81\t" << format_number(k2 / k1, 5) << "\n";
        o << "Q\xE2\x82\x81\xE2\x82\x80 (per 10 K near T\xE2\x82\x81)\t" << format_number(fc::arrhenius(m_prefactor, m_ea, m_t1 + 10.0) / k1, 4) << "\n";
        o << "Fraction of collisions with E \xE2\x89\xA5 E\xE2\x82\x90 at T\xE2\x82\x81\t" << format_number(std::exp(-m_ea / (R_GAS * m_t1)), 4) << "\n";
        o << "Raising T makes many more collisions energetic enough, so k grows exponentially.\n";
        m_readouts = o.str();
    }

    void build_profile() {
        const double ea = m_ea / 1000.0, dh = m_dh / 1000.0;
        const auto curve = [&](double barrier) {
            std::vector<std::pair<double, double>> pts;
            for (int i = 0; i <= 300; ++i) {
                const double x = i / 300.0;
                const double base = dh * (0.5 - 0.5 * std::cos(3.14159265358979 * std::min(1.0, std::max(0.0, (x - 0.1) / 0.8))));
                const double bump = (barrier - 0.5 * dh) * std::exp(-std::pow((x - 0.5) / 0.13, 2.0));
                pts.emplace_back(x, base + bump);
            }
            return pts;
        };
        m_series.push_back(Series{ "uncatalyzed", palette(4), curve(ea), false, false });
        if (finite(m_ea_cat) && m_ea_cat > 0) m_series.push_back(Series{ "catalyzed", palette(1), curve(m_ea_cat / 1000.0), true, false });
        m_markers.push_back(Marker{ 0.05, 0.0, "reactants", palette(0) });
        m_markers.push_back(Marker{ 0.95, dh, "products", palette(2) });
        m_markers.push_back(Marker{ 0.5, ea + 0.0 * dh, "transition state", palette(5) });
        m_title = dh < 0 ? "exothermic reaction profile" : "endothermic reaction profile";
        m_xlabel = "reaction progress";
        m_ylabel = "energy (kJ/mol, relative to reactants)";
        const UnitPrefs& u = m_shared.units;
        std::ostringstream o;
        o << "#Energy profile\n";
        o << "\xCE\x94H\t" << u.format(m_dh, "molar energy", 4) << (m_dh < 0 ? " (exothermic)" : " (endothermic)") << "\n";
        o << "E\xE2\x82\x90 forward\t" << u.format(m_ea, "molar energy", 4) << "\n";
        o << "E\xE2\x82\x90 reverse\t" << u.format(m_ea - m_dh, "molar energy", 4) << "\n";
        if (finite(m_ea_cat) && m_ea_cat > 0) {
            o << "E\xE2\x82\x90 with catalyst\t" << u.format(m_ea_cat, "molar energy", 4) << "\n";
            o << "Speed-up at 25 \xC2\xB0" "C\t" << format_number(std::exp((m_ea - m_ea_cat) / (R_GAS * 298.15)), 4) << "\xC3\x97\n";
            o << "A catalyst lowers the barrier but not \xCE\x94H, so it speeds up both directions and doesn't move the equilibrium.\n";
        }
        if (m_ea < std::max(0.0, m_dh)) o << "!E\xE2\x82\x90 can't be lower than \xCE\x94H for an endothermic reaction\n";
        m_readouts = o.str();
    }

    void build_boltzmann() {
        const double m = m_molar / 1000.0 / N_A;
        if (!(m > 0)) throw std::invalid_argument("give a molar mass");
        const auto f = [&](double v, double t) { return 4.0 * 3.14159265358979 * std::pow(m / (2.0 * 3.14159265358979 * K_B * t), 1.5) * v * v * std::exp(-m * v * v / (2.0 * K_B * t)); };
        const double vmax = 4.0 * std::sqrt(2.0 * K_B * std::max(m_t1, m_t2) / m);
        Series a{ "T\xE2\x82\x81", palette(0), {}, false, false }, b{ "T\xE2\x82\x82", palette(4), {}, false, false };
        for (int i = 0; i <= 400; ++i) {
            const double v = vmax * i / 400.0;
            a.points.emplace_back(v, f(v, m_t1) * 1000.0);
            b.points.emplace_back(v, f(v, m_t2) * 1000.0);
        }
        m_series = { a, b };
        const double vea = std::sqrt(2.0 * m_ea / N_A / m);
        if (vea < vmax) m_series.push_back(Series{ "speed for E\xE2\x82\x90", palette(5), { { vea, 0.0 }, { vea, f(std::sqrt(2.0 * K_B * m_t1 / m), m_t1) * 1000.0 } }, true, false });
        for (double t : { m_t1, m_t2 }) m_markers.push_back(Marker{ std::sqrt(2.0 * K_B * t / m), f(std::sqrt(2.0 * K_B * t / m), t) * 1000.0, "most probable", t == m_t1 ? palette(0) : palette(4) });
        m_title = "Maxwell\xE2\x80\x93" "Boltzmann speeds";
        m_xlabel = "speed (m/s)";
        m_ylabel = "probability density (per km/s)";
        const UnitPrefs& u = m_shared.units;
        std::ostringstream o;
        o << "#Molecular speeds\n";
        for (double t : { m_t1, m_t2 }) {
            o << "#At " << u.format(t, "temperature", 4) << "\n";
            o << "Most probable\t" << format_number(std::sqrt(2.0 * K_B * t / m), 4) << " m/s\n";
            o << "Average\t" << format_number(std::sqrt(8.0 * K_B * t / (3.14159265358979 * m)), 4) << " m/s\n";
            o << "RMS\t" << format_number(std::sqrt(3.0 * K_B * t / m), 4) << " m/s\n";
            o << "Fraction above E\xE2\x82\x90\t" << format_number(std::exp(-m_ea / (R_GAS * t)) * (1.0 + 0.0), 4) << "\n";
        }
        o << "Heating spreads the curve out and moves it right, so far more molecules clear the activation energy.\n";
        m_readouts = o.str();
    }

    void build_fit() {
        std::vector<std::pair<double, double>> data;
        for (const std::string& pair : pieces(m_data_text, ";\n")) {
            const auto xy = pieces(pair, ", \t");
            if (xy.size() < 2) continue;
            const double t = num(xy[0], fc::NO_VALUE), c = num(xy[1], fc::NO_VALUE);
            if (finite(t) && finite(c) && c > 0) data.emplace_back(t, c);
        }
        if (data.size() < 3) throw std::invalid_argument("give at least three t, [A] pairs, like 0, 1.0; 10, 0.74; 20, 0.55");
        struct FitResult {
            double slope, intercept, r2;
        };
        const auto fit = [&](int order) {
            double sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
            const double n = static_cast<double>(data.size());
            for (const auto& p : data) {
                const double x = p.first, y = order == 0 ? p.second : order == 1 ? std::log(p.second) : 1.0 / p.second;
                sx += x;
                sy += y;
                sxx += x * x;
                sxy += x * y;
                syy += y * y;
            }
            const double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
            const double inter = (sy - slope * sx) / n;
            const double r = (n * sxy - sx * sy) / std::sqrt((n * sxx - sx * sx) * (n * syy - sy * sy));
            return FitResult{ slope, inter, r * r };
        };
        FitResult fits[3] = { fit(0), fit(1), fit(2) };
        int best = 0;
        for (int i = 1; i < 3; ++i) {
            if (fits[i].r2 > fits[best].r2) best = i;
        }
        Series d{ "data", palette(5), data, false, true };
        m_series.push_back(d);
        const double tmax = data.back().first;
        static const char* names[] = { "zero order", "first order", "second order" };
        for (int order = 0; order < 3; ++order) {
            Series s{ names[order], palette(order), {}, order != best, false };
            for (int i = 0; i <= 200; ++i) {
                const double t = tmax * i / 200.0;
                double c = 0.0;
                if (order == 0) c = fits[0].intercept + fits[0].slope * t;
                else if (order == 1) c = std::exp(fits[1].intercept + fits[1].slope * t);
                else c = 1.0 / (fits[2].intercept + fits[2].slope * t);
                s.points.emplace_back(t, c);
            }
            m_series.push_back(s);
        }
        m_title = "best fit: " + std::string(names[best]);
        m_xlabel = "time (s)";
        m_ylabel = "concentration (M)";
        std::ostringstream o;
        o << "#Which order fits?\n";
        for (int i = 0; i < 3; ++i) {
            const double k = i == 0 ? -fits[i].slope : i == 1 ? -fits[i].slope : fits[i].slope;
            o << names[i] << "\tR\xC2\xB2 = " << format_number(fits[i].r2, 5) << ", k = " << format_number(k, 4) << (i == best ? "  (best)" : "") << "\n";
        }
        o << "The order whose straight-line plot ([A], ln[A] or 1/[A] against t) is straightest wins.\n";
        m_readouts = o.str();
    }

    ChemShared& m_shared;
    std::string m_mode = "rate";
    double m_order = 1.0, m_k = 0.1, m_a0 = 1.0, m_tmax = 0.0;
    std::string m_plot = "conc";
    std::string m_steps_text = "A -> B : 0.3; B -> C : 0.1";
    std::string m_init_text = "A = 1";
    double m_ea = 50000.0, m_prefactor = 1e10, m_t1 = 298.15, m_t2 = 348.15;
    double m_dh = -50000.0, m_ea_cat = 50000.0, m_molar = 28.0;
    std::string m_data_text = "0, 1.0; 10, 0.74; 20, 0.55; 30, 0.41; 40, 0.30; 60, 0.17";
    std::vector<Series> m_series;
    std::vector<Marker> m_markers;
    std::string m_title, m_xlabel, m_ylabel, m_readouts, m_error;
    PlotFrame m_frame;
    Vec2 m_mouse{ -1.0, -1.0 };
};

struct AcidSpec {
    std::string name;
    std::vector<double> pka;
    bool base = false;
    bool strong = false;
};

bool acid_spec(const std::string& text, AcidSpec& out) {
    const std::string t = trim(text);
    if (t.empty()) return false;
    if (t.rfind("pKa", 0) == 0 || t.rfind("pka", 0) == 0) {
        const std::size_t eq = t.find('=');
        out.name = "custom acid";
        for (const std::string& p : pieces(eq == std::string::npos ? t.substr(3) : t.substr(eq + 1), ", ")) out.pka.push_back(num(p, 4.76));
        return !out.pka.empty();
    }
    if (const fc::AcidConstantEntry* a = fc::find_acid(t)) {
        out.name = a->name;
        out.strong = a->strong;
        for (double p : { a->pka1, a->pka2, a->pka3 }) {
            if (p == p) out.pka.push_back(p);
        }
        if (out.pka.empty()) out.pka.push_back(-6.0);
        if (a->strong) out.pka[0] = std::min(out.pka[0], -3.0);
        return true;
    }
    if (const fc::BaseConstantEntry* b = fc::find_base(t)) {
        out.name = b->name;
        out.base = true;
        out.pka = { 14.0 - b->pkb };
        return true;
    }
    try {
        if (fc::is_strong_acid(t)) {
            const fc::Formula f = fc::Formula::parse(t);
            out.name = t;
            out.strong = true;
            const int h = std::max(1, static_cast<int>(std::lround(f.count(1))));
            out.pka.assign(1, -6.0);
            if (h > 1) out.pka.push_back(t == "H2SO4" ? 1.99 : -4.0);
            return true;
        }
        if (fc::is_strong_base(t)) {
            const fc::Formula f = fc::Formula::parse(t);
            out.name = t;
            out.strong = true;
            out.base = true;
            const int oh = std::max(1, static_cast<int>(std::lround(f.count(8))));
            out.pka.assign(static_cast<std::size_t>(oh), 20.0);
            return true;
        }
    } catch (...) {
    }
    return false;
}

class AcidBaseView : public ChemView {
public:
    explicit AcidBaseView(ChemShared& shared) : m_shared(shared) { rebuild(); }

    void draw(Canvas& canvas, const Theme& theme) override {
        canvas.clear(theme.background);
        if (!m_error.empty()) {
            draw_hud(canvas, theme, { "Acids and bases", m_error }, false);
            return;
        }
        if (m_mode == "titration" || m_mode == "species") {
            m_frame.xmin = m_xmin;
            m_frame.xmax = m_xmax;
            m_frame.ymin = m_ymin;
            m_frame.ymax = m_ymax;
            plot_series(canvas, theme, m_frame, m_series, m_markers, m_xlabel, m_ylabel, false, 58.0);
            if (m_mode == "titration") {
                for (const auto& band : m_bands) {
                    const Vec2 a = m_frame.to_screen(m_xmin, band.second.first), b = m_frame.to_screen(m_xmax, band.second.second);
                    canvas.fill_rect_px(a.x, b.y, b.x - a.x, a.y - b.y, rgb(band.first.r, band.first.g, band.first.b, 28));
                }
            }
        } else {
            draw_scale(canvas, theme);
        }
        draw_hud(canvas, theme, { m_title }, false);
        if (m_frame.inside(m_mouse) && (m_mode == "titration" || m_mode == "species")) {
            const double x = m_frame.x_of(m_mouse.x);
            std::string line = (m_mode == "titration" ? "V = " + format_number(x, 4) + " " + m_vunit : (m_poh ? "pOH " : "pH ") + format_number(x, 3) + " (" + (m_poh ? "pH " : "pOH ") + format_number(pkw() - x, 3) + ")") + ":";
            for (const Series& s : m_series) {
                if (s.points.empty()) continue;
                double best = 0.0, bd = 1e300;
                for (const auto& p : s.points) {
                    if (std::fabs(p.first - x) < bd) {
                        bd = std::fabs(p.first - x);
                        best = p.second;
                    }
                }
                line += "  " + (s.name.empty() ? std::string(m_poh ? "pOH" : "pH") : s.name) + " " + format_number(best, 4);
                if (s.name.empty()) line += std::string(m_poh ? " (pH " : " (pOH ") + format_number(pkw() - best, 4) + ")";
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
        if (c == "acids") {
            std::string out;
            for (const fc::AcidConstantEntry& e : fc::acid_constants()) out += std::string(e.name) + "\t" + pretty_formula(e.formula) + "  " + e.name + (e.strong ? " (strong)" : "") + "\n";
            for (const fc::BaseConstantEntry& e : fc::base_constants()) out += std::string(e.name) + "\t" + pretty_formula(e.formula) + "  " + e.name + " (base)\n";
            for (const char* s : { "NaOH", "KOH", "Ca(OH)2", "Ba(OH)2" }) out += std::string(s) + "\t" + pretty_formula(s) + " (strong base)\n";
            return out;
        }
        if (c == "readouts") return m_readouts;
        if (c == "scale") {
            m_poh = arg(a, 1) == "poh";
            rebuild();
            return m_poh ? "poh" : "ph";
        }
        if (c == "mode") m_mode = arg(a, 1);
        else if (c == "titration") {
            m_mode = "titration";
            m_analyte = arg(a, 1);
            m_ca = num(arg(a, 2), 0.1);
            m_va = num(arg(a, 3), 0.025);
            m_ct = num(arg(a, 4), 0.1);
            m_vmax = num(arg(a, 5), 0.0);
            m_titrant = arg(a, 6);
            m_flask_role = arg(a, 7);
            m_vadd = arg(a, 8).empty() ? fc::NO_VALUE : num(arg(a, 8), fc::NO_VALUE);
            m_vunit = arg(a, 9).empty() ? "mL" : arg(a, 9);
        } else if (c == "species") {
            m_mode = "species";
            m_analyte = arg(a, 1);
        } else if (c == "mixture") {
            m_mode = "mixture";
            m_mixture = arg(a, 1);
            m_t = num(arg(a, 2), 298.15);
        } else if (c == "mix") {
            m_mode = "mix";
            m_mix_spec = arg(a, 1);
            m_mix_water = std::max(0.0, num(arg(a, 2), 0.0));
            m_t = num(arg(a, 3), 298.15);
        } else if (c == "buffer") {
            m_mode = "buffer";
            m_analyte = arg(a, 1);
            m_ca = num(arg(a, 2), 0.1);
            m_cb = num(arg(a, 3), 0.1);
            m_added = num(arg(a, 4), 0.0);
        } else {
            return "error\tunknown command " + c;
        }
        rebuild();
        return m_error.empty() ? "ok\t" + format_number(m_ph, 4) + "\t" + (finite(m_ph) ? format_number(pkw() - m_ph, 4) : std::string()) : "error\t" + m_error;
    }

private:
    void draw_scale(Canvas& canvas, const Theme& theme) {
        const double x0 = 60.0, x1 = canvas.width() - 60.0, y = canvas.height() * 0.55;
        for (int i = 0; i < 140; ++i) {
            const double p = i / 10.0;
            const double t = p / 14.0;
            const Color c = t < 0.5 ? mix(rgb(230, 70, 70), rgb(120, 200, 90), t * 2.0) : mix(rgb(120, 200, 90), rgb(110, 90, 210), (t - 0.5) * 2.0);
            canvas.fill_rect_px(x0 + (x1 - x0) * i / 140.0, y - 22.0, (x1 - x0) / 140.0 + 1.0, 44.0, c);
        }
        for (int p = 0; p <= 14; ++p) {
            const double x = x0 + (x1 - x0) * p / 14.0;
            canvas.line_px(Vec2{ x, y + 22.0 }, Vec2{ x, y + 30.0 }, theme.text, 1.2);
            canvas.text(Vec2{ x, y + 44.0 }, std::to_string(p), 12.0, theme.text, TextAlign::Center);
        }
        for (int p = 0; p <= 14; ++p) {
            const double x = x0 + (x1 - x0) * (14 - p) / 14.0;
            canvas.text(Vec2{ x, y + 66.0 }, std::to_string(p), 11.0, theme.muted, TextAlign::Center);
        }
        canvas.text(Vec2{ x0 - 18.0, y + 44.0 }, "pH", 12.0, theme.text, TextAlign::Right);
        canvas.text(Vec2{ x0 - 18.0, y + 66.0 }, "pOH", 11.0, theme.muted, TextAlign::Right);
        canvas.text(Vec2{ x0, y - 34.0 }, "acidic", 12.0, theme.red);
        canvas.text(Vec2{ x1, y - 34.0 }, "basic", 12.0, theme.accent, TextAlign::Right);
        canvas.text(Vec2{ 0.5 * (x0 + x1), y - 34.0 }, "neutral", 12.0, theme.green, TextAlign::Center);
        if (finite(m_ph)) {
            const double x = x0 + (x1 - x0) * std::max(0.0, std::min(14.0, m_ph)) / 14.0;
            canvas.fill_triangle(Vec2{ x, y - 26.0 }, Vec2{ x - 9.0, y - 44.0 }, Vec2{ x + 9.0, y - 44.0 }, theme.text);
            canvas.text(Vec2{ x, y - 52.0 }, "pH " + format_number(m_ph, 4) + "  \xC2\xB7  pOH " + format_number(pkw() - m_ph, 4), 16.0, theme.text, TextAlign::Center, 1.2);
            canvas.fill_triangle(Vec2{ x, y + 76.0 }, Vec2{ x - 7.0, y + 90.0 }, Vec2{ x + 7.0, y + 90.0 }, theme.muted);
        }
    }

    double pkw() const {
        return -std::log10(fc::water_ion_product(m_t));
    }

    void to_poh_scale() {
        const double k = pkw();
        if (m_mode == "titration") {
            for (Series& se : m_series) {
                for (auto& pt : se.points) pt.second = k - pt.second;
            }
            for (Marker& mk : m_markers) mk.y = k - mk.y;
            for (auto& band : m_bands) band.second = { k - band.second.second, k - band.second.first };
            m_ylabel = "pOH";
        } else if (m_mode == "species") {
            for (Series& se : m_series) {
                for (auto& pt : se.points) pt.first = k - pt.first;
                std::reverse(se.points.begin(), se.points.end());
            }
            for (Marker& mk : m_markers) {
                mk.x = k - mk.x;
                mk.label = "pOH " + format_number(mk.x, 3);
            }
            m_xlabel = "pOH";
        }
    }

    void rebuild() {
        m_error.clear();
        m_series.clear();
        m_markers.clear();
        m_bands.clear();
        m_ph = fc::NO_VALUE;
        try {
            if (m_mode == "titration") build_titration();
            else if (m_mode == "species") build_species();
            else if (m_mode == "buffer") build_buffer();
            else if (m_mode == "mix") build_mix();
            else build_mixture();
            if (m_poh) to_poh_scale();
        } catch (const std::exception& e) {
            m_error = e.what();
            m_readouts = std::string("#Acids and bases\n!") + e.what() + "\n";
        }
    }

    void build_titration() {
        const int role = m_flask_role == "acid" ? 1 : m_flask_role == "base" ? -1 : 0;
        const fc::AcidBaseSpecies* analyte = fc::find_acid_base_species(m_analyte, role);
        if (!analyte) {
            build_titration_legacy();
            return;
        }
        const fc::AcidBaseSpecies* titrant = nullptr;
        if (!m_titrant.empty()) {
            titrant = fc::find_acid_base_species(m_titrant, analyte->is_acid() ? -1 : 1);
            if (!titrant) {
                const fc::AcidBaseSpecies* other = fc::find_acid_base_species(m_titrant, 0);
                if (other && other->role == analyte->role) throw std::invalid_argument(std::string("both are ") + (analyte->is_acid() ? "acids" : "bases") + ": titrate an acid with a base or a base with an acid");
                throw std::invalid_argument("'" + m_titrant + "' isn't a" + (analyte->is_acid() ? " base" : "n acid") + " I know");
            }
        } else {
            titrant = fc::find_acid_base_species(analyte->is_acid() ? "sodium hydroxide" : "hydrochloric acid", analyte->is_acid() ? -1 : 1);
        }
        const double kw = fc::water_ion_product(m_t);
        double vmax = m_vmax;
        if (!(vmax > 0.0) && finite(m_vadd) && m_vadd > 0.0) {
            const double last = m_ca * m_va * std::max(1, analyte->equivalents) / (m_ct * std::max(1, titrant->equivalents));
            vmax = std::max(1.5 * last, 1.1 * m_vadd);
        }
        const fc::TitrationResult r = fc::titrate(*analyte, m_ca, m_va, *titrant, m_ct, vmax, 500, kw);
        if (!r.ok) throw std::invalid_argument(r.error);
        Series curve{ "", palette(0), {}, false, false };
        for (const auto& pt : r.curve) curve.points.emplace_back(vshow(pt.first), pt.second);
        m_series.push_back(curve);
        m_xmin = 0.0;
        m_xmax = vshow(r.curve.back().first);
        m_ymin = 0.0;
        m_ymax = 14.0;
        std::ostringstream o;
        o << "#" << analyte->name << " titrated with " << titrant->name << "\n";
        o << "Kind\t" << r.kind << "\n";
        const auto describe = [&](const fc::AcidBaseSpecies& sp) {
            std::string d = sp.kind;
            if (finite(sp.pk)) d += sp.is_acid() ? ", pKa " + format_number(sp.pk, 4) : ", pKb " + format_number(sp.pk, 4);
            return d;
        };
        o << "In the flask\t" << analyte->name << " (" << describe(*analyte) << ")\n";
        o << "Titrant\t" << titrant->name << " (" << describe(*titrant) << ")\n";
        m_ph = r.start_ph;
        o << "Starting pH\t" << format_number(r.start_ph, 4) << "\n";
        o << "Starting pOH\t" << format_number(pkw() - r.start_ph, 4) << "\n";
        double first_eq_ph = fc::NO_VALUE;
        for (std::size_t k = 0; k < r.equivalence_volumes.size(); ++k) {
            const double v = r.equivalence_volumes[k];
            const double phk = fc::titration_point_ph(*analyte, m_ca, m_va, *titrant, m_ct, v, kw);
            if (k == 0) first_eq_ph = phk;
            m_markers.push_back(Marker{ vshow(v), phk, "equivalence " + std::to_string(k + 1), palette(4) });
            o << "Equivalence point " << (k + 1) << "\t" << format_number(vshow(v), 4) << " " << m_vunit << " at pH " << format_number(phk, 4) << " (pOH " << format_number(pkw() - phk, 4) << ")\n";
        }
        for (std::size_t k = 0; k < r.half_equivalence_volumes.size(); ++k) {
            const double v = r.half_equivalence_volumes[k];
            if (vshow(v) > m_xmax) continue;
            const double phh = fc::titration_point_ph(*analyte, m_ca, m_va, *titrant, m_ct, v, kw);
            const std::string what = analyte->is_acid() ? "pH \xE2\x89\x88 pKa" : "pOH \xE2\x89\x88 pKb";
            m_markers.push_back(Marker{ vshow(v), phh, "half-equivalence: " + what, palette(1) });
            o << "Half-equivalence " << (k + 1) << "\tpH " << format_number(phh, 4) << ", pOH " << format_number(pkw() - phh, 4) << "\n";
            m_bands.push_back({ palette(1), { phh - 1.0, phh + 1.0 } });
        }
        double steepest = 0.0;
        for (std::size_t i = 1; i < r.curve.size(); ++i) {
            const double dv = vshow(r.curve[i].first) - vshow(r.curve[i - 1].first);
            if (dv > 0) steepest = std::max(steepest, std::fabs(r.curve[i].second - r.curve[i - 1].second) / dv);
        }
        o << "Steepest rise\t" << format_number(steepest, 3) << " pH units per " << m_vunit << "\n";
        o << "#Amounts\n";
        const double n_analyte = m_ca * m_va;
        o << "In the flask\t" << format_number(vshow(m_va), 4) << " " << m_vunit << " of " << format_number(m_ca, 4) << " M = " << format_number(n_analyte * 1000.0, 4) << " mmol " << analyte->name << "\n";
        o << "Titrant\t" << format_number(m_ct, 4) << " M " << titrant->name << " (" << format_number(m_ct * 1.0, 4) << " mmol per mL)\n";
        if (finite(m_vadd) && m_vadd >= 0.0) {
            const double v = m_vadd;
            const double phv = fc::titration_point_ph(*analyte, m_ca, m_va, *titrant, m_ct, v, kw);
            const double n_tit = m_ct * v * std::max(1, titrant->equivalents);
            const double n_need = n_analyte * std::max(1, analyte->equivalents);
            m_markers.push_back(Marker{ vshow(v), phv, "now", theme_marker() });
            m_ph = phv;
            o << "#After adding " << format_number(vshow(v), 4) << " " << m_vunit << "\n";
            o << "pH\t" << format_number(phv, 4) << "\n";
            o << "pOH\t" << format_number(pkw() - phv, 4) << "\n";
            o << "Total volume\t" << format_number(vshow(m_va + v), 4) << " " << m_vunit << "\n";
            o << "Titrant added\t" << format_number(m_ct * v * 1000.0, 4) << " mmol (" << format_number(n_tit * 1000.0, 4) << " mmol " << (titrant->is_base() ? "of base capacity" : "H\xE2\x81\xBA") << ")\n";
            o << "Neutralized\t" << format_number(std::min(100.0, 100.0 * n_tit / n_need), 4) << "% of the way to the last equivalence point\n";
            const double left = n_need - n_tit;
            if (std::fabs(left) < 1e-12) o << "Exactly at the last equivalence point\n";
            else if (left > 0) o << (analyte->is_acid() ? "Acid" : "Base") << " still to neutralize\t" << format_number(left * 1000.0, 4) << " mmol (" << format_number(vshow(left / (m_ct * std::max(1, titrant->equivalents))), 4) << " " << m_vunit << " more titrant)\n";
            else o << "Excess titrant\t" << format_number(-left * 1000.0, 4) << " mmol\n";
        }
        o << "#Indicator\n";
        struct Ind {
            const char* name;
            double lo, hi;
        };
        static const Ind inds[] = { { "methyl orange", 3.1, 4.4 }, { "methyl red", 4.4, 6.2 }, { "bromothymol blue", 6.0, 7.6 }, { "phenol red", 6.8, 8.4 }, { "phenolphthalein", 8.2, 10.0 }, { "thymolphthalein", 9.3, 10.5 } };
        const Ind* best = &inds[0];
        double bd = 1e300;
        for (const Ind& i : inds) {
            const double d = std::fabs(0.5 * (i.lo + i.hi) - first_eq_ph);
            if (d < bd) {
                bd = d;
                best = &i;
            }
        }
        o << "Best indicator\t" << best->name << " (changes " << format_number(best->lo, 2) << "\xE2\x80\x93" << format_number(best->hi, 2) << ")\n";
        const bool sa = analyte->strong, st = titrant->strong;
        if (sa && st) o << "Strong with strong: the salt left at the equivalence point doesn't react with water, so it's neutral (pH 7 at 25 \xC2\xB0" "C) and the jump is big.\n";
        else if (sa != st) {
            const bool weak_acid_side = (analyte->is_acid() && !sa) || (titrant->is_acid() && !st);
            o << (weak_acid_side ? "A weak acid neutralized by a strong base leaves its conjugate base, so the equivalence point is basic." : "A weak base neutralized by a strong acid leaves its conjugate acid, so the equivalence point is acidic.") << "\n";
            if (!sa) o << "Before the equivalence point the flask is a buffer; at half-equivalence " << (analyte->is_acid() ? "pH = pKa" : "pOH = pKb") << ".\n";
        } else {
            o << "Weak with weak: both conjugates react with water, so the equivalence pH is about \xC2\xBD(pKw + pKa \xE2\x88\x92 pKb) and the jump is small and gradual. Indicators can't show it well; a pH meter or conductivity works better.\n";
            o << "!no sharp end point\n";
        }
        m_title = analyte->name + " + " + titrant->name;
        m_xlabel = "titrant added (" + m_vunit + ")";
        m_ylabel = "pH";
        m_readouts = o.str();
    }

    void build_titration_legacy() {
        AcidSpec s;
        if (!acid_spec(m_analyte, s)) throw std::invalid_argument("pick an acid or base (or type pKa=4.76)");
        const int protons = static_cast<int>(s.pka.size());
        const double moles = m_ca * m_va;
        const double veq = fc::equivalence_volume(moles, m_ct, protons);
        const double vmax = m_vmax > 0 ? m_vmax : 1.5 * veq;
        Series curve{ "", palette(0), {}, false, false };
        for (int i = 0; i <= 500; ++i) {
            const double v = vmax * i / 500.0;
            curve.points.emplace_back(vshow(v), fc::titration_ph(moles, s.pka, m_ct * v, m_va + v, !s.base));
        }
        m_series.push_back(curve);
        m_xmin = 0.0;
        m_xmax = vshow(vmax);
        m_ymin = 0.0;
        m_ymax = 14.0;
        const UnitPrefs& u = m_shared.units;
        std::ostringstream o;
        o << "#" << s.name << " titrated with a strong " << (s.base ? "acid" : "base") << "\n";
        const double ph0 = fc::titration_ph(moles, s.pka, 0.0, m_va, !s.base);
        m_ph = ph0;
        o << "Starting pH\t" << format_number(ph0, 4) << "\n";
        o << "Starting pOH\t" << format_number(pkw() - ph0, 4) << "\n";
        for (std::size_t i = 0; i < s.pka.size() && !s.strong; ++i) {
            const std::string idx = s.pka.size() > 1 ? std::to_string(i + 1) : std::string();
            if (s.base) o << "pKb" << idx << "\t" << format_number(pkw() - s.pka[i], 4) << " (Kb " << format_number(std::pow(10.0, s.pka[i] - pkw()), 3) << "; pKa of BH\xE2\x81\xBA " << format_number(s.pka[i], 4) << ")\n";
            else o << "pKa" << idx << "\t" << format_number(s.pka[i], 4) << " (Ka " << format_number(std::pow(10.0, -s.pka[i]), 3) << "; pKb of its base " << format_number(pkw() - s.pka[i], 4) << ")\n";
        }
        for (int k = 1; k <= protons; ++k) {
            const double v = veq * k / protons;
            const double phk = fc::titration_ph(moles, s.pka, m_ct * v, m_va + v, !s.base);
            m_markers.push_back(Marker{ vshow(v), phk, "equivalence " + std::to_string(k), palette(4) });
            o << "Equivalence point " << k << "\t" << u.format(v * 1e-3, "volume", 4) << " at pH " << format_number(phk, 4) << " (pOH " << format_number(pkw() - phk, 4) << ")\n";
            if (!s.strong) {
                const double vh = veq * (k - 0.5) / protons;
                const double phh = fc::titration_ph(moles, s.pka, m_ct * vh, m_va + vh, !s.base);
                m_markers.push_back(Marker{ vshow(vh), phh, "half-equivalence: pH = pK" + std::string(s.base ? "a of BH\xE2\x81\xBA" : "a"), palette(1) });
                o << "Half-equivalence " << k << "\tpH " << format_number(phh, 4) << ", pOH " << format_number(pkw() - phh, 4) << " \xE2\x89\x88 pKa " << format_number(s.pka[static_cast<std::size_t>(k - 1)], 4) << "\n";
                m_bands.push_back({ palette(1), { s.pka[static_cast<std::size_t>(k - 1)] - 1.0, s.pka[static_cast<std::size_t>(k - 1)] + 1.0 } });
            }
        }
        const double ph_eq = fc::titration_ph(moles, s.pka, m_ct * veq, m_va + veq, !s.base);
        o << "#Indicator\n";
        struct Ind {
            const char* name;
            double lo, hi;
        };
        static const Ind inds[] = { { "methyl orange", 3.1, 4.4 }, { "methyl red", 4.4, 6.2 }, { "bromothymol blue", 6.0, 7.6 }, { "phenol red", 6.8, 8.4 }, { "phenolphthalein", 8.2, 10.0 }, { "thymolphthalein", 9.3, 10.5 } };
        const Ind* best = &inds[0];
        double bd = 1e300;
        for (const Ind& i : inds) {
            const double d = std::fabs(0.5 * (i.lo + i.hi) - ph_eq);
            if (d < bd) {
                bd = d;
                best = &i;
            }
        }
        o << "Best indicator\t" << best->name << " (changes " << format_number(best->lo, 2) << "\xE2\x80\x93" << format_number(best->hi, 2) << ")\n";
        o << (s.strong ? "A strong acid with a strong base gives a neutral equivalence point (pH 7)." : s.base ? "A weak base gives an acidic equivalence point, because BH\xE2\x81\xBA is a weak acid." : "A weak acid gives a basic equivalence point, because its conjugate base takes protons from water.") << "\n";
        m_title = s.name + " titration";
        m_xlabel = "titrant added (" + m_vunit + ")";
        m_ylabel = "pH";
        m_readouts = o.str();
    }

    void build_species() {
        AcidSpec s;
        if (!acid_spec(m_analyte, s)) throw std::invalid_argument("pick an acid (or type pKa=2.15, 7.20, 12.35)");
        const std::size_t n = s.pka.size() + 1;
        std::vector<Series> series(n);
        for (std::size_t i = 0; i < n; ++i) {
            std::string label = s.base ? (i == 0 ? "BH\xE2\x81\xBA" : "B") : (i == 0 ? "H" + std::string(s.pka.size() > 1 ? subscript_digits(std::to_string(s.pka.size())) : "") + "A" : "A" + superscript_charge(-static_cast<int>(i)));
            if (!s.base && i > 0 && i < n - 1) label = "H" + (n - 1 - i > 1 ? subscript_digits(std::to_string(n - 1 - i)) : std::string()) + "A" + superscript_charge(-static_cast<int>(i));
            series[i] = Series{ label, palette(static_cast<int>(i)), {}, false, false };
        }
        for (int k = 0; k <= 280; ++k) {
            const double p = k / 20.0;
            const std::vector<double> fr = fc::species_fractions(p, s.pka);
            for (std::size_t i = 0; i < n && i < fr.size(); ++i) series[i].points.emplace_back(p, fr[i]);
        }
        m_series = series;
        for (double pk : s.pka) m_markers.push_back(Marker{ pk, 0.5, "pKa " + format_number(pk, 3), palette(5) });
        m_xmin = 0.0;
        m_xmax = 14.0;
        m_ymin = 0.0;
        m_ymax = 1.05;
        m_title = s.name + ": which form is present at each pH";
        m_xlabel = "pH";
        m_ylabel = "fraction";
        std::ostringstream o;
        o << "#" << s.name << "\n";
        for (std::size_t i = 0; i < s.pka.size(); ++i) o << "pKa" << (s.pka.size() > 1 ? std::to_string(i + 1) : std::string()) << "\t" << format_number(s.pka[i], 4) << " (Ka " << format_number(std::pow(10.0, -s.pka[i]), 3) << ")\n";
        for (std::size_t i = 0; i < s.pka.size(); ++i) o << "pKb of " << (s.base ? "B" : "the conjugate base") << (s.pka.size() > 1 ? " " + std::to_string(i + 1) : std::string()) << "\t" << format_number(pkw() - s.pka[i], 4) << " (Kb " << format_number(std::pow(10.0, s.pka[i] - pkw()), 3) << ")\n";
        o << "pKa + pKb = pKw\t" << format_number(pkw(), 4) << "\n";
        o << "At pH = pKa the two neighboring forms are 50/50; one pH unit away it's about 91/9.\n";
        m_readouts = o.str();
    }

    void build_mix() {
        struct Part {
            const fc::AcidBaseSpecies* species = nullptr;
            double moles = 0.0;
            double volume = 0.0;
            std::string how;
        };
        std::vector<Part> parts;
        for (const std::string& item : pieces(m_mix_spec, ";")) {
            const std::vector<std::string> f = pieces(item, "|");
            if (f.size() < 4) continue;
            Part p;
            const std::string role = trim(f[0]);
            const std::string name = trim(f[1]);
            p.how = trim(f[2]);
            p.species = fc::find_acid_base_species(name, role == "acid" ? 1 : role == "base" ? -1 : 0);
            if (!p.species) throw std::invalid_argument("'" + name + "' isn't an acid or base I know");
            const double x = num(f[3], 0.0);
            const double y = f.size() > 4 ? num(f[4], 0.0) : 0.0;
            if (p.how == "sol") {
                p.moles = std::max(0.0, x) * std::max(0.0, y);
                p.volume = std::max(0.0, y);
            } else if (p.how == "g") {
                p.moles = std::max(0.0, x) / salt_molar_mass(*p.species);
            } else {
                p.moles = std::max(0.0, x);
            }
            parts.push_back(p);
        }
        if (parts.empty()) throw std::invalid_argument("add at least one acid or base");
        double volume = m_mix_water;
        for (const Part& p : parts) volume += p.volume;
        if (!(volume > 0.0)) throw std::invalid_argument("the mixture needs some volume: give a solution's volume or add water");
        std::vector<fc::AcidBaseComponent> comps;
        double strong_acid = 0.0, strong_base = 0.0, acid_eq = 0.0, base_eq = 0.0;
        for (const Part& p : parts) {
            fc::add_acid_base_species(*p.species, p.moles / volume, comps, strong_acid, strong_base);
            (p.species->is_acid() ? acid_eq : base_eq) += p.moles * std::max(1, p.species->equivalents);
        }
        const double kw = fc::water_ion_product(m_t);
        m_ph = fc::equilibrium_ph(comps, strong_acid, strong_base, kw);
        const double h = std::pow(10.0, -m_ph);
        const UnitPrefs& u = m_shared.units;
        std::ostringstream o;
        o << "#What went in\n";
        for (const Part& p : parts) {
            std::string how = p.how == "sol" ? format_number(p.volume * 1000.0, 4) + " mL of solution" : p.how == "g" ? format_number(p.moles * salt_molar_mass(*p.species), 4) + " g of solid" : "added as moles";
            o << p.species->name << "\t" << format_number(p.moles * 1000.0, 4) << " mmol (" << how << "), now " << format_number(p.moles / volume, 4) << " M\n";
        }
        if (m_mix_water > 0.0) o << "Water\t" << format_number(m_mix_water * 1000.0, 4) << " mL\n";
        o << "Total volume\t" << u.format(volume * 1e-3, "volume", 5) << "\n";
        o << "#Result at " << u.format(m_t, "temperature", 4) << "\n";
        o << "pH\t" << format_number(m_ph, 4) << "\n";
        o << "pOH\t" << format_number(-std::log10(kw) - m_ph, 4) << "\n";
        o << "[H\xE2\x82\x83O\xE2\x81\xBA]\t" << format_number(h, 4) << " M\n";
        o << "[OH\xE2\x81\xBB]\t" << format_number(kw / h, 4) << " M\n";
        if (acid_eq > 0.0 && base_eq > 0.0) {
            const double diff = acid_eq - base_eq;
            if (std::fabs(diff) <= 1e-9 * std::max(acid_eq, base_eq)) o << "Neutralization\tthe acid and base exactly cancel (an equivalence point)\n";
            else o << "Neutralization\t" << format_number(std::fabs(diff) * 1000.0, 4) << " mmol of " << (diff > 0 ? "acid" : "base") << " capacity left over\n";
        }
        bool weak = false, strong_left = false;
        for (const Part& p : parts) weak = weak || !p.species->strong;
        strong_left = (acid_eq > base_eq && strong_acid > strong_base) || (base_eq > acid_eq && strong_base > strong_acid + 1e-12);
        if (weak && !strong_left && acid_eq > 0.0 && base_eq > 0.0 && std::fabs(acid_eq - base_eq) > 1e-12) o << "A weak acid or base with some of its partner is a buffer: it resists changes in pH.\n";
        o << "The pH comes from solving the charge balance with every equilibrium at once, after diluting everything into the total volume.\n";
        m_title = "pH " + format_number(m_ph, 4) + "  \xC2\xB7  pOH " + format_number(-std::log10(kw) - m_ph, 4);
        m_readouts = o.str();
    }

    static double salt_molar_mass(const fc::AcidBaseSpecies& s) {
        double m = 0.0;
        try {
            std::string f = s.formula;
            if (!f.empty() && (f.back() == '+' || f.back() == '-')) f.pop_back();
            m = fc::Formula::parse(f).molar_mass();
        } catch (...) {
        }
        if (s.kind.find("anion") != std::string::npos) m += s.strong_base_per_mole * 22.98977;
        if (s.kind.find("cation") != std::string::npos) m += 35.453;
        return m > 0.0 ? m : 1.0;
    }

    void build_mixture() {
        std::vector<fc::AcidBaseComponent> comps;
        double strong_acid = 0.0, strong_base = 0.0;
        std::ostringstream o;
        o << "#Solution\n";
        for (const std::string& entry : pieces(m_mixture, ";\n")) {
            const auto parts = pieces(entry, " \t");
            if (parts.empty()) continue;
            const double c = parts.size() > 1 ? num(parts.back(), 0.1) : 0.1;
            std::string name = entry.substr(0, entry.rfind(parts.back()));
            if (parts.size() == 1) name = parts[0];
            name = trim(name);
            AcidSpec s;
            if (acid_spec(name, s)) {
                if (s.strong && s.base) {
                    strong_base += c * static_cast<double>(s.pka.size());
                } else if (s.strong) {
                    comps.push_back(fc::AcidBaseComponent{ c, s.pka, 0 });
                } else {
                    comps.push_back(fc::AcidBaseComponent{ c, s.pka, s.base ? 1 : 0 });
                }
                o << s.name << "\t" << format_number(c, 4) << " M" << (s.strong ? " (strong)" : "") << "\n";
                continue;
            }
            const fc::PolyatomicIonEntry* ion = fc::find_polyatomic_ion(name);
            if (ion && ion->charge < 0 && ion->acid && *ion->acid && acid_spec(ion->acid, s)) {
                comps.push_back(fc::AcidBaseComponent{ c, s.pka, 0 });
                strong_base += c * static_cast<double>(-ion->charge);
                o << ion->name << " (as a sodium salt)\t" << format_number(c, 4) << " M\n";
                continue;
            }
            o << "!don't know " << name << "\n";
        }
        const double kw = fc::water_ion_product(m_t);
        m_ph = fc::equilibrium_ph(comps, strong_acid, strong_base, kw);
        const double h = std::pow(10.0, -m_ph);
        o << "#Result at " << m_shared.units.format(m_t, "temperature", 4) << "\n";
        o << "pH\t" << format_number(m_ph, 4) << "\n";
        o << "pOH\t" << format_number(-std::log10(kw) - m_ph, 4) << "\n";
        o << "[H\xE2\x82\x83O\xE2\x81\xBA]\t" << format_number(h, 4) << " M\n";
        o << "[OH\xE2\x81\xBB]\t" << format_number(kw / h, 4) << " M\n";
        o << "K_w\t" << format_number(kw, 4) << " (pK_w " << format_number(-std::log10(kw), 4) << ")\n";
        o << "Neutral pH here\t" << format_number(-0.5 * std::log10(kw), 4) << "\n";
        o << "The pH comes from solving the charge balance with every equilibrium at once, so mixtures and buffers come out right.\n";
        m_title = "pH " + format_number(m_ph, 4) + "  \xC2\xB7  pOH " + format_number(-std::log10(kw) - m_ph, 4);
        o << "pOH comes from [OH\xE2\x81\xBB] the way pH comes from [H\xE2\x82\x83O\xE2\x81\xBA]; they always add up to pK_w.\n";
        m_readouts = o.str();
    }

    void build_buffer() {
        AcidSpec s;
        if (!acid_spec(m_analyte, s) || s.strong) throw std::invalid_argument("buffers need a weak acid (or weak base)");
        const double pka = s.pka[0];
        double acid = m_ca, base = m_cb;
        acid += m_added;
        base -= m_added;
        std::ostringstream o;
        o << "#" << s.name << " buffer\n";
        const double ph_start = fc::henderson_hasselbalch(pka, m_cb, m_ca);
        o << "pKa\t" << format_number(pka, 4) << " (Ka " << format_number(std::pow(10.0, -pka), 3) << ")\n";
        o << "pKb of the base\t" << format_number(pkw() - pka, 4) << " (Kb " << format_number(std::pow(10.0, pka - pkw()), 3) << ")\n";
        o << "pH (Henderson\xE2\x80\x93Hasselbalch)\t" << format_number(ph_start, 4) << "\n";
        if (std::fabs(m_added) > 0) {
            if (acid <= 0 || base <= 0) {
                o << "!the added " << (m_added > 0 ? "acid" : "base") << " used up the buffer\n";
                m_ph = m_added > 0 ? fc::strong_acid_ph(-base) : fc::strong_base_ph(-acid);
            } else {
                m_ph = fc::henderson_hasselbalch(pka, base, acid);
            }
            o << "After adding " << format_number(std::fabs(m_added), 4) << " M strong " << (m_added > 0 ? "acid" : "base") << "\tpH " << format_number(m_ph, 4) << " (change " << format_number(m_ph - ph_start, 3) << ")\n";
            const double water = m_added > 0 ? fc::strong_acid_ph(m_added) : fc::strong_base_ph(-m_added);
            o << "Same amount in pure water\tpH " << format_number(water, 4) << "\n";
        } else {
            m_ph = ph_start;
        }
        const double h = std::pow(10.0, -m_ph);
        o << "pH\t" << format_number(m_ph, 4) << "\n";
        o << "pOH\t" << format_number(pkw() - m_ph, 4) << " (pOH = pKb + log([HA]/[A\xE2\x81\xBB]))\n";
        o << "[H\xE2\x82\x83O\xE2\x81\xBA]\t" << format_number(h, 4) << " M\n";
        o << "[OH\xE2\x81\xBB]\t" << format_number(std::pow(10.0, -(pkw() - m_ph)), 4) << " M\n";
        o << "Buffer capacity\t" << format_number(fc::buffer_capacity(m_ca + m_cb, std::pow(10.0, -pka), h), 4) << " mol/(L\xC2\xB7pH)\n";
        o << "Best range\tpH " << format_number(pka - 1.0, 3) << " to " << format_number(pka + 1.0, 3) << "\n";
        m_title = s.name + " buffer: pH " + format_number(m_ph, 4) + ", pOH " + format_number(pkw() - m_ph, 4);
        m_readouts = o.str();
    }

    ChemShared& m_shared;
    std::string m_mode = "titration";
    std::string m_analyte = "acetic acid";
    std::string m_titrant;
    std::string m_flask_role;
    std::string m_vunit = "mL";
    std::string m_mix_spec;
    double m_mix_water = 0.0;
    double m_vadd = fc::NO_VALUE;

    static Color theme_marker() {
        return rgb(250, 250, 255);
    }

    double vshow(double liters) const {
        return from_si_or(liters * 1e-3, m_vunit, fc::Dimension::Volume);
    }
    double m_ca = 0.1, m_va = 0.025, m_ct = 0.1, m_vmax = 0.0, m_cb = 0.1, m_added = 0.0, m_t = 298.15;
    std::string m_mixture = "acetic acid 0.1; sodium hydroxide 0.0";
    double m_ph = fc::NO_VALUE;
    bool m_poh = false;
    std::vector<Series> m_series;
    std::vector<Marker> m_markers;
    std::vector<std::pair<Color, std::pair<double, double>>> m_bands;
    double m_xmin = 0, m_xmax = 1, m_ymin = 0, m_ymax = 14;
    std::string m_title, m_xlabel, m_ylabel, m_readouts, m_error;
    PlotFrame m_frame;
    Vec2 m_mouse{ -1.0, -1.0 };
};

} // namespace

std::unique_ptr<ChemView> make_kinetics_view(ChemShared& shared) {
    return std::unique_ptr<ChemView>(new KineticsView(shared));
}

std::unique_ptr<ChemView> make_acidbase_view(ChemShared& shared) {
    return std::unique_ptr<ChemView>(new AcidBaseView(shared));
}

} // namespace sapling

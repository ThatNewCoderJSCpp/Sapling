#include "simulations.hpp"

#include <algorithm>
#include <cmath>

namespace sapling {

namespace {

const double k_pi = 3.141592653589793;

num::PdeField constant(double c) {
    return [c](const num::PdeInputs& in, double* out) {
        for (std::size_t i = 0; i < in.count; ++i) out[i] = c;
    };
}

struct Coefficient {
    int term = 0;
    num::PdeField value;
};

num::PdeEquation linear_equation(const std::vector<Coefficient>& parts, const num::PdeField& source = nullptr) {
    num::PdeEquation e;
    e.linear = true;
    e.constant = true;
    for (const Coefficient& p : parts) {
        e.terms.push_back(p.term);
        e.partials.push_back({ p.term, p.value });
    }
    e.residual = [parts, source](const num::PdeInputs& in, double* out) {
        std::vector<double> c(in.count);
        for (std::size_t i = 0; i < in.count; ++i) out[i] = 0.0;
        for (const Coefficient& p : parts) {
            p.value(in, c.data());
            const double* v = in.term(static_cast<std::size_t>(p.term));
            for (std::size_t i = 0; i < in.count; ++i) out[i] += c[i] * v[i];
        }
        if (source) {
            source(in, c.data());
            for (std::size_t i = 0; i < in.count; ++i) out[i] -= c[i];
        }
    };
    return e;
}

num::PdeCondition fixed(int function, int term, double value) {
    num::PdeCondition c;
    c.function = function;
    c.label = "fixed value";
    c.equation = linear_equation({ { term, constant(1.0) } }, constant(value));
    return c;
}

class Tunneling : public Simulation {
public:
    Tunneling() {
        add_parameter("k", "wave number", 2.0, 0.5, 4.0, 0.05);
        add_parameter("V0", "barrier height", 2.5, -3.0, 8.0, 0.05);
        add_parameter("width", "barrier width", 1.0, 0.1, 4.0, 0.05);
        add_parameter("sigma", "packet width", 2.0, 0.5, 5.0, 0.1);
    }

    void reset() override {
        m_time = 0.0;
        const double V0 = parameter("V0");
        const double width = parameter("width");
        const double k = parameter("k");
        const double sigma = parameter("sigma");
        num::PdeField potential = [V0, width](const num::PdeInputs& in, double* out) {
            for (std::size_t i = 0; i < in.count; ++i) out[i] = (in.x()[i] >= 0.0 && in.x()[i] <= width) ? V0 : 0.0;
        };
        num::PdeField minus_potential = [potential](const num::PdeInputs& in, double* out) {
            potential(in, out);
            for (std::size_t i = 0; i < in.count; ++i) out[i] = -out[i];
        };
        num::PdeProblem p;
        p.functions = 2;
        p.terms = { num::PdeTerm(0, 0, 0, 0, 1), num::PdeTerm(1, 2, 0, 0, 0), num::PdeTerm(1, 0, 0, 0, 0), num::PdeTerm(1, 0, 0, 0, 1), num::PdeTerm(0, 2, 0, 0, 0), num::PdeTerm(0, 0, 0, 0, 0) };
        p.equations = { linear_equation({ { 0, constant(1.0) }, { 1, constant(0.5) }, { 2, minus_potential } }), linear_equation({ { 3, constant(1.0) }, { 4, constant(-0.5) }, { 5, potential } }) };
        p.lo[0] = -40.0;
        p.hi[0] = 40.0;
        p.points[0] = 1201;
        for (int s = 0; s < 2; ++s) p.edges[0][s].conditions = { fixed(0, 5, 0.0), fixed(1, 2, 0.0) };
        const double x0 = -12.0;
        p.initial = { { [=](const num::PdeInputs& in, double* out) {
            for (std::size_t i = 0; i < in.count; ++i) {
                const double x = in.x()[i];
                out[i] = std::exp(-(x - x0) * (x - x0) / (4.0 * sigma * sigma)) * std::cos(k * x);
            }
        } }, { [=](const num::PdeInputs& in, double* out) {
            for (std::size_t i = 0; i < in.count; ++i) {
                const double x = in.x()[i];
                out[i] = std::exp(-(x - x0) * (x - x0) / (4.0 * sigma * sigma)) * std::sin(k * x);
            }
        } } };
        p.t0 = 0.0;
        p.t1 = m_end = std::min(40.0, 34.0 / std::max(0.5, k));
        p.steps = 1600;
        p.frames = 401;
        p.method = num::PdeMethod::Trapezoid;
        m_solution = num::PdeSolution();
        std::string err;
        m_ok = num::solve_pde(p, m_solution, err);
        m_error = err;
        m_norm0 = m_ok ? probability(0.0, -1e300, 1e300) : 1.0;
    }

    void step(double dt) override {
        m_time = std::min(m_end, m_time + dt);
    }

    void draw(Canvas& c, const Theme& th) const override {
        c.fit_stretched(-30.0, 30.0, -0.35, 0.75, 0.03);
        c.clear(th.background);
        c.grid(5.0, th.grid);
        c.line(Vec2{ -40.0, 0.0 }, Vec2{ 40.0, 0.0 }, th.axes, 1.2);
        const double V0 = parameter("V0");
        const double width = parameter("width");
        const double vs = 0.6 / std::max(4.0, std::fabs(V0) * 1.2);
        c.rect(Vec2{ 0.0, 0.0 }, Vec2{ width, V0 * vs }, mix(th.background, th.warm, 0.35));
        const double E = 0.5 * parameter("k") * parameter("k");
        c.line(Vec2{ -40.0, E * vs }, Vec2{ 40.0, E * vs }, mix(th.background, th.green, 0.5), 1.0);
        if (!m_ok) return;
        const std::vector<double> u = m_solution.at_time(m_time, 0);
        const std::vector<double> v = m_solution.at_time(m_time, 1);
        std::vector<Vec2> re, density;
        for (std::size_t i = 0; i < u.size(); ++i) {
            const double x = m_solution.coord(0, static_cast<int>(i));
            re.push_back(Vec2{ x, 0.3 * u[i] });
            density.push_back(Vec2{ x, 0.6 * (u[i] * u[i] + v[i] * v[i]) });
        }
        c.polyline(re, mix(th.background, th.accent, 0.45), 1.2);
        c.polyline(density, th.accent, 2.2);
    }

    std::vector<Readout> readouts() const override {
        if (!m_ok) return { { "problem", m_error } };
        const double width = parameter("width");
        const double total = probability(m_time, -1e300, 1e300);
        const double through = probability(m_time, width, 1e300);
        const double back = probability(m_time, -1e300, 0.0);
        const double k = parameter("k");
        const double sigma = parameter("sigma");
        return {
            { "time", format_number(m_time, 4) },
            { "energy", format_number(0.5 * k * k + 1.0 / (8.0 * sigma * sigma), 5) + " (barrier " + format_number(parameter("V0"), 4) + ")" },
            { "got through", format_number(100.0 * through / total, 4) + "%" },
            { "bounced back", format_number(100.0 * back / total, 4) + "%" },
            { "total probability", format_number(total / m_norm0, 8) }
        };
    }

    double suggested_speed() const override { return 2.5; }

    std::string note() const override {
        return "Schr\xC3\xB6" "dinger's equation with \xC4\xA7 = m = 1, solved for \xCF\x88 = u + iv as a pair of real PDEs with Crank\xE2\x80\x93Nicolson steps, which keep the total probability at 1. Bright: |\xCF\x88|\xC2\xB2. Faint: the real part. Green line: the packet's energy.";
    }

private:
    double probability(double t, double from, double to) const {
        std::size_t frame = 0;
        for (std::size_t k = 1; k < m_solution.times.size(); ++k) {
            if (std::fabs(m_solution.times[k] - t) < std::fabs(m_solution.times[frame] - t)) frame = k;
        }
        const std::size_t n = m_solution.points();
        const std::vector<double>& f = m_solution.frames[frame];
        const double h = m_solution.spacing(0);
        double s = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double x = m_solution.coord(0, static_cast<int>(i));
            if (x < from || x > to) continue;
            s += (f[i] * f[i] + f[n + i] * f[n + i]) * h;
        }
        return s;
    }

    num::PdeSolution m_solution;
    bool m_ok = false;
    std::string m_error;
    double m_end = 10.0;
    double m_norm0 = 1.0;
};

class HeatRod : public Simulation {
public:
    HeatRod() {
        add_parameter("alpha", "diffusivity", 0.05, 0.005, 0.5, 0.005, "m\xC2\xB2/s");
        add_parameter("left", "left end", 100.0, -50.0, 300.0, 1.0, "\xC2\xB0" "C");
        add_parameter("right", "right end", 0.0, -50.0, 300.0, 1.0, "\xC2\xB0" "C");
        add_parameter("start", "starting temperature", 20.0, -50.0, 300.0, 1.0, "\xC2\xB0" "C");
    }

    void reset() override {
        m_time = 0.0;
        const double alpha = parameter("alpha");
        const double start = parameter("start");
        num::PdeProblem p;
        p.terms = { num::PdeTerm(0, 0, 0, 0, 1), num::PdeTerm(0, 2, 0, 0, 0), num::PdeTerm(0, 0, 0, 0, 0) };
        p.equations = { linear_equation({ { 0, constant(1.0) }, { 1, constant(-alpha) } }) };
        p.lo[0] = 0.0;
        p.hi[0] = 1.0;
        p.points[0] = 201;
        p.edges[0][0].conditions = { fixed(0, 2, parameter("left")) };
        p.edges[0][1].conditions = { fixed(0, 2, parameter("right")) };
        p.initial = { { constant(start) } };
        p.t1 = m_end = 0.8 / alpha;
        p.steps = 800;
        p.frames = 401;
        p.method = num::PdeMethod::Bdf2;
        m_solution = num::PdeSolution();
        std::string err;
        m_ok = num::solve_pde(p, m_solution, err);
        m_error = err;
    }

    void step(double dt) override {
        m_time = std::min(m_end, m_time + dt);
    }

    void draw(Canvas& c, const Theme& th) const override {
        c.fit_stretched(-0.05, 1.05, -0.25, 1.1, 0.04);
        c.clear(th.background);
        if (!m_ok) return;
        const std::vector<double> u = m_solution.at_time(m_time, 0);
        const double lo = std::min({ parameter("left"), parameter("right"), parameter("start") });
        const double hi = std::max({ parameter("left"), parameter("right"), parameter("start") });
        const double span = std::max(1e-9, hi - lo);
        for (int k = 0; k <= 4; ++k) {
            const double y = 0.2 + 0.8 * k / 4.0;
            c.line(Vec2{ 0.0, y }, Vec2{ 1.0, y }, th.grid, 1.0);
        }
        const double h = m_solution.spacing(0);
        for (std::size_t i = 0; i < u.size(); ++i) {
            const double x = m_solution.coord(0, static_cast<int>(i));
            c.rect(Vec2{ x - 0.5 * h, -0.15 }, Vec2{ x + 0.5 * h, 0.05 }, to_color(num::colormap("turbo", (u[i] - lo) / span)));
        }
        std::vector<Vec2> curve;
        std::vector<Vec2> steady;
        for (std::size_t i = 0; i < u.size(); ++i) {
            const double x = m_solution.coord(0, static_cast<int>(i));
            curve.push_back(Vec2{ x, 0.2 + 0.8 * (u[i] - lo) / span });
            const double s = parameter("left") + (parameter("right") - parameter("left")) * x;
            steady.push_back(Vec2{ x, 0.2 + 0.8 * (s - lo) / span });
        }
        c.polyline(steady, mix(th.background, th.muted, 0.5), 1.2);
        c.polyline(curve, th.warm, 2.4);
    }

    std::vector<Readout> readouts() const override {
        if (!m_ok) return { { "problem", m_error } };
        const std::vector<double> u = m_solution.at_time(m_time, 0);
        const double mid = u[u.size() / 2];
        double mean = 0.0;
        for (double v : u) mean += v;
        mean /= static_cast<double>(u.size());
        const double target = 0.5 * (parameter("left") + parameter("right"));
        return {
            { "time", format_number(m_time, 4) + " s" },
            { "middle", format_number(mid, 5) + " \xC2\xB0" "C" },
            { "average", format_number(mean, 5) + " \xC2\xB0" "C" },
            { "settles to", format_number(target, 5) + " \xC2\xB0" "C in the middle" },
            { "time scale L\xC2\xB2/(\xCF\x80\xC2\xB2\xCE\xB1)", format_number(1.0 / (k_pi * k_pi * parameter("alpha")), 4) + " s" }
        };
    }

    double suggested_speed() const override { return 0.08 / parameter("alpha"); }

    std::string note() const override {
        return "The heat equation u_t = \xCE\xB1u_xx on a 1 m rod, solved with BDF2 steps. The bar shows the temperature as color; the grey line is where it settles.";
    }

private:
    num::PdeSolution m_solution;
    bool m_ok = false;
    std::string m_error;
    double m_end = 10.0;
};

class CylinderFlow : public Simulation {
public:
    CylinderFlow() {
        add_parameter("U", "wind speed", 1.0, 0.1, 3.0, 0.05, "m/s");
        add_parameter("R", "cylinder radius", 1.0, 0.3, 2.0, 0.05, "m");
        add_parameter("gamma", "circulation", -4.0, -15.0, 15.0, 0.1, "m\xC2\xB2/s");
        add_parameter("tracers", "tracers", 120.0, 0.0, 400.0, 10.0);
    }

    void reset() override {
        m_time = 0.0;
        m_tracers.clear();
        const int n = static_cast<int>(parameter("tracers"));
        for (int i = 0; i < n; ++i) m_tracers.push_back(spawn(i, n));
        m_age.assign(m_tracers.size(), 0.0);
    }

    void step(double dt) override {
        if (!(dt > 0.0)) return;
        m_time += dt;
        const num::OdeSystem<double> f = [this](const double&, const std::vector<double>& y, std::vector<double>& d) {
            const Vec2 v = velocity(Vec2{ y[0], y[1] });
            d = { v.x, v.y };
        };
        const int n = static_cast<int>(m_tracers.size());
        for (int i = 0; i < n; ++i) {
            Vec2& p = m_tracers[static_cast<std::size_t>(i)];
            std::vector<double> y = { p.x, p.y }, out;
            const int sub = 4;
            for (int s = 0; s < sub; ++s) {
                num::rk4_step<double>(f, 0.0, y, dt / sub, out);
                y = out;
            }
            p = Vec2{ y[0], y[1] };
            m_age[static_cast<std::size_t>(i)] += dt;
            if (p.x > 9.0 || std::fabs(p.y) > 6.0 || std::hypot(p.x, p.y) < 0.98 * parameter("R") || !std::isfinite(p.x)) {
                p = spawn(i, n);
                p.x = -9.0;
                m_age[static_cast<std::size_t>(i)] = 0.0;
            }
        }
    }

    void draw(Canvas& c, const Theme& th) const override {
        c.fit(-6.0, 6.0, -4.0, 4.0, 0.02);
        c.clear(th.background);
        double x0, x1, y0, y1;
        c.visible(x0, x1, y0, y1);
        const double R = parameter("R");
        const num::ScalarGrid2 g = num::sample_grid([this, R](double x, double y) {
            return std::hypot(x, y) < R ? std::nan("") : stream(Vec2{ x, y });
        }, x0, x1, y0, y1, 241, 161);
        const double U = parameter("U");
        for (int k = -14; k <= 14; ++k) {
            const double level = 0.32 * U * k;
            const std::vector<num::Point2> seg = num::marching_squares(g, level);
            const Color col = k == 0 ? mix(th.background, th.accent, 0.8) : mix(th.background, th.accent, 0.35);
            for (std::size_t i = 0; i + 1 < seg.size(); i += 2) c.line(Vec2{ seg[i].x, seg[i].y }, Vec2{ seg[i + 1].x, seg[i + 1].y }, col, 1.2);
        }
        for (std::size_t i = 0; i < m_tracers.size(); ++i) c.disc(m_tracers[i], 2.2, th.yellow);
        c.disc_world(Vec2{ 0.0, 0.0 }, R, mix(th.background, th.muted, 0.5));
        c.circle(Vec2{ 0.0, 0.0 }, R, th.text, 1.5);
        const double lift = -1.225 * U * parameter("gamma");
        if (std::fabs(lift) > 1e-9) c.arrow(Vec2{ 0.0, 0.0 }, Vec2{ 0.0, std::max(-3.0, std::min(3.0, lift * 0.4)) }, th.green, 2.5);
    }

    std::vector<Readout> readouts() const override {
        const double U = parameter("U"), R = parameter("R"), G = parameter("gamma");
        const double s = G / (4.0 * k_pi * U * R);
        std::string stag;
        if (std::fabs(s) <= 1.0) stag = "on the cylinder at " + format_number(std::asin(s) * 180.0 / k_pi, 4) + "\xC2\xB0 and " + format_number(180.0 - std::asin(s) * 180.0 / k_pi, 4) + "\xC2\xB0";
        else stag = "off the cylinder (one point in the flow)";
        return {
            { "lift per meter (air)", format_number(-1.225 * U * G, 5) + " N/m" },
            { "Kutta\xE2\x80\x93Joukowski", "L = -\xCF\x81U\xCE\x93 (clockwise circulation lifts)" },
            { "fastest speed", format_number(max_speed(), 4) + " m/s" },
            { "stagnation points", stag }
        };
    }

    std::string note() const override {
        return "Ideal (potential) flow: w(z) = U(z + R\xC2\xB2/z) - i\xCE\x93/(2\xCF\x80) ln z. Lines are streamlines from marching squares; yellow tracers are carried by the flow with RK4. The green arrow is the lift.";
    }

private:
    double stream(Vec2 p) const {
        const double U = parameter("U"), R = parameter("R"), G = parameter("gamma");
        const double r2 = p.x * p.x + p.y * p.y;
        return U * p.y * (1.0 - R * R / r2) - G / (4.0 * k_pi) * std::log(r2);
    }

    Vec2 velocity(Vec2 p) const {
        const double U = parameter("U"), R = parameter("R"), G = parameter("gamma");
        const double r2 = p.x * p.x + p.y * p.y;
        if (r2 < 1e-12) return Vec2{ 0.0, 0.0 };
        const double r4 = r2 * r2;
        const double u = U * (1.0 - R * R * (p.x * p.x - p.y * p.y) / r4) - G / (2.0 * k_pi) * p.y / r2;
        const double v = -U * R * R * 2.0 * p.x * p.y / r4 + G / (2.0 * k_pi) * p.x / r2;
        return Vec2{ u, v };
    }

    double max_speed() const {
        double best = 0.0;
        const double R = parameter("R");
        for (int i = 0; i < 360; ++i) {
            const double a = 2.0 * k_pi * i / 360.0;
            best = std::max(best, length(velocity(Vec2{ R * 1.0001 * std::cos(a), R * 1.0001 * std::sin(a) })));
        }
        return best;
    }

    Vec2 spawn(int i, int n) const {
        const double y = -3.8 + 7.6 * (i + 0.5) / std::max(1, n);
        const double x = -6.0 - 3.0 * std::fmod(i * 0.618033988749895, 1.0);
        return Vec2{ x, y };
    }

    std::vector<Vec2> m_tracers;
    std::vector<double> m_age;
};

class FieldLines : public Simulation {
public:
    FieldLines() {
        add_parameter("q1", "left charge", 1.0, -3.0, 3.0, 0.1, "nC");
        add_parameter("q2", "right charge", -1.0, -3.0, 3.0, 0.1, "nC");
        add_parameter("q3", "moving charge", 0.5, -3.0, 3.0, 0.1, "nC");
        add_parameter("lines", "lines per unit charge", 14.0, 4.0, 40.0, 1.0);
        add_parameter("omega", "orbit speed", 0.4, 0.0, 3.0, 0.05, "rad/s");
    }

    void reset() override {
        m_time = 0.0;
        rebuild();
    }

    void step(double dt) override {
        if (!(dt > 0.0) || parameter("omega") == 0.0) return;
        m_time += dt;
        rebuild();
    }

    void draw(Canvas& c, const Theme& th) const override {
        c.fit(-3.0, 3.0, -2.2, 2.2, 0.02);
        c.clear(th.background);
        double x0, x1, y0, y1;
        c.visible(x0, x1, y0, y1);
        const num::ScalarGrid2 g = num::sample_grid([this](double x, double y) { return potential(Vec2{ x, y }); }, x0, x1, y0, y1, 181, 121);
        for (int k = -8; k <= 8; ++k) {
            if (k == 0) continue;
            const double level = 0.6 * k * std::fabs(k) * 0.25;
            const std::vector<num::Point2> seg = num::marching_squares(g, level);
            const Color col = k > 0 ? mix(th.background, th.red, 0.3) : mix(th.background, th.accent, 0.3);
            for (std::size_t i = 0; i + 1 < seg.size(); i += 2) c.line(Vec2{ seg[i].x, seg[i].y }, Vec2{ seg[i + 1].x, seg[i + 1].y }, col, 1.0);
        }
        for (const auto& line : m_lines) {
            c.polyline(line, th.text, 1.3);
            if (line.size() > 30) {
                const std::size_t m = line.size() / 2;
                c.arrow(line[m - 1], line[m + 1], th.text, 1.3);
            }
        }
        for (const Charge& q : m_charges) {
            c.disc(q.at, 8.0 + 3.0 * std::fabs(q.q), q.q > 0 ? th.red : q.q < 0 ? th.accent : th.muted);
        }
    }

    std::vector<Readout> readouts() const override {
        double total = 0.0;
        for (const Charge& q : m_charges) total += q.q;
        const Vec2 e = field(Vec2{ 0.0, 0.0 });
        return {
            { "total charge", format_number(total, 4) + " nC" },
            { "field at the origin", format_number(length(e) * 8.98755, 5) + " V/m" },
            { "potential at the origin", format_number(potential(Vec2{ 0.0, 0.0 }) * 8.98755, 5) + " V" },
            { "lines drawn", std::to_string(m_lines.size()) }
        };
    }

    std::string note() const override {
        return "Positions in meters. Each line follows dr/ds = E/|E| from a positive charge (or into a negative one), traced with an adaptive ODE solver that stops when it reaches a charge. Faint curves are equipotentials.";
    }

private:
    struct Charge {
        Vec2 at;
        double q = 0.0;
    };

    void rebuild() {
        const double a = parameter("omega") * m_time;
        m_charges = { { Vec2{ -1.2, 0.0 }, parameter("q1") }, { Vec2{ 1.2, 0.0 }, parameter("q2") }, { Vec2{ 1.4 * std::cos(a + 1.2), 1.4 * std::sin(a + 1.2) }, parameter("q3") } };
        m_lines.clear();
        double positive = 0.0;
        for (const Charge& q : m_charges) positive += std::max(0.0, q.q);
        const bool outward = positive > 0.0;
        for (const Charge& q : m_charges) {
            if ((outward && q.q <= 0.0) || (!outward && q.q >= 0.0)) continue;
            const int n = std::max(4, static_cast<int>(std::lround(std::fabs(q.q) * parameter("lines"))));
            for (int i = 0; i < n; ++i) {
                const double ang = 2.0 * k_pi * (i + 0.5) / n;
                trace(q.at + Vec2{ 0.08 * std::cos(ang), 0.08 * std::sin(ang) }, outward ? 1.0 : -1.0);
            }
        }
    }

    void trace(Vec2 start, double sign) {
        const num::OdeSystem<double> f = [this, sign](const double&, const std::vector<double>& y, std::vector<double>& d) {
            const Vec2 e = field(Vec2{ y[0], y[1] });
            const double l = std::max(1e-12, length(e));
            d = { sign * e.x / l, sign * e.y / l };
        };
        num::OdeOptions<double> opt;
        opt.relative_tolerance = 1e-6;
        opt.absolute_tolerance = 1e-8;
        opt.max_step = 0.05;
        for (const Charge& q : m_charges) {
            const Vec2 at = q.at;
            opt.events.push_back({ [at](const double&, const std::vector<double>& y) { return std::hypot(y[0] - at.x, y[1] - at.y) - 0.06; }, -1, true });
        }
        opt.events.push_back({ [](const double&, const std::vector<double>& y) { return 6.0 - std::hypot(y[0], y[1]); }, -1, true });
        const num::OdeSolution<double> s = num::solve_adaptive<double>(f, 0.0, { start.x, start.y }, 14.0, opt);
        std::vector<Vec2> line;
        for (const auto& y : s.y) line.push_back(Vec2{ y[0], y[1] });
        if (line.size() > 1) m_lines.push_back(std::move(line));
    }

    Vec2 field(Vec2 p) const {
        Vec2 e{ 0.0, 0.0 };
        for (const Charge& q : m_charges) {
            const Vec2 d = p - q.at;
            const double r2 = d.x * d.x + d.y * d.y;
            const double r = std::sqrt(r2);
            if (r < 1e-9) continue;
            e = e + d * (q.q / (r2 * r));
        }
        return e;
    }

    double potential(Vec2 p) const {
        double v = 0.0;
        for (const Charge& q : m_charges) v += q.q / std::max(1e-3, length(p - q.at));
        return v;
    }

    std::vector<Charge> m_charges;
    std::vector<std::vector<Vec2>> m_lines;
};

} // namespace

std::unique_ptr<Simulation> make_tunneling() {
    return std::make_unique<Tunneling>();
}

std::unique_ptr<Simulation> make_heat_rod() {
    return std::make_unique<HeatRod>();
}

std::unique_ptr<Simulation> make_cylinder_flow() {
    return std::make_unique<CylinderFlow>();
}

std::unique_ptr<Simulation> make_field_lines() {
    return std::make_unique<FieldLines>();
}

} // namespace sapling

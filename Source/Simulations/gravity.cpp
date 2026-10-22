#include "simulations.hpp"

#include <cmath>

namespace sapling {

namespace {

const double k_pi = 3.141592653589793;

class Kepler : public Simulation {
public:
    Kepler() {
        add_parameter("r0", "start distance", 1.0, 0.2, 3.0, 0.05, "AU");
        add_parameter("v0", "start speed", 1.2, 0.0, 2.0, 0.01, "\xC3\x97 circular");
        add_parameter("angle", "start direction", 90.0, 0.0, 180.0, 1.0, "\xC2\xB0");
        add_parameter("substeps", "steps per frame", 40.0, 1.0, 400.0, 1.0);
    }

    void reset() override {
        m_time = 0.0;
        const double r0 = parameter("r0");
        const double v = parameter("v0") / std::sqrt(r0);
        const double a = parameter("angle") * k_pi / 180.0;
        m_q = { r0, 0.0 };
        m_v = { v * std::cos(a), v * std::sin(a) };
        m_e0 = energy();
        m_trail.clear();
        m_trail.push_back(Vec2{ m_q[0], m_q[1] });
    }

    void step(double dt) override {
        if (!(dt > 0.0)) return;
        const num::AccelerationFunction<double> acc = [](const double&, const std::vector<double>& q, const std::vector<double>&, std::vector<double>& a) {
            const double r = std::sqrt(q[0] * q[0] + q[1] * q[1]);
            const double r3 = r * r * r;
            a = { -q[0] / r3, -q[1] / r3 };
        };
        const std::size_t steps = static_cast<std::size_t>(parameter("substeps"));
        const num::OdeSolution<double> s = num::solve_symplectic<double>(acc, m_time, m_q, m_v, m_time + dt, steps, num::SymplecticMethod::Yoshida4);
        if (!s.ok) return;
        const std::vector<double>& y = s.back();
        m_q = { y[0], y[1] };
        m_v = { y[2], y[3] };
        m_time += dt;
        m_trail.push_back(Vec2{ m_q[0], m_q[1] });
        trim_trail(m_trail, 2000);
    }

    void draw(Canvas& c, const Theme& th) const override {
        const double r0 = parameter("r0");
        const double e = energy();
        double reach = 2.5 * r0;
        if (e < 0.0) reach = std::min(8.0, std::max(reach, -1.0 / e * 1.05 + 0.2));
        c.fit(-reach, reach, -reach, reach, 0.04);
        c.clear(th.background);
        c.grid(reach > 4 ? 1.0 : 0.5, th.grid);
        c.fading_polyline(m_trail, th.accent, 2.0);
        c.disc(Vec2{ 0.0, 0.0 }, 11.0, th.yellow);
        const Vec2 p{ m_q[0], m_q[1] };
        c.arrow(p, p + Vec2{ m_v[0], m_v[1] } * 0.3, th.green, 2.0);
        c.disc(p, 6.0, th.accent);
    }

    std::vector<Readout> readouts() const override {
        const double e = energy();
        const double l = m_q[0] * m_v[1] - m_q[1] * m_v[0];
        const double ecc = std::sqrt(std::max(0.0, 1.0 + 2.0 * e * l * l));
        std::vector<Readout> out = {
            { "time", format_number(m_time / (2.0 * k_pi), 5) + " years" },
            { "distance", format_number(std::hypot(m_q[0], m_q[1]), 6) + " AU" },
            { "energy", format_number(e, 10) },
            { "energy drift", format_number(std::fabs(e - m_e0), 3) },
            { "angular momentum", format_number(l, 10) },
            { "eccentricity", format_number(ecc, 6) }
        };
        if (e < 0.0) out.push_back({ "period", format_number(std::pow(-1.0 / (2.0 * e), 1.5), 6) + " years" });
        else out.push_back({ "orbit", "escapes (not bound)" });
        return out;
    }

    double suggested_speed() const override { return 2.0; }

    std::string note() const override {
        return "Distances are in AU with the Sun's GM set to 1, so a circular orbit at 1 AU takes one year. Stepped with Yoshida's 4th-order symplectic method, so the energy never drifts away.";
    }

private:
    double energy() const {
        return 0.5 * (m_v[0] * m_v[0] + m_v[1] * m_v[1]) - 1.0 / std::hypot(m_q[0], m_q[1]);
    }

    std::vector<double> m_q;
    std::vector<double> m_v;
    std::vector<Vec2> m_trail;
    double m_e0 = 0.0;
};

class BlackHole : public Simulation {
public:
    BlackHole() {
        add_parameter("r0", "start radius", 20.0, 3.0, 40.0, 0.1, "M");
        add_parameter("L", "angular momentum", 4.4, 0.0, 10.0, 0.01, "M");
        add_parameter("vr", "start radial speed", 0.0, -0.5, 0.5, 0.01, "c");
        add_parameter("newton", "compare with Newton", 1.0, 0.0, 1.0, 1.0);
    }

    void reset() override {
        m_time = 0.0;
        m_coord_time = 0.0;
        const double r0 = parameter("r0");
        m_gr = { r0, parameter("vr"), 0.0 };
        m_nt = m_gr;
        m_inside = false;
        m_trail.clear();
        m_newton_trail.clear();
        m_trail.push_back(point(m_gr));
        m_newton_trail.push_back(point(m_nt));
        const double L = parameter("L");
        m_E = std::sqrt((1.0 - 2.0 / r0) * (1.0 + L * L / (r0 * r0)) + m_gr[1] * m_gr[1]);
        m_last_pr = m_gr[1];
        m_periapses.clear();
    }

    void step(double dt) override {
        if (!(dt > 0.0) || m_inside) return;
        const double L = parameter("L");
        const num::OdeSystem<double> gr = [L](const double&, const std::vector<double>& y, std::vector<double>& d) {
            const double r = y[0];
            d = { y[1], -1.0 / (r * r) + L * L / (r * r * r) - 3.0 * L * L / (r * r * r * r), L / (r * r) };
        };
        const num::OdeSystem<double> nt = [L](const double&, const std::vector<double>& y, std::vector<double>& d) {
            const double r = y[0];
            d = { y[1], -1.0 / (r * r) + L * L / (r * r * r), L / (r * r) };
        };
        num::OdeOptions<double> opt;
        opt.relative_tolerance = 1e-10;
        opt.absolute_tolerance = 1e-12;
        opt.events.push_back({ [](const double&, const std::vector<double>& y) { return y[0] - 2.0; }, -1, true });
        const num::OdeSolution<double> s = num::solve_adaptive<double>(gr, m_time, m_gr, m_time + dt, opt);
        if (!s.ok && !s.stopped_by_event) return;
        for (std::size_t i = 1; i < s.size(); ++i) {
            const double dtau = s.t[i] - s.t[i - 1];
            const double rm = 0.5 * (s.y[i][0] + s.y[i - 1][0]);
            m_coord_time += dtau * m_E / std::max(1e-9, 1.0 - 2.0 / rm);
            const double pr = s.y[i][1];
            if (m_last_pr < 0.0 && pr >= 0.0) m_periapses.push_back(s.y[i][2]);
            m_last_pr = pr;
        }
        m_gr = s.back();
        if (s.stopped_by_event) m_inside = true;
        num::OdeOptions<double> nopt;
        nopt.relative_tolerance = 1e-10;
        nopt.absolute_tolerance = 1e-12;
        if (parameter("newton") > 0.5 && m_nt[0] > 0.05) {
            const num::OdeSolution<double> n = num::solve_adaptive<double>(nt, m_time, m_nt, m_time + dt, nopt);
            if (n.ok) m_nt = n.back();
        }
        m_time = s.t.back();
        m_trail.push_back(point(m_gr));
        m_newton_trail.push_back(point(m_nt));
        trim_trail(m_trail, 4000);
        trim_trail(m_newton_trail, 4000);
    }

    void draw(Canvas& c, const Theme& th) const override {
        const double reach = std::max(10.0, 1.3 * parameter("r0"));
        c.fit(-reach, reach, -reach, reach, 0.04);
        c.clear(th.background);
        c.grid(reach > 30 ? 10.0 : 5.0, th.grid);
        c.circle(Vec2{ 0.0, 0.0 }, 6.0, mix(th.background, th.green, 0.45), 1.2);
        c.circle(Vec2{ 0.0, 0.0 }, 3.0, mix(th.background, th.warm, 0.6), 1.2);
        if (parameter("newton") > 0.5) c.fading_polyline(m_newton_trail, mix(th.background, th.muted, 0.6), 1.4);
        c.fading_polyline(m_trail, th.accent, 2.0);
        c.disc_world(Vec2{ 0.0, 0.0 }, 2.0, rgb(0, 0, 0));
        c.circle(Vec2{ 0.0, 0.0 }, 2.0, th.red, 1.5);
        if (!m_inside) c.disc(point(m_gr), 5.0, th.accent);
    }

    std::vector<Readout> readouts() const override {
        std::vector<Readout> out = {
            { "proper time \xCF\x84", format_number(m_time, 6) + " M" },
            { "far-away time t", m_inside ? std::string("\xE2\x88\x9E (frozen at the horizon)") : format_number(m_coord_time, 6) + " M" },
            { "radius", format_number(m_gr[0], 6) + " M" },
            { "laps", format_number(m_gr[2] / (2.0 * k_pi), 4) },
            { "energy per mass", format_number(m_E, 8) }
        };
        if (m_periapses.size() >= 2) {
            const double shift = m_periapses.back() - m_periapses[m_periapses.size() - 2] - 2.0 * k_pi;
            out.push_back({ "periapsis shift per lap", format_number(shift * 180.0 / k_pi, 5) + "\xC2\xB0" });
        }
        if (m_inside) out.push_back({ "status", "crossed the event horizon" });
        return out;
    }

    double suggested_speed() const override { return 40.0; }

    std::string note() const override {
        return "Units with G = c = M = 1. The red circle is the horizon (r = 2M), orange the photon sphere (3M), green the innermost stable circular orbit (6M). The grey path is the Newtonian orbit from the same start.";
    }

private:
    static Vec2 point(const std::vector<double>& y) {
        return Vec2{ y[0] * std::cos(y[2]), y[0] * std::sin(y[2]) };
    }

    std::vector<double> m_gr;
    std::vector<double> m_nt;
    std::vector<Vec2> m_trail;
    std::vector<Vec2> m_newton_trail;
    std::vector<double> m_periapses;
    double m_coord_time = 0.0;
    double m_E = 1.0;
    double m_last_pr = 0.0;
    bool m_inside = false;
};

} // namespace

std::unique_ptr<Simulation> make_kepler() {
    return std::make_unique<Kepler>();
}

std::unique_ptr<Simulation> make_black_hole() {
    return std::make_unique<BlackHole>();
}

} // namespace sapling

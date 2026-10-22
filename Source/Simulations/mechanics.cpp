#include "simulations.hpp"

#include <cmath>

namespace sapling {

Color to_color(const num::Rgba8& c) {
    return Color{ c.r, c.g, c.b, c.a };
}

void trim_trail(std::vector<Vec2>& trail, std::size_t most) {
    if (trail.size() > most) trail.erase(trail.begin(), trail.begin() + static_cast<long>(trail.size() - most));
}

namespace {

const double k_pi = 3.141592653589793;

class DoublePendulum : public Simulation {
public:
    DoublePendulum() {
        add_parameter("theta1", "first angle", 120.0, -180.0, 180.0, 1.0, "\xC2\xB0");
        add_parameter("theta2", "second angle", -10.0, -180.0, 180.0, 1.0, "\xC2\xB0");
        add_parameter("l1", "first rod", 1.0, 0.2, 2.0, 0.05, "m");
        add_parameter("l2", "second rod", 1.0, 0.2, 2.0, 0.05, "m");
        add_parameter("m1", "first mass", 1.0, 0.1, 5.0, 0.1, "kg");
        add_parameter("m2", "second mass", 1.0, 0.1, 5.0, 0.1, "kg");
        add_parameter("g", "gravity", 9.81, 0.0, 30.0, 0.01, "m/s\xC2\xB2");
        add_parameter("ghost", "twin starts off by", 0.001, 0.0, 5.0, 0.001, "\xC2\xB0");
    }

    void reset() override {
        m_time = 0.0;
        const double t1 = parameter("theta1") * k_pi / 180.0;
        const double t2 = parameter("theta2") * k_pi / 180.0;
        m_state = { t1, t2, 0.0, 0.0 };
        m_twin = { t1 + parameter("ghost") * k_pi / 180.0, t2, 0.0, 0.0 };
        m_energy0 = energy(m_state);
        m_trail.clear();
        m_twin_trail.clear();
    }

    void step(double dt) override {
        if (!(dt > 0.0)) return;
        num::OdeOptions<double> opt;
        opt.relative_tolerance = 1e-10;
        opt.absolute_tolerance = 1e-12;
        const num::OdeSystem<double> f = [this](const double&, const std::vector<double>& y, std::vector<double>& d) { rhs(y, d); };
        const num::OdeSolution<double> a = num::solve_adaptive<double>(f, m_time, m_state, m_time + dt, opt);
        const num::OdeSolution<double> b = num::solve_adaptive<double>(f, m_time, m_twin, m_time + dt, opt);
        if (!a.ok || !b.ok) return;
        m_state = a.back();
        m_twin = b.back();
        m_time += dt;
        m_trail.push_back(bob2(m_state));
        m_twin_trail.push_back(bob2(m_twin));
        trim_trail(m_trail, 700);
        trim_trail(m_twin_trail, 700);
    }

    void draw(Canvas& c, const Theme& th) const override {
        const double reach = parameter("l1") + parameter("l2");
        c.fit(-reach, reach, -reach, reach, 0.06);
        c.clear(th.background);
        c.grid(0.5, th.grid);
        if (parameter("ghost") > 0.0) {
            c.fading_polyline(m_twin_trail, mix(th.background, th.warm, 0.6), 1.5);
            draw_pendulum(c, m_twin, mix(th.background, th.muted, 0.5), mix(th.background, th.warm, 0.6));
        }
        c.fading_polyline(m_trail, th.accent, 2.0);
        draw_pendulum(c, m_state, th.text, th.accent);
        c.disc(Vec2{ 0.0, 0.0 }, 4.0, th.muted);
    }

    std::vector<Readout> readouts() const override {
        const double e = energy(m_state);
        const Vec2 a = bob2(m_state);
        const Vec2 b = bob2(m_twin);
        return {
            { "time", format_number(m_time, 5) + " s" },
            { "angles", format_number(m_state[0] * 180.0 / k_pi, 4) + "\xC2\xB0, " + format_number(m_state[1] * 180.0 / k_pi, 4) + "\xC2\xB0" },
            { "energy", format_number(e, 8) + " J" },
            { "energy drift", format_number(std::fabs(e - m_energy0), 3) + " J" },
            { "twin's distance", format_number(length(a - b), 4) + " m" }
        };
    }

    std::string note() const override {
        return "Solved with an adaptive Dormand\xE2\x80\x93Prince method. The faint twin starts a tiny bit off; watch how soon the two disagree.";
    }

private:
    void rhs(const std::vector<double>& y, std::vector<double>& d) const {
        const double m1 = parameter("m1"), m2 = parameter("m2"), l1 = parameter("l1"), l2 = parameter("l2"), g = parameter("g");
        const double t1 = y[0], t2 = y[1], w1 = y[2], w2 = y[3];
        const double delta = t1 - t2;
        const double den = 2.0 * m1 + m2 - m2 * std::cos(2.0 * delta);
        d.resize(4);
        d[0] = w1;
        d[1] = w2;
        d[2] = (-g * (2.0 * m1 + m2) * std::sin(t1) - m2 * g * std::sin(t1 - 2.0 * t2) - 2.0 * std::sin(delta) * m2 * (w2 * w2 * l2 + w1 * w1 * l1 * std::cos(delta))) / (l1 * den);
        d[3] = (2.0 * std::sin(delta) * (w1 * w1 * l1 * (m1 + m2) + g * (m1 + m2) * std::cos(t1) + w2 * w2 * l2 * m2 * std::cos(delta))) / (l2 * den);
    }

    double energy(const std::vector<double>& y) const {
        const double m1 = parameter("m1"), m2 = parameter("m2"), l1 = parameter("l1"), l2 = parameter("l2"), g = parameter("g");
        const double kin = 0.5 * m1 * l1 * l1 * y[2] * y[2] + 0.5 * m2 * (l1 * l1 * y[2] * y[2] + l2 * l2 * y[3] * y[3] + 2.0 * l1 * l2 * y[2] * y[3] * std::cos(y[0] - y[1]));
        const double pot = -(m1 + m2) * g * l1 * std::cos(y[0]) - m2 * g * l2 * std::cos(y[1]);
        return kin + pot;
    }

    Vec2 bob1(const std::vector<double>& y) const {
        return Vec2{ parameter("l1") * std::sin(y[0]), -parameter("l1") * std::cos(y[0]) };
    }

    Vec2 bob2(const std::vector<double>& y) const {
        return bob1(y) + Vec2{ parameter("l2") * std::sin(y[1]), -parameter("l2") * std::cos(y[1]) };
    }

    void draw_pendulum(Canvas& c, const std::vector<double>& y, Color rod, Color bob) const {
        const Vec2 a = bob1(y);
        const Vec2 b = bob2(y);
        c.line(Vec2{ 0.0, 0.0 }, a, rod, 3.0);
        c.line(a, b, rod, 3.0);
        c.disc(a, 6.0 + 3.0 * std::sqrt(parameter("m1")), bob);
        c.disc(b, 6.0 + 3.0 * std::sqrt(parameter("m2")), bob);
    }

    std::vector<double> m_state;
    std::vector<double> m_twin;
    std::vector<Vec2> m_trail;
    std::vector<Vec2> m_twin_trail;
    double m_energy0 = 0.0;
};

} // namespace

std::unique_ptr<Simulation> make_double_pendulum() {
    return std::make_unique<DoublePendulum>();
}

} // namespace sapling

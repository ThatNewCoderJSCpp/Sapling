#include "sapling/canvas.hpp"
#include "sapling/simulation.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace sapling;

namespace {

int g_passed = 0;
int g_failed = 0;

void check(bool ok, const std::string& what) {
    if (ok) {
        ++g_passed;
    } else {
        ++g_failed;
        std::printf("FAIL: %s\n", what.c_str());
    }
}

std::string readout(const Simulation& s, const std::string& label) {
    for (const Readout& r : s.readouts()) {
        if (r.label == label) return r.value;
    }
    return std::string();
}

double number(const Simulation& s, const std::string& label) {
    return std::atof(readout(s, label).c_str());
}

void run(Simulation& s, double seconds, double dt = 1.0 / 60.0) {
    for (double t = 0.0; t < seconds; t += dt) s.step(dt);
}

void test_every_simulation_runs() {
    Canvas c;
    c.resize(320, 200);
    Theme th;
    for (const SimulationInfo& info : simulations()) {
        auto s = make_simulation(info.id);
        check(s != nullptr, info.id + " can be made");
        if (!s) continue;
        bool known_section = false;
        for (const Section& sec : sections()) known_section = known_section || sec.id == info.section;
        check(known_section, info.id + " is in a known section");
        for (int i = 0; i < 20; ++i) s->step(s->suggested_speed() / 60.0);
        s->draw(c, th);
        check(!s->readouts().empty(), info.id + " has readouts");
        bool painted = false;
        for (std::size_t k = 0; k + 3 < c.pixels().size(); k += 4) painted = painted || c.pixels()[k] != th.background.r;
        check(painted, info.id + " draws something");
        for (const Parameter& p : s->parameters()) {
            check(p.min <= p.value && p.value <= p.max, info.id + ": " + p.id + " starts inside its range");
        }
    }
}

void test_physics() {
    auto pendulum = make_simulation("double_pendulum");
    run(*pendulum, 10.0);
    check(number(*pendulum, "energy drift") < 1e-6, "the double pendulum keeps its energy: " + readout(*pendulum, "energy drift"));

    auto kepler = make_simulation("kepler");
    kepler->set_parameter("v0", 1.0);
    for (int i = 0; i < 100; ++i) kepler->step(2.0 * 3.141592653589793 / 100.0);
    check(std::fabs(number(*kepler, "distance") - 1.0) < 1e-9, "a circular orbit stays at 1 AU: " + readout(*kepler, "distance"));
    check(std::fabs(number(*kepler, "time") - 1.0) < 1e-9, "one lap takes a year: " + readout(*kepler, "time"));
    check(number(*kepler, "energy drift") < 1e-10, "the symplectic orbit keeps its energy");

    auto hole = make_simulation("black_hole");
    hole->set_parameter("r0", 20.0);
    hole->set_parameter("L", 4.6);
    run(*hole, 3000.0, 2.0);
    const std::string shift = readout(*hole, "periapsis shift per lap");
    check(!shift.empty() && std::atof(shift.c_str()) > 5.0, "the black hole orbit's periapsis moves forward: " + shift);
    auto plunge = make_simulation("black_hole");
    plunge->set_parameter("r0", 10.0);
    plunge->set_parameter("L", 3.0);
    run(*plunge, 500.0, 1.0);
    check(readout(*plunge, "status") == "crossed the event horizon", "a low angular momentum orbit plunges in");

    auto tunnel = make_simulation("tunneling");
    run(*tunnel, 12.0, 0.1);
    check(std::fabs(number(*tunnel, "total probability") - 1.0) < 1e-6, "tunneling keeps the total probability at 1: " + readout(*tunnel, "total probability"));
    const double through = number(*tunnel, "got through");
    check(through > 0.5 && through < 60.0, "part of the packet tunnels: " + readout(*tunnel, "got through"));

    auto heat = make_simulation("heat_rod");
    run(*heat, 1000.0, 1.0);
    check(std::fabs(number(*heat, "middle") - 50.0) < 0.5, "the rod settles to 50 in the middle: " + readout(*heat, "middle"));

    auto flow = make_simulation("cylinder_flow");
    check(std::fabs(number(*flow, "lift per meter (air)") - 4.9) < 1e-9, "Kutta-Joukowski lift: " + readout(*flow, "lift per meter (air)"));

    auto field = make_simulation("field_lines");
    check(std::atoi(readout(*field, "lines drawn").c_str()) > 10, "field lines are traced: " + readout(*field, "lines drawn"));
}

} // namespace

int main() {
    test_every_simulation_runs();
    test_physics();
    std::printf("%d passed, %d failed\n", g_passed, g_failed);
    return g_failed == 0 ? 0 : 1;
}

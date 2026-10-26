#ifndef SAPLING_SIMULATION_HPP
#define SAPLING_SIMULATION_HPP

#include "sapling/canvas.hpp"

#include <string>
#include <vector>

namespace sapling {

struct Parameter {
    std::string id;
    std::string label;
    double value = 0.0;
    double min = 0.0;
    double max = 1.0;
    double step = 0.01;
    std::string unit;
};

struct Readout {
    std::string label;
    std::string value;
};

struct Theme {
    Color background = rgb(30, 30, 46);
    Color grid = rgb(49, 50, 68);
    Color axes = rgb(88, 91, 112);
    Color text = rgb(205, 214, 244);
    Color accent = rgb(137, 180, 250);
    Color warm = rgb(250, 179, 135);
    Color red = rgb(243, 139, 168);
    Color green = rgb(166, 227, 161);
    Color yellow = rgb(249, 226, 175);
    Color muted = rgb(147, 153, 178);
};

class Simulation {
public:
    virtual ~Simulation() = default;

    virtual void reset() = 0;
    virtual void step(double dt) = 0;
    virtual void draw(Canvas& canvas, const Theme& theme) const = 0;
    virtual std::vector<Readout> readouts() const = 0;
    virtual double suggested_speed() const { return 1.0; }
    virtual std::string note() const { return std::string(); }

    const std::vector<Parameter>& parameters() const { return m_params; }
    bool set_parameter(const std::string& id, double value);
    double parameter(const std::string& id) const;
    double time() const { return m_time; }

protected:
    void add_parameter(const std::string& id, const std::string& label, double value, double min, double max, double step, const std::string& unit = std::string());
    virtual void parameter_changed(const std::string& id);

    std::vector<Parameter> m_params;
    double m_time = 0.0;
};

std::string format_number(double v, int digits = 6);

} // namespace sapling

#endif // SAPLING_SIMULATION_HPP

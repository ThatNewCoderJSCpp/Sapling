#ifndef SAPLING_SIMULATIONS_HPP
#define SAPLING_SIMULATIONS_HPP

#include "sapling/simulation.hpp"

#ifndef FIZMO
#define FIZMO
#endif
#ifndef FIZMO_NUMERICS
#define FIZMO_NUMERICS
#endif
#include <fizmo/includes.hpp>

#include <memory>

namespace sapling {

namespace num = fizmo::numerics;

std::unique_ptr<Simulation> make_double_pendulum();
std::unique_ptr<Simulation> make_kepler();
std::unique_ptr<Simulation> make_black_hole();
std::unique_ptr<Simulation> make_tunneling();
std::unique_ptr<Simulation> make_heat_rod();
std::unique_ptr<Simulation> make_cylinder_flow();
std::unique_ptr<Simulation> make_field_lines();

Color to_color(const num::Rgba8& c);
void trim_trail(std::vector<Vec2>& trail, std::size_t most);

} // namespace sapling

#endif // SAPLING_SIMULATIONS_HPP

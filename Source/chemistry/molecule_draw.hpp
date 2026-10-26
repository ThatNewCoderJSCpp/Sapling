#ifndef SAPLING_CHEMISTRY_MOLECULE_DRAW_HPP
#define SAPLING_CHEMISTRY_MOLECULE_DRAW_HPP

#include "sapling/chemistry.hpp"

#include <fizmo/../chemistry.hpp>

#include <string>
#include <vector>

namespace sapling {

namespace fc = fizmo::chemistry;

struct DrawAtom {
    Vec3 p;
    int z = 6;
    int formal_charge = 0;
    int lone_pairs = 0;
    int unpaired = 0;
    double partial = 0.0;
    bool selected = false;
    bool highlight = false;
    bool problem = false;
};

struct DrawBond {
    int a = 0;
    int b = 0;
    double order = 1.0;
    bool selected = false;
    bool preview = false;
};

struct DrawOptions {
    bool flat = false;
    bool lone_pairs = true;
    bool formal_charges = true;
    bool partial_charges = false;
    bool bond_dipoles = false;
    bool molecular_dipole = false;
    bool labels = true;
    double atom_scale = 1.0;
    Vec3 dipole;
    double dipole_debye = 0.0;
};

double draw_radius_pm(int z, bool flat);
void draw_molecule(Canvas& canvas, const Camera3D& camera, const Theme& theme, const std::vector<DrawAtom>& atoms, const std::vector<DrawBond>& bonds, const DrawOptions& options);
int pick_atom(const Camera3D& camera, const std::vector<DrawAtom>& atoms, Vec2 screen, bool flat, double atom_scale);
int pick_bond(const Camera3D& camera, const std::vector<DrawAtom>& atoms, const std::vector<DrawBond>& bonds, Vec2 screen);
std::vector<double> lone_pair_angles(const std::vector<double>& bond_angles, int count);
double bond_dipole_debye(int za, int zb, double length_pm);
Vec3 aligned_vector(const std::vector<Vec3>& from, const std::vector<Vec3>& to, Vec3 v);
void draw_hud(Canvas& canvas, const Theme& theme, const std::vector<std::string>& lines, bool bottom);

} // namespace sapling

#endif // SAPLING_CHEMISTRY_MOLECULE_DRAW_HPP

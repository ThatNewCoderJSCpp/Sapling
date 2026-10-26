#ifndef SAPLING_CHEMISTRY_BONDING_VIEW_HPP
#define SAPLING_CHEMISTRY_BONDING_VIEW_HPP

#include "molecule_draw.hpp"

#include <cstdint>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace sapling {

enum class BondingMode { Covalent, Ionic, Metallic };
enum class SandboxTool { Atom, Bond, Move, Erase, Inspect };

struct SandboxAtom {
    int z = 6;
    Vec3 p;
    int neutrons = 0;
};

int sandbox_mass_number(const SandboxAtom& a);

struct SandboxBond {
    int a = 0;
    int b = 0;
    int order = 1;
};

struct CovalentState {
    std::vector<SandboxAtom> atoms;
    std::vector<SandboxBond> bonds;
    int charge = 0;
    std::string name;
};

struct BondReport {
    double avg_order = 1.0;
    double table_length = 0.0;
    double drawn_length = 0.0;
    double energy = 0.0;
    fc::BondThermo thermo;
    double dipole = 0.0;
    double en_difference = 0.0;
    double ionic_percent = 0.0;
};

struct CovalentAnalysis {
    bool ok = false;
    std::string error;
    fc::Molecule molecule;
    fc::LewisStructure lewis;
    fc::LewisStructure best;
    bool better_exists = false;
    fc::ResonanceSet resonance;
    std::vector<double> avg_orders;
    std::vector<double> partial;
    std::vector<int> lone_pairs;
    std::vector<int> unpaired;
    std::vector<int> formal;
    std::vector<bool> problem;
    double mu_shape = 0.0;
    double mu_drawn = 0.0;
    Vec3 dipole_vector;
    bool polar = false;
    std::vector<Vec3> relaxed;
    std::vector<BondReport> bonds;
    double intact_fraction = 1.0;
    std::string verdict;
    std::vector<std::string> issues;
};

struct LatticeAtom {
    Vec3 site;
    Vec3 pos;
    Vec3 start;
    Vec3 vel;
    int species = 0;
    double phase = 0.0;
    bool interstitial = false;
};

class BondingView : public ChemView {
public:
    explicit BondingView(ChemShared& shared);

    void draw(Canvas& canvas, const Theme& theme) override;
    bool mouse_down(Vec2 p, int button, unsigned mods) override;
    bool mouse_move(Vec2 p, bool down) override;
    bool mouse_up(Vec2 p, int button) override;
    bool wheel(Vec2 p, double dy) override;
    bool double_click(Vec2 p) override;
    bool key(const std::string& key, unsigned mods) override;
    bool step(double dt) override;
    std::string command(const std::vector<std::string>& args) override;

private:
    void changed(bool structure = true);
    void push_undo();
    void ensure_analysis();
    std::vector<DrawAtom> draw_atoms() const;
    std::vector<DrawBond> draw_bonds() const;
    bool flat() const { return !m_3d; }
    double current_scale() const;

    int add_atom(int z, Vec3 p);
    void remove_atom(int index);
    int find_bond(int a, int b) const;
    void connect(int a, int b, bool cycle);
    Vec3 place_point(Vec2 screen, double depth) const;
    double ideal_length(int za, int zb, int order) const;
    void fill_hydrogens();
    void tidy();
    void apply_best_lewis();
    void apply_resonance(int index);
    bool load_smiles(const std::string& smiles, int charge, const std::string& name);
    void frame_view();

    std::string covalent_readouts();
    std::string selection_readouts();
    std::string bond_readout(int bond);
    std::string atom_readout(int atom);

    void set_ionic(const std::string& cation, const std::string& anion);
    void set_metal(const std::string& spec);
    void rebuild_lattice();
    void restart_assembly();
    void step_lattice(double dt);
    std::string lattice_state() const;
    void draw_lattice(Canvas& canvas, const Theme& theme);
    std::string ionic_readouts();
    std::string metallic_readouts();
    std::string lattice_selection();
    int pick_lattice(Vec2 screen) const;
    Color species_color(int species) const;
    double species_radius(int species) const;
    std::string species_label(int species) const;

    std::string to_json() const;
    bool from_json(const std::string& json, std::string* error);
    std::string kind_name() const;
    std::string state_text() const;
    std::string format_t() const;

    ChemShared& m_shared;
    BondingMode m_mode = BondingMode::Covalent;
    SandboxTool m_tool = SandboxTool::Atom;
    bool m_3d = false;
    int m_element = 6;
    std::map<std::string, bool> m_options;
    double m_temperature = 298.15;
    double m_pressure = 101325.0;
    double m_frequency = 60.0;
    std::uint64_t m_version = 1;

    CovalentState m_cov;
    std::vector<CovalentState> m_undo;
    std::vector<CovalentState> m_redo;
    CovalentAnalysis m_analysis;
    bool m_analysis_dirty = true;
    int m_sel_atom = -1;
    int m_sel_bond = -1;
    int m_hover_atom = -1;
    int m_hover_bond = -1;

    Camera3D m_cam;
    bool m_framed = false;
    bool m_down = false;
    int m_button = 0;
    unsigned m_mods = 0;
    Vec2 m_down_at;
    Vec2 m_last;
    Vec2 m_mouse;
    bool m_dragged = false;
    int m_drag_atom = -1;
    int m_drag_bond = -1;
    double m_drag_depth = 0.0;
    bool m_drag_pushed = false;
    std::string m_hint;

    std::string m_cation = "Na+";
    std::string m_anion = "Cl-";
    fc::IonicCrystal m_ionic;
    std::vector<std::pair<int, double>> m_components;
    fc::MetallicCrystal m_metal;
    std::string m_lattice_error;
    int m_cells = 2;
    std::vector<LatticeAtom> m_lattice;
    std::vector<std::pair<std::size_t, std::size_t>> m_neighbors;
    std::vector<LatticeAtom> m_electrons;
    double m_assembly = 1.0;
    int m_lattice_state = 0;
    double m_nn_distance = 250.0;
    double m_time = 0.0;
    int m_sel_lattice = -1;
    Vec3 m_origin;
    Vec3 m_box_lo;
    Vec3 m_box_hi;
    std::mt19937 m_rng{ 7u };
};

} // namespace sapling

#endif // SAPLING_CHEMISTRY_BONDING_VIEW_HPP

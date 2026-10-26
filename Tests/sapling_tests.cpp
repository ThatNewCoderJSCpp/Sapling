#include "sapling/canvas.hpp"
#include "sapling/chemistry.hpp"
#include "sapling/simulation.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

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

class Drift : public Simulation {
public:
    Drift() {
        add_parameter("v", "speed", 2.0, -5.0, 5.0, 0.1, "m/s");
    }

    void reset() override {
        m_time = 0.0;
        m_x = 0.0;
        ++resets;
    }

    void step(double dt) override {
        m_time += dt;
        m_x += parameter("v") * dt;
    }

    void draw(Canvas& canvas, const Theme& theme) const override {
        canvas.clear(theme.background);
        canvas.fit(-10.0, 10.0, -10.0, 10.0);
        canvas.disc_world(Vec2{ m_x, 0.0 }, 1.0, theme.accent);
    }

    std::vector<Readout> readouts() const override {
        return { Readout{ "x", format_number(m_x) } };
    }

    double x() const { return m_x; }
    int resets = 0;

private:
    double m_x = 0.0;
};

bool painted(const Canvas& c, Color background) {
    const auto& px = c.pixels();
    for (std::size_t k = 0; k + 3 < px.size(); k += 4) {
        if (px[k] != background.r || px[k + 1] != background.g || px[k + 2] != background.b) return true;
    }
    return false;
}

void test_canvas() {
    Canvas c;
    c.resize(200, 100);
    check(c.width() == 200 && c.height() == 100, "the canvas resizes");
    check(c.pixels().size() == 200u * 100u * 4u, "the canvas holds RGBA pixels");

    const Color bg = rgb(10, 20, 30);
    c.clear(bg);
    check(!painted(c, bg), "clear fills every pixel");

    c.set_view(0.0, 0.0, 10.0);
    const Vec2 s = c.to_screen(Vec2{ 1.0, 2.0 });
    const Vec2 w = c.to_world(s);
    check(std::fabs(w.x - 1.0) < 1e-9 && std::fabs(w.y - 2.0) < 1e-9, "to_world undoes to_screen");

    c.line(Vec2{ -3.0, -3.0 }, Vec2{ 3.0, 3.0 }, rgb(255, 255, 255));
    check(painted(c, bg), "a line draws something");

    const std::vector<std::uint8_t> taken = c.take_pixels();
    check(taken.size() == 200u * 100u * 4u, "take_pixels hands over the image");
    check(c.restore(taken), "restore takes the pixels back");
    check(painted(c, bg), "the restored image is the same one");

    const Color m = mix(rgb(0, 0, 0), rgb(200, 100, 50), 0.5);
    check(m.r == 100 && m.g == 50 && m.b == 25, "mix blends halfway");
    check(std::fabs(length(Vec2{ 3.0, 4.0 }) - 5.0) < 1e-12, "length of (3, 4) is 5");
}

void test_simulation_base() {
    Drift d;
    d.reset();
    check(d.parameters().size() == 1, "parameters are listed");
    check(d.parameter("v") == 2.0, "a parameter starts at its value");
    check(d.parameter("nope") == 0.0, "an unknown parameter reads as 0");

    for (int i = 0; i < 60; ++i) d.step(1.0 / 60.0);
    check(std::fabs(d.x() - 2.0) < 1e-9, "stepping a second at 2 m/s moves 2 m");
    check(std::fabs(d.time() - 1.0) < 1e-9, "time adds up");

    const int before = d.resets;
    check(d.set_parameter("v", 99.0), "setting a known parameter works");
    check(d.parameter("v") == 5.0, "a parameter is clamped to its range");
    check(d.resets == before + 1, "changing a parameter resets by default");
    check(d.set_parameter("v", 5.0) && d.resets == before + 1, "setting the same value doesn't reset");
    check(!d.set_parameter("v", std::nan("")), "a NaN parameter is refused");
    check(!d.set_parameter("nope", 1.0), "an unknown parameter is refused");

    Canvas c;
    c.resize(64, 64);
    Theme th;
    d.draw(c, th);
    check(painted(c, th.background), "the simulation draws");
    check(!d.readouts().empty() && d.readouts()[0].label == "x", "readouts are listed");
}

void test_format_number() {
    check(format_number(1.5) == "1.5", "format_number prints plainly");
    check(format_number(1e-320) == "0", "tiny numbers print as 0");
    check(format_number(std::nan("")) == "undefined", "NaN prints as undefined");
}

bool contains(const std::string& s, const std::string& part) {
    return s.find(part) != std::string::npos;
}

void set_env(const char* name, const std::string& value) {
#if defined(_WIN32)
    _putenv_s(name, value.c_str());
#else
    setenv(name, value.c_str(), 1);
#endif
}

void test_canvas_extras() {
    Canvas c;
    c.resize(120, 80);
    Theme th;
    c.clear(th.background);
    check(c.text_width("H2O", 12.0) > 10.0, "text has a width");
    check(c.text_width("", 12.0) == 0.0, "empty text has no width");
    c.text(Vec2{ 10.0, 40.0 }, "\xCE\x94H = -286 kJ/mol \xE2\x86\x92 H\xE2\x82\x82O", 14.0, th.text);
    check(painted(c, th.background), "text draws, Greek letters and subscripts included");
    c.clear(th.background);
    c.fill_polygon({ Vec2{ 10, 10 }, Vec2{ 100, 20 }, Vec2{ 50, 70 } }, th.accent);
    check(painted(c, th.background), "polygons fill");
    c.clear(th.background);
    c.sphere_px(Vec2{ 60, 40 }, 20.0, th.warm);
    check(painted(c, th.background), "spheres draw");
}

void test_render_helpers() {
    Camera3D cam;
    cam.width = 200;
    cam.height = 100;
    cam.zoom = 2.0;
    cam.yaw = 0.7;
    cam.pitch = 0.3;
    const Vec3 p{ 30.0, -12.0, 5.0 };
    double depth = 0.0;
    const Vec2 s = cam.project(p, &depth);
    const Vec3 back = cam.unproject(s, depth);
    check(length(back - p) < 1e-9, "unproject undoes project");
    check(superscript_charge(-2) == "\xC2\xB2\xE2\x81\xBB", "charges become superscripts");
    check(subscript_digits("H2O") == "H\xE2\x82\x82O", "digits become subscripts");
    const std::vector<double> t = nice_ticks(0.0, 10.0, 5);
    check(!t.empty() && t.front() == 0.0 && t.back() == 10.0, "ticks cover the range");
    UnitPrefs u;
    u.set("temperature", "C");
    check(u.get("temperature") == "\xC2\xB0" "C", "unit aliases are stored by symbol");
    check(contains(u.format(373.15, "temperature"), "100"), "temperatures convert for display");
}

void test_json_and_storage() {
    JsonValue v;
    check(parse_json("{\"a\": [1, 2.5, \"x\\u00e9\"], \"b\": true, \"c\": null}", v), "JSON parses");
    check(v.find("a") && v.find("a")->array.size() == 3 && v.find("a")->array[1].number == 2.5, "JSON arrays and numbers");
    check(v.find("a")->array[2].string == "x\xC3\xA9", "JSON unicode escapes");
    check(!parse_json("{\"a\": }", v), "bad JSON is refused");
    check(json_quote("a\"b\n") == "\"a\\\"b\\n\"", "JSON quoting escapes");
    check(safe_file_name("../evil/name") == "_evil_name", "file names can't leave the folder");

    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "sapling-tests-save";
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
    set_env("SAPLING_SAVE_DIR", dir.string());
    check(chemistry_save_root() == dir / "chemistry", "SAPLING_SAVE_DIR moves the save folder");
    check(save_chemistry("covalent", "water", "{\"type\": \"covalent\"}"), "saving works");
    check(!save_chemistry("plasma", "x", "{}"), "unknown kinds are refused");
    const std::vector<std::string> names = list_saved_chemistry("covalent");
    check(names.size() == 1 && names[0] == "water", "saved files are listed");
    std::string json;
    check(load_chemistry("covalent", "water", json) && contains(json, "covalent"), "saved files load");
    check(delete_saved_chemistry("covalent", "water") && list_saved_chemistry("covalent").empty(), "saved files delete");
    std::filesystem::remove_all(dir, ec);
}

void test_bonding_view() {
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "sapling-tests-bond";
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
    set_env("SAPLING_SAVE_DIR", dir.string());
    ChemShared shared;
    std::unique_ptr<ChemView> view = make_bonding_view(shared);
    Canvas c;
    c.resize(400, 300);
    Theme th;
    view->draw(c, th);
    check(painted(c, th.background), "the bonding sandbox draws");
    check(contains(view->command({ "readouts" }), "water"), "the sandbox starts with water and names it");
    check(contains(view->command({ "readouts" }), "polar"), "water is analysed for polarity");
    check(contains(view->command({ "load_ion", "nitrate" }), "ok"), "polyatomic ions load as templates");
    const std::string nitrate = view->command({ "readouts" });
    check(contains(nitrate, "1.33"), "nitrate shows its resonance bond order");
    check(contains(nitrate, "nonpolar"), "nitrate is nonpolar");
    check(contains(view->command({ "clear" }), "ok") && contains(view->command({ "readouts" }), "Pick an element"), "clearing empties the sandbox");
    view->command({ "element", "C" });
    view->mouse_down(Vec2{ 200, 150 }, 1, 0);
    view->mouse_up(Vec2{ 200, 150 }, 1);
    view->command({ "fill_h" });
    check(contains(view->command({ "readouts" }), "methane"), "a clicked carbon plus Fill H makes methane");
    check(contains(view->command({ "save", "my methane" }), "ok"), "covalent saves work");
    check(contains(view->command({ "saved" }), "my methane"), "saves are listed by kind");
    view->command({ "clear" });
    check(contains(view->command({ "load", "covalent", "my methane" }), "ok") && contains(view->command({ "readouts" }), "methane"), "saves load back");
    check(contains(view->command({ "mode", "ionic" }), "ionic"), "ionic mode");
    check(contains(view->command({ "ionic", "Mg", "O" }), "MgO"), "magnesium oxide builds");
    check(contains(view->command({ "readouts" }), "Lattice energy"), "ionic readouts show the lattice energy");
    for (int i = 0; i < 50; ++i) view->step(0.05);
    view->draw(c, th);
    check(painted(c, th.background), "the lattice draws");
    check(contains(view->command({ "mode", "metallic" }), "metallic"), "metallic mode");
    check(contains(view->command({ "metal", "Cu=70;Zn=30" }), "brass"), "brass is recognized");
    check(contains(view->command({ "tp", "1500", "K", "1", "atm" }), "1500"), "temperature and pressure are set");
    check(contains(view->command({ "readouts" }), "liquid"), "brass melts at 1500 K");
    std::filesystem::remove_all(dir, ec);
}

void test_other_views_and_tools() {
    ChemShared shared;
    Canvas c;
    c.resize(320, 240);
    Theme th;
    const char* names[] = { "structure", "phase", "orbital", "cloud", "element" };
    std::unique_ptr<ChemView> views[] = { make_structure_view(shared), make_phase_view(shared), make_orbital_view(shared), make_cloud_view(shared), make_element_view(shared) };
    for (int i = 0; i < 5; ++i) {
        views[i]->draw(c, th);
        check(painted(c, th.background), std::string("the ") + names[i] + " view draws");
        check(!views[i]->command({ "readouts" }).empty(), std::string("the ") + names[i] + " view has readouts");
    }
    check(contains(views[1]->command({ "readouts" }), "373"), "water boils at 373 K");
    check(contains(views[1]->command({ "point", "-10", "C", "1", "atm" }), "solid"), "water is solid at -10 C");
    views[1]->command({ "view", "3d" });
    views[1]->draw(c, th);
    check(painted(c, th.background), "the 3D P-v-T surface draws");
    check(contains(chemistry_tool(shared, { "name", "Fe2(SO4)3" }), "iron(III) sulfate"), "the namer names ionic compounds");
    check(contains(chemistry_tool(shared, { "name", "sulfurous acid" }), "H\xE2\x82\x82SO\xE2\x82\x83"), "names turn into formulas");
    check(contains(chemistry_tool(shared, { "balance", "Al + O2", "->", "Al2O3" }), "4 Al + 3 O\xE2\x82\x82 \xE2\x86\x92 2 Al\xE2\x82\x82O\xE2\x82\x83"), "the balancer balances");
    check(contains(chemistry_tool(shared, { "predict", "C2H6 + O2" }), "combustion"), "the predictor finds combustion");
    check(contains(chemistry_tool(shared, { "solve", "molarity", "n=0.5|mol;V=250|mL", "M", "M" }), "2"), "the solver finds molarity");
    check(contains(chemistry_tool(shared, { "convert", "100", "C", "F", "temperature" }), "212"), "units convert");
    check(contains(chemistry_tool(shared, { "stoich", "2 H2 + O2 -> 2 H2O", "0|2|g;1|32|g", "25", "C", "1", "atm", "g" }), "Limiting"), "stoichiometry finds the limiting reactant");
    check(contains(chemistry_tool(shared, { "ions", "", "phosph" }), "phosphate"), "ions can be searched");
}

void test_new_chemistry_views() {
    ChemShared shared;
    Canvas c;
    c.resize(400, 300);
    Theme th;
    std::unique_ptr<ChemView> kin = make_kinetics_view(shared);
    kin->draw(c, th);
    check(painted(c, th.background), "the kinetics view draws");
    check(contains(kin->command({ "rate", "1", "0.1", "1", "0", "conc" }), "ok"), "first-order rate law runs");
    check(contains(kin->command({ "readouts" }), "6.93"), "first-order half-life is ln2/k");
    check(contains(kin->command({ "mechanism", "A -> B : 0.3; B -> C : 0.1", "A = 1", "40" }), "ok"), "mechanisms integrate");
    check(contains(kin->command({ "fit", "0, 1.0; 10, 0.74; 20, 0.55; 30, 0.41; 40, 0.30; 60, 0.17" }), "ok") && contains(kin->command({ "readouts" }), "first order"), "data fitting picks first order");
    check(contains(kin->command({ "boltzmann", "N2", "298", "800", "50000" }), "ok"), "Maxwell-Boltzmann takes a formula");
    check(contains(kin->command({ "arrhenius", "50000", "1e10", "298", "348", "ln" }), "ok"), "Arrhenius runs");
    kin->draw(c, th);
    check(painted(c, th.background), "the Arrhenius plot draws");
    std::unique_ptr<ChemView> ab = make_acidbase_view(shared);
    ab->draw(c, th);
    check(painted(c, th.background), "the acid-base view draws");
    check(contains(ab->command({ "acids" }), "acetic acid"), "acids are listed");
    check(contains(ab->command({ "mixture", "HCl 0.01", "298.15" }), "ok\t2"), "0.01 M HCl has pH 2");
    check(contains(ab->command({ "mixture", "NaOH 0.01", "298.15" }), "ok\t12") || contains(ab->command({ "mixture", "NaOH 0.01", "298.15" }), "ok\t11.99"), "0.01 M NaOH has pH 12");
    check(contains(ab->command({ "buffer", "acetic acid", "0.1", "0.1", "0" }), "ok\t4.7"), "an equimolar acetate buffer sits at the pKa");
    check(contains(ab->command({ "titration", "phosphoric acid", "0.1", "0.025", "0.1", "" }), "ok"), "polyprotic titrations run");
    check(contains(ab->command({ "readouts" }), "quivalence"), "titrations report equivalence points");
    check(contains(ab->command({ "readouts" }), "pOH"), "titrations report pOH");
    check(contains(ab->command({ "titration", "ammonia", "0.1", "0.025", "0.1", "", "hydrochloric acid", "base" }), "ok\t11.1") && contains(ab->command({ "readouts" }), "weak base titrated with a strong acid"), "a weak base in the flask titrated with a strong acid");
    check(contains(ab->command({ "titration", "acetic acid", "0.1", "0.025", "0.1", "", "ammonia", "acid" }), "ok") && contains(ab->command({ "readouts" }), "no sharp end point"), "weak with weak has no sharp end point");
    check(contains(ab->command({ "titration", "carbonate", "0.1", "0.025", "0.1", "", "hydrochloric acid", "base" }), "ok") && contains(ab->command({ "readouts" }), "Equivalence point 2\t50 mL"), "carbonate shows two equivalence points");
    check(contains(ab->command({ "titration", "acetic acid", "0.1", "0.025", "0.1", "", "nitric acid", "acid" }), "both are acids"), "two acids are refused");
    check(contains(ab->command({ "mix", "acid|acetic acid|sol|0.1|0.05;base|sodium hydroxide|sol|0.1|0.025;", "0", "298.15" }), "ok\t4.75"), "half-neutralized acetic acid mixed by volume is at its pKa");
    check(contains(ab->command({ "mix", "base|sodium hydroxide|g|0.39997", "0.1", "298.15" }), "ok\t12.99") && contains(ab->command({ "readouts" }), "Total volume"), "a weighed solid dissolved in water");
    check(contains(ab->command({ "mix", "acid|hydrochloric acid|mol|0.001;base|ammonia|sol|0.1|0.01", "0.09", "298.15" }), "ok") && contains(ab->command({ "readouts" }), "exactly cancel"), "moles and solutions mix to an equivalence point");
    check(contains(ab->command({ "titration", "acetic acid", "0.1", "0.025", "0.1", "", "sodium hydroxide", "acid", "0.0125", "L" }), "ok") && contains(ab->command({ "readouts" }), "After adding 0.0125 L") && contains(ab->command({ "readouts" }), "pH\t4.75") && contains(ab->command({ "readouts" }), "Equivalence point 1\t0.025 L"), "volumes follow the chosen unit and the added amount is reported");
    check(contains(ab->command({ "scale", "poh" }), "poh") && contains(ab->command({ "mixture", "HCl 0.01", "298.15" }), "\t11.99"), "pOH comes back alongside pH");
    ab->command({ "scale", "ph" });
    std::unique_ptr<ChemView> el = make_element_view(shared);
    el->command({ "element", "C", "0", "2" });
    check(contains(el->command({ "readouts" }), "14"), "extra neutrons make carbon-14");
    el->command({ "style", "nucleus" });
    el->draw(c, th);
    check(painted(c, th.background), "the nucleus picture draws");
    el->command({ "element", "Cu", "0", "0" });
    const std::string cu = el->command({ "readouts" });
    check(contains(cu, "Specific heat") || contains(cu, "specific heat"), "element readouts show the specific heat");
    check(contains(cu, "esistivity"), "element readouts show the resistivity");
    check(contains(cu, "[Ar]"), "the noble-gas configuration is shown");
    check(contains(cu, "Work function") && contains(cu, "Young's modulus") && contains(cu, "Fermi energy") && contains(cu, "susceptibility"), "element readouts cover electrons, stiffness and magnetism");
    check(contains(el->command({ "nucleus", "2d" }), "2d") && contains(el->command({ "nucleus", "3d" }), "3d") && contains(el->command({ "nucleus", "circle" }), "circle"), "the nucleus can be one circle, flat or 3D");
    el->command({ "style", "both" });
    el->command({ "nucleus", "2d" });
    el->draw(c, th);
    check(painted(c, th.background), "the shells draw with every nucleon");
    std::unique_ptr<ChemView> ph = make_phase_view(shared);
    ph->command({ "axes", "0", "700", "0", "3e7" });
    check(contains(ph->command({ "axis_values" }), "700"), "the axis range is kept");
    ph->command({ "extreme", "1" });
    check(contains(ph->command({ "readouts" }), "ioniz"), "plasma readouts appear");
    check(!ph->command({ "selection" }).empty(), "the selection panel has text");
    ph->command({ "substance", "H2O2" });
    check(!contains(ph->command({ "props" }), "vaporization\t\t"), "hydrogen peroxide gets an enthalpy of vaporization");
    std::unique_ptr<ChemView> orb = make_orbital_view(shared);
    orb->command({ "element", "Cs" });
    orb->command({ "show", "1s" });
    check(std::strtod(orb->command({ "rmax", "auto" }).c_str(), nullptr) > 0.0, "the radial range resets to auto");
    std::unique_ptr<ChemView> bond = make_bonding_view(shared);
    check(contains(bond->command({ "load_compound", "XeF4" }), "ok"), "XeF4 loads");
    check(!contains(bond->command({ "readouts" }), "No bond energy"), "Xe-F has a bond energy");
    bond->command({ "mode", "ionic" });
    check(contains(bond->command({ "ionic", "ammonium", "sulfate" }), "ok"), "polyatomic ions build ionic lattices");
    check(contains(bond->command({ "ionic", "Cu2+", "O2-" }), "CuO") && contains(bond->command({ "ionic", "Cu+", "Cl-" }), "CuCl"), "metal ions of either charge build lattices");
    check(contains(bond->command({ "readouts" }), "eat capacity") || contains(bond->command({ "readouts" }), "Specific heat"), "ionic readouts show heat capacity");
    bond->command({ "mode", "metallic" });
    bond->command({ "metal", "Cu=70;Zn=30" });
    check(contains(bond->command({ "readouts" }), "Skin depth"), "metal readouts show the skin depth");
    bond->command({ "mode", "covalent" });
    bond->command({ "load_compound", "water" });
    check(contains(bond->command({ "readouts" }), "78.36"), "water's permittivity is shown");
    bond->command({ "mode", "ionic" });
    bond->command({ "ionic", "Na", "Cl" });
    check(contains(bond->command({ "readouts" }), "Relative permittivity") && contains(bond->command({ "readouts" }), "Bulk modulus"), "salt shows its permittivity and stiffness");
}

} // namespace

int main() {
    test_canvas();
    test_simulation_base();
    test_format_number();
    test_canvas_extras();
    test_render_helpers();
    test_json_and_storage();
    test_bonding_view();
    test_other_views_and_tools();
    test_new_chemistry_views();
    std::printf("%d passed, %d failed\n", g_passed, g_failed);
    return g_failed == 0 ? 0 : 1;
}

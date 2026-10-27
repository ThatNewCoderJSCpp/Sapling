#ifndef SAPLING_CHEMISTRY_HPP
#define SAPLING_CHEMISTRY_HPP

#include "sapling/canvas.hpp"
#include "sapling/simulation.hpp"

#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace sapling {

constexpr unsigned k_shift = 1u;
constexpr unsigned k_control = 2u;
constexpr unsigned k_alt = 4u;

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vec3 operator+(Vec3 a, Vec3 b);
Vec3 operator-(Vec3 a, Vec3 b);
Vec3 operator*(Vec3 a, double s);
double dot(Vec3 a, Vec3 b);
Vec3 cross(Vec3 a, Vec3 b);
double length(Vec3 a);
Vec3 normalized(Vec3 a);

struct Camera3D {
    double yaw = 0.55;
    double pitch = 0.35;
    double zoom = 1.0;
    double perspective = 0.0;
    Vec3 target;
    Vec2 pan;
    int width = 1;
    int height = 1;

    Vec3 to_view(Vec3 world) const;
    Vec3 from_view(Vec3 view) const;
    Vec2 project(Vec3 world, double* depth = nullptr, double* scale = nullptr) const;
    Vec3 unproject(Vec2 screen, double view_depth) const;
    void orbit(double dx_px, double dy_px);
    void zoom_by(double factor);
    void frame(const std::vector<Vec3>& points, double padding);
    void face_front();
};

struct SceneItem {
    enum Kind { Sphere, Line, Disc, Ring, Text, Arrow, Dashed, Cylinder };
    Kind kind = Sphere;
    double depth = 0.0;
    Vec2 a;
    Vec2 b;
    double size = 1.0;
    double width = 1.0;
    Color color;
    Color color2;
    std::string text;
    TextAlign align = TextAlign::Center;
};

class Scene3D {
public:
    explicit Scene3D(const Camera3D& camera) : m_camera(camera) {}

    void sphere(Vec3 p, double radius_world, Color c, double depth_bias = 0.0);
    void sphere_px(Vec3 p, double radius_px, Color c, double depth_bias = 0.0);
    void cylinder(Vec3 a, Vec3 b, double width_world, Color ca, Color cb, double depth_bias = 0.0);
    void line(Vec3 a, Vec3 b, double width_px, Color c, double depth_bias = 0.0);
    void dashed(Vec3 a, Vec3 b, double width_px, Color c, double depth_bias = 0.0);
    void arrow(Vec3 a, Vec3 b, double width_px, Color c, double depth_bias = 0.0);
    void disc(Vec3 p, double radius_px, Color c, double depth_bias = 0.0);
    void ring(Vec3 p, double radius_px, double width_px, Color c, double depth_bias = 0.0);
    void text(Vec3 p, const std::string& s, double size_px, Color c, double depth_bias = 0.0, Vec2 offset_px = Vec2{});
    void add(SceneItem item) { m_items.push_back(std::move(item)); }
    void draw(Canvas& canvas);
    const Camera3D& camera() const { return m_camera; }

private:
    const Camera3D& m_camera;
    std::vector<SceneItem> m_items;
};

struct PlotFrame {
    double x0 = 0.0;
    double y0 = 0.0;
    double w = 1.0;
    double h = 1.0;
    double xmin = 0.0;
    double xmax = 1.0;
    double ymin = 0.0;
    double ymax = 1.0;
    bool log_x = false;
    bool log_y = false;

    Vec2 to_screen(double x, double y) const;
    double x_of(double sx) const;
    double y_of(double sy) const;
    bool inside(Vec2 s) const;
};

std::vector<double> nice_ticks(double lo, double hi, int target = 6);
void draw_plot_axes(Canvas& canvas, const PlotFrame& frame, const Theme& theme, const std::string& x_label, const std::string& y_label);

Color element_color(int atomic_number);
std::string superscript_charge(int charge);
std::string subscript_digits(const std::string& s);
std::string superscript_digits(const std::string& s);
std::string pretty_formula(const std::string& formula);
std::string pretty_text(const std::string& text);

struct UnitPrefs {
    std::map<std::string, std::string> units;

    UnitPrefs();
    std::string get(const std::string& dimension) const;
    void set(const std::string& dimension, const std::string& unit);
    std::string format(double value_si, const std::string& dimension, int significant = 5) const;
    std::string format_kj(double kj_per_mol, int significant = 5) const;
    std::string unit_of(const std::string& dimension) const { return get(dimension); }
};

struct ChemShared {
    UnitPrefs units;
    std::string last_prediction;
    std::string last_prediction_sides;
    std::string last_balanced;
};

class ChemView {
public:
    virtual ~ChemView() = default;

    virtual void draw(Canvas& canvas, const Theme& theme) = 0;
    virtual bool mouse_down(Vec2, int, unsigned) { return false; }
    virtual bool mouse_move(Vec2, bool) { return false; }
    virtual bool mouse_up(Vec2, int) { return false; }
    virtual bool wheel(Vec2, double) { return false; }
    virtual bool double_click(Vec2) { return false; }
    virtual bool key(const std::string&, unsigned) { return false; }
    virtual bool step(double) { return false; }
    virtual std::string command(const std::vector<std::string>& args) = 0;
    virtual std::string hover_text() const { return std::string(); }
};

std::unique_ptr<ChemView> make_bonding_view(ChemShared& shared);
std::unique_ptr<ChemView> make_structure_view(ChemShared& shared);
std::unique_ptr<ChemView> make_phase_view(ChemShared& shared);
std::unique_ptr<ChemView> make_orbital_view(ChemShared& shared);
std::unique_ptr<ChemView> make_cloud_view(ChemShared& shared);
std::unique_ptr<ChemView> make_element_view(ChemShared& shared);
std::unique_ptr<ChemView> make_kinetics_view(ChemShared& shared);
std::unique_ptr<ChemView> make_acidbase_view(ChemShared& shared);
std::unique_ptr<ChemView> make_spectra_view(ChemShared& shared);
std::unique_ptr<ChemView> make_network_view(ChemShared& shared);

std::string chemistry_tool(ChemShared& shared, const std::vector<std::string>& args);

struct JsonValue {
    enum Type { Null, Boolean, Number, String, Array, Object };
    Type type = Null;
    bool boolean = false;
    double number = 0.0;
    std::string string;
    std::vector<JsonValue> array;
    std::vector<std::pair<std::string, JsonValue>> object;

    const JsonValue* find(const std::string& key) const;
    double number_or(const std::string& key, double fallback) const;
    std::string string_or(const std::string& key, const std::string& fallback) const;
};

bool parse_json(const std::string& text, JsonValue& out, std::string* error = nullptr);
std::string json_quote(const std::string& s);

std::filesystem::path chemistry_save_root();
std::vector<std::string> saved_chemistry_kinds();
std::string safe_file_name(const std::string& name);
bool save_chemistry(const std::string& kind, const std::string& name, const std::string& json, std::string* error = nullptr);
std::vector<std::string> list_saved_chemistry(const std::string& kind);
bool load_chemistry(const std::string& kind, const std::string& name, std::string& json);
bool delete_saved_chemistry(const std::string& kind, const std::string& name);

} // namespace sapling

#endif // SAPLING_CHEMISTRY_HPP

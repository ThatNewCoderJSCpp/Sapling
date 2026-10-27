#include "sapling/chemistry.hpp"

#include <fizmo/../chemistry.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>

namespace sapling {

namespace fc = fizmo::chemistry;

Vec3 operator+(Vec3 a, Vec3 b) { return Vec3{ a.x + b.x, a.y + b.y, a.z + b.z }; }
Vec3 operator-(Vec3 a, Vec3 b) { return Vec3{ a.x - b.x, a.y - b.y, a.z - b.z }; }
Vec3 operator*(Vec3 a, double s) { return Vec3{ a.x * s, a.y * s, a.z * s }; }
double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 cross(Vec3 a, Vec3 b) { return Vec3{ a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }
double length(Vec3 a) { return std::sqrt(dot(a, a)); }

Vec3 normalized(Vec3 a) {
    const double l = length(a);
    return l > 1e-12 ? a * (1.0 / l) : Vec3{ 0.0, 0.0, 0.0 };
}

Vec3 Camera3D::to_view(Vec3 p) const {
    const Vec3 q = p - target;
    const double cy = std::cos(yaw), sy = std::sin(yaw);
    const double cp = std::cos(pitch), sp = std::sin(pitch);
    const double x1 = cy * q.x + sy * q.z;
    const double z1 = -sy * q.x + cy * q.z;
    const double y2 = cp * q.y - sp * z1;
    const double z2 = sp * q.y + cp * z1;
    return Vec3{ x1, y2, z2 };
}

Vec3 Camera3D::from_view(Vec3 v) const {
    const double cy = std::cos(yaw), sy = std::sin(yaw);
    const double cp = std::cos(pitch), sp = std::sin(pitch);
    const double y1 = cp * v.y + sp * v.z;
    const double z1 = -sp * v.y + cp * v.z;
    const double qx = cy * v.x - sy * z1;
    const double qz = sy * v.x + cy * z1;
    return Vec3{ qx, y1, qz } + target;
}

namespace {

double perspective_factor(double perspective, double depth) {
    if (perspective <= 0.0) return 1.0;
    return 1.0 / std::max(0.2, 1.0 - depth * perspective);
}

} // namespace

Vec2 Camera3D::project(Vec3 world, double* depth, double* scale) const {
    const Vec3 v = to_view(world);
    const double s = perspective_factor(perspective, v.z) * zoom;
    if (depth) *depth = v.z;
    if (scale) *scale = s;
    return Vec2{ 0.5 * width + pan.x + v.x * s, 0.5 * height + pan.y - v.y * s };
}

Vec3 Camera3D::unproject(Vec2 screen, double view_depth) const {
    const double s = perspective_factor(perspective, view_depth) * zoom;
    const double x = (screen.x - 0.5 * width - pan.x) / s;
    const double y = -(screen.y - 0.5 * height - pan.y) / s;
    return from_view(Vec3{ x, y, view_depth });
}

void Camera3D::orbit(double dx, double dy) {
    yaw += dx * 0.008;
    pitch += dy * 0.008;
    pitch = std::max(-1.55, std::min(1.55, pitch));
}

void Camera3D::zoom_by(double factor) {
    zoom = std::max(1e-6, std::min(1e6, zoom * factor));
}

void Camera3D::frame(const std::vector<Vec3>& points, double padding) {
    pan = Vec2{};
    if (points.empty()) return;
    Vec3 lo = points[0], hi = points[0];
    for (const Vec3& p : points) {
        lo = Vec3{ std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z) };
        hi = Vec3{ std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z) };
    }
    target = (lo + hi) * 0.5;
    double r = 0.0;
    for (const Vec3& p : points) r = std::max(r, length(p - target));
    r += padding;
    if (r <= 0.0) r = 100.0;
    zoom = 0.44 * std::min(width, height) / r;
}

void Camera3D::face_front() {
    yaw = 0.0;
    pitch = 0.0;
}

void Scene3D::sphere(Vec3 p, double radius_world, Color c, double depth_bias) {
    double d = 0.0, s = 1.0;
    const Vec2 sp = m_camera.project(p, &d, &s);
    SceneItem it;
    it.kind = SceneItem::Sphere;
    it.depth = d + depth_bias;
    it.a = sp;
    it.size = radius_world * s;
    it.color = c;
    m_items.push_back(std::move(it));
}

void Scene3D::sphere_px(Vec3 p, double radius_px, Color c, double depth_bias) {
    double d = 0.0;
    const Vec2 sp = m_camera.project(p, &d);
    SceneItem it;
    it.kind = SceneItem::Sphere;
    it.depth = d + depth_bias;
    it.a = sp;
    it.size = radius_px;
    it.color = c;
    m_items.push_back(std::move(it));
}

void Scene3D::cylinder(Vec3 a, Vec3 b, double width_world, Color ca, Color cb, double depth_bias) {
    double da = 0.0, db = 0.0, s = 1.0;
    const Vec2 pa = m_camera.project(a, &da, &s);
    const Vec2 pb = m_camera.project(b, &db);
    SceneItem it;
    it.kind = SceneItem::Cylinder;
    it.depth = 0.5 * (da + db) + depth_bias;
    it.a = pa;
    it.b = pb;
    it.width = std::max(1.0, width_world * s);
    it.color = ca;
    it.color2 = cb;
    m_items.push_back(std::move(it));
}

void Scene3D::line(Vec3 a, Vec3 b, double width_px, Color c, double depth_bias) {
    double da = 0.0, db = 0.0;
    SceneItem it;
    it.kind = SceneItem::Line;
    it.a = m_camera.project(a, &da);
    it.b = m_camera.project(b, &db);
    it.depth = 0.5 * (da + db) + depth_bias;
    it.width = width_px;
    it.color = c;
    m_items.push_back(std::move(it));
}

void Scene3D::dashed(Vec3 a, Vec3 b, double width_px, Color c, double depth_bias) {
    line(a, b, width_px, c, depth_bias);
    m_items.back().kind = SceneItem::Dashed;
}

void Scene3D::arrow(Vec3 a, Vec3 b, double width_px, Color c, double depth_bias) {
    line(a, b, width_px, c, depth_bias);
    m_items.back().kind = SceneItem::Arrow;
}

void Scene3D::disc(Vec3 p, double radius_px, Color c, double depth_bias) {
    double d = 0.0;
    SceneItem it;
    it.kind = SceneItem::Disc;
    it.a = m_camera.project(p, &d);
    it.depth = d + depth_bias;
    it.size = radius_px;
    it.color = c;
    m_items.push_back(std::move(it));
}

void Scene3D::ring(Vec3 p, double radius_px, double width_px, Color c, double depth_bias) {
    disc(p, radius_px, c, depth_bias);
    m_items.back().kind = SceneItem::Ring;
    m_items.back().width = width_px;
}

void Scene3D::text(Vec3 p, const std::string& s, double size_px, Color c, double depth_bias, Vec2 offset_px) {
    double d = 0.0;
    SceneItem it;
    it.kind = SceneItem::Text;
    it.a = m_camera.project(p, &d) + offset_px;
    it.depth = d + depth_bias;
    it.size = size_px;
    it.color = c;
    it.text = s;
    m_items.push_back(std::move(it));
}

namespace {

Color darker(Color c, double k) {
    return rgb(static_cast<int>(c.r * k), static_cast<int>(c.g * k), static_cast<int>(c.b * k), c.a);
}

Color lighter(Color c, double k) {
    return mix(c, rgb(255, 255, 255, c.a), k);
}

void draw_cylinder_half(Canvas& canvas, Vec2 a, Vec2 b, double w, Color c) {
    canvas.line_px(a, b, darker(c, 0.55), w);
    canvas.line_px(a, b, c, w * 0.72);
    const Vec2 d = b - a;
    const double l = length(d);
    if (l < 1e-6) return;
    const Vec2 n{ -d.y / l, d.x / l };
    const double side = n.y > 0.0 ? -1.0 : 1.0;
    const Vec2 off = n * (side * w * 0.16);
    canvas.line_px(a + off, b + off, lighter(c, 0.35), w * 0.25);
}

} // namespace

void Scene3D::draw(Canvas& canvas) {
    std::stable_sort(m_items.begin(), m_items.end(), [](const SceneItem& x, const SceneItem& y) { return x.depth < y.depth; });
    for (const SceneItem& it : m_items) {
        switch (it.kind) {
            case SceneItem::Sphere:
                canvas.sphere_px(it.a, it.size, it.color, 1.0);
                break;
            case SceneItem::Cylinder: {
                const Vec2 mid = (it.a + it.b) * 0.5;
                draw_cylinder_half(canvas, it.a, mid, it.width, it.color);
                draw_cylinder_half(canvas, mid, it.b, it.width, it.color2);
                break;
            }
            case SceneItem::Line:
                canvas.line_px(it.a, it.b, it.color, it.width);
                break;
            case SceneItem::Dashed:
                canvas.dashed_line_px(it.a, it.b, it.color, it.width, 5.0);
                break;
            case SceneItem::Arrow:
                canvas.arrow_px(it.a, it.b, it.color, it.width);
                break;
            case SceneItem::Disc:
                canvas.disc_px(it.a, it.size, it.color);
                break;
            case SceneItem::Ring:
                canvas.circle_px(it.a, it.size, it.color, it.width);
                break;
            case SceneItem::Text:
                canvas.text(it.a, it.text, it.size, it.color, it.align);
                break;
        }
    }
    m_items.clear();
}

Vec2 PlotFrame::to_screen(double x, double y) const {
    const auto tx = [this](double v) { return log_x ? std::log10(std::max(v, 1e-300)) : v; };
    const auto ty = [this](double v) { return log_y ? std::log10(std::max(v, 1e-300)) : v; };
    const double ax = tx(xmin), bx = tx(xmax), ay = ty(ymin), by = ty(ymax);
    const double fx = (tx(x) - ax) / (bx - ax);
    const double fy = (ty(y) - ay) / (by - ay);
    return Vec2{ x0 + fx * w, y0 + h - fy * h };
}

double PlotFrame::x_of(double sx) const {
    const double f = (sx - x0) / w;
    if (log_x) return std::pow(10.0, std::log10(xmin) + f * (std::log10(xmax) - std::log10(xmin)));
    return xmin + f * (xmax - xmin);
}

double PlotFrame::y_of(double sy) const {
    const double f = (y0 + h - sy) / h;
    if (log_y) return std::pow(10.0, std::log10(ymin) + f * (std::log10(ymax) - std::log10(ymin)));
    return ymin + f * (ymax - ymin);
}

bool PlotFrame::inside(Vec2 s) const {
    return s.x >= x0 && s.x <= x0 + w && s.y >= y0 && s.y <= y0 + h;
}

std::vector<double> nice_ticks(double lo, double hi, int target) {
    std::vector<double> out;
    if (!(hi > lo) || !std::isfinite(lo) || !std::isfinite(hi)) return out;
    const double raw = (hi - lo) / std::max(1, target);
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double r = raw / mag;
    const double step = (r < 1.5 ? 1.0 : r < 3.0 ? 2.0 : r < 7.0 ? 5.0 : 10.0) * mag;
    for (double v = std::ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) {
        out.push_back(std::fabs(v) < step * 1e-9 ? 0.0 : v);
        if (out.size() > 60) break;
    }
    return out;
}

namespace {

std::string tick_label(double v) {
    char buf[32];
    const double a = std::fabs(v);
    if (a != 0.0 && (a >= 1e5 || a < 1e-3)) {
        std::snprintf(buf, sizeof buf, "%.0e", v);
    } else {
        std::snprintf(buf, sizeof buf, "%.6g", v);
    }
    return buf;
}

std::string power_label(int e) {
    std::string out = "10";
    std::string digits = std::to_string(e);
    for (char c : digits) {
        if (c == '-') out += "\xE2\x81\xBB";
        else {
            static const char* sup[] = { "\xE2\x81\xB0", "\xC2\xB9", "\xC2\xB2", "\xC2\xB3", "\xE2\x81\xB4", "\xE2\x81\xB5", "\xE2\x81\xB6", "\xE2\x81\xB7", "\xE2\x81\xB8", "\xE2\x81\xB9" };
            out += sup[c - '0'];
        }
    }
    return out;
}

} // namespace

void draw_plot_axes(Canvas& canvas, const PlotFrame& f, const Theme& theme, const std::string& x_label, const std::string& y_label) {
    const Color grid = mix(theme.background, theme.grid, 0.8);
    canvas.fill_rect_px(f.x0, f.y0, f.w, f.h, mix(theme.background, rgb(0, 0, 0), 0.15));
    const auto ticks_for = [](double lo, double hi, bool log) {
        std::vector<double> t;
        if (log) {
            const int a = static_cast<int>(std::ceil(std::log10(lo) - 1e-9));
            const int b = static_cast<int>(std::floor(std::log10(hi) + 1e-9));
            const int stride = std::max(1, (b - a) / 8 + 1);
            for (int e = a; e <= b; e += stride) t.push_back(std::pow(10.0, e));
        } else {
            t = nice_ticks(lo, hi, 7);
        }
        return t;
    };
    for (double x : ticks_for(f.xmin, f.xmax, f.log_x)) {
        const Vec2 s = f.to_screen(x, f.log_y ? f.ymin : f.ymin);
        canvas.line_px(Vec2{ s.x, f.y0 }, Vec2{ s.x, f.y0 + f.h }, grid, 1.0);
        const std::string lbl = f.log_x ? power_label(static_cast<int>(std::lround(std::log10(x)))) : tick_label(x);
        canvas.text(Vec2{ s.x, f.y0 + f.h + 15.0 }, lbl, 10.0, theme.muted, TextAlign::Center);
    }
    for (double y : ticks_for(f.ymin, f.ymax, f.log_y)) {
        const Vec2 s = f.to_screen(f.xmin, y);
        canvas.line_px(Vec2{ f.x0, s.y }, Vec2{ f.x0 + f.w, s.y }, grid, 1.0);
        const std::string lbl = f.log_y ? power_label(static_cast<int>(std::lround(std::log10(y)))) : tick_label(y);
        canvas.text(Vec2{ f.x0 - 6.0, s.y + 4.0 }, lbl, 10.0, theme.muted, TextAlign::Right);
    }
    canvas.stroke_rect_px(f.x0, f.y0, f.w, f.h, theme.axes, 1.2);
    canvas.text(Vec2{ f.x0 + 0.5 * f.w, f.y0 + f.h + 34.0 }, x_label, 12.0, theme.text, TextAlign::Center);
    canvas.text(Vec2{ f.x0 + f.w, f.y0 - 10.0 }, y_label, 12.0, theme.text, TextAlign::Right);
}

Color element_color(int z) {
    const fc::Element* e = fc::find_element(z);
    if (e && e->cpk_color != 0) {
        Color c = rgb(e->cpk_red(), e->cpk_green(), static_cast<int>(e->cpk_color & 0xFF));
        if (c.r + c.g + c.b < 120) c = mix(c, rgb(120, 120, 130), 0.5);
        return c;
    }
    const unsigned h = static_cast<unsigned>(z) * 2654435761u;
    return rgb(100 + static_cast<int>(h % 140), 100 + static_cast<int>((h >> 8) % 140), 100 + static_cast<int>((h >> 16) % 140));
}

namespace {

const char* const SUPER[] = { "\xE2\x81\xB0", "\xC2\xB9", "\xC2\xB2", "\xC2\xB3", "\xE2\x81\xB4", "\xE2\x81\xB5", "\xE2\x81\xB6", "\xE2\x81\xB7", "\xE2\x81\xB8", "\xE2\x81\xB9" };
const char* const SUB[] = { "\xE2\x82\x80", "\xE2\x82\x81", "\xE2\x82\x82", "\xE2\x82\x83", "\xE2\x82\x84", "\xE2\x82\x85", "\xE2\x82\x86", "\xE2\x82\x87", "\xE2\x82\x88", "\xE2\x82\x89" };

} // namespace

std::string superscript_charge(int charge) {
    if (charge == 0) return std::string();
    std::string out;
    const int a = std::abs(charge);
    if (a > 1) {
        for (char c : std::to_string(a)) out += SUPER[c - '0'];
    }
    out += charge > 0 ? "\xE2\x81\xBA" : "\xE2\x81\xBB";
    return out;
}

std::string subscript_digits(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (c >= '0' && c <= '9') out += SUB[c - '0'];
        else out += c;
    }
    return out;
}

std::string superscript_digits(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (c >= '0' && c <= '9') out += SUPER[c - '0'];
        else if (c == '-') out += "\xE2\x81\xBB";
        else if (c == '+') out += "\xE2\x81\xBA";
        else if (c == '.') out += "\xC2\xB7";
        else out += c;
    }
    return out;
}

std::string pretty_formula(const std::string& formula) {
    try {
        const fc::Formula f = fc::Formula::parse(formula);
        const std::string d = f.display_string(true);
        if (!d.empty()) return d;
    } catch (...) {
    }
    std::string out;
    char prev = 0;
    for (char c : formula) {
        if (c >= '0' && c <= '9' && (std::isalpha(static_cast<unsigned char>(prev)) || prev == ')' || prev == ']' || (prev & 0x80))) {
            out += SUB[c - '0'];
        } else {
            out += c;
            prev = c;
            continue;
        }
    }
    return out;
}

namespace {

bool formula_like(const std::string& core) {
    if (core.size() < 2 || !(core[0] >= 'A' && core[0] <= 'Z')) return false;
    bool digit = false;
    bool charge = false;
    for (std::size_t i = 0; i < core.size(); ++i) {
        const char c = core[i];
        if (c >= '0' && c <= '9') digit = true;
        else if (c == '+' || c == '-') {
            if (i + 1 != core.size() && !(i + 2 == core.size() && (core[i + 1] == '+' || core[i + 1] == '-'))) {
                const std::string tail = core.substr(i + 1);
                bool digits_after = !tail.empty();
                for (char t : tail) if (!(t >= '0' && t <= '9')) digits_after = false;
                if (!digits_after) return false;
            }
            charge = true;
        } else if (c == '^' || c == '{' || c == '}' || c == '_' || c == '(' || c == ')' || c == '[' || c == ']' || c == '*') {
        } else if (!std::isalpha(static_cast<unsigned char>(c))) {
            return false;
        }
    }
    if (!digit && !charge) return false;
    try {
        const fc::Formula f = fc::Formula::parse(core);
        if (f.molar_mass() <= 0.0) return false;
        return !f.display_string(true).empty();
    } catch (...) {
        return false;
    }
}

}

std::string pretty_text(const std::string& text) {
    std::string out;
    std::size_t i = 0;
    while (i < text.size()) {
        if (text[i] == ' ' || text[i] == '\t' || text[i] == '\n') {
            out += text[i++];
            continue;
        }
        std::size_t j = i;
        while (j < text.size() && text[j] != ' ' && text[j] != '\t' && text[j] != '\n') ++j;
        std::string word = text.substr(i, j - i);
        i = j;
        if (word == "->" || word == "-->") { out += "\xE2\x86\x92"; continue; }
        if (word == "<=>" || word == "<->") { out += "\xE2\x87\x8C"; continue; }
        std::string lead;
        while (!word.empty() && (word[0] == '(' || word[0] == '"' || word[0] == '\'')) {
            if (word[0] == '(' && word.size() > 1 && word[1] >= 'A' && word[1] <= 'Z' && word.find(')') != std::string::npos && word.find(')') + 1 < word.size()) break;
            lead += word[0];
            word.erase(0, 1);
        }
        std::string trail;
        while (!word.empty() && std::string(".,;:!?\"')").find(word.back()) != std::string::npos) {
            if (word.back() == ')' && std::count(word.begin(), word.end(), '(') >= std::count(word.begin(), word.end(), ')')) break;
            trail.insert(trail.begin(), word.back());
            word.pop_back();
        }
        std::string phase;
        for (const char* p : { "(aq)", "(s)", "(l)", "(g)" }) {
            const std::string ps = p;
            if (word.size() > ps.size() && word.compare(word.size() - ps.size(), ps.size(), ps) == 0) {
                phase = ps;
                word.erase(word.size() - ps.size());
                break;
            }
        }
        std::string coefficient;
        std::size_t k = 0;
        while (k < word.size() && word[k] >= '0' && word[k] <= '9') ++k;
        if (k > 0 && k < word.size() && word[k] >= 'A' && word[k] <= 'Z') {
            coefficient = word.substr(0, k);
            word.erase(0, k);
        }
        if (formula_like(word)) out += lead + coefficient + pretty_formula(word) + phase + trail;
        else out += lead + coefficient + word + phase + trail;
    }
    return out;
}

UnitPrefs::UnitPrefs() {
    units["temperature"] = "K";
    units["temperature difference"] = "K";
    units["pressure"] = "atm";
    units["molar energy"] = "kJ/mol";
    units["energy"] = "kJ";
    units["length"] = "pm";
    units["mass"] = "g";
    units["volume"] = "L";
    units["amount"] = "mol";
    units["density"] = "g/cm\xC2\xB3";
    units["specific volume"] = "m\xC2\xB3/kg";
    units["molar volume"] = "L/mol";
    units["molar mass"] = "g/mol";
    units["dipole moment"] = "D";
    units["concentration"] = "M";
    units["molar entropy"] = "J/(mol\xC2\xB7K)";
    units["time"] = "s";
}

std::string UnitPrefs::get(const std::string& dimension) const {
    const auto it = units.find(dimension);
    if (it != units.end()) return it->second;
    return fc::si_unit(fc::dimension_named(dimension)).symbol;
}

void UnitPrefs::set(const std::string& dimension, const std::string& unit) {
    const fc::Dimension d = fc::dimension_named(dimension);
    if (fc::find_unit(unit, d)) units[dimension] = fc::find_unit(unit, d)->symbol;
}

std::string UnitPrefs::format(double value_si, const std::string& dimension, int significant) const {
    if (!std::isfinite(value_si)) return "n/a";
    return fc::format_quantity(value_si, fc::dimension_named(dimension), get(dimension), significant);
}

std::string UnitPrefs::format_kj(double kj_per_mol, int significant) const {
    return format(kj_per_mol * 1000.0, "molar energy", significant);
}

} // namespace sapling

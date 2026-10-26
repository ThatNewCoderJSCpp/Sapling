#ifndef SAPLING_CANVAS_HPP
#define SAPLING_CANVAS_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace sapling {

struct Color {
    std::uint8_t r = 0;
    std::uint8_t g = 0;
    std::uint8_t b = 0;
    std::uint8_t a = 255;
};

Color rgb(int r, int g, int b, int a = 255);
Color mix(Color a, Color b, double t);

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

Vec2 operator+(Vec2 a, Vec2 b);
Vec2 operator-(Vec2 a, Vec2 b);
Vec2 operator*(Vec2 a, double s);
double length(Vec2 a);

enum class TextAlign {
    Left,
    Center,
    Right
};

class Canvas {
public:
    void resize(int width, int height);
    int width() const { return m_w; }
    int height() const { return m_h; }

    void clear(Color c);
    void set_view(double cx, double cy, double pixels_per_unit);
    void fit(double x0, double x1, double y0, double y1, double margin = 0.08);
    void fit_stretched(double x0, double x1, double y0, double y1, double margin = 0.04);
    double scale() const { return m_scale; }
    double scale_y() const { return m_scale_y; }
    Vec2 to_screen(Vec2 world) const;
    Vec2 to_world(Vec2 screen) const;
    void visible(double& x0, double& x1, double& y0, double& y1) const;

    void line(Vec2 a, Vec2 b, Color c, double width_px = 2.0);
    void polyline(const std::vector<Vec2>& points, Color c, double width_px = 2.0);
    void fading_polyline(const std::vector<Vec2>& points, Color c, double width_px = 2.0);
    void disc(Vec2 center, double radius_px, Color c);
    void disc_world(Vec2 center, double radius, Color c);
    void circle(Vec2 center, double radius, Color c, double width_px = 1.5);
    void rect(Vec2 a, Vec2 b, Color c);
    void arrow(Vec2 from, Vec2 to, Color c, double width_px = 2.0);
    void grid(double step, Color c);
    void axes(Color c);

    void line_px(Vec2 a, Vec2 b, Color c, double width_px = 1.5);
    void dashed_line_px(Vec2 a, Vec2 b, Color c, double width_px = 1.5, double dash_px = 6.0);
    void disc_px(Vec2 center, double radius_px, Color c);
    void circle_px(Vec2 center, double radius_px, Color c, double width_px = 1.5);
    void arrow_px(Vec2 from, Vec2 to, Color c, double width_px = 2.0, double head_px = -1.0);
    void fill_rect_px(double x, double y, double w, double h, Color c);
    void stroke_rect_px(double x, double y, double w, double h, Color c, double width_px = 1.0);
    void fill_polygon(const std::vector<Vec2>& screen_points, Color c);
    void fill_triangle(Vec2 a, Vec2 b, Vec2 c, Color color);
    void sphere_px(Vec2 center, double radius_px, Color c, double light = 1.0);
    void put_pixel(int x, int y, Color c);
    void text(Vec2 screen, const std::string& s, double size_px, Color c, TextAlign align = TextAlign::Left, double weight = 1.0);
    double text_width(const std::string& s, double size_px) const;

    const std::vector<std::uint8_t>& pixels() const { return m_px; }
    bool restore(const std::vector<std::uint8_t>& rgba);
    std::vector<std::uint8_t> take_pixels();

private:
    void blend(int x, int y, Color c, double coverage);
    void capsule(Vec2 a, Vec2 b, double radius, Color c);

    int m_w = 1;
    int m_h = 1;
    double m_cx = 0.0;
    double m_cy = 0.0;
    double m_scale = 100.0;
    double m_scale_y = 100.0;
    std::vector<std::uint8_t> m_px = std::vector<std::uint8_t>(4, 0);
};

} // namespace sapling

#endif // SAPLING_CANVAS_HPP

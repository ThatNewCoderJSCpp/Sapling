#include "sapling/canvas.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace sapling {

Color rgb(int r, int g, int b, int a) {
    auto c = [](int v) { return static_cast<std::uint8_t>(std::min(255, std::max(0, v))); };
    return Color{ c(r), c(g), c(b), c(a) };
}

Color mix(Color a, Color b, double t) {
    t = std::min(1.0, std::max(0.0, t));
    auto m = [t](std::uint8_t x, std::uint8_t y) { return static_cast<int>(std::lround(x + (y - x) * t)); };
    return rgb(m(a.r, b.r), m(a.g, b.g), m(a.b, b.b), m(a.a, b.a));
}

Vec2 operator+(Vec2 a, Vec2 b) { return Vec2{ a.x + b.x, a.y + b.y }; }
Vec2 operator-(Vec2 a, Vec2 b) { return Vec2{ a.x - b.x, a.y - b.y }; }
Vec2 operator*(Vec2 a, double s) { return Vec2{ a.x * s, a.y * s }; }
double length(Vec2 a) { return std::hypot(a.x, a.y); }

void Canvas::resize(int width, int height) {
    m_w = std::max(1, width);
    m_h = std::max(1, height);
    m_px.assign(static_cast<std::size_t>(m_w) * m_h * 4, 0);
}

void Canvas::clear(Color c) {
    for (std::size_t i = 0; i + 3 < m_px.size(); i += 4) {
        m_px[i] = c.r;
        m_px[i + 1] = c.g;
        m_px[i + 2] = c.b;
        m_px[i + 3] = 255;
    }
}

void Canvas::set_view(double cx, double cy, double pixels_per_unit) {
    m_cx = cx;
    m_cy = cy;
    m_scale = pixels_per_unit > 0.0 ? pixels_per_unit : 1.0;
    m_scale_y = m_scale;
}

void Canvas::fit(double x0, double x1, double y0, double y1, double margin) {
    const double w = std::max(1e-12, x1 - x0);
    const double h = std::max(1e-12, y1 - y0);
    const double sx = m_w * (1.0 - 2.0 * margin) / w;
    const double sy = m_h * (1.0 - 2.0 * margin) / h;
    set_view(0.5 * (x0 + x1), 0.5 * (y0 + y1), std::min(sx, sy));
}

void Canvas::fit_stretched(double x0, double x1, double y0, double y1, double margin) {
    const double w = std::max(1e-12, x1 - x0);
    const double h = std::max(1e-12, y1 - y0);
    m_cx = 0.5 * (x0 + x1);
    m_cy = 0.5 * (y0 + y1);
    m_scale = m_w * (1.0 - 2.0 * margin) / w;
    m_scale_y = m_h * (1.0 - 2.0 * margin) / h;
}

Vec2 Canvas::to_screen(Vec2 p) const {
    return Vec2{ 0.5 * m_w + (p.x - m_cx) * m_scale, 0.5 * m_h - (p.y - m_cy) * m_scale_y };
}

Vec2 Canvas::to_world(Vec2 s) const {
    return Vec2{ m_cx + (s.x - 0.5 * m_w) / m_scale, m_cy - (s.y - 0.5 * m_h) / m_scale_y };
}

void Canvas::visible(double& x0, double& x1, double& y0, double& y1) const {
    const Vec2 a = to_world(Vec2{ 0.0, static_cast<double>(m_h) });
    const Vec2 b = to_world(Vec2{ static_cast<double>(m_w), 0.0 });
    x0 = a.x;
    y0 = a.y;
    x1 = b.x;
    y1 = b.y;
}

void Canvas::blend(int x, int y, Color c, double coverage) {
    if (x < 0 || y < 0 || x >= m_w || y >= m_h) return;
    const double a = std::min(1.0, std::max(0.0, coverage)) * c.a / 255.0;
    if (a <= 0.0) return;
    std::uint8_t* p = &m_px[(static_cast<std::size_t>(y) * m_w + x) * 4];
    p[0] = static_cast<std::uint8_t>(p[0] + (c.r - p[0]) * a + 0.5);
    p[1] = static_cast<std::uint8_t>(p[1] + (c.g - p[1]) * a + 0.5);
    p[2] = static_cast<std::uint8_t>(p[2] + (c.b - p[2]) * a + 0.5);
    p[3] = 255;
}

void Canvas::capsule(Vec2 a, Vec2 b, double r, Color c) {
    if (!std::isfinite(a.x) || !std::isfinite(a.y) || !std::isfinite(b.x) || !std::isfinite(b.y)) return;
    const int x0 = std::max(0, static_cast<int>(std::floor(std::min(a.x, b.x) - r - 1.0)));
    const int x1 = std::min(m_w - 1, static_cast<int>(std::ceil(std::max(a.x, b.x) + r + 1.0)));
    const int y0 = std::max(0, static_cast<int>(std::floor(std::min(a.y, b.y) - r - 1.0)));
    const int y1 = std::min(m_h - 1, static_cast<int>(std::ceil(std::max(a.y, b.y) + r + 1.0)));
    if (x0 > x1 || y0 > y1) return;
    const Vec2 d = b - a;
    const double dd = d.x * d.x + d.y * d.y;
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            const Vec2 p{ x + 0.5, y + 0.5 };
            double t = dd > 0.0 ? ((p.x - a.x) * d.x + (p.y - a.y) * d.y) / dd : 0.0;
            t = std::min(1.0, std::max(0.0, t));
            const double dist = length(p - (a + d * t));
            const double cov = r + 0.5 - dist;
            if (cov > 0.0) blend(x, y, c, cov);
        }
    }
}

void Canvas::line(Vec2 a, Vec2 b, Color c, double width_px) {
    capsule(to_screen(a), to_screen(b), 0.5 * width_px, c);
}

void Canvas::polyline(const std::vector<Vec2>& points, Color c, double width_px) {
    for (std::size_t i = 0; i + 1 < points.size(); ++i) line(points[i], points[i + 1], c, width_px);
}

void Canvas::fading_polyline(const std::vector<Vec2>& points, Color c, double width_px) {
    const std::size_t n = points.size();
    for (std::size_t i = 0; i + 1 < n; ++i) {
        Color f = c;
        f.a = static_cast<std::uint8_t>(c.a * (0.1 + 0.9 * static_cast<double>(i + 1) / static_cast<double>(n)));
        line(points[i], points[i + 1], f, width_px);
    }
}

void Canvas::disc(Vec2 center, double radius_px, Color c) {
    const Vec2 s = to_screen(center);
    capsule(s, s, radius_px, c);
}

void Canvas::disc_world(Vec2 center, double radius, Color c) {
    disc(center, radius * m_scale, c);
}

void Canvas::circle(Vec2 center, double radius, Color c, double width_px) {
    const int n = std::max(24, static_cast<int>(radius * m_scale * 0.6));
    Vec2 prev{ center.x + radius, center.y };
    for (int i = 1; i <= n; ++i) {
        const double t = 2.0 * 3.141592653589793 * i / n;
        const Vec2 p{ center.x + radius * std::cos(t), center.y + radius * std::sin(t) };
        line(prev, p, c, width_px);
        prev = p;
    }
}

void Canvas::rect(Vec2 a, Vec2 b, Color c) {
    Vec2 sa = to_screen(a);
    Vec2 sb = to_screen(b);
    const int x0 = std::max(0, static_cast<int>(std::floor(std::min(sa.x, sb.x))));
    const int x1 = std::min(m_w - 1, static_cast<int>(std::ceil(std::max(sa.x, sb.x))) - 1);
    const int y0 = std::max(0, static_cast<int>(std::floor(std::min(sa.y, sb.y))));
    const int y1 = std::min(m_h - 1, static_cast<int>(std::ceil(std::max(sa.y, sb.y))) - 1);
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) blend(x, y, c, 1.0);
    }
}

void Canvas::arrow(Vec2 from, Vec2 to, Color c, double width_px) {
    line(from, to, c, width_px);
    const Vec2 s0 = to_screen(from);
    const Vec2 s1 = to_screen(to);
    const Vec2 d = s1 - s0;
    const double l = length(d);
    if (l < 1e-9) return;
    const Vec2 u = d * (1.0 / l);
    const Vec2 nrm{ -u.y, u.x };
    const double head = std::min(10.0 + 2.0 * width_px, 0.45 * l);
    capsule(s1, s1 - u * head + nrm * (0.5 * head), 0.5 * width_px, c);
    capsule(s1, s1 - u * head - nrm * (0.5 * head), 0.5 * width_px, c);
}

void Canvas::grid(double step, Color c) {
    if (!(step > 0.0)) return;
    double x0, x1, y0, y1;
    visible(x0, x1, y0, y1);
    if ((x1 - x0) / step > 400 || (y1 - y0) / step > 400) return;
    for (double x = std::ceil(x0 / step) * step; x <= x1; x += step) line(Vec2{ x, y0 }, Vec2{ x, y1 }, c, 1.0);
    for (double y = std::ceil(y0 / step) * step; y <= y1; y += step) line(Vec2{ x0, y }, Vec2{ x1, y }, c, 1.0);
}

void Canvas::axes(Color c) {
    double x0, x1, y0, y1;
    visible(x0, x1, y0, y1);
    line(Vec2{ x0, 0.0 }, Vec2{ x1, 0.0 }, c, 1.4);
    line(Vec2{ 0.0, y0 }, Vec2{ 0.0, y1 }, c, 1.4);
}

std::vector<std::uint8_t> Canvas::take_pixels() {
    std::vector<std::uint8_t> out = std::move(m_px);
    m_px.assign(static_cast<std::size_t>(m_w) * m_h * 4, 0);
    return out;
}

} // namespace sapling

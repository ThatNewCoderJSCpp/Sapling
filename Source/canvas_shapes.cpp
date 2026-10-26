#include "sapling/canvas.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace sapling {

void Canvas::put_pixel(int x, int y, Color c) {
    blend(x, y, c, 1.0);
}

void Canvas::line_px(Vec2 a, Vec2 b, Color c, double width_px) {
    capsule(a, b, 0.5 * width_px, c);
}

void Canvas::dashed_line_px(Vec2 a, Vec2 b, Color c, double width_px, double dash_px) {
    const Vec2 d = b - a;
    const double l = length(d);
    if (l < 1e-9) return;
    const double dash = std::max(1.0, dash_px);
    const Vec2 u = d * (1.0 / l);
    for (double s = 0.0; s < l; s += 2.0 * dash) {
        const double e = std::min(l, s + dash);
        capsule(a + u * s, a + u * e, 0.5 * width_px, c);
    }
}

void Canvas::disc_px(Vec2 center, double radius_px, Color c) {
    capsule(center, center, radius_px, c);
}

void Canvas::circle_px(Vec2 center, double radius_px, Color c, double width_px) {
    const int n = std::max(16, static_cast<int>(radius_px * 0.8));
    Vec2 prev{ center.x + radius_px, center.y };
    for (int i = 1; i <= n; ++i) {
        const double t = 2.0 * 3.141592653589793 * i / n;
        const Vec2 p{ center.x + radius_px * std::cos(t), center.y + radius_px * std::sin(t) };
        capsule(prev, p, 0.5 * width_px, c);
        prev = p;
    }
}

void Canvas::arrow_px(Vec2 from, Vec2 to, Color c, double width_px, double head_px) {
    capsule(from, to, 0.5 * width_px, c);
    const Vec2 d = to - from;
    const double l = length(d);
    if (l < 1e-9) return;
    const Vec2 u = d * (1.0 / l);
    const Vec2 nrm{ -u.y, u.x };
    const double head = head_px > 0.0 ? std::min(head_px, 0.6 * l) : std::min(8.0 + 2.0 * width_px, 0.45 * l);
    fill_triangle(to, to - u * head + nrm * (0.45 * head), to - u * head - nrm * (0.45 * head), c);
}

void Canvas::fill_rect_px(double x, double y, double w, double h, Color c) {
    if (w < 0.0) {
        x += w;
        w = -w;
    }
    if (h < 0.0) {
        y += h;
        h = -h;
    }
    const int x0 = std::max(0, static_cast<int>(std::floor(x)));
    const int x1 = std::min(m_w - 1, static_cast<int>(std::ceil(x + w)) - 1);
    const int y0 = std::max(0, static_cast<int>(std::floor(y)));
    const int y1 = std::min(m_h - 1, static_cast<int>(std::ceil(y + h)) - 1);
    for (int py = y0; py <= y1; ++py) {
        const double cy = std::min<double>(py + 1, y + h) - std::max<double>(py, y);
        for (int px = x0; px <= x1; ++px) {
            const double cx = std::min<double>(px + 1, x + w) - std::max<double>(px, x);
            blend(px, py, c, cx * cy);
        }
    }
}

void Canvas::stroke_rect_px(double x, double y, double w, double h, Color c, double width_px) {
    const Vec2 a{ x, y }, b{ x + w, y }, d{ x + w, y + h }, e{ x, y + h };
    capsule(a, b, 0.5 * width_px, c);
    capsule(b, d, 0.5 * width_px, c);
    capsule(d, e, 0.5 * width_px, c);
    capsule(e, a, 0.5 * width_px, c);
}

void Canvas::fill_polygon(const std::vector<Vec2>& pts, Color c) {
    if (pts.size() < 3) return;
    double ymin = pts[0].y, ymax = pts[0].y;
    for (const Vec2& p : pts) {
        if (!std::isfinite(p.x) || !std::isfinite(p.y)) return;
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
    }
    const int y0 = std::max(0, static_cast<int>(std::floor(ymin)));
    const int y1 = std::min(m_h - 1, static_cast<int>(std::ceil(ymax)));
    const int sub = 4;
    std::vector<double> cover;
    std::vector<double> xs;
    for (int y = y0; y <= y1; ++y) {
        cover.assign(static_cast<std::size_t>(m_w) + 1, 0.0);
        bool any = false;
        int lo = m_w, hi = -1;
        for (int s = 0; s < sub; ++s) {
            const double sy = y + (s + 0.5) / sub;
            xs.clear();
            for (std::size_t i = 0; i < pts.size(); ++i) {
                const Vec2 a = pts[i];
                const Vec2 b = pts[(i + 1) % pts.size()];
                if ((a.y <= sy && b.y > sy) || (b.y <= sy && a.y > sy)) {
                    xs.push_back(a.x + (sy - a.y) * (b.x - a.x) / (b.y - a.y));
                }
            }
            std::sort(xs.begin(), xs.end());
            for (std::size_t k = 0; k + 1 < xs.size(); k += 2) {
                const double xa = std::max(0.0, xs[k]);
                const double xb = std::min(static_cast<double>(m_w), xs[k + 1]);
                if (xb <= xa) continue;
                any = true;
                const int ia = static_cast<int>(std::floor(xa));
                const int ib = std::min(m_w - 1, static_cast<int>(std::floor(xb)));
                lo = std::min(lo, ia);
                hi = std::max(hi, ib);
                if (ia == ib) {
                    cover[static_cast<std::size_t>(ia)] += (xb - xa) / sub;
                } else {
                    cover[static_cast<std::size_t>(ia)] += (ia + 1 - xa) / sub;
                    for (int x = ia + 1; x < ib; ++x) cover[static_cast<std::size_t>(x)] += 1.0 / sub;
                    cover[static_cast<std::size_t>(ib)] += (xb - ib) / sub;
                }
            }
        }
        if (!any) continue;
        for (int x = lo; x <= hi; ++x) {
            const double cv = cover[static_cast<std::size_t>(x)];
            if (cv > 0.0) blend(x, y, c, cv);
        }
    }
}

void Canvas::fill_triangle(Vec2 a, Vec2 b, Vec2 c, Color color) {
    fill_polygon(std::vector<Vec2>{ a, b, c }, color);
}

void Canvas::sphere_px(Vec2 center, double r, Color c, double light) {
    if (!(r > 0.0) || !std::isfinite(center.x) || !std::isfinite(center.y)) return;
    const int x0 = std::max(0, static_cast<int>(std::floor(center.x - r - 1.0)));
    const int x1 = std::min(m_w - 1, static_cast<int>(std::ceil(center.x + r + 1.0)));
    const int y0 = std::max(0, static_cast<int>(std::floor(center.y - r - 1.0)));
    const int y1 = std::min(m_h - 1, static_cast<int>(std::ceil(center.y + r + 1.0)));
    const double lx = -0.45, ly = -0.55, lz = 0.70;
    const double ll = std::sqrt(lx * lx + ly * ly + lz * lz);
    const double strength = std::min(1.5, std::max(0.0, light));
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            const double dx = (x + 0.5 - center.x) / r;
            const double dy = (y + 0.5 - center.y) / r;
            const double d2 = dx * dx + dy * dy;
            const double edge = (1.0 - std::sqrt(d2)) * r + 0.5;
            if (edge <= 0.0) continue;
            const double nz = std::sqrt(std::max(0.0, 1.0 - std::min(1.0, d2)));
            const double diff = std::max(0.0, (dx * lx + dy * ly + nz * lz) / ll);
            const double hx = lx / ll, hy = ly / ll, hz = lz / ll + 1.0;
            const double hl = std::sqrt(hx * hx + hy * hy + hz * hz);
            const double spec = std::pow(std::max(0.0, (dx * hx + dy * hy + nz * hz) / hl), 40.0);
            const double shade = (0.32 + 0.68 * diff) * (1.0 - 0.35 * strength) + 0.35 * strength * diff;
            Color px = rgb(static_cast<int>(c.r * shade), static_cast<int>(c.g * shade), static_cast<int>(c.b * shade), c.a);
            px = mix(px, rgb(255, 255, 255, c.a), 0.55 * spec * strength);
            blend(x, y, px, std::min(1.0, edge));
        }
    }
}

} // namespace sapling

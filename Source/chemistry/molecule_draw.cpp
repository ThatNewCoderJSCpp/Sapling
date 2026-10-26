#include "molecule_draw.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace sapling {

namespace {

constexpr double PI = 3.14159265358979323846;

double angle_diff(double a, double b) {
    double d = std::fmod(std::fabs(a - b), 2.0 * PI);
    return d > PI ? 2.0 * PI - d : d;
}

Color label_color(int z, const Theme& theme) {
    if (z == 6 || z == 1) return theme.text;
    return mix(element_color(z), theme.text, 0.25);
}

std::string charge_text(int q) {
    if (q == 0) return std::string();
    std::string out = std::abs(q) > 1 ? std::to_string(std::abs(q)) : std::string();
    out += q > 0 ? "+" : "\xE2\x88\x92";
    return out;
}

const fc::Element* elem(int z) {
    return fc::find_element(z);
}

double electronegativity(int z) {
    const fc::Element* e = elem(z);
    return e && e->electronegativity_pauling == e->electronegativity_pauling ? e->electronegativity_pauling : 0.0;
}

Vec3 view_direction(const Camera3D& cam) {
    return normalized(cam.from_view(Vec3{ 0.0, 0.0, 1.0 }) - cam.from_view(Vec3{ 0.0, 0.0, 0.0 }));
}

void dipole_arrow_2d(Canvas& canvas, Vec2 from, Vec2 to, Color c) {
    canvas.arrow_px(from, to, c, 1.8, 8.0);
    const Vec2 d = to - from;
    const double l = length(d);
    if (l < 1e-6) return;
    const Vec2 u = d * (1.0 / l);
    const Vec2 n{ -u.y, u.x };
    const Vec2 t = from + u * std::min(9.0, 0.25 * l);
    canvas.line_px(t + n * 5.0, t - n * 5.0, c, 1.8);
}

struct Placed {
    Vec2 s;
    double depth = 0.0;
    double r_px = 10.0;
};

std::vector<double> neighbor_angles(const std::vector<Placed>& placed, const std::vector<DrawBond>& bonds, int atom) {
    std::vector<double> out;
    for (const DrawBond& b : bonds) {
        if (b.a != atom && b.b != atom) continue;
        const int o = b.a == atom ? b.b : b.a;
        if (o < 0 || o >= static_cast<int>(placed.size())) continue;
        const Vec2 d = placed[static_cast<std::size_t>(o)].s - placed[static_cast<std::size_t>(atom)].s;
        if (length(d) < 1e-6) continue;
        out.push_back(std::atan2(d.y, d.x));
    }
    return out;
}

} // namespace

double draw_radius_pm(int z, bool flat) {
    if (flat) return 24.0;
    const fc::Element* e = elem(z);
    const double cov = e && e->covalent_radius == e->covalent_radius ? e->covalent_radius : 75.0;
    return std::max(15.0, 0.42 * cov);
}

std::vector<double> lone_pair_angles(const std::vector<double>& bond_angles, int count) {
    std::vector<double> used = bond_angles;
    std::vector<double> out;
    for (int k = 0; k < count; ++k) {
        if (used.empty()) {
            out.push_back(-0.5 * PI);
            used.push_back(-0.5 * PI);
            continue;
        }
        double best = 0.0, best_score = -1.0;
        for (int i = 0; i < 72; ++i) {
            const double a = -0.5 * PI + i * (2.0 * PI / 72.0);
            double m = 10.0;
            for (double u : used) m = std::min(m, angle_diff(a, u));
            if (m > best_score + 1e-9) {
                best_score = m;
                best = a;
            }
        }
        out.push_back(best);
        used.push_back(best);
    }
    return out;
}

double bond_dipole_debye(int za, int zb, double length_pm) {
    const double dx = std::fabs(electronegativity(za) - electronegativity(zb));
    const double ionic = fc::percent_ionic_character(dx) / 100.0;
    return fc::bond_dipole(ionic, length_pm);
}

void draw_hud(Canvas& canvas, const Theme& theme, const std::vector<std::string>& lines, bool bottom) {
    if (lines.empty()) return;
    double w = 0.0;
    for (const std::string& l : lines) w = std::max(w, canvas.text_width(l, 12.0));
    const double h = 18.0 * static_cast<double>(lines.size()) + 8.0;
    const double x = 10.0;
    const double y = bottom ? canvas.height() - h - 10.0 : 10.0;
    canvas.fill_rect_px(x, y, w + 16.0, h, rgb(theme.background.r, theme.background.g, theme.background.b, 200));
    canvas.stroke_rect_px(x, y, w + 16.0, h, theme.grid, 1.0);
    for (std::size_t i = 0; i < lines.size(); ++i) {
        canvas.text(Vec2{ x + 8.0, y + 18.0 * static_cast<double>(i + 1) }, lines[i], 12.0, i == 0 ? theme.text : theme.muted);
    }
}

void draw_molecule(Canvas& canvas, const Camera3D& cam, const Theme& theme, const std::vector<DrawAtom>& atoms, const std::vector<DrawBond>& bonds, const DrawOptions& opt) {
    std::vector<Placed> placed(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        double s = 1.0;
        placed[i].s = cam.project(atoms[i].p, &placed[i].depth, &s);
        placed[i].r_px = draw_radius_pm(atoms[i].z, opt.flat) * opt.atom_scale * s;
    }
    const auto valid_bond = [&](const DrawBond& b) {
        return b.a >= 0 && b.b >= 0 && b.a < static_cast<int>(atoms.size()) && b.b < static_cast<int>(atoms.size()) && b.a != b.b;
    };

    if (opt.flat) {
        for (const DrawBond& b : bonds) {
            if (!valid_bond(b)) continue;
            const Placed& pa = placed[static_cast<std::size_t>(b.a)];
            const Placed& pb = placed[static_cast<std::size_t>(b.b)];
            Vec2 d = pb.s - pa.s;
            const double l = length(d);
            if (l < 1e-6) continue;
            const Vec2 u = d * (1.0 / l);
            const Vec2 n{ -u.y, u.x };
            const double ra = std::min(0.45 * l, pa.r_px * 0.8);
            const double rb = std::min(0.45 * l, pb.r_px * 0.8);
            const Vec2 a = pa.s + u * ra;
            const Vec2 e = pb.s - u * rb;
            Color c = b.selected ? theme.accent : (b.preview ? theme.muted : theme.text);
            const double gap = std::max(3.0, pa.r_px * 0.22);
            const int whole = static_cast<int>(std::floor(b.order + 1e-6));
            const bool partial = b.order - whole > 0.05;
            const int lines = whole + (partial ? 1 : 0);
            for (int k = 0; k < lines; ++k) {
                const double off = (k - 0.5 * (lines - 1)) * gap;
                if (partial && k == lines - 1) canvas.dashed_line_px(a + n * off, e + n * off, c, 1.6, 4.0);
                else canvas.line_px(a + n * off, e + n * off, c, 1.8);
            }
            if (b.selected) canvas.line_px(a, e, rgb(theme.accent.r, theme.accent.g, theme.accent.b, 70), gap * lines + 8.0);
        }
        if (opt.bond_dipoles) {
            for (const DrawBond& b : bonds) {
                if (!valid_bond(b)) continue;
                const int za = atoms[static_cast<std::size_t>(b.a)].z, zb = atoms[static_cast<std::size_t>(b.b)].z;
                const double xa = electronegativity(za), xb = electronegativity(zb);
                if (std::fabs(xa - xb) < 0.05) continue;
                Vec2 from = placed[static_cast<std::size_t>(xa < xb ? b.a : b.b)].s;
                Vec2 to = placed[static_cast<std::size_t>(xa < xb ? b.b : b.a)].s;
                const Vec2 d = to - from;
                const double l = length(d);
                if (l < 1e-6) continue;
                const Vec2 u = d * (1.0 / l);
                const Vec2 n{ -u.y, u.x };
                const double mu = bond_dipole_debye(za, zb, length(atoms[static_cast<std::size_t>(b.a)].p - atoms[static_cast<std::size_t>(b.b)].p));
                const double frac = std::min(0.8, 0.25 + 0.35 * mu);
                const Vec2 mid = (from + to) * 0.5 + n * 12.0;
                dipole_arrow_2d(canvas, mid - u * (0.5 * frac * l), mid + u * (0.5 * frac * l), theme.warm);
            }
        }
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            const DrawAtom& at = atoms[i];
            const Placed& p = placed[i];
            const fc::Element* e = elem(at.z);
            canvas.disc_px(p.s, p.r_px * 0.78, theme.background);
            if (at.highlight) canvas.disc_px(p.s, p.r_px * 0.95, rgb(theme.yellow.r, theme.yellow.g, theme.yellow.b, 60));
            if (at.selected) canvas.circle_px(p.s, p.r_px * 0.95, theme.accent, 2.0);
            if (at.problem) canvas.circle_px(p.s, p.r_px * 1.05, rgb(theme.red.r, theme.red.g, theme.red.b, 170), 1.4);
            const double size = std::max(10.0, std::min(34.0, p.r_px * 1.05));
            canvas.text(Vec2{ p.s.x, p.s.y + 0.5 * size }, e ? e->symbol : "?", size, label_color(at.z, theme), TextAlign::Center, 1.15);
            if (opt.lone_pairs && (at.lone_pairs > 0 || at.unpaired > 0)) {
                const std::vector<double> used = neighbor_angles(placed, bonds, static_cast<int>(i));
                const std::vector<double> lp = lone_pair_angles(used, at.lone_pairs + at.unpaired);
                const double rr = p.r_px * 0.95;
                const double dot = std::max(1.6, p.r_px * 0.09);
                for (std::size_t k = 0; k < lp.size(); ++k) {
                    const Vec2 u{ std::cos(lp[k]), std::sin(lp[k]) };
                    const Vec2 n{ -u.y, u.x };
                    const Vec2 c = p.s + u * rr;
                    if (static_cast<int>(k) < at.lone_pairs) {
                        canvas.disc_px(c + n * (dot * 1.9), dot, theme.text);
                        canvas.disc_px(c - n * (dot * 1.9), dot, theme.text);
                    } else {
                        canvas.disc_px(c, dot, theme.yellow);
                    }
                }
            }
            if (opt.formal_charges && at.formal_charge != 0) {
                const Vec2 c = p.s + Vec2{ p.r_px * 0.85, -p.r_px * 0.85 };
                const double fs = std::max(9.0, size * 0.5);
                canvas.circle_px(c, fs * 0.75, at.formal_charge > 0 ? theme.accent : theme.red, 1.2);
                canvas.text(Vec2{ c.x, c.y + fs * 0.42 }, charge_text(at.formal_charge), fs * 0.85, at.formal_charge > 0 ? theme.accent : theme.red, TextAlign::Center);
            }
            if (opt.partial_charges && std::fabs(at.partial) > 0.04) {
                const Vec2 c = p.s + Vec2{ p.r_px * 0.9, p.r_px * 1.2 };
                canvas.text(c, at.partial > 0 ? "\xCE\xB4+" : "\xCE\xB4\xE2\x88\x92", std::max(9.0, size * 0.5), at.partial > 0 ? theme.accent : theme.red, TextAlign::Left);
            }
        }
    } else {
        Scene3D scene(cam);
        const Vec3 vd = view_direction(cam);
        for (const DrawBond& b : bonds) {
            if (!valid_bond(b)) continue;
            const DrawAtom& A = atoms[static_cast<std::size_t>(b.a)];
            const DrawAtom& B = atoms[static_cast<std::size_t>(b.b)];
            const Vec3 dir = normalized(B.p - A.p);
            Vec3 perp = normalized(cross(dir, vd));
            if (length(perp) < 0.5) perp = normalized(cross(dir, Vec3{ 0.0, 1.0, 0.0 }));
            const int whole = std::max(1, static_cast<int>(std::floor(b.order + 1e-6)));
            const bool partial = b.order - whole > 0.05;
            const int lines = whole + (partial ? 1 : 0);
            const double w = (lines == 1 ? 13.0 : lines == 2 ? 8.0 : 6.0) * opt.atom_scale;
            const double gap = (lines == 2 ? 11.0 : 10.0) * opt.atom_scale;
            Color ca = b.preview ? theme.muted : mix(element_color(A.z), rgb(200, 200, 200), 0.25);
            Color cb = b.preview ? theme.muted : mix(element_color(B.z), rgb(200, 200, 200), 0.25);
            if (b.selected) {
                ca = theme.accent;
                cb = theme.accent;
            }
            for (int k = 0; k < lines; ++k) {
                const Vec3 off = perp * ((k - 0.5 * (lines - 1)) * gap);
                if (partial && k == lines - 1) scene.dashed(A.p + off, B.p + off, 2.0, theme.text);
                else scene.cylinder(A.p + off, B.p + off, w, ca, cb);
            }
        }
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            const DrawAtom& at = atoms[i];
            const double r = draw_radius_pm(at.z, false) * opt.atom_scale;
            scene.sphere(at.p, r, element_color(at.z));
            const Placed& p = placed[i];
            if (at.selected) scene.ring(at.p, p.r_px + 3.0, 2.2, theme.accent, r + 1.0);
            if (at.highlight) scene.ring(at.p, p.r_px + 3.0, 2.0, theme.yellow, r + 1.0);
            if (at.problem) scene.ring(at.p, p.r_px + 5.0, 1.4, rgb(theme.red.r, theme.red.g, theme.red.b, 170), r + 1.0);
            const fc::Element* e = elem(at.z);
            if (opt.labels && e) {
                const double size = std::max(9.0, std::min(18.0, p.r_px * 0.9));
                const Color tc = (element_color(at.z).r + element_color(at.z).g + element_color(at.z).b) > 450 ? rgb(30, 30, 40) : rgb(240, 240, 245);
                scene.text(at.p, e->symbol, size, tc, r + 0.5, Vec2{ 0.0, size * 0.42 });
            }
            if (opt.lone_pairs && (at.lone_pairs > 0 || at.unpaired > 0)) {
                const std::vector<double> used = neighbor_angles(placed, bonds, static_cast<int>(i));
                const std::vector<double> lp = lone_pair_angles(used, at.lone_pairs + at.unpaired);
                const double rr = p.r_px + 6.0;
                for (std::size_t k = 0; k < lp.size(); ++k) {
                    const Vec2 u{ std::cos(lp[k]), std::sin(lp[k]) };
                    const Vec2 n{ -u.y, u.x };
                    SceneItem it;
                    it.kind = SceneItem::Disc;
                    it.depth = p.depth + r + 0.4;
                    it.size = 2.0;
                    it.color = static_cast<int>(k) < at.lone_pairs ? theme.text : theme.yellow;
                    if (static_cast<int>(k) < at.lone_pairs) {
                        it.a = p.s + u * rr + n * 3.6;
                        scene.add(it);
                        it.a = p.s + u * rr - n * 3.6;
                        scene.add(it);
                    } else {
                        it.a = p.s + u * rr;
                        scene.add(it);
                    }
                }
            }
            if (opt.formal_charges && at.formal_charge != 0) {
                SceneItem it;
                it.kind = SceneItem::Text;
                it.depth = p.depth + r + 0.6;
                it.a = p.s + Vec2{ p.r_px * 0.8 + 4.0, -p.r_px * 0.8 };
                it.text = charge_text(at.formal_charge);
                it.size = 12.0;
                it.color = at.formal_charge > 0 ? theme.accent : theme.red;
                scene.add(it);
            }
            if (opt.partial_charges && std::fabs(at.partial) > 0.04) {
                SceneItem it;
                it.kind = SceneItem::Text;
                it.depth = p.depth + r + 0.6;
                it.a = p.s + Vec2{ p.r_px * 0.8 + 4.0, p.r_px + 10.0 };
                it.text = at.partial > 0 ? "\xCE\xB4+" : "\xCE\xB4\xE2\x88\x92";
                it.size = 11.0;
                it.align = TextAlign::Left;
                it.color = at.partial > 0 ? theme.accent : theme.red;
                scene.add(it);
            }
        }
        scene.draw(canvas);
        if (opt.bond_dipoles) {
            for (const DrawBond& b : bonds) {
                if (!valid_bond(b)) continue;
                const int za = atoms[static_cast<std::size_t>(b.a)].z, zb = atoms[static_cast<std::size_t>(b.b)].z;
                const double xa = electronegativity(za), xb = electronegativity(zb);
                if (std::fabs(xa - xb) < 0.05) continue;
                const Vec2 from = placed[static_cast<std::size_t>(xa < xb ? b.a : b.b)].s;
                const Vec2 to = placed[static_cast<std::size_t>(xa < xb ? b.b : b.a)].s;
                const Vec2 d = to - from;
                const double l = length(d);
                if (l < 1e-6) continue;
                const Vec2 u = d * (1.0 / l);
                const Vec2 n{ -u.y, u.x };
                const double mu = bond_dipole_debye(za, zb, length(atoms[static_cast<std::size_t>(b.a)].p - atoms[static_cast<std::size_t>(b.b)].p));
                const double frac = std::min(0.8, 0.25 + 0.35 * mu);
                const Vec2 mid = (from + to) * 0.5 + n * 12.0;
                dipole_arrow_2d(canvas, mid - u * (0.5 * frac * l), mid + u * (0.5 * frac * l), theme.warm);
            }
        }
    }

    if (opt.molecular_dipole && !atoms.empty()) {
        Vec3 c;
        for (const DrawAtom& a : atoms) c = c + a.p;
        c = c * (1.0 / static_cast<double>(atoms.size()));
        const Vec2 sc = cam.project(c);
        if (opt.dipole_debye < 0.05) {
            double low = sc.y;
            for (const DrawAtom& a : atoms) low = std::max(low, cam.project(a.p).y);
            canvas.text(Vec2{ sc.x, low + 40.0 }, "\xCE\xBC \xE2\x89\x88 0 (nonpolar)", 12.0, theme.green, TextAlign::Center);
        } else {
            const Vec3 dir = normalized(opt.dipole);
            const Vec2 tip = cam.project(c + dir * 100.0);
            Vec2 d = tip - sc;
            const double l = length(d);
            const double len = std::min(170.0, std::max(40.0, 45.0 * opt.dipole_debye));
            if (l < 4.0) {
                canvas.circle_px(sc, 12.0, theme.yellow, 2.0);
                canvas.disc_px(sc, 3.0, theme.yellow);
            } else {
                const Vec2 u = d * (1.0 / l);
                const double proj = std::min(1.0, l / (100.0 * cam.zoom));
                const Vec2 from = sc - u * (0.5 * len * proj);
                const Vec2 to = sc + u * (0.5 * len * proj);
                canvas.arrow_px(from, to, theme.yellow, 3.0, 13.0);
                const Vec2 n{ -u.y, u.x };
                const Vec2 t = from + u * 10.0;
                canvas.line_px(t + n * 8.0, t - n * 8.0, theme.yellow, 3.0);
            }
            char buf[48];
            std::snprintf(buf, sizeof buf, "\xCE\xBC = %.2f D", opt.dipole_debye);
            canvas.text(Vec2{ sc.x + 14.0, sc.y - 14.0 }, buf, 12.0, theme.yellow, TextAlign::Left);
        }
    }
}

int pick_atom(const Camera3D& cam, const std::vector<DrawAtom>& atoms, Vec2 screen, bool flat, double atom_scale) {
    int best = -1;
    double best_depth = -1e300;
    double best_d = 1e300;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        double depth = 0.0, s = 1.0;
        const Vec2 p = cam.project(atoms[i].p, &depth, &s);
        const double r = std::max(8.0, draw_radius_pm(atoms[i].z, flat) * atom_scale * s) + 3.0;
        const double d = length(p - screen);
        if (d > r) continue;
        if (flat ? d < best_d : depth > best_depth) {
            best = static_cast<int>(i);
            best_depth = depth;
            best_d = d;
        }
    }
    return best;
}

int pick_bond(const Camera3D& cam, const std::vector<DrawAtom>& atoms, const std::vector<DrawBond>& bonds, Vec2 screen) {
    int best = -1;
    double best_d = 7.0;
    for (std::size_t i = 0; i < bonds.size(); ++i) {
        const DrawBond& b = bonds[i];
        if (b.a < 0 || b.b < 0 || b.a >= static_cast<int>(atoms.size()) || b.b >= static_cast<int>(atoms.size())) continue;
        const Vec2 a = cam.project(atoms[static_cast<std::size_t>(b.a)].p);
        const Vec2 c = cam.project(atoms[static_cast<std::size_t>(b.b)].p);
        const Vec2 d = c - a;
        const double dd = d.x * d.x + d.y * d.y;
        if (dd < 1e-9) continue;
        double t = ((screen.x - a.x) * d.x + (screen.y - a.y) * d.y) / dd;
        if (t < 0.12 || t > 0.88) continue;
        const Vec2 q = a + d * t;
        const double dist = length(screen - q);
        if (dist < best_d) {
            best_d = dist;
            best = static_cast<int>(i);
        }
    }
    return best;
}

namespace {

void jacobi4(double a[4][4], double v[4][4]) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) v[i][j] = i == j ? 1.0 : 0.0;
    }
    for (int sweep = 0; sweep < 60; ++sweep) {
        double off = 0.0;
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) off += a[i][j] * a[i][j];
        }
        if (off < 1e-22) break;
        for (int p = 0; p < 4; ++p) {
            for (int q = p + 1; q < 4; ++q) {
                if (std::fabs(a[p][q]) < 1e-30) continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;
                for (int k = 0; k < 4; ++k) {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 4; ++k) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 4; ++k) {
                    const double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
}

} // namespace

Vec3 aligned_vector(const std::vector<Vec3>& from, const std::vector<Vec3>& to, Vec3 v) {
    const std::size_t n = std::min(from.size(), to.size());
    if (n < 2) return v;
    Vec3 cf, ct;
    for (std::size_t i = 0; i < n; ++i) {
        cf = cf + from[i];
        ct = ct + to[i];
    }
    cf = cf * (1.0 / static_cast<double>(n));
    ct = ct * (1.0 / static_cast<double>(n));
    double S[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    double spread = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 a = from[i] - cf;
        const Vec3 b = to[i] - ct;
        spread += dot(b, b);
        const double av[3] = { a.x, a.y, a.z };
        const double bv[3] = { b.x, b.y, b.z };
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) S[r][c] += av[r] * bv[c];
        }
    }
    if (spread < 1e-9) return v;
    const double Sxx = S[0][0], Sxy = S[0][1], Sxz = S[0][2];
    const double Syx = S[1][0], Syy = S[1][1], Syz = S[1][2];
    const double Szx = S[2][0], Szy = S[2][1], Szz = S[2][2];
    double N[4][4] = {
        { Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx },
        { Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz },
        { Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy },
        { Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz }
    };
    double V[4][4];
    jacobi4(N, V);
    int best = 0;
    for (int i = 1; i < 4; ++i) {
        if (N[i][i] > N[best][best]) best = i;
    }
    const double qw = V[0][best], qx = V[1][best], qy = V[2][best], qz = V[3][best];
    const Vec3 u{ qx, qy, qz };
    const Vec3 t = cross(u, v) * 2.0;
    return v + t * qw + cross(u, t);
}

} // namespace sapling

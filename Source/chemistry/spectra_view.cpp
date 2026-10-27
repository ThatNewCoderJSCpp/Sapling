#include "molecule_draw.hpp"
#include "view_util.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace sapling {

using namespace chem_util;

namespace {

Color to_color(const fc::Srgb8& c, std::uint8_t a = 255) {
    return rgb(c.r, c.g, c.b, a);
}

std::string symbol_of(int z) {
    return z >= 1 && z <= 118 ? std::string(fc::element(z).symbol) : std::string("?");
}

std::string roman(int ion) {
    static const char* const r[] = { "", "I", "II", "III", "IV" };
    return ion >= 1 && ion <= 4 ? r[ion] : std::to_string(ion);
}

std::string band_name(double nm) {
    if (nm < 380.0) return "ultraviolet";
    if (nm < 450.0) return "violet";
    if (nm < 485.0) return "blue";
    if (nm < 500.0) return "cyan";
    if (nm < 565.0) return "green";
    if (nm < 590.0) return "yellow";
    if (nm < 625.0) return "orange";
    if (nm <= 750.0) return "red";
    return "infrared";
}

std::vector<std::pair<double, double>> parse_pairs(const std::string& text, bool& had_values) {
    std::vector<std::pair<double, double>> out;
    had_values = false;
    std::string cleaned;
    for (char c : text) cleaned += (c == ';' || c == '\r') ? '\n' : c;
    for (const std::string& raw : split(cleaned, '\n')) {
        std::string line = trim(raw);
        if (line.empty() || line[0] == '#') continue;
        for (char& c : line) {
            if (c == ',' || c == '\t' || c == ':') c = ' ';
        }
        std::istringstream in(line);
        double nm = 0.0, v = 1.0;
        if (!(in >> nm)) continue;
        if (in >> v) had_values = true;
        else v = 1.0;
        out.emplace_back(nm, v);
    }
    return out;
}

class SpectraView : public ChemView {
public:
    explicit SpectraView(ChemShared& shared) : m_shared(shared) { rebuild(); }

    void draw(Canvas& canvas, const Theme& theme) override {
        canvas.clear(theme.background);
        const double w = canvas.width(), h = canvas.height();
        const double left = 56.0, right = 120.0;
        const double strip_y = 80.0, strip_h = std::max(40.0, h * 0.16);
        m_frame.x0 = left;
        m_frame.w = std::max(50.0, w - left - right);
        m_frame.y0 = strip_y + strip_h + 34.0;
        m_frame.h = std::max(50.0, h - m_frame.y0 - 48.0);
        m_frame.xmin = m_start;
        m_frame.xmax = m_end;
        double top = 0.0;
        for (double v : m_shown.intensity) top = std::max(top, v);
        for (double v : m_observed.intensity) top = std::max(top, v);
        for (double v : m_fit.fitted.intensity) top = std::max(top, v);
        m_frame.ymin = 0.0;
        m_frame.ymax = top > 0.0 ? top * 1.1 : 1.0;
        draw_strip(canvas, theme, left, strip_y, m_frame.w, strip_h);
        draw_plot_axes(canvas, m_frame, theme, "wavelength (nm)", m_kind == "absorption" ? "transmitted light" : "relative intensity");
        const auto curve = [&](const fc::Spectrum& s, Color c, double width, bool fill) {
            if (s.size() < 2) return;
            Vec2 prev{};
            for (std::size_t i = 0; i < s.size(); ++i) {
                const double nm = s.wavelength(i);
                if (nm < m_start || nm > m_end) continue;
                const Vec2 q = m_frame.to_screen(nm, s.intensity[i]);
                if (fill) {
                    const Vec2 base = m_frame.to_screen(nm, 0.0);
                    canvas.line_px(base, q, to_color(fc::wavelength_color(nm), 150), std::max(1.0, m_frame.w / std::max<double>(1, s.size()) + 0.5));
                }
                if (i > 0) canvas.line_px(prev, q, c, width);
                prev = q;
            }
        };
        if (m_mode == "identify" && m_observed.size() > 1) {
            curve(m_observed, theme.text, 1.6, true);
            if (m_fit.ok) curve(m_fit.fitted, theme.yellow, 1.4, false);
            canvas.text(Vec2{ m_frame.x0 + m_frame.w - 6.0, m_frame.y0 + 52.0 }, "white: measured   yellow: best fit", 10.0, theme.muted, TextAlign::Right);
        } else {
            curve(m_shown, theme.text, 1.4, m_kind != "absorption");
        }
        draw_line_labels(canvas, theme);
        draw_swatches(canvas, theme, w - right + 14.0, strip_y, right - 24.0);
        std::string head = m_mode == "identify" ? "Identify elements from a spectrum" : m_mode == "color" ? "Match a color" : (m_kind == "absorption" ? "Absorption spectrum" : "Emission spectrum");
        draw_hud(canvas, theme, { head, m_title }, false);
        if (m_frame.inside(m_mouse) || (m_mouse.y >= strip_y && m_mouse.y <= strip_y + strip_h && m_mouse.x >= left && m_mouse.x <= left + m_frame.w)) {
            const double nm = m_frame.x_of(m_mouse.x);
            std::string line = format_number(nm, 5) + " nm (" + band_name(nm) + ")";
            const fc::EmissionLine* best = nullptr;
            double bd = 1e300;
            for (const fc::EmissionLine& l : fc::emission_lines()) {
                if (!shown_element(l.atomic_number)) continue;
                const double d = std::fabs(l.wavelength_nm - nm);
                if (d < bd) {
                    bd = d;
                    best = &l;
                }
            }
            if (best && bd < 3.0) line += "  \xC2\xB7  " + symbol_of(best->atomic_number) + " " + roman(best->ion) + " " + format_number(best->wavelength_nm, 6) + " nm";
            canvas.line_px(Vec2{ m_mouse.x, strip_y }, Vec2{ m_mouse.x, m_frame.y0 + m_frame.h }, rgb(255, 255, 255, 70), 1.0);
            draw_hud(canvas, theme, { line }, true);
        }
    }

    bool mouse_move(Vec2 p, bool) override {
        m_mouse = p;
        return true;
    }

    std::string command(const std::vector<std::string>& a) override {
        const std::string c = arg(a, 0);
        if (c == "readouts") return m_readouts;
        if (c == "elements_list") {
            std::string out;
            for (int z : fc::elements_with_emission_lines(300.0, 1000.0)) out += symbol_of(z) + "\t" + fc::element(z).name + "\n";
            return out;
        }
        if (c == "mix") {
            m_mode = "mix";
            m_mix_text = arg(a, 1);
        } else if (c == "kind") {
            m_kind = arg(a, 1) == "absorption" ? "absorption" : "emission";
        } else if (c == "range") {
            const double lo = num(arg(a, 1), 380.0), hi = num(arg(a, 2), 780.0);
            if (hi > lo + 5.0) {
                m_start = std::max(100.0, lo);
                m_end = std::min(2000.0, hi);
            }
        } else if (c == "fwhm") {
            m_fwhm = std::max(0.05, std::min(50.0, num(arg(a, 1), 1.0)));
        } else if (c == "candidates") {
            m_candidates_text = arg(a, 1);
        } else if (c == "observed") {
            m_mode = "identify";
            m_observed_text = arg(a, 1);
            m_profile.clear();
        } else if (c == "profile") {
            m_mode = "identify";
            m_profile_start = num(arg(a, 1), 380.0);
            m_profile_end = num(arg(a, 2), 780.0);
            m_profile.clear();
            for (const std::string& v : split(arg(a, 3), ',')) {
                const double x = num(trim(v), fc::NO_VALUE);
                if (finite(x)) m_profile.push_back(x);
            }
            m_observed_text.clear();
        } else if (c == "identify") {
            m_mode = "identify";
        } else if (c == "color") {
            m_mode = "color";
            m_color_text = arg(a, 1);
            m_color_flame = arg(a, 2) != "emission";
        } else if (c == "mode") {
            m_mode = arg(a, 1);
        } else {
            return "error\tunknown command " + c;
        }
        rebuild();
        return m_error.empty() ? "ok\t" + m_summary : "error\t" + m_error;
    }

private:
    bool shown_element(int z) const {
        if (m_mode == "identify") {
            for (const auto& p : m_fit.percentages) {
                if (p.first == z) return true;
            }
            return false;
        }
        for (const auto& f : m_fractions) {
            if (f.first == z) return true;
        }
        return false;
    }

    std::vector<int> candidate_list() const {
        std::vector<int> out;
        for (const std::string& w : split(m_candidates_text, ' ')) {
            for (const std::string& part : split(w, ',')) {
                const std::string t = trim(part);
                if (t.empty()) continue;
                if (const fc::Element* e = element_of(t)) out.push_back(e->atomic_number);
            }
        }
        return out;
    }

    void parse_mix() {
        m_fractions.clear();
        std::string text = m_mix_text;
        for (char& c : text) {
            if (c == ',' || c == '\n') c = ';';
        }
        for (const std::string& part : split(text, ';')) {
            std::string t = trim(part);
            if (t.empty()) continue;
            for (char& c : t) {
                if (c == '=' || c == ':' || c == '%') c = ' ';
            }
            std::istringstream in(t);
            std::string sym;
            double pct = 100.0;
            in >> sym;
            if (!(in >> pct)) pct = 100.0;
            const fc::Element* e = element_of(sym);
            if (!e) throw std::invalid_argument("'" + sym + "' isn't an element");
            if (pct > 0.0) m_fractions.emplace_back(e->atomic_number, pct);
        }
        if (m_fractions.empty()) throw std::invalid_argument("give at least one element, like Na 30, Li 70");
    }

    void rebuild() {
        m_error.clear();
        m_summary.clear();
        m_fit = fc::SpectrumFit{};
        m_color_fit = fc::ColorFit{};
        m_shown = fc::Spectrum{};
        m_observed = fc::Spectrum{};
        std::ostringstream o;
        const double step = std::max(0.05, std::min(1.0, (m_end - m_start) / 1600.0));
        try {
            if (m_mode == "identify") {
                build_identify(o, step);
            } else if (m_mode == "color") {
                build_color(o);
            } else {
                parse_mix();
                build_mix(o, step);
            }
        } catch (const std::exception& e) {
            m_error = e.what();
            o << "#Spectra\n!" << e.what() << "\n";
        }
        m_readouts = o.str();
    }

    void element_lines(std::ostringstream& o, int z, double weight) {
        const fc::Element& e = fc::element(z);
        o << "#" << e.name << " (" << e.symbol << ")" << (weight > 0 ? "  " + format_number(weight, 4) + "%" : std::string()) << "\n";
        std::vector<fc::EmissionLine> lines = fc::emission_lines_of(z, m_start, m_end);
        std::sort(lines.begin(), lines.end(), [](const fc::EmissionLine& a, const fc::EmissionLine& b) { return a.intensity > b.intensity; });
        if (lines.size() > 6) lines.resize(6);
        if (lines.empty()) o << "no strong lines between " << format_number(m_start, 4) << " and " << format_number(m_end, 4) << " nm\n";
        for (const fc::EmissionLine& l : lines) o << e.symbol << " " << roman(l.ion) << " " << format_number(l.wavelength_nm, 6) << " nm\t" << band_name(l.wavelength_nm) << ", strength " << format_number(l.intensity, 6) << "\n";
        const fc::Srgb8 own = fc::spectrum_color(fc::emission_spectrum({ { z, 1.0 } }, 1.0, 380.0, 780.0, 1.0));
        o << "Glow color\t" << fc::to_hex(own) << "\n";
        if (fc::has_flame_color(z)) o << "Flame test\t" << fc::flame_color_name(z) << " (" << fc::to_hex(fc::flame_color(z)) << ")\n";
    }

    void build_mix(std::ostringstream& o, double step) {
        const fc::Spectrum emit = fc::emission_spectrum(m_fractions, m_fwhm, m_start, m_end, step);
        m_shown = m_kind == "absorption" ? fc::absorption_spectrum(m_fractions, m_fwhm, 0.9, m_start, m_end, step) : emit;
        const fc::Spectrum visible = fc::emission_spectrum(m_fractions, m_fwhm, 380.0, 780.0, 0.5);
        m_glow = fc::spectrum_color(visible);
        double total = 0.0;
        for (const auto& f : m_fractions) total += f.second;
        o << "#Mixture\n";
        for (const auto& f : m_fractions) o << fc::element(f.first).name << "\t" << format_number(100.0 * f.second / total, 4) << "% of the atoms\n";
        o << "Color of the light\t" << fc::to_hex(m_glow) << "\n";
        o << "Each element's lines are scaled by its share; the intensities assume every element glows equally strongly per atom, which real flames and lamps don't quite do.\n";
        m_title = "";
        for (const auto& f : m_fractions) m_title += (m_title.empty() ? "" : " + ") + symbol_of(f.first) + " " + format_number(100.0 * f.second / total, 3) + "%";
        m_summary = "Color of the light: " + fc::to_hex(m_glow);
        for (const auto& f : m_fractions) element_lines(o, f.first, 100.0 * f.second / total);
        if (m_kind == "absorption") o << "#Absorption\nA cool gas of the same atoms in front of white light takes out exactly these wavelengths, leaving dark lines (how the Sun's composition was found).\n";
    }

    void build_identify(std::ostringstream& o, double step) {
        if (!m_profile.empty()) {
            std::vector<std::pair<double, double>> samples;
            const std::size_t n = m_profile.size();
            for (std::size_t i = 0; i < n; ++i) samples.emplace_back(m_profile_start + (m_profile_end - m_profile_start) * (n > 1 ? static_cast<double>(i) / (n - 1) : 0.0), m_profile[i]);
            const double lo = std::min(m_profile_start, m_profile_end), hi = std::max(m_profile_start, m_profile_end);
            m_start = std::max(m_start, lo);
            m_end = std::min(m_end, hi);
            if (m_end <= m_start + 5.0) {
                m_start = lo;
                m_end = hi;
            }
            m_observed = fc::spectrum_from_samples(samples, m_start, m_end, step);
            double floor_v = 1e300;
            for (double v : m_observed.intensity) floor_v = std::min(floor_v, v);
            if (floor_v < 1e300 && floor_v > 0.0) {
                for (double& v : m_observed.intensity) v -= floor_v;
            }
            o << "#Measured spectrum\nFrom an image\t" << n << " columns, " << format_number(m_profile_start, 4) << " to " << format_number(m_profile_end, 4) << " nm\n";
        } else {
            bool values = false;
            const std::vector<std::pair<double, double>> pts = parse_pairs(m_observed_text, values);
            if (pts.empty()) throw std::invalid_argument("paste wavelengths (nm), one per line, or pairs of wavelength and intensity");
            bool dense = pts.size() >= 25;
            if (dense) {
                std::vector<double> xs;
                for (const auto& p : pts) xs.push_back(p.first);
                std::sort(xs.begin(), xs.end());
                double mx = 0.0;
                for (std::size_t i = 1; i < xs.size(); ++i) mx = std::max(mx, xs[i] - xs[i - 1]);
                dense = mx < 10.0;
            }
            if (dense) {
                m_observed = fc::spectrum_from_samples(pts, m_start, m_end, step);
                o << "#Measured spectrum\nSamples\t" << pts.size() << " points\n";
            } else {
                m_observed = fc::spectrum_from_lines(pts, m_fwhm, m_start, m_end, step);
                o << "#Measured lines\n";
                for (const auto& p : pts) o << format_number(p.first, 6) << " nm\t" << band_name(p.first) << (values ? ", intensity " + format_number(p.second, 4) : std::string()) << "\n";
            }
        }
        double area = 0.0;
        for (double v : m_observed.intensity) area += v * m_observed.step_nm;
        if (area > 0.0) {
            for (double& v : m_observed.intensity) v /= area;
        }
        m_fit = fc::fit_emission_spectrum(m_observed, candidate_list(), m_fwhm, 1.0);
        if (!m_fit.ok) throw std::invalid_argument(m_fit.message);
        m_glow = fc::spectrum_color(fc::resample_spectrum(m_observed, 380.0, 780.0, 0.5));
        o << "#Best match\n";
        for (const auto& p : m_fit.percentages) o << fc::element(p.first).name << " (" << symbol_of(p.first) << ")\t" << format_number(p.second, 3) << "%\n";
        if (finite(m_fit.r_squared)) o << "Fit quality R\xC2\xB2\t" << format_number(m_fit.r_squared, 4) << (m_fit.r_squared > 0.9 ? " (good)" : m_fit.r_squared > 0.6 ? " (fair: some lines aren't explained)" : " (poor: lines from something not in the list, or a bad calibration)") << "\n";
        o << "Candidates tried\t" << m_fit.candidates.size() << " elements\n";
        o << "Color of the light\t" << fc::to_hex(m_glow) << "\n";
        o << "The percentages are each element's share of the atoms if they all glowed equally per atom; for real samples, calibrate against known mixtures. Elements with lines at nearly the same wavelengths can trade places, so narrow the candidates if you know what might be there.\n";
        m_title = "";
        for (std::size_t i = 0; i < m_fit.percentages.size() && i < 5; ++i) m_title += (i ? ", " : "") + symbol_of(m_fit.percentages[i].first) + " " + format_number(m_fit.percentages[i].second, 3) + "%";
        m_summary = m_title;
        for (const auto& p : m_fit.percentages) element_lines(o, p.first, p.second);
        std::vector<std::pair<int, double>> fr;
        for (const auto& p : m_fit.percentages) fr.push_back(p);
        m_fractions = fr;
    }

    void build_color(std::ostringstream& o) {
        fc::Srgb8 target;
        if (!fc::parse_hex_color(m_color_text, target)) throw std::invalid_argument("give a color like #ff8040");
        std::vector<int> cands = candidate_list();
        if (cands.empty()) {
            for (int z = 1; z <= 118; ++z) {
                if (fc::has_flame_color(z)) cands.push_back(z);
            }
        }
        m_color_fit = fc::fit_color_mixture(target, cands, m_color_flame);
        if (!m_color_fit.ok) throw std::invalid_argument(m_color_fit.message);
        m_target = target;
        m_glow = m_color_fit.matched;
        o << "#Color match\n";
        o << "Target\t" << fc::to_hex(target) << "\n";
        o << "Best mix\t" << fc::to_hex(m_color_fit.matched) << " (off by " << format_number(m_color_fit.distance, 3) << " of 441)\n";
        for (const auto& p : m_color_fit.percentages) o << fc::element(p.first).name << " (" << symbol_of(p.first) << ")\t" << format_number(p.second, 3) << "% of the light\n";
        if (!m_color_fit.unique) o << "!a color has only three numbers (red, green, blue), so with more than three candidates many mixes look the same; this is one of them. A spectrum pins it down; or list just the elements you suspect.\n";
        o << "Matched against each element's " << (m_color_flame ? "flame-test color" : "glow color worked out from its lines") << ".\n";
        m_title = fc::to_hex(target);
        m_summary = "";
        for (const auto& p : m_color_fit.percentages) m_summary += (m_summary.empty() ? "" : ", ") + symbol_of(p.first) + " " + format_number(p.second, 3) + "%";
        m_fractions = m_color_fit.percentages;
        const fc::Spectrum s = fc::emission_spectrum(m_fractions, m_fwhm, m_start, m_end, std::max(0.05, (m_end - m_start) / 1600.0));
        m_shown = s;
        for (const auto& p : m_color_fit.percentages) element_lines(o, p.first, p.second);
    }

    void draw_strip(Canvas& canvas, const Theme& theme, double x0, double y0, double width, double height) {
        canvas.fill_rect_px(x0, y0, width, height, rgb(0, 0, 0));
        const fc::Spectrum& s = m_mode == "identify" ? m_observed : m_shown;
        double peak = 0.0;
        for (double v : s.intensity) peak = std::max(peak, v);
        const bool absorb = m_kind == "absorption" && m_mode == "mix";
        const int cols = static_cast<int>(width);
        for (int i = 0; i < cols; ++i) {
            const double nm = m_start + (m_end - m_start) * (i + 0.5) / std::max(1, cols);
            const fc::Srgb8 c = fc::wavelength_color(nm);
            double v = s.size() > 1 ? s.at(nm) : 0.0;
            if (!absorb) v = peak > 0.0 ? std::min(1.0, std::pow(std::max(0.0, v) / peak, 0.5)) : 0.0;
            const double fade = nm < 380.0 || nm > 750.0 ? 0.35 : 1.0;
            const double k = std::max(0.0, std::min(1.0, v)) * fade;
            canvas.fill_rect_px(x0 + i, y0, 1.0, height, rgb(static_cast<std::uint8_t>(c.r * k), static_cast<std::uint8_t>(c.g * k), static_cast<std::uint8_t>(c.b * k)));
        }
        canvas.stroke_rect_px(x0, y0, width, height, theme.axes, 1.0);
        canvas.text(Vec2{ x0, y0 - 6.0 }, absorb ? "what you'd see through a spectroscope (dark lines on a rainbow)" : "what you'd see through a spectroscope", 11.0, theme.muted);
    }

    void draw_line_labels(Canvas& canvas, const Theme& theme) {
        std::vector<std::pair<double, std::string>> labels;
        double total = 0.0;
        for (const auto& f : m_fractions) total += f.second;
        for (const auto& f : m_fractions) {
            if (total > 0.0 && f.second / total < 0.05) continue;
            std::vector<fc::EmissionLine> lines = fc::emission_lines_of(f.first, m_start, m_end);
            std::sort(lines.begin(), lines.end(), [](const fc::EmissionLine& a, const fc::EmissionLine& b) { return a.intensity > b.intensity; });
            for (std::size_t i = 0; i < lines.size() && i < 3; ++i) labels.emplace_back(lines[i].wavelength_nm, symbol_of(f.first) + " " + format_number(lines[i].wavelength_nm, 4));
        }
        std::sort(labels.begin(), labels.end());
        double row_end[3] = { -1e300, -1e300, -1e300 };
        for (const auto& l : labels) {
            const Vec2 q = m_frame.to_screen(l.first, m_frame.ymax);
            const double half = 0.5 * canvas.text_width(l.second, 10.0) + 3.0;
            for (int row = 0; row < 3; ++row) {
                if (q.x - half < row_end[row]) continue;
                row_end[row] = q.x + half;
                canvas.text(Vec2{ q.x, m_frame.y0 + 12.0 + row * 12.0 }, l.second, 10.0, to_color(fc::wavelength_color(l.first)), TextAlign::Center);
                break;
            }
        }
        (void)theme;
    }

    void draw_swatches(Canvas& canvas, const Theme& theme, double x, double y, double size) {
        const double s = std::min(size, 80.0);
        if (m_mode == "color") {
            canvas.text(Vec2{ x, y - 6.0 }, "target", 11.0, theme.muted);
            canvas.fill_rect_px(x, y, s, s * 0.6, to_color(m_target));
            canvas.text(Vec2{ x, y + s * 0.6 + 16.0 }, "best mix", 11.0, theme.muted);
            canvas.fill_rect_px(x, y + s * 0.6 + 22.0, s, s * 0.6, to_color(m_glow));
            return;
        }
        canvas.text(Vec2{ x, y - 6.0 }, "glow", 11.0, theme.muted);
        canvas.fill_rect_px(x, y, s, s * 0.6, to_color(m_glow));
        canvas.text(Vec2{ x, y + s * 0.6 + 14.0 }, fc::to_hex(m_glow), 10.0, theme.muted);
        double yy = y + s * 0.6 + 34.0;
        for (const auto& f : m_fractions) {
            if (!fc::has_flame_color(f.first) || yy > canvas.height() - 30.0) continue;
            canvas.fill_rect_px(x, yy, 16.0, 16.0, to_color(fc::flame_color(f.first)));
            canvas.text(Vec2{ x + 22.0, yy + 12.0 }, symbol_of(f.first) + " flame", 10.0, theme.muted);
            yy += 22.0;
        }
    }

    ChemShared& m_shared;
    std::string m_mode = "mix";
    std::string m_kind = "emission";
    std::string m_mix_text = "Na 50, Sr 30, Cu 20";
    std::string m_observed_text;
    std::string m_candidates_text;
    std::string m_color_text = "#ff6040";
    bool m_color_flame = true;
    std::vector<double> m_profile;
    double m_profile_start = 380.0, m_profile_end = 780.0;
    double m_start = 380.0, m_end = 780.0, m_fwhm = 1.0;
    std::vector<std::pair<int, double>> m_fractions;
    fc::Spectrum m_shown, m_observed;
    fc::SpectrumFit m_fit;
    fc::ColorFit m_color_fit;
    fc::Srgb8 m_glow, m_target;
    std::string m_title, m_readouts, m_error, m_summary;
    PlotFrame m_frame;
    Vec2 m_mouse{ -1.0, -1.0 };
};

} // namespace

std::unique_ptr<ChemView> make_spectra_view(ChemShared& shared) {
    return std::make_unique<SpectraView>(shared);
}

} // namespace sapling

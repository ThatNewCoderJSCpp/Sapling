#include "sapling/canvas.hpp"
#include "sapling/simulation.hpp"

#include <grove/canvas.hpp>
#include <grove/dom_tools.hpp>
#include <grove/grove.hpp>
#include <grove/runtime/runtime.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#ifndef SAPLING_VERSION_STRING
#define SAPLING_VERSION_STRING "0.1.0"
#endif

namespace fs = std::filesystem;
using namespace sapling;

namespace {

std::string env(const char* name) {
    const char* v = std::getenv(name);
    return v ? std::string(v) : std::string();
}

std::string clean(std::string s) {
    for (char& c : s) {
        if (c == '\t' || c == '\n' || c == '\r') c = ' ';
    }
    return s;
}

std::string arg(const std::vector<std::string>& args, std::size_t i) {
    return i < args.size() ? args[i] : std::string();
}

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == sep) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    out.push_back(cur);
    return out;
}

std::vector<std::string> fields_of(const std::vector<std::string>& args) {
    if (args.size() == 1) return split(args[0], '\t');
    return args;
}

fs::path data_dir() {
    const std::string over = env("SAPLING_DATA_DIR");
    if (!over.empty()) return over;
#if defined(_WIN32)
    const std::string appdata = env("APPDATA");
    if (!appdata.empty()) return fs::path(appdata) / "Sapling";
#else
    const std::string xdg = env("XDG_DATA_HOME");
    if (!xdg.empty()) return fs::path(xdg) / "sapling";
    const std::string home = env("HOME");
    if (!home.empty()) return fs::path(home) / ".local" / "share" / "sapling";
#endif
    return fs::current_path() / "sapling-data";
}

std::string find_ui() {
    std::vector<fs::path> candidates;
    const std::string over = env("SAPLING_UI_DIR");
    if (!over.empty()) candidates.push_back(over);
    const fs::path exe = fs::path(grove::executable_path()).parent_path();
    candidates.push_back(exe / "ui");
    candidates.push_back(exe.parent_path() / "share" / "sapling" / "ui");
#ifdef SAPLING_SOURCE_UI_DIR
    candidates.push_back(SAPLING_SOURCE_UI_DIR);
#endif
    for (const fs::path& c : candidates) {
        std::error_code ec;
        if (fs::is_regular_file(c / "sapling.grove", ec)) return (c / "sapling.grove").string();
    }
    return std::string();
}

class App {
public:
    explicit App(grove::Project& project) : m_doc(project.document) {
        register_natives();
    }

    void install() {
        m_canvas = m_doc.find_by_id("view");
        if (m_canvas) {
            grove::canvas_of(m_doc, *m_canvas);
            m_doc.events().add_listener(m_canvas, "canvas_resize", [this](grove::Event&) { m_dirty = true; });
        }
        m_last = std::chrono::steady_clock::now();
        m_doc.events().set_interval(16.0, [this] { tick(); });
    }

private:
    void native(const std::string& name, std::function<std::string(const std::vector<std::string>&)> fn) {
        grove::register_native(m_doc, "sapling." + name, std::move(fn));
    }

    void register_natives() {
        native("version", [](const std::vector<std::string>&) { return std::string(SAPLING_VERSION_STRING); });
        native("sections", [](const std::vector<std::string>&) {
            std::string out;
            for (const Section& s : sections()) out += s.id + "\t" + clean(s.title) + "\t" + clean(s.description) + "\n";
            return out;
        });
        native("simulations", [](const std::vector<std::string>& a) {
            std::string out;
            for (const SimulationInfo& s : simulations()) {
                if (!arg(a, 0).empty() && s.section != arg(a, 0)) continue;
                out += s.id + "\t" + clean(s.title) + "\t" + clean(s.description) + "\n";
            }
            return out;
        });
        native("open", [this](const std::vector<std::string>& a) {
            const SimulationInfo* info = find_simulation(arg(a, 0));
            if (!info) return std::string();
            m_sim = make_simulation(info->id);
            if (!m_sim) return std::string();
            m_info = info;
            m_playing = true;
            m_speed = 1.0;
            m_dirty = true;
            return clean(info->title) + "\t" + clean(info->description) + "\t" + clean(m_sim->note());
        });
        native("close", [this](const std::vector<std::string>&) {
            m_sim.reset();
            m_info = nullptr;
            return std::string("ok");
        });
        native("params", [this](const std::vector<std::string>&) {
            std::string out;
            if (!m_sim) return out;
            for (const Parameter& p : m_sim->parameters()) {
                out += p.id + "\t" + clean(p.label) + "\t" + format_number(p.value, 10) + "\t" + format_number(p.min, 10) + "\t" + format_number(p.max, 10) + "\t" + format_number(p.step, 10) + "\t" + clean(p.unit) + "\n";
            }
            return out;
        });
        native("set_param", [this](const std::vector<std::string>& args) {
            const auto f = fields_of(args);
            if (!m_sim) return std::string();
            char* end = nullptr;
            const std::string text = arg(f, 1);
            const double v = std::strtod(text.c_str(), &end);
            if (end == text.c_str()) return std::string();
            m_sim->set_parameter(arg(f, 0), v);
            m_dirty = true;
            return format_number(m_sim->parameter(arg(f, 0)), 10);
        });
        native("play", [this](const std::vector<std::string>& a) {
            m_playing = arg(a, 0) == "1";
            m_last = std::chrono::steady_clock::now();
            return std::string(m_playing ? "1" : "0");
        });
        native("playing", [this](const std::vector<std::string>&) { return std::string(m_playing ? "1" : "0"); });
        native("reset", [this](const std::vector<std::string>&) {
            if (m_sim) m_sim->reset();
            m_dirty = true;
            return std::string("ok");
        });
        native("speed", [this](const std::vector<std::string>& a) {
            const double v = std::strtod(arg(a, 0).c_str(), nullptr);
            if (v > 0.0 && std::isfinite(v)) m_speed = v;
            return format_number(m_speed, 6);
        });
        native("step_once", [this](const std::vector<std::string>&) {
            if (m_sim) m_sim->step(m_sim->suggested_speed() * m_speed / 60.0);
            m_dirty = true;
            return std::string("ok");
        });
        native("readouts", [this](const std::vector<std::string>&) {
            std::string out;
            if (!m_sim) return out;
            for (const Readout& r : m_sim->readouts()) out += clean(r.label) + "\t" + clean(r.value) + "\n";
            return out;
        });
    }

    void tick() {
        const auto now = std::chrono::steady_clock::now();
        double dt = std::chrono::duration<double>(now - m_last).count();
        m_last = now;
        if (!m_sim || !m_canvas) return;
        dt = std::min(dt, 0.05);
        if (m_playing) {
            m_sim->step(dt * m_speed * m_sim->suggested_speed());
            m_dirty = true;
        }
        if (m_dirty) draw();
    }

    void draw() {
        unsigned w = 0, h = 0;
        if (!grove::canvas_size(m_doc, *m_canvas, w, h) || w == 0 || h == 0) return;
        m_dirty = false;
        if (static_cast<int>(w) != m_image.width() || static_cast<int>(h) != m_image.height()) m_image.resize(static_cast<int>(w), static_cast<int>(h));
        m_sim->draw(m_image, m_theme);
        grove::set_canvas_pixels(m_doc, *m_canvas, w, h, m_image.take_pixels());
    }

    grove::Document& m_doc;
    grove::Node* m_canvas = nullptr;
    std::unique_ptr<Simulation> m_sim;
    const SimulationInfo* m_info = nullptr;
    Canvas m_image;
    Theme m_theme;
    bool m_playing = true;
    bool m_dirty = true;
    double m_speed = 1.0;
    std::chrono::steady_clock::time_point m_last;
};

int bundle_ui(const std::string& output, bool portable) {
    const std::string ui = find_ui();
    if (ui.empty()) {
        std::cerr << "sapling: can't find the UI (sapling.grove) to bundle; set SAPLING_UI_DIR to the folder that holds it\n";
        return 1;
    }
    grove::AppBuildOptions options;
    options.output = output;
    options.stub = grove::executable_path();
    options.tools = grove::default_script_tools(ui);
    options.portable = portable;
    grove::Diagnostics diagnostics;
    const grove::AppBuildResult result = grove::build_app(ui, options, diagnostics, &std::cout);
    if (!result.ok) {
        diagnostics.print(std::cerr);
        std::cerr << "sapling: couldn't bundle the UI into " << output << "\n";
        return 1;
    }
    std::cout << "sapling: wrote " << result.output << " (" << result.files << " files, " << result.modules << " compiled scripts)\n";
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    bool check = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--version") {
            std::cout << "sapling " << SAPLING_VERSION_STRING << "\n";
            return 0;
        }
        if (a == "--help") {
            std::cout << "usage: sapling [--check] [--bundle OUTPUT [--portable]]\n\n"
                         "  --check          load the window and its script without opening it, then exit\n"
                         "  --bundle OUTPUT  write a copy of sapling with the window and its compiled script built in\n"
                         "  --portable       with --bundle, link the script with the static C++ runtime (when the compiler has one)\n\n"
                         "Set SAPLING_UI_DIR to use a different UI folder and SAPLING_DATA_DIR to keep data somewhere else.\n";
            return 0;
        }
        if (a == "--bundle") {
            if (i + 1 >= argc) {
                std::cerr << "sapling: --bundle needs an output path\n";
                return 2;
            }
            bool portable = false;
            for (int j = 1; j < argc; ++j) portable = portable || std::string(argv[j]) == "--portable";
            return bundle_ui(argv[i + 1], portable);
        }
        if (a == "--check") check = true;
    }

    grove::LoadedApp loaded;
    grove::Project from_files;
    grove::Project* project = nullptr;
    grove::runtime::RunOptions options;
    std::optional<grove::Bundle> bundle;
    if (env("SAPLING_UI_DIR").empty()) bundle = grove::read_embedded_bundle(grove::executable_path());
    if (bundle) {
        std::string error;
        if (!grove::load_app(*bundle, loaded, error)) {
            std::cerr << "sapling: " << error << "\n";
            return 1;
        }
        project = &loaded.project;
        options.window_id = loaded.window_id;
        options.prebuilt = &loaded.scripts;
        if (check) {
            std::cout << "sapling " << SAPLING_VERSION_STRING << ": built-in window ok, " << loaded.scripts.runnable_count() << " compiled script(s)\n";
            return loaded.scripts.ok() ? 0 : 1;
        }
    } else {
        const std::string ui = find_ui();
        if (ui.empty()) {
            std::cerr << "sapling: can't find the UI (sapling.grove); set SAPLING_UI_DIR to the folder that holds it\n";
            return 1;
        }
        from_files = grove::load_project(ui);
        if (!from_files.diagnostics.ok()) {
            from_files.diagnostics.print(std::cerr);
            return 1;
        }
        if (check) {
            std::cout << "sapling " << SAPLING_VERSION_STRING << ": window ok (" << ui << ")\n";
            return 0;
        }
        project = &from_files;
        options.tools = grove::default_script_tools(ui);
        options.tools.cache_dir = (data_dir() / "script-cache").string();
    }
    auto app = std::make_unique<App>(*project);
    app->install();
    return grove::runtime::run(*project, options);
}

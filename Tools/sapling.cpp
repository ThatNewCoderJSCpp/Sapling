#include "sapling/canvas.hpp"
#include "sapling/chemistry.hpp"
#include "sapling/simulation.hpp"

#include <grove/canvas.hpp>
#include <grove/dom_tools.hpp>
#include <grove/grove.hpp>
#include <grove/runtime/runtime.hpp>

#if defined(__ANDROID__)
#ifndef FIZMO
#define FIZMO
#endif
#ifndef FIZMO_GRAPHICS
#define FIZMO_GRAPHICS
#endif
#include <fizmo/includes.hpp>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
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

std::string panel_lines(const std::string& text) {
    std::string out;
    for (const std::string& line : split(text, '\n')) {
        if (line.empty()) continue;
        if (line[0] == '#') out += "h\t" + line.substr(1) + "\n";
        else if (line[0] == '!') out += "w\t" + line.substr(1) + "\n";
        else if (line.find('\t') != std::string::npos) out += "r\t" + line + "\n";
        else out += "n\t" + line + "\n";
    }
    return out;
}

std::vector<std::string> fields_of(const std::vector<std::string>& args) {
    if (args.size() == 1) return split(args[0], '\t');
    return args;
}

[[maybe_unused]] fs::path data_dir() {
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

[[maybe_unused]] std::string find_ui() {
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
        install_chemistry();
        m_last = std::chrono::steady_clock::now();
        m_doc.events().set_interval(16.0, [this] { tick(); });
    }

private:
    struct ChemCanvas {
        std::string name;
        grove::Node* node = nullptr;
        std::unique_ptr<ChemView> view;
        Canvas image;
        bool dirty = true;
    };

    void add_chem_view(const std::string& name, const std::string& id, std::unique_ptr<ChemView> view) {
        auto c = std::make_unique<ChemCanvas>();
        c->name = name;
        c->view = std::move(view);
        c->node = m_doc.find_by_id(id);
        m_chem_views[name] = std::move(c);
    }

    ChemCanvas* chem_at(double x, double y) {
        for (auto& kv : m_chem_views) {
            ChemCanvas& c = *kv.second;
            if (!c.node) continue;
            const grove::CanvasRect r = grove::canvas_rect(m_doc, *c.node);
            if (r.shown && x >= r.x && y >= r.y && x < r.x + r.width && y < r.y + r.height) return &c;
        }
        return nullptr;
    }

    Vec2 local(const ChemCanvas& c, double x, double y) const {
        const grove::CanvasRect r = grove::canvas_rect(m_doc, *c.node);
        unsigned w = 0, h = 0;
        grove::canvas_size(m_doc, *c.node, w, h);
        const double sx = r.width > 0.0 && w > 0 ? w / r.width : 1.0;
        const double sy = r.height > 0.0 && h > 0 ? h / r.height : 1.0;
        return Vec2{ (x - r.x) * sx, (y - r.y) * sy };
    }

    void install_chemistry() {
        add_chem_view("bond", "chem-bond-view", make_bonding_view(m_chem));
        add_chem_view("structure", "chem-ion-view", make_structure_view(m_chem));
        add_chem_view("phase", "chem-phase-view", make_phase_view(m_chem));
        add_chem_view("orbital", "chem-orbital-view", make_orbital_view(m_chem));
        add_chem_view("cloud", "chem-cloud-view", make_cloud_view(m_chem));
        add_chem_view("element", "chem-element-view", make_element_view(m_chem));
        for (auto& kv : m_chem_views) {
            ChemCanvas* c = kv.second.get();
            if (!c->node) continue;
            grove::canvas_of(m_doc, *c->node);
            m_doc.events().add_listener(c->node, "canvas_resize", [c](grove::Event&) { c->dirty = true; });
            m_doc.events().add_listener(c->node, "mouse_down", [this, c](grove::Event& e) {
                m_chem_drag = c;
                if (c->view->mouse_down(local(*c, e.x, e.y), e.button, e.modifiers)) c->dirty = true;
            });
            m_doc.events().add_listener(c->node, "wheel", [this, c](grove::Event& e) {
                if (c->view->wheel(local(*c, e.x, e.y), e.delta_y)) {
                    c->dirty = true;
                    e.default_prevented = true;
                }
            });
            m_doc.events().add_listener(c->node, "double_click", [this, c](grove::Event& e) {
                if (c->view->double_click(local(*c, e.x, e.y))) c->dirty = true;
            });
            m_doc.events().add_listener(c->node, "key_down", [c](grove::Event& e) {
                if (c->view->key(e.key, e.modifiers)) {
                    c->dirty = true;
                    e.default_prevented = true;
                }
            });
        }
        m_doc.events().add_listener(&m_doc.root(), "mouse_move", [this](grove::Event& e) {
            ChemCanvas* c = m_chem_drag ? m_chem_drag : chem_at(e.x, e.y);
            if (c != m_chem_hover && m_chem_hover && !m_chem_drag) {
                if (m_chem_hover->view->mouse_move(Vec2{ -1e6, -1e6 }, false)) m_chem_hover->dirty = true;
            }
            m_chem_hover = c;
            if (c && c->view->mouse_move(local(*c, e.x, e.y), m_chem_drag != nullptr)) c->dirty = true;
        });
        m_doc.events().add_listener(&m_doc.root(), "mouse_up", [this](grove::Event& e) {
            ChemCanvas* c = m_chem_drag;
            m_chem_drag = nullptr;
            if (c && c->view->mouse_up(local(*c, e.x, e.y), e.button)) c->dirty = true;
        });
    }

    void chem_dirty_all() {
        for (auto& kv : m_chem_views) kv.second->dirty = true;
    }

    void tick_chemistry(double dt) {
        for (auto& kv : m_chem_views) {
            ChemCanvas& c = *kv.second;
            if (!c.node) continue;
            const grove::CanvasRect r = grove::canvas_rect(m_doc, *c.node);
            if (!r.shown) continue;
            if (c.view->step(dt)) c.dirty = true;
            if (!c.dirty) continue;
            unsigned w = 0, h = 0;
            if (!grove::canvas_size(m_doc, *c.node, w, h) || w == 0 || h == 0) continue;
            c.dirty = false;
            if (static_cast<int>(w) != c.image.width() || static_cast<int>(h) != c.image.height()) c.image.resize(static_cast<int>(w), static_cast<int>(h));
            c.view->draw(c.image, m_theme);
            grove::set_canvas_pixels(m_doc, *c.node, w, h, c.image.take_pixels());
        }
    }

    void native(const std::string& name, std::function<std::string(const std::vector<std::string>&)> fn) {
        grove::register_native(m_doc, "sapling." + name, std::move(fn));
    }

    void register_natives() {
        native("version", [](const std::vector<std::string>&) { return std::string(SAPLING_VERSION_STRING); });
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
        native("chem.cmd", [this](const std::vector<std::string>& args) {
            const auto f = fields_of(args);
            const auto it = m_chem_views.find(arg(f, 0));
            if (it == m_chem_views.end()) return std::string("error\tno chemistry view called ") + arg(f, 0);
            const std::vector<std::string> rest(f.begin() + 1, f.end());
            std::string out;
            try {
                out = it->second->view->command(rest);
            } catch (const std::exception& e) {
                out = std::string("error\t") + e.what();
            }
            const std::string c = arg(rest, 0);
            if (c != "readouts" && c != "selection" && c != "version" && c != "state" && c != "props" && c != "saved") it->second->dirty = true;
            return out;
        });
        native("chem.panel", [this](const std::vector<std::string>& args) {
            const auto f = fields_of(args);
            const auto it = m_chem_views.find(arg(f, 0));
            if (it == m_chem_views.end()) return std::string("w\tno chemistry view called ") + arg(f, 0) + "\n";
            try {
                return panel_lines(it->second->view->command(std::vector<std::string>(f.begin() + 1, f.end())));
            } catch (const std::exception& e) {
                return std::string("w\t") + e.what() + "\n";
            }
        });
        native("chem.toolpanel", [this](const std::vector<std::string>& args) {
            try {
                return panel_lines(chemistry_tool(m_chem, fields_of(args)));
            } catch (const std::exception& e) {
                return std::string("w\t") + e.what() + "\n";
            }
        });
        native("chem.tool", [this](const std::vector<std::string>& args) {
            const auto f = fields_of(args);
            std::string out;
            try {
                out = chemistry_tool(m_chem, f);
            } catch (const std::exception& e) {
                out = std::string("error\t") + e.what();
            }
            if (arg(f, 0) == "units") chem_dirty_all();
            return out;
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
        dt = std::min(dt, 0.1);
        tick_chemistry(dt);
        if (!m_canvas) return;
        if (m_sim && m_playing) {
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
        if (m_sim) {
            m_sim->draw(m_image, m_theme);
        } else {
            draw_empty();
        }
        grove::set_canvas_pixels(m_doc, *m_canvas, w, h, m_image.take_pixels());
    }

    void draw_empty() {
        m_image.clear(m_theme.background);
        m_image.set_view(0.0, 0.0, 60.0);
        m_image.grid(1.0, m_theme.grid);
        m_image.axes(m_theme.axes);
    }

    grove::Document& m_doc;
    grove::Node* m_canvas = nullptr;
    ChemShared m_chem;
    std::map<std::string, std::unique_ptr<ChemCanvas>> m_chem_views;
    ChemCanvas* m_chem_drag = nullptr;
    ChemCanvas* m_chem_hover = nullptr;
    std::unique_ptr<Simulation> m_sim;
    Canvas m_image;
    Theme m_theme;
    bool m_playing = true;
    bool m_dirty = true;
    double m_speed = 1.0;
    std::chrono::steady_clock::time_point m_last;
};

#if !defined(__ANDROID__)
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

#endif

} // namespace

#if !defined(__ANDROID__)
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
#endif

#if defined(__ANDROID__)
namespace {

int show_start_error(fizmo::android::App& app, const std::string& message) {
    fizmo::android::log(fizmo::android::LogLevel::Error, "sapling", message);
    fizmo::graphics::Renderer renderer(app.surface());
    renderer.bind();
    fizmo::android::Event e;
    while (app.running()) {
        bool redraw = false;
        while (app.poll(e, 250)) {
            if (e.type == fizmo::android::EventType::Back) app.finish();
            redraw = true;
        }
        if (!redraw || !renderer.is_bound()) continue;
        const float k = app.density();
        const fizmo::android::Insets in = app.content_insets();
        const int x = in.left + static_cast<int>(24 * k);
        const int y = in.top + static_cast<int>(32 * k);
        const unsigned w = app.surface().width() > static_cast<unsigned>(x * 2) ? app.surface().width() - static_cast<unsigned>(x * 2) : app.surface().width();
        renderer.begin_frame();
        renderer.clear(fizmo::colors::Color(30, 30, 36));
        renderer.draw_text(x, y, w, static_cast<unsigned>(60 * k), std::string("Sapling couldn't start"), fizmo::text::TextStyle(22.0 * k, fizmo::colors::Color(243, 139, 168)));
        renderer.draw_text(x, y + static_cast<int>(48 * k), w, app.surface().height(), message, fizmo::text::TextStyle(15.0 * k, fizmo::colors::Color(205, 214, 244)));
        renderer.present();
    }
    return 1;
}

} // namespace

int fizmo_android_main(fizmo::android::App& app) {
    fizmo::android::redirect_stdio("sapling");
    const std::string data = app.internal_data_path().empty() ? std::string("/data/local/tmp") : app.internal_data_path();
    std::error_code ec;
    fs::create_directories(fs::path(data) / "data", ec);
    setenv("HOME", data.c_str(), 1);
    setenv("XDG_DATA_HOME", (fs::path(data) / "share").string().c_str(), 1);
    setenv("XDG_CONFIG_HOME", (fs::path(data) / "config").string().c_str(), 1);
    setenv("XDG_CACHE_HOME", (fs::path(data) / "cache").string().c_str(), 1);
    setenv("TMPDIR", app.cache_path().c_str(), 1);
    setenv("SAPLING_DATA_DIR", (fs::path(data) / "data").string().c_str(), 1);

    const std::optional<std::string> text = app.assets().read("grove/app.bundle");
    if (!text) return show_start_error(app, "this APK has no grove/app.bundle asset; build it again with grove build");
    std::optional<grove::Bundle> bundle = grove::parse_bundle(*text);
    if (!bundle) return show_start_error(app, "grove/app.bundle is damaged; build Sapling again");

    grove::AppLoadOptions options;
    options.root = data;
    options.cache_dir = app.cache_path();
    options.module_dir = app.native_library_dir();
    grove::LoadedApp loaded;
    std::string error;
    if (!grove::load_app(*bundle, loaded, error, options)) return show_start_error(app, error);

    auto shell = std::make_unique<App>(loaded.project);
    shell->install();
    grove::runtime::RunOptions run;
    run.window_id = loaded.window_id;
    run.scripts = true;
    return grove::runtime::run_android(loaded.project, app, &loaded.scripts, run);
}
#endif // __ANDROID__

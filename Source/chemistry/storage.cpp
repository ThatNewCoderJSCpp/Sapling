#include "sapling/chemistry.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <system_error>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

namespace sapling {

namespace fs = std::filesystem;

const JsonValue* JsonValue::find(const std::string& key) const {
    for (const auto& kv : object) {
        if (kv.first == key) return &kv.second;
    }
    return nullptr;
}

double JsonValue::number_or(const std::string& key, double fallback) const {
    const JsonValue* v = find(key);
    return v && v->type == Number ? v->number : fallback;
}

std::string JsonValue::string_or(const std::string& key, const std::string& fallback) const {
    const JsonValue* v = find(key);
    return v && v->type == String ? v->string : fallback;
}

namespace {

class JsonReader {
public:
    explicit JsonReader(const std::string& s) : m_s(s) {}

    bool value(JsonValue& out, int depth = 0) {
        if (depth > 64) return fail("nested too deeply");
        skip();
        if (m_i >= m_s.size()) return fail("unexpected end");
        const char c = m_s[m_i];
        if (c == '{') return object(out, depth);
        if (c == '[') return array(out, depth);
        if (c == '"') {
            out.type = JsonValue::String;
            return string(out.string);
        }
        if (word("true")) {
            out.type = JsonValue::Boolean;
            out.boolean = true;
            return true;
        }
        if (word("false")) {
            out.type = JsonValue::Boolean;
            out.boolean = false;
            return true;
        }
        if (word("null")) {
            out.type = JsonValue::Null;
            return true;
        }
        return number(out);
    }

    bool finish() {
        skip();
        return m_i == m_s.size() || fail("text after the value");
    }

    std::string error;

private:
    bool fail(const std::string& why) {
        if (error.empty()) error = why + " at " + std::to_string(m_i);
        return false;
    }

    void skip() {
        while (m_i < m_s.size() && std::isspace(static_cast<unsigned char>(m_s[m_i]))) ++m_i;
    }

    bool word(const char* w) {
        const std::size_t n = std::char_traits<char>::length(w);
        if (m_s.compare(m_i, n, w) != 0) return false;
        m_i += n;
        return true;
    }

    bool number(JsonValue& out) {
        const char* begin = m_s.c_str() + m_i;
        char* end = nullptr;
        const double v = std::strtod(begin, &end);
        if (end == begin) return fail("expected a value");
        m_i += static_cast<std::size_t>(end - begin);
        out.type = JsonValue::Number;
        out.number = v;
        return true;
    }

    static void put_utf8(std::string& out, unsigned cp) {
        if (cp < 0x80) {
            out += static_cast<char>(cp);
        } else if (cp < 0x800) {
            out += static_cast<char>(0xC0 | (cp >> 6));
            out += static_cast<char>(0x80 | (cp & 0x3F));
        } else if (cp < 0x10000) {
            out += static_cast<char>(0xE0 | (cp >> 12));
            out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
            out += static_cast<char>(0x80 | (cp & 0x3F));
        } else {
            out += static_cast<char>(0xF0 | (cp >> 18));
            out += static_cast<char>(0x80 | ((cp >> 12) & 0x3F));
            out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
            out += static_cast<char>(0x80 | (cp & 0x3F));
        }
    }

    bool hex4(unsigned& cp) {
        if (m_i + 4 > m_s.size()) return fail("bad escape");
        cp = 0;
        for (int k = 0; k < 4; ++k) {
            const char h = m_s[m_i++];
            cp <<= 4;
            if (h >= '0' && h <= '9') cp |= static_cast<unsigned>(h - '0');
            else if (h >= 'a' && h <= 'f') cp |= static_cast<unsigned>(h - 'a' + 10);
            else if (h >= 'A' && h <= 'F') cp |= static_cast<unsigned>(h - 'A' + 10);
            else return fail("bad escape");
        }
        return true;
    }

    bool string(std::string& out) {
        ++m_i;
        out.clear();
        while (m_i < m_s.size()) {
            const char c = m_s[m_i++];
            if (c == '"') return true;
            if (c != '\\') {
                out += c;
                continue;
            }
            if (m_i >= m_s.size()) break;
            const char e = m_s[m_i++];
            switch (e) {
                case 'n': out += '\n'; break;
                case 't': out += '\t'; break;
                case 'r': out += '\r'; break;
                case 'b': out += '\b'; break;
                case 'f': out += '\f'; break;
                case 'u': {
                    unsigned cp = 0;
                    if (!hex4(cp)) return false;
                    if (cp >= 0xD800 && cp < 0xDC00 && m_i + 1 < m_s.size() && m_s[m_i] == '\\' && m_s[m_i + 1] == 'u') {
                        m_i += 2;
                        unsigned lo = 0;
                        if (!hex4(lo)) return false;
                        cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                    }
                    put_utf8(out, cp);
                    break;
                }
                default: out += e; break;
            }
        }
        return fail("unterminated string");
    }

    bool array(JsonValue& out, int depth) {
        ++m_i;
        out.type = JsonValue::Array;
        skip();
        if (m_i < m_s.size() && m_s[m_i] == ']') {
            ++m_i;
            return true;
        }
        while (true) {
            JsonValue v;
            if (!value(v, depth + 1)) return false;
            out.array.push_back(std::move(v));
            skip();
            if (m_i < m_s.size() && m_s[m_i] == ',') {
                ++m_i;
                continue;
            }
            if (m_i < m_s.size() && m_s[m_i] == ']') {
                ++m_i;
                return true;
            }
            return fail("expected , or ]");
        }
    }

    bool object(JsonValue& out, int depth) {
        ++m_i;
        out.type = JsonValue::Object;
        skip();
        if (m_i < m_s.size() && m_s[m_i] == '}') {
            ++m_i;
            return true;
        }
        while (true) {
            skip();
            if (m_i >= m_s.size() || m_s[m_i] != '"') return fail("expected a key");
            std::string key;
            if (!string(key)) return false;
            skip();
            if (m_i >= m_s.size() || m_s[m_i] != ':') return fail("expected :");
            ++m_i;
            JsonValue v;
            if (!value(v, depth + 1)) return false;
            out.object.emplace_back(std::move(key), std::move(v));
            skip();
            if (m_i < m_s.size() && m_s[m_i] == ',') {
                ++m_i;
                continue;
            }
            if (m_i < m_s.size() && m_s[m_i] == '}') {
                ++m_i;
                return true;
            }
            return fail("expected , or }");
        }
    }

    const std::string& m_s;
    std::size_t m_i = 0;
};

std::string env_value(const char* name) {
    const char* v = std::getenv(name);
    return v ? std::string(v) : std::string();
}

bool valid_kind(const std::string& kind) {
    const auto kinds = saved_chemistry_kinds();
    return std::find(kinds.begin(), kinds.end(), kind) != kinds.end();
}

void mark_hidden(const fs::path& p) {
#if defined(_WIN32)
    const DWORD attrs = GetFileAttributesW(p.wstring().c_str());
    if (attrs != INVALID_FILE_ATTRIBUTES) SetFileAttributesW(p.wstring().c_str(), attrs | FILE_ATTRIBUTE_HIDDEN);
#else
    (void)p;
#endif
}

} // namespace

bool parse_json(const std::string& text, JsonValue& out, std::string* error) {
    JsonReader r(text);
    out = JsonValue{};
    const bool ok = r.value(out) && r.finish();
    if (!ok && error) *error = r.error;
    return ok;
}

std::string json_quote(const std::string& s) {
    std::string out = "\"";
    for (char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\t': out += "\\t"; break;
            case '\r': out += "\\r"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof buf, "\\u%04x", static_cast<unsigned>(static_cast<unsigned char>(c)));
                    out += buf;
                } else {
                    out += c;
                }
        }
    }
    out += '"';
    return out;
}

fs::path chemistry_save_root() {
    const std::string over = env_value("SAPLING_SAVE_DIR");
    if (!over.empty()) return fs::path(over) / "chemistry";
#if defined(__ANDROID__)
    const std::string data = env_value("SAPLING_DATA_DIR");
    if (!data.empty()) return fs::path(data) / ".sapling" / "chemistry";
#endif
    std::error_code ec;
    fs::path base = fs::current_path(ec);
    if (ec) base = fs::path(".");
    return base / ".sapling" / "chemistry";
}

std::vector<std::string> saved_chemistry_kinds() {
    return { "covalent", "ionic", "metallic" };
}

std::string safe_file_name(const std::string& name) {
    std::string out;
    for (unsigned char c : name) {
        if (c >= 0x80 || std::isalnum(c) || c == '-' || c == '_' || c == '.' || c == ' ' || c == '(' || c == ')' || c == '+') out += static_cast<char>(c);
        else out += '_';
    }
    while (!out.empty() && (out.back() == ' ' || out.back() == '.')) out.pop_back();
    while (!out.empty() && (out.front() == ' ' || out.front() == '.')) out.erase(out.begin());
    if (out.size() > 80) out.resize(80);
    return out;
}

bool save_chemistry(const std::string& kind, const std::string& name, const std::string& json, std::string* error) {
    if (!valid_kind(kind)) {
        if (error) *error = "unknown kind '" + kind + "'";
        return false;
    }
    const std::string file = safe_file_name(name);
    if (file.empty()) {
        if (error) *error = "give it a name first";
        return false;
    }
    const fs::path root = chemistry_save_root();
    const fs::path dir = root / kind;
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) {
        if (error) *error = "can't make " + dir.string() + ": " + ec.message();
        return false;
    }
    mark_hidden(root.parent_path());
    std::ofstream out(dir / (file + ".json"), std::ios::binary | std::ios::trunc);
    if (!out) {
        if (error) *error = "can't write " + (dir / (file + ".json")).string();
        return false;
    }
    out << json << "\n";
    return static_cast<bool>(out);
}

std::vector<std::string> list_saved_chemistry(const std::string& kind) {
    std::vector<std::string> out;
    if (!valid_kind(kind)) return out;
    std::error_code ec;
    const fs::path dir = chemistry_save_root() / kind;
    if (!fs::is_directory(dir, ec)) return out;
    for (const auto& entry : fs::directory_iterator(dir, ec)) {
        if (!entry.is_regular_file(ec) || entry.path().extension() != ".json") continue;
        const auto u = entry.path().stem().u8string();
        out.push_back(std::string(u.begin(), u.end()));
    }
    std::sort(out.begin(), out.end());
    return out;
}

bool load_chemistry(const std::string& kind, const std::string& name, std::string& json) {
    if (!valid_kind(kind)) return false;
    std::ifstream in(chemistry_save_root() / kind / (safe_file_name(name) + ".json"), std::ios::binary);
    if (!in) return false;
    std::ostringstream ss;
    ss << in.rdbuf();
    json = ss.str();
    return true;
}

bool delete_saved_chemistry(const std::string& kind, const std::string& name) {
    if (!valid_kind(kind)) return false;
    std::error_code ec;
    return fs::remove(chemistry_save_root() / kind / (safe_file_name(name) + ".json"), ec);
}

} // namespace sapling

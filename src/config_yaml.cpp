#include "config_yaml.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

// Generic YAML tree

struct Node {
    enum class Kind { Scalar, Map, List } kind = Kind::Scalar;
    int line = 0; // 1-based source line, for error messages
    std::string scalar;
    std::vector<std::pair<std::string, Node>> map; // keeps file order
    std::vector<Node> list;
};

[[noreturn]] void fail(int line, const std::string& msg)
{
    throw std::invalid_argument("yaml line " + std::to_string(line) + ": " + msg);
}

std::string trim(const std::string& s)
{
    const auto b = s.find_first_not_of(" \t");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t");
    return s.substr(b, e - b + 1);
}

struct SourceLine {
    int number;
    int indent;
    std::string text; // without indentation and comment
};

// Strips comments and blank lines, measures indentation
std::vector<SourceLine> split_lines(const std::string& text)
{
    std::vector<SourceLine> lines;
    std::size_t pos = 0;
    int number = 0;
    while (pos <= text.size()) {
        std::size_t end = text.find('\n', pos);
        if (end == std::string::npos) end = text.size();
        std::string raw = text.substr(pos, end - pos);
        pos = end + 1;
        ++number;

        if (!raw.empty() && raw.back() == '\r') raw.pop_back();
        // A comment starts at '#' at the beginning of the line or after whitespace
        for (std::size_t i = 0; i < raw.size(); ++i) {
            if (raw[i] == '#' && (i == 0 || raw[i - 1] == ' ' || raw[i - 1] == '\t')) {
                raw.erase(i);
                break;
            }
        }
        const auto first = raw.find_first_not_of(' ');
        if (first == std::string::npos || trim(raw).empty()) continue;
        if (raw[first] == '\t') fail(number, "tabs are not allowed for indentation");
        lines.push_back({number, static_cast<int>(first), trim(raw)});
    }
    return lines;
}

Node make_scalar(const std::string& s, int line)
{
    if (s.empty()) fail(line, "missing value");
    if (std::string("\"'&*|>[]{}!%@`").find(s.front()) != std::string::npos)
        fail(line, "unsupported YAML syntax: " + s);
    Node n;
    n.line = line;
    n.scalar = s;
    return n;
}

// Splits "key: rest" (or "key:") into key and trimmed rest
std::pair<std::string, std::string> split_key(const std::string& s, int line)
{
    const auto colon = s.find(':');
    if (colon == std::string::npos || (colon + 1 < s.size() && s[colon + 1] != ' '))
        fail(line, "expected 'key: value', got: " + s);
    std::string key = trim(s.substr(0, colon));
    if (key.empty()) fail(line, "empty key");
    return {key, trim(s.substr(colon + 1))};
}

void add_entry(Node& map, std::string key, Node value, int line)
{
    for (const auto& entry : map.map)
        if (entry.first == key) fail(line, "duplicate key '" + key + "'");
    map.map.emplace_back(std::move(key), std::move(value));
}

// Parses a single-line flow map "{a: 1, b: x}"
Node parse_flow_map(const std::string& s, int line)
{
    if (s.size() < 2 || s.front() != '{' || s.back() != '}')
        fail(line, "unsupported flow map: " + s);
    Node n;
    n.kind = Node::Kind::Map;
    n.line = line;
    const std::string inner = s.substr(1, s.size() - 2);
    if (trim(inner).empty()) return n;

    std::size_t pos = 0;
    while (pos <= inner.size()) {
        std::size_t comma = inner.find(',', pos);
        if (comma == std::string::npos) comma = inner.size();
        auto [key, value] = split_key(trim(inner.substr(pos, comma - pos)), line);
        add_entry(n, key, make_scalar(value, line), line);
        pos = comma + 1;
    }
    return n;
}

Node parse_value(const std::string& s, int line)
{
    if (!s.empty() && s.front() == '{') return parse_flow_map(s, line);
    return make_scalar(s, line);
}

bool is_list_item(const std::string& s) { return s == "-" || s.rfind("- ", 0) == 0; }

// Parses the block starting at lines[i] with the given indentation
Node parse_block(const std::vector<SourceLine>& lines, std::size_t& i, int indent)
{
    Node n;
    n.line = lines[i].number;
    const bool is_list = is_list_item(lines[i].text);
    n.kind = is_list ? Node::Kind::List : Node::Kind::Map;

    while (i < lines.size() && lines[i].indent >= indent) {
        const SourceLine& l = lines[i];
        if (l.indent > indent) fail(l.number, "unexpected indentation");
        if (is_list_item(l.text) != is_list) fail(l.number, "cannot mix list items and keys");

        if (is_list) {
            n.list.push_back(parse_value(trim(l.text.substr(1)), l.number));
            ++i;
            continue;
        }

        auto [key, rest] = split_key(l.text, l.number);
        ++i;
        if (!rest.empty()) {
            add_entry(n, key, parse_value(rest, l.number), l.number);
        } else {
            // Nested block: lists may sit at the same indentation as their key
            const bool nested = i < lines.size() && (lines[i].indent > indent ||
                                                     (lines[i].indent == indent &&
                                                      is_list_item(lines[i].text)));
            if (!nested) fail(l.number, "missing value for '" + key + "'");
            Node child = parse_block(lines, i, lines[i].indent);
            add_entry(n, key, std::move(child), l.number);
        }
    }
    return n;
}

Node parse_yaml(const std::string& text)
{
    const std::vector<SourceLine> lines = split_lines(text);
    if (lines.empty()) throw std::invalid_argument("yaml: document is empty");
    if (lines[0].indent != 0) fail(lines[0].number, "unexpected indentation");
    std::size_t i = 0;
    Node root = parse_block(lines, i, 0);
    if (i != lines.size()) fail(lines[i].number, "unexpected content");
    return root;
}

// Typed access

const Node* find(const Node& map, const std::string& key)
{
    for (const auto& entry : map.map)
        if (entry.first == key) return &entry.second;
    return nullptr;
}

void expect_kind(const Node& n, Node::Kind kind, const std::string& what)
{
    static const char* names[] = {"a value", "a map", "a list"};
    if (n.kind != kind) fail(n.line, what + " must be " + names[static_cast<int>(kind)]);
}

const Node& require(const Node& map, const std::string& key, const std::string& where)
{
    const Node* n = find(map, key);
    if (!n) fail(map.line, where + ": missing key '" + key + "'");
    return *n;
}

// Rejects keys not in allowed, so typos fail loudly
void check_keys(const Node& map, const std::vector<std::string>& allowed, const std::string& where)
{
    for (const auto& entry : map.map)
        if (std::find(allowed.begin(), allowed.end(), entry.first) == allowed.end())
            fail(entry.second.line, where + ": unknown key '" + entry.first + "'");
}

double to_double(const Node& n, const std::string& what)
{
    expect_kind(n, Node::Kind::Scalar, what);
    try {
        std::size_t used = 0;
        const double v = std::stod(n.scalar, &used);
        if (used == n.scalar.size()) return v;
    } catch (const std::exception&) {
    }
    fail(n.line, what + " must be a number, got '" + n.scalar + "'");
}

int to_int(const Node& n, const std::string& what)
{
    expect_kind(n, Node::Kind::Scalar, what);
    try {
        std::size_t used = 0;
        const int v = std::stoi(n.scalar, &used);
        if (used == n.scalar.size()) return v;
    } catch (const std::exception&) {
    }
    fail(n.line, what + " must be an integer, got '" + n.scalar + "'");
}

bool to_bool(const Node& n, const std::string& what)
{
    expect_kind(n, Node::Kind::Scalar, what);
    if (n.scalar == "true") return true;
    if (n.scalar == "false") return false;
    fail(n.line, what + " must be true or false, got '" + n.scalar + "'");
}

std::string to_string(const Node& n, const std::string& what)
{
    expect_kind(n, Node::Kind::Scalar, what);
    return n.scalar;
}

double get_double(const Node& map, const std::string& key, const std::string& where)
{
    return to_double(require(map, key, where), key);
}

SlacklinerConfig parse_slackliner(const Node& n)
{
    expect_kind(n, Node::Kind::Map, "slackliner");
    check_keys(n, {"m", "l_leg", "l_leash", "x_coor"}, "slackliner");
    SlacklinerConfig s;
    s.m = get_double(n, "m", "slackliner");
    s.l_leg = get_double(n, "l_leg", "slackliner");
    s.l_leash = get_double(n, "l_leash", "slackliner");
    s.x_coor = get_double(n, "x_coor", "slackliner");
    return s;
}

SegmentConfig parse_segment(const Node& n, const std::string& where)
{
    expect_kind(n, Node::Kind::Map, where);
    check_keys(n, {"main", "backup", "L_main", "L_backup", "break_mainline"}, where);
    SegmentConfig s;
    s.main = to_string(require(n, "main", where), "main");
    s.backup = to_string(require(n, "backup", where), "backup");
    s.L_main = get_double(n, "L_main", where);
    s.L_backup = get_double(n, "L_backup", where);
    if (const Node* b = find(n, "break_mainline")) s.break_mainline = to_bool(*b, "break_mainline");
    return s;
}

} // namespace

WebbingCatalog parse_webbings(const std::string& yaml_text)
{
    const Node root = parse_yaml(yaml_text);
    expect_kind(root, Node::Kind::Map, "webbing catalogue");

    WebbingCatalog webbings;
    for (const auto& [name, n] : root.map) {
        expect_kind(n, Node::Kind::Map, "webbing '" + name + "'");
        check_keys(n, {"stretch_pct", "tension_kN", "weight_g_m"}, name);
        webbings[name] = Webbing{get_double(n, "stretch_pct", name), get_double(n, "tension_kN", name),
                                 get_double(n, "weight_g_m", name)};
    }
    return webbings;
}

SimConfig parse_config(const std::string& yaml_text)
{
    const Node root = parse_yaml(yaml_text);
    expect_kind(root, Node::Kind::Map, "config");
    check_keys(root,
               {"plots", "gif", "L", "N", "T", "tension_kN", "pull_webbing", "pull_side",
                "slackliner", "segments"},
               "config");

    SimConfig cfg;
    if (const Node* n = find(root, "plots")) cfg.plots = to_bool(*n, "plots");
    if (const Node* n = find(root, "gif")) cfg.gif = to_bool(*n, "gif");
    cfg.L = get_double(root, "L", "config");
    cfg.N = to_int(require(root, "N", "config"), "N");
    cfg.T = get_double(root, "T", "config");
    if (const Node* n = find(root, "tension_kN")) cfg.tension_kN = to_double(*n, "tension_kN");
    if (const Node* n = find(root, "pull_webbing")) cfg.pull_webbing = to_double(*n, "pull_webbing");
    if (const Node* n = find(root, "pull_side")) cfg.pull_side = to_string(*n, "pull_side");
    cfg.slackliner = parse_slackliner(require(root, "slackliner", "config"));

    const Node& segments = require(root, "segments", "config");
    expect_kind(segments, Node::Kind::List, "segments");
    for (std::size_t i = 0; i < segments.list.size(); ++i)
        cfg.segments.push_back(
            parse_segment(segments.list[i], "segments[" + std::to_string(i) + "]"));

    validate(cfg);
    return cfg;
}

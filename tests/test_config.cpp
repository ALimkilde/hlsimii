// Tests for validate, build_line and the YAML front end.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "config.h"
#include "config_yaml.h"

namespace {

int failures = 0;

#define CHECK(cond)                                                               \
    do {                                                                          \
        if (!(cond)) {                                                            \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            ++failures;                                                           \
        }                                                                         \
    } while (0)

bool close(double a, double b) { return std::abs(a - b) <= 1e-9 * std::max(1.0, std::abs(b)); }

bool throws(const std::function<void()>& f)
{
    try {
        f();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

std::string read_file(const std::string& relative_path)
{
    std::ifstream in(std::string(HLSIM_SOURCE_DIR) + "/" + relative_path);
    if (!in) throw std::runtime_error("cannot open " + relative_path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// A small valid config, modified by the error tests
const std::string base_config = R"(
L: 10
N: 5
T: 1
tension_kN: 1
slackliner:
  m: 80
  l_leg: 1
  l_leash: 1
  x_coor: 5
segments:
  - {main: stretchy, backup: stiff, L_main: 10, L_backup: 11}
)";

std::string replace(std::string s, const std::string& from, const std::string& to)
{
    s.replace(s.find(from), from.size(), to);
    return s;
}

// Fixed test data lives in tests/data/; examples/ and data/ are free to edit.
void test_webbings()
{
    const WebbingCatalog w = parse_webbings(read_file("tests/data/webbings.yaml"));
    CHECK(w.size() == 2);
    CHECK(close(w.at("stretchy").stretch_pct, 10.0));
    CHECK(close(w.at("stiff").tension_kN, 10.0));
    CHECK(close(w.at("stiff").weight_g_m, 40.0));
}

void test_fixtures()
{
    const SimConfig t = parse_config(read_file("tests/data/tension.yaml"));
    CHECK(t.plots && !t.gif);
    CHECK(close(t.L, 100.0));
    CHECK(t.N == 35);
    CHECK(close(t.T, 15.0));
    CHECK(t.tension_kN && close(*t.tension_kN, 1.35));
    CHECK(!t.pull_webbing);
    CHECK(t.pull_side == "left");
    CHECK(close(t.slackliner.m, 89.0));
    CHECK(close(t.slackliner.x_coor, 50.0));
    CHECK(t.segments.size() == 3);
    CHECK(t.segments[0].break_mainline);
    CHECK(!t.segments[1].break_mainline);
    CHECK(!t.segments[2].break_mainline); // omitted: defaults to false
    CHECK(t.segments[2].main == "stiff" && t.segments[2].backup == "stretchy");
    CHECK(close(t.segments[2].L_main, 40.0) && close(t.segments[2].L_backup, 43.0));

    const SimConfig p = parse_config(read_file("tests/data/pull.yaml"));
    CHECK(!p.plots && p.gif);
    CHECK(!p.tension_kN);
    CHECK(p.pull_webbing && close(*p.pull_webbing, 1.5));
    CHECK(p.pull_side == "right"); // omitted: default
    CHECK(close(p.slackliner.m, 75.0));
    CHECK(p.segments.size() == 1);
}

// The user-facing files must keep parsing, whatever values they hold
void test_shipped_files_parse()
{
    const WebbingCatalog w = parse_webbings(read_file("data/webbings.yaml"));
    CHECK(!w.empty());
    for (const char* path : {"examples/example.yaml", "examples/1km.yaml", "examples/short_and_static.yaml"}) {
        try {
            build_line(parse_config(read_file(path)), w);
        } catch (const std::exception& e) {
            std::cerr << path << ": " << e.what() << "\n";
            ++failures;
        }
    }
}

void test_build_line()
{
    const WebbingCatalog w = parse_webbings(read_file("tests/data/webbings.yaml"));
    const Line line = build_line(parse_config(read_file("tests/data/tension.yaml")), w);
    CHECK(line.num_segments() == 3);
    const Segment& s = line.segments()[0];
    CHECK(close(s.ea_main, w.at("stretchy").ea()));
    CHECK(close(s.ea_backup, w.at("stiff").ea()));
    CHECK(close(s.rho_main, w.at("stretchy").rho()));
    CHECK(close(s.rho_backup, w.at("stiff").rho()));
    CHECK(close(s.L_main, 30.0) && close(s.L_backup, 32.0));

    SimConfig cfg = parse_config(base_config);
    CHECK(!throws([&] { build_line(cfg, w); }));
    cfg.segments[0].backup = "rope";
    CHECK(throws([&] { build_line(cfg, w); }));
}

// SimConfig filled without YAML, as the WASM front end will do
void test_validate()
{
    SimConfig cfg;
    cfg.L = 10.0;
    cfg.N = 5;
    cfg.T = 1.0;
    cfg.pull_webbing = 0.5;
    cfg.segments.push_back({"stretchy", "stiff", 10.0, 11.0, false});
    validate(cfg);

    SimConfig both = cfg;
    both.tension_kN = 1.0;
    CHECK(throws([&] { validate(both); }));
    SimConfig neither = cfg;
    neither.pull_webbing.reset();
    CHECK(throws([&] { validate(neither); }));
    SimConfig side = cfg;
    side.pull_side = "middle";
    CHECK(throws([&] { validate(side); }));
    SimConfig empty = cfg;
    empty.segments.clear();
    CHECK(throws([&] { validate(empty); }));
    SimConfig few = cfg;
    few.N = 1;
    CHECK(throws([&] { validate(few); }));
}

void test_yaml_errors()
{
    parse_config(base_config);
    CHECK(throws([] { parse_config(replace(base_config, "N: 5\n", "")); }));        // missing key
    CHECK(throws([] { parse_config(replace(base_config, "N: 5", "N: 5.5")); }));     // not an int
    CHECK(throws([] { parse_config(replace(base_config, "T: 1", "T: one")); }));     // not a number
    CHECK(throws([] { parse_config(replace(base_config, "T: 1", "T: 1\nTT: 2")); })); // unknown key
    CHECK(throws([] { parse_config(replace(base_config, "T: 1", "T: 1\nT: 2")); }));  // duplicate
    CHECK(throws([] { parse_config(replace(base_config, "T: 1", "T: 1\npull_webbing: 1")); }));
    CHECK(throws([] { parse_config(replace(base_config, "T: 1", "T: 1\npull_side: up")); }));
    CHECK(throws([] { parse_config(replace(base_config, "m: 80", "m: \"80\"")); }));  // quotes
    CHECK(throws([] { parse_config(replace(base_config, "  m: 80", "   m: 80")); })); // indentation
    CHECK(throws([] { parse_config(replace(base_config, "L_backup: 11}", "L_backup: 11")); }));
    CHECK(throws([] { parse_config(replace(base_config, "L_main: 10", "L_mian: 10")); }));
    CHECK(throws([] { parse_config(""); }));
    CHECK(throws([] { parse_webbings("joker: {stretch_pct: 3.6, tension_kN: 5}"); }));
}

} // namespace

int main()
{
    test_webbings();
    test_fixtures();
    test_shipped_files_parse();
    test_build_line();
    test_validate();
    test_yaml_errors();

    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all tests passed\n";
    return 0;
}

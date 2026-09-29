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
  - {main: joker, backup: solid, L_main: 10, L_backup: 11}
)";

std::string replace(std::string s, const std::string& from, const std::string& to)
{
    s.replace(s.find(from), from.size(), to);
    return s;
}

void test_webbings()
{
    const WebbingCatalog w = parse_webbings(read_file("data/webbings.yaml"));
    CHECK(w.size() == 6);
    CHECK(close(w.at("joker").stretch_pct, 3.6));
    CHECK(close(w.at("hmpe").tension_kN, 10.0));
    CHECK(close(w.at("nylove").weight_g_m, 62.0));
}

void test_examples()
{
    const SimConfig ex = parse_config(read_file("examples/example.yaml"));
    CHECK(ex.plots && !ex.gif);
    CHECK(close(ex.L, 100.0));
    CHECK(ex.N == 35);
    CHECK(close(ex.T, 15.0));
    CHECK(ex.tension_kN && close(*ex.tension_kN, 1.35));
    CHECK(!ex.pull_webbing);
    CHECK(ex.pull_side == "right");
    CHECK(close(ex.slackliner.m, 89.0));
    CHECK(close(ex.slackliner.x_coor, 50.0));
    CHECK(ex.segments.size() == 3);
    CHECK(ex.segments[0].break_mainline);
    CHECK(!ex.segments[1].break_mainline);
    CHECK(ex.segments[2].main == "joker" && ex.segments[2].backup == "solid");
    CHECK(close(ex.segments[2].L_main, 40.0) && close(ex.segments[2].L_backup, 43.0));

    const SimConfig km = parse_config(read_file("examples/1km.yaml"));
    CHECK(km.segments.size() == 20);
    CHECK(km.N == 100);

    const SimConfig ss = parse_config(read_file("examples/short_and_static.yaml"));
    CHECK(!ss.tension_kN);
    CHECK(ss.pull_webbing && close(*ss.pull_webbing, 1.0));
    CHECK(!ss.plots && ss.gif);
    CHECK(ss.segments.size() == 1);
}

void test_build_line()
{
    const WebbingCatalog w = parse_webbings(read_file("data/webbings.yaml"));
    const Line line = build_line(parse_config(read_file("examples/example.yaml")), w);
    CHECK(line.num_segments() == 3);
    const Segment& s = line.segments()[0];
    CHECK(close(s.ea_main, w.at("joker").ea()));
    CHECK(close(s.ea_backup, w.at("solid").ea()));
    CHECK(close(s.rho_main, w.at("joker").rho()));
    CHECK(close(s.rho_backup, w.at("solid").rho()));
    CHECK(close(s.L_main, 30.0) && close(s.L_backup, 32.0));

    SimConfig cfg = parse_config(base_config);
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
    cfg.segments.push_back({"joker", "solid", 10.0, 11.0, false});
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
    test_examples();
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

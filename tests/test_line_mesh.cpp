// Tests for pull_webbing, elements_per_segment and discretize.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

#include "line.h"
#include "mesh.h"

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

double total_mass(const Line& line)
{
    double m = 0.0;
    for (const Segment& s : line.segments()) m += s.mass();
    return m;
}

bool same_segment(const Segment& a, const Segment& b)
{
    return a.L_main == b.L_main && a.L_backup == b.L_backup && a.ea_main == b.ea_main &&
           a.ea_backup == b.ea_backup && a.rho_main == b.rho_main && a.rho_backup == b.rho_backup;
}

const Webbing w1{5.0, 10.0, 50.0};
const Webbing w2{8.0, 10.0, 30.0};

// Segments with different backup/main ratios
Line three_segment_line()
{
    return Line({Segment(w1, w2, 10.0, 12.0), Segment(w1, w1, 0.4, 0.6), Segment(w2, w1, 7.3, 8.0)});
}

void test_pull_webbing()
{
    const Line line = three_segment_line();
    const Line pulled = pull_webbing(line, 1.0, 0.5);
    const auto& a = line.segments();
    const auto& b = pulled.segments();

    CHECK(b.size() == a.size());
    CHECK(same_segment(a[1], b[1])); // middle untouched

    CHECK(close(b[0].L_main, a[0].L_main - 1.0));
    CHECK(close(b[2].L_main, a[2].L_main - 0.5));
    CHECK(close(b[0].L_backup / b[0].L_main, a[0].L_backup / a[0].L_main));
    CHECK(close(b[2].L_backup / b[2].L_main, a[2].L_backup / a[2].L_main));
    CHECK(close(b[0].L_backup, 12.0 - 1.2)); // backup pulled 1.2x the main

    const double removed = a[0].rho_main * 1.0 + a[0].rho_backup * 1.0 * (12.0 / 10.0) +
                           a[2].rho_main * 0.5 + a[2].rho_backup * 0.5 * (8.0 / 7.3);
    CHECK(close(total_mass(pulled), total_mass(line) - removed));

    // No pull leaves the line unchanged
    const Line same = pull_webbing(line, 0.0, 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) CHECK(same_segment(a[i], same.segments()[i]));

    // Single segment: both pulls come off the same segment
    const Line single({Segment(w1, w2, 10.0, 11.0)});
    const Segment s = pull_webbing(single, 1.0, 2.0).segments()[0];
    CHECK(close(s.L_main, 7.0));
    CHECK(close(s.L_backup, 11.0 * 0.7));

    // Invalid pulls
    CHECK(throws([&] { pull_webbing(line, -0.1, 0.0); }));
    CHECK(throws([&] { pull_webbing(line, 0.0, 7.3); }));
    CHECK(throws([&] { pull_webbing(line, 10.5, 0.0); }));
    CHECK(throws([&] { pull_webbing(single, 5.0, 5.0); }));
    CHECK(throws([&] { pull_webbing(Line({}), 0.0, 0.0); }));
}

void test_elements_per_segment()
{
    const std::vector<int> n = elements_per_segment(three_segment_line(), 1.0);
    CHECK(n.size() == 3);
    CHECK(n[0] == 10); // exact multiple
    CHECK(n[1] == 1);  // shorter than h
    CHECK(n[2] == 7);  // rounded

    CHECK(elements_per_segment(three_segment_line(), 100.0) == std::vector<int>({1, 1, 1}));
    CHECK(throws([] { elements_per_segment(three_segment_line(), 0.0); }));
}

void test_discretize()
{
    const Line line = pull_webbing(three_segment_line(), 1.0, 0.5);
    const std::vector<int> n = elements_per_segment(line, 0.7);
    const DiscreteLine d = discretize(line, n);

    int num_el = 0;
    for (int ni : n) num_el += ni;
    CHECK(d.num_elements() == num_el);
    CHECK(d.num_nodes() == num_el + 1);
    CHECK(static_cast<int>(d.k_backup.size()) == num_el);
    CHECK(static_cast<int>(d.l_main.size()) == num_el);
    CHECK(static_cast<int>(d.l_backup.size()) == num_el);
    CHECK(static_cast<int>(d.element_segment.size()) == num_el);
    CHECK(!d.break_mainline);

    double node_mass_sum = 0.0;
    for (double m : d.node_mass) node_mass_sum += m;
    CHECK(close(node_mass_sum, total_mass(line)));

    for (int s = 0; s < line.num_segments(); ++s) {
        const Segment& seg = line.segments()[s];
        double lm = 0.0, lb = 0.0, compliance_main = 0.0, compliance_backup = 0.0;
        int count = 0;
        for (int e = 0; e < d.num_elements(); ++e) {
            if (d.element_segment[e] != s) continue;
            ++count;
            lm += d.l_main[e];
            lb += d.l_backup[e];
            compliance_main += 1.0 / d.k_main[e];
            compliance_backup += 1.0 / d.k_backup[e];
        }
        CHECK(count == n[s]);
        CHECK(close(lm, seg.L_main));
        CHECK(close(lb, seg.L_backup));
        // Springs in series recover the segment stiffness EA/L
        CHECK(close(1.0 / compliance_main, seg.ea_main / seg.L_main));
        CHECK(close(1.0 / compliance_backup, seg.ea_backup / seg.L_backup));
    }

    CHECK(throws([&] { discretize(line, {1, 1}); }));
    CHECK(throws([&] { discretize(line, {1, 0, 1}); }));
}

} // namespace

int main()
{
    test_pull_webbing();
    test_elements_per_segment();
    test_discretize();

    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all tests passed\n";
    return 0;
}

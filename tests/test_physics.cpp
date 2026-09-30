// Tests for LineModel::edge_force and LineModel::edge_forces.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <algorithm>
#include <iostream>

#include "physics.h"

namespace {

int failures = 0;

#define CHECK(cond)                                                               \
    do {                                                                          \
        if (!(cond)) {                                                            \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            ++failures;                                                           \
        }                                                                         \
    } while (0)

bool close(const Vec2& a, const Vec2& b) { return (a - b).norm() <= 1e-9 * std::max(1.0, b.norm()); }

// Three edges, main rest length 1.0 and backup rest length 1.2.
// The main stiffness differs per edge, so a mixed-up edge index gives a wrong force.
DiscreteLine three_edge_line()
{
    DiscreteLine d;
    d.k_main = {100.0, 200.0, 300.0};
    d.k_backup = {50.0, 50.0, 50.0};
    d.l_main = {1.0, 1.0, 1.0};
    d.l_backup = {1.2, 1.2, 1.2};
    d.node_mass = {0.5, 1.0, 1.0, 0.5};
    d.element_segment = {0, 0, 0};
    return d;
}

LineModel model_with_anchors(Vec2 left, Vec2 right) { return LineModel(three_edge_line(), left, right, Params{}); }

void test_edge_force()
{
    const LineModel m = model_with_anchors({0.0, 0.0}, {3.0, 0.0});
    const Vec2 origin(0.0, 0.0);

    // Shorter than both rest lengths: slack
    CHECK(close(m.edge_force(0, origin, {0.9, 0.0}), Vec2(0.0, 0.0)));
    // Exactly at the main rest length: still zero
    CHECK(close(m.edge_force(0, origin, {1.0, 0.0}), Vec2(0.0, 0.0)));
    // Main stretched 0.1, backup slack: 100 * 0.1
    CHECK(close(m.edge_force(0, origin, {1.1, 0.0}), Vec2(10.0, 0.0)));
    // Exactly at the backup rest length: main only, 100 * 0.2
    CHECK(close(m.edge_force(0, origin, {1.2, 0.0}), Vec2(20.0, 0.0)));
    // Both stretched: 100 * 0.5 + 50 * 0.3
    CHECK(close(m.edge_force(0, origin, {1.5, 0.0}), Vec2(65.0, 0.0)));

    // The direction follows the edge; the magnitude does not change
    CHECK(close(m.edge_force(0, origin, {0.0, 1.1}), Vec2(0.0, 10.0)));
    const Vec2 dir(0.6, 0.8);
    CHECK(close(m.edge_force(0, origin, 1.5 * dir), 65.0 * dir));
    // Only the relative position matters
    const Vec2 shift(-4.0, 7.0);
    CHECK(close(m.edge_force(0, shift, shift + 1.5 * dir), 65.0 * dir));
    // Swapping the endpoints flips the force
    CHECK(close(m.edge_force(0, 1.5 * dir, origin), -65.0 * dir));

    // The edge index selects the stiffness: 200 * 0.1 and 300 * 0.1
    CHECK(close(m.edge_force(1, origin, {1.1, 0.0}), Vec2(20.0, 0.0)));
    CHECK(close(m.edge_force(2, origin, {1.1, 0.0}), Vec2(30.0, 0.0)));
}

void test_edge_forces()
{
    const Mat2X v = Mat2X::Zero(2, 2);

    // Every edge at its main rest length: no force anywhere
    {
        LineModel m = model_with_anchors({0.0, 0.0}, {3.0, 0.0});
        Mat2X q(2, 2);
        q << 1.0, 2.0,
             0.0, 0.0;
        m.edge_forces(q, v);
        const Mat2X Fe = m.get_Fe();
        CHECK(Fe.cols() == 3);
        CHECK(Fe.isZero());
    }

    // Every edge stretched 0.1: the force scales with each edge's own stiffness
    {
        LineModel m = model_with_anchors({0.0, 0.0}, {3.3, 0.0});
        Mat2X q(2, 2);
        q << 1.1, 2.2,
             0.0, 0.0;
        m.edge_forces(q, v);
        const Mat2X Fe = m.get_Fe();
        CHECK(close(Fe.col(0), Vec2(10.0, 0.0)));
        CHECK(close(Fe.col(1), Vec2(20.0, 0.0)));
        CHECK(close(Fe.col(2), Vec2(30.0, 0.0)));
    }

    // Edge lengths 1.1, 1.5 and 1.0: checks that each column uses the right endpoints,
    // including the anchors on the first and last edge
    {
        LineModel m = model_with_anchors({0.0, 0.0}, {3.6, 0.0});
        Mat2X q(2, 2);
        q << 1.1, 2.6,
             0.0, 0.0;
        m.edge_forces(q, v);
        const Mat2X Fe = m.get_Fe();
        CHECK(close(Fe.col(0), Vec2(10.0, 0.0)));                 // 100 * 0.1
        CHECK(close(Fe.col(1), Vec2(200.0 * 0.5 + 50.0 * 0.3, 0.0)));
        CHECK(close(Fe.col(2), Vec2(0.0, 0.0)));                  // at rest length
    }

    // U shape: up from the left anchor, across, down to the right anchor.
    // Each Fe column points from the edge's left node to its right node.
    {
        LineModel m = model_with_anchors({0.0, 0.0}, {1.0, 0.0});
        Mat2X q(2, 2);
        q << 0.0, 1.0,
             1.1, 1.1;
        m.edge_forces(q, v);
        const Mat2X Fe = m.get_Fe();
        CHECK(close(Fe.col(0), Vec2(0.0, 10.0)));    // up, 100 * 0.1
        CHECK(close(Fe.col(1), Vec2(0.0, 0.0)));     // across, at rest length
        CHECK(close(Fe.col(2), Vec2(0.0, -30.0)));   // down, 300 * 0.1
    }
}

} // namespace

int main()
{
    test_edge_force();
    test_edge_forces();

    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all tests passed\n";
    return 0;
}

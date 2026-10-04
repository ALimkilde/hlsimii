// Tests for LineModel::edge_force and LineModel::assemble.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>

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

// Gravity is switched off, so the net forces below are pure spring forces
LineModel model_with_anchors(Vec2 left, Vec2 right)
{
    Params p;
    p.gravity.setZero();
    return LineModel(three_edge_line(), left, right, p);
}

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

// Finite-difference check of the edge Jacobian K = dF/dqb (and dF/dqa = -K).
// The central-difference error should shrink like h^2, so the observed order should be close to 2.
// Prints a convergence table for each case.
void check_edge_jacobian(const LineModel& m, Eigen::Index e, const Vec2& qa, const Vec2& qb, const char* name)
{
    Mat2 K;
    m.edge_force(e, qa, qb, &K);

    // A generic direction: not along the edge, so the geometric term is exercised
    const Vec2 d = Vec2(0.3, -0.7).normalized();
    const Vec2 Kd = K * d;

    std::cout << "edge Jacobian FD check: " << name << "\n";
    std::cout << "        h        err(qb)    order      err(qa)    order\n";

    double h = 1e-1;
    double prev_err_b = 0.0, prev_err_a = 0.0;
    double order_b = 0.0, order_a = 0.0;
    for (int i = 0; i < 6; ++i, h /= 2.0) {
        const Vec2 fd_b = (m.edge_force(e, qa, qb + h * d) - m.edge_force(e, qa, qb - h * d)) / (2.0 * h);
        const Vec2 fd_a = (m.edge_force(e, qa + h * d, qb) - m.edge_force(e, qa - h * d, qb)) / (2.0 * h);
        const double err_b = (fd_b - Kd).norm();
        const double err_a = (fd_a + Kd).norm();

        std::printf("  %9.2e  %10.3e", h, err_b);
        if (i > 0) {
            order_b = std::log2(prev_err_b / err_b);
            order_a = std::log2(prev_err_a / err_a);
            std::printf("  %6.3f", order_b);
        } else {
            std::printf("        ");
        }
        std::printf("   %10.3e", err_a);
        if (i > 0) std::printf("  %6.3f", order_a);
        std::printf("\n");

        prev_err_b = err_b;
        prev_err_a = err_a;
    }

    CHECK(std::abs(order_b - 2.0) < 0.1);
    CHECK(std::abs(order_a - 2.0) < 0.1);
}

void test_edge_jacobian()
{
    const LineModel m = model_with_anchors({0.0, 0.0}, {3.0, 0.0});
    const Vec2 shift(-0.4, 0.7);
    const Vec2 dir(0.6, 0.8);

    // Lengths are kept away from the rest lengths 1.0 and 1.2, where the force has a kink
    check_edge_jacobian(m, 0, shift, shift + 1.1 * dir, "main only (len 1.1)");
    check_edge_jacobian(m, 0, shift, shift + 1.5 * dir, "main + backup (len 1.5)");
    check_edge_jacobian(m, 2, shift, shift + 1.5 * dir, "edge 2, main + backup (len 1.5)");

    // Slack edge: no force and no stiffness
    Mat2 K = Mat2::Ones();
    m.edge_force(0, shift, shift + 0.9 * dir, &K);
    CHECK(K.isZero());
}

// Net force on each interior node from assemble(). Edge e pulls its left node towards its right node
// with f_e = edge_force(e, ...) and the right node back with -f_e, so interior node i gets
// f_{i+1} - f_i (nodes numbered from the left anchor, which is node 0).
void test_assemble()
{
    const Mat2X v = Mat2X::Zero(2, 2);

    // F starts as NaN: assemble must overwrite it, not add to stale values
    const double nan = std::numeric_limits<double>::quiet_NaN();
    auto fresh_F = [&] { return Mat2X::Constant(2, 2, nan); };

    // Every edge at its main rest length: no force anywhere
    {
        const LineModel m = model_with_anchors({0.0, 0.0}, {3.0, 0.0});
        Mat2X q(2, 2);
        q << 1.0, 2.0,
             0.0, 0.0;
        Mat2X F = fresh_F();
        m.assemble(q, v, F);
        CHECK(F.isZero());
    }

    // Every edge stretched 0.1: edge forces 10, 20, 30 along +x
    {
        const LineModel m = model_with_anchors({0.0, 0.0}, {3.3, 0.0});
        Mat2X q(2, 2);
        q << 1.1, 2.2,
             0.0, 0.0;
        Mat2X F = fresh_F();
        m.assemble(q, v, F);
        CHECK(close(F.col(0), Vec2(10.0, 0.0)));   // 20 - 10
        CHECK(close(F.col(1), Vec2(10.0, 0.0)));   // 30 - 20
    }

    // Edge lengths 1.1, 1.5 and 1.0: edge forces 10, 200 * 0.5 + 50 * 0.3 = 115 and 0.
    // Checks that each node uses the right edges, including the anchors on the first and last edge
    {
        const LineModel m = model_with_anchors({0.0, 0.0}, {3.6, 0.0});
        Mat2X q(2, 2);
        q << 1.1, 2.6,
             0.0, 0.0;
        Mat2X F = fresh_F();
        m.assemble(q, v, F);
        CHECK(close(F.col(0), Vec2(105.0, 0.0)));    // 115 - 10
        CHECK(close(F.col(1), Vec2(-115.0, 0.0)));   // 0 - 115
    }

    // U shape: up from the left anchor, across, down to the right anchor.
    // Edge forces (0, 10), 0 and (0, -30): both nodes are pulled down towards their anchor
    {
        const LineModel m = model_with_anchors({0.0, 0.0}, {1.0, 0.0});
        Mat2X q(2, 2);
        q << 0.0, 1.0,
             1.1, 1.1;
        Mat2X F = fresh_F();
        m.assemble(q, v, F);
        CHECK(close(F.col(0), Vec2(0.0, -10.0)));   // 0 - (0, 10)
        CHECK(close(F.col(1), Vec2(0.0, -30.0)));   // (0, -30) - 0
    }
}

} // namespace

int main()
{
    test_edge_force();
    test_edge_jacobian();
    test_assemble();

    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all tests passed\n";
    return 0;
}

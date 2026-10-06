// Tests for LineModel::edge_vector, edge_force, assemble, residual, both_neighboring_edges_slack
// and place/remove_slackliner.
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

// Edge e runs from node e-1 to node e, with the anchors standing in for node -1 and node q.cols()
void test_edge_vector()
{
    const LineModel m = model_with_anchors({0.0, 0.0}, {3.0, 0.5});
    Mat2X q(2, 2);
    q << 1.1, 2.3,
         0.2, -0.4;

    CHECK(close(m.edge_vector(q, 0), Vec2(1.1, 0.2)));    // left anchor -> node 0
    CHECK(close(m.edge_vector(q, 1), Vec2(1.2, -0.6)));   // node 0 -> node 1
    CHECK(close(m.edge_vector(q, 2), Vec2(0.7, 0.9)));    // node 1 -> right anchor
}

// edge_force takes the edge vector dq = qb - qa and returns the force pulling qa towards qb
void test_edge_force()
{
    const LineModel m = model_with_anchors({0.0, 0.0}, {3.0, 0.0});

    // Shorter than both rest lengths: slack
    CHECK(close(m.edge_force(0, {0.9, 0.0}), Vec2(0.0, 0.0)));
    // Exactly at the main rest length: still zero
    CHECK(close(m.edge_force(0, {1.0, 0.0}), Vec2(0.0, 0.0)));
    // Main stretched 0.1, backup slack: 100 * 0.1
    CHECK(close(m.edge_force(0, {1.1, 0.0}), Vec2(10.0, 0.0)));
    // Exactly at the backup rest length: main only, 100 * 0.2
    CHECK(close(m.edge_force(0, {1.2, 0.0}), Vec2(20.0, 0.0)));
    // Both stretched: 100 * 0.5 + 50 * 0.3
    CHECK(close(m.edge_force(0, {1.5, 0.0}), Vec2(65.0, 0.0)));

    // The direction follows the edge; the magnitude does not change
    CHECK(close(m.edge_force(0, {0.0, 1.1}), Vec2(0.0, 10.0)));
    const Vec2 dir(0.6, 0.8);
    CHECK(close(m.edge_force(0, 1.5 * dir), 65.0 * dir));
    // Reversing the edge flips the force
    CHECK(close(m.edge_force(0, -1.5 * dir), -65.0 * dir));

    // The edge index selects the stiffness: 200 * 0.1 and 300 * 0.1
    CHECK(close(m.edge_force(1, {1.1, 0.0}), Vec2(20.0, 0.0)));
    CHECK(close(m.edge_force(2, {1.1, 0.0}), Vec2(30.0, 0.0)));
}

// Finite-difference check of the edge Jacobian K = dF/d(dq).
// The central-difference error should shrink like h^2, so the observed order should be close to 2.
// Prints a convergence table for each case.
void check_edge_jacobian(const LineModel& m, Eigen::Index e, const Vec2& dq, const char* name)
{
    Mat2 K;
    m.edge_force(e, dq, &K);

    // A generic direction: not along the edge, so the geometric term is exercised
    const Vec2 d = Vec2(0.3, -0.7).normalized();
    const Vec2 Kd = K * d;

    std::cout << "edge Jacobian FD check: " << name << "\n";
    std::cout << "        h        err        order\n";

    double h = 1e-1;
    double prev_err = 0.0, order = 0.0;
    for (int i = 0; i < 6; ++i, h /= 2.0) {
        const Vec2 fd = (m.edge_force(e, dq + h * d) - m.edge_force(e, dq - h * d)) / (2.0 * h);
        const double err = (fd - Kd).norm();

        std::printf("  %9.2e  %10.3e", h, err);
        if (i > 0) {
            order = std::log2(prev_err / err);
            std::printf("  %6.3f", order);
        }
        std::printf("\n");

        prev_err = err;
    }

    CHECK(std::abs(order - 2.0) < 0.1);
}

void test_edge_jacobian()
{
    const LineModel m = model_with_anchors({0.0, 0.0}, {3.0, 0.0});
    const Vec2 dir(0.6, 0.8);

    // Lengths are kept away from the rest lengths 1.0 and 1.2, where the force has a kink
    check_edge_jacobian(m, 0, 1.1 * dir, "main only (len 1.1)");
    check_edge_jacobian(m, 0, 1.5 * dir, "main + backup (len 1.5)");
    check_edge_jacobian(m, 2, 1.5 * dir, "edge 2, main + backup (len 1.5)");

    // Slack edge: no force and no stiffness
    Mat2 K = Mat2::Ones();
    m.edge_force(0, 0.9 * dir, &K);
    CHECK(K.isZero());
}

// Net force on each interior node from assemble(). Edge e pulls its left node towards its right node
// with f_e = edge_force(e, ...) and the right node back with -f_e, so interior node i gets
// f_{i+1} - f_i (nodes numbered from the left anchor, which is node 0).
void test_assemble()
{

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
        m.assemble(q, F);
        CHECK(F.isZero());
    }

    // Every edge stretched 0.1: edge forces 10, 20, 30 along +x
    {
        const LineModel m = model_with_anchors({0.0, 0.0}, {3.3, 0.0});
        Mat2X q(2, 2);
        q << 1.1, 2.2,
             0.0, 0.0;
        Mat2X F = fresh_F();
        m.assemble(q, F);
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
        m.assemble(q, F);
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
        m.assemble(q, F);
        CHECK(close(F.col(0), Vec2(0.0, -10.0)));   // 0 - (0, 10)
        CHECK(close(F.col(1), Vec2(0.0, -30.0)));   // (0, -30) - 0
    }
}

// Five edges with distinct stiffnesses and distinct node masses, so a shifted edge or mass index
// changes the result. Rest lengths: main 1.0, backup 1.2.
DiscreteLine five_edge_line()
{
    DiscreteLine d;
    d.k_main = {100.0, 200.0, 300.0, 400.0, 500.0};
    d.k_backup = {50.0, 60.0, 70.0, 80.0, 90.0};
    d.l_main = {1.0, 1.0, 1.0, 1.0, 1.0};
    d.l_backup = {1.2, 1.2, 1.2, 1.2, 1.2};
    d.node_mass = {0.3, 1.1, 1.3, 1.7, 1.9, 0.7};
    d.element_segment = {0, 0, 0, 0, 0};
    return d;
}

// A zigzag through the plane with edge lengths 1.1, 0.9, 1.5, 1.3, 1.05, i.e.
// main only, slack, main + backup, main + backup, main only. All lengths stay away from the kinks.
struct Zigzag {
    Vec2 left, right;
    Mat2X q{2, 4};
};

Zigzag zigzag()
{
    const double len[5] = {1.1, 0.9, 1.5, 1.3, 1.05};
    const double angle[5] = {0.3, -0.5, 0.1, 0.8, -0.2};   // edge directions [rad]

    Zigzag z;
    z.left = Vec2(0.2, -0.1);
    Vec2 p = z.left;
    for (int e = 0; e < 5; ++e) {
        p += len[e] * Vec2(std::cos(angle[e]), std::sin(angle[e]));
        if (e < 4) z.q.col(e) = p;
    }
    z.right = p;
    return z;
}

// Fills every block with NaN, so a block that assemble forgets to write shows up in the check
void poison(SymmBlockTriMat& K)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < K.num_blocks(); ++i) K.diag(i) = Mat2::Constant(nan);
    for (int i = 0; i + 1 < K.num_blocks(); ++i) K.upper(i) = Mat2::Constant(nan);
}

// Finite-difference check of the assembled Jacobian K = dF/dq on the zigzag, with gravity on.
// Each column j of the dense K is compared to the central difference of F when the flat
// coordinate j of q is perturbed. The error should shrink like h^2.
void test_assemble_jacobian()
{
    const Zigzag z = zigzag();
    const LineModel m(five_edge_line(), z.left, z.right, Params{});
    const Eigen::Index n = z.q.size();   // 8 coordinates, ordered [x0 y0 x1 y1 ...] like the dense K

    SymmBlockTriMat K(4);
    poison(K);
    Mat2X F(2, 4);
    m.assemble(z.q, F, &K);
    const Mat J = K.to_dense();
    CHECK(J.allFinite());

    // Forces must not depend on whether the Jacobian was requested
    Mat2X F_no_K(2, 4);
    m.assemble(z.q, F_no_K);
    CHECK(F_no_K == F);

    std::cout << "assembled Jacobian FD check (zigzag, 4 free nodes)\n";
    std::cout << "        h        max err    order\n";

    double h = 1e-2;
    double prev_err = 0.0, order = 0.0;
    for (int step = 0; step < 5; ++step, h /= 2.0) {
        Mat J_fd(n, n);
        for (Eigen::Index j = 0; j < n; ++j) {
            Mat2X qp = z.q, qm = z.q;
            qp.data()[j] += h;
            qm.data()[j] -= h;
            Mat2X Fp(2, 4), Fm(2, 4);
            m.assemble(qp, Fp);
            m.assemble(qm, Fm);
            J_fd.col(j) = Eigen::Map<const Vec>((Fp - Fm).eval().data(), n) / (2.0 * h);
        }
        const double err = (J_fd - J).cwiseAbs().maxCoeff();

        std::printf("  %9.2e  %10.3e", h, err);
        if (step > 0) {
            order = std::log2(prev_err / err);
            std::printf("  %6.3f", order);
        }
        std::printf("\n");
        prev_err = err;
    }
    CHECK(std::abs(order - 2.0) < 0.1);
}

// With every edge slack the only force is gravity, so node i feels m_{i+1} * g.
// The node masses differ, which catches an off-by-one in the mass index (including the last node).
void test_assemble_gravity()
{
    Params p;
    p.gravity = Vec2(0.0, -10.0);
    const LineModel m(five_edge_line(), {0.0, 0.0}, {4.0, 0.0}, p);

    Mat2X q(2, 4);
    q << 0.8, 1.6, 2.4, 3.2,
         0.0, 0.0, 0.0, 0.0;   // edges 0.8 long: all slack
    Mat2X F(2, 4);
    SymmBlockTriMat K(4);
    poison(K);
    m.assemble(q, F, &K);

    CHECK(close(F.col(0), Vec2(0.0, -11.0)));
    CHECK(close(F.col(1), Vec2(0.0, -13.0)));
    CHECK(close(F.col(2), Vec2(0.0, -17.0)));
    CHECK(close(F.col(3), Vec2(0.0, -19.0)));

    // Slack everywhere: no stiffness
    CHECK(K.to_dense().isZero());
}

// Smallest valid line: two edges and one free node, so the loop runs once, the node is both
// first and last, and K has no upper block to write.
void test_assemble_single_node()
{
    DiscreteLine d;
    d.k_main = {100.0, 200.0};
    d.k_backup = {50.0, 50.0};
    d.l_main = {1.0, 1.0};
    d.l_backup = {1.2, 1.2};
    d.node_mass = {0.5, 1.0, 0.5};
    d.element_segment = {0, 0};

    Params p;
    p.gravity = Vec2(0.0, -10.0);
    const LineModel m(d, {0.0, 0.0}, {2.2, 0.0}, p);

    Mat2X q(2, 1);
    q << 1.1, 0.0;   // both edges stretched 0.1: forces 10 and 20 along +x
    Mat2X F(2, 1);
    SymmBlockTriMat K(1);
    poison(K);
    m.assemble(q, F, &K);

    CHECK(close(F.col(0), Vec2(10.0, -10.0)));   // (20 - 10, 1.0 * -10)

    // Collinear along x, main only: each edge has K = diag(k, T/len), and dF/dq = -(K_0 + K_1)
    Mat2 expected;
    expected << -(100.0 + 200.0), 0.0,
                0.0, -(10.0 / 1.1 + 20.0 / 1.1);
    CHECK(K.diag(0).isApprox(expected));
}

// residual() is assemble() on flat vectors ordered [x0 y0 x1 y1 ...].
// The entries are compared one by one, so a wrong layout (e.g. all x first, then all y) would fail.
void test_residual()
{
    const Zigzag z = zigzag();
    const LineModel m(five_edge_line(), z.left, z.right, Params{});

    Mat2X F(2, 4);
    SymmBlockTriMat K(4);
    m.assemble(z.q, F, &K);

    Vec q(8);
    for (int i = 0; i < 4; ++i) {
        q[2 * i] = z.q(0, i);
        q[2 * i + 1] = z.q(1, i);
    }

    Vec r(8);
    SymmBlockTriMat jac(4);
    poison(jac);
    m.residual(q, r, &jac);
    for (int i = 0; i < 4; ++i) {
        CHECK(r[2 * i] == F(0, i));
        CHECK(r[2 * i + 1] == F(1, i));
    }
    CHECK(jac.to_dense() == K.to_dense());

    // Without a Jacobian
    Vec r_no_jac(8);
    m.residual(q, r_no_jac);
    CHECK(r_no_jac == r);
}

// Lays the five edges of five_edge_line() out along x with the given lengths, and checks that
// the Mat2X and flat versions agree. Main rest length is 1.0, so 0.9 is slack and 1.1 is taut.
bool neighbors_slack(const double (&len)[5])
{
    Mat2X q(2, 4);
    double x = 0.0;
    for (int e = 0; e < 4; ++e) {
        x += len[e];
        q.col(e) = Vec2(x, 0.0);
    }
    const LineModel m(five_edge_line(), {0.0, 0.0}, {x + len[4], 0.0}, Params{});

    const bool result = m.both_neighboring_edges_slack(q);
    const Vec q_flat = Eigen::Map<const Vec>(q.data(), q.size());
    CHECK(m.both_neighboring_edges_slack_flat(q_flat) == result);
    return result;
}

// Only two slack edges next to each other count, including the pairs at either anchor
void test_both_neighboring_edges_slack()
{
    CHECK(!neighbors_slack({1.1, 1.1, 1.1, 1.1, 1.1}));   // all taut
    CHECK(!neighbors_slack({0.9, 1.1, 0.9, 1.1, 0.9}));   // slack, but never two in a row
    CHECK(neighbors_slack({0.9, 0.9, 1.1, 1.1, 1.1}));    // pair at the left anchor
    CHECK(neighbors_slack({1.1, 0.9, 0.9, 1.1, 1.1}));    // pair in the middle
    CHECK(neighbors_slack({1.1, 1.1, 1.1, 0.9, 0.9}));    // pair at the right anchor: last edge must be checked
}


// A slackliner placed on free node i adds mass * g to F.col(i) and nothing else. The line is all
// slack, so F is pure gravity. Every free node is tried, including the first and the last.
void test_place_slackliner()
{
    Params p;
    p.gravity = Vec2(0.0, -10.0);
    const LineModel m(five_edge_line(), {0.0, 0.0}, {4.0, 0.0}, p);
    const double mass = 80.0;

    Mat2X q(2, 4);
    q << 0.8, 1.6, 2.4, 3.2,
         0.0, 0.0, 0.0, 0.0;   // edges 0.8 long: all slack
    Mat2X F_bare(2, 4), F(2, 4);
    m.assemble(q, F_bare);

    for (int i = 0; i < 4; ++i) {
        m.place_slackliner(i, mass);
        m.assemble(q, F);
        for (int j = 0; j < 4; ++j) {
            const Vec2 extra = (j == i) ? Vec2(mass * p.gravity) : Vec2::Zero();
            if (!close(F.col(j), Vec2(F_bare.col(j) + extra)))
                std::cerr << "  slackliner on node " << i << ": wrong force on node " << j << "\n";
            CHECK(close(F.col(j), Vec2(F_bare.col(j) + extra)));
        }
    }

    // Removing restores the bare line, and a second remove is harmless
    m.remove_slackliner();
    m.assemble(q, F);
    CHECK(F == F_bare);
    m.remove_slackliner();
    m.assemble(q, F);
    CHECK(F == F_bare);
}

} // namespace

int main()
{
    test_edge_vector();
    test_edge_force();
    test_edge_jacobian();
    test_assemble();
    test_assemble_jacobian();
    test_assemble_gravity();
    test_assemble_single_node();
    test_residual();
    test_both_neighboring_edges_slack();
    test_place_slackliner();

    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all tests passed\n";
    return 0;
}

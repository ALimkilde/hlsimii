// Tests for LineModel::edge_vector, edge_force (incl. Kelvin-Voigt damping), edge_drag, assemble, residual, both_neighboring_edges_slack
// place/remove_slackliner and rhs.
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

// Zero velocities switch Kelvin-Voigt damping off, so the forces are the pure spring forces
const Vec2 dv_rest = Vec2::Zero();
Mat2X at_rest(const Mat2X& q) { return Mat2X::Zero(2, q.cols()); }

// Three edges, main rest length 1.0 and backup rest length 1.2.
// The main stiffness differs per edge, so a mixed-up edge index gives a wrong force.
DiscreteLine three_edge_line()
{
    DiscreteLine d;
    d.k_main = {100.0, 200.0, 300.0};
    d.k_backup = {50.0, 50.0, 50.0};
    d.l_main = {1.0, 1.0, 1.0};
    d.l_backup = {1.2, 1.2, 1.2};
    d.node_mass = Vec::Constant(2, 1.0);
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
    CHECK(close(m.edge_force(0, {0.9, 0.0}, dv_rest), Vec2(0.0, 0.0)));
    // Exactly at the main rest length: still zero
    CHECK(close(m.edge_force(0, {1.0, 0.0}, dv_rest), Vec2(0.0, 0.0)));
    // Main stretched 0.1, backup slack: 100 * 0.1
    CHECK(close(m.edge_force(0, {1.1, 0.0}, dv_rest), Vec2(10.0, 0.0)));
    // Exactly at the backup rest length: main only, 100 * 0.2
    CHECK(close(m.edge_force(0, {1.2, 0.0}, dv_rest), Vec2(20.0, 0.0)));
    // Both stretched: 100 * 0.5 + 50 * 0.3
    CHECK(close(m.edge_force(0, {1.5, 0.0}, dv_rest), Vec2(65.0, 0.0)));

    // The direction follows the edge; the magnitude does not change
    CHECK(close(m.edge_force(0, {0.0, 1.1}, dv_rest), Vec2(0.0, 10.0)));
    const Vec2 dir(0.6, 0.8);
    CHECK(close(m.edge_force(0, 1.5 * dir, dv_rest), 65.0 * dir));
    // Reversing the edge flips the force
    CHECK(close(m.edge_force(0, -1.5 * dir, dv_rest), -65.0 * dir));

    // The edge index selects the stiffness: 200 * 0.1 and 300 * 0.1
    CHECK(close(m.edge_force(1, {1.1, 0.0}, dv_rest), Vec2(20.0, 0.0)));
    CHECK(close(m.edge_force(2, {1.1, 0.0}, dv_rest), Vec2(30.0, 0.0)));
}

// Finite-difference check of the edge Jacobian K = dF/d(dq).
// The central-difference error should shrink like h^2, so the observed order should be close to 2.
// Prints a convergence table for each case.
void check_edge_jacobian(const LineModel& m, Eigen::Index e, const Vec2& dq, const char* name)
{
    Mat2 K;
    m.edge_force(e, dq, dv_rest, &K);

    // A generic direction: not along the edge, so the geometric term is exercised
    const Vec2 d = Vec2(0.3, -0.7).normalized();
    const Vec2 Kd = K * d;

    std::cout << "edge Jacobian FD check: " << name << "\n";
    std::cout << "        h        err        order\n";

    double h = 1e-1;
    double prev_err = 0.0, order = 0.0;
    for (int i = 0; i < 6; ++i, h /= 2.0) {
        const Vec2 fd = (m.edge_force(e, dq + h * d, dv_rest) - m.edge_force(e, dq - h * d, dv_rest)) / (2.0 * h);
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
    m.edge_force(0, 0.9 * dir, dv_rest, &K);
    CHECK(K.isZero());
}

// === Kelvin-Voigt damping ===
// Damping only acts while an edge retracts: the tension becomes T + c * dlen/dt with c = tau_KV * k,
// clipped so it never goes below 0. dlen/dt = dv . u, with dv = vb - va.

// Same as model_with_anchors, with a round tau_KV so the damping forces below are easy to check
LineModel damped_model(Vec2 left, Vec2 right)
{
    Params p;
    p.gravity.setZero();
    p.tau_KV = 0.1;
    return LineModel(three_edge_line(), left, right, p);
}

// Edge 0 at length 1.5: T = 100 * 0.5 + 50 * 0.3 = 65 and c = 0.1 * (100 + 50) = 15
void test_edge_force_damping()
{
    const LineModel m = damped_model({0.0, 0.0}, {3.0, 0.0});
    const Vec2 u(0.6, 0.8);
    const Vec2 dq = 1.5 * u;

    // Stretching: no damping, whatever the speed
    CHECK(close(m.edge_force(0, dq, 2.0 * u), 65.0 * u));
    CHECK(close(m.edge_force(0, dq, 100.0 * u), 65.0 * u));
    // Retracting at 2 m/s: 65 - 15 * 2
    CHECK(close(m.edge_force(0, dq, -2.0 * u), 35.0 * u));
    // Retracting fast: 65 - 15 * 10 < 0, clipped to no tension (webbing cannot push)
    CHECK(close(m.edge_force(0, dq, -10.0 * u), Vec2(0.0, 0.0)));
    // Sideways motion does not change the length, so no damping
    const Vec2 u_perp(-0.8, 0.6);
    CHECK(close(m.edge_force(0, dq, 3.0 * u_perp), 65.0 * u));
    // Slack edge: no force, even when retracting
    CHECK(close(m.edge_force(0, 0.9 * u, -2.0 * u), Vec2(0.0, 0.0)));
    // Main only, edge 2 at length 1.1: T = 300 * 0.1 = 30, c = 0.1 * 300 = 30, so 30 - 30 * 0.5
    CHECK(close(m.edge_force(2, 1.1 * u, -0.5 * u), 15.0 * u));
}

// For many velocities dv, compares the damped force f(dv) to the spring force f(0) on one edge:
//  - damping never adds energy: the power it puts into the two nodes, -(f(dv) - f(0)) . dv, is <= 0
//  - the tension never goes negative: f(dv) . u >= 0
//  - stretching is never damped: dv . u > 0 gives f(dv) == f(0)
void test_damping_dissipates()
{
    const LineModel m = damped_model({0.0, 0.0}, {3.0, 0.0});
    int checked = 0;

    for (const double len : {1.1, 1.5}) {               // main only, main + backup
        for (int a = 0; a < 8; ++a) {                     // edge direction
            const double th = 0.8 * a;
            const Vec2 u(std::cos(th), std::sin(th));
            const Vec2 f0 = m.edge_force(0, len * u, dv_rest);

            for (int b = 0; b < 16; ++b) {                // velocity direction
                for (const double speed : {0.01, 0.3, 2.0, 50.0}) {
                    const double phi = 0.4 * b;
                    const Vec2 dv = speed * Vec2(std::cos(phi), std::sin(phi));
                    const Vec2 f = m.edge_force(0, len * u, dv);

                    CHECK(-(f - f0).dot(dv) <= 1e-12);
                    CHECK(f.dot(u) >= -1e-12);
                    if (dv.dot(u) > 0.0) CHECK(close(f, f0));
                    ++checked;
                }
            }
        }
    }
    std::cout << "damping dissipation check: " << checked << " velocities\n";
}

// The force is continuous where an edge turns from stretching to retracting (dlen/dt = 0):
// a small retraction speed eps changes the force by only c * eps
void test_damping_continuous()
{
    const LineModel m = damped_model({0.0, 0.0}, {3.0, 0.0});
    const Vec2 u(0.6, 0.8);
    const Vec2 dq = 1.5 * u;   // T = 65, c = 15
    const Vec2 f0 = m.edge_force(0, dq, dv_rest);

    for (const double eps : {1e-3, 1e-6, 1e-9}) {
        CHECK(close(m.edge_force(0, dq, eps * u), f0));
        CHECK(std::abs((f0 - m.edge_force(0, dq, -eps * u)).norm() - 15.0 * eps) <= 1e-9);
    }
}

// Finite-difference check of both edge Jacobians with damping: Kq = df/d(dq) and Kv = df/d(dv).
// The cases cover every branch of edge_force, and dv always has a sideways part so the terms from
// the turning edge direction are exercised. Each Jacobian is also requested on its own, and must
// come out the same as when both are requested.
void test_edge_jacobian_damped()
{
    const LineModel m = damped_model({0.0, 0.0}, {3.0, 0.0});
    const Vec2 u(0.6, 0.8), u_perp(-0.8, 0.6);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double h = 1e-6;

    // Edge 0 at len 1.5: T = 65, c = 15, so damping clips at dlen/dt = -65/15 = -4.33
    // Edge 0 at len 1.1: T = 10, c = 10, so damping clips at dlen/dt = -1
    struct Case { const char* name; double len; Vec2 dv; };
    const Case cases[] = {
        {"stretching",              1.5,   2.0 * u + 1.0 * u_perp},
        {"damped, main + backup",   1.5,  -2.0 * u + 1.5 * u_perp},
        {"damped, mostly sideways", 1.5,  -0.7 * u - 3.0 * u_perp},
        {"damped, main only",       1.1,  -0.5 * u + 0.8 * u_perp},
        {"clipped",                 1.5, -10.0 * u + 1.0 * u_perp},
        {"slack, retracting",       0.9,  -2.0 * u + 1.0 * u_perp},
    };

    std::cout << "damped edge Jacobian FD check\n";
    std::cout << "  case                       max err Kq   max err Kv\n";

    for (const Case& c : cases) {
        const Vec2 dq = c.len * u;

        Mat2 Kq = Mat2::Constant(nan), Kv = Mat2::Constant(nan);
        m.edge_force(0, dq, c.dv, &Kq, &Kv);

        Mat2 Kq_fd, Kv_fd;
        for (int j = 0; j < 2; ++j) {
            const Vec2 d = h * Vec2::Unit(j);
            Kq_fd.col(j) = (m.edge_force(0, dq + d, c.dv) - m.edge_force(0, dq - d, c.dv)) / (2.0 * h);
            Kv_fd.col(j) = (m.edge_force(0, dq, c.dv + d) - m.edge_force(0, dq, c.dv - d)) / (2.0 * h);
        }
        const double err_q = (Kq - Kq_fd).cwiseAbs().maxCoeff();
        const double err_v = (Kv - Kv_fd).cwiseAbs().maxCoeff();
        std::printf("  %-25s  %10.3e   %10.3e\n", c.name, err_q, err_v);

        // NaN left in a Jacobian makes the error NaN, which fails these checks too
        CHECK(err_q <= 1e-6 * std::max(1.0, Kq_fd.norm()));
        CHECK(err_v <= 1e-6 * std::max(1.0, Kv_fd.norm()));

        Mat2 Kq_only = Mat2::Constant(nan), Kv_only = Mat2::Constant(nan);
        m.edge_force(0, dq, c.dv, &Kq_only);
        m.edge_force(0, dq, c.dv, nullptr, &Kv_only);
        CHECK(Kq_only == Kq);
        CHECK(Kv_only == Kv);
    }
}

// assemble passes each edge its velocity difference, with the anchors standing still.
// Every edge stretched 0.1 (spring forces 10, 20, 30, c = 10, 20, 30).
// Node 0 moves left and node 1 moves right:
//   edge 0: dv = v0 - 0  = -0.5 -> retracting, 10 - 10 * 0.5 = 5
//   edge 1: dv = v1 - v0 = +1.0 -> stretching, 20
//   edge 2: dv = 0 - v1  = -0.5 -> retracting, 30 - 30 * 0.5 = 15
void test_assemble_damping()
{
    const LineModel m = damped_model({0.0, 0.0}, {3.3, 0.0});
    Mat2X q(2, 2), v(2, 2), F(2, 2);
    q << 1.1, 2.2,
         0.0, 0.0;
    v << -0.5, 0.5,
          0.0, 0.0;
    m.assemble(q, v, F);
    CHECK(close(F.col(0), Vec2(15.0, 0.0)));    // 20 - 5
    CHECK(close(F.col(1), Vec2(-5.0, 0.0)));    // 15 - 20
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
        m.assemble(q, at_rest(q), F);
        CHECK(F.isZero());
    }

    // Every edge stretched 0.1: edge forces 10, 20, 30 along +x
    {
        const LineModel m = model_with_anchors({0.0, 0.0}, {3.3, 0.0});
        Mat2X q(2, 2);
        q << 1.1, 2.2,
             0.0, 0.0;
        Mat2X F = fresh_F();
        m.assemble(q, at_rest(q), F);
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
        m.assemble(q, at_rest(q), F);
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
        m.assemble(q, at_rest(q), F);
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
    d.node_mass = (Vec(4) << 1.1, 1.3, 1.7, 1.9).finished();
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
void poison(BlockTriMat& K)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < K.num_blocks(); ++i) K.diag(i) = Mat2::Constant(nan);
    for (int i = 0; i + 1 < K.num_blocks(); ++i) {
        K.upper(i) = Mat2::Constant(nan);
        K.lower(i) = Mat2::Constant(nan);
    }
}

// Finite-difference check of the assembled Jacobian K = dF/dq on the zigzag, with gravity on.
// Each column j of the dense K is compared to the central difference of F when the flat
// coordinate j of q is perturbed. The error should shrink like h^2.
void test_assemble_jacobian()
{
    const Zigzag z = zigzag();
    const LineModel m(five_edge_line(), z.left, z.right, Params{});
    const Eigen::Index n = z.q.size();   // 8 coordinates, ordered [x0 y0 x1 y1 ...] like the dense K

    BlockTriMat K(4);
    poison(K);
    Mat2X F(2, 4);
    m.assemble(z.q, at_rest(z.q), F, &K);
    const Mat J = K.to_dense();
    CHECK(J.allFinite());

    // Forces must not depend on whether the Jacobian was requested
    Mat2X F_no_K(2, 4);
    m.assemble(z.q, at_rest(z.q), F_no_K);
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
            m.assemble(qp, at_rest(qp), Fp);
            m.assemble(qm, at_rest(qm), Fm);
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

// Node velocities for the zigzag, built from a chosen stretch rate dlen/dt per edge
// (plus a sideways part), so with damping on the edges hit every branch of edge_force:
//   edge 0: main only,     rate -0.5: damped     (T = 10,  c = 10, clips below -1)
//   edge 1: slack,         rate -0.3
//   edge 2: main + backup, rate -0.3: damped     (T = 171, c = 37, clips below -4.6)
//   edge 3: main + backup, rate -4.0: clipped    (T = 128, c = 48, clips below -2.67)
//   edge 4: main only,     rate +3.8: stretching (not chosen: fixed by the right anchor standing still)
// The node velocities are also mostly sideways to the edges, so drag (when on) is far from zero.
Mat2X zigzag_velocities(const LineModel& m, const Zigzag& z)
{
    const double rate[4] = {-0.5, -0.3, -0.3, -4.0};
    const double side[4] = {0.4, 0.7, 0.5, 0.6};
    Mat2X v(2, 4);
    Vec2 v_prev = Vec2::Zero();   // left anchor
    for (Eigen::Index e = 0; e < 4; ++e) {
        const Vec2 u = m.edge_vector(z.q, e).normalized();
        const Vec2 u_perp(-u.y(), u.x());
        v.col(e) = v_prev + rate[e] * u + side[e] * u_perp;
        v_prev = v.col(e);
    }
    return v;
}

// Finite-difference check of both assembled Jacobians, Kq = dF/dq and Kv = dF/dv, on the zigzag
// with damping on (and drag, if p has it), using zigzag_velocities.
void check_assemble_jacobians(const Params& p, const char* name)
{
    const Zigzag z = zigzag();
    const LineModel m(five_edge_line(), z.left, z.right, p);
    const Eigen::Index n = z.q.size();

    const Mat2X v = zigzag_velocities(m, z);
    const Vec2 u4 = m.edge_vector(z.q, 4).normalized();
    CHECK(m.vel_vector(v, 4).dot(u4) > 1.0);   // edge 4 really is stretching, away from the kink

    BlockTriMat Kq(4), Kv(4);
    poison(Kq);
    poison(Kv);
    Mat2X F(2, 4);
    m.assemble(z.q, v, F, &Kq, &Kv);
    const Mat Jq = Kq.to_dense();
    const Mat Jv = Kv.to_dense();
    CHECK(Jq.allFinite());
    CHECK(Jv.allFinite());

    // Forces must not depend on whether the Jacobians were requested
    Mat2X F_no_K(2, 4);
    m.assemble(z.q, v, F_no_K);
    CHECK(F_no_K == F);

    // Column j: central difference of F when flat coordinate j of q (or of v) is perturbed
    const double h = 1e-6;
    Mat Jq_fd(n, n), Jv_fd(n, n);
    for (Eigen::Index j = 0; j < n; ++j) {
        Mat2X qp = z.q, qm = z.q, vp = v, vm = v;
        qp.data()[j] += h;
        qm.data()[j] -= h;
        vp.data()[j] += h;
        vm.data()[j] -= h;
        Mat2X Fp(2, 4), Fm(2, 4);
        m.assemble(qp, v, Fp);
        m.assemble(qm, v, Fm);
        Jq_fd.col(j) = Eigen::Map<const Vec>((Fp - Fm).eval().data(), n) / (2.0 * h);
        m.assemble(z.q, vp, Fp);
        m.assemble(z.q, vm, Fm);
        Jv_fd.col(j) = Eigen::Map<const Vec>((Fp - Fm).eval().data(), n) / (2.0 * h);
    }
    const double err_q = (Jq - Jq_fd).cwiseAbs().maxCoeff();
    const double err_v = (Jv - Jv_fd).cwiseAbs().maxCoeff();
    std::printf("assembled Jacobian FD check (zigzag, %s): max err Kq %.3e, Kv %.3e\n", name, err_q, err_v);

    CHECK(err_q <= 1e-6 * std::max(1.0, Jq_fd.cwiseAbs().maxCoeff()));
    CHECK(err_v <= 1e-6 * std::max(1.0, Jv_fd.cwiseAbs().maxCoeff()));
}

void test_assemble_jacobian_damped()
{
    Params p;
    p.tau_KV = 0.1;

    // Drag off: only the edge forces depend on q and v
    p.c_drag = 0.0;
    check_assemble_jacobians(p, "damping only");

    // Drag on, with a wide webbing so the drag terms are as large as the damping terms.
    // Drag on node a depends on its own velocity and on dq = q_b - q_a, so it fills diag(a)
    // and upper(a) of Kq but only diag(a) of Kv (and likewise for node b)
    p.c_drag = 1.15;
    p.webbing_width = 1.0;
    check_assemble_jacobians(p, "damping + drag");
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
    BlockTriMat K(4);
    poison(K);
    m.assemble(q, at_rest(q), F, &K);

    CHECK(close(F.col(0), Vec2(0.0, -11.0)));
    CHECK(close(F.col(1), Vec2(0.0, -13.0)));
    CHECK(close(F.col(2), Vec2(0.0, -17.0)));
    CHECK(close(F.col(3), Vec2(0.0, -19.0)));

    // Slack everywhere: no stiffness
    CHECK(K.to_dense().isZero());
}

// Smallest valid line: two edges and one free node, so the loop runs once, the node is both
// first and last, and K has no off-diagonal blocks to write.
void test_assemble_single_node()
{
    DiscreteLine d;
    d.k_main = {100.0, 200.0};
    d.k_backup = {50.0, 50.0};
    d.l_main = {1.0, 1.0};
    d.l_backup = {1.2, 1.2};
    d.node_mass = Vec::Constant(1, 1.0);
    d.element_segment = {0, 0};

    Params p;
    p.gravity = Vec2(0.0, -10.0);
    const LineModel m(d, {0.0, 0.0}, {2.2, 0.0}, p);

    Mat2X q(2, 1);
    q << 1.1, 0.0;   // both edges stretched 0.1: forces 10 and 20 along +x
    Mat2X F(2, 1);
    BlockTriMat K(1);
    poison(K);
    m.assemble(q, at_rest(q), F, &K);

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
    BlockTriMat K(4);
    m.assemble(z.q, at_rest(z.q), F, &K);

    Vec q(8);
    for (int i = 0; i < 4; ++i) {
        q[2 * i] = z.q(0, i);
        q[2 * i + 1] = z.q(1, i);
    }

    Vec r(8);
    BlockTriMat jac(4);
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
    m.assemble(q, at_rest(q), F_bare);

    for (int i = 0; i < 4; ++i) {
        m.place_slackliner(i, mass);
        m.assemble(q, at_rest(q), F);
        for (int j = 0; j < 4; ++j) {
            const Vec2 extra = (j == i) ? Vec2(mass * p.gravity) : Vec2::Zero();
            if (!close(F.col(j), Vec2(F_bare.col(j) + extra)))
                std::cerr << "  slackliner on node " << i << ": wrong force on node " << j << "\n";
            CHECK(close(F.col(j), Vec2(F_bare.col(j) + extra)));
        }
    }

    // Removing restores the bare line, and a second remove is harmless
    m.remove_slackliner();
    m.assemble(q, at_rest(q), F);
    CHECK(F == F_bare);
    m.remove_slackliner();
    m.assemble(q, at_rest(q), F);
    CHECK(F == F_bare);
}

// Packs positions and velocities into the ODE state y = [x0 y0 x1 y1 ... | vx0 vy0 vx1 vy1 ...]
Vec pack_state(const Mat2X& q, const Mat2X& v)
{
    Vec y(q.size() + v.size());
    y << Eigen::Map<const Vec>(q.data(), q.size()), Eigen::Map<const Vec>(v.data(), v.size());
    return y;
}

// rhs gives dy = [v | F/m]. On the zigzag with gravity on, F comes from assemble and each node is
// divided by its own mass. The flat entries are compared one by one, so a wrong layout fails.
void test_rhs_layout()
{
    const Zigzag z = zigzag();
    const LineModel m(five_edge_line(), z.left, z.right, Params{});
    const double node_mass[4] = {1.1, 1.3, 1.7, 1.9};

    Mat2X v(2, 4);
    v << 0.1, -0.2, 0.3, -0.4,
         0.5, -0.6, 0.7, -0.8;
    const Vec y = pack_state(z.q, v);

    // dy starts as NaN: rhs must overwrite every entry
    Vec dy = Vec::Constant(16, std::numeric_limits<double>::quiet_NaN());
    m.rhs(0.0, y, dy);
    CHECK(dy.allFinite());

    Mat2X F(2, 4);
    m.assemble(z.q, v, F);   // same velocities as rhs, so damping matches

    for (int i = 0; i < 4; ++i) {
        CHECK(dy[2 * i] == v(0, i));
        CHECK(dy[2 * i + 1] == v(1, i));
        CHECK(std::abs(dy[8 + 2 * i] - F(0, i) / node_mass[i]) <= 1e-12 * std::abs(F(0, i)));
        CHECK(std::abs(dy[8 + 2 * i + 1] - F(1, i) / node_mass[i]) <= 1e-12 * std::abs(F(1, i)));
    }
}

// With every edge slack only gravity acts, and m * g / m = g on every node, whatever its mass.
// That holds with a slackliner too: the extra mass enters both the weight and the division.
void test_rhs_free_fall()
{
    Params p;
    p.gravity = Vec2(0.0, -10.0);
    const LineModel m(five_edge_line(), {0.0, 0.0}, {4.0, 0.0}, p);

    Mat2X q(2, 4);
    q << 0.8, 1.6, 2.4, 3.2,
         0.0, 0.0, 0.0, 0.0;   // edges 0.8 long: all slack
    const Vec y = pack_state(q, Mat2X::Zero(2, 4));
    Vec dy(16);

    auto all_nodes_fall_with_g = [&] {
        Eigen::Map<const Mat2X> dv(dy.data() + 8, 2, 4);
        for (int i = 0; i < 4; ++i)
            if (!close(Vec2(dv.col(i)), p.gravity)) return false;
        return true;
    };

    m.rhs(0.0, y, dy);
    CHECK(all_nodes_fall_with_g());

    m.place_slackliner(2, 80.0);
    m.rhs(0.0, y, dy);
    CHECK(all_nodes_fall_with_g());
    m.remove_slackliner();
}

// At the static solution the net force is zero, so a line at rest stays at rest: dy = 0.
// Checked without and with a slackliner, which moves the equilibrium.
void test_rhs_static_equilibrium()
{
    const LineModel m(five_edge_line(), {0.0, 0.0}, {5.5, 0.0}, Params{});
    const double tol = 1e-9;

    auto at_rest_after_static_solve = [&] {
        Vec q = m.static_solver_initial_guess();
        if (!m.static_solver(q, tol)) return false;

        const Mat2X Q = Eigen::Map<const Mat2X>(q.data(), 2, 4);
        Vec dy(16);
        m.rhs(0.0, pack_state(Q, Mat2X::Zero(2, 4)), dy);
        // |F| < tol and every node mass is above 1, so |F/m| < tol
        return dy.head(8).isZero() && dy.tail(8).norm() < tol;
    };

    CHECK(at_rest_after_static_solve());

    m.place_slackliner(1, 80.0);
    CHECK(at_rest_after_static_solve());
    m.remove_slackliner();
}

// Finite-difference check of the full ODE Jacobian J = d(dy)/dy that rhs returns through Kq and Kv.
// With y = [q | v] and dy = [v | F/m], J has the block structure
//     J = [ 0             | I             ]
//         [ diag(1/m) Kq  | diag(1/m) Kv  ]
// The dense J is built from the returned Kq, Kv and compared column by column to the central
// difference of rhs. Damping and drag are on and the zigzag velocities hit every edge_force branch.
// Checked without and with a slackliner, so the per-node mass division (incl. the extra mass) is tested.
void check_rhs_jacobian(const LineModel& m, const Zigzag& z, const char* name)
{
    const Eigen::Index n = z.q.size();   // 8 position coordinates, 16 state entries
    const Vec y = pack_state(z.q, zigzag_velocities(m, z));

    BlockTriMat Kq(4), Kv(4);
    poison(Kq);
    poison(Kv);
    Vec dy(2 * n);
    m.rhs(0.0, y, dy, &Kq, &Kv);
    CHECK(Kq.to_dense().allFinite());
    CHECK(Kv.to_dense().allFinite());

    // dy must not depend on whether the Jacobians were requested
    Vec dy_no_K(2 * n);
    m.rhs(0.0, y, dy_no_K);
    CHECK(dy_no_K == dy);

    // Requesting only one Jacobian gives the same blocks
    BlockTriMat Kq_only(4), Kv_only(4);
    poison(Kq_only);
    poison(Kv_only);
    m.rhs(0.0, y, dy_no_K, &Kq_only);
    m.rhs(0.0, y, dy_no_K, nullptr, &Kv_only);
    CHECK(Kq_only.to_dense() == Kq.to_dense());
    CHECK(Kv_only.to_dense() == Kv.to_dense());

    Mat J = Mat::Zero(2 * n, 2 * n);
    J.topRightCorner(n, n).setIdentity();
    J.bottomLeftCorner(n, n) = Kq.to_dense();
    J.bottomRightCorner(n, n) = Kv.to_dense();

    const double h = 1e-6;
    Mat J_fd(2 * n, 2 * n);
    Vec dy_p(2 * n), dy_m(2 * n);
    for (Eigen::Index j = 0; j < 2 * n; ++j) {
        Vec yp = y, ym = y;
        yp[j] += h;
        ym[j] -= h;
        m.rhs(0.0, yp, dy_p);
        m.rhs(0.0, ym, dy_m);
        J_fd.col(j) = (dy_p - dy_m) / (2.0 * h);
    }

    const double err = (J - J_fd).cwiseAbs().maxCoeff();
    std::printf("rhs Jacobian FD check (zigzag, %s): max err %.3e (max |J| %.3e)\n",
                name, err, J_fd.cwiseAbs().maxCoeff());
    CHECK(err <= 1e-6 * std::max(1.0, J_fd.cwiseAbs().maxCoeff()));
}

void test_rhs_jacobian()
{
    Params p;
    p.tau_KV = 0.1;
    p.c_drag = 1.15;
    p.webbing_width = 1.0;   // wide webbing: drag terms as large as the damping terms

    const Zigzag z = zigzag();
    const LineModel m(five_edge_line(), z.left, z.right, p);

    check_rhs_jacobian(m, z, "no slackliner");

    m.place_slackliner(2, 80.0);
    check_rhs_jacobian(m, z, "slackliner on node 2");
    m.remove_slackliner();
}

// Round drag parameters and no gravity, so the drag forces below can be checked by hand.
// For three_edge_line, every edge has c = 0.5 * rho * C_d * w * (l_main + l_backup) = 0.5 * 0.1 * 2.2 = 0.11.
// edge_drag returns the share of one end node, i.e. with c/2 = 0.055.
Params drag_params()
{
    Params p;
    p.gravity.setZero();
    p.rho_air = 1.0;
    p.c_drag = 1.0;
    p.webbing_width = 0.1;
    return p;
}

// edge_drag(e, dq, v) = -(c/2) |v_perp| v_perp, where v_perp is the part of v normal to the edge
void test_edge_drag()
{
    const LineModel m(three_edge_line(), {0.0, 0.0}, {3.0, 0.0}, drag_params());

    // Horizontal edge moving down at 2 m/s: (c/2) * 2 * 2 = 0.22 upwards
    CHECK(close(m.edge_drag(0, {1.0, 0.0}, {0.0, -2.0}), Vec2(0.0, 0.22)));
    // Sliding along the edge adds nothing
    CHECK(close(m.edge_drag(0, {1.0, 0.0}, {3.0, -2.0}), Vec2(0.0, 0.22)));
    CHECK(close(m.edge_drag(0, {1.0, 0.0}, {3.0, 0.0}), Vec2(0.0, 0.0)));
    // Diagonal edge, v = (2, 0): v_perp = (1, -1), |v_perp| = sqrt(2)
    const double s = 0.055 * std::sqrt(2.0);
    CHECK(close(m.edge_drag(0, {1.0, 1.0}, {2.0, 0.0}), Vec2(-s, s)));
}

// Finite-difference check of both drag Jacobians: Kq = d(drag)/d(dq) and Kv = d(drag)/dv.
// At rest and when sliding along the edge v_perp = 0, where |v_perp| v_perp has zero derivative:
// both Jacobians must come out finite and zero. There the central difference itself is only
// accurate to about c * h, hence the absolute part of the tolerance.
// Each Jacobian is also requested on its own, and must come out the same as when both are requested.
void test_edge_drag_jacobian()
{
    const LineModel m(three_edge_line(), {0.0, 0.0}, {3.0, 0.0}, drag_params());
    const Vec2 u(0.6, 0.8), u_perp(-0.8, 0.6);
    const Vec2 dq = 1.5 * u;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double h = 1e-6;

    struct Case { const char* name; Vec2 v; };
    const Case cases[] = {
        {"generic",           1.3 * u + 2.0 * u_perp},
        {"mostly sideways",  -0.4 * u - 3.0 * u_perp},
        {"mostly along edge", 3.0 * u + 0.2 * u_perp},
        {"along edge",        2.0 * u},
        {"at rest",           Vec2::Zero()},
    };

    std::cout << "edge drag Jacobian FD check\n";
    std::cout << "  case                 max err Kq   max err Kv\n";

    for (const Case& c : cases) {
        Mat2 Kq = Mat2::Constant(nan), Kv = Mat2::Constant(nan);
        m.edge_drag(0, dq, c.v, &Kq, &Kv);

        Mat2 Kq_fd, Kv_fd;
        for (int j = 0; j < 2; ++j) {
            const Vec2 d = h * Vec2::Unit(j);
            Kq_fd.col(j) = (m.edge_drag(0, dq + d, c.v) - m.edge_drag(0, dq - d, c.v)) / (2.0 * h);
            Kv_fd.col(j) = (m.edge_drag(0, dq, c.v + d) - m.edge_drag(0, dq, c.v - d)) / (2.0 * h);
        }
        const double err_q = (Kq - Kq_fd).cwiseAbs().maxCoeff();
        const double err_v = (Kv - Kv_fd).cwiseAbs().maxCoeff();
        std::printf("  %-19s  %10.3e   %10.3e\n", c.name, err_q, err_v);

        // NaN left in a Jacobian makes the error NaN, which fails these checks too
        CHECK(err_q <= 1e-6 * std::max(1.0, Kq_fd.norm()));
        CHECK(err_v <= 1e-6 * std::max(1.0, Kv_fd.norm()));

        Mat2 Kq_only = Mat2::Constant(nan), Kv_only = Mat2::Constant(nan);
        m.edge_drag(0, dq, c.v, &Kq_only);
        m.edge_drag(0, dq, c.v, nullptr, &Kv_only);
        CHECK(Kq_only == Kq);
        CHECK(Kv_only == Kv);
    }

    // Exactly at rest the Jacobians are exactly zero
    Mat2 Kq = Mat2::Constant(nan), Kv = Mat2::Constant(nan);
    m.edge_drag(0, dq, Vec2::Zero(), &Kq, &Kv);
    CHECK(Kq.isZero(0.0));
    CHECK(Kv.isZero(0.0));
}

// A straight line falling broadside. The edges are at rest length, so only gravity and drag act.
// Each node gets half the drag of both its edges: 0.11 v^2 in total. With m = 1 and g = 10 the
// forces balance at the terminal velocity v_t = sqrt(10 / 0.11).
void test_assemble_drag_terminal_velocity()
{
    Params p = drag_params();
    p.gravity = Vec2(0.0, -10.0);
    const LineModel m(three_edge_line(), {0.0, 0.0}, {3.0, 0.0}, p);

    Mat2X q(2, 2);
    q << 1.0, 2.0,
         0.0, 0.0;
    Mat2X F(2, 2);
    const double v_t = std::sqrt(10.0 / 0.11);

    auto falling_at = [](double speed) {
        Mat2X v(2, 2);
        v << 0.0, 0.0,
             -speed, -speed;
        return v;
    };

    m.assemble(q, falling_at(v_t), F);
    CHECK(close(F.col(0), Vec2(0.0, 0.0)));
    CHECK(close(F.col(1), Vec2(0.0, 0.0)));

    // Slower than v_t it still accelerates down, faster it is braked
    m.assemble(q, falling_at(0.5 * v_t), F);
    CHECK(F(1, 0) < 0.0 && F(1, 1) < 0.0);
    m.assemble(q, falling_at(2.0 * v_t), F);
    CHECK(F(1, 0) > 0.0 && F(1, 1) > 0.0);
}

// The same 3 m line split into N equal edges, all at rest length
DiscreteLine uniform_line(int N)
{
    DiscreteLine d;
    d.k_main.assign(N, 100.0);
    d.k_backup.assign(N, 50.0);
    d.l_main.assign(N, 3.0 / N);
    d.l_backup.assign(N, 3.6 / N);
    d.node_mass = Vec::Constant(N - 1, 1.0);
    d.element_segment.assign(N, 0);
    return d;
}

// The total drag on a line must not depend on the mesh. The whole line moving down at 2 m/s
// feels 0.5 * rho * C_d * w * (3.0 + 3.6) * 2^2 = 1.32 N. The outer half of each anchor edge
// pushes on the anchor, not on a free node, so the free nodes get (1 - 1/N) of it.
void test_drag_mesh_refinement()
{
    const double full_drag = 1.32;

    for (int N : {3, 6, 12, 24}) {
        const LineModel m(uniform_line(N), {0.0, 0.0}, {3.0, 0.0}, drag_params());

        Mat2X q = Mat2X::Zero(2, N - 1);
        for (int i = 0; i < N - 1; ++i) q(0, i) = 3.0 * (i + 1) / N;
        Mat2X v = Mat2X::Zero(2, N - 1);
        v.row(1).setConstant(-2.0);

        Mat2X F(2, N - 1);
        m.assemble(q, v, F);

        const double total = F.row(1).sum();
        CHECK(std::abs(total - (1.0 - 1.0 / N) * full_drag) <= 1e-12);
    }
}

} // namespace

// node_force(i) is column i of assemble(), with gravity, Kelvin-Voigt damping and drag all on.
// The line is stretched and the velocities make some edges retract, so the damping is active.
void test_node_force_matches_assemble()
{
    const LineModel m(three_edge_line(), {0.0, 0.0}, {3.6, 0.3}, Params{});

    Mat2X q(2, 2);
    q << 1.3, 2.4,
         -0.5, -0.4;
    Mat2X v(2, 2);
    v << 0.7, -1.1,
         -2.0, 0.4;

    Mat2X F(2, 2);
    m.assemble(q, v, F);
    for (Eigen::Index i = 0; i < q.cols(); ++i) {
        CHECK(close(m.node_force(q, v, i), F.col(i)));
    }
}

int main()
{
    test_edge_vector();
    test_edge_force();
    test_edge_jacobian();
    test_edge_force_damping();
    test_damping_dissipates();
    test_damping_continuous();
    test_edge_jacobian_damped();
    test_assemble_damping();
    test_assemble();
    test_assemble_jacobian();
    test_assemble_jacobian_damped();
    test_assemble_gravity();
    test_assemble_single_node();
    test_residual();
    test_both_neighboring_edges_slack();
    test_place_slackliner();
    test_rhs_layout();
    test_rhs_free_fall();
    test_rhs_static_equilibrium();
    test_rhs_jacobian();
    test_edge_drag();
    test_edge_drag_jacobian();
    test_assemble_drag_terminal_velocity();
    test_drag_mesh_refinement();
    test_node_force_matches_assemble();

    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all tests passed\n";
    return 0;
}

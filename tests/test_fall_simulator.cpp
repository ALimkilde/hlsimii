// Tests for FallSimulator::event_leash_loose, event_leash_loose_flat and free_slackliner.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.
//
// Sign convention tested: the event function is negative while the leash is taut,
// zero at the threshold, and positive once the leash would have to push (loose).

#include <cmath>
#include <iostream>

#include "fall_simulator.h"

namespace {

int failures = 0;

#define CHECK(cond)                                                               \
    do {                                                                          \
        if (!(cond)) {                                                            \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            ++failures;                                                           \
        }                                                                         \
    } while (0)

bool near_zero(double a) { return std::abs(a) <= 1e-9; }

const double gy = -9.81;  // y points up; FallSimulator takes its free fall gravity from Params{}

// Three edges between anchors (0,0) and (3,0), main rest length 1.0 and backup 1.2.
// Two free nodes of 1 kg each.
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

SlacklinerConfig slackliner()
{
    SlacklinerConfig s;
    s.m = 80.0;
    s.l_leash = 1.0;
    s.l_leg = 1.5;  // differs from l_leash so mixing them up shows
    return s;
}

// FallSimulator copies the LineModel, so the slackliner is placed before construction
FallSimulator simulator_with_slackliner_at(Eigen::Index node, const Params& p = Params{})
{
    const LineModel m(three_edge_line(), {0.0, 0.0}, {3.0, 0.0}, p);
    m.place_slackliner(static_cast<int>(node), slackliner().m);
    return FallSimulator(m, slackliner(), node, SolverType::RungeKutta45);
}

// z = [q | v] with q given column-wise as (x0, y0, x1, y1) and the nodes at rest.
// Zero velocities keep Kelvin-Voigt damping out of the force checks; reading v instead of q
// would see zeros and fail them.
Vec make_z(double x0, double y0, double x1, double y1)
{
    Vec z(8);
    z << x0, y0, x1, y1,
         0.0, 0.0, 0.0, 0.0;
    return z;
}

// Same, with nonzero junk velocities so that reading the wrong half of z shows up
Vec make_z_moving(double x0, double y0, double x1, double y1)
{
    Vec z(8);
    z << x0, y0, x1, y1,
         7.0, -3.0, 5.0, 11.0;
    return z;
}

// ====== event_leash_loose(y_acc) ======

void test_loose_threshold_is_free_fall()
{
    const FallSimulator sim = simulator_with_slackliner_at(0);

    // Node falls with exactly gravity: leash carries no force
    CHECK(near_zero(sim.event_leash_loose(gy)));
}

void test_loose_sign_when_taut()
{
    const FallSimulator sim = simulator_with_slackliner_at(0);

    CHECK(sim.event_leash_loose(0.0) < 0.0);         // hanging at rest
    CHECK(sim.event_leash_loose(20.0) < 0.0);        // line decelerating the fall hard
    CHECK(sim.event_leash_loose(gy + 0.1) < 0.0);    // falling slightly slower than gravity
}

void test_loose_sign_when_loose()
{
    const FallSimulator sim = simulator_with_slackliner_at(0);

    CHECK(sim.event_leash_loose(gy - 0.1) > 0.0);    // line pulls node down faster than gravity
    CHECK(sim.event_leash_loose(-50.0) > 0.0);       // rebound whipping the node down
}

// ====== event_leash_loose_flat(t, z) ======

// Both edges at node 0 are slack, so only gravity acts on node 0: threshold.
// Edge 2 is stretched, so node 1 has a large x-force; reading F with swapped
// indices (x of node 1 instead of y of node 0) gives a nonzero value.
void test_flat_slack_node_is_threshold_node0()
{
    const FallSimulator sim = simulator_with_slackliner_at(0);
    const Vec z = make_z(0.5, -0.1, 1.0, -0.1);
    CHECK(near_zero(sim.event_leash_loose_flat(0.0, z)));
}

// Same as above with the slackliner on node 1 (edges 1 and 2 slack, edge 0 stretched)
void test_flat_slack_node_is_threshold_node1()
{
    const FallSimulator sim = simulator_with_slackliner_at(1);
    const Vec z = make_z(2.0, -0.1, 2.5, -0.1);
    CHECK(near_zero(sim.event_leash_loose_flat(0.0, z)));
}

// Node 0 sags below the anchor: stretched edge 0 pulls it up (edge 1 is at rest length) -> taut
void test_flat_line_pulls_up_node0()
{
    const FallSimulator sim = simulator_with_slackliner_at(0);
    const Vec z = make_z(1.0, -0.5, 2.0, -0.5);
    CHECK(sim.event_leash_loose_flat(0.0, z) < 0.0);
}

// Node 0 above the anchor: stretched edge 0 pulls it down on top of gravity -> loose
void test_flat_line_pulls_down_node0()
{
    const FallSimulator sim = simulator_with_slackliner_at(0);
    const Vec z = make_z(1.0, 0.5, 2.0, 0.5);
    CHECK(sim.event_leash_loose_flat(0.0, z) > 0.0);
}

// Mirror cases for node 1, where edge 2 to the right anchor is the stretched one
void test_flat_line_pulls_up_node1()
{
    const FallSimulator sim = simulator_with_slackliner_at(1);
    const Vec z = make_z(1.0, -0.5, 2.0, -0.5);
    CHECK(sim.event_leash_loose_flat(0.0, z) < 0.0);
}

void test_flat_line_pulls_down_node1()
{
    const FallSimulator sim = simulator_with_slackliner_at(1);
    const Vec z = make_z(1.0, 0.5, 2.0, 0.5);
    CHECK(sim.event_leash_loose_flat(0.0, z) > 0.0);
}

// The event sees the node velocities through Kelvin-Voigt damping. Node 1 at (2, -0.5) as in
// test_flat_line_pulls_up_node1; edge 1 is at rest length, so only edge 2 holds node 1 up.
// Moving node 1 up makes edge 2 retract, which lowers its tension: the leash is closer to loose.
// Moving it up fast clips edge 2 to no tension, so node 1 falls freely: exactly the threshold.
// Air drag is switched off so that only the damping is tested.
void test_flat_sees_damping_node1()
{
    Params no_drag;
    no_drag.rho_air = 0.0;
    const FallSimulator sim = simulator_with_slackliner_at(1, no_drag);
    Vec z = make_z(1.0, -0.5, 2.0, -0.5);
    const double at_rest = sim.event_leash_loose_flat(0.0, z);

    z[7] = 0.1;    // vy of node 1
    const double slow = sim.event_leash_loose_flat(0.0, z);
    CHECK(at_rest < slow);
    CHECK(slow < 0.0);

    z[7] = 50.0;
    CHECK(near_zero(sim.event_leash_loose_flat(0.0, z)));
}

// The flat wrapper agrees with the core function fed the acceleration computed by hand.
// Node 0 at (1, -0.5): edge 0 has length sqrt(1.25), stretch sqrt(1.25) - 1, k = 100,
// unit vector towards the left anchor (-1, 0.5)/sqrt(1.25). Edge 1 is at rest length.
void test_flat_matches_core()
{
    const FallSimulator sim = simulator_with_slackliner_at(0);
    const Vec z = make_z(1.0, -0.5, 2.0, -0.5);

    const double len = std::sqrt(1.25);
    const double Fy_line = 100.0 * (len - 1.0) * 0.5 / len;
    const double m_total = 1.0 + slackliner().m;
    const double ay = (Fy_line + m_total * gy) / m_total;

    CHECK(std::abs(sim.event_leash_loose_flat(0.0, z) - sim.event_leash_loose(ay)) <= 1e-9);
}

// ====== free_slackliner(t, q, v) ======

// Right after release the slackliner hangs exactly one leash length below the ring,
// so the taut event starts at its threshold. t0 is nonzero so an ignored t shows up.
void test_free_starts_at_taut_threshold()
{
    FallSimulator sim = simulator_with_slackliner_at(1);
    const Vec z = make_z_moving(1.0, -0.5, 2.0, -0.7);
    const ConstStateView s = const_state_view(z);

    const double t0 = 3.0;
    sim.free_slackliner(t0, s.q, s.v);
    CHECK(near_zero(sim.event_leash_taut(t0, s.q)));
}

// The slackliner leaves with the node's y-velocity. With the ring held still, after dt
// the slackliner has moved vy*dt + 0.5*gy*dt^2, and the taut event equals minus that.
// Node 1 has vy = 11 in make_z_moving, so the slackliner is still rising: event negative.
void test_free_takes_node_velocity()
{
    FallSimulator sim = simulator_with_slackliner_at(1);
    const Vec z = make_z_moving(1.0, -0.5, 2.0, -0.7);
    const ConstStateView s = const_state_view(z);

    const double t0 = 3.0;
    const double dt = 0.1;
    const double vy = 11.0;
    sim.free_slackliner(t0, s.q, s.v);

    const double expected = -(vy * dt + 0.5 * gy * dt * dt);
    CHECK(std::abs(sim.event_leash_taut(t0 + dt, s.q) - expected) <= 1e-9);
}

// After release the node has only its own 1 kg mass. Node 0 at (1, -0.5) as in
// test_flat_matches_core: a = Fy_line / 1 + gy, so the loose event is gy - a = -Fy_line.
void test_free_removes_mass()
{
    FallSimulator sim = simulator_with_slackliner_at(0);
    const Vec z = make_z(1.0, -0.5, 2.0, -0.5);
    const ConstStateView s = const_state_view(z);

    sim.free_slackliner(0.0, s.q, s.v);

    const double len = std::sqrt(1.25);
    const double Fy_line = 100.0 * (len - 1.0) * 0.5 / len;
    CHECK(std::abs(sim.event_leash_loose_flat(0.0, z) - (-Fy_line)) <= 1e-9);
}

}  // namespace

int main()
{
    test_loose_threshold_is_free_fall();
    test_loose_sign_when_taut();
    test_loose_sign_when_loose();
    test_flat_slack_node_is_threshold_node0();
    test_flat_slack_node_is_threshold_node1();
    test_flat_line_pulls_up_node0();
    test_flat_line_pulls_down_node0();
    test_flat_line_pulls_up_node1();
    test_flat_line_pulls_down_node1();
    test_flat_sees_damping_node1();
    test_flat_matches_core();
    test_free_starts_at_taut_threshold();
    test_free_takes_node_velocity();
    test_free_removes_mass();

    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all tests passed\n";
    return 0;
}

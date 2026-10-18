// Tests for the DP5 stepper and the generic integrate() driver, through integrate_rk45().
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "rk45.h"

namespace {

using rk45::DP5;
using rk45::RKOptions;
using rk45::RKResult;
using rk45::StopReason;
using rk45::integrate_rk45;

int failures = 0;

#define CHECK(cond)                                                               \
    do {                                                                          \
        if (!(cond)) {                                                            \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            ++failures;                                                           \
        }                                                                         \
    } while (0)

const double pi = std::acos(-1.0);

// Harmonic oscillator x'' = -x as z = [x, v]. From x(0) = 1, v(0) = 0: x = cos t, v = -sin t
auto oscillator = [](double, CRef<Vec> z, Ref<Vec> dz) { dz[0] = z[1]; dz[1] = -z[0]; };
auto no_events  = [](double, CRef<Vec>, Ref<Vec>) {};
auto x_zero     = [](double, CRef<Vec> z, Ref<Vec> g) { g[0] = z[0]; };   // fires at pi/2 + k pi
auto keep_going = [](double, CRef<Vec>) { return true; };

Vec oscillator_start()
{
    Vec z(2);
    z << 1.0, 0.0;
    return z;
}

RKOptions tolerances(double tol)
{
    RKOptions opt;
    opt.rtol = tol;
    opt.atol = tol;
    return opt;
}

// Records every observer call
struct Recorder {
    std::vector<double> t;
    std::vector<Vec> z;
    bool operator()(double ti, CRef<Vec> zi) { t.push_back(ti); z.emplace_back(zi); return true; }
};

// Global error at tend follows the tolerance
void test_accuracy_follows_tolerance()
{
    for (double tol : {1e-4, 1e-6, 1e-8, 1e-10}) {
        Vec z = oscillator_start();
        const RKResult r = integrate_rk45(oscillator, no_events, 0, keep_going, keep_going, 0.0, 10.0, z, tolerances(tol));
        CHECK(r.reason == StopReason::ReachedEnd);
        CHECK(std::abs(z[0] - std::cos(10.0)) < 10.0 * tol);
        CHECK(std::abs(z[1] + std::sin(10.0)) < 10.0 * tol);
    }
}

// 5th order: 10^4 tighter tolerance needs about 10^(4/5) ~ 6.3 times more steps
void test_step_count_scales_with_order()
{
    auto count_steps = [](double tol) {
        int n = 0;
        auto on_step = [&](double, CRef<Vec>) { ++n; return true; };
        Vec z = oscillator_start();
        integrate_rk45(oscillator, no_events, 0, keep_going, on_step, 0.0, 10.0, z, tolerances(tol));
        return n;
    };
    const double ratio = static_cast<double>(count_steps(1e-10)) / count_steps(1e-6);
    CHECK(ratio > 4.0 && ratio < 9.0);
}

// Ends exactly at tend, with on_step called on increasing times up to tend
void test_reaches_tend_exactly()
{
    Recorder steps;
    Vec z = oscillator_start();
    const RKResult r = integrate_rk45(oscillator, no_events, 0, keep_going, std::ref(steps), 0.0, 3.0, z, tolerances(1e-8));
    CHECK(r.reason == StopReason::ReachedEnd);
    CHECK(r.t == 3.0);
    CHECK(r.event_index == -1);
    CHECK(steps.t.back() == 3.0);
    for (size_t i = 1; i < steps.t.size(); ++i) CHECK(steps.t[i] > steps.t[i - 1]);
}

// Output on the global grid k * dt_out in (t0, tend], interpolated to the step accuracy
void test_output_grid()
{
    RKOptions opt = tolerances(1e-8);
    opt.dt_out = 0.1;
    Recorder out;
    Vec z = oscillator_start();
    integrate_rk45(oscillator, no_events, 0, std::ref(out), keep_going, 0.0, 10.0, z, opt);

    CHECK(out.t.size() == 100);
    for (size_t i = 0; i < out.t.size(); ++i) {
        CHECK(std::abs(out.t[i] - (i + 1) * 0.1) < 1e-12);
        CHECK(std::abs(out.z[i][0] - std::cos(out.t[i])) < 1e-7);
    }
}

// Starting off the grid: first output is the next grid point, not t0 + dt_out
void test_output_grid_from_offgrid_start()
{
    RKOptions opt = tolerances(1e-8);
    opt.dt_out = 0.1;
    Recorder out;
    Vec z(2);
    z << std::cos(0.25), -std::sin(0.25);
    integrate_rk45(oscillator, no_events, 0, std::ref(out), keep_going, 0.25, 1.0, z, opt);

    CHECK(out.t.size() == 8);   // 0.3, ..., 1.0
    CHECK(std::abs(out.t.front() - 0.3) < 1e-12);
    CHECK(std::abs(out.t.back() - 1.0) < 1e-12);
}

// No output when dt_out <= 0
void test_output_disabled()
{
    Recorder out;
    Vec z = oscillator_start();
    integrate_rk45(oscillator, no_events, 0, std::ref(out), keep_going, 0.0, 1.0, z, tolerances(1e-8));
    CHECK(out.t.empty());
}

// x = cos t crosses 0 at pi/2; z is the state at the event, no output after it
void test_event_located()
{
    RKOptions opt = tolerances(1e-10);
    opt.dt_out = 0.5;
    Recorder out, steps;
    Vec z = oscillator_start();
    const RKResult r = integrate_rk45(oscillator, x_zero, 1, std::ref(out), std::ref(steps), 0.0, 10.0, z, opt);

    CHECK(r.reason == StopReason::Event);
    CHECK(r.event_index == 0);
    CHECK(std::abs(r.t - pi / 2) < 1e-8);
    CHECK(std::abs(z[0]) < 1e-12);       // Root found to machine precision on the interpolant
    CHECK(std::abs(z[1] + 1.0) < 1e-8);
    CHECK(out.t.back() == 1.5);          // 2.0 is after the event
    CHECK(steps.t.back() == r.t);
}

// After a restart at the event, g starts at ~0 on the new side and does not fire again.
// The next event is at 3 pi/2 and the output grid continues without duplicates.
void test_restart_after_event()
{
    RKOptions opt = tolerances(1e-10);
    opt.dt_out = 0.5;
    Recorder out;
    Vec z = oscillator_start();
    RKResult r = integrate_rk45(oscillator, x_zero, 1, std::ref(out), keep_going, 0.0, 10.0, z, opt);
    r = integrate_rk45(oscillator, x_zero, 1, std::ref(out), keep_going, r.t, 10.0, z, opt);

    CHECK(r.reason == StopReason::Event);
    CHECK(std::abs(r.t - 3 * pi / 2) < 1e-8);
    CHECK(out.t.size() == 9);   // 0.5, ..., 4.5
    for (size_t i = 0; i < out.t.size(); ++i) CHECK(std::abs(out.t[i] - (i + 1) * 0.5) < 1e-12);
}

// With several events the earliest one wins: x - 0.5 crosses at pi/3, before x at pi/2
void test_earliest_event_wins()
{
    auto two_events = [](double, CRef<Vec> z, Ref<Vec> g) { g[0] = z[0]; g[1] = z[0] - 0.5; };
    RKOptions opt = tolerances(1e-10);
    opt.h0 = 2.0;   // First trial step spans both roots
    Vec z = oscillator_start();
    const RKResult r = integrate_rk45(oscillator, two_events, 2, keep_going, keep_going, 0.0, 10.0, z, opt);

    CHECK(r.reason == StopReason::Event);
    CHECK(r.event_index == 1);
    CHECK(std::abs(r.t - pi / 3) < 1e-8);
}

// An observer returning false stops the integration at its time, with z at that time
void test_observer_stop()
{
    RKOptions opt = tolerances(1e-10);
    opt.dt_out = 0.1;
    auto stop_after_1 = [](double t, CRef<Vec>) { return t < 1.0 - 1e-12; };
    Vec z = oscillator_start();
    RKResult r = integrate_rk45(oscillator, no_events, 0, stop_after_1, keep_going, 0.0, 10.0, z, opt);
    CHECK(r.reason == StopReason::ObserverStop);
    CHECK(std::abs(r.t - 1.0) < 1e-12);
    CHECK(std::abs(z[0] - std::cos(1.0)) < 1e-8);

    z = oscillator_start();
    r = integrate_rk45(oscillator, no_events, 0, keep_going, stop_after_1, 0.0, 10.0, z, opt);
    CHECK(r.reason == StopReason::ObserverStop);
    CHECK(r.t >= 1.0 - 1e-12 && r.t < 10.0);
    CHECK(std::abs(z[0] - std::cos(r.t)) < 1e-8);
}

// A rhs that always gives NaN is rejected until the step size underflows
void test_step_too_small()
{
    auto nan_rhs = [](double, CRef<Vec>, Ref<Vec> dz) { dz.setConstant(std::nan("")); };
    Vec z = oscillator_start();
    const RKResult r = integrate_rk45(nan_rhs, no_events, 0, keep_going, keep_going, 0.0, 1.0, z, tolerances(1e-8));
    CHECK(r.reason == StopReason::StepTooSmall);
    CHECK(r.t == 0.0);
    CHECK(z[0] == 1.0 && z[1] == 0.0);   // State untouched
}

// Stops after max_steps attempts
void test_max_steps()
{
    RKOptions opt = tolerances(1e-10);
    opt.max_steps = 3;
    Vec z = oscillator_start();
    const RKResult r = integrate_rk45(oscillator, no_events, 0, keep_going, keep_going, 0.0, 10.0, z, opt);
    CHECK(r.reason == StopReason::MaxSteps);
    CHECK(r.t > 0.0 && r.t < 10.0);
}

// The dense output reproduces the step endpoints
void test_interpolant_endpoints()
{
    DP5 stepper(2);
    auto rhs = oscillator;
    stepper.init(rhs, 0.0, oscillator_start());
    stepper.try_step(rhs, 0.3, 1e-8, 1e-8);
    stepper.accept();

    Vec z(2);
    stepper.interpolate(stepper.t_prev(), z);
    CHECK((z - oscillator_start()).norm() < 1e-15);
    stepper.interpolate(stepper.t(), z);
    CHECK((z - stepper.state()).norm() < 1e-15);
}

} // namespace

int main()
{
    test_accuracy_follows_tolerance();
    test_step_count_scales_with_order();
    test_reaches_tend_exactly();
    test_output_grid();
    test_output_grid_from_offgrid_start();
    test_output_disabled();
    test_event_located();
    test_restart_after_event();
    test_earliest_event_wins();
    test_observer_stop();
    test_step_too_small();
    test_max_steps();
    test_interpolant_endpoints();

    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all tests passed\n";
    return 0;
}

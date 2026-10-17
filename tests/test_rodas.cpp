// Tests for the Rodas4 stepper and its integrate() driver.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "linalg.h"
#include "rodas4.h"

namespace {

using DenseRodas4 = Rodas4<RodasDenseLinSys>;

int failures = 0;

#define CHECK(cond)                                                               \
    do {                                                                          \
        if (!(cond)) {                                                            \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            ++failures;                                                           \
        }                                                                         \
    } while (0)

// Harmonic oscillator x'' = -x as z = [x, v]. From x(0) = 1, v(0) = 0: x = cos t, v = -sin t
struct Oscillator {
    int jac_calls = 0;   // One per attempted step
    void operator()(double, CRef<Vec> z, Ref<Vec> dz) const { dz[0] = z[1]; dz[1] = -z[0]; }
    void operator()(double t, CRef<Vec> z, Ref<Vec> dz, Ref<Mat> J)
    {
        ++jac_calls;
        (*this)(t, z, dz);
        J << 0.0, 1.0,
            -1.0, 0.0;
    }
};

// Stiff linear system: y0' = -y0, y1' = -k (y1 - y0) with k = 1000.
// From y(0) = [1, 0]: y0 = e^-t, y1 = k/(k-1) (e^-t - e^-kt)
struct StiffLinear {
    static constexpr double k = 1000.0;
    int jac_calls = 0;
    void operator()(double, CRef<Vec> y, Ref<Vec> dy) const { dy[0] = -y[0]; dy[1] = -k * (y[1] - y[0]); }
    void operator()(double t, CRef<Vec> y, Ref<Vec> dy, Ref<Mat> J)
    {
        ++jac_calls;
        (*this)(t, y, dy);
        J << -1.0, 0.0,
               k,  -k;
    }
    static double y1_exact(double t) { return k / (k - 1.0) * (std::exp(-t) - std::exp(-k * t)); }
};

// Always NaN, so every step is rejected
struct NanRhs {
    void operator()(double, CRef<Vec>, Ref<Vec> dz) const { dz.setConstant(std::numeric_limits<double>::quiet_NaN()); }
    void operator()(double t, CRef<Vec> z, Ref<Vec> dz, Ref<Mat> J) const { (*this)(t, z, dz); J.setZero(); }
};

Vec oscillator_start()
{
    Vec z(2);
    z << 1.0, 0.0;
    return z;
}

Rodas4Options tolerances(double tol)
{
    Rodas4Options opt;
    opt.rtol = tol;
    opt.atol = tol;
    return opt;
}

// Event function with no components: never fires
void no_events(double, CRef<Vec>, Ref<Vec>) {}

// integrate() without events, with observers that never stop
template <class Rhs>
Rodas4Result integrate_quiet(DenseRodas4& stepper, Rhs& rhs, double t0, double tend, Ref<Vec> z, const Rodas4Options& opt)
{
    auto keep_going = [](double, CRef<Vec>) { return true; };
    return integrate(stepper, rhs, no_events, 0, keep_going, keep_going, t0, tend, z, opt);
}

// Fixed steps through the stepper directly: halving h divides the error by ~2^4
void test_fixed_step_order()
{
    auto error_with_steps = [](int n) {
        Oscillator rhs;
        DenseRodas4 stepper(2);
        stepper.init(0.0, oscillator_start());
        const double h = 2.0 / n;
        for (int i = 0; i < n; ++i) {
            stepper.try_step(rhs, h, 1e-6, 1e-6);
            stepper.accept();
        }
        return std::abs(stepper.state()[0] - std::cos(2.0));
    };
    const double e1 = error_with_steps(20), e2 = error_with_steps(40), e3 = error_with_steps(80);
    const double order12 = std::log2(e1 / e2), order23 = std::log2(e2 / e3);
    std::cout << "  observed order: " << order12 << ", " << order23 << "\n";
    CHECK(order12 > 3.7 && order12 < 4.5);
    CHECK(order23 > 3.7 && order23 < 4.5);
}

// Global error at tend follows the tolerance
void test_accuracy_follows_tolerance()
{
    for (double tol : {1e-4, 1e-6, 1e-8}) {
        Oscillator rhs;
        DenseRodas4 stepper(2);
        Vec z = oscillator_start();
        const Rodas4Result r = integrate_quiet(stepper, rhs, 0.0, 10.0, z, tolerances(tol));
        std::cout << "  tol " << tol << ": error " << std::abs(z[0] - std::cos(10.0))
                  << ", attempts " << rhs.jac_calls << "\n";
        CHECK(r.reason == StopReason::ReachedEnd);
        CHECK(std::abs(z[0] - std::cos(10.0)) < 100.0 * tol);
        CHECK(std::abs(z[1] + std::sin(10.0)) < 100.0 * tol);
    }
}

// Ends exactly at tend
void test_reaches_tend_exactly()
{
    Oscillator rhs;
    DenseRodas4 stepper(2);
    Vec z = oscillator_start();
    const Rodas4Result r = integrate_quiet(stepper, rhs, 0.0, 3.0, z, tolerances(1e-8));
    CHECK(r.reason == StopReason::ReachedEnd);
    CHECK(r.t == 3.0);
    CHECK(r.event_index == -1);
}

// Stiff problem: an explicit method needs h < ~2/k, i.e. thousands of steps on [0, 10].
// Rodas4 should resolve the e^-kt transient and then take large steps.
void test_stiff_takes_few_steps()
{
    StiffLinear rhs;
    DenseRodas4 stepper(2);
    Vec y(2);
    y << 1.0, 0.0;
    const Rodas4Result r = integrate_quiet(stepper, rhs, 0.0, 10.0, y, tolerances(1e-6));
    std::cout << "  stiff: attempts " << rhs.jac_calls << ", error " << std::abs(y[1] - StiffLinear::y1_exact(10.0))
              << "\n";
    CHECK(r.reason == StopReason::ReachedEnd);
    CHECK(rhs.jac_calls < 200);
    CHECK(std::abs(y[0] - std::exp(-10.0)) < 1e-5);
    CHECK(std::abs(y[1] - StiffLinear::y1_exact(10.0)) < 1e-5);
}

// A far too large first step is rejected, and the step size still recovers afterwards
void test_recovers_after_reject()
{
    Oscillator rhs;
    DenseRodas4 stepper(2);
    Rodas4Options opt = tolerances(1e-8);
    opt.h0 = 5.0;
    Vec z = oscillator_start();
    const Rodas4Result r = integrate_quiet(stepper, rhs, 0.0, 10.0, z, opt);
    CHECK(r.reason == StopReason::ReachedEnd);
    CHECK(std::abs(z[0] - std::cos(10.0)) < 1e-6);

    // Same work as from a sensible h0, within a few attempts
    Oscillator rhs_ref;
    DenseRodas4 stepper_ref(2);
    Vec z_ref = oscillator_start();
    integrate_quiet(stepper_ref, rhs_ref, 0.0, 10.0, z_ref, tolerances(1e-8));
    std::cout << "  h0 = 5: attempts " << rhs.jac_calls << " (h0 default: " << rhs_ref.jac_calls << ")\n";
    CHECK(rhs.jac_calls < rhs_ref.jac_calls + 15);
}

// A tiny h0 must grow: without growth [0, 10] would take 10 / 1e-6 steps
void test_step_grows_from_tiny_h0()
{
    Oscillator rhs;
    DenseRodas4 stepper(2);
    Rodas4Options opt = tolerances(1e-6);
    opt.h0 = 1e-6;
    Vec z = oscillator_start();
    const Rodas4Result r = integrate_quiet(stepper, rhs, 0.0, 10.0, z, opt);
    CHECK(r.reason == StopReason::ReachedEnd);
    CHECK(rhs.jac_calls < 200);
}

// A rhs that always gives NaN is rejected until the step size underflows
void test_step_too_small()
{
    NanRhs rhs;
    DenseRodas4 stepper(2);
    Vec z = oscillator_start();
    const Rodas4Result r = integrate_quiet(stepper, rhs, 0.0, 10.0, z, tolerances(1e-6));
    CHECK(r.reason == StopReason::StepTooSmall);
    CHECK(r.t == 0.0);
    CHECK(z[0] == 1.0 && z[1] == 0.0);   // State untouched
}

// Running out of steps stops early with the state at the last accepted step
void test_max_steps()
{
    Oscillator rhs;
    DenseRodas4 stepper(2);
    Rodas4Options opt = tolerances(1e-8);
    opt.max_steps = 3;
    Vec z = oscillator_start();
    const Rodas4Result r = integrate_quiet(stepper, rhs, 0.0, 10.0, z, opt);
    CHECK(r.reason == StopReason::MaxSteps);
    CHECK(r.t > 0.0 && r.t < 10.0);
    CHECK(std::abs(z[0] - std::cos(r.t)) < 1e-6);
}

// ---------------------------------------------------------------------------
// Dense output: Rodas4::interpolate
// ---------------------------------------------------------------------------

// Before any accepted step the interpolant is the initial state
void test_interpolate_before_first_accept()
{
    DenseRodas4 stepper(2);
    stepper.init(0.0, oscillator_start());
    Vec out(2);
    stepper.interpolate(0.0, out);
    CHECK(out == oscillator_start());
}

// s = 0 and s = 1 reproduce the step's endpoints bit for bit
void test_interpolate_endpoints()
{
    Oscillator rhs;
    DenseRodas4 stepper(2);
    stepper.init(0.0, oscillator_start());
    stepper.try_step(rhs, 0.3, 1e-6, 1e-6);
    stepper.accept();
    const Vec y_start = oscillator_start(), y_end = stepper.state();

    Vec out(2);
    stepper.interpolate(0.0, out);
    CHECK(out == y_start);
    stepper.interpolate(stepper.t(), out);
    CHECK(out == y_end);
}

// One step from t = 0: the interpolant is 3rd order, so its error inside the step is O(h^4).
// Uses the max over interior points: the error at a single point can pass through zero.
void test_interpolate_order()
{
    auto max_interior_error = [](double h) {
        Oscillator rhs;
        DenseRodas4 stepper(2);
        stepper.init(0.0, oscillator_start());
        stepper.try_step(rhs, h, 1e-6, 1e-6);
        stepper.accept();
        Vec out(2);
        double err = 0.0;
        for (int i = 1; i < 10; ++i) {
            const double t = 0.1 * i * h;
            stepper.interpolate(t, out);
            err = std::max(err, std::abs(out[0] - std::cos(t)) + std::abs(out[1] + std::sin(t)));
        }
        return err;
    };
    const double e1 = max_interior_error(0.2), e2 = max_interior_error(0.1), e3 = max_interior_error(0.05);
    const double order12 = std::log2(e1 / e2), order23 = std::log2(e2 / e3);
    std::cout << "  dense output local order: " << order12 << ", " << order23 << "\n";
    CHECK(order12 > 3.7 && order12 < 4.5);
    CHECK(order23 > 3.7 && order23 < 4.5);
}

// A rejected trial step after an accept must not change the interpolant of the accepted step
void test_interpolate_after_rejected_step()
{
    Oscillator rhs;
    DenseRodas4 stepper(2);
    stepper.init(0.0, oscillator_start());
    stepper.try_step(rhs, 0.3, 1e-6, 1e-6);
    stepper.accept();

    Vec before(2), after(2);
    stepper.interpolate(0.17, before);
    stepper.try_step(rhs, 0.05, 1e-6, 1e-6);   // Trial step with a different h, never accepted
    stepper.interpolate(0.17, after);
    CHECK(before == after);
}

// ---------------------------------------------------------------------------
// Dense output: the output grid in integrate()
// ---------------------------------------------------------------------------

struct Call {
    bool is_output;   // on_output, otherwise on_step
    double t;
    Vec z;
};

// Runs integrate() and records every observer call in order
std::vector<Call> record_calls(double t0, double tend, const Rodas4Options& opt, Rodas4Result* result = nullptr)
{
    Oscillator rhs;
    DenseRodas4 stepper(2);
    Vec z(2);
    z << std::cos(t0), -std::sin(t0);
    std::vector<Call> calls;
    auto on_output = [&](double t, CRef<Vec> zt) { calls.push_back({true, t, zt}); return true; };
    auto on_step = [&](double t, CRef<Vec> zt) { calls.push_back({false, t, zt}); return true; };
    const Rodas4Result r = integrate(stepper, rhs, no_events, 0, on_output, on_step, t0, tend, z, opt);
    if (result) *result = r;
    return calls;
}

std::vector<double> output_times(const std::vector<Call>& calls)
{
    std::vector<double> ts;
    for (const Call& c : calls)
        if (c.is_output) ts.push_back(c.t);
    return ts;
}

// Output lands exactly on t_k = k * dt_out for k = 1..10 on (0, 1], and is accurate there
void test_output_on_grid()
{
    Rodas4Options opt = tolerances(1e-8);
    opt.dt_out = 0.1;
    const std::vector<Call> calls = record_calls(0.0, 1.0, opt);
    const std::vector<double> ts = output_times(calls);
    CHECK(ts.size() == 10);
    for (std::size_t i = 0; i < ts.size(); ++i) CHECK(ts[i] == static_cast<double>(i + 1) * 0.1);
    for (const Call& c : calls) {
        if (!c.is_output) continue;
        CHECK(std::abs(c.z[0] - std::cos(c.t)) < 1e-6);
        CHECK(std::abs(c.z[1] + std::sin(c.t)) < 1e-6);
    }
}

// t0 off the grid: the first output is the first grid point after t0
void test_output_t0_off_grid()
{
    Rodas4Options opt = tolerances(1e-8);
    opt.dt_out = 0.1;
    const std::vector<double> ts = output_times(record_calls(0.25, 1.0, opt));
    CHECK(!ts.empty() && ts.front() == 3 * 0.1);
    CHECK(!ts.empty() && ts.back() == 10 * 0.1);
    CHECK(ts.size() == 8);
}

// The grid point at tend is emitted even when k * dt_out rounds to just above tend
// (3 * 0.1 = 0.30000000000000004 > 0.3)
void test_output_at_tend_with_rounding()
{
    Rodas4Options opt = tolerances(1e-8);
    opt.dt_out = 0.1;
    const std::pair<int, double> cases[] = {{3, 0.3}, {6, 0.6}, {7, 0.7}, {10, 1.0}};
    for (const auto& [k_end, tend] : cases) {
        const std::vector<Call> calls = record_calls(0.0, tend, opt);
        const std::vector<double> ts = output_times(calls);
        std::cout << "  tend " << tend << ": " << ts.size() << " outputs (expected " << k_end << ")\n";
        CHECK(static_cast<int>(ts.size()) == k_end);

        // The last output is at tend (up to rounding of k * dt_out) with the state there
        const Call* last = nullptr;
        for (const Call& c : calls)
            if (c.is_output) last = &c;
        CHECK(last && std::abs(last->t - tend) < 1e-12);
        CHECK(last && std::abs(last->z[0] - std::cos(tend)) < 1e-6);
        // Output times strictly increase: no output is repeated
        for (std::size_t i = 1; i < ts.size(); ++i) CHECK(ts[i - 1] < ts[i]);
    }
}

// dt_out <= 0 disables output; a negative dt_out must not hang
void test_output_disabled()
{
    for (double dt : {0.0, -0.1}) {
        Rodas4Options opt = tolerances(1e-6);
        opt.dt_out = dt;
        Rodas4Result r{};
        const std::vector<double> ts = output_times(record_calls(0.0, 1.0, opt, &r));
        CHECK(ts.empty());
        CHECK(r.reason == StopReason::ReachedEnd);
    }
}

// All observer calls arrive in time order, and every step's outputs come before its on_step
void test_observer_calls_in_time_order()
{
    Rodas4Options opt = tolerances(1e-6);
    opt.dt_out = 0.05;   // Several outputs per step
    const std::vector<Call> calls = record_calls(0.0, 5.0, opt);
    int n_steps = 0;
    for (std::size_t i = 1; i < calls.size(); ++i) CHECK(calls[i - 1].t <= calls[i].t);
    for (const Call& c : calls) n_steps += !c.is_output;
    CHECK(n_steps > 0);
    CHECK(!calls.empty() && !calls.back().is_output);   // The final step's on_step is the last call
}

// on_output returning false stops at that output time, with z set to the output state
void test_output_observer_stops()
{
    Oscillator rhs;
    DenseRodas4 stepper(2);
    Rodas4Options opt = tolerances(1e-8);
    opt.dt_out = 0.1;
    Vec z = oscillator_start();
    Vec z_seen(2);
    auto on_output = [&](double t, CRef<Vec> zt) { z_seen = zt; return t < 0.5 - 1e-12; };
    auto on_step = [](double, CRef<Vec>) { return true; };
    const Rodas4Result r = integrate(stepper, rhs, no_events, 0, on_output, on_step, 0.0, 10.0, z, opt);
    CHECK(r.reason == StopReason::ObserverStop);
    CHECK(r.t == 5 * 0.1);
    CHECK(z == z_seen);
    CHECK(std::abs(z[0] - std::cos(r.t)) < 1e-6);
}

// on_step returning false stops at the end of that step, after its outputs
void test_step_observer_stops()
{
    Oscillator rhs;
    DenseRodas4 stepper(2);
    Rodas4Options opt = tolerances(1e-8);
    opt.dt_out = 1e-3;
    Vec z = oscillator_start();
    int n_outputs = 0;
    auto on_output = [&](double, CRef<Vec>) { ++n_outputs; return true; };
    auto on_step = [](double, CRef<Vec>) { return false; };
    const Rodas4Result r = integrate(stepper, rhs, no_events, 0, on_output, on_step, 0.0, 10.0, z, opt);
    CHECK(r.reason == StopReason::ObserverStop);
    CHECK(r.t > 0.0);
    CHECK(z == stepper.state());
    CHECK(n_outputs == static_cast<int>(std::floor(r.t / 1e-3)));   // All grid points inside the first step
}

// ---------------------------------------------------------------------------
// Events
// ---------------------------------------------------------------------------

// Runs integrate() with the given events and records every observer call
template <class Events>
std::vector<Call> run_with_events(Events&& events, int n_events, double t0, Vec z0, double tend,
                                  const Rodas4Options& opt, Rodas4Result& r, Vec& z_end)
{
    Oscillator rhs;
    DenseRodas4 stepper(2);
    z_end = z0;
    std::vector<Call> calls;
    auto on_output = [&](double t, CRef<Vec> zt) { calls.push_back({true, t, zt}); return true; };
    auto on_step = [&](double t, CRef<Vec> zt) { calls.push_back({false, t, zt}); return true; };
    r = integrate(stepper, rhs, events, n_events, on_output, on_step, t0, tend, z_end, opt);
    return calls;
}

// x = cos t falls through 0 at pi/2: stop there with the state on the new side of the root
void test_event_located()
{
    auto events = [](double, CRef<Vec> z, Ref<Vec> g) { g[0] = z[0]; };
    Rodas4Result r{};
    Vec z;
    run_with_events(events, 1, 0.0, oscillator_start(), 10.0, tolerances(1e-8), r, z);
    std::cout << "  event at t = " << r.t << " (pi/2 = " << M_PI / 2 << "), x = " << z[0] << "\n";
    CHECK(r.reason == StopReason::Event);
    CHECK(r.event_index == 0);
    CHECK(std::abs(r.t - M_PI / 2) < 1e-6);
    CHECK(z[0] <= 0.0);                        // Just past the root, so a restart does not re-fire
    CHECK(std::abs(z[1] + 1.0) < 1e-6);        // v = -sin(pi/2)
}

// An event that only depends on t is located to rounding: the bisection is not limited by the interpolant
void test_event_time_only()
{
    auto events = [](double t, CRef<Vec>, Ref<Vec> g) { g[0] = t - 0.7; };
    Rodas4Result r{};
    Vec z;
    run_with_events(events, 1, 0.0, oscillator_start(), 10.0, tolerances(1e-6), r, z);
    CHECK(r.reason == StopReason::Event);
    CHECK(r.t >= 0.7);
    CHECK(r.t - 0.7 < 1e-14);
}

// With several events, the earliest in time wins, not the lowest index
void test_event_earliest_wins()
{
    // g0 = x crosses at pi/2, g1 = x - cos(1) crosses at t = 1
    auto events = [](double, CRef<Vec> z, Ref<Vec> g) { g[0] = z[0]; g[1] = z[0] - std::cos(1.0); };
    Rodas4Options opt = tolerances(1e-8);
    opt.h0 = 0.9;   // First steps are large, so both are likely in play in the same step
    Rodas4Result r{};
    Vec z;
    run_with_events(events, 2, 0.0, oscillator_start(), 10.0, opt, r, z);
    CHECK(r.reason == StopReason::Event);
    CHECK(r.event_index == 1);
    CHECK(std::abs(r.t - 1.0) < 1e-6);
}

// g(t0) == 0 is not an event (restart convention); the next real crossing still fires
void test_event_zero_at_start_does_not_fire()
{
    // g = x - 1 starts at 0, goes negative, and returns to 0 at 2 pi
    auto events = [](double, CRef<Vec> z, Ref<Vec> g) { g[0] = z[0] - 1.0; };
    Rodas4Result r{};
    Vec z;
    run_with_events(events, 1, 0.0, oscillator_start(), 5.0, tolerances(1e-8), r, z);
    CHECK(r.reason == StopReason::ReachedEnd);
}

// Restarting from an event (as FallSimulator does) does not re-fire immediately: next stop is 3 pi / 2
void test_event_restart()
{
    auto events = [](double, CRef<Vec> z, Ref<Vec> g) { g[0] = z[0]; };
    Rodas4Result r{};
    Vec z;
    run_with_events(events, 1, 0.0, oscillator_start(), 10.0, tolerances(1e-8), r, z);
    CHECK(r.reason == StopReason::Event);

    Rodas4Result r2{};
    Vec z2;
    run_with_events(events, 1, r.t, z, 10.0, tolerances(1e-8), r2, z2);
    std::cout << "  restart: next event at t = " << r2.t << " (3 pi/2 = " << 1.5 * M_PI << ")\n";
    CHECK(r2.reason == StopReason::Event);
    CHECK(std::abs(r2.t - 1.5 * M_PI) < 1e-6);
}

// Observers see nothing past the event; the last call is on_step at the event time with the event state
void test_event_observers()
{
    auto events = [](double, CRef<Vec> z, Ref<Vec> g) { g[0] = z[0]; };
    Rodas4Options opt = tolerances(1e-8);
    opt.dt_out = 0.01;
    Rodas4Result r{};
    Vec z;
    const std::vector<Call> calls = run_with_events(events, 1, 0.0, oscillator_start(), 10.0, opt, r, z);
    CHECK(r.reason == StopReason::Event);
    for (const Call& c : calls) CHECK(c.t <= r.t);
    CHECK(!calls.empty() && !calls.back().is_output);
    CHECK(!calls.empty() && calls.back().t == r.t);
    CHECK(!calls.empty() && calls.back().z == z);
    CHECK(output_times(calls).size() == 157);   // k * 0.01 <= pi/2 = 1.5708 for k = 1..157
}

} // namespace

int main()
{
    test_fixed_step_order();
    test_accuracy_follows_tolerance();
    test_reaches_tend_exactly();
    test_stiff_takes_few_steps();
    test_recovers_after_reject();
    test_step_grows_from_tiny_h0();
    test_step_too_small();
    test_max_steps();

    test_interpolate_before_first_accept();
    test_interpolate_endpoints();
    test_interpolate_order();
    test_interpolate_after_rejected_step();

    test_output_on_grid();
    test_output_t0_off_grid();
    test_output_at_tend_with_rounding();
    test_output_disabled();
    test_observer_calls_in_time_order();
    test_output_observer_stops();
    test_step_observer_stops();

    test_event_located();
    test_event_time_only();
    test_event_earliest_wins();
    test_event_zero_at_start_does_not_fire();
    test_event_restart();
    test_event_observers();

    if (failures == 0) std::cout << "All rodas tests passed\n";
    return failures == 0 ? 0 : 1;
}

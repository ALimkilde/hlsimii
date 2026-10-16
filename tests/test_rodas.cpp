// Tests for the Rodas4 stepper and its integrate() driver.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <cmath>
#include <iostream>
#include <limits>

#include "rodas4.h"

namespace {

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

// Fixed steps through the stepper directly: halving h divides the error by ~2^4
void test_fixed_step_order()
{
    auto error_with_steps = [](int n) {
        Oscillator rhs;
        Rodas4 stepper(2);
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
        Rodas4 stepper(2);
        Vec z = oscillator_start();
        const Rodas4Result r = integrate(stepper, rhs, 0.0, 10.0, z, tolerances(tol));
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
    Rodas4 stepper(2);
    Vec z = oscillator_start();
    const Rodas4Result r = integrate(stepper, rhs, 0.0, 3.0, z, tolerances(1e-8));
    CHECK(r.reason == StopReason::ReachedEnd);
    CHECK(r.t == 3.0);
    CHECK(r.event_index == -1);
}

// Stiff problem: an explicit method needs h < ~2/k, i.e. thousands of steps on [0, 10].
// Rodas4 should resolve the e^-kt transient and then take large steps.
void test_stiff_takes_few_steps()
{
    StiffLinear rhs;
    Rodas4 stepper(2);
    Vec y(2);
    y << 1.0, 0.0;
    const Rodas4Result r = integrate(stepper, rhs, 0.0, 10.0, y, tolerances(1e-6));
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
    Rodas4 stepper(2);
    Rodas4Options opt = tolerances(1e-8);
    opt.h0 = 5.0;
    Vec z = oscillator_start();
    const Rodas4Result r = integrate(stepper, rhs, 0.0, 10.0, z, opt);
    CHECK(r.reason == StopReason::ReachedEnd);
    CHECK(std::abs(z[0] - std::cos(10.0)) < 1e-6);

    // Same work as from a sensible h0, within a few attempts
    Oscillator rhs_ref;
    Rodas4 stepper_ref(2);
    Vec z_ref = oscillator_start();
    integrate(stepper_ref, rhs_ref, 0.0, 10.0, z_ref, tolerances(1e-8));
    std::cout << "  h0 = 5: attempts " << rhs.jac_calls << " (h0 default: " << rhs_ref.jac_calls << ")\n";
    CHECK(rhs.jac_calls < rhs_ref.jac_calls + 15);
}

// A tiny h0 must grow: without growth [0, 10] would take 10 / 1e-6 steps
void test_step_grows_from_tiny_h0()
{
    Oscillator rhs;
    Rodas4 stepper(2);
    Rodas4Options opt = tolerances(1e-6);
    opt.h0 = 1e-6;
    Vec z = oscillator_start();
    const Rodas4Result r = integrate(stepper, rhs, 0.0, 10.0, z, opt);
    CHECK(r.reason == StopReason::ReachedEnd);
    CHECK(rhs.jac_calls < 200);
}

// A rhs that always gives NaN is rejected until the step size underflows
void test_step_too_small()
{
    NanRhs rhs;
    Rodas4 stepper(2);
    Vec z = oscillator_start();
    const Rodas4Result r = integrate(stepper, rhs, 0.0, 10.0, z, tolerances(1e-6));
    CHECK(r.reason == StopReason::StepTooSmall);
    CHECK(r.t == 0.0);
    CHECK(z[0] == 1.0 && z[1] == 0.0);   // State untouched
}

// Running out of steps stops early with the state at the last accepted step
void test_max_steps()
{
    Oscillator rhs;
    Rodas4 stepper(2);
    Rodas4Options opt = tolerances(1e-8);
    opt.max_steps = 3;
    Vec z = oscillator_start();
    const Rodas4Result r = integrate(stepper, rhs, 0.0, 10.0, z, opt);
    CHECK(r.reason == StopReason::MaxSteps);
    CHECK(r.t > 0.0 && r.t < 10.0);
    CHECK(std::abs(z[0] - std::cos(r.t)) < 1e-6);
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

    if (failures == 0) std::cout << "All rodas tests passed\n";
    return failures == 0 ? 0 : 1;
}

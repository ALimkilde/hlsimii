// Tests for FreeFall::get_state and set_initial_conditions.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <algorithm>
#include <cmath>
#include <iostream>

#include "free_fall.h"

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

const double gy = -9.81;  // y points up

// At t = t0 the state is exactly the initial condition
void test_state_at_t0()
{
    const FreeFall f(gy, 2.0, 10.0, 3.0);
    const FallState s = f.get_state(2.0);
    CHECK(close(s.y, 10.0));
    CHECK(close(s.vy, 3.0));
}

// Dropped from rest: after 1 s it moves down and has fallen 4.905 m
void test_drop_from_rest()
{
    const FreeFall f(gy, 0.0, 0.0, 0.0);
    const FallState s = f.get_state(1.0);
    CHECK(close(s.vy, -9.81));
    CHECK(close(s.y, -4.905));
}

// Thrown upward at a nonzero t0: dt = 1.5, vy = 3 - 9.81*1.5, y = 10 + 3*1.5 - 0.5*9.81*1.5^2
void test_nonzero_initial_conditions()
{
    const FreeFall f(gy, 2.0, 10.0, 3.0);
    const FallState s = f.get_state(3.5);
    CHECK(close(s.vy, -11.715));
    CHECK(close(s.y, 3.46375));
}

// The gravity passed in is used, not a hard-coded value
void test_uses_given_gravity()
{
    const FreeFall f(-2.0, 0.0, 0.0, 0.0);
    const FallState s = f.get_state(3.0);
    CHECK(close(s.vy, -6.0));
    CHECK(close(s.y, -9.0));
}

// After set_initial_conditions the trajectory restarts from the new state
void test_set_initial_conditions()
{
    FreeFall f(gy, 0.0, 0.0, 0.0);
    f.set_initial_conditions(1.0, 4.0, -1.0);
    const FallState s = f.get_state(2.0);
    CHECK(close(s.vy, -1.0 - 9.81));
    CHECK(close(s.y, 4.0 - 1.0 - 0.5 * 9.81));
}

} // namespace

int main()
{
    test_state_at_t0();
    test_drop_from_rest();
    test_nonzero_initial_conditions();
    test_uses_given_gravity();
    test_set_initial_conditions();

    if (failures > 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all tests passed\n";
    return 0;
}

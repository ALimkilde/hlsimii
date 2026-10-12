// Tests for ProgressBar: exact drawn text, throttling, clamping and the final newline.
// Uses a CHECK macro instead of assert so the tests also run in Release builds.

#include <cmath>
#include <iostream>
#include <sstream>

#include "progress_bar.h"

namespace {

int failures = 0;

#define CHECK(cond)                                                                    \
    do {                                                                               \
        if (!(cond)) {                                                                 \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            ++failures;                                                                \
        }                                                                              \
    } while (0)

// Half done with width 10: five '#', five spaces
void test_draw_half()
{
    std::ostringstream os;
    ProgressBar bar(os, true, 10);
    bar.update(0.5);
    CHECK(os.str() == "\r[#####     ] 50%");
}

// Updates within the same whole percent draw nothing new
void test_same_percent_not_redrawn()
{
    std::ostringstream os;
    ProgressBar bar(os, true, 10);
    bar.update(0.501);
    const std::string once = os.str();
    bar.update(0.509);
    CHECK(os.str() == once);
}

// Disabled bar writes nothing, not even the final newline
void test_disabled_writes_nothing()
{
    std::ostringstream os;
    {
        ProgressBar bar(os, false, 10);
        bar.update(0.5);
    }
    CHECK(os.str().empty());
}

// The destructor ends the line after a draw
void test_destructor_newline()
{
    std::ostringstream os;
    {
        ProgressBar bar(os, true, 10);
        bar.update(1.0);
    }
    CHECK(os.str() == "\r[##########] 100%\n");
}

// Out-of-range fractions are clamped, NaN counts as 0
void test_clamping()
{
    std::ostringstream high;
    ProgressBar bar_high(high, true, 10);
    bar_high.update(1.5);
    CHECK(high.str() == "\r[##########] 100%");

    std::ostringstream nan;
    ProgressBar bar_nan(nan, true, 10);
    bar_nan.update(std::nan(""));
    CHECK(nan.str() == "\r[          ] 0%");
}

} // namespace

int main()
{
    test_draw_half();
    test_same_percent_not_redrawn();
    test_disabled_writes_nothing();
    test_destructor_newline();
    test_clamping();

    if (failures == 0) std::cout << "all progress bar tests passed\n";
    return failures == 0 ? 0 : 1;
}

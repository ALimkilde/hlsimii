#pragma once

// One-line text progress bar, redrawn in place with '\r'.
// Redraws only when the whole percent changes; prints the final newline when destroyed.

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

#ifndef __EMSCRIPTEN__
#include <unistd.h>
#endif

class ProgressBar {
public:
    // By default only draws when stderr is a terminal, so redirected logs and tests stay clean
    explicit ProgressBar(std::ostream& os = std::cerr, bool enabled = stderr_is_terminal(), int width = 40)
        : os_(os), enabled_(enabled), width_(width)
    {
    }

    ~ProgressBar()
    {
        if (drawn_) os_ << '\n' << std::flush;
    }

    ProgressBar(const ProgressBar&) = delete;
    ProgressBar& operator=(const ProgressBar&) = delete;

    // fraction in [0, 1]; values outside are clamped and NaN counts as 0
    void update(double fraction)
    {
        if (!enabled_) return;
        if (!(fraction >= 0.0)) fraction = 0.0; // also catches NaN
        fraction = std::min(fraction, 1.0);

        const int pct = static_cast<int>(100.0 * fraction);
        if (pct == last_pct_) return;
        last_pct_ = pct;

        const int filled = pct * width_ / 100;
        os_ << "\r[" << std::string(filled, '#') << std::string(width_ - filled, ' ') << "] " << pct << '%'
            << std::flush;
        drawn_ = true;
    }

private:
    static bool stderr_is_terminal()
    {
#ifdef __EMSCRIPTEN__
        return false; // no terminal in the browser
#else
        return isatty(fileno(stderr));
#endif
    }

    std::ostream& os_;
    bool enabled_;
    int width_;
    int last_pct_ = -1;
    bool drawn_ = false;
};

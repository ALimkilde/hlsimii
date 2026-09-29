#include "line.h"

#include <stdexcept>
#include <utility>

// Webbing

double Webbing::rho() const { return weight_g_m / 1000.0; }

double Webbing::ea() const { return 100.0 * tension_kN / stretch_pct * 1e3; }

// Segment

Segment::Segment(Webbing webbing_main, Webbing webbing_backup, double l_main, double l_backup)
    : ea_main(webbing_main.ea()),
      ea_backup(webbing_backup.ea()),
      rho_main(webbing_main.rho()),
      rho_backup(webbing_backup.rho()),
      L_main(l_main),
      L_backup(l_backup)
{
}

double Segment::mass() const { return rho_main * L_main + rho_backup * L_backup; }

// Line 

Line::Line(std::vector<Segment> segments)
    : segments_(std::move(segments)) {}

int Line::num_segments() const { return segments_.size(); }
const std::vector<Segment>& Line::segments() const { return segments_; }

// Pulling

namespace {

// Shortens the main of a segment by pull, scaling the backup by the same factor
void shorten(Segment& s, double pull)
{
    if (pull >= s.L_main)
        throw std::invalid_argument("pull_webbing: pull exceeds the end segment's main length");
    const double factor = (s.L_main - pull) / s.L_main;
    s.L_main *= factor;
    s.L_backup *= factor;
}

} // namespace

Line pull_webbing(const Line& line, double pull_left, double pull_right)
{
    if (line.segments().empty())
        throw std::invalid_argument("pull_webbing: line has no segments");
    if (pull_left < 0.0 || pull_right < 0.0)
        throw std::invalid_argument("pull_webbing: pulls must be non-negative");

    std::vector<Segment> segments = line.segments();
    if (segments.size() == 1) {
        shorten(segments.front(), pull_left + pull_right);
    } else {
        shorten(segments.front(), pull_left);
        shorten(segments.back(), pull_right);
    }
    return Line(std::move(segments));
}

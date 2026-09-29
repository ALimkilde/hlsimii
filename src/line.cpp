#include "line.h"

#include <utility>

// Webbing

double Webbing::rho() const { return weight_g_m / 1000.0; }

double Webbing::ea() const { return 100.0 * tension_kN / stretch_pct * 1e3; }

// Segment

Segment::Segment(Webbing webbing_main, Webbing webbing_backup, double l_main, double l_backup)
    : ea_main(webbing_main.ea()),
      ea_backup(webbing_backup.ea()),
      L_main(l_main),
      L_backup(l_backup),
      weight(webbing_main.rho() * l_main + webbing_backup.rho() * l_backup)
{
}

// Line 

Line::Line(std::vector<Segment> segments)
    : segments_(std::move(segments)) {}

int Line::num_segments() const { return segments_.size(); }
const std::vector<Segment>& Line::segments() const { return segments_; }

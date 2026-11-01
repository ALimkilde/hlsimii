#include "summary.h"

#include <iomanip>
#include <limits>

namespace {

// Writes "key": value pairs; the caller writes the braces and the separating commas
void field(std::ostream& out, const char* key, double value) { out << "\"" << key << "\": " << value; }

void static_state(std::ostream& out, const char* height_key, double height, const StaticData& s)
{
    out << "{";
    field(out, height_key, height);
    out << ", ";
    field(out, "max_tension", s.max_tension);
    out << ", ";
    field(out, "tension_left", s.anchor_tension.left);
    out << ", ";
    field(out, "tension_right", s.anchor_tension.right);
    out << "}";
}

void fall(std::ostream& out, const SimData& d)
{
    field(out, "y_harness_start", d.y_slackliner.front());
    out << ", ";
    field(out, "y_harness_min", d.lowest_point_slackliner);
    out << ", ";
    field(out, "peak_tension_left", d.peak_force_anchors[0]);
    out << ", ";
    field(out, "peak_tension_right", d.peak_force_anchors[1]);
    out << ", ";
    field(out, "peak_tension_leash", d.peak_force_leash);
}

} // namespace

void write_summary_json(std::ostream& out, const SimConfig& cfg, const ExperimentOutput& leash,
                        const std::optional<ExperimentOutput>& backup)
{
    const Eigen::Index ring = leash.dynamics.node_slackliner;
    const auto ring_y = [ring](const StaticData& s) { return as_const_mat2x(s.q)(1, ring); };

    // Enough digits that a script reading the file gets the doubles back exactly
    out << std::setprecision(std::numeric_limits<double>::max_digits10);

    out << "{\n  ";
    field(out, "setup_mass", leash.setup_mass);
    out << ",\n  \"empty_line\": ";
    static_state(out, "sag", leash.empty_line.sag, leash.empty_line);
    out << ",\n  \"walking\": ";
    static_state(out, "y_ring", ring_y(leash.walking_line), leash.walking_line);
    out << ",\n  \"leash_fall\": {";
    fall(out, leash.dynamics);
    out << "},\n  \"backup_fall\": ";

    if (backup) {
        const StaticData& settled = backup->settled_line;
        out << "{";
        fall(out, backup->dynamics);
        out << ", \"settled\": {";
        field(out, "y_harness", ring_y(settled) - cfg.slackliner.l_leash);
        out << ", ";
        field(out, "tension_left", settled.anchor_tension.left);
        out << ", ";
        field(out, "tension_right", settled.anchor_tension.right);
        out << "}}";
    } else {
        out << "null";
    }
    out << "\n}\n";
}

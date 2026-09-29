#include "config.h"

#include <stdexcept>

void validate(const SimConfig& cfg)
{
    if (cfg.tension_kN.has_value() == cfg.pull_webbing.has_value())
        throw std::invalid_argument("config: set exactly one of tension_kN and pull_webbing");
    if (cfg.pull_side != "left" && cfg.pull_side != "right")
        throw std::invalid_argument("config: pull_side must be 'left' or 'right', got '" +
                                    cfg.pull_side + "'");
    if (cfg.L <= 0.0 || cfg.T <= 0.0)
        throw std::invalid_argument("config: L and T must be positive");
    if (cfg.N < 2)
        throw std::invalid_argument("config: N must be at least 2");
    if (cfg.segments.empty())
        throw std::invalid_argument("config: at least one segment is required");
}

namespace {

const Webbing& lookup(const WebbingCatalog& webbings, const std::string& name)
{
    auto it = webbings.find(name);
    if (it != webbings.end()) return it->second;

    std::string known;
    for (const auto& [n, w] : webbings) known += (known.empty() ? "" : ", ") + n;
    throw std::invalid_argument("config: unknown webbing '" + name + "' (known: " + known + ")");
}

} // namespace

Line build_line(const SimConfig& cfg, const WebbingCatalog& webbings)
{
    std::vector<Segment> segments;
    segments.reserve(cfg.segments.size());
    for (const SegmentConfig& s : cfg.segments)
        segments.emplace_back(lookup(webbings, s.main), lookup(webbings, s.backup), s.L_main,
                              s.L_backup);
    return Line(std::move(segments));
}

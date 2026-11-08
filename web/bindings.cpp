// embind glue between the browser and hlsim_core. The page calls simulateLeashFall and
// simulateBackupFall (in two workers, in parallel); everything here only converts their
// ExperimentOutput to JS. parseConfig and listWebbings feed the setup form.

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <exception>
#include <string>
#include <vector>

#include "config.h"
#include "config_yaml.h"
#include "experiments.h"

using emscripten::val;

namespace {

// Copies a C++ vector into a fresh JS Float64Array (the view into WASM memory is not kept)
val to_f64(const std::vector<double>& v)
{
    return val::global("Float64Array").new_(emscripten::typed_memory_view(v.size(), v.data()));
}

val error(const std::string& msg)
{
    val out = val::object();
    out.set("error", msg);
    return out;
}

val config_to_js(const SimConfig& cfg)
{
    val out = val::object();
    out.set("L", cfg.L);
    out.set("N", cfg.N);
    out.set("T", cfg.T);
    out.set("tension_kN", cfg.tension_kN ? val(*cfg.tension_kN) : val::null());
    out.set("pull_webbing", cfg.pull_webbing ? val(*cfg.pull_webbing) : val::null());
    out.set("pull_side", cfg.pull_side);
    out.set("kv_scale", cfg.kv_scale);
    out.set("drag_scale", cfg.drag_scale);
    out.set("anchor_height_left", cfg.anchor_height_left);
    out.set("anchor_height_right", cfg.anchor_height_right);

    val sl = val::object();
    sl.set("m", cfg.slackliner.m);
    sl.set("l_leg", cfg.slackliner.l_leg);
    sl.set("l_leash", cfg.slackliner.l_leash);
    sl.set("x_coor", cfg.slackliner.x_coor);
    out.set("slackliner", sl);

    val segs = val::array();
    for (const SegmentConfig& s : cfg.segments) {
        val seg = val::object();
        seg.set("main", s.main);
        seg.set("backup", s.backup);
        seg.set("L_main", s.L_main);
        seg.set("L_backup", s.L_backup);
        seg.set("break_mainline", s.break_mainline);
        segs.call<void>("push", seg);
    }
    out.set("segments", segs);
    return out;
}

// Line shape with the anchors included: left anchor, free nodes, right anchor
void push_shape(CRef<Mat2X> q, const Anchors& anchors, std::vector<double>& x, std::vector<double>& y)
{
    x.push_back(anchors.left.x());
    y.push_back(anchors.left.y());
    for (Eigen::Index i = 0; i < q.cols(); ++i) {
        x.push_back(q(0, i));
        y.push_back(q(1, i));
    }
    x.push_back(anchors.right.x());
    y.push_back(anchors.right.y());
}

// One StaticData as computed by Experiments. ring is the slackliner's node.
val static_to_js(const StaticData& s, const Anchors& anchors, Eigen::Index ring)
{
    const Eigen::Map<const Mat2X> q = as_const_mat2x(s.q);

    std::vector<double> x, y;
    push_shape(q, anchors, x, y);

    val out = val::object();
    out.set("x", to_f64(x));
    out.set("y", to_f64(y));
    out.set("sag", s.sag);
    out.set("ringX", q(0, ring));
    out.set("ringY", q(1, ring));
    out.set("anchorLeft", s.anchor_tension.left);
    out.set("anchorRight", s.anchor_tension.right);
    out.set("tensionMax", s.max_tension);
    return out;
}

// Per-frame series of one fall, for the charts and the animation. The forces are the ones
// FallSimulator::on_output recorded at rest, so the damping and drag scales don't enter here.
val dynamics_to_js(const SimData& data, const std::vector<int>& element_segment, const Anchors& anchors, const SimConfig& cfg)
{
    const double l_leash = cfg.slackliner.l_leash;

    const std::size_t frames = data.time.size();
    const Eigen::Index num_nodes = data.num_nodes;
    const Eigen::Index ring = data.node_slackliner;
    const std::size_t stride = static_cast<std::size_t>(num_nodes) + 2; // anchors included

    std::vector<double> x, y;
    x.reserve(frames * stride);
    y.reserve(frames * stride);
    std::vector<double> ring_x(frames), ring_y(frames), taut(frames);

    for (std::size_t f = 0; f < frames; ++f) {
        const Eigen::Map<const Mat2X> q = as_const_mat2x(data.frame(f));
        push_shape(q, anchors, x, y);

        ring_x[f] = q(0, ring);
        ring_y[f] = q(1, ring);
        const bool is_taut = ring_y[f] - data.y_slackliner[f] >= l_leash * (1.0 - 1e-9);
        taut[f] = is_taut ? 1.0 : 0.0;
    }

    val out = val::object();
    out.set("frames", static_cast<double>(frames));
    out.set("stride", static_cast<double>(stride));
    out.set("time", to_f64(data.time));
    out.set("x", to_f64(x));
    out.set("y", to_f64(y));
    out.set("ringX", to_f64(ring_x));
    out.set("ringY", to_f64(ring_y));
    out.set("yHarness", to_f64(data.y_slackliner));
    out.set("taut", to_f64(taut));
    out.set("tensionMax", to_f64(data.max_tension));
    out.set("anchorLeft", to_f64(data.force_anchor_left));
    out.set("anchorRight", to_f64(data.force_anchor_right));
    out.set("leash", to_f64(data.leash_tension));
    out.set("backupActive", to_f64(data.backup_activated)); // frames x edges; edge e joins frame-row points e and e+1

    // Free nodes where one rig segment meets the next, as indices into a frame row
    val junctions = val::array();
    for (std::size_t e = 1; e < element_segment.size(); ++e)
        if (element_segment[e] != element_segment[e - 1]) junctions.call<void>("push", static_cast<int>(e));
    out.set("junctions", junctions);

    // Peaks as reported by the simulator (sampled every integration step)
    val peaks = val::object();
    peaks.set("webbing", data.peak_force_webbing);
    peaks.set("leash", data.peak_force_leash);
    peaks.set("anchorLeft", data.peak_force_anchors[0]);
    peaks.set("anchorRight", data.peak_force_anchors[1]);
    peaks.set("lowestWebbing", data.lowest_point_webbing);
    peaks.set("lowestSlackliner", data.lowest_point_slackliner);
    peaks.set("stepsAccepted", data.n_accepts);
    peaks.set("stepsRejected", data.n_rejects);
    out.set("summary", peaks);
    return out;
}

// One fall as returned by Experiments. The anchors are the ones the simulation used
// (exp.anchors), so the drawn line ends where the simulated one does.
val experiment_to_js(const ExperimentOutput& exp, const SimConfig& cfg)
{
    const Eigen::Index ring = exp.dynamics.node_slackliner;
    const Anchors& anchors = exp.anchors;

    val out = val::object();
    out.set("dynamics", dynamics_to_js(exp.dynamics, exp.element_segment, anchors, cfg));
    out.set("emptyLine", static_to_js(exp.empty_line, anchors, ring));
    out.set("walkingLine", static_to_js(exp.walking_line, anchors, ring));
    out.set("settledLine", static_to_js(exp.settled_line, anchors, ring));
    out.set("setupMass", exp.setup_mass); // webbing in the span after tensioning [kg]
    return out;
}

// {config, fall}: what the page needs from one worker. fall is null for a backup
// fall when no segment breaks.
val simulate(const std::string& config_yaml, const std::string& webbings_yaml, bool backup)
{
    try {
        const SimConfig cfg = parse_config(config_yaml);
        const WebbingCatalog webbings = parse_webbings(webbings_yaml);

        val out = val::object();
        out.set("config", config_to_js(cfg));

        const RunOptions opt; // Rodas4, the same tolerances and frame rate as main.cpp
        out.set("fall", val::null());
        if (!backup) {
            const ExperimentOutput exp = Experiments::simulateLeashFall(cfg, webbings, opt);
            out.set("fall", experiment_to_js(exp, cfg));
        } else if (has_backup_fall(cfg)) {
            const ExperimentOutput exp = Experiments::simulateBackupFall(cfg, webbings, opt);
            out.set("fall", experiment_to_js(exp, cfg));
        }
        return out;
    } catch (const std::exception& e) {
        return error(e.what());
    }
}

} // namespace

val parseConfig(const std::string& yaml)
{
    try {
        return config_to_js(parse_config(yaml));
    } catch (const std::exception& e) {
        return error(e.what());
    }
}

val listWebbings(const std::string& yaml)
{
    try {
        val out = val::object();
        for (const auto& [name, w] : parse_webbings(yaml)) {
            val o = val::object();
            o.set("stretch_pct", w.stretch_pct);
            o.set("tension_kN", w.tension_kN);
            o.set("weight_g_m", w.weight_g_m);
            o.set("width_mm", w.width * 1000.0);
            out.set(name, o);
        }
        return out;
    } catch (const std::exception& e) {
        return error(e.what());
    }
}

val simulateLeashFall(const std::string& config_yaml, const std::string& webbings_yaml)
{
    return simulate(config_yaml, webbings_yaml, false);
}

val simulateBackupFall(const std::string& config_yaml, const std::string& webbings_yaml)
{
    return simulate(config_yaml, webbings_yaml, true);
}

EMSCRIPTEN_BINDINGS(hlsim)
{
    emscripten::function("parseConfig", &parseConfig);
    emscripten::function("listWebbings", &listWebbings);
    emscripten::function("simulateLeashFall", &simulateLeashFall);
    emscripten::function("simulateBackupFall", &simulateBackupFall);
}

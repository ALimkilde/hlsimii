#pragma once

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "line.h"

// Simulation input as plain data. Front ends (the YAML parser in config_yaml.h, or JS via
// embind in the WASM build) fill these structs; the simulation only ever sees them.

struct SlacklinerConfig {
    double m = 0.0;       // mass [kg]
    double l_leg = 0.0;   // length from harness to feet [m]
    double l_leash = 0.0; // length of leash [m]
    double x_coor = 0.0;  // x-coordinate of slackliner [m]
};

struct SegmentConfig {
    std::string main;   // webbing name of the main line
    std::string backup; // webbing name of the backup
    double L_main = 0.0;
    double L_backup = 0.0;
    bool break_mainline = false;
};

struct SimConfig {
    bool plots = false;
    bool gif = false;
    double L = 0.0; // length of highline spot [m]
    int N = 0;      // number of discretization vertices
    double T = 0.0; // length of simulation [s]
    std::optional<double> tension_kN;   // standing tension; exactly one of these two is set
    std::optional<double> pull_webbing; // webbing pulled out of the span [m]
    std::string pull_side = "right";    // anchor side to pull from: "left" or "right"
    SlacklinerConfig slackliner;
    std::vector<SegmentConfig> segments;
};

using WebbingCatalog = std::map<std::string, Webbing>;

// Semantic checks shared by every front end. Throws std::invalid_argument.
void validate(const SimConfig& cfg);

// Looks up the webbing names of each segment and builds the line
Line build_line(const SimConfig& cfg, const WebbingCatalog& webbings);

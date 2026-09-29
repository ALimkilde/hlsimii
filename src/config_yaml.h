#pragma once

#include <string>

#include "config.h"

// YAML front end for SimConfig. Supports only the small YAML subset used by examples/*.yaml
// and data/webbings.yaml: comments, `key: scalar`, nested block maps, block lists and
// single-line flow maps `{k: v, ...}`. Anything else throws std::invalid_argument.

WebbingCatalog parse_webbings(const std::string& yaml_text);

// Parses a simulation config and runs validate() on it
SimConfig parse_config(const std::string& yaml_text);

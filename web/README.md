# hlsimii web

Runs the hlsimii fall simulator in the browser (WebAssembly) and plots the leash and backup falls.

Compiles the C++ in `../src` unchanged; everything web-specific lives in this folder:

- `bindings.cpp`: the embind glue. The simulation interface is two functions, `simulateLeashFall` and
  `simulateBackupFall` (config YAML, webbings YAML), which call the `Experiments` functions of the same name with
  the default `RunOptions` (Rodas4, same tolerances and 1/60 s frame rate as `main.cpp`). Each returns
  `{config, fall}`; `fall` is the fall's `ExperimentOutput`, or null from `simulateBackupFall` when
  no segment breaks:
  - `setupMass`: `ExperimentOutput::setup_mass`, the webbing in the span after tensioning
  - the three `StaticData` (empty, walking and settled line). The summary table is built from these and from
    the simulator's peaks in `SimData` (as in `main.cpp`), nothing is re-solved in the glue
  - per-frame series for the charts and animation, as recorded by `FallSimulator::on_output` at rest: max
    webbing tension, anchor forces, leash force while attached, and which edges have a loaded backup
    (`SimData::backup_activated`). The nodes where rig segments connect come from
    `ExperimentOutput::element_segment`. The glue doesn't build any line model of its own

  `parseConfig` and `listWebbings` feed the setup form.
  `kv_scale` and `drag_scale` (sliders in the Damping box, log scale 0.01–10, far left = 0) are written to the
  YAML and read by the simulator's own parser.
- `CMakeLists.txt`: a separate emscripten project, independent of the native build in `../CMakeLists.txt`.
  It lists the simulator sources directly. Eigen comes from `/usr/include/eigen3`.
- `site/`: the hand-written page. Each run starts a fresh pair of module workers, each with its own WASM instance, that run
  the leash and the backup fall in parallel and are thrown away afterwards, so a crashed or hung simulator never outlives
  its run (a spare pair is warmed up in the background; Stop and a 60 s timeout terminate a run). A separate helper worker
  parses the form and the webbings and is replaced if it crashes. uPlot draws the charts, and a canvas animates the falls
  (mesh nodes as dots, segment connections as diamonds, edges with a loaded backup in orange).
  The Webbings tab edits the catalogue: a table synced with a YAML box, validated by `listWebbings`.
  The page sends that catalogue with every simulation, so `data/webbings.yaml` is only the starting point;
  edits are kept in the browser's localStorage until "reset to shipped catalogue".

## Build and run

Needs emsdk (installed in `~/emsdk`; `build.sh` sources it automatically).

```
./build.sh          # -> dist/: site/ + hlsim.mjs, hlsim.wasm + webbings.yaml, examples and docs/hlsim_numerical_model.pdf from ..
./build.sh serve    # build, then serve dist/ at http://localhost:8000
./build.sh test     # node build + headless smoke test (test/smoke.mjs) on examples/110m.yaml
```

## Publishing

`dist/` is the complete static site and the only thing to upload; it is gitignored and rebuilt by
`./build.sh`. The examples it ships are `HLSIM_EXAMPLES` in `CMakeLists.txt` and must match the
Example options in `site/index.html`. The host must serve `.wasm` as `application/wasm`.

The page must be served over http (not file://) because it uses a module worker and fetch.

## Notes

- The setup form and the YAML box stay in sync. Run always sends the YAML to the simulator's own
  parser, so you get the same validation and error messages as the native program.
- The simulator has no progress callback, so the page only shows elapsed time.
  With Rodas4, `examples/110m.yaml` takes about 0.1 s in WASM, about the same as the native release build.
- Rebuild after changing anything in `../src`, `../data/webbings.yaml`, the examples or `site/`.

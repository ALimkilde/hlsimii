// Runs the WASM simulator off the main thread. The page starts two of these and runs one fall in each.
// Messages: {id, type: 'parse'|'webbings'|'leash'|'backup', yaml, webbings}
//   parse:        yaml is a config    -> parsed config
//   webbings:     yaml is a catalogue -> {name: {stretch_pct, ...}} or {error}
//   leash/backup: yaml is a config, webbings the catalogue to simulate with
// The shipped webbings.yaml is sent with the ready message; the page owns the catalogue after that.
import createHlsim from "./hlsim.mjs";

const ready = (async () => {
  const [hlsim, webbingsYaml] = await Promise.all([
    createHlsim({ printErr: (s) => postMessage({ type: "log", text: s }) }),
    fetch("webbings.yaml").then((r) => r.text()),
  ]);
  return { hlsim, webbingsYaml };
})();

ready.then(({ hlsim, webbingsYaml }) => {
  postMessage({ type: "ready", webbingsYaml, webbings: hlsim.listWebbings(webbingsYaml) });
}, (err) => postMessage({ type: "ready", error: String(err) }));

onmessage = async ({ data: { id, type, yaml, webbings } }) => {
  const { hlsim } = await ready;
  let result;
  try {
    if (type === "parse") result = hlsim.parseConfig(yaml);
    else if (type === "webbings") result = hlsim.listWebbings(yaml);
    else if (type === "leash" || type === "backup") {
      const t0 = performance.now();
      result = type === "leash" ? hlsim.simulateLeashFall(yaml, webbings) : hlsim.simulateBackupFall(yaml, webbings);
      result.seconds = (performance.now() - t0) / 1000;
    }
  } catch (err) {
    // The bindings turn C++ exceptions into {error}, so anything thrown here is a trap or abort:
    // this WASM instance can't be trusted any more, and the page replaces the worker
    result = { error: String(err), crashed: true };
  }
  // Typed arrays are copied by structured clone; transfer them instead
  const transfer = [];
  const collect = (o) => {
    if (o instanceof Float64Array) transfer.push(o.buffer);
    else if (o && typeof o === "object") Object.values(o).forEach(collect);
  };
  collect(result);
  postMessage({ id, type: "result", result }, transfer);
};

// Headless check of the WASM build: node test/smoke.mjs [config.yaml]
// Needs the node build: ./build.sh test
import { readFileSync } from "node:fs";
import createHlsim from "../build-node/out/hlsim.mjs";

const cfgPath = process.argv[2] ?? new URL("../../examples/110m.yaml", import.meta.url);
const webbings = readFileSync(new URL("../../data/webbings.yaml", import.meta.url), "utf8");
const config = readFileSync(cfgPath, "utf8");

// The config with some top-level keys set, replacing them if the file already has them
const withKeys = (yaml, keys) =>
  yaml.split("\n").filter((l) => !Object.keys(keys).some((k) => l.startsWith(`${k}:`))).join("\n") +
  "\n" + Object.entries(keys).map(([k, v]) => `${k}: ${v}`).join("\n") + "\n";

const hlsim = await createHlsim();

// Both falls, combined the way app.js combines the two workers' results
const run = (yaml, web) => {
  const leash = hlsim.simulateLeashFall(yaml, web), backup = hlsim.simulateBackupFall(yaml, web);
  const failed = [leash, backup].find((x) => x.error);
  if (failed) return failed;
  return { config: leash.config, setupMass: leash.fall.setupMass, leash: leash.fall, backup: backup.fall };
};

const t0 = performance.now();
const r = run(config, webbings);
const secs = (performance.now() - t0) / 1000;
if (r.error) {
  console.error("error:", r.error);
  process.exit(1);
}

const kN = (n) => (n / 1e3).toFixed(5);
const max = (a) => a.reduce((m, x) => Math.max(m, x), -Infinity);
const staticLine = (label, s) =>
  console.log(`  ${label}: anchors ${kN(s.anchorLeft)} / ${kN(s.anchorRight)} kN, sag ${s.sag.toFixed(3)} m, ring (${s.ringX.toFixed(2)}, ${s.ringY.toFixed(3)})`);

console.log(`scales: kv ${r.config.kv_scale}, drag ${r.config.drag_scale}`);
// Setup mass comes from ExperimentOutput; both falls run on the same rig, so it must agree
console.log(`setup mass ${r.setupMass.toFixed(3)} kg`);
if (!(r.setupMass > 0) || (r.backup && r.backup.setupMass !== r.leash.setupMass)) {
  console.error("setup mass missing, or differs between the leash and backup fall");
  process.exit(1);
}
for (const [label, exp] of [["Leash fall", r.leash], ["Backup fall", r.backup]]) {
  if (!exp) continue;
  const fall = exp.dynamics, s = fall.summary;
  console.log(`${label}: ${fall.frames} frames, t_end = ${fall.time[fall.frames - 1]}`);
  staticLine("empty line  ", exp.emptyLine);
  staticLine("walking line", exp.walkingLine);
  staticLine("settled line", exp.settledLine);
  console.log(`  peak webbing tension ${kN(s.webbing)} kN (frames: ${kN(max(fall.tensionMax))})`);
  console.log(`  peak anchor forces   ${kN(s.anchorLeft)} / ${kN(s.anchorRight)} kN`);
  console.log(`  peak leash force     ${kN(s.leash)} kN (frames: ${kN(max(fall.leash))})`);
  console.log(`  lowest webbing ${s.lowestWebbing.toFixed(5)} m, lowest slackliner ${s.lowestSlackliner.toFixed(5)} m`);
  console.log(`  solver steps ${s.stepsAccepted} accepted, ${s.stepsRejected} rejected; junction nodes [${fall.junctions}]`);
  const edges = fall.stride - 1;
  console.log(`  edges with backup loaded, last frame: ${fall.backupActive.slice((fall.frames - 1) * edges, fall.frames * edges).reduce((a, b) => a + b, 0)} of ${edges}`);
}
console.log(`wasm run took ${secs.toFixed(1)} s`);

// The damping scales must reach the dynamics: no damping and no drag gives a different leash peak
const undamped = run(withKeys(config, { kv_scale: 0, drag_scale: 0 }), webbings);
if (undamped.error) {
  console.error("undamped run error:", undamped.error);
  process.exit(1);
}
const peak = (res) => res.leash.dynamics.summary.leash;
console.log(`leash peak, scales 0/0: ${kN(peak(undamped))} kN vs ${kN(peak(r))} kN`);
if (peak(undamped) === peak(r)) {
  console.error("kv_scale/drag_scale had no effect on the leash fall");
  process.exit(1);
}

// Anchor heights must reach the simulation: the drawn line starts and ends at the anchors
const sloped = run(withKeys(config, { anchor_height_left: 3, anchor_height_right: -2 }), webbings);
if (sloped.error) {
  console.error("sloped run error:", sloped.error);
  process.exit(1);
}
const shape = sloped.leash.emptyLine, last = shape.y.length - 1;
console.log(`anchors 3 / -2 m: line ends at y = ${shape.y[0]} / ${shape.y[last]} m, sag ${shape.sag.toFixed(3)} m`);
const frameEnds = sloped.leash.dynamics.y;
if (shape.y[0] !== 3 || shape.y[last] !== -2 || frameEnds[0] !== 3 || frameEnds[sloped.leash.dynamics.stride - 1] !== -2) {
  console.error("anchor heights did not reach the line shape");
  process.exit(1);
}

// The page edits the catalogue: listWebbings reports width, and an edited catalogue reaches the dynamics
const catalogue = hlsim.listWebbings(webbings);
console.log(`catalogue: ${Object.keys(catalogue).length} webbings, solid width ${catalogue.solid?.width_mm} mm`);
if (catalogue.solid?.width_mm !== 20 || catalogue.joker?.width_mm !== 25) {
  console.error("listWebbings did not report width_mm");
  process.exit(1);
}
const stretchier = run(config, webbings.replace(/^joker:.*$/m, "joker: {stretch_pct: 7.2, tension_kN: 5, weight_g_m: 54}"));
if (stretchier.error) {
  console.error("edited catalogue run error:", stretchier.error);
  process.exit(1);
}
console.log(`joker stretch 3.6 -> 7.2 %: leash peak ${kN(peak(stretchier))} kN vs ${kN(peak(r))} kN`);
if (peak(stretchier) === peak(r)) {
  console.error("the edited catalogue had no effect on the leash fall");
  process.exit(1);
}
console.log("bad catalogue ->", hlsim.listWebbings("joker: {stretch_pct: 3.6}").error);

// Error paths
const bad = run(config.replace(/main: joker/, "main: nosuchwebbing"), webbings);
console.log("bad webbing ->", bad.error);
const negative = run(withKeys(config, { kv_scale: -1 }), webbings);
console.log("negative kv_scale ->", negative.error);

// No mainline break: the backup worker returns fall: null instead of simulating
const unbroken = hlsim.simulateBackupFall(config.replaceAll("break_mainline: true", "break_mainline: false"), webbings);
console.log("no break -> backup fall", unbroken.error ?? unbroken.fall);
if (unbroken.error || unbroken.fall !== null) {
  console.error("expected a null backup fall without a mainline break");
  process.exit(1);
}

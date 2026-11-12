// Front end: setup form <-> YAML, runs the simulator in a worker, plots and animates the falls.

const $ = (sel, root = document) => root.querySelector(sel);
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

// ---------------------------------------------------------------- worker

// Two workers, each with its own copy of the WASM module, so the leash and backup falls run in
// parallel. The first one also answers the form's parse requests.
function startWorker() {
  const worker = new Worker("worker.js", { type: "module" });
  const pending = new Map();
  let nextId = 0;
  const ready = new Promise((resolve) => {
    worker.onmessage = ({ data }) => {
      if (data.type === "ready") resolve(data);
      else if (data.type === "result") {
        pending.get(data.id)(data.result);
        pending.delete(data.id);
      } else if (data.type === "log") console.log("[hlsim]", data.text);
    };
  });
  const call = (type, yaml, webbings) =>
    new Promise((resolve) => {
      const id = nextId++;
      pending.set(id, resolve);
      worker.postMessage({ id, type, yaml, webbings });
    });
  return { ready, call };
}

const workers = [startWorker(), startWorker()];
const call = workers[0].call;
const workerReady = Promise.all(workers.map((w) => w.ready)).then((r) => r.find((x) => x.error) ?? r[0]);
let webbingNames = [];
let webbingsYaml = ""; // the catalogue the simulations use, owned by the Webbings tab

// Both falls at once; the result has the shape showResults expects
async function simulate(yaml) {
  const [leash, backup] = await Promise.all([
    workers[0].call("leash", yaml, webbingsYaml),
    workers[1].call("backup", yaml, webbingsYaml),
  ]);
  const failed = [leash, backup].find((r) => r.error);
  if (failed) return failed;
  return { config: leash.config, setupMass: leash.fall.setupMass, leash: leash.fall, backup: backup.fall };
}

// ---------------------------------------------------------------- form <-> YAML

const form = $("#form");
const yamlBox = $("#yaml");
const errorBox = $("#error");
const runButton = $("#run");

function showError(msg) {
  errorBox.hidden = !msg;
  errorBox.textContent = msg || "";
}

// Options are the catalogue plus the current value, so a value no longer in the catalogue stays
// visible and the simulator reports the unknown webbing instead of the form silently switching it
function setWebbingOptions(sel, value) {
  sel.replaceChildren();
  for (const name of new Set([...webbingNames, value].filter(Boolean))) sel.add(new Option(name, name));
  sel.value = value;
}

function webbingSelect(value) {
  const sel = document.createElement("select");
  setWebbingOptions(sel, value);
  return sel;
}

function addSegmentRow(seg) {
  const tr = document.createElement("tr");
  const cell = (el) => { tr.insertCell().append(el); return el; };
  const num = (v) => Object.assign(document.createElement("input"), { type: "number", step: "any", min: "0", value: v });

  const main = cell(webbingSelect(seg.main));
  const backup = cell(webbingSelect(seg.backup));
  const lMain = cell(num(seg.L_main));
  const lBackup = cell(num(seg.L_backup));
  const brk = cell(Object.assign(document.createElement("input"), { type: "checkbox", checked: seg.break_mainline }));
  const del = Object.assign(document.createElement("button"), { type: "button", className: "secondary", textContent: "×", title: "remove segment" });
  del.onclick = () => { tr.remove(); formChanged(); };
  cell(del);

  tr.read = () => ({ main: main.value, backup: backup.value, L_main: lMain.value, L_backup: lBackup.value, break_mainline: brk.checked });
  $("#segments tbody").append(tr);
}

// Damping/drag scale sliders: position 0 is off, positions 1..300 span 10^-2 .. 10^1 (1 at position 200).
// The exact value lives in dataset.value, so a YAML value like 2.5 survives a round trip through the slider.
const SCALES = ["kv_scale", "drag_scale"];
const sliderToScale = (s) => (s <= 0 ? 0 : 10 ** ((s - 200) / 100));
const scaleToSlider = (v) => (v <= 0 ? 0 : Math.min(300, Math.max(1, Math.round(200 + 100 * Math.log10(v)))));
const formatScale = (v) => (v === 0 ? "0 (off)" : `${+v.toPrecision(v < 1 ? 2 : 3)}×`);

function setScale(name, value) {
  const slider = form.elements[name];
  slider.dataset.value = value;
  slider.value = scaleToSlider(value);
  $(`#${name}-value`).textContent = formatScale(value);
  $(`button[data-reset="${name}"]`).disabled = value === 1;
}

function fillForm(cfg) {
  const f = form.elements;
  for (const k of ["L", "N", "T", "pull_side", "anchor_height_left", "anchor_height_right"]) f[k].value = cfg[k];
  f.pull_webbing.value = cfg.pull_webbing ?? "";
  f.tension_kN.value = cfg.tension_kN ?? "";
  for (const k of ["m", "l_leg", "l_leash", "x_coor"]) f[k].value = cfg.slackliner[k];
  for (const k of SCALES) setScale(k, cfg[k]);
  $("#segments tbody").replaceChildren();
  cfg.segments.forEach(addSegmentRow);
}

// Writes the subset of YAML that hlsimii's config_yaml.cpp reads
function formToYaml() {
  const f = form.elements;
  const v = (name) => f[name].value.trim();
  const lines = [
    "# Generated from the form. Webbing names refer to webbings.yaml.",
    `L: ${v("L")}              # length of highline spot [m]`,
    `N: ${v("N")}              # number of discretization vertices`,
    `T: ${v("T")}              # length of simulation [s]`,
    `anchor_height_left: ${v("anchor_height_left") || 0}     # [m]`,
    `anchor_height_right: ${v("anchor_height_right") || 0}    # [m]`,
  ];
  // The config takes exactly one of the two; an empty form falls back to pull_webbing so the parser reports it
  if (v("tension_kN") !== "") lines.push(`tension_kN: ${v("tension_kN")}       # standing tension [kN]`);
  else lines.push(`pull_webbing: ${v("pull_webbing")}       # webbing pulled [m]`);
  lines.push(`pull_side: ${v("pull_side")}`);
  lines.push(`kv_scale: ${+(+(f.kv_scale.dataset.value ?? 1)).toPrecision(3)}        # multiplies tau_KV`);
  lines.push(`drag_scale: ${+(+(f.drag_scale.dataset.value ?? 1)).toPrecision(3)}      # multiplies c_drag`);
  lines.push("", "slackliner:");
  lines.push(`  m: ${v("m")}             # mass [kg]`, `  l_leg: ${v("l_leg")}`, `  l_leash: ${v("l_leash")}`, `  x_coor: ${v("x_coor")}`);
  lines.push("", "segments:");
  for (const tr of $("#segments tbody").rows) {
    const s = tr.read();
    lines.push(`  - {main: ${s.main}, backup: ${s.backup}, L_main: ${s.L_main}, L_backup: ${s.L_backup}, break_mainline: ${s.break_mainline}}`);
  }
  return lines.join("\n") + "\n";
}

const debounce = (fn, ms) => {
  let h;
  return (...a) => { clearTimeout(h); h = setTimeout(() => fn(...a), ms); };
};

// Validate with the simulator's own parser, so errors match the native program
const validateYaml = debounce(async (fill) => {
  const cfg = await call("parse", yamlBox.value);
  if (cfg.error) return showError(cfg.error);
  showError("");
  if (fill) fillForm(cfg);
}, 300);

function formChanged(event) {
  // Typing in one of tension_kN / pull_webbing clears the other
  const f = form.elements;
  if (event?.target === f.pull_webbing && f.pull_webbing.value.trim() !== "") f.tension_kN.value = "";
  if (event?.target === f.tension_kN && f.tension_kN.value.trim() !== "") f.pull_webbing.value = "";
  for (const k of SCALES) if (event?.target === f[k]) setScale(k, sliderToScale(+f[k].value));
  yamlBox.value = formToYaml();
  validateYaml(false);
}

// Reset buttons put a scale back to the simulator's default of 1
for (const button of form.querySelectorAll("button[data-reset]"))
  button.onclick = () => { setScale(button.dataset.reset, 1); formChanged(); };

form.addEventListener("input", formChanged);
form.addEventListener("change", formChanged);
yamlBox.addEventListener("input", () => validateYaml(true));
$("#add-segment").onclick = () => {
  const rows = $("#segments tbody").rows;
  const last = rows.length ? rows[rows.length - 1].read() : { main: webbingNames[0], backup: webbingNames[0], L_main: 10, L_backup: 11, break_mainline: false };
  addSegmentRow({ ...last, break_mainline: false });
  formChanged();
};

// Per-viewer storage; can be unavailable (private windows, blocked site data)
const store = {
  get: (k) => { try { return localStorage.getItem(k); } catch { return null; } },
  set: (k, v) => { try { localStorage.setItem(k, v); return true; } catch { return false; } },
  remove: (k) => { try { localStorage.removeItem(k); } catch {} },
};

// ---------------------------------------------------------------- setups (examples + saved)

// Saved setups are the config YAML by name, kept in this browser: { name: yaml }
const SETUPS_KEY = "hlsim.setups";
const SAVED_PREFIX = "saved:"; // option values of saved setups; examples are file names
const setupSelect = $("#example");
const savedGroup = $("#saved-setups");

function savedSetups() {
  try { return JSON.parse(store.get(SETUPS_KEY)) || {}; } catch { return {}; }
}

function fillSavedSetups(selected) {
  const setups = savedSetups();
  savedGroup.replaceChildren(...Object.keys(setups).sort((a, b) => a.localeCompare(b)).map((name) => new Option(name, SAVED_PREFIX + name)));
  savedGroup.hidden = savedGroup.children.length === 0;
  if (selected !== undefined) setupSelect.value = SAVED_PREFIX + selected;
  $("#delete-setup").hidden = !setupSelect.value.startsWith(SAVED_PREFIX);
}

async function loadSetup(value) {
  $("#delete-setup").hidden = !value.startsWith(SAVED_PREFIX);
  yamlBox.value = value.startsWith(SAVED_PREFIX)
    ? savedSetups()[value.slice(SAVED_PREFIX.length)] ?? ""
    : await fetch(value).then((r) => r.text());
  const cfg = await call("parse", yamlBox.value);
  if (cfg.error) return showError(cfg.error);
  showError("");
  fillForm(cfg);
}
setupSelect.onchange = (e) => loadSetup(e.target.value);

$("#save-setup").onclick = () => {
  const current = setupSelect.value.startsWith(SAVED_PREFIX) ? setupSelect.value.slice(SAVED_PREFIX.length) : "";
  const name = prompt("Save the current setup as:", current)?.trim();
  if (!name) return;
  const setups = savedSetups();
  if (name in setups && name !== current && !confirm(`Replace the saved setup "${name}"?`)) return;
  setups[name] = yamlBox.value;
  if (!store.set(SETUPS_KEY, JSON.stringify(setups)))
    return alert("Couldn't save: this browser blocks site storage (private window or blocked site data).");
  fillSavedSetups(name);
};

$("#delete-setup").onclick = () => {
  const name = setupSelect.value.slice(SAVED_PREFIX.length);
  if (!confirm(`Delete the saved setup "${name}"?`)) return;
  const setups = savedSetups();
  delete setups[name];
  store.set(SETUPS_KEY, JSON.stringify(setups));
  // Back to the first example, so the dropdown names what the form shows
  setupSelect.selectedIndex = 0;
  fillSavedSetups();
  loadSetup(setupSelect.value);
};

// ---------------------------------------------------------------- webbings tab
const WEBBINGS_KEY = "hlsim.webbings";
const webbingsBox = $("#webbings-yaml");
const webbingsError = $("#webbings-error");
let shippedWebbingsYaml = "";

function showWebbingsError(msg) {
  webbingsError.hidden = !msg;
  webbingsError.textContent = msg || "";
}

function addWebbingRow(name, w = {}) {
  const tr = document.createElement("tr");
  const cell = (el) => { tr.insertCell().append(el); return el; };
  const input = (props) => cell(Object.assign(document.createElement("input"), props));
  const num = (v, placeholder = "") => input({ type: "number", step: "any", min: "0", value: v ?? "", placeholder });

  const nameIn = input({ type: "text", className: "name", value: name, spellcheck: false });
  const stretch = num(w.stretch_pct);
  const tension = num(w.tension_kN);
  const weight = num(w.weight_g_m);
  const width = num(w.width_mm === 25 ? "" : w.width_mm, "25");
  const del = cell(Object.assign(document.createElement("button"), { type: "button", className: "secondary", textContent: "×", title: "remove webbing" }));
  del.onclick = () => { tr.remove(); webbingsTableChanged(); };

  // Renaming a webbing renames it in the segments that use it
  tr.dataset.name = name;
  nameIn.addEventListener("change", () => {
    const from = tr.dataset.name, to = nameIn.value.trim();
    tr.dataset.name = to;
    if (from && to && from !== to) renameInSegments(from, to);
  });

  tr.read = () => ({ name: nameIn.value.trim(), stretch_pct: stretch.value.trim(), tension_kN: tension.value.trim(), weight_g_m: weight.value.trim(), width_mm: width.value.trim() });
  $("#webbings tbody").append(tr);
  return nameIn;
}

function fillWebbingsTable(webbings) {
  $("#webbings tbody").replaceChildren();
  for (const [name, w] of Object.entries(webbings)) addWebbingRow(name, w);
}

// The same format as the shipped data/webbings.yaml
function webbingsTableToYaml() {
  const lines = [
    "# Webbing catalogue, edited in the web page.",
    "#   stretch_pct  stretch at tension_kN [%]",
    "#   tension_kN   tension the stretch is measured at [kN]",
    "#   weight_g_m   weight [g/m]",
    "#   width_mm     webbing width [mm], optional (default 25), only used for air drag",
  ];
  for (const tr of $("#webbings tbody").rows) {
    const w = tr.read();
    const width = w.width_mm !== "" ? `, width_mm: ${w.width_mm}` : "";
    lines.push(`${w.name}: {stretch_pct: ${w.stretch_pct}, tension_kN: ${w.tension_kN}, weight_g_m: ${w.weight_g_m}${width}}`);
  }
  return lines.join("\n") + "\n";
}

function refreshSegmentSelects() {
  for (const sel of $("#segments tbody").querySelectorAll("select")) setWebbingOptions(sel, sel.value);
}

function renameInSegments(from, to) {
  for (const sel of $("#segments tbody").querySelectorAll("select"))
    if (sel.value === from) setWebbingOptions(sel, to);
  formChanged();
}

// Validates a catalogue with the simulator's own reader and makes it the active one
async function setWebbings(yaml, { fillTable }) {
  const webbings = await call("webbings", yaml);
  if (webbings.error) return showWebbingsError(webbings.error);
  showWebbingsError("");
  webbingsYaml = yaml;
  webbingNames = Object.keys(webbings);
  if (yaml === shippedWebbingsYaml) store.remove(WEBBINGS_KEY);
  else store.set(WEBBINGS_KEY, yaml);
  if (fillTable) fillWebbingsTable(webbings);
  refreshSegmentSelects();
}

const webbingsTableChanged = debounce(() => {
  webbingsBox.value = webbingsTableToYaml();
  // The parser would only say "missing value" on some line; name the cell instead
  for (const tr of $("#webbings tbody").rows) {
    const w = tr.read();
    const empty = ["name", "stretch_pct", "tension_kN", "weight_g_m"].find((k) => w[k] === "");
    if (empty) return showWebbingsError(`webbing '${w.name || "(unnamed)"}': ${empty} is empty`);
  }
  setWebbings(webbingsBox.value, { fillTable: false });
}, 300);

$("#webbings").addEventListener("input", webbingsTableChanged);
webbingsBox.addEventListener("input", debounce(() => setWebbings(webbingsBox.value, { fillTable: true }), 300));
$("#add-webbing").onclick = () => {
  const taken = new Set(webbingNames);
  let name = "new_webbing";
  for (let i = 2; taken.has(name); i++) name = `new_webbing_${i}`;
  addWebbingRow(name, { stretch_pct: 3, tension_kN: 5, weight_g_m: 50 }).select();
  webbingsTableChanged();
};
$("#reset-webbings").onclick = async () => {
  webbingsBox.value = shippedWebbingsYaml;
  await setWebbings(shippedWebbingsYaml, { fillTable: true });
};

// ---------------------------------------------------------------- tabs

const TAB_KEY = "hlsim.tab";
function selectTab(id) {
  for (const tab of document.querySelectorAll('[role="tab"]')) {
    const on = tab.id === id;
    tab.setAttribute("aria-selected", on);
    $(`#${tab.getAttribute("aria-controls")}`).hidden = !on;
  }
  store.set(TAB_KEY, id);
}
for (const tab of document.querySelectorAll('[role="tab"]')) tab.onclick = () => selectTab(tab.id);
if (store.get(TAB_KEY) === "tab-webbings") selectTab("tab-webbings");

// ---------------------------------------------------------------- run

runButton.onclick = async () => {
  runButton.disabled = true;
  runButton.textContent = "Simulating…";
  showError("");
  const t0 = performance.now();
  const timer = setInterval(() => ($("#status").textContent = `${((performance.now() - t0) / 1000).toFixed(0)} s`), 250);
  try {
    const result = await simulate(yamlBox.value);
    if (result.error) showError(result.error);
    else {
      showResults(result);
      $("#status").textContent = `done in ${((performance.now() - t0) / 1000).toFixed(1)} s`;
    }
  } finally {
    clearInterval(timer);
    if (errorBox.hidden === false) $("#status").textContent = "";
    runButton.disabled = false;
    runButton.textContent = "Run simulation";
  }
};

// Ctrl+R (Cmd+R on macOS) runs the simulation instead of reloading the page
document.addEventListener("keydown", (e) => {
  if ((e.ctrlKey || e.metaKey) && !e.shiftKey && !e.altKey && e.key.toLowerCase() === "r") {
    e.preventDefault();
    if (!runButton.disabled) runButton.click();
  }
});

// ---------------------------------------------------------------- results

const view = { result: null, falls: [], frame: 0, playing: false };
const kN = (n) => (n / 1000).toFixed(2);
const max = (a) => a.reduce((m, x) => Math.max(m, x), -Infinity);
const min = (a) => a.reduce((m, x) => Math.min(m, x), Infinity);

function showResults(result) {
  view.result = result;
  view.frame = 0;
  // Each fall is an ExperimentOutput: per-frame dynamics and three StaticData.
  // The charts and the animation draw the dynamics.
  const falls = [["Leash fall", result.leash], ["Backup fall", result.backup]]
    .filter(([, exp]) => exp)
    .map(([label, exp]) => [label, exp.dynamics]);

  // Summary table, same rows and columns as the ISA backup fall calculator.
  // Like ISA, the height is the leash ring for "Walking" and the harness for the falls
  // (harness = ring - l_leash once the leash is taut).
  const m = (v) => v.toFixed(2);
  const l_leash = result.config.slackliner.l_leash;
  const staticRow = (label, height, s) => [label, m(height), kN(s.anchorLeft), kN(s.anchorRight), "-"];
  // A fall's row: the simulator's peaks
  const fallRow = (label, exp) => {
    const S = exp.dynamics.summary;
    return [label, m(S.lowestSlackliner), kN(S.anchorLeft), kN(S.anchorRight), kN(S.leash)];
  };
  const rows = [];
  const L = result.leash;
  rows.push(staticRow("Standing, no one on line (sag)", L.emptyLine.sag, L.emptyLine));
  rows.push(staticRow("Walking (ring)", L.walkingLine.ringY, L.walkingLine));
  rows.push(fallRow("Leash fall (harness)", L));
  if (result.backup) {
    const B = result.backup, S = B.settledLine;
    rows.push(fallRow("Backup fall - impact (harness)", B));
    rows.push(staticRow("Backup fall - settled (harness)", S.ringY - l_leash, S));
  } else {
    rows.push(["Backup fall", "no mainline break in this setup", "", "", ""]);
  }
  $("#summary tbody").replaceChildren(
    ...rows.map((cells) => {
      const tr = document.createElement("tr");
      for (const text of cells) tr.insertCell().textContent = text;
      return tr;
    }),
  );
  $("#setup-weight").textContent = `The setup weights ${result.setupMass.toFixed(2)} kg.`;

  // Shared axis ranges so the two falls compare directly
  const all = falls.map(([, f]) => f);
  const forceMax = max(all.map((f) => Math.max(max(f.tensionMax), max(f.leash)))) / 1000;
  const yLo = min(all.map((f) => Math.min(min(f.yHarness), min(f.y))));
  const yHi = max(all.map((f) => Math.max(max(f.yHarness), max(f.y))));
  const ranges = { force: [0, forceMax * 1.08], yLo, yHi };

  for (const f of view.falls) f.forces.destroy();
  const container = $("#falls");
  container.replaceChildren();
  view.falls = falls.map(([label, fall]) => {
    const el = $("#fall-template").content.firstElementChild.cloneNode(true);
    $("h2", el).textContent = label;
    container.append(el);
    return {
      fall,
      el,
      whole: $("canvas.whole", el),
      zoom: $("canvas.zoom", el),
      forces: forceChart($(".forces", el), fall, ranges.force),
      ranges,
    };
  });

  const slider = $("#time");
  slider.max = all[0].frames - 1;
  slider.value = 0;
  $("#results").hidden = false;
  resizeCharts();
  setFrame(0);
}

// --- charts (uPlot)

const seriesStyle = (label, color, extra = {}) => ({ label, stroke: color, width: 2, points: { show: false }, ...extra });

function axes(unit) {
  const grid = { stroke: css("--grid"), width: 1 };
  const common = { stroke: css("--text-secondary"), grid, ticks: { stroke: css("--grid"), width: 1, size: 4 } };
  return [{ ...common, label: "t [s]", labelSize: 18 }, { ...common, values: (u, v) => v.map((x) => `${+x.toFixed(6)} ${unit}`), size: 64 }];
}

// Solid baseline at y = 0, so the zero level stands out from the grid
const zeroLinePlugin = () => ({
  hooks: {
    draw: (u) => {
      const { min, max } = u.scales.y;
      if (!(min <= 0 && 0 <= max)) return;
      const y = Math.round(u.valToPos(0, "y", true));
      const { ctx } = u;
      ctx.save();
      ctx.strokeStyle = css("--zero");
      ctx.lineWidth = 1.5 * devicePixelRatio;
      ctx.beginPath();
      ctx.moveTo(u.bbox.left, y);
      ctx.lineTo(u.bbox.left + u.bbox.width, y);
      ctx.stroke();
      ctx.restore();
    },
  },
});

// Marks the current frame on a chart; also used by the time slider
const cursorPlugin = () => ({
  hooks: {
    draw: (u) => {
      const t = view.result && view.falls.length ? view.falls[0].fall.time[view.frame] : null;
      if (t == null) return;
      const x = u.valToPos(t, "x", true);
      const { ctx } = u;
      ctx.save();
      ctx.strokeStyle = css("--text-muted");
      ctx.setLineDash([3, 3]);
      ctx.lineWidth = devicePixelRatio;
      ctx.beginPath();
      ctx.moveTo(x, u.bbox.top);
      ctx.lineTo(x, u.bbox.top + u.bbox.height);
      ctx.stroke();
      ctx.restore();
    },
  },
});

// Direct label at the peak of selected series (values in ink, not series color)
const peakLabelPlugin = (seriesIdx) => ({
  hooks: {
    draw: (u) => {
      const { ctx } = u;
      ctx.save();
      ctx.font = `${12 * devicePixelRatio}px system-ui, sans-serif`;
      ctx.textBaseline = "bottom";
      for (const si of seriesIdx) {
        if (!u.series[si].show) continue;
        const ys = u.data[si];
        let i = 0;
        for (let k = 1; k < ys.length; k++) if (ys[k] > ys[i]) i = k;
        const x = u.valToPos(u.data[0][i], "x", true);
        const y = u.valToPos(ys[i], "y", true);
        ctx.fillStyle = u.series[si].stroke();
        ctx.beginPath();
        ctx.arc(x, y, 4 * devicePixelRatio, 0, 2 * Math.PI);
        ctx.fill();
        ctx.strokeStyle = css("--surface-1");
        ctx.lineWidth = 2 * devicePixelRatio;
        ctx.stroke();
        ctx.fillStyle = css("--text-primary");
        const text = `${u.series[si].label} ${ys[i].toFixed(2)} kN`;
        const w = ctx.measureText(text).width;
        const tx = Math.min(Math.max(x + 8 * devicePixelRatio, u.bbox.left), u.bbox.left + u.bbox.width - w);
        ctx.fillText(text, tx, y - 4 * devicePixelRatio);
      }
      ctx.restore();
    },
  },
});

function forceChart(el, fall, range) {
  const toKN = (a) => Array.from(a, (v) => v / 1000);
  const data = [Array.from(fall.time), toKN(fall.tensionMax), toKN(fall.leash), toKN(fall.anchorLeft), toKN(fall.anchorRight)];
  return new uPlot(
    {
      width: 400, height: 240,
      scales: { x: { time: false }, y: { range: () => range } },
      axes: axes("kN"),
      series: [
        { label: "t [s]", value: (u, v) => (v == null ? "–" : v.toFixed(3)) },
        seriesStyle("max webbing", css("--series-1")),
        seriesStyle("leash", css("--series-2")),
        seriesStyle("left anchor", css("--series-3"), { dash: [6, 4] }),
        seriesStyle("right anchor", css("--series-4"), { dash: [2, 3] }),
      ].map((s, i) => (i ? { ...s, value: (u, v) => (v == null ? "–" : v.toFixed(2)) } : s)),
      plugins: [zeroLinePlugin(), cursorPlugin(), peakLabelPlugin([1, 2]), clickToFramePlugin()],
      cursor: { drag: { x: true, y: false } },
    },
    data,
    el,
  );
}

function resizeCharts() {
  for (const f of view.falls) {
    const w = f.forces.root.parentElement.clientWidth;
    if (w > 0) f.forces.setSize({ width: w, height: 240 });
  }
  drawAnimations();
}
new ResizeObserver(debounce(resizeCharts, 100)).observe($("#falls"));

// --- line animation (canvas, true scale)

function setupCanvas(canvas) {
  const w = canvas.clientWidth, h = canvas.clientHeight;
  const dpr = devicePixelRatio;
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
  }
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.fillStyle = css("--surface-1");
  ctx.fillRect(0, 0, w, h);
  return { ctx, w, h };
}

// Largest of 1, 2, 5 x 10^k that is at most target and at most cap
function niceStep(target, cap) {
  if (!(target > 0)) return cap;
  const p = 10 ** Math.floor(Math.log10(target));
  const m = target / p;
  return Math.min(cap, (m >= 5 ? 5 : m >= 2 ? 2 : 1) * p);
}

// Maps world [x0,x1]x[y0,y1] into the canvas with equal scale on both axes
function viewport(w, h, x0, x1, y0, y1, pad = 14) {
  const s = Math.min((w - 2 * pad) / (x1 - x0), (h - 2 * pad) / (y1 - y0));
  const ox = (w - s * (x1 - x0)) / 2 - s * x0;
  const oy = (h + s * (y1 + y0)) / 2;
  return { s, X: (x) => ox + s * x, Y: (y) => oy - s * y };
}

function drawScene(canvas, fall, f, vp, cfg, { gridStep }) {
  const { ctx, w, h } = setupCanvas(canvas);
  const { X, Y } = vp;
  const n = fall.stride, row = f * n;

  // Faint metre grid, for scale. Horizontal lines are denser when few would fit, and labelled with their height.
  ctx.strokeStyle = css("--grid");
  ctx.lineWidth = 1;
  ctx.beginPath();
  const xMin = (0 - X(0)) / vp.s, xMax = (w - X(0)) / vp.s;
  const yMax = (Y(0) - 0) / vp.s, yMin = (Y(0) - h) / vp.s;
  const yStep = niceStep((yMax - yMin) / 4, gridStep);
  const yLines = [];
  for (let k = Math.ceil(yMin / yStep); k * yStep <= yMax; k++) yLines.push(k * yStep);
  for (let x = Math.ceil(xMin / gridStep) * gridStep; x <= xMax; x += gridStep) { ctx.moveTo(X(x), 0); ctx.lineTo(X(x), h); }
  for (const y of yLines) if (y !== 0) { ctx.moveTo(0, Y(y)); ctx.lineTo(w, Y(y)); }
  ctx.stroke();

  // Height zero (anchor level reference) as a solid line
  if (yMin <= 0 && 0 <= yMax) {
    ctx.strokeStyle = css("--zero");
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(0, Math.round(Y(0)) + 0.25);
    ctx.lineTo(w, Math.round(Y(0)) + 0.25);
    ctx.stroke();
  }

  // Height labels at the left edge, on a surface-coloured backing so the line doesn't cross the text
  ctx.font = "11px system-ui, sans-serif";
  ctx.textBaseline = "middle";
  const decimals = yStep < 1 ? 1 : 0;
  for (const y of yLines) {
    const py = Y(y);
    if (py < 8 || py > h - 20) continue; // keep clear of the edges and the grid note
    const zero = Math.abs(y) < 1e-9;
    const text = `${zero ? 0 : y.toFixed(decimals)} m`;
    ctx.font = `${zero ? 600 : 400} 11px system-ui, sans-serif`;
    const tw = ctx.measureText(text).width;
    ctx.fillStyle = css("--surface-1");
    ctx.fillRect(3, py - 7, tw + 6, 14);
    ctx.fillStyle = css(zero ? "--text-primary" : "--text-muted");
    ctx.fillText(text, 6, py);
  }
  ctx.textBaseline = "alphabetic";
  ctx.font = "11px system-ui, sans-serif";
  ctx.fillStyle = css("--text-muted");
  ctx.fillText(`grid ${gridStep} m`, 6, h - 6);

  // Webbing, one path per run of edges whose backup is slack (blue) or loaded (orange)
  const px = (i) => X(fall.x[row + i]), py = (i) => Y(fall.y[row + i]);
  const edges = n - 1, active = (e) => fall.backupActive[f * edges + e] > 0;
  ctx.lineJoin = "round";
  ctx.lineCap = "round";
  for (let e = 0; e < edges; ) {
    const on = active(e);
    let j = e;
    while (j + 1 < edges && active(j + 1) === on) j++;
    ctx.strokeStyle = css(on ? "--backup" : "--series-1");
    ctx.lineWidth = on ? 3 : 2;
    ctx.beginPath();
    ctx.moveTo(px(e), py(e));
    for (let i = e + 1; i <= j + 1; i++) ctx.lineTo(px(i), py(i));
    ctx.stroke();
    e = j + 1;
  }
  ctx.lineCap = "butt";

  // Mesh nodes as tiny dots
  ctx.fillStyle = css("--text-primary");
  for (let i = 1; i < n - 1; i++) {
    ctx.beginPath();
    ctx.arc(px(i), py(i), 1.3, 0, 2 * Math.PI);
    ctx.fill();
  }

  // Nodes where one segment connects to the next: hollow diamonds
  ctx.fillStyle = css("--surface-1");
  ctx.strokeStyle = css("--text-primary");
  ctx.lineWidth = 1.5;
  for (const i of fall.junctions) {
    const x = px(i), y = py(i), r = 4.5;
    ctx.beginPath();
    ctx.moveTo(x, y - r);
    ctx.lineTo(x + r, y);
    ctx.lineTo(x, y + r);
    ctx.lineTo(x - r, y);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();
  }

  // Anchors
  ctx.fillStyle = css("--text-secondary");
  for (const i of [0, n - 1]) ctx.fillRect(X(fall.x[row + i]) - 4, Y(fall.y[row + i]) - 4, 8, 8);

  // Leash (solid when taut) and slackliner
  const rx = fall.ringX[f], ry = fall.ringY[f], hy = fall.yHarness[f];
  const slack = cfg.slackliner;
  const taut = fall.taut[f] > 0;
  const color = css("--series-2");
  ctx.strokeStyle = color;
  ctx.lineWidth = 2;
  ctx.setLineDash(taut ? [] : [4, 3]);
  ctx.beginPath();
  ctx.moveTo(X(rx), Y(ry));
  ctx.lineTo(X(rx), Y(hy));
  ctx.stroke();
  ctx.setLineDash([]);
  // Slackliner as harness and feet points, l_leg apart (true scale)
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(X(rx), Y(hy));
  ctx.lineTo(X(rx), Y(hy - slack.l_leg));
  ctx.stroke();
  ctx.fillStyle = color;
  for (const y of [hy, hy - slack.l_leg]) {
    ctx.beginPath();
    ctx.arc(X(rx), Y(y), 3, 0, 2 * Math.PI);
    ctx.fill();
  }
  // ring
  ctx.fillStyle = css("--surface-1");
  ctx.strokeStyle = css("--text-primary");
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.arc(X(rx), Y(ry), 3.5, 0, 2 * Math.PI);
  ctx.fill();
  ctx.stroke();
}

// The figure is drawn true to scale, but never smaller than this [px/m], so it stays readable in the whole-line view
const MIN_FIGURE_SCALE = 22;

// Cartoon slackliner hanging from the harness at canvas point (hx, hy), k px per metre.
// Feet are l_leg below the harness, as in the simulation. Smiles when the leash holds, "o" mouth in free fall.
function drawSlackliner(ctx, hx, hy, k, l_leg, taut, color) {
  const P = (u, v) => [hx + u * k, hy - v * k]; // body coordinates in metres, v up
  const ink = css("--text-primary");
  const outline = Math.max(1, 0.025 * k);
  ctx.save();
  ctx.lineCap = "round";
  ctx.lineJoin = "round";

  // Legs and shoes (slightly apart, knees a little bent)
  ctx.strokeStyle = ink;
  ctx.lineWidth = Math.max(2.5, 0.1 * k);
  for (const side of [-1, 1]) {
    ctx.beginPath();
    ctx.moveTo(...P(0.17 + 0.05 * side, 0.02));
    ctx.quadraticCurveTo(...P(0.17 + 0.14 * side, -0.55 * l_leg), ...P(0.17 + 0.07 * side, -l_leg + 0.06));
    ctx.stroke();
  }
  ctx.fillStyle = color;
  for (const side of [-1, 1]) {
    ctx.beginPath();
    ctx.ellipse(...P(0.17 + 0.12 * side, -l_leg + 0.04), 0.11 * k, 0.065 * k, 0, 0, 2 * Math.PI);
    ctx.fill();
  }

  // Arms: one holds the leash, the other waves
  ctx.strokeStyle = color;
  ctx.lineWidth = Math.max(2, 0.07 * k);
  ctx.beginPath();
  ctx.moveTo(...P(0.1, 0.4));
  ctx.quadraticCurveTo(...P(-0.02, 0.5), ...P(0, 0.62));
  ctx.moveTo(...P(0.27, 0.4));
  ctx.quadraticCurveTo(...P(0.45, 0.45), ...P(0.5, taut ? 0.62 : 0.72));
  ctx.stroke();
  ctx.fillStyle = "#ffe2c6";
  for (const [u, v] of [[0, 0.62], [0.5, taut ? 0.62 : 0.72]]) {
    ctx.beginPath();
    ctx.arc(...P(u, v), 0.045 * k, 0, 2 * Math.PI);
    ctx.fill();
  }

  // Torso: a rounded orange bean with the harness belt at the bottom
  ctx.fillStyle = color;
  ctx.beginPath();
  ctx.ellipse(...P(0.18, 0.24), 0.21 * k, 0.28 * k, 0, 0, 2 * Math.PI);
  ctx.fill();
  ctx.strokeStyle = ink;
  ctx.lineWidth = Math.max(1.5, 0.04 * k);
  ctx.beginPath();
  ctx.moveTo(...P(0.03, 0.04));
  ctx.lineTo(...P(0.33, 0.04));
  ctx.stroke();

  // Head: big and round, with a helmet
  const [cx, cy] = P(0.2, 0.76);
  const r = 0.27 * k;
  ctx.fillStyle = "#ffe2c6";
  ctx.strokeStyle = ink;
  ctx.lineWidth = outline;
  ctx.beginPath();
  ctx.arc(cx, cy, r, 0, 2 * Math.PI);
  ctx.fill();
  ctx.stroke();
  ctx.fillStyle = color;
  ctx.beginPath();
  ctx.arc(cx, cy - 0.02 * k, r * 1.06, Math.PI * 1.05, Math.PI * 1.95);
  ctx.closePath();
  ctx.fill();
  ctx.stroke();

  // Face (skipped when tiny)
  if (r >= 5) {
    ctx.fillStyle = ink;
    for (const side of [-1, 1]) {
      ctx.beginPath();
      ctx.arc(cx + side * 0.35 * r, cy + 0.15 * r, Math.max(1, 0.1 * r), 0, 2 * Math.PI);
      ctx.fill();
    }
    ctx.fillStyle = "rgba(235, 104, 52, 0.35)"; // blush
    for (const side of [-1, 1]) {
      ctx.beginPath();
      ctx.ellipse(cx + side * 0.6 * r, cy + 0.42 * r, 0.16 * r, 0.1 * r, 0, 0, 2 * Math.PI);
      ctx.fill();
    }
    ctx.strokeStyle = ink;
    ctx.lineWidth = Math.max(1, 0.08 * r);
    ctx.beginPath();
    if (taut) ctx.arc(cx, cy + 0.38 * r, 0.25 * r, 0.15 * Math.PI, 0.85 * Math.PI);
    else ctx.ellipse(cx, cy + 0.55 * r, 0.11 * r, 0.15 * r, 0, 0, 2 * Math.PI);
    ctx.stroke();
  }
  ctx.restore();
}

const ZOOM_HALF_WIDTH = 4; // [m], as in tools/animate_fall.py

function drawAnimations() {
  if (!view.result) return;
  const cfg = view.result.config;
  for (const v of view.falls) {
    const { fall, ranges } = v;
    const f = Math.min(view.frame, fall.frames - 1);
    const top = Math.max(ranges.yHi, 0) + 2;
    let c = v.whole;
    drawScene(c, fall, f, viewport(c.clientWidth, c.clientHeight, 0, cfg.L, ranges.yLo - 1, top), cfg, { gridStep: 10 });

    c = v.zoom;
    const rx = fall.ringX[f], ry = fall.ringY[f];
    const zoomH = cfg.slackliner.l_leash + cfg.slackliner.l_leg + 3;
    drawScene(c, fall, f, viewport(c.clientWidth, c.clientHeight, rx - ZOOM_HALF_WIDTH, rx + ZOOM_HALF_WIDTH,
      ry - zoomH + 1.5, ry + 1.5), cfg, { gridStep: 1 });
  }
}

function setFrame(f) {
  view.frame = f;
  $("#time").value = f;
  const t = view.falls[0].fall.time[f];
  $("#time-label").textContent = `t = ${t.toFixed(3)} s`;
  for (const { forces: u } of view.falls) {
    // Park the chart cursor on the frame so the legend shows its values
    if (!u.hovering) u.setCursor({ left: u.valToPos(t, "x"), top: u.valToPos(u.data[1][f], "y") }, false);
    u.redraw(false);
  }
  drawAnimations();
}

// Clicking a chart moves the time cursor there
const clickToFramePlugin = () => ({
  hooks: {
    ready: (u) => {
      u.over.addEventListener("mouseenter", () => (u.hovering = true));
      u.over.addEventListener("mouseleave", () => (u.hovering = false));
      u.over.addEventListener("click", () => {
        if (u.cursor.idx == null) return;
        stop();
        setFrame(Math.min(u.cursor.idx, view.falls[0].fall.frames - 1));
      });
    },
  },
});

$("#time").oninput = (e) => { stop(); setFrame(+e.target.value); };

let raf = 0, lastTs = 0, simTime = 0;
function stop() {
  view.playing = false;
  cancelAnimationFrame(raf);
  $("#play").textContent = "▶";
}
function tick(ts) {
  const time = view.falls[0].fall.time;
  simTime += ((ts - lastTs) / 1000) * +$("#speed").value;
  lastTs = ts;
  let f = view.frame;
  while (f + 1 < time.length && time[f + 1] <= simTime) f++;
  if (f !== view.frame) setFrame(f);
  if (f >= time.length - 1) return stop();
  raf = requestAnimationFrame(tick);
}
$("#play").onclick = () => {
  if (view.playing) return stop();
  if (!view.falls.length) return;
  if (view.frame >= view.falls[0].fall.frames - 1) setFrame(0);
  view.playing = true;
  $("#play").textContent = "❚❚";
  simTime = view.falls[0].fall.time[view.frame];
  lastTs = performance.now();
  raf = requestAnimationFrame(tick);
};

// Redraw canvases and charts in the new palette when the color scheme flips
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => view.result && showResults(view.result));

// ---------------------------------------------------------------- start

const ready = await workerReady;
if (ready.error) showError(`Failed to load the simulator: ${ready.error}`);
else {
  // A catalogue edited in an earlier visit wins over the shipped one, unless it no longer parses
  shippedWebbingsYaml = ready.webbingsYaml;
  const saved = store.get(WEBBINGS_KEY);
  webbingsBox.value = saved ?? shippedWebbingsYaml;
  await setWebbings(webbingsBox.value, { fillTable: true });
  if (webbingsYaml === "") {
    webbingsBox.value = shippedWebbingsYaml;
    await setWebbings(shippedWebbingsYaml, { fillTable: true });
  }
  fillSavedSetups();
  await loadSetup(setupSelect.value);
  runButton.disabled = false;
  runButton.textContent = "Run simulation";
}

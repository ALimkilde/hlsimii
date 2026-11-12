# Adding or removing a YAML example on the website

The website's example dropdown is not generated from `examples/`: each example is listed in two
places, and both must agree. Paths are relative to the repo root.

## Add an example

1. Put the file in `examples/` (e.g. `examples/1km.yaml`).
2. Add the file name to `HLSIM_EXAMPLES` in `web/CMakeLists.txt`. This list decides what the build
   copies into `dist/`:
   ```cmake
   set(HLSIM_EXAMPLES 110m.yaml short_and_static.yaml 250m.yaml 1km.yaml)
   ```
3. Add an option to the "Examples" group of `<select id="example">` in `web/site/index.html`.
   `value` is the file name, the text is what the dropdown shows:
   ```html
   <option value="1km.yaml">1km</option>
   ```
4. Run `web/build.sh` and reload the page.

## Remove an example

1. Remove its name from `HLSIM_EXAMPLES` in `web/CMakeLists.txt`.
2. Remove its `<option>` from `web/site/index.html`.
3. Optionally delete the file from `examples/`.
4. Rebuild.

## Things to know

- **The first option is the default.** The page loads the first example in the dropdown when it
  opens (currently `110m`), so reorder the options to change the default.
- **Removed files stay in `dist/`.** The build copies into `dist/` but never deletes, so a removed
  example stays there (unused) until you delete it or clear `dist/`.
- **`110m.yaml` is also used elsewhere.** It is the default config of the smoke test
  (`web/test/smoke.mjs`, line 6, used by `./build.sh test`) and is mentioned twice in
  `web/README.md`. Update those if you remove or rename it.
- **Webbing names must exist.** The segments of an example must use webbing names from
  `data/webbings.yaml`, or the simulator reports an error when the example is loaded.
- **Saved setups are separate.** Setups saved with "save current setup" live in the browser's
  localStorage and need none of this.

## Possible improvement

Generate the dropdown options from `HLSIM_EXAMPLES` at build time (or from a small JSON list the
page fetches), so `web/CMakeLists.txt` is the only place to edit.

# Cell Quant Studio Universal - handoff

Claude Code reads this file automatically at the start of a session in this repository.

Repository: https://github.com/KiViKatt/ImageJ-SingleCellQuantificationPlugin. The 2.3.0 package is on `main`; the earlier macros are in `legacy/`.

## What this is

A Fiji/ImageJ application, written in Jython 2.7, for quantifying multichannel fluorescence images of human brain organoid sections:

- per-channel preparation and independent segmentation;
- ROI-based object counts, marker-positive area and intensity;
- folder batches, from the window and headless;
- validation against manual annotation;
- six optional analysis modules (tissue outline, depth profile, threshold spectrum, tuning, corrected counts, group robustness).

The owner is preparing a short first-author software paper.

## Status (2026-11-19)

Version 2.3.0 has been run and tested in the owner's real Fiji: Jython 2.7.4, ImageJ 1.54p, Java 21.

**Automated tests:**

- `tests/run_tests.py`: 126 tests, about 4300 checks, all pass.
- `tests/gui_smoke.py`: 30 window steps, all pass.
- Both also pass on macOS (Apple silicon, the same Fiji versions; 4302 checks there), run on 2026-11-19. No real image was used on macOS.
- `tools/test_cellquant_figures.py`: 6 tests, all pass.

**Real data (the owner's Zeiss Airyscan CZI files, 4 channels x 4 Z, 4084 x 4084, 16-bit, 0.2076 micron per pixel; read-only, unpublished):**

- A 2.3 run reproduces the owner's 2.2 run of 2026-11-18 exactly for channels 2-4, on every summary column and all object rows.
- That 2.2 run's channel 1 had been measured on an image edited in memory; one Backspace on C1 Z1 with the ROI active reproduces it exactly.
- A two-file batch, the window, and the headless script through the Jaunch launcher all worked.

**Not verified:**

- Linux;
- the last engine fixes (holed-object intensity, overlay-based ROI collection, path shortening) on a real CZI after they were applied (synthetic tests pass);
- real multi-series and time-lapse files;
- a person using the window by hand after the last changes;
- validation against manual annotation of real images (the paper's main gap; see `docs/VALIDATION_STUDY.md`).

## Files

| File | Role |
|---|---|
| `Cell_Quant_Studio_final.py` | Launcher. `execfile`s the engine, then `load_engine_modules(globals(), package_dir)`, then the window, into one namespace. Holds `LAUNCHER_VERSION`. |
| `cellquantengine.py` | Analysis engine. No Swing. Holds `APP_VERSION`. Also the optional-module framework (`RunExtension`, `register_*`, `load_engine_modules`). |
| `cellquantcore.py` | Swing window (about 166 KB). Holds `CORE_VERSION`. Analyses tab built from `SETTINGS_FORMS`; background jobs via `start_background`. |
| `cellquanttissue.py`, `cellquantdepth.py`, `cellquantspectrum.py`, `cellquanttune.py`, `cellquantcorrect.py`, `cellquantrobust.py` | Optional modules, executed into the engine namespace. Each ends with `register_engine_module(<file>, <version>)`. |
| `cellquant_batch_headless.py` | Folder batch from a saved profile, no window. SciJava `#@` parameters. Loads the modules. |
| `tests/run_tests.py` | Test runner and core tests; loads every `tests/test_*.py` (each ends with `TESTS.append(...)`). |
| `tests/test_*.py` | Module tests, engine-fix regression tests, Bio-Formats tests, ROI-file tests. |
| `tests/gui_smoke.py` | Drives the real window with synthetic images; auto-dismisses dialogs. |
| `tests/run_tests.sh`, `tests/cli_runner.py` | Run Jython scripts with Fiji's Java and classpath from a terminal; exit status 0/1. |
| `tools/cellquant_figures.py` | CPython 3, no dependencies: SVG figures from result tables. Tests in `tools/test_cellquant_figures.py`. |
| `docs/VALIDATION.md`, `docs/VALIDATION_STUDY.md` | How to validate; the planned study for the paper. |
| `paper/manuscript_draft.md`, `paper/references.md` | Manuscript draft (F1000 software-tool format) and verified references. |
| `README.md`, `CHANGELOG.md`, `CITATION.cff`, `RELEASE_CHECKLIST.md` | User docs, release notes, citation metadata, the owner's remaining steps. |
| `legacy/` | The earlier macros and launcher, unchanged. |

All version constants and every module's registered version must match. The launcher, the window and `load_engine_modules` refuse to run otherwise. Bump them all together.

## Rules for the code (Jython 2.7 inside Fiji)

1. **Python 2.7 syntax that is also valid Python 3.** No f-strings, no `nonlocal`, no keyword-only arguments, no print statements. Check with `python -I -m py_compile <file>` (then delete `__pycache__`).
2. **Plain ASCII source, no encoding declaration line.** Fiji's script engine rejects a coding line in some launch paths. The only non-ASCII escapes are the two unicode-escape literals in `is_micrometer_unit`. Tests build the micro sign with `b"\xc2\xb5".decode("utf-8")`.
3. **Catch errors with `except ANY_ERROR`, never `except Exception`.** Under Jython `except Exception` misses Java exceptions. `ANY_ERROR = (Exception, JavaThrowable)`.
4. **Use `to_text()`, never `str()`, on anything that may hold text.** `str()` on numbers is fine. `to_text(float)` uses `repr` so settings keep every digit.
5. **Jython compiles each method into a JVM method capped at 64 KB.** Big files load fine (the window is about 166 KB), but keep functions small and avoid huge literal tables.
6. **Python 2 semantics:**
   - integer division floors;
   - list comprehensions leak their variable into the enclosing scope (in the shared launcher namespace too, so give them unusual names);
   - no `nonlocal` (use a dict).
7. **The engine and the modules never touch Swing.** Loading them with `CELLQUANT_NO_GUI = True` opens nothing.
8. **Restore every ImageJ global you change in `finally`:** `Prefs`, the Analyze measurement redirect, and the like. The test runner checks Prefs, the redirect, the ROI Manager instance and the open-image count.

## Invariants found the hard way - do not regress

- **Objects are confined to the true outline of the parent ROI:** `analyze_objects` ANDs the mask with the ROI mask first (`ParticleAnalyzer` only honours the bounding rectangle).
- **Particle ROIs are collected with `ParticleAnalyzer.OVERLAY` from the work image's overlay, not through a `RoiManager`.** `RoiManager(False)` throws `HeadlessException` under `--headless`, and closing a hidden manager disconnects the user's ROI Manager. Do not use `ADD_TO_OVERLAY` (it returned no ROIs).
- **Per-object intensity is measured on the object's own pixels**, through a ParticleAnalyzer pass redirected to the prepared image. The redirect is saved and restored. Bit-identical to the old per-ROI measurement for objects without holes.
- **An empty threshold range gives an all-background mask** with `mask_data["empty_range"] = True`, not an error. The tune and robust modules rely on that flag.
- **Edge exclusion is done by the engine** with a one-pixel band of the ROI mask.
- **ROI sets are read and written directly** (`save_roi_array`, `load_rois_from_file`). The file type is sniffed from content: ROI Manager saves a one-ROI set as a plain `.roi` even under a `.zip` name.
- **Binary operations go through `binary_morphology()`** with explicit iterations and neighbour count, never `IJ.run`.
- **8-bit conversion before Local Phansalkar uses `convertToByte(True)` after `setMinAndMax`.** Local Phansalkar and CLAHE are called through their classes (pixel-identical to the menu commands), so they work without menus and without a display.
- **`IJ.run` does not raise when a plugin fails.** Verify the effect.
- **The appearance is a global Swing look.** `look_install(name)` switches Fiji to FlatLaf dark or light with the palette of the chosen look (`LOOK_PALETTES`: Dark and Light, each with its FlatLaf colors under `"extras"` and the colors of the painted parts) and the fonts (`fonts/`, loaded by `look_fonts_load`) when the window opens, and `look_restore()` hands back the previous look object when it closes (`LookWindowListener`) or when Fiji default is chosen. Keep the look object itself: under Jython `obj.getClass().getName()` raises. A variable font loads at its thinnest weight in Java, so only fixed-weight font files are bundled (IBM Plex Sans from the IBM/plex repository, Chakra Petch from google/fonts). Styling uses FlatLaf client properties and UI defaults only, so other looks ignore it. Painted parts take every color from `look_palette()` at painting time, never from a constant, so a change of look recolors what is already shown; a new painted color goes into both palettes. Custom-painted text goes through `paint_fitted()` and fonts through `look_font()`: nothing may overlap. Tests press buttons by their text: keep button texts unique and unchanged.
- **One layout pass shapes every tab.** `tidy_tab()` runs once after `build_ui`: one right-aligned label column per tab (long labels wrap), buttons behind fields of one width, button rows of one size, forms in scroll areas pinned to the top. Tabs other than Session and Results sit in a centered column (`in_column`). Build new rows with `add_form_row` / `add_grid_row` and plain `FlowLayout` button rows and they are shaped the same way; do not size widgets by hand.
- **CLAHE runs in strips on worker threads** (`clahe_in_strips`), overlapping by the block radius, and must stay pixel-identical to the one-piece call (test "CLAHE in strips on several cores is pixel-identical"). An image without a display range, a stack, or a float or RGB plane goes through the one-piece path.
- **The Results tab only displays.** `results_from_summary()` reads the summary rows of the finished run; it must never measure or write anything, and a failure in `show_results()` must not fail the analysis.
- **Window threading:**
  - Status messages before the window is displayable are appended directly; afterwards through `run_on_ui`.
  - Never use `invokeAndWait` from a worker thread.
  - Script-thread calls to window actions go through `@ui_action` / `call_on_ui`.
  - Long jobs run in `start_background` workers.
- **The source image may be edited in memory.** `confirm_source_changes` asks before measuring it, and edited prepared images are re-prepared, not measured.
- **Modules load before the window** (the Batch tab's ROI-policy combo is built from `BATCH_ROI_POLICIES`). A refused module leaves nothing registered.
- **Batches refuse Display-minimum thresholds.** Images with a calibration function are refused for segmentation.
- **Long result paths are shortened** (`short_name_base`) only when the longest path would exceed 259 characters. Short names stay byte-identical to 2.2.
- **Results for rectangular ROIs and the whole image stay identical to 2.2** at default Fiji options, for objects without holes.

## Compatibility policy

The owner's standing instruction: do not remove features, only improve and add, and keep things compatible with the previous version.

- Profiles saved by 2.2.x must keep loading; new keys get defaults.
- Every table keeps its 2.2 columns in the same order; new columns are appended (module columns last).
- Property keys, output file names and folder names stay as they are (except the documented long-name shortening).
- Any change that alters numbers goes under "Results that change" in `CHANGELOG.md`.

## Running things

**Fiji on the owner's PC:** ask the owner (or see Claude's memory) for its folder. An old 2.2 copy of the app is installed in its `plugins/Scripts/Cell Quant Studio/`. Do not overwrite it without asking.

**Tests from a terminal** (Git Bash). These use Fiji's bundled Java, with `Prefs.load` first:

- `bash tests/run_tests.sh <Fiji folder>` runs the full suite, about 1 minute.
- `bash tests/run_tests.sh <Fiji folder> "<absolute path>/tests/gui_smoke.py"` runs the window test. It needs an absolute path, because the runner changes into the Fiji folder, and it opens windows for about 30 s.
- `CELLQUANT_TEST_FILTER=word,word` selects tests.
- `CELLQUANT_JAVA_OPTS=-Xmx16g` gives Java more memory for large images; `-Djava.awt.headless=true` runs without AWT.
- Always wrap runs in `timeout`, and check that no `java.exe` is left running.

**Headless through the real launcher** (works; ends after about 8 s, but prints nothing to Git Bash or PowerShell, so read `batch_log.txt`):

```text
fiji-windows-x64.exe --jaunch-skip-console-check --headless --console --run "C:/path/script.py" 'name="value",...'
```

**Things to know about Fiji script parameters:**

- Fiji persists omitted `#@` parameters from the last run, so give them all (`roi_policy` uses `persist=false`).
- Fiji exits with status 0 even when a script fails.

**Script Editor equivalent:** `org.scijava.script.ScriptService.run(File, True, inputs).get()` from a Jython script.

**Viewing SVG figures:** Fiji's Batik jars (`batik-*`) can render SVG to PNG through `GVTBuilder` (no transcoder jar).

## Roadmap status

Roadmap items 1-5 are implemented as optional modules with tests (depth profile, threshold spectrum, auto-tune, corrected counts, group robustness), plus the automatic tissue outline.

Parked: recording-site tissue composition around a probe track (no suitable organoids yet).

For new features: engine or module function first with a known-answer test, then the window control, then README and CHANGELOG. Outputs go in new files or appended columns.

## Toward a paper

- Literature check (2026-11-18): every module has prior art, and `paper/references.md` cites it. The most defensible contribution is the combination in one auditable Fiji tool of validation against manual points, calibrated counts, full provenance and a routine robustness (multiverse) report.
- The substance must come from validation on real organoid images (`docs/VALIDATION_STUDY.md`): 40-60 fields, two blind annotators, StarDist and Cellpose comparison, Bland-Altman, ICC, F1.
- Venues: F1000Research Software Tool Article or JORS after a Zenodo DOI. eNeuro "Open Source Tools and Methods" after validation. JOSS needs about 6 months of public history.

## Open decisions for the owner

- **Licence: decided 2026-11-20, GPL-3.0-or-later** (`LICENSE`, `license:` in `CITATION.cff`, an SPDX header in every source file). Chosen to match the Fiji ecosystem and keep the option of selling separate commercial licences later, which needs the owner to hold all copyright (so outside contributions would need a contributor agreement). Every new source file gets the same two-line header.
- **Demo profile.** The owner's own demo profile was never added to `examples/`.
- **MAP2 threshold.** The main preset uses a by-eye Display minimum. Capture a fixed Manual value, or validate between two raters.
- **Old data.** Re-run 2.2.x counts made with non-rectangular ROIs, and the 2026-11-18 org1-ROI1 run.
- **Tissue outline.** It sits about 1.6 smoothing widths outside a sharp edge on synthetic data. Compare with hand-drawn outlines before relying on depth near the surface.

## Working with the owner

- A researcher, not a software engineer: explain findings in plain terms, lead with what changed in the results.
- Sometimes prompts from a phone; keep replies short and self-contained.
- Say plainly what was and was not tested.

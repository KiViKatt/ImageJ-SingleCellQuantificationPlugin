# Remaining steps before a public release

A working note for the maintainer. Delete it once the steps are done.

## What has been tested (2026-11-19)

- In the real Fiji on the development PC (Jython 2.7.4, ImageJ 1.54p, Java 21):
  - 125 engine and module tests, about 4300 checks, all pass;
  - the 29-step automated window test passes.
- On macOS (Apple silicon, same Fiji versions): both automated tests pass. No real image was used there.
- On the author's Zeiss CZI files (read-only):
  - A 2.3 run reproduces the 2.2 run of 2026-11-18 exactly for channels 2-4. That covers every summary column and all 36,453 object rows.
  - Channel 1 of that 2.2 run was measured on an image edited in memory; see CHANGELOG.
  - A two-file batch took about 65 s per image and about 4.6 GB of memory at peak.
  - The tissue outline and the depth profile were run on one CZI.
- Headless batches through the Fiji launcher, and the launcher through Fiji's script service.

Not tested yet:

- Linux;
- a person using the window by hand on real data after the last round of changes;
- time-lapse CZI files and multi-series files with label images (only synthetic ones);
- Excel or File Explorer on result paths over 259 characters (they are now shortened);
- the final engine fixes on a real CZI (the holed-object intensity change was measured on a CZI before it was applied; all synthetic tests pass after it).

## 1. Checks to do yourself in Fiji (about 30 minutes)

1. **Install.** Copy the launcher, the engine, the window file and the six `cellquant*.py` modules into `Fiji.app/plugins/Scripts/Cell Quant Studio/`, replacing the 2.2 files there. Run `Cell_Quant_Studio_final.py`. The window should open with seven tabs (Session, Channels, ROIs and output, Batch, Validation and QC, Analyses, Help).
2. **Real image, rectangle.** Open one CZI and load your usual profile. Prepare, preview, and run with a rectangular ROI. Compare `roi_summary.csv` with an older run: channels without holed objects should be identical.
3. **Freehand ROI.** Object counts can be lower than in 2.2; positive area and intensities should be identical.
4. **Second image** without reloading ROI Manager.
5. **Analyses tab.** Click **Detect tissue outline**, check the outline, switch on the depth profile and run. Look at `depth_profile.csv` or make a figure with `tools/cellquant_figures.py depth`.
6. **Small batch** of 2-3 CZIs, first with "Whole image", then with "Automatic tissue outline".
7. **One validation** with a few multi-point marks.

## 2. Publishing

1. ~~Choose a licence.~~ Done: GPL-3.0-or-later (`LICENSE`, `CITATION.cff`, source headers). Before selling any commercial licence, check your institution's IP policy.
2. **Tag a release** on GitHub (for example `v2.3.0`) and connect the repository to Zenodo to get a DOI. Add the DOI to `CITATION.cff` and the README.
3. **Optional: a Fiji update site**, so others install with Help > Update.

## 3. The validation study (the substance of the paper)

Follow `docs/VALIDATION_STUDY.md`:

1. Choose 40-60 fields from at least 3 batches.
2. Two blind annotators mark them.
3. Compare with StarDist and Cellpose.
4. Report Bland-Altman, ICC and F1.

Then fill the bracketed places in `paper/manuscript_draft.md`.

## 4. Decisions still open

- **MAP2 threshold.** The main preset uses a by-eye Display minimum. For batches and for the paper, capture it as a fixed Manual value, or validate it between two raters.
- **Old data.** Re-run 2.2.x object counts made with non-rectangular ROIs, and the 2026-11-18 org1-ROI1 run.
- **Demo profile.** Add your own `CellQuant_4Channel_Demo_SOXCaspaseMAP2DAPI.properties` to `examples/`.

# Changelog

## 2.3.0 (2026-11-19)

Nothing was removed. Profiles saved by 2.2.x load unchanged, and every table keeps its 2.2 columns in the same order; new columns are added after them. This release was run and tested in a real Fiji (Jython 2.7.4, ImageJ 1.54p, Java 21): 125 automated engine and module tests, a 30-step automated window test, and comparisons with a 2.2 run on the author's Zeiss CZI images.

### Results that change compared with 2.2.x

Read this section before comparing new numbers with old ones.

1. **Object counts inside non-rectangular parent ROIs.**
   - **Before:** objects were detected inside the *bounding rectangle* of a parent ROI, not inside its outline. With a freehand, polygon or oval ROI, objects outside the ROI but inside its bounding box were counted, and the count was divided by the true ROI area.
   - **Affected:** object count, density, total and mean object area, and the object table. Positive area, intensity, rectangular ROIs and **Whole image** were not affected.
   - **Now:** objects are detected inside the true outline; an object crossing the outline is clipped there. **Counts made in 2.2.x with non-rectangular ROIs should be re-run.**
2. **Edge exclusion with non-rectangular ROIs.** It used to test the bounding rectangle only. It now removes every object with a pixel on the ROI outline or the image border.
3. **Open and close with more than one iteration.** "Open iterations = N" used to repeat a one-step opening N times, which gives the same mask as N = 1. It is now N erosions followed by N dilations (ImageJ's meaning of iterations); close is the reverse. Iterations of 0 or 1 are unaffected.
4. **Installations with non-default Fiji options.** Results no longer depend on **Process > Binary > Options** or on **Edit > Options > Conversions > Scale when converting**. With default options the masks are identical to before.
5. **Per-object intensity of objects with holes.**
   - **Before:** an object's Mean, StdDev, Min, Max, Median, Integrated and Raw Integrated Density (including the cross-channel columns) were measured inside its outer outline, so they included its holes and any object lying inside them.
   - **Now:** they are measured on the object's own pixels.
   - **Size of the change:** on the author's CZI images 0.2-3.5% of objects were affected, with mean changes of 0.7-16% (largest for ring-like or net-like staining). Objects without holes give bit-identical values.
6. **Images in which no pixel passes the threshold** (negative controls, empty fields, a threshold above the brightest pixel) now give zero positive area and zero objects. Before, the image failed (a batch reported it as FAILED).
7. **"Current Z" after moving the Z slider** now measures the slice shown. Before, the window could reuse the prepared image of the previous slice.
8. **A prepared channel image edited by hand** (Brightness/Contrast "Apply", Edit > Clear, a filter) is no longer measured; the channel is prepared again from the source.
9. **Long file names.** When a result path would exceed 259 characters, the run folder and the file names inside it use the start of the image name plus an 8-character code. Names that fit are unchanged. `Run_ID` then holds the shortened folder name; `Source_Image` keeps the full name.
10. **Line endings.** Tables and logs are written with LF line endings on every system. On Windows, 2.2.x wrote CRLF. Columns and values are unaffected.

**Note on the run of 2026-11-18 (org1-ROI1, version 2.2.0).** Its channel 1 was measured on an image that had been edited in memory. Its numbers are reproduced exactly when slice 1 of channel 1 is cleared inside the parent ROI, which one Backspace with that ROI active does. Channels 2-4 matched the unedited data exactly. Re-run that image. The window now warns before measuring an edited image.

### Added

- **Results tab.** After an analysis the window shows headline numbers (object counts, positive-area share), a chart of the object counts per parent ROI, a merged color picture of all channels with the parent ROIs and detected objects outlined, and the object overlay image. It reads the summary table and the saved files of the run; no table, file or number changes.
- **A new appearance** in the macOS design language (FlatLaf's macOS theme, dark or light): neutral surfaces, one blue accent, one filled button for the main action, rounded group outlines, and text that is shortened instead of overlapping. It applies to all Fiji windows while Cell Quant Studio is open; the previous look is put back when the window closes. **Appearance > Fiji default** at the top right switches it off, and the choice is remembered. Without FlatLaf the window keeps Fiji's look.
- **Consistent layout on every tab:** one aligned label column per tab, buttons of one size within a row and behind fields, content in a centered column instead of stretched across the window, group titles above their outline, and forms that start at the top of their area. The Session tab's actions are two rows: getting an image, then the steps of an analysis. No button text changed.
- **Session tab:** a drop zone for the image file and a small picture of the source image. The workflow text of the old drop area is now its tooltip.
- **Optional analysis modules**, each in its own file, with settings saved in profiles and a new **Analyses** tab. With a module off, every table is byte-for-byte unchanged (tested).
  - `cellquanttissue.py`: automatic tissue outline (Gaussian smoothing in micrometers, automatic threshold on a log scale, hole filling, region-size filter); window button and the batch region policy "Automatic tissue outline"; `tissue_outline.csv`, outline ROIs and a check image; `Tissue_Outlines/` in batches.
  - `cellquantdepth.py`: depth profile from the ROI outline (`depth_profile.csv`; object depth columns in `objects.csv`).
  - `cellquantspectrum.py`: threshold spectrum of marker coverage (`threshold_spectrum.csv`; coverage, sensitivity, elasticity, stable range and normalized mean in `roi_summary.csv`), after Johnson and Walker (2015).
  - `cellquanttune.py`: tuning of one channel's settings against manual point marks, with leave-one-field-out held-out scoring and a tuned profile.
  - `cellquantcorrect.py`: error-corrected counts from validation data, with 95% prediction intervals, Bland-Altman statistics and a leave-one-out check of whether correction helps.
  - `cellquantrobust.py`: group-level robustness over a grid of settings (specification curve, Hedges g, Welch and Mann-Whitney tests, unit averaging).
- **Folder batch processing**:
  - the Batch tab and `cellquant_batch_headless.py`;
  - region choices: whole image, an ROI file per image, one ROI Manager set for all images, or module policies;
  - combined tables written image by image, including `batch_<table>.csv` for module tables;
  - a per-file report;
  - background execution with progress and cancel;
  - Bio-Formats opening and multi-series support.
- **Headless batches** load the optional modules, accept module region policies, and count objects without a display.
- **Validation and QC tab**: object detection scored against manual point marks (precision, recall, F1, count error, colour-coded overlay); positive-area mask against a reference mask (Dice, Jaccard); parameter sensitivity sweeps.
- **Window**:
  - analyses, validation and sensitivity run in the background, so the window stays usable, with a progress bar and Cancel;
  - a warning before preparing or measuring a source image with unsaved in-memory edits, recorded in `run_log.txt`.
- **Binary neighbour count** per channel (1-8), stored in the profile.
- **New output columns**:
  - `objects.csv`: `Touches_ROI_Edge`, `Object_MinFeret`, `Object_Aspect_Ratio`, `Object_Roundness`, `Object_Solidity`;
  - `roi_summary.csv`: edge-touching and excluded object counts;
  - `thresholds.csv`: neighbour count, pad-edges state, 8-bit rescaling range and ImageJ version;
  - `batch_report.csv`: series name, image size, channels, slices, frames and calibration.
- **Tests**:
  - `tests/run_tests.py` (engine and modules; every `tests/test_*.py` is loaded);
  - `tests/gui_smoke.py` (the window);
  - `tests/run_tests.sh` and `tests/cli_runner.py` (run from a terminal); a broken test file is reported and the rest still run.
- `tools/cellquant_figures.py`: SVG figures from result tables (Bland-Altman, depth profiles, threshold spectra, specification curves), with tests.
- An example profile, a validation guide (`docs/VALIDATION.md`), a validation-study plan (`docs/VALIDATION_STUDY.md`), a manuscript draft and verified references (`paper/`), and citation metadata.

### Fixed

- **Channel settings erased by an image with fewer channels.** Opening a 1- to 3-channel image, or clicking "Use active image" on a prepared window, switched off every higher channel and set it to Ignore. The settings did not come back with the next 4-channel image, and a profile saved afterwards stored the erased values. The tabs of missing channels are still greyed out, but their settings are kept and are simply not used for that image. Folder jobs (batch, tuning, robustness) use every enabled channel whatever image is open.
- **Validation, sensitivity and "Capture display minimum" with a Display-minimum threshold.** When the prepared image was missing or out of date, these rebuilt it and used its default display range instead of the minimum chosen in Brightness/Contrast, without any warning. They now stop with the same message as "Run analysis".
- **A sensitivity or tuning sweep stopped by one ImageJ error.** An error raised inside ImageJ or a plugin for one factor or one combination of settings ended the whole sweep, because only Python errors were caught there. It is now noted for that factor or combination, like any other failure, and the sweep goes on. Cancelling still stops it.
- **Blended pixels around the label of the tissue check image on macOS.** macOS always smooths text, which left pixels between yellow and gray around the letters. The label is now drawn in plain yellow on every system. On Windows the picture is unchanged.
- **ROI Manager disconnected after an analysis.** Hidden ROI Manager objects were closed after each run, so ImageJ lost the visible manager. The next run reported "ROI Manager contains no parent ROIs". Objects are now found without any ROI Manager.
- **A measurement redirect set in Analyze > Set Measurements** could stop object detection; it is now saved, cleared and restored around each analysis.
- **Non-ASCII text.** A micrometer sign or an accented file or channel name could stop loading or export. Text handling is now Unicode-safe and tables are UTF-8.
- **CLAHE and Auto Local Threshold without an image window or without a display.** Both are now called through their libraries, with pixel-identical results. A failed local threshold is detected instead of passing silently.
- **Errors inside ImageJ or a plugin.** Java errors were not caught under Jython and skipped the error dialog. Both kinds are now caught, and a batch continues with the next file.
- **A one-ROI set that 2.2 saved under a `.zip` name** (ROI Manager writes such a set as a plain `.roi`) now loads.
- **Batch image names.** Images opened by Bio-Formats are named after their file, not after an abbreviated full path.
- **"Analyze every series"** really opens every series under Jython (before, only the first). Series open one at a time; label, macro and thumbnail images are skipped and reported. Unreadable files are released, and Bio-Formats console noise is reduced.
- **Images with a calibration function** (for example signed 16-bit) were thresholded on raw values against calibrated limits. They are now refused for segmentation with a clear message.
- **Other images:**
  - An RGB "Intensity only" channel no longer triggers a lookup-table error.
  - Auto thresholds on images with an inverted lookup table pick the right phase.
- **Profiles:**
  - Settings are saved with every digit (before, floats were cut to 12 digits).
  - An unknown frame mode or automatic method is refused when the profile loads.
  - A bad value in a switched-off channel no longer blocks a profile.
  - Channel settings missing from an old profile get their defaults in the window.
- **Window:**
  - Window actions called from a script thread run on the event thread.
  - A prepared image or mask can no longer be chosen as the source.
  - The Log no longer shows "None" lines.
- **Batches:**
  - Each image's table rows are freed before the next image opens (about 0.6 GB less held for 65,000 objects).
  - A batch writes tables faster, with byte-identical output.
- The version written to every table matches the package. The launcher, engine, window and modules refuse to run together if their versions differ, and a refused module leaves nothing registered.
- A failed run no longer leaves an empty run folder when a display-minimum threshold is rejected.

### Changed

- The application is a launcher, an engine, a window, six optional modules and a headless batch script. Each file stays within what Fiji's Jython can compile.
- The automatic threshold is computed from the full data range explicitly.
- ROI sets are read and written directly instead of through a hidden ROI Manager, so a damaged file raises an error message instead of opening a dialog.
- The window's default output is still the image folder; the README recommends an output folder near the drive root.

# Changelog

## 2.3.0

Nothing was removed. Profiles saved by 2.2.x load unchanged, and every table keeps its 2.2 columns in the same order; new columns are added after them.

### Results that change compared with 2.2.x

Read this section before comparing new numbers with old ones.

1. **Object counts inside non-rectangular parent ROIs.** In 2.2.x, objects were detected inside the *bounding rectangle* of a parent ROI, not inside its outline. With a freehand, polygon, or oval ROI, objects outside the ROI but inside its bounding box were counted, and the count was then divided by the true ROI area. Object count, object density, total object area, mean object area, and the object table were affected. Positive area and intensity were not affected. Rectangular ROIs and **Whole image** were not affected. In 2.3.0 objects are detected inside the true ROI outline; an object crossing the outline is clipped there. **Counts made in 2.2.x with non-rectangular ROIs should be re-run.**
2. **Edge exclusion with non-rectangular ROIs.** "Exclude objects touching ROI/image edges" used to test the bounding rectangle only. It now removes every object with a pixel on the ROI outline or the image border.
3. **Open and close with more than one iteration.** "Open iterations = N" used to repeat a one-step opening N times, which gives the same mask as N = 1. It is now N erosions followed by N dilations (ImageJ's own meaning of iterations); close is the reverse. Profiles with open or close iterations of 0 or 1 are unaffected. Erode and dilate iterations are unchanged.
4. **Installations with non-default Fiji options.** Results no longer depend on **Process > Binary > Options** (iterations, count, pad edges) or on **Edit > Options > Conversions > Scale when converting**. On a Fiji with default options the masks are identical to before. With "Scale when converting" switched off, 2.2.x turned the whole local-threshold copy white; that no longer happens.

### Fixed

- **ROI Manager disconnected after an analysis.** Each run closed hidden ROI Manager objects, which made ImageJ forget the visible ROI Manager window. The next run reported "ROI Manager contains no parent ROIs", and "Open ROI Manager" or "Load ROI set" opened a second, empty manager. The ROI Manager now stays connected.
- **Non-ASCII text.** A calibration unit such as the micrometer sign, or an accented file or channel name, could stop image loading or table export. All text handling is now Unicode-safe and tables are written as UTF-8.
- **CLAHE without an open image window.** The CLAHE menu command fails when no image window exists, and Fiji reports that only in the Log. CLAHE is now called through its library, which gives pixel-identical results and works in batches. Local thresholding is checked to have produced a two-level mask.
- **Errors inside ImageJ or a plugin.** Under Jython the application's error handling caught Python errors only, so a failure raised by ImageJ, Bio-Formats, or a plugin skipped the error dialog. Both kinds are now caught and reported, and a batch continues with the next file.
- **A blank or very dim channel.** When the threshold lay entirely above (or below) the pixel values of an image, for example a negative control thresholded at a fixed Manual value, the whole image failed with "Threshold upper value is below the lower value" and none of its channels or ROIs were written; in a batch the file was reported as FAILED. Such a channel now gives an empty mask: positive area 0 and object count 0, with a note in the run log. Images that were analyzed successfully before give the same numbers.
- **Batch without a display.** `cellquant_batch_headless.py` run with `--headless` failed on every image that counts objects, because ImageJ's ROI Manager cannot be created without a display. Object outlines are collected through an image overlay in that case. A test compares both routes outline by outline and measurement by measurement; with a display the previous route is still used.
- **Channel settings erased by an image with fewer channels.** Opening a 1- to 3-channel image, or clicking "Use active image" on a prepared window, switched off every higher channel and set it to Ignore; the settings did not come back with the next 4-channel image, and a profile saved afterwards stored the erased values. The channel tabs are still greyed out, but the settings are kept and are simply not used for an image that lacks the channel.
- **Older profiles loaded in the window.** A key missing from a profile (for example the binary neighbor count, which 2.2.x did not save) kept whatever the window showed before, while the headless batch used the default. Both now use the default, so one profile gives one set of settings.
- **Validation, sensitivity, and "Capture display minimum" with a Display-minimum threshold.** If the prepared image was missing or out of date, these silently rebuilt it and used its default display range instead of the minimum chosen in Brightness/Contrast. They now stop with the same message as "Run analysis" and ask for the channels to be prepared first.
- A particle analysis started by hand after a run no longer adds its particles to the ROI Manager unexpectedly.
- The version written to every table and settings file now matches the package (2.2.2 packages reported 2.2.0). The launcher, engine, and window refuse to run together if their versions differ.
- A failed run no longer leaves an empty run folder behind when a display-minimum threshold is rejected.

### Added

- **Folder batch processing** (Batch tab and `cellquant_batch_headless.py`): whole image, an ROI file per image, or one ROI Manager set for all images; combined tables written image by image; a per-file report; background execution with progress and cancel; optional Bio-Formats opening and multi-series support.
- **Validation and QC tab**: object detection scored against manual point marks (precision, recall, F1, count error, color-coded overlay); positive-area mask compared with a reference mask (Dice, Jaccard); parameter sensitivity sweeps.
- **Binary neighbor count** per channel (1-8), stored in the profile.
- New output columns: `Touches_ROI_Edge`, `Object_MinFeret`, `Object_Aspect_Ratio`, `Object_Roundness`, `Object_Solidity` in `objects.csv`; edge-touching and excluded object counts in `roi_summary.csv`; neighbor count, pad-edges state, 8-bit rescaling range, and ImageJ version in `thresholds.csv`.
- **Analysis engine separated from the window** (`cellquantengine.py`), usable from scripts and tests without any user interface.
- **Automated tests** (`tests/run_tests.py`) on synthetic images with known answers.
- An example four-channel profile, a validation guide, and citation metadata.

### Changed

- The application is now three files: launcher, engine, and window. This keeps each source file well inside the size Fiji's Jython can compile.
- The automatic threshold is explicitly computed from the full data range. (On ImageJ 1.54 this was already the case; the change guards older versions.)
- ROI sets are read and written directly instead of through a hidden ROI Manager, so a damaged or unwritable file raises an error message instead of opening a dialog mid-run.

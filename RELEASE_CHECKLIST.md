# Before publishing 2.3.0

This file is a working note for the maintainer. Delete it once the steps are done.

## What has and has not been tested

The engine tests and an end-to-end exercise of the window (open image, prepare, preview, analyze with a polygon ROI, profile save and load, validation, sensitivity, batch with cancel and error cases, headless batch) were run against **ImageJ 1.54u** with the Auto Local Threshold and CLAHE plugins, all built from their official sources.

Not tested, because Fiji and Jython could not be installed in the build environment:

- **Fiji's Jython 2.7 interpreter.** The code was run under CPython through a Java bridge, with stand-ins for the two Jython features the window relies on (Python classes extending Java classes, and listener keyword arguments such as `actionPerformed=`). Syntax was checked against the Python 2.7 grammar. Jython-specific behavior is therefore the main thing to confirm.
  Two Jython behaviors that the bridge cannot show were checked by reading the Jython 2.7.3 source instead: how `except` matches Java exceptions (which led to a fix), and the per-file size limit of its compiler (each application file is estimated at about a quarter of the limit, the same as the 2.2 core).
- **Bio-Formats.** Opening CZI and other proprietary formats, and the "every series" option, were not exercised.
- **A real microscope image.** Only synthetic images were used.

## Checks to run in your Fiji

1. Run `tests/run_tests.py` from the Script Editor (Python/Jython). Every test should report PASS; plugin-dependent tests report SKIPPED only if a plugin is missing.
2. Run `Cell_Quant_Studio_final.py`. The window should open with six tabs.
3. Open one of your CZI files, load your usual profile, prepare, preview, and run with a **rectangular** ROI. Compare `roi_summary.csv` with a 2.2.x run of the same image, ROI, and settings: the numbers should be identical (unless the profile uses open or close iterations above 1).
4. Repeat with a **freehand** ROI. The object count should now be lower than in 2.2.x wherever objects lie in the bounding box but outside the ROI; the positive area and intensities should be identical. Check the overlay: no outlined object should lie outside the ROI.
5. Run a second image without reloading ROI Manager. In 2.2.x this failed with "ROI Manager contains no parent ROIs".
6. Run a small batch on a folder of CZI files with the whole-image policy, then with ROI files.
7. Run one validation with a few multi-point marks.

If step 1 fails with a syntax or import error, note the exact message: it will point at a Jython difference that the build environment could not reveal.

## Publishing steps

- Choose a license (MIT and BSD-3-Clause are common for Fiji scripts), add a `LICENSE` file, and add the `license:` line to `CITATION.cff`.
- Replace the legacy macro files in the repository with this package, or move them to a `legacy/` folder.
- Add your own demo profile and a small shareable test image to `examples/`.
- Commit, tag `v2.3.0`, and archive the tag (for example on Zenodo) to obtain a DOI; add it to `CITATION.cff` and the README.

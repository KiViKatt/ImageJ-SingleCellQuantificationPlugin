# Cell Quant Studio Universal

**Configurable multichannel fluorescence microscopy analysis in Fiji/ImageJ.**

Cell Quant Studio combines image preparation, independent channel segmentation, and ROI-based quantification in a graphical application. Researchers configure channel assignments, processing steps, thresholds, and analysis modes through the interface, save those settings as reusable profiles, and apply the workflow to single images or to whole folders.

The application supports nuclear-object counting, marker-positive area measurements, and fluorescence intensity analysis. DAPI and markers such as MAP2 are analyzed independently, so a nuclear segmentation does not have to define the area measured in a non-nuclear marker channel.

> **Version scope:** This documentation describes Universal **2.3.0**, checked against the packaged source. The earlier `CellQuant_01_Prep.ijm` and `CellQuant_02_CellBasedQuant.ijm` macros implement a separate four-channel workflow and are not covered here. If you are upgrading from 2.2.x, read [CHANGELOG.md](CHANGELOG.md) first: some results change, on purpose.

## Capabilities

| Area | Available workflow |
|---|---|
| Channel configuration | Assign marker names and choose a separate analysis mode for each of up to eight channels. |
| Image preparation | Prepare channels using Z-projection or a selected Z slice; optionally subtract background or apply median/Gaussian filtering. |
| Segmentation | Choose automatic global thresholds, manual intensity thresholds, a captured display minimum, or local Phansalkar thresholding. |
| Local contrast | Apply optional CLAHE to the local-threshold segmentation copy only. |
| Object separation | Apply configurable binary cleanup, hole filling, and watershed separation. |
| Quantification | Measure segmented objects, marker-positive area, and fluorescence intensity within analysis ROIs of any shape. |
| Object filtering | Set particle-area and circularity limits and exclude objects on the ROI outline. |
| Batch processing | Analyze every image of a folder with one configuration, from the window or without it. |
| Validation and QC | Score object detection against manual point marks, compare a mask with a reference mask, and measure parameter sensitivity. |
| Export | Save measurement tables and the exact settings used, with optional prepared images, masks, ROIs, and overlays. |

## Analysis modes

Choose a mode according to the biological endpoint and the appearance of the marker.

| Mode | Purpose | Example |
|---|---|---|
| Ignore | Exclude the channel from analysis. | An unused imaging channel. |
| Intensity only | Measure fluorescence within the analysis ROI without identifying individual objects. | Regional Caspase-3 signal. |
| Positive area | Measure the area occupied by threshold-positive signal. | MAP2-positive tissue or neurite coverage. |
| Objects only | Detect and measure segmented objects after applying the selected particle filters. | Discrete SOX2-stained objects. |
| Objects + area | Combine object analysis with positive-area measurements. | DAPI nuclear-object count and DAPI-positive area. |

An object is a segmented image component. Its biological identity must be established for the experiment. A connected MAP2 component can contain processes from multiple neurons; it is not automatically a single neuron.

## Requirements

- Fiji with ImageJ 1.53 or later and Python/Jython scripting support.
- Image files readable by the installed Fiji environment. Bio-Formats is used for formats such as CZI.
- The **Auto Local Threshold** plugin for local Phansalkar thresholding (included in Fiji).
- The **Enhance Local Contrast (CLAHE)** plugin when CLAHE is enabled (included in Fiji).
- Correct channel, Z, time, and spatial-calibration metadata for the images being analyzed.

## Package contents

| File | Purpose |
|---|---|
| `Cell_Quant_Studio_final.py` | Launcher. Run this file. |
| `cellquantengine.py` | Analysis engine: preparation, segmentation, measurement, batch, validation. No user interface. |
| `cellquantcore.py` | The application window. |
| `cellquant_batch_headless.py` | Runs a folder batch from a saved profile without the window. |
| `tests/run_tests.py` | Automated tests of the engine on synthetic images. |
| `examples/` | An example four-channel settings profile to adapt. |
| `docs/VALIDATION.md` | How to produce validation data with the built-in tools. |
| `CHANGELOG.md`, `CITATION.cff` | Release notes and citation metadata. |

The three application files check that they come from the same version and stop with a clear message if they do not.

## Installation and launch

1. Extract the release and keep its folder structure. A suitable location is `Fiji.app/plugins/Scripts/Cell Quant Studio/`.
2. Keep `Cell_Quant_Studio_final.py`, `cellquantengine.py`, and `cellquantcore.py` together.
3. Open `Cell_Quant_Studio_final.py` in Fiji's Script Editor.
4. Select **Language > Python (Jython)** if necessary and click **Run**.

The launcher loads the engine and then the window. The application is split into several files because Fiji's Jython has a size limit per source file; keeping the files separate keeps each one well inside it.

## Quick start

1. Open a representative microscopy image in Fiji and confirm its channel order, Z dimensions, and calibration.
2. Launch Cell Quant Studio and select the input image (drag and drop, **Browse**, or **Use active image**).
3. Assign channel names and analysis modes on the **Channels** tab, or load an existing settings profile.
4. Select the projection or slice preparation for each enabled channel.
5. Configure optional background correction, filtering, thresholding, and object separation.
6. Define the analysis regions on **ROIs and output**: the whole image, the current image selection, or parent area ROIs from ROI Manager.
7. Click **1. Prepare enabled channels**, then **2. Preview all masks**, and inspect the segmentation before processing an experiment.
8. Choose the output folder and export options, then click **3. Run analysis**.
9. Review the measurement tables, masks, and settings together.

Save the settings profile after checking it on representative images. A reusable profile records a procedure; it does not establish that the segmentation is valid for every acquisition or biological sample.

## Batch processing

The **Batch** tab applies the current channel settings, calibration choice, and saved-output choices to every image in a folder.

1. Choose the input folder and the file extensions to include. Tick **Include subfolders** if needed.
2. Choose the regions to quantify:
   - **Whole image**.
   - **ROI file per image**: a ROI Manager set saved as `<image name>_ROIs.zip` (also accepted: `<image name>.zip`, `<image name>_RoiSet.zip`, `<image name>.roi`), beside the image or in a separate ROI folder. Images without a file are either skipped or analyzed whole, as selected. **Save ROI Manager set for the current image** writes the file with the right name.
   - **Current ROI Manager set for every image**.
3. Click **Scan folder** to see which files would be analyzed and which have an ROI file.
4. Click **Run batch**. The batch runs in the background with a progress bar and can be cancelled; images already finished are kept.

Results are written to a new `CellQuant_Batch_<timestamp>` folder inside the output folder chosen on the Session tab (or inside the input folder if none is set):

| File | Content |
|---|---|
| One run folder per image | Identical to a single-image run. |
| `batch_roi_summary.csv`, `batch_objects.csv`, `batch_thresholds.csv` | The tables of all images combined, with `Batch_ID` and `Image_Relative_Path` columns. |
| `batch_report.csv` | One line per file: `OK`, `SKIPPED`, or `FAILED`, with the reason, the ROI source, and the time taken. |
| `settings.properties`, `batch_log.txt` | The configuration and the log of the batch. |

Things to know:

- A **Display minimum** threshold is chosen by eye on one image, so a batch refuses it. Click **Capture display minimum** to turn it into a fixed Manual threshold, or use Auto or Local Phansalkar.
- With **Current frame only**, a batch analyzes the first time point of each image. **All frames** analyzes every time point.
- A file that cannot be opened, or that lacks an enabled channel, is reported as `FAILED` and the batch continues.
- Result folders written by the application are never picked up as input, so the output folder can be the input folder.
- While a batch runs, analysis actions in the window are disabled. Avoid running Analyze Particles or binary commands by hand in Fiji until it finishes.
- Multi-series files: by default the first series is analyzed. **Analyze every series** processes each one as its own image.

### Batch without the window

`cellquant_batch_headless.py` runs the same batch from a saved profile. Open it in the Script Editor and click Run to get a dialog, or start it without a display, for example:

```text
ImageJ-linux64 --ij2 --headless --console --run cellquant_batch_headless.py \
  'input_dir="/data/images",profile="/data/profile.properties",output_dir="/data/results",roi_policy="Whole image",extensions="tif, czi",recursive=false,force_bioformats=false,all_series=false'
```

The launcher executable name depends on the Fiji build. Keep the script beside `cellquantengine.py`.

## Validation and QC

The **Validation and QC** tab tests one channel with its current settings, on the current frame, inside the parent ROIs chosen on **ROIs and output**. See [docs/VALIDATION.md](docs/VALIDATION.md) for a complete procedure.

**Object detection against manual marks.** Mark one point per object with Fiji's Multi-point tool on the source image (or load saved point ROIs), then click **Run validation**. Inside each parent ROI:

| Outcome | Meaning |
|---|---|
| True positive | A detected object containing at least one mark. |
| False positive | A detected object containing no mark. |
| Missed | A mark outside every detected object. |
| Merged | Extra marks inside one object: the detection joined objects counted separately by the annotator. |

Precision, recall, F1, and count error are reported per ROI and pooled, with per-object and per-mark tables and an overlay (green = matched, orange = merged, red = false positive, magenta = missed mark).

**Positive area against a reference mask.** Load a binary reference image of the same size (nonzero = positive) to get Dice, Jaccard, pixel precision and recall, and the area difference for the channel's positive-area mask.

**Parameter sensitivity.** Choose a setting (threshold level, local radius, Phansalkar parameter 1, CLAHE slope, or minimum object area) and a list of factors. The channel is re-segmented with the setting multiplied by each factor, and positive area and object count are reported with their percent change from the configured value.

## Example: SOX2, Caspase-3, MAP2, and DAPI

| Channel | Marker | Analysis mode |
|---|---|---|
| 1 | SOX2 | Objects only |
| 2 | Caspase-3 | Intensity only |
| 3 | MAP2 | Positive area |
| 4 | DAPI | Objects + area |

`examples/Example_4Channel_SOX2_Casp3_MAP2_DAPI.properties` holds these assignments with automatic and local thresholds as a starting point. Load it with **ROIs and output > Profiles and output access > Load settings profile**, then check every threshold and filter on your own images. The assignments are an example, not a requirement: confirm the channel order for each dataset. Profiles saved by 2.2.x load unchanged.

### DAPI segmentation

The DAPI workflow can use local Phansalkar thresholding with optional CLAHE, followed by configurable binary cleanup and watershed separation. The local workflow uses an 8-bit segmentation copy. Pre-watershed and watershed masks are saved for review when **Save raw and final masks** is on.

The pre-watershed ("raw") mask already includes the configured morphology and hole filling. It should not be described as an untouched threshold result.

Review both missed dim nuclei and artificial splits of bright nuclei. Watershed separates connected regions in a binary image; its output still needs to be checked against the source channel.

### Independent MAP2 analysis

Use **Positive area** when the endpoint is MAP2 coverage within a tissue ROI. MAP2 is segmented from its own channel rather than restricted to DAPI nuclear objects. This retains marker signal in non-nuclear structures.

Use an object mode only when connected-component measurements answer the intended question. MAP2 object count, MAP2-positive area, and neuron count are different endpoints.

## Brightness, contrast, and measurement intensities

Brightness/contrast display limits control how an image is shown. A marker display range of **1200-2000** is not automatically a segmentation threshold. Only the **Display minimum** threshold mode uses the lower display limit, and the value used is recorded in `thresholds.csv`. Automatic thresholds are computed from the full data range of the image and do not depend on the display range.

CLAHE changes pixel intensities to enhance local contrast. It operates on a separate segmentation copy, never on the prepared measurement image. Optional background subtraction, filtering, and Z-projection do change the prepared measurement values. Describe measurements as coming from the configured prepared image rather than claiming they always represent unmodified raw pixels.

Using identical local-threshold settings across images applies the same algorithm and parameters. It does not impose one shared absolute intensity cutoff.

## Interpreting the outputs

| Quantity | Meaning |
|---|---|
| Object count | Number of segmented objects inside the parent ROI that pass the particle filters. |
| Object count density | Object count divided by the area of the parent ROI. |
| Positive area | Area of the threshold-positive ("raw") mask inside the parent ROI. |
| Positive-area fraction | Positive area divided by ROI area, also given as a percentage. |
| Mean intensity | Average pixel value of the prepared image inside the measurement region. |
| Integrated density | Mean multiplied by calibrated area; raw integrated density is the sum of pixel values. |

**Objects and ROI outlines.** Objects are detected inside the true outline of each parent ROI, whatever its shape. An object that crosses the outline is clipped there and still counted; `Touches_ROI_Edge` in `objects.csv` flags it. **Exclude objects touching ROI/image edges** removes every object with a pixel on the ROI outline or the image border. `roi_summary.csv` reports both how many edge-touching objects were found and how many were excluded.

Positive-area measurements and the sum of object areas can differ because watershed lines, particle filters, and edge exclusion change which pixels belong to counted objects.

Report area and count-density units from the actual image calibration. `Object_Density_Per_mm2` is filled only when the calibration unit is micrometers. Measurements on projected or selected-slice images describe two-dimensional regions; they are not three-dimensional cell densities.

## Exported results

Each run creates a timestamped folder:

| File or output group | Purpose |
|---|---|
| `roi_summary.csv` | One row per parent ROI and frame: intensity, positive area, object counts and densities for every enabled channel. |
| `objects.csv` | One row per detected object: geometry, shape descriptors, edge flag, and intensity. |
| `thresholds.csv` | The threshold values and every segmentation setting actually used, per channel and frame, with the application and ImageJ versions. |
| `settings.properties` | The complete configuration; it can be loaded again as a profile. |
| `run_log.txt` | Processing messages. |
| `Prepared/`, `Masks/`, `ROIs/`, `Overlays/` | Optional prepared images, raw and final masks, parent and object ROI sets, and object-outline overlays. |

Tables are UTF-8 text without a byte-order mark. In Excel, import them through **Data > From Text/CSV** and choose UTF-8 so that the micrometer sign displays correctly. Columns introduced in 2.3 are placed after the 2.2 columns.

Keep each run's tables, settings, and spatial outputs together.

## Reproducibility

A saved profile is intended to give the same masks on any Fiji installation. To that end:

- Binary open, close, erode, and dilate take their iteration and neighbor counts from the profile, not from **Process > Binary > Options**. Open with N iterations is N erosions followed by N dilations.
- The 8-bit conversion before local thresholding rescales from the image minimum to its maximum regardless of **Edit > Options > Conversions**.
- CLAHE is called through its library with explicit parameters (precise method, no mask).
- Fiji's black-background and pad-edges options are set as needed during processing and restored afterwards.

The automated tests check these points. To run them, open `tests/run_tests.py` in the Script Editor, choose Python (Jython), and click Run. They build synthetic images with known answers, open no window, and print a pass/fail summary.

## Scope and limitations

- The workflow analyzes two-dimensional projections or selected slices. Projection can merge structures that are separated along Z.
- Thresholds, CLAHE settings, watershed, morphology, and particle filters alter object counts and apparent sizes. Inspect their effects on representative images and use the sensitivity tool to quantify them.
- DAPI identifies nuclear signal. A detected DAPI object is not automatically an identified neuron or a validated one-cell segmentation.
- MAP2-positive coverage and connected-object measurements do not directly establish individual-neuron size or neurite length.
- Objects from the same image or organoid share a biological sample. Statistical comparisons should preserve the experiment's replicate structure.
- Segmentation is threshold-based. Densely packed or touching nuclei that watershed cannot separate are better served by learned segmentation tools; the validation tab measures how well the threshold approach works on your images.
- Group statistics, plotting, and reporting are outside this application.

## Troubleshooting

**The launcher cannot find the application files**

Re-extract the release and keep `Cell_Quant_Studio_final.py`, `cellquantengine.py`, and `cellquantcore.py` together. Launch the first one through Fiji's Script Editor.

**A message says the files are from different versions**

A file from an older package is mixed in. Replace all three application files with those of one release.

**An encoding or parser error appears before the interface opens**

Run the launcher, not `cellquantcore.py`. Record the exact error and the Fiji, Java, and package versions when reporting the issue.

**"ROI Manager contains no parent ROIs"**

Open ROI Manager and add at least one area ROI, or switch the quantification region to Whole image.

**DAPI nuclei merge or dim nuclei are missed**

Inspect the source channel and both mask stages. Review projection choice, threshold method, local radius, CLAHE settings, cleanup, and particle-area filters. Change settings on a representative validation subset before applying them across a group comparison.

**MAP2 results appear restricted to nuclei**

Confirm that MAP2 is configured as its own channel with **Positive area**, and that cross-channel object measurement is not being read as total MAP2 coverage.

**A saved preset gives unexpected results**

Check channel assignments, calibration, threshold units, projection settings, and analysis modes. Presets from earlier workflows can use different channel assignments. See the changelog for results that differ between versions.

**A batch reports FAILED files**

Open `batch_report.csv`; the `Message` column gives the reason for each file.

## Reporting an issue

Include the package version, Fiji and Java versions, operating system, image dimensions, channel/Z/time layout, bit depth, calibration, the `settings.properties` of the run, and the complete error message. Where possible, provide a small representative image that may be shared.

## Citation

Citation metadata is in [CITATION.cff](CITATION.cff). Cite the version used, which is recorded in every output table.

## References

- [Project repository](https://github.com/KiViKatt/ImageJ-SingleCellQuantificationPlugin)
- [ImageJ: Jython scripting](https://imagej.net/scripting/jython/)
- [ImageJ: Enhance Local Contrast (CLAHE)](https://imagej.net/plugins/clahe)
- [ImageJ: Auto Local Threshold](https://imagej.net/plugins/auto-local-threshold)

## Acknowledgments

Cell Quant Studio builds on Fiji, ImageJ, Jython, Java Swing, Bio-Formats, and the Auto Local Threshold and CLAHE plugins.

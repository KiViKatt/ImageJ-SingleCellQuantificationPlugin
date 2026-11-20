# Validating a Cell Quant Studio configuration

A settings profile is a procedure, not evidence that the procedure is right for your images. This guide describes how to produce that evidence with the **Validation and QC** tab, and what to report.

## 1. Object counts against manual annotation

Use this for any channel in an Objects mode (for example DAPI nuclei).

**Choose the fields.** Pick fields that span the conditions of the experiment: every experimental group, more than one organoid or sample per group, and both easy and difficult regions (sparse and dense, bright and dim). Five to ten fields with a few hundred objects in total is a common minimum. Choose them before looking at the segmentation.

**Annotate without seeing the result.**

1. Open the image and prepare the channels. Do not preview the mask yet.
2. On the source image, select Fiji's **Multi-point** tool and click once on every object inside the parent ROI. For a Z-projection, annotate what is visible in the projection, because that is what is being segmented.
3. Keep the annotation: press **t** to add the point selection to ROI Manager, select it there, and use **More > Save** to write a `.roi` file. Then delete it from ROI Manager so that only parent area ROIs remain.

A second annotator marking the same fields shows how much two people differ, which is the yardstick for judging the software. Run the validation once with each annotator's marks and report both, together with the two manual counts per field.

**Score.**

1. On **Validation and QC**, choose the channel.
2. Choose **Multi-point selection on the source image** (with the marks active on the image) or **Point ROI file** (with the saved `.roi`).
3. Click **Run validation**.

The run folder holds `validation_summary.csv` (one row per parent ROI plus a pooled row), `validation_objects.csv`, `validation_points.csv`, the mask, the thresholds used, the complete settings, and an overlay:

| Overlay color | Meaning |
|---|---|
| Green outline | Detected object with one mark (true positive). |
| Orange outline | Detected object with several marks (merged objects). |
| Red outline | Detected object without a mark (false positive). |
| Cyan circle | Mark inside a detected object. |
| Magenta circle | Mark outside every detected object (missed). |

**Report.** Per field and pooled: manual count, detected count, precision, recall, and F1. Across fields, show detected count against manual count (a scatter plot with the identity line, and the difference against the mean) rather than a correlation coefficient alone, since two counts can correlate strongly while one is consistently higher. State the number of fields, objects, samples, and annotators.

**Read the errors.** Many orange outlines mean touching objects are not being separated: review watershed, the local radius, and the maximum object area. Many red outlines mean debris or background is being segmented: review the threshold and the minimum object area. Many magenta marks mean dim objects are missed: review the threshold method, CLAHE, and the projection.

**Do not tune on the validation set.** If settings are adjusted after looking at validation results, score the final settings on fields that were not used for tuning.

## 2. Positive area against a reference mask

Use this for a Positive-area channel (for example MAP2 coverage).

1. Make a reference mask of the same pixel size as the image: a binary image in which every nonzero pixel is positive. It can be drawn by hand on the prepared image (for example with the brush tool, then **Edit > Selection > Create Mask**) or produced by an expert-set threshold on that image. Save it as TIFF.
2. On **Validation and QC**, choose the channel, choose **Reference binary mask image (file)**, select the file, and click **Run validation**.

`mask_validation_summary.csv` reports, per parent ROI, the positive pixels of both masks, their overlap, Dice, Jaccard, pixel precision and recall, the two area fractions, and the area difference in percent. For a coverage endpoint the area difference is the most direct figure; Dice and Jaccard show whether the same pixels are being selected.

Hand-drawn masks of fine processes are themselves variable. Where possible, have two people produce a reference and report their agreement with each other alongside the software's agreement with each.

## 3. Sensitivity to settings

A result that changes sharply when a threshold moves slightly is fragile, whatever its validation score.

1. On **Validation and QC**, choose the channel and the parameter to vary.
2. Enter the factors, for example `0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3`. Factor 1.0 is the configured value.
3. Click **Run sensitivity analysis**.

`sensitivity.csv` holds positive area, positive-area fraction, object count, and density for every factor and parent ROI, with the percent change from the configured value. Plot the measured quantity against the parameter value for several fields. A flat region around the configured value supports the choice; a steep slope means the reported numbers depend on it and that dependence should be stated.

Which parameter applies depends on the threshold mode:

| Parameter | Applies to |
|---|---|
| Threshold level | Auto, Manual, and Display minimum thresholds. |
| Local threshold radius, Phansalkar parameter 1 | Local Phansalkar thresholds. |
| CLAHE maximum slope | Channels with segmentation-only CLAHE. |
| Minimum object area | Objects modes. |

## 4. Manual and display-minimum thresholds

A threshold chosen by eye for each image cannot be validated as a procedure, and it can differ between people and between groups. If a channel uses **Display minimum**:

- Click **Capture display minimum** to record the value as a Manual threshold, and apply one value (or one rule) to all images of an experiment; or
- have two people set the threshold independently on the same images, blind to group, and report how far the resulting measurements differ.

`thresholds.csv` records the value actually used for every image either way.

## 5. Going further with the optional modules

- **Tuning with held-out scoring.** Section 1 warns against tuning on the validation set. The tuning module (Analyses tab, **Run tuning...**) does it properly: it reads a folder of images with `<image name>_points.zip` (or `_points.roi`) mark files, searches the chosen settings, and reports the F1 on fields that were not used for tuning. Quote that held-out F1, not the optimistic score on all fields.
- **Corrected counts.** The `validation_summary.csv` files of several validation runs are calibration data for the count-correction module (**Run count correction...**). It fits manual against detected counts and gives corrected counts with 95% prediction intervals, plus a check of whether correction actually reduces the error.
- **Robustness of a group difference.** The robustness module (**Run robustness...**) reports whether a group difference keeps its direction across a grid of defensible settings.
- **Figures.** `tools/cellquant_figures.py validation <folder>` draws the detected-against-manual and Bland-Altman plots from your validation runs.
- **A full study for publication** is planned in [VALIDATION_STUDY.md](VALIDATION_STUDY.md).

## 6. What goes into a methods section

- The application version and the ImageJ version (both are in `thresholds.csv`).
- Per channel: projection, preprocessing, threshold mode and parameters, binary cleanup, watershed, and particle filters. `settings.properties` holds all of them and can be deposited as supplementary material.
- How parent ROIs were defined and by whom.
- The validation figures above, and the sensitivity of the main endpoints.
- The unit of replication used for statistics (for example organoid or batch, not individual objects).

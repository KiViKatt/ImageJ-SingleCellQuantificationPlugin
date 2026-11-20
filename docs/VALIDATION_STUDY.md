# Validation study for the paper

`VALIDATION.md` explains how to validate one settings profile. This file is the plan for the study a reviewer will expect before the software can be published. It is sized to be at or above the published precedents in this niche (Quanty-cFOS 2023, Courtney et al. 2021, AxoNet 2020; see `paper/references.md`) without needing a large team.

Everything here uses tools already in the application. Steps that need a planned module are marked with the module.

## 1. Images

| Item | Minimum | Better |
|---|---|---|
| Fields | 40 | 60 |
| Independent differentiation batches | 3 | 4 or more |
| Conditions or time points | 2 (for example H9 and mutant) | 3 |
| Markers | DAPI (nuclei), one nuclear marker (for example SOX2, PAX6 or TBR2), one non-nuclear area marker (for example MAP2) | add one dim or sparse marker |
| Annotated objects | about 2,000 in total | 4,000 |

- Pick fields before looking at any segmentation, stratified by density (sparse, medium, dense) and by image quality.
- Record for every field the batch, organoid, section, condition and acquisition date. Keep this in a table (`fields.csv`). It is needed for the splits and for the statistics.
- Split by **organoid or batch, never by field**: about 60% for tuning and calibration, 40% as an untouched test set. Fields from one organoid must all fall on the same side.

## 2. Annotation

- At least **two independent annotators**, blind to condition and to the segmentation, mark every test field (Multi-point tool, one click per object; see `VALIDATION.md` section 1). Three is better.
- One annotator re-marks at least 10 fields after at least two weeks, for intra-rater agreement.
- File naming for the auto-tune module: `<image name>_points.zip` beside the image, plus the parent-ROI file `<image name>_ROIs.zip`.
- Area marker: 10 to 15 fields with hand-drawn reference masks from two annotators (`VALIDATION.md` section 2).
- **Human against human is the ceiling.** Report it: annotator A scored against annotator B with the same matching rule.

## 3. Methods compared, on the same test fields and the same matching rule

1. Cell Quant Studio with the expert's settings (as used before the study).
2. Cell Quant Studio with settings from the auto-tune module (tuned on the tuning fields only).
3. A fixed global automatic threshold (Otsu, no CLAHE, no watershed) as a naive baseline.
4. StarDist 2D, pretrained `2D_versatile_fluo` model, in Fiji with default probability and overlap thresholds.
5. Cellpose, pretrained `nuclei` (or `cyto3`) model, default settings.
6. Optional: QuPath cell detection with default settings.

Export each comparator's labels or ROIs and score them with Cell Quant Studio's validation (load them as a reference or convert the labels to object ROIs) so every method is scored by the same code. If a deep-learning method wins on nuclei, say so; the case for this tool then rests on non-nuclear markers, provenance and the robustness report.

## 4. Metrics

Detection (per field and pooled; the validation tab already writes these):

- Precision, recall and F1, matching a mark to the object that contains it.
- Report the sensitivity of F1 to the matching rule. Count merged marks separately (they are the typical failure in dense organoid regions).
- Follow the terminology of Metrics Reloaded (Maier-Hein et al. 2024).

Counts per field:

- Bland-Altman bias and 95% limits of agreement, with a check for proportional bias (difference against mean).
- ICC(2,1), absolute agreement, with its 95% CI (Koo and Li 2016).
- Mean absolute error and mean absolute percentage error.
- Deming regression slope and intercept. A Pearson correlation alone is not enough.

Area marker:

- Dice and Jaccard against each annotator, and between the annotators.
- Area-fraction difference.

## 5. One analysis per planned module

| Module | What to show |
|---|---|
| Auto-tune (`cellquanttune.py`) | Held-out F1 of tuned settings against the expert's and the default settings. Learning curve: held-out F1 with 2, 4, 8 and all tuning fields. |
| Corrected counts (`cellquantcorrect.py`) | On the test fields: bias and limits of agreement before and after correction, and the share of manual counts that fall inside the stated 95% prediction intervals (should be close to 95%). |
| Group robustness (`cellquantrobust.py`) | The real H9 against mutant contrast as a specification curve (effect size with CI for every setting combination). A null control: split one group at random into two halves; the curve should straddle zero. |
| Depth profile (`cellquantdepth.py`) | A synthetic image with a known gradient (exact by construction, from the tests), plus a biological sanity check: progenitor markers (SOX2/PAX6) should peak at a different depth than neuronal markers. |
| Threshold spectrum (`cellquantspectrum.py`) | The coverage curves of both groups for MAP2, the reported value with its plateau width, and its relation to mean intensity (the area under the full curve is the normalized mean intensity). |

## 6. Statistics and reporting

- The statistical unit is the organoid (or batch), not the field and never the object. Average fields within an organoid first.
- Deposit with a DOI (for example Zenodo or the BioImage Archive): the test images or crops of them, all annotations, every `settings.properties`, and the result tables.
- Fill in the community checklist for image analysis reporting (Schmied et al. 2024) as a supplementary table.
- State software versions: the Cell Quant Studio version and the ImageJ version are written to every `thresholds.csv`.

## 7. Order of work

1. Choose and list the fields (`fields.csv`), and make the organoid-level split.
2. Annotators mark the fields, blind.
3. Run all methods on the test fields, and auto-tune on the tuning fields only.
4. Compute the metrics; make the Bland-Altman plots and the specification curve.
5. Write the results into `paper/manuscript_draft.md`.

Expected time: annotation is the bulk. At roughly 1 to 2 seconds per mark, 2,000 marks by two annotators is about 2 to 3 hours each, plus the reference masks.

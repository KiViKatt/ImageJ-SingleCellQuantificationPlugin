# Manuscript draft

Working title: **Cell Quant Studio: validated, reproducible quantification of multichannel fluorescence images of brain organoid sections in Fiji**

Format: F1000Research Software Tool Article (Introduction; Methods: Implementation and Operation; Use case; Discussion; availability). It converts with little change to JORS, SoftwareX or eNeuro "Open Source Tools and Methods". References are numbered by their first appearance below; full details are in `references.md`.

Text in **[brackets]** is still to be filled in from the validation study (`docs/VALIDATION_STUDY.md`). Do not submit before every bracket is gone.

---

## Abstract (about 250 words)

Brain organoid sections are usually quantified with thresholds and counting rules chosen by eye, image by image, and these choices are rarely reported in a form someone else can rerun. We present Cell Quant Studio, an open-source application for Fiji [1, 2] that turns those choices into a saved, reusable settings profile and applies it to single images or whole folders. Each channel is prepared and segmented independently, by global, manual or local (Phansalkar [3]) thresholding with optional contrast enhancement, binary cleanup and watershed. The application then measures object counts, marker-positive area and intensity inside parent regions of any shape.

Every run records the exact settings, threshold values and software versions next to its results. Built-in tools score detected objects against manual point annotations (precision, recall, F1), compare positive-area masks with reference masks (Dice), and measure how strongly results depend on each setting. Optional modules add four functions:

- depth profiles from the organoid surface;
- marker coverage across the full range of thresholds [4];
- tuning of segmentation settings against manual annotations, scored on held-out fields;
- a report of whether a group difference holds across a grid of defensible settings, in the manner of a specification-curve analysis [5, 6].

On **[N] fields from [N] organoids in [N] batches**, annotated independently by **[N]** people, detection reached F1 **[x]** against a human-versus-human ceiling of **[x]**, compared with **[x]** for StarDist and **[x]** for Cellpose. **[One sentence on the robustness result for the H9 versus mutant contrast.]** Cell Quant Studio runs in a standard Fiji installation without programming and is tested automatically.

## Introduction

- Organoid sections are quantified widely. Common endpoints are progenitor and neuron counts, marker coverage such as MAP2, and layering from the surface inward.
- The analysis usually depends on choices made by eye: a threshold, a display range, how touching nuclei are split, which region counts. Two people, or one person on two days, can reach different numbers. Choices made this way are seldom reported in rerunnable form [7, 8].
- General tools exist: CellProfiler [9], QuPath [10], ilastik [11], and deep-learning nucleus segmentation (StarDist [12], Cellpose [13]). Several Fiji macros and plugins count marker-positive cells in sections [14-16].
- What is missing, to our knowledge, is one auditable workflow for organoid sections in Fiji that:
  - keeps nuclear and non-nuclear markers independent;
  - confines counts to hand-drawn or detected organoid outlines;
  - records every setting with the results;
  - includes the means to validate itself and to show how much the conclusions depend on the settings.
- The contribution is therefore not a new segmentation algorithm. It is transparency, validation and robustness checks for a standard threshold-based workflow, packaged for researchers who do not program.

## Methods

### Implementation

- Language and platform: Jython 2.7 inside Fiji (ImageJ **[1.54p]**). No installation beyond copying one folder. Reads proprietary microscope formats through Bio-Formats [17].
- Architecture:
  - an analysis engine with no user interface;
  - a Swing window;
  - optional analysis modules that add columns and tables without changing existing ones;
  - a headless batch script for unattended runs.
- Reproducibility measures (each checked by an automated test):
  - binary operations take their iteration and neighbour counts from the profile, not from Fiji's global options;
  - the 8-bit conversion before local thresholding rescales explicitly;
  - CLAHE [18] is called through its library with explicit parameters;
  - global options touched during processing are restored afterwards;
  - objects are detected inside the true outline of each parent region, not its bounding box.
- Tests: **[N]** automated tests with **[N]** checks run on synthetic images whose correct answers are known by construction. They run in a standard Fiji in about **[N]** seconds.

### Operation

1. Open an image and assign each channel a name and an analysis mode: intensity only, positive area, objects, or objects plus area.
2. Choose preparation (Z projection or slice, background subtraction, filtering) and segmentation settings for each channel.
3. Define parent regions: whole image, ROI Manager outlines, or automatic tissue outlines.
4. Prepare, preview the masks, run. Save the profile and apply it to a folder.
5. Validate on annotated fields; review the sensitivity and robustness reports.

Outputs (UTF-8 CSV):

- per-region summary;
- per-object table;
- thresholds and settings actually used;
- optional masks, overlays and ROI sets;
- module tables (depth profile, threshold spectrum, tuning, calibration, robustness).

### Optional analyses

- **Depth profile.** Euclidean distance from the region outline; pixels and objects binned by depth. The image border is not treated as surface by default.
- **Threshold spectrum.** Coverage at every threshold level, following Johnson and Walker [4]. The area under the full curve equals the normalized mean intensity, so it is reported as such. Summary values are the coverage at the configured threshold, its local sensitivity and the plateau width.
- **Auto-tune.** Grid search of chosen settings for the best pooled F1 against point annotations [19, 20]. Leave-one-field-out held-out scoring; the optimistic score on the tuning fields is labelled as such.
- **Corrected counts.** Regression of manual on detected counts from validated fields, with 95% prediction intervals and a leave-one-out check of whether correction helps, following [21, 22].
- **Group robustness.** Every combination of a settings grid is applied to all images. The experimental unit (organoid) is averaged first. Per combination: effect size (Hedges g), Welch and Mann-Whitney tests. Overall: the share of combinations in which the reference direction holds, reported descriptively [5, 6].

## Use case: **[H9 versus mutant]** organoids

- Data: **[N]** sections of 5-month organoids stained for TBR2, PAX6, HOPX and DAPI, imaged on an Airyscan microscope at about 0.21 µm per pixel, 4 Z slices.
- Validation (table and Bland-Altman figure):
  - Detection: **[precision / recall / F1 per method]**.
  - Count agreement: **[ICC, bias, limits of agreement]**.
- Depth profile: **[figure: progenitor and neuronal markers by depth]**.
- Robustness: **[specification-curve figure; share of combinations with the same sign]**.
- Time: **[minutes per image; minutes for the whole batch]** on a **[CPU, RAM]** workstation.

## Discussion

- Strengths:
  - every number traceable to its settings;
  - validation built in;
  - robustness reported rather than assumed;
  - no programming needed.
- Limitations:
  - threshold-based segmentation is outperformed by deep learning on dense nuclei **[if the comparison shows it]**;
  - measurements are 2-D on projections or single slices, not 3-D densities;
  - corrected counts are only valid for the settings and imaging conditions they were calibrated on;
  - the robustness share is descriptive, not a corrected statistical test.
- Future work: importing StarDist or Cellpose labels as the object source, so the validation and robustness reports can be applied to learned segmentations.

## Data and software availability

- Software: **[GitHub URL]**, archived release **[Zenodo DOI]**, licence GPL-3.0-or-later.
- Data: test images, annotations, settings profiles and result tables at **[repository and DOI]**.

## Use of AI tools

**[Required by JOSS and expected by bioRxiv.]** Describe which parts of the code, tests and text were written with AI assistance, and that the author reviewed and validated them.

## Author contributions, competing interests, funding

**[To fill.]**

---

### Reference numbers used above

1 Schindelin 2012; 2 Schneider 2012; 3 Phansalkar 2011; 4 Johnson and Walker 2015; 5 Steegen 2016; 6 Simonsohn 2020; 7 Botvinik-Nezer 2020; 8 Gomes 2019; 9 Stirling 2021; 10 Bankhead 2017; 11 Berg 2019; 12 Schmidt 2018; 13 Stringer 2021; 14 Beretta 2023; 15 Unger 2025; 16 Waisman 2025; 17 Linkert 2010; 18 Zuiderveld 1994; 19 Held 2013; 20 Teodoro 2017; 21 Reynaud 2012; 22 Ritch 2020.

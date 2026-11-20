# Figures from result tables

`cellquant_figures.py` turns Cell Quant Studio result tables into figures for a paper or a lab meeting. It writes SVG files, which open in a web browser, Inkscape, Illustrator, PowerPoint or Word, and next to each figure a CSV holding exactly the numbers it shows.

It runs with any Python 3 and needs no extra packages. Run it with your normal Python, not with Fiji.

## The four figures

| Command | Reads | Shows |
|---|---|---|
| `validation` | `validation_summary.csv` (one or many validation runs) | Detected against manual count per field with the identity line, and a Bland-Altman plot (difference against mean, with bias and 95% limits of agreement). The pooled row is left out. |
| `depth` | `depth_profile.csv` or `batch_depth_profile.csv` | Per channel: object density, positive-area fraction or mean intensity against depth from the outline. Thin lines: one per image and parent ROI. Thick lines: group means. |
| `spectrum` | `threshold_spectrum.csv` or `batch_threshold_spectrum.csv` | Per channel: coverage at every threshold, one line per image and ROI. The dot marks the threshold the run used. |
| `robustness` | a `CellQuant_Robustness_...` folder | The specification curve: the group difference with its 95% confidence interval for every combination of settings, sorted, with the configured settings marked. Below it, the settings each combination used. |

## Examples

Open a terminal (PowerShell) in the repository folder:

```powershell
python tools\cellquant_figures.py validation "D:\results" -o figures
```

```powershell
python tools\cellquant_figures.py depth "D:\results\CellQuant_Batch_20261119_120000" -o figures --metric Positive_Area_Fraction_Percent
```

```powershell
python tools\cellquant_figures.py depth "D:\results\CellQuant_Batch_20261119_120000" -o figures --relative
```

```powershell
python tools\cellquant_figures.py robustness "D:\results\CellQuant_Robustness_20261119_130000" -o figures
```

A folder is searched with its subfolders. Pointed at a batch folder, the tool uses the combined `batch_...` table and does not count the copies in the run folders again.

## Options

| Option | Meaning |
|---|---|
| `-o FOLDER` | Where the figures go (default `figures`). |
| `--group-regex RULE` | How each image gets its group, from its file name. Default: the start of the name up to the first `-`, `_` or space, so `Control-org1...` is `Control` and `Mutant-org1...` is `Mutant`. Example: `--group-regex "^(Control|Mutant)"`. At most eight groups. |
| `--metric NAME` | Depth: the column to plot, for example `Positive_Area_Fraction_Percent`, `Mean_Intensity` or `Object_Density_Per_mm2`. A channel without that column falls back to its own best measure. Robustness: `g` for Hedges g instead of the difference. |
| `--relative` | Depth: 0 at the outline, 1 at the deepest point, so organoids of different sizes can be compared. |
| `--frame N` | Only this time point. |
| `--title TEXT` | Figure title. |

## What the figures do and do not show

- Group means in the depth figure average the profiles of images and parent ROIs equally. For statistics, average per organoid first; the CSV beside the figure has every value needed.
- The robustness figure's count of combinations that agree in direction is descriptive. The combinations share data and are not independent tests.
- Colours follow a colour-blind-checked palette. Every figure also states its numbers in the CSV, so it can be read without colour.

## Tests

```powershell
python tools\test_cellquant_figures.py
```

The tests build small tables whose plotted numbers are known and check each figure.

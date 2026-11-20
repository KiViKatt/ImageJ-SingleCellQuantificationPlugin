# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach
# Copyright (C) 2026 Maksym Taran

"""Tests of tools/cellquant_figures.py, with small tables whose plotted numbers
are known. Run with any Python 3:

    python tools/test_cellquant_figures.py

The source is kept plain ASCII with no encoding declaration on purpose.
"""

import csv
import math
import os
import shutil
import sys
import tempfile
import unittest
import xml.dom.minidom

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellquant_figures as F  # noqa: E402


def write(path, headers, rows):
    folder = os.path.dirname(path)
    if not os.path.isdir(folder):
        os.makedirs(folder)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, quoting=csv.QUOTE_ALL)
        writer.writerow(headers)
        for row in rows:
            writer.writerow([row.get(h, "") for h in headers])


def read(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class FigureTests(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="cq_fig_")
        self.out = os.path.join(self.tmp, "out")

    def tearDown(self):
        shutil.rmtree(self.tmp, True)

    def assertValidSvg(self, path):
        xml.dom.minidom.parse(path)  # raises when the SVG is not well-formed

    def test_validation_bias_and_limits(self):
        # Differences detected - manual: -10, 0, 10, 20 -> bias 5,
        # sample SD sqrt(((-15)^2 + 5^2 + 5^2 + 15^2) / 3) = sqrt(500/3).
        rows = [
            {"Source_Image": "H9-a.czi", "Parent_ROI_Name": "r1", "Manual_Count": "100", "Detected_Count": "90"},
            {"Source_Image": "H9-b.czi", "Parent_ROI_Name": "r1", "Manual_Count": "200", "Detected_Count": "200"},
            {"Source_Image": "Mut-a.czi", "Parent_ROI_Name": "r1", "Manual_Count": "300", "Detected_Count": "310"},
            {"Source_Image": "Mut-b.czi", "Parent_ROI_Name": "r1", "Manual_Count": "400", "Detected_Count": "420"},
            {"Source_Image": "Mut-b.czi", "Parent_ROI_Name": "ALL_ROIS_POOLED", "Manual_Count": "1000", "Detected_Count": "1020"},
        ]
        headers = ["Source_Image", "Parent_ROI_Name", "Channel_Name", "Manual_Count", "Detected_Count", "F1"]
        write(os.path.join(self.tmp, "runA", "validation_summary.csv"), headers, rows[:2])
        write(os.path.join(self.tmp, "runB", "validation_summary.csv"), headers, rows[2:])
        F.main(["validation", self.tmp, "-o", self.out])
        table = read(os.path.join(self.out, "validation_agreement.csv"))
        self.assertEqual(len(table), 4, "the pooled row is left out")
        self.assertEqual(sorted(set(r["Group"] for r in table)), ["H9", "Mut"])
        diffs = [float(r["Difference"]) for r in table]
        self.assertAlmostEqual(sum(diffs) / 4.0, 5.0)
        self.assertValidSvg(os.path.join(self.out, "validation_agreement.svg"))
        with open(os.path.join(self.out, "validation_agreement.svg"), encoding="utf-8") as handle:
            text = handle.read()
        sd = math.sqrt(500.0 / 3.0)
        self.assertIn("bias 5", text)
        self.assertIn("%.2f" % (5 + 1.96 * sd), text)

    def test_depth_uses_batch_table_once(self):
        headers = ["Source_Image", "Parent_ROI_Name", "Frame", "Channel", "Channel_Name", "Bin_Index",
                   "Depth_Mid", "Depth_Unit", "Relative_Depth_Start", "Relative_Depth_End",
                   "Positive_Area_Fraction_Percent", "Mean_Intensity", "Object_Density_Per_mm2"]
        rows = []
        for image, scale in [("H9-x.czi", 1.0), ("Mut-x.czi", 2.0)]:
            for b in range(3):
                rows.append({"Source_Image": image, "Parent_ROI_Name": "Tissue_001", "Frame": "1", "Channel": "2",
                             "Channel_Name": "DAPI", "Bin_Index": str(b + 1), "Depth_Mid": str(12.5 + 25 * b),
                             "Depth_Unit": "micron", "Relative_Depth_Start": str(b / 3.0),
                             "Relative_Depth_End": str((b + 1) / 3.0), "Mean_Intensity": str(100 * scale * (b + 1)),
                             "Object_Density_Per_mm2": str(1000 * scale)})
        batch = os.path.join(self.tmp, "CellQuant_Batch_1")
        write(os.path.join(batch, "batch_depth_profile.csv"), headers, rows)
        write(os.path.join(batch, "H9-x_CellQuant_1", "depth_profile.csv"), headers, rows[:3])
        F.main(["depth", self.tmp, "-o", self.out])
        table = read(os.path.join(self.out, "depth_profile.csv"))
        self.assertEqual(len(table), 6, "the per-run copy inside a batch is not counted again")
        self.assertEqual(set(r["Measure"] for r in table), set(["Object_Density_Per_mm2"]))
        F.main(["depth", self.tmp, "-o", self.out, "--metric", "Mean_Intensity", "--relative"])
        rel = read(os.path.join(self.out, "depth_profile_relative.csv"))
        self.assertAlmostEqual(float(rel[0]["Depth"]), 1 / 6.0)
        self.assertEqual(rel[2]["Value"], "300.0")
        self.assertValidSvg(os.path.join(self.out, "depth_profile_relative.svg"))

    def test_spectrum_marks_configured_threshold(self):
        headers = ["Source_Image", "Parent_ROI_Name", "Frame", "Channel", "Channel_Name", "Level_Index",
                   "Threshold", "Coverage_Percent"]
        rows = [{"Source_Image": "H9-x.czi", "Parent_ROI_Name": "r", "Frame": "1", "Channel": "3",
                 "Channel_Name": "MAP2", "Level_Index": "0", "Threshold": "500", "Coverage_Percent": "40"}]
        for i, (t, c) in enumerate([(0, 100), (400, 60), (800, 20), (1200, 0)]):
            rows.append(dict(rows[0], Level_Index=str(i + 1), Threshold=str(t), Coverage_Percent=str(c)))
        write(os.path.join(self.tmp, "threshold_spectrum.csv"), headers, rows)
        F.main(["spectrum", self.tmp, "-o", self.out])
        self.assertEqual(len(read(os.path.join(self.out, "threshold_spectrum.csv"))), 5)
        path = os.path.join(self.out, "threshold_spectrum.svg")
        self.assertValidSvg(path)
        with open(path, encoding="utf-8") as handle:
            self.assertIn("threshold 500, coverage 40.00%", handle.read())

    def test_robustness_curve_order(self):
        headers = ["Curve_Rank", "Specification_ID", "Specification", "Is_Reference", "Channel", "Channel_Name",
                   "Metric", "Metric_Unit", "Group_A", "Group_B", "Difference_B_Minus_A", "Difference_CI95_Low",
                   "Difference_CI95_High", "Hedges_g", "Welch_P", "Same_Sign_As_Reference"]
        rows = [
            {"Curve_Rank": "2", "Specification_ID": "S2", "Is_Reference": "True", "Group_A": "H9", "Group_B": "Mut",
             "Difference_B_Minus_A": "5", "Difference_CI95_Low": "1", "Difference_CI95_High": "9",
             "Same_Sign_As_Reference": "True"},
            {"Curve_Rank": "1", "Specification_ID": "S1", "Is_Reference": "False", "Group_A": "H9", "Group_B": "Mut",
             "Difference_B_Minus_A": "-2", "Difference_CI95_Low": "-6", "Difference_CI95_High": "2",
             "Same_Sign_As_Reference": "False"},
        ]
        folder = os.path.join(self.tmp, "CellQuant_Robustness_1")
        write(os.path.join(folder, "robustness_specifications.csv"), headers, rows)
        write(os.path.join(folder, "robustness_specification_parameters.csv"),
              ["Specification_ID", "Parameter", "Value", "Applies"],
              [{"Specification_ID": "S1", "Parameter": "threshold_factor", "Value": "1.2", "Applies": "True"},
               {"Specification_ID": "S2", "Parameter": "threshold_factor", "Value": "1.0", "Applies": "True"}])
        F.main(["robustness", folder, "-o", self.out])
        table = read(os.path.join(self.out, "specification_curve.csv"))
        self.assertEqual([r["Specification_ID"] for r in table], ["S1", "S2"], "sorted by Curve_Rank")
        path = os.path.join(self.out, "specification_curve.svg")
        self.assertValidSvg(path)
        with open(path, encoding="utf-8") as handle:
            text = handle.read()
        self.assertIn("threshold_factor = 1.0", text)
        self.assertIn("configured", text)

    def test_too_many_groups_is_refused(self):
        with self.assertRaises(SystemExit):
            F.group_colours(["g%d" % i for i in range(9)])

    def test_nice_ticks_cover_range(self):
        ticks = F.nice_ticks(0.0, 0.0153)
        self.assertLessEqual(ticks[0], 0.0)
        self.assertGreaterEqual(ticks[-1], 0.0153)


if __name__ == "__main__":
    unittest.main()

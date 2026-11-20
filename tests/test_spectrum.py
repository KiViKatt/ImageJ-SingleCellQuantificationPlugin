# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach

# Cell Quant Studio Universal - tests of the threshold spectrum (cellquantspectrum.py)
#
# Loaded by run_tests.py, which provides the helpers used here.
# The source is kept plain ASCII with no encoding declaration on purpose.
#
# The test image ("plateaus"), 100 x 80 = 8000 pixels, single channel:
#   background 100 everywhere,
#   a 40 x 30 block of 500 at (10, 10)   -> 1200 pixels at 500 or more,
#   inside it a 20 x 10 block of 2000 at (20, 15) -> 200 pixels at 2000.
# So 6800 pixels hold 100, 1000 hold 500, and 200 hold 2000, and the number
# of pixels at or above a threshold t is
#   8000 for t <= 100, 1200 for 100 < t <= 500, 200 for 500 < t <= 2000, 0 above.
# The 8-bit version uses 10, 100, 200 for the three values, and the 32-bit
# version 100.5, 500.5, 2000.5 (all exact in single precision).

from java.util import Random as JavaRandom
from ij.measure import Calibration as SpCalibration
from jarray import array as jarray_array


def sp_single(values):
    """values as a Java float[]: rounded to single precision."""
    return jarray_array(values, "f")

SP_W = 100
SP_H = 80
SP_LEVELS = (100, 500, 2000)
SP_LEVELS_8BIT = (10, 100, 200)
SP_SAVE_KEYS = ["save_prepared", "save_masks", "save_overlays", "save_object_rois", "save_parent_rois"]


def sp_plateau_processor(values=SP_LEVELS, eight_bit=False, add=0, like=None):
    """The plateau image. like, when given, is a blank processor of the
    wanted kind (used for the 32-bit version)."""
    ip = ByteProcessor(SP_W, SP_H) if eight_bit else ShortProcessor(SP_W, SP_H)
    if like is not None:
        ip = like
    ip.setValue(float(values[0] + add))
    ip.fill()
    ip.setValue(float(values[1] + add))
    ip.fill(Roi(10, 10, 40, 30))
    ip.setValue(float(values[2] + add))
    ip.fill(Roi(20, 15, 20, 10))
    ip.resetRoi()
    return ip


def sp_image(title, planes):
    """A single-channel image; several planes become Z slices."""
    stack = ImageStack(planes[0].getWidth(), planes[0].getHeight())
    for plane in planes:
        stack.addSlice("", plane)
    imp = ImagePlus(title, stack)
    if len(planes) > 1:
        imp.setDimensions(1, len(planes), 1)
        imp.setOpenAsHyperStack(True)
    return imp


def sp_config(E, overrides=None):
    values = {
        "enabled": True, "name": "MAP2", "analysis_mode": "Positive area",
        "threshold_mode": "Manual", "manual_low": 300, "manual_high": -1
    }
    if overrides:
        values.update(overrides)
    return E.normalize_channel_config(1, values)


def sp_triangle():
    polygon = Polygon()
    polygon.addPoint(50, 5)
    polygon.addPoint(95, 75)
    polygon.addPoint(5, 75)
    roi = PolygonRoi(polygon, Roi.POLYGON)
    roi.setName("triangle")
    return roi


def sp_options(E, overrides=None, enabled=True):
    options = E.default_global_options()
    options["spectrum_enabled"] = enabled
    for key in SP_SAVE_KEYS:
        options[key] = False
    if overrides:
        options.update(overrides)
    return options


def sp_run(E, tmp, name, imp, configs, rois, overrides=None, enabled=True):
    options = sp_options(E, overrides, enabled)
    run_dir = os.path.join(tmp, name)
    os.makedirs(run_dir)
    return E.AnalysisEngine().analyze_image(
        imp, "", configs, [E.clone_roi(roi) for roi in rois], E.engine_options_from_global(options),
        run_dir, [1], settings=(options, configs)
    )


def sp_rows(result, roi_index, channel=1):
    """(row of the configured threshold or None, rows of the axis levels)."""
    configured = None
    levels = []
    for row in result["extra_tables"]["threshold_spectrum.csv"]["rows"]:
        if row["Parent_ROI_Index"] != roi_index + 1 or row["Channel"] != channel:
            continue
        if row["Level_Index"] == 0:
            configured = row
        else:
            levels.append(row)
    return (configured, levels)


def sp_walk(ip, roi, low, high):
    """Independent answer, pixel by pixel through the ROI's own mask:
    (ROI pixels, pixels with low <= value <= high, mean of the values after
    limiting each to [low, high])."""
    bounds = roi.getBounds()
    mask = roi.getMask()
    pixels = 0
    hits = 0
    total = 0.0
    for y in range(bounds.height):
        for x in range(bounds.width):
            if mask is not None and mask.get(x, y) == 0:
                continue
            value = ip.getf(bounds.x + x, bounds.y + y)
            pixels += 1
            if low <= value <= high:
                hits += 1
            total += min(max(value, low), high)
    return (pixels, hits, total / pixels if pixels > 0 else float("nan"))


def sp_at_least(threshold, counts=(8000, 1200, 200), values=SP_LEVELS):
    if threshold <= values[0]:
        return counts[0]
    if threshold <= values[1]:
        return counts[1]
    if threshold <= values[2]:
        return counts[2]
    return 0


def sp_non_increasing(values):
    for index in range(1, len(values)):
        if values[index] > values[index - 1]:
            return False
    return True


def test_spectrum_known_plateaus(E, t, tmp):
    """Exact pixel counts at every level for a whole-image ROI, a ROI on the
    image border, and a polygon ROI; and the summary values."""
    plane = sp_plateau_processor()
    imp = sp_image("sp_plateaus.tif", [plane])
    rois = [Roi(0, 0, SP_W, SP_H), Roi(0, 0, 30, 25), sp_triangle()]
    # Axis 0, 100, ... 2500: 26 levels, so every level is a whole hundred.
    result = sp_run(E, tmp, "sp_plateaus", imp, [sp_config(E)], rois,
                    {"spectrum_low": "0", "spectrum_high": "2500", "spectrum_steps": "26"})
    imp.close()
    t.equal(result["extra_table_names"], ["threshold_spectrum.csv"], "the spectrum table is reported")
    text = read_utf8(os.path.join(result["run_dir"], "threshold_spectrum.csv"))
    lines = [line for line in text.split("\n") if line.strip() != ""]
    t.equal(lines[0], ",".join(['"' + name + '"' for name in E.SPECTRUM_HEADERS]), "the table header is the fixed list")
    t.equal(len(lines), 1 + 3 * 27, "one row per ROI and level, plus the configured threshold of each ROI")

    # ROI 1, the whole image.
    configured, levels = sp_rows(result, 0)
    t.equal([row["Threshold"] for row in levels], [100.0 * index for index in range(26)], "the axis is evenly spaced from low to high")
    t.equal([row["Level_Index"] for row in levels], list(range(1, 27)), "levels are numbered from 1")
    t.equal([row["Positive_Pixels"] for row in levels], [sp_at_least(100.0 * index) for index in range(26)],
            "whole image: exact pixel count at every level")
    t.check(sp_non_increasing([row["Positive_Pixels"] for row in levels]), "coverage never rises with the threshold")
    t.equal(levels[0]["ROI_Pixels"], 8000, "whole image: ROI pixel count")
    t.close(levels[2]["Coverage_Percent"], 15.0, "whole image: 1200 of 8000 pixels at threshold 200")
    t.equal((levels[0]["Axis_Low"], levels[0]["Axis_High"]), (0.0, 2500.0), "the axis limits are recorded")
    t.equal(levels[0]["Direction"], "above", "a dark-background channel counts pixels at or above the threshold")
    t.equal(levels[0]["Run_ID"], result["run_id"], "rows carry the run identifier")
    t.equal(levels[0]["Application_Version"], E.APP_VERSION, "rows carry the application version")
    t.equal(configured["Threshold"], 300.0, "the configured threshold has its own row")
    t.equal(configured["Positive_Pixels"], 1200, "configured threshold 300: the 500 block and the 2000 block")
    t.check("configured threshold" in configured["Note"], "the configured row is labelled")

    summary = result["summary_rows"][0]
    t.equal(summary["C1_MAP2_Positive_Area"], 1200.0, "the engine's own positive area is 1200 pixels")
    t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], summary["C1_MAP2_Positive_Area_Fraction_Percent"],
            "coverage at the configured threshold equals the positive-area fraction", 1e-12)
    t.equal(summary["C1_MAP2_Spectrum_Threshold"], 300.0, "the configured threshold is in the summary")
    # Within 10 percent of 1200 pixels means 1080 to 1320 pixels: true for
    # every threshold from 101 to 500 and for no other.
    t.equal(summary["C1_MAP2_Spectrum_Stable_From"], 101, "stable range starts just above the background value")
    t.equal(summary["C1_MAP2_Spectrum_Stable_To"], 500, "stable range ends at the block value")
    t.close(summary["C1_MAP2_Spectrum_Stable_Width_Percent"], 100.0 * (500 - 101) / 300.0, "stable width relative to the threshold")
    # No pixel value lies between 297 and 303, so a 1 percent change of the
    # threshold changes nothing.
    t.equal(summary["C1_MAP2_Spectrum_Sensitivity_PP_Per_Percent"], 0.0, "flat spectrum: no sensitivity")
    t.equal(summary["C1_MAP2_Spectrum_Elasticity"], 0.0, "flat spectrum: no elasticity")
    # Mean intensity (6800*100 + 1000*500 + 200*2000) / 8000 = 197.5; on the
    # axis 0..2500 that is 7.9 percent.
    t.close(summary["C1_MAP2_Spectrum_Normalized_Mean_Percent"], 7.9, "normalized mean on the fixed axis")
    t.close(summary["C1_MAP2_Spectrum_Normalized_Mean_Percent"], 100.0 * summary["C1_MAP2_Mean"] / 2500.0,
            "the area under the spectrum is the mean intensity rescaled to the axis", 1e-12)
    t.equal(summary["C1_MAP2_Spectrum_Note"], "", "nothing to note for a plain manual threshold")

    # ROI 2 touches the image border: 30 x 25 = 750 pixels, of which the 500
    # block covers x 10..29, y 10..24 (300) and the 2000 block x 20..29,
    # y 15..24 (100).
    configured, levels = sp_rows(result, 1)
    t.equal(levels[0]["ROI_Pixels"], 750, "border ROI: pixel count")
    t.equal([row["Positive_Pixels"] for row in levels],
            [sp_at_least(100.0 * index, (750, 300, 100)) for index in range(26)], "border ROI: exact pixel count at every level")
    t.equal(configured["Positive_Pixels"], 300, "border ROI: configured threshold")
    t.close(configured["Coverage_Percent"], 40.0, "border ROI: 300 of 750 pixels")
    t.close(result["summary_rows"][1]["C1_MAP2_Spectrum_Coverage_Percent"],
            result["summary_rows"][1]["C1_MAP2_Positive_Area_Fraction_Percent"], "border ROI: equals the positive-area fraction", 1e-12)

    # ROI 3, a triangle: the expected counts come from walking its pixels.
    configured, levels = sp_rows(result, 2)
    walked = [sp_walk(plane, rois[2], 100.0 * index, 1e9) for index in range(26)]
    t.equal(levels[0]["ROI_Pixels"], walked[0][0], "polygon ROI: pixel count")
    t.check(0 < walked[0][0] < 70 * 90, "polygon ROI: fewer pixels than its bounding box")
    t.equal([row["Positive_Pixels"] for row in levels], [item[1] for item in walked], "polygon ROI: exact pixel count at every level")
    t.equal(configured["Positive_Pixels"], sp_walk(plane, rois[2], 300, 1e9)[1], "polygon ROI: configured threshold")
    t.close(result["summary_rows"][2]["C1_MAP2_Spectrum_Coverage_Percent"],
            result["summary_rows"][2]["C1_MAP2_Positive_Area_Fraction_Percent"], "polygon ROI: equals the positive-area fraction", 1e-12)
    t.close(result["summary_rows"][2]["C1_MAP2_Spectrum_Normalized_Mean_Percent"],
            100.0 * sp_walk(plane, rois[2], 0, 2500)[2] / 2500.0, "polygon ROI: normalized mean", 1e-9)


def test_spectrum_ramp_summaries(E, t, tmp):
    """Sensitivity, stable range and the area-under-the-curve identity on a
    ramp, where the density of pixel values is known."""
    # 200 x 50 pixels, value 1000 + 10 * x: each of the 200 values 1000,
    # 1010, ... 2990 is held by exactly 50 pixels.
    ip = ShortProcessor(200, 50)
    for x in range(200):
        ip.setValue(float(1000 + 10 * x))
        ip.fill(Roi(x, 0, 1, 50))
    ip.resetRoi()
    imp = sp_image("sp_ramp.tif", [ip])
    # No axis limits: the axis runs from the image minimum to its maximum.
    result = sp_run(E, tmp, "sp_ramp", imp, [sp_config(E, {"manual_low": 2000})], [Roi(0, 0, 200, 50)],
                    {"spectrum_steps": "200"})
    imp.close()
    configured, levels = sp_rows(result, 0)
    summary = result["summary_rows"][0]
    t.equal((levels[0]["Axis_Low"], levels[0]["Axis_High"]), (1000.0, 2990.0), "the axis is the image minimum and maximum")
    t.equal(len(levels), 200, "one row per level")
    t.equal([row["Applied_Threshold"] for row in levels], [1000.0 + 10 * index for index in range(200)],
            "16-bit levels are applied as whole values")
    t.equal([row["Positive_Pixels"] for row in levels], [50 * (200 - index) for index in range(200)],
            "ramp: exact pixel count at every level")
    t.close(levels[0]["Coverage_Percent"], 100.0, "at the image minimum everything is positive")
    t.close(levels[199]["Coverage_Percent"], 0.5, "at the image maximum only the brightest column is positive")
    # Threshold 2000 selects columns 100..199: 5000 of 10000 pixels.
    t.equal(configured["Positive_Pixels"], 5000, "ramp: configured threshold")
    t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], 50.0, "ramp: coverage at the configured threshold")
    t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], summary["C1_MAP2_Positive_Area_Fraction_Percent"],
            "ramp: equals the positive-area fraction", 1e-12)
    # 50 pixels per 10 intensity units is 0.05 percent of the ROI per unit.
    # One percent of the threshold is 20 units: 1.0 percentage point, which
    # is 2 percent of the coverage of 50 percent.
    t.close(summary["C1_MAP2_Spectrum_Sensitivity_PP_Per_Percent"], -1.0, "ramp: percentage points per percent of threshold", 1e-12)
    t.close(summary["C1_MAP2_Spectrum_Elasticity"], -2.0, "ramp: elasticity", 1e-12)
    # Within 10 percent of 5000 means 4500 to 5500 pixels, that is 90 to 110
    # columns: thresholds 1891 (first to leave out column 89) to 2100.
    t.equal(summary["C1_MAP2_Spectrum_Stable_From"], 1891, "ramp: start of the stable range")
    t.equal(summary["C1_MAP2_Spectrum_Stable_To"], 2100, "ramp: end of the stable range")
    t.close(summary["C1_MAP2_Spectrum_Stable_Width_Percent"], 10.45, "ramp: stable width")
    # Area under the curve = (mean - min) / (max - min); the mean is 1995.
    t.close(summary["C1_MAP2_Spectrum_Normalized_Mean_Percent"], 50.0, "ramp: normalized mean", 1e-12)
    t.close(summary["C1_MAP2_Spectrum_Normalized_Mean_Percent"],
            100.0 * (summary["C1_MAP2_Mean"] - summary["C1_MAP2_Min"]) / (summary["C1_MAP2_Max"] - summary["C1_MAP2_Min"]),
            "area under the spectrum equals (mean - min) / (max - min)", 1e-12)
    area = 0.0
    for index in range(1, 200):
        area += 0.5 * (levels[index]["Coverage_Percent"] + levels[index - 1]["Coverage_Percent"]) * \
            (levels[index]["Threshold"] - levels[index - 1]["Threshold"])
    t.check(abs(area / 1990.0 - summary["C1_MAP2_Spectrum_Normalized_Mean_Percent"]) < 0.5,
            "the trapezoid area under the exported curve agrees with the normalized mean")

    # A wider tolerance gives a wider stable range: 25 percent means 3750 to
    # 6250 pixels, 75 to 125 columns, thresholds 1741 to 2250.
    imp = sp_image("sp_ramp2.tif", [ip])
    wide = sp_run(E, tmp, "sp_ramp_wide", imp, [sp_config(E, {"manual_low": 2000})], [Roi(0, 0, 200, 50)],
                  {"spectrum_steps": "2", "spectrum_tolerance_percent": "25"})
    imp.close()
    t.equal((wide["summary_rows"][0]["C1_MAP2_Spectrum_Stable_From"], wide["summary_rows"][0]["C1_MAP2_Spectrum_Stable_To"]),
            (1741, 2250), "the tolerance setting widens the stable range")
    t.equal(len(sp_rows(wide, 0)[1]), 2, "two levels is the smallest axis")

    # The sensitivity is the slope of the curve between the two thresholds
    # ImageJ applies for threshold - 1 % and + 1 %. Threshold 2010: 1989.9
    # and 2030.1 are applied as 1990 and 2030, which drops 4 columns (1990 to
    # 2020), 2 points of coverage over 40 units: 0.05 points per unit, the
    # ramp's own slope (50 of 10000 pixels every 10 units). One percent of
    # 2010 is 20.1 units, so -1.005 points per percent; dividing by the typed
    # width of 40.2 would give -1.0 and understate the slope of the line.
    # Coverage at 2010 is 99 columns, 49.5 percent, so the elasticity is
    # 100 * -1.005 / 49.5.
    imp = sp_image("sp_ramp3.tif", [ip])
    odd = sp_run(E, tmp, "sp_ramp_odd", imp, [sp_config(E, {"manual_low": 2010})], [Roi(0, 0, 200, 50)],
                 {"spectrum_steps": "2"})
    imp.close()
    summary = odd["summary_rows"][0]
    t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], 49.5, "ramp, threshold 2010: coverage")
    t.close(summary["C1_MAP2_Spectrum_Sensitivity_PP_Per_Percent"], -1.005, "ramp, threshold 2010: slope of the line times 1 percent", 1e-12)
    t.close(summary["C1_MAP2_Spectrum_Elasticity"], 100.0 * -1.005 / 49.5, "ramp, threshold 2010: elasticity", 1e-12)


def test_spectrum_matches_engine_thresholds(E, t, tmp):
    """The configured row reproduces the engine's mask for automatic
    thresholds, a light background, an inverting lookup table, an upper
    threshold, and fractional thresholds on 8-, 16- and 32-bit images."""
    whole = Roi(0, 0, SP_W, SP_H)
    axis = {"spectrum_low": "0", "spectrum_high": "2500", "spectrum_steps": "26"}
    # A third ROI, a polygon that reaches past the top and right image
    # borders and cuts through both blocks: the ROI mask must be clipped to
    # the image the same way the engine clips it.
    outside = Polygon()
    outside.addPoint(30, -20)
    outside.addPoint(130, 20)
    outside.addPoint(30, 60)
    beyond = PolygonRoi(outside, Roi.POLYGON)

    def run(name, planes, overrides, options=None):
        imp = sp_image(name + ".tif", planes)
        try:
            return sp_run(E, tmp, name, imp, [sp_config(E, overrides)], [whole, sp_triangle(), beyond], options or axis)
        finally:
            imp.close()

    def agrees(result, label):
        t.equal(len(result["summary_rows"]), 3, label + ": three ROIs are measured")
        for index in range(len(result["summary_rows"])):
            summary = result["summary_rows"][index]
            configured = sp_rows(result, index)[0]
            # The image is not calibrated, so the engine's positive area is a pixel count.
            t.equal(float(configured["Positive_Pixels"]), summary["C1_MAP2_Positive_Area"],
                    label + ", ROI " + str(index + 1) + ": the configured row has the engine's positive pixel count")
            t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], summary["C1_MAP2_Positive_Area_Fraction_Percent"],
                    label + ", ROI " + str(index + 1) + ": coverage equals the positive-area fraction", 1e-12)
            t.check("CHECK FAILED" not in summary["C1_MAP2_Spectrum_Note"],
                    label + ", ROI " + str(index + 1) + ": the module's own check against the mask passes")

    # Automatic threshold (Otsu), dark background.
    result = run("sp_otsu", [sp_plateau_processor()], {"threshold_mode": "Auto", "auto_method": "Otsu"})
    agrees(result, "Otsu")
    configured, levels = sp_rows(result, 0)
    t.equal(configured["Threshold"], result["threshold_rows"][0]["Threshold_Low"], "Otsu: the configured row uses the recorded threshold")
    t.equal(configured["Direction"], "above", "Otsu on a dark background counts at or above")
    t.check(configured["Positive_Pixels"] in [1200, 200], "Otsu separates the blocks from the background")

    # Automatic threshold, light background: positive is at or BELOW.
    result = run("sp_light", [sp_plateau_processor()], {"threshold_mode": "Auto", "auto_method": "Otsu", "dark_background": False})
    agrees(result, "light background")
    configured, levels = sp_rows(result, 0)
    t.equal(configured["Direction"], "below", "a light-background channel counts at or below the threshold")
    t.equal(configured["Threshold"], result["threshold_rows"][0]["Threshold_High"], "light background: the moving threshold is the upper one")
    # At or below t: 0 below 100, 6800 from 100, 7800 from 500, 8000 from 2000.
    expected = []
    for index in range(26):
        expected.append(8000 - sp_at_least(100.0 * index + 1))
    t.equal([row["Positive_Pixels"] for row in levels], expected, "light background: exact pixel count at every level")
    t.check(sp_non_increasing([-row["Positive_Pixels"] for row in levels]), "light background: coverage never falls with the threshold")
    t.check(configured["Positive_Pixels"] in [6800, 7800], "light background: the dark side is selected")

    # An inverting lookup table swaps the side ImageJ selects. Checked
    # against the mask build_masks makes from the same image.
    engine = E.AnalysisEngine()
    plane = sp_plateau_processor()
    plane.invertLut()
    prepared = ImagePlus("sp_inverted", plane)
    config = sp_config(E, {"threshold_mode": "Auto", "auto_method": "Otsu"})
    data = engine.build_masks(prepared, config, "sp_inverted")
    plan, measured = E.spectrum_for_roi(prepared, sp_triangle(), config, data, sp_options(E, axis))
    t.equal(plan["direction"], "below", "inverting lookup table with 'dark background': at or below")
    t.equal(measured["count"], E.count_positive_pixels(data["raw"].getProcessor(), sp_triangle()),
            "inverting lookup table: the count is the mask's count")
    t.check(prepared.getRoi() is None and prepared.getProcessor().getMask() is None, "the prepared image is left without a ROI")
    E.safe_close_image(data["raw"])
    E.safe_close_image(data["final"])

    # Manual threshold with an upper limit: 300..1000 selects only the 1000
    # pixels at 500; the 200 pixels at 2000 are left out at every level.
    result = run("sp_band", [sp_plateau_processor()], {"manual_low": 300, "manual_high": 1000})
    agrees(result, "upper limit")
    configured, levels = sp_rows(result, 0)
    t.equal(configured["Positive_Pixels"], 1000, "upper limit: only the 500 block")
    expected = []
    for index in range(26):
        expected.append(max(0, sp_at_least(100.0 * index) - 200))
    t.equal([row["Positive_Pixels"] for row in levels], expected, "upper limit: every level leaves out the pixels above it")
    t.check("upper threshold" in configured["Note"], "upper limit: the table says so")
    # Within 10 percent of 1000 pixels: thresholds 101 to 500.
    t.equal((result["summary_rows"][0]["C1_MAP2_Spectrum_Stable_From"], result["summary_rows"][0]["C1_MAP2_Spectrum_Stable_To"]),
            (101, 500), "upper limit: stable range")

    # Fractional thresholds. ImageJ rounds on 16-bit images: 500.4 keeps the
    # pixels at 500, 500.5 drops them.
    result = run("sp_round_down", [sp_plateau_processor()], {"manual_low": 500.4})
    agrees(result, "16-bit 500.4")
    t.equal(sp_rows(result, 0)[0]["Positive_Pixels"], 1200, "16-bit: 500.4 is applied as 500")
    t.equal(sp_rows(result, 0)[0]["Applied_Threshold"], 500.0, "16-bit: the applied threshold is recorded")
    result = run("sp_round_up", [sp_plateau_processor()], {"manual_low": 500.5})
    agrees(result, "16-bit 500.5")
    t.equal(sp_rows(result, 0)[0]["Positive_Pixels"], 200, "16-bit: 500.5 is applied as 501")

    # 8-bit images: ImageJ cuts the fraction off, so 100.7 keeps the pixels at 100.
    result = run("sp_8bit", [sp_plateau_processor(SP_LEVELS_8BIT, True)], {"manual_low": 100.7}, {"spectrum_steps": "20"})
    agrees(result, "8-bit 100.7")
    configured, levels = sp_rows(result, 0)
    t.equal(configured["Positive_Pixels"], 1200, "8-bit: 100.7 is applied as 100")
    t.equal((levels[0]["Axis_Low"], levels[0]["Axis_High"]), (10.0, 200.0), "8-bit: the axis is the image minimum and maximum")
    # Levels 10, 20, ... 200.
    expected = []
    for index in range(20):
        expected.append(sp_at_least(10.0 + 10 * index, (8000, 1200, 200), SP_LEVELS_8BIT))
    t.equal([row["Positive_Pixels"] for row in levels], expected, "8-bit: exact pixel count at every level")
    t.close(result["summary_rows"][0]["C1_MAP2_Spectrum_Normalized_Mean_Percent"],
            100.0 * (result["summary_rows"][0]["C1_MAP2_Mean"] - 10.0) / 190.0, "8-bit: normalized mean identity", 1e-12)

    # 32-bit: the image holds 100.5, 500.5 and 2000.5. A 32-bit threshold is
    # not rounded to a whole value, so 500.5 keeps the pixels at 500.5 and
    # 500.51 drops them.
    def float_planes():
        return [sp_plateau_processor(SP_LEVELS, False, 0.5, E.FloatProcessor(SP_W, SP_H))]

    planes = float_planes()
    average = {"manual_low": 500.5}
    probe = sp_image("sp_probe", planes)
    prepared = engine.prepare_single_channel(probe, sp_config(E, average), 1)
    t.equal(prepared.getBitDepth(), 32, "a 32-bit image stays 32-bit when it is prepared")
    walked_low = sp_walk(prepared.getProcessor(), sp_triangle(), 500.5, 1e9)
    walked_axis = sp_walk(prepared.getProcessor(), sp_triangle(), 200.0, 1000.0)
    t.check(walked_low[1] > 0, "32-bit polygon ROI: the triangle holds positive pixels")
    prepared.close()
    probe.close()
    # The axis 200..1000 is narrower than the data, so the normalized mean
    # limits 100.5 to 200 and 2000.5 to 1000:
    # (6800*200 + 1000*500.5 + 200*1000) / 8000 = 257.5625 -> 7.1953125 percent.
    narrow = {"spectrum_low": "200", "spectrum_high": "1000", "spectrum_steps": "9"}
    result = run("sp_float", planes, average, narrow)
    agrees(result, "32-bit 500.5")
    configured, levels = sp_rows(result, 0)
    t.equal(configured["Positive_Pixels"], 1200, "32-bit: 500.5 keeps the pixels at 500.5")
    t.equal([row["Positive_Pixels"] for row in levels], [1200, 1200, 1200, 1200, 200, 200, 200, 200, 200],
            "32-bit: exact pixel count at 200, 300, ... 1000")
    t.close(result["summary_rows"][0]["C1_MAP2_Spectrum_Normalized_Mean_Percent"], 7.1953125, "32-bit: normalized mean on a narrow axis", 1e-9)
    # 1080..1320 pixels: thresholds above 100.5 up to 500.5.
    t.equal((result["summary_rows"][0]["C1_MAP2_Spectrum_Stable_From"], result["summary_rows"][0]["C1_MAP2_Spectrum_Stable_To"]),
            (100.5, 500.5), "32-bit: stable range")
    # No pixel value lies within 1 percent of 500.5 except 500.5 itself:
    # 1000 of 8000 pixels (12.5 points) are lost between 495.495 and 505.505,
    # each applied in single precision, so the slope is -12.5 points over
    # their single-precision distance (10.01 up to rounding), and the
    # sensitivity is that slope times 5.005 units (1 percent of 500.5): about
    # -6.25 points per percent. The coverage is 15 percent, so the elasticity
    # is 100 * sensitivity / 15 (about -41.67), defined as on the ramp
    # (-1 point of 50 percent = -2).
    applied = float(sp_single([505.505])[0]) - float(sp_single([495.495])[0])
    sensitivity = -12.5 / applied * 5.005
    t.close(sensitivity, -6.25, "32-bit: the expected sensitivity is -6.25 up to single-precision rounding", 1e-5)
    t.close(result["summary_rows"][0]["C1_MAP2_Spectrum_Sensitivity_PP_Per_Percent"], sensitivity, "32-bit: sensitivity", 1e-12)
    t.close(result["summary_rows"][0]["C1_MAP2_Spectrum_Elasticity"], 100.0 * sensitivity / 15.0, "32-bit: elasticity", 1e-12)
    configured, levels = sp_rows(result, 1)
    t.equal(configured["ROI_Pixels"], walked_low[0], "32-bit polygon ROI: pixel count")
    t.equal(configured["Positive_Pixels"], walked_low[1], "32-bit polygon ROI: configured threshold")
    t.close(result["summary_rows"][1]["C1_MAP2_Spectrum_Normalized_Mean_Percent"], 100.0 * (walked_axis[2] - 200.0) / 800.0,
            "32-bit polygon ROI: normalized mean", 1e-9)
    result = run("sp_float_up", float_planes(), {"manual_low": 500.51}, narrow)
    agrees(result, "32-bit 500.51")
    t.equal(sp_rows(result, 0)[0]["Positive_Pixels"], 200, "32-bit: 500.51 drops the pixels at 500.5")
    # Automatic threshold on a 32-bit image, with no fixed axis: the axis
    # is the image minimum and maximum, 100.5 and 2000.5.
    result = run("sp_float_otsu", float_planes(), {"threshold_mode": "Auto", "auto_method": "Otsu"}, {"spectrum_steps": "20"})
    agrees(result, "32-bit Otsu")
    configured, levels = sp_rows(result, 0)
    t.equal((levels[0]["Axis_Low"], levels[0]["Axis_High"]), (100.5, 2000.5), "32-bit: the axis is the image minimum and maximum")
    t.equal((levels[0]["Positive_Pixels"], levels[19]["Positive_Pixels"]), (8000, 200), "32-bit: both ends of the axis are counted exactly")
    t.close(result["summary_rows"][0]["C1_MAP2_Spectrum_Normalized_Mean_Percent"],
            100.0 * (result["summary_rows"][0]["C1_MAP2_Mean"] - 100.5) / 1900.0, "32-bit: normalized mean identity", 1e-9)


def test_spectrum_degenerate_cases(E, t, tmp):
    """Binary cleanup, a constant image, a one-pixel ROI, a ROI with no
    positive pixel, and an empty distribution."""
    axis = {"spectrum_low": "0", "spectrum_high": "2500", "spectrum_steps": "26"}

    # One isolated bright pixel is removed by a binary opening, so the mask
    # holds 1200 pixels while the plain threshold selects 1201.
    plane = sp_plateau_processor()
    plane.set(90, 70, 2000)
    imp = sp_image("sp_cleanup.tif", [plane])
    result = sp_run(E, tmp, "sp_cleanup", imp, [sp_config(E, {"open_iterations": 1})], [Roi(0, 0, SP_W, SP_H)], axis)
    imp.close()
    summary = result["summary_rows"][0]
    t.equal(summary["C1_MAP2_Positive_Area"], 1200.0, "cleanup: the mask loses the isolated pixel")
    t.equal(sp_rows(result, 0)[0]["Positive_Pixels"], 1201, "cleanup: the spectrum is the plain threshold count")
    t.check("before cleanup" in summary["C1_MAP2_Spectrum_Note"], "cleanup: the summary says the two differ")
    t.check("before cleanup" in sp_rows(result, 0)[1][0]["Note"], "cleanup: the table says the two differ")

    # A constant image has no range: one level, everything positive.
    flat = ShortProcessor(40, 30)
    flat.setValue(300.0)
    flat.fill()
    imp = sp_image("sp_flat.tif", [flat])
    result = sp_run(E, tmp, "sp_flat", imp, [sp_config(E, {"manual_low": 100})], [Roi(0, 0, 40, 30)])
    imp.close()
    configured, levels = sp_rows(result, 0)
    summary = result["summary_rows"][0]
    t.equal(len(levels), 1, "constant image: a single level")
    # build_masks raises a manual threshold below the image minimum to the
    # minimum (it selects the same pixels), so the run's threshold is 300.
    t.equal(summary["C1_MAP2_Spectrum_Threshold"], 300, "constant image: the threshold the run used")
    t.equal((levels[0]["Threshold"], levels[0]["Positive_Pixels"]), (300.0, 1200), "constant image: every pixel is at the level")
    t.close(configured["Coverage_Percent"], 100.0, "constant image: full coverage")
    t.equal(summary["C1_MAP2_Spectrum_Normalized_Mean_Percent"], "", "constant image: no normalized mean")
    t.check("axis is empty" in summary["C1_MAP2_Spectrum_Note"], "constant image: the note explains the single level")
    t.equal((summary["C1_MAP2_Spectrum_Stable_From"], summary["C1_MAP2_Spectrum_Stable_To"]), (300, 300),
            "constant image: the stable range is limited to the one value present")
    t.equal(summary["C1_MAP2_Spectrum_Sensitivity_PP_Per_Percent"], -50.0,
            "constant image: all pixels sit at the threshold (100 points lost between 297 and 303, times 3 units per percent)")

    # A one-pixel ROI on the 2000 block and a ROI on background only.
    imp = sp_image("sp_small.tif", [sp_plateau_processor()])
    result = sp_run(E, tmp, "sp_small", imp, [sp_config(E)], [Roi(25, 20, 1, 1), Roi(60, 50, 30, 20)], axis)
    imp.close()
    configured, levels = sp_rows(result, 0)
    summary = result["summary_rows"][0]
    t.equal((configured["ROI_Pixels"], configured["Positive_Pixels"]), (1, 1), "one-pixel ROI: one positive pixel")
    t.equal([row["Positive_Pixels"] for row in levels], [1] * 21 + [0] * 5, "one-pixel ROI: positive up to its value of 2000")
    # Every threshold up to 2000 selects the one pixel, so the range has no
    # low end. It must still hold the configured threshold 300 (a range of
    # 2000 to 2000 would wrongly say the coverage is fragile), so it is cut
    # at 300, below the ROI's only value, and the note says it is open.
    t.equal((summary["C1_MAP2_Spectrum_Stable_From"], summary["C1_MAP2_Spectrum_Stable_To"]),
            (300, 2000), "one-pixel ROI: stable range from the configured threshold to the value present")
    t.check("open at its low end" in summary["C1_MAP2_Spectrum_Note"], "one-pixel ROI: the stable range is noted as open")
    t.close(summary["C1_MAP2_Spectrum_Stable_Width_Percent"], 100.0 * (2000 - 300) / 300.0, "one-pixel ROI: stable width")
    configured, levels = sp_rows(result, 1)
    summary = result["summary_rows"][1]
    t.equal((configured["ROI_Pixels"], configured["Positive_Pixels"]), (600, 0), "background ROI: nothing positive")
    t.equal(summary["C1_MAP2_Spectrum_Coverage_Percent"], 0.0, "background ROI: zero coverage")
    t.equal(summary["C1_MAP2_Positive_Area_Fraction_Percent"], 0.0, "background ROI: the engine agrees")
    t.equal(summary["C1_MAP2_Spectrum_Elasticity"], "", "background ROI: no elasticity of zero coverage")
    t.equal((summary["C1_MAP2_Spectrum_Stable_From"], summary["C1_MAP2_Spectrum_Stable_To"]), ("", ""),
            "background ROI: no stable range of zero coverage")
    t.equal([row["Positive_Pixels"] for row in levels], [600, 600] + [0] * 24, "background ROI: every pixel holds 100")

    # No pixels at all: nothing to divide by, and no error.
    settings = E.spectrum_settings(sp_options(E, axis))
    probe = ImagePlus("sp_probe", sp_plateau_processor())
    config = sp_config(E)
    data = E.AnalysisEngine().build_masks(probe, config, "sp_probe")
    plan = E.spectrum_plan(probe, config, data, settings)
    E.safe_close_image(data["raw"])
    E.safe_close_image(data["final"])
    empty = E.spectrum_measure(E.SpectrumCounts([0] * 65536), plan, settings["tolerance"])
    t.equal(empty["roi_pixels"], 0, "empty distribution: no pixels")
    t.equal((empty["coverage"], empty["normalized_mean"], empty["stable_from"]), ("", "", ""), "empty distribution: values stay empty")
    t.check("no pixels" in empty["note"], "empty distribution: the note says why")
    rows = E.spectrum_table_rows({}, plan, empty)
    t.equal(len(rows), 26, "empty distribution: the level rows are still written")
    t.equal(rows[0]["Coverage_Percent"], "", "empty distribution: coverage is left empty")

    # The distribution classes against direct counting, including the ends.
    counts = E.SpectrumCounts([3, 0, 2, 5] + [0] * 252)
    t.equal([counts.at_least(v) for v in [-1, 0, 0.5, 1, 2, 2.5, 3, 3.5, 300]], [10, 10, 7, 7, 7, 5, 5, 0, 0], "counts at or above a value")
    t.equal([counts.greater(v) for v in [-1, 0, 0.5, 2, 3, 300]], [10, 7, 7, 5, 0, 0], "counts above a value")
    t.equal([counts.kth_largest(k) for k in [1, 5, 6, 7, 8, 10]], [3, 3, 2, 2, 0, 0], "k-th largest value")
    t.equal([counts.kth_smallest(k) for k in [1, 3, 4, 5, 6, 10]], [0, 0, 2, 2, 3, 3], "k-th smallest value")
    # Values 0,0,0,2,2,3,3,3,3,3 limited to [1, 2.5]: 1,1,1,2,2,2.5 x 5 -> 19.5 / 10.
    t.close(counts.clipped_mean(1.0, 2.5), 1.95, "mean of the limited values", 1e-12)
    t.close(counts.clipped_mean(-2.0, 400.0), 1.9, "limits outside the data leave the mean alone", 1e-12)


def test_spectrum_off_settings_and_batch(E, t, tmp):
    """Switched off nothing changes; the settings are saved and checked; a
    batch on a fixed axis writes one combined table."""
    # Channel 1 counts objects only, so it has no positive area and no
    # spectrum; channel 2 is a positive-area channel.
    configs = [nuclei_config(E, {"analysis_mode": "Objects only"}), marker_config(E)]
    rois = [Roi(0, 0, WIDTH, HEIGHT), triangle_roi()]
    engine = E.AnalysisEngine()
    base_headers = engine.build_summary_headers(configs)
    t.equal(E.default_global_options()["spectrum_enabled"], False, "the spectrum is off by default")
    t.equal(E.ENGINE_MODULE_VERSIONS.get("cellquantspectrum.py"), E.APP_VERSION, "the module reports the application version")
    # "Objects + area" also measures positive area, so it gets a spectrum too,
    # and the columns follow the order of the channels.
    both = [nuclei_config(E), marker_config(E)]
    t.equal(E.ThresholdSpectrumExtension().summary_headers(both, sp_options(E)),
            [E.channel_column_prefix(config) + suffix for config in both for suffix in E.SPECTRUM_SUMMARY_SUFFIXES],
            "an 'Objects + area' channel has spectrum columns, in channel order")

    imp = make_image("sp_switch.tif")
    off = sp_run(E, tmp, "sp_off", imp, configs, rois, None, False)
    on = sp_run(E, tmp, "sp_on", imp, configs, rois)
    imp.close()
    t.equal(off["summary_headers"], base_headers, "switched off: summary columns are unchanged")
    t.equal(off["extra_table_names"], [], "switched off: no extra table")
    t.check(not os.path.exists(os.path.join(off["run_dir"], "threshold_spectrum.csv")), "switched off: no spectrum file")
    added = [E.channel_column_prefix(configs[1]) + suffix for suffix in E.SPECTRUM_SUMMARY_SUFFIXES]
    t.equal(on["summary_headers"], base_headers + added, "switched on: columns are appended for the positive-area channels in order")
    t.equal(on["object_headers"], off["object_headers"], "the object table gets no new column")
    for name in ["objects.csv", "thresholds.csv"]:
        t.equal(read_utf8(os.path.join(on["run_dir"], name)).replace("sp_on", "sp_off"),
                read_utf8(os.path.join(off["run_dir"], name)), name + " is the same with the spectrum on or off")
    off_lines = read_utf8(os.path.join(off["run_dir"], "roi_summary.csv")).split("\n")
    on_lines = read_utf8(os.path.join(on["run_dir"], "roi_summary.csv")).replace("sp_on", "sp_off").split("\n")
    t.equal(len(on_lines), len(off_lines), "roi_summary.csv keeps its rows")
    same = True
    for index in range(len(off_lines)):
        if off_lines[index] != "" and not on_lines[index].startswith(off_lines[index] + ","):
            same = False
    t.check(same, "every roi_summary.csv line with the spectrum on starts with the line written with it off")
    # The standard marker rectangle: 100 x 80 = 8000 pixels at 900 or more.
    t.equal(sp_rows(on, 0, 2)[0]["Positive_Pixels"], 8000, "two-channel image: the marker channel is measured")
    t.equal(len(sp_rows(on, 0, 1)[1]), 0, "a channel without positive area has no spectrum")
    t.close(on["summary_rows"][1]["C2_MAP2_Spectrum_Coverage_Percent"], on["summary_rows"][1]["C2_MAP2_Positive_Area_Fraction_Percent"],
            "calibrated image, polygon ROI: coverage equals the positive-area fraction", 1e-12)
    t.check(not E.ThresholdSpectrumExtension().is_enabled(sp_options(E), [nuclei_config(E, {"analysis_mode": "Objects only"})]),
            "nothing is added when no channel is analysed for positive area")

    # Settings: saved in a profile, defaults for an older profile, and checked.
    options = sp_options(E, {"spectrum_low": "0", "spectrum_high": "4095", "spectrum_steps": "64", "spectrum_tolerance_percent": "5"})
    path = os.path.join(tmp, "sp_profile.properties")
    E.write_settings_properties(path, options, configs)
    loaded = E.load_profile_settings(path)[0]
    for key in ["spectrum_enabled", "spectrum_low", "spectrum_high", "spectrum_steps", "spectrum_tolerance_percent"]:
        t.equal(loaded[key], options[key], "setting '" + key + "' survives a save and load")
    legacy = E.read_properties_file(path)
    for spec in E.SPECTRUM_SETTING_SPECS:
        legacy.remove(spec[0])
    legacy_options = E.global_options_from_properties(legacy)
    t.equal(legacy_options["spectrum_enabled"], False, "a profile without the keys loads with the spectrum off")
    parsed = E.spectrum_settings(legacy_options)
    t.equal((parsed["steps"], parsed["low"], parsed["high"], parsed["tolerance"]), (256, None, None, 10.0), "defaults of the settings")
    for bad in [{"spectrum_steps": "1"}, {"spectrum_steps": "many"}, {"spectrum_low": "500", "spectrum_high": "500"},
                {"spectrum_low": "x"}, {"spectrum_tolerance_percent": "100"}, {"spectrum_tolerance_percent": "-1"}]:
        raised = False
        try:
            E.spectrum_settings(sp_options(E, bad))
        except ValueError:
            raised = True
        t.check(raised, "a setting that cannot be used is refused: " + repr(bad))

    # Batch: two images of different brightness on one fixed axis.
    folder = os.path.join(tmp, "sp_batch")
    os.makedirs(folder)
    for index in range(2):
        imp = sp_image("sp_b" + str(index) + ".tif", [sp_plateau_processor(SP_LEVELS, False, 50 * index)])
        FileSaver(imp).saveAsTiff(os.path.join(folder, "sp_b" + str(index) + ".tif"))
        imp.close()
    batch_configs = [sp_config(E)]

    def batch(overrides, out_name):
        result = E.run_batch(folder, os.path.join(tmp, out_name), batch_configs, batch_configs, sp_options(E, overrides),
                             {"extensions": "tif", "roi_policy": E.BATCH_ROI_WHOLE})
        t.equal(result["analyzed"], 2, out_name + ": the batch analyzed both images")
        path = os.path.join(result["batch_dir"], "batch_threshold_spectrum.csv")
        t.check(os.path.isfile(path), out_name + ": the batch writes the combined spectrum table")
        lines = [line for line in read_utf8(path).split("\n") if line.strip() != ""]
        t.check(lines[0].endswith('"Note","Batch_ID","Image_Relative_Path"'), out_name + ": the combined table carries the batch columns")
        names = lines[0].replace('"', "").split(",")
        return [dict(zip(names, line.replace('"', "").split(","))) for line in lines[1:]]

    rows = batch({"spectrum_low": "0", "spectrum_high": "2500", "spectrum_steps": "26"}, "sp_batch_fixed")
    t.equal(len(rows), 2 * 27, "fixed axis: one row per image and level, plus the configured threshold")
    t.equal(set([(row["Axis_Low"], row["Axis_High"]) for row in rows]), set([("0.0", "2500.0")]), "fixed axis: every image shares the axis")
    first = [row for row in rows if row["Source_Image"] == "sp_b0.tif" and row["Level_Index"] != "0"]
    second = [row for row in rows if row["Source_Image"] == "sp_b1.tif" and row["Level_Index"] != "0"]
    t.equal([row["Threshold"] for row in first], [row["Threshold"] for row in second], "fixed axis: the same levels in both images")
    # The second image is 50 brighter: 150, 550, 2050. At level 100 * i the
    # first image has the standard counts and the second those of t - 50.
    t.equal([int(row["Positive_Pixels"]) for row in first], [sp_at_least(100.0 * index) for index in range(26)], "fixed axis: first image")
    t.equal([int(row["Positive_Pixels"]) for row in second], [sp_at_least(100.0 * index - 50) for index in range(26)], "fixed axis: second image")
    summary_header = read_utf8(os.path.join(os.path.dirname(os.path.join(tmp, "sp_batch_fixed", "x")),
                                            os.listdir(os.path.join(tmp, "sp_batch_fixed"))[0], "batch_roi_summary.csv")).split("\n")[0]
    t.check('"C1_MAP2_Spectrum_Note","Batch_ID","Image_Relative_Path"' in summary_header, "the combined summary carries the spectrum columns")

    rows = batch({"spectrum_steps": "20"}, "sp_batch_own")
    t.equal(set([(row["Source_Image"], row["Axis_Low"], row["Axis_High"]) for row in rows]),
            set([("sp_b0.tif", "100.0", "2000.0"), ("sp_b1.tif", "150.0", "2050.0")]), "no limits: each image has its own axis")


def test_spectrum_local_threshold(E, t, tmp):
    """A Local Phansalkar channel has no single threshold: the global
    spectrum is still exported and the values at the threshold stay empty."""
    if not local_threshold_available(E):
        raise SkipTest("Auto Local Threshold plugin is not installed")
    imp = sp_image("sp_local.tif", [sp_plateau_processor()])
    result = sp_run(E, tmp, "sp_local", imp, [sp_config(E, {"threshold_mode": "Local Phansalkar", "local_radius": 15})],
                    [Roi(0, 0, SP_W, SP_H)], {"spectrum_low": "0", "spectrum_high": "2500", "spectrum_steps": "26"})
    imp.close()
    configured, levels = sp_rows(result, 0)
    summary = result["summary_rows"][0]
    t.check(configured is None, "local threshold: no configured-threshold row")
    t.equal([row["Positive_Pixels"] for row in levels], [sp_at_least(100.0 * index) for index in range(26)],
            "local threshold: the global spectrum of the prepared image")
    t.equal(levels[0]["Threshold_Mode"], "Local Phansalkar", "local threshold: the mode is recorded")
    t.check("no single threshold" in levels[0]["Note"], "local threshold: the table says why values are missing")
    for suffix in ["_Spectrum_Threshold", "_Spectrum_Coverage_Percent", "_Spectrum_Sensitivity_PP_Per_Percent",
                   "_Spectrum_Elasticity", "_Spectrum_Stable_From", "_Spectrum_Stable_To", "_Spectrum_Stable_Width_Percent"]:
        t.equal(summary["C1_MAP2" + suffix], "", "local threshold: " + suffix + " is empty")
    t.close(summary["C1_MAP2_Spectrum_Normalized_Mean_Percent"], 7.9, "local threshold: the normalized mean does not need a threshold")
    t.check("no single threshold" in summary["C1_MAP2_Spectrum_Note"], "local threshold: the summary says why values are missing")


def sp_uniform_processor(eight_bit, offset=0):
    """256 x 20 pixels, column x holds offset + x: each of the 256 values is
    held by exactly 20 pixels, so the coverage curve is a straight line."""
    ip = ByteProcessor(256, 20) if eight_bit else ShortProcessor(256, 20)
    for x in range(256):
        ip.setValue(float(offset + x))
        ip.fill(Roi(x, 0, 1, 20))
    ip.resetRoi()
    return ip


def sp_micron(imp, size):
    cal = SpCalibration()
    cal.pixelWidth = size
    cal.pixelHeight = size
    cal.setUnit("micron")
    imp.setCalibration(cal)
    return imp


def sp_walk_clipped(ip, roi, low):
    """Independent answer for a ROI that may reach past the image border:
    (pixels of the ROI inside the image, those with value >= low)."""
    bounds = roi.getBounds()
    mask = roi.getMask()
    pixels = 0
    hits = 0
    for y in range(bounds.height):
        for x in range(bounds.width):
            px = bounds.x + x
            py = bounds.y + y
            if px < 0 or py < 0 or px >= ip.getWidth() or py >= ip.getHeight():
                continue
            if mask is not None and mask.get(x, y) == 0:
                continue
            pixels += 1
            if ip.getf(px, py) >= low:
                hits += 1
    return (pixels, hits)


def test_spectrum_verifier_cases(E, t, tmp):
    """Known answers derived independently: a straight-line spectrum at 8 and
    16 bits (sensitivity, elasticity, stable range, both directions), a
    border ROI with an open range, not-a-number pixels, two frames with a
    calibrated oval ROI cut by the border, scaled automatic and Display
    minimum thresholds, and a profile without any spectrum key."""
    # A. 8-bit straight line, calibrated 0.25 micron, Manual threshold 125.
    # Whole image: 5120 pixels, 20 * (256 - 125) = 2620 at or above 125,
    # coverage 51.171875 percent. The curve falls 100 / 256 points per unit,
    # so the sensitivity is -(100 / 256) * 1.25 = -0.48828125 points per
    # percent at every threshold (ImageJ applies 123.75 and 126.25 as 123
    # and 126; the typed width 2.5 would give -0.5859375). Stable range at 10
    # percent: 2358 to 2882 pixels, 20 * (256 - t) in that band for t = 112
    # to 138, width 100 * 26 / 125 = 20.8 percent. Mean 127.5 on the axis
    # 0..255: normalized mean 50 percent.
    # Second ROI: 100 x 40 at (200, -5) is cut by the top, bottom and right
    # borders to x 200..255, y 0..19: 1120 pixels, all at 200 or more, so
    # 100 percent. No threshold below can add pixels: open at the low end, cut
    # at the configured 125; 1008 pixels (90 percent) last up to t = 205.
    imp = sp_micron(sp_image("sp_line8.tif", [sp_uniform_processor(True)]), 0.25)
    rois = [Roi(0, 0, 256, 20), Roi(200, -5, 100, 40)]
    result = sp_run(E, tmp, "sp_line8", imp, [sp_config(E, {"manual_low": 125})], rois, {"spectrum_steps": "256"})
    imp.close()
    configured, levels = sp_rows(result, 0)
    summary = result["summary_rows"][0]
    t.equal(summary["Spatial_Unit"], "micron", "8-bit line: the image is calibrated")
    t.equal([row["Threshold"] for row in levels], [float(index) for index in range(256)], "8-bit line: axis 0..255 in 256 levels")
    t.equal([row["Positive_Pixels"] for row in levels], [20 * (256 - index) for index in range(256)], "8-bit line: exact count at every level")
    t.equal((configured["ROI_Pixels"], configured["Positive_Pixels"]), (5120, 2620), "8-bit line: configured threshold 125")
    t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], 51.171875, "8-bit line: coverage", 1e-12)
    t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], summary["C1_MAP2_Positive_Area_Fraction_Percent"],
            "8-bit line, calibrated: coverage equals the positive-area fraction", 1e-12)
    t.close(summary["C1_MAP2_Spectrum_Sensitivity_PP_Per_Percent"], -0.48828125, "8-bit line: sensitivity is the line's slope", 1e-12)
    t.close(summary["C1_MAP2_Spectrum_Elasticity"], 100.0 * -0.48828125 / 51.171875, "8-bit line: elasticity", 1e-12)
    t.equal((summary["C1_MAP2_Spectrum_Stable_From"], summary["C1_MAP2_Spectrum_Stable_To"]), (112.0, 138.0), "8-bit line: stable range")
    t.close(summary["C1_MAP2_Spectrum_Stable_Width_Percent"], 20.8, "8-bit line: stable width", 1e-12)
    t.close(summary["C1_MAP2_Spectrum_Normalized_Mean_Percent"], 50.0, "8-bit line: normalized mean", 1e-12)
    configured, levels = sp_rows(result, 1)
    summary = result["summary_rows"][1]
    t.equal((configured["ROI_Pixels"], configured["Positive_Pixels"]), (1120, 1120), "border ROI: clipped to 56 x 20, all positive")
    t.equal([row["Positive_Pixels"] for row in levels], [min(1120, 20 * (256 - index)) for index in range(256)],
            "border ROI: exact count at every level")
    t.equal((summary["C1_MAP2_Spectrum_Stable_From"], summary["C1_MAP2_Spectrum_Stable_To"]), (125.0, 205.0),
            "border ROI: open range cut at the configured threshold")
    t.check("open at its low end" in summary["C1_MAP2_Spectrum_Note"], "border ROI: the range is noted as open")
    t.close(summary["C1_MAP2_Spectrum_Stable_Width_Percent"], 64.0, "border ROI: stable width 100 * 80 / 125", 1e-12)
    t.equal(summary["C1_MAP2_Spectrum_Sensitivity_PP_Per_Percent"], 0.0, "border ROI: no pixel near the threshold")

    # B. 16-bit straight line, values 1000..1255, Manual threshold 1125: one
    # percent is 11.25, and 1113.75 and 1136.25 are applied as 1114 and
    # 1136, 22 values or 440 of 5120 pixels: 8.59375 points over 22 units,
    # 100 / 256 points per unit as the line has, so -(100 / 256) * 11.25 =
    # -4.39453125 (the typed width 22.5 would give -4.296875).
    imp = sp_image("sp_line16.tif", [sp_uniform_processor(False, 1000)])
    result = sp_run(E, tmp, "sp_line16", imp, [sp_config(E, {"manual_low": 1125})], [Roi(0, 0, 256, 20)], {"spectrum_steps": "2"})
    imp.close()
    summary = result["summary_rows"][0]
    t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], 100.0 * 20 * 131 / 5120.0, "16-bit line: coverage", 1e-12)
    t.close(summary["C1_MAP2_Spectrum_Sensitivity_PP_Per_Percent"], -4.39453125, "16-bit line: sensitivity is the line's slope", 1e-12)

    # C. Light background, at or below: Otsu on the 8-bit line with a dark
    # foreground. With T the recorded upper threshold (lower = 0), the
    # count is 20 * (T + 1), the curve rises 100 / 256 points per unit, so
    # the sensitivity is +(100 / 256) * T / 100, and the stable range runs
    # from the smallest t with 20 * (t + 1) >= least to the largest with
    # 20 * (t + 1) <= most.
    imp = sp_image("sp_line_light.tif", [sp_uniform_processor(True)])
    result = sp_run(E, tmp, "sp_line_light", imp, [sp_config(E, {"threshold_mode": "Auto", "auto_method": "Otsu", "dark_background": False})],
                    [Roi(0, 0, 256, 20)], {"spectrum_steps": "256"})
    imp.close()
    configured, levels = sp_rows(result, 0)
    summary = result["summary_rows"][0]
    top = result["threshold_rows"][0]["Threshold_High"]
    t.equal(result["threshold_rows"][0]["Threshold_Low"], 0.0, "light line: the lower threshold is the image minimum")
    t.check(20 < top < 235 and top == int(top), "light line: Otsu gives a whole threshold inside the range")
    count = 20 * (int(top) + 1)
    least = int(math.ceil(0.9 * count - 1e-9))
    most = int(math.floor(1.1 * count + 1e-9))
    t.equal(configured["Direction"], "below", "light line: counted at or below")
    t.equal(configured["Positive_Pixels"], count, "light line: values 0..T are selected")
    t.equal([row["Positive_Pixels"] for row in levels], [20 * (index + 1) for index in range(256)], "light line: exact count at every level")
    t.close(summary["C1_MAP2_Spectrum_Sensitivity_PP_Per_Percent"], (100.0 / 256) * top / 100.0, "light line: sensitivity is positive", 1e-12)
    t.close(summary["C1_MAP2_Spectrum_Elasticity"], 100.0 * (100.0 / 256) * top / 100.0 / (100.0 * count / 5120.0),
            "light line: elasticity", 1e-12)
    t.equal((summary["C1_MAP2_Spectrum_Stable_From"], summary["C1_MAP2_Spectrum_Stable_To"]),
            (float(int(math.ceil(least / 20.0)) - 1), float(most // 20 - 1)), "light line: stable range")
    t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], summary["C1_MAP2_Positive_Area_Fraction_Percent"],
            "light line: coverage equals the positive-area fraction", 1e-12)

    # D. A 32-bit image with not-a-number pixels inside the ROI: 40 x 30 at
    # 10.0, a 10 x 10 block of 50.0, and a 5 x 30 strip of NaN. ImageJ
    # leaves NaN out of the ROI area, so threshold 20 gives 100 of 1050
    # pixels. Stable range at 10 percent: 90 to 110 pixels, every threshold
    # above 10.0 (exclusive, the 32-bit convention) up to 50.0. A ROI that
    # holds only NaN has nothing to measure.
    plane = E.FloatProcessor(40, 30)
    plane.setValue(10.0)
    plane.fill()
    plane.setValue(50.0)
    plane.fill(Roi(20, 10, 10, 10))
    plane.setValue(float("nan"))
    plane.fill(Roi(0, 0, 5, 30))
    plane.resetRoi()
    imp = sp_image("sp_nan.tif", [plane])
    result = sp_run(E, tmp, "sp_nan", imp, [sp_config(E, {"manual_low": 20})], [Roi(0, 0, 40, 30), Roi(0, 0, 5, 30)],
                    {"spectrum_low": "0", "spectrum_high": "100", "spectrum_steps": "11"})
    imp.close()
    configured, levels = sp_rows(result, 0)
    summary = result["summary_rows"][0]
    t.equal((configured["ROI_Pixels"], configured["Positive_Pixels"]), (1050, 100), "NaN: they are left out of the ROI pixels")
    t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], 100.0 * 100 / 1050.0, "NaN: coverage of the pixels with a number", 1e-12)
    t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], summary["C1_MAP2_Positive_Area_Fraction_Percent"],
            "NaN: coverage equals the positive-area fraction", 1e-12)
    t.check("150 pixel(s)" in summary["C1_MAP2_Spectrum_Note"], "NaN: the note gives how many were left out")
    t.equal([row["Positive_Pixels"] for row in levels], [1050, 1050, 100, 100, 100, 100, 0, 0, 0, 0, 0], "NaN: exact count at 0, 10, ... 100")
    t.close(levels[1]["Coverage_Percent"], 100.0, "NaN: at 10 every pixel with a number is positive")
    t.equal((summary["C1_MAP2_Spectrum_Stable_From"], summary["C1_MAP2_Spectrum_Stable_To"]), (10.0, 50.0), "NaN: stable range")
    # (950 * 10 + 100 * 50) / 1050 on the axis 0..100.
    t.close(summary["C1_MAP2_Spectrum_Normalized_Mean_Percent"], (950 * 10.0 + 100 * 50.0) / 1050.0, "NaN: normalized mean", 1e-9)
    configured, levels = sp_rows(result, 1)
    summary = result["summary_rows"][1]
    t.equal(levels[0]["ROI_Pixels"], 0, "only NaN: no pixel to measure")
    t.check(configured is None, "only NaN: no configured-threshold row")
    t.equal(summary["C1_MAP2_Spectrum_Coverage_Percent"], "", "only NaN: coverage is empty")
    t.check("150 pixel(s)" in summary["C1_MAP2_Spectrum_Note"] and "no pixels" in summary["C1_MAP2_Spectrum_Note"],
            "only NaN: the note says why")

    # E. Two frames, calibrated 0.5 micron, an oval cut by the top and left
    # borders. Frame 2 has its 500 block lowered to 250, below the threshold
    # 300. Expected counts come from walking the oval's own mask.
    planes = [sp_plateau_processor(), sp_plateau_processor((100, 250, 2000))]
    stack = ImageStack(SP_W, SP_H)
    for plane in planes:
        stack.addSlice("", plane)
    imp = sp_micron(ImagePlus("sp_frames.tif", stack), 0.5)
    imp.setDimensions(1, 1, 2)
    imp.setOpenAsHyperStack(True)
    oval = OvalRoi(-20, -15, 60, 50)
    options = sp_options(E, {"spectrum_low": "0", "spectrum_high": "2500", "spectrum_steps": "26"})
    run_dir = os.path.join(tmp, "sp_frames")
    os.makedirs(run_dir)
    config = sp_config(E)
    result = E.AnalysisEngine().analyze_image(imp, "", [config], [E.clone_roi(oval)], E.engine_options_from_global(options),
                                              run_dir, [1, 2], settings=(options, [config]))
    imp.close()
    rows = result["extra_tables"]["threshold_spectrum.csv"]["rows"]
    t.equal(len(rows), 2 * 27, "two frames: one set of rows per frame")
    t.equal(len(result["summary_rows"]), 2, "two frames: one summary row per frame")
    for frame in [1, 2]:
        plane = planes[frame - 1]
        frame_rows = [row for row in rows if row["Frame"] == frame]
        configured = [row for row in frame_rows if row["Level_Index"] == 0][0]
        walked = [sp_walk_clipped(plane, oval, 100.0 * index) for index in range(26)]
        label = "frame " + str(frame)
        t.check(0 < walked[0][0] < 60 * 50, label + ": the oval is cut by the border")
        t.equal(configured["ROI_Pixels"], walked[0][0], label + ": ROI pixels of the clipped oval")
        t.equal(configured["Positive_Pixels"], sp_walk_clipped(plane, oval, 300)[1], label + ": configured threshold")
        t.equal([row["Positive_Pixels"] for row in frame_rows if row["Level_Index"] > 0], [item[1] for item in walked],
                label + ": exact count at every level")
        summary = result["summary_rows"][frame - 1]
        t.equal(summary["Frame"], frame, label + ": summary row order")
        t.close(summary["C1_MAP2_Spectrum_Coverage_Percent"], summary["C1_MAP2_Positive_Area_Fraction_Percent"],
                label + ": coverage equals the positive-area fraction", 1e-12)
    t.check(rows[0]["Positive_Pixels"] > [row for row in rows if row["Frame"] == 2][0]["Positive_Pixels"],
            "two frames: lowering the block in frame 2 lowers its coverage")

    # F. Thresholds that are not whole values on 8- and 16-bit images: a
    # scaled automatic threshold on a light background (upper end x 0.973),
    # a scaled automatic threshold on 8 bits (lower end x 1.0137), and a
    # Display-minimum threshold. Each must count exactly what the mask holds.
    cases = [
        ("sp_scaled_light", sp_uniform_processor(False, 1000),
         {"threshold_mode": "Auto", "auto_method": "Otsu", "dark_background": False, "auto_scale": 0.973}, None),
        ("sp_scaled_8bit", sp_uniform_processor(True),
         {"threshold_mode": "Auto", "auto_method": "Otsu", "auto_scale": 1.0137}, None),
        ("sp_display_min", sp_uniform_processor(False, 1000), {"threshold_mode": "Display minimum"}, (1100.6, 1255.0))
    ]
    for name, plane, overrides, display in cases:
        imp = sp_image(name + ".tif", [plane])
        if display is not None:
            imp.setDisplayRange(display[0], display[1])
        result = sp_run(E, tmp, name, imp, [sp_config(E, overrides)], [Roi(0, 0, 256, 20), sp_triangle()], {"spectrum_steps": "2"})
        imp.close()
        for index in range(2):
            configured = sp_rows(result, index)[0]
            summary = result["summary_rows"][index]
            t.equal(float(configured["Positive_Pixels"]), summary["C1_MAP2_Positive_Area"], name + ", ROI " + str(index + 1) + ": the mask's count")
            t.check("CHECK FAILED" not in summary["C1_MAP2_Spectrum_Note"], name + ", ROI " + str(index + 1) + ": the module's own check passes")
        moving = "Threshold_High" if overrides.get("dark_background") is False else "Threshold_Low"
        t.equal(sp_rows(result, 0)[0]["Threshold"], result["threshold_rows"][0][moving], name + ": the recorded threshold is used")

    # G. A profile from before this module (no spectrum key at all) gives
    # the same tables as the spectrum switched off, byte for byte.
    configs = [nuclei_config(E), marker_config(E)]
    rois = [Roi(0, 0, WIDTH, HEIGHT), triangle_roi()]
    imp = make_image("sp_legacy.tif")
    off = sp_run(E, tmp, "sp_legacy_off", imp, configs, rois, None, False)
    legacy = sp_options(E, None, False)
    for spec in E.SPECTRUM_SETTING_SPECS:
        del legacy[spec[1]]
    run_dir = os.path.join(tmp, "sp_legacy_none")
    os.makedirs(run_dir)
    none = E.AnalysisEngine().analyze_image(imp, "", configs, [E.clone_roi(roi) for roi in rois],
                                            E.engine_options_from_global(legacy), run_dir, [1], settings=(legacy, configs))
    imp.close()
    for name in ["roi_summary.csv", "objects.csv", "thresholds.csv"]:
        t.equal(read_utf8(os.path.join(none["run_dir"], name)).replace("sp_legacy_none", "sp_legacy_off"),
                read_utf8(os.path.join(off["run_dir"], name)), name + ": a profile without spectrum keys writes the same table")
    t.equal(none["extra_table_names"], [], "a profile without spectrum keys adds no table")


TESTS.append(("threshold spectrum: exact counts for whole, border and polygon ROIs", test_spectrum_known_plateaus))
TESTS.append(("threshold spectrum: sensitivity, stable range and area identity on a ramp", test_spectrum_ramp_summaries))
TESTS.append(("threshold spectrum: agrees with the engine's own thresholds at every bit depth", test_spectrum_matches_engine_thresholds))
TESTS.append(("threshold spectrum: cleanup note and degenerate cases", test_spectrum_degenerate_cases))
TESTS.append(("threshold spectrum: off by default, settings, and batch on a shared axis", test_spectrum_off_settings_and_batch))
TESTS.append(("threshold spectrum: local threshold channel", test_spectrum_local_threshold))
TESTS.append(("threshold spectrum: straight-line, NaN, two-frame and legacy-profile known answers", test_spectrum_verifier_cases))

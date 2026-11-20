# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach

# Cell Quant Studio Universal - regression tests for engine fixes of 2.3.0
#
# Loaded by tests/run_tests.py, which executes this file in its own
# namespace (so the helpers defined there, such as make_image, nuclei_config
# and Recorder, are available here). Every test has a known answer by
# construction. Each states which defect it covers and how it failed before
# the fix ("Before the fix: ...").
#
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij import IJ, ImagePlus, ImageStack, CompositeImage
from ij.process import ByteProcessor, ShortProcessor, FloatProcessor, ColorProcessor, ImageProcessor
from ij.gui import Roi, OvalRoi
from ij.io import FileSaver, RoiEncoder
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import Analyzer, ParticleAnalyzer
from java.awt import GraphicsEnvironment
from java.io import File, PrintStream, ByteArrayOutputStream
from java.lang import System
from java.util import Properties, Random

import os
import sys
import codecs


def fix_raises(function, word=None):
    """True when function() raises ValueError (whose text holds word)."""
    try:
        function()
    except ValueError as exc:
        return word is None or word in E_TEXT(exc)
    return False


def E_TEXT(value):
    try:
        return unicode(value)
    except NameError:
        return str(value)


def fix_csv_rows(path):
    """Rows of a CSV written by write_csv (every field quoted) as dicts."""
    handle = codecs.open(path, "r", "utf-8")
    try:
        lines = [line.rstrip(u"\n") for line in handle.read().split(u"\n") if line != u""]
    finally:
        handle.close()

    def split(line):
        fields = []
        position = 0
        while position < len(line):
            end = position + 1
            text = []
            while True:
                quote = line.index(u'"', end)
                if quote + 1 < len(line) and line[quote + 1] == u'"':
                    text.append(line[end:quote + 1])
                    end = quote + 2
                    continue
                text.append(line[end:quote])
                break
            fields.append(u"".join(text))
            position = quote + 2
        return fields

    headers = split(lines[0])
    return [dict(zip(headers, split(line))) for line in lines[1:]]


def fix_read_bytes(path):
    handle = open(path, "rb")
    try:
        return handle.read()
    finally:
        handle.close()


def fix_options(E, **changes):
    options = E.default_global_options()
    options.update(changes)
    return options


def fix_analyze(E, imp, configs, rois, run_dir, options=None, log=None):
    options = options if options is not None else fix_options(E)
    engine = E.AnalysisEngine(log=log)
    return engine.analyze_image(imp, "", configs, rois, E.engine_options_from_global(options), run_dir,
                                [1], settings=(options, configs))


def fix_stack_image(title, planes, calibrate=False):
    """A one-slice, one-frame image with one channel per processor."""
    stack = ImageStack(planes[0].getWidth(), planes[0].getHeight())
    for plane in planes:
        stack.addSlice("", plane)
    imp = ImagePlus(title, stack)
    imp.setDimensions(len(planes), 1, 1)
    if len(planes) > 1:
        imp = CompositeImage(imp, CompositeImage.GRAYSCALE)
    if calibrate:
        cal = imp.getCalibration()
        cal.pixelWidth = PIXEL
        cal.pixelHeight = PIXEL
        cal.setUnit("um")
    return imp


# -----------------------------------------------------------------------------
# A. An empty threshold range gives zero instead of an error
# -----------------------------------------------------------------------------

def test_fix_empty_threshold_range(E, t, tmp):
    """A field whose brightest pixel is below the threshold (a negative
    control) has zero positive area and zero objects.
    Before the fix: build_masks raised "Threshold upper value is below the
    lower value", so analyze_image failed and a batch reported the image as
    FAILED."""
    engine = E.AnalysisEngine()
    flat = ImagePlus("flat", blank_plane(100))
    data = engine.build_masks(flat, nuclei_config(E, {"manual_low": 1000}), "flat")
    t.equal(data["empty_range"], True, "Manual 1000 on a field of 100 is an empty range")
    t.equal(data["raw"].getProcessor().getHistogram()[255], 0, "the raw mask is all background")
    t.equal(data["final"].getProcessor().getHistogram()[255], 0, "the final mask is all background")
    t.equal(data["threshold_low"], 1000.0, "the configured lower value is recorded")
    t.equal(data["threshold_high"], 100.0, "the upper value is the data maximum")
    for key in ["raw", "final"]:
        E.safe_close_image(data[key])
    two_level = ImagePlus("two", nuclei_plane(DISKS))
    normal = engine.build_masks(two_level, nuclei_config(E), "two")
    t.equal(normal["empty_range"], False, "an ordinary threshold is not an empty range")
    E.safe_close_image(normal["raw"])
    E.safe_close_image(normal["final"])
    low_band = engine.build_masks(flat, nuclei_config(E, {"manual_low": 10, "manual_high": 50}), "band")
    t.equal(low_band["empty_range"], True, "Manual 10 to 50 on data of 100 is an empty range")
    E.safe_close_image(low_band["raw"])
    E.safe_close_image(low_band["final"])
    scaled = engine.build_masks(two_level, nuclei_config(E, {"threshold_mode": "Auto", "auto_scale": 30.0}), "otsu")
    t.equal(scaled["empty_range"], True, "Otsu times 30 above the brightest pixel is an empty range")
    E.safe_close_image(scaled["raw"])
    E.safe_close_image(scaled["final"])

    # The whole run, with every optional module loaded and the depth profile
    # and threshold spectrum switched on.
    imp = make_image("negative_control.tif", slices=1, disks_per_frame=[[]])
    configs = [nuclei_config(E, {"watershed": True, "fill_holes": True}), marker_config(E, {"manual_low": 5000})]
    options = fix_options(E, depth_enabled=True, spectrum_enabled=True)
    lines = []
    result = fix_analyze(E, imp, configs, [E.whole_image_roi(imp), OvalRoi(20, 20, 100, 80)],
                         E.unique_run_dir(tmp, "negative_control"), options, lines.append)
    for row in result["summary_rows"]:
        t.equal(row["C1_DAPI_Object_Count"], 0, "no objects in " + E.to_text(row["Parent_ROI_Name"]))
        t.equal(row["C1_DAPI_Positive_Area"], 0.0, "no positive nuclear area")
        t.equal(row["C2_MAP2_Positive_Area"], 0.0, "no positive marker area")
        t.equal(row["C2_MAP2_Positive_Area_Fraction_Percent"], 0.0, "a positive fraction of 0 percent")
    t.equal(len(result["object_rows"]), 0, "objects.csv has no rows")
    t.equal([row["Threshold_Low"] for row in result["threshold_rows"]], [1000.0, 5000.0],
            "thresholds.csv records the configured lower values")
    t.check(len([line for line in lines if "no pixel value lies in the threshold range" in line]) == 2,
            "the run log says why nothing is positive, once per channel")
    imp.close()

    # A batch reports the image as analyzed, with zero objects.
    folder = os.path.join(tmp, "fix_empty_batch")
    os.makedirs(folder)
    control = make_image("control.tif", slices=1, disks_per_frame=[[]])
    FileSaver(control).saveAsTiff(os.path.join(folder, "control.tif"))
    control.close()
    batch = E.run_batch(folder, os.path.join(tmp, "fix_empty_batch_out"), configs, configs,
                        fix_options(E), {"extensions": "tif", "roi_policy": E.BATCH_ROI_WHOLE})
    t.equal((batch["analyzed"], batch["failed"]), (1, 0), "the batch analyzes the negative control")
    report = fix_csv_rows(os.path.join(batch["batch_dir"], "batch_report.csv"))
    t.equal(report[0]["Status"], "OK", "batch_report.csv says OK")
    summary = fix_csv_rows(os.path.join(batch["batch_dir"], "batch_roi_summary.csv"))
    t.equal(summary[0]["C1_DAPI_Object_Count"], "0", "batch_roi_summary.csv counts 0 objects")


def test_fix_sensitivity_empty_range(E, t, tmp):
    """A sensitivity factor that lifts the threshold above every pixel gives
    zero with a note. Before the fix: that factor was skipped with the error
    text as its note and no values."""
    imp = ImagePlus("sens", nuclei_plane(DISKS))
    config = nuclei_config(E)
    result = E.run_sensitivity(E.AnalysisEngine(), imp, config, [E.whole_image_roi(imp)], E.SENS_THRESHOLD,
                               [1.0, 3.0], E.unique_dir(tmp, "fix_sensitivity"), "sens.tif", 1)
    rows = dict([(row["Factor"], row) for row in result["rows"]])
    t.equal(rows[3.0]["Parameter_Value"], 3000.0, "factor 3 asks for a threshold of 3000")
    t.equal(rows[3.0]["Positive_Area"], 0.0, "nothing is positive above the brightest pixel (2000)")
    t.equal(rows[3.0]["Object_Count"], 0, "and there are no objects")
    t.equal(rows[3.0]["Positive_Area_Change_Percent"], -100.0, "a change of -100 percent")
    t.equal(rows[3.0]["Note"], E.SENSITIVITY_EMPTY_NOTE, "the note says why")
    t.equal(rows[1.0]["Object_Count"], len(DISKS), "the configured threshold finds every disk")
    t.equal(rows[1.0]["Note"], "", "the configured threshold has no note")
    imp.close()


def test_fix_sensitivity_java_error(E, t, tmp):
    """An error raised inside ImageJ or a plugin (a Java error) for one
    factor is noted for that factor and the sweep goes on. Before the fix:
    only Python errors were caught, so the whole sweep stopped. A cancel
    still stops it."""
    from java.lang import IllegalStateException

    class FailingEngine(E.AnalysisEngine):
        failure = None

        def build_masks(self, prepared_imp, config, title_prefix):
            if config["manual_low"] > 1500:
                raise self.failure
            return E.AnalysisEngine.build_masks(self, prepared_imp, config, title_prefix)

    imp = ImagePlus("sens java", nuclei_plane(DISKS))
    config = nuclei_config(E)
    engine = FailingEngine()
    engine.failure = IllegalStateException("plugin failure")
    result = E.run_sensitivity(engine, imp, config, [E.whole_image_roi(imp)], E.SENS_THRESHOLD,
                               [1.0, 1.2, 1.8], E.unique_dir(tmp, "fix_sensitivity_java"), "sens.tif", 1)
    rows = dict([(row["Factor"], row) for row in result["rows"]])
    t.equal(sorted(rows.keys()), [1.0, 1.2, 1.8], "every factor has a row")
    t.check("plugin failure" in E_TEXT(rows[1.8]["Note"]), "the failed factor carries the error as its note")
    t.equal(rows[1.2]["Object_Count"], len(DISKS), "the factors after and before it are still measured")
    t.equal(rows[1.2]["Note"], "", "a working factor has no note")

    engine.failure = E.AnalysisCancelled("Cancelled by user.")
    cancelled = False
    try:
        E.run_sensitivity(engine, imp, config, [E.whole_image_roi(imp)], E.SENS_THRESHOLD,
                          [1.0, 1.8], E.unique_dir(tmp, "fix_sensitivity_cancel"), "sens.tif", 1)
    except E.AnalysisCancelled:
        cancelled = True
    t.check(cancelled, "a cancel still ends the sweep")
    imp.close()


# -----------------------------------------------------------------------------
# B. Images with a calibration function are refused for segmentation
# -----------------------------------------------------------------------------

def test_fix_calibration_function_refused(E, t, tmp):
    """Signed 16-bit data: thresholds would act on stored values while the
    statistics are calibrated. Before the fix: Manual 1000 on a disk at 2000
    (calibrated) found 0 positive pixels without any message."""
    ip = ShortProcessor(60, 60)
    ip.setValue(32768 + 100)
    ip.fill()
    ip.setValue(32768 + 2000)
    ip.fill(OvalRoi(10, 10, 30, 30))
    ip.resetRoi()
    imp = ImagePlus("signed.tif", ip)
    imp.getCalibration().setSigned16BitCalibration()
    t.check(E.has_calibration_function(imp), "the image has a calibration function")
    t.close(E.AnalysisEngine().measure_intensity(imp, Roi(0, 0, 5, 5))["mean"], 100.0, "its values are calibrated")
    segmented = nuclei_config(E, {"analysis_mode": "Positive area"})
    t.check(fix_raises(lambda: E.validate_channel_configs(imp, [segmented]), "calibration function"),
            "a segmented channel is refused with a message naming the calibration function")
    t.check(fix_raises(lambda: E.validate_channel_configs(imp, [segmented]), "32-bit"),
            "and the message says to convert to 32-bit")
    intensity_only = nuclei_config(E, {"analysis_mode": "Intensity only"})
    E.validate_channel_configs(imp, [intensity_only])
    t.check(True, "an intensity-only channel is accepted")
    plain = ImagePlus("plain.tif", nuclei_plane(DISKS))
    t.check(not E.has_calibration_function(plain), "an ordinary 16-bit image has no calibration function")
    E.validate_channel_configs(plain, [segmented])

    folder = os.path.join(tmp, "fix_signed_batch")
    os.makedirs(folder)
    FileSaver(imp).saveAsTiff(os.path.join(folder, "signed.tif"))
    batch = E.run_batch(folder, os.path.join(tmp, "fix_signed_out"), [segmented], [segmented], fix_options(E),
                        {"extensions": "tif", "roi_policy": E.BATCH_ROI_WHOLE})
    report = fix_csv_rows(os.path.join(batch["batch_dir"], "batch_report.csv"))
    if report[0]["Status"] == "FAILED":
        t.check("calibration function" in report[0]["Message"], "a batch reports the reason for the image")
    else:
        # A TIFF reader that drops the calibration function gives plain
        # unsigned data, which is then segmented normally.
        t.equal(report[0]["Status"], "OK", "the reopened image has no calibration function")
    imp.close()
    plain.close()


# -----------------------------------------------------------------------------
# C. Per-object intensity on the object's own pixels
# -----------------------------------------------------------------------------

def fix_ring_image():
    """80 x 80 pixels, two channels, uncalibrated. A square ring of 500
    pixels (30 x 30 minus a 20 x 20 hole) with a 6 x 6 object (36 pixels)
    in the middle of the hole, 7 pixels from the ring on every side.
    C1: background 100, ring 2000, hole 300, inner object 5000.
    C2: background 10, ring 700, hole 50, inner object 900."""
    planes = []
    for values in [(100, 2000, 300, 5000), (10, 700, 50, 900)]:
        ip = ShortProcessor(80, 80)
        ip.setValue(values[0])
        ip.fill()
        ip.setValue(values[1])
        ip.fill(Roi(20, 20, 30, 30))
        ip.setValue(values[2])
        ip.fill(Roi(25, 25, 20, 20))
        ip.setValue(values[3])
        ip.fill(Roi(32, 32, 6, 6))
        ip.resetRoi()
        planes.append(ip)
    return fix_stack_image("ring.tif", planes)


def test_fix_object_intensity_ring(E, t, tmp):
    """Before the fix: the ring row reported Object_Area 500 but Mean
    1432.44 and Raw Integrated Density 1289200, because its outline also
    covers the 364 hole pixels at 300 and the 36 pixels of the inner object
    at 5000 ((500 x 2000 + 364 x 300 + 36 x 5000) / 900)."""
    imp = fix_ring_image()
    c1 = E.normalize_channel_config(1, {"enabled": True, "name": "Nuc", "analysis_mode": "Objects only",
                                        "threshold_mode": "Manual", "manual_low": 1000, "cross_measure": True})
    c2 = E.normalize_channel_config(2, {"enabled": True, "name": "Mark", "analysis_mode": "Intensity only"})
    result = fix_analyze(E, imp, [c1, c2], [E.whole_image_roi(imp)], E.unique_run_dir(tmp, "fix_ring"))
    rows = sorted(result["object_rows"], key=lambda row: -row["Object_Area"])
    t.equal(len(rows), 2, "the ring and the object inside its hole are two objects")
    ring = rows[0]
    inner = rows[1]
    t.equal(ring["Object_Area"], 500.0, "the ring's area excludes its hole")
    expected = [("Nuc", 2000.0, 1000000.0), ("Mark", 700.0, 350000.0)]
    for name, mean, density in expected:
        prefix = ("C1_" if name == "Nuc" else "C2_") + name
        t.equal(ring[prefix + "_Object_Mean"], mean, name + ": ring mean on its own pixels")
        t.equal(ring[prefix + "_Object_StdDev"], 0.0, name + ": ring standard deviation")
        t.equal((ring[prefix + "_Object_Min"], ring[prefix + "_Object_Max"], ring[prefix + "_Object_Median"]),
                (mean, mean, mean), name + ": ring min, max and median")
        t.equal(ring[prefix + "_Object_Integrated_Density"], density, name + ": ring integrated density")
        t.equal(ring[prefix + "_Object_Raw_Integrated_Density"], density, name + ": ring raw integrated density")
    t.equal(inner["Object_Area"], 36.0, "the inner object has 36 pixels")
    t.equal(inner["C1_Nuc_Object_Mean"], 5000.0, "the inner object's own mean")
    t.equal(inner["C1_Nuc_Object_Raw_Integrated_Density"], 180000.0, "and raw integrated density")
    t.equal(inner["C2_Mark_Object_Mean"], 900.0, "the inner object's cross-channel mean")
    written = fix_csv_rows(os.path.join(result["run_dir"], "objects.csv"))
    t.check("2000.0" in [row["C1_Nuc_Object_Mean"] for row in written], "objects.csv holds the ring's true mean")

    engine = E.AnalysisEngine()
    prepared = {1: fix_channel(imp, 1), 2: fix_channel(imp, 2)}
    mask = engine.build_masks(prepared[1], c1, "ring")
    whole = E.whole_image_roi(imp)
    # Include holes: the ring is one filled object of 900 pixels; its pixels
    # and its outline are the same, so the value is that of 2.2 exactly.
    filled, info = engine.analyze_objects(mask["final"], whole, E.normalize_channel_config(1, dict(c1, include_holes=True)), None, prepared)
    t.equal(len(filled), 1, "with holes included the inner object is part of the ring")
    t.equal(filled[0]["area"], 900.0, "which covers 900 pixels")
    old = engine.measure_intensity(prepared[1], filled[0]["roi"])
    t.equal(filled[0]["intensities"][1][0], old["mean"], "with holes included the mean equals the outline's, bit for bit")
    t.close(filled[0]["intensities"][1][0], 1289200.0 / 900.0, "which is (500 x 2000 + 364 x 300 + 36 x 5000) / 900")
    t.equal(filled[0]["intensities"][1][6], old["raw_integrated_density"], "and so does the raw integrated density")
    # Composite ROIs: the ring's outline excludes the hole; same values.
    composite, info = engine.analyze_objects(mask["final"], whole, E.normalize_channel_config(1, dict(c1, composite_rois=True)), None, prepared)
    by_area = dict([(obj["area"], obj) for obj in composite])
    t.equal(sorted(by_area.keys()), [36.0, 500.0], "with composite ROIs the ring and the inner object are found")
    t.equal(by_area[500.0]["intensities"][1][0], 2000.0, "and the ring's mean is 2000")
    t.equal(by_area[500.0]["intensities"][1][0], engine.measure_intensity(prepared[1], by_area[500.0]["roi"])["mean"],
            "equal to the composite outline's mean")
    for imp_closing in [mask["raw"], mask["final"], prepared[1], prepared[2]]:
        E.safe_close_image(imp_closing)
    imp.close()


def fix_channel(imp, channel):
    from ij.plugin import Duplicator
    return Duplicator().run(imp, channel, channel, 1, 1, 1, 1)


def fix_random_planes(width, height, bit_depth, seed):
    """A mask of random disks (some touching) and three intensity planes of
    random values: (mask processor, [plane, plane, plane])."""
    random = Random(seed)
    mask = ByteProcessor(width, height)
    mask.setValue(255.0)
    for index in range(90):
        radius = 2 + random.nextInt(6)
        mask.fill(OvalRoi(random.nextInt(width - 16), random.nextInt(height - 16), 2 * radius + 1, 2 * radius + 3))
    mask.resetRoi()
    planes = []
    for channel in range(3):
        if bit_depth == 8:
            ip = ByteProcessor(width, height)
        elif bit_depth == 16:
            ip = ShortProcessor(width, height)
        else:
            ip = FloatProcessor(width, height)
        top = 250 if bit_depth == 8 else 60000
        for y in range(height):
            for x in range(width):
                value = random.nextInt(top)
                if bit_depth == 32:
                    ip.setf(x, y, value + random.nextFloat())
                else:
                    ip.set(x, y, value)
        planes.append(ip)
    return (mask, planes)


def fix_components(mask_ip, four_connected):
    """Connected components of the 255 pixels of a mask, found independently
    of ImageJ: a list holding the component number of every pixel (0 for
    background)."""
    width = mask_ip.getWidth()
    height = mask_ip.getHeight()
    labels = [0] * (width * height)
    steps = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    if not four_connected:
        steps += [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    current = 0
    for start in range(width * height):
        if labels[start] != 0 or mask_ip.get(start % width, start // width) == 0:
            continue
        current += 1
        labels[start] = current
        pending = [start]
        while len(pending) > 0:
            position = pending.pop()
            x = position % width
            y = position // width
            for dx, dy in steps:
                nx = x + dx
                ny = y + dy
                if 0 <= nx < width and 0 <= ny < height:
                    neighbor = ny * width + nx
                    if labels[neighbor] == 0 and mask_ip.get(nx, ny) != 0:
                        labels[neighbor] = current
                        pending.append(neighbor)
    return labels


def fix_own_pixels(labels, width, roi):
    """The pixels of the component that fills most of an object outline."""
    counts = {}
    inside = [(int(point.x), int(point.y)) for point in roi.getContainedPoints()]
    for x, y in inside:
        label = labels[y * width + x]
        if label != 0:
            counts[label] = counts.get(label, 0) + 1
    best = max(counts.keys(), key=lambda label: counts[label])
    return ([(x, y) for x, y in inside if labels[y * width + x] == best], len(inside))


def test_fix_object_intensity_bit_identical(E, t, tmp):
    """(1) For every object whose outline covers exactly its own pixels, the
    new per-object values equal those of 2.2 (measure_intensity on the
    outline) bit for bit, on 8-, 16- and 32-bit images, for every particle
    setting and for every channel (cross-channel columns). (2) For every
    object, with or without holes, the new mean, minimum, maximum and pixel
    count are those of its own pixels, found here by an independent
    connected-component search. Before the fix: objects with holes, and
    four-connected objects whose outline also covers a diagonal neighbour,
    included foreign pixels."""
    engine = E.AnalysisEngine()
    variants = [
        ("default", {}), ("include holes", {"include_holes": True}), ("composite ROIs", {"composite_rois": True}),
        ("four-connected", {"four_connected": True}),
        ("size and shape limits", {"min_size": 12, "max_size": 150, "min_circularity": 0.4, "size_units": "Pixel^2"})
    ]
    width = 160
    for bit_depth in [8, 16, 32]:
        mask_ip, planes = fix_random_planes(width, 120, bit_depth, 17 + bit_depth)
        mask = ImagePlus("mask", mask_ip)
        images = dict([(index + 1, ImagePlus("C" + str(index + 1), planes[index])) for index in range(3)])
        components = {False: fix_components(mask_ip, False), True: fix_components(mask_ip, True)}
        for label, change in variants:
            values = {"enabled": True, "name": "Nuc", "analysis_mode": "Objects only", "threshold_mode": "Manual"}
            values.update(change)
            config = E.normalize_channel_config(1, values)
            for parent in [Roi(0, 0, width, 120), OvalRoi(10, 8, 140, 104)]:
                objects, info = engine.analyze_objects(mask, parent, config, None, images)
                compared = 0
                mismatches = 0
                wrong = 0
                phantoms = 0
                inconsistent = 0
                first_wrong = ""
                for obj in objects:
                    outline = engine.measure_intensity(images[1], obj["roi"])["pixel_count"]
                    if outline == int(round(obj["area"])):
                        compared += 1
                        for channel in [1, 2, 3]:
                            old = engine.measure_intensity(images[channel], obj["roi"])
                            if list(obj["intensities"][channel]) != [old[key] for key in E.OBJECT_INTENSITY_KEYS]:
                                mismatches += 1
                    if obj["area"] == 0:
                        # Four-connected analysis in ImageJ reports an
                        # object that touches an earlier one only at a
                        # corner once per pixel, with area 0 (as in 2.2).
                        phantoms += 1
                        continue
                    # The measured pixels are exactly Object_Area's pixels.
                    measured = obj["intensities"][1][6] / obj["intensities"][1][0]
                    if abs(measured - obj["area"]) > 1e-6 * obj["area"]:
                        inconsistent += 1
                    # ImageJ's four-connected particles can differ by a pixel
                    # from a plain four-connected search, so the independent
                    # check is made for eight-connected objects.
                    if config["include_holes"] or config["four_connected"] or parent.getType() != Roi.RECTANGLE:
                        continue
                    own, outline_pixels = fix_own_pixels(components[config["four_connected"]], width, obj["roi"])
                    for channel in [1, 2, 3]:
                        pixel_values = [planes[channel - 1].getf(x, y) for x, y in own]
                        new = obj["intensities"][channel]
                        mean = sum(pixel_values) / float(len(pixel_values))
                        if (len(own) != int(round(obj["area"])) or new[2] != min(pixel_values) or new[3] != max(pixel_values) or
                                abs(new[0] - mean) > 1e-9 * max(1.0, abs(mean))):
                            wrong += 1
                            first_wrong = first_wrong or (
                                " (first: area " + repr(obj["area"]) + ", own pixels " + str(len(own)) + ", outline pixels " +
                                str(outline_pixels) + ", bounds " + E.to_text(obj["roi"].getBounds()) + ")")
                where = str(bit_depth) + "-bit, " + label + ", " + parent.getClass().getSimpleName()
                t.check(compared >= 5, where + ": objects whose outline holds only their own pixels are compared (" + str(compared) + ")")
                t.equal(mismatches, 0, where + ": for those every value of every channel equals 2.2's bit for bit")
                t.equal(wrong, 0, where + ": every object is measured on exactly its own pixels" + first_wrong)
                t.equal(inconsistent, 0, where + ": every object is measured on as many pixels as its Object_Area")
                if not config["four_connected"]:
                    t.equal(phantoms, 0, where + ": no object has zero area")
        mask.close()


def test_fix_object_rois_without_roi_manager(E, t, tmp):
    """Objects are collected from an overlay, not a hidden ROI Manager (which
    cannot be made without a display). The outlines, their names, colours
    and order are those of the ROI Manager route of 2.2 (compared byte for
    byte where a display exists). Before the fix: under
    -Djava.awt.headless=true every object channel failed with
    java.awt.HeadlessException."""
    engine = E.AnalysisEngine()
    mask_ip, planes = fix_random_planes(200, 150, 8, 5)
    mask = ImagePlus("mask", mask_ip)
    headless = GraphicsEnvironment.isHeadless()
    variants = [{}, {"include_holes": True}, {"composite_rois": True}, {"four_connected": True},
                {"min_size": 20, "size_units": "Pixel^2", "min_circularity": 0.5}]
    for change in variants:
        values = {"enabled": True, "name": "Nuc", "analysis_mode": "Objects only", "threshold_mode": "Manual"}
        values.update(change)
        config = E.normalize_channel_config(1, values)
        for parent in [Roi(0, 0, 200, 150), OvalRoi(15, 10, 170, 130)]:
            objects, info = engine.analyze_objects(mask, parent, config)
            label = repr(sorted(change.keys())) + " in " + parent.getClass().getSimpleName()
            t.check(len(objects) > 5, label + ": objects are found")
            if headless:
                continue
            reference = fix_manager_rois(E, mask, parent, config)
            t.equal(len(objects), len(reference), label + ": as many objects as through a ROI Manager")
            same = 0
            for index in range(min(len(objects), len(reference))):
                if list(RoiEncoder.saveAsByteArray(objects[index]["roi"])) == list(RoiEncoder.saveAsByteArray(reference[index])):
                    same += 1
            t.equal(same, len(reference), label + ": every outline, name and order is identical")
    mask.close()


def fix_manager_rois(E, mask, parent, config):
    """The object outlines exactly as 2.2 collected them: through a hidden ROI
    Manager (only possible with a display)."""
    from ij.plugin.frame import RoiManager
    inside, band = E.roi_masks(parent, mask.getWidth(), mask.getHeight())
    work_ip = mask.getProcessor().duplicate()
    if parent.getMask() is not None:
        from ij.process import Blitter
        work_ip.copyBits(inside, 0, 0, Blitter.AND)
    work = ImagePlus("work", work_ip)
    options = ParticleAnalyzer.SHOW_NONE | ParticleAnalyzer.ADD_TO_MANAGER
    if config["include_holes"]:
        options |= ParticleAnalyzer.INCLUDE_HOLES
    if config["composite_rois"]:
        options |= ParticleAnalyzer.COMPOSITE_ROIS
    if config["four_connected"]:
        options |= ParticleAnalyzer.FOUR_CONNECTED
    minimum, maximum = E.AnalysisEngine().particle_size_pixels(config, mask.getCalibration())
    bounds = parent.getBounds()
    work.setRoi(Roi(bounds.x, bounds.y, bounds.width, bounds.height))
    work_ip.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
    manager = RoiManager(False)
    try:
        ParticleAnalyzer.setRoiManager(manager)
        ParticleAnalyzer(options, Measurements.AREA, ResultsTable(), minimum, maximum,
                         config["min_circularity"], config["max_circularity"]).analyze(work)
        return [roi.clone() for roi in manager.getRoisAsArray()]
    finally:
        ParticleAnalyzer.setRoiManager(None)
        manager.reset()
        manager.dispose()


def test_fix_measurement_redirect(E, t, tmp):
    """A 'Redirect to' image set in Analyze > Set Measurements neither breaks
    object detection nor is changed by it. Before the fix: a redirect image
    of another size made analyze_objects raise "Particle analysis failed"."""
    engine = E.AnalysisEngine()
    mask = mask_with_disks([(40, 40), (100, 100), (150, 60)])
    other = ImagePlus("other", ByteProcessor(10, 10))
    t.check(E.save_measurement_redirect() is None, "no redirect is set at the start")
    Analyzer.setRedirectImage(other)
    try:
        saved = E.save_measurement_redirect()
        values = ImagePlus("values", mask.getProcessor().duplicate())
        values.setDisplayRange(5, 300)
        pixels = list(values.getProcessor().getPixels())
        objects, info = engine.analyze_objects(mask, Roi(0, 0, 200, 200), nuclei_config(E), None, {1: values})
        t.equal(len(objects), 3, "every object is found with a redirect of another size set")
        t.equal((values.getDisplayRangeMin(), values.getDisplayRangeMax()), (5.0, 300.0),
                "the measured image keeps its display range")
        t.check(list(values.getProcessor().getPixels()) == pixels and values.getProcessor().getRoi().width == 200,
                "and its pixels, with no selection left on it")
        t.check(Analyzer.isRedirectImage(), "the redirect is still set afterwards")
        after = E.save_measurement_redirect()
        t.equal(after[0], saved[0], "the redirect target is unchanged")
        t.equal(after[1], saved[1], "the redirect title is unchanged")
        t.check(after[2] is saved[2], "the redirect image is unchanged")
    finally:
        Analyzer.setRedirectImage(None)
    objects, info = engine.analyze_objects(mask, Roi(0, 0, 200, 200), nuclei_config(E))
    t.check(not Analyzer.isRedirectImage(), "without a redirect none is left set")
    mask.close()


# -----------------------------------------------------------------------------
# D. Local Phansalkar without the menu command
# -----------------------------------------------------------------------------

def test_fix_local_threshold_direct(E, t, tmp):
    """The plugin is called directly (it then also works headless, where the
    menu command is unknown) with pixel-identical output."""
    if E.local_threshold_plugin() is None:
        raise SkipTest("Auto Local Threshold plugin is not installed")
    source = nuclei_plane(DISKS).convertToByte(True)
    for x in range(0, WIDTH, 3):
        source.set(x, 100, 180)
    for parameters in [(15, 0.25, 0.5), (35, 1.0 / 3.0, 0.5), (7, -0.1, 0.6)]:
        direct = ImagePlus("direct", source.duplicate())
        E.run_local_phansalkar(direct, parameters[0], parameters[1], parameters[2])
        t.check(E.is_two_level_mask(direct.getProcessor()), "the direct call gives a mask " + repr(parameters))
        menu = ImagePlus("menu", source.duplicate())
        IJ.run(menu, "Auto Local Threshold",
               "method=Phansalkar radius=" + str(parameters[0]) + " parameter_1=" + str(parameters[1]) +
               " parameter_2=" + str(parameters[2]) + " white")
        if E.is_two_level_mask(menu.getProcessor()):
            t.equal(list(direct.getProcessor().getPixels()), list(menu.getProcessor().getPixels()),
                    "the direct call equals the menu command pixel for pixel " + repr(parameters))
        else:
            t.check(GraphicsEnvironment.isHeadless(), "only a headless Fiji lacks the menu command")
        direct.close()
        menu.close()


# -----------------------------------------------------------------------------
# E. RGB lookup table, inverted LUT, float text
# -----------------------------------------------------------------------------

def test_fix_rgb_lut_and_inverted_auto(E, t, tmp):
    """Before the fix: (1) apply_lut on an RGB image fell back to IJ.run with a
    LUT command, which shows the error 'LUTs cannot be assigned to RGB
    Images' on every run; (2) Auto Otsu 'dark' on an image with an inverted
    LUT selected the 5684 background pixels instead of the disk."""
    rgb = ImagePlus("rgb", ColorProcessor(20, 20))
    namespace = getattr(E.apply_lut, "__globals__", None) or E.apply_lut.func_globals
    calls = []

    class FixIJ(object):
        @staticmethod
        def run(*arguments):
            calls.append(arguments)

    saved = namespace["IJ"]
    namespace["IJ"] = FixIJ
    try:
        E.apply_lut(rgb, "Green")
        E.apply_lut(rgb, "Grays")
    finally:
        namespace["IJ"] = saved
    t.equal(len(calls), 0, "no LUT command is run on an RGB image")
    t.equal(rgb.getBitDepth(), 24, "the RGB image is left as it is")
    rgb.close()

    ip = ByteProcessor(100, 100)
    ip.setValue(20)
    ip.fill()
    ip.setValue(200)
    ip.fill(OvalRoi(35, 35, 31, 31))
    ip.resetRoi()
    disk = ip.getHistogram()[200]
    ip.invertLut()
    imp = ImagePlus("inverted", ip)
    engine = E.AnalysisEngine()
    for mode in ["Auto", "Manual"]:
        config = E.normalize_channel_config(1, {"enabled": True, "name": "X", "analysis_mode": "Positive area",
                                                "threshold_mode": mode, "manual_low": 100, "lut": "Original"})
        data = engine.build_masks(imp, config, "inv")
        t.equal(data["raw"].getProcessor().getHistogram()[255], disk, mode + ": the bright disk is positive")
        E.safe_close_image(data["raw"])
        E.safe_close_image(data["final"])
    t.check(imp.getProcessor().isInvertedLut(), "the prepared image keeps its inverted LUT")
    imp.close()


def test_fix_unknown_command_refused(E, t, tmp):
    """A filter or binary command that this Fiji does not know stops the run
    with a message. Before: IJ.run left the image unchanged (for example no
    watershed split) and only wrote "Unrecognized command" to the Log. The
    unknown command is given to a stand-in for IJ, so no ImageJ error
    dialog can open."""
    plane = ByteProcessor(20, 20)
    plane.setValue(100)
    plane.fill()
    imp = ImagePlus("known", plane)
    E.run_imagej_command(imp, "Invert", "")
    t.equal(imp.getProcessor().get(0, 0), 155, "a known command is applied")
    namespace = getattr(E.run_imagej_command, "__globals__", None) or E.run_imagej_command.func_globals

    class FixQuietIJ(object):
        @staticmethod
        def run(*arguments):
            pass

    saved = namespace["IJ"]
    namespace["IJ"] = FixQuietIJ
    refused = False
    try:
        try:
            E.run_imagej_command(imp, "No Such Command For The Test", "")
        except RuntimeError as exc:
            refused = "No Such Command For The Test" in E.to_text(exc)
    finally:
        namespace["IJ"] = saved
    t.check(refused, "an unknown command is refused with its name")
    imp.close()


def test_fix_float_settings_text(E, t, tmp):
    """Before the fix: to_text(1.0 / 3) gave '0.333333333333' (12 digits), so
    a saved profile did not give back the value it was saved with."""
    t.equal(E.to_text(1.0 / 3.0), u"0.3333333333333333", "every digit of 1/3 is kept")
    for value, text in [(0.25, u"0.25"), (35.0, u"35.0"), (1e-05, u"1e-05"), (500.0, u"500.0"), (-1.0, u"-1.0")]:
        t.equal(E.to_text(value), text, "short values are written as before: " + text)
    config = nuclei_config(E, {"manual_low": 1234.56789012345, "auto_scale": 1.0 / 3.0, "min_size": 0.1 + 0.2})
    path = os.path.join(tmp, "fix_float.properties")
    E.write_settings_properties(path, E.default_global_options(), [config])
    loaded = E.load_profile_settings(path)[1][0]
    for key in ["manual_low", "auto_scale", "min_size"]:
        t.equal(loaded[key], config[key], key + " comes back exactly")
    example = os.path.join(E.CELLQUANT_PACKAGE_DIR, "examples", "Example_4Channel_SOX2_Casp3_MAP2_DAPI.properties")
    if os.path.isfile(example):
        options, enabled, everything = E.load_profile_settings(example)
        t.equal(len(enabled), 4, "the example profile of 2.2 still loads with its four channels")
        props = E.read_properties_file(example)
        for config in everything:
            for key in ["manual_low", "auto_scale", "local_parameter_1", "clahe_max_slope"]:
                raw = props.getProperty("channel." + str(config["channel"]) + "." + key)
                if raw is not None:
                    t.equal(config[key], float(raw), "example C" + str(config["channel"]) + " " + key + " is read exactly")


# -----------------------------------------------------------------------------
# E. The test runner
# -----------------------------------------------------------------------------

def test_fix_runner_reports_broken_files(E, t, tmp):
    """Before the fix: one test file with a syntax error stopped the whole
    suite with a traceback and no PASS or FAIL lines; and a test that left
    Analyze's redirect set was not reported."""
    folder = os.path.join(tmp, "fix_runner")
    os.makedirs(folder)
    here = find_tests_dir()
    copy = os.path.join(folder, "run_tests.py")
    handle = open(os.path.join(here, "run_tests.py"), "rb")
    try:
        source = handle.read()
    finally:
        handle.close()
    handle = open(copy, "wb")
    try:
        handle.write(source)
    finally:
        handle.close()
    probes = {
        "test_aaa_probe.py": (
            "def fix_probe_good(E, t, tmp):\n"
            "    t.check(True, 'fine')\n"
            "TESTS.append(('runner probe passes', fix_probe_good))\n"
            "def fix_probe_leaky(E, t, tmp):\n"
            "    from ij.plugin.filter import Analyzer\n"
            "    Analyzer.setRedirectImage(ImagePlus('probe', ByteProcessor(4, 4)))\n"
            "TESTS.append(('runner probe leaves a redirect', fix_probe_leaky))\n"
        ),
        "test_bbb_probe.py": "def broken(:\n    pass\n"
    }
    for name in probes:
        handle = open(os.path.join(folder, name), "w")
        try:
            handle.write(probes[name])
        finally:
            handle.close()
    namespace = {"__file__": copy, "__name__": "__main__", "CELLQUANT_ENGINE_PATH": find_engine()}
    saved_filter = os.environ.get("CELLQUANT_TEST_FILTER")
    saved_out = sys.stdout
    saved_java_out = System.out
    captured = FixCapture()
    java_bytes = ByteArrayOutputStream()
    raised = False
    os.environ["CELLQUANT_TEST_FILTER"] = "runner probe"
    sys.stdout = captured
    System.setOut(PrintStream(java_bytes, True))
    try:
        try:
            run_source_file(copy, namespace)
        except RuntimeError:
            raised = True
    finally:
        sys.stdout = saved_out
        System.setOut(saved_java_out)
        Analyzer.setRedirectImage(None)
        if saved_filter is None:
            del os.environ["CELLQUANT_TEST_FILTER"]
        else:
            os.environ["CELLQUANT_TEST_FILTER"] = saved_filter
    report = captured.text()
    t.check("FAIL     could not load test_bbb_probe.py" in report, "a broken test file is reported as FAIL")
    t.check("PASS     runner probe passes" in report, "the tests of the other files still run")
    t.check("FAIL     runner probe leaves a redirect" in report and "measurement redirect" in report,
            "a test that leaves Analyze's redirect set fails and the report says so")
    t.check("1 passed, 3 failed" in report, "the totals count the broken file and the state left changed at the end")
    t.check(raised, "the run ends with an error")


class FixCapture(object):
    def __init__(self):
        self.parts = []

    def write(self, text):
        self.parts.append(text)

    def flush(self):
        pass

    def text(self):
        return "".join([E_TEXT(part) for part in self.parts])


# -----------------------------------------------------------------------------
# F. The faster CSV writer gives the same bytes
# -----------------------------------------------------------------------------

def fix_reference_csv(E, path, headers, rows):
    """write_csv as it was before the fast path (one csv_escape per cell)."""
    handle = codecs.open(path, "w", "utf-8")
    try:
        handle.write(u",".join([E.csv_escape(h) for h in headers]) + u"\n")
        for row in rows:
            handle.write(u",".join([E.csv_escape(row.get(h, "")) for h in headers]) + u"\n")
    finally:
        handle.close()


def fix_table(count):
    random = Random(11)
    micro = b"\xc2\xb5".decode("utf-8")
    specials = [float("nan"), float("inf"), -float("inf"), -0.0, 1e16, 1e-05, 0.1, 2.5e-310, 1.0 / 3.0, 35.0]
    texts = [u"plain", u'with "quotes"', micro + u"m", u"line\nbreak", u"", u"comma, inside", b"bytes".decode("ascii")]
    headers = ["F" + str(index) for index in range(12)] + ["I", "L", "B", "N", "T", "S", "R", "Missing"]
    rows = []
    for index in range(count):
        row = {}
        for column in range(12):
            row["F" + str(column)] = random.nextGaussian() * 1000.0 if column > 1 else specials[(index + column) % len(specials)]
        row["I"] = random.nextInt(100000) - 50000
        row["L"] = 2 ** 40 + index
        row["B"] = (index % 3 == 0)
        row["N"] = None
        row["T"] = texts[index % len(texts)]
        row["S"] = "byte text " + str(index)
        row["R"] = index * 0.5
        rows.append(row)
    return (headers, rows)


def test_fix_csv_fast_path(E, t, tmp):
    """write_csv and CsvAppender give byte for byte the file of the plain
    writer, on a large table of every kind of value."""
    headers, rows = fix_table(12000)
    reference = os.path.join(tmp, "fix_reference.csv")
    fast = os.path.join(tmp, "fix_fast.csv")
    appended = os.path.join(tmp, "fix_appended.csv")
    fix_reference_csv(E, reference, headers, rows)
    E.write_csv(fast, headers, rows)
    appender = E.CsvAppender(appended, headers)
    appender.append(rows[:5000])
    appender.append(rows[5000:5001])
    appender.append(rows[5001:])
    appender.close()
    expected = fix_read_bytes(reference)
    t.check(len(expected) > 1000000, "the table is large (" + str(len(expected)) + " bytes)")
    t.check(fix_read_bytes(fast) == expected, "write_csv writes the same bytes")
    t.check(fix_read_bytes(appended) == expected, "CsvAppender writes the same bytes")
    t.equal(appender.count, len(rows), "CsvAppender counts every row")
    for value in [float("nan"), -0.0, 1e16, 7, 2 ** 70, True, None, u'a"b', "x"]:
        t.equal(E.csv_cell(value), E.csv_escape(value), "csv_cell equals csv_escape for " + repr(value))


# -----------------------------------------------------------------------------
# G. Long image names and Windows path lengths
# -----------------------------------------------------------------------------

FIX_LONG_STEM = (
    "Day60_organoid-B-SOX2-green-PAX6-red-HOPX-magenta-DAPI-blue-zstack-zoom1-10x-org1-ROI1_Airyscan_Proc"
)[:98]


def fix_deep_root(tmp, length, name):
    """A folder whose absolute path has about length characters."""
    root = os.path.join(os.path.abspath(tmp), name)
    while len(root) < length - 1:
        root = os.path.join(root, "d" * min(60, max(1, length - len(root) - 2)))
    return root


def fix_all_paths(folder):
    found = []
    for root, folders, names in os.walk(folder):
        for name in names:
            found.append(File(os.path.join(root, name)).getAbsolutePath())
    return found


def test_fix_short_names_unchanged(E, t, tmp):
    """A name that fits keeps exactly the folder and file names of 2.2."""
    imp = make_image("short name.tif", slices=1)
    configs = [nuclei_config(E), marker_config(E)]
    plan = E.plan_run_folder(tmp, "short name", "", configs, multi_frame=False)
    t.check(plan["name"].startswith("short_name_CellQuant_") and not plan["shortened"], "the folder name is the full name")
    t.equal(len(plan["name"]), len("short_name_CellQuant_") + 15, "followed by the time stamp only")
    result = fix_analyze(E, imp, configs, [E.whole_image_roi(imp)], E.unique_run_dir(tmp, "short name", configs))
    t.equal(result["file_base"], "short_name", "files are named with the full name")
    names = sorted([os.path.relpath(path, result["run_dir"]).replace("\\", "/") for path in fix_all_paths(result["run_dir"])])
    t.equal(names, sorted([
        "Masks/short_name_C1_DAPI_FINAL_MASK.tif", "Masks/short_name_C1_DAPI_RAW_MASK.tif",
        "Masks/short_name_C2_MAP2_FINAL_MASK.tif", "Masks/short_name_C2_MAP2_RAW_MASK.tif",
        "Overlays/short_name_F0001_C1_DAPI_Objects_Overlay.tif",
        "Prepared/short_name_C1_DAPI_PREP.tif", "Prepared/short_name_C2_MAP2_PREP.tif",
        "ROIs/Parent_ROIs_Input.zip", "ROIs/short_name_F0001_C1_DAPI_Objects.zip",
        "objects.csv", "roi_summary.csv", "settings.properties", "thresholds.csv"
    ]), "every file has the name 2.2 gives it")
    t.equal(E.plan_run_folder(tmp, "x", "Validation")["name"][:22], "x_CellQuant_Validation", "validation folders keep their pattern")
    imp.close()


def fix_long_run(E, t, root, imp, configs, label):
    """Analyze imp into a new run folder below root; check that every path
    fits and that the image files use the folder's short form. Returns
    (short form, log lines, result)."""
    lines = []
    run_dir = E.unique_run_dir(root, FIX_LONG_STEM, configs, multi_frame=False)
    folder = os.path.basename(run_dir)
    short = folder[:folder.index("_CellQuant_")]
    t.check(short != FIX_LONG_STEM and short.endswith("_" + E.name_code(FIX_LONG_STEM)) and len(short) <= 48,
            label + ": the run folder uses a short form ending in the code of the full name (" + short + ")")
    result = fix_analyze(E, imp, configs, [E.whole_image_roi(imp)], run_dir, None, lines.append)
    paths = fix_all_paths(run_dir)
    t.equal(len(paths), 13, label + ": every result file is written")
    t.check(max([len(path) for path in paths]) <= 259, label + ": every path has at most 259 characters (longest " +
            str(max([len(path) for path in paths])) + ")")
    inside = [os.path.basename(path) for path in paths if os.path.basename(path).startswith(short + "_")]
    t.equal(len(inside), 8, label + ": every image file in the run folder starts with the short form")
    t.check(E.shortened_name_note(FIX_LONG_STEM, short) in lines, label + ": the log pairs the short form with the full name")
    t.check(len([line for line in lines if line.startswith("NOTE - the longest result path")]) == 0, label + ": no long-path note")
    t.equal(result["summary_rows"][0]["Source_Image"], FIX_LONG_STEM + ".tif", label + ": Source_Image keeps the full name")
    t.equal(result["file_base"], short, label + ": the run reports its short form")
    return (short, lines, result)


def test_fix_long_names_shortened(E, t, tmp):
    """A 98-character image name: the run folder and the files in it use the
    first part of the name and a code, every path fits in 259 characters,
    and the log records the full name. Before the fix: below a
    100-character folder the longest path had about 300 characters."""
    t.equal(len(FIX_LONG_STEM), 98, "the test name has 98 characters")
    short = E.short_name_base(FIX_LONG_STEM)
    t.equal(len(short), 48, "the usual short form has 48 characters")
    t.check(short.endswith("_" + E.name_code(FIX_LONG_STEM)), "and ends with the code of the full name")
    imp = make_image(FIX_LONG_STEM + ".tif", slices=1)
    configs = [nuclei_config(E), marker_config(E)]
    plain = fix_long_run(E, t, fix_deep_root(tmp, 60, "fix_long_60"), imp, configs, "60-character folder")
    t.equal(plain[0], short, "below a 60-character folder the usual 48-character form is enough")
    deeper = fix_long_run(E, t, fix_deep_root(tmp, 100, "fix_long_100"), imp, configs, "100-character folder")
    t.check(len(deeper[0]) < 48, "below a 100-character folder a shorter form is used (" + str(len(deeper[0])) + ")")
    root = fix_deep_root(tmp, 100, "fix_long_100")

    # Only files that are written count: without saved images the full name fits.
    nothing_saved = fix_options(E, save_prepared=False, save_masks=False, save_object_rois=False, save_overlays=False,
                                tissue_save_qc=False)
    shallow = fix_deep_root(tmp, 60, "fix_long_60")
    t.check(E.plan_run_folder(shallow, FIX_LONG_STEM, "", configs, multi_frame=False)["shortened"],
            "with images saved the name is shortened below a 60-character folder")
    kept = E.plan_run_folder(shallow, FIX_LONG_STEM, "", configs, multi_frame=False, options=nothing_saved)
    t.check(not kept["shortened"] and kept["name"].startswith(FIX_LONG_STEM + "_CellQuant_"),
            "without saved images the full name is kept there")
    t.check(kept["longest"] <= 259, "and the longest expected path fits (" + str(kept["longest"]) + ")")

    # Two long names that share their first 60 characters stay apart.
    other = FIX_LONG_STEM[:60] + "-second-section-of-the-same-organoid-imaged-later"
    first = E.plan_run_folder(root, FIX_LONG_STEM, "", configs)
    second = E.plan_run_folder(root, other, "", configs)
    t.check(second["shortened"] and second["file_base"] != first["file_base"], "a second long name gets another short form")

    # Validation in such a folder names its files the same way.
    lines = []
    validation_dir = E.unique_run_dir(root, FIX_LONG_STEM, [configs[0]], "Validation", multi_frame=False)
    folder = os.path.basename(validation_dir)
    validation_short = folder[:folder.index("_CellQuant_Validation_")]
    t.check(validation_short.endswith("_" + E.name_code(FIX_LONG_STEM)), "validation folders are shortened too")
    prepared = E.AnalysisEngine().prepare_single_channel(imp, configs[0], 1)
    E.validate_objects_against_points(E.AnalysisEngine(log=lines.append), prepared, configs[0], [E.whole_image_roi(imp)],
                                      [(30, 30), (5, 5)], validation_dir, imp.getTitle(), 1)
    names = [os.path.basename(path) for path in fix_all_paths(validation_dir)]
    t.equal(len([name for name in names if name.startswith(validation_short + "_")]), 2,
            "the validation overlay and mask use the short form")
    t.check(max([len(path) for path in fix_all_paths(validation_dir)]) <= 259, "and fit")
    t.check(E.shortened_name_note(FIX_LONG_STEM, validation_short) in lines, "and the validation log records the full name")
    prepared.close()

    # Deeper still: even the shortest form does not fit, so the run says so.
    if os.sep == "\\":
        deep = fix_deep_root(tmp, 236, "fix_deep")
        deep_lines = []
        fix_analyze(E, imp, configs, [E.whole_image_roi(imp)], E.unique_run_dir(deep, FIX_LONG_STEM, configs),
                    None, deep_lines.append)
        notes = [line for line in deep_lines if line.startswith("NOTE - the longest result path")]
        t.equal(len(notes), 1, "a single run logs the long-path note")
        t.check(len(notes) == 1 and "Python or R scripts" in notes[0], "which names Python and R scripts")
    imp.close()


def test_fix_long_names_in_a_batch(E, t, tmp):
    """In a batch the run folder is shortened, the input-side ROI file keeps
    the full name and is found, and batch_log.txt records the mapping."""
    root = fix_deep_root(tmp, 90, "fix_long_batch")
    folder = os.path.join(root, "in")
    os.makedirs(folder)
    imp = make_image(FIX_LONG_STEM + ".tif", slices=1)
    FileSaver(imp).saveAsTiff(os.path.join(folder, FIX_LONG_STEM + ".tif"))
    oval = OvalRoi(20, 20, 150, 120)
    oval.setName("organoid")
    E.save_roi_array([oval], os.path.join(folder, FIX_LONG_STEM + "_ROIs.zip"))
    configs = [nuclei_config(E), marker_config(E)]
    batch = E.run_batch(folder, os.path.join(root, "out"), configs, configs, fix_options(E),
                        {"extensions": "tif", "roi_policy": E.BATCH_ROI_FILE_SKIP})
    report = fix_csv_rows(os.path.join(batch["batch_dir"], "batch_report.csv"))[0]
    t.equal(report["Status"], "OK", "the image is analyzed")
    t.check(report["ROI_Source"].endswith(FIX_LONG_STEM + "_ROIs.zip") or FIX_LONG_STEM in report["ROI_Source"],
            "with the ROI file of the full name")
    t.equal(report["Image_File"], FIX_LONG_STEM + ".tif", "batch_report.csv names the image in full")
    folder = report["Run_Folder"]
    short = folder[:folder.index("_CellQuant_")]
    t.check(short.endswith("_" + E.name_code(FIX_LONG_STEM)) and len(short) <= 48, "and its run folder in short (" + short + ")")
    log = read_utf8(os.path.join(batch["batch_dir"], "batch_log.txt"))
    t.check(E.shortened_name_note(FIX_LONG_STEM, short) in log, "batch_log.txt records the full name")
    t.check(batch["longest_path"] <= 259, "every batch path fits (" + str(batch["longest_path"]) + ")")
    run_log = read_utf8(os.path.join(batch["batch_dir"], report["Run_Folder"], "run_log.txt"))
    t.check("Long image name shortened" in run_log, "and so does the image's run_log.txt")
    imp.close()


# -----------------------------------------------------------------------------
# H. Profile checks
# -----------------------------------------------------------------------------

def fix_profile(E, tmp, name, changes, removals=None):
    path = os.path.join(tmp, name)
    configs = [nuclei_config(E), marker_config(E)]
    E.write_settings_properties(path, E.default_global_options(), configs)
    props = E.read_properties_file(path)
    for key in changes:
        props.setProperty(key, changes[key])
    for key in (removals or []):
        props.remove(key)
    from java.io import FileOutputStream
    stream = FileOutputStream(path)
    try:
        props.store(stream, "test")
    finally:
        stream.close()
    return path


def test_fix_profile_checks(E, t, tmp):
    """Before the fixes: an unknown frame mode silently analyzed frame 1 only;
    a bad value in a switched-off channel blocked the whole profile; an
    unknown automatic method was found only in the middle of a run."""
    path = fix_profile(E, tmp, "fix_frames.properties", {"global.frame_mode": "Every frame"})
    t.check(fix_raises(lambda: E.load_profile_settings(path), "frame mode"), "an unknown frame mode is refused at load")
    t.check(fix_raises(lambda: E.run_batch(tmp, tmp, [nuclei_config(E)], [nuclei_config(E)],
                                           fix_options(E, frame_mode="Every frame"), {"extensions": "tif"}), "frame mode"),
            "and by a batch")
    path = fix_profile(E, tmp, "fix_off.properties", {
        "channel.6.enabled": "false", "channel.6.analysis_mode": "Ignore", "channel.6.manual_low": "",
        "channel.6.name": "Spare", "channel.6.auto_method": "Otsuu"
    })
    options, enabled, everything = E.load_profile_settings(path)
    spare = [config for config in everything if config["channel"] == 6]
    t.equal(len(enabled), 2, "a bad value in a switched-off channel does not block the profile")
    t.equal((spare[0]["manual_low"], spare[0]["auto_method"], spare[0]["name"]), (500.0, "Otsu", "Spare"),
            "the bad values get their defaults, the others are kept")
    path = fix_profile(E, tmp, "fix_on.properties", {"channel.1.manual_low": ""})
    t.check(fix_raises(lambda: E.load_profile_settings(path), "numeric"), "a bad value in a switched-on channel is refused")
    path = fix_profile(E, tmp, "fix_method.properties", {"channel.1.threshold_mode": "Auto", "channel.1.auto_method": "Otsuu"})
    t.check(fix_raises(lambda: E.load_profile_settings(path), "Otsuu"), "an unknown automatic method is refused at load")
    path = fix_profile(E, tmp, "fix_case.properties", {"channel.1.threshold_mode": "Auto", "channel.1.auto_method": "otsu"})
    t.check(fix_raises(lambda: E.load_profile_settings(path), "spelled exactly"), "the method must be spelled exactly")
    for name in ["CellQuant_Tuning_20261119_101500", "CellQuant_Robustness_20261119_101500_2",
                 "roi_summary_CellQuant_CountCorrection_20261119_101500", "CellQuant_Batch_20261119_101500",
                 "image_CellQuant_20261119_101500", "CellQuant_Anything_20261119_101500"]:
        t.check(E.is_cellquant_output_dir(name), name + " is a result folder")
    for name in ["CellQuant data", "organoids", "CellQuant_notes", "CellQuant_Tuning_old"]:
        t.check(not E.is_cellquant_output_dir(name), name + " is not a result folder")


# -----------------------------------------------------------------------------
# I. Module registration, batch ROI providers, positive fraction
# -----------------------------------------------------------------------------

FIX_MODULE_SOURCE = (
    "register_global_settings([('probe.key', 'probe_key', 'text', 'x')])\n"
    "class FixProbeExtension(RunExtension):\n"
    "    name = 'fix_probe'\n"
    "register_run_extension(FixProbeExtension())\n"
    "register_batch_roi_policy('Fix probe policy', lambda request: (None, 'probe'))\n"
    "register_settings_form('Fix probe', 'note', [('probe_key', 'Probe', 'text', None)])\n"
    "FIX_PROBE_NAME = 1\n"
)


def fix_registries(namespace):
    return [list(namespace[name]) if isinstance(namespace[name], list) else dict(namespace[name])
            for name in ["GLOBAL_SETTING_SPECS", "RUN_EXTENSIONS", "BATCH_ROI_POLICIES",
                         "BATCH_ROI_PROVIDERS", "SETTINGS_FORMS", "ENGINE_MODULE_VERSIONS"]]


def test_fix_refused_module_leaves_nothing(E, t, tmp):
    """A module of another release (or one that fails while loading) is
    refused and everything it registered is taken back. Before the fix: its
    settings, extension, ROI policy and form stayed registered."""
    cases = [
        ("another release", FIX_MODULE_SOURCE + "register_engine_module('cellquanttissue.py', '0.0.1')\n", RuntimeError),
        ("a failure while loading", FIX_MODULE_SOURCE + "raise ValueError('broken module')\n", ValueError)
    ]
    for label, source, error in cases:
        namespace = {"CELLQUANT_NO_GUI": True}
        run_source_file(find_engine(), namespace)
        folder = os.path.join(tmp, "fix_module_" + label.replace(" ", "_"))
        os.makedirs(folder)
        handle = open(os.path.join(folder, "cellquanttissue.py"), "w")
        try:
            handle.write(source)
        finally:
            handle.close()
        before = fix_registries(namespace)
        names = set(namespace.keys())
        refused = False
        try:
            namespace["load_engine_modules"](namespace, folder)
        except error:
            refused = True
        t.check(refused, label + ": the module is refused")
        t.check(fix_registries(namespace) == before, label + ": every registry is as before")
        t.check(set(namespace.keys()) == names, label + ": no name of the module is left")
    namespace = {"CELLQUANT_NO_GUI": True}
    run_source_file(find_engine(), namespace)
    folder = os.path.join(tmp, "fix_module_good")
    os.makedirs(folder)
    handle = open(os.path.join(folder, "cellquanttissue.py"), "w")
    try:
        handle.write(FIX_MODULE_SOURCE + "register_engine_module('cellquanttissue.py', APP_VERSION)\n")
    finally:
        handle.close()
    t.equal(namespace["load_engine_modules"](namespace, folder), ["cellquanttissue.py"], "a module of this release loads")
    t.check("Fix probe policy" in namespace["BATCH_ROI_POLICIES"] and "FIX_PROBE_NAME" in namespace,
            "and stays registered")


def test_fix_roi_provider_can_cancel(E, t, tmp):
    """A batch ROI policy of a module receives the batch's cancel check.
    Before the fix: the request had no "is_cancelled"."""
    folder = os.path.join(tmp, "fix_provider")
    os.makedirs(folder)
    imp = make_image("one.tif", slices=1)
    FileSaver(imp).saveAsTiff(os.path.join(folder, "one.tif"))
    imp.close()
    seen = {}

    def provider(request):
        seen.update(request)
        return (None, "probe: skipped")

    def never():
        return False

    name = "Fix probe ROI policy"
    E.register_batch_roi_policy(name, provider)
    try:
        configs = [nuclei_config(E)]
        batch = E.run_batch(folder, os.path.join(tmp, "fix_provider_out"), configs, configs, fix_options(E),
                            {"extensions": "tif", "roi_policy": name}, is_cancelled=never)
    finally:
        E.BATCH_ROI_POLICIES.remove(name)
        del E.BATCH_ROI_PROVIDERS[name]
    t.check(seen.get("is_cancelled") is never, "the policy gets the batch's cancel check")
    t.equal(batch["skipped"], 1, "and its answer is used")


def test_fix_positive_fraction_own_area(E, t, tmp):
    """The positive-area fraction of a channel is taken over that channel's
    own ROI area. Identical when no pixel is NaN; on a 32-bit image whose
    second channel has NaN in 20 of 100 ROI pixels, 40 positive pixels are
    50 percent (40 of 80), not 40 percent as before the fix."""
    first = FloatProcessor(40, 40)
    first.setValue(100.0)
    first.fill()
    second = FloatProcessor(40, 40)
    for index in range(40):
        second.setf(10 + index % 10, 10 + index // 10, 1000.0)
    for index in range(20):
        second.setf(10 + index % 10, 16 + index // 10, float("nan"))
    imp = fix_stack_image("nan.tif", [first, second])
    c1 = E.normalize_channel_config(1, {"enabled": True, "name": "A", "analysis_mode": "Intensity only"})
    c2 = E.normalize_channel_config(2, {"enabled": True, "name": "B", "analysis_mode": "Positive area",
                                        "threshold_mode": "Manual", "manual_low": 500})
    roi = Roi(10, 10, 10, 10)
    row = fix_analyze(E, imp, [c1, c2], [roi], E.unique_run_dir(tmp, "fix_nan"))["summary_rows"][0]
    t.equal(row["ROI_Area"], 100.0, "ROI_Area is still that of the first channel")
    t.equal(row["C2_B_Positive_Area"], 40.0, "40 pixels are positive")
    t.equal(row["C2_B_Positive_Area_Fraction_Percent"], 50.0, "which is half of the 80 pixels that hold a number")
    imp.close()

    plain = make_image("plain.tif", slices=1)
    configs = [nuclei_config(E), marker_config(E)]
    result = fix_analyze(E, plain, configs, [E.whole_image_roi(plain), OvalRoi(20, 30, 150, 120)], E.unique_run_dir(tmp, "fix_plain"))
    for row in result["summary_rows"]:
        for prefix in ["C1_DAPI", "C2_MAP2"]:
            t.equal(row[prefix + "_Positive_Area_Fraction_Percent"],
                    row[prefix + "_Positive_Area"] / row["ROI_Area"] * 100.0,
                    prefix + ": without NaN the fraction is exactly the one of 2.2 (positive area / ROI_Area)")
    plain.close()


TESTS.append(("engine fix: an empty threshold range gives zero, not an error", test_fix_empty_threshold_range))
TESTS.append(("engine fix: sensitivity factor above every pixel gives zero with a note", test_fix_sensitivity_empty_range))
TESTS.append(("engine fix: a Java error in one sensitivity factor is noted, not fatal", test_fix_sensitivity_java_error))
TESTS.append(("engine fix: images with a calibration function are refused for segmentation", test_fix_calibration_function_refused))
TESTS.append(("engine fix: per-object intensity of a ring and the object in its hole", test_fix_object_intensity_ring))
TESTS.append(("engine fix: per-object intensity bit-identical without holes", test_fix_object_intensity_bit_identical))
TESTS.append(("engine fix: object outlines without a ROI Manager, identical", test_fix_object_rois_without_roi_manager))
TESTS.append(("engine fix: a measurement redirect neither breaks nor changes", test_fix_measurement_redirect))
TESTS.append(("engine fix: Local Phansalkar called directly, pixel-identical", test_fix_local_threshold_direct))
TESTS.append(("engine fix: RGB lookup table and inverted LUT with Auto", test_fix_rgb_lut_and_inverted_auto))
TESTS.append(("engine fix: an unknown ImageJ command is refused", test_fix_unknown_command_refused))
TESTS.append(("engine fix: float settings keep every digit", test_fix_float_settings_text))
TESTS.append(("engine fix: the test runner reports broken files and changed state", test_fix_runner_reports_broken_files))
TESTS.append(("engine fix: fast CSV writer gives the same bytes", test_fix_csv_fast_path))
TESTS.append(("engine fix: short names stay as in 2.2", test_fix_short_names_unchanged))
TESTS.append(("engine fix: long image names are shortened to fit Windows paths", test_fix_long_names_shortened))
TESTS.append(("engine fix: long image names in a batch", test_fix_long_names_in_a_batch))
TESTS.append(("engine fix: profile checks (frame mode, switched-off channels, method, folders)", test_fix_profile_checks))
TESTS.append(("engine fix: a refused module leaves nothing registered", test_fix_refused_module_leaves_nothing))
TESTS.append(("engine fix: batch ROI policies receive the cancel check", test_fix_roi_provider_can_cancel))
TESTS.append(("engine fix: positive fraction over the channel's own area", test_fix_positive_fraction_own_area))

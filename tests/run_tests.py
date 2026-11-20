# Cell Quant Studio Universal - automated tests for the analysis engine
#
# What this checks: the measurement code in cellquantengine.py, on small
# synthetic images whose correct answers are known by construction (how many
# objects, how many positive pixels, which objects lie inside an ROI).
#
# How to run it in Fiji:
#   Open this file in the Script Editor, choose Language > Python (Jython),
#   and click Run. A summary is printed to the editor console and the Log.
# Headless:
#   ImageJ-linux64 --ij2 --headless --console --run tests/run_tests.py
#
# The tests open no window, need no data files, write only to a temporary
# folder that is removed afterwards, and restore every Fiji option they touch.
# Tests that need the Auto Local Threshold or CLAHE plugin are reported as
# SKIPPED when the plugin is not installed.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij import IJ, ImagePlus, ImageStack, CompositeImage, Prefs
from ij.process import ShortProcessor, ByteProcessor, ImageProcessor, ImageConverter
from ij.gui import Roi, OvalRoi, PolygonRoi, PointRoi, ShapeRoi
from ij.io import FileSaver
from ij.plugin.frame import RoiManager
from java.awt import Polygon
from java.lang import Throwable as JavaThrowable

import os
import math
import shutil
import tempfile
import traceback
import codecs


# Under Jython, "except ANY_ERROR" alone does not catch Java exceptions.
ANY_ERROR = (Exception, JavaThrowable)

MICRO = b"\xc2\xb5".decode("utf-8")
WIDTH = 240
HEIGHT = 200
PIXEL = 0.5
DISK_RADIUS = 6
DISKS = [
    (30, 30), (70, 30), (110, 30), (150, 30), (190, 30),
    (30, 80), (70, 80), (110, 80), (150, 80), (190, 80),
    (50, 140), (100, 140), (150, 140), (200, 140), (120, 175)
]
MARKER_OUTER = (40, 50, 100, 80)
MARKER_INNER = (60, 70, 60, 40)


class SkipTest(Exception):
    pass


class Namespace(object):
    def __init__(self, values):
        self.__dict__.update(values)


class Recorder(object):
    def __init__(self, name):
        self.name = name
        self.checks = 0
        self.failures = []

    def check(self, condition, message):
        self.checks += 1
        if not condition:
            self.failures.append(message)

    def equal(self, actual, expected, message):
        self.check(actual == expected, message + " (expected " + repr(expected) + ", got " + repr(actual) + ")")

    def close(self, actual, expected, message, tolerance=1e-6):
        ok = False
        try:
            ok = abs(float(actual) - float(expected)) <= tolerance * max(1.0, abs(float(expected)))
        except ANY_ERROR:
            ok = False
        self.check(ok, message + " (expected " + repr(expected) + ", got " + repr(actual) + ")")


# -----------------------------------------------------------------------------
# Loading the engine
# -----------------------------------------------------------------------------

def find_engine():
    candidates = []
    if "CELLQUANT_ENGINE_PATH" in globals():
        candidates.append(globals()["CELLQUANT_ENGINE_PATH"])
    if os.environ.get("CELLQUANT_ENGINE_PATH"):
        candidates.append(os.environ.get("CELLQUANT_ENGINE_PATH"))
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        candidates.append(os.path.join(os.path.dirname(here), "cellquantengine.py"))
        candidates.append(os.path.join(here, "cellquantengine.py"))
    except ANY_ERROR:
        pass
    try:
        root = IJ.getDirectory("imagej")
        if root:
            candidates.append(os.path.join(root, "plugins", "Scripts", "Cell Quant Studio", "cellquantengine.py"))
    except ANY_ERROR:
        pass
    for candidate in candidates:
        if candidate and os.path.isfile(candidate):
            return candidate
    raise IOError(
        "Could not find cellquantengine.py. Keep the tests folder inside the Cell Quant Studio folder, "
        "or install the application in Fiji.app/plugins/Scripts/Cell Quant Studio/."
    )


def load_engine():
    path = find_engine()
    namespace = {"CELLQUANT_NO_GUI": True}
    try:
        runner = execfile
    except NameError:
        runner = None
    if runner is not None:
        runner(path, namespace)
    else:
        handle = open(path, "r")
        try:
            source = handle.read()
        finally:
            handle.close()
        exec(compile(source, path, "exec"), namespace)
    return Namespace(namespace)


# -----------------------------------------------------------------------------
# Synthetic images with known answers
# -----------------------------------------------------------------------------

def blank_plane(value):
    ip = ShortProcessor(WIDTH, HEIGHT)
    ip.setValue(float(value))
    ip.fill()
    return ip


def draw_disk(ip, cx, cy, radius, value):
    ip.setValue(float(value))
    ip.fill(OvalRoi(cx - radius, cy - radius, 2 * radius + 1, 2 * radius + 1))
    ip.resetRoi()


def draw_rect(ip, rect, value):
    ip.setValue(float(value))
    ip.fill(Roi(rect[0], rect[1], rect[2], rect[3]))
    ip.resetRoi()


def nuclei_plane(disks):
    ip = blank_plane(100)
    for disk in disks:
        draw_disk(ip, disk[0], disk[1], DISK_RADIUS, 2000)
    return ip


def marker_plane():
    ip = blank_plane(100)
    draw_rect(ip, MARKER_OUTER, 900)
    draw_rect(ip, MARKER_INNER, 1500)
    return ip


def disk_pixel_count():
    ip = ByteProcessor(40, 40)
    ip.setValue(255.0)
    ip.fill(OvalRoi(20 - DISK_RADIUS, 20 - DISK_RADIUS, 2 * DISK_RADIUS + 1, 2 * DISK_RADIUS + 1))
    return ip.getHistogram()[255]


def make_image(title, slices=2, frames=1, disks_per_frame=None, channels=2):
    """Two channels: C1 holds bright disks ('nuclei'), C2 a two-level
    rectangle ('marker'). Calibrated at 0.5 micrometer per pixel."""
    stack = ImageStack(WIDTH, HEIGHT)
    for t in range(frames):
        disks = DISKS if disks_per_frame is None else disks_per_frame[t]
        for z in range(slices):
            for c in range(channels):
                stack.addSlice("", nuclei_plane(disks) if c == 0 else marker_plane())
    imp = ImagePlus(title, stack)
    imp.setDimensions(channels, slices, frames)
    if channels > 1:
        imp = CompositeImage(imp, CompositeImage.GRAYSCALE)
    elif slices * frames > 1:
        imp.setOpenAsHyperStack(True)
    cal = imp.getCalibration()
    cal.pixelWidth = PIXEL
    cal.pixelHeight = PIXEL
    cal.setUnit("um")
    return imp


def nuclei_config(E, overrides=None):
    values = {
        "enabled": True, "name": "DAPI", "analysis_mode": "Objects + area",
        "threshold_mode": "Manual", "manual_low": 1000, "manual_high": -1
    }
    if overrides:
        values.update(overrides)
    return E.normalize_channel_config(1, values)


def marker_config(E, overrides=None):
    values = {
        "enabled": True, "name": "MAP2", "analysis_mode": "Positive area",
        "threshold_mode": "Manual", "manual_low": 800, "manual_high": -1
    }
    if overrides:
        values.update(overrides)
    return E.normalize_channel_config(2, values)


def triangle_roi():
    polygon = Polygon()
    polygon.addPoint(100, 30)
    polygon.addPoint(170, 170)
    polygon.addPoint(30, 170)
    roi = PolygonRoi(polygon, Roi.POLYGON)
    roi.setName("triangle")
    return roi


def mask_with_disks(centers):
    ip = ByteProcessor(200, 200)
    ip.setValue(255.0)
    for center in centers:
        ip.fill(OvalRoi(center[0] - DISK_RADIUS, center[1] - DISK_RADIUS, 2 * DISK_RADIUS + 1, 2 * DISK_RADIUS + 1))
    ip.resetRoi()
    return ImagePlus("mask", ip)


def read_utf8(path):
    handle = codecs.open(path, "r", "utf-8")
    try:
        return handle.read()
    finally:
        handle.close()


def run_image(E, imp, configs, rois, out_root, frames=None, options_override=None):
    global_options = E.default_global_options()
    if options_override:
        global_options.update(options_override)
    run_dir = E.unique_run_dir(out_root, E.strip_extension(imp.getTitle()))
    engine = E.AnalysisEngine()
    return engine.analyze_image(
        imp, "", configs, rois, E.engine_options_from_global(global_options), run_dir,
        frames if frames is not None else [1], settings=(global_options, configs)
    )


def local_threshold_available(E):
    probe = ImagePlus("probe", nuclei_plane(DISKS))
    config = nuclei_config(E, {"threshold_mode": "Local Phansalkar"})
    try:
        data = E.AnalysisEngine().build_masks(probe, config, "probe")
    except RuntimeError:
        return False
    E.safe_close_image(data["raw"])
    E.safe_close_image(data["final"])
    return True


def clahe_available():
    try:
        from mpicbg.ij.clahe import Flat
        return Flat is not None
    except ANY_ERROR:
        return False


# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

def test_text_and_csv(E, t, tmp):
    name = "caf" + b"\xc3\xa9".decode("utf-8") + " " + MICRO + "m"
    t.equal(E.to_text(None), "", "to_text(None)")
    t.equal(E.to_text(name), name, "to_text keeps non-ASCII text")
    t.check(MICRO in E.sanitize_name(name), "sanitize_name keeps non-ASCII characters")
    t.check(" " not in E.sanitize_name(name), "sanitize_name replaces spaces")
    t.check(E.is_micrometer_unit(MICRO + "m"), "micro sign unit recognised")
    t.check(E.is_micrometer_unit("micron"), "'micron' recognised")
    path = os.path.join(tmp, "unicode.csv")
    E.write_csv(path, ["Name", "Value", "Missing"], [{"Name": name, "Value": 1.5}, {"Name": 'quote"d', "Value": float("nan")}])
    text = read_utf8(path)
    t.check(name in text, "CSV holds the non-ASCII text as UTF-8")
    t.check('"1.5"' in text, "CSV holds the number")
    t.check('"quote""d"' in text, "CSV doubles embedded quotes")
    t.equal(len(text.strip().split("\n")), 3, "CSV has a header and two rows")
    t.equal(E.parse_factor_list("0.8, 1.2"), [0.8, 1.0, 1.2], "factor list always includes the baseline 1.0")


def test_settings_roundtrip(E, t, tmp):
    c1 = nuclei_config(E, {"open_iterations": 2, "morph_count": 3, "fill_holes": True, "name": "DAPI " + MICRO})
    c2 = marker_config(E)
    options = E.default_global_options()
    options["spatial_unit"] = MICRO + "m"
    options["batch_recursive"] = True
    path = os.path.join(tmp, "profile.properties")
    E.write_settings_properties(path, options, [c1, c2])
    loaded_options, enabled, everything = E.load_profile_settings(path)
    t.equal(len(enabled), 2, "both channels come back enabled")
    for key in c1:
        t.equal(enabled[0][key], c1[key], "C1 setting '" + key + "' survives a save and load")
    t.equal(loaded_options["spatial_unit"], MICRO + "m", "non-ASCII unit survives a save and load")
    t.equal(loaded_options["batch_recursive"], True, "batch option survives a save and load")

    # A profile written by 2.2.x has no neighbor count and no batch keys.
    legacy = E.read_properties_file(path)
    legacy.remove("channel.1.morph_count")
    legacy.remove("batch.recursive")
    t.equal(E.normalize_channel_config(1, E.channel_values_from_properties(legacy, 1))["morph_count"], 1,
            "missing neighbor count defaults to 1")
    t.equal(E.global_options_from_properties(legacy)["batch_recursive"], False, "missing batch key uses its default")


def test_config_validation(E, t, tmp):
    def raises(function):
        try:
            function()
        except ValueError:
            return True
        return False

    t.check(raises(lambda: E.normalize_channel_config(1, {"manual_low": "abc"})), "non-numeric threshold is rejected")
    t.check(raises(lambda: E.normalize_channel_config(1, {"morph_count": 9})), "neighbor count above 8 is rejected")
    t.check(raises(lambda: E.normalize_channel_config(1, {"analysis_mode": "Nonsense"})), "unknown analysis mode is rejected")
    imp = make_image("validate.tif", slices=1)
    good = [nuclei_config(E), marker_config(E)]
    E.validate_channel_configs(imp, good)
    t.check(raises(lambda: E.validate_channel_configs(imp, [nuclei_config(E), marker_config(E, {"name": "dapi"})])),
            "duplicate channel names are rejected")
    t.check(raises(lambda: E.validate_channel_configs(imp, [E.normalize_channel_config(3, {"enabled": True, "analysis_mode": "Intensity only"})])),
            "a channel the image does not have is rejected")
    t.check(raises(lambda: E.validate_channel_configs(imp, [nuclei_config(E, {"clahe": True})])),
            "CLAHE without a local threshold is rejected")
    t.check(raises(lambda: E.validate_unattended_configs([marker_config(E, {"threshold_mode": "Display minimum"})], "Batch processing")),
            "a display-minimum threshold is rejected for unattended runs")
    line = Roi(0, 0, 10, 10)
    t.check(raises(lambda: E.validate_parent_rois([Roi(200, 150, 100, 100)], WIDTH, HEIGHT)), "an ROI outside the image is rejected")
    t.equal(E.validate_parent_rois([line], WIDTH, HEIGHT)[0].getName(), "ROI_001", "an unnamed ROI gets a name")
    imp.close()


def test_full_run(E, t, tmp):
    imp = make_image("full " + MICRO + ".tif")
    configs = [nuclei_config(E), marker_config(E)]
    result = run_image(E, imp, configs, [E.whole_image_roi(imp)], tmp)
    row = result["summary_rows"][0]
    disk_area = disk_pixel_count() * PIXEL * PIXEL
    outer = MARKER_OUTER[2] * MARKER_OUTER[3]
    inner = MARKER_INNER[2] * MARKER_INNER[3]
    total = WIDTH * HEIGHT

    t.close(row["ROI_Area"], total * PIXEL * PIXEL, "ROI area is calibrated")
    t.equal(row["Spatial_Unit"], MICRO + "m", "unit is reported")
    t.equal(row["C1_DAPI_Object_Count"], len(DISKS), "every disk is counted once")
    t.close(row["C1_DAPI_Total_Object_Area"], len(DISKS) * disk_area, "total object area")
    t.close(row["C1_DAPI_Object_Density_Per_mm2"], len(DISKS) * 1000000.0 / (total * PIXEL * PIXEL), "density per mm2")
    t.close(row["C1_DAPI_Positive_Area"], len(DISKS) * disk_area, "nuclear positive area equals the disks")
    t.close(row["C2_MAP2_Positive_Area"], outer * PIXEL * PIXEL, "marker positive area")
    t.close(row["C2_MAP2_Positive_Area_Fraction_Percent"], 100.0 * outer / total, "marker positive-area fraction")
    expected_mean = (inner * 1500.0 + (outer - inner) * 900.0 + (total - outer) * 100.0) / total
    t.close(row["C2_MAP2_Mean"], expected_mean, "marker mean intensity")
    t.close(row["C2_MAP2_Raw_Integrated_Density"], expected_mean * total, "raw integrated density")
    t.equal(row["C1_DAPI_Edge_Touching_Object_Count"], 0, "no disk touches the image border")

    objects = result["object_rows"]
    t.equal(len(objects), len(DISKS), "one object row per disk")
    t.close(objects[0]["Object_Area"], disk_area, "object area")
    t.close(objects[0]["C1_DAPI_Object_Mean"], 2000.0, "object intensity on its own channel")
    t.equal(objects[0]["C2_MAP2_Object_Mean"], "", "no cross-channel measurement unless requested")
    # A disk centred on pixel (x, y) has its centroid at (x + 0.5, y + 0.5).
    centroids = sorted([(int(math.floor(o["Centroid_X"] / PIXEL + 0.001)), int(math.floor(o["Centroid_Y"] / PIXEL + 0.001))) for o in objects])
    t.equal(centroids, sorted(DISKS), "object centroids are the disk centres")

    run_dir = result["run_dir"]
    for name in ["roi_summary.csv", "objects.csv", "thresholds.csv", "settings.properties"]:
        t.check(os.path.isfile(os.path.join(run_dir, name)), name + " is written")
    t.check(os.path.isfile(os.path.join(run_dir, "ROIs", "Parent_ROIs_Input.zip")), "parent ROIs are saved")
    t.equal(len(os.listdir(os.path.join(run_dir, "Masks"))), 4, "raw and final masks are saved for both channels")
    t.equal(len(os.listdir(os.path.join(run_dir, "Prepared"))), 2, "prepared images are saved")
    t.equal(len(os.listdir(os.path.join(run_dir, "Overlays"))), 1, "object overlay is saved")
    saved = E.load_rois_from_file(os.path.join(run_dir, "ROIs", [n for n in os.listdir(os.path.join(run_dir, "ROIs")) if "Objects" in n][0]))
    t.equal(len(saved), len(DISKS), "every object ROI is in the saved ROI set")
    text = read_utf8(os.path.join(run_dir, "roi_summary.csv"))
    t.check(MICRO + "m" in text, "summary CSV is UTF-8 and holds the unit")
    thresholds = read_utf8(os.path.join(run_dir, "thresholds.csv"))
    t.check("Binary_Neighbor_Count" in thresholds and "ImageJ_Version" in thresholds, "thresholds.csv records the new provenance columns")

    crossed = run_image(E, imp, [nuclei_config(E, {"cross_measure": True}), marker_config(E)], [E.whole_image_roi(imp)], tmp)
    values = [o["C2_MAP2_Object_Mean"] for o in crossed["object_rows"]]
    t.check(len([v for v in values if v == ""]) == 0, "cross-channel measurement fills the other channel")
    imp.close()


def test_column_order_compatible(E, t, tmp):
    engine = E.AnalysisEngine()
    configs = [nuclei_config(E), marker_config(E)]
    summary = engine.build_summary_headers(configs)
    t.equal(summary[:12], [
        "Run_ID", "Application_Version", "Source_Image", "Source_Path", "Frame", "Parent_ROI_Index",
        "Parent_ROI_Name", "ROI_Area", "ROI_Area_Unit", "Pixel_Width", "Pixel_Height", "Spatial_Unit"
    ], "summary keeps the 2.2 leading columns")
    t.equal(summary[-2:], ["C1_DAPI_Edge_Touching_Object_Count", "C1_DAPI_Excluded_Edge_Object_Count"],
            "columns added in 2.3 come last in the summary")
    objects = engine.build_object_headers(configs)
    t.equal(objects[9:15], ["Object_Area", "Object_Perimeter", "Object_Circularity", "Object_Feret", "Centroid_X", "Centroid_Y"],
            "object table keeps the 2.2 geometry columns in place")
    t.equal(objects[-5:], ["Touches_ROI_Edge", "Object_MinFeret", "Object_Aspect_Ratio", "Object_Roundness", "Object_Solidity"],
            "columns added in 2.3 come last in the object table")
    thresholds = engine.build_threshold_headers()
    t.equal(thresholds[41], "Cross_Channel_Object_Measurement", "threshold table keeps the 2.2 columns in place")


def test_objects_confined_to_roi(E, t, tmp):
    """Objects must be confined to the true outline of a non-rectangular ROI."""
    inside = (100, 100)
    straddling = (60, 110)
    outside_in_box = [(40, 50), (160, 60)]
    far_outside = (10, 10)
    mask = mask_with_disks([inside, straddling, far_outside] + outside_in_box)
    triangle = triangle_roi()
    engine = E.AnalysisEngine()
    full = disk_pixel_count()

    objects, info = engine.analyze_objects(mask, triangle, nuclei_config(E))
    t.equal(len(objects), 2, "only objects inside the triangle are counted")
    areas = sorted([o["area"] for o in objects])
    t.close(areas[1], full, "the object inside the ROI keeps its full area")
    t.check(0 < areas[0] < full, "the object crossing the outline is clipped at the outline")
    expected = E.count_positive_pixels(mask.getProcessor().duplicate(), triangle)
    t.close(sum(areas), expected, "object area inside the ROI equals positive area inside the ROI")
    t.equal(info["edge_touching"], 1, "the clipped object is flagged as touching the ROI edge")
    t.equal(sorted([o["touches_edge"] for o in objects]), [False, True], "edge flag is per object")

    kept, info = engine.analyze_objects(mask, triangle, nuclei_config(E, {"exclude_edges": True}))
    t.equal(len(kept), 1, "edge exclusion removes the object on the outline")
    t.equal(info["excluded"], 1, "the excluded object is counted")
    t.close(kept[0]["area"], full, "the interior object remains")

    whole = Roi(0, 0, 200, 200)
    objects, info = engine.analyze_objects(mask, whole, nuclei_config(E))
    t.equal(len(objects), 5, "all objects are counted in the whole image")
    t.equal(info["edge_touching"], 0, "no object touches the image border here")

    border = mask_with_disks([(3, 100), (100, 100)])
    objects, info = engine.analyze_objects(border, whole, nuclei_config(E, {"exclude_edges": True}))
    t.equal(len(objects), 1, "an object on the image border is excluded")

    # For any ROI shape, the detected objects must cover exactly the mask
    # pixels that ImageJ measures as positive inside that ROI.
    field = mask_with_disks([(20 + 23 * (i % 8), 18 + 27 * (i // 8)) for i in range(48)])
    ring = ShapeRoi(OvalRoi(15, 10, 170, 175)).xor(ShapeRoi(OvalRoi(70, 60, 60, 70)))
    slanted = Polygon()
    for corner in [(12, 40), (150, 8), (192, 120), (95, 190), (30, 150)]:
        slanted.addPoint(corner[0], corner[1])
    shapes = [("oval", OvalRoi(25, 30, 150, 120)), ("ring with a hole", ring), ("irregular polygon", PolygonRoi(slanted, Roi.POLYGON))]
    for name, shape in shapes:
        found, info = engine.analyze_objects(field, shape, nuclei_config(E))
        positive = E.count_positive_pixels(field.getProcessor().duplicate(), shape)
        t.close(sum([o["area"] for o in found]), positive, name + ": object area equals positive area inside the ROI")
        t.check(0 < len(found) < 48, name + ": only part of the field is counted")
        centres_inside = len([o for o in found if not o["touches_edge"]])
        kept, info = engine.analyze_objects(field, shape, nuclei_config(E, {"exclude_edges": True}))
        t.equal(len(kept), centres_inside, name + ": edge exclusion keeps exactly the objects not on the outline")
        t.close(sum([o["area"] for o in kept]), len(kept) * full, name + ": every object kept after edge exclusion is whole")

    box = Roi(90, 90, 80, 80)
    manager_before = RoiManager.getInstance()
    objects, info = engine.analyze_objects(mask, box, nuclei_config(E))
    t.equal(len(objects), 1, "a rectangular ROI still confines objects to the rectangle")
    t.check(RoiManager.getInstance() == manager_before, "the ROI Manager window stays connected after an analysis")


def test_morphology_is_explicit(E, t, tmp):
    def shapes():
        ip = ByteProcessor(120, 80)
        ip.setValue(255.0)
        ip.fill(Roi(10, 10, 20, 20))
        ip.fill(Roi(50, 20, 4, 30))
        ip.resetRoi()
        return ip

    def foreground(ip):
        return ip.getHistogram()[255]

    square = 20 * 20
    bar = 4 * 30
    once = shapes()
    E.binary_morphology(once, "open", 1, 1)
    t.equal(foreground(once), square + bar, "open x1 keeps a 4-pixel-wide bar")
    twice = shapes()
    E.binary_morphology(twice, "open", 2, 1)
    t.equal(foreground(twice), square, "open x2 removes the 4-pixel-wide bar and restores the square")
    closed = shapes()
    E.binary_morphology(closed, "close", 2, 1)
    t.check(foreground(closed) >= square + bar, "close never removes foreground")
    eroded = shapes()
    E.binary_morphology(eroded, "erode", 1, 1)
    t.equal(foreground(eroded), 18 * 18 + 2 * 28, "erode x1 removes one pixel from every side")
    gentle = shapes()
    E.binary_morphology(gentle, "erode", 1, 8)
    t.equal(foreground(gentle), square + bar, "neighbor count 8 erodes only fully isolated pixels")

    old_pad = Prefs.padEdges
    try:
        results = []
        for pad in [False, True]:
            Prefs.padEdges = pad
            ip = ByteProcessor(40, 40)
            ip.setValue(255.0)
            ip.fill(Roi(0, 10, 12, 12))
            ip.resetRoi()
            E.binary_morphology(ip, "erode", 1, 1)
            results.append(foreground(ip))
            t.equal(Prefs.padEdges, pad, "Fiji's pad-edges option is restored")
        t.equal(results[0], results[1], "erosion does not depend on Fiji's pad-edges option")
    finally:
        Prefs.padEdges = old_pad

    plane = blank_plane(100)
    draw_rect(plane, (10, 10, 20, 20), 2000)
    draw_rect(plane, (50, 20, 4, 30), 2000)
    data = E.AnalysisEngine().build_masks(ImagePlus("shapes", plane), nuclei_config(E, {"open_iterations": 2}), "shapes")
    t.equal(foreground(data["raw"].getProcessor()), square, "open iterations from the configuration reach the mask")
    E.safe_close_image(data["raw"])
    E.safe_close_image(data["final"])


def test_fill_holes_and_watershed(E, t, tmp):
    plane = blank_plane(100)
    draw_disk(plane, 60, 60, 20, 2000)
    draw_disk(plane, 60, 60, 5, 100)
    draw_disk(plane, 150, 100, 14, 2000)
    draw_disk(plane, 172, 100, 14, 2000)
    imp = ImagePlus("holes", plane)
    engine = E.AnalysisEngine()
    whole = Roi(0, 0, WIDTH, HEIGHT)
    old_black = Prefs.blackBackground
    try:
        for black in [False, True]:
            Prefs.blackBackground = black
            plain = engine.build_masks(imp, nuclei_config(E), "holes")
            filled = engine.build_masks(imp, nuclei_config(E, {"fill_holes": True}), "holes")
            split = engine.build_masks(imp, nuclei_config(E, {"fill_holes": True, "watershed": True}), "holes")
            hole = plain["raw"].getProcessor().getHistogram()[255]
            t.check(filled["raw"].getProcessor().getHistogram()[255] > hole, "fill holes fills the hole (black background option " + str(black) + ")")
            before, info = engine.analyze_objects(filled["final"], whole, nuclei_config(E))
            after, info = engine.analyze_objects(split["final"], whole, nuclei_config(E))
            t.equal(len(before), 2, "two touching disks are one object before watershed")
            t.equal(len(after), 3, "watershed separates the touching disks (black background option " + str(black) + ")")
            t.equal(Prefs.blackBackground, black, "Fiji's black-background option is restored")
            for data in [plain, filled, split]:
                E.safe_close_image(data["raw"])
                E.safe_close_image(data["final"])
    finally:
        Prefs.blackBackground = old_black


def test_auto_threshold_ignores_display_range(E, t, tmp):
    engine = E.AnalysisEngine()
    config = marker_config(E, {"threshold_mode": "Auto", "auto_method": "Otsu"})
    lows = []
    for display in [(0, 65535), (120, 400), (1000, 1400)]:
        imp = ImagePlus("auto", marker_plane())
        imp.setDisplayRange(display[0], display[1])
        data = engine.build_masks(imp, config, "auto")
        lows.append(data["threshold_low"])
        E.safe_close_image(data["raw"])
        E.safe_close_image(data["final"])
    t.check(lows[0] == lows[1] == lows[2], "automatic threshold is the same for every display range: " + repr(lows))
    t.check(100 < lows[0] <= 1500, "automatic threshold separates background from signal")


def test_local_threshold_ignores_conversion_option(E, t, tmp):
    if not local_threshold_available(E):
        raise SkipTest("Auto Local Threshold plugin is not installed")
    engine = E.AnalysisEngine()
    config = nuclei_config(E, {"threshold_mode": "Local Phansalkar", "local_radius": 15})
    old_scaling = ImageConverter.getDoScaling()
    try:
        counts = []
        for scaling in [True, False]:
            ImageConverter.setDoScaling(scaling)
            data = engine.build_masks(ImagePlus("local", nuclei_plane(DISKS)), config, "local")
            counts.append(data["raw"].getProcessor().getHistogram()[255])
            t.close(data["scale_max"], 2000.0, "the 8-bit rescaling range is recorded")
            E.safe_close_image(data["raw"])
            E.safe_close_image(data["final"])
            t.equal(ImageConverter.getDoScaling(), scaling, "Fiji's conversion option is left alone")
        t.equal(counts[0], counts[1], "local threshold does not depend on 'Scale when converting'")
        t.check(0 < counts[0] < WIDTH * HEIGHT / 2, "local threshold finds the disks, not the whole image")
    finally:
        ImageConverter.setDoScaling(old_scaling)


def test_clahe_without_image_window(E, t, tmp):
    if not clahe_available():
        raise SkipTest("CLAHE plugin is not installed")
    # CLAHE leaves perfectly flat regions alone, so test it on a ramp with disks.
    ramp = ByteProcessor(WIDTH, HEIGHT)
    for x in range(WIDTH):
        ramp.setValue(float(20 + x / 2))
        ramp.fill(Roi(x, 0, 1, HEIGHT))
    for disk in DISKS:
        ramp.setValue(230.0)
        ramp.fill(OvalRoi(disk[0] - DISK_RADIUS, disk[1] - DISK_RADIUS, 2 * DISK_RADIUS + 1, 2 * DISK_RADIUS + 1))
    ramp.resetRoi()
    imp = ImagePlus("clahe", ramp)
    before = E.pixel_fingerprint(imp)
    E.apply_clahe(imp, 31, 256, 3.0)
    t.check(E.pixel_fingerprint(imp) != before, "CLAHE changes the segmentation copy with no image window open")
    t.check(not imp.isLocked(), "CLAHE leaves the image unlocked")
    if not local_threshold_available(E):
        return
    engine = E.AnalysisEngine()
    prepared = ImagePlus("clahe_prepared", nuclei_plane(DISKS))
    reference = E.pixel_fingerprint(prepared)
    config = nuclei_config(E, {"threshold_mode": "Local Phansalkar", "local_radius": 15, "clahe": True, "clahe_block_size": 63})
    data = engine.build_masks(prepared, config, "clahe")
    t.check(E.is_two_level_mask(data["raw"].getProcessor()), "CLAHE + local threshold yields a mask")
    t.equal(E.pixel_fingerprint(prepared), reference, "CLAHE never touches the prepared measurement image")
    E.safe_close_image(data["raw"])
    E.safe_close_image(data["final"])


def test_all_frames(E, t, tmp):
    imp = make_image("frames.tif", slices=1, frames=2, disks_per_frame=[DISKS, DISKS[:7]])
    result = run_image(E, imp, [nuclei_config(E), marker_config(E)], [E.whole_image_roi(imp)], tmp, frames=[1, 2])
    t.equal([r["Frame"] for r in result["summary_rows"]], [1, 2], "one summary row per frame")
    t.equal([r["C1_DAPI_Object_Count"] for r in result["summary_rows"]], [len(DISKS), 7], "each frame is counted on its own")
    t.equal(len(result["threshold_rows"]), 4, "thresholds are recorded per channel and frame")
    imp.close()


def test_batch(E, t, tmp):
    folder = os.path.join(tmp, "batch_in")
    os.makedirs(os.path.join(folder, "sub"))
    for name in ["a.tif", "b.tif", os.path.join("sub", "c.tif")]:
        imp = make_image(os.path.basename(name), slices=1)
        FileSaver(imp).saveAsTiff(os.path.join(folder, name))
        imp.close()
    single = make_image("bad.tif", slices=1, channels=1)
    FileSaver(single).saveAsTiff(os.path.join(folder, "bad.tif"))
    single.close()
    for name in ["notes.txt", "corrupt.tif"]:
        handle = open(os.path.join(folder, name), "w")
        handle.write("not an image")
        handle.close()
    five = Roi(10, 10, 205, 45)
    five.setName("top_row")
    E.save_roi_array([five], os.path.join(folder, "a_ROIs.zip"))

    configs = [nuclei_config(E), marker_config(E)]
    options = E.default_global_options()
    options["save_prepared"] = False
    options["save_overlays"] = False
    base = {"extensions": "tif", "recursive": True}

    t.equal(len(E.discover_batch_images(folder, ["tif"], False)), 4, "non-recursive scan finds the top-level images")
    t.equal(len(E.discover_batch_images(folder, ["tif"], True)), 5, "recursive scan finds the subfolder image")
    t.equal(len(E.open_image_series(os.path.join(folder, "corrupt.tif"))), 0, "an unreadable file opens as nothing, without a dialog")
    t.check(E.find_roi_file(os.path.join(folder, "a.tif")) is not None, "the ROI file beside an image is found")
    t.check(E.find_roi_file(os.path.join(folder, "b.tif")) is None, "no ROI file is reported as missing")

    batch_options = dict(base)
    batch_options["roi_policy"] = E.BATCH_ROI_FILE_SKIP
    result = E.run_batch(folder, folder, configs, configs, options, batch_options)
    t.equal((result["analyzed"], result["skipped"], result["failed"]), (1, 2, 2), "ROI-file policy: one analyzed, two skipped, two failed")
    report = read_utf8(os.path.join(result["batch_dir"], "batch_report.csv"))
    t.check('"OK"' in report and '"SKIPPED"' in report and '"FAILED"' in report, "the report lists every status")
    summary = read_utf8(os.path.join(result["batch_dir"], "batch_roi_summary.csv")).strip().split("\n")
    t.equal(len(summary), 2, "combined summary holds the analyzed image only")
    t.check('"top_row"' in summary[1], "the saved ROI is used and named")
    objects = read_utf8(os.path.join(result["batch_dir"], "batch_objects.csv")).strip().split("\n")
    t.equal(len(objects) - 1, 5, "five disks lie in the saved ROI")
    for name in ["batch_thresholds.csv", "batch_log.txt", "settings.properties"]:
        t.check(os.path.isfile(os.path.join(result["batch_dir"], name)), name + " is written")

    # Results were written into the input folder; a second scan must ignore them.
    t.equal(len(E.discover_batch_images(folder, ["tif"], True)), 5, "saved masks in result folders are not picked up as input")

    batch_options["roi_policy"] = E.BATCH_ROI_FILE_WHOLE
    result = E.run_batch(folder, folder, configs, configs, options, batch_options)
    t.equal((result["analyzed"], result["skipped"], result["failed"]), (3, 0, 2), "whole-image fallback: three analyzed, two failed")
    objects = read_utf8(os.path.join(result["batch_dir"], "batch_objects.csv")).strip().split("\n")
    t.equal(len(objects) - 1, 5 + 2 * len(DISKS), "objects of all images are combined")

    batch_options["roi_policy"] = E.BATCH_ROI_MANAGER
    result = E.run_batch(folder, folder, configs, configs, options, batch_options, manager_rois=[five])
    t.equal(result["analyzed"], 3, "one ROI set can be applied to every image")

    calls = {"n": 0}

    def cancel_after_first():
        calls["n"] += 1
        return calls["n"] > 3

    batch_options["roi_policy"] = E.BATCH_ROI_WHOLE
    result = E.run_batch(folder, folder, configs, configs, options, batch_options, is_cancelled=cancel_after_first)
    t.check(result["cancelled"], "a batch can be cancelled")
    t.check(result["analyzed"] < 3, "a cancelled batch stops early")
    t.check(os.path.isfile(os.path.join(result["batch_dir"], "batch_report.csv")), "a cancelled batch still writes its report")

    rejected = False
    try:
        E.run_batch(folder, folder, [marker_config(E, {"threshold_mode": "Display minimum"})], configs, options, batch_options)
    except ValueError:
        rejected = True
    t.check(rejected, "a batch refuses display-minimum thresholds")


def test_point_validation(E, t, tmp):
    imp = make_image("validate_points.tif", slices=1)
    engine = E.AnalysisEngine()
    config = nuclei_config(E)
    prepared = engine.prepare_single_channel(imp, config, 1)
    points = [DISKS[i] for i in range(12)]
    points.append(DISKS[12])
    points.append((DISKS[12][0] + 2, DISKS[12][1] + 1))
    points.append((5, 5))

    marks = PointRoi(points[0][0], points[0][1])
    for point in points[1:]:
        marks.addPoint(point[0], point[1])
    recovered = E.annotation_points_from_rois([marks, Roi(0, 0, 5, 5), None])
    t.equal(sorted(recovered), sorted(points), "point marks are read back at their pixel positions")

    run_dir = E.unique_dir(tmp, "validation_points")
    result = E.validate_objects_against_points(engine, prepared, config, [E.whole_image_roi(imp)], recovered, run_dir, "validate_points.tif", 1)
    row = result["rows"][0]
    t.equal(row["Manual_Count"], 15, "manual count")
    t.equal(row["Detected_Count"], 15, "detected count")
    t.equal(row["True_Positive_Objects"], 13, "objects holding a mark")
    t.equal(row["False_Positive_Objects"], 2, "objects without a mark")
    t.equal(row["Missed_Annotations"], 1, "mark outside every object")
    t.equal(row["Merged_Annotations"], 1, "second mark inside one object")
    t.close(row["Precision"], 13.0 / 15.0, "precision")
    t.close(row["Recall"], 13.0 / 15.0, "recall")
    t.close(row["F1"], 13.0 / 15.0, "F1")
    for name in ["validation_summary.csv", "validation_objects.csv", "validation_points.csv", "thresholds.csv"]:
        t.check(os.path.isfile(os.path.join(run_dir, name)), name + " is written")
    t.check(len([n for n in os.listdir(run_dir) if n.endswith("_Validation_Overlay.tif")]) == 1, "validation overlay is written")

    left = Roi(0, 0, 120, HEIGHT)
    left.setName("left")
    right = Roi(120, 0, 120, HEIGHT)
    right.setName("right")
    run_dir = E.unique_dir(tmp, "validation_points_two")
    result = E.validate_objects_against_points(engine, prepared, config, [left, right], recovered, run_dir, "validate_points.tif", 1)
    t.equal(len(result["rows"]), 3, "two ROIs give two rows and a pooled row")
    t.equal(result["rows"][2]["Parent_ROI_Name"], "ALL_ROIS_POOLED", "pooled row is labelled")
    t.equal(result["rows"][0]["Manual_Count"] + result["rows"][1]["Manual_Count"], result["rows"][2]["Manual_Count"], "pooled manual count is the sum")
    prepared.close()
    imp.close()


def test_mask_validation(E, t, tmp):
    imp = make_image("validate_mask.tif", slices=1)
    engine = E.AnalysisEngine()
    config = marker_config(E)
    prepared = engine.prepare_single_channel(imp, config, 1)
    whole = [E.whole_image_roi(imp)]

    same = ByteProcessor(WIDTH, HEIGHT)
    same.setValue(255.0)
    same.fill(Roi(MARKER_OUTER[0], MARKER_OUTER[1], MARKER_OUTER[2], MARKER_OUTER[3]))
    result = E.validate_mask_against_reference(engine, prepared, config, whole, ImagePlus("same", same), E.unique_dir(tmp, "mask_same"), "validate_mask.tif", 1)
    t.close(result["rows"][0]["Dice"], 1.0, "identical masks have Dice 1")
    t.close(result["rows"][0]["Area_Difference_Percent"], 0.0, "identical masks have no area difference")

    shifted = ShortProcessor(WIDTH, HEIGHT)
    shifted.setValue(7.0)
    shifted.fill(Roi(MARKER_OUTER[0] + 20, MARKER_OUTER[1], MARKER_OUTER[2], MARKER_OUTER[3]))
    result = E.validate_mask_against_reference(engine, prepared, config, whole, ImagePlus("shifted", shifted), E.unique_dir(tmp, "mask_shifted"), "validate_mask.tif", 1)
    row = result["rows"][0]
    t.equal(row["Overlap_Pixels"], 80 * 80, "overlap of a mask shifted by 20 pixels")
    t.close(row["Dice"], 0.8, "Dice of the shifted mask")
    t.close(row["Jaccard"], 6400.0 / 9600.0, "Jaccard of the shifted mask")
    t.close(row["Precision"], 0.8, "pixel precision")

    wrong_size = False
    try:
        E.validate_mask_against_reference(engine, prepared, config, whole, ImagePlus("small", ByteProcessor(10, 10)), E.unique_dir(tmp, "mask_small"), "x", 1)
    except ValueError:
        wrong_size = True
    t.check(wrong_size, "a reference mask of another size is rejected")
    prepared.close()
    imp.close()


def test_sensitivity(E, t, tmp):
    imp = make_image("sensitivity.tif", slices=1)
    engine = E.AnalysisEngine()
    whole = [E.whole_image_roi(imp)]
    marker = marker_config(E)
    prepared = engine.prepare_single_channel(imp, marker, 1)
    result = E.run_sensitivity(engine, prepared, marker, whole, E.SENS_THRESHOLD, [0.7, 1.0, 1.2], E.unique_dir(tmp, "sens_threshold"), "sensitivity.tif", 1)
    rows = result["rows"]
    outer = MARKER_OUTER[2] * MARKER_OUTER[3] * PIXEL * PIXEL
    inner = MARKER_INNER[2] * MARKER_INNER[3] * PIXEL * PIXEL
    t.equal([r["Factor"] for r in rows], [0.7, 1.0, 1.2], "one row per factor")
    t.close(rows[1]["Parameter_Value"], 800.0, "factor 1.0 is the configured threshold")
    t.close(rows[0]["Positive_Area"], outer, "a lower threshold keeps the whole rectangle")
    t.close(rows[2]["Positive_Area"], inner, "a higher threshold keeps only the bright core")
    t.close(rows[2]["Positive_Area_Change_Percent"], (inner - outer) / outer * 100.0, "change is relative to the configured value")
    t.close(rows[1]["Positive_Area_Change_Percent"], 0.0, "the baseline has zero change")
    prepared.close()

    nuclei = nuclei_config(E, {"min_size": 10})
    prepared = engine.prepare_single_channel(imp, nuclei, 1)
    result = E.run_sensitivity(engine, prepared, nuclei, whole, E.SENS_MIN_SIZE, [1.0, 5.0], E.unique_dir(tmp, "sens_size"), "sensitivity.tif", 1)
    t.equal([r["Object_Count"] for r in result["rows"]], [len(DISKS), 0], "raising the minimum size above the disk area removes every object")
    refused = False
    try:
        E.run_sensitivity(engine, prepared, nuclei, whole, E.SENS_LOCAL_RADIUS, [1.0], E.unique_dir(tmp, "sens_bad"), "x", 1)
    except ValueError:
        refused = True
    t.check(refused, "a parameter that does not apply to the channel is refused")
    prepared.close()
    imp.close()


def test_projection_and_calibration(E, t, tmp):
    stack = ImageStack(WIDTH, HEIGHT)
    for value in [100, 300, 200]:
        stack.addSlice("", blank_plane(value))
    imp = ImagePlus("z.tif", stack)
    imp.setDimensions(1, 3, 1)
    engine = E.AnalysisEngine()
    expected = {"Max Z": 300.0, "Min Z": 100.0, "Average Z": 200.0, "Sum Z": 600.0, "Median Z": 200.0}
    for mode in sorted(expected.keys()):
        config = E.normalize_channel_config(1, {"enabled": True, "analysis_mode": "Intensity only", "projection": mode})
        prepared = engine.prepare_single_channel(imp, config, 1)
        t.close(engine.measure_intensity(prepared, Roi(0, 0, 10, 10))["mean"], expected[mode], mode + " projection")
        prepared.close()
    config = E.normalize_channel_config(1, {"enabled": True, "analysis_mode": "Intensity only", "projection": "Specific Z", "z_slice": 2})
    prepared = engine.prepare_single_channel(imp, config, 1, (2.0, 2.0, "nm"))
    t.close(engine.measure_intensity(prepared, Roi(0, 0, 10, 10))["mean"], 300.0, "specific Z slice")
    t.close(engine.measure_intensity(prepared, Roi(0, 0, 10, 10))["area"], 400.0, "calibration override is applied")
    prepared.close()
    imp.close()


TESTS = [
    ("text and CSV handle non-ASCII", test_text_and_csv),
    ("settings profile round trip", test_settings_roundtrip),
    ("configuration validation", test_config_validation),
    ("full run on a known image", test_full_run),
    ("table columns stay compatible with 2.2", test_column_order_compatible),
    ("objects are confined to the parent ROI outline", test_objects_confined_to_roi),
    ("binary morphology is explicit", test_morphology_is_explicit),
    ("fill holes and watershed", test_fill_holes_and_watershed),
    ("automatic threshold ignores the display range", test_auto_threshold_ignores_display_range),
    ("local threshold ignores the conversion option", test_local_threshold_ignores_conversion_option),
    ("CLAHE works without an image window", test_clahe_without_image_window),
    ("all frames", test_all_frames),
    ("folder batch", test_batch),
    ("validation against point marks", test_point_validation),
    ("validation against a reference mask", test_mask_validation),
    ("parameter sensitivity", test_sensitivity),
    ("Z projection and calibration override", test_projection_and_calibration)
]


def run_all():
    E = load_engine()
    lines = [E.APP_NAME + " " + E.APP_VERSION + " tests on ImageJ " + E.IMAGEJ_VERSION]
    saved = (Prefs.blackBackground, Prefs.padEdges, ImageConverter.getDoScaling())
    tmp = tempfile.mkdtemp(prefix="cellquant_tests_")
    passed = 0
    failed = 0
    skipped = 0
    total_checks = 0
    try:
        for name, function in TESTS:
            recorder = Recorder(name)
            status = "PASS"
            try:
                function(E, recorder, tmp)
            except SkipTest as skip:
                status = "SKIPPED (" + str(skip) + ")"
            except:
                recorder.failures.append("unexpected error:\n" + traceback.format_exc())
            total_checks += recorder.checks
            if status.startswith("SKIPPED"):
                skipped += 1
            elif len(recorder.failures) > 0:
                status = "FAIL"
                failed += 1
            else:
                passed += 1
            lines.append(status.ljust(8) + " " + name + "  [" + str(recorder.checks) + " checks]")
            for failure in recorder.failures:
                lines.append("         - " + failure)
    finally:
        shutil.rmtree(tmp, True)
    restored = (Prefs.blackBackground, Prefs.padEdges, ImageConverter.getDoScaling()) == saved
    if not restored:
        failed += 1
        lines.append("FAIL     Fiji options were not restored after the tests")
    lines.append("")
    lines.append(
        str(passed) + " passed, " + str(failed) + " failed, " + str(skipped) + " skipped; " +
        str(total_checks) + " checks."
    )
    report = "\n".join(lines)
    print(report)
    try:
        IJ.log(report)
    except ANY_ERROR:
        pass
    return failed


CELLQUANT_TEST_FAILURES = run_all()
if CELLQUANT_TEST_FAILURES > 0:
    raise RuntimeError("Cell Quant Studio tests failed: " + str(CELLQUANT_TEST_FAILURES) + " test(s). See the report above.")

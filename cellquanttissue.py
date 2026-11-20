# Cell Quant Studio Universal - automatic tissue outline (optional module)
#
# Finds the outline of the tissue (for example an organoid section) in an
# image, so that a depth profile has a parent ROI and a folder of images does
# not need one hand-drawn outline per image.
#
# The method, in the order it is applied (every parameter is a setting that
# is saved in the profile and in each run record):
#   1. Tissue image: one channel, or the average of all analysed channels
#      after each has been rescaled to its own range. Only the channel's Z
#      projection is used, not its background subtraction or filters, so the
#      outline does not move when segmentation settings are tuned.
#   2. Gaussian smoothing, sigma given in micrometers (pixels when the image
#      has no calibration), so that nuclei and the gaps between them merge
#      into one region.
#   3. A named ImageJ automatic threshold, computed from a 256-bin histogram
#      of the full range of the smoothed image, by default on a logarithmic
#      intensity scale. Fluorescence intensities are strongly skewed: on a
#      linear scale the dim background is squeezed into the first few bins,
#      and when the tissue fills most of the field (the usual case for a
#      zoomed organoid section) Triangle, Li and Otsu then split the tissue
#      itself. On a log scale background and tissue form two separate peaks
#      whichever of them is larger, and Triangle puts the threshold at the
#      foot of the taller peak. A consequence to keep in mind: at a sharp
#      tissue edge the outline lies up to about two smoothing widths
#      (sigma) outside the true edge when the background peak is the taller
#      one, and up to about two sigma inside it when the tissue peak is.
#      Which peak is taller depends on the image (on the log scale a peak
#      with a small relative spread is tall), so the record of every outline
#      keeps the threshold and the mean intensity on both sides of it.
#   4. Optional closing of gaps, filling of holes, and an inward shift.
#   5. Connected regions above a minimum share of the image are kept (or
#      only the largest) and traced as ROIs named Tissue_001, Tissue_002 ...
#      by decreasing area.
#
# Nothing here reads Process > Binary > Options or the black-background
# preference, so the same settings give the same outline on any Fiji.
#
# This file is loaded by load_engine_modules() into the namespace of
# cellquantengine.py: the engine's names are globals here. It has no user
# interface and opens no window.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij import ImagePlus
from ij.gui import Roi
from ij.measure import Measurements
from ij.plugin import RoiScaler
from ij.plugin.filter import GaussianBlur, RankFilters
from ij.process import (
    AutoThresholder, ByteProcessor, FloatProcessor, FloodFiller,
    ImageProcessor, ImageStatistics, Blitter, LUT
)
from java.awt import Color

import os
import time
import math


TISSUE_POLICY = "Automatic tissue outline"
TISSUE_TABLE_NAME = "tissue_outline.csv"
TISSUE_COLLECTION_DIR = "Tissue_Outlines"
# Every ROI made here carries its detection record as ROI properties with
# this prefix. They are stored in ROI files, so the record travels with the
# outline and a later run can report how its parent ROI was made.
TISSUE_TAG = "CellQuantTissue_"
TISSUE_BLUR_ACCURACY = 0.002
TISSUE_QC_MAX_SIZE = 1024
# On the log scale, intensities below this fraction of the maximum are
# raised to it, so that zero has a logarithm. Such pixels are always
# background, because the threshold lies above the lowest histogram bin.
TISSUE_LOG_FLOOR = 0.001
# An image whose smoothed range is below this fraction of its largest
# absolute value is treated as uniform (no tissue can be told apart).
TISSUE_UNIFORM_TOLERANCE = 0.00001

# (property key, option key, kind, default)
TISSUE_SETTING_SPECS = [
    ("tissue.channel", "tissue_channel", "text", "0"),
    ("tissue.sigma", "tissue_sigma", "text", "5"),
    ("tissue.method", "tissue_method", "text", "Triangle"),
    ("tissue.log_scale", "tissue_log_scale", "bool", True),
    ("tissue.min_area_percent", "tissue_min_area_percent", "text", "1"),
    ("tissue.largest_only", "tissue_largest_only", "bool", False),
    ("tissue.fill_holes", "tissue_fill_holes", "bool", True),
    ("tissue.close_radius", "tissue_close_radius", "text", "0"),
    ("tissue.shrink", "tissue_shrink", "text", "0"),
    ("tissue.save_qc", "tissue_save_qc", "bool", True),
    ("tissue.save_rois", "tissue_save_rois", "bool", True)
]

# Micrometers per calibration unit, for images not calibrated in micrometers
# (spellings that is_micrometer_unit() does not know are listed here too).
TISSUE_UNIT_MICROMETERS = {
    "micrometre": 1.0, "micrometres": 1.0,
    "nm": 0.001, "nanometer": 0.001, "nanometers": 0.001, "nanometre": 0.001, "nanometres": 0.001,
    "mm": 1000.0, "millimeter": 1000.0, "millimeters": 1000.0, "millimetre": 1000.0, "millimetres": 1000.0,
    "cm": 10000.0, "centimeter": 10000.0, "centimeters": 10000.0, "centimetre": 10000.0, "centimetres": 10000.0,
    "m": 1000000.0, "meter": 1000000.0, "meters": 1000000.0, "metre": 1000000.0, "metres": 1000000.0,
    "inch": 25400.0, "inches": 25400.0, "in": 25400.0
}
# Units that mean "not calibrated".
TISSUE_PIXEL_UNITS = ["", "pixel", "pixels", "px"]

# (ROI property, table column) of the detection record, in table order.
TISSUE_RECORD = [
    ("Version", "Detection_Version"),
    ("Image", "Detected_On_Image"),
    ("Frame", "Detection_Frame"),
    ("Channels", "Tissue_Channels"),
    ("Method", "Threshold_Method"),
    ("Scale", "Threshold_Scale"),
    ("Threshold", "Threshold"),
    ("Mean_Above", "Mean_Above_Threshold"),
    ("Mean_Below", "Mean_Below_Threshold"),
    ("Sigma", "Smoothing_Sigma"),
    ("Sigma_Pixels", "Smoothing_Sigma_Pixels"),
    ("Length_Unit", "Length_Unit"),
    ("Min_Area_Percent", "Min_Area_Percent"),
    ("Largest_Only", "Largest_Only"),
    ("Fill_Holes", "Fill_Holes"),
    ("Close_Radius", "Close_Radius"),
    ("Shrink", "Shrink"),
    ("Area_Pixels", "Detected_Area_Pixels"),
    ("Filled_Hole_Pixels", "Filled_Hole_Pixels"),
    ("Centroid_X", "Centroid_X_Pixels"),
    ("Centroid_Y", "Centroid_Y_Pixels"),
    ("Border_Sides", "Border_Sides"),
    ("Border_Contact_Pixels", "Border_Contact_Pixels"),
    ("Regions_Found", "Regions_Found"),
    ("Regions_Kept", "Regions_Kept"),
    ("Rejected_Area_Percent", "Rejected_Area_Percent"),
    ("Sigma_Pixels_Y", "Smoothing_Sigma_Pixels_Y")
]

TISSUE_TABLE_HEADERS = [
    "Run_ID", "Application_Version", "Source_Image", "Source_Path",
    "Parent_ROI_Index", "Parent_ROI_Name", "Outline_Status",
    "ROI_Area_Pixels", "ROI_Area", "ROI_Area_Unit", "Image_Area_Percent",
    "Touches_Image_Border"
] + [pair[1] for pair in TISSUE_RECORD]


# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

def tissue_settings(global_options):
    """Read and check the tissue-outline settings. A missing key uses its
    default, so options saved by an older version keep working."""
    options = global_options if global_options is not None else {}
    raw = {}
    for spec in TISSUE_SETTING_SPECS:
        value = options.get(spec[1])
        raw[spec[1]] = spec[3] if value is None else value
    method = to_text(raw["tissue_method"]).strip()
    if method not in AUTO_METHODS:
        raise ValueError(
            "Tissue outline: unknown threshold method '" + method + "'. Choose one of: " + ", ".join(AUTO_METHODS) + "."
        )
    settings = {
        "channel": parse_int_value(raw["tissue_channel"], "Tissue channel", 0),
        "sigma": parse_float_value(raw["tissue_sigma"], "Tissue smoothing width", 0),
        "method": method,
        "log_scale": parse_bool(raw["tissue_log_scale"], True),
        "min_area_percent": parse_float_value(raw["tissue_min_area_percent"], "Smallest tissue region (percent of the image)", 0),
        "largest_only": parse_bool(raw["tissue_largest_only"], False),
        "fill_holes": parse_bool(raw["tissue_fill_holes"], True),
        "close_radius": parse_float_value(raw["tissue_close_radius"], "Tissue gap-closing radius", 0),
        "shrink": parse_float_value(raw["tissue_shrink"], "Tissue outline inward shift", 0),
        "save_qc": parse_bool(raw["tissue_save_qc"], True),
        "save_rois": parse_bool(raw["tissue_save_rois"], True)
    }
    if settings["min_area_percent"] > 100:
        raise ValueError("Smallest tissue region (percent of the image) must be at most 100.")
    return settings


def tissue_length_scale(imp):
    """How lengths in the settings turn into pixels for this image.

    Returns (unit, size_x, size_y). For an image calibrated in a physical
    unit the settings are in micrometers and size is micrometers per pixel.
    For an uncalibrated image the settings are in pixels and size is 1."""
    cal = imp.getCalibration()
    unit = to_text(cal.getUnit()).strip().lower()
    factor = None
    if is_micrometer_unit(unit):
        factor = 1.0
    elif unit in TISSUE_UNIT_MICROMETERS:
        factor = TISSUE_UNIT_MICROMETERS[unit]
    if factor is None or not (cal.pixelWidth > 0) or not (cal.pixelHeight > 0):
        return ("pixel", 1.0, 1.0)
    return ("micrometer", cal.pixelWidth * factor, cal.pixelHeight * factor)


def tissue_scale_warnings(imp, scale, settings):
    """Plain-language warnings about how lengths were turned into pixels."""
    warnings = []
    unit = to_text(imp.getCalibration().getUnit()).strip()
    if scale[0] == "pixel" and unit.lower() not in TISSUE_PIXEL_UNITS:
        warnings.append(
            "the image unit '" + unit + "' is not a known length unit, so lengths are taken as pixels"
        )
    uses_length = settings["close_radius"] > 0 or settings["shrink"] > 0
    if uses_length and abs(scale[1] - scale[2]) > 0.001 * max(scale[1], scale[2]):
        warnings.append(
            "pixels are not square: gap closing and the inward shift use their mean size " +
            "(the smoothing uses each direction)"
        )
    return warnings


# -----------------------------------------------------------------------------
# Step 1 and 2: the smoothed tissue image
# -----------------------------------------------------------------------------

def tissue_plane_configs(source, configs, channel):
    """The channels that make the tissue image, each as a channel
    configuration that keeps only the Z projection of that channel."""
    total = source.getNChannels()
    by_channel = {}
    for config in (configs or []):
        by_channel[config["channel"]] = config
    if channel == 0:
        numbers = sorted(by_channel.keys())
        if len(numbers) == 0:
            numbers = list(range(1, total + 1))
    else:
        numbers = [channel]
    planes = []
    for number in numbers:
        if number < 1 or number > total:
            raise ValueError(
                "Tissue outline: channel C" + str(number) + " does not exist in this image, which has " +
                str(total) + " channel(s)."
            )
        values = {"name": "Tissue", "lut": "Original"}
        base = by_channel.get(number)
        if base is not None:
            values["projection"] = base["projection"]
            values["z_slice"] = base["z_slice"]
        plane = normalize_channel_config(number, values)
        if plane["projection"] == "Specific Z" and plane["z_slice"] > source.getNSlices():
            raise ValueError(
                "Tissue outline: C" + str(number) + " requests Z" + str(plane["z_slice"]) +
                " but the image has only " + str(source.getNSlices()) + " Z slices."
            )
        planes.append(plane)
    return planes


def tissue_rescale_unit(ip):
    """Rescale a float image in place to 0-1 over its own range."""
    ip.resetMinAndMax()
    low = ip.getMin()
    high = ip.getMax()
    if high > low:
        ip.subtract(low)
        ip.multiply(1.0 / (high - low))
    else:
        ip.multiply(0.0)


def tissue_build_image(engine, source, configs, settings, frame, calibration):
    """Steps 1 and 2. Returns a dictionary with the smoothed tissue image
    ("ip", 32-bit) and how it was made."""
    planes = tissue_plane_configs(source, configs, settings["channel"])
    built = {"channels": "+".join(["C" + str(plane["channel"]) for plane in planes])}
    total = None
    for plane in planes:
        engine.check_cancelled()
        prepared = engine.prepare_single_channel(source, plane, frame, calibration)
        try:
            if "unit" not in built:
                scale = tissue_length_scale(prepared)
                built["unit"] = scale[0]
                built["pixel_x"] = scale[1]
                built["pixel_y"] = scale[2]
                built["sigma_x"] = settings["sigma"] / scale[1]
                built["sigma_y"] = settings["sigma"] / scale[2]
                built["warnings"] = tissue_scale_warnings(prepared, scale, settings)
                limit = min(prepared.getWidth(), prepared.getHeight()) / 4.0
                if max(built["sigma_x"], built["sigma_y"]) > limit:
                    raise ValueError(
                        "Tissue outline: a smoothing width of " + repr(settings["sigma"]) + " " + scale[0] +
                        " is " + ("%.1f" % max(built["sigma_x"], built["sigma_y"])) +
                        " pixels, more than a quarter of the image. Check the setting and the image calibration."
                    )
            ip = prepared.getProcessor().duplicate().convertToFloat()
        finally:
            safe_close_image(prepared)
        ip.resetRoi()
        if built["sigma_x"] > 0 and built["sigma_y"] > 0:
            GaussianBlur().blurGaussian(ip, built["sigma_x"], built["sigma_y"], TISSUE_BLUR_ACCURACY)
        if len(planes) == 1:
            total = ip
        else:
            # Several channels: each is brought to 0-1 so that a bright
            # channel does not outweigh a dim one, then they are averaged.
            tissue_rescale_unit(ip)
            if total is None:
                total = ip
            else:
                total.copyBits(ip, 0, 0, Blitter.ADD)
    if len(planes) > 1:
        total.multiply(1.0 / len(planes))
    total.resetMinAndMax()
    built["ip"] = total
    return built


# -----------------------------------------------------------------------------
# Step 3 and 4: threshold and clean-up of the mask
# -----------------------------------------------------------------------------

def tissue_histogram(stats):
    """The histogram of an ImageStatistics as a list of whole numbers.
    (Under Jython, stats.histogram is hidden by a method of the same name.)"""
    return [int(count) for count in stats.getHistogram()]


def tissue_threshold(ip, method, log_scale=True):
    """Automatic threshold of a float image from its full data range.

    The range [minimum, maximum] of the data (log_scale off) or of its
    natural logarithm (log_scale on; values below TISSUE_LOG_FLOOR times
    the maximum are first raised to that floor) is divided into 256 equal
    bins, and the named method picks a bin t from that histogram. Tissue is
    every pixel at or above the lower edge of bin t + 1, that is minimum +
    (t + 1) * bin width, turned back into an intensity on the log scale.
    (ImageJ's own Image > Adjust > Threshold maps t to minimum + (t + 1) *
    range / 255 for 16- and 32-bit images; the value used here is the true
    bin edge of the histogram the method saw.) Returns None when the image
    is uniform or the method finds no threshold."""
    ip.resetRoi()
    ip.resetThreshold()
    stats = ImageStatistics.getStatistics(ip, Measurements.MIN_MAX | Measurements.MEAN, None)
    # Smoothing a constant image leaves differences of float rounding only
    # (about 1e-7 of the value); they must not count as contrast.
    if not (stats.max - stats.min > TISSUE_UNIFORM_TOLERANCE * max(abs(stats.max), abs(stats.min))):
        return None
    scale = "linear"
    hist_stats = stats
    if log_scale and stats.max > 0:
        scale = "log"
        work = ip.duplicate()
        work.min(stats.max * TISSUE_LOG_FLOOR)
        work.log()
        hist_stats = ImageStatistics.getStatistics(work, Measurements.MIN_MAX, None)
        work = None
        if not (hist_stats.max > hist_stats.min):
            return None
    index = AutoThresholder().getThreshold(method, tissue_histogram(hist_stats))
    if index < 0 or index >= 255:
        return None
    level = hist_stats.histMin + (index + 1) * hist_stats.binSize
    if scale == "log":
        level = math.exp(level)
    ip.setThreshold(level, stats.max, ImageProcessor.NO_LUT_UPDATE)
    above = ImageStatistics.getStatistics(ip, Measurements.MEAN | Measurements.LIMIT, None)
    below_count = stats.pixelCount - above.pixelCount
    mean_below = float("nan")
    if below_count > 0:
        mean_below = (stats.mean * stats.pixelCount - above.mean * above.pixelCount) / float(below_count)
    return {
        "level": level, "scale": scale, "data_min": stats.min, "data_max": stats.max,
        "mean_above": above.mean if above.pixelCount > 0 else float("nan"),
        "mean_below": mean_below
    }


def tissue_fill_holes(ip):
    """Fill the holes of a 0/255 mask in place: background that is not
    connected (through edge-sharing pixels) to the image border becomes
    foreground. Same result as Process > Binary > Fill Holes, written out so
    that it does not depend on the black-background preference."""
    width = ip.getWidth()
    height = ip.getHeight()
    ip.resetRoi()
    ip.setValue(127.0)
    filler = FloodFiller(ip)
    for y in range(height):
        if ip.get(0, y) == 0:
            filler.fill(0, y)
        if ip.get(width - 1, y) == 0:
            filler.fill(width - 1, y)
    for x in range(width):
        if ip.get(x, 0) == 0:
            filler.fill(x, 0)
        if ip.get(x, height - 1) == 0:
            filler.fill(x, height - 1)
    # Now: outside background 127, holes 0, tissue 255.
    # XOR with 127 gives 0, 127, 128; everything above 0 is tissue.
    ip.xor(127)
    ip.threshold(0)


def tissue_clean_mask(mask, settings, built):
    """Step 4, in place. Returns the number of hole pixels that were filled.

    Closing and the inward shift use a circular neighbourhood (the kernel of
    Process > Filters > Maximum / Minimum, so the radius is rounded the way
    those commands round it: along x and y the outline moves by a whole
    number of pixels, 2 for a radius of 1.5 to 2.4 pixels, 3 for 2.5 to
    about 3.8).
    With non-square pixels the radius uses the geometric mean pixel size.
    Pixels beyond the image border count as copies of the nearest border
    pixel, so tissue that is cut by the border is not pulled away from it."""
    pixel = math.sqrt(built["pixel_x"] * built["pixel_y"])
    close_pixels = settings["close_radius"] / pixel
    shrink_pixels = settings["shrink"] / pixel
    filters = RankFilters()
    if close_pixels >= 0.5:
        filters.rank(mask, close_pixels, RankFilters.MAX)
        filters.rank(mask, close_pixels, RankFilters.MIN)
    filled = 0
    if settings["fill_holes"]:
        before = mask.getHistogram()[255]
        tissue_fill_holes(mask)
        filled = mask.getHistogram()[255] - before
    if shrink_pixels >= 0.5:
        filters.rank(mask, shrink_pixels, RankFilters.MIN)
    return filled


# -----------------------------------------------------------------------------
# Step 5: regions and their ROIs
# -----------------------------------------------------------------------------

def tissue_roi_pixels(roi):
    """Number of pixels ImageJ measures inside an ROI."""
    mask = roi.getMask()
    if mask is None:
        bounds = roi.getBounds()
        return bounds.width * bounds.height
    mask.resetRoi()
    return mask.getHistogram()[255]


def tissue_roi_centroid(roi):
    """Centroid (x, y) in pixels of the pixels inside an ROI, pixel centres
    at +0.5 as in ImageJ's Measure. Computed here rather than taken from the
    particle analysis, whose X and Y follow Analyze > Set Scale > Global."""
    bounds = roi.getBounds()
    mask = roi.getMask()
    if mask is None:
        return (bounds.x + bounds.width / 2.0, bounds.y + bounds.height / 2.0)
    work = mask.duplicate()
    work.resetRoi()
    work.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
    stats = ImageStatistics.getStatistics(work, Measurements.CENTROID | Measurements.LIMIT, None)
    return (bounds.x + stats.xCentroid, bounds.y + stats.yCentroid)


def tissue_border_contact(roi, width, height):
    """Which image borders an ROI touches and with how many pixels.
    Returns (sides, pixels); a corner pixel is counted once."""
    bounds = roi.getBounds()
    mask = roi.getMask()
    if mask is None:
        mask = ByteProcessor(bounds.width, bounds.height)
        mask.setValue(255.0)
        mask.fill()
    left = bounds.x <= 0
    top = bounds.y <= 0
    right = bounds.x + bounds.width >= width
    bottom = bounds.y + bounds.height >= height
    sides = []
    for name, touched in [("left", left), ("top", top), ("right", right), ("bottom", bottom)]:
        if touched:
            sides.append(name)

    def count(x, y, w, h):
        mask.setRoi(x, y, w, h)
        value = mask.getHistogram()[255]
        mask.resetRoi()
        return value

    pixels = 0
    if top:
        pixels += count(0, 0, bounds.width, 1)
    if bottom and not (top and bounds.height == 1):
        pixels += count(0, bounds.height - 1, bounds.width, 1)
    first = 1 if top else 0
    last = bounds.height - 1 if bottom else bounds.height
    if last > first:
        if left:
            pixels += count(0, first, 1, last - first)
        if right and not (left and bounds.width == 1):
            pixels += count(bounds.width - 1, first, 1, last - first)
    return ("+".join(sides), pixels)


def tissue_find_regions(engine, mask, settings):
    """Connected regions of the mask that reach the minimum area, largest
    first. Each is a dictionary with its ROI and pixel measurements. The
    engine's own object detection is used, on an uncalibrated copy."""
    width = mask.getWidth()
    height = mask.getHeight()
    minimum = settings["min_area_percent"] / 100.0 * width * height
    config = normalize_channel_config(1, {
        "name": "Tissue", "analysis_mode": "Objects only", "size_units": "Pixel^2",
        "min_size": minimum, "max_size": -1, "composite_rois": True
    })
    objects = engine.analyze_objects(ImagePlus("CellQuant_tissue_mask", mask), Roi(0, 0, width, height), config)[0]
    regions = []
    for obj in objects:
        contact = tissue_border_contact(obj["roi"], width, height)
        centroid = tissue_roi_centroid(obj["roi"])
        regions.append({
            "roi": obj["roi"], "pixels": tissue_roi_pixels(obj["roi"]),
            "x": centroid[0], "y": centroid[1], "sides": contact[0], "contact": contact[1]
        })
    regions.sort(key=lambda region: -region["pixels"])
    return regions


def tissue_tag(roi, key):
    """One value of the detection record of an ROI, or "" when absent."""
    try:
        value = roi.getProperty(TISSUE_TAG + key)
    except ANY_ERROR:
        value = None
    return u"" if value is None else to_text(value)


def tissue_number_text(value):
    if isinstance(value, bool):
        return bool_text(value)
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return u""
        return to_text(repr(value))
    return to_text(value)


def tissue_tag_roi(roi, record):
    for key in record:
        roi.setProperty(TISSUE_TAG + key, tissue_number_text(record[key]))


def tissue_description(info):
    """One line saying what was found, for the batch report and the log."""
    head = TISSUE_POLICY + ": "
    if info["regions_kept"] == 0:
        return head + "no tissue found (" + info["reason"] + ")"
    text = (
        head + str(info["regions_kept"]) + " region(s), " + ("%.1f" % info["area_percent"]) +
        "% of the image, " + info["method"] + " threshold " + ("%.5g" % info["threshold"]) +
        " (" + info["scale"] + " scale) on " + info["channels"]
    )
    if info["regions_found"] > info["regions_kept"]:
        text += ", largest of " + str(info["regions_found"])
    if info["touches_border"]:
        text += ", touches the image border (" + info["border_sides"].replace("+", ", ") + ")"
    for warning in info["warnings"]:
        text += "; " + warning
    return text


def detect_tissue_outline(engine, source, configs, global_options, frame=1, keep_image=False):
    """Find the tissue outline of one image.

    engine: an AnalysisEngine (None makes one). source: the image, which is
    left unchanged. configs: the enabled channel configurations (may be
    empty). global_options: the settings dictionary; missing tissue keys
    use their defaults. Returns (rois, info): rois is a list of area ROIs
    named Tissue_001 ... by decreasing area, empty when nothing was found;
    info is a dictionary of plain values describing what was done, with
    the one-line summary under "description". With keep_image, info also
    holds the smoothed tissue image under "tissue_image".
    """
    started = time.time()
    if engine is None:
        engine = AnalysisEngine()
    settings = tissue_settings(global_options)
    calibration = None
    if global_options is not None:
        calibration = calibration_override_from_options(global_options)
    width = source.getWidth()
    height = source.getHeight()
    built = tissue_build_image(engine, source, configs, settings, frame, calibration)
    info = {
        "image": to_text(source.getTitle()), "width": width, "height": height, "frame": frame,
        "channels": built["channels"], "method": settings["method"],
        "scale": "log" if settings["log_scale"] else "linear",
        "sigma": settings["sigma"], "sigma_pixels": built["sigma_x"], "sigma_pixels_y": built["sigma_y"],
        "length_unit": built["unit"],
        "threshold": float("nan"), "mean_above": float("nan"), "mean_below": float("nan"),
        "regions_found": 0, "regions_kept": 0, "area_pixels": 0, "area_percent": 0.0,
        "rejected_area_percent": 0.0, "filled_hole_pixels": 0,
        "touches_border": False, "border_sides": "", "regions": [], "warnings": list(built["warnings"]), "reason": ""
    }
    if 0 < built["sigma_x"] < 0.5:
        info["warnings"].append("smoothing is under half a pixel - check the image calibration")
    rois = []
    threshold = tissue_threshold(built["ip"], settings["method"], settings["log_scale"])
    if threshold is None:
        info["reason"] = "the image is uniform or the " + settings["method"] + " method found no threshold"
    else:
        info["scale"] = threshold["scale"]
        info["threshold"] = threshold["level"]
        info["mean_above"] = threshold["mean_above"]
        info["mean_below"] = threshold["mean_below"]
        mask = ImagePlus("CellQuant_tissue", built["ip"]).createThresholdMask()
        if mask.isInvertedLut():
            mask.invertLut()
        info["filled_hole_pixels"] = tissue_clean_mask(mask, settings, built)
        engine.check_cancelled()
        regions = tissue_find_regions(engine, mask, settings)
        info["regions_found"] = len(regions)
        if settings["largest_only"]:
            regions = regions[:1]
        rois = tissue_finish_regions(regions, mask, settings, built, info)
        if len(rois) == 0:
            info["reason"] = (
                ("%.2f" % (100.0 * mask.getHistogram()[255] / float(width * height))) +
                "% of the image is above the " + settings["method"] + " threshold, but no connected region reaches " +
                repr(settings["min_area_percent"]) + "% of the image"
            )
    built["ip"].resetThreshold()
    if keep_image:
        info["tissue_image"] = ImagePlus("Tissue image", built["ip"])
    info["seconds"] = time.time() - started
    info["description"] = tissue_description(info)
    engine.log(info["description"])
    return (rois, info)


def tissue_finish_regions(regions, mask, settings, built, info):
    """Name the kept regions, measure them, attach the detection record to
    each ROI, and complete info. Returns the ROIs."""
    width = info["width"]
    height = info["height"]
    image_pixels = float(width * height)
    kept_pixels = sum([region["pixels"] for region in regions])
    info["regions_kept"] = len(regions)
    info["area_pixels"] = kept_pixels
    info["area_percent"] = 100.0 * kept_pixels / image_pixels
    info["rejected_area_percent"] = 100.0 * (mask.getHistogram()[255] - kept_pixels) / image_pixels
    sides = []
    rois = []
    for index in range(len(regions)):
        region = regions[index]
        roi = region["roi"]
        roi.setName("Tissue_" + str(index + 1).zfill(3))
        bounds = roi.getBounds()
        for side in region["sides"].split("+"):
            if side != "" and side not in sides:
                sides.append(side)
        record = {
            "Version": APP_VERSION, "Image": info["image"], "Width": width, "Height": height,
            "Bounds": ",".join([str(v) for v in [bounds.x, bounds.y, bounds.width, bounds.height]]),
            "Frame": info["frame"], "Channels": info["channels"], "Method": info["method"],
            "Scale": info["scale"], "Threshold": info["threshold"], "Mean_Above": info["mean_above"], "Mean_Below": info["mean_below"],
            "Sigma": info["sigma"], "Sigma_Pixels": info["sigma_pixels"], "Length_Unit": info["length_unit"],
            "Min_Area_Percent": settings["min_area_percent"], "Largest_Only": settings["largest_only"],
            "Fill_Holes": settings["fill_holes"], "Close_Radius": settings["close_radius"],
            "Shrink": settings["shrink"], "Area_Pixels": region["pixels"],
            "Filled_Hole_Pixels": info["filled_hole_pixels"],
            "Centroid_X": region["x"], "Centroid_Y": region["y"],
            "Border_Sides": region["sides"], "Border_Contact_Pixels": region["contact"],
            "Regions_Found": info["regions_found"], "Regions_Kept": len(regions),
            "Rejected_Area_Percent": info["rejected_area_percent"],
            "Sigma_Pixels_Y": info["sigma_pixels_y"]
        }
        tissue_tag_roi(roi, record)
        info["regions"].append({
            "name": to_text(roi.getName()), "area_pixels": region["pixels"],
            "area": region["pixels"] * built["pixel_x"] * built["pixel_y"],
            "area_unit": "micrometer^2" if built["unit"] == "micrometer" else "pixel^2",
            "area_percent": 100.0 * region["pixels"] / image_pixels,
            "centroid_x": region["x"], "centroid_y": region["y"],
            "border_sides": region["sides"], "border_contact_pixels": region["contact"]
        })
        rois.append(roi)
    order = ["left", "top", "right", "bottom"]
    info["border_sides"] = "+".join([side for side in order if side in sides])
    info["touches_border"] = len(sides) > 0
    return rois


# -----------------------------------------------------------------------------
# Check image: the outline drawn on the tissue image
# -----------------------------------------------------------------------------

def tissue_qc_image(planes, rois, label=""):
    """Small RGB picture for checking an outline by eye.

    planes: one or more same-sized image processors (they are not changed);
    several are averaged after rescaling each to its own range. The picture
    is reduced so that its longer side is at most TISSUE_QC_MAX_SIZE pixels,
    shown from its minimum to its 99.8th percentile, with the ROIs in yellow."""
    width = planes[0].getWidth()
    height = planes[0].getHeight()
    factor = min(1.0, TISSUE_QC_MAX_SIZE / float(max(width, height)))
    small_width = max(1, int(round(width * factor)))
    small_height = max(1, int(round(height * factor)))
    total = None
    for plane in planes:
        small = plane.duplicate()
        small.resetRoi()
        if small_width != width or small_height != height:
            small = small.resize(small_width, small_height, True)
        small = small.convertToFloat()
        tissue_rescale_unit(small)
        if total is None:
            total = small
        else:
            total.copyBits(small, 0, 0, Blitter.ADD)
    stats = ImageStatistics.getStatistics(total, Measurements.MIN_MAX, None)
    high = stats.max
    if stats.max > stats.min:
        wanted = 0.998 * stats.pixelCount
        seen = 0
        counts = tissue_histogram(stats)
        for index in range(len(counts)):
            seen += counts[index]
            if seen >= wanted:
                high = stats.histMin + (index + 1) * stats.binSize
                break
    total.setMinAndMax(stats.min, high)
    gray = total.convertToByte(True)
    # Shown in gray whatever colour table the channel has, so that the
    # yellow outline stands out the same way on every image.
    gray.setColorModel(LUT.createLutFromColor(Color.white))
    picture = gray.convertToRGB()
    picture.setColor(Color.yellow)
    picture.setLineWidth(2)
    scale_x = small_width / float(width)
    scale_y = small_height / float(height)
    for roi in rois:
        picture.draw(RoiScaler.scale(clone_roi(roi), scale_x, scale_y, False))
    if label != "":
        try:
            picture.setLineWidth(1)
            picture.drawString(to_text(label), 6, 18)
        except ANY_ERROR:
            pass
    return ImagePlus("Tissue outline check", picture)


def save_tissue_qc(path, planes, rois, label=""):
    """Write the check image of tissue_qc_image() as a TIFF file."""
    picture = tissue_qc_image(planes, rois, label)
    try:
        save_tiff(picture, path)
    finally:
        safe_close_image(picture)
    return path


# -----------------------------------------------------------------------------
# Batch ROI policy
# -----------------------------------------------------------------------------

def tissue_batch_rois(request):
    """Batch ROI policy "Automatic tissue outline": the parent ROIs of each
    image are its detected tissue regions. An image without tissue is
    skipped, and the reason goes into batch_report.csv."""
    # The batch does not hand over its cancel check yet; it is used as soon
    # as the request carries one under "is_cancelled".
    engine = AnalysisEngine(log=request.get("log"), is_cancelled=request.get("is_cancelled"))
    rois, info = detect_tissue_outline(
        engine, request["imp"], request.get("configs"), request.get("global_options"), 1
    )
    if len(rois) == 0:
        return (None, info["description"])
    return (rois, info["description"])


# -----------------------------------------------------------------------------
# Record of the outline in every run that uses one
# -----------------------------------------------------------------------------

def tissue_file_stem(source_title):
    """Image name without extension, usable as a file name. Kept as it is
    when possible so that <stem>_ROIs.zip is found by the batch policy
    "ROI file per image"."""
    stem = strip_extension(source_title)
    for token in ["\\", "/", ":", "*", "?", '"', "<", ">", "|"]:
        if token in stem:
            return sanitize_name(stem)
    return stem if stem.strip() != "" else "image"


def tissue_collection_folder(run_dir):
    """The folder that collects the outlines of a whole batch, or None when
    the run is not part of a batch. A batch folder is recognised by its
    name, CellQuant_Batch_<time>, which holds one run folder per image."""
    parent = os.path.dirname(os.path.abspath(to_text(run_dir)))
    if to_text(os.path.basename(parent)).startswith("CellQuant_Batch_"):
        return os.path.join(parent, TISSUE_COLLECTION_DIR)
    return None


def tissue_outline_row(context, index, roi):
    """The tissue_outline.csv row of one parent ROI that carries a
    detection record."""
    width = context["width"]
    height = context["height"]
    pixels = tissue_roi_pixels(roi)
    bounds = roi.getBounds()
    bounds_text = ",".join([str(v) for v in [bounds.x, bounds.y, bounds.width, bounds.height]])
    unchanged = (
        tissue_tag(roi, "Image") == to_text(context["source_title"]) and
        tissue_tag(roi, "Width") == str(width) and tissue_tag(roi, "Height") == str(height) and
        tissue_tag(roi, "Bounds") == bounds_text and tissue_tag(roi, "Area_Pixels") == str(pixels)
    )
    summary = context["summary_rows"][index]
    row = {
        "Run_ID": context["run_id"], "Application_Version": APP_VERSION,
        "Source_Image": context["source_title"], "Source_Path": context["source_path"],
        "Parent_ROI_Index": index + 1, "Parent_ROI_Name": summary.get("Parent_ROI_Name", ""),
        "Outline_Status": "as detected" if unchanged else "edited or reused",
        "ROI_Area_Pixels": pixels, "ROI_Area": summary.get("ROI_Area", ""),
        "ROI_Area_Unit": summary.get("ROI_Area_Unit", ""),
        "Image_Area_Percent": 100.0 * pixels / float(width * height),
        "Touches_Image_Border": bool_text(tissue_border_contact(roi, width, height)[0] != "")
    }
    for pair in TISSUE_RECORD:
        row[pair[1]] = tissue_tag(roi, pair[0])
    return row


def tissue_save_run_files(context, tagged, settings):
    """Write the check image and the ROI set of the detected outlines into
    the run folder and, in a batch, also into one folder for all images."""
    rois = [context["parent_rois"][index] for index in tagged]
    run_dir = context["run_dir"]
    stem = tissue_file_stem(context["source_title"])
    collection = tissue_collection_folder(run_dir)
    shared = stem
    if collection is not None:
        # Two images of one batch (in different subfolders) can have the
        # same name. The second must not overwrite the first one's outline.
        number = 1
        while (os.path.exists(os.path.join(collection, shared + "_ROIs.zip")) or
               os.path.exists(os.path.join(collection, shared + "_Tissue_Outline_QC.tif"))):
            number += 1
            shared = stem + "_" + str(number)
        if shared != stem:
            context["engine"].log(
                "Another image of this batch is also named " + stem + ": its outline is collected as " +
                shared + " in " + TISSUE_COLLECTION_DIR + ". Rename one of the images before reusing these outlines."
            )
    if settings["save_rois"]:
        ensure_dir(os.path.join(run_dir, "ROIs"))
        save_roi_array(rois, os.path.join(run_dir, "ROIs", "Tissue_Outline_ROIs.zip"))
        if collection is not None:
            ensure_dir(collection)
            save_roi_array(rois, os.path.join(collection, shared + "_ROIs.zip"))
    if settings["save_qc"]:
        prepared = context["prepared"]
        used = [text for text in tissue_tag(rois[0], "Channels").replace("C", "").split("+") if text.isdigit()]
        channels = [int(text) for text in used if int(text) in prepared]
        if len(channels) == 0:
            channels = sorted(prepared.keys())
        planes = [prepared[channel].getProcessor() for channel in channels]
        label = (
            stem + "  " + tissue_tag(rois[0], "Method") + "  shown: " +
            "+".join(["C" + str(channel) for channel in channels])
        )
        picture = tissue_qc_image(planes, rois, label)
        try:
            # Inside the run folder the image name is shortened like every
            # other file name there when it had to be (long paths).
            qc_name = shortened_file_name(
                sanitize_name(stem) + "_Tissue_Outline_QC.tif",
                context.get("base", u""), context.get("file_base", context.get("base", u""))
            )
            save_tiff(picture, os.path.join(run_dir, "Overlays", qc_name))
            if collection is not None:
                save_tiff(picture, os.path.join(collection, shared + "_Tissue_Outline_QC.tif"))
        finally:
            safe_close_image(picture)


class TissueOutlineExtension(RunExtension):
    """Records how an automatically detected parent ROI was made.

    Active in every run but does nothing unless a parent ROI carries the
    detection record that detect_tissue_outline() attaches. Then it writes
    tissue_outline.csv and, when switched on, a check image and the ROI
    set. It adds no column to any existing table."""

    name = "tissue_outline"

    def is_enabled(self, global_options, configs):
        return True

    def begin_run(self, context):
        state = context["state"][self.name]
        state["rows"] = []
        state["done"] = False
        state["tagged"] = []
        rois = context["parent_rois"]
        for index in range(len(rois)):
            if tissue_tag(rois[index], "Method") != "":
                state["tagged"].append(index)

    def frame_done(self, context):
        state = context["state"][self.name]
        if state["done"] or len(state["tagged"]) == 0:
            return
        state["done"] = True
        for index in state["tagged"]:
            state["rows"].append(tissue_outline_row(context, index, context["parent_rois"][index]))
        # A problem with these optional files (or with a tissue setting, when
        # the outlines were made earlier and only reused here) is logged and
        # never stops the analysis itself.
        try:
            settings = tissue_settings(context["global_options"])
            if settings["save_qc"] or settings["save_rois"]:
                tissue_save_run_files(context, state["tagged"], settings)
        except ANY_ERROR as exc:
            context["engine"].log("Could not save the tissue outline files: " + to_text(exc))

    def end_run(self, context):
        rows = context["state"][self.name].get("rows", [])
        if len(rows) == 0:
            return []
        return [(TISSUE_TABLE_NAME, TISSUE_TABLE_HEADERS, rows)]


register_global_settings(TISSUE_SETTING_SPECS)
register_settings_form(
    "Tissue outline",
    "Finds the outline of the tissue for use as parent ROI. Lengths are in micrometers "
    "(in pixels when the image is not calibrated).",
    [
        ("tissue_channel", "Tissue channel (0 = all analysed channels, averaged)", "text", None),
        ("tissue_sigma", "Smoothing width (sigma)", "text", None),
        ("tissue_method", "Threshold method", "choice", list(AUTO_METHODS)),
        ("tissue_log_scale", "Threshold on a log intensity scale (recommended)", "bool", None),
        ("tissue_min_area_percent", "Smallest region kept (percent of the image)", "text", None),
        ("tissue_largest_only", "Keep only the largest region", "bool", None),
        ("tissue_fill_holes", "Fill holes inside the tissue", "bool", None),
        ("tissue_close_radius", "Close gaps up to this radius (0 = off)", "text", None),
        ("tissue_shrink", "Move the outline inward by (0 = off)", "text", None),
        ("tissue_save_qc", "Save a check image of each outline", "bool", None),
        ("tissue_save_rois", "Save each detected outline as a ROI file", "bool", None)
    ]
)
register_batch_roi_policy(TISSUE_POLICY, tissue_batch_rois)
register_run_extension(TissueOutlineExtension())
register_engine_module("cellquanttissue.py", "2.3.0")

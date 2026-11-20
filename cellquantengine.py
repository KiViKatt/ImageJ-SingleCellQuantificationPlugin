# Cell Quant Studio Universal - analysis engine
#
# Image preparation, independent channel segmentation, ROI-based measurement,
# object counting, folder batch processing, validation against manual
# annotation, and parameter sensitivity for Fiji/ImageJ (Jython).
#
# This file has no user interface. It never imports or touches a Swing
# component, so the same code runs from the window (cellquantcore.py), from
# the headless batch script, and from the automated tests.
#
# Installation:
#   Keep this file beside the launcher (Cell_Quant_Studio_final.py) and
#   cellquantcore.py, for example in
#   Fiji.app/plugins/Scripts/Cell Quant Studio/
#   The launcher loads this file first and then the window.
#
# Compatibility target: Fiji with ImageJ 1.53+ and Jython 2.7.
# The source is kept plain ASCII with no encoding declaration on purpose:
# Fiji's script engine rejects an encoding line in some launch paths.

from ij import IJ, WindowManager, ImagePlus, Prefs
from ij.plugin import Duplicator, ZProjector
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable, Measurements
from ij.gui import Roi, Overlay, OvalRoi
from ij.process import ImageProcessor, LUT, AutoThresholder, ByteProcessor, FloatProcessor, Blitter
from ij.io import FileSaver, RoiDecoder, RoiEncoder

from java.awt import Color, GraphicsEnvironment
from java.io import (
    FileInputStream, FileOutputStream, BufferedInputStream, BufferedOutputStream,
    ByteArrayOutputStream, DataOutputStream
)
from java.lang import Throwable as JavaThrowable
from java.util import Properties
from java.util.zip import ZipInputStream, ZipOutputStream, ZipEntry

import os
import time
import traceback
import math
import codecs


# Under Jython, "except ANY_ERROR" does not catch an exception thrown by Java
# code (ImageJ, Bio-Formats, a plugin). Error handling in this application
# uses this pair so that both Python and Java failures are reported.
ANY_ERROR = (Exception, JavaThrowable)


APP_NAME = "Cell Quant Studio Universal"
APP_VERSION = "2.3.0"
MAX_CHANNELS = 8

ANALYSIS_MODES = [
    "Ignore",
    "Intensity only",
    "Positive area",
    "Objects only",
    "Objects + area"
]

PROJECTION_MODES = [
    "Max Z",
    "Average Z",
    "Sum Z",
    "Min Z",
    "Standard deviation Z",
    "Median Z",
    "Current Z",
    "Specific Z"
]

THRESHOLD_MODES = [
    "Auto",
    "Manual",
    "Display minimum",
    "Local Phansalkar"
]

LUT_NAMES = [
    "Original", "Grays", "Red", "Green", "Blue", "Magenta", "Cyan", "Yellow"
]

SIZE_UNITS = ["Calibrated area", "Pixel^2"]
FRAME_MODES = ["Current frame only", "All frames"]
ROI_SOURCES = ["Whole image", "Current image ROI", "ROI Manager"]

PRESET_NAMES = [
    "MAP2 + DAPI (C3 marker, C2 nuclei, CLAHE)",
    "MAP2 + DAPI (C3 marker, C4 nuclei)",
    "Nuclear objects + independent marker",
    "Four-channel cell-based",
    "Positive-area quantification",
    "Intensity-only quantification",
    "Custom"
]

BATCH_ROI_WHOLE = "Whole image"
BATCH_ROI_FILE_SKIP = "ROI file per image (skip image if missing)"
BATCH_ROI_FILE_WHOLE = "ROI file per image (whole image if missing)"
BATCH_ROI_MANAGER = "Current ROI Manager set for every image"
BATCH_ROI_POLICIES = [
    BATCH_ROI_WHOLE, BATCH_ROI_FILE_SKIP, BATCH_ROI_FILE_WHOLE, BATCH_ROI_MANAGER
]
DEFAULT_BATCH_EXTENSIONS = "tif, tiff, czi, lif, nd2, lsm, oib, oir, ics, dv, vsi"
NATIVE_IMAGE_EXTENSIONS = ["tif", "tiff", "png", "jpg", "jpeg", "gif", "bmp"]
ROI_FILE_SUFFIXES = ["_ROIs.zip", "_RoiSet.zip", "_rois.zip", ".zip", "_ROI.roi", ".roi"]

ANNOTATION_POINTS_IMAGE = "Multi-point selection on the source image"
ANNOTATION_POINTS_FILE = "Point ROI file (.roi or .zip)"
ANNOTATION_MASK_FILE = "Reference binary mask image (file)"
ANNOTATION_SOURCES = [ANNOTATION_POINTS_IMAGE, ANNOTATION_POINTS_FILE, ANNOTATION_MASK_FILE]

SENS_THRESHOLD = "Threshold level"
SENS_LOCAL_RADIUS = "Local threshold radius"
SENS_LOCAL_K = "Phansalkar parameter 1"
SENS_CLAHE_SLOPE = "CLAHE maximum slope"
SENS_MIN_SIZE = "Minimum object area"
SENSITIVITY_PARAMETERS = [
    SENS_THRESHOLD, SENS_LOCAL_RADIUS, SENS_LOCAL_K, SENS_CLAHE_SLOPE, SENS_MIN_SIZE
]
DEFAULT_SENSITIVITY_FACTORS = "0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3"

try:
    AUTO_METHODS = list(AutoThresholder.getMethods())
except ANY_ERROR:
    AUTO_METHODS = [
        "Default", "Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li",
        "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu",
        "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"
    ]

try:
    IMAGEJ_VERSION = str(IJ.getFullVersion())
except ANY_ERROR:
    IMAGEJ_VERSION = "unknown"

# Python 2 (Jython) and Python 3 name differences, kept in one place.
try:
    text_type = unicode
except NameError:
    text_type = str
try:
    number_types = (int, long, float)
except NameError:
    number_types = (int, float)


# -----------------------------------------------------------------------------
# General utilities
# -----------------------------------------------------------------------------

def to_text(value):
    """Return value as text. Never raises on non-ASCII characters such as the
    micro sign in a calibration unit or an accented file name."""
    if value is None:
        return u""
    if isinstance(value, text_type):
        return value
    if isinstance(value, bytes):
        try:
            return value.decode("utf-8")
        except ANY_ERROR:
            return value.decode("latin-1", "replace")
    try:
        return text_type(value)
    except ANY_ERROR:
        try:
            return text_type(repr(value))
        except ANY_ERROR:
            return u"?"


def bool_text(value):
    return "true" if bool(value) else "false"


def parse_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return to_text(value).strip().lower() in ["1", "true", "yes", "on"]


def is_number(value):
    return isinstance(value, number_types) and not isinstance(value, bool)


def parse_float_value(raw, label, minimum=None, allow_negative_one=False):
    raw_text = to_text(raw).strip()
    try:
        value = float(raw_text)
    except ANY_ERROR:
        raise ValueError(label + " must be numeric. Current value: " + raw_text)
    if math.isnan(value) or math.isinf(value):
        raise ValueError(label + " must be a finite number. Current value: " + raw_text)
    if allow_negative_one and value == -1:
        return value
    if minimum is not None and value < minimum:
        raise ValueError(label + " must be at least " + str(minimum) + ".")
    return value


def parse_int_value(raw, label, minimum=None, maximum=None):
    raw_text = to_text(raw).strip()
    try:
        value = int(float(raw_text))
    except ANY_ERROR:
        raise ValueError(label + " must be an integer. Current value: " + raw_text)
    if minimum is not None and value < minimum:
        raise ValueError(label + " must be at least " + str(minimum) + ".")
    if maximum is not None and value > maximum:
        raise ValueError(label + " must be at most " + str(maximum) + ".")
    return value


def parse_factor_list(text):
    """Parse '0.8, 0.9, 1.0' into a list of positive floats. 1.0 is always
    included because it is the baseline every other value is compared with."""
    factors = []
    for part in to_text(text).replace(";", ",").split(","):
        part = part.strip()
        if part == "":
            continue
        value = parse_float_value(part, "Sensitivity factor", 0.000001)
        if value not in factors:
            factors.append(value)
    if 1.0 not in factors:
        factors.append(1.0)
    factors.sort()
    return factors


def sanitize_name(value):
    value = to_text(value).strip()
    replacements = ["\\", "/", ":", "*", "?", '"', "<", ">", "|"]
    for token in replacements:
        value = value.replace(token, "_")
    value = value.replace(" ", "_")
    while "__" in value:
        value = value.replace("__", "_")
    value = value.strip("_.")
    return value if value else "unnamed"


def strip_extension(value):
    name = os.path.basename(to_text(value))
    dot = name.rfind(".")
    if dot > 0:
        name = name[:dot]
    return name


def timestamp_text():
    return time.strftime("%Y%m%d_%H%M%S")


def csv_escape(value):
    if value is None:
        text = u""
    elif isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            text = u""
        else:
            text = to_text(repr(value))
    else:
        text = to_text(value)
    return u'"' + text.replace(u'"', u'""') + u'"'


def csv_line(headers, row):
    return u",".join([csv_escape(row.get(h, "")) for h in headers]) + u"\n"


def write_csv(path, headers, rows):
    """Write a UTF-8 CSV file (no byte-order mark)."""
    handle = codecs.open(path, "w", "utf-8")
    try:
        handle.write(u",".join([csv_escape(h) for h in headers]) + u"\n")
        for row in rows:
            handle.write(csv_line(headers, row))
    finally:
        handle.close()


class CsvAppender(object):
    """Writes a CSV file row by row so that a long batch keeps its results on
    disk as it goes instead of holding every row in memory."""

    def __init__(self, path, headers):
        self.path = path
        self.headers = list(headers)
        self.count = 0
        self.handle = codecs.open(path, "w", "utf-8")
        self.handle.write(u",".join([csv_escape(h) for h in self.headers]) + u"\n")
        self.handle.flush()

    def append(self, rows):
        for row in rows:
            self.handle.write(csv_line(self.headers, row))
            self.count += 1
        self.handle.flush()

    def close(self):
        try:
            self.handle.close()
        except ANY_ERROR:
            pass


def write_text_file(path, lines):
    handle = codecs.open(path, "w", "utf-8")
    try:
        for line in lines:
            handle.write(to_text(line) + u"\n")
    finally:
        handle.close()


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def unique_dir(root, name):
    ensure_dir(root)
    initial = os.path.join(root, name)
    candidate = initial
    counter = 2
    while os.path.exists(candidate):
        candidate = initial + "_" + str(counter)
        counter += 1
    os.makedirs(candidate)
    return candidate


def unique_run_dir(root, base):
    return unique_dir(root, sanitize_name(base) + "_CellQuant_" + timestamp_text())


def safe_close_image(imp):
    if imp is None:
        return
    try:
        imp.changes = False
    except ANY_ERROR:
        pass
    try:
        imp.close()
    except ANY_ERROR:
        try:
            win = imp.getWindow()
            if win is not None:
                win.close()
        except ANY_ERROR:
            pass


def is_headless():
    try:
        return bool(GraphicsEnvironment.isHeadless())
    except ANY_ERROR:
        return False


def clone_roi(roi):
    if roi is None:
        return None
    try:
        return roi.clone()
    except ANY_ERROR:
        return roi


def roi_name(roi, fallback):
    try:
        name = roi.getName()
        if name is not None and to_text(name).strip() != "":
            return to_text(name)
    except ANY_ERROR:
        pass
    return fallback


def rt_value(rt, row, names, default_value=float("nan")):
    for name in names:
        try:
            value = rt.getValue(name, row)
            if not math.isnan(value):
                return value
        except ANY_ERROR:
            pass
    return default_value


def is_micrometer_unit(unit):
    normalized = to_text(unit).strip().lower()
    return normalized in [
        u"um", u"micron", u"microns", u"micrometer", u"micrometers",
        u"\u00b5m", u"\u03bcm"
    ]


def safe_ratio(numerator, denominator):
    if denominator is None or denominator == 0:
        return float("nan")
    return float(numerator) / float(denominator)


# -----------------------------------------------------------------------------
# Settings model: one definition shared by the window, profiles, the headless
# batch script, and the tests
# -----------------------------------------------------------------------------

CHANNEL_DEFAULTS = {
    "enabled": False,
    "name": "",
    "analysis_mode": "Ignore",
    "projection": "Max Z",
    "z_slice": 1,
    "subtract_background": False,
    "background_radius": 30,
    "light_background": False,
    "median_radius": 0,
    "gaussian_sigma": 0,
    "lut": "Grays",
    "threshold_mode": "Auto",
    "auto_method": "Otsu",
    "dark_background": True,
    "auto_scale": 1.0,
    "manual_low": 500,
    "manual_high": -1,
    "local_radius": 35,
    "local_parameter_1": 0.25,
    "local_parameter_2": 0.5,
    "clahe": False,
    "clahe_block_size": 127,
    "clahe_histogram_bins": 256,
    "clahe_max_slope": 3.0,
    "open_iterations": 0,
    "close_iterations": 0,
    "erode_iterations": 0,
    "dilate_iterations": 0,
    "morph_count": 1,
    "fill_holes": False,
    "watershed": False,
    "size_units": "Calibrated area",
    "min_size": 0,
    "max_size": -1,
    "min_circularity": 0.0,
    "max_circularity": 1.0,
    "exclude_edges": False,
    "include_holes": False,
    "composite_rois": False,
    "four_connected": False,
    "cross_measure": False
}

CHANNEL_BOOLEAN_KEYS = [
    "enabled", "subtract_background", "light_background", "dark_background", "clahe",
    "fill_holes", "watershed", "exclude_edges", "include_holes",
    "composite_rois", "four_connected", "cross_measure"
]
CHANNEL_STRING_KEYS = [
    "name", "analysis_mode", "projection", "lut", "threshold_mode",
    "auto_method", "size_units"
]
CHANNEL_NUMERIC_KEYS = [
    "z_slice", "background_radius", "median_radius", "gaussian_sigma",
    "auto_scale", "manual_low", "manual_high", "local_radius",
    "local_parameter_1", "local_parameter_2", "clahe_block_size",
    "clahe_histogram_bins", "clahe_max_slope", "open_iterations",
    "close_iterations", "erode_iterations", "dilate_iterations", "morph_count",
    "min_size", "max_size", "min_circularity", "max_circularity"
]

# (property key, option key, kind, default). Property keys written by 2.2.x
# keep their names, so older profiles load unchanged.
GLOBAL_SETTING_SPECS = [
    ("global.frame_mode", "frame_mode", "text", "Current frame only"),
    ("global.roi_source", "roi_source", "text", "ROI Manager"),
    ("global.preset", "preset", "text", "Custom"),
    ("global.override_calibration", "override_calibration", "bool", False),
    ("global.pixel_width", "pixel_width", "text", "1.0"),
    ("global.pixel_height", "pixel_height", "text", "1.0"),
    ("global.spatial_unit", "spatial_unit", "text", "pixel"),
    ("output.save_prepared", "save_prepared", "bool", True),
    ("output.save_masks", "save_masks", "bool", True),
    ("output.save_parent_rois", "save_parent_rois", "bool", True),
    ("output.save_object_rois", "save_object_rois", "bool", True),
    ("output.save_overlays", "save_overlays", "bool", True),
    ("output.show_summary", "show_summary", "bool", True),
    ("lifecycle.keep_prepared_open", "keep_prepared_open", "bool", True),
    ("lifecycle.keep_masks_open", "keep_masks_open", "bool", True),
    ("lifecycle.close_source", "close_source", "bool", False),
    ("lifecycle.close_generated_after", "close_generated_after", "bool", False),
    ("lifecycle.tile_prepared", "tile_prepared", "bool", False),
    ("batch.extensions", "batch_extensions", "text", DEFAULT_BATCH_EXTENSIONS),
    ("batch.recursive", "batch_recursive", "bool", False),
    ("batch.roi_policy", "batch_roi_policy", "text", BATCH_ROI_WHOLE),
    ("batch.force_bioformats", "batch_force_bioformats", "bool", False),
    ("batch.all_series", "batch_all_series", "bool", False)
]


def default_global_options():
    options = {}
    for spec in GLOBAL_SETTING_SPECS:
        options[spec[1]] = spec[3]
    return options


def normalize_channel_config(channel_number, values):
    """Turn raw values (text from the window or a profile, or typed values)
    into one validated channel configuration. Missing keys use the defaults."""
    raw = dict(CHANNEL_DEFAULTS)
    raw["name"] = "Channel_" + str(channel_number)
    if values is not None:
        for key in values:
            raw[key] = values[key]
    label = "C" + str(channel_number)
    analysis_mode = to_text(raw["analysis_mode"])
    name = to_text(raw["name"]).strip()
    config = {
        "channel": channel_number,
        "enabled": parse_bool(raw["enabled"]) and analysis_mode != "Ignore",
        "name": name if name else "Channel_" + str(channel_number),
        "analysis_mode": analysis_mode,
        "projection": to_text(raw["projection"]),
        "z_slice": parse_int_value(raw["z_slice"], label + " Z slice", 1),
        "subtract_background": parse_bool(raw["subtract_background"]),
        "background_radius": parse_float_value(raw["background_radius"], label + " background radius", 0),
        "light_background": parse_bool(raw["light_background"]),
        "median_radius": parse_float_value(raw["median_radius"], label + " median radius", 0),
        "gaussian_sigma": parse_float_value(raw["gaussian_sigma"], label + " Gaussian sigma", 0),
        "lut": to_text(raw["lut"]),
        "threshold_mode": to_text(raw["threshold_mode"]),
        "auto_method": to_text(raw["auto_method"]),
        "dark_background": parse_bool(raw["dark_background"]),
        "auto_scale": parse_float_value(raw["auto_scale"], label + " automatic threshold multiplier", 0),
        "manual_low": parse_float_value(raw["manual_low"], label + " manual lower threshold"),
        "manual_high": parse_float_value(raw["manual_high"], label + " manual upper threshold", allow_negative_one=True),
        "local_radius": parse_int_value(raw["local_radius"], label + " local threshold radius", 1),
        "local_parameter_1": parse_float_value(raw["local_parameter_1"], label + " Phansalkar parameter 1"),
        "local_parameter_2": parse_float_value(raw["local_parameter_2"], label + " Phansalkar parameter 2"),
        "clahe": parse_bool(raw["clahe"]),
        "clahe_block_size": parse_int_value(raw["clahe_block_size"], label + " CLAHE block size", 3),
        "clahe_histogram_bins": parse_int_value(raw["clahe_histogram_bins"], label + " CLAHE histogram bins", 2),
        "clahe_max_slope": parse_float_value(raw["clahe_max_slope"], label + " CLAHE maximum slope", 0.000001),
        "open_iterations": parse_int_value(raw["open_iterations"], label + " open iterations", 0),
        "close_iterations": parse_int_value(raw["close_iterations"], label + " close iterations", 0),
        "erode_iterations": parse_int_value(raw["erode_iterations"], label + " erode iterations", 0),
        "dilate_iterations": parse_int_value(raw["dilate_iterations"], label + " dilate iterations", 0),
        "morph_count": parse_int_value(raw["morph_count"], label + " binary neighbor count", 1, 8),
        "fill_holes": parse_bool(raw["fill_holes"]),
        "watershed": parse_bool(raw["watershed"]),
        "size_units": to_text(raw["size_units"]),
        "min_size": parse_float_value(raw["min_size"], label + " minimum object area", 0),
        "max_size": parse_float_value(raw["max_size"], label + " maximum object area", allow_negative_one=True),
        "min_circularity": parse_float_value(raw["min_circularity"], label + " minimum circularity", 0),
        "max_circularity": parse_float_value(raw["max_circularity"], label + " maximum circularity", 0),
        "exclude_edges": parse_bool(raw["exclude_edges"]),
        "include_holes": parse_bool(raw["include_holes"]),
        "composite_rois": parse_bool(raw["composite_rois"]),
        "four_connected": parse_bool(raw["four_connected"]),
        "cross_measure": parse_bool(raw["cross_measure"])
    }
    if config["analysis_mode"] not in ANALYSIS_MODES:
        raise ValueError(label + " has an unknown analysis mode: " + config["analysis_mode"])
    if config["projection"] not in PROJECTION_MODES:
        raise ValueError(label + " has an unknown Z processing mode: " + config["projection"])
    if config["threshold_mode"] not in THRESHOLD_MODES:
        raise ValueError(label + " has an unknown threshold mode: " + config["threshold_mode"])
    if config["size_units"] not in SIZE_UNITS:
        raise ValueError(label + " has unknown particle-size units: " + config["size_units"])
    return config


def needs_mask(config):
    return config["analysis_mode"] in ["Positive area", "Objects only", "Objects + area"]


def needs_objects(config):
    return config["analysis_mode"] in ["Objects only", "Objects + area"]


def needs_positive_area(config):
    return config["analysis_mode"] in ["Positive area", "Objects + area"]


def channel_column_prefix(config):
    return "C" + str(config["channel"]) + "_" + sanitize_name(config["name"])


def validate_channel_configs(source, configs):
    """Check a set of enabled channel configurations against one image."""
    if source is None:
        raise ValueError("No source image is selected.")
    if len(configs) == 0:
        raise ValueError("No channels are enabled for analysis.")
    if source.getBitDepth() == 24:
        for config in configs:
            if config["analysis_mode"] != "Intensity only":
                raise ValueError("RGB images are not supported for segmentation. Open the channels as grayscale using Bio-Formats.")

    source_channels = source.getNChannels()
    source_slices = source.getNSlices()
    used_names = set()
    for config in configs:
        ch = config["channel"]
        if ch < 1 or ch > source_channels:
            raise ValueError("Enabled channel C" + str(ch) + " does not exist in the source image.")
        if config["projection"] == "Specific Z" and config["z_slice"] > source_slices:
            raise ValueError("C" + str(ch) + " requests Z" + str(config["z_slice"]) + " but the source has only " + str(source_slices) + " Z slices.")
        clean_name = sanitize_name(config["name"]).lower()
        if clean_name in used_names:
            raise ValueError("Channel names must be unique. Duplicate name: " + config["name"])
        used_names.add(clean_name)

        if config["manual_high"] != -1 and config["manual_high"] < config["manual_low"]:
            raise ValueError("C" + str(ch) + " manual upper threshold is below its lower threshold.")
        if config["max_size"] != -1 and config["max_size"] < config["min_size"]:
            raise ValueError("C" + str(ch) + " maximum object area is below its minimum object area.")
        if config["min_circularity"] > 1 or config["max_circularity"] > 1:
            raise ValueError("Circularity must be between 0 and 1 for C" + str(ch) + ".")
        if config["max_circularity"] < config["min_circularity"]:
            raise ValueError("C" + str(ch) + " maximum circularity is below its minimum circularity.")
        if config["threshold_mode"] == "Local Phansalkar" and config["local_radius"] < 1:
            raise ValueError("C" + str(ch) + " local threshold radius must be at least 1 pixel.")
        if config["clahe"] and config["threshold_mode"] != "Local Phansalkar":
            raise ValueError(
                "C" + str(ch) + " has segmentation-only CLAHE enabled, but CLAHE is intentionally restricted "
                "to Local Phansalkar mode so Manual/Display thresholds retain raw-intensity meaning."
            )


def validate_unattended_configs(configs, context):
    """Batch and headless runs cannot stop for a Brightness/Contrast choice."""
    for config in configs:
        if needs_mask(config) and config["threshold_mode"] == "Display minimum":
            raise ValueError(
                context + " cannot use a Display-minimum threshold because that value is chosen by eye on one image. "
                "For C" + str(config["channel"]) + " (" + config["name"] + "), click 'Capture display minimum' to turn it "
                "into a fixed Manual threshold, or choose Auto or Local Phansalkar."
            )


def validate_parent_rois(rois, width, height):
    validated = []
    for i in range(len(rois)):
        roi = rois[i]
        if roi is None or not roi.isArea():
            raise ValueError("Parent ROI " + str(i + 1) + " is not a closed area ROI.")
        bounds = roi.getBounds()
        if bounds.width <= 0 or bounds.height <= 0:
            raise ValueError("Parent ROI " + str(i + 1) + " has zero area.")
        if bounds.x < 0 or bounds.y < 0 or bounds.x + bounds.width > width or bounds.y + bounds.height > height:
            raise ValueError("Parent ROI " + str(i + 1) + " extends outside the image bounds.")
        if roi.getName() is None or to_text(roi.getName()).strip() == "":
            roi.setName("ROI_" + str(i + 1).zfill(3))
        validated.append(roi)
    return validated


def whole_image_roi(imp):
    roi = Roi(0, 0, imp.getWidth(), imp.getHeight())
    roi.setName("Whole_Image")
    return roi


def calibration_override_from_options(global_options):
    if not parse_bool(global_options.get("override_calibration"), False):
        return None
    pixel_width = parse_float_value(global_options.get("pixel_width"), "Pixel width", 0.0000000001)
    pixel_height = parse_float_value(global_options.get("pixel_height"), "Pixel height", 0.0000000001)
    unit = to_text(global_options.get("spatial_unit")).strip()
    if unit == "":
        unit = "pixel"
    return (pixel_width, pixel_height, unit)


def settings_to_properties(global_options, channel_configs):
    """Build the .properties content for a profile or a run record.
    channel_configs holds every channel that should be recorded."""
    props = Properties()
    props.setProperty("application.name", APP_NAME)
    props.setProperty("application.version", APP_VERSION)
    props.setProperty("application.imagej_version", IMAGEJ_VERSION)
    for spec in GLOBAL_SETTING_SPECS:
        value = global_options.get(spec[1], spec[3])
        if spec[2] == "bool":
            props.setProperty(spec[0], bool_text(parse_bool(value, spec[3])))
        else:
            props.setProperty(spec[0], to_text(value))
    for config in channel_configs:
        prefix = "channel." + str(config["channel"]) + "."
        for key in config:
            if key == "channel":
                continue
            value = config[key]
            props.setProperty(prefix + key, bool_text(value) if isinstance(value, bool) else to_text(value))
    return props


def write_settings_properties(path, global_options, channel_configs):
    props = settings_to_properties(global_options, channel_configs)
    stream = FileOutputStream(path)
    try:
        props.store(stream, APP_NAME + " profile")
    finally:
        stream.close()


def read_properties_file(path):
    props = Properties()
    stream = FileInputStream(path)
    try:
        props.load(stream)
    finally:
        stream.close()
    return props


def global_options_from_properties(props):
    options = default_global_options()
    for spec in GLOBAL_SETTING_SPECS:
        value = props.getProperty(spec[0])
        if value is None:
            continue
        if spec[2] == "bool":
            options[spec[1]] = parse_bool(value, spec[3])
        else:
            options[spec[1]] = to_text(value)
    return options


def channel_values_from_properties(props, channel_number):
    prefix = "channel." + str(channel_number) + "."
    values = {}
    for key in CHANNEL_BOOLEAN_KEYS:
        value = props.getProperty(prefix + key)
        if value is not None:
            values[key] = parse_bool(value)
    for key in CHANNEL_STRING_KEYS + CHANNEL_NUMERIC_KEYS:
        value = props.getProperty(prefix + key)
        if value is not None:
            values[key] = to_text(value)
    return values


def load_profile_settings(path):
    """Read a settings profile without the window.
    Returns (global_options, enabled_configs, all_configs)."""
    props = read_properties_file(path)
    global_options = global_options_from_properties(props)
    all_configs = []
    enabled = []
    for channel_number in range(1, MAX_CHANNELS + 1):
        values = channel_values_from_properties(props, channel_number)
        if len(values) == 0:
            continue
        config = normalize_channel_config(channel_number, values)
        all_configs.append(config)
        if config["enabled"]:
            enabled.append(config)
    return (global_options, enabled, all_configs)


# -----------------------------------------------------------------------------
# Binary-image helpers that do not depend on Fiji's global options
# -----------------------------------------------------------------------------

def binary_morphology(ip, operation, iterations, count):
    """Erode, dilate, open, or close a 0/255 mask (255 = foreground).

    ImageJ's Process > Binary commands read their iteration count, neighbor
    count, and edge padding from Process > Binary > Options, which differ
    between installations. This applies the same pixel operations with every
    parameter given explicitly, so a saved profile yields the same mask
    everywhere. Open is erosions followed by dilations; close is the reverse.
    """
    if iterations <= 0:
        return
    old_pad = Prefs.padEdges
    Prefs.padEdges = False
    try:
        if operation == "erode":
            for i in range(iterations):
                ip.erode(count, 0)
        elif operation == "dilate":
            for i in range(iterations):
                ip.dilate(count, 0)
        elif operation == "open":
            for i in range(iterations):
                ip.erode(count, 0)
            for i in range(iterations):
                ip.dilate(count, 0)
        elif operation == "close":
            for i in range(iterations):
                ip.dilate(count, 0)
            for i in range(iterations):
                ip.erode(count, 0)
        else:
            raise ValueError("Unsupported binary operation: " + to_text(operation))
    finally:
        Prefs.padEdges = old_pad


def is_two_level_mask(ip):
    """True when an 8-bit image holds only the values 0 and 255."""
    if ip.getBitDepth() != 8:
        return False
    histogram = ip.getHistogram()
    for value in range(1, 255):
        if histogram[value] != 0:
            return False
    return True


def pixel_fingerprint(imp):
    stats = imp.getStatistics(Measurements.MEAN | Measurements.STD_DEV | Measurements.MIN_MAX)
    return (stats.mean, stats.stdDev, stats.min, stats.max)


def apply_clahe(imp, block_size, histogram_bins, max_slope):
    """Apply CLAHE to an image in place (precise method, no mask).

    The CLAHE library is called directly with the same parameter mapping as
    the Process > Enhance Local Contrast (CLAHE) dialog. The menu command is
    not used because it fails when no image window is open, which is the
    normal situation in a batch, and IJ.run would not report that failure.
    The menu command remains as a fallback for installations where the
    library class cannot be imported.
    """
    block_radius = int((block_size - 1) / 2)
    bins = int(histogram_bins) - 1
    flat = None
    try:
        from mpicbg.ij.clahe import Flat
        flat = Flat
    except ANY_ERROR:
        flat = None
    imp.deleteRoi()
    if flat is not None:
        flat.getInstance().run(imp, block_radius, bins, float(max_slope), None, False)
        return
    if WindowManager.getIDList() is None:
        raise RuntimeError(
            "the CLAHE menu command needs at least one open image window and the CLAHE library could not be loaded directly"
        )
    IJ.run(
        imp, "Enhance Local Contrast (CLAHE)",
        "blocksize=" + str(block_size) + " histogram=" + str(histogram_bins) +
        " maximum=" + str(max_slope) + " mask=*None*"
    )
    # The command locks the image and unlocks it when it finishes. A lock that
    # is still held means it stopped part-way; IJ.run does not raise for that.
    if imp.isLocked():
        imp.unlock()
        raise RuntimeError("the CLAHE command did not finish; see the Log window")


def roi_inside_mask(roi, width, height):
    """Full-size 0/255 image that is 255 exactly on the pixels ImageJ treats
    as inside the ROI when it measures."""
    inside = ByteProcessor(width, height)
    inside.setValue(255.0)
    inside.fill(roi)
    inside.resetRoi()
    return inside


def roi_masks(roi, width, height):
    """Return (inside, edge band) for a parent ROI.

    inside is the full-size 0/255 image of the ROI. The edge band is the
    one-pixel band of ROI pixels that touch the outside of the ROI or the
    image border; an object touches the ROI edge when it overlaps this band.
    Both are built from the same mask ImageJ uses to measure inside the ROI,
    so object detection and area measurement agree on every boundary pixel.
    """
    inside = roi_inside_mask(roi, width, height)
    eroded = inside.duplicate()
    binary_morphology(eroded, "erode", 1, 1)
    band = inside.duplicate()
    band.copyBits(eroded, 0, 0, Blitter.SUBTRACT)
    return (inside, ImagePlus("CellQuant_edge_band", band))


def roi_edge_band(roi, width, height):
    return roi_masks(roi, width, height)[1]


def count_positive_pixels(ip, roi):
    """Number of 255-valued pixels of a mask inside an ROI."""
    imp = ImagePlus("CellQuant_count", ip)
    imp.setRoi(clone_roi(roi))
    ip.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
    stats = imp.getStatistics(Measurements.AREA | Measurements.LIMIT)
    imp.deleteRoi()
    return int(stats.pixelCount)


def to_binary_processor(imp):
    """Convert any grayscale image to a 0/255 mask where nonzero is positive."""
    ip = imp.getProcessor().duplicate()
    ip.resetRoi()
    if ip.getBitDepth() == 24:
        raise ValueError("The reference mask must be a grayscale image, not RGB.")
    if ip.getBitDepth() != 8:
        ip = ip.convertToByte(False)
    if ip.isInvertedLut():
        ip.invertLut()
    ip.threshold(0)
    return ip


def save_tiff(imp, path):
    ensure_dir(os.path.dirname(path))
    saver = FileSaver(imp)
    if imp.getStackSize() > 1:
        ok = saver.saveAsTiffStack(path)
    else:
        ok = saver.saveAsTiff(path)
    if not ok:
        raise RuntimeError("Could not save TIFF: " + to_text(path))


def save_roi_array(rois, path):
    """Save ROIs as a ROI Manager .zip set (readable by ROI Manager).

    The archive is written directly instead of through a ROI Manager, so a
    write error raises here rather than opening a dialog in the middle of a
    batch."""
    if len(rois) == 0:
        return
    used = set()
    stream = ZipOutputStream(BufferedOutputStream(FileOutputStream(path)))
    try:
        data = DataOutputStream(BufferedOutputStream(stream))
        encoder = RoiEncoder(data)
        for index in range(len(rois)):
            roi = rois[index]
            if roi is None:
                continue
            label = roi_name(roi, "ROI_" + str(index + 1).zfill(4))
            label = label.replace("/", "_").replace("\\", "_")
            unique = label
            suffix = 1
            while unique in used:
                unique = label + "-" + str(suffix)
                suffix += 1
            used.add(unique)
            stream.putNextEntry(ZipEntry(unique + ".roi"))
            encoder.write(roi)
            data.flush()
            stream.closeEntry()
        data.flush()
    finally:
        stream.close()


def load_rois_from_file(path):
    """Read every ROI from a ROI Manager .zip set or a single .roi file.
    Read directly, so a damaged file raises instead of opening a dialog."""
    path = to_text(path)
    if not os.path.isfile(path):
        raise ValueError("ROI file does not exist: " + path)
    rois = []
    if path.lower().endswith(".zip"):
        stream = ZipInputStream(BufferedInputStream(FileInputStream(path)))
        try:
            entry = stream.getNextEntry()
            while entry is not None:
                name = to_text(entry.getName())
                if name.lower().endswith(".roi"):
                    buffer = ByteArrayOutputStream()
                    value = stream.read()
                    while value != -1:
                        buffer.write(value)
                        value = stream.read()
                    roi = RoiDecoder(buffer.toByteArray(), name[:-4]).getRoi()
                    if roi is not None:
                        if roi.getName() is None or to_text(roi.getName()).strip() == "":
                            roi.setName(name[:-4])
                        rois.append(roi)
                entry = stream.getNextEntry()
        finally:
            stream.close()
        if len(rois) == 0:
            raise ValueError("The archive contains no .roi entries: " + path)
    else:
        roi = RoiDecoder(path).getRoi()
        if roi is None:
            raise ValueError("Could not read ROI file: " + path)
        rois.append(roi)
    return rois


def apply_lut(imp, lut_name):
    if lut_name == "Original":
        return
    colors = {
        "Grays": Color.white,
        "Red": Color.red,
        "Green": Color.green,
        "Blue": Color.blue,
        "Magenta": Color.magenta,
        "Cyan": Color.cyan,
        "Yellow": Color.yellow
    }
    color = colors.get(lut_name, Color.white)
    try:
        imp.setLut(LUT.createLutFromColor(color))
        imp.updateAndDraw()
    except ANY_ERROR:
        try:
            IJ.run(imp, lut_name, "")
        except ANY_ERROR:
            pass


# -----------------------------------------------------------------------------
# Analysis engine (no user interface)
# -----------------------------------------------------------------------------

class AnalysisCancelled(Exception):
    pass


class AnalysisEngine(object):
    """Image preparation, segmentation, and measurement.

    Every setting arrives as an argument; nothing here reads a window. The
    window, the batch runner, and the tests all call this one class.
    """

    def __init__(self, log=None, is_cancelled=None):
        self._log = log
        self._is_cancelled = is_cancelled
        # Object outlines normally come back through a hidden ROI Manager.
        # A ROI Manager is a window class and cannot be created when Fiji
        # runs without a display, so there the outlines are collected from
        # an image overlay instead. Both routes give the same outlines.
        self.collect_rois_from_overlay = is_headless()

    def log(self, message):
        if self._log is not None:
            self._log(message)

    def show_status(self, message):
        try:
            IJ.showStatus(message)
        except ANY_ERROR:
            pass

    def check_cancelled(self):
        if self._is_cancelled is not None and self._is_cancelled():
            raise AnalysisCancelled("Cancelled by user.")

    # ----- preparation --------------------------------------------------------

    def prepare_single_channel(self, source, config, frame, calibration=None):
        channel = config["channel"]
        total_z = source.getNSlices()
        projection = config["projection"]

        if projection == "Current Z":
            z1 = source.getZ()
            z2 = z1
        elif projection == "Specific Z":
            z1 = config["z_slice"]
            z2 = z1
        else:
            z1 = 1
            z2 = total_z

        source_roi = clone_roi(source.getRoi())
        try:
            source.deleteRoi()
            duplicate = Duplicator().run(source, channel, channel, z1, z2, frame, frame)
        finally:
            if source_roi is not None:
                source.setRoi(source_roi)

        if duplicate is None:
            raise RuntimeError("Could not duplicate C" + str(channel) + ".")

        result = duplicate
        method_map = {
            "Max Z": "max",
            "Average Z": "avg",
            "Sum Z": "sum",
            "Min Z": "min",
            "Standard deviation Z": "sd",
            "Median Z": "median"
        }
        if projection in method_map and duplicate.getStackSize() > 1:
            projected = ZProjector.run(duplicate, method_map[projection])
            if projected is None:
                safe_close_image(duplicate)
                raise RuntimeError(projection + " failed for C" + str(channel) + ".")
            safe_close_image(duplicate)
            result = projected

        base = sanitize_name(strip_extension(source.getTitle()))
        frame_tag = "" if source.getNFrames() == 1 else "_F" + str(frame).zfill(4)
        result.setTitle(base + frame_tag + "_C" + str(channel) + "_" + sanitize_name(config["name"]) + "_PREP")

        if calibration is not None:
            cal = result.getCalibration().copy()
            cal.pixelWidth = calibration[0]
            cal.pixelHeight = calibration[1]
            cal.setUnit(calibration[2])
            result.setCalibration(cal)

        if config["subtract_background"]:
            options = "rolling=" + str(config["background_radius"])
            if config["light_background"]:
                options += " light"
            IJ.run(result, "Subtract Background...", options)
        if config["median_radius"] > 0:
            IJ.run(result, "Median...", "radius=" + str(config["median_radius"]))
        if config["gaussian_sigma"] > 0:
            IJ.run(result, "Gaussian Blur...", "sigma=" + str(config["gaussian_sigma"]))

        apply_lut(result, config["lut"])
        result.deleteRoi()
        return result

    def prepare_channels(self, source, configs, frame, calibration=None):
        prepared = {}
        for config in configs:
            self.show_status("Preparing C" + str(config["channel"]) + ": " + config["name"])
            prepared[config["channel"]] = self.prepare_single_channel(source, config, frame, calibration)
        return prepared

    # ----- segmentation -------------------------------------------------------

    def build_masks(self, prepared_imp, config, title_prefix):
        # Work only on a duplicate. The prepared image remains untouched and is
        # still used for all intensity measurements. This is especially important
        # for segmentation-only CLAHE.
        work = prepared_imp.duplicate()
        work.deleteRoi()

        original_stats = prepared_imp.getStatistics(Measurements.MIN_MAX)
        data_min = original_stats.min
        data_max = original_stats.max
        display_min = prepared_imp.getDisplayRangeMin()
        display_max = prepared_imp.getDisplayRangeMax()
        threshold_mode = config["threshold_mode"]

        threshold_low = ""
        threshold_high = ""
        segmentation_min = data_min
        segmentation_max = data_max
        scale_min = ""
        scale_max = ""

        if threshold_mode == "Local Phansalkar":
            # The local threshold works on an 8-bit copy. A deeper image is
            # rescaled from its own minimum and maximum to 0-255, then CLAHE
            # is optionally applied, then the local threshold. The rescaling
            # is done explicitly so that it does not depend on Fiji's
            # Edit > Options > Conversions "Scale when converting" setting.
            work_stats = work.getStatistics(Measurements.MIN_MAX)
            low = work_stats.min
            high = work_stats.max
            if high <= low:
                high = low + 1.0
            if work.getBitDepth() != 8:
                source_ip = work.getProcessor()
                source_ip.setMinAndMax(low, high)
                work.setProcessor(source_ip.convertToByte(True))
                scale_min = low
                scale_max = high

            if config["clahe"]:
                try:
                    apply_clahe(
                        work, config["clahe_block_size"], config["clahe_histogram_bins"], config["clahe_max_slope"]
                    )
                except ANY_ERROR as exc:
                    safe_close_image(work)
                    raise RuntimeError(
                        "CLAHE failed for C" + str(config["channel"]) + " (" + to_text(exc) +
                        "). Verify that Fiji's Enhance Local Contrast (CLAHE) plugin is available."
                    )

            seg_stats = work.getStatistics(Measurements.MIN_MAX)
            segmentation_min = seg_stats.min
            segmentation_max = seg_stats.max

            local_options = (
                "method=Phansalkar" +
                " radius=" + str(config["local_radius"]) +
                " parameter_1=" + str(config["local_parameter_1"]) +
                " parameter_2=" + str(config["local_parameter_2"]) +
                " white"
            )
            try:
                IJ.run(work, "Auto Local Threshold", local_options)
            except ANY_ERROR:
                safe_close_image(work)
                raise RuntimeError(
                    "Local Phansalkar thresholding failed for C" + str(config["channel"]) +
                    ". Verify that Fiji's Auto Local Threshold command is available."
                )
            # IJ.run reports a plugin failure in the Log window without raising,
            # so confirm that the result really is a two-level mask.
            if not is_two_level_mask(work.getProcessor()):
                safe_close_image(work)
                raise RuntimeError(
                    "Local Phansalkar thresholding did not produce a mask for C" + str(config["channel"]) +
                    ". Verify that Fiji's Auto Local Threshold command is available and check the Log window."
                )

            raw = ImagePlus(title_prefix + "_RAW_MASK", work.getProcessor().duplicate())
            raw.setCalibration(prepared_imp.getCalibration().copy())
            safe_close_image(work)

        else:
            stats = work.getStatistics(Measurements.MIN_MAX)
            work_data_min = stats.min
            work_data_max = stats.max

            if threshold_mode == "Display minimum":
                threshold_low = display_min
                threshold_high = work_data_max
            elif threshold_mode == "Manual":
                threshold_low = config["manual_low"]
                threshold_high = work_data_max if config["manual_high"] == -1 else config["manual_high"]
            elif threshold_mode == "Auto":
                method = config["auto_method"]
                method_text = method + (" dark" if config["dark_background"] else "")
                # The automatic threshold is computed from the full data range
                # of the image, whatever Brightness/Contrast range is showing.
                work.getProcessor().resetMinAndMax()
                work.setAutoThreshold(method_text)
                threshold_low = work.getProcessor().getMinThreshold()
                threshold_high = work.getProcessor().getMaxThreshold()
                if threshold_low == ImageProcessor.NO_THRESHOLD or threshold_low != threshold_low:
                    safe_close_image(work)
                    raise RuntimeError("Automatic threshold failed for C" + str(config["channel"]) + " using " + method + ".")
                if config["dark_background"]:
                    threshold_low = threshold_low * config["auto_scale"]
                else:
                    threshold_high = threshold_high * config["auto_scale"]
            else:
                safe_close_image(work)
                raise ValueError("Unsupported threshold mode: " + threshold_mode)

            if threshold_low > work_data_max or threshold_high < work_data_min:
                # The whole threshold range lies outside the pixel values of
                # this image (a blank or very dim channel, a negative control).
                # Nothing is positive: the mask is empty, not an error.
                mask_processor = ByteProcessor(work.getWidth(), work.getHeight())
                self.log(
                    "C" + str(config["channel"]) + " (" + config["name"] + "): the threshold range " +
                    str(threshold_low) + " to " + str(threshold_high) + " lies outside the image values " +
                    str(work_data_min) + " to " + str(work_data_max) + "; nothing is positive."
                )
            else:
                threshold_low = max(work_data_min, threshold_low)
                threshold_high = min(work_data_max, threshold_high)
                if threshold_high < threshold_low:
                    safe_close_image(work)
                    raise ValueError("Threshold upper value is below the lower value for C" + str(config["channel"]) + ".")
                work.getProcessor().setThreshold(threshold_low, threshold_high, ImageProcessor.NO_LUT_UPDATE)
                mask_processor = work.createThresholdMask()
            raw = ImagePlus(title_prefix + "_RAW_MASK", mask_processor)
            raw.setCalibration(prepared_imp.getCalibration().copy())
            safe_close_image(work)

        # From here on a mask is defined by pixel value: 255 is foreground.
        raw_ip = raw.getProcessor()
        if raw_ip.isInvertedLut():
            raw_ip.invertLut()
        raw_ip.resetThreshold()

        count = config.get("morph_count", 1)
        binary_morphology(raw_ip, "open", config["open_iterations"], count)
        binary_morphology(raw_ip, "close", config["close_iterations"], count)
        binary_morphology(raw_ip, "erode", config["erode_iterations"], count)
        binary_morphology(raw_ip, "dilate", config["dilate_iterations"], count)

        old_black_background = Prefs.blackBackground
        Prefs.blackBackground = True
        try:
            if config["fill_holes"]:
                IJ.run(raw, "Fill Holes", "")

            final = raw.duplicate()
            if config["watershed"]:
                final.setTitle(title_prefix + "_WATERSHED_MASK")
                IJ.run(final, "Watershed", "")
            else:
                final.setTitle(title_prefix + "_FINAL_MASK")
        finally:
            Prefs.blackBackground = old_black_background

        raw.getProcessor().setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
        final.getProcessor().setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)

        return {
            "raw": raw,
            "final": final,
            "threshold_mode": threshold_mode,
            "threshold_low": threshold_low,
            "threshold_high": threshold_high,
            "data_min": data_min,
            "data_max": data_max,
            "display_min": display_min,
            "display_max": display_max,
            "segmentation_min": segmentation_min,
            "segmentation_max": segmentation_max,
            "scale_min": scale_min,
            "scale_max": scale_max
        }

    # ----- measurement --------------------------------------------------------

    def measure_intensity(self, imp, roi):
        imp.setRoi(clone_roi(roi))
        flags = (
            Measurements.AREA | Measurements.MEAN | Measurements.STD_DEV |
            Measurements.MIN_MAX | Measurements.MEDIAN | Measurements.CENTROID
        )
        stats = imp.getStatistics(flags)
        raw_stats = imp.getRawStatistics()
        imp.deleteRoi()
        return {
            "area": stats.area,
            "pixel_count": raw_stats.pixelCount,
            "mean": stats.mean,
            "std_dev": stats.stdDev,
            "min": stats.min,
            "max": stats.max,
            "median": stats.median,
            "x_centroid": stats.xCentroid,
            "y_centroid": stats.yCentroid,
            "integrated_density": stats.mean * stats.area,
            "raw_integrated_density": raw_stats.mean * raw_stats.pixelCount
        }

    def measure_positive_area(self, mask_imp, roi):
        mask_imp.setRoi(clone_roi(roi))
        mask_imp.getProcessor().setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
        stats = mask_imp.getStatistics(Measurements.AREA | Measurements.LIMIT)
        mask_imp.deleteRoi()
        return {
            "positive_area": stats.area,
            "positive_pixels": stats.pixelCount
        }

    def particle_size_pixels(self, config, calibration):
        min_size = config["min_size"]
        max_size = float("inf") if config["max_size"] == -1 else config["max_size"]
        if config["size_units"] == "Pixel^2":
            return (min_size, max_size)
        pixel_area = calibration.pixelWidth * calibration.pixelHeight
        if pixel_area <= 0:
            pixel_area = 1.0
        min_px = min_size / pixel_area
        max_px = float("inf") if math.isinf(max_size) else max_size / pixel_area
        return (min_px, max_px)

    def analyze_objects(self, mask_imp, parent_roi, config, masks=None):
        """Detect the objects of one mask inside one parent ROI.

        The mask is first cleared outside the true outline of the parent ROI,
        so objects are confined to the ROI itself and an object crossing the
        outline is clipped there. (Handing ParticleAnalyzer an image with a
        non-rectangular ROI is not enough: called this way it only honours the
        ROI's bounding rectangle.) Edge exclusion is then decided exactly:
        an object is on the edge when one of its pixels touches the outside
        of the ROI or the image border.

        masks is the result of roi_masks() for this parent ROI; pass it when
        several channels or frames use the same ROI, to avoid rebuilding it.

        Returns (objects, info). Each object is a dictionary holding its ROI
        and geometry; info counts detected, edge-touching, and excluded objects.
        """
        width = mask_imp.getWidth()
        height = mask_imp.getHeight()
        bounds = parent_roi.getBounds()
        if masks is None:
            masks = roi_masks(parent_roi, width, height)
        inside_ip = masks[0]
        edge_band = masks[1]

        work_ip = mask_imp.getProcessor().duplicate()
        work_ip.resetRoi()
        if parent_roi.getMask() is not None:
            work_ip.copyBits(inside_ip, 0, 0, Blitter.AND)
        work_imp = ImagePlus("CellQuant_particles", work_ip)
        work_imp.setCalibration(mask_imp.getCalibration().copy())

        result_table = ResultsTable()
        use_overlay = self.collect_rois_from_overlay
        temp_manager = None
        if use_overlay:
            options = ParticleAnalyzer.SHOW_OVERLAY_OUTLINES
        else:
            temp_manager = RoiManager(False)
            options = ParticleAnalyzer.SHOW_NONE | ParticleAnalyzer.ADD_TO_MANAGER
        if config["include_holes"]:
            options |= ParticleAnalyzer.INCLUDE_HOLES
        if config["composite_rois"]:
            options |= ParticleAnalyzer.COMPOSITE_ROIS
        if config["four_connected"]:
            options |= ParticleAnalyzer.FOUR_CONNECTED

        measurements = (
            Measurements.AREA | Measurements.MEAN | Measurements.STD_DEV |
            Measurements.MIN_MAX | Measurements.CENTROID | Measurements.PERIMETER |
            Measurements.SHAPE_DESCRIPTORS | Measurements.FERET |
            Measurements.INTEGRATED_DENSITY | Measurements.MEDIAN | Measurements.RECT
        )

        min_px, max_px = self.particle_size_pixels(config, mask_imp.getCalibration())
        work_imp.setRoi(Roi(bounds.x, bounds.y, bounds.width, bounds.height))
        work_ip.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)

        try:
            if not use_overlay:
                ParticleAnalyzer.setRoiManager(temp_manager)
            analyzer = ParticleAnalyzer(
                options, measurements, result_table,
                min_px, max_px, config["min_circularity"], config["max_circularity"]
            )
            ok = analyzer.analyze(work_imp)
            if not ok:
                raise RuntimeError("Particle analysis failed for C" + str(config["channel"]) + ".")
            if use_overlay:
                rois = []
                particle_overlay = work_imp.getOverlay()
                if particle_overlay is not None:
                    for roi in particle_overlay.toArray():
                        roi = clone_roi(roi)
                        # The overlay colors its outlines; a ROI Manager does not.
                        roi.setStrokeColor(None)
                        roi.setFillColor(None)
                        rois.append(roi)
                if len(rois) != result_table.size():
                    raise RuntimeError(
                        "Particle analysis returned " + str(len(rois)) + " outline(s) for " +
                        str(result_table.size()) + " measured object(s) for C" + str(config["channel"]) + "."
                    )
            else:
                rois = [clone_roi(roi) for roi in temp_manager.getRoisAsArray()]
        finally:
            work_imp.deleteRoi()
            # Never leave a manager registered: ParticleAnalyzer would force
            # the next particle analysis anywhere in Fiji to add to it.
            try:
                ParticleAnalyzer.setRoiManager(None)
            except ANY_ERROR:
                pass
            # Release the hidden manager with dispose(), never close():
            # RoiManager.close() clears ImageJ's reference to the ROI Manager
            # window, which would disconnect the user's own ROI Manager.
            if temp_manager is not None:
                try:
                    temp_manager.reset()
                    temp_manager.dispose()
                except ANY_ERROR:
                    pass

        objects = []
        edge_touching = 0
        excluded = 0
        for index in range(len(rois)):
            object_roi = rois[index]
            edge_band.setRoi(clone_roi(object_roi))
            touches = edge_band.getStatistics(Measurements.MIN_MAX).max > 0
            edge_band.deleteRoi()
            if touches:
                edge_touching += 1
                if config["exclude_edges"]:
                    excluded += 1
                    continue
            object_bounds = object_roi.getBounds()
            objects.append({
                "roi": object_roi,
                "area": rt_value(result_table, index, ["Area"]),
                "perimeter": rt_value(result_table, index, ["Perim.", "Perimeter"]),
                "circularity": rt_value(result_table, index, ["Circ.", "Circularity"]),
                "feret": rt_value(result_table, index, ["Feret"]),
                "min_feret": rt_value(result_table, index, ["MinFeret"]),
                "aspect_ratio": rt_value(result_table, index, ["AR"]),
                "roundness": rt_value(result_table, index, ["Round"]),
                "solidity": rt_value(result_table, index, ["Solidity"]),
                "x": rt_value(result_table, index, ["X", "XM"], object_bounds.x + object_bounds.width / 2.0),
                "y": rt_value(result_table, index, ["Y", "YM"], object_bounds.y + object_bounds.height / 2.0),
                "touches_edge": touches
            })

        info = {
            "detected": len(rois),
            "edge_touching": edge_touching,
            "excluded": excluded
        }
        return (objects, info)

    # ----- table layout -------------------------------------------------------

    def build_summary_headers(self, configs):
        headers = [
            "Run_ID", "Application_Version", "Source_Image", "Source_Path", "Frame", "Parent_ROI_Index",
            "Parent_ROI_Name", "ROI_Area", "ROI_Area_Unit", "Pixel_Width",
            "Pixel_Height", "Spatial_Unit"
        ]
        for config in configs:
            prefix = channel_column_prefix(config)
            headers += [
                prefix + "_Mean", prefix + "_StdDev", prefix + "_Min", prefix + "_Max",
                prefix + "_Median", prefix + "_Integrated_Density", prefix + "_Raw_Integrated_Density"
            ]
            if needs_positive_area(config):
                headers += [
                    prefix + "_Positive_Area", prefix + "_Positive_Area_Density",
                    prefix + "_Positive_Area_Fraction_Percent"
                ]
            if needs_objects(config):
                headers += [
                    prefix + "_Object_Count", prefix + "_Object_Density_Per_Area",
                    prefix + "_Object_Density_Per_1e6_Area", prefix + "_Object_Density_Per_mm2",
                    prefix + "_Total_Object_Area", prefix + "_Mean_Object_Area"
                ]
        # Columns added in 2.3 follow the 2.2 columns so older scripts that
        # read by position keep working.
        for config in configs:
            if needs_objects(config):
                prefix = channel_column_prefix(config)
                headers += [
                    prefix + "_Edge_Touching_Object_Count",
                    prefix + "_Excluded_Edge_Object_Count"
                ]
        return headers

    def build_object_headers(self, configs):
        headers = [
            "Run_ID", "Application_Version", "Frame", "Segmentation_Channel", "Segmentation_Name",
            "Parent_ROI_Index", "Parent_ROI_Name", "Object_Index_In_Parent",
            "Global_Object_ID", "Object_Area", "Object_Perimeter", "Object_Circularity",
            "Object_Feret", "Centroid_X", "Centroid_Y", "Bounds_X", "Bounds_Y",
            "Bounds_Width", "Bounds_Height"
        ]
        for config in configs:
            prefix = channel_column_prefix(config)
            headers += [
                prefix + "_Object_Mean", prefix + "_Object_StdDev", prefix + "_Object_Min",
                prefix + "_Object_Max", prefix + "_Object_Median",
                prefix + "_Object_Integrated_Density", prefix + "_Object_Raw_Integrated_Density"
            ]
        headers += [
            "Touches_ROI_Edge", "Object_MinFeret", "Object_Aspect_Ratio",
            "Object_Roundness", "Object_Solidity"
        ]
        return headers

    def build_threshold_headers(self):
        return [
            "Run_ID", "Application_Version", "Source_Image", "Frame", "Channel", "Channel_Name", "Analysis_Mode",
            "Threshold_Mode", "Automatic_Method", "Dark_Background", "Automatic_Multiplier",
            "Threshold_Low", "Threshold_High", "Data_Min", "Data_Max", "Display_Min", "Display_Max",
            "Local_Radius", "Local_Parameter_1", "Local_Parameter_2",
            "CLAHE_Enabled", "CLAHE_Block_Size", "CLAHE_Histogram_Bins", "CLAHE_Maximum_Slope",
            "Segmentation_Min", "Segmentation_Max",
            "Open_Iterations", "Close_Iterations", "Erode_Iterations", "Dilate_Iterations",
            "Fill_Holes", "Watershed", "Object_Size_Units", "Minimum_Object_Area",
            "Maximum_Object_Area", "Minimum_Circularity", "Maximum_Circularity",
            "Exclude_Edges", "Include_Holes", "Composite_ROIs", "Four_Connected",
            "Cross_Channel_Object_Measurement",
            "Binary_Neighbor_Count", "Binary_Pad_Edges",
            "Segmentation_8bit_Scale_Min", "Segmentation_8bit_Scale_Max", "ImageJ_Version"
        ]

    def threshold_row(self, run_id, source_title, frame, config, mask_data):
        return {
            "Run_ID": run_id,
            "Application_Version": APP_VERSION,
            "Source_Image": source_title,
            "Frame": frame,
            "Channel": config["channel"],
            "Channel_Name": config["name"],
            "Analysis_Mode": config["analysis_mode"],
            "Threshold_Mode": mask_data["threshold_mode"],
            "Automatic_Method": "Phansalkar" if config["threshold_mode"] == "Local Phansalkar" else config["auto_method"],
            "Dark_Background": config["dark_background"],
            "Automatic_Multiplier": config["auto_scale"],
            "Threshold_Low": mask_data["threshold_low"],
            "Threshold_High": mask_data["threshold_high"],
            "Data_Min": mask_data["data_min"],
            "Data_Max": mask_data["data_max"],
            "Display_Min": mask_data["display_min"],
            "Display_Max": mask_data["display_max"],
            "Local_Radius": config["local_radius"],
            "Local_Parameter_1": config["local_parameter_1"],
            "Local_Parameter_2": config["local_parameter_2"],
            "CLAHE_Enabled": config["clahe"],
            "CLAHE_Block_Size": config["clahe_block_size"],
            "CLAHE_Histogram_Bins": config["clahe_histogram_bins"],
            "CLAHE_Maximum_Slope": config["clahe_max_slope"],
            "Segmentation_Min": mask_data["segmentation_min"],
            "Segmentation_Max": mask_data["segmentation_max"],
            "Open_Iterations": config["open_iterations"],
            "Close_Iterations": config["close_iterations"],
            "Erode_Iterations": config["erode_iterations"],
            "Dilate_Iterations": config["dilate_iterations"],
            "Fill_Holes": config["fill_holes"],
            "Watershed": config["watershed"],
            "Object_Size_Units": config["size_units"],
            "Minimum_Object_Area": config["min_size"],
            "Maximum_Object_Area": config["max_size"],
            "Minimum_Circularity": config["min_circularity"],
            "Maximum_Circularity": config["max_circularity"],
            "Exclude_Edges": config["exclude_edges"],
            "Include_Holes": config["include_holes"],
            "Composite_ROIs": config["composite_rois"],
            "Four_Connected": config["four_connected"],
            "Cross_Channel_Object_Measurement": config["cross_measure"],
            "Binary_Neighbor_Count": config.get("morph_count", 1),
            "Binary_Pad_Edges": False,
            "Segmentation_8bit_Scale_Min": mask_data.get("scale_min", ""),
            "Segmentation_8bit_Scale_Max": mask_data.get("scale_max", ""),
            "ImageJ_Version": IMAGEJ_VERSION
        }

    def build_object_row(self, run_id, frame, roi_index, parent_name, global_object_id,
                         object_index, segmentation_config, obj, prepared, configs):
        object_roi = obj["roi"]
        bounds = object_roi.getBounds()
        area = obj["area"]
        centroid_x = obj["x"]
        centroid_y = obj["y"]

        # ResultsTable headings vary slightly among ImageJ releases. Recover
        # missing essential geometry directly from the object ROI rather than
        # exporting NaN silently. Valid ParticleAnalyzer area values are kept so
        # hole/composite-particle behavior remains consistent with its settings.
        if area != area or centroid_x != centroid_x or centroid_y != centroid_y:
            object_stats = self.measure_intensity(
                prepared[segmentation_config["channel"]], object_roi
            )
            if area != area:
                area = object_stats["area"]
            if centroid_x != centroid_x:
                centroid_x = object_stats["x_centroid"]
            if centroid_y != centroid_y:
                centroid_y = object_stats["y_centroid"]

        row = {
            "Run_ID": run_id,
            "Application_Version": APP_VERSION,
            "Frame": frame,
            "Segmentation_Channel": segmentation_config["channel"],
            "Segmentation_Name": segmentation_config["name"],
            "Parent_ROI_Index": roi_index + 1,
            "Parent_ROI_Name": parent_name,
            "Object_Index_In_Parent": object_index + 1,
            "Global_Object_ID": global_object_id,
            "Object_Area": area,
            "Object_Perimeter": obj["perimeter"],
            "Object_Circularity": obj["circularity"],
            "Object_Feret": obj["feret"],
            "Centroid_X": centroid_x,
            "Centroid_Y": centroid_y,
            "Bounds_X": bounds.x,
            "Bounds_Y": bounds.y,
            "Bounds_Width": bounds.width,
            "Bounds_Height": bounds.height,
            "Touches_ROI_Edge": 1 if obj["touches_edge"] else 0,
            "Object_MinFeret": obj["min_feret"],
            "Object_Aspect_Ratio": obj["aspect_ratio"],
            "Object_Roundness": obj["roundness"],
            "Object_Solidity": obj["solidity"]
        }

        for config in configs:
            prefix = channel_column_prefix(config)
            should_measure = (
                config["channel"] == segmentation_config["channel"] or
                segmentation_config["cross_measure"]
            )
            if should_measure:
                intensity = self.measure_intensity(prepared[config["channel"]], object_roi)
                row[prefix + "_Object_Mean"] = intensity["mean"]
                row[prefix + "_Object_StdDev"] = intensity["std_dev"]
                row[prefix + "_Object_Min"] = intensity["min"]
                row[prefix + "_Object_Max"] = intensity["max"]
                row[prefix + "_Object_Median"] = intensity["median"]
                row[prefix + "_Object_Integrated_Density"] = intensity["integrated_density"]
                row[prefix + "_Object_Raw_Integrated_Density"] = intensity["raw_integrated_density"]
            else:
                row[prefix + "_Object_Mean"] = ""
                row[prefix + "_Object_StdDev"] = ""
                row[prefix + "_Object_Min"] = ""
                row[prefix + "_Object_Max"] = ""
                row[prefix + "_Object_Median"] = ""
                row[prefix + "_Object_Integrated_Density"] = ""
                row[prefix + "_Object_Raw_Integrated_Density"] = ""
        return row

    # ----- one complete image -------------------------------------------------

    def analyze_image(self, source, source_path, configs, parent_rois, options, run_dir, frames,
                      prepared_cache=None, settings=None, retain_masks=False):
        """Run the full measurement of one image and write its run folder.

        options: save_prepared, save_masks, save_parent_rois, save_object_rois,
                 save_overlays (booleans) and calibration (None or
                 (pixel_width, pixel_height, unit)).
        prepared_cache: prepared images to reuse for the first frame, or None.
        settings: (global_options, channel_configs) recorded as settings.properties.
        Returns a dictionary with the table rows, headers, and file paths.
        """
        source_title = to_text(source.getTitle())
        base = sanitize_name(strip_extension(source_title))
        prepared_dir = os.path.join(run_dir, "Prepared")
        masks_dir = os.path.join(run_dir, "Masks")
        roi_dir = os.path.join(run_dir, "ROIs")
        overlay_dir = os.path.join(run_dir, "Overlays")
        ensure_dir(prepared_dir)
        ensure_dir(masks_dir)
        ensure_dir(roi_dir)
        ensure_dir(overlay_dir)

        run_id = os.path.basename(run_dir)
        run_start = time.time()
        self.log("Started analysis run: " + run_id)
        self.log(APP_NAME + " " + APP_VERSION + " on ImageJ " + IMAGEJ_VERSION)
        if settings is not None:
            write_settings_properties(os.path.join(run_dir, "settings.properties"), settings[0], settings[1])
        if options.get("save_parent_rois", True):
            save_roi_array(parent_rois, os.path.join(roi_dir, "Parent_ROIs_Input.zip"))

        calibration = options.get("calibration")
        width = source.getWidth()
        height = source.getHeight()

        summary_headers = self.build_summary_headers(configs)
        object_headers = self.build_object_headers(configs)
        threshold_headers = self.build_threshold_headers()
        summary_rows = []
        object_rows = []
        threshold_rows = []
        global_object_id = 0
        retained_masks = []
        edge_bands = {}
        frames = list(frames)

        for frame_position in range(len(frames)):
            frame = frames[frame_position]
            self.check_cancelled()
            self.show_status("Cell Quant Studio: frame " + str(frame_position + 1) + " of " + str(len(frames)))
            self.log("Processing frame " + str(frame) + " (" + str(frame_position + 1) + " of " + str(len(frames)) + ").")

            use_cache = prepared_cache is not None and frame_position == 0
            if use_cache:
                prepared = prepared_cache
            else:
                prepared = self.prepare_channels(source, configs, frame, calibration)
            frame_masks = []
            frame_ok = False

            try:
                if options.get("save_prepared", True):
                    for config in configs:
                        imp = prepared[config["channel"]]
                        save_tiff(imp, os.path.join(prepared_dir, sanitize_name(imp.getTitle()) + ".tif"))

                masks = {}
                for config in configs:
                    if not needs_mask(config):
                        continue
                    prefix = to_text(prepared[config["channel"]].getTitle()).replace("_PREP", "")
                    mask_data = self.build_masks(prepared[config["channel"]], config, prefix)
                    masks[config["channel"]] = mask_data
                    raw = mask_data["raw"]
                    final = mask_data["final"]
                    frame_masks.append(raw)
                    frame_masks.append(final)

                    if options.get("save_masks", True):
                        save_tiff(raw, os.path.join(masks_dir, sanitize_name(raw.getTitle()) + ".tif"))
                        save_tiff(final, os.path.join(masks_dir, sanitize_name(final.getTitle()) + ".tif"))

                    threshold_rows.append(self.threshold_row(run_id, source_title, frame, config, mask_data))

                object_roi_lists = {}
                object_overlays = {}
                for config in configs:
                    if needs_objects(config):
                        object_roi_lists[config["channel"]] = []
                        object_overlays[config["channel"]] = Overlay()
                        for parent_roi in parent_rois:
                            outline = clone_roi(parent_roi)
                            outline.setStrokeColor(Color.yellow)
                            outline.setStrokeWidth(1.5)
                            object_overlays[config["channel"]].add(outline)

                for roi_index in range(len(parent_rois)):
                    self.check_cancelled()
                    parent_roi = parent_rois[roi_index]
                    parent_name = roi_name(parent_roi, "ROI_" + str(roi_index + 1).zfill(3))
                    reference_stats = self.measure_intensity(prepared[configs[0]["channel"]], parent_roi)
                    roi_area = reference_stats["area"]
                    cal = prepared[configs[0]["channel"]].getCalibration()
                    unit = to_text(cal.getUnit()) if cal.getUnit() else u"pixel"

                    summary = {
                        "Run_ID": run_id,
                        "Application_Version": APP_VERSION,
                        "Source_Image": source_title,
                        "Source_Path": source_path,
                        "Frame": frame,
                        "Parent_ROI_Index": roi_index + 1,
                        "Parent_ROI_Name": parent_name,
                        "ROI_Area": roi_area,
                        "ROI_Area_Unit": unit + u"^2",
                        "Pixel_Width": cal.pixelWidth,
                        "Pixel_Height": cal.pixelHeight,
                        "Spatial_Unit": unit
                    }

                    for config in configs:
                        channel = config["channel"]
                        prefix = channel_column_prefix(config)
                        intensity = self.measure_intensity(prepared[channel], parent_roi)
                        summary[prefix + "_Mean"] = intensity["mean"]
                        summary[prefix + "_StdDev"] = intensity["std_dev"]
                        summary[prefix + "_Min"] = intensity["min"]
                        summary[prefix + "_Max"] = intensity["max"]
                        summary[prefix + "_Median"] = intensity["median"]
                        summary[prefix + "_Integrated_Density"] = intensity["integrated_density"]
                        summary[prefix + "_Raw_Integrated_Density"] = intensity["raw_integrated_density"]

                        if needs_positive_area(config):
                            positive = self.measure_positive_area(masks[channel]["raw"], parent_roi)
                            positive_area = positive["positive_area"]
                            fraction = positive_area / roi_area if roi_area > 0 else float("nan")
                            summary[prefix + "_Positive_Area"] = positive_area
                            summary[prefix + "_Positive_Area_Density"] = fraction
                            summary[prefix + "_Positive_Area_Fraction_Percent"] = fraction * 100.0

                        if needs_objects(config):
                            if roi_index not in edge_bands:
                                edge_bands[roi_index] = roi_masks(parent_roi, width, height)
                            objects, info = self.analyze_objects(
                                masks[channel]["final"], parent_roi, config, edge_bands[roi_index]
                            )
                            object_count = len(objects)
                            object_areas = []
                            for object_index in range(object_count):
                                global_object_id += 1
                                obj = objects[object_index]
                                object_roi = obj["roi"]
                                object_label = (
                                    "F" + str(frame).zfill(4) + "_C" + str(channel) + "_ROI" +
                                    str(roi_index + 1).zfill(3) + "_OBJ" + str(object_index + 1).zfill(5)
                                )
                                object_roi.setName(object_label)
                                object_roi_lists[channel].append(clone_roi(object_roi))

                                outline = clone_roi(object_roi)
                                outline.setStrokeColor(Color.red)
                                outline.setStrokeWidth(1.0)
                                object_overlays[channel].add(outline)

                                row = self.build_object_row(
                                    run_id, frame, roi_index, parent_name, global_object_id,
                                    object_index, config, obj, prepared, configs
                                )
                                object_rows.append(row)
                                if row.get("Object_Area", "") != "":
                                    object_areas.append(row["Object_Area"])

                            total_object_area = sum(object_areas) if len(object_areas) > 0 else 0.0
                            summary[prefix + "_Object_Count"] = object_count
                            summary[prefix + "_Object_Density_Per_Area"] = object_count / roi_area if roi_area > 0 else float("nan")
                            summary[prefix + "_Object_Density_Per_1e6_Area"] = object_count * 1000000.0 / roi_area if roi_area > 0 else float("nan")
                            summary[prefix + "_Object_Density_Per_mm2"] = (
                                object_count * 1000000.0 / roi_area if roi_area > 0 and is_micrometer_unit(unit) else ""
                            )
                            summary[prefix + "_Total_Object_Area"] = total_object_area
                            summary[prefix + "_Mean_Object_Area"] = total_object_area / object_count if object_count > 0 else 0.0
                            summary[prefix + "_Edge_Touching_Object_Count"] = info["edge_touching"]
                            summary[prefix + "_Excluded_Edge_Object_Count"] = info["excluded"]

                    summary_rows.append(summary)

                for config in configs:
                    channel = config["channel"]
                    if not needs_objects(config):
                        continue
                    if options.get("save_object_rois", True) and len(object_roi_lists[channel]) > 0:
                        roi_path = os.path.join(
                            roi_dir,
                            base + "_F" + str(frame).zfill(4) + "_C" + str(channel) + "_" + sanitize_name(config["name"]) + "_Objects.zip"
                        )
                        save_roi_array(object_roi_lists[channel], roi_path)
                    if options.get("save_overlays", True):
                        overlay_source = prepared[channel].duplicate()
                        overlay_source.setOverlay(object_overlays[channel])
                        flattened = overlay_source.flatten()
                        overlay_path = os.path.join(
                            overlay_dir,
                            base + "_F" + str(frame).zfill(4) + "_C" + str(channel) + "_" + sanitize_name(config["name"]) + "_Objects_Overlay.tif"
                        )
                        save_tiff(flattened, overlay_path)
                        safe_close_image(flattened)
                        safe_close_image(overlay_source)
                frame_ok = True
            finally:
                if retain_masks and frame_ok:
                    retained_masks += frame_masks
                else:
                    for imp in frame_masks + ([] if frame_ok else retained_masks):
                        safe_close_image(imp)
                if not use_cache:
                    for imp in prepared.values():
                        safe_close_image(imp)

        summary_path = os.path.join(run_dir, "roi_summary.csv")
        object_path = os.path.join(run_dir, "objects.csv")
        threshold_path = os.path.join(run_dir, "thresholds.csv")
        write_csv(summary_path, summary_headers, summary_rows)
        write_csv(object_path, object_headers, object_rows)
        write_csv(threshold_path, threshold_headers, threshold_rows)

        duration = time.time() - run_start
        self.log("Analysis completed in " + IJ.d2s(duration, 2) + " seconds.")
        self.log("Saved run: " + to_text(run_dir))

        return {
            "run_dir": run_dir,
            "run_id": run_id,
            "base": base,
            "summary_headers": summary_headers,
            "summary_rows": summary_rows,
            "object_headers": object_headers,
            "object_rows": object_rows,
            "threshold_headers": threshold_headers,
            "threshold_rows": threshold_rows,
            "mask_images": retained_masks,
            "frames": frames,
            "duration": duration
        }


def engine_options_from_global(global_options):
    return {
        "save_prepared": parse_bool(global_options.get("save_prepared"), True),
        "save_masks": parse_bool(global_options.get("save_masks"), True),
        "save_parent_rois": parse_bool(global_options.get("save_parent_rois"), True),
        "save_object_rois": parse_bool(global_options.get("save_object_rois"), True),
        "save_overlays": parse_bool(global_options.get("save_overlays"), True),
        "calibration": calibration_override_from_options(global_options)
    }


# -----------------------------------------------------------------------------
# Folder batch processing (no user interface)
# -----------------------------------------------------------------------------

def parse_extension_list(text):
    result = []
    for part in to_text(text).replace(";", ",").replace(" ", ",").split(","):
        part = part.strip().lower().lstrip("*").lstrip(".")
        if part != "" and part not in result:
            result.append(part)
    return result


def is_cellquant_output_dir(name):
    name = to_text(name)
    return "_CellQuant_" in name or name.startswith("CellQuant_Batch_")


def discover_batch_images(input_dir, extensions, recursive):
    """List image files in a folder, sorted by name. Result folders written by
    this application are skipped so saved masks are never re-analyzed."""
    found = []
    pending = [input_dir]
    while len(pending) > 0:
        folder = pending.pop(0)
        try:
            names = sorted(os.listdir(folder))
        except ANY_ERROR:
            continue
        subfolders = []
        for name in names:
            path = os.path.join(folder, name)
            lower = to_text(name).lower()
            if lower.startswith("."):
                continue
            if os.path.isdir(path):
                if recursive and not is_cellquant_output_dir(name):
                    subfolders.append(path)
            elif os.path.isfile(path):
                for extension in extensions:
                    if lower.endswith("." + extension):
                        found.append(path)
                        break
        pending = subfolders + pending
    return found


def roi_file_candidates(image_path, roi_dir=None):
    folder = roi_dir if roi_dir else os.path.dirname(image_path)
    name = os.path.basename(to_text(image_path))
    stems = [strip_extension(name)]
    if stems[0].lower().endswith(".ome"):
        stems.append(stems[0][:-4])
    stems.append(name)
    candidates = []
    for stem in stems:
        for suffix in ROI_FILE_SUFFIXES:
            candidates.append(os.path.join(folder, stem + suffix))
    return candidates


def find_roi_file(image_path, roi_dir=None):
    for candidate in roi_file_candidates(image_path, roi_dir):
        if os.path.isfile(candidate):
            return candidate
    return None


def is_native_image_file(path):
    lower = to_text(path).lower()
    for extension in NATIVE_IMAGE_EXTENSIONS:
        if lower.endswith("." + extension):
            return True
    return False


def open_with_bioformats_api(path, all_series=False, log=None):
    """Open a file through the Bio-Formats API without showing a dialog.
    Returns a list of images (possibly empty) or None if Bio-Formats is missing."""
    try:
        from loci.plugins import BF
    except ANY_ERROR:
        return None
    images = None
    if all_series:
        try:
            # 'in' is a Python keyword, so the package is reached by name.
            import loci.plugins as loci_plugins
            importer_options = getattr(getattr(loci_plugins, "in"), "ImporterOptions")()
            importer_options.setId(path)
            importer_options.setOpenAllSeries(True)
            images = BF.openImagePlus(importer_options)
        except ANY_ERROR as exc:
            if log is not None:
                log("Could not request every series (" + to_text(exc) + "); opening the first series only.")
            images = None
    if images is None:
        images = BF.openImagePlus(path)
    if images is None:
        return []
    return [images[i] for i in range(len(images))]


def open_image_series(path, force_bioformats=False, all_series=False, log=None):
    """Open one file for unattended processing.
    Returns a list of (series_label, ImagePlus); the label is empty for a
    single image."""
    use_bioformats_first = force_bioformats or all_series or not is_native_image_file(path)
    images = None
    if use_bioformats_first:
        try:
            images = open_with_bioformats_api(path, all_series, log)
        except ANY_ERROR as exc:
            if log is not None:
                log("Bio-Formats could not open " + to_text(os.path.basename(path)) + ": " + to_text(exc))
            images = None
    if images is None or len(images) == 0:
        imp = None
        # Send ImageJ's "cannot open" message to the Log window instead of a
        # dialog, so one unreadable file does not hold up an unattended run.
        try:
            IJ.redirectErrorMessages()
            imp = IJ.openImage(path)
        except ANY_ERROR:
            imp = None
        try:
            IJ.redirectErrorMessages(False)
        except ANY_ERROR:
            pass
        if imp is None and not use_bioformats_first:
            try:
                images = open_with_bioformats_api(path, all_series, log)
            except ANY_ERROR:
                images = None
        elif imp is not None:
            images = [imp]
    if images is None or len(images) == 0:
        return []
    if len(images) == 1:
        return [("", images[0])]
    result = []
    for index in range(len(images)):
        result.append(("series_" + str(index + 1).zfill(3), images[index]))
    return result


def batch_parent_rois(imp, image_path, policy, roi_dir, manager_rois):
    """Resolve the parent ROIs for one batch image.
    Returns (rois, description). rois is None when the image must be skipped."""
    if policy == BATCH_ROI_WHOLE:
        return ([whole_image_roi(imp)], "Whole image")
    if policy == BATCH_ROI_MANAGER:
        if manager_rois is None or len(manager_rois) == 0:
            raise ValueError("The ROI Manager set is empty.")
        return ([clone_roi(roi) for roi in manager_rois], "ROI Manager set")
    if policy in [BATCH_ROI_FILE_SKIP, BATCH_ROI_FILE_WHOLE]:
        roi_path = find_roi_file(image_path, roi_dir)
        if roi_path is not None:
            rois = [roi for roi in load_rois_from_file(roi_path) if roi is not None and roi.isArea()]
            if len(rois) > 0:
                return (rois, "ROI file " + to_text(os.path.basename(roi_path)))
            if policy == BATCH_ROI_FILE_SKIP:
                return (None, "ROI file " + to_text(os.path.basename(roi_path)) + " contains no area ROI")
        if policy == BATCH_ROI_FILE_WHOLE:
            return ([whole_image_roi(imp)], "Whole image (no ROI file)")
        if roi_path is None:
            return (None, "No ROI file found")
    raise ValueError("Unsupported batch ROI policy: " + to_text(policy))


BATCH_REPORT_HEADERS = [
    "Batch_ID", "Image_Index", "Image_File", "Image_Relative_Path", "Series", "Status",
    "Message", "Run_Folder", "ROI_Source", "Parent_ROI_Count", "Frames",
    "Summary_Rows", "Object_Rows", "Seconds"
]


def run_batch(input_dir, output_root, configs, all_configs, global_options, batch_options,
              manager_rois=None, log=None, progress=None, is_cancelled=None, file_list=None):
    """Analyze every image of a folder with one fixed configuration.

    batch_options: extensions (text), recursive, roi_policy, roi_dir,
                   force_bioformats, all_series.
    Each image gets its own run folder, identical to a single-image run. The
    batch folder also holds combined tables for all images, a per-image report,
    the settings used, and a log. Tables are written image by image, so an
    interrupted batch keeps everything finished so far.
    """
    log_lines = []

    def batch_log(message):
        line = "[" + time.strftime("%H:%M:%S") + "] " + to_text(message)
        log_lines.append(line)
        if log is not None:
            log(message)

    validate_unattended_configs(configs, "Batch processing")
    if len(configs) == 0:
        raise ValueError("No channels are enabled for analysis.")
    input_dir = os.path.abspath(to_text(input_dir))
    if not os.path.isdir(input_dir):
        raise ValueError("Batch input folder does not exist: " + input_dir)
    extensions = parse_extension_list(batch_options.get("extensions", DEFAULT_BATCH_EXTENSIONS))
    if len(extensions) == 0:
        raise ValueError("List at least one file extension for the batch.")
    recursive = parse_bool(batch_options.get("recursive"), False)
    policy = to_text(batch_options.get("roi_policy", BATCH_ROI_WHOLE))
    if policy not in BATCH_ROI_POLICIES:
        raise ValueError("Unsupported batch ROI policy: " + policy)
    roi_dir = to_text(batch_options.get("roi_dir", "")).strip()
    if roi_dir != "" and not os.path.isdir(roi_dir):
        raise ValueError("ROI folder does not exist: " + roi_dir)
    force_bioformats = parse_bool(batch_options.get("force_bioformats"), False)
    all_series = parse_bool(batch_options.get("all_series"), False)

    if file_list is None:
        files = discover_batch_images(input_dir, extensions, recursive)
    else:
        files = list(file_list)
    if len(files) == 0:
        raise ValueError("No image files with the listed extensions were found in " + input_dir)

    batch_dir = unique_dir(os.path.abspath(to_text(output_root)), "CellQuant_Batch_" + timestamp_text())
    batch_id = os.path.basename(batch_dir)
    write_settings_properties(os.path.join(batch_dir, "settings.properties"), global_options, all_configs)

    header_engine = AnalysisEngine()
    extra = ["Batch_ID", "Image_Relative_Path"]
    summary_out = CsvAppender(os.path.join(batch_dir, "batch_roi_summary.csv"), header_engine.build_summary_headers(configs) + extra)
    object_out = CsvAppender(os.path.join(batch_dir, "batch_objects.csv"), header_engine.build_object_headers(configs) + extra)
    threshold_out = CsvAppender(os.path.join(batch_dir, "batch_thresholds.csv"), header_engine.build_threshold_headers() + extra)
    report_out = CsvAppender(os.path.join(batch_dir, "batch_report.csv"), BATCH_REPORT_HEADERS)

    engine_options = engine_options_from_global(global_options)
    all_frames = to_text(global_options.get("frame_mode")) == "All frames"
    counts = {"ok": 0, "failed": 0, "skipped": 0}
    cancelled = False
    batch_start = time.time()
    batch_log("Batch " + batch_id + ": " + str(len(files)) + " file(s) from " + input_dir)
    batch_log(APP_NAME + " " + APP_VERSION + " on ImageJ " + IMAGEJ_VERSION)

    try:
        for file_index in range(len(files)):
            if is_cancelled is not None and is_cancelled():
                cancelled = True
                break
            path = files[file_index]
            try:
                relative = os.path.relpath(path, input_dir)
            except ANY_ERROR:
                relative = os.path.basename(path)
            if progress is not None:
                progress(file_index, len(files), os.path.basename(path))
            batch_log("Image " + str(file_index + 1) + " of " + str(len(files)) + ": " + to_text(relative))

            series = []
            open_error = ""
            try:
                series = open_image_series(path, force_bioformats, all_series, batch_log)
            except ANY_ERROR as exc:
                open_error = to_text(exc)
            if len(series) == 0:
                counts["failed"] += 1
                message = "Could not open the file" + (": " + open_error if open_error else ".")
                batch_log("  FAILED - " + message)
                report_out.append([{
                    "Batch_ID": batch_id, "Image_Index": file_index + 1,
                    "Image_File": os.path.basename(path), "Image_Relative_Path": relative,
                    "Status": "FAILED", "Message": message
                }])
                continue

            for series_label, imp in series:
                report = {
                    "Batch_ID": batch_id, "Image_Index": file_index + 1,
                    "Image_File": os.path.basename(path), "Image_Relative_Path": relative,
                    "Series": series_label
                }
                image_start = time.time()
                image_log = []

                def image_logger(message, _lines=image_log):
                    _lines.append("[" + time.strftime("%H:%M:%S") + "] " + to_text(message))

                try:
                    if is_cancelled is not None and is_cancelled():
                        raise AnalysisCancelled("Cancelled by user.")
                    if series_label != "":
                        imp.setTitle(strip_extension(os.path.basename(path)) + "_" + series_label)
                    validate_channel_configs(imp, configs)
                    rois, roi_note = batch_parent_rois(imp, path, policy, roi_dir, manager_rois)
                    report["ROI_Source"] = roi_note
                    if rois is None:
                        counts["skipped"] += 1
                        report["Status"] = "SKIPPED"
                        report["Message"] = roi_note
                        batch_log("  SKIPPED - " + roi_note)
                    else:
                        rois = validate_parent_rois(rois, imp.getWidth(), imp.getHeight())
                        frames = range(1, imp.getNFrames() + 1) if all_frames else [1]
                        run_dir = unique_run_dir(batch_dir, strip_extension(imp.getTitle()))
                        engine = AnalysisEngine(log=image_logger, is_cancelled=is_cancelled)
                        result = engine.analyze_image(
                            imp, path, configs, rois, engine_options, run_dir, frames,
                            prepared_cache=None, settings=(global_options, all_configs), retain_masks=False
                        )
                        write_text_file(os.path.join(run_dir, "run_log.txt"), [APP_NAME + " " + APP_VERSION] + image_log)
                        for rows in [result["summary_rows"], result["object_rows"], result["threshold_rows"]]:
                            for row in rows:
                                row["Batch_ID"] = batch_id
                                row["Image_Relative_Path"] = relative
                        summary_out.append(result["summary_rows"])
                        object_out.append(result["object_rows"])
                        threshold_out.append(result["threshold_rows"])
                        counts["ok"] += 1
                        report["Status"] = "OK"
                        report["Message"] = ""
                        report["Run_Folder"] = os.path.basename(run_dir)
                        report["Parent_ROI_Count"] = len(rois)
                        report["Frames"] = len(result["frames"])
                        report["Summary_Rows"] = len(result["summary_rows"])
                        report["Object_Rows"] = len(result["object_rows"])
                        batch_log(
                            "  OK - " + str(len(rois)) + " ROI(s), " + str(len(result["object_rows"])) +
                            " object(s), " + roi_note
                        )
                except AnalysisCancelled:
                    cancelled = True
                    report["Status"] = "CANCELLED"
                    report["Message"] = "Cancelled by user."
                except ANY_ERROR as exc:
                    counts["failed"] += 1
                    report["Status"] = "FAILED"
                    report["Message"] = to_text(exc)
                    batch_log("  FAILED - " + to_text(exc))
                    try:
                        IJ.log(traceback.format_exc())
                    except ANY_ERROR:
                        pass
                finally:
                    report["Seconds"] = round(time.time() - image_start, 2)
                    report_out.append([report])
                    safe_close_image(imp)
                    try:
                        imp.flush()
                    except ANY_ERROR:
                        pass
                if cancelled:
                    break
            if cancelled:
                # Close any series of this file that were opened but not reached.
                for series_label, imp in series:
                    safe_close_image(imp)
                break
    finally:
        summary_out.close()
        object_out.close()
        threshold_out.close()
        report_out.close()
        duration = time.time() - batch_start
        if cancelled:
            batch_log("Batch cancelled by user.")
        batch_log(
            "Batch finished in " + IJ.d2s(duration, 1) + " s: " + str(counts["ok"]) + " analyzed, " +
            str(counts["skipped"]) + " skipped, " + str(counts["failed"]) + " failed."
        )
        batch_log("Saved batch: " + to_text(batch_dir))
        write_text_file(os.path.join(batch_dir, "batch_log.txt"), [APP_NAME + " " + APP_VERSION] + log_lines)
        if progress is not None:
            progress(len(files), len(files), "")

    return {
        "batch_dir": batch_dir,
        "batch_id": batch_id,
        "files": len(files),
        "analyzed": counts["ok"],
        "skipped": counts["skipped"],
        "failed": counts["failed"],
        "cancelled": cancelled,
        "duration": duration
    }


# -----------------------------------------------------------------------------
# Validation against manual annotation (no user interface)
# -----------------------------------------------------------------------------

def annotation_points_from_rois(rois):
    """Collect the pixel positions of every point ROI in a list of ROIs."""
    points = []
    for roi in rois:
        if roi is None or roi.getType() != Roi.POINT:
            continue
        for point in roi.getContainedPoints():
            points.append((int(point.x), int(point.y)))
    return points


def match_points_to_objects(object_rois, points, width, height):
    """Compare detected objects with manually marked positions.

    An object holding at least one mark is a true positive. An object holding
    none is a false positive. A mark outside every object is missed. When one
    object holds several marks, the extra marks are counted as merged: the
    detection joined objects the annotator counted separately.
    """
    labels = FloatProcessor(width, height)
    for index in range(len(object_rois)):
        labels.setValue(float(index + 1))
        labels.fill(object_rois[index])
    labels.resetRoi()

    points_per_object = [0] * len(object_rois)
    point_labels = []
    for point in points:
        label = 0
        if 0 <= point[0] < width and 0 <= point[1] < height:
            label = int(round(labels.getf(point[0], point[1])))
        point_labels.append(label)
        if label > 0:
            points_per_object[label - 1] += 1

    true_positive = 0
    false_positive = 0
    merged = 0
    for count in points_per_object:
        if count == 0:
            false_positive += 1
        else:
            true_positive += 1
            merged += count - 1
    missed = 0
    for label in point_labels:
        if label == 0:
            missed += 1

    manual = len(points)
    detected = len(object_rois)
    precision = safe_ratio(true_positive, true_positive + false_positive)
    recall = safe_ratio(true_positive, manual)
    f1 = float("nan")
    if precision == precision and recall == recall and (precision + recall) > 0:
        f1 = 2.0 * precision * recall / (precision + recall)
    return {
        "manual": manual,
        "detected": detected,
        "true_positive": true_positive,
        "false_positive": false_positive,
        "missed": missed,
        "merged": merged,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "points_per_object": points_per_object,
        "point_labels": point_labels
    }


VALIDATION_SUMMARY_HEADERS = [
    "Run_ID", "Application_Version", "Source_Image", "Frame", "Channel", "Channel_Name",
    "Parent_ROI_Index", "Parent_ROI_Name", "Manual_Count", "Detected_Count",
    "True_Positive_Objects", "False_Positive_Objects", "Missed_Annotations", "Merged_Annotations",
    "Precision", "Recall", "F1", "Count_Difference", "Count_Error_Percent"
]
VALIDATION_OBJECT_HEADERS = [
    "Run_ID", "Parent_ROI_Index", "Parent_ROI_Name", "Object_Index_In_Parent",
    "Centroid_X", "Centroid_Y", "Object_Area", "Annotations_Inside", "Classification"
]
VALIDATION_POINT_HEADERS = [
    "Run_ID", "Parent_ROI_Index", "Parent_ROI_Name", "Point_X_Pixel", "Point_Y_Pixel",
    "Matched_Object_Index", "Classification"
]
MASK_VALIDATION_HEADERS = [
    "Run_ID", "Application_Version", "Source_Image", "Frame", "Channel", "Channel_Name",
    "Parent_ROI_Index", "Parent_ROI_Name", "ROI_Pixels", "Test_Positive_Pixels",
    "Reference_Positive_Pixels", "Overlap_Pixels", "Dice", "Jaccard", "Precision", "Recall",
    "Test_Area_Fraction_Percent", "Reference_Area_Fraction_Percent", "Area_Difference_Percent"
]


def validation_summary_row(run_id, source_title, frame, config, roi_index, parent_name, match):
    return {
        "Run_ID": run_id,
        "Application_Version": APP_VERSION,
        "Source_Image": source_title,
        "Frame": frame,
        "Channel": config["channel"],
        "Channel_Name": config["name"],
        "Parent_ROI_Index": roi_index,
        "Parent_ROI_Name": parent_name,
        "Manual_Count": match["manual"],
        "Detected_Count": match["detected"],
        "True_Positive_Objects": match["true_positive"],
        "False_Positive_Objects": match["false_positive"],
        "Missed_Annotations": match["missed"],
        "Merged_Annotations": match["merged"],
        "Precision": match["precision"],
        "Recall": match["recall"],
        "F1": match["f1"],
        "Count_Difference": match["detected"] - match["manual"],
        "Count_Error_Percent": safe_ratio(match["detected"] - match["manual"], match["manual"]) * 100.0
    }


def validate_objects_against_points(engine, prepared_imp, config, parent_rois, points,
                                    run_dir, source_title, frame):
    """Segment one channel with the current settings and score the detected
    objects against manual point annotations, parent ROI by parent ROI.
    Writes validation_summary.csv, validation_objects.csv,
    validation_points.csv, and a color-coded overlay."""
    if not needs_objects(config):
        raise ValueError(
            "Point annotations validate object detection. Set C" + str(config["channel"]) +
            " to 'Objects only' or 'Objects + area'."
        )
    run_id = os.path.basename(run_dir)
    width = prepared_imp.getWidth()
    height = prepared_imp.getHeight()
    prefix = to_text(prepared_imp.getTitle()).replace("_PREP", "") + "_VALIDATION"
    mask_data = engine.build_masks(prepared_imp, config, prefix)
    summary_rows = []
    object_rows = []
    point_rows = []
    overlay = Overlay()
    totals = {"manual": 0, "detected": 0, "true_positive": 0, "false_positive": 0, "missed": 0, "merged": 0}
    try:
        for roi_index in range(len(parent_rois)):
            parent_roi = parent_rois[roi_index]
            parent_name = roi_name(parent_roi, "ROI_" + str(roi_index + 1).zfill(3))
            objects, info = engine.analyze_objects(mask_data["final"], parent_roi, config)
            object_rois = [obj["roi"] for obj in objects]
            parent_points = [point for point in points if parent_roi.contains(point[0], point[1])]
            match = match_points_to_objects(object_rois, parent_points, width, height)
            summary_rows.append(validation_summary_row(run_id, source_title, frame, config, roi_index + 1, parent_name, match))
            for key in totals:
                totals[key] += match[key]

            outline = clone_roi(parent_roi)
            outline.setStrokeColor(Color.yellow)
            outline.setStrokeWidth(1.5)
            overlay.add(outline)
            for object_index in range(len(objects)):
                inside = match["points_per_object"][object_index]
                if inside == 0:
                    classification = "FALSE_POSITIVE"
                    color = Color.red
                elif inside == 1:
                    classification = "TRUE_POSITIVE"
                    color = Color.green
                else:
                    classification = "TRUE_POSITIVE_MERGED"
                    color = Color.orange
                object_outline = clone_roi(objects[object_index]["roi"])
                object_outline.setStrokeColor(color)
                object_outline.setStrokeWidth(1.0)
                overlay.add(object_outline)
                object_rows.append({
                    "Run_ID": run_id,
                    "Parent_ROI_Index": roi_index + 1,
                    "Parent_ROI_Name": parent_name,
                    "Object_Index_In_Parent": object_index + 1,
                    "Centroid_X": objects[object_index]["x"],
                    "Centroid_Y": objects[object_index]["y"],
                    "Object_Area": objects[object_index]["area"],
                    "Annotations_Inside": inside,
                    "Classification": classification
                })
            for point_index in range(len(parent_points)):
                point = parent_points[point_index]
                label = match["point_labels"][point_index]
                marker = OvalRoi(point[0] - 2, point[1] - 2, 5, 5)
                marker.setStrokeColor(Color.cyan if label > 0 else Color.magenta)
                marker.setStrokeWidth(1.0)
                overlay.add(marker)
                point_rows.append({
                    "Run_ID": run_id,
                    "Parent_ROI_Index": roi_index + 1,
                    "Parent_ROI_Name": parent_name,
                    "Point_X_Pixel": point[0],
                    "Point_Y_Pixel": point[1],
                    "Matched_Object_Index": label if label > 0 else "",
                    "Classification": "MATCHED" if label > 0 else "MISSED"
                })

        if len(parent_rois) > 1:
            precision = safe_ratio(totals["true_positive"], totals["true_positive"] + totals["false_positive"])
            recall = safe_ratio(totals["true_positive"], totals["manual"])
            f1 = float("nan")
            if precision == precision and recall == recall and (precision + recall) > 0:
                f1 = 2.0 * precision * recall / (precision + recall)
            pooled = dict(totals)
            pooled["precision"] = precision
            pooled["recall"] = recall
            pooled["f1"] = f1
            summary_rows.append(validation_summary_row(run_id, source_title, frame, config, "", "ALL_ROIS_POOLED", pooled))

        base = sanitize_name(to_text(prepared_imp.getTitle()).replace("_PREP", ""))
        overlay_source = prepared_imp.duplicate()
        overlay_source.setOverlay(overlay)
        flattened = overlay_source.flatten()
        save_tiff(flattened, os.path.join(run_dir, base + "_Validation_Overlay.tif"))
        safe_close_image(flattened)
        safe_close_image(overlay_source)
        save_tiff(mask_data["final"], os.path.join(run_dir, sanitize_name(mask_data["final"].getTitle()) + ".tif"))
    finally:
        safe_close_image(mask_data["raw"])
        safe_close_image(mask_data["final"])

    write_csv(os.path.join(run_dir, "validation_summary.csv"), VALIDATION_SUMMARY_HEADERS, summary_rows)
    write_csv(os.path.join(run_dir, "validation_objects.csv"), VALIDATION_OBJECT_HEADERS, object_rows)
    write_csv(os.path.join(run_dir, "validation_points.csv"), VALIDATION_POINT_HEADERS, point_rows)
    write_csv(
        os.path.join(run_dir, "thresholds.csv"), engine.build_threshold_headers(),
        [engine.threshold_row(run_id, source_title, frame, config, mask_data)]
    )
    return {"headers": VALIDATION_SUMMARY_HEADERS, "rows": summary_rows, "run_dir": run_dir}


def validate_mask_against_reference(engine, prepared_imp, config, parent_rois, reference_imp,
                                    run_dir, source_title, frame):
    """Compare this channel's positive-area mask with a reference binary mask
    (nonzero = positive), pixel by pixel inside each parent ROI.
    Writes mask_validation_summary.csv."""
    if not needs_mask(config):
        raise ValueError(
            "A reference mask validates a segmented channel. Set C" + str(config["channel"]) +
            " to 'Positive area' or an Objects mode."
        )
    width = prepared_imp.getWidth()
    height = prepared_imp.getHeight()
    if reference_imp.getWidth() != width or reference_imp.getHeight() != height:
        raise ValueError(
            "The reference mask is " + str(reference_imp.getWidth()) + " x " + str(reference_imp.getHeight()) +
            " pixels but the image is " + str(width) + " x " + str(height) + "."
        )
    run_id = os.path.basename(run_dir)
    reference_ip = to_binary_processor(reference_imp)
    prefix = to_text(prepared_imp.getTitle()).replace("_PREP", "") + "_VALIDATION"
    mask_data = engine.build_masks(prepared_imp, config, prefix)
    rows = []
    try:
        test_ip = mask_data["raw"].getProcessor().duplicate()
        test_ip.resetRoi()
        overlap_ip = test_ip.duplicate()
        overlap_ip.copyBits(reference_ip, 0, 0, Blitter.AND)
        for roi_index in range(len(parent_rois)):
            parent_roi = parent_rois[roi_index]
            roi_pixels = engine.measure_intensity(prepared_imp, parent_roi)["pixel_count"]
            test_pixels = count_positive_pixels(test_ip, parent_roi)
            reference_pixels = count_positive_pixels(reference_ip, parent_roi)
            overlap_pixels = count_positive_pixels(overlap_ip, parent_roi)
            rows.append({
                "Run_ID": run_id,
                "Application_Version": APP_VERSION,
                "Source_Image": source_title,
                "Frame": frame,
                "Channel": config["channel"],
                "Channel_Name": config["name"],
                "Parent_ROI_Index": roi_index + 1,
                "Parent_ROI_Name": roi_name(parent_roi, "ROI_" + str(roi_index + 1).zfill(3)),
                "ROI_Pixels": roi_pixels,
                "Test_Positive_Pixels": test_pixels,
                "Reference_Positive_Pixels": reference_pixels,
                "Overlap_Pixels": overlap_pixels,
                "Dice": safe_ratio(2.0 * overlap_pixels, test_pixels + reference_pixels),
                "Jaccard": safe_ratio(overlap_pixels, test_pixels + reference_pixels - overlap_pixels),
                "Precision": safe_ratio(overlap_pixels, test_pixels),
                "Recall": safe_ratio(overlap_pixels, reference_pixels),
                "Test_Area_Fraction_Percent": safe_ratio(test_pixels, roi_pixels) * 100.0,
                "Reference_Area_Fraction_Percent": safe_ratio(reference_pixels, roi_pixels) * 100.0,
                "Area_Difference_Percent": safe_ratio(test_pixels - reference_pixels, reference_pixels) * 100.0
            })
        save_tiff(mask_data["raw"], os.path.join(run_dir, sanitize_name(mask_data["raw"].getTitle()) + ".tif"))
    finally:
        safe_close_image(mask_data["raw"])
        safe_close_image(mask_data["final"])
    write_csv(os.path.join(run_dir, "mask_validation_summary.csv"), MASK_VALIDATION_HEADERS, rows)
    write_csv(
        os.path.join(run_dir, "thresholds.csv"), engine.build_threshold_headers(),
        [engine.threshold_row(run_id, source_title, frame, config, mask_data)]
    )
    return {"headers": MASK_VALIDATION_HEADERS, "rows": rows, "run_dir": run_dir}


# -----------------------------------------------------------------------------
# Parameter sensitivity (no user interface)
# -----------------------------------------------------------------------------

SENSITIVITY_HEADERS = [
    "Run_ID", "Application_Version", "Source_Image", "Frame", "Channel", "Channel_Name",
    "Parameter", "Factor", "Parameter_Value", "Parent_ROI_Index", "Parent_ROI_Name", "ROI_Area",
    "Positive_Area", "Positive_Area_Fraction_Percent", "Positive_Area_Change_Percent",
    "Object_Count", "Object_Density_Per_Area", "Object_Count_Change_Percent", "Note"
]


def sensitivity_variant(config, parameter, factor, baseline):
    """Return (parameter value, modified configuration) for one sweep step.
    baseline is the mask data of the unmodified configuration."""
    variant = dict(config)
    label = "C" + str(config["channel"])
    if parameter == SENS_THRESHOLD:
        if config["threshold_mode"] == "Local Phansalkar":
            raise ValueError(
                label + " uses a local threshold, which has no single threshold level. "
                "Sweep the local radius or Phansalkar parameter 1 instead."
            )
        low = float(baseline["threshold_low"])
        high = float(baseline["threshold_high"])
        variant["threshold_mode"] = "Manual"
        if config["threshold_mode"] == "Auto" and not config["dark_background"]:
            variant["manual_low"] = low
            variant["manual_high"] = high * factor
            return (high * factor, variant)
        variant["manual_low"] = low * factor
        variant["manual_high"] = high
        return (low * factor, variant)
    if parameter == SENS_LOCAL_RADIUS:
        if config["threshold_mode"] != "Local Phansalkar":
            raise ValueError(label + " does not use Local Phansalkar thresholding.")
        value = max(1, int(round(config["local_radius"] * factor)))
        variant["local_radius"] = value
        return (value, variant)
    if parameter == SENS_LOCAL_K:
        if config["threshold_mode"] != "Local Phansalkar":
            raise ValueError(label + " does not use Local Phansalkar thresholding.")
        value = config["local_parameter_1"] * factor
        variant["local_parameter_1"] = value
        return (value, variant)
    if parameter == SENS_CLAHE_SLOPE:
        if not config["clahe"]:
            raise ValueError(label + " does not have segmentation-only CLAHE enabled.")
        value = config["clahe_max_slope"] * factor
        variant["clahe_max_slope"] = value
        return (value, variant)
    if parameter == SENS_MIN_SIZE:
        if not needs_objects(config):
            raise ValueError(label + " is not set to an Objects analysis mode.")
        value = config["min_size"] * factor
        variant["min_size"] = value
        return (value, variant)
    raise ValueError("Unsupported sensitivity parameter: " + to_text(parameter))


def run_sensitivity(engine, prepared_imp, config, parent_rois, parameter, factors,
                    run_dir, source_title, frame):
    """Repeat the segmentation of one channel while one parameter is scaled by
    each factor, and record positive area and object count per parent ROI.
    Factor 1.0 is the configured value; changes are reported relative to it.
    Writes sensitivity.csv."""
    if not needs_mask(config):
        raise ValueError("C" + str(config["channel"]) + " is not segmented; sensitivity needs a Positive area or Objects mode.")
    run_id = os.path.basename(run_dir)
    prefix = to_text(prepared_imp.getTitle()).replace("_PREP", "") + "_SENSITIVITY"
    width = prepared_imp.getWidth()
    height = prepared_imp.getHeight()

    baseline = engine.build_masks(prepared_imp, config, prefix)
    safe_close_image(baseline["raw"])
    safe_close_image(baseline["final"])

    edge_bands = {}
    results = {}
    notes = {}
    values = {}
    for factor in factors:
        engine.check_cancelled()
        try:
            value, variant = sensitivity_variant(config, parameter, factor, baseline)
            values[factor] = value
            mask_data = engine.build_masks(prepared_imp, variant, prefix)
        except (ValueError, RuntimeError) as exc:
            if factor == 1.0:
                raise
            notes[factor] = to_text(exc)
            engine.log("Sensitivity factor " + str(factor) + " skipped: " + to_text(exc))
            continue
        try:
            per_roi = []
            for roi_index in range(len(parent_rois)):
                parent_roi = parent_rois[roi_index]
                roi_area = engine.measure_intensity(prepared_imp, parent_roi)["area"]
                positive_area = engine.measure_positive_area(mask_data["raw"], parent_roi)["positive_area"]
                object_count = ""
                if needs_objects(variant):
                    if roi_index not in edge_bands:
                        edge_bands[roi_index] = roi_masks(parent_roi, width, height)
                    objects, info = engine.analyze_objects(mask_data["final"], parent_roi, variant, edge_bands[roi_index])
                    object_count = len(objects)
                per_roi.append((roi_area, positive_area, object_count))
            results[factor] = per_roi
        finally:
            safe_close_image(mask_data["raw"])
            safe_close_image(mask_data["final"])

    rows = []
    reference = results.get(1.0)
    for factor in factors:
        for roi_index in range(len(parent_rois)):
            parent_roi = parent_rois[roi_index]
            row = {
                "Run_ID": run_id,
                "Application_Version": APP_VERSION,
                "Source_Image": source_title,
                "Frame": frame,
                "Channel": config["channel"],
                "Channel_Name": config["name"],
                "Parameter": parameter,
                "Factor": factor,
                "Parameter_Value": values.get(factor, ""),
                "Parent_ROI_Index": roi_index + 1,
                "Parent_ROI_Name": roi_name(parent_roi, "ROI_" + str(roi_index + 1).zfill(3)),
                "Note": notes.get(factor, "")
            }
            if factor in results:
                roi_area, positive_area, object_count = results[factor][roi_index]
                row["ROI_Area"] = roi_area
                row["Positive_Area"] = positive_area
                row["Positive_Area_Fraction_Percent"] = safe_ratio(positive_area, roi_area) * 100.0
                if reference is not None:
                    row["Positive_Area_Change_Percent"] = safe_ratio(positive_area - reference[roi_index][1], reference[roi_index][1]) * 100.0
                if object_count != "":
                    row["Object_Count"] = object_count
                    row["Object_Density_Per_Area"] = safe_ratio(object_count, roi_area)
                    if reference is not None and reference[roi_index][2] != "":
                        row["Object_Count_Change_Percent"] = safe_ratio(object_count - reference[roi_index][2], reference[roi_index][2]) * 100.0
            rows.append(row)

    write_csv(os.path.join(run_dir, "sensitivity.csv"), SENSITIVITY_HEADERS, rows)
    return {"headers": SENSITIVITY_HEADERS, "rows": rows, "run_dir": run_dir}

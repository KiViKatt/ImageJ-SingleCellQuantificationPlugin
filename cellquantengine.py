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
from ij.plugin.filter import ParticleAnalyzer, Analyzer
from ij.measure import ResultsTable, Measurements
from ij.gui import Roi, Overlay, OvalRoi
from ij.process import ImageProcessor, LUT, AutoThresholder, ByteProcessor, FloatProcessor, Blitter
from ij.io import FileSaver, RoiDecoder, RoiEncoder

from java.awt import Color
from java.io import (
    FileInputStream, FileOutputStream, BufferedInputStream, BufferedOutputStream,
    ByteArrayOutputStream, DataOutputStream
)
from java.lang import Throwable as JavaThrowable
from java.lang import Runtime
from java.util import Properties
from java.util.zip import ZipInputStream, ZipOutputStream, ZipEntry

import os
import time
import threading
import traceback
import math
import codecs
import binascii
import re


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
    long_type = long
except NameError:
    number_types = (int, float)
    long_type = int

# Windows programs (Explorer, Excel, Python and R on a default Windows
# installation) cannot open a path longer than this many characters.
WINDOWS_PATH_LIMIT = 259


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
    if isinstance(value, float):
        # repr keeps every digit (1/3 -> 0.3333333333333333). Python 2's
        # unicode() keeps only 12 significant digits, so a setting written to
        # a profile would not come back exactly. Short values such as 0.25,
        # 35.0 or 1e-05 are written exactly as before.
        return text_type(repr(value))
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


CSV_INFINITY = float("inf")
CSV_CHUNK_ROWS = 2000


def csv_cell(value):
    """csv_escape(value), with a shortcut for the common types. The text is
    character for character what csv_escape gives (tested); the shortcut
    only avoids its calls, which dominate the time of a large table."""
    kind = type(value)
    if kind is float:
        if value != value or value == CSV_INFINITY or value == -CSV_INFINITY:
            return u'""'
        return u'"' + repr(value) + u'"'
    if kind is int or kind is long_type:
        return u'"' + str(value) + u'"'
    if kind is not text_type:
        return csv_escape(value)
    if u'"' in value:
        value = value.replace(u'"', u'""')
    return u'"' + value + u'"'


def csv_line(headers, row):
    get = row.get
    return u",".join([csv_cell(get(h, "")) for h in headers]) + u"\n"


def write_csv_rows(handle, headers, rows):
    """Write rows in chunks: one large write is much faster than many small
    ones through a codecs writer under Jython. Returns the number of rows."""
    chunk = []
    count = 0
    for row in rows:
        chunk.append(csv_line(headers, row))
        if len(chunk) >= CSV_CHUNK_ROWS:
            handle.write(u"".join(chunk))
            count += len(chunk)
            chunk = []
    if len(chunk) > 0:
        handle.write(u"".join(chunk))
        count += len(chunk)
    return count


def write_csv(path, headers, rows):
    """Write a UTF-8 CSV file (no byte-order mark)."""
    handle = codecs.open(path, "w", "utf-8")
    try:
        handle.write(u",".join([csv_escape(h) for h in headers]) + u"\n")
        write_csv_rows(handle, headers, rows)
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
        self.count += write_csv_rows(self.handle, self.headers, rows)
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


# Long image names. A run folder and the files in it are named after the
# image. With a long image name in a deep folder the longest path can exceed
# the Windows limit; then (and only then) the image-name part of the folder
# name and of the file names inside it is cut to its first characters plus
# an 8-character code of the whole name, for example
# "Control-day150-org1-SOX2-MAP2-DAPI-zst_1a2b3c4d". Every name that fits
# stays exactly as before. Names that are read back as input
# (<image>_ROIs.zip, <image>_points.zip) are never shortened.
SHORT_NAME_ROOM = 48
SHORT_NAME_MIN_ROOM = 16
# Longest name inside a run folder that does not contain the image name
# (for example ROIs/Tissue_Outline_ROIs.zip or mask_validation_summary.csv).
RUN_FOLDER_FIXED_ROOM = 40
# Room kept for the "_2", "_3" that unique_dir adds to a folder name in use.
RUN_FOLDER_COUNTER_ROOM = 2


def name_code(base):
    """Eight hexadecimal digits identifying a name (CRC-32 of its UTF-8)."""
    return "%08x" % (binascii.crc32(to_text(base).encode("utf-8")) & 0xffffffff)


def short_name_base(base, room=SHORT_NAME_ROOM):
    """base itself when it has at most room characters; otherwise its first
    characters and "_" + name_code(base), at most room characters in all."""
    base = to_text(base)
    if len(base) <= room:
        return base
    head = base[:max(7, room - 9)].rstrip("_.-")
    return head + "_" + name_code(base)


def run_folder_file_suffixes(configs=None, kind="", multi_frame=True, options=None):
    """(subfolder, suffix) of every file a run of these channel settings may
    write whose name starts with the image name: subfolder/<image name><suffix>.
    Prepared images and masks carry a frame tag only when the image has
    several frames (multi_frame); ROI sets and overlays always do. options
    (the global options) leaves out the files an analysis run is set not to
    save. Without configs, every kind of file of a channel C8 named
    Channel_8 is counted."""
    if not configs:
        configs = [{"channel": MAX_CHANNELS, "name": "Channel_" + str(MAX_CHANNELS),
                    "analysis_mode": "Objects + area", "watershed": True}]

    # An analysis run (kind "") writes prepared images, masks, ROI sets and
    # overlays; a validation run writes its own overlay and masks; a
    # sensitivity run writes no image file. Any other kind counts all.
    known = kind in ["", "Validation", "Sensitivity"]

    def saved(key):
        if not known:
            return True
        return kind == "" and (options is None or parse_bool(options.get(key), True))

    validation = kind == "Validation" or not known
    pairs = []
    # The check image of an automatic tissue outline (cellquanttissue.py).
    if (kind == "" or not known) and "cellquanttissue.py" in ENGINE_MODULE_VERSIONS and (
            options is None or parse_bool(options.get("tissue_save_qc"), True)):
        pairs.append(("Overlays", "_Tissue_Outline_QC.tif"))
    for config in configs:
        channel_part = "_C" + str(config["channel"]) + "_" + sanitize_name(config["name"])
        part = ("_F0001" if multi_frame else "") + channel_part
        mask = "_WATERSHED_MASK.tif" if config.get("watershed", True) else "_FINAL_MASK.tif"
        if saved("save_prepared"):
            pairs.append(("Prepared", part + "_PREP.tif"))
        if needs_mask(config):
            if saved("save_masks"):
                pairs.append(("Masks", part + mask))
            if validation:
                pairs += [("", part + "_VALIDATION" + mask), ("", part + "_VALIDATION_RAW_MASK.tif")]
        if needs_objects(config):
            if saved("save_object_rois"):
                pairs.append(("ROIs", "_F0001" + channel_part + "_Objects.zip"))
            if saved("save_overlays"):
                pairs.append(("Overlays", "_F0001" + channel_part + "_Objects_Overlay.tif"))
            if validation:
                pairs.append(("", part + "_Validation_Overlay.tif"))
    return pairs


def run_folder_inside_length(file_base, configs=None, fixed_names=None, inside=None, kind="", multi_frame=True,
                             options=None):
    """Length of the longest path below a run folder (from its separator on).
    inside replaces the (subfolder, suffix) pairs of an image run folder for
    a folder that holds other files (only fixed_names when it is [])."""
    longest = max([RUN_FOLDER_FIXED_ROOM] + [len(to_text(name)) for name in (fixed_names or [])])
    pairs = run_folder_file_suffixes(configs, kind, multi_frame, options) if inside is None else inside
    for sub, suffix in pairs:
        longest = max(longest, len(sub) + (1 if sub else 0) + len(file_base) + len(suffix))
    return 1 + longest


def plan_run_folder(root, base, kind="", configs=None, fixed_names=None, inside=None, multi_frame=True,
                    options=None):
    """Name a run folder for an image. Returns a dictionary: name (the
    folder name), file_base (the image-name part of the file names inside),
    base (the full sanitized image name), shortened, and longest (the
    expected length of the longest path). kind is "" for an analysis run,
    otherwise for example "Validation": <base>_CellQuant_<kind>_<time>.
    configs are the enabled channels, multi_frame tells whether the image
    has several frames, and options (the global options) which files an
    analysis saves, so that only files the run can write are counted.
    This is the one place where every run folder named after an image or a
    table is named: single runs, batches, validation, sensitivity, the
    window, and the modules."""
    base = sanitize_name(base)
    middle = "_CellQuant_" + (to_text(kind) + "_" if kind else "") + timestamp_text()
    root_length = len(to_text(os.path.abspath(to_text(root))))

    def longest(file_base):
        return (root_length + 1 + len(file_base) + len(middle) + RUN_FOLDER_COUNTER_ROOM +
                run_folder_inside_length(file_base, configs, fixed_names, inside, kind, multi_frame, options))

    file_base = base
    if longest(base) > WINDOWS_PATH_LIMIT:
        file_base = short_name_base(base, SHORT_NAME_ROOM)
        if longest(file_base) > WINDOWS_PATH_LIMIT:
            # The image name occurs twice in the longest path: once in the
            # folder name and once in the file name.
            room = max(SHORT_NAME_MIN_ROOM, (WINDOWS_PATH_LIMIT - longest(u"")) // 2)
            file_base = short_name_base(base, min(room, SHORT_NAME_ROOM))
    return {
        "name": file_base + middle, "file_base": file_base, "base": base,
        "shortened": file_base != base, "longest": longest(file_base)
    }


def unique_run_dir(root, base, configs=None, kind="", fixed_names=None, inside=None, multi_frame=True,
                   options=None):
    """Create the run folder of an image, <base>_CellQuant_<time>, below
    root (see plan_run_folder for long names) and return its path."""
    plan = plan_run_folder(root, base, kind, configs, fixed_names, inside, multi_frame, options)
    return unique_dir(root, plan["name"])


def run_file_base(run_dir, base, configs=None):
    """The image-name part of the file names inside run_dir. It is the
    full name in a folder that plan_run_folder named with the full name and
    the shortened name in one it shortened. In a folder named otherwise
    (by a caller), the full name unless the file names would not fit."""
    base = sanitize_name(base)
    folder = to_text(os.path.basename(os.path.normpath(to_text(run_dir))))
    if folder.startswith(base + "_CellQuant_"):
        return base
    marker = "_" + name_code(base)
    index = folder.find(marker + "_CellQuant_")
    if index > 0:
        return folder[:index + len(marker)]
    used = len(to_text(os.path.abspath(to_text(run_dir))))
    if used + run_folder_inside_length(base, configs) <= WINDOWS_PATH_LIMIT:
        return base
    shorter = short_name_base(base, SHORT_NAME_ROOM)
    if used + run_folder_inside_length(shorter, configs) > WINDOWS_PATH_LIMIT:
        room = max(SHORT_NAME_MIN_ROOM, WINDOWS_PATH_LIMIT - used - run_folder_inside_length(u"", configs))
        shorter = short_name_base(base, min(room, SHORT_NAME_ROOM))
    return shorter


def shortened_file_name(name, base, file_base):
    """A file name that starts with the image name, with that part replaced
    by file_base (unchanged when file_base is the full name)."""
    name = to_text(name)
    if file_base != base and name.startswith(base):
        return file_base + name[len(base):]
    return name


def shortened_name_note(base, file_base):
    return "Long image name shortened in file names: " + to_text(file_base) + " = " + to_text(base)


def run_folder_file_name(run_dir, source_title, config, name):
    """sanitize_name(name), shortened like the run folder of source_title
    when that folder had to be shortened (validation runs)."""
    base = sanitize_name(strip_extension(source_title))
    return shortened_file_name(sanitize_name(name), base, run_file_base(run_dir, base, [config]))


def log_run_folder_notes(engine, run_dir, source_title, config):
    """Record in the run's log a shortened image name and a result path that
    is too long for Windows programs (validation and sensitivity runs)."""
    base = sanitize_name(strip_extension(source_title))
    file_base = run_file_base(run_dir, base, [config])
    if file_base != base:
        engine.log(shortened_name_note(base, file_base))
    note = long_path_note(longest_path_length(run_dir))
    if note != "":
        engine.log(note)


def longest_path_length(folder):
    """Length of the longest file path below a folder."""
    longest = 0
    try:
        for root, folders, names in os.walk(folder):
            for name in names:
                longest = max(longest, len(os.path.join(root, name)))
    except ANY_ERROR:
        pass
    return longest


def long_path_note(longest):
    """The warning logged when a result path is too long for Windows
    programs, or "" when it is not (or the system is not Windows)."""
    if longest <= WINDOWS_PATH_LIMIT or os.sep != "\\":
        return ""
    return (
        "NOTE - the longest result path has " + str(longest) + " characters. Windows programs such as "
        "Excel, File Explorer, and Python or R scripts may refuse to open a path longer than " +
        str(WINDOWS_PATH_LIMIT) + " characters. If a result file will not open, move the results folder nearer "
        "to the drive root, or use a shorter output folder or shorter image file names."
    )


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
    if config["threshold_mode"] == "Auto" and config["auto_method"] not in AUTO_METHODS:
        # Found here rather than as an ImageJ error in the middle of a run.
        raise ValueError(
            label + " has an unknown automatic threshold method: " + config["auto_method"] +
            ". Choose one of: " + ", ".join(AUTO_METHODS) + " (spelled exactly so)."
        )
    return config


def needs_mask(config):
    return config["analysis_mode"] in ["Positive area", "Objects only", "Objects + area"]


def needs_objects(config):
    return config["analysis_mode"] in ["Objects only", "Objects + area"]


def needs_positive_area(config):
    return config["analysis_mode"] in ["Positive area", "Objects + area"]


def channel_column_prefix(config):
    return "C" + str(config["channel"]) + "_" + sanitize_name(config["name"])


def has_calibration_function(imp):
    """True when the image's pixel values are mapped to other values by a
    calibration function (for example signed 16-bit data, as in many
    DeltaVision .dv files). Statistics are then calibrated values, but a
    threshold acts on the stored values."""
    try:
        return bool(imp.getCalibration().calibrated())
    except ANY_ERROR:
        return False


def calibration_function_message(imp, config):
    return (
        "C" + str(config["channel"]) + " (" + to_text(config["name"]) + ") cannot be segmented: the image '" +
        to_text(imp.getTitle()) + "' has a calibration function (for example signed 16-bit data), so a threshold "
        "would act on the stored values while measurements report calibrated values. Convert the image to "
        "32-bit first (Image > Type > 32-bit, which applies the calibration), save it, and analyze that copy. "
        "Channels set to 'Intensity only' are not affected."
    )


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
    segmented = [config for config in configs if needs_mask(config)]
    if len(segmented) > 0 and has_calibration_function(source):
        raise ValueError(calibration_function_message(source, segmented[0]))

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
    check_frame_mode(options.get("frame_mode"))
    return options


def check_frame_mode(value):
    """An unknown frame mode would silently analyze the first frame only."""
    if to_text(value) not in FRAME_MODES:
        raise ValueError(
            "Unknown frame mode '" + to_text(value) + "' (setting global.frame_mode). Choose one of: " +
            ", ".join(FRAME_MODES) + "."
        )


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


def channel_is_switched_on(values):
    return parse_bool(values.get("enabled"), False) and to_text(values.get("analysis_mode", "Ignore")) != "Ignore"


def usable_channel_values(channel_number, values):
    """The values of a switched-off channel without those that cannot be
    read; each of those falls back to its default."""
    kept = {}
    for key in values:
        trial = {key: values[key]}
        if key == "auto_method":
            trial["threshold_mode"] = "Auto"
        try:
            normalize_channel_config(channel_number, trial)
            kept[key] = values[key]
        except ValueError:
            pass
    try:
        normalize_channel_config(channel_number, kept)
    except ValueError:
        # Values that are valid alone but not together (an unknown method
        # with an Auto threshold): keep only the name and the switch.
        kept = dict([(key, kept[key]) for key in ["name", "enabled"] if key in kept])
    return kept


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
        try:
            config = normalize_channel_config(channel_number, values)
        except ValueError:
            if channel_is_switched_on(values):
                raise
            # A channel that is switched off is not analyzed, so a bad value
            # there must not block the profile: that value gets its default.
            config = normalize_channel_config(channel_number, usable_channel_values(channel_number, values))
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


# CLAHE is by far the slowest step on a large image (about 20 s for 4084 x
# 4084 pixels on one core). The value of a pixel depends only on the pixels
# within the block radius around it, so the image can be cut into horizontal
# strips that overlap by that radius and are processed on several cores.
# The result is pixel-identical to processing the image in one piece.
CLAHE_STRIP_ROWS = 256
CLAHE_STRIP_LIMIT = 16


def clahe_strip_count(height):
    """In how many strips an image of this height is processed."""
    try:
        cores = int(Runtime.getRuntime().availableProcessors())
    except ANY_ERROR:
        cores = 1
    return max(1, min(cores, CLAHE_STRIP_LIMIT, int(height) // CLAHE_STRIP_ROWS))


def clahe_in_strips(flat, imp, block_radius, bins, max_slope, strips):
    """CLAHE on one 8-bit or 16-bit plane, strip by strip on worker
    threads. Returns False (and leaves the image alone) when the image is
    not of a kind this handles; raises what a strip raised."""
    source = imp.getProcessor()
    if imp.getStackSize() != 1 or source.getBitDepth() not in (8, 16):
        return False
    width = source.getWidth()
    height = source.getHeight()
    low = source.getMin()
    high = source.getMax()
    if not high > low:
        # No display range to hand to the strips (setMinAndMax(0, 0) would
        # mean "reset" there): leave such an image to the one-piece path.
        return False
    bounds = [int(round(height * index / float(strips))) for index in range(strips + 1)]
    pieces = [None] * strips
    errors = []

    def work(index):
        try:
            first = bounds[index]
            last = bounds[index + 1]
            top = max(0, first - block_radius)
            bottom = min(height, last + block_radius)
            part = source.duplicate()
            part.setRoi(0, top, width, bottom - top)
            part = part.crop()
            # A 16-bit plane is scaled by its display range inside CLAHE.
            part.setMinAndMax(low, high)
            part_imp = ImagePlus("CellQuant_clahe_strip", part)
            flat.getInstance().run(part_imp, block_radius, bins, float(max_slope), None, False)
            done = part_imp.getProcessor()
            done.setRoi(0, first - top, width, last - first)
            pieces[index] = (first, done.crop())
        except ANY_ERROR as exc:
            errors.append(exc)

    workers = [threading.Thread(target=work, args=(index,)) for index in range(strips)]
    for worker in workers:
        worker.start()
    for worker in workers:
        worker.join()
    if len(errors) > 0:
        raise errors[0]
    for piece in pieces:
        source.insert(piece[1], 0, piece[0])
    source.setMinAndMax(low, high)
    return True


def apply_clahe(imp, block_size, histogram_bins, max_slope, strips=None):
    """Apply CLAHE to an image in place (precise method, no mask).

    strips: None lets the image height and the number of processor cores
    decide; 1 processes the image in one piece. The pixels are the same
    either way.

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
        if strips is None:
            strips = clahe_strip_count(imp.getHeight())
        if strips > 1 and clahe_in_strips(flat, imp, block_radius, bins, max_slope, int(strips)):
            return
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


def run_imagej_command(imp, command, options):
    """IJ.run(imp, command, options), but refused with a clear message where
    this Fiji does not know the command (for example ImageJ started from
    plain Java without a display, where its menus are not set up). IJ.run
    would then leave the image unchanged and only write "Unrecognized
    command" to the Log, which would silently change the results."""
    IJ.run(imp, command, options)
    try:
        from ij import Menus
        commands = Menus.getCommands()
    except ANY_ERROR:
        commands = None
    if commands is not None and not commands.containsKey(command):
        raise RuntimeError(
            "The ImageJ command '" + command + "' is not available in this Fiji, so it could not be applied. "
            "Start Fiji normally (or with its own --headless option) and run again."
        )


def local_threshold_plugin():
    """An instance of Fiji's Auto Local Threshold plugin, or None."""
    try:
        from fiji.threshold import Auto_Local_Threshold
        return Auto_Local_Threshold()
    except ANY_ERROR:
        return None


def run_local_phansalkar(imp, radius, parameter_1, parameter_2):
    """Phansalkar local threshold of an 8-bit image in place (white objects).

    The plugin is called directly, as its menu command does after reading
    its dialog, so it also works where Fiji's menus are not loaded (a
    headless run). The numbers go through the same text the menu command
    would read ("parameter_1=0.25"), so both ways give the same pixels. The
    menu command is the fallback when the plugin class cannot be loaded."""
    radius_text = str(radius)
    first_text = str(parameter_1)
    second_text = str(parameter_2)
    plugin = local_threshold_plugin()
    if plugin is not None:
        plugin.exec(imp, "Phansalkar", int(float(radius_text)), float(first_text), float(second_text), True)
        return
    IJ.run(
        imp, "Auto Local Threshold",
        "method=Phansalkar radius=" + radius_text + " parameter_1=" + first_text +
        " parameter_2=" + second_text + " white"
    )


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


def roi_file_kind(path):
    """Return "zip", "roi", or "" from the first bytes of a file.

    The extension cannot be trusted: ROI Manager's Save command writes a
    plain .roi file when the set holds a single ROI, whatever name it is
    given, so a one-ROI "Parent_ROIs_Input.zip" written by version 2.2 is
    not an archive."""
    stream = FileInputStream(path)
    try:
        head = [stream.read(), stream.read(), stream.read(), stream.read()]
    finally:
        stream.close()
    if head[0] == 0x50 and head[1] == 0x4B:
        return "zip"
    if head == [0x49, 0x6F, 0x75, 0x74]:
        return "roi"
    return ""


def load_rois_from_file(path):
    """Read every ROI from a ROI Manager .zip set or a single .roi file.
    Read directly, so a damaged file raises instead of opening a dialog."""
    path = to_text(path)
    if not os.path.isfile(path):
        raise ValueError("ROI file does not exist: " + path)
    rois = []
    kind = roi_file_kind(path)
    if kind == "":
        kind = "zip" if path.lower().endswith(".zip") else "roi"
    if kind == "zip":
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
        # A file without a stored name is named after the file; drop the
        # extension of a one-ROI file that was saved with a .zip name.
        name = to_text(roi.getName()) if roi.getName() is not None else u""
        if name.strip() == "" or name == to_text(os.path.basename(path)):
            roi.setName(strip_extension(path))
        rois.append(roi)
    return rois


def apply_lut(imp, lut_name):
    # An RGB image has no lookup table; trying to set one only produced an
    # ImageJ error message (a dialog in the window) on every run.
    if lut_name == "Original" or imp.getBitDepth() == 24:
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
# Optional analysis modules
#
# Additional measurements (for example a depth profile from the ROI outline)
# live in their own source files beside this one. A module registers its
# settings and a RunExtension; analyze_image and run_batch then call it at
# fixed points. With no module loaded, or with a module's setting switched
# off, every table is exactly what it was before.
# -----------------------------------------------------------------------------

ENGINE_MODULE_FILES = [
    "cellquanttissue.py", "cellquantdepth.py", "cellquantspectrum.py",
    "cellquanttune.py", "cellquantcorrect.py", "cellquantrobust.py"
]
ENGINE_MODULE_VERSIONS = {}
RUN_EXTENSIONS = []
BATCH_ROI_PROVIDERS = {}
SETTINGS_FORMS = []


class RunExtension(object):
    """Base class of an optional per-run measurement.

    One instance is registered when its module loads, so it must keep no
    state of its own between runs; per-run state belongs in
    context["state"][self.name], a dictionary the engine creates.

    context keys available in every call: engine, run_dir, run_id,
    source_title, source_path, configs, parent_rois, options,
    global_options, width, height, frames, state, and get_roi_masks
    (a function: parent ROI index -> the roi_masks() pair, cached).
    Added for frame_done: frame, prepared (channel -> prepared image),
    masks (channel -> build_masks result), objects ((channel, parent ROI
    index) -> list of object dictionaries, each with its table row under
    "row"), and summary_rows (the roi_summary rows of this frame, in parent
    ROI order). Images in the context are closed after frame_done returns.
    """

    name = "extension"

    def is_enabled(self, global_options, configs):
        return False

    def summary_headers(self, configs, global_options):
        """Columns to append to roi_summary.csv; fill them in frame_done."""
        return []

    def object_headers(self, configs, global_options):
        """Columns to append to objects.csv; fill them in frame_done."""
        return []

    def begin_run(self, context):
        pass

    def frame_done(self, context):
        pass

    def end_run(self, context):
        """Return a list of (file name, headers, rows) tables. The engine
        writes each into the run folder, and a batch combines them into
        batch_<file name>."""
        return []


def register_run_extension(extension):
    """Register a RunExtension, replacing one of the same name."""
    for index in range(len(RUN_EXTENSIONS)):
        if RUN_EXTENSIONS[index].name == extension.name:
            RUN_EXTENSIONS[index] = extension
            return extension
    RUN_EXTENSIONS.append(extension)
    return extension


def register_global_settings(specs):
    """Add (property key, option key, kind, default) settings so they are
    saved in profiles and run records. Keys that already exist are kept."""
    known = set([spec[0] for spec in GLOBAL_SETTING_SPECS])
    for spec in specs:
        if spec[0] not in known:
            GLOBAL_SETTING_SPECS.append(tuple(spec))
            known.add(spec[0])


def register_engine_module(file_name, version):
    ENGINE_MODULE_VERSIONS[file_name] = version


def register_batch_roi_policy(name, provider):
    """Add a way of choosing the parent ROIs of each batch image.

    provider(request) returns (rois, description); rois is None when the
    image must be skipped. request is a dictionary with imp (the opened
    image), path, configs, global_options, and log (a function)."""
    if name not in BATCH_ROI_POLICIES:
        BATCH_ROI_POLICIES.append(name)
    BATCH_ROI_PROVIDERS[name] = provider


def register_settings_form(title, note, fields):
    """Describe the settings of a module so the window can show them.

    Plain data only; nothing here touches the user interface. fields is a
    list of (option key, label, kind, choices): kind is "bool", "text", or
    "choice" (choices is then the list of allowed values, otherwise None).
    Every option key must also be registered with register_global_settings."""
    form = {"title": title, "note": note, "fields": [tuple(field) for field in fields]}
    for index in range(len(SETTINGS_FORMS)):
        if SETTINGS_FORMS[index]["title"] == title:
            SETTINGS_FORMS[index] = form
            return form
    SETTINGS_FORMS.append(form)
    return form


def active_run_extensions(global_options, configs):
    if global_options is None:
        return []
    return [extension for extension in RUN_EXTENSIONS if extension.is_enabled(global_options, configs)]


def run_source_file(path, namespace):
    """Execute a source file in a namespace (Jython 2.7 and Python 3)."""
    try:
        runner = execfile
    except NameError:
        runner = None
    if runner is not None:
        runner(path, namespace)
        return
    handle = open(path, "r")
    try:
        source = handle.read()
    finally:
        handle.close()
    exec(compile(source, path, "exec"), namespace)


def load_engine_modules(namespace, folder):
    """Load every optional analysis module found in folder into namespace,
    which must be the namespace this engine was loaded into. Returns the
    names of the files loaded. A module from another release is refused."""
    loaded = []
    for file_name in ENGINE_MODULE_FILES:
        path = os.path.join(to_text(folder), file_name)
        if not os.path.isfile(path):
            continue
        # A module registers itself while it runs. When it is refused (or
        # fails), everything it registered and every name it defined is
        # taken back, so a caller that catches the error does not go on with
        # a half-registered module.
        snapshot = registry_snapshot(namespace)
        try:
            run_source_file(path, namespace)
            version = ENGINE_MODULE_VERSIONS.get(file_name)
            if version != APP_VERSION:
                raise RuntimeError(
                    "Cell Quant Studio: " + file_name + " is version " + to_text(version) +
                    " but cellquantengine.py is version " + APP_VERSION +
                    ". Re-extract the release so that all files come from one package."
                )
        except ANY_ERROR:
            registry_restore(namespace, snapshot)
            raise
        loaded.append(file_name)
    return loaded


ENGINE_REGISTRY_NAMES = [
    "GLOBAL_SETTING_SPECS", "RUN_EXTENSIONS", "BATCH_ROI_POLICIES", "BATCH_ROI_PROVIDERS",
    "SETTINGS_FORMS", "ENGINE_MODULE_VERSIONS"
]


def registry_snapshot(namespace):
    """Copies of the registries and of the namespace, for registry_restore."""
    contents = {}
    for name in ENGINE_REGISTRY_NAMES:
        registry = namespace.get(name)
        if isinstance(registry, dict):
            contents[name] = dict(registry)
        elif registry is not None:
            contents[name] = list(registry)
    return {"names": dict(namespace), "contents": contents}


def registry_restore(namespace, snapshot):
    """Put the registries (in place, so every reference sees it) and the
    names of namespace back to what registry_snapshot recorded."""
    names = snapshot["names"]
    for name in list(namespace.keys()):
        if name not in names:
            del namespace[name]
    for name in names:
        if namespace.get(name) is not names[name]:
            namespace[name] = names[name]
    for name in snapshot["contents"]:
        registry = namespace[name]
        saved = snapshot["contents"][name]
        if isinstance(registry, dict):
            registry.clear()
            registry.update(saved)
        else:
            registry[:] = saved


# -----------------------------------------------------------------------------
# Particle analysis helpers
# -----------------------------------------------------------------------------

# ImageJ keeps "Analyze > Set Measurements > Redirect to" in three private
# static fields of Analyzer. They are read and written back exactly, so a
# redirect the user chose survives an analysis unchanged.
REDIRECT_FIELDS = ["redirectTarget", "redirectTitle", "redirectImage"]


def analyzer_field(name):
    from java.lang import Class as JavaClass
    owner = JavaClass.forName("ij.plugin.filter.Analyzer", False, ImagePlus().getClass().getClassLoader())
    field = owner.getDeclaredField(name)
    field.setAccessible(True)
    return field


def save_measurement_redirect():
    """The current 'Redirect to' setting, for restore_measurement_redirect.
    None means no redirect is set."""
    if not Analyzer.isRedirectImage():
        return None
    try:
        return [analyzer_field(name).get(None) for name in REDIRECT_FIELDS]
    except ANY_ERROR:
        # Cannot happen with ImageJ 1.54; then the redirect is kept as it
        # is, which is what earlier versions did.
        return "unreadable"


def restore_measurement_redirect(saved):
    if saved == "unreadable":
        return
    if saved is None:
        Analyzer.setRedirectImage(None)
        return
    for index in range(len(REDIRECT_FIELDS)):
        analyzer_field(REDIRECT_FIELDS[index]).set(None, saved[index])


def run_particle_pass(settings, extra_options, measurements, table):
    """One ParticleAnalyzer pass over the bounding rectangle of the parent
    ROI of a prepared work mask (see analyze_objects). Returns its result."""
    work_imp, bounds, options, min_px, max_px, min_circularity, max_circularity = settings
    work_imp.setRoi(Roi(bounds.x, bounds.y, bounds.width, bounds.height))
    try:
        analyzer = ParticleAnalyzer(
            options | extra_options, measurements, table, min_px, max_px, min_circularity, max_circularity
        )
        return analyzer.analyze(work_imp)
    finally:
        work_imp.deleteRoi()


def particle_roi(roi, index, image_height):
    """An object outline taken from the analysis overlay, named and coloured
    exactly as when ParticleAnalyzer adds it to a ROI Manager (as 2.2 did):
    name <object number>-<y of the bounds centre> with at least 4 digits
    (RoiManager's label), and no stroke colour."""
    roi = clone_roi(roi)
    bounds = roi.getBounds()
    center = max(0, bounds.y + bounds.height // 2)
    number = index + 1
    digits = max(4, len(str(center)), len(str(number)))
    if digits == 4 and image_height >= 10000:
        digits = 5
    roi.setName(("000000" + str(number))[-digits:] + "-" + ("000000" + str(center))[-digits:])
    roi.setStrokeColor(None)
    return roi


OBJECT_INTENSITY_MEASUREMENTS = (
    Measurements.AREA | Measurements.MEAN | Measurements.STD_DEV | Measurements.MIN_MAX |
    Measurements.MEDIAN | Measurements.INTEGRATED_DENSITY | Measurements.RECT
)


OBJECT_INTENSITY_KEYS = [
    "mean", "std_dev", "min", "max", "median", "integrated_density", "raw_integrated_density"
]
OBJECT_INTENSITY_COLUMNS = [
    "_Object_Mean", "_Object_StdDev", "_Object_Min", "_Object_Max", "_Object_Median",
    "_Object_Integrated_Density", "_Object_Raw_Integrated_Density"
]


def same_pixel_geometry(first, second):
    a = first.getCalibration()
    b = second.getCalibration()
    return (a.pixelWidth == b.pixelWidth and a.pixelHeight == b.pixelHeight and
            a.xOrigin == b.xOrigin and a.yOrigin == b.yOrigin)


def particle_intensities(settings, image, objects_table, segmentation_channel, channel):
    """Per-object intensities of one prepared image, measured on exactly the
    pixels of each object.

    The objects are found again on the same work mask with the same
    settings, with "Redirect to" pointing at the prepared image, so
    ParticleAnalyzer measures each particle's own (flood-filled) pixels.
    Measuring the traced outline instead, as 2.2 did, also counts a hole in
    the object and any object inside that hole. Returns one tuple per
    object, in the order of objects_table: (mean, std_dev, min, max,
    median, integrated_density, raw_integrated_density). For an object
    without holes the values are bit for bit those of measure_intensity on
    its outline (tested)."""
    work_imp = settings[0]
    if image.getWidth() != work_imp.getWidth() or image.getHeight() != work_imp.getHeight():
        raise ValueError("The image of C" + str(channel) + " is not the size of the mask of C" + str(segmentation_channel) + ".")
    table = ResultsTable()
    # The redirect points at an undisplayed image sharing the pixels, so
    # ImageJ finds it whatever happens to a window showing the image.
    target = ImagePlus("CellQuant_measured", image.getProcessor())
    target.setCalibration(image.getCalibration())
    Analyzer.setRedirectImage(target)
    try:
        ok = run_particle_pass(settings, 0, OBJECT_INTENSITY_MEASUREMENTS, table)
    finally:
        Analyzer.setRedirectImage(None)
        image.getProcessor().resetRoi()
    count = objects_table.size()
    if not ok or table.size() != count:
        raise RuntimeError(
            "Per-object measurement of C" + str(channel) + " in the objects of C" + str(segmentation_channel) +
            " found " + str(table.size()) + " objects instead of " + str(count) + "."
        )
    # The second pass must have found the same objects in the same order:
    # their bounding rectangles (which do not depend on pixel values) agree.
    check = same_pixel_geometry(image, work_imp)
    columns = ["Mean", "StdDev", "Min", "Max", "Median", "IntDen", "RawIntDen"]
    values = []
    for row in range(count):
        if check:
            for name in ["BX", "BY", "Width", "Height"]:
                if table.getValue(name, row) != objects_table.getValue(name, row):
                    raise RuntimeError(
                        "Per-object measurement of C" + str(channel) + " did not line up with the objects of C" +
                        str(segmentation_channel) + " (object " + str(row + 1) + ")."
                    )
        values.append(tuple([rt_value(table, row, [name]) for name in columns]))
    return values


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
            run_imagej_command(result, "Subtract Background...", options)
        if config["median_radius"] > 0:
            run_imagej_command(result, "Median...", "radius=" + str(config["median_radius"]))
        if config["gaussian_sigma"] > 0:
            run_imagej_command(result, "Gaussian Blur...", "sigma=" + str(config["gaussian_sigma"]))

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
        empty_range = False

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

            try:
                run_local_phansalkar(
                    work, config["local_radius"], config["local_parameter_1"], config["local_parameter_2"]
                )
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
                # ImageJ reads "dark" relative to the lookup table: with an
                # inverted LUT it would select the dark pixels. "dark" here
                # always means bright objects, as in Manual mode, so the copy
                # gets a normal LUT (the pixel values do not change).
                if work.getProcessor().isInvertedLut():
                    work.getProcessor().invertLut()
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

            threshold_low = max(work_data_min, threshold_low)
            threshold_high = min(work_data_max, threshold_high)
            if threshold_high < threshold_low:
                # No pixel value lies in the range, for example a threshold
                # above the brightest pixel of a negative control. That is a
                # valid answer: nothing is positive and there are no objects.
                empty_range = True
                mask_processor = ByteProcessor(work.getWidth(), work.getHeight())
                self.log(
                    "C" + str(config["channel"]) + ": no pixel value lies in the threshold range (" +
                    to_text(threshold_low) + " to " + to_text(threshold_high) + "), so nothing is positive."
                )
            else:
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
                run_imagej_command(raw, "Fill Holes", "")

            final = raw.duplicate()
            if config["watershed"]:
                final.setTitle(title_prefix + "_WATERSHED_MASK")
                run_imagej_command(final, "Watershed", "")
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
            "scale_max": scale_max,
            # True when the threshold range held no pixel value of the image:
            # both masks are then empty (all background).
            "empty_range": empty_range
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

    def analyze_objects(self, mask_imp, parent_roi, config, masks=None, measure_images=None):
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

        measure_images maps channel numbers to prepared images. Each object
        then also gets "intensities": channel -> the per-object intensity
        values (measure_intensity's keys mean, std_dev, min, max, median,
        integrated_density, raw_integrated_density), measured on exactly the
        object's own pixels: not on a hole inside it, nor on an object
        inside such a hole (see particle_intensities).

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

        options = ParticleAnalyzer.SHOW_NONE
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
        work_ip.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
        pass_settings = (
            work_imp, bounds, options, min_px, max_px, config["min_circularity"], config["max_circularity"]
        )

        # "Redirect to" (Analyze > Set Measurements) is a global setting that
        # ParticleAnalyzer obeys: a redirect image of another size makes the
        # analysis fail. It is cleared for the analysis and put back after.
        redirect = save_measurement_redirect()
        intensities = {}
        try:
            Analyzer.setRedirectImage(None)
            # The objects are collected from an overlay of the work image,
            # not through a hidden ROI Manager, which cannot be created
            # without a display (a headless run). Same ROIs, same order.
            ok = run_particle_pass(pass_settings, ParticleAnalyzer.OVERLAY, measurements, result_table)
            if not ok:
                raise RuntimeError("Particle analysis failed for C" + str(config["channel"]) + ".")
            overlay = work_imp.getOverlay()
            rois = [] if overlay is None else [particle_roi(overlay.get(index), index, height) for index in range(overlay.size())]
            work_imp.setOverlay(None)
            if len(rois) != result_table.size():
                raise RuntimeError("Particle analysis of C" + str(config["channel"]) + " returned " +
                                   str(len(rois)) + " outlines for " + str(result_table.size()) + " objects.")
            for channel in sorted((measure_images or {}).keys()):
                self.check_cancelled()
                intensities[channel] = particle_intensities(
                    pass_settings, measure_images[channel], result_table, config["channel"], channel
                )
        finally:
            work_imp.deleteRoi()
            restore_measurement_redirect(redirect)

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
                "touches_edge": touches,
                "intensities": dict([(channel, intensities[channel][index]) for channel in intensities])
            })

        info = {
            "detected": len(rois),
            "edge_touching": edge_touching,
            "excluded": excluded
        }
        return (objects, info)

    # ----- table layout -------------------------------------------------------

    def build_summary_headers(self, configs, global_options=None):
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
        # Columns of optional analysis modules come last of all.
        for extension in active_run_extensions(global_options, configs):
            headers += list(extension.summary_headers(configs, global_options))
        return headers

    def build_object_headers(self, configs, global_options=None):
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
        for extension in active_run_extensions(global_options, configs):
            headers += list(extension.object_headers(configs, global_options))
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
                # Measured by analyze_objects on the object's own pixels.
                # Without that measurement (an object from elsewhere), the
                # outline is measured as before, holes included.
                values = obj.get("intensities", {}).get(config["channel"])
                if values is None:
                    intensity = self.measure_intensity(prepared[config["channel"]], object_roi)
                    values = tuple([intensity[key] for key in OBJECT_INTENSITY_KEYS])
                for position in range(len(OBJECT_INTENSITY_COLUMNS)):
                    row[prefix + OBJECT_INTENSITY_COLUMNS[position]] = values[position]
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
        # The image-name part of the file names in the run folder: the full
        # name unless it would make a path too long for Windows programs.
        file_base = run_file_base(run_dir, base, configs)

        def file_name(name):
            return shortened_file_name(sanitize_name(name), base, file_base)

        self.log("Started analysis run: " + run_id)
        self.log(APP_NAME + " " + APP_VERSION + " on ImageJ " + IMAGEJ_VERSION)
        if file_base != base:
            self.log(shortened_name_note(base, file_base))
        if settings is not None:
            write_settings_properties(os.path.join(run_dir, "settings.properties"), settings[0], settings[1])
        if options.get("save_parent_rois", True):
            save_roi_array(parent_rois, os.path.join(roi_dir, "Parent_ROIs_Input.zip"))

        calibration = options.get("calibration")
        width = source.getWidth()
        height = source.getHeight()

        global_options = options.get("global_options")
        if global_options is None and settings is not None:
            global_options = settings[0]
        extensions = active_run_extensions(global_options, configs)

        summary_headers = self.build_summary_headers(configs, global_options)
        object_headers = self.build_object_headers(configs, global_options)
        threshold_headers = self.build_threshold_headers()
        summary_rows = []
        object_rows = []
        threshold_rows = []
        global_object_id = 0
        retained_masks = []
        edge_bands = {}
        frames = list(frames)

        def get_roi_masks(roi_index):
            if roi_index not in edge_bands:
                edge_bands[roi_index] = roi_masks(parent_rois[roi_index], width, height)
            return edge_bands[roi_index]

        context = {
            "engine": self, "run_dir": run_dir, "run_id": run_id,
            "source_title": source_title, "source_path": source_path,
            "configs": configs, "parent_rois": parent_rois, "options": options,
            "global_options": global_options, "width": width, "height": height,
            "frames": frames, "state": {}, "get_roi_masks": get_roi_masks,
            # base: the image name as used in file names; file_base: what
            # file names in the run folder start with instead (the same
            # unless the name had to be shortened, see shortened_file_name).
            "base": base, "file_base": file_base
        }
        for extension in extensions:
            context["state"][extension.name] = {}
            extension.begin_run(context)

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
                        save_tiff(imp, os.path.join(prepared_dir, file_name(imp.getTitle()) + ".tif"))

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
                        save_tiff(raw, os.path.join(masks_dir, file_name(raw.getTitle()) + ".tif"))
                        save_tiff(final, os.path.join(masks_dir, file_name(final.getTitle()) + ".tif"))

                    threshold_rows.append(self.threshold_row(run_id, source_title, frame, config, mask_data))

                object_roi_lists = {}
                object_overlays = {}
                frame_objects = {}
                frame_summary_rows = []
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
                            # The fraction of this channel's own ROI area. It
                            # equals ROI_Area (measured on the first channel)
                            # unless a 32-bit image holds NaN (no number) in
                            # different pixels of different channels.
                            channel_area = intensity["area"]
                            fraction = positive_area / channel_area if channel_area > 0 else float("nan")
                            summary[prefix + "_Positive_Area"] = positive_area
                            summary[prefix + "_Positive_Area_Density"] = fraction
                            summary[prefix + "_Positive_Area_Fraction_Percent"] = fraction * 100.0

                        if needs_objects(config):
                            measured = [config] + ([other for other in configs if other is not config] if config["cross_measure"] else [])
                            objects, info = self.analyze_objects(
                                masks[channel]["final"], parent_roi, config, get_roi_masks(roi_index),
                                dict([(other["channel"], prepared[other["channel"]]) for other in measured])
                            )
                            frame_objects[(channel, roi_index)] = objects
                            object_count = len(objects)
                            object_areas = []
                            for object_index in range(object_count):
                                if object_index % 2000 == 1999:
                                    self.check_cancelled()
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
                                obj["row"] = row
                                # The values are in the row now.
                                obj.pop("intensities", None)
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
                    frame_summary_rows.append(summary)

                if len(extensions) > 0:
                    context["frame"] = frame
                    context["prepared"] = prepared
                    context["masks"] = masks
                    context["objects"] = frame_objects
                    context["summary_rows"] = frame_summary_rows
                    for extension in extensions:
                        self.check_cancelled()
                        extension.frame_done(context)

                for config in configs:
                    channel = config["channel"]
                    if not needs_objects(config):
                        continue
                    if options.get("save_object_rois", True) and len(object_roi_lists[channel]) > 0:
                        roi_path = os.path.join(
                            roi_dir,
                            file_base + "_F" + str(frame).zfill(4) + "_C" + str(channel) + "_" + sanitize_name(config["name"]) + "_Objects.zip"
                        )
                        save_roi_array(object_roi_lists[channel], roi_path)
                    if options.get("save_overlays", True):
                        overlay_source = prepared[channel].duplicate()
                        overlay_source.setOverlay(object_overlays[channel])
                        flattened = overlay_source.flatten()
                        overlay_path = os.path.join(
                            overlay_dir,
                            file_base + "_F" + str(frame).zfill(4) + "_C" + str(channel) + "_" + sanitize_name(config["name"]) + "_Objects_Overlay.tif"
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

        extra_tables = {}
        extra_table_names = []
        for key in ["frame", "prepared", "masks", "objects", "summary_rows"]:
            context.pop(key, None)
        for extension in extensions:
            for table in (extension.end_run(context) or []):
                write_csv(os.path.join(run_dir, table[0]), table[1], table[2])
                extra_tables[table[0]] = {"headers": list(table[1]), "rows": table[2]}
                extra_table_names.append(table[0])

        duration = time.time() - run_start
        self.log("Analysis completed in " + IJ.d2s(duration, 2) + " seconds.")
        path_note = long_path_note(longest_path_length(run_dir))
        if path_note != "":
            self.log(path_note)
        self.log("Saved run: " + to_text(run_dir))

        return {
            "run_dir": run_dir,
            "run_id": run_id,
            "base": base,
            "file_base": file_base,
            "summary_headers": summary_headers,
            "summary_rows": summary_rows,
            "object_headers": object_headers,
            "object_rows": object_rows,
            "threshold_headers": threshold_headers,
            "threshold_rows": threshold_rows,
            "extra_tables": extra_tables,
            "extra_table_names": extra_table_names,
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
        "calibration": calibration_override_from_options(global_options),
        "global_options": global_options
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


# Folders written by this application or its modules: CellQuant_Batch_<time>,
# CellQuant_Tuning_<time>, CellQuant_Robustness_<time>, and any other
# CellQuant_<Word>_<yyyymmdd_hhmmss> (with the _2, _3 of a folder in use).
CELLQUANT_FOLDER_PATTERN = re.compile(r"^CellQuant_[A-Za-z]+_[0-9]{8}_[0-9]{6}(_[0-9]+)?$")


def is_cellquant_output_dir(name):
    """True for a folder this application wrote: a run folder
    (<image>_CellQuant_<time>, also _Validation_, _Sensitivity_, or a
    table's _CellQuant_CountCorrection_ folder) or a folder such as
    CellQuant_Batch_<time>. A batch never looks inside these."""
    name = to_text(name)
    return "_CellQuant_" in name or name.startswith("CellQuant_Batch_") or CELLQUANT_FOLDER_PATTERN.match(name) is not None


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


def roi_file_candidates(image_path, roi_dir=None, series_label=""):
    """ROI file names tried for one image, best match first. A series of a
    multi-series file may have its own file, named after the image and the
    series (image_series_002_ROIs.zip); otherwise the file of the image is
    used for every series."""
    folder = roi_dir if roi_dir else os.path.dirname(image_path)
    name = os.path.basename(to_text(image_path))
    stems = [strip_extension(name)]
    if stems[0].lower().endswith(".ome"):
        stems.append(stems[0][:-4])
    stems.append(name)
    if series_label:
        stems = [stem + "_" + series_label for stem in stems] + stems
    candidates = []
    for stem in stems:
        for suffix in ROI_FILE_SUFFIXES:
            candidates.append(os.path.join(folder, stem + suffix))
    return candidates


def find_roi_file(image_path, roi_dir=None, series_label=""):
    for candidate in roi_file_candidates(image_path, roi_dir, series_label):
        if os.path.isfile(candidate):
            return candidate
    return None


def is_native_image_file(path):
    lower = to_text(path).lower()
    # An OME-TIFF may hold several series, which only Bio-Formats can tell
    # apart; ImageJ hands such a file to Bio-Formats anyway.
    if lower.endswith(".ome.tif") or lower.endswith(".ome.tiff"):
        return False
    for extension in NATIVE_IMAGE_EXTENSIONS:
        if lower.endswith("." + extension):
            return True
    return False


# Series that slide scanners and some microscopes store beside the data: a
# photograph of the slide label or a low-resolution overview.
THUMBNAIL_SERIES_NAMES = ["label image", "macro image", "thumbnail image", "thumbnail"]


def bioformats_importer_class(name):
    """Return a class of Bio-Formats' importer package, loci.plugins.in.
    'in' is a Python keyword, so the package cannot be written in an import
    statement, and under Jython the loci.plugins package has no attribute
    'in' until something else has loaded it. The class is therefore loaded
    by name through Java."""
    from java.lang import Class
    from loci.plugins import BF
    full_name = "loci.plugins.in." + name
    try:
        return Class.forName(full_name, True, BF.getClassLoader())
    except ANY_ERROR:
        return IJ.getClassLoader().loadClass(full_name)


def quiet_bioformats_log():
    """Bio-Formats prints progress lines ("Reading IFDs", ...) to the console
    through its logger. Limit it to warnings while a file is opened. Returns
    what restore_bioformats_log() needs to put the previous level back."""
    try:
        from org.slf4j import LoggerFactory
        from ch.qos.logback.classic import Level
        logger = LoggerFactory.getLogger("loci")
        if to_text(logger.getEffectiveLevel()) != "INFO":
            return None
        previous = logger.getLevel()
        logger.setLevel(Level.WARN)
        return (logger, previous)
    except ANY_ERROR:
        return None


def restore_bioformats_log(state):
    if state is not None:
        try:
            state[0].setLevel(state[1])
        except ANY_ERROR:
            pass


def bioformats_import(path, series_index=None, read_pixels=True):
    """Run the Bio-Formats importer on one file without any dialog.
    Returns (images, series). series describes every series of the file
    (index, name, width, height, channels, slices, frames, thumbnail);
    images holds the requested series, the first one when series_index is
    None, or nothing when read_pixels is False. The file is always released,
    also after a failed read: Bio-Formats alone keeps a damaged file open,
    which on Windows stops the user from moving or deleting it."""
    options = bioformats_importer_class("ImporterOptions")()
    options.setId(path)
    if series_index is not None:
        options.clearSeries()
        options.setSeriesOn(series_index, True)
    process = bioformats_importer_class("ImportProcess")(options)
    images = []
    series = []
    log_state = quiet_bioformats_log()
    try:
        if process.execute():
            reader = process.getReader()
            metadata = process.getOMEMetadata()
            for index in range(process.getSeriesCount()):
                reader.setSeries(index)
                name = ""
                try:
                    name = to_text(metadata.getImageName(index)).strip()
                except ANY_ERROR:
                    pass
                series.append({
                    "index": index, "name": name,
                    "width": reader.getSizeX(), "height": reader.getSizeY(),
                    "channels": reader.getSizeC(), "slices": reader.getSizeZ(),
                    "frames": reader.getSizeT(), "thumbnail": bool(reader.isThumbnailSeries())
                })
            reader.setSeries(0)
            if read_pixels:
                opened = bioformats_importer_class("ImagePlusReader")(process).openImagePlus()
                if opened is not None:
                    images = [opened[i] for i in range(len(opened))]
    finally:
        restore_bioformats_log(log_state)
        for getter in ["getReader", "getImageReader"]:
            try:
                getattr(process, getter)().close()
            except ANY_ERROR:
                pass
    return (images, series)


def is_thumbnail_series(info):
    return bool(info.get("thumbnail")) or to_text(info.get("name")).strip().lower() in THUMBNAIL_SERIES_NAMES


class BatchSeries(object):
    """One image of a batch file. A series of a multi-series file is read
    only when its turn comes, so the series of one file are never all in
    memory together."""

    def __init__(self, label, name="", image=None, opener=None, skip_reason=""):
        self.label = label
        self.name = name
        self.image = image
        self.opener = opener
        self.skip_reason = skip_reason

    def open(self):
        if self.image is None and self.opener is not None:
            opener = self.opener
            self.opener = None
            self.image = opener()
        return self.image

    def release(self):
        image = self.image
        self.image = None
        self.opener = None
        if image is not None:
            safe_close_image(image)
            try:
                image.flush()
            except ANY_ERROR:
                pass


def bioformats_series_opener(path, index):
    def opener():
        images, series = bioformats_import(path, index)
        return images[0] if len(images) > 0 else None
    return opener


def plan_bioformats_series(path, all_series=False, log=None, skip_thumbnails=True):
    """List what Bio-Formats would open from one file as BatchSeries entries.
    Returns None when Bio-Formats is not installed and raises when it cannot
    read the file."""
    try:
        from loci.plugins import BF
    except ANY_ERROR:
        return None
    if not all_series:
        images, series = bioformats_import(path)
        if len(series) > 1 and log is not None:
            log(
                "  NOTE - " + to_text(os.path.basename(path)) + " holds " + str(len(series)) +
                " series and only the first is analyzed. Turn on 'Analyze every series' to analyze them all."
            )
        first_name = series[0]["name"] if len(series) > 0 else ""
        if len(images) == 1:
            return [BatchSeries("", first_name, image=images[0])]
        return [
            BatchSeries("series_" + str(i + 1).zfill(3), first_name, image=images[i])
            for i in range(len(images))
        ]
    images, series = bioformats_import(path, read_pixels=False)
    entries = []
    for info in series:
        label = ""
        reason = ""
        if len(series) > 1:
            label = "series_" + str(info["index"] + 1).zfill(3)
            if skip_thumbnails and is_thumbnail_series(info):
                reason = (
                    "Not analyzed: '" + info["name"] + "' (" + str(info["width"]) + " x " + str(info["height"]) +
                    ") is a label, overview or thumbnail series, not microscope data."
                )
        entries.append(BatchSeries(
            label, info["name"], opener=bioformats_series_opener(path, info["index"]), skip_reason=reason
        ))
    return entries


def plan_image_series(path, force_bioformats=False, all_series=False, log=None, skip_thumbnails=True):
    """Find out what one file holds, for unattended processing.
    Returns a list of BatchSeries, empty when the file cannot be opened. A
    single image is opened right away; the series of a multi-series file
    are opened one at a time by BatchSeries.open()."""
    use_bioformats_first = force_bioformats or all_series or not is_native_image_file(path)
    entries = None
    if use_bioformats_first:
        try:
            entries = plan_bioformats_series(path, all_series, log, skip_thumbnails)
        except ANY_ERROR as exc:
            if log is not None:
                log("Bio-Formats could not open " + to_text(os.path.basename(path)) + ": " + to_text(exc))
            entries = None
    if entries is None or len(entries) == 0:
        imp = None
        # Send ImageJ's "cannot open" message to the Log window instead of a
        # dialog, so one unreadable file does not hold up an unattended run.
        log_state = quiet_bioformats_log()
        try:
            IJ.redirectErrorMessages()
            imp = IJ.openImage(path)
        except ANY_ERROR:
            imp = None
        restore_bioformats_log(log_state)
        try:
            IJ.redirectErrorMessages(False)
        except ANY_ERROR:
            pass
        if imp is None:
            # ImageJ keeps a file it could not read open until Java next
            # collects garbage, which on Windows stops the user from moving
            # or deleting the damaged file. This is rare, so collect now.
            try:
                from java.lang import System as JavaSystem
                JavaSystem.gc()
            except ANY_ERROR:
                pass
        if imp is None and not use_bioformats_first:
            try:
                entries = plan_bioformats_series(path, all_series, log, skip_thumbnails)
            except ANY_ERROR as exc:
                if log is not None:
                    log("Bio-Formats could not open " + to_text(os.path.basename(path)) + ": " + to_text(exc))
                entries = None
        elif imp is not None:
            entries = [BatchSeries("", "", image=imp)]
    if entries is None:
        return []
    return entries


def open_with_bioformats_api(path, all_series=False, log=None):
    """Open a file through the Bio-Formats API without showing a dialog.
    Returns a list of images (possibly empty) or None if Bio-Formats is missing."""
    entries = plan_bioformats_series(path, all_series, log, False)
    if entries is None:
        return None
    images = [entry.open() for entry in entries]
    return [imp for imp in images if imp is not None]


def open_image_series(path, force_bioformats=False, all_series=False, log=None):
    """Open one file for unattended processing.
    Returns a list of (series_label, ImagePlus); the label is empty for a
    single image. Label and overview series are left out. Every image is in
    memory at once; a batch uses plan_image_series() to open them in turn."""
    result = []
    for entry in plan_image_series(path, force_bioformats, all_series, log):
        if entry.skip_reason != "":
            continue
        try:
            imp = entry.open()
        except ANY_ERROR as exc:
            if log is not None:
                log("Could not open " + (entry.label or to_text(os.path.basename(path))) + ": " + to_text(exc))
            imp = None
        if imp is not None:
            result.append((entry.label, imp))
    return result


def batch_parent_rois(imp, image_path, policy, roi_dir, manager_rois, series_label=""):
    """Resolve the parent ROIs for one batch image.
    Returns (rois, description). rois is None when the image must be skipped."""
    if policy == BATCH_ROI_WHOLE:
        return ([whole_image_roi(imp)], "Whole image")
    if policy == BATCH_ROI_MANAGER:
        if manager_rois is None or len(manager_rois) == 0:
            raise ValueError("The ROI Manager set is empty.")
        return ([clone_roi(roi) for roi in manager_rois], "ROI Manager set")
    if policy in [BATCH_ROI_FILE_SKIP, BATCH_ROI_FILE_WHOLE]:
        roi_path = find_roi_file(image_path, roi_dir, series_label)
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
    "Summary_Rows", "Object_Rows", "Seconds",
    "Series_Name", "Image_Width", "Image_Height", "Image_Channels", "Image_Slices", "Image_Frames",
    "Pixel_Width", "Pixel_Height", "Spatial_Unit"
]


def run_batch(input_dir, output_root, configs, all_configs, global_options, batch_options,
              manager_rois=None, log=None, progress=None, is_cancelled=None, file_list=None):
    """Analyze every image of a folder with one fixed configuration.

    batch_options: extensions (text), recursive, roi_policy, roi_dir,
                   force_bioformats, all_series, skip_thumbnail_series
                   (default true: with all_series, label and overview
                   series are reported as skipped instead of analyzed).
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
    if global_options.get("frame_mode") is not None:
        check_frame_mode(global_options.get("frame_mode"))
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
    skip_thumbnails = parse_bool(batch_options.get("skip_thumbnail_series"), True)
    first_channels = [None]

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
    summary_out = CsvAppender(os.path.join(batch_dir, "batch_roi_summary.csv"), header_engine.build_summary_headers(configs, global_options) + extra)
    object_out = CsvAppender(os.path.join(batch_dir, "batch_objects.csv"), header_engine.build_object_headers(configs, global_options) + extra)
    extra_out = {}
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
            open_notes = []
            file_start = time.time()

            def open_log(message, _notes=open_notes):
                _notes.append(to_text(message).strip())
                batch_log(message)

            try:
                series = plan_image_series(path, force_bioformats, all_series, open_log, skip_thumbnails)
            except ANY_ERROR as exc:
                open_error = to_text(exc)
            if len(series) == 0:
                counts["failed"] += 1
                if open_error == "" and len(open_notes) > 0:
                    open_error = open_notes[-1]
                message = "Could not open the file" + (": " + open_error if open_error else ".")
                batch_log("  FAILED - " + message)
                report_out.append([{
                    "Batch_ID": batch_id, "Image_Index": file_index + 1,
                    "Image_File": os.path.basename(path), "Image_Relative_Path": relative,
                    "Status": "FAILED", "Message": message,
                    "Seconds": round(time.time() - file_start, 2)
                }])
                continue

            for entry in series:
                series_label = entry.label
                report = {
                    "Batch_ID": batch_id, "Image_Index": file_index + 1,
                    "Image_File": os.path.basename(path), "Image_Relative_Path": relative,
                    "Series": series_label, "Series_Name": entry.name
                }
                if entry.skip_reason != "":
                    counts["skipped"] += 1
                    report["Status"] = "SKIPPED"
                    report["Message"] = entry.skip_reason
                    batch_log("  SKIPPED - " + series_label + ": " + entry.skip_reason)
                    report_out.append([report])
                    continue
                # The time of an image includes opening it. A single image
                # was opened while the file was examined.
                image_start = file_start if entry.image is not None else time.time()
                image_log = []
                imp = None

                def image_logger(message, _lines=image_log):
                    _lines.append("[" + time.strftime("%H:%M:%S") + "] " + to_text(message))

                try:
                    if is_cancelled is not None and is_cancelled():
                        raise AnalysisCancelled("Cancelled by user.")
                    imp = entry.open()
                    if imp is None:
                        raise ValueError("Bio-Formats returned no image for " + (series_label or "this file") + ".")
                    # Name the image after its file. Bio-Formats titles an
                    # image with its whole path and replaces the middle of a
                    # long one with "...", which would make Source_Image and
                    # the run folder name ambiguous. A series keeps the file
                    # extension after its label, so that a name with dots in
                    # it is not cut short when the extension is removed.
                    file_name = to_text(os.path.basename(path))
                    if series_label != "":
                        stem = strip_extension(file_name)
                        imp.setTitle(stem + "_" + series_label + file_name[len(stem):])
                    else:
                        imp.setTitle(file_name)
                    calibration = imp.getCalibration()
                    report["Image_Width"] = imp.getWidth()
                    report["Image_Height"] = imp.getHeight()
                    report["Image_Channels"] = imp.getNChannels()
                    report["Image_Slices"] = imp.getNSlices()
                    report["Image_Frames"] = imp.getNFrames()
                    report["Pixel_Width"] = calibration.pixelWidth
                    report["Pixel_Height"] = calibration.pixelHeight
                    report["Spatial_Unit"] = to_text(calibration.getUnit())
                    image_note = ""
                    if first_channels[0] is None:
                        first_channels[0] = imp.getNChannels()
                    elif imp.getNChannels() != first_channels[0]:
                        image_note = (
                            "Note: this image has " + str(imp.getNChannels()) + " channel(s); the first image of the batch had " +
                            str(first_channels[0]) + ". Check that the channel numbers of the settings mean the same stain here."
                        )
                        batch_log("  " + image_note)
                    validate_channel_configs(imp, configs)
                    if policy in BATCH_ROI_PROVIDERS:
                        rois, roi_note = BATCH_ROI_PROVIDERS[policy]({
                            "imp": imp, "path": path, "configs": configs,
                            "global_options": global_options, "log": batch_log,
                            "is_cancelled": is_cancelled
                        })
                    else:
                        rois, roi_note = batch_parent_rois(imp, path, policy, roi_dir, manager_rois, series_label)
                    report["ROI_Source"] = roi_note
                    if rois is None:
                        counts["skipped"] += 1
                        report["Status"] = "SKIPPED"
                        report["Message"] = roi_note
                        batch_log("  SKIPPED - " + roi_note)
                    else:
                        rois = validate_parent_rois(rois, imp.getWidth(), imp.getHeight())
                        frames = range(1, imp.getNFrames() + 1) if all_frames else [1]
                        plan = plan_run_folder(
                            batch_dir, strip_extension(imp.getTitle()), "", configs,
                            multi_frame=imp.getNFrames() > 1, options=global_options
                        )
                        run_dir = unique_dir(batch_dir, plan["name"])
                        if plan["shortened"]:
                            batch_log("  " + shortened_name_note(plan["base"], plan["file_base"]))
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
                        # Tables of optional analysis modules, combined as batch_<name>.
                        for table_name in result.get("extra_table_names", []):
                            table = result["extra_tables"][table_name]
                            for row in table["rows"]:
                                row["Batch_ID"] = batch_id
                                row["Image_Relative_Path"] = relative
                            if table_name not in extra_out:
                                extra_out[table_name] = CsvAppender(
                                    os.path.join(batch_dir, "batch_" + table_name), table["headers"] + extra
                                )
                            extra_out[table_name].append(table["rows"])
                        counts["ok"] += 1
                        report["Status"] = "OK"
                        report["Message"] = image_note
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
                    # Close the image and give its pixels back before the
                    # next one is opened. The table rows of this image are
                    # already written; let them go too (about 0.6 GB for
                    # 65000 objects).
                    entry.release()
                    imp = None
                    result = None
                if cancelled:
                    break
            # Close whatever of this file was opened but not reached.
            for entry in series:
                entry.release()
            if cancelled:
                break
    finally:
        summary_out.close()
        object_out.close()
        threshold_out.close()
        report_out.close()
        for table_name in extra_out:
            extra_out[table_name].close()
        duration = time.time() - batch_start
        if cancelled:
            batch_log("Batch cancelled by user.")
        batch_log(
            "Batch finished in " + IJ.d2s(duration, 1) + " s: " + str(counts["ok"]) + " analyzed, " +
            str(counts["skipped"]) + " skipped, " + str(counts["failed"]) + " failed."
        )
        longest_path = longest_path_length(batch_dir)
        path_note = long_path_note(longest_path)
        if path_note != "":
            batch_log(path_note)
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
        "duration": duration,
        "longest_path": longest_path
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
        save_tiff(flattened, os.path.join(run_dir, run_folder_file_name(run_dir, source_title, config, base + "_Validation_Overlay.tif")))
        safe_close_image(flattened)
        safe_close_image(overlay_source)
        save_tiff(mask_data["final"], os.path.join(
            run_dir, run_folder_file_name(run_dir, source_title, config, mask_data["final"].getTitle()) + ".tif"
        ))
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
    log_run_folder_notes(engine, run_dir, source_title, config)
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
        save_tiff(mask_data["raw"], os.path.join(
            run_dir, run_folder_file_name(run_dir, source_title, config, mask_data["raw"].getTitle()) + ".tif"
        ))
    finally:
        safe_close_image(mask_data["raw"])
        safe_close_image(mask_data["final"])
    write_csv(os.path.join(run_dir, "mask_validation_summary.csv"), MASK_VALIDATION_HEADERS, rows)
    write_csv(
        os.path.join(run_dir, "thresholds.csv"), engine.build_threshold_headers(),
        [engine.threshold_row(run_id, source_title, frame, config, mask_data)]
    )
    log_run_folder_notes(engine, run_dir, source_title, config)
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


SENSITIVITY_EMPTY_NOTE = "No pixel value lies in the threshold range: nothing is positive, no objects."


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
        except ANY_ERROR as exc:
            # ANY_ERROR: an error raised inside ImageJ or a plugin is a Java
            # error, and must be noted for this factor like any other.
            if factor == 1.0 or isinstance(exc, AnalysisCancelled):
                raise
            notes[factor] = to_text(exc)
            engine.log("Sensitivity factor " + str(factor) + " skipped: " + to_text(exc))
            continue
        if mask_data.get("empty_range", False):
            notes[factor] = SENSITIVITY_EMPTY_NOTE
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
    log_run_folder_notes(engine, run_dir, source_title, config)
    return {"headers": SENSITIVITY_HEADERS, "rows": rows, "run_dir": run_dir}

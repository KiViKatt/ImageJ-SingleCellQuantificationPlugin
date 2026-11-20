# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach
# Copyright (C) 2026 Maksym Taran

# Cell Quant Studio Universal - threshold spectrum of marker coverage (optional module)
#
# The positive-area fraction of a marker depends on one threshold chosen by
# the analyst. This module shows how much: for every parent ROI and every
# channel analysed for positive area it exports the coverage (percent of the
# ROI pixels that would be positive) at every level of a threshold axis,
# plus a few numbers that say how robust the reported coverage is.
# Cumulative threshold spectra are a known idea from the immunolabelling
# literature; this is an implementation of it, not a new measure.
#
# Definitions (n = pixels of the ROI that hold a number, all counts are exact
# pixel counts taken from the full histogram of the prepared image inside the
# ROI; only a 32-bit image can hold pixels without a number, NaN, and ImageJ
# leaves those out of every measurement, the ROI area included, so they are
# left out of n too and the note says how many there were):
#   coverage(t)   100 * (pixels a threshold t would select) / n. A pixel is
#                 selected at or above t; for an Auto threshold on a light
#                 background at or below t, as build_masks does. The other
#                 threshold of the run (upper or lower limit) stays as it was
#                 used, so the curve passes exactly through the run's own
#                 value. A level is applied the way ImageJ applies a typed
#                 threshold: rounded to the nearest whole value on a 16-bit
#                 image, cut to a whole value on an 8-bit image, and rounded
#                 to single precision on a 32-bit image.
#   Coverage at the configured threshold
#                 the plain threshold coverage, before any binary cleanup
#                 (open, close, erode, dilate, fill holes). Without cleanup it
#                 equals the positive-area fraction of roi_summary.csv.
#   Sensitivity   average slope of coverage(t) over the configured threshold
#                 plus and minus 1 percent (at least one intensity level on 8-
#                 and 16-bit images), given as percentage points of coverage
#                 per +1 percent of threshold (slope * threshold / 100), and
#                 as elasticity (percent change of coverage per +1 percent of
#                 threshold, that is 100 * sensitivity / coverage). The slope
#                 is the change of coverage between the two ends of the
#                 interval divided by the distance between the two ends as
#                 ImageJ applies them, so it is the exact slope of the curve
#                 between two of its points.
#   Stable range  the thresholds, around the configured one, for which the
#                 coverage stays within the tolerance (default 10 percent,
#                 relative) of the reported coverage. It always holds the
#                 configured threshold. When no threshold on the side of more
#                 coverage can leave the tolerance, the range has no end on
#                 that side: it is then cut at the last pixel value of the
#                 ROI (or at the configured threshold when that lies beyond
#                 it) and the note says that the range is open. On 8- and
#                 16-bit images both ends are whole values inside the range;
#                 on a 32-bit image a closed end on the side of more coverage
#                 is the pixel value just outside the range (thresholds above
#                 it, or below it for at-or-below counting, are inside).
#   Normalized mean
#                 the area under the at-or-above curve divided by the axis
#                 width. This is mathematically the mean intensity (clipped
#                 to the axis) rescaled to 0-100 percent of the axis, so it
#                 says nothing that mean intensity does not say. It is named
#                 for what it is.
#
# This file is loaded by load_engine_modules() into the namespace of
# cellquantengine.py: the engine's names are globals here. It has no user
# interface, opens no window, keeps no state between runs, and changes no
# Fiji option.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij.measure import Measurements
from ij.process import ImageStatistics
from java.util import Arrays
from jarray import array as jarray_array

import math


SPECTRUM_TABLE_NAME = "threshold_spectrum.csv"
SPECTRUM_MAX_STEPS = 4096
# Half-width of the interval used for the sensitivity, as a share of the
# configured threshold.
SPECTRUM_SLOPE_WINDOW = 0.01
SPECTRUM_ABOVE = "above"
SPECTRUM_BELOW = "below"

# (property key, option key, kind, default)
SPECTRUM_SETTING_SPECS = [
    ("spectrum.enabled", "spectrum_enabled", "bool", False),
    ("spectrum.steps", "spectrum_steps", "text", "256"),
    ("spectrum.low", "spectrum_low", "text", ""),
    ("spectrum.high", "spectrum_high", "text", ""),
    ("spectrum.tolerance_percent", "spectrum_tolerance_percent", "text", "10")
]

SPECTRUM_HEADERS = [
    "Run_ID", "Application_Version", "Source_Image", "Frame", "Parent_ROI_Index",
    "Parent_ROI_Name", "Channel", "Channel_Name", "Threshold_Mode", "Direction",
    "Level_Index", "Threshold", "Applied_Threshold", "Positive_Pixels", "ROI_Pixels",
    "Coverage_Percent", "Axis_Low", "Axis_High", "Note"
]

# Appended to roi_summary.csv for every positive-area channel, after the
# channel's column prefix.
SPECTRUM_SUMMARY_SUFFIXES = [
    "_Spectrum_Threshold", "_Spectrum_Coverage_Percent",
    "_Spectrum_Sensitivity_PP_Per_Percent", "_Spectrum_Elasticity",
    "_Spectrum_Stable_From", "_Spectrum_Stable_To", "_Spectrum_Stable_Width_Percent",
    "_Spectrum_Normalized_Mean_Percent", "_Spectrum_Note"
]

SPECTRUM_NOTE_LOCAL = (
    "Local Phansalkar has no single threshold: this is the global spectrum of the prepared image "
    "and the values at the configured threshold are not available"
)
SPECTRUM_NOTE_CLEANUP = (
    "binary cleanup is configured: the spectrum is the plain threshold coverage before cleanup "
    "and differs from the positive-area fraction"
)
SPECTRUM_NOTE_AXIS = "the threshold axis is empty for this image (its lowest level is not below its highest)"
SPECTRUM_NOTE_EMPTY_ROI = "the ROI holds no pixels that can be measured"
SPECTRUM_NOTE_CONFIGURED = "configured threshold"
SPECTRUM_NOTE_OPEN_LOW = (
    "the stable range is open at its low end: no lower threshold changes the coverage by more than the tolerance"
)
SPECTRUM_NOTE_OPEN_HIGH = (
    "the stable range is open at its high end: no higher threshold changes the coverage by more than the tolerance"
)
SPECTRUM_NOTE_NAN = (
    " pixel(s) of the ROI hold no number (NaN) and are left out, as ImageJ leaves them out of the ROI area"
)
SPECTRUM_NOTE_MISMATCH = (
    "CHECK FAILED: the mask of this run holds a different number of positive pixels in this ROI, namely "
)


def spectrum_join(parts):
    """Join the notes that are not empty."""
    kept = []
    for part in parts:
        if part != "":
            kept.append(part)
    return "; ".join(kept)


# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

def spectrum_settings(global_options):
    """Read and check the settings of this module. Raises ValueError with a
    plain message when a value cannot be used."""
    options = global_options if global_options is not None else {}

    def text_of(key, default):
        text = to_text(options.get(key)).strip()
        return default if text == "" else text

    def bound(key, label):
        text = text_of(key, "")
        if text == "":
            return None
        return parse_float_value(text, label)

    settings = {
        "enabled": parse_bool(options.get("spectrum_enabled"), False),
        "steps": parse_int_value(
            text_of("spectrum_steps", "256"), "Threshold spectrum: number of levels", 2, SPECTRUM_MAX_STEPS
        ),
        "low": bound("spectrum_low", "Threshold spectrum: lowest threshold"),
        "high": bound("spectrum_high", "Threshold spectrum: highest threshold"),
        "tolerance": parse_float_value(
            text_of("spectrum_tolerance_percent", "10"), "Threshold spectrum: stable-range tolerance", 0
        )
    }
    if settings["tolerance"] >= 100.0:
        raise ValueError("Threshold spectrum: stable-range tolerance must be below 100 percent.")
    if settings["low"] is not None and settings["high"] is not None and not settings["high"] > settings["low"]:
        raise ValueError("Threshold spectrum: the highest threshold must be above the lowest threshold.")
    return settings


# -----------------------------------------------------------------------------
# Exact intensity distribution of the pixels of one ROI
# -----------------------------------------------------------------------------

def spectrum_cumulate(counts, first, last):
    """(above, total) for the values first..last of a histogram:
    above[v] = pixels with a value from v to last, and total = above[first]."""
    above = [0] * (len(counts) + 1)
    running = 0
    value = last
    while value >= first:
        running += counts[value]
        above[value] = running
        value -= 1
    return (above, running)


class SpectrumCounts(object):
    """Pixels of an 8-bit or 16-bit image inside one ROI, as cumulative
    counts built from the full histogram (256 or 65536 whole values).

    first and last may name the lowest and highest value of the whole image.
    They only save time (a 12-bit image in a 16-bit file uses a sixteenth of
    the values); when they do not hold every pixel they are ignored."""

    whole = True

    def __init__(self, histogram, first=None, last=None):
        counts = list(histogram)
        size = len(counts)
        total = sum(counts)
        low = 0 if first is None else max(0, min(size - 1, int(first)))
        high = size - 1 if last is None else max(low, min(size - 1, int(last)))
        above, running = spectrum_cumulate(counts, low, high)
        if running != total:
            low = 0
            high = size - 1
            above, running = spectrum_cumulate(counts, low, high)
        if low > 0:
            above[0:low] = [running] * low
        self.above = above      # above[v] = pixels with value >= v
        self.size = size
        self.first = low        # no pixel is below first or above last
        self.last = high
        self.n = running        # pixels of the ROI
        self.valid = running    # the same: every pixel holds a number
        self.missing = 0        # pixels without a number (never on 8/16-bit)

    def at_least(self, threshold):
        """Pixels with value >= threshold."""
        if threshold <= 0:
            return self.n
        if threshold > self.size - 1:
            return 0
        return self.above[int(math.ceil(threshold))]

    def greater(self, threshold):
        """Pixels with value > threshold."""
        if threshold < 0:
            return self.n
        if threshold >= self.size - 1:
            return 0
        return self.above[int(math.floor(threshold)) + 1]

    def kth_largest(self, k):
        """The k-th largest pixel value (k from 1 to n)."""
        low = 0
        high = self.size
        while high - low > 1:
            middle = (low + high) // 2
            if self.above[middle] >= k:
                low = middle
            else:
                high = middle
        return low

    def kth_smallest(self, k):
        """The k-th smallest pixel value (k from 1 to n)."""
        low = 0
        high = self.size
        limit = self.n - k
        while high - low > 1:
            middle = (low + high) // 2
            if self.above[middle] <= limit:
                high = middle
            else:
                low = middle
        return high - 1

    def clipped_mean(self, low, high):
        """Mean of the pixel values after limiting each to [low, high].

        Uses the identity mean(clip(X)) - low = integral over [low, high] of
        the share of pixels at or above t, summed exactly over whole values."""
        total = 0.0
        if low < self.first:
            # every pixel is at or above each threshold up to first
            total += self.n * (min(high, float(self.first)) - low)
        value = max(self.first + 1, int(math.floor(low)) + 1)
        last = min(self.last, int(math.ceil(high)))
        while value <= last:
            left = max(low, value - 1.0)
            right = min(high, float(value))
            if right > left:
                total += self.above[value] * (right - left)
            value += 1
        return low + total / float(self.n)


class SpectrumSorted(object):
    """Pixels of a 32-bit image inside one ROI, sorted. processor holds the
    ROI's bounding rectangle with every pixel outside the ROI set to
    not-a-number; sorting puts those after all real values.

    Pixels of the ROI that hold no number are left out of n, exactly as
    ImageJ leaves them out of the ROI area and of every statistic, so that
    the coverage has the same denominator as the positive-area fraction.
    missing says how many there were."""

    whole = False

    def __init__(self, processor, roi_pixels):
        self.processor = processor
        self.pixels = processor.getPixels()
        Arrays.sort(self.pixels)
        low = 0
        high = len(self.pixels)
        while low < high:
            middle = (low + high) // 2
            value = self.pixels[middle]
            if value != value:
                high = middle
            else:
                low = middle + 1
        self.valid = low        # pixels that hold a number
        self.n = low
        self.missing = roi_pixels - low

    def first_index(self, threshold, strict):
        """Index of the first sorted value >= threshold (> when strict)."""
        low = 0
        high = self.valid
        while low < high:
            middle = (low + high) // 2
            value = self.pixels[middle]
            if value < threshold or (strict and value == threshold):
                low = middle + 1
            else:
                high = middle
        return low

    def at_least(self, threshold):
        return self.valid - self.first_index(threshold, False)

    def greater(self, threshold):
        return self.valid - self.first_index(threshold, True)

    def kth_largest(self, k):
        return self.pixels[self.valid - k]

    def kth_smallest(self, k):
        return self.pixels[k - 1]

    def clipped_mean(self, low, high):
        processor = self.processor
        if low > self.pixels[0] or high < self.pixels[self.valid - 1]:
            processor = processor.duplicate()
            processor.min(low)
            processor.max(high)
        # Not-a-number pixels (outside the ROI) are left out of the mean.
        return ImageStatistics.getStatistics(processor, Measurements.MEAN, None).mean


def spectrum_distribution(imp, roi, value_range=None):
    """Exact distribution of the pixel values of a single-plane grayscale
    image inside roi (None means the whole image). The ROI is applied through
    the same mask ImageJ uses when it measures; the image is not changed.
    value_range, the (lowest, highest) value of the whole image, is optional
    and only saves time on 8- and 16-bit images."""
    ip = imp.getProcessor()
    depth = imp.getBitDepth()
    if depth not in [8, 16, 32]:
        raise ValueError("The threshold spectrum needs a grayscale image, not RGB.")
    if roi is None:
        ip.resetRoi()
    else:
        ip.setRoi(clone_roi(roi))
    try:
        if depth != 32:
            if value_range is None:
                return SpectrumCounts(ip.getHistogram())
            return SpectrumCounts(ip.getHistogram(), value_range[0], value_range[1])
        mask = ip.getMask()
        crop = ip.crop()
    finally:
        ip.resetRoi()
    crop.resetRoi()
    roi_pixels = crop.getWidth() * crop.getHeight()
    if mask is not None:
        outside = mask.duplicate()
        roi_pixels -= outside.getHistogram()[0]
        outside.invert()
        crop.setValue(float("nan"))
        crop.fill(outside)
    return SpectrumSorted(crop, roi_pixels)


# -----------------------------------------------------------------------------
# The measurement
# -----------------------------------------------------------------------------

def spectrum_effective(value, bit_depth):
    """The value ImageJ really compares pixels with when a threshold is
    typed: rounded to the nearest whole value on a 16-bit image, cut to a
    whole value on an 8-bit image, and rounded to single precision (the
    precision of the pixels) on a 32-bit image."""
    if value != value or math.isinf(value):
        return value
    if bit_depth == 16:
        return float(math.floor(value + 0.5))
    if bit_depth == 8:
        return float(int(value))
    # A one-element Java float[] rounds to single precision exactly as a
    # cast does (the Float(double) constructor is marked for removal).
    return float(jarray_array([value], "f")[0])


def spectrum_band_count(dist, low, high, bit_depth):
    """Pixels the thresholds [low, high] select, as createThresholdMask does."""
    count = dist.at_least(spectrum_effective(low, bit_depth)) - dist.greater(spectrum_effective(high, bit_depth))
    return max(0, count)


def spectrum_level_count(dist, plan, threshold):
    """Pixels selected when the run's threshold is moved to threshold."""
    if plan["direction"] == SPECTRUM_BELOW:
        return spectrum_band_count(dist, plan["fixed"], threshold, plan["bit_depth"])
    return spectrum_band_count(dist, threshold, plan["fixed"], plan["bit_depth"])


def spectrum_levels(axis_low, axis_high, steps):
    """Evenly spaced thresholds from axis_low to axis_high, both included."""
    if not axis_high > axis_low:
        return [float(axis_low)]
    span = float(axis_high) - float(axis_low)
    levels = []
    for index in range(steps):
        levels.append(axis_low + span * index / float(steps - 1))
    levels[-1] = float(axis_high)
    return levels


def spectrum_plan(imp, config, mask_data, settings):
    """Everything about one channel of one frame that is the same for every
    parent ROI: the threshold axis, the direction, and the thresholds the
    run really used (taken from the build_masks result)."""
    depth = imp.getBitDepth()
    mode = to_text(mask_data.get("threshold_mode", config["threshold_mode"]))
    data_min = float(mask_data["data_min"])
    data_max = float(mask_data["data_max"])
    axis_low = data_min if settings["low"] is None else settings["low"]
    axis_high = data_max if settings["high"] is None else settings["high"]
    notes = []
    plan = {
        "bit_depth": depth, "mode": mode, "direction": SPECTRUM_ABOVE,
        "axis_low": axis_low, "axis_high": axis_high,
        "levels": spectrum_levels(axis_low, axis_high, settings["steps"]),
        "low": None, "high": None, "free": None, "fixed": float("inf"),
        "cleanup": False, "value_range": (data_min, data_max)
    }
    if not axis_high > axis_low:
        notes.append(SPECTRUM_NOTE_AXIS)

    low = mask_data.get("threshold_low")
    high = mask_data.get("threshold_high")
    if mode == "Local Phansalkar" or not is_number(low) or not is_number(high):
        notes.append(SPECTRUM_NOTE_LOCAL)
    else:
        low = float(low)
        high = float(high)
        plan["low"] = low
        plan["high"] = high
        # ImageJ's automatic threshold selects the dark side when the
        # background is light, and swaps sides under an inverting lookup
        # table; build_masks follows it, and so does this.
        inverted = bool(imp.getProcessor().isInvertedLut())
        if mode == "Auto" and bool(config["dark_background"]) == inverted:
            plan["direction"] = SPECTRUM_BELOW
            plan["free"] = high
            plan["fixed"] = low
            if spectrum_effective(low, depth) > data_min:
                notes.append("pixels below the lower threshold " + to_text(low) + " are not counted at any level")
        else:
            plan["free"] = low
            plan["fixed"] = high
            if spectrum_effective(high, depth) < data_max:
                notes.append("pixels above the upper threshold " + to_text(high) + " are not counted at any level")
        cleanup = config["fill_holes"]
        for key in ["open_iterations", "close_iterations", "erode_iterations", "dilate_iterations"]:
            if config[key] > 0:
                cleanup = True
        if cleanup:
            plan["cleanup"] = True
            notes.append(SPECTRUM_NOTE_CLEANUP)
    plan["note"] = "; ".join(notes)
    return plan


def spectrum_slope(dist, plan, threshold):
    """Average slope of the coverage curve (percentage points per intensity
    unit) over threshold plus and minus 1 percent, from exact pixel counts.
    Returns None when it cannot be computed.

    The slope is the change of coverage between the two ends of the interval
    divided by the distance between the two ends as ImageJ applies them
    (whole values on 8- and 16-bit images). That is the exact slope of the
    coverage curve between two of its points: on a curve that is a straight
    line it is the slope of that line, whatever the threshold. Dividing by
    the width as typed instead would be off by up to a third on 8- and
    16-bit images, because rounding changes the width of the interval but
    the pixels counted are those between the applied ends (8-bit threshold
    125: 123.75 and 126.25 are applied as 123 and 126, a width of 3 rather
    than 2.5; 16-bit threshold 125: as 124 and 126, a width of 2)."""
    width = abs(threshold) * SPECTRUM_SLOPE_WINDOW
    if dist.whole:
        width = max(width, 1.0)
    if dist.n == 0 or not width > 0:
        return None
    depth = plan["bit_depth"]
    low = threshold - width
    high = threshold + width
    applied = spectrum_effective(high, depth) - spectrum_effective(low, depth)
    if not applied > 0:
        return None
    change = spectrum_level_count(dist, plan, high) - spectrum_level_count(dist, plan, low)
    return 100.0 * change / float(dist.n) / applied


def spectrum_stable_range(dist, plan, count, tolerance):
    """(from, to, note): the thresholds for which the number of selected
    pixels stays within tolerance percent of count, the number selected at
    the configured threshold. The range always holds the configured
    threshold (as ImageJ applies it).

    The end toward less coverage is always a pixel value of the ROI. Toward
    more coverage the range may have no end, because selecting every pixel
    that can be selected is still within the tolerance; it is then cut at the
    last pixel value of the ROI on that side, or at the configured threshold
    when that lies beyond it, and note says that the range is open. On a
    32-bit image a closed end toward more coverage is a limit value that
    itself is just outside the range. (None, None, "") when count is zero."""
    if count <= 0:
        return (None, None, "")
    fraction = tolerance / 100.0
    most = int(math.floor(count * (1.0 + fraction) + 1e-9))
    least = max(1, int(math.ceil(count * (1.0 - fraction) - 1e-9)))
    step = 1 if dist.whole else 0
    depth = plan["bit_depth"]
    fixed = spectrum_effective(plan["fixed"], depth)
    used = spectrum_effective(plan["free"], depth)
    note = ""
    if plan["direction"] == SPECTRUM_BELOW:
        # selected(t) = (pixels <= t) - skipped, rising with t
        skipped = dist.valid - dist.at_least(fixed)
        start = dist.kth_smallest(least + skipped)
        if most + skipped >= dist.valid:
            end = max(dist.kth_largest(1), used)
            note = SPECTRUM_NOTE_OPEN_HIGH
        else:
            end = dist.kth_smallest(most + skipped + 1) - step
        return (float(start), float(end), note)
    # selected(t) = (pixels >= t) - skipped, falling with t
    skipped = dist.greater(fixed)
    end = dist.kth_largest(least + skipped)
    if most + skipped >= dist.valid:
        start = min(dist.kth_smallest(1), used)
        note = SPECTRUM_NOTE_OPEN_LOW
    else:
        start = dist.kth_largest(most + skipped + 1) + step
    return (float(start), float(end), note)


def spectrum_measure(dist, plan, tolerance):
    """All numbers of one parent ROI and channel. Values that do not exist
    are returned as empty text, which is what the tables show."""
    depth = plan["bit_depth"]
    n = dist.n
    counts = []
    for threshold in plan["levels"]:
        counts.append(spectrum_level_count(dist, plan, threshold))
    result = {
        "roi_pixels": n, "counts": counts, "threshold": "", "count": "", "coverage": "",
        "sensitivity": "", "elasticity": "", "stable_from": "", "stable_to": "",
        "stable_width": "", "normalized_mean": "", "note": plan["note"]
    }
    notes = [plan["note"]]
    if dist.missing > 0:
        # Only a 32-bit image can hold pixels without a number.
        notes.append(str(dist.missing) + SPECTRUM_NOTE_NAN)
    if n == 0:
        notes.append(SPECTRUM_NOTE_EMPTY_ROI)
    result["note"] = spectrum_join(notes)
    if n == 0:
        return result
    if plan["axis_high"] > plan["axis_low"]:
        mean = dist.clipped_mean(plan["axis_low"], plan["axis_high"])
        result["normalized_mean"] = 100.0 * (mean - plan["axis_low"]) / (plan["axis_high"] - plan["axis_low"])
    free = plan["free"]
    if free is None:
        return result

    count = spectrum_band_count(dist, plan["low"], plan["high"], depth)
    coverage = 100.0 * count / float(n)
    result["threshold"] = free
    result["count"] = count
    result["coverage"] = coverage
    slope = spectrum_slope(dist, plan, free)
    if slope is not None and free > 0:
        result["sensitivity"] = slope * free / 100.0
        if coverage > 0:
            result["elasticity"] = slope * free / coverage
    stable = spectrum_stable_range(dist, plan, count, tolerance)
    if stable[0] is not None:
        result["stable_from"] = stable[0]
        result["stable_to"] = stable[1]
        if free > 0:
            result["stable_width"] = 100.0 * (stable[1] - stable[0]) / free
        result["note"] = spectrum_join([result["note"], stable[2]])
    return result


def spectrum_check_mask(engine, mask_data, roi, plan, result):
    """Compare the count at the configured threshold with the mask the run
    really made, whenever the two must agree (a single threshold and no
    binary cleanup). A difference is written into the note instead of being
    passed over: roi_summary.csv then holds the value of the mask."""
    if result["count"] == "" or plan["cleanup"]:
        return
    in_mask = int(engine.measure_positive_area(mask_data["raw"], roi)["positive_pixels"])
    if in_mask != result["count"]:
        result["note"] = spectrum_join([result["note"], SPECTRUM_NOTE_MISMATCH + str(in_mask)])


def spectrum_for_roi(imp, roi, config, mask_data, global_options):
    """Threshold spectrum of one prepared image inside one ROI.

    mask_data is the build_masks result of the same image and channel
    settings. Returns (plan, result): plan["levels"] holds the thresholds and
    result["counts"] the selected pixels at each of them."""
    settings = spectrum_settings(global_options)
    plan = spectrum_plan(imp, config, mask_data, settings)
    return (plan, spectrum_measure(spectrum_distribution(imp, roi), plan, settings["tolerance"]))


# -----------------------------------------------------------------------------
# Tables
# -----------------------------------------------------------------------------

def spectrum_table_rows(base, plan, result):
    """Rows of threshold_spectrum.csv for one frame, parent ROI, and channel:
    level 0 is the threshold the run used (when there is one), levels 1 and
    up are the axis."""
    depth = plan["bit_depth"]
    n = result["roi_pixels"]
    rows = []
    if result["count"] != "":
        row = dict(base)
        row["Level_Index"] = 0
        row["Threshold"] = plan["free"]
        row["Applied_Threshold"] = spectrum_effective(plan["free"], depth)
        row["Positive_Pixels"] = result["count"]
        row["Coverage_Percent"] = result["coverage"]
        row["Note"] = spectrum_join([SPECTRUM_NOTE_CONFIGURED, result["note"]])
        rows.append(row)
    for index in range(len(plan["levels"])):
        row = dict(base)
        row["Level_Index"] = index + 1
        row["Threshold"] = plan["levels"][index]
        row["Applied_Threshold"] = spectrum_effective(plan["levels"][index], depth)
        row["Positive_Pixels"] = result["counts"][index]
        row["Coverage_Percent"] = 100.0 * result["counts"][index] / float(n) if n > 0 else ""
        row["Note"] = result["note"]
        rows.append(row)
    return rows


def spectrum_fill_summary(summary, prefix, result):
    summary[prefix + "_Spectrum_Threshold"] = result["threshold"]
    summary[prefix + "_Spectrum_Coverage_Percent"] = result["coverage"]
    summary[prefix + "_Spectrum_Sensitivity_PP_Per_Percent"] = result["sensitivity"]
    summary[prefix + "_Spectrum_Elasticity"] = result["elasticity"]
    summary[prefix + "_Spectrum_Stable_From"] = result["stable_from"]
    summary[prefix + "_Spectrum_Stable_To"] = result["stable_to"]
    summary[prefix + "_Spectrum_Stable_Width_Percent"] = result["stable_width"]
    summary[prefix + "_Spectrum_Normalized_Mean_Percent"] = result["normalized_mean"]
    summary[prefix + "_Spectrum_Note"] = result["note"]


class ThresholdSpectrumExtension(RunExtension):
    """Adds threshold_spectrum.csv and the _Spectrum_ columns of
    roi_summary.csv to a run when spectrum.enabled is set."""

    name = "threshold_spectrum"

    def is_enabled(self, global_options, configs):
        if not parse_bool(global_options.get("spectrum_enabled"), False):
            return False
        for config in configs:
            if needs_positive_area(config):
                return True
        return False

    def summary_headers(self, configs, global_options):
        headers = []
        for config in configs:
            if needs_positive_area(config):
                prefix = channel_column_prefix(config)
                for suffix in SPECTRUM_SUMMARY_SUFFIXES:
                    headers.append(prefix + suffix)
        return headers

    def begin_run(self, context):
        # Checked here so that a setting that cannot be used stops the run
        # with a plain message before any work is done.
        spectrum_settings(context["global_options"])
        context["state"][self.name]["rows"] = []

    def frame_done(self, context):
        settings = spectrum_settings(context["global_options"])
        rows = context["state"][self.name]["rows"]
        engine = context["engine"]
        for config in context["configs"]:
            if not needs_positive_area(config):
                continue
            channel = config["channel"]
            imp = context["prepared"][channel]
            mask_data = context["masks"][channel]
            plan = spectrum_plan(imp, config, mask_data, settings)
            prefix = channel_column_prefix(config)
            engine.log(
                "Threshold spectrum of C" + str(channel) + ": " + str(len(plan["levels"])) +
                " level(s) from " + to_text(plan["axis_low"]) + " to " + to_text(plan["axis_high"]) + "."
            )
            for roi_index in range(len(context["parent_rois"])):
                engine.check_cancelled()
                summary = context["summary_rows"][roi_index]
                roi = context["parent_rois"][roi_index]
                dist = spectrum_distribution(imp, roi, plan["value_range"])
                result = spectrum_measure(dist, plan, settings["tolerance"])
                spectrum_check_mask(engine, mask_data, roi, plan, result)
                spectrum_fill_summary(summary, prefix, result)
                base = {
                    "Run_ID": context["run_id"],
                    "Application_Version": APP_VERSION,
                    "Source_Image": context["source_title"],
                    "Frame": context["frame"],
                    "Parent_ROI_Index": roi_index + 1,
                    "Parent_ROI_Name": summary.get("Parent_ROI_Name", ""),
                    "Channel": channel,
                    "Channel_Name": config["name"],
                    "Threshold_Mode": plan["mode"],
                    "Direction": plan["direction"],
                    "ROI_Pixels": result["roi_pixels"],
                    "Axis_Low": plan["axis_low"],
                    "Axis_High": plan["axis_high"]
                }
                rows.extend(spectrum_table_rows(base, plan, result))

    def end_run(self, context):
        return [(SPECTRUM_TABLE_NAME, SPECTRUM_HEADERS, context["state"][self.name].get("rows", []))]


register_global_settings(SPECTRUM_SETTING_SPECS)
register_settings_form(
    "Threshold spectrum",
    "For every positive-area channel and parent ROI, exports the coverage at every level of a threshold "
    "axis (threshold_spectrum.csv) and how robust the reported coverage is to the threshold. Leave both "
    "axis limits empty to use each image's own minimum and maximum; type both to put a whole batch on "
    "one axis.",
    [
        ("spectrum_enabled", "Export the threshold spectrum", "bool", None),
        ("spectrum_steps", "Number of threshold levels", "text", None),
        ("spectrum_low", "Lowest threshold of the axis (empty = image minimum)", "text", None),
        ("spectrum_high", "Highest threshold of the axis (empty = image maximum)", "text", None),
        ("spectrum_tolerance_percent", "Stable-range tolerance (percent of the reported coverage)", "text", None)
    ]
)
register_run_extension(ThresholdSpectrumExtension())
register_engine_module("cellquantspectrum.py", "2.3.0")

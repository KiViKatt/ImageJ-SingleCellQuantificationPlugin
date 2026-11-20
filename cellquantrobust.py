# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach

# Cell Quant Studio Universal - robustness of a group difference to the
# analysis settings (optional module)
#
# A difference between two groups of images (for example control and mutant
# organoids) in an object density or a marker-positive area fraction should
# not depend on one arbitrary threshold. This module repeats the measurement
# of one channel over a grid of defensible settings ("specifications") and
# reports, for every specification, the group means and the size of the
# difference, and over all of them in how many the direction of the
# difference is the same as with the configured settings.
#
# This is an application of multiverse / specification-curve analysis
# (Steegen et al. 2016; Simonsohn et al. 2020) to segmentation settings. It
# is a description of how much a result depends on the settings, not a
# statistical test: the specifications reuse the same images and are
# strongly correlated.
#
# How it works:
#   1. Every image is opened once. Only the chosen channel is prepared and
#      segmented again for each specification, with the same engine
#      functions a normal run uses, inside the same parent ROIs a batch
#      would use.
#   2. Each image belongs to a group and to an experimental unit (for
#      example the organoid). Measurements are first combined within each
#      unit; group statistics are computed across units only, because
#      several images or ROIs of one organoid are not independent.
#   3. Per specification: number of units, mean and SD per group, difference
#      and ratio of the means, Hedges g, a Welch test with a 95 percent
#      interval of the difference, and a Mann-Whitney test.
#   4. Over all specifications: in how many the difference has the direction
#      it has with the configured settings, the range of the effect size,
#      and which values of which setting go with a reversed direction.
#
# This file is loaded by load_engine_modules() into the namespace of
# cellquantengine.py: the engine's names are globals here. It has no user
# interface, opens no window, changes no Fiji option, and adds nothing to
# the tables of a normal run or batch.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij import IJ, ImagePlus, ImageStack
from ij.gui import Roi
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
from ij.process import ImageProcessor, ByteProcessor, Blitter

import os
import re
import time
import math
import codecs
import itertools
import traceback

# The Welch test and its interval need Student's t distribution, taken from
# the Apache Commons Math library that ships with Fiji. Without it the Welch
# p value and the interval are left empty; everything else is computed here.
ROBUST_MATH = {"available": False}
try:
    from org.apache.commons.math3.stat.inference import TTest as RobustTTest
    from org.apache.commons.math3.distribution import TDistribution as RobustTDistribution
    from jarray import array as robust_java_array
    ROBUST_MATH["available"] = True
except ANY_ERROR:
    pass


ROBUST_METRIC_DENSITY = "Object count density"
ROBUST_METRIC_AREA_FRACTION = "Positive-area fraction (percent)"
ROBUST_METRIC_OBJECT_AREA = "Mean object area"
ROBUST_METRICS = [ROBUST_METRIC_DENSITY, ROBUST_METRIC_AREA_FRACTION, ROBUST_METRIC_OBJECT_AREA]
# Asked for by name so that the refusal can explain itself.
ROBUST_METRIC_INTENSITY = "Mean intensity"
# value = numerator / denominator * scale (see robust_metric_value)
ROBUST_METRIC_SCALE = {
    ROBUST_METRIC_DENSITY: 1.0, ROBUST_METRIC_AREA_FRACTION: 100.0, ROBUST_METRIC_OBJECT_AREA: 1.0
}

ROBUST_UNIT_MEAN = "Mean of the measurements of the unit"
ROBUST_UNIT_POOLED = "Pooled over the unit (sum / sum)"
ROBUST_UNIT_AGGREGATIONS = [ROBUST_UNIT_MEAN, ROBUST_UNIT_POOLED]

ROBUST_GROUP_PREFIX = "File name before the first separator"
ROBUST_GROUP_REGEX = "Regular expression on the relative path"
ROBUST_GROUP_CSV = "Table (CSV file)"
ROBUST_GROUP_SOURCES = [ROBUST_GROUP_PREFIX, ROBUST_GROUP_REGEX, ROBUST_GROUP_CSV]

ROBUST_DEFAULT_MAX_SPECIFICATIONS = 60
ROBUST_RAW_GRID_LIMIT = 5000     # combinations examined before duplicates are removed
ROBUST_MIN_UNITS = 3             # fewer units in a group: means only
ROBUST_EXACT_LIMIT = 20000       # arrangements enumerated for an exact Mann-Whitney p
ROBUST_ALPHA = 0.05
ROBUST_DEFAULT_GRID = "threshold_factor=0.8,0.9,1.0,1.1,1.2"

# Grid vocabulary. threshold_factor multiplies the threshold level of the
# configured setting (as the sensitivity sweep does); every other key is a
# channel setting and takes the values that setting takes in a profile.
ROBUST_FACTOR_KEY = "threshold_factor"
ROBUST_PREPARATION_KEYS = ["subtract_background", "background_radius", "median_radius", "gaussian_sigma"]
ROBUST_OBJECT_KEYS = [
    "watershed", "min_size", "max_size", "min_circularity", "max_circularity",
    "exclude_edges", "include_holes", "four_connected"
]
ROBUST_CLAHE_KEYS = ["clahe_block_size", "clahe_histogram_bins", "clahe_max_slope"]
ROBUST_LOCAL_KEYS = ["local_radius", "local_parameter_1", "local_parameter_2", "clahe"] + ROBUST_CLAHE_KEYS
ROBUST_MORPHOLOGY_KEYS = ["open_iterations", "close_iterations", "erode_iterations", "dilate_iterations"]
ROBUST_GRID_KEYS = (
    [ROBUST_FACTOR_KEY] + ROBUST_PREPARATION_KEYS +
    ["threshold_mode", "auto_method", "auto_scale", "manual_low", "manual_high"] +
    ROBUST_LOCAL_KEYS + ROBUST_MORPHOLOGY_KEYS + ["morph_count", "fill_holes"] + ROBUST_OBJECT_KEYS
)

# Earlier engines reported a threshold range that holds no pixel value with
# this error text; the engine now returns an empty mask with
# mask_data["empty_range"] True instead. Kept for scripts that refer to it.
ROBUST_EMPTY_RANGE_TEXT = "Threshold upper value is below the lower value"

ROBUST_STATISTICS_NOTE = (
    "Descriptive only. The specifications reuse the same images and are strongly correlated, so the "
    "fractions below are not significance tests and are not corrected for multiple comparisons."
)

# (property key, option key, kind, default)
ROBUST_SETTING_SPECS = [
    ("robust.channel", "robust_channel", "text", "1"),
    ("robust.metric", "robust_metric", "text", ROBUST_METRIC_DENSITY),
    ("robust.grid", "robust_grid", "text", ROBUST_DEFAULT_GRID),
    ("robust.max_specifications", "robust_max_specifications", "text", str(ROBUST_DEFAULT_MAX_SPECIFICATIONS)),
    ("robust.group_source", "robust_group_source", "text", ROBUST_GROUP_PREFIX),
    ("robust.separator", "robust_separator", "text", "_"),
    ("robust.group_regex", "robust_group_regex", "text", ""),
    ("robust.unit_regex", "robust_unit_regex", "text", ""),
    ("robust.group_csv", "robust_group_csv", "text", ""),
    ("robust.control_group", "robust_control_group", "text", ""),
    ("robust.unit_aggregation", "robust_unit_aggregation", "text", ROBUST_UNIT_MEAN)
]

ROBUST_ID_HEADERS = ["Run_ID", "Application_Version"]
ROBUST_SPEC_ID_HEADERS = ROBUST_ID_HEADERS + ["Specification_ID", "Specification", "Is_Reference"]
ROBUST_METRIC_HEADERS = ["Channel", "Channel_Name", "Metric", "Metric_Unit"]

ROBUST_MEASUREMENT_HEADERS = ROBUST_SPEC_ID_HEADERS + [
    "Image_File", "Image_Relative_Path", "Series", "Frame", "Group", "Unit",
    "Parent_ROI_Index", "Parent_ROI_Name"
] + ROBUST_METRIC_HEADERS + [
    "Value", "ROI_Area", "ROI_Area_Unit", "Object_Count", "Total_Object_Area", "Positive_Area",
    "Threshold_Low", "Threshold_High", "Status", "Note"
]
ROBUST_UNIT_HEADERS = ROBUST_SPEC_ID_HEADERS + [
    "Group", "Unit", "N_Images", "N_Measurements", "N_Measurements_Used", "Unit_Aggregation"
] + ROBUST_METRIC_HEADERS + ["Value"]
ROBUST_GROUP_HEADERS = ROBUST_SPEC_ID_HEADERS + ["Group", "N_Units"] + ROBUST_METRIC_HEADERS + [
    "Mean", "SD", "Minimum", "Maximum"
]
ROBUST_SPECIFICATION_HEADERS = ROBUST_ID_HEADERS + [
    "Curve_Rank", "Specification_ID", "Specification", "Is_Reference"
] + ROBUST_METRIC_HEADERS + [
    "N_Groups", "Group_A", "Group_B", "N_Units_A", "Mean_A", "SD_A", "N_Units_B", "Mean_B", "SD_B",
    "Difference_B_Minus_A", "Ratio_B_Over_A", "Hedges_g", "Welch_t", "Welch_df", "Welch_P",
    "Difference_CI95_Low", "Difference_CI95_High", "Mann_Whitney_U", "Mann_Whitney_P",
    "Mann_Whitney_Method", "Sign", "Same_Sign_As_Reference", "N_Measurements",
    "N_Failed_Measurements", "Complete", "Note"
]
ROBUST_PARAMETER_HEADERS = ROBUST_SPEC_ID_HEADERS + ["Parameter", "Value", "Reference_Value", "Applies"]
ROBUST_VALUE_HEADERS = ROBUST_ID_HEADERS + [
    "Parameter", "Value", "Is_Reference_Value", "N_Specifications", "N_Same_Sign_As_Reference",
    "N_Opposite_Sign", "N_Zero_Or_Undefined", "Fraction_Opposite_Sign", "More_Often_Reversed_Than_Overall"
]
ROBUST_SUMMARY_HEADERS = ROBUST_ID_HEADERS + ROBUST_METRIC_HEADERS + [
    "Group_A", "Group_B", "N_Groups", "Unit_Aggregation", "Group_Rule", "Grid",
    "N_Images_Analyzed", "N_Images_Failed", "N_Images_Skipped", "N_Images_Not_Assigned",
    "N_Specifications", "N_Specifications_Incomplete", "N_Specifications_Compared",
    "Reference_Specification_ID", "Reference_N_Units_A", "Reference_N_Units_B",
    "Reference_Difference_B_Minus_A", "Reference_Hedges_g", "Reference_Sign",
    "N_Same_Sign_As_Reference", "N_Opposite_Sign", "N_Zero_Or_Undefined",
    "Fraction_Same_Sign_As_Reference",
    "Hedges_g_Minimum", "Hedges_g_Median", "Hedges_g_Maximum",
    "Difference_Minimum", "Difference_Median", "Difference_Maximum",
    "Ratio_Minimum", "Ratio_Median", "Ratio_Maximum",
    "Descriptive_N_Welch_P_Below_0.05_Same_Sign", "Descriptive_Fraction_Welch_P_Below_0.05_Same_Sign",
    "Descriptive_N_Mann_Whitney_P_Below_0.05_Same_Sign",
    "Descriptive_Fraction_Mann_Whitney_P_Below_0.05_Same_Sign",
    "Sign_Flipping_Values", "Curve_Sorted_By", "Statistics_Note", "Note"
]
ROBUST_ASSIGNMENT_HEADERS = ROBUST_ID_HEADERS + [
    "Image_File", "Image_Relative_Path", "Group", "Unit", "Status", "Message", "Rule"
]
ROBUST_REPORT_HEADERS = ROBUST_ID_HEADERS + [
    "Image_Index", "Image_File", "Image_Relative_Path", "Series", "Group", "Unit", "Status",
    "Message", "ROI_Source", "Parent_ROI_Count", "Frames", "Measurements", "Seconds"
]


# -----------------------------------------------------------------------------
# Small numeric helpers
# -----------------------------------------------------------------------------

def robust_mean(values):
    if len(values) == 0:
        return None
    return sum(values) / float(len(values))


def robust_variance(values):
    """Sample variance (divisor n - 1); None with fewer than two values."""
    if len(values) < 2:
        return None
    mean = robust_mean(values)
    return sum([(value - mean) * (value - mean) for value in values]) / float(len(values) - 1)


def robust_sd(values):
    variance = robust_variance(values)
    return None if variance is None else math.sqrt(variance)


def robust_median(values):
    if len(values) == 0:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2 == 1:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0


def robust_cell(value):
    """A table cell: empty for a value that could not be computed."""
    return "" if value is None else value


def robust_value_text(value):
    if isinstance(value, bool):
        return bool_text(value)
    if isinstance(value, float):
        return to_text(repr(value))
    return to_text(value)


def robust_sign(difference, scale):
    """1, -1, or 0. A difference far below rounding error of the means
    counts as zero."""
    if difference is None:
        return None
    if abs(difference) <= 1e-12 * abs(scale):
        return 0
    return 1 if difference > 0 else -1


def robust_sign_text(sign):
    if sign is None:
        return ""
    return {1: "positive", -1: "negative", 0: "zero"}[sign]


# -----------------------------------------------------------------------------
# Statistics of two groups of unit values
# -----------------------------------------------------------------------------

def robust_midranks(values):
    """Ranks 1..n of the values; tied values share the mean of their ranks."""
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start
        while stop + 1 < len(order) and values[order[stop + 1]] == values[order[start]]:
            stop += 1
        shared = (start + stop) / 2.0 + 1.0
        for position in range(start, stop + 1):
            ranks[order[position]] = shared
        start = stop + 1
    return ranks


def robust_combinations(total, chosen):
    count = 1
    for index in range(1, chosen + 1):
        count = count * (total - chosen + index) // index
    return count


def robust_mann_whitney(a, b):
    """Two-sided Mann-Whitney (Wilcoxon rank-sum) test of two samples.

    Returns (U of b, p, method). U of b is the number of pairs in which the
    value of b is the larger one, ties counting one half. For small samples
    p is exact: the share of all ways of dividing the pooled ranks into two
    groups of these sizes whose rank sum lies at least as far from its
    expectation as the observed one (ties are handled by mid-ranks). For
    larger samples (more than ROBUST_EXACT_LIMIT ways) the usual normal
    approximation is used, with the correction for ties and the continuity
    correction: z = (|U - nA nB / 2| - 1/2) / sd, where
    sd^2 = nA nB / 12 * ((N + 1) - sum(t^3 - t) / (N (N - 1))), N = nA + nB,
    and t runs over the sizes of the groups of tied values."""
    n_a = len(a)
    n_b = len(b)
    ranks = robust_midranks(list(a) + list(b))
    rank_sum_b = sum(ranks[n_a:])
    u_b = rank_sum_b - n_b * (n_b + 1) / 2.0
    arrangements = robust_combinations(n_a + n_b, min(n_a, n_b))
    if arrangements <= ROBUST_EXACT_LIMIT:
        # Enumerate the smaller group; the two rank sums add to a constant,
        # so both give the same distance from expectation.
        small = n_b if n_b <= n_a else n_a
        observed = rank_sum_b if n_b <= n_a else sum(ranks[:n_a])
        expected = small * (n_a + n_b + 1) / 2.0
        limit = abs(observed - expected) - 1e-9
        extreme = 0
        for chosen in itertools.combinations(ranks, small):
            if abs(sum(chosen) - expected) >= limit:
                extreme += 1
        return (u_b, extreme / float(arrangements), "exact permutation")
    total = n_a + n_b
    tie_term = 0.0
    tied = {}
    for rank in ranks:
        tied[rank] = tied.get(rank, 0) + 1
    for rank in tied:
        tie_term += tied[rank] ** 3 - tied[rank]
    variance = n_a * n_b / 12.0 * ((total + 1) - tie_term / float(total * (total - 1)))
    method = "normal approximation (tie and continuity corrected)"
    if variance <= 0:
        return (u_b, 1.0, method)
    z = max(0.0, abs(u_b - n_a * n_b / 2.0) - 0.5) / math.sqrt(variance)
    return (u_b, math.erfc(z / math.sqrt(2.0)), method)


def robust_two_group_statistics(a, b):
    """Compare the unit values of group A and group B.

    difference = mean(B) - mean(A); ratio = mean(B) / mean(A).
    Hedges g = J * difference / pooled SD, with the pooled SD
    sqrt(((nA-1) varA + (nB-1) varB) / (nA+nB-2)) and the small-sample
    correction J = 1 - 3 / (4 (nA+nB-2) - 1).
    Welch test: t = difference / sqrt(varA/nA + varB/nB) with the
    Welch-Satterthwaite degrees of freedom; two-sided p; the 95 percent
    interval of the difference uses the same standard error and degrees of
    freedom. With fewer than ROBUST_MIN_UNITS units in a group only the
    means, their difference, and their ratio are returned."""
    a = [float(value) for value in a]
    b = [float(value) for value in b]
    result = {
        "n_a": len(a), "n_b": len(b), "mean_a": robust_mean(a), "mean_b": robust_mean(b),
        "sd_a": robust_sd(a), "sd_b": robust_sd(b), "difference": None, "ratio": None,
        "hedges_g": None, "welch_t": None, "welch_df": None, "welch_p": None,
        "ci_low": None, "ci_high": None, "mw_u": None, "mw_p": None, "mw_method": "",
        "sign": None, "note": ""
    }
    if len(a) == 0 or len(b) == 0:
        result["note"] = "A group has no unit with a value."
        return result
    difference = result["mean_b"] - result["mean_a"]
    scale = max(abs(result["mean_a"]), abs(result["mean_b"]))
    result["difference"] = difference
    result["sign"] = robust_sign(difference, scale)
    if result["mean_a"] != 0:
        result["ratio"] = result["mean_b"] / result["mean_a"]
    if min(len(a), len(b)) < ROBUST_MIN_UNITS:
        result["note"] = (
            "Fewer than " + str(ROBUST_MIN_UNITS) + " units in a group: means only, "
            "no effect size, interval, or test."
        )
        return result
    notes = []
    var_a = robust_variance(a)
    var_b = robust_variance(b)
    freedom = len(a) + len(b) - 2
    pooled = math.sqrt(((len(a) - 1) * var_a + (len(b) - 1) * var_b) / float(freedom))
    part_a = var_a / len(a)
    part_b = var_b / len(b)
    error = math.sqrt(part_a + part_b)
    if pooled <= 1e-12 * scale or error <= 1e-12 * scale:
        notes.append("No variation between the units of the groups: effect size and Welch test are undefined.")
    else:
        result["hedges_g"] = (1.0 - 3.0 / (4.0 * freedom - 1.0)) * difference / pooled
        result["welch_t"] = difference / error
        result["welch_df"] = (part_a + part_b) ** 2 / (
            part_a * part_a / (len(a) - 1) + part_b * part_b / (len(b) - 1)
        )
        if ROBUST_MATH["available"]:
            result["welch_p"] = RobustTTest().tTest(robust_java_array(a, "d"), robust_java_array(b, "d"))
            quantile = RobustTDistribution(result["welch_df"]).inverseCumulativeProbability(0.975)
            result["ci_low"] = difference - quantile * error
            result["ci_high"] = difference + quantile * error
        else:
            notes.append("Commons Math is missing: no Welch p value and no interval.")
    u_b, p, method = robust_mann_whitney(a, b)
    result["mw_u"] = u_b
    result["mw_p"] = p
    result["mw_method"] = method
    result["note"] = " ".join(notes)
    return result


# -----------------------------------------------------------------------------
# The grid of specifications
# -----------------------------------------------------------------------------

def parse_robust_grid(text):
    """Parse 'threshold_factor=0.8,1.2; min_size=10,20; watershed=true,false'
    into a dictionary of setting -> list of values (as text)."""
    grid = {}
    for part in to_text(text).replace("\n", ";").split(";"):
        part = part.strip()
        if part == "":
            continue
        if "=" not in part:
            raise ValueError("Robustness grid: write each entry as setting=value,value (found '" + part + "').")
        key = part.split("=", 1)[0].strip()
        values = [item.strip() for item in part.split("=", 1)[1].split(",") if item.strip() != ""]
        if key in grid:
            raise ValueError("Robustness grid: '" + key + "' is listed twice.")
        if len(values) == 0:
            raise ValueError("Robustness grid: '" + key + "' has no values.")
        grid[key] = values
    return grid


def robust_grid_text(axes):
    """The checked grid as one line of text that parse_robust_grid reads back."""
    return "; ".join([axis[0] + "=" + ",".join([robust_value_text(value) for value in axis[1]]) for axis in axes])


def robust_parse_value(config, key, raw):
    """One grid value, checked and typed exactly as the profile would."""
    if key == ROBUST_FACTOR_KEY:
        return parse_float_value(raw, "Robustness grid threshold_factor", 0.000001)
    if key in CHANNEL_BOOLEAN_KEYS and not isinstance(raw, bool):
        if to_text(raw).strip().lower() not in ["1", "true", "yes", "on", "0", "false", "no", "off"]:
            raise ValueError("Robustness grid: " + key + " takes true or false, not '" + to_text(raw) + "'.")
    merged = dict(config)
    merged[key] = raw.strip() if isinstance(raw, text_type) else raw
    value = normalize_channel_config(config["channel"], merged)[key]
    if key == "auto_method" and value not in AUTO_METHODS:
        raise ValueError("Robustness grid: unknown automatic threshold method '" + value + "'.")
    return value


def robust_grid_values(config, key, raw_values):
    """The values of one grid setting: typed, without duplicates, with the
    configured value always included, numbers in ascending order."""
    if isinstance(raw_values, (list, tuple)):
        items = list(raw_values)
    elif isinstance(raw_values, bool) or is_number(raw_values):
        items = [raw_values]
    else:
        items = [part for part in to_text(raw_values).split(",") if part.strip() != ""]
    values = []
    for item in items:
        value = robust_parse_value(config, key, item)
        if value not in values:
            values.append(value)
    reference = 1.0 if key == ROBUST_FACTOR_KEY else config[key]
    if reference not in values:
        values.insert(0, reference)
    if key not in CHANNEL_STRING_KEYS:
        values.sort()
    return values


def robust_variant(config, chosen, metric):
    """Apply one combination of grid values to the channel configuration.

    Returns (configuration, threshold factor, names of the grid settings
    that have no effect in this combination). A setting without effect (for
    example the local radius while the threshold is automatic, or an object
    filter while the metric is the positive-area fraction) is put back to
    its configured value, so that combinations that differ only in such a
    setting are recognised as one specification."""
    variant = dict(config)
    factor = 1.0
    for key in chosen:
        if key == ROBUST_FACTOR_KEY:
            factor = chosen[key]
        else:
            variant[key] = chosen[key]
    unused = []
    mode = variant["threshold_mode"]
    if mode != "Auto":
        unused += ["auto_method", "auto_scale"]
    if mode != "Manual":
        unused += ["manual_low", "manual_high"]
    elif variant["manual_low"] == 0:
        # The factor multiplies the lower level, and 0 times any factor is
        # still 0: every factor would segment the same pixels.
        unused.append(ROBUST_FACTOR_KEY)
    elif variant["manual_low"] < 0 and factor != 1.0:
        raise ValueError(
            "Robustness grid: threshold_factor multiplies the Manual lower threshold, which is negative here (" +
            robust_value_text(variant["manual_low"]) + "), so a factor below 1 would raise it. "
            "Vary manual_low directly instead."
        )
    if mode != "Local Phansalkar":
        unused += ROBUST_LOCAL_KEYS
    else:
        unused.append(ROBUST_FACTOR_KEY)
        if not variant["clahe"]:
            unused += ROBUST_CLAHE_KEYS
    if not variant["subtract_background"]:
        unused.append("background_radius")
    if sum([variant[key] for key in ROBUST_MORPHOLOGY_KEYS]) == 0:
        unused.append("morph_count")
    if metric == ROBUST_METRIC_AREA_FRACTION or not needs_objects(variant):
        unused += ROBUST_OBJECT_KEYS
    for key in unused:
        if key == ROBUST_FACTOR_KEY:
            factor = 1.0
        else:
            variant[key] = config[key]
    if mode != "Local Phansalkar":
        # CLAHE exists only for the local threshold; the engine refuses it elsewhere.
        variant["clahe"] = False
    return (variant, factor, [key for key in unused if key in chosen])


def robust_check_image():
    """A one-pixel image with every channel, used only to run the engine's
    own consistency checks on each specification before any file is read."""
    stack = ImageStack(1, 1)
    for index in range(MAX_CHANNELS):
        stack.addSlice("", ByteProcessor(1, 1))
    imp = ImagePlus("CellQuant_robust_check", stack)
    imp.setDimensions(MAX_CHANNELS, 1, 1)
    return imp


def robust_check_specification(spec, probe):
    try:
        validate_unattended_configs([spec["config"]], "A robustness analysis")
        check = dict(spec["config"])
        check["projection"] = "Max Z"
        validate_channel_configs(probe, [check])
    except ValueError as exc:
        raise ValueError("Robustness specification " + spec["label"] + " is not usable: " + to_text(exc))


def robust_build_specifications(config, grid, metric, max_specifications=ROBUST_DEFAULT_MAX_SPECIFICATIONS):
    """Turn a grid into the list of specifications to evaluate.

    Every combination of the listed values (full factorial) is one
    specification; the configured setting is always among them and is the
    reference. Combinations that segment identically are kept once. More
    specifications than max_specifications are refused rather than cut,
    because dropping some would be an arbitrary choice."""
    if not isinstance(grid, dict):
        grid = parse_robust_grid(grid)
    for key in grid:
        if key not in ROBUST_GRID_KEYS:
            raise ValueError(
                "Robustness grid: '" + to_text(key) + "' is not a setting that can be varied. "
                "Choose from: " + ", ".join(ROBUST_GRID_KEYS) + "."
            )
    axes = []
    combinations = 1
    for key in ROBUST_GRID_KEYS:
        if key in grid:
            values = robust_grid_values(config, key, grid[key])
            axes.append((key, values))
            combinations *= len(values)
    if combinations > ROBUST_RAW_GRID_LIMIT:
        raise ValueError(
            "The robustness grid holds " + str(combinations) + " combinations. Use fewer settings or fewer values."
        )
    reference = [1.0 if axis[0] == ROBUST_FACTOR_KEY else config[axis[0]] for axis in axes]
    specs = []
    seen = set()
    prep_groups = {}
    for combo in itertools.product(*[axis[1] for axis in axes]):
        chosen = {}
        for position in range(len(axes)):
            chosen[axes[position][0]] = combo[position]
        variant, factor, unused = robust_variant(config, chosen, metric)
        values = []
        for axis in axes:
            key = axis[0]
            values.append((key, factor if key == ROBUST_FACTOR_KEY else variant[key], key not in unused))
        signature = tuple([robust_value_text(item[1]) for item in values])
        if signature in seen:
            continue
        seen.add(signature)
        prep = tuple([robust_value_text(variant[key]) for key in ROBUST_PREPARATION_KEYS])
        if prep not in prep_groups:
            prep_groups[prep] = len(prep_groups)
        label = "; ".join([item[0] + "=" + robust_value_text(item[1]) for item in values if item[2]])
        specs.append({
            "config": variant, "factor": factor, "values": values, "prep_group": prep_groups[prep],
            "label": label if label != "" else "configured settings",
            "is_reference": [item[1] for item in values] == reference
        })
    for axis in axes:
        if len(axis[1]) < 2:
            continue
        used = False
        for spec in specs:
            for item in spec["values"]:
                if item[0] == axis[0] and item[2]:
                    used = True
        if not used:
            raise ValueError(
                "Robustness grid: '" + axis[0] + "' has no effect on any specification with these settings and "
                "this metric (for example a local-threshold setting while the channel uses an automatic "
                "threshold, an object filter while the metric is the positive-area fraction, or threshold_factor "
                "while the Manual lower threshold is 0). Remove it."
            )
    if len(specs) < 2:
        raise ValueError("The robustness grid holds only the configured setting. List at least one other value.")
    if len(specs) > int(max_specifications):
        raise ValueError(
            "The robustness grid gives " + str(len(specs)) + " specifications; the limit is " +
            str(int(max_specifications)) + ". Use fewer settings or values, or raise the limit on purpose."
        )
    probe = robust_check_image()
    for index in range(len(specs)):
        specs[index]["index"] = index
        specs[index]["id"] = "S" + str(index + 1).zfill(3)
        robust_check_specification(specs[index], probe)
    return {"specs": specs, "axes": axes, "reference_values": reference}


# -----------------------------------------------------------------------------
# Groups and experimental units
# -----------------------------------------------------------------------------

def robust_group_rule(groups):
    """Check and complete the description of how images get their group."""
    raw = groups if groups is not None else {}
    rule = {
        "source": to_text(raw.get("source", ROBUST_GROUP_PREFIX)),
        "separator": to_text(raw.get("separator", "_")),
        "group_regex": to_text(raw.get("group_regex", "")),
        "unit_regex": to_text(raw.get("unit_regex", "")),
        "csv_path": to_text(raw.get("csv_path", "")).strip(),
        "control_group": to_text(raw.get("control_group", "")).strip()
    }
    if rule["source"] not in ROBUST_GROUP_SOURCES:
        raise ValueError("Unknown way of assigning groups: " + rule["source"])
    if rule["source"] == ROBUST_GROUP_PREFIX and rule["separator"] == "":
        raise ValueError("Give the separator that ends the group name in the file name (for example _).")
    used = ["group_regex"] if rule["source"] == ROBUST_GROUP_REGEX else []
    if rule["unit_regex"] != "":
        used.append("unit_regex")
    for key in used:
        try:
            compiled = re.compile(rule[key], re.UNICODE)
        except ANY_ERROR as exc:
            raise ValueError("The " + key.replace("_", " ") + " is not a valid regular expression: " + to_text(exc))
        if compiled.groups != 1:
            raise ValueError(
                "The " + key.replace("_", " ") + " must contain exactly one pair of parentheses around "
                "the part of the path that names the " + key.split("_")[0] + "."
            )
    if rule["source"] == ROBUST_GROUP_CSV and not os.path.isfile(rule["csv_path"]):
        raise ValueError("The group table does not exist: " + rule["csv_path"])
    return rule


def robust_rule_text(rule):
    text = rule["source"]
    if rule["source"] == ROBUST_GROUP_PREFIX:
        text += " '" + rule["separator"] + "'"
    elif rule["source"] == ROBUST_GROUP_REGEX:
        text += ": " + rule["group_regex"]
    else:
        text += ": " + to_text(os.path.basename(rule["csv_path"]))
    if rule["unit_regex"] != "":
        text += "; unit = regular expression on the relative path: " + rule["unit_regex"]
    elif rule["source"] == ROBUST_GROUP_CSV:
        text += "; unit = column Unit, or the image file where it is empty"
    else:
        text += "; unit = the image file"
    return text


def robust_parse_table(text):
    """Rows of a small CSV table (comma, semicolon, or tab; quotes allowed)."""
    first_line = text.split("\n")[0]
    delimiter = ","
    for candidate in [";", "\t"]:
        if first_line.count(candidate) > first_line.count(delimiter):
            delimiter = candidate
    rows = []
    row = []
    cell = []
    quoted = False
    position = 0
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    while position < len(text):
        char = text[position]
        if quoted:
            if char == '"' and text[position + 1:position + 2] == '"':
                cell.append('"')
                position += 1
            elif char == '"':
                quoted = False
            else:
                cell.append(char)
        elif char == '"':
            quoted = True
        elif char == delimiter:
            row.append(u"".join(cell).strip())
            cell = []
        elif char == "\n":
            row.append(u"".join(cell).strip())
            cell = []
            rows.append(row)
            row = []
        else:
            cell.append(char)
        position += 1
    if len(cell) > 0 or len(row) > 0:
        row.append(u"".join(cell).strip())
        rows.append(row)
    return [item for item in rows if len([value for value in item if value != ""]) > 0]


def robust_read_group_table(path):
    """Read the group table: a column Image (file name or path relative to
    the input folder), a column Group, and optionally a column Unit."""
    text = None
    for encoding in ["utf-8-sig", "cp1252", "latin-1"]:
        try:
            handle = codecs.open(path, "r", encoding)
            try:
                text = handle.read()
            finally:
                handle.close()
            break
        except (UnicodeError, LookupError):
            text = None
    if text is None:
        raise ValueError("Could not read the group table: " + to_text(path))
    rows = robust_parse_table(text)
    if len(rows) == 0:
        raise ValueError("The group table is empty: " + to_text(path))
    headers = [value.strip().lower().replace(" ", "_") for value in rows[0]]
    columns = {}
    for name in ["image", "image_relative_path", "image_file", "file", "file_name", "filename", "path"]:
        if name in headers and "image" not in columns:
            columns["image"] = headers.index(name)
    for name in ["group", "unit"]:
        if name in headers:
            columns[name] = headers.index(name)
    if "image" not in columns or "group" not in columns:
        raise ValueError(
            "The group table needs a column named Image and a column named Group (and may have a column "
            "Unit). Found: " + ", ".join(rows[0])
        )
    entries = []
    for row in rows[1:]:
        entry = {}
        for name in columns:
            entry[name] = row[columns[name]] if columns[name] < len(row) else u""
        entries.append(entry)
    return entries


def robust_path_keys(relative):
    """The names under which a table may refer to an image, lower case."""
    relative = to_text(relative).replace("\\", "/").lower()
    name = relative.split("/")[-1]
    keys = [relative, name]
    for value in [relative, name]:
        dot = value.rfind(".")
        if dot > value.rfind("/") + 1:
            keys.append(value[:dot])
    return keys


def robust_assign_from_table(relatives, rule):
    entries = robust_read_group_table(rule["csv_path"])
    index = {}
    exact = {}
    for position in range(len(relatives)):
        keys = robust_path_keys(relatives[position])
        # The full relative path names one image even when a file of the
        # same name lies in a subfolder, so it is looked up first.
        full = keys[0]
        exact.setdefault(full, []).append(position)
        dot = full.rfind(".")
        if dot > full.rfind("/") + 1:
            exact.setdefault(full[:dot], []).append(position)
        for key in set(keys):
            index.setdefault(key, []).append(position)
    rows = []
    for relative in relatives:
        rows.append({"relative": relative, "group": u"", "unit": u"", "status": "NOT ASSIGNED",
                     "message": "Not listed in the group table."})
    extra = []
    for entry in entries:
        key = entry["image"].replace("\\", "/").lower()
        matches = exact.get(key, [])
        if len(matches) != 1:
            matches = index.get(key, [])
        problem = ""
        if len(matches) == 0:
            problem = "No image with this name was found."
        elif len(matches) > 1:
            problem = "Several images have this name; write the path relative to the input folder."
        elif entry["group"] == "":
            problem = "The table gives no group for this image."
        if problem != "":
            if len(matches) == 1:
                rows[matches[0]]["message"] = problem
            else:
                extra.append({"relative": entry["image"], "group": entry["group"], "unit": entry.get("unit", u""),
                              "status": "NO IMAGE FOR TABLE ROW", "message": problem})
            continue
        row = rows[matches[0]]
        unit = entry.get("unit", u"")
        listed = row["status"] == "ASSIGNED"
        if listed and (row["group"] != entry["group"] or (unit != "" and row["unit"] != "" and row["unit"] != unit)):
            row["status"] = "NOT ASSIGNED"
            row["message"] = "Listed more than once with different groups or units."
            row["conflict"] = True
        elif not row.get("conflict", False):
            row["status"] = "ASSIGNED"
            row["message"] = ""
            row["group"] = entry["group"]
            # A repeated row with an empty Unit keeps the unit given before.
            if unit != "" or not listed:
                row["unit"] = unit
    return rows + extra


def robust_assign_groups(relatives, rule):
    """Give every image (by its path relative to the input folder) a group
    and an experimental unit. Returns one dictionary per image, plus one per
    table row that matches no image, with status and reason."""
    if rule["source"] == ROBUST_GROUP_CSV:
        rows = robust_assign_from_table(relatives, rule)
    else:
        rows = []
        pattern = re.compile(rule["group_regex"], re.UNICODE) if rule["source"] == ROBUST_GROUP_REGEX else None
        for relative in relatives:
            row = {"relative": relative, "group": u"", "unit": u"", "status": "NOT ASSIGNED", "message": ""}
            forward = to_text(relative).replace("\\", "/")
            if pattern is None:
                name = forward.split("/")[-1]
                cut = name.find(rule["separator"])
                if cut > 0:
                    row["group"] = name[:cut]
                else:
                    row["message"] = "The file name has no text before a '" + rule["separator"] + "'."
            else:
                match = pattern.search(forward)
                if match is not None and match.group(1) is not None and match.group(1) != "":
                    row["group"] = match.group(1)
                else:
                    row["message"] = "The group expression does not match this path."
            if row["group"] != "":
                row["status"] = "ASSIGNED"
            rows.append(row)
    unit_pattern = re.compile(rule["unit_regex"], re.UNICODE) if rule["unit_regex"] != "" else None
    for row in rows:
        if row["status"] != "ASSIGNED":
            continue
        forward = to_text(row["relative"]).replace("\\", "/")
        if unit_pattern is not None:
            match = unit_pattern.search(forward)
            if match is not None and match.group(1) is not None and match.group(1) != "":
                row["unit"] = match.group(1)
            else:
                row["status"] = "NOT ASSIGNED"
                row["message"] = "The unit expression does not match this path."
        elif row["unit"] == "":
            row["unit"] = forward
    return rows


def robust_group_order(names, control_group):
    """Group names in alphabetical order, the control group first. Group A
    of the statistics is the first name, group B the second."""
    ordered = sorted(names, key=lambda name: (to_text(name).lower(), to_text(name)))
    if control_group != "":
        if control_group not in ordered:
            raise ValueError(
                "The control group '" + control_group + "' is not among the groups found: " + ", ".join(ordered)
            )
        ordered.remove(control_group)
        ordered.insert(0, control_group)
    return ordered


# -----------------------------------------------------------------------------
# Measuring one image under every specification
# -----------------------------------------------------------------------------

def robust_metric_value(metric, roi_area, count, total_area, positive_area):
    """(value, numerator, denominator) of one measurement.

    Object count density: objects / ROI area. Positive-area fraction:
    positive area / ROI area * 100. Mean object area: total object area /
    objects. The arithmetic follows roi_summary.csv exactly. The value is
    None where the denominator is zero."""
    if metric == ROBUST_METRIC_AREA_FRACTION:
        numerator, denominator = positive_area, roi_area
    elif metric == ROBUST_METRIC_DENSITY:
        numerator, denominator = count, roi_area
    else:
        numerator, denominator = total_area, count
    if numerator is None or denominator is None or denominator <= 0:
        return (None, numerator, denominator)
    return (float(numerator) / float(denominator) * ROBUST_METRIC_SCALE[metric], numerator, denominator)


def robust_count_objects(engine, mask_imp, parent_roi, config, get_masks):
    """Number and total area of the objects of one mask inside one ROI.

    The answer is the one analyze_objects gives. Without edge exclusion it
    is obtained faster: the same particle analysis on the same mask, cleared
    outside the true outline of the ROI, but collecting only the areas
    instead of one ROI and a dozen measurements per object. With edge
    exclusion analyze_objects itself is used."""
    if config["exclude_edges"]:
        objects, info = engine.analyze_objects(mask_imp, parent_roi, config, get_masks())
        areas = [obj["area"] for obj in objects]
        return (len(objects), sum(areas) if len(areas) > 0 else 0.0)
    bounds = parent_roi.getBounds()
    work_ip = mask_imp.getProcessor().duplicate()
    work_ip.resetRoi()
    if parent_roi.getMask() is not None:
        work_ip.copyBits(get_masks()[0], 0, 0, Blitter.AND)
    work_imp = ImagePlus("CellQuant_robust_particles", work_ip)
    work_imp.setCalibration(mask_imp.getCalibration().copy())
    options = ParticleAnalyzer.SHOW_NONE
    if config["include_holes"]:
        options |= ParticleAnalyzer.INCLUDE_HOLES
    if config["composite_rois"]:
        options |= ParticleAnalyzer.COMPOSITE_ROIS
    if config["four_connected"]:
        options |= ParticleAnalyzer.FOUR_CONNECTED
    min_px, max_px = engine.particle_size_pixels(config, mask_imp.getCalibration())
    table = ResultsTable()
    work_imp.setRoi(Roi(bounds.x, bounds.y, bounds.width, bounds.height))
    work_ip.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
    analyzer = ParticleAnalyzer(
        options, Measurements.AREA, table, min_px, max_px, config["min_circularity"], config["max_circularity"]
    )
    if not analyzer.analyze(work_imp):
        raise RuntimeError("Particle analysis failed for C" + str(config["channel"]) + ".")
    work_imp.deleteRoi()
    count = table.size()
    total = 0.0
    if count > 0:
        for area in table.getColumnAsDoubles(table.getColumnIndex("Area")):
            total += area
    return (count, total)


def robust_resolve_threshold(engine, state, spec):
    """The configuration that is segmented on this image: the specification
    with its threshold factor turned into a threshold level.

    Manual threshold: the configured lower level times the factor. Automatic
    threshold: the level the method finds on this image, times the configured
    multiplier, times the factor; the result is used as a Manual threshold,
    exactly as the sensitivity sweep does (sensitivity_variant). With a dark
    background the factor moves the lower level, with a light background
    the upper one. Factor 1.0 changes nothing: the reference specification
    is segmented by the engine exactly as in a normal run."""
    variant = spec["config"]
    factor = spec["factor"]
    if factor == 1.0:
        return variant
    scaled = dict(variant)
    if variant["threshold_mode"] == "Manual":
        scaled["manual_low"] = variant["manual_low"] * factor
        return scaled
    # The level of the method alone (multiplier 1) always lies inside the
    # range of the image, so this never stops on an empty range, and one
    # look at the image serves every multiplier in the grid.
    key = (variant["auto_method"], variant["dark_background"])
    baselines = state["baselines"]
    if key not in baselines:
        plain = dict(variant)
        for name in ROBUST_MORPHOLOGY_KEYS:
            plain[name] = 0
        plain["fill_holes"] = False
        plain["watershed"] = False
        plain["auto_scale"] = 1.0
        data = engine.build_masks(state["prepared"], plain, "CellQuant_robust_baseline")
        safe_close_image(data["raw"])
        safe_close_image(data["final"])
        baselines[key] = (float(data["threshold_low"]), float(data["threshold_high"]))
    low, high = baselines[key]
    data_min, data_max = state["range"]
    scaled["threshold_mode"] = "Manual"
    if variant["dark_background"]:
        # build_masks never reports a level below the darkest pixel.
        level = low * variant["auto_scale"]
        if level <= data_max:
            level = max(data_min, level)
        scaled["manual_low"] = level * factor
        scaled["manual_high"] = high
    else:
        level = high * variant["auto_scale"]
        if level >= data_min:
            level = min(data_max, level)
        scaled["manual_low"] = low
        scaled["manual_high"] = level * factor
    return scaled


def robust_threshold_is_empty(variant, data_min, data_max):
    """True when a Manual threshold range holds no pixel value of the image
    (the same condition for which build_masks stops with an error)."""
    if variant["threshold_mode"] != "Manual":
        return False
    high = data_max if variant["manual_high"] == -1 else variant["manual_high"]
    return min(data_max, high) < max(data_min, variant["manual_low"])


def robust_measure_specification(engine, state, spec, rois, roi_info, metric, get_masks):
    """Segment the prepared channel with one specification and measure the
    metric in every parent ROI. Returns one result dictionary per ROI. An
    error of this specification on this image is recorded, not raised."""
    prepared = state["prepared"]
    mask_data = None
    status = "OK"
    note = ""
    thresholds = ("", "")
    measured = []
    try:
        variant = robust_resolve_threshold(engine, state, spec)
        empty = robust_threshold_is_empty(variant, state["range"][0], state["range"][1])
        if empty:
            thresholds = (variant["manual_low"], variant["manual_high"] if variant["manual_high"] != -1 else state["range"][1])
        else:
            mask_data = engine.build_masks(prepared, variant, "CellQuant_robust")
            thresholds = (mask_data["threshold_low"], mask_data["threshold_high"])
            # The engine gives an empty mask, marked so, when the threshold
            # range holds no pixel value of the image.
            if mask_data.get("empty_range", False):
                empty = True
        if empty:
            status = "EMPTY"
            note = "The threshold range holds no pixel value of this image: nothing is positive."
        for roi_index in range(len(rois)):
            count = None
            total_area = None
            positive_area = None
            if metric == ROBUST_METRIC_AREA_FRACTION:
                positive_area = 0.0
                if not empty:
                    positive_area = engine.measure_positive_area(mask_data["raw"], rois[roi_index])["positive_area"]
            elif empty:
                count, total_area = 0, 0.0
            else:
                count, total_area = robust_count_objects(
                    engine, mask_data["final"], rois[roi_index], variant,
                    lambda _index=roi_index: get_masks(_index)
                )
            measured.append((count, total_area, positive_area))
    except AnalysisCancelled:
        raise
    except ANY_ERROR as exc:
        status = "FAILED"
        note = to_text(exc)
        measured = [(None, None, None)] * len(rois)
    finally:
        if mask_data is not None:
            safe_close_image(mask_data["raw"])
            safe_close_image(mask_data["final"])
    results = []
    for roi_index in range(len(rois)):
        count, total_area, positive_area = measured[roi_index]
        value, numerator, denominator = (None, None, None)
        row_status = status
        row_note = note
        if status != "FAILED":
            value, numerator, denominator = robust_metric_value(
                metric, roi_info[roi_index]["area"], count, total_area, positive_area
            )
            if value is None:
                row_status = "UNDEFINED"
                row_note = (note + " " if note != "" else "") + (
                    "No object: the mean object area is undefined." if metric == ROBUST_METRIC_OBJECT_AREA
                    else "The ROI has no area."
                )
        results.append({
            "spec": spec["index"], "roi": roi_index, "value": value, "numerator": numerator,
            "denominator": denominator, "count": count, "total_area": total_area,
            "positive_area": positive_area, "low": thresholds[0], "high": thresholds[1],
            "status": row_status, "note": row_note
        })
    return results


def robust_measure_frame(engine, imp, frame, specs, rois, calibration, metric):
    """Measure every specification on one frame of an opened image.

    The channel is prepared once for each distinct set of preparation
    settings in the grid (once in all when the grid varies segmentation
    settings only). Returns (results, ROI information)."""
    width = imp.getWidth()
    height = imp.getHeight()
    mask_cache = {}

    def get_masks(roi_index):
        if roi_index not in mask_cache:
            mask_cache[roi_index] = roi_masks(rois[roi_index], width, height)
        return mask_cache[roi_index]

    state = {"prepared": None, "group": None, "range": None, "baselines": {}}
    roi_info = None
    results = []
    try:
        for spec in sorted(specs, key=lambda item: (item["prep_group"], item["index"])):
            engine.check_cancelled()
            if state["prepared"] is None or state["group"] != spec["prep_group"]:
                safe_close_image(state["prepared"])
                state["prepared"] = None
                state["prepared"] = engine.prepare_single_channel(imp, spec["config"], frame, calibration)
                state["group"] = spec["prep_group"]
                stats = state["prepared"].getStatistics(Measurements.MIN_MAX)
                state["range"] = (stats.min, stats.max)
                state["baselines"] = {}
                if roi_info is None:
                    roi_info = []
                    cal = state["prepared"].getCalibration()
                    unit = to_text(cal.getUnit()) if cal.getUnit() else u"pixel"
                    for roi in rois:
                        roi_info.append({"area": engine.measure_intensity(state["prepared"], roi)["area"], "unit": unit})
            engine.show_status("Cell Quant Studio robustness: specification " + spec["id"])
            results += robust_measure_specification(engine, state, spec, rois, roi_info, metric, get_masks)
    finally:
        safe_close_image(state["prepared"])
    results.sort(key=lambda item: (item["spec"], item["roi"]))
    return (results, roi_info)


def robust_metric_unit(metric, units):
    """Unit of the metric, from the spatial units of the measured images."""
    if metric == ROBUST_METRIC_AREA_FRACTION:
        return u"percent"
    if len(units) != 1:
        return u"mixed units" if len(units) > 1 else u""
    unit = list(units)[0]
    if metric == ROBUST_METRIC_DENSITY:
        return u"objects per " + unit + u"^2"
    return unit + u"^2"


# -----------------------------------------------------------------------------
# From measurements to units, groups, and the summary
# -----------------------------------------------------------------------------

def robust_unit_values(records, aggregation, metric):
    """Combine the measurements of each experimental unit.

    records: dictionaries with spec, group, unit, image, value, numerator,
    denominator, status. Returns spec index -> {"units": {(group, unit):
    entry}, "measurements": n, "failed": n}; entry["value"] is the mean of
    the unit's measurement values, or with pooling the sum of numerators
    over the sum of denominators (an area-weighted mean for the two
    ROI-area metrics). A unit belongs to one group; the same unit name in
    two groups is two units."""
    by_spec = {}
    for record in records:
        info = by_spec.setdefault(record["spec"], {"units": {}, "measurements": 0, "failed": 0})
        entry = info["units"].setdefault((record["group"], record["unit"]), {
            "values": [], "numerator": 0.0, "denominator": 0.0, "images": set(), "count": 0, "value": None
        })
        info["measurements"] += 1
        entry["count"] += 1
        entry["images"].add(record["image"])
        if record["status"] == "FAILED":
            info["failed"] += 1
            continue
        if record["value"] is not None:
            entry["values"].append(record["value"])
        if record["numerator"] is not None and record["denominator"] is not None:
            entry["numerator"] += record["numerator"]
            entry["denominator"] += record["denominator"]
    for spec_index in by_spec:
        for key in by_spec[spec_index]["units"]:
            entry = by_spec[spec_index]["units"][key]
            if aggregation == ROBUST_UNIT_POOLED:
                if entry["denominator"] > 0:
                    entry["value"] = entry["numerator"] / entry["denominator"] * ROBUST_METRIC_SCALE[metric]
            else:
                entry["value"] = robust_mean(entry["values"])
    return by_spec


def robust_specification_statistics(spec, info, group_names):
    """Group values and, with exactly two groups, the comparison of B with A
    for one specification. Returns (per-group list, statistics or None)."""
    groups = []
    for name in group_names:
        values = []
        for key in sorted(info["units"].keys()):
            if key[0] == name and info["units"][key]["value"] is not None:
                values.append(info["units"][key]["value"])
        groups.append({"group": name, "values": values})
    stats = None
    if len(group_names) == 2:
        stats = robust_two_group_statistics(groups[0]["values"], groups[1]["values"])
    return (groups, stats)


def robust_range(values):
    if len(values) == 0:
        return ("", "", "")
    return (min(values), robust_median(values), max(values))


def robust_summarize(entries, axes):
    """Summarise the comparison over all specifications.

    entries: one dictionary per specification with spec (from
    robust_build_specifications), stats (from robust_two_group_statistics,
    or None), and complete (no failed measurement). Only complete
    specifications with a defined difference are compared. Returns
    (summary dictionary, per-parameter-value rows)."""
    reference = None
    for entry in entries:
        if entry["spec"]["is_reference"]:
            reference = entry
    compared = [entry for entry in entries
                if entry["complete"] and entry["stats"] is not None and entry["stats"]["sign"] is not None]
    reference_sign = None
    if reference is not None and reference["complete"] and reference["stats"] is not None:
        reference_sign = reference["stats"]["sign"]
    summary = {
        "N_Specifications": len(entries),
        "N_Specifications_Incomplete": len([entry for entry in entries if not entry["complete"]]),
        "N_Specifications_Compared": len(compared),
        "Reference_Specification_ID": reference["spec"]["id"] if reference is not None else "",
        "Reference_Sign": robust_sign_text(reference_sign)
    }
    notes = []
    if reference is not None and reference["stats"] is not None:
        summary["Reference_N_Units_A"] = reference["stats"]["n_a"]
        summary["Reference_N_Units_B"] = reference["stats"]["n_b"]
        summary["Reference_Difference_B_Minus_A"] = robust_cell(reference["stats"]["difference"])
        summary["Reference_Hedges_g"] = robust_cell(reference["stats"]["hedges_g"])
    same = 0
    opposite = 0
    welch = 0
    mann_whitney = 0
    for entry in compared:
        sign = entry["stats"]["sign"]
        entry["same_sign"] = None
        if reference_sign is None or reference_sign == 0:
            continue
        entry["same_sign"] = sign == reference_sign
        if sign == reference_sign:
            same += 1
            if entry["stats"]["welch_p"] is not None and entry["stats"]["welch_p"] < ROBUST_ALPHA:
                welch += 1
            if entry["stats"]["mw_p"] is not None and entry["stats"]["mw_p"] < ROBUST_ALPHA:
                mann_whitney += 1
        elif sign == -reference_sign:
            opposite += 1
    if reference_sign is None:
        notes.append("The reference specification has no defined difference, so directions are not compared.")
    elif reference_sign == 0:
        notes.append("The difference is zero at the reference specification, so directions are not compared.")
    elif len(compared) > 0:
        summary["N_Same_Sign_As_Reference"] = same
        summary["N_Opposite_Sign"] = opposite
        summary["N_Zero_Or_Undefined"] = len(compared) - same - opposite
        summary["Fraction_Same_Sign_As_Reference"] = same / float(len(compared))
        # Left empty when no compared specification has that p value (fewer
        # than 3 units in a group, or no Commons Math): a 0 would read as
        # "tested and never below 0.05".
        for name, key, count in [("Welch", "welch_p", welch), ("Mann_Whitney", "mw_p", mann_whitney)]:
            if len([candidate for candidate in compared if candidate["stats"][key] is not None]) > 0:
                summary["Descriptive_N_" + name + "_P_Below_0.05_Same_Sign"] = count
                summary["Descriptive_Fraction_" + name + "_P_Below_0.05_Same_Sign"] = count / float(len(compared))
    for name, key in [("Hedges_g", "hedges_g"), ("Difference", "difference"), ("Ratio", "ratio")]:
        low, middle, high = robust_range(
            [entry["stats"][key] for entry in compared if entry["stats"][key] is not None]
        )
        summary[name + "_Minimum"] = low
        summary[name + "_Median"] = middle
        summary[name + "_Maximum"] = high

    # Which values of which setting go with a reversed direction. A value is
    # named when the specifications that use it are reversed more often than
    # the specifications as a whole; in a grid of several settings this picks
    # out the setting that matters and leaves out the ones that only happen
    # to be combined with it.
    overall = opposite / float(len(compared)) if len(compared) > 0 else 0.0
    value_rows = []
    flips = []
    for position in range(len(axes)):
        key = axes[position][0]
        for value in axes[position][1]:
            counts = [0, 0, 0]
            for entry in compared:
                item = entry["spec"]["values"][position]
                if not item[2] or item[1] != value:
                    continue
                counts[0] += 1
                if entry.get("same_sign") is True:
                    counts[1] += 1
                elif entry.get("same_sign") is False and entry["stats"]["sign"] != 0:
                    counts[2] += 1
            named = counts[2] > 0 and counts[2] / float(counts[0]) > overall + 1e-12
            value_rows.append({
                "Parameter": key, "Value": robust_value_text(value),
                "Is_Reference_Value": value == (1.0 if key == ROBUST_FACTOR_KEY else reference["spec"]["config"][key])
                if reference is not None else "",
                "N_Specifications": counts[0], "N_Same_Sign_As_Reference": counts[1],
                "N_Opposite_Sign": counts[2], "N_Zero_Or_Undefined": counts[0] - counts[1] - counts[2],
                "Fraction_Opposite_Sign": counts[2] / float(counts[0]) if counts[0] > 0 else "",
                "More_Often_Reversed_Than_Overall": named
            })
            if named:
                flips.append(
                    key + "=" + robust_value_text(value) + " (" + str(counts[2]) + " of " + str(counts[0]) + ")"
                )
    if reference_sign not in [1, -1]:
        summary["Sign_Flipping_Values"] = ""
    elif opposite == 0:
        summary["Sign_Flipping_Values"] = "none"
    elif len(flips) == 0:
        summary["Sign_Flipping_Values"] = (
            "no single value (the reversals are spread evenly over the values of every setting; "
            "see robustness_specifications.csv)"
        )
    else:
        summary["Sign_Flipping_Values"] = "; ".join(flips)
    summary["Note"] = " ".join(notes)
    return (summary, value_rows)


def robust_curve_order(entries):
    """Order the specifications for a specification curve: by Hedges g when
    every compared specification has one, otherwise by the difference of
    the means; specifications without a value come last."""
    with_stats = [entry for entry in entries if entry["stats"] is not None and entry["stats"]["difference"] is not None]
    key = "hedges_g"
    for entry in with_stats:
        if entry["stats"]["hedges_g"] is None:
            key = "difference"
    if len(with_stats) == 0:
        key = ""

    def rank(entry):
        if key == "" or entry["stats"] is None or entry["stats"][key] is None:
            return (1, 0.0, entry["spec"]["index"])
        return (0, entry["stats"][key], entry["spec"]["index"])

    label = {"hedges_g": "Hedges_g", "difference": "Difference_B_Minus_A", "": "Specification_ID"}[key]
    return (sorted(entries, key=rank), label)


# -----------------------------------------------------------------------------
# The complete analysis
# -----------------------------------------------------------------------------

def robust_plan(configs, channel, grid, groups, metric, unit_aggregation, max_specifications):
    """Check everything that can be checked before a file is opened."""
    if len(configs) == 0:
        raise ValueError("No channels are enabled for analysis.")
    validate_unattended_configs(configs, "A robustness analysis")
    channel = parse_int_value(channel, "Robustness channel", 1, MAX_CHANNELS)
    config = None
    for candidate in configs:
        if candidate["channel"] == channel:
            config = candidate
    if config is None:
        raise ValueError("C" + str(channel) + " is not enabled for analysis, so its robustness cannot be studied.")
    metric = to_text(metric)
    if metric == ROBUST_METRIC_INTENSITY:
        raise ValueError(
            "Mean intensity is measured on the prepared image and does not depend on the segmentation "
            "settings, so there is nothing to study. Choose: " + ", ".join(ROBUST_METRICS) + "."
        )
    if metric not in ROBUST_METRICS:
        raise ValueError("Unknown robustness metric: " + metric + ". Choose: " + ", ".join(ROBUST_METRICS) + ".")
    if metric == ROBUST_METRIC_AREA_FRACTION and not needs_positive_area(config):
        raise ValueError("C" + str(channel) + " is not set to a Positive area mode, so it has no positive-area fraction.")
    if metric != ROBUST_METRIC_AREA_FRACTION and not needs_objects(config):
        raise ValueError("C" + str(channel) + " is not set to an Objects mode, so it has no objects to count.")
    if to_text(unit_aggregation) not in ROBUST_UNIT_AGGREGATIONS:
        raise ValueError("Unknown way of combining the measurements of a unit: " + to_text(unit_aggregation))
    plan = robust_build_specifications(config, grid, metric, max_specifications)
    plan["config"] = config
    plan["metric"] = metric
    plan["aggregation"] = to_text(unit_aggregation)
    plan["rule"] = robust_group_rule(groups)
    plan["grid_text"] = robust_grid_text(plan["axes"])
    return plan


def robust_statistics_tables(plan, records, group_names, base, metric_unit):
    """Build every table that follows from the measurements.
    Returns a dictionary of table name -> rows, plus the summary row."""
    specs = plan["specs"]
    by_spec = robust_unit_values(records, plan["aggregation"], plan["metric"])
    metric_cells = {
        "Channel": plan["config"]["channel"], "Channel_Name": plan["config"]["name"],
        "Metric": plan["metric"], "Metric_Unit": metric_unit
    }

    def spec_row(spec, extra):
        row = dict(base)
        row.update(metric_cells)
        row["Specification_ID"] = spec["id"]
        row["Specification"] = spec["label"]
        row["Is_Reference"] = spec["is_reference"]
        row.update(extra)
        return row

    unit_rows = []
    group_rows = []
    entries = []
    for spec in specs:
        info = by_spec.get(spec["index"], {"units": {}, "measurements": 0, "failed": 0})
        for key in sorted(info["units"].keys(), key=lambda item: (group_names.index(item[0]), item[1])):
            entry = info["units"][key]
            unit_rows.append(spec_row(spec, {
                "Group": key[0], "Unit": key[1], "N_Images": len(entry["images"]),
                "N_Measurements": entry["count"], "N_Measurements_Used": len(entry["values"]),
                "Unit_Aggregation": plan["aggregation"], "Value": robust_cell(entry["value"])
            }))
        groups, stats = robust_specification_statistics(spec, info, group_names)
        for group in groups:
            values = group["values"]
            group_rows.append(spec_row(spec, {
                "Group": group["group"], "N_Units": len(values), "Mean": robust_cell(robust_mean(values)),
                "SD": robust_cell(robust_sd(values)),
                "Minimum": min(values) if len(values) > 0 else "",
                "Maximum": max(values) if len(values) > 0 else ""
            }))
        entries.append({
            "spec": spec, "stats": stats, "complete": info["failed"] == 0 and info["measurements"] > 0,
            "measurements": info["measurements"], "failed": info["failed"]
        })

    summary, value_rows = robust_summarize(entries, plan["axes"])
    ordered, sorted_by = robust_curve_order(entries)
    specification_rows = []
    for rank in range(len(ordered)):
        entry = ordered[rank]
        stats = entry["stats"]
        notes = []
        row = spec_row(entry["spec"], {
            "Curve_Rank": rank + 1, "N_Groups": len(group_names),
            "Group_A": group_names[0], "Group_B": group_names[1] if len(group_names) > 1 else "",
            "N_Measurements": entry["measurements"], "N_Failed_Measurements": entry["failed"],
            "Complete": entry["complete"]
        })
        if not entry["complete"]:
            notes.append("A measurement failed for this specification; it is left out of the summary.")
        if stats is None:
            row["Group_B"] = ""
            notes.append("More than two groups: see robustness_groups.csv for the group means.")
        else:
            pairs = [
                ("N_Units_A", "n_a"), ("Mean_A", "mean_a"), ("SD_A", "sd_a"), ("N_Units_B", "n_b"),
                ("Mean_B", "mean_b"), ("SD_B", "sd_b"), ("Difference_B_Minus_A", "difference"),
                ("Ratio_B_Over_A", "ratio"), ("Hedges_g", "hedges_g"), ("Welch_t", "welch_t"),
                ("Welch_df", "welch_df"), ("Welch_P", "welch_p"), ("Difference_CI95_Low", "ci_low"),
                ("Difference_CI95_High", "ci_high"), ("Mann_Whitney_U", "mw_u"), ("Mann_Whitney_P", "mw_p"),
                ("Mann_Whitney_Method", "mw_method")
            ]
            for pair in pairs:
                row[pair[0]] = robust_cell(stats[pair[1]])
            row["Sign"] = robust_sign_text(stats["sign"])
            row["Same_Sign_As_Reference"] = robust_cell(entry.get("same_sign"))
            if stats["note"] != "":
                notes.append(stats["note"])
        row["Note"] = " ".join(notes)
        specification_rows.append(row)

    summary_row = dict(base)
    summary_row.update(metric_cells)
    summary_row.update(summary)
    summary_row.update({
        "Group_A": group_names[0], "Group_B": group_names[1] if len(group_names) == 2 else "",
        "N_Groups": len(group_names), "Unit_Aggregation": plan["aggregation"],
        "Group_Rule": robust_rule_text(plan["rule"]), "Grid": plan["grid_text"],
        "Curve_Sorted_By": sorted_by, "Statistics_Note": ROBUST_STATISTICS_NOTE
    })
    if len(group_names) != 2:
        summary_row["Note"] = (
            str(len(group_names)) + " groups: the comparison needs exactly two, so only the group means "
            "are reported (robustness_groups.csv)."
        )
    for row in value_rows:
        row.update(base)
    return {
        "units": unit_rows, "groups": group_rows, "specifications": specification_rows,
        "values": value_rows, "summary": summary_row
    }


def robust_parameter_rows(plan, base):
    rows = []
    for spec in plan["specs"]:
        for position in range(len(spec["values"])):
            item = spec["values"][position]
            row = dict(base)
            row.update({
                "Specification_ID": spec["id"], "Specification": spec["label"],
                "Is_Reference": spec["is_reference"], "Parameter": item[0],
                "Value": robust_value_text(item[1]),
                "Reference_Value": robust_value_text(plan["reference_values"][position]),
                "Applies": item[2]
            })
            rows.append(row)
    return rows


def robust_process_image(engine, imp, path, item, plan, run, series_label):
    """Measure one opened image; returns the report cells of its row.
    Raises on a failure of the whole image; rows are only kept when every
    frame was measured, so that a broken image leaves nothing half-counted."""
    configs = run["configs"]
    # The image is named after its file, as in a batch.
    file_name = to_text(os.path.basename(path))
    if series_label != "":
        stem = strip_extension(file_name)
        imp.setTitle(stem + "_" + series_label + file_name[len(stem):])
    else:
        imp.setTitle(file_name)
    validate_channel_configs(imp, configs)
    if run["policy"] in BATCH_ROI_PROVIDERS:
        rois, roi_note = BATCH_ROI_PROVIDERS[run["policy"]]({
            "imp": imp, "path": path, "configs": configs,
            "global_options": run["global_options"], "log": run["log"]
        })
    else:
        rois, roi_note = batch_parent_rois(
            imp, path, run["policy"], run["roi_dir"], run["manager_rois"], series_label
        )
    report = {"ROI_Source": roi_note}
    if rois is None:
        report["Status"] = "SKIPPED"
        report["Message"] = roi_note
        return report
    rois = validate_parent_rois(rois, imp.getWidth(), imp.getHeight())
    frames = range(1, imp.getNFrames() + 1) if run["all_frames"] else [1]
    specs = plan["specs"]
    config = plan["config"]
    rows = []
    records = []
    units = set()
    for frame in frames:
        results, roi_info = robust_measure_frame(engine, imp, frame, specs, rois, run["calibration"], plan["metric"])
        for result in results:
            spec = specs[result["spec"]]
            info = roi_info[result["roi"]]
            units.add(info["unit"])
            row = dict(run["base"])
            row.update({
                "Specification_ID": spec["id"], "Specification": spec["label"], "Is_Reference": spec["is_reference"],
                "Image_File": os.path.basename(path), "Image_Relative_Path": item["relative"],
                "Series": series_label, "Frame": frame, "Group": item["group"], "Unit": item["unit"],
                "Parent_ROI_Index": result["roi"] + 1,
                "Parent_ROI_Name": roi_name(rois[result["roi"]], "ROI_" + str(result["roi"] + 1).zfill(3)),
                "Channel": config["channel"], "Channel_Name": config["name"], "Metric": plan["metric"],
                "Metric_Unit": robust_metric_unit(plan["metric"], [info["unit"]]),
                "Value": robust_cell(result["value"]), "ROI_Area": info["area"],
                "ROI_Area_Unit": info["unit"] + u"^2", "Object_Count": robust_cell(result["count"]),
                "Total_Object_Area": robust_cell(result["total_area"]),
                "Positive_Area": robust_cell(result["positive_area"]),
                "Threshold_Low": result["low"], "Threshold_High": result["high"],
                "Status": result["status"], "Note": result["note"]
            })
            rows.append(row)
            records.append({
                "spec": result["spec"], "group": item["group"], "unit": item["unit"],
                "image": item["relative"] + "|" + series_label, "value": result["value"],
                "numerator": result["numerator"], "denominator": result["denominator"],
                "status": result["status"]
            })
    run["measurements"].append(rows)
    run["records"] += records
    run["units"].update(units)
    report.update({
        "Status": "OK", "Message": "", "Parent_ROI_Count": len(rois), "Frames": len(frames),
        "Measurements": len(rows)
    })
    failed = len([record for record in records if record["status"] == "FAILED"])
    if failed > 0:
        report["Message"] = str(failed) + " measurement(s) failed; see robustness_measurements.csv."
    return report


def run_robustness(input_dir, output_root, configs, all_configs, global_options, batch_options,
                   channel, grid, groups, metric=ROBUST_METRIC_DENSITY,
                   unit_aggregation=ROBUST_UNIT_MEAN, max_specifications=ROBUST_DEFAULT_MAX_SPECIFICATIONS,
                   manager_rois=None, log=None, progress=None, is_cancelled=None, file_list=None):
    """Measure a set of images under a grid of settings of one channel and
    report how a group difference depends on them.

    input_dir, output_root, configs, all_configs, global_options,
    batch_options, manager_rois, log, progress, is_cancelled, file_list:
    as for run_batch (the same files are found, opened, and given the same
    parent ROIs as in a batch).
    channel: number of the channel whose settings are varied.
    grid: dictionary of setting -> list of values (see ROBUST_GRID_KEYS), or
          the same as text, for example "threshold_factor=0.8,1.2; min_size=10,20".
    groups: dictionary describing how images get a group and a unit: source
            (one of ROBUST_GROUP_SOURCES), separator, group_regex,
            unit_regex, csv_path, control_group.
    metric: one of ROBUST_METRICS.

    Writes the folder CellQuant_Robustness_<time> under output_root and
    returns a dictionary with its path, the counts of images, and the rows
    of the statistics tables. When cancelled, the measurements made so far
    are kept but no statistics are written."""
    log_lines = []

    def run_log(message):
        log_lines.append("[" + time.strftime("%H:%M:%S") + "] " + to_text(message))
        if log is not None:
            log(message)

    plan = robust_plan(configs, channel, grid, groups, metric, unit_aggregation, max_specifications)
    input_dir = os.path.abspath(to_text(input_dir))
    if not os.path.isdir(input_dir):
        raise ValueError("Input folder does not exist: " + input_dir)
    extensions = parse_extension_list(batch_options.get("extensions", DEFAULT_BATCH_EXTENSIONS))
    if len(extensions) == 0:
        raise ValueError("List at least one file extension.")
    policy = to_text(batch_options.get("roi_policy", BATCH_ROI_WHOLE))
    if policy not in BATCH_ROI_POLICIES:
        raise ValueError("Unsupported batch ROI policy: " + policy)
    roi_dir = to_text(batch_options.get("roi_dir", "")).strip()
    if roi_dir != "" and not os.path.isdir(roi_dir):
        raise ValueError("ROI folder does not exist: " + roi_dir)
    force_bioformats = parse_bool(batch_options.get("force_bioformats"), False)
    all_series = parse_bool(batch_options.get("all_series"), False)
    skip_thumbnails = parse_bool(batch_options.get("skip_thumbnail_series"), True)
    recursive = parse_bool(batch_options.get("recursive"), False)
    if file_list is None:
        files = discover_batch_images(input_dir, extensions, recursive)
    else:
        files = list(file_list)
    if len(files) == 0:
        raise ValueError("No image files with the listed extensions were found in " + input_dir)
    relatives = []
    for path in files:
        try:
            relatives.append(to_text(os.path.relpath(path, input_dir)))
        except ANY_ERROR:
            relatives.append(to_text(os.path.basename(path)))
    rule = plan["rule"]
    assignment = robust_assign_groups(relatives, rule)
    items = []
    for position in range(len(files)):
        if assignment[position]["status"] == "ASSIGNED":
            items.append({"path": files[position], "relative": relatives[position],
                          "group": assignment[position]["group"], "unit": assignment[position]["unit"]})
    group_names = robust_group_order(list(set([item["group"] for item in items])), rule["control_group"])
    if len(group_names) < 2:
        raise ValueError(
            "A comparison needs images of at least two groups, but the rule (" + robust_rule_text(rule) +
            ") gives " + (("only the group '" + group_names[0] + "'") if len(group_names) == 1 else "no group") +
            " for the " + str(len(files)) + " image file(s) found, for example " + relatives[0] + "."
        )

    run_dir = unique_dir(os.path.abspath(to_text(output_root)), "CellQuant_Robustness_" + timestamp_text())
    run_id = os.path.basename(run_dir)
    base = {"Run_ID": run_id, "Application_Version": APP_VERSION}
    # The record holds the settings as they were used, so that the saved
    # file reproduces this analysis when it is loaded as a profile.
    recorded = dict(global_options)
    recorded.update({
        "batch_extensions": ", ".join(extensions), "batch_recursive": recursive, "batch_roi_policy": policy,
        "batch_force_bioformats": force_bioformats, "batch_all_series": all_series,
        "robust_channel": str(plan["config"]["channel"]), "robust_metric": plan["metric"],
        "robust_grid": plan["grid_text"], "robust_max_specifications": str(int(max_specifications)),
        "robust_group_source": rule["source"], "robust_separator": rule["separator"],
        "robust_group_regex": rule["group_regex"], "robust_unit_regex": rule["unit_regex"],
        "robust_group_csv": rule["csv_path"], "robust_control_group": rule["control_group"],
        "robust_unit_aggregation": plan["aggregation"]
    })
    write_settings_properties(os.path.join(run_dir, "settings.properties"), recorded, all_configs)
    assignment_rows = []
    for row in assignment:
        cells = dict(base)
        cells.update({
            "Image_File": to_text(row["relative"]).replace("\\", "/").split("/")[-1],
            "Image_Relative_Path": row["relative"], "Group": row["group"], "Unit": row["unit"],
            "Status": row["status"], "Message": row["message"], "Rule": robust_rule_text(rule)
        })
        assignment_rows.append(cells)
    write_csv(os.path.join(run_dir, "group_assignment.csv"), ROBUST_ASSIGNMENT_HEADERS, assignment_rows)
    write_csv(os.path.join(run_dir, "robustness_specification_parameters.csv"),
              ROBUST_PARAMETER_HEADERS, robust_parameter_rows(plan, base))

    run = {
        "configs": configs, "global_options": global_options, "policy": policy, "roi_dir": roi_dir,
        "manager_rois": manager_rois, "log": run_log, "base": base, "records": [], "measurements": [],
        "units": set(), "calibration": calibration_override_from_options(global_options),
        "all_frames": to_text(global_options.get("frame_mode")) == "All frames"
    }
    counts = {"ok": 0, "failed": 0, "skipped": 0}
    not_assigned = len([row for row in assignment[:len(files)] if row["status"] != "ASSIGNED"])
    cancelled = False
    tables = None
    start = time.time()
    run_log("Robustness " + run_id + ": " + str(len(items)) + " assigned image file(s) of " + str(len(files)) +
            " from " + input_dir)
    run_log(APP_NAME + " " + APP_VERSION + " on ImageJ " + IMAGEJ_VERSION)
    run_log("C" + str(plan["config"]["channel"]) + " (" + plan["config"]["name"] + "), " + plan["metric"] + ", " +
            str(len(plan["specs"])) + " specifications: " + plan["grid_text"])
    run_log("Groups: " + ", ".join(group_names) + " by " + robust_rule_text(rule))
    run_log("Units are combined as: " + plan["aggregation"] + ". Parent ROIs: " + policy +
            (" (ROI folder " + roi_dir + ")" if roi_dir != "" else "") + ". Frames: " +
            ("all" if run["all_frames"] else "the first") + ".")
    for row in assignment:
        if row["status"] != "ASSIGNED":
            run_log("  " + row["status"] + " - " + to_text(row["relative"]) + ": " + row["message"])
    labels = {}
    for item in items:
        labels.setdefault(item["unit"], set()).add(item["group"])
    for label in sorted(labels.keys()):
        if len(labels[label]) > 1:
            run_log("  Note: unit '" + label + "' occurs in more than one group and is treated as one unit per group.")

    measurement_out = CsvAppender(os.path.join(run_dir, "robustness_measurements.csv"), ROBUST_MEASUREMENT_HEADERS)
    report_out = CsvAppender(os.path.join(run_dir, "robustness_image_report.csv"), ROBUST_REPORT_HEADERS)
    try:
        for file_index in range(len(items)):
            if is_cancelled is not None and is_cancelled():
                cancelled = True
                break
            item = items[file_index]
            path = item["path"]
            if progress is not None:
                progress(file_index, len(items), os.path.basename(path))
            run_log("Image " + str(file_index + 1) + " of " + str(len(items)) + ": " + item["relative"] +
                    " (" + item["group"] + ", unit " + item["unit"] + ")")
            # Files are examined and opened as a batch does: one series of
            # a multi-series file at a time, each released when it is done.
            series = []
            open_error = ""
            open_notes = []
            file_start = time.time()

            def open_log(message, _notes=open_notes):
                _notes.append(to_text(message).strip())
                run_log(message)

            try:
                series = plan_image_series(path, force_bioformats, all_series, open_log, skip_thumbnails)
            except ANY_ERROR as exc:
                open_error = to_text(exc)
            report_base = dict(base)
            report_base.update({
                "Image_Index": file_index + 1, "Image_File": os.path.basename(path),
                "Image_Relative_Path": item["relative"], "Group": item["group"], "Unit": item["unit"]
            })
            if len(series) == 0:
                counts["failed"] += 1
                if open_error == "" and len(open_notes) > 0:
                    open_error = open_notes[-1]
                report_base["Status"] = "FAILED"
                report_base["Message"] = "Could not open the file" + (": " + open_error if open_error else ".")
                report_base["Seconds"] = round(time.time() - file_start, 2)
                run_log("  FAILED - " + report_base["Message"])
                report_out.append([report_base])
                continue
            for entry in series:
                series_label = entry.label
                report = dict(report_base)
                report["Series"] = series_label
                if entry.skip_reason != "":
                    counts["skipped"] += 1
                    report["Status"] = "SKIPPED"
                    report["Message"] = entry.skip_reason
                    run_log("  SKIPPED - " + series_label + ": " + entry.skip_reason)
                    report_out.append([report])
                    continue
                # The time of an image includes opening it. A single image
                # was opened while the file was examined.
                image_start = file_start if entry.image is not None else time.time()
                try:
                    if is_cancelled is not None and is_cancelled():
                        raise AnalysisCancelled("Cancelled by user.")
                    imp = entry.open()
                    if imp is None:
                        raise ValueError("Bio-Formats returned no image for " + (series_label or "this file") + ".")
                    engine = AnalysisEngine(log=run_log, is_cancelled=is_cancelled)
                    report.update(robust_process_image(engine, imp, path, item, plan, run, series_label))
                    if report["Status"] == "OK":
                        counts["ok"] += 1
                        measurement_out.append(run["measurements"].pop())
                        run_log("  OK - " + str(report["Parent_ROI_Count"]) + " ROI(s), " + to_text(report["ROI_Source"]))
                    else:
                        counts["skipped"] += 1
                        run_log("  SKIPPED - " + to_text(report["Message"]))
                except AnalysisCancelled:
                    cancelled = True
                    report["Status"] = "CANCELLED"
                    report["Message"] = "Cancelled by user."
                except ANY_ERROR as exc:
                    counts["failed"] += 1
                    report["Status"] = "FAILED"
                    report["Message"] = to_text(exc)
                    run_log("  FAILED - " + to_text(exc))
                    try:
                        run_log(traceback.format_exc())
                    except ANY_ERROR:
                        pass
                finally:
                    report["Seconds"] = round(time.time() - image_start, 2)
                    report_out.append([report])
                    # Close the image and give its pixels back before the
                    # next one is opened.
                    entry.release()
                if cancelled:
                    break
            # Close whatever of this file was opened but not reached.
            for entry in series:
                entry.release()
            if cancelled:
                break

        if not cancelled:
            metric_unit = robust_metric_unit(plan["metric"], run["units"])
            tables = robust_statistics_tables(plan, run["records"], group_names, base, metric_unit)
            summary = tables["summary"]
            summary.update({
                "N_Images_Analyzed": counts["ok"], "N_Images_Failed": counts["failed"],
                "N_Images_Skipped": counts["skipped"], "N_Images_Not_Assigned": not_assigned
            })
            if len(run["units"]) > 1 and plan["metric"] != ROBUST_METRIC_AREA_FRACTION:
                summary["Note"] = (to_text(summary.get("Note", "")) + " The images have different spatial units (" +
                                   ", ".join(sorted(run["units"])) + "); their values are not comparable.").strip()
            write_csv(os.path.join(run_dir, "robustness_units.csv"), ROBUST_UNIT_HEADERS, tables["units"])
            write_csv(os.path.join(run_dir, "robustness_groups.csv"), ROBUST_GROUP_HEADERS, tables["groups"])
            write_csv(os.path.join(run_dir, "robustness_specifications.csv"),
                      ROBUST_SPECIFICATION_HEADERS, tables["specifications"])
            write_csv(os.path.join(run_dir, "robustness_parameter_values.csv"), ROBUST_VALUE_HEADERS, tables["values"])
            write_csv(os.path.join(run_dir, "robustness_summary.csv"), ROBUST_SUMMARY_HEADERS, [summary])
            if to_text(summary.get("Fraction_Same_Sign_As_Reference", "")) != "":
                run_log(
                    "Direction of the difference (" + group_names[1] + " minus " + group_names[0] + ") is the same "
                    "as with the configured settings in " + str(summary["N_Same_Sign_As_Reference"]) + " of " +
                    str(summary["N_Specifications_Compared"]) + " specifications. Reversed with: " +
                    to_text(summary["Sign_Flipping_Values"]) + "."
                )
            if to_text(summary.get("Note", "")) != "":
                run_log("Note: " + to_text(summary["Note"]))
            run_log(ROBUST_STATISTICS_NOTE)
    finally:
        measurement_out.close()
        report_out.close()
        duration = time.time() - start
        if cancelled:
            run_log("Cancelled by user. The measurements made so far are kept; no statistics were written.")
        run_log(
            "Robustness analysis finished in " + IJ.d2s(duration, 1) + " s: " + str(counts["ok"]) + " analyzed, " +
            str(counts["skipped"]) + " skipped, " + str(counts["failed"]) + " failed, " + str(not_assigned) +
            " without a group."
        )
        run_log("Saved: " + to_text(run_dir))
        write_text_file(os.path.join(run_dir, "robustness_log.txt"), [APP_NAME + " " + APP_VERSION] + log_lines)
        if progress is not None:
            progress(len(items), len(items), "")

    return {
        "run_dir": run_dir, "run_id": run_id, "files": len(files), "assigned": len(items),
        "not_assigned": not_assigned, "analyzed": counts["ok"], "skipped": counts["skipped"],
        "failed": counts["failed"], "cancelled": cancelled, "duration": duration,
        "specifications": len(plan["specs"]), "groups": group_names, "tables": tables,
        "summary": tables["summary"] if tables is not None else None
    }


def robust_request_from_options(global_options):
    """The arguments of run_robustness as stored in the settings (a profile
    or the window). A missing key uses its default."""
    options = global_options if global_options is not None else {}
    raw = {}
    for spec in ROBUST_SETTING_SPECS:
        value = options.get(spec[1])
        raw[spec[1]] = spec[3] if value is None else to_text(value)
    return {
        "channel": parse_int_value(raw["robust_channel"], "Robustness channel", 1, MAX_CHANNELS),
        "metric": raw["robust_metric"].strip(),
        "grid": parse_robust_grid(raw["robust_grid"]),
        "max_specifications": parse_int_value(raw["robust_max_specifications"], "Largest number of specifications", 2),
        "unit_aggregation": raw["robust_unit_aggregation"].strip(),
        "groups": {
            "source": raw["robust_group_source"].strip(), "separator": raw["robust_separator"],
            "group_regex": raw["robust_group_regex"], "unit_regex": raw["robust_unit_regex"],
            "csv_path": raw["robust_group_csv"], "control_group": raw["robust_control_group"]
        }
    }


def run_robustness_from_options(input_dir, output_root, configs, all_configs, global_options, batch_options,
                                manager_rois=None, log=None, progress=None, is_cancelled=None, file_list=None):
    """run_robustness with its own arguments read from the settings."""
    request = robust_request_from_options(global_options)
    return run_robustness(
        input_dir, output_root, configs, all_configs, global_options, batch_options,
        request["channel"], request["grid"], request["groups"], request["metric"],
        request["unit_aggregation"], request["max_specifications"],
        manager_rois=manager_rois, log=log, progress=progress, is_cancelled=is_cancelled, file_list=file_list
    )


register_global_settings(ROBUST_SETTING_SPECS)
register_settings_form(
    "Robustness to settings",
    "Repeats the measurement of one channel over a grid of settings and reports whether the difference "
    "between two groups of images keeps its direction. Grid example: " + ROBUST_DEFAULT_GRID +
    "; min_size=10,20; watershed=true,false",
    [
        ("robust_channel", "Channel whose settings are varied", "text", None),
        ("robust_metric", "Measurement to compare", "choice", list(ROBUST_METRICS)),
        ("robust_grid", "Grid of settings (setting=value,value; ...)", "text", None),
        ("robust_max_specifications", "Largest number of specifications", "text", None),
        ("robust_group_source", "Groups come from", "choice", list(ROBUST_GROUP_SOURCES)),
        ("robust_separator", "Separator that ends the group name", "text", None),
        ("robust_group_regex", "Group expression (one pair of parentheses)", "text", None),
        ("robust_unit_regex", "Unit expression, for example the organoid (empty = each image file)", "text", None),
        ("robust_group_csv", "Group table (CSV with columns Image, Group, Unit)", "text", None),
        ("robust_control_group", "Control group (empty = first in alphabetical order)", "text", None),
        ("robust_unit_aggregation", "Combine the measurements of a unit as", "choice", list(ROBUST_UNIT_AGGREGATIONS))
    ]
)
register_engine_module("cellquantrobust.py", "2.3.0")

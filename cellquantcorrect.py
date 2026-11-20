# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach
# Copyright (C) 2026 Maksym Taran

# Cell Quant Studio Universal - error-corrected object counts
#
# Automated counts are biased: touching nuclei that are detected as one object
# are undercounted, debris is overcounted. This module uses validated fields,
# where both the detected count and a manual count are known, to fit a linear
# calibration of manual against detected counts, and applies it to the object
# counts of a result table. Every corrected count comes with a 95 percent
# prediction interval.
#
# It is a calculation on tables only. It opens no image and no window, changes
# no Fiji option, and never modifies an existing file: the results are three
# new files (count_calibration.csv, count_calibration_fit.csv and
# <table name>_corrected.csv).
#
# The statistics, in short (the README has the full description):
#   - one calibration pair per validated parent ROI: (detected D, manual M)
#   - M is regressed on D, because D is what will be observed later and M is
#     what has to be predicted (direct prediction, also called inverse
#     calibration)
#   - two models are always fitted and reported: M = b D (through the origin,
#     the default) and M = a + b D
#   - the interval of a corrected count is the prediction interval of a new
#     field, with the t distribution on n - 1 or n - 2 degrees of freedom
#   - the leave-one-out error of the corrected counts is compared with the
#     error of the raw counts, and the correction is recommended only when it
#     is smaller
#
# This file is loaded into the namespace of cellquantengine.py by
# load_engine_modules, so every engine name used here is a global.
# The source is kept plain ASCII with no encoding declaration on purpose.

import os
import math
import re
import hashlib

# Student t quantiles come from Apache Commons Math, which ships with Fiji.
# Without it a built-in table and a series approximation are used instead.
try:
    from org.apache.commons.math3.distribution import TDistribution as CorrectTDistribution
except ANY_ERROR:
    CorrectTDistribution = None


CORRECT_MODEL_ORIGIN = "Proportional (through the origin)"
CORRECT_MODEL_LINE = "Line with intercept"
CORRECT_MODELS = [CORRECT_MODEL_ORIGIN, CORRECT_MODEL_LINE]

CORRECT_WEIGHT_NONE = "None (ordinary least squares)"
CORRECT_WEIGHT_COUNT = "1 / detected count"
CORRECT_WEIGHTINGS = [CORRECT_WEIGHT_NONE, CORRECT_WEIGHT_COUNT]

CORRECT_MIN_POINTS = 3            # fewer validated fields: refuse to fit
CORRECT_FEW_POINTS = 5            # fewer: strong warning
CORRECT_RECOMMENDED_POINTS = 8    # fewer: mild warning
CORRECT_RECOMMENDED_IMAGES = 3    # fields should come from at least this many images
CORRECT_POOLED_ROW = "ALL_ROIS_POOLED"
CORRECT_BY_EYE_MODE = "Display minimum"   # threshold set by hand per image
CORRECT_NAN = float("nan")
CORRECT_PATH_LIMIT = 259          # longest path that Windows programs open

COUNT_CALIBRATION_HEADERS = [
    "Run_ID", "Application_Version", "Calibration_ID", "Pair_Index", "Channel", "Channel_Name",
    "Source_Image", "Frame", "Parent_ROI_Name", "Validation_Run_ID", "Validation_Application_Version",
    "Calibration_Source_File", "Detected_Count", "Manual_Count", "Difference_Detected_Minus_Manual",
    "Percent_Error", "Model", "Weighting", "Fitted_Manual_Count", "Residual_Manual_Minus_Fitted",
    "LeaveOneOut_Corrected_Count", "LeaveOneOut_Error", "Settings_Fingerprint"
]
COUNT_CALIBRATION_FIT_HEADERS = [
    "Run_ID", "Application_Version", "Calibration_ID", "Channel", "Channel_Name", "Model", "Selected",
    "Weighting", "Fit_Status", "N_Fields", "N_Images", "Detected_Min", "Detected_Max",
    "Slope", "Slope_SE", "Slope_Lower95", "Slope_Upper95",
    "Intercept", "Intercept_SE", "Intercept_Lower95", "Intercept_Upper95", "Intercept_CI_Includes_Zero",
    "Degrees_Of_Freedom", "T_Quantile_975", "T_Quantile_Source", "Residual_SD",
    "R_Squared", "R_Squared_Uncentered", "AICc",
    "Mean_Difference_Detected_Minus_Manual", "SD_Difference",
    "Limits_Of_Agreement_Lower", "Limits_Of_Agreement_Upper",
    "MAPE_Percent", "MAPE_Fields_Used", "Raw_RMSE", "Raw_MAE",
    "LeaveOneOut_RMSE_Corrected", "LeaveOneOut_MAE_Corrected", "Correction_Recommended", "Recommendation",
    "Settings_Fingerprint", "Settings_Fingerprint_Source", "Validation_Application_Versions",
    "Corrected_Table", "Corrected_Count_Column", "Table_Settings_Fingerprint", "Settings_Match",
    "Rows_Corrected", "Rows_Extrapolated", "Rows_Clipped", "Warnings", "Interval_Note"
]

CORRECT_INTERVAL_NOTE = (
    "Lower95 and Upper95 are a 95 percent prediction interval for the manual count of one new field with "
    "this detected count. It covers counting error only, not biological or sampling variability, and holds "
    "only for the settings profile, stain and imaging conditions of the calibration fields."
)

# Channel settings that thresholds.csv also records, as (settings key, column).
# They allow a comparison of settings when a folder has no settings.properties.
CORRECT_THRESHOLD_COLUMNS = [
    ("threshold_mode", "Threshold_Mode"), ("auto_method", "Automatic_Method"),
    ("dark_background", "Dark_Background"), ("auto_scale", "Automatic_Multiplier"),
    ("local_radius", "Local_Radius"), ("local_parameter_1", "Local_Parameter_1"),
    ("local_parameter_2", "Local_Parameter_2"), ("clahe", "CLAHE_Enabled"),
    ("clahe_block_size", "CLAHE_Block_Size"), ("clahe_histogram_bins", "CLAHE_Histogram_Bins"),
    ("clahe_max_slope", "CLAHE_Maximum_Slope"), ("open_iterations", "Open_Iterations"),
    ("close_iterations", "Close_Iterations"), ("erode_iterations", "Erode_Iterations"),
    ("dilate_iterations", "Dilate_Iterations"), ("fill_holes", "Fill_Holes"),
    ("watershed", "Watershed"), ("size_units", "Object_Size_Units"),
    ("min_size", "Minimum_Object_Area"), ("max_size", "Maximum_Object_Area"),
    ("min_circularity", "Minimum_Circularity"), ("max_circularity", "Maximum_Circularity"),
    ("exclude_edges", "Exclude_Edges"), ("include_holes", "Include_Holes"),
    ("composite_rois", "Composite_ROIs"), ("four_connected", "Four_Connected"),
    ("morph_count", "Binary_Neighbor_Count")
]


# -----------------------------------------------------------------------------
# Reading CSV files (Jython 2.7 and Python 3, UTF-8, quoted or plain fields)
# -----------------------------------------------------------------------------

def correct_csv_delimiter(text):
    """Comma unless the header line clearly uses semicolons or tabs, as
    spreadsheet programs write in some regional settings."""
    first = to_text(text).split(u"\n", 1)[0]
    plain = re.sub(u'"[^"]*"', u"", first)
    best = u","
    for candidate in [u";", u"\t"]:
        if plain.count(candidate) > plain.count(best):
            best = candidate
    return best


def correct_parse_line(line, delimiter, carry=None):
    """Split one line of a CSV text into fields.

    carry is None at the start of a record. When a quoted field is still
    open at the end of the line (it holds a line break), the result is
    (None, carry) and the next line must be passed with that carry;
    otherwise it is (fields, None). A Windows line end is removed, except
    inside a quoted field, where the text is kept as it was."""
    quote = u'"'
    ending = u""
    if line.endswith(u"\r"):
        line = line[:-1]
        ending = u"\r"
    size = len(line)
    pos = 0
    if carry is None:
        fields = []
        parts = None
    else:
        fields = carry[0]
        parts = carry[1]
        parts.append(u"\n")
    while True:
        value = u""
        if parts is None and pos < size and line[pos] == quote:
            parts = []
            pos += 1
        if parts is not None:
            while True:
                end = line.find(quote, pos)
                if end < 0:
                    parts.append(line[pos:] + ending)
                    return (None, (fields, parts))
                parts.append(line[pos:end])
                if end + 1 < size and line[end + 1] == quote:
                    parts.append(quote)
                    pos = end + 2
                else:
                    pos = end + 1
                    break
            value = u"".join(parts)
            parts = None
        stop = line.find(delimiter, pos)
        if stop < 0:
            stop = size
        value += line[pos:stop]
        fields.append(value)
        if stop >= size:
            return (fields, None)
        pos = stop + 1


def correct_parse_csv(text, delimiter=None):
    """Split CSV text into records, each a list of text fields.

    Understands what write_csv produces (every field quoted, quotes doubled)
    and also plain unquoted fields, line breaks inside quoted fields, a
    byte-order mark, and Windows or Unix line ends. Blank lines are dropped.

    The text is split into lines first and each line is parsed on its own:
    under Jython, searching a long text from a position becomes very slow
    once the text is larger than a few hundred kilobytes (a 750 kB table took
    47 s when it was searched as one text)."""
    text = to_text(text)
    bom = b"\xef\xbb\xbf".decode("utf-8")
    if text.startswith(bom):
        text = text[len(bom):]
    if delimiter is None:
        delimiter = correct_csv_delimiter(text)
    if u"\n" not in text and u"\r" in text:
        text = text.replace(u"\r", u"\n")    # old Mac line ends
    records = []
    carry = None
    for line in text.split(u"\n"):
        fields, carry = correct_parse_line(line, delimiter, carry)
        if fields is not None:
            records.append(fields)
    if carry is not None:
        # A quote that is never closed: the field runs to the end of the text.
        carry[0].append(u"".join(carry[1]))
        records.append(carry[0])
    return [item for item in records if not (len(item) == 1 and item[0].strip() == u"")]


def correct_read_csv(path):
    """Read a CSV file into {"headers": [...], "rows": [[...], ...]}.
    Every value is text. Short rows are padded with empty fields."""
    handle = open(path, "rb")
    try:
        data = handle.read()
    finally:
        handle.close()
    text = to_text(data)
    delimiter = correct_csv_delimiter(text)
    records = correct_parse_csv(text, delimiter)
    if len(records) == 0:
        raise ValueError("The table is empty: " + to_text(path))
    headers = [to_text(name) for name in records[0]]
    rows = []
    for index in range(1, len(records)):
        row = records[index]
        if len(row) > len(headers):
            extra = [value for value in row[len(headers):] if value.strip() != u""]
            if len(extra) > 0:
                raise ValueError(
                    "Row " + str(index + 1) + " of " + to_text(os.path.basename(path)) +
                    " has more fields than the header line."
                )
            row = row[:len(headers)]
        while len(row) < len(headers):
            row.append(u"")
        rows.append(row)
    return {"headers": headers, "rows": rows, "delimiter": delimiter, "path": path}


def correct_number(value, decimal_comma=False):
    """A finite number, or None when the field is empty or not a number."""
    raw = to_text(value).strip()
    if raw == u"":
        return None
    number = None
    try:
        number = float(raw)
    except ANY_ERROR:
        if decimal_comma and raw.count(u",") == 1 and u"." not in raw:
            try:
                number = float(raw.replace(u",", u"."))
            except ANY_ERROR:
                number = None
    if number is None or number != number or math.isinf(number):
        return None
    return number


def correct_header_lookup(headers):
    """Map a simplified header name (lower case, underscores) to its index."""
    lookup = {}
    for index in range(len(headers)):
        key = to_text(headers[index]).strip().lower().replace(u" ", u"_")
        if key not in lookup:
            lookup[key] = index
    return lookup


def correct_table_value(row, lookup, name):
    index = lookup.get(name.lower())
    if index is None:
        return u""
    return to_text(row[index])


# -----------------------------------------------------------------------------
# Student t quantile
# -----------------------------------------------------------------------------

CORRECT_T975_TABLE = [
    12.706204736, 4.302652730, 3.182446305, 2.776445105, 2.570581836,
    2.446911851, 2.364624252, 2.306004135, 2.262157163, 2.228138852,
    2.200985160, 2.178812830, 2.160368656, 2.144786688, 2.131449546,
    2.119905299, 2.109815578, 2.100922040, 2.093024054, 2.085963447,
    2.079613845, 2.073873068, 2.068657610, 2.063898562, 2.059538553,
    2.055529439, 2.051830516, 2.048407142, 2.045229642, 2.042272456
]


def correct_t975_builtin(df):
    """The 97.5 percent point of Student's t without any library: a table up
    to 30 degrees of freedom, then the Cornish-Fisher series around the
    normal quantile (accurate to better than 0.000001 from 31 upward)."""
    df = int(df)
    if df < 1:
        return CORRECT_NAN
    if df <= len(CORRECT_T975_TABLE):
        return CORRECT_T975_TABLE[df - 1]
    z = 1.959963984540054
    n = float(df)
    g1 = (z ** 3 + z) / 4.0
    g2 = (5.0 * z ** 5 + 16.0 * z ** 3 + 3.0 * z) / 96.0
    g3 = (3.0 * z ** 7 + 19.0 * z ** 5 + 17.0 * z ** 3 - 15.0 * z) / 384.0
    g4 = (79.0 * z ** 9 + 776.0 * z ** 7 + 1482.0 * z ** 5 - 1920.0 * z ** 3 - 945.0 * z) / 92160.0
    return z + g1 / n + g2 / n ** 2 + g3 / n ** 3 + g4 / n ** 4


def correct_t975(df):
    """(two-sided 95 percent t multiplier, where the value came from)."""
    df = int(df)
    if df < 1:
        return (CORRECT_NAN, "not defined")
    if CorrectTDistribution is not None:
        try:
            value = float(CorrectTDistribution(float(df)).inverseCumulativeProbability(0.975))
            return (value, "Apache Commons Math TDistribution")
        except ANY_ERROR:
            pass
    return (correct_t975_builtin(df), "built-in table and series approximation")


# -----------------------------------------------------------------------------
# Fitting: manual count regressed on detected count
# -----------------------------------------------------------------------------

def correct_weights(detected, weighting):
    if weighting == CORRECT_WEIGHT_NONE:
        return [1.0] * len(detected)
    if weighting == CORRECT_WEIGHT_COUNT:
        for value in detected:
            if value <= 0:
                raise ValueError(
                    "The weighting '" + CORRECT_WEIGHT_COUNT + "' needs at least one detected object in "
                    "every calibration field. Remove the empty field or use '" + CORRECT_WEIGHT_NONE + "'."
                )
        return [1.0 / value for value in detected]
    raise ValueError("Unknown calibration weighting: " + to_text(weighting))


def correct_coefficients(detected, manual, model, weighting):
    """Weighted least-squares coefficients of manual = intercept + slope *
    detected, without any inference. Raises ValueError when the data cannot
    determine the model."""
    weights = correct_weights(detected, weighting)
    count = len(detected)
    sum_w = sum(weights)
    if model == CORRECT_MODEL_ORIGIN:
        sdd = sum([weights[i] * detected[i] * detected[i] for i in range(count)])
        sdm = sum([weights[i] * detected[i] * manual[i] for i in range(count)])
        if count < 1 or sdd <= 0:
            raise ValueError("Every calibration field has zero detected objects, so no correction factor can be fitted.")
        return {
            "intercept": 0.0, "slope": sdm / sdd, "weights": weights, "sum_w": sum_w,
            "d_mean": 0.0, "spread": sdd, "parameters": 1
        }
    if model == CORRECT_MODEL_LINE:
        if count < 2:
            raise ValueError("A line needs at least two calibration fields.")
        d_mean = sum([weights[i] * detected[i] for i in range(count)]) / sum_w
        m_mean = sum([weights[i] * manual[i] for i in range(count)]) / sum_w
        sxx = sum([weights[i] * (detected[i] - d_mean) ** 2 for i in range(count)])
        sxy = sum([weights[i] * (detected[i] - d_mean) * (manual[i] - m_mean) for i in range(count)])
        if sxx <= 1e-12 * sum_w * max(1.0, d_mean * d_mean):
            raise ValueError(
                "All calibration fields have the same detected count, so a line with an intercept cannot be fitted."
            )
        slope = sxy / sxx
        return {
            "intercept": m_mean - slope * d_mean, "slope": slope, "weights": weights, "sum_w": sum_w,
            "d_mean": d_mean, "spread": sxx, "parameters": 2
        }
    raise ValueError("Unknown calibration model: " + to_text(model))


def correct_checked_counts(values, label):
    result = []
    for value in values:
        number = float(value)
        if number != number or math.isinf(number) or number < 0:
            raise ValueError(label + " counts must be finite numbers that are not negative.")
        result.append(number)
    return result


def fit_count_calibration(detected, manual, model=CORRECT_MODEL_ORIGIN, weighting=CORRECT_WEIGHT_NONE):
    """Fit manual = intercept + slope * detected and return the coefficients
    with their standard errors, 95 percent confidence intervals, the residual
    standard deviation, R squared and AICc.

    model: CORRECT_MODEL_ORIGIN (intercept fixed at zero) or CORRECT_MODEL_LINE.
    weighting: CORRECT_WEIGHT_NONE (constant error variance) or
               CORRECT_WEIGHT_COUNT (error variance proportional to the
               detected count; weights 1 / detected).
    Raises ValueError with fewer than three fields or when the model cannot
    be determined."""
    detected = correct_checked_counts(detected, "Detected")
    manual = correct_checked_counts(manual, "Manual")
    count = len(detected)
    if len(manual) != count:
        raise ValueError("The detected and manual counts are not paired one to one.")
    if count < CORRECT_MIN_POINTS:
        raise ValueError(
            "A count calibration needs at least " + str(CORRECT_MIN_POINTS) + " validated fields; " +
            str(count) + " found. Validate more fields (8 or more are recommended)."
        )
    core = correct_coefficients(detected, manual, model, weighting)
    weights = core["weights"]
    slope = core["slope"]
    intercept = core["intercept"]
    df = count - core["parameters"]
    fitted = [intercept + slope * value for value in detected]
    residuals = [manual[i] - fitted[i] for i in range(count)]
    sse = sum([weights[i] * residuals[i] ** 2 for i in range(count)])
    residual_sd = math.sqrt(sse / float(df))
    t_value, t_source = correct_t975(df)

    slope_se = residual_sd / math.sqrt(core["spread"])
    intercept_se = CORRECT_NAN
    intercept_low = CORRECT_NAN
    intercept_high = CORRECT_NAN
    if model == CORRECT_MODEL_LINE:
        intercept_se = residual_sd * math.sqrt(1.0 / core["sum_w"] + core["d_mean"] ** 2 / core["spread"])
        intercept_low = intercept - t_value * intercept_se
        intercept_high = intercept + t_value * intercept_se

    # R squared is 1 - SSE / SS(manual about its mean) for both models, so
    # the two can be compared. The uncentered form (1 - SSE / sum of squared
    # manual counts) is what most programs print for a fit through the
    # origin; it is reported separately for that model.
    m_mean = sum([weights[i] * manual[i] for i in range(count)]) / core["sum_w"]
    syy = sum([weights[i] * (manual[i] - m_mean) ** 2 for i in range(count)])
    smm = sum([weights[i] * manual[i] ** 2 for i in range(count)])
    r_squared = 1.0 - sse / syy if syy > 0 else CORRECT_NAN
    r_uncentered = CORRECT_NAN
    if model == CORRECT_MODEL_ORIGIN and smm > 0:
        r_uncentered = 1.0 - sse / smm

    # AICc with k = coefficients + 1 (the error variance). Only differences
    # between the two models, fitted with the same weighting, mean anything.
    k = core["parameters"] + 1
    aicc = CORRECT_NAN
    if sse > 1e-12 * max(smm, 1e-300) and count - k - 1 > 0:
        aicc = count * math.log(sse / count) + 2.0 * k + 2.0 * k * (k + 1) / float(count - k - 1)

    return {
        "model": model, "weighting": weighting, "n": count, "df": df,
        "slope": slope, "intercept": intercept,
        "slope_se": slope_se, "slope_low": slope - t_value * slope_se, "slope_high": slope + t_value * slope_se,
        "intercept_se": intercept_se, "intercept_low": intercept_low, "intercept_high": intercept_high,
        "t": t_value, "t_source": t_source, "residual_sd": residual_sd, "sse": sse,
        "r_squared": r_squared, "r_squared_uncentered": r_uncentered, "aicc": aicc,
        "d_min": min(detected), "d_max": max(detected),
        "sum_w": core["sum_w"], "d_mean": core["d_mean"], "spread": core["spread"],
        "fitted": fitted, "residuals": residuals
    }


def predict_corrected_count(fit, detected):
    """Corrected count for one detected count, with its 95 percent prediction
    interval: the range expected to hold the manual count of ONE new field
    with this detected count. It is wider than the confidence interval of
    the fitted line because it includes the scatter of single fields.

    Returns None when detected is not a usable count. A corrected value or
    limit below zero is set to zero; "clipped" tells when the corrected value
    itself was negative."""
    if detected is None:
        return None
    d0 = float(detected)
    if d0 != d0 or math.isinf(d0) or d0 < 0:
        return None
    estimate = fit["intercept"] + fit["slope"] * d0
    # Variance of a new observation in units of the residual variance:
    # 1 for constant variance, D for variance proportional to the count.
    unit = d0 if fit["weighting"] == CORRECT_WEIGHT_COUNT else 1.0
    if fit["model"] == CORRECT_MODEL_LINE:
        unit += 1.0 / fit["sum_w"]
    unit += (d0 - fit["d_mean"]) ** 2 / fit["spread"]
    half = fit["t"] * fit["residual_sd"] * math.sqrt(unit)
    clipped = estimate < 0
    return {
        "detected": d0,
        "estimate": max(0.0, estimate),
        "lower": max(0.0, estimate - half),
        "upper": max(0.0, estimate + half),
        "clipped": clipped,
        "extrapolated": d0 < fit["d_min"] or d0 > fit["d_max"]
    }


def count_agreement(detected, manual):
    """Agreement of the raw detected counts with the manual counts: mean
    difference (bias), Bland-Altman 95 percent limits of agreement
    (bias +/- 1.96 SD of the differences), mean absolute percentage error
    over the fields with a manual count above zero, RMSE and MAE."""
    count = len(detected)
    diffs = [float(detected[i]) - float(manual[i]) for i in range(count)]
    mean_diff = sum(diffs) / float(count) if count > 0 else CORRECT_NAN
    sd_diff = CORRECT_NAN
    if count > 1:
        sd_diff = math.sqrt(sum([(value - mean_diff) ** 2 for value in diffs]) / float(count - 1))
    percents = [abs(diffs[i]) / float(manual[i]) * 100.0 for i in range(count) if manual[i] > 0]
    return {
        "mean_difference": mean_diff,
        "sd_difference": sd_diff,
        "loa_low": mean_diff - 1.96 * sd_diff,
        "loa_high": mean_diff + 1.96 * sd_diff,
        "mape": sum(percents) / float(len(percents)) if len(percents) > 0 else CORRECT_NAN,
        "mape_fields": len(percents),
        "rmse": math.sqrt(sum([value ** 2 for value in diffs]) / float(count)) if count > 0 else CORRECT_NAN,
        "mae": sum([abs(value) for value in diffs]) / float(count) if count > 0 else CORRECT_NAN
    }


def leave_one_out_correction(detected, manual, model, weighting):
    """Correct every calibration field with a fit that did not see it.
    Returns the corrected counts (None where the remaining fields cannot
    determine the model) with their RMSE and MAE against the manual counts."""
    count = len(detected)
    predictions = []
    errors = []
    for index in range(count):
        rest_d = [detected[i] for i in range(count) if i != index]
        rest_m = [manual[i] for i in range(count) if i != index]
        try:
            core = correct_coefficients(rest_d, rest_m, model, weighting)
            value = max(0.0, core["intercept"] + core["slope"] * detected[index])
            predictions.append(value)
            errors.append(value - manual[index])
        except ValueError:
            predictions.append(None)
    rmse = CORRECT_NAN
    mae = CORRECT_NAN
    if count > 0 and len(errors) == count:
        rmse = math.sqrt(sum([value ** 2 for value in errors]) / float(count))
        mae = sum([abs(value) for value in errors]) / float(count)
    return {"predictions": predictions, "rmse": rmse, "mae": mae}


# -----------------------------------------------------------------------------
# Settings fingerprint: was the table made with the calibrated settings?
# -----------------------------------------------------------------------------

def correct_setting_keys(config):
    """The channel settings that can change which objects are detected, given
    the modes that are switched on. Display settings are left out."""
    keys = [
        "projection", "subtract_background", "median_radius", "gaussian_sigma", "threshold_mode",
        "open_iterations", "close_iterations", "erode_iterations", "dilate_iterations",
        "fill_holes", "watershed", "size_units", "min_size", "max_size",
        "min_circularity", "max_circularity", "exclude_edges", "include_holes",
        "composite_rois", "four_connected"
    ]
    if to_text(config.get("projection")) == "Specific Z":
        keys.append("z_slice")
    if parse_bool(config.get("subtract_background"), False):
        keys += ["background_radius", "light_background"]
    mode = to_text(config.get("threshold_mode"))
    if mode == "Auto":
        keys += ["auto_method", "dark_background", "auto_scale"]
    elif mode == "Manual":
        keys += ["manual_low", "manual_high"]
    elif mode == "Local Phansalkar":
        keys += ["local_radius", "local_parameter_1", "local_parameter_2", "clahe"]
        if parse_bool(config.get("clahe"), False):
            keys += ["clahe_block_size", "clahe_histogram_bins", "clahe_max_slope"]
    for name in ["open_iterations", "close_iterations", "erode_iterations", "dilate_iterations"]:
        if (correct_number(config.get(name)) or 0) > 0:
            keys.append("morph_count")
            break
    return keys


def correct_setting_value(key, value):
    if key in CHANNEL_BOOLEAN_KEYS:
        return bool_text(parse_bool(value, False))
    if key in CHANNEL_NUMERIC_KEYS:
        number = correct_number(value)
        if number is not None:
            return repr(number)
    return to_text(value).strip()


def correct_settings_fingerprint(config, partial=False, extra=None):
    """A short code that is the same for two channel configurations exactly
    when their detection settings are the same. partial=True uses only the
    settings that thresholds.csv records."""
    recorded = set([pair[0] for pair in CORRECT_THRESHOLD_COLUMNS])
    parts = []
    for key in correct_setting_keys(config):
        if partial and key not in recorded:
            continue
        parts.append(key + u"=" + correct_setting_value(key, config.get(key, CHANNEL_DEFAULTS.get(key))))
    if extra and not partial:
        parts += [to_text(item) for item in extra]
    return hashlib.sha1(u";".join(parts).encode("utf-8")).hexdigest()[:12]


def correct_folder_fingerprints(folder, channel):
    """(full, partial, by eye) for one channel's settings in a run,
    validation or batch folder. The fingerprints are an empty text where
    not available: the full one needs settings.properties, the partial one
    can also be built from thresholds.csv or batch_thresholds.csv.
    "by eye" is True when the channel uses the 'Display minimum' threshold,
    whose value is set by hand for each image and is not part of a profile."""
    full = u""
    partial = u""
    by_eye = False
    try:
        channel = int(float(to_text(channel)))
    except ANY_ERROR:
        return (full, partial, by_eye)
    path = os.path.join(folder, "settings.properties")
    if os.path.isfile(path):
        try:
            props = read_properties_file(path)
            values = channel_values_from_properties(props, channel)
            if len(values) > 0:
                config = normalize_channel_config(channel, values)
                by_eye = config["threshold_mode"] == CORRECT_BY_EYE_MODE
                extra = []
                if parse_bool(props.getProperty("global.override_calibration"), False):
                    extra = [
                        u"pixel_width=" + correct_setting_value("min_size", props.getProperty("global.pixel_width")),
                        u"pixel_height=" + correct_setting_value("min_size", props.getProperty("global.pixel_height"))
                    ]
                full = correct_settings_fingerprint(config, False, extra)
                partial = correct_settings_fingerprint(config, True)
        except ANY_ERROR:
            full = u""
            partial = u""
    if partial == u"":
        for name in ["thresholds.csv", "batch_thresholds.csv"]:
            path = os.path.join(folder, name)
            if not os.path.isfile(path):
                continue
            try:
                table = correct_read_csv(path)
                lookup = correct_header_lookup(table["headers"])
                for row in table["rows"]:
                    if correct_number(correct_table_value(row, lookup, "Channel")) != channel:
                        continue
                    config = {}
                    for key, column in CORRECT_THRESHOLD_COLUMNS:
                        value = correct_table_value(row, lookup, column)
                        config[key] = value if value.strip() != u"" else CHANNEL_DEFAULTS.get(key)
                    partial = correct_settings_fingerprint(config, True)
                    by_eye = to_text(config.get("threshold_mode")).strip() == CORRECT_BY_EYE_MODE
                    break
            except ANY_ERROR:
                partial = u""
            if partial != u"":
                break
    return (full, partial, by_eye)


# -----------------------------------------------------------------------------
# Collecting calibration pairs
# -----------------------------------------------------------------------------

def correct_channel_number(value):
    """The channel number as text: '2' from 2, 2.0, '2', 'C2' or 'C2_DAPI'.
    An empty text when the value holds no channel number."""
    match = re.match(u"^[cC]?([0-9]+)(\\.0+)?(_.*)?$", to_text(value).strip())
    if match is None:
        return u""
    return to_text(str(int(match.group(1))))


def correct_channel_matches(wanted, channel, channel_name):
    """wanted may be a channel number (2), a label (C2), a channel name
    (DAPI), or a column prefix (C2_DAPI). An empty wanted matches all."""
    wanted = to_text(wanted).strip().lower()
    if wanted == u"":
        return True
    number = correct_channel_number(channel)
    name = to_text(channel_name).strip().lower()
    accepted = []
    if number != u"":
        accepted += [number, u"c" + number]
        if name != u"":
            accepted.append((u"c" + number + u"_" + sanitize_name(channel_name)).lower())
    if name != u"":
        accepted += [name, sanitize_name(channel_name).lower()]
    return wanted in accepted


def correct_find_validation_summaries(folder):
    """Every validation_summary.csv below a folder, in a fixed order."""
    found = []
    for root, dirs, files in os.walk(to_text(folder)):
        dirs.sort()
        for name in sorted(files):
            if name == "validation_summary.csv":
                found.append(os.path.join(root, name))
    return found


def correct_pairs_from_table(table, source_label, channel, folder, cache, skipped):
    """Calibration pairs from one table that has Manual_Count and
    Detected_Count columns. Other columns are used when present."""
    lookup = correct_header_lookup(table["headers"])
    if "manual_count" not in lookup or "detected_count" not in lookup:
        raise ValueError(
            to_text(source_label) + " has no Manual_Count and Detected_Count columns."
        )
    decimal_comma = table["delimiter"] != u","
    pairs = []
    for row in table["rows"]:
        if correct_table_value(row, lookup, "Parent_ROI_Name") == CORRECT_POOLED_ROW:
            continue
        row_channel = correct_channel_number(correct_table_value(row, lookup, "Channel"))
        row_name = correct_table_value(row, lookup, "Channel_Name")
        if (row_channel != u"" or row_name != u"") and not correct_channel_matches(channel, row_channel, row_name):
            continue
        manual = correct_number(correct_table_value(row, lookup, "Manual_Count"), decimal_comma)
        detected = correct_number(correct_table_value(row, lookup, "Detected_Count"), decimal_comma)
        if manual is None or detected is None or manual < 0 or detected < 0:
            skipped.append(to_text(source_label))
            continue
        full = u""
        partial = u""
        by_eye = False
        if folder is not None and row_channel != u"":
            key = row_channel
            if key not in cache:
                cache[key] = correct_folder_fingerprints(folder, row_channel)
            full, partial, by_eye = cache[key]
        pairs.append({
            "detected": detected, "manual": manual, "by_eye": by_eye,
            "channel": row_channel, "channel_name": row_name,
            "source_image": correct_table_value(row, lookup, "Source_Image"),
            "frame": correct_table_value(row, lookup, "Frame"),
            "parent_roi": correct_table_value(row, lookup, "Parent_ROI_Name"),
            "run_id": correct_table_value(row, lookup, "Run_ID"),
            "app_version": correct_table_value(row, lookup, "Application_Version"),
            "source_file": to_text(source_label),
            "fp_full": full, "fp_partial": partial
        })
    return pairs


def collect_calibration_pairs(folder=None, csv_path=None, channel=None):
    """Gather the calibration pairs of ONE channel.

    folder: scanned, with its subfolders, for validation_summary.csv files
            written by the object validation (the pooled row is skipped).
    csv_path: a table with Manual_Count and Detected_Count columns, one row
            per validated field (other columns are optional).
    channel: number, C-label, channel name or column prefix. May be left
            empty when the data hold a single channel.
    Returns {"pairs", "channel", "channel_name", "files", "skipped", "notes"};
    notes are warnings about the data, in plain language."""
    pairs = []
    skipped = []
    files = 0
    wanted = u"" if channel is None else to_text(channel)
    if folder is not None and to_text(folder).strip() != u"":
        folder = os.path.abspath(to_text(folder))
        if not os.path.isdir(folder):
            raise ValueError("The calibration folder does not exist: " + folder)
        for path in correct_find_validation_summaries(folder):
            label = os.path.relpath(path, folder).replace(os.sep, "/")
            pairs += correct_pairs_from_table(correct_read_csv(path), label, wanted, os.path.dirname(path), {}, skipped)
            files += 1
    if csv_path is not None and to_text(csv_path).strip() != u"":
        csv_path = to_text(csv_path)
        if not os.path.isfile(csv_path):
            raise ValueError("The calibration table does not exist: " + csv_path)
        pairs += correct_pairs_from_table(correct_read_csv(csv_path), os.path.basename(csv_path), wanted, None, {}, skipped)
        files += 1
    if files == 0:
        raise ValueError(
            "No calibration data: no validation_summary.csv was found in the folder and no calibration table was given."
        )
    numbers = correct_distinct(pairs, "channel")
    names = correct_distinct(pairs, "channel_name")
    if len(numbers) > 1 or (len(numbers) == 0 and len(names) > 1):
        found = u", ".join([u"C" + item for item in numbers] if len(numbers) > 1 else names)
        raise ValueError(
            "The calibration data hold more than one channel (" + found +
            "). A calibration is fitted for one channel at a time: choose the channel."
        )
    if len(pairs) == 0:
        raise ValueError(
            "No calibration pairs were found" + (u" for channel " + wanted if wanted != u"" else u"") +
            ". Each validated parent ROI gives one pair of a manual and a detected count."
        )
    notes = []
    if len(names) > 1:
        notes.append(
            "The calibration fields carry different channel names (" + u", ".join(names) +
            ") for the same channel number. A calibration is valid for one stain only: check the validation runs."
        )
    if len(skipped) > 0:
        notes.append(str(len(skipped)) + " calibration row(s) without a usable manual and detected count were skipped.")
    return {
        "pairs": pairs,
        "channel": numbers[0] if len(numbers) == 1 else u"",
        "channel_name": names[0] if len(names) == 1 else u"",
        "files": files, "skipped": skipped, "notes": notes
    }


def correct_pairs_from_counts(detected, manual):
    """Calibration pairs from two plain lists (for scripts and tests)."""
    if len(detected) != len(manual):
        raise ValueError("The detected and manual counts are not paired one to one.")
    return [{"detected": float(detected[i]), "manual": float(manual[i])} for i in range(len(detected))]


# -----------------------------------------------------------------------------
# The calibration: both fits, agreement, cross-validation, warnings
# -----------------------------------------------------------------------------

def correct_number_text(value):
    try:
        if value != value:
            return "not available"
        return "%.4g" % value
    except ANY_ERROR:
        return to_text(value)


def correct_distinct(pairs, key):
    return sorted(set([to_text(pair.get(key, u"")) for pair in pairs if to_text(pair.get(key, u"")) != u""]))


def correct_calibration_fingerprint(pairs):
    """(fingerprint, where it came from, True when the pairs disagree)."""
    full = correct_distinct(pairs, "fp_full")
    partial = correct_distinct(pairs, "fp_partial")
    with_full = len([pair for pair in pairs if to_text(pair.get("fp_full", u"")) != u""])
    with_partial = len([pair for pair in pairs if to_text(pair.get("fp_partial", u"")) != u""])
    if len(full) > 1 or len(partial) > 1:
        return (u"MIXED", "validation runs with different settings", True)
    if len(full) == 1 and with_full == len(pairs):
        return (full[0], "settings.properties", False)
    if len(partial) == 1 and with_partial == len(pairs):
        return (partial[0], "thresholds.csv (segmentation settings only)", False)
    return (u"", "not available", False)


def correct_calibration_warnings(pairs, fits, model, weighting):
    """Plain-language warnings about the calibration data and the model."""
    warnings = []
    count = len(pairs)
    if count < CORRECT_FEW_POINTS:
        warnings.append(
            "Only " + str(count) + " validated fields: the calibration is barely determined and the "
            "intervals are very wide. Validate at least " + str(CORRECT_RECOMMENDED_POINTS) + " fields."
        )
    elif count < CORRECT_RECOMMENDED_POINTS:
        warnings.append(
            "Only " + str(count) + " validated fields; " + str(CORRECT_RECOMMENDED_POINTS) +
            " or more, spanning sparse and dense regions, are recommended."
        )
    images = correct_distinct(pairs, "source_image")
    with_image = len([pair for pair in pairs if to_text(pair.get("source_image", u"")) != u""])
    if with_image < count:
        warnings.append("The number of images behind the calibration fields is not known (no Source_Image column).")
    elif len(images) < CORRECT_RECOMMENDED_IMAGES:
        warnings.append(
            "The calibration fields come from " + str(len(images)) + " image(s). Fields of one image are not "
            "independent; use at least " + str(CORRECT_RECOMMENDED_IMAGES) + " images or sections."
        )
    seen = set()
    repeated = 0
    for pair in pairs:
        key = (to_text(pair.get("source_image", u"")), to_text(pair.get("frame", u"")), to_text(pair.get("parent_roi", u"")))
        if key[0] == u"" and key[2] == u"":
            continue
        if key in seen:
            repeated += 1
        seen.add(key)
    if repeated > 0:
        warnings.append(
            str(repeated) + " calibration pair(s) repeat a field (same image, frame and parent ROI) that is "
            "already in the calibration. Repeated fields make the intervals too narrow; keep one validation per field."
        )
    fingerprint = correct_calibration_fingerprint(pairs)
    if fingerprint[2]:
        warnings.append(
            "The validation runs were made with DIFFERENT detection settings. A calibration is valid for one "
            "settings profile only: keep the validation runs of one profile in the calibration folder."
        )
    elif fingerprint[0] == u"":
        warnings.append(
            "The detection settings of the calibration fields are not recorded, so it cannot be checked that "
            "the corrected table was made with the same settings."
        )
    if len([pair for pair in pairs if pair.get("by_eye")]) > 0:
        warnings.append(
            "Calibration fields were thresholded with '" + CORRECT_BY_EYE_MODE + "', a value set by eye for each "
            "image. It is not part of the settings profile, so the fields may not share one threshold. Validate "
            "with a fixed threshold (Manual, Auto or Local Phansalkar) before correcting counts."
        )
    versions = correct_distinct(pairs, "app_version")
    if len(versions) > 1:
        warnings.append("The validation runs were made with different application versions: " + u", ".join(versions) + ".")
    line = fits.get(CORRECT_MODEL_LINE)
    if model == CORRECT_MODEL_LINE:
        warnings.append(
            "A line with an intercept adds the same number of objects to every ROI whatever its size. "
            "Use it only for ROIs of about the same area as the calibration fields."
        )
        if line["intercept_low"] <= 0 <= line["intercept_high"]:
            warnings.append(
                "The intercept cannot be told apart from zero (95 percent CI " +
                correct_number_text(line["intercept_low"]) + " to " + correct_number_text(line["intercept_high"]) +
                "), so the simpler proportional model is adequate for these data."
            )
    if model == CORRECT_MODEL_ORIGIN and line is not None and "error" not in line:
        if line["intercept_low"] > 0 or line["intercept_high"] < 0:
            warnings.append(
                "The intercept of the line model differs from zero (95 percent CI " +
                correct_number_text(line["intercept_low"]) + " to " + correct_number_text(line["intercept_high"]) +
                "), so a purely proportional correction may not describe these data. Compare both models in "
                "count_calibration_fit.csv."
            )
    for name in CORRECT_MODELS:
        if name != model and "error" in fits.get(name, {}):
            warnings.append("The model '" + name + "' could not be fitted for comparison: " + fits[name]["error"])
    return warnings


def build_count_calibration(pairs, model=CORRECT_MODEL_ORIGIN, weighting=CORRECT_WEIGHT_NONE,
                            channel=u"", channel_name=u""):
    """Fit both models to the calibration pairs, score them, and decide
    whether the selected model's correction is worth applying.

    pairs: dictionaries with "detected" and "manual" (see
           collect_calibration_pairs and correct_pairs_from_counts).
    Returns the calibration used by correct_count_table and
    write_count_calibration. Raises ValueError when the selected model
    cannot be fitted (fewer than three fields, no spread, unknown name)."""
    if model not in CORRECT_MODELS:
        raise ValueError("Unknown calibration model: " + to_text(model))
    if weighting not in CORRECT_WEIGHTINGS:
        raise ValueError("Unknown calibration weighting: " + to_text(weighting))
    detected = [float(pair["detected"]) for pair in pairs]
    manual = [float(pair["manual"]) for pair in pairs]
    fits = {}
    loo = {}
    for name in CORRECT_MODELS:
        try:
            fits[name] = fit_count_calibration(detected, manual, name, weighting)
            loo[name] = leave_one_out_correction(detected, manual, name, weighting)
        except ValueError as exc:
            if name == model:
                raise
            fits[name] = {"model": name, "error": to_text(exc)}
    agreement = count_agreement(detected, manual)
    warnings = correct_calibration_warnings(pairs, fits, model, weighting)

    # The correction is worth applying only if counts corrected by a fit
    # that never saw the field are closer to the manual counts than the raw
    # detected counts are.
    recommended = {}
    texts = {}
    for name in CORRECT_MODELS:
        if name not in loo:
            recommended[name] = None
            texts[name] = "Not fitted: " + fits[name]["error"]
            continue
        rmse = loo[name]["rmse"]
        if rmse != rmse:
            recommended[name] = None
            texts[name] = "The leave-one-out error could not be computed; judge the correction from the residuals."
        elif rmse < agreement["rmse"]:
            recommended[name] = True
            texts[name] = (
                "Correction recommended: the leave-one-out RMSE of the corrected counts (" + correct_number_text(rmse) +
                ") is smaller than the RMSE of the raw detected counts (" + correct_number_text(agreement["rmse"]) + ")."
            )
        else:
            recommended[name] = False
            texts[name] = (
                "Correction NOT recommended: the leave-one-out RMSE of the corrected counts (" + correct_number_text(rmse) +
                ") is not smaller than the RMSE of the raw detected counts (" + correct_number_text(agreement["rmse"]) +
                "). Report the raw counts together with the agreement statistics."
            )
    if recommended[model] is False:
        warnings.insert(0, texts[model])

    fingerprint = correct_calibration_fingerprint(pairs)
    identity = [to_text(model), to_text(weighting), to_text(channel), to_text(channel_name)]
    identity += sorted([repr(detected[i]) + u":" + repr(manual[i]) for i in range(len(pairs))])
    calibration_id = u"CAL_" + hashlib.sha1(u"|".join(identity).encode("utf-8")).hexdigest()[:12].upper()
    return {
        "id": calibration_id, "model": model, "weighting": weighting,
        "channel": to_text(channel), "channel_name": to_text(channel_name),
        "pairs": pairs, "detected": detected, "manual": manual,
        "fits": fits, "fit": fits[model], "loo": loo, "agreement": agreement,
        "recommended": recommended, "recommendation_texts": texts,
        "warnings": warnings,
        "n_images": len(correct_distinct(pairs, "source_image")),
        "versions": correct_distinct(pairs, "app_version"),
        "fingerprint": fingerprint[0], "fingerprint_source": fingerprint[1],
        "fp_full": correct_distinct(pairs, "fp_full"), "fp_partial": correct_distinct(pairs, "fp_partial")
    }


# -----------------------------------------------------------------------------
# Applying the calibration to a result table
# -----------------------------------------------------------------------------

def correct_find_count_column(headers, channel=u"", count_column=None, channel_name=u""):
    """The object-count column to correct: the one named, or the single
    '<prefix>_Object_Count' column of the channel (found by its number, or
    by its name when the calibration data carry no channel number)."""
    if count_column is not None and to_text(count_column).strip() != u"":
        if to_text(count_column) not in headers:
            raise ValueError("The table has no column named " + to_text(count_column) + ".")
        return to_text(count_column)
    candidates = []
    for name in headers:
        if not name.endswith("_Object_Count"):
            continue
        if name.endswith("_Edge_Touching_Object_Count") or name.endswith("_Excluded_Edge_Object_Count"):
            continue
        candidates.append(name)
    number = correct_channel_number(channel)
    chosen = candidates
    if number != u"":
        # The engine's own column when number and name are both known, so
        # that further count columns of the channel cannot be confused with it.
        if to_text(channel_name).strip() != u"":
            exact = u"C" + number + u"_" + sanitize_name(channel_name) + u"_Object_Count"
            if exact in candidates:
                return exact
        chosen = [name for name in candidates if name.startswith(u"C" + number + u"_")]
    elif to_text(channel_name).strip() != u"":
        ending = re.escape(u"_" + sanitize_name(channel_name) + u"_Object_Count") + u"$"
        named = [name for name in candidates if re.match(u"^C[0-9]+" + ending, name) is not None]
        if len(named) > 0:
            chosen = named
    if len(chosen) == 1:
        return chosen[0]
    if len(chosen) == 0:
        raise ValueError(
            "The table has no object-count column" + (u" for channel C" + number if number != u"" else u"") +
            ". Columns found: " + (u", ".join(candidates) if len(candidates) > 0 else u"none") + "."
        )
    raise ValueError("Choose the column to correct; the table has several: " + u", ".join(chosen) + ".")


def corrected_column_names(count_column):
    """The appended columns, in order, for one object-count column."""
    stem = count_column[:-len("_Object_Count")] if count_column.endswith("_Object_Count") else count_column
    density = stem + "_Object_Density_Per_Area_Corrected"
    return {
        "count": count_column + "_Corrected",
        "count_low": count_column + "_Corrected_Lower95",
        "count_high": count_column + "_Corrected_Upper95",
        "density": density,
        "density_low": density + "_Lower95",
        "density_high": density + "_Upper95",
        "extrapolated": count_column + "_Correction_Extrapolated",
        "clipped": count_column + "_Correction_Clipped",
        "order": [
            count_column + "_Corrected", count_column + "_Corrected_Lower95", count_column + "_Corrected_Upper95",
            density, density + "_Lower95", density + "_Upper95",
            count_column + "_Correction_Extrapolated", count_column + "_Correction_Clipped", "Calibration_ID",
            "Correction_Run_ID", "Correction_Application_Version"
        ]
    }


def correct_larger_roi_note(calibration):
    """What to keep in mind when ROIs hold more objects than any calibration
    field, which is the usual case when small fields were counted by hand."""
    if calibration["model"] == CORRECT_MODEL_LINE:
        return (
            "ROIs with more objects than any calibration field: a line with an intercept was only checked "
            "inside the calibrated range and should not be used for them. Use the proportional model."
        )
    text = (
        "ROIs with more objects than any calibration field: the proportional model assumes that the same "
        "fraction of objects is missed or merged in the larger ROI, which holds only if its cell density and "
        "staining resemble the calibration fields. This cannot be checked from the counts. For a much larger "
        "ROI the interval mostly shows how well the average correction factor is known; it does not include "
        "differences of that factor between sections or regions, and is too narrow when the factor varies."
    )
    if calibration["weighting"] == CORRECT_WEIGHT_NONE:
        text += (
            " Ordinary least squares also assumes the same scatter at every count. When the scatter grows "
            "with the count, as it usually does, the interval is too narrow from the upper end of the "
            "calibrated range upward; the weighting '" + CORRECT_WEIGHT_COUNT + "' lets the scatter grow "
            "with the count and suits large ROIs better."
        )
    return text


def correct_table_settings_check(table, table_dir, calibration, count_column):
    """Compare the settings and version of the table to be corrected with
    those of the calibration. Returns (table fingerprint, verdict, warnings)."""
    warnings = []
    channel = calibration["channel"]
    if channel == u"":
        match = re.match(u"^C([0-9]+)_", count_column)
        if match is not None:
            channel = match.group(1)
    full, partial, by_eye = correct_folder_fingerprints(table_dir, channel)
    verdict = "NOT_CHECKED"
    shown = full if full != u"" else partial
    cal_full = calibration["fp_full"]
    cal_partial = calibration["fp_partial"]
    all_full = len([pair for pair in calibration["pairs"] if to_text(pair.get("fp_full", u"")) != u""]) == len(calibration["pairs"])
    if calibration["fingerprint"] == u"MIXED":
        verdict = "NOT_CHECKED"
    elif full != u"" and len(cal_full) == 1 and all_full:
        verdict = "SAME" if full == cal_full[0] else "DIFFERENT"
    elif partial != u"" and len(cal_partial) == 1 and calibration["fingerprint"] != u"":
        verdict = "SAME_SEGMENTATION_SETTINGS_ONLY" if partial == cal_partial[0] else "DIFFERENT"
        shown = partial
    if by_eye and verdict != "DIFFERENT":
        # The profile agrees, but the threshold itself was set by hand for
        # each image, so equal settings do not mean equal thresholds.
        verdict = "THRESHOLD_SET_BY_EYE"
        warnings.append(
            "The table was made with the '" + CORRECT_BY_EYE_MODE + "' threshold, a value set by eye for each "
            "image. It cannot be checked that it equals the threshold of the calibration fields."
        )
    elif verdict == "DIFFERENT":
        warnings.append(
            "The table was made with DIFFERENT detection settings than the calibration fields. The "
            "calibration does not apply to it: validate fields with the settings of this table."
        )
    elif verdict == "SAME_SEGMENTATION_SETTINGS_ONLY":
        warnings.append(
            "Only the segmentation settings recorded in thresholds.csv could be compared with the "
            "calibration (they agree). Projection, preprocessing and manual threshold values were not checked."
        )
    elif verdict == "NOT_CHECKED":
        warnings.append(
            "It could not be checked that the table was made with the settings of the calibration fields. "
            "A calibration is valid only for the same settings profile, stain and imaging conditions."
        )
    lookup = correct_header_lookup(table["headers"])
    if "application_version" in lookup and len(calibration["versions"]) > 0:
        versions = sorted(set([correct_table_value(row, lookup, "Application_Version") for row in table["rows"]]))
        other = [value for value in versions if value != u"" and value not in calibration["versions"]]
        if len(other) > 0:
            warnings.append(
                "The table was written by application version " + u", ".join(other) +
                " but the calibration fields by " + u", ".join(calibration["versions"]) + "."
            )
    return (shown, verdict, warnings)


def correct_count_table(table_path, calibration, out_dir=None, count_column=None):
    """Write <table name>_corrected.csv: the table at table_path (a
    roi_summary.csv or batch_roi_summary.csv) with every original column
    unchanged and the corrected count, its 95 percent prediction interval,
    the corrected density and the flags appended at the end.

    The original file is only read. Returns a dictionary with the new
    file's path, headers and rows, the counts of flagged rows, and warnings."""
    table_path = to_text(table_path)
    if not os.path.isfile(table_path):
        raise ValueError("The table to correct does not exist: " + table_path)
    table = correct_read_csv(table_path)
    headers = table["headers"]
    if len(set(headers)) != len(headers):
        raise ValueError("The table has two columns with the same name and cannot be corrected safely.")
    column = correct_find_count_column(headers, calibration["channel"], count_column, calibration["channel_name"])
    names = corrected_column_names(column)
    for name in names["order"]:
        if name in headers:
            raise ValueError(
                "The table already has the column " + name + ". Correct the original table, not a corrected copy."
            )
    fit = calibration["fit"]
    decimal_comma = table["delimiter"] != u","
    count_index = headers.index(column)
    area_index = headers.index("ROI_Area") if "ROI_Area" in headers else -1
    fingerprint, verdict, warnings = correct_table_settings_check(table, os.path.dirname(os.path.abspath(table_path)), calibration, column)
    if area_index < 0:
        warnings.append("The table has no ROI_Area column, so the corrected density is left empty.")
    if calibration["channel_name"] != u"":
        ending = u"_" + sanitize_name(calibration["channel_name"]) + u"_Object_Count"
        if column.endswith("_Object_Count") and not column.endswith(ending):
            warnings.append(
                "The corrected column " + column + " does not carry the channel name of the calibration fields (" +
                calibration["channel_name"] + "). A calibration is valid for one stain only."
            )

    if out_dir is None:
        out_dir = os.path.dirname(os.path.abspath(table_path))
    out_path = os.path.join(out_dir, strip_extension(table_path) + "_corrected.csv")
    if os.path.exists(out_path):
        raise ValueError(
            "The file " + to_text(out_path) + " already exists and is not overwritten. "
            "Choose another output folder, or move the earlier corrected table away."
        )
    # Like every table of the application, the new columns name the folder
    # they were written into and the version that computed them.
    correction_run = os.path.basename(os.path.abspath(out_dir))

    rows = []
    totals = {"corrected": 0, "extrapolated": 0, "clipped": 0, "blank": 0, "above": 0}
    for values in table["rows"]:
        row = {}
        for index in range(len(headers)):
            row[headers[index]] = values[index]
        row["Calibration_ID"] = calibration["id"]
        row["Correction_Run_ID"] = correction_run
        row["Correction_Application_Version"] = APP_VERSION
        result = predict_corrected_count(fit, correct_number(values[count_index], decimal_comma))
        if result is None:
            totals["blank"] += 1
            rows.append(row)
            continue
        totals["corrected"] += 1
        row[names["count"]] = result["estimate"]
        row[names["count_low"]] = result["lower"]
        row[names["count_high"]] = result["upper"]
        row[names["extrapolated"]] = 1 if result["extrapolated"] else 0
        row[names["clipped"]] = 1 if result["clipped"] else 0
        if result["extrapolated"]:
            totals["extrapolated"] += 1
            if result["detected"] > fit["d_max"]:
                totals["above"] += 1
        if result["clipped"]:
            totals["clipped"] += 1
        area = correct_number(values[area_index], decimal_comma) if area_index >= 0 else None
        if area is not None and area > 0:
            row[names["density"]] = result["estimate"] / area
            row[names["density_low"]] = result["lower"] / area
            row[names["density_high"]] = result["upper"] / area
        rows.append(row)

    if totals["extrapolated"] > 0:
        warnings.append(
            str(totals["extrapolated"]) + " of " + str(totals["corrected"]) + " rows have a detected count OUTSIDE the "
            "calibrated range (" + correct_number_text(fit["d_min"]) + " to " + correct_number_text(fit["d_max"]) +
            "). Their corrected values are extrapolations and are flagged in the table; treat them with caution."
        )
        if totals["above"] > 0:
            warnings.append(correct_larger_roi_note(calibration))
    if totals["clipped"] > 0:
        warnings.append(str(totals["clipped"]) + " corrected count(s) were below zero and were set to zero (flagged).")
    if totals["blank"] > 0:
        warnings.append(str(totals["blank"]) + " row(s) have no usable count in " + column + " and were left uncorrected.")

    ensure_dir(out_dir)
    out_headers = headers + names["order"]
    write_csv(out_path, out_headers, rows)
    return {
        "path": out_path, "headers": out_headers, "rows": rows, "count_column": column, "columns": names,
        "table_name": os.path.basename(table_path), "table_fingerprint": fingerprint, "settings_match": verdict,
        "warnings": warnings, "corrected": totals["corrected"], "extrapolated": totals["extrapolated"],
        "clipped": totals["clipped"], "blank": totals["blank"]
    }


# -----------------------------------------------------------------------------
# Writing the calibration tables
# -----------------------------------------------------------------------------

def correct_count_value(value):
    """A whole number is written as 7, not 7.0; anything else unchanged."""
    if value == int(value):
        return int(value)
    return value


def correct_pair_rows(calibration, run_id):
    fit = calibration["fit"]
    loo = calibration["loo"][calibration["model"]]["predictions"]
    rows = []
    for index in range(len(calibration["pairs"])):
        pair = calibration["pairs"][index]
        detected = calibration["detected"][index]
        manual = calibration["manual"][index]
        fingerprint = to_text(pair.get("fp_full", u""))
        if fingerprint == u"":
            fingerprint = to_text(pair.get("fp_partial", u""))
        rows.append({
            "Run_ID": run_id, "Application_Version": APP_VERSION, "Calibration_ID": calibration["id"],
            "Pair_Index": index + 1,
            "Channel": pair.get("channel", u"") or calibration["channel"],
            "Channel_Name": pair.get("channel_name", u"") or calibration["channel_name"],
            "Source_Image": pair.get("source_image", u""), "Frame": pair.get("frame", u""),
            "Parent_ROI_Name": pair.get("parent_roi", u""),
            "Validation_Run_ID": pair.get("run_id", u""),
            "Validation_Application_Version": pair.get("app_version", u""),
            "Calibration_Source_File": pair.get("source_file", u""),
            "Detected_Count": correct_count_value(detected), "Manual_Count": correct_count_value(manual),
            "Difference_Detected_Minus_Manual": correct_count_value(detected - manual),
            "Percent_Error": (detected - manual) / manual * 100.0 if manual > 0 else u"",
            "Model": calibration["model"], "Weighting": calibration["weighting"],
            "Fitted_Manual_Count": fit["fitted"][index],
            "Residual_Manual_Minus_Fitted": fit["residuals"][index],
            "LeaveOneOut_Corrected_Count": loo[index] if loo[index] is not None else u"",
            "LeaveOneOut_Error": loo[index] - manual if loo[index] is not None else u"",
            "Settings_Fingerprint": fingerprint
        })
    return rows


def correct_fit_rows(calibration, run_id, corrected=None):
    """One row per model; "Selected" marks the model that was applied."""
    agreement = calibration["agreement"]
    warnings = list(calibration["warnings"])
    if corrected is not None:
        warnings += corrected["warnings"]
    rows = []
    for name in CORRECT_MODELS:
        fit = calibration["fits"][name]
        row = {
            "Run_ID": run_id, "Application_Version": APP_VERSION, "Calibration_ID": calibration["id"],
            "Channel": calibration["channel"], "Channel_Name": calibration["channel_name"],
            "Model": name, "Selected": 1 if name == calibration["model"] else 0,
            "Weighting": calibration["weighting"],
            "N_Fields": len(calibration["pairs"]), "N_Images": calibration["n_images"] if calibration["n_images"] > 0 else u"",
            "Detected_Min": correct_count_value(min(calibration["detected"])),
            "Detected_Max": correct_count_value(max(calibration["detected"])),
            "Mean_Difference_Detected_Minus_Manual": agreement["mean_difference"],
            "SD_Difference": agreement["sd_difference"],
            "Limits_Of_Agreement_Lower": agreement["loa_low"], "Limits_Of_Agreement_Upper": agreement["loa_high"],
            "MAPE_Percent": agreement["mape"], "MAPE_Fields_Used": agreement["mape_fields"],
            "Raw_RMSE": agreement["rmse"], "Raw_MAE": agreement["mae"],
            "Recommendation": calibration["recommendation_texts"][name],
            "Settings_Fingerprint": calibration["fingerprint"],
            "Settings_Fingerprint_Source": calibration["fingerprint_source"],
            "Validation_Application_Versions": u"; ".join(calibration["versions"]),
            "Warnings": u" | ".join(warnings), "Interval_Note": CORRECT_INTERVAL_NOTE
        }
        if calibration["recommended"][name] is not None:
            row["Correction_Recommended"] = 1 if calibration["recommended"][name] else 0
        if "error" in fit:
            row["Fit_Status"] = "NOT_FITTED: " + fit["error"]
        else:
            row["Fit_Status"] = "OK"
            row["Slope"] = fit["slope"]
            row["Slope_SE"] = fit["slope_se"]
            row["Slope_Lower95"] = fit["slope_low"]
            row["Slope_Upper95"] = fit["slope_high"]
            row["Intercept"] = fit["intercept"]
            row["Degrees_Of_Freedom"] = fit["df"]
            row["T_Quantile_975"] = fit["t"]
            row["T_Quantile_Source"] = fit["t_source"]
            row["Residual_SD"] = fit["residual_sd"]
            row["R_Squared"] = fit["r_squared"]
            row["R_Squared_Uncentered"] = fit["r_squared_uncentered"]
            row["AICc"] = fit["aicc"]
            row["LeaveOneOut_RMSE_Corrected"] = calibration["loo"][name]["rmse"]
            row["LeaveOneOut_MAE_Corrected"] = calibration["loo"][name]["mae"]
            if name == CORRECT_MODEL_LINE:
                row["Intercept_SE"] = fit["intercept_se"]
                row["Intercept_Lower95"] = fit["intercept_low"]
                row["Intercept_Upper95"] = fit["intercept_high"]
                row["Intercept_CI_Includes_Zero"] = 1 if fit["intercept_low"] <= 0 <= fit["intercept_high"] else 0
        if corrected is not None:
            row["Corrected_Table"] = corrected["table_name"]
            row["Corrected_Count_Column"] = corrected["count_column"]
            row["Table_Settings_Fingerprint"] = corrected["table_fingerprint"]
            row["Settings_Match"] = corrected["settings_match"]
            if name == calibration["model"]:
                row["Rows_Corrected"] = corrected["corrected"]
                row["Rows_Extrapolated"] = corrected["extrapolated"]
                row["Rows_Clipped"] = corrected["clipped"]
        rows.append(row)
    return rows


def write_count_calibration(calibration, out_dir, corrected=None):
    """Write count_calibration.csv (one row per calibration pair) and
    count_calibration_fit.csv (one row per model) into out_dir.
    corrected: the result of correct_count_table, when a table was corrected
    with this calibration, so that the fit table records it."""
    ensure_dir(out_dir)
    run_id = os.path.basename(os.path.abspath(out_dir))
    pair_rows = correct_pair_rows(calibration, run_id)
    fit_rows = correct_fit_rows(calibration, run_id, corrected)
    pairs_path = os.path.join(out_dir, "count_calibration.csv")
    fit_path = os.path.join(out_dir, "count_calibration_fit.csv")
    write_csv(pairs_path, COUNT_CALIBRATION_HEADERS, pair_rows)
    write_csv(fit_path, COUNT_CALIBRATION_FIT_HEADERS, fit_rows)
    return {"pairs_path": pairs_path, "fit_path": fit_path, "pair_rows": pair_rows, "fit_rows": fit_rows}


def correct_remove_if_empty(folder):
    try:
        if len(os.listdir(folder)) == 0:
            os.rmdir(folder)
    except ANY_ERROR:
        pass


def run_count_correction(calibration_folder=None, calibration_csv=None, table_path=None, channel=None,
                         output_root=None, model=None, weighting=None, count_column=None,
                         global_options=None, log=None):
    """The whole procedure: collect the calibration pairs of one channel, fit
    and score the calibration, correct one result table (optional), and write
    everything into a new folder below output_root.

    model and weighting default to the profile settings correct_model and
    correct_weighting of global_options, then to the proportional model with
    ordinary least squares. Nothing existing is modified.
    Returns a dictionary with out_dir, calibration, corrected (or None),
    the written paths, and every warning."""
    def say(message):
        if log is not None:
            log(message)

    options = global_options if global_options is not None else {}
    if model is None or to_text(model).strip() == u"":
        model = to_text(options.get("correct_model", CORRECT_MODEL_ORIGIN))
    if weighting is None or to_text(weighting).strip() == u"":
        weighting = to_text(options.get("correct_weighting", CORRECT_WEIGHT_NONE))
    has_table = table_path is not None and to_text(table_path).strip() != u""
    if has_table and not os.path.isfile(to_text(table_path)):
        # Checked before any folder is made, so a mistyped path leaves nothing behind.
        raise ValueError("The table to correct does not exist: " + to_text(table_path))
    collected = collect_calibration_pairs(calibration_folder, calibration_csv, channel)
    channel_id = collected["channel"]
    if channel_id == u"":
        channel_id = correct_channel_number(channel)
    calibration = build_count_calibration(collected["pairs"], model, weighting, channel_id, collected["channel_name"])
    calibration["warnings"] += collected["notes"]

    if output_root is None or to_text(output_root).strip() == u"":
        if has_table:
            output_root = os.path.dirname(os.path.abspath(to_text(table_path)))
        elif calibration_folder is not None and to_text(calibration_folder).strip() != u"":
            output_root = to_text(calibration_folder)
        else:
            output_root = os.path.dirname(os.path.abspath(to_text(calibration_csv)))
    base = sanitize_name(strip_extension(table_path)) if has_table else "Calibration"
    # Named by the engine's run-folder helper: <base>_CellQuant_CountCorrection_<time>,
    # with base shortened only when a written path would be too long.
    out_dir = unique_run_dir(
        os.path.abspath(to_text(output_root)), base, None, "CountCorrection",
        ["count_calibration_fit.csv", (strip_extension(table_path) + "_corrected.csv") if has_table else u""], []
    )

    corrected = None
    try:
        if has_table:
            corrected = correct_count_table(table_path, calibration, out_dir, count_column)
        written = write_count_calibration(calibration, out_dir, corrected)
    except ANY_ERROR as exc:
        # A table that cannot be corrected leaves no empty folder behind.
        correct_remove_if_empty(out_dir)
        raise exc

    fit = calibration["fit"]
    say("Count calibration " + calibration["id"] + ": " + str(fit["n"]) + " fields, model '" + model + "'.")
    say(
        "Corrected count = " + correct_number_text(fit["intercept"]) + " + " + correct_number_text(fit["slope"]) +
        " x detected count (slope 95 percent CI " + correct_number_text(fit["slope_low"]) + " to " +
        correct_number_text(fit["slope_high"]) + ", residual SD " + correct_number_text(fit["residual_sd"]) + ")."
    )
    say(calibration["recommendation_texts"][model])
    warnings = list(calibration["warnings"])
    if corrected is not None:
        warnings += corrected["warnings"]
        say("Corrected table: " + to_text(corrected["path"]))
    written_paths = [written["pairs_path"], written["fit_path"]]
    if corrected is not None:
        written_paths.append(corrected["path"])
    longest = max([len(to_text(path)) for path in written_paths])
    if longest > CORRECT_PATH_LIMIT:
        # Java writes such paths, but Excel, Explorer and Python on Windows
        # cannot open a path longer than 259 characters.
        warnings.append(
            "A result path is " + str(longest) + " characters long, more than the " + str(CORRECT_PATH_LIMIT) +
            " that Excel, Explorer and Python on Windows can open. Choose a shorter output folder."
        )
    for warning in warnings:
        say("WARNING: " + warning)
    say(CORRECT_INTERVAL_NOTE)
    return {
        "out_dir": out_dir, "calibration": calibration, "corrected": corrected,
        "pairs_path": written["pairs_path"], "fit_path": written["fit_path"],
        "corrected_path": corrected["path"] if corrected is not None else None,
        "warnings": warnings
    }


# -----------------------------------------------------------------------------
# Registration: two profile settings, shown by the window; no per-run step
# -----------------------------------------------------------------------------

register_global_settings([
    ("correct.model", "correct_model", "text", CORRECT_MODEL_ORIGIN),
    ("correct.weighting", "correct_weighting", "text", CORRECT_WEIGHT_NONE)
])
register_settings_form(
    "Corrected counts (calibration against manual counts)",
    "Fits manual against detected counts of validated fields and adds corrected counts with 95 percent "
    "prediction intervals to a copy of a result table. Valid only for the settings, stain and imaging "
    "conditions of the validated fields; the intervals cover counting error only.",
    [
        ("correct_model", "Calibration model", "choice", CORRECT_MODELS),
        ("correct_weighting", "Weighting of the calibration fields", "choice", CORRECT_WEIGHTINGS)
    ]
)
register_engine_module("cellquantcorrect.py", "2.3.0")

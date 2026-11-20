# Cell Quant Studio Universal - tuning of segmentation settings from manual
# point marks, with held-out scoring (optional module)
#
# Purpose: instead of adjusting a threshold or a size limit by eye, search
# chosen settings of ONE channel for the combination that agrees best with
# manual marks (one point per object, made with Fiji's multi-point tool) on
# several annotated fields, and report how well settings found this way work
# on fields that were NOT used to find them. This applies established
# practice (tuning segmentation parameters against annotation: Held et al.
# 2013; Teodoro et al. 2017) to this application; it is not a new method.
#
# Definitions used everywhere in this file:
#   - A field is one image (one frame) with optional parent ROIs (the whole
#     image when none are given) and a list of marks.
#   - Scoring is the same as in validate_objects_against_points: per parent
#     ROI, objects come from analyze_objects on the final mask, a mark belongs
#     to the parent ROI that contains it, and match_points_to_objects counts
#     true positives (objects holding at least one mark), false positives
#     (objects holding none), missed marks, and merged marks (the second and
#     later marks inside one object). Marks outside every parent ROI are
#     ignored and reported.
#   - Counts are pooled (summed) over parent ROIs and fields before a score
#     is computed: precision = TP / (TP + FP), recall = TP / marks,
#     F1 = 2 * precision * recall / (precision + recall), which equals
#     2 TP / (2 TP + FP + missed + merged). The absolute count error is
#     the sum, over every parent ROI, of |detected - marks| in that ROI.
#   - A threshold range that holds no pixel value of a field (for example a
#     Manual threshold above the brightest pixel of a dim field) detects
#     nothing on that field: every mark there is missed. That is a result,
#     not a failure; such rows say "OK: nothing detected ...".
#   - The search is a full grid over the listed values. The best combination
#     has the highest objective; ties go to the combination with the fewest
#     settings changed from the starting configuration, then the smallest
#     summed relative change, then the earliest position in the grid.
#   - Held-out scoring: every group of fields (by default every field) is
#     left out once; the settings are tuned on the remaining fields and
#     scored on the left-out fields. The pooled held-out score is the number
#     to quote. The score of the final settings (tuned on all fields) on the
#     same fields is optimistic and is labelled so in every output.
#
# Each field is segmented once per combination; cross-validation then only
# re-adds the stored counts, so it costs no extra image processing.
#
# This file is loaded by load_engine_modules() into the namespace of
# cellquantengine.py: the engine's names are globals here. It has no user
# interface, opens no window, and changes no Fiji option.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij import IJ
from ij.gui import Roi

import math
import os
import time
import traceback


TUNE_OBJECTIVE_F1 = "F1"
TUNE_OBJECTIVE_COUNT = "Absolute count error"
TUNE_OBJECTIVES = [TUNE_OBJECTIVE_F1, TUNE_OBJECTIVE_COUNT]
TUNE_DEFAULT_FACTORS = "0.7, 0.85, 1.0, 1.15, 1.3"
TUNE_MAX_COMBINATIONS = 2000
TUNE_POINT_SUFFIXES = ["_points.zip", "_points.roi"]
TUNE_GROUP_IMAGE = "Image"
TUNE_GROUP_FOLDER = "Subfolder"
TUNE_GROUPINGS = [TUNE_GROUP_IMAGE, TUNE_GROUP_FOLDER]
TUNE_NO_HELD_OUT = "none: tuned and scored on the same field (optimistic)"
TUNE_NO_HELD_OUT_GROUP = "none: tuned and scored on the same fields (optimistic)"
# Text of the engine's error for a threshold range that holds no pixel value
# of the image (build_masks). Such a combination detects nothing; it has not
# failed. An engine that returns an empty mask instead gives the same counts.
TUNE_EMPTY_RANGE_TEXT = "Threshold upper value is below the lower value"
TUNE_EMPTY_RANGE_STATUS = "OK: no pixel value lies in the threshold range, nothing detected"

# Settings that can be searched: (key, level, threshold mode it needs).
# The level says which intermediate result a change invalidates:
#   prep     - the prepared image (filters applied during preparation)
#   mask     - the mask (threshold and binary operations)
#   particle - the particle analysis of an unchanged mask
#   size     - only the size limits, applied to an unchanged object list
# The order of this list is the order of the grid and of table columns.
TUNE_PARAMETERS = [
    ("gaussian_sigma", "prep", ""),
    ("median_radius", "prep", ""),
    ("manual_low", "mask", "Manual"),
    ("auto_method", "mask", "Auto"),
    ("auto_scale", "mask", "Auto"),
    ("local_radius", "mask", "Local Phansalkar"),
    ("local_parameter_1", "mask", "Local Phansalkar"),
    ("clahe", "mask", "Local Phansalkar"),
    ("clahe_max_slope", "mask", "Local Phansalkar"),
    ("open_iterations", "mask", ""),
    ("fill_holes", "mask", ""),
    ("watershed", "mask", ""),
    ("min_circularity", "particle", ""),
    ("min_size", "size", ""),
    ("max_size", "size", "")
]


def tune_prefixed(prefix, names):
    """prefix + each name. A function, so that building the header lists
    below leaves no loop variable behind in the shared namespace."""
    return [prefix + name for name in names]


def tune_parameter_keys():
    return [spec[0] for spec in TUNE_PARAMETERS]


TUNE_PARAMETER_KEYS = tune_parameter_keys()

TUNE_COUNT_KEYS = ["manual", "detected", "true_positive", "false_positive", "missed", "merged"]
TUNE_METRIC_HEADERS = [
    "Manual_Count", "Detected_Count", "True_Positive_Objects", "False_Positive_Objects",
    "Missed_Annotations", "Merged_Annotations", "Precision", "Recall", "F1",
    "Count_Difference", "Count_Error_Percent", "Abs_Count_Error", "Abs_Count_Error_Percent"
]
TUNE_VALUE_HEADERS = tune_prefixed("P_", TUNE_PARAMETER_KEYS)
TUNE_ID_HEADERS = ["Run_ID", "Application_Version", "Channel", "Channel_Name"]

TUNE_FIELD_HEADERS = [
    "Run_ID", "Application_Version", "Field_Index", "Field", "Group", "Image_Path",
    "ROI_Source", "Points_Source", "Status", "Message", "Width", "Height",
    "Parent_ROI_Count", "Marks_Total", "Marks_Inside_ROIs", "Marks_Outside_ROIs", "Seconds"
]
TUNE_GRID_HEADERS = (
    TUNE_ID_HEADERS + ["Combination", "Settings"] + TUNE_VALUE_HEADERS +
    ["Changed_Settings", "Relative_Change", "Is_Start", "Status", "Fields"] + TUNE_METRIC_HEADERS +
    ["Objective", "Objective_Value", "Rank", "Selected_Final", "Folds_Selected"]
)
TUNE_SCORE_HEADERS = (
    ["Run_ID", "Application_Version", "Combination", "Field_Index", "Field", "Group", "Status"] +
    TUNE_METRIC_HEADERS
)
TUNE_FOLD_HEADERS = (
    TUNE_ID_HEADERS + ["Fold", "Held_Out_Group", "Held_Out_Field", "Training_Fields", "Training_Marks",
                       "Selected_Combination", "Selected_Settings"] + TUNE_VALUE_HEADERS +
    ["Training_F1_Optimistic"] + TUNE_METRIC_HEADERS +
    tune_prefixed("Baseline_", TUNE_METRIC_HEADERS) + ["Note"]
)
TUNE_STABILITY_HEADERS = TUNE_ID_HEADERS + [
    "Parameter", "Value", "Is_Start_Value", "Folds_Selected", "Folds_Scored",
    "Fraction_Of_Folds", "Selected_Final"
]
TUNE_SUMMARY_HEADERS = TUNE_ID_HEADERS + [
    "Objective", "Validation_Design", "Fields", "Fields_Failed", "Images_Without_Annotation",
    "Groups", "Folds", "Folds_Scored", "Marks_Inside_ROIs", "Marks_Outside_ROIs",
    "Search_Space", "Combinations", "Combinations_Usable", "Start_Settings",
    "Tuned_Combination", "Tuned_Settings", "Tuned_Profile",
    "HeldOut_Fields", "HeldOut_F1_Field_Min", "HeldOut_F1_Field_Max"
] + (
    tune_prefixed("HeldOut_", TUNE_METRIC_HEADERS) +
    tune_prefixed("HeldOut_Baseline_", TUNE_METRIC_HEADERS) +
    tune_prefixed("Resubstitution_Optimistic_", TUNE_METRIC_HEADERS) +
    tune_prefixed("AllFields_Baseline_", TUNE_METRIC_HEADERS)
) + ["Evaluations", "Seconds", "Evaluations_Per_Second", "Note"]

# (property key, option key, kind, default): saved in profiles so that the
# window and the headless script can start a tuning run from a profile.
TUNE_SETTING_SPECS = [
    ("tune.channel", "tune_channel", "text", "0"),
    ("tune.objective", "tune_objective", "text", TUNE_OBJECTIVE_F1),
    ("tune.search_space", "tune_search_space", "text", ""),
    ("tune.factors", "tune_factors", "text", TUNE_DEFAULT_FACTORS),
    ("tune.max_combinations", "tune_max_combinations", "text", str(TUNE_MAX_COMBINATIONS)),
    ("tune.group_by", "tune_group_by", "text", TUNE_GROUP_IMAGE)
]


# -----------------------------------------------------------------------------
# Search space
# -----------------------------------------------------------------------------

def tune_text_value(value):
    """One setting value as short text."""
    if isinstance(value, bool):
        return bool_text(value)
    if isinstance(value, float):
        return to_text("%.12g" % value)
    return to_text(value)


def tune_settings_text(space, values):
    return "; ".join([entry[0] + "=" + tune_text_value(values[entry[0]]) for entry in space])


def describe_tune_search_space(space):
    """Text form of a normalized search space; parse_tune_search_space reads it back."""
    parts = []
    for entry in space:
        parts.append(entry[0] + "=" + ",".join([tune_text_value(value) for value in entry[1]]))
    return "; ".join(parts)


def parse_tune_search_space(text):
    """Turn 'manual_low = 800, 1000; watershed = true, false' into a
    dictionary setting -> list of values (still text; they are checked and
    converted by normalize_tune_search_space)."""
    space = {}
    for part in to_text(text).replace("\n", ";").split(";"):
        part = part.strip()
        if part == "":
            continue
        if "=" not in part:
            raise ValueError("Search space entry '" + part + "' must look like: setting = value, value")
        pieces = part.split("=", 1)
        key = pieces[0].strip().lower()
        if key in space:
            raise ValueError("The search space lists " + key + " twice.")
        space[key] = [item.strip() for item in pieces[1].split(",") if item.strip() != ""]
    return space


def default_tune_search_space(config, factors=None):
    """A modest grid around the current settings of one channel: the
    threshold setting of its threshold mode and the minimum object area,
    each multiplied by the factors, plus watershed and fill holes on/off.
    A setting whose current value is 0 is left out, because multiplying it
    changes nothing; list explicit values for it instead."""
    validate_unattended_configs([config], "Tuning")
    if factors is None or to_text(factors).strip() == "":
        factors = TUNE_DEFAULT_FACTORS
    if not isinstance(factors, (list, tuple)):
        factors = parse_factor_list(factors)
    space = {}
    numeric = ["min_size"]
    if config["threshold_mode"] == "Manual":
        numeric.append("manual_low")
    elif config["threshold_mode"] == "Auto":
        numeric.append("auto_scale")
    elif config["threshold_mode"] == "Local Phansalkar":
        numeric.append("local_radius")
        numeric.append("local_parameter_1")
    for key in numeric:
        if config[key] == 0:
            continue
        values = []
        for factor in factors:
            # Six significant digits keep the saved profile readable; the
            # factor 1 is the current value itself, never a rounded copy.
            value = config[key] if factor == 1.0 else float("%.6g" % (config[key] * factor))
            if key == "local_radius":
                # Half up, the same under Jython 2.7 and Python 3.
                value = max(1, int(math.floor(value + 0.5)))
            if key == "manual_low" and config["manual_high"] != -1 and value > config["manual_high"]:
                continue
            if value not in values:
                values.append(value)
        space[key] = values
    space["watershed"] = [config["watershed"], not config["watershed"]]
    space["fill_holes"] = [config["fill_holes"], not config["fill_holes"]]
    return space


TUNE_BOOLEAN_WORDS = ["true", "false", "1", "0", "yes", "no", "on", "off"]


def tune_check_raw(key, raw, label):
    """A switch must be written as a recognised word: parse_bool would read
    a typing error such as 'ture' silently as false."""
    if key in CHANNEL_BOOLEAN_KEYS and not isinstance(raw, bool):
        if to_text(raw).strip().lower() not in TUNE_BOOLEAN_WORDS:
            raise ValueError(
                label + ": '" + to_text(raw) + "' is not a value for " + key + ". Use true or false."
            )


def tune_check_value(config, key, value, label):
    """Checks that normalize_channel_config does not make for one value."""
    if key == "auto_method" and value not in AUTO_METHODS:
        raise ValueError(label + ": '" + to_text(value) + "' is not an automatic threshold method.")
    if key == "min_circularity" and (value > 1 or value > config["max_circularity"]):
        raise ValueError(label + ": minimum circularity " + str(value) + " is above the maximum circularity.")
    if key == "manual_low" and config["manual_high"] != -1 and value > config["manual_high"]:
        raise ValueError(label + ": manual lower threshold " + str(value) + " is above the upper threshold.")


def normalize_tune_search_space(config, search_space, max_combinations=None):
    """Check a search space against one channel configuration.

    Returns a list of (setting, values) in the fixed order of
    TUNE_PARAMETERS, each value converted exactly as a profile would hold
    it. Raises ValueError for a setting that cannot be tuned, a setting
    without effect in this threshold mode, a bad value, or a grid above the
    limit."""
    label = "C" + str(config["channel"])
    if not needs_objects(config):
        raise ValueError(
            "Point marks tune object detection. Set " + label + " to 'Objects only' or 'Objects + area'."
        )
    validate_unattended_configs([config], "Tuning")
    for key in search_space:
        if key not in TUNE_PARAMETER_KEYS:
            raise ValueError(
                "'" + to_text(key) + "' cannot be tuned. Settings that can be tuned: " +
                ", ".join(TUNE_PARAMETER_KEYS) + "."
            )
    space = []
    for spec in TUNE_PARAMETERS:
        key = spec[0]
        if key not in search_space:
            continue
        if spec[2] != "" and config["threshold_mode"] != spec[2]:
            raise ValueError(
                key + " cannot be tuned for " + label + ": it only has an effect in the '" + spec[2] +
                "' threshold mode, and " + label + " uses '" + config["threshold_mode"] + "'."
            )
        raw_values = search_space[key]
        if not isinstance(raw_values, (list, tuple)):
            raw_values = [raw_values]
        values = []
        for raw in raw_values:
            tune_check_raw(key, raw, label)
            candidate = dict(config)
            candidate[key] = raw
            value = normalize_channel_config(config["channel"], candidate)[key]
            tune_check_value(config, key, value, label)
            if value not in values:
                values.append(value)
        if len(values) == 0:
            raise ValueError("The search space lists no value for " + key + ".")
        space.append((key, values))
    if len(space) == 0:
        raise ValueError("The search space is empty: list at least one setting to tune.")
    searched = dict(space)
    if "clahe_max_slope" in searched and not config["clahe"] and True not in searched.get("clahe", []):
        raise ValueError(
            "clahe_max_slope cannot be tuned for " + label + " while CLAHE is off. "
            "Switch CLAHE on or add 'clahe = false, true' to the search space."
        )
    limit = TUNE_MAX_COMBINATIONS
    if max_combinations is not None and to_text(max_combinations).strip() != "":
        limit = parse_int_value(max_combinations, "Largest number of combinations", 1)
    total = 1
    for entry in space:
        total = total * len(entry[1])
    if total > limit:
        raise ValueError(
            "The search grid has " + str(total) + " combinations (" +
            " x ".join([str(len(entry[1])) + " " + entry[0] for entry in space]) +
            "); the limit is " + str(limit) + ". Shrink the grid: list fewer values or fewer settings."
        )
    return space


def tune_change_size(config, space, values):
    """How far a combination is from the starting configuration:
    (number of settings changed, summed relative change of the numbers).
    A number that starts at 0 is scaled by the largest value in its list."""
    changes = 0
    distance = 0.0
    for entry in space:
        key = entry[0]
        start = config[key]
        value = values[key]
        if value == start:
            continue
        changes += 1
        if not is_number(value) or not is_number(start):
            continue
        if key == "max_size" and (value == -1 or start == -1):
            # -1 means "no upper limit": a different kind of setting, not a distance.
            continue
        scale = abs(start)
        if scale == 0:
            scale = max([abs(option) for option in entry[1] if is_number(option)] + [0])
        if scale == 0:
            scale = 1.0
        distance += abs(value - start) / float(scale)
    return (changes, distance)


def tune_make_combination(config, space, values, number):
    variant = dict(config)
    variant.update(values)
    changes, distance = tune_change_size(config, space, values)
    status = ""
    if variant["max_size"] != -1 and variant["max_size"] < variant["min_size"]:
        status = "INVALID: maximum object area is below minimum object area"
    return {
        "number": number, "values": dict(values), "config": variant,
        "changes": changes, "distance": distance, "status": status, "baseline": False
    }


def tune_build_combinations(config, space):
    """Every combination of the grid, in a fixed order: the first setting of
    the space changes slowest. The starting configuration is added at the
    end as the baseline; it is scored like the others but never selected."""
    total = 1
    for entry in space:
        total = total * len(entry[1])
    combos = []
    for number in range(total):
        remainder = number
        values = {}
        for position in range(len(space) - 1, -1, -1):
            options = space[position][1]
            values[space[position][0]] = options[remainder % len(options)]
            remainder = remainder // len(options)
        combos.append(tune_make_combination(config, space, values, number + 1))
    start_values = {}
    for entry in space:
        start_values[entry[0]] = config[entry[0]]
    baseline = tune_make_combination(config, space, start_values, 0)
    baseline["baseline"] = True
    baseline["status"] = ""
    combos.append(baseline)
    return combos


# -----------------------------------------------------------------------------
# Annotated folders
# -----------------------------------------------------------------------------

def points_file_candidates(image_path, points_dir=None):
    folder = points_dir if points_dir else os.path.dirname(image_path)
    name = os.path.basename(to_text(image_path))
    stems = [strip_extension(name)]
    if stems[0].lower().endswith(".ome"):
        stems.append(stems[0][:-4])
    stems.append(name)
    candidates = []
    for stem in stems:
        for suffix in TUNE_POINT_SUFFIXES:
            candidates.append(os.path.join(folder, stem + suffix))
    return candidates


def find_points_file(image_path, points_dir=None):
    for candidate in points_file_candidates(image_path, points_dir):
        if os.path.isfile(candidate):
            return candidate
    return None


def discover_annotated_fields(folder, extensions=None, recursive=False, roi_dir=None,
                              points_dir=None, group_by=None):
    """List the annotated fields of a folder.

    For every image file: its parent ROIs come from the ROI file found by
    find_roi_file (the whole image when there is none), its marks from
    <image name>_points.zip or <image name>_points.roi (a ROI Manager set or
    one multi-point ROI). Without a points file, point ROIs stored in the
    parent-ROI file are used. A points file that holds no point ROI means
    "annotated, zero marks". An image with no marks from either source is
    not annotated and is reported, not used.

    Returns a dictionary: fields (list of field dictionaries for
    run_tuning), report (one row per image: Status ANNOTATED, NO_ANNOTATION
    or FAILED), and missing (paths of the images without annotation)."""
    folder = os.path.abspath(to_text(folder))
    if not os.path.isdir(folder):
        raise ValueError("The annotated folder does not exist: " + folder)
    if extensions is None or to_text(extensions).strip() == "":
        extensions = DEFAULT_BATCH_EXTENSIONS
    if not isinstance(extensions, (list, tuple)):
        extensions = parse_extension_list(extensions)
    group_by = TUNE_GROUP_IMAGE if group_by is None else to_text(group_by)
    if group_by not in TUNE_GROUPINGS:
        raise ValueError("Unknown grouping of fields: " + group_by)
    fields = []
    report = []
    missing = []
    for path in discover_batch_images(folder, extensions, parse_bool(recursive)):
        try:
            name = to_text(os.path.relpath(path, folder))
        except ANY_ERROR:
            name = to_text(os.path.basename(path))
        group = name
        if group_by == TUNE_GROUP_FOLDER:
            group = to_text(os.path.dirname(name))
            if group == "":
                group = to_text(os.path.basename(folder))
        row = {"Field": name, "Group": group, "Image_Path": to_text(path), "Status": "ANNOTATED", "Message": ""}
        report.append(row)
        try:
            rois = None
            stored_points = []
            row["ROI_Source"] = "Whole image (no ROI file)"
            roi_path = find_roi_file(path, roi_dir if roi_dir else None)
            if roi_path is not None:
                loaded = load_rois_from_file(roi_path)
                stored_points = [roi for roi in loaded if roi is not None and roi.getType() == Roi.POINT]
                area = [roi for roi in loaded if roi is not None and roi.isArea()]
                if len(area) > 0:
                    rois = area
                    row["ROI_Source"] = "ROI file " + to_text(os.path.basename(roi_path))
                else:
                    row["ROI_Source"] = "Whole image (no area ROI in " + to_text(os.path.basename(roi_path)) + ")"
            points_path = find_points_file(path, points_dir if points_dir else None)
            if points_path is not None:
                points = annotation_points_from_rois(load_rois_from_file(points_path))
                row["Points_Source"] = "Points file " + to_text(os.path.basename(points_path))
                if len(points) == 0:
                    row["Message"] = "The points file holds no point ROI: used as a field with zero marks."
            elif len(stored_points) > 0:
                points = annotation_points_from_rois(stored_points)
                row["Points_Source"] = "Point ROIs in " + to_text(os.path.basename(roi_path))
            else:
                row["Status"] = "NO_ANNOTATION"
                row["Message"] = "No points file (" + strip_extension(path) + "_points.zip or _points.roi) was found."
                missing.append(path)
                continue
            row["Marks_Total"] = len(points)
            fields.append({
                "name": name, "group": group, "path": path, "rois": rois, "points": points,
                "roi_source": row["ROI_Source"], "points_source": row["Points_Source"]
            })
        except ANY_ERROR as exc:
            row["Status"] = "FAILED"
            row["Message"] = "Could not read the annotation: " + to_text(exc)
    return {"fields": fields, "report": report, "missing": missing}


# -----------------------------------------------------------------------------
# Scoring one field for every combination
# -----------------------------------------------------------------------------

def tune_apply_filters(base, config):
    """The last steps of AnalysisEngine.prepare_single_channel (median, then
    Gaussian), applied to a copy of an image prepared without them. The
    earlier steps do not depend on these two settings, so the result is the
    image a full preparation with this configuration gives."""
    result = base.duplicate()
    result.setTitle(base.getTitle())
    if config["median_radius"] > 0:
        IJ.run(result, "Median...", "radius=" + str(config["median_radius"]))
    if config["gaussian_sigma"] > 0:
        IJ.run(result, "Gaussian Blur...", "sigma=" + str(config["gaussian_sigma"]))
    result.deleteRoi()
    return result


def tune_open_field(engine, field, config, needs_base, calibration, force_bioformats, log):
    """Get one field ready for scoring.

    A field dictionary holds name, points (list of (x, y) pixel positions)
    or point_rois, optional rois (parent ROIs), optional group and frame,
    and ONE of: prepared (the prepared image of the tuned channel; it is
    left open), source (an open multichannel image; left open), or path (an
    image file; opened and closed here).

    Returns a dictionary with the prepared image (or, when a preparation
    filter is searched, the image prepared without median and Gaussian
    filters), the validated parent ROIs, and the marks per parent ROI."""
    work = {"prepared": None, "base": None, "owned": []}
    opened = None
    try:
        if field.get("prepared") is not None:
            if needs_base:
                raise ValueError(
                    "gaussian_sigma and median_radius are applied while an image is prepared, so they cannot be "
                    "tuned on a field given as an already prepared image. Give the source image or the image file."
                )
            if field["prepared"].getBitDepth() == 24:
                raise ValueError("RGB images are not supported for segmentation.")
            work["prepared"] = field["prepared"]
        else:
            source = field.get("source")
            if source is None:
                path = to_text(field.get("path", ""))
                if path == "":
                    raise ValueError("The field has no image: give 'prepared', 'source' or 'path'.")
                series = open_image_series(path, force_bioformats, False, log)
                if len(series) == 0:
                    raise ValueError("Could not open the image file.")
                opened = series[0][1]
                opened.setTitle(to_text(os.path.basename(path)))
                source = opened
            validate_channel_configs(source, [config])
            prepare_config = dict(config)
            if needs_base:
                prepare_config["median_radius"] = 0
                prepare_config["gaussian_sigma"] = 0
            image = engine.prepare_single_channel(source, prepare_config, int(field.get("frame", 1)), calibration)
            work["owned"].append(image)
            work["base" if needs_base else "prepared"] = image
    finally:
        if opened is not None:
            safe_close_image(opened)
            try:
                opened.flush()
            except ANY_ERROR:
                pass
    try:
        reference = work["prepared"] if work["base"] is None else work["base"]
        width = reference.getWidth()
        height = reference.getHeight()
        rois = field.get("rois")
        if rois is None or len(rois) == 0:
            rois = [whole_image_roi(reference)]
        else:
            rois = validate_parent_rois([clone_roi(roi) for roi in rois], width, height)
        points = field.get("points")
        if points is None:
            points = annotation_points_from_rois(field.get("point_rois") or [])
        points = [(int(point[0]), int(point[1])) for point in points]
        outside = 0
        for point in points:
            inside_any = False
            for roi in rois:
                if roi.contains(point[0], point[1]):
                    inside_any = True
                    break
            if not inside_any:
                outside += 1
        work.update({
            "width": width, "height": height, "rois": rois, "points": points,
            # Exactly the rule of validate_objects_against_points.
            "roi_points": [[point for point in points if roi.contains(point[0], point[1])] for roi in rois],
            "outside": outside, "roi_masks": {}
        })
        work["inside"] = sum([len(part) for part in work["roi_points"]])
    except ANY_ERROR:
        tune_close_field(work)
        raise
    return work


def tune_close_field(work):
    if work is None:
        return
    for imp in work.get("owned", []):
        safe_close_image(imp)
    work["owned"] = []
    work["roi_masks"] = {}


def tune_level_keys(combos):
    """For every combination, the values that decide which intermediate
    result it can share with other combinations, per level."""
    keys = {"prep": [], "mask": [], "particle": [], "size": []}
    for combo in combos:
        variant = combo["config"]
        for level in keys:
            parts = []
            for spec in TUNE_PARAMETERS:
                if spec[1] != level:
                    continue
                value = variant[spec[0]]
                if spec[0] == "clahe_max_slope" and not variant["clahe"]:
                    value = None
                parts.append(value)
            keys[level].append(tuple(parts))
    return keys


def tune_group_by_key(indices, keys):
    """Split indices into groups of equal key, in order of first appearance."""
    order = []
    members = {}
    for index in indices:
        key = keys[index]
        if key not in members:
            members[key] = []
            order.append(key)
        members[key].append(index)
    return [(key, members[key]) for key in order]


def tune_field_roi_masks(work, roi_index):
    if roi_index not in work["roi_masks"]:
        work["roi_masks"][roi_index] = roi_masks(work["rois"][roi_index], work["width"], work["height"])
    return work["roi_masks"][roi_index]


def tune_match_counts(work, object_roi_lists):
    """Pool the match of every parent ROI of one field. abs_count_error is
    summed ROI by ROI, |detected - marks| of each parent ROI, because the
    parent ROI is the unit whose count is reported: a count too high in one
    ROI and too low in another of the same image must not cancel."""
    totals = {"abs_count_error": 0}
    for key in TUNE_COUNT_KEYS:
        totals[key] = 0
    for roi_index in range(len(work["rois"])):
        match = match_points_to_objects(
            object_roi_lists[roi_index], work["roi_points"][roi_index], work["width"], work["height"]
        )
        for key in TUNE_COUNT_KEYS:
            totals[key] += match[key]
        totals["abs_count_error"] += abs(match["detected"] - match["manual"])
    return totals


def tune_direct_counts(engine, work, final_mask, variant):
    """Counts of one configuration, straight through the engine."""
    lists = []
    for roi_index in range(len(work["rois"])):
        objects, info = engine.analyze_objects(
            final_mask, work["rois"][roi_index], variant, tune_field_roi_masks(work, roi_index)
        )
        lists.append([obj["roi"] for obj in objects])
    return tune_match_counts(work, lists)


def tune_score_sizes(engine, work, final_mask, combos, size_groups, state):
    """Score combinations that share one mask and differ only in the size
    limits. Returns one count dictionary per size group.

    The particle analysis accepts an object when its pixel count lies
    between the two limits, object by object. So the objects for strict
    limits are a subset of the objects for the loosest limits of the group:
    the analysis runs once with the loosest limits and each combination
    keeps the objects whose pixel count passes its own limits. The first
    time this shortcut is used on a field it is checked against a direct
    analysis; if the two ever disagree the shortcut is switched off."""
    first = combos[size_groups[0][1][0]]["config"]
    if first["composite_rois"]:
        # One combined ROI for all objects: there is no object list to filter.
        state["fast"] = False
    if not state["fast"] or len(size_groups) < 2:
        return [tune_direct_counts(engine, work, final_mask, combos[group[1][0]]["config"]) for group in size_groups]
    calibration = final_mask.getCalibration()
    pixel_area = calibration.pixelWidth * calibration.pixelHeight
    minimums = [combos[group[1][0]]["config"]["min_size"] for group in size_groups]
    maximums = [combos[group[1][0]]["config"]["max_size"] for group in size_groups]
    loose = dict(first)
    loose["min_size"] = min(minimums)
    loose["max_size"] = -1 if -1 in maximums else max(maximums)
    object_lists = []
    pixel_lists = []
    for roi_index in range(len(work["rois"])):
        objects, info = engine.analyze_objects(
            final_mask, work["rois"][roi_index], loose, tune_field_roi_masks(work, roi_index)
        )
        pixels = [obj["area"] / pixel_area for obj in objects]
        for value in pixels:
            if value != value:
                state["fast"] = False
                return tune_score_sizes(engine, work, final_mask, combos, size_groups, state)
        object_lists.append(objects)
        pixel_lists.append([int(round(value)) for value in pixels])
    results = []
    for group in size_groups:
        engine.check_cancelled()
        variant = combos[group[1][0]]["config"]
        min_px, max_px = engine.particle_size_pixels(variant, calibration)
        lists = []
        for roi_index in range(len(work["rois"])):
            objects = object_lists[roi_index]
            pixels = pixel_lists[roi_index]
            lists.append([objects[k]["roi"] for k in range(len(objects)) if min_px <= pixels[k] <= max_px])
        counts = tune_match_counts(work, lists)
        is_loose = variant["min_size"] == loose["min_size"] and variant["max_size"] == loose["max_size"]
        if not state["checked"] and not is_loose:
            state["checked"] = True
            if tune_direct_counts(engine, work, final_mask, variant) != counts:
                state["fast"] = False
                engine.log("Tuning: the size-limit shortcut disagreed with a direct analysis; it is switched off.")
                return tune_score_sizes(engine, work, final_mask, combos, size_groups, state)
        results.append(counts)
    return results


def tune_score_field(engine, work, combos, fast_sizes=True, tick=None):
    """Score every combination on one field.

    Returns one entry per combination: a dictionary of pooled counts
    (manual, detected, true_positive, false_positive, missed, merged), or
    {"error": text} when the combination cannot be run on this field (for
    example a manual threshold above the brightest pixel).

    Combinations are grouped so that each prepared image, each mask and
    each particle analysis is computed once and shared."""
    results = [None] * len(combos)
    usable = []
    for index in range(len(combos)):
        if combos[index]["status"] != "":
            results[index] = {"error": combos[index]["status"]}
        else:
            usable.append(index)
    keys = tune_level_keys(combos)
    state = {"fast": bool(fast_sizes), "checked": False}

    def assign(members, value):
        for index in members:
            results[index] = value
        if tick is not None:
            tick(len(members))

    if tick is not None:
        tick(len(combos) - len(usable))
    for prep_group in tune_group_by_key(usable, keys["prep"]):
        engine.check_cancelled()
        variant = combos[prep_group[1][0]]["config"]
        temporary = None
        if work["base"] is None:
            prepared = work["prepared"]
        elif variant["median_radius"] > 0 or variant["gaussian_sigma"] > 0:
            temporary = tune_apply_filters(work["base"], variant)
            prepared = temporary
        else:
            prepared = work["base"]
        try:
            for mask_group in tune_group_by_key(prep_group[1], keys["mask"]):
                engine.check_cancelled()
                variant = combos[mask_group[1][0]]["config"]
                try:
                    mask_data = engine.build_masks(prepared, variant, "CellQuant_tuning")
                except (ValueError, RuntimeError) as exc:
                    if TUNE_EMPTY_RANGE_TEXT in to_text(exc):
                        # No pixel value of this field lies in the threshold
                        # range (for example a threshold above the brightest
                        # pixel of a negative control): the true answer is
                        # an empty mask, so nothing is detected.
                        empty = tune_match_counts(work, [[] for roi in work["rois"]])
                        empty["empty_range"] = True
                        assign(mask_group[1], empty)
                    else:
                        assign(mask_group[1], {"error": to_text(exc)})
                    continue
                safe_close_image(mask_data["raw"])
                try:
                    for particle_group in tune_group_by_key(mask_group[1], keys["particle"]):
                        engine.check_cancelled()
                        size_groups = tune_group_by_key(particle_group[1], keys["size"])
                        counts = tune_score_sizes(engine, work, mask_data["final"], combos, size_groups, state)
                        for position in range(len(size_groups)):
                            if mask_data.get("empty_range"):
                                # An engine that returns an empty mask for an
                                # empty threshold range says so; keep the note.
                                counts[position]["empty_range"] = True
                            assign(size_groups[position][1], counts[position])
                finally:
                    safe_close_image(mask_data["final"])
        finally:
            safe_close_image(temporary)
    return results


# -----------------------------------------------------------------------------
# Scores, selection, and cross-validation (arithmetic on stored counts)
# -----------------------------------------------------------------------------

def tune_pool(scores):
    """Sum the counts of several fields. Returns None when one of them could
    not be scored. abs_count_error is the sum over every parent ROI of
    every field of |detected - manual| (see tune_match_counts), so errors
    of opposite sign do not cancel."""
    totals = {"abs_count_error": 0, "fields": 0}
    for key in TUNE_COUNT_KEYS:
        totals[key] = 0
    for score in scores:
        if score is None or "error" in score:
            return None
        for key in TUNE_COUNT_KEYS:
            totals[key] += score[key]
        totals["abs_count_error"] += score["abs_count_error"]
        totals["fields"] += 1
    return totals


def tune_all_missed(marks):
    """Counts of a field on which the chosen settings cannot run at all:
    nothing is detected, every mark is missed."""
    return {"manual": marks, "detected": 0, "true_positive": 0, "false_positive": 0, "missed": marks,
            "merged": 0, "abs_count_error": marks}


def tune_metrics(totals):
    """Precision, recall and F1 from pooled counts, with the formulas of
    validate_objects_against_points, and the count errors."""
    precision = safe_ratio(totals["true_positive"], totals["true_positive"] + totals["false_positive"])
    recall = safe_ratio(totals["true_positive"], totals["manual"])
    f1 = float("nan")
    if precision == precision and recall == recall and (precision + recall) > 0:
        f1 = 2.0 * precision * recall / (precision + recall)
    difference = totals["detected"] - totals["manual"]
    return {
        "precision": precision, "recall": recall, "f1": f1,
        "count_difference": difference,
        "count_error_percent": safe_ratio(difference, totals["manual"]) * 100.0,
        "abs_count_error_percent": safe_ratio(totals["abs_count_error"], totals["manual"]) * 100.0
    }


def tune_metric_cells(totals, prefix=""):
    """Table cells for pooled counts; empty when there are none."""
    if totals is None:
        return {}
    metrics = tune_metrics(totals)
    return {
        prefix + "Manual_Count": totals["manual"],
        prefix + "Detected_Count": totals["detected"],
        prefix + "True_Positive_Objects": totals["true_positive"],
        prefix + "False_Positive_Objects": totals["false_positive"],
        prefix + "Missed_Annotations": totals["missed"],
        prefix + "Merged_Annotations": totals["merged"],
        prefix + "Precision": metrics["precision"],
        prefix + "Recall": metrics["recall"],
        prefix + "F1": metrics["f1"],
        prefix + "Count_Difference": metrics["count_difference"],
        prefix + "Count_Error_Percent": metrics["count_error_percent"],
        prefix + "Abs_Count_Error": totals["abs_count_error"],
        prefix + "Abs_Count_Error_Percent": metrics["abs_count_error_percent"]
    }


def tune_objective_value(totals, objective):
    """The number that is maximized. F1 is computed from the counts with a
    single division, 2 TP / (TP + FP + marks), so that equal scores compare
    as exactly equal and ties are broken by the stated rule, not by rounding."""
    if objective == TUNE_OBJECTIVE_COUNT:
        return -float(totals["abs_count_error"])
    denominator = totals["true_positive"] + totals["false_positive"] + totals["manual"]
    if denominator == 0:
        return 0.0
    return 2.0 * totals["true_positive"] / float(denominator)


def tune_rank(combos, field_scores, field_indices, objective):
    """Order the combinations that could be scored on every listed field,
    best first. Each entry is (sort key, combination index, pooled counts).
    The key encodes the whole selection rule: objective, then fewest
    settings changed, then smallest relative change, then grid order."""
    ranked = []
    for index in range(len(combos)):
        combo = combos[index]
        if combo["baseline"] or combo["status"] != "":
            continue
        totals = tune_pool([field_scores[field_index][index] for field_index in field_indices])
        if totals is None:
            continue
        key = (tune_objective_value(totals, objective), -combo["changes"], -combo["distance"], -index)
        ranked.append((key, index, totals))
    ranked.sort(reverse=True)
    return ranked


def tune_validation_design(groups):
    labels = []
    for group in groups:
        if group not in labels:
            labels.append(group)
    if len(labels) < 2:
        return TUNE_NO_HELD_OUT if len(groups) == 1 else TUNE_NO_HELD_OUT_GROUP
    if len(labels) < len(groups):
        return "leave-one-group-out cross-validation (" + str(len(labels)) + " groups)"
    if len(labels) == 2:
        return "two-way split (each field held out once)"
    return "leave-one-field-out cross-validation"


def tune_cross_validate(combos, field_scores, groups, marks, objective):
    """Leave each group of fields out once: tune on the others, score the
    left-out fields with the settings found. Pure arithmetic on the stored
    counts. Returns one dictionary per fold; selected is None (with a note)
    when the fold could not be tuned."""
    labels = []
    for group in groups:
        if group not in labels:
            labels.append(group)
    folds = []
    if len(labels) < 2:
        return folds
    baseline_index = len(combos) - 1
    for label in labels:
        test = [index for index in range(len(groups)) if groups[index] == label]
        train = [index for index in range(len(groups)) if groups[index] != label]
        fold = {
            "group": label, "test": test, "train": train, "selected": None, "note": "",
            "train_marks": sum([marks[index] for index in train]), "train_totals": None,
            "held_out": [], "baseline": []
        }
        folds.append(fold)
        if fold["train_marks"] == 0:
            fold["note"] = "Not tuned: the training fields hold no marks."
            continue
        ranked = tune_rank(combos, field_scores, train, objective)
        if len(ranked) == 0:
            fold["note"] = "Not tuned: no combination could be run on every training field."
            continue
        fold["selected"] = ranked[0][1]
        fold["train_totals"] = ranked[0][2]
        notes = []
        for field_index in test:
            score = field_scores[field_index][fold["selected"]]
            if "error" in score:
                notes.append("The selected settings could not be run on this field (" + score["error"] + "); every mark is counted as missed.")
                score = tune_all_missed(marks[field_index])
            fold["held_out"].append(score)
            score = field_scores[field_index][baseline_index]
            if "error" in score:
                notes.append("The starting settings could not be run on this field; every mark is counted as missed.")
                score = tune_all_missed(marks[field_index])
            fold["baseline"].append(score)
        fold["note"] = " ".join(notes)
    return folds


# -----------------------------------------------------------------------------
# A tuning run
# -----------------------------------------------------------------------------

def tune_value_cells(config):
    cells = {}
    for key in TUNE_PARAMETER_KEYS:
        value = config[key]
        cells["P_" + key] = bool_text(value) if isinstance(value, bool) else value
    return cells


def tune_tables(run_id, config, space, combos, names, groups, marks, field_scores, folds, objective):
    """Build the rows of tuning_grid.csv, tuning_field_scores.csv,
    tuning_folds.csv and tuning_stability.csv, and the selection results."""
    ident = {
        "Run_ID": run_id, "Application_Version": APP_VERSION,
        "Channel": config["channel"], "Channel_Name": config["name"]
    }
    everything = range(len(names))
    ranked = tune_rank(combos, field_scores, everything, objective)
    final = ranked[0][1] if len(ranked) > 0 else None
    rank_of = {}
    for position in range(len(ranked)):
        rank_of[ranked[position][1]] = position + 1
    scored_folds = [fold for fold in folds if fold["selected"] is not None]

    grid_rows = []
    score_rows = []
    for index in range(len(combos) - 1):
        combo = combos[index]
        totals = tune_pool([field_scores[field_index][index] for field_index in everything])
        row = dict(ident)
        row.update(tune_value_cells(combo["config"]))
        row.update(tune_metric_cells(totals))
        status = "OK"
        if combo["status"] != "":
            status = combo["status"]
        elif totals is None:
            failed = [field_index for field_index in everything if "error" in field_scores[field_index][index]]
            status = (
                "FAILED on " + str(len(failed)) + " field(s), first " + names[failed[0]] + ": " +
                field_scores[failed[0]][index]["error"]
            )
        else:
            empty = [field_index for field_index in everything if field_scores[field_index][index].get("empty_range")]
            if len(empty) > 0:
                status = (
                    "OK: nothing detected on " + str(len(empty)) + " field(s), first " + names[empty[0]] +
                    ", because no pixel value lies in the threshold range"
                )
        row.update({
            "Combination": combo["number"], "Settings": tune_settings_text(space, combo["values"]),
            "Changed_Settings": combo["changes"], "Relative_Change": combo["distance"],
            "Is_Start": bool_text(combo["changes"] == 0), "Status": status,
            "Fields": len(names), "Objective": objective,
            "Objective_Value": tune_objective_value(totals, objective) if totals is not None else "",
            "Rank": rank_of.get(index, ""), "Selected_Final": bool_text(index == final),
            "Folds_Selected": len([fold for fold in scored_folds if fold["selected"] == index])
        })
        grid_rows.append(row)
        for field_index in everything:
            score = field_scores[field_index][index]
            field_status = "OK"
            if "error" in score:
                field_status = "FAILED: " + score["error"]
            elif score.get("empty_range"):
                field_status = TUNE_EMPTY_RANGE_STATUS
            score_row = {
                "Run_ID": run_id, "Application_Version": APP_VERSION, "Combination": combo["number"],
                "Field_Index": field_index + 1, "Field": names[field_index], "Group": groups[field_index],
                "Status": field_status
            }
            if "error" not in score:
                score_row.update(tune_metric_cells(tune_pool([score])))
            score_rows.append(score_row)

    fold_rows = []
    for fold_index in range(len(folds)):
        fold = folds[fold_index]
        for position in range(len(fold["test"])):
            field_index = fold["test"][position]
            row = dict(ident)
            row.update({
                "Fold": fold_index + 1, "Held_Out_Group": fold["group"], "Held_Out_Field": names[field_index],
                "Training_Fields": len(fold["train"]), "Training_Marks": fold["train_marks"], "Note": fold["note"]
            })
            if fold["selected"] is not None:
                combo = combos[fold["selected"]]
                row.update(tune_value_cells(combo["config"]))
                row.update(tune_metric_cells(tune_pool([fold["held_out"][position]])))
                row.update(tune_metric_cells(tune_pool([fold["baseline"][position]]), "Baseline_"))
                row["Selected_Combination"] = combo["number"]
                row["Selected_Settings"] = tune_settings_text(space, combo["values"])
                row["Training_F1_Optimistic"] = tune_metrics(fold["train_totals"])["f1"]
            fold_rows.append(row)

    stability_rows = []
    for entry in space:
        for value in entry[1]:
            chosen = len([fold for fold in scored_folds if combos[fold["selected"]]["values"][entry[0]] == value])
            row = dict(ident)
            row.update({
                "Parameter": entry[0], "Value": tune_text_value(value),
                "Is_Start_Value": bool_text(value == config[entry[0]]),
                "Folds_Selected": chosen, "Folds_Scored": len(scored_folds),
                "Fraction_Of_Folds": safe_ratio(chosen, len(scored_folds)),
                "Selected_Final": bool_text(final is not None and combos[final]["values"][entry[0]] == value)
            })
            stability_rows.append(row)

    return {
        "ident": ident, "ranked": ranked, "final": final, "grid_rows": grid_rows,
        "score_rows": score_rows, "fold_rows": fold_rows, "stability_rows": stability_rows
    }


def tune_summary_row(tables, combos, field_scores, folds, marks, objective):
    """The one row of tuning_summary.csv (without the run bookkeeping)."""
    row = dict(tables["ident"])
    baseline_index = len(combos) - 1
    everything = range(len(marks))
    held_out = []
    held_out_baseline = []
    field_f1 = []
    for fold in folds:
        if fold["selected"] is None:
            continue
        held_out += fold["held_out"]
        held_out_baseline += fold["baseline"]
        for score in fold["held_out"]:
            value = tune_metrics(tune_pool([score]))["f1"]
            if value == value:
                field_f1.append(value)
    held_out_totals = tune_pool(held_out) if len(held_out) > 0 else None
    resubstitution = tables["ranked"][0][2]
    all_baseline = tune_pool([
        tune_all_missed(marks[index]) if "error" in field_scores[index][baseline_index]
        else field_scores[index][baseline_index] for index in everything
    ])
    row.update(tune_metric_cells(held_out_totals, "HeldOut_"))
    row.update(tune_metric_cells(tune_pool(held_out_baseline) if len(held_out) > 0 else None, "HeldOut_Baseline_"))
    row.update(tune_metric_cells(resubstitution, "Resubstitution_Optimistic_"))
    row.update(tune_metric_cells(all_baseline, "AllFields_Baseline_"))
    row.update({
        "Objective": objective, "Folds": len(folds),
        "Folds_Scored": len([fold for fold in folds if fold["selected"] is not None]),
        "HeldOut_Fields": len(held_out),
        "HeldOut_F1_Field_Min": min(field_f1) if len(field_f1) > 0 else "",
        "HeldOut_F1_Field_Max": max(field_f1) if len(field_f1) > 0 else ""
    })
    return {
        "row": row, "held_out": held_out_totals, "resubstitution": resubstitution,
        "baseline": all_baseline,
        "held_out_baseline": tune_pool(held_out_baseline) if len(held_out) > 0 else None
    }


def tune_score_all_fields(engine, fields, config, combos, run, options, progress):
    """Open, score and close every field in turn. Fills run["field_rows"]
    and returns the lists that describe the fields that could be used."""
    needs_base = len([entry for entry in run["space"] if entry[0] in ["gaussian_sigma", "median_radius"]]) > 0
    used = {"names": [], "groups": [], "marks": [], "outside": 0, "scores": []}
    for field_index in range(len(fields)):
        field = fields[field_index]
        name = to_text(field.get("name", "Field_" + str(field_index + 1).zfill(3)))
        row = {
            "Run_ID": run["run_id"], "Application_Version": APP_VERSION, "Field_Index": "", "Field": name,
            "Group": to_text(field.get("group", name)), "Image_Path": to_text(field.get("path", "")),
            "ROI_Source": to_text(field.get("roi_source", "")),
            "Points_Source": to_text(field.get("points_source", "")), "Message": ""
        }
        run["field_rows"].append(row)
        started = time.time()
        offset = field_index * len(combos)
        local = {"done": 0}

        def tick(count, _offset=offset, _local=local, _name=name):
            _local["done"] += count
            if progress is not None:
                progress(_offset + _local["done"], len(fields) * len(combos), _name)

        work = None
        try:
            engine.check_cancelled()
            engine.show_status("Tuning: field " + str(field_index + 1) + " of " + str(len(fields)))
            work = tune_open_field(
                engine, field, config, needs_base, run["calibration"],
                parse_bool(options.get("force_bioformats"), False), run["log"]
            )
            row.update({
                "Width": work["width"], "Height": work["height"], "Parent_ROI_Count": len(work["rois"]),
                "Marks_Total": len(work["points"]), "Marks_Inside_ROIs": work["inside"],
                "Marks_Outside_ROIs": work["outside"]
            })
            scores = tune_score_field(engine, work, combos, parse_bool(options.get("fast_size_filter"), True), tick)
            used["names"].append(name)
            used["groups"].append(row["Group"])
            used["marks"].append(work["inside"])
            used["outside"] += work["outside"]
            used["scores"].append(scores)
            row["Status"] = "USED"
            row["Field_Index"] = len(used["names"])
            run["log"](
                "Field " + name + ": " + str(work["inside"]) + " mark(s) inside " + str(len(work["rois"])) +
                " parent ROI(s)" + (", " + str(work["outside"]) + " outside every ROI (ignored)" if work["outside"] > 0 else "") +
                ", " + IJ.d2s(time.time() - started, 1) + " s."
            )
        except AnalysisCancelled:
            row["Status"] = "CANCELLED"
            raise
        except ANY_ERROR as exc:
            row["Status"] = "FAILED"
            row["Message"] = to_text(exc)
            run["log"]("Field " + name + " FAILED and is not used: " + to_text(exc))
            if not isinstance(exc, (ValueError, RuntimeError, IOError)):
                # An unexpected error: keep the details for a bug report.
                try:
                    IJ.log(traceback.format_exc())
                except ANY_ERROR:
                    pass
        finally:
            row["Seconds"] = round(time.time() - started, 2)
            tune_close_field(work)
    return used


def run_tuning(engine, fields, config, all_configs, global_options, search_space, output_root,
               options=None, progress=None):
    """Tune settings of one channel against manual marks and write a run
    folder CellQuant_Tuning_<timestamp> under output_root.

    engine: an AnalysisEngine (its log and cancel functions are used).
    fields: list of field dictionaries (see tune_open_field), for example
            from discover_annotated_fields.
    config: the starting configuration of the channel to tune.
    all_configs: every channel of the starting profile; the tuned profile
            holds them unchanged except the tuned channel.
    search_space: dictionary setting -> list of values, or None for
            default_tune_search_space(config).
    options: objective (TUNE_OBJECTIVES), max_combinations, factors (for
            the default search space), fast_size_filter, force_bioformats,
            unused (report rows of images without annotation).
    progress: function (done, total, field name), counted in field x
            combination evaluations.

    Returns a dictionary; "cancelled" is True (and no profile is written)
    when the run was cancelled."""
    options = options or {}
    objective = to_text(options.get("objective", TUNE_OBJECTIVE_F1))
    if objective not in TUNE_OBJECTIVES:
        raise ValueError("Unknown tuning objective: " + objective)
    if fields is None or len(fields) == 0:
        raise ValueError("No annotated field was given.")
    # Read the configuration the way a saved profile is read, so that the
    # settings scored here are exactly the settings the tuned profile holds.
    config = normalize_channel_config(config["channel"], config)
    if search_space is None:
        search_space = default_tune_search_space(config, options.get("factors"))
    space = normalize_tune_search_space(config, search_space, options.get("max_combinations"))
    combos = tune_build_combinations(config, space)

    run_dir = unique_dir(os.path.abspath(to_text(output_root)), "CellQuant_Tuning_" + timestamp_text())
    log_lines = []

    def tune_log(message):
        log_lines.append("[" + time.strftime("%H:%M:%S") + "] " + to_text(message))
        engine.log(message)

    run = {
        "run_id": os.path.basename(run_dir), "space": space, "log": tune_log, "field_rows": [],
        "calibration": calibration_override_from_options(global_options)
    }
    result = {"run_dir": run_dir, "run_id": run["run_id"], "cancelled": False, "profile_path": None}
    unused = [dict(row) for row in options.get("unused", [])]
    for row in unused:
        row["Run_ID"] = run["run_id"]
        row["Application_Version"] = APP_VERSION
    started = time.time()
    try:
        write_settings_properties(os.path.join(run_dir, "settings.properties"), global_options, all_configs)
        tune_log("Tuning run " + run["run_id"] + ": C" + str(config["channel"]) + " (" + config["name"] + "), objective " + objective)
        tune_log(APP_NAME + " " + APP_VERSION + " on ImageJ " + IMAGEJ_VERSION)
        tune_log("Search space (" + str(len(combos) - 1) + " combinations): " + describe_tune_search_space(space))
        used = tune_score_all_fields(engine, fields, config, combos, run, options, progress)
        seconds = time.time() - started
        names = used["names"]
        if len(names) == 0:
            raise ValueError("No annotated field could be used. " + to_text(run["field_rows"][0].get("Message", "")))
        if sum(used["marks"]) == 0:
            raise ValueError("The annotated fields hold no mark inside a parent ROI, so there is nothing to tune against.")
        folds = tune_cross_validate(combos, used["scores"], used["groups"], used["marks"], objective)
        tables = tune_tables(run["run_id"], config, space, combos, names, used["groups"], used["marks"], used["scores"], folds, objective)
        if tables["final"] is None:
            raise ValueError(
                "No combination could be run on every field. First problem: " + tables["grid_rows"][0]["Status"]
            )
        unusable = (len(combos) - 1) - len(tables["ranked"])
        if unusable > 0:
            tune_log(
                str(unusable) + " of " + str(len(combos) - 1) + " combination(s) could not be run on every field "
                "and were left out of the selection (see the Status column of tuning_grid.csv)."
            )
        final = combos[tables["final"]]
        tuned_config = dict(final["config"])
        tuned_all = []
        replaced = False
        for other in all_configs:
            if other["channel"] == config["channel"]:
                tuned_all.append(tuned_config)
                replaced = True
            else:
                tuned_all.append(other)
        if not replaced:
            tuned_all.append(tuned_config)
        profile_path = os.path.join(run_dir, "tuned_profile.properties")
        write_settings_properties(profile_path, global_options, tuned_all)

        summary = tune_summary_row(tables, combos, used["scores"], folds, used["marks"], objective)
        design = tune_validation_design(used["groups"])
        evaluations = len(names) * len(combos)
        note = "HeldOut_ columns are the validation result. Resubstitution_Optimistic_ columns score the final settings on the fields they were tuned on."
        why = design.split(": ", 1)[-1]
        if summary["held_out"] is None and len(set(used["groups"])) < 2 and len(names) > 1:
            note = (
                "No held-out estimate: all " + str(len(names)) + " fields belong to one group, so none can be held out. "
                "Put the fields of each organoid or section in its own subfolder and hold out by subfolder."
            )
        elif summary["held_out"] is None:
            note = "No held-out estimate: " + why + ". Annotate at least two fields (three or more recommended) for a validation result."
        summary["row"].update({
            "Validation_Design": design, "Fields": len(names),
            "Fields_Failed": len([row for row in run["field_rows"] if row.get("Status") == "FAILED"]),
            "Images_Without_Annotation": len(unused),
            "Groups": len(set(used["groups"])), "Marks_Inside_ROIs": sum(used["marks"]),
            "Marks_Outside_ROIs": used["outside"], "Search_Space": describe_tune_search_space(space),
            "Combinations": len(combos) - 1, "Combinations_Usable": len(tables["ranked"]),
            "Start_Settings": tune_settings_text(space, combos[-1]["values"]),
            "Tuned_Combination": final["number"], "Tuned_Settings": tune_settings_text(space, final["values"]),
            "Tuned_Profile": "tuned_profile.properties", "Evaluations": evaluations,
            "Seconds": round(seconds, 2), "Evaluations_Per_Second": safe_ratio(evaluations, seconds), "Note": note
        })
        write_csv(os.path.join(run_dir, "tuning_grid.csv"), TUNE_GRID_HEADERS, tables["grid_rows"])
        write_csv(os.path.join(run_dir, "tuning_field_scores.csv"), TUNE_SCORE_HEADERS, tables["score_rows"])
        write_csv(os.path.join(run_dir, "tuning_folds.csv"), TUNE_FOLD_HEADERS, tables["fold_rows"])
        write_csv(os.path.join(run_dir, "tuning_stability.csv"), TUNE_STABILITY_HEADERS, tables["stability_rows"])
        write_csv(os.path.join(run_dir, "tuning_summary.csv"), TUNE_SUMMARY_HEADERS, [summary["row"]])

        tune_log(
            str(evaluations) + " evaluations (field x combination) in " + IJ.d2s(seconds, 1) + " s: " +
            IJ.d2s(safe_ratio(evaluations, seconds), 2) + " per second."
        )
        tune_log("Tuned settings (from all " + str(len(names)) + " field(s)): " + tune_settings_text(space, final["values"]))
        if summary["held_out"] is not None:
            tune_log(
                "HELD-OUT result, " + design + ": F1 " + IJ.d2s(tune_metrics(summary["held_out"])["f1"], 3) +
                ", precision " + IJ.d2s(tune_metrics(summary["held_out"])["precision"], 3) +
                ", recall " + IJ.d2s(tune_metrics(summary["held_out"])["recall"], 3) +
                ". Starting settings on the same fields: F1 " +
                IJ.d2s(tune_metrics(summary["held_out_baseline"])["f1"], 3) + "."
            )
        else:
            tune_log("NO HELD-OUT RESULT: " + why + ".")
        tune_log(
            "Final settings on the fields they were tuned on (optimistic, not a validation result): F1 " +
            IJ.d2s(tune_metrics(summary["resubstitution"])["f1"], 3) + "."
        )
        tune_log("Saved tuning run: " + to_text(run_dir))
        result.update({
            "profile_path": profile_path, "tuned_config": tuned_config, "tuned_configs": tuned_all,
            "tuned_values": dict(final["values"]), "design": design, "objective": objective,
            "space": space, "folds": folds, "field_names": names, "field_scores": used["scores"],
            "combinations": combos, "held_out": summary["held_out"],
            "held_out_baseline": summary["held_out_baseline"],
            "resubstitution": summary["resubstitution"], "baseline": summary["baseline"],
            "summary_row": summary["row"], "grid_rows": tables["grid_rows"], "fold_rows": tables["fold_rows"],
            "stability_rows": tables["stability_rows"], "score_rows": tables["score_rows"],
            "seconds": seconds, "evaluations_per_second": safe_ratio(evaluations, seconds)
        })
    except AnalysisCancelled:
        result["cancelled"] = True
        tune_log("Tuning cancelled by user. No tuned profile was written.")
    except ANY_ERROR as exc:
        tune_log("Tuning stopped: " + to_text(exc))
        raise
    finally:
        result["field_rows"] = run["field_rows"] + unused
        write_csv(os.path.join(run_dir, "tuning_fields.csv"), TUNE_FIELD_HEADERS, result["field_rows"])
        write_text_file(os.path.join(run_dir, "tuning_log.txt"), [APP_NAME + " " + APP_VERSION] + log_lines)
        if progress is not None:
            progress(len(fields) * len(combos), len(fields) * len(combos), "")
    return result


def tune_config_from_options(global_options, configs):
    """The channel to tune: the channel number in the tune_channel setting,
    or with 0 the first enabled channel that detects objects."""
    wanted = parse_int_value(global_options.get("tune_channel", "0"), "Channel to tune", 0, MAX_CHANNELS)
    for config in configs:
        if (wanted == 0 and needs_objects(config)) or config["channel"] == wanted:
            return config
    if wanted == 0:
        raise ValueError("No enabled channel detects objects, so there is nothing to tune.")
    raise ValueError("C" + str(wanted) + " is not an enabled channel.")


def run_tuning_folder(engine, folder, output_root, configs, all_configs, global_options,
                      folder_options=None, progress=None):
    """Tune from an annotated folder using the tune_* settings of a profile.

    folder_options (all optional): extensions, recursive (default: the batch
    settings of the profile), roi_dir, points_dir, force_bioformats, and
    search_space (a dictionary that replaces the tune_search_space text).
    With tune_group_by = Subfolder, subfolders are always searched."""
    folder_options = folder_options or {}
    config = tune_config_from_options(global_options, configs)
    group_by = to_text(global_options.get("tune_group_by", TUNE_GROUP_IMAGE))
    recursive = parse_bool(folder_options.get("recursive", global_options.get("batch_recursive")), False)
    if group_by == TUNE_GROUP_FOLDER and not recursive:
        # Holding out by subfolder only makes sense when the subfolders are read.
        recursive = True
        engine.log("Tuning: fields are held out by subfolder, so subfolders are searched too.")
    found = discover_annotated_fields(
        folder,
        folder_options.get("extensions", global_options.get("batch_extensions")),
        recursive, folder_options.get("roi_dir"), folder_options.get("points_dir"), group_by
    )
    for row in found["report"]:
        if row["Status"] != "ANNOTATED":
            engine.log("Not used for tuning: " + row["Field"] + " - " + row["Message"])
    if len(found["fields"]) == 0:
        raise ValueError(
            "No annotated image was found in " + to_text(folder) + ". Save the marks of each image as "
            "<image name>_points.zip or <image name>_points.roi beside it."
        )
    search_space = folder_options.get("search_space")
    if search_space is None and to_text(global_options.get("tune_search_space", "")).strip() != "":
        search_space = parse_tune_search_space(global_options.get("tune_search_space"))
    options = {
        "objective": global_options.get("tune_objective", TUNE_OBJECTIVE_F1),
        "max_combinations": parse_int_value(
            global_options.get("tune_max_combinations", TUNE_MAX_COMBINATIONS), "Largest number of combinations", 1
        ),
        "factors": global_options.get("tune_factors"),
        "force_bioformats": folder_options.get("force_bioformats", global_options.get("batch_force_bioformats")),
        "unused": [row for row in found["report"] if row["Status"] != "ANNOTATED"]
    }
    return run_tuning(engine, found["fields"], config, all_configs, global_options, search_space,
                      output_root, options, progress)


register_global_settings(TUNE_SETTING_SPECS)
register_settings_form(
    "Tune settings from point marks",
    "Searches settings of one channel for the best agreement with manual point marks "
    "(<image name>_points.zip beside each image) and scores the result on fields left out of the search.",
    [
        ("tune_channel", "Channel to tune (0 = first channel that detects objects)", "text", None),
        ("tune_objective", "What to optimise", "choice", list(TUNE_OBJECTIVES)),
        ("tune_search_space", "Settings to search, e.g. manual_low = 800, 1000; min_size = 20, 40 (empty = automatic)", "text", None),
        ("tune_factors", "Factors of the automatic search", "text", None),
        ("tune_max_combinations", "Largest number of combinations", "text", None),
        ("tune_group_by", "Fields that are held out together", "choice", list(TUNE_GROUPINGS))
    ]
)
register_engine_module("cellquanttune.py", "2.3.0")

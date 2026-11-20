# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach

# Cell Quant Studio Universal - depth profile from the ROI outline (optional module)
#
# Brain organoids are layered from the surface inward. When the parent ROI is
# the outline of the tissue, this module sorts every pixel and every detected
# object of the ROI by its distance from that outline ("depth") and reports,
# for each depth bin and each analysed channel: the area of the bin, the mean
# prepared intensity, the marker-positive area fraction, and the number and
# density of objects. It needs no annotation. Measuring along the distance
# from the tissue edge is an established idea (for example SCOUT for cleared
# whole organoids); this is a convenience for 2-D sections, not a new method.
#
# Exact definitions (all of them are checked by tests/test_depth.py):
#   - Depth of a pixel = distance from the centre of the pixel to the outline
#     of the ROI = (ImageJ's Euclidean distance map of the ROI mask) minus
#     half a pixel. The distance map gives 1 for a pixel that lies directly
#     beside a pixel outside the ROI; the outline runs between the two, half
#     a pixel from each centre. So the outermost row of pixels has depth 0.5
#     pixel and covers depths 0 to 1 pixel. Along straight horizontal and
#     vertical outlines this is exact; along slanted and curved outlines the
#     distance is measured to the centre of the nearest outside pixel, which
#     is within about half a pixel of the true distance. ImageJ's distance
#     map itself is exact in almost every pixel; on irregular shapes a few
#     pixels in a thousand come out up to about 0.05 pixel too large
#     (measured against a brute-force search in tests/test_depth.py).
#   - The ROI mask is the one ImageJ uses to measure inside the ROI, so the
#     bin areas add up to ROI_Area exactly.
#   - A hole inside the ROI (a composite ROI) is an outline too: pixels next
#     to the hole are shallow. Separate parts of a composite ROI each have
#     their own outline.
#   - The image border. Where the ROI is cut by the edge of the image, the
#     real outline is outside the picture. By default the image border does
#     NOT count as outline; depth is then the distance to the nearest outline
#     that is visible, which can be larger than the true depth near the cut
#     side. With depth.image_border_is_surface = true the border counts as
#     outline, and depth near the cut side can be smaller than the true
#     depth. The true depth lies between the two. Such ROIs are flagged in
#     ROI_Touches_Image_Border. A ROI that covers the whole image has no
#     visible outline at all: no numbers are produced, only a note.
#   - Bins: bin n (counted from 1) holds depths from (n - 1) x bin width up
#     to, but not including, n x bin width (a depth exactly on an edge goes
#     into the deeper bin; see DEPTH_EDGE_FACTOR). The bin width is in the
#     length unit of the image (pixels when the image is not calibrated). With a
#     maximum depth set, the limit is rounded up to a whole number of bins
#     and everything at or beyond it goes into one last open-ended bin, so
#     nothing is dropped.
#   - Relative depth = depth divided by the largest depth found in the ROI.
#   - An object belongs to the bin of the pixel that contains its centroid.
#     When the centroid lies outside the object's own outline (a crescent),
#     the mean depth of the object's pixels is used instead, and the object
#     row says so.
#   - Pixels must be square (width and height within 0.1 percent); otherwise
#     no numbers are produced, only a note, because the distance map counts
#     in pixels.
#
# This file is loaded by load_engine_modules() into the namespace of
# cellquantengine.py: the engine's names are globals here. It has no user
# interface, opens no window, and changes no Fiji option.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij.measure import Measurements
from ij.plugin.filter import EDM
from ij.process import ImageProcessor, ImageStatistics, Blitter

import math


DEPTH_TABLE_NAME = "depth_profile.csv"
# More bins than this for one ROI is taken to be a mistake in the bin width.
DEPTH_MAX_BINS = 2000
# Pixel width and height may differ by this fraction at most.
DEPTH_SQUARE_TOLERANCE = 0.001
# The distance map holds single-precision numbers. Without help, a depth that
# lies exactly on a bin edge (for example depth 34.5 pixel with bins of 2.3
# pixel) can come out a hair below the edge and land in the shallower bin.
# Depths are therefore stretched by half a millionth before they are sorted
# into bins, which is far below anything the distance map can resolve, so a
# depth exactly on an edge always goes into the deeper bin.
DEPTH_EDGE_FACTOR = 1.0000005

# (property key, option key, kind, default)
DEPTH_SETTING_SPECS = [
    ("depth.enabled", "depth_enabled", "bool", False),
    ("depth.bin_width", "depth_bin_width", "text", "25"),
    ("depth.max_depth", "depth_max_depth", "text", "0"),
    ("depth.image_border_is_surface", "depth_image_border_is_surface", "bool", False)
]

DEPTH_PROFILE_HEADERS = [
    "Run_ID", "Application_Version", "Source_Image", "Frame",
    "Parent_ROI_Index", "Parent_ROI_Name", "ROI_Touches_Image_Border",
    "ROI_Border_Contact_Pixels", "Image_Border_Is_Surface", "ROI_Max_Depth",
    "Channel", "Channel_Name", "Analysis_Mode", "Bin_Index",
    "Depth_Start", "Depth_End", "Depth_Mid", "Depth_Unit",
    "Relative_Depth_Start", "Relative_Depth_End",
    "Bin_Pixels", "Bin_Area", "Mean_Intensity",
    "Positive_Pixels", "Positive_Area_Fraction_Percent",
    "Object_Count", "Object_Density_Per_Area", "Object_Density_Per_mm2", "Note"
]

DEPTH_OBJECT_HEADERS = [
    "Object_Depth_From_Outline", "Object_Relative_Depth",
    "Object_Depth_Bin", "Object_Depth_Measured_At"
]

DEPTH_AT_CENTROID = "centroid"
DEPTH_AT_OBJECT_MEAN = "object mean (centroid outside the object)"
DEPTH_AT_NOTHING = "not available"

DEPTH_NOTE_EMPTY = "No depth profile: the ROI contains no pixel."
DEPTH_NOTE_NO_SURFACE = (
    "No depth profile: the ROI covers the whole image, so no part of its outline is visible. "
    "Draw the tissue outline, or set depth.image_border_is_surface to true to measure from the image border."
)
DEPTH_NOTE_NOT_SQUARE = (
    "No depth profile: the pixels are not square (width and height differ by more than 0.1 percent)."
)
DEPTH_NOTE_OPEN_BIN = "Open-ended last bin: every depth at or beyond Depth_Start."


# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

def depth_settings(global_options):
    """Read and check the depth-profile settings. A missing key uses its
    default, so options saved by an older version keep working."""
    options = global_options if global_options is not None else {}
    raw = {}
    for spec in DEPTH_SETTING_SPECS:
        value = options.get(spec[1])
        raw[spec[1]] = spec[3] if value is None else value
    bin_width = parse_float_value(raw["depth_bin_width"], "Depth bin width", 0)
    if bin_width <= 0:
        raise ValueError("Depth bin width must be greater than zero.")
    return {
        "enabled": parse_bool(raw["depth_enabled"], False),
        "bin_width": bin_width,
        "max_depth": parse_float_value(raw["depth_max_depth"], "Maximum depth of the depth profile", 0),
        "border_is_surface": parse_bool(raw["depth_image_border_is_surface"], False)
    }


def depth_summary_headers(configs):
    """Columns appended to roi_summary.csv."""
    headers = ["ROI_Max_Depth", "ROI_Mean_Depth", "ROI_Touches_Image_Border", "Depth_Profile_Note"]
    for config in configs:
        prefix = channel_column_prefix(config)
        if needs_objects(config):
            headers.append(prefix + "_Median_Object_Depth")
        if needs_positive_area(config):
            headers.append(prefix + "_Positive_Area_Mean_Depth")
    return headers


# -----------------------------------------------------------------------------
# The distance map of one ROI
# -----------------------------------------------------------------------------

def depth_crop(ip, geometry):
    """Copy of the part of a full-size image that the ROI's distance map
    covers. The image itself is left as it was."""
    ip.setRoi(geometry["x"], geometry["y"], geometry["width"], geometry["height"])
    try:
        return ip.crop()
    finally:
        ip.resetRoi()


def depth_border_pixels(inside, width, height):
    """Number of ROI pixels on the outermost rows and columns of the image
    (a corner pixel is counted once). inside is the full-size ROI mask."""
    strips = [(0, 0, width, 1)]
    if height > 1:
        strips.append((0, height - 1, width, 1))
    if height > 2:
        strips.append((0, 1, 1, height - 2))
        if width > 1:
            strips.append((width - 1, 1, 1, height - 2))
    total = 0
    try:
        for strip in strips:
            inside.setRoi(strip[0], strip[1], strip[2], strip[3])
            total += inside.getHistogram()[255]
    finally:
        inside.resetRoi()
    return total


def depth_roi_geometry(inside, bounds, width, height):
    """Cut the ROI mask down to the ROI's bounding box plus one pixel of
    margin on every side where the image has one. The margin holds pixels
    outside the ROI, which is what the distance map measures to; a side
    without margin is a side where the ROI reaches the image border."""
    x0 = max(0, bounds.x - 1)
    y0 = max(0, bounds.y - 1)
    x1 = min(width, bounds.x + bounds.width + 1)
    y1 = min(height, bounds.y + bounds.height + 1)
    geometry = {"x": x0, "y": y0, "width": x1 - x0, "height": y1 - y0}
    if x1 - x0 < 1 or y1 - y0 < 1:
        # The ROI lies entirely outside the image: it holds no pixel, and
        # an empty crop cannot be made.
        geometry.update({"x": 0, "y": 0, "width": 0, "height": 0, "mask": None,
                         "roi_pixels": 0, "outside_pixels": 0, "border_pixels": 0})
        return geometry
    geometry["mask"] = depth_crop(inside, geometry)
    histogram = geometry["mask"].getHistogram()
    geometry["roi_pixels"] = histogram[255]
    geometry["outside_pixels"] = histogram[0]
    geometry["border_pixels"] = depth_border_pixels(inside, width, height)
    return geometry


def depth_distance_map(geometry, border_is_surface):
    """Add the distance map to a geometry: "edm" holds, for every ROI pixel,
    the distance in pixels from its centre to the centre of the nearest pixel
    outside the ROI (0 outside the ROI). Pixels beyond the image border count
    as outside the ROI only when border_is_surface is true. Also stores the
    largest and the mean value inside the ROI."""
    edm = EDM().makeFloatEDM(geometry["mask"], 0, bool(border_is_surface))
    if edm is None:
        raise RuntimeError("Depth profile: the distance map was interrupted.")
    edm.setMask(geometry["mask"])
    stats = ImageStatistics.getStatistics(edm, Measurements.MEAN | Measurements.MIN_MAX, None)
    edm.resetRoi()
    geometry["edm"] = edm
    geometry["edm_pixels"] = edm.getPixels()
    geometry["max_edm"] = stats.max
    geometry["mean_edm"] = stats.mean
    return geometry


def depth_bin_codes(edm, mask, bins_per_pixel, regular_bins):
    """Turn the distance map into a 16-bit image of bin numbers: 0 outside
    the ROI, n for a pixel of bin n (counted from 1). The distance map is
    used up. regular_bins > 0 puts everything beyond that many bins into bin
    regular_bins + 1. Returns (codes, pixels per bin as a list)."""
    edm.resetRoi()
    # bin number = floor((distance - 0.5) x bins per pixel) + 1. The
    # conversion to 16-bit below rounds, so 0.5 is added here instead of 1.
    edm.subtract(0.5)
    edm.multiply(bins_per_pixel * DEPTH_EDGE_FACTOR)
    edm.add(0.5)
    outside = mask.duplicate()
    outside.invert()
    edm.setValue(-1.0)
    edm.fill(outside)
    codes = edm.convertToShort(False)
    if regular_bins > 0 and regular_bins < 65535:
        codes.max(float(regular_bins + 1))
    histogram = codes.getHistogram()
    last = len(histogram) - 1
    while last > 0 and histogram[last] == 0:
        last -= 1
    return (codes, [int(histogram[code]) for code in range(1, last + 1)])


# -----------------------------------------------------------------------------
# Measurements per bin
# -----------------------------------------------------------------------------

def depth_shift_image(codes, base, bin_count):
    """Bin numbers base+1 .. base+255 moved into the high byte of a 16-bit
    image (bin number x 256); every other pixel is 0. Adding a value of
    0-255 to it gives one number per (bin, value) pair, so one histogram of
    the sum counts all pairs at once."""
    shift = codes.duplicate()
    if bin_count > 255:
        codes.setThreshold(base + 1, base + 255, ImageProcessor.NO_LUT_UPDATE)
        outside = codes.createMask()
        codes.resetThreshold()
        outside.invert()
        shift.subtract(float(base))
        shift.setValue(0.0)
        shift.fill(outside)
    shift.multiply(256.0)
    return shift


def depth_add_bin_sums(crop, shift, base, bin_count, sums):
    """Add the sum of the pixel values of each bin base+1 .. base+255 to
    sums (a list, one entry per bin). crop is an 8-bit or 16-bit copy and is
    used up. The value is split into its high and low byte; each byte is
    counted against the bin number in one histogram, which is exact."""
    if crop.getBitDepth() == 8:
        crop = crop.convertToShort(False)
    low = crop.duplicate()
    getattr(low, "and")(255)
    crop.copyBits(low, 0, 0, Blitter.SUBTRACT)
    crop.multiply(1.0 / 256.0)
    crop.copyBits(shift, 0, 0, Blitter.ADD)
    low.copyBits(shift, 0, 0, Blitter.ADD)
    high_histogram = crop.getHistogram()
    low_histogram = low.getHistogram()
    for code in range(1, min(255, bin_count - base) + 1):
        offset = code * 256
        high_sum = 0
        low_sum = 0
        for value in range(1, 256):
            high_sum += value * high_histogram[offset + value]
            low_sum += value * low_histogram[offset + value]
        sums[base + code - 1] += 256 * high_sum + low_sum


def depth_masked_bin_means(crop, calibration, codes, bin_count):
    """Mean value of each bin, one bin at a time, with ImageJ's own
    statistics. Used for images the histogram method cannot handle exactly:
    32-bit, RGB, and images with a density calibration."""
    means = []
    try:
        for index in range(bin_count):
            codes.setThreshold(index + 1, index + 1, ImageProcessor.NO_LUT_UPDATE)
            crop.setMask(codes.createMask())
            stats = ImageStatistics.getStatistics(crop, Measurements.MEAN, calibration)
            means.append(stats.mean if stats.pixelCount > 0 else float("nan"))
    finally:
        codes.resetThreshold()
        crop.resetRoi()
    return means


def depth_bin_means(context, geometry, codes, bin_pixels):
    """Mean prepared intensity of every bin for every analysed channel.
    Returns channel -> list of means (NaN for a bin without pixels)."""
    bin_count = len(bin_pixels)
    means = {}
    exact = []
    for config in context["configs"]:
        channel = config["channel"]
        imp = context["prepared"][channel]
        if imp.getBitDepth() in [8, 16] and not imp.getCalibration().calibrated():
            exact.append((channel, imp.getProcessor(), [0] * bin_count))
        else:
            means[channel] = depth_masked_bin_means(
                depth_crop(imp.getProcessor(), geometry), imp.getCalibration(), codes, bin_count
            )
    base = 0
    while base < bin_count and len(exact) > 0:
        shift = depth_shift_image(codes, base, bin_count)
        for entry in exact:
            depth_add_bin_sums(depth_crop(entry[1], geometry), shift, base, bin_count, entry[2])
        base += 255
    for entry in exact:
        means[entry[0]] = [safe_ratio(entry[2][index], bin_pixels[index]) for index in range(bin_count)]
    return means


def depth_positive_mask(raw_ip, geometry):
    """The marker-positive pixels of the ROI: value 255 in the raw mask (the
    same pixels roi_summary counts) and inside the ROI."""
    positive = depth_crop(raw_ip, geometry)
    positive.resetThreshold()
    positive.threshold(254)
    positive.copyBits(geometry["mask"], 0, 0, Blitter.AND)
    return positive


def depth_positive_counts(positive, codes, bin_count):
    """Number of positive pixels in each bin."""
    negative = positive.duplicate()
    negative.invert()
    keyed = codes.duplicate()
    keyed.setValue(0.0)
    keyed.fill(negative)
    histogram = keyed.getHistogram()
    return [int(histogram[code]) for code in range(1, bin_count + 1)]


# -----------------------------------------------------------------------------
# Objects
# -----------------------------------------------------------------------------

def depth_place_object(obj, calibration, image_height, geometry):
    """Where the depth of one object is read. Returns (position, distance,
    how): position is the index, in the distance map, of the pixel holding
    the object's centroid, with the distance-map value there. When the
    centroid is not inside the object's own outline or not inside the parent
    ROI, position is -1 and distance is the mean distance-map value of the
    object's pixels inside the parent ROI."""
    edm = geometry["edm"]
    object_roi = obj["roi"]
    row = obj.get("row") or {}
    try:
        # Centroids are stored in calibrated units; convert back to pixels.
        x = int(math.floor(calibration.getRawX(float(row.get("Centroid_X", obj.get("x")))) + 0.0000001))
        y = int(math.floor(calibration.getRawY(float(row.get("Centroid_Y", obj.get("y"))), image_height) + 0.0000001))
        local_x = x - geometry["x"]
        local_y = y - geometry["y"]
        if (0 <= local_x < geometry["width"] and 0 <= local_y < geometry["height"] and
                object_roi.contains(x, y)):
            position = local_y * geometry["width"] + local_x
            distance = geometry["edm_pixels"][position]
            if distance > 0:
                return (position, distance, DEPTH_AT_CENTROID)
    except ANY_ERROR:
        pass
    bounds = object_roi.getBounds()
    moved = clone_roi(object_roi)
    moved.setLocation(bounds.x - geometry["x"], bounds.y - geometry["y"])
    try:
        edm.setRoi(moved)
        edm.setThreshold(0.5, 1.0e30, ImageProcessor.NO_LUT_UPDATE)
        stats = ImageStatistics.getStatistics(edm, Measurements.MEAN | Measurements.LIMIT, None)
    finally:
        edm.resetThreshold()
        edm.resetRoi()
    if stats.pixelCount > 0:
        return (-1, stats.mean, DEPTH_AT_OBJECT_MEAN)
    return (-1, float("nan"), DEPTH_AT_NOTHING)


def depth_median(values):
    ordered = sorted(values)
    count = len(ordered)
    if count == 0:
        return float("nan")
    if count % 2 == 1:
        return ordered[count // 2]
    return (ordered[count // 2 - 1] + ordered[count // 2]) / 2.0


def depth_assign_objects(objects, placed, codes, plan):
    """Write the depth columns of every object row and count the objects of
    each bin. Returns (objects per bin, list of object depths)."""
    code_pixels = codes.getPixels()
    bin_count = plan["bin_count"]
    counts = [0] * bin_count
    depths = []
    for index in range(len(objects)):
        position, distance, how = placed[index]
        row = objects[index].get("row")
        if how == DEPTH_AT_NOTHING:
            if row is not None:
                row["Object_Depth_Measured_At"] = how
            continue
        if position >= 0:
            bin_index = (code_pixels[position] & 0xffff) - 1
        else:
            bin_index = int(math.floor((distance - 0.5) * plan["bins_per_pixel"] * DEPTH_EDGE_FACTOR))
        bin_index = max(0, min(bin_count - 1, bin_index))
        depth = (distance - 0.5) * plan["pixel"]
        counts[bin_index] += 1
        depths.append(depth)
        if row is not None:
            row["Object_Depth_From_Outline"] = depth
            row["Object_Relative_Depth"] = safe_ratio(depth, plan["max_depth"])
            row["Object_Depth_Bin"] = bin_index + 1
            row["Object_Depth_Measured_At"] = how
    return (counts, depths)


# -----------------------------------------------------------------------------
# One parent ROI of one frame
# -----------------------------------------------------------------------------

def depth_bin_plan(geometry, settings, pixel):
    """How depths turn into bins for this ROI. Returns a dictionary; its
    "note" is not empty when the bin width gives too many bins."""
    bin_width = settings["bin_width"]
    bins_per_pixel = pixel / bin_width
    natural = int(math.floor((geometry["max_edm"] - 0.5) * bins_per_pixel * DEPTH_EDGE_FACTOR)) + 1
    regular_bins = 0
    if settings["max_depth"] > 0:
        regular_bins = max(1, int(math.ceil(settings["max_depth"] / bin_width - 0.000000001)))
    expected = natural if regular_bins == 0 else min(natural, regular_bins + 1)
    note = ""
    if expected > DEPTH_MAX_BINS:
        note = (
            "No depth profile: a bin width of " + str(bin_width) + " would give " + str(expected) +
            " bins (the limit is " + str(DEPTH_MAX_BINS) + "). Use a larger depth.bin_width or set depth.max_depth."
        )
    return {
        "pixel": pixel, "bin_width": bin_width, "bins_per_pixel": bins_per_pixel,
        "regular_bins": regular_bins, "max_depth": (geometry["max_edm"] - 0.5) * pixel,
        "mean_depth": (geometry["mean_edm"] - 0.5) * pixel, "note": note
    }


def depth_base_row(context, roi_index, summary, geometry, settings):
    return {
        "Run_ID": context["run_id"],
        "Application_Version": APP_VERSION,
        "Source_Image": context["source_title"],
        "Frame": context["frame"],
        "Parent_ROI_Index": roi_index + 1,
        "Parent_ROI_Name": summary.get("Parent_ROI_Name", ""),
        "ROI_Touches_Image_Border": 1 if geometry["border_pixels"] > 0 else 0,
        "ROI_Border_Contact_Pixels": geometry["border_pixels"],
        "Image_Border_Is_Surface": 1 if settings["border_is_surface"] else 0,
        "Depth_Unit": summary.get("Spatial_Unit", "")
    }


def depth_channel_row(base, config):
    row = dict(base)
    row["Channel"] = config["channel"]
    row["Channel_Name"] = config["name"]
    row["Analysis_Mode"] = config["analysis_mode"]
    return row


def depth_bin_edges(plan, index):
    """(start, end, note) of bin index (counted from 0), in image units."""
    start = index * plan["bin_width"]
    if plan["regular_bins"] > 0 and index >= plan["regular_bins"]:
        return (start, max(start, plan["max_depth"]), DEPTH_NOTE_OPEN_BIN)
    return (start, (index + 1) * plan["bin_width"], "")


def depth_profile_rows(base, configs, plan, bin_pixels, measured, pixel_area):
    """The table rows of one ROI: one per channel and bin."""
    rows = []
    unit = base["Depth_Unit"]
    max_depth = plan["max_depth"]
    for config in configs:
        channel = config["channel"]
        for index in range(len(bin_pixels)):
            start, end, note = depth_bin_edges(plan, index)
            area = bin_pixels[index] * pixel_area
            row = depth_channel_row(base, config)
            row["ROI_Max_Depth"] = max_depth
            row["Bin_Index"] = index + 1
            row["Depth_Start"] = start
            row["Depth_End"] = end
            row["Depth_Mid"] = (start + end) / 2.0
            row["Relative_Depth_Start"] = safe_ratio(start, max_depth)
            row["Relative_Depth_End"] = safe_ratio(min(end, max_depth), max_depth)
            row["Bin_Pixels"] = bin_pixels[index]
            row["Bin_Area"] = area
            row["Mean_Intensity"] = measured["means"][channel][index]
            row["Note"] = note
            if channel in measured["positive"]:
                positive = measured["positive"][channel][index]
                row["Positive_Pixels"] = positive
                row["Positive_Area_Fraction_Percent"] = safe_ratio(positive, bin_pixels[index]) * 100.0
            if channel in measured["objects"]:
                count = measured["objects"][channel][index]
                row["Object_Count"] = count
                row["Object_Density_Per_Area"] = safe_ratio(count, area)
                row["Object_Density_Per_mm2"] = (
                    safe_ratio(count * 1000000.0, area) if is_micrometer_unit(unit) else ""
                )
            rows.append(row)
    return rows


def depth_measure_roi(context, settings, roi_index):
    """Depth profile of one parent ROI in the current frame. Fills the
    appended columns of its roi_summary row and of its object rows, and
    returns its rows of the depth table."""
    configs = context["configs"]
    summary = context["summary_rows"][roi_index]
    inside = context["get_roi_masks"](roi_index)[0]
    geometry = depth_roi_geometry(
        inside, context["parent_rois"][roi_index].getBounds(), context["width"], context["height"]
    )
    base = depth_base_row(context, roi_index, summary, geometry, settings)
    summary["ROI_Touches_Image_Border"] = base["ROI_Touches_Image_Border"]
    pixel_width = float(summary.get("Pixel_Width", 1.0))
    pixel_height = float(summary.get("Pixel_Height", 1.0))

    note = ""
    plan = None
    if geometry["roi_pixels"] == 0:
        note = DEPTH_NOTE_EMPTY
    elif abs(pixel_width - pixel_height) > DEPTH_SQUARE_TOLERANCE * max(pixel_width, pixel_height):
        note = DEPTH_NOTE_NOT_SQUARE
    elif geometry["outside_pixels"] == 0 and not settings["border_is_surface"]:
        note = DEPTH_NOTE_NO_SURFACE
    else:
        depth_distance_map(geometry, settings["border_is_surface"])
        pixel = (pixel_width + pixel_height) / 2.0
        plan = depth_bin_plan(geometry, settings, pixel if pixel > 0 else 1.0)
        note = plan["note"]
    if note != "":
        summary["Depth_Profile_Note"] = note
        context["engine"].log("Depth profile, " + to_text(base["Parent_ROI_Name"]) + ": " + note)
        rows = []
        for config in configs:
            row = depth_channel_row(base, config)
            row["Note"] = note
            rows.append(row)
        return rows

    summary["ROI_Max_Depth"] = plan["max_depth"]
    summary["ROI_Mean_Depth"] = plan["mean_depth"]
    edm = geometry["edm"]

    # Everything that reads the distance map itself comes first, because the
    # map is turned into bin numbers in place afterwards.
    placed = {}
    positive_masks = {}
    for config in configs:
        channel = config["channel"]
        prefix = channel_column_prefix(config)
        if needs_objects(config):
            calibration = context["prepared"][channel].getCalibration()
            placed[channel] = [
                depth_place_object(obj, calibration, context["height"], geometry)
                for obj in context["objects"].get((channel, roi_index), [])
            ]
        if needs_positive_area(config):
            positive = depth_positive_mask(context["masks"][channel]["raw"].getProcessor(), geometry)
            edm.setMask(positive)
            stats = ImageStatistics.getStatistics(edm, Measurements.MEAN, None)
            edm.resetRoi()
            positive_masks[channel] = positive
            summary[prefix + "_Positive_Area_Mean_Depth"] = (
                (stats.mean - 0.5) * plan["pixel"] if stats.pixelCount > 0 else float("nan")
            )

    codes, bin_pixels = depth_bin_codes(edm, geometry["mask"], plan["bins_per_pixel"], plan["regular_bins"])
    geometry["edm"] = None
    geometry["edm_pixels"] = None
    plan["bin_count"] = len(bin_pixels)

    measured = {"means": depth_bin_means(context, geometry, codes, bin_pixels), "positive": {}, "objects": {}}
    for config in configs:
        channel = config["channel"]
        if channel in positive_masks:
            measured["positive"][channel] = depth_positive_counts(positive_masks[channel], codes, plan["bin_count"])
        if channel in placed:
            counts, depths = depth_assign_objects(
                context["objects"].get((channel, roi_index), []), placed[channel], codes, plan
            )
            measured["objects"][channel] = counts
            summary[channel_column_prefix(config) + "_Median_Object_Depth"] = depth_median(depths)
    return depth_profile_rows(base, configs, plan, bin_pixels, measured, pixel_width * pixel_height)


# -----------------------------------------------------------------------------
# The run extension
# -----------------------------------------------------------------------------

class DepthProfileExtension(RunExtension):
    """Writes depth_profile.csv and appends the depth columns to the summary
    and object tables when depth.enabled is true."""

    name = "depth_profile"

    def is_enabled(self, global_options, configs):
        return parse_bool((global_options or {}).get("depth_enabled"), False)

    def summary_headers(self, configs, global_options):
        return depth_summary_headers(configs)

    def object_headers(self, configs, global_options):
        return list(DEPTH_OBJECT_HEADERS)

    def begin_run(self, context):
        state = context["state"][self.name]
        # Checked here, before any image work, so a mistyped bin width stops
        # the run at once with a clear message.
        state["settings"] = depth_settings(context["global_options"])
        state["rows"] = []

    def frame_done(self, context):
        state = context["state"][self.name]
        for roi_index in range(len(context["parent_rois"])):
            context["engine"].check_cancelled()
            state["rows"] += depth_measure_roi(context, state["settings"], roi_index)

    def end_run(self, context):
        return [(DEPTH_TABLE_NAME, list(DEPTH_PROFILE_HEADERS), context["state"][self.name]["rows"])]


register_global_settings(DEPTH_SETTING_SPECS)
register_settings_form(
    "Depth profile",
    "Measures every parent ROI in bins of distance from its outline. The parent ROI must be the outline "
    "of the tissue. Lengths are in the unit of the image (pixels when the image is not calibrated).",
    [
        ("depth_enabled", "Measure a depth profile from the ROI outline", "bool", None),
        ("depth_bin_width", "Bin width", "text", None),
        ("depth_max_depth", "Maximum depth (0 = to the centre; deeper goes into one last bin)", "text", None),
        ("depth_image_border_is_surface", "Count the image border as outline where it cuts the ROI", "bool", None)
    ]
)
register_run_extension(DepthProfileExtension())
register_engine_module("cellquantdepth.py", "2.3.0")

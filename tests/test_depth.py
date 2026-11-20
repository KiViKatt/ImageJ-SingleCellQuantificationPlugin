# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach
# Copyright (C) 2026 Maksym Taran

# Cell Quant Studio Universal - tests of the depth profile (cellquantdepth.py)
#
# Loaded by run_tests.py, which provides the helpers used here.
# The source is kept plain ASCII with no encoding declaration on purpose.
#
# Every expected number below follows from how the test image is drawn.
# The image is 240 x 200 pixels at 0.5 micrometer per pixel, four channels:
#   C1 "DAPI"  bright disks of radius 6 (13 pixels across) at chosen places.
#              A disk drawn around pixel (cx, cy) is symmetric about that
#              pixel, so its centroid lies in pixel (cx, cy).
#   C2 "MAP2"  bright (900) everywhere from 10 pixels outside the main test
#              ROI to 20 pixels inside it, dark (100) further in.
#   C3 "Flat"  500 everywhere.
#   C4 "Ramp"  200 x column + 7: values need both bytes of a 16-bit number.
# The main test ROI is the rectangle x 20..219, y 20..179 (200 x 160 pixels).
# Depth of a pixel = distance of its centre from the outline, so the pixels
# of "ring" k (k = 0 on the outline) have depth k + 0.5 pixel, and ring k of
# a W x H rectangle holds (W - 2k)(H - 2k) - (W - 2k - 2)(H - 2k - 2) pixels.

DEPTH_RECT = (20, 20, 200, 160)
DEPTH_DISKS = [(27, 100), (45, 100), (120, 100), (120, 32), (205, 150), (60, 60), (8, 8)]
DEPTH_SAVE_KEYS = ["save_prepared", "save_masks", "save_overlays", "save_object_rois", "save_parent_rois"]


def depth_ramp_plane():
    ip = ShortProcessor(WIDTH, HEIGHT)
    for y in range(HEIGHT):
        for x in range(WIDTH):
            ip.set(x, y, 200 * x + 7)
    return ip


def depth_band_plane():
    ip = blank_plane(100)
    draw_rect(ip, (10, 10, 220, 180), 900)
    draw_rect(ip, (40, 40, 160, 120), 100)
    return ip


def depth_image(title, disks=None, frames=1, slices=1, calibrated=True, first_plane=None):
    """The four-channel test image described at the top of this file.
    disks is one list, or one list per frame."""
    stack = ImageStack(WIDTH, HEIGHT)
    for frame in range(frames):
        frame_disks = DEPTH_DISKS if disks is None else disks
        if len(frame_disks) > 0 and isinstance(frame_disks[0], list):
            frame_disks = frame_disks[frame]
        for z in range(slices):
            stack.addSlice("", first_plane.duplicate() if first_plane is not None else nuclei_plane(frame_disks))
            stack.addSlice("", depth_band_plane())
            stack.addSlice("", blank_plane(500))
            stack.addSlice("", depth_ramp_plane())
    imp = ImagePlus(title, stack)
    imp.setDimensions(4, slices, frames)
    imp = CompositeImage(imp, CompositeImage.GRAYSCALE)
    if calibrated:
        cal = imp.getCalibration()
        cal.pixelWidth = PIXEL
        cal.pixelHeight = PIXEL
        cal.setUnit("um")
    return imp


def depth_configs(E, map2_overrides=None, ramp_projection="Max Z"):
    return [
        nuclei_config(E),
        marker_config(E, map2_overrides),
        E.normalize_channel_config(3, {"enabled": True, "name": "Flat", "analysis_mode": "Intensity only"}),
        E.normalize_channel_config(4, {
            "enabled": True, "name": "Ramp", "analysis_mode": "Intensity only", "projection": ramp_projection
        })
    ]


def depth_run(E, imp, configs, rois, tmp, name, overrides=None, frames=None):
    """One analysis with the depth profile switched on (unless overridden).
    The run folder is always called "run", so the Run_ID is the same in
    every test run and whole files can be compared."""
    options = E.default_global_options()
    options["depth_enabled"] = True
    options["depth_bin_width"] = "5"
    if overrides:
        options.update(overrides)
    engine_options = E.engine_options_from_global(options)
    for key in DEPTH_SAVE_KEYS:
        engine_options[key] = False
    run_dir = os.path.join(tmp, name, "run")
    os.makedirs(run_dir)
    return E.AnalysisEngine().analyze_image(
        imp, "", configs, [E.clone_roi(roi) for roi in rois], engine_options, run_dir,
        frames if frames is not None else [1], settings=(options, configs)
    )


def depth_rows(result, channel, roi_index=1, frame=1):
    rows = []
    for row in result["extra_tables"]["depth_profile.csv"]["rows"]:
        if row["Channel"] == channel and row["Parent_ROI_Index"] == roi_index and row["Frame"] == frame:
            rows.append(row)
    return rows


def depth_column(rows, key):
    return [row.get(key, "") for row in rows]


def depth_rect_ring(width, height, k):
    """Pixels of ring k of a width x height rectangle (see the top of the file)."""
    outer = max(0, width - 2 * k) * max(0, height - 2 * k)
    inner = max(0, width - 2 * k - 2) * max(0, height - 2 * k - 2)
    return outer - inner


def depth_rect_bins(width, height, rings_per_bin):
    counts = []
    k = 0
    while depth_rect_ring(width, height, k) > 0:
        if k % rings_per_bin == 0:
            counts.append(0)
        counts[-1] += depth_rect_ring(width, height, k)
        k += 1
    return counts


def depth_object_row(result, cx, cy, channel=1):
    """The object row of the disk drawn around pixel (cx, cy)."""
    for row in result["object_rows"]:
        if (row["Segmentation_Channel"] == channel and abs(row["Centroid_X"] - (cx + 0.5) * PIXEL) < 0.01 and
                abs(row["Centroid_Y"] - (cy + 0.5) * PIXEL) < 0.01):
            return row
    return None


def depth_object_row_pixels(result, cx, cy, channel=1):
    """The same for an image without calibration (centroids in pixels)."""
    for row in result["object_rows"]:
        if (row["Segmentation_Channel"] == channel and abs(row["Centroid_X"] - (cx + 0.5)) < 0.01 and
                abs(row["Centroid_Y"] - (cy + 0.5)) < 0.01):
            return row
    return None


def depth_brute_force(inside, box, border_is_surface):
    """Depth in pixels of every ROI pixel inside box = (x0, y0, x1, y1), by
    trying every pixel that is not in the ROI: distance between the two
    pixel centres, minus half a pixel. Positions one step beyond the image
    count as "not in the ROI" only when border_is_surface is true. Returns
    {(x, y): depth}. box must hold the ROI with one pixel to spare on every
    side (a pixel further away is never the nearest one)."""
    members = []
    others = []
    for y in range(box[1], box[3] + 1):
        for x in range(box[0], box[2] + 1):
            if x < 0 or y < 0 or x >= WIDTH or y >= HEIGHT:
                if border_is_surface:
                    others.append((x, y))
            elif inside.get(x, y) != 0:
                members.append((x, y))
            else:
                others.append((x, y))
    depths = {}
    for member in members:
        best = None
        for other in others:
            dx = member[0] - other[0]
            dy = member[1] - other[1]
            squared = dx * dx + dy * dy
            if best is None or squared < best:
                best = squared
        depths[member] = math.sqrt(best) - 0.5
    return depths


def test_depth_rectangle(E, t, tmp):
    """A rectangle away from the image border: every bin is an exact
    rectangular frame, so every number is known."""
    imp = depth_image("depth_rect.tif")
    configs = depth_configs(E)
    result = depth_run(E, imp, configs, [Roi(*DEPTH_RECT)], tmp, "rect")
    imp.close()
    summary = result["summary_rows"][0]

    # 10 rings of one pixel per 5 um bin; 80 rings in all -> 8 bins.
    expected = depth_rect_bins(200, 160, 10)
    t.equal(len(expected), 8, "the construction gives 8 bins")
    flat = depth_rows(result, 3)
    t.equal(depth_column(flat, "Bin_Index"), list(range(1, 9)), "bins are numbered from 1 to 8")
    t.equal(depth_column(flat, "Bin_Pixels"), expected, "every bin holds exactly its rectangular frame")
    t.equal(expected[0], 200 * 160 - 180 * 140, "bin 1 is the outer 10-pixel frame")
    t.equal(expected[7], 60 * 20, "bin 8 is the innermost 60 x 20 block")
    t.equal(depth_column(flat, "Depth_Start"), [5.0 * j for j in range(8)], "bin starts are multiples of the bin width")
    t.equal(depth_column(flat, "Depth_End"), [5.0 * (j + 1) for j in range(8)], "bin ends are multiples of the bin width")
    t.equal(depth_column(flat, "Depth_Mid"), [5.0 * j + 2.5 for j in range(8)], "bin middles")
    # ImageJ stores the unit "um" as the micro sign followed by "m".
    t.equal(set(depth_column(flat, "Depth_Unit")), set([summary["Spatial_Unit"]]), "depth is in the unit of the image")
    t.check(E.is_micrometer_unit(flat[0]["Depth_Unit"]), "and that unit is the micrometer")
    for index in range(8):
        t.close(flat[index]["Bin_Area"], expected[index] * PIXEL * PIXEL, "bin " + str(index + 1) + " area is calibrated")
        t.close(flat[index]["Relative_Depth_Start"], 5.0 * index / 39.75, "bin " + str(index + 1) + " relative start")
    t.close(flat[7]["Relative_Depth_End"], 1.0, "the last bin ends at relative depth 1")
    t.close(sum(depth_column(flat, "Bin_Area")), summary["ROI_Area"], "bin areas add up to ROI_Area")
    t.equal(sum(depth_column(flat, "Bin_Pixels")), 200 * 160, "bin pixels add up to the ROI")
    # The deepest pixels are ring 79: depth 79.5 pixel = 39.75 um.
    t.close(summary["ROI_Max_Depth"], 39.75, "ROI_Max_Depth is the depth of the innermost ring")
    t.equal(set(depth_column(flat, "ROI_Max_Depth")), set([summary["ROI_Max_Depth"]]), "the table repeats ROI_Max_Depth")
    rings = [depth_rect_ring(200, 160, k) for k in range(80)]
    mean_depth = sum([(k + 0.5) * rings[k] for k in range(80)]) * PIXEL / float(sum(rings))
    t.close(summary["ROI_Mean_Depth"], mean_depth, "ROI_Mean_Depth is the mean over all rings")
    t.equal(summary["ROI_Touches_Image_Border"], 0, "the ROI does not touch the image border")
    t.equal(set(depth_column(flat, "ROI_Touches_Image_Border")), set([0]), "the table says so too")
    t.equal(summary.get("Depth_Profile_Note", ""), "", "no note for a normal ROI")

    # Intensity: uniform image -> same mean everywhere; the ramp is symmetric
    # about the middle column 119.5 of every frame -> 200 x 119.5 + 7.
    for row in flat:
        t.close(row["Mean_Intensity"], 500.0, "uniform channel, bin " + str(row["Bin_Index"]))
        t.equal(row.get("Positive_Pixels", ""), "", "an intensity-only channel has no positive area")
        t.equal(row.get("Object_Count", ""), "", "an intensity-only channel has no object count")
    ramp = depth_rows(result, 4)
    for row in ramp:
        t.close(row["Mean_Intensity"], 23907.0, "ramp channel, bin " + str(row["Bin_Index"]))
    weighted = sum([row["Mean_Intensity"] * row["Bin_Pixels"] for row in ramp]) / float(200 * 160)
    t.close(weighted, summary["C4_Ramp_Mean"], "bin means weighted by bin pixels give the ROI mean")

    # Marker: positive exactly in rings 0-19 = bins 1 and 2.
    band = depth_rows(result, 2)
    t.equal(depth_column(band, "Positive_Pixels"), expected[:2] + [0] * 6, "positive pixels are in the first two bins only")
    t.equal(depth_column(band, "Positive_Area_Fraction_Percent"), [100.0, 100.0] + [0.0] * 6, "positive fraction per bin")
    t.equal(depth_column(band, "Mean_Intensity"), [900.0, 900.0] + [100.0] * 6, "marker mean per bin")
    t.close(sum(depth_column(band, "Positive_Pixels")) * PIXEL * PIXEL, summary["C2_MAP2_Positive_Area"],
            "positive pixels add up to the Positive_Area of roi_summary")
    band_depth = sum([(k + 0.5) * rings[k] for k in range(20)]) * PIXEL / float(sum(rings[:20]))
    t.close(summary["C2_MAP2_Positive_Area_Mean_Depth"], band_depth, "mean depth of the positive pixels")

    # Objects: ring of the centroid pixel = distance to the nearest side.
    #   (27,100) ring 7 -> bin 1    (120,32) ring 12 -> bin 2    (205,150) ring 14 -> bin 2
    #   (45,100) ring 25 -> bin 3   (60,60) ring 40 -> bin 5     (120,100) ring 79 -> bin 8
    #   (8,8) lies outside the ROI.
    nuclei = depth_rows(result, 1)
    t.equal(depth_column(nuclei, "Object_Count"), [1, 2, 1, 0, 1, 0, 0, 1], "objects per bin")
    t.equal(sum(depth_column(nuclei, "Object_Count")), summary["C1_DAPI_Object_Count"], "object counts add up to the ROI count")
    t.close(nuclei[1]["Object_Density_Per_Area"], 2.0 / (expected[1] * PIXEL * PIXEL), "object density of bin 2")
    t.close(nuclei[1]["Object_Density_Per_mm2"], 2.0e6 / (expected[1] * PIXEL * PIXEL), "object density per mm2 of bin 2")
    t.equal(sum(depth_column(nuclei, "Positive_Pixels")), 6 * disk_pixel_count(), "positive pixels of an objects + area channel")
    known = [((27, 100), 7, 1), ((120, 32), 12, 2), ((205, 150), 14, 2), ((45, 100), 25, 3), ((60, 60), 40, 5), ((120, 100), 79, 8)]
    for center, ring, bin_index in known:
        row = depth_object_row(result, center[0], center[1])
        t.check(row is not None, "the disk at " + str(center) + " has an object row")
        if row is None:
            continue
        t.close(row["Object_Depth_From_Outline"], (ring + 0.5) * PIXEL, "depth of the disk at " + str(center))
        t.close(row["Object_Relative_Depth"], (ring + 0.5) * PIXEL / 39.75, "relative depth of the disk at " + str(center))
        t.equal(row["Object_Depth_Bin"], bin_index, "bin of the disk at " + str(center))
        t.equal(row["Object_Depth_Measured_At"], E.DEPTH_AT_CENTROID, "the disk at " + str(center) + " is read at its centroid")
    # depths 3.75, 6.25, 7.25, 12.75, 20.25, 39.75 -> median (7.25 + 12.75) / 2
    t.close(summary["C1_DAPI_Median_Object_Depth"], 10.0, "median object depth")

    # Files.
    t.equal(result["extra_tables"]["depth_profile.csv"]["headers"], E.DEPTH_PROFILE_HEADERS, "the depth table has its fixed columns")
    text = read_utf8(os.path.join(result["run_dir"], "depth_profile.csv"))
    lines = [line for line in text.split("\n") if line.strip() != ""]
    t.equal(len(lines), 1 + 4 * 8, "depth_profile.csv has one row per channel and bin")
    t.check(lines[0].startswith('"Run_ID","Application_Version","Source_Image"'), "the depth table starts with the run columns")
    t.check(lines[1].startswith('"run","' + E.APP_VERSION + '"'), "every row carries Run_ID and Application_Version")
    summary_header = read_utf8(os.path.join(result["run_dir"], "roi_summary.csv")).split("\n")[0]
    t.check(summary_header.endswith('"C1_DAPI_Median_Object_Depth","C1_DAPI_Positive_Area_Mean_Depth","C2_MAP2_Positive_Area_Mean_Depth"'),
            "the summary depth columns come last")
    object_header = read_utf8(os.path.join(result["run_dir"], "objects.csv")).split("\n")[0]
    t.check(object_header.endswith(",".join(['"' + name + '"' for name in E.DEPTH_OBJECT_HEADERS])),
            "the object depth columns come last")


def test_depth_switched_off(E, t, tmp):
    """Switched off (the default), every standard table is byte for byte
    what it is without the module; switched on, no existing value changes."""
    defaults = E.default_global_options()
    t.equal(defaults["depth_enabled"], False, "the depth profile is off by default")
    t.equal(defaults["depth_bin_width"], "25", "default bin width")
    t.equal(defaults["depth_max_depth"], "0", "default maximum depth")
    t.equal(defaults["depth_image_border_is_surface"], False, "the image border is not outline by default")
    t.equal(E.ENGINE_MODULE_VERSIONS.get("cellquantdepth.py"), E.APP_VERSION, "the module reports the application version")

    configs = depth_configs(E)
    rois = [Roi(*DEPTH_RECT), triangle_roi()]
    imp = depth_image("depth_off.tif")
    engine = E.AnalysisEngine()
    base_summary = engine.build_summary_headers(configs)
    base_objects = engine.build_object_headers(configs)
    off = depth_run(E, imp, configs, rois, tmp, "off", {"depth_enabled": False})
    on = depth_run(E, imp, configs, rois, tmp, "on")
    saved = list(E.RUN_EXTENSIONS)
    try:
        E.RUN_EXTENSIONS[:] = [extension for extension in saved if extension.name != "depth_profile"]
        absent = depth_run(E, imp, configs, rois, tmp, "absent", {"depth_enabled": False})
    finally:
        E.RUN_EXTENSIONS[:] = saved
    imp.close()

    t.equal(off["summary_headers"], base_summary, "switched off: summary columns are unchanged")
    t.equal(off["object_headers"], base_objects, "switched off: object columns are unchanged")
    t.check("depth_profile.csv" not in off["extra_table_names"], "switched off: no depth table")
    t.check(not os.path.exists(os.path.join(off["run_dir"], "depth_profile.csv")), "switched off: no depth file")
    for name in ["roi_summary.csv", "objects.csv", "thresholds.csv"]:
        t.check(read_utf8(os.path.join(off["run_dir"], name)) == read_utf8(os.path.join(absent["run_dir"], name)),
                "switched off: " + name + " is identical to a run without the module")

    t.equal(on["summary_headers"], base_summary + E.depth_summary_headers(configs), "switched on: summary columns are appended")
    t.equal(on["object_headers"], base_objects + E.DEPTH_OBJECT_HEADERS, "switched on: object columns are appended")
    t.equal(len(on["summary_rows"]), len(off["summary_rows"]), "same number of summary rows")
    t.equal(len(on["object_rows"]), len(off["object_rows"]), "same number of object rows")
    changed = []
    for index in range(len(off["summary_rows"])):
        for key in base_summary:
            if off["summary_rows"][index].get(key) != on["summary_rows"][index].get(key):
                changed.append(key)
    for index in range(len(off["object_rows"])):
        for key in base_objects:
            if off["object_rows"][index].get(key) != on["object_rows"][index].get(key):
                changed.append(key)
    t.equal(changed, [], "switched on: no existing summary or object value changes")
    t.check(read_utf8(os.path.join(off["run_dir"], "thresholds.csv")) == read_utf8(os.path.join(on["run_dir"], "thresholds.csv")),
            "switched on: thresholds.csv is unchanged")

    # A profile saved before this module existed has no depth keys.
    path = os.path.join(tmp, "depth_profile_settings.properties")
    options = E.default_global_options()
    options["depth_enabled"] = True
    options["depth_bin_width"] = "12.5"
    E.write_settings_properties(path, options, configs)
    loaded = E.load_profile_settings(path)[0]
    t.equal(loaded["depth_enabled"], True, "the switch survives a save and load")
    t.equal(E.depth_settings(loaded)["bin_width"], 12.5, "the bin width survives a save and load")
    legacy = E.read_properties_file(path)
    for key in ["depth.enabled", "depth.bin_width", "depth.max_depth", "depth.image_border_is_surface"]:
        legacy.remove(key)
    old = E.depth_settings(E.global_options_from_properties(legacy))
    t.equal(old, {"enabled": False, "bin_width": 25.0, "max_depth": 0.0, "border_is_surface": False},
            "a profile without depth keys loads with the defaults")

    def rejected(values):
        try:
            E.depth_settings(values)
        except ValueError:
            return True
        return False

    t.check(rejected({"depth_bin_width": "0"}), "a bin width of zero is rejected")
    t.check(rejected({"depth_bin_width": "wide"}), "a bin width that is not a number is rejected")
    t.check(rejected({"depth_max_depth": "-5"}), "a negative maximum depth is rejected")


def test_depth_disk_roi(E, t, tmp):
    """A circular ROI: the outline is curved, so depths are accurate to
    about half a pixel and bin areas follow ring areas within that."""
    imp = depth_image("depth_disk.tif")
    configs = depth_configs(E)
    oval = OvalRoi(20, 10, 181, 181)
    result = depth_run(E, imp, configs, [oval], tmp, "disk")
    imp.close()
    summary = result["summary_rows"][0]
    flat = depth_rows(result, 3)
    pixels = depth_column(flat, "Bin_Pixels")
    t.close(sum(depth_column(flat, "Bin_Area")), summary["ROI_Area"], "bin areas add up to ROI_Area for a circle")
    t.equal(sum(pixels), oval.getMask().getHistogram()[255], "bin pixels add up to the pixels of the circle")
    # Radius 90.5 pixel: depths reach about 90 pixel = 45 um -> 9 or 10 bins.
    t.check(len(flat) in [9, 10], "a circle of radius 90.5 pixel gives 9 or 10 bins of 10 pixel (got " + str(len(flat)) + ")")
    t.check(abs(summary["ROI_Max_Depth"] - 45.0) <= 0.5, "ROI_Max_Depth is the radius within one pixel (got " + repr(summary["ROI_Max_Depth"]) + ")")
    # Bin j (from 0) is the ring between radius 90.5 - 10 j and 80.5 - 10 j.
    # Each of its two edges may lie up to half a pixel off (the outer edge of
    # bin 1 is the ROI outline itself, which is exact).
    for j in range(8):
        outer = 90.5 - 10 * j
        inner = outer - 10
        area = math.pi * (outer * outer - inner * inner)
        slack = 0.5 * 2 * math.pi * inner + (0.5 * 2 * math.pi * outer if j > 0 else 0.01 * area)
        t.check(abs(pixels[j] - area) <= slack,
                "ring " + str(j + 1) + " of the circle holds " + str(pixels[j]) + " pixels, expected " + str(int(area)) + " +/- " + str(int(slack)))
    for row in flat:
        t.close(row["Mean_Intensity"], 500.0, "uniform channel in the circle, bin " + str(row["Bin_Index"]))

    # Pixel by pixel against geometry: true depth = 90.5 - distance of the
    # pixel centre from the centre of the circle (110.5, 100.5).
    inside = E.roi_inside_mask(oval, WIDTH, HEIGHT)
    geometry = E.depth_roi_geometry(inside, oval.getBounds(), WIDTH, HEIGHT)
    E.depth_distance_map(geometry, False)
    t.equal((geometry["x"], geometry["y"], geometry["width"], geometry["height"]), (19, 9, 183, 183),
            "the distance map covers the bounding box plus one pixel of margin")
    worst = 0.0
    total = 0.0
    count = 0
    for y in range(geometry["height"]):
        for x in range(geometry["width"]):
            if geometry["mask"].get(x, y) == 0:
                continue
            true_depth = 90.5 - math.sqrt((x + 19 + 0.5 - 110.5) ** 2 + (y + 9 + 0.5 - 100.5) ** 2)
            error = geometry["edm"].getf(x, y) - 0.5 - true_depth
            worst = max(worst, abs(error))
            total += error
            count += 1
    t.check(worst <= 0.75, "depth of every pixel of the circle is within 0.75 pixel of geometry (worst " + repr(worst) + ")")
    t.check(abs(total / count) <= 0.5, "mean depth error over the circle is within half a pixel (got " + repr(total / count) + ")")


def test_depth_image_border(E, t, tmp):
    """A ROI cut by the image border, measured both ways, and a ROI that is
    the whole image."""
    disks = [(10, 100), (60, 60)]
    configs = depth_configs(E)
    cut = Roi(0, 20, 100, 160)

    # Border not outline: depth = distance to the right, top or bottom side.
    # Pixels deeper than 10 j pixel form a (100 - 10 j) x (160 - 20 j) block.
    imp = depth_image("depth_cut.tif", disks)
    result = depth_run(E, imp, configs, [cut], tmp, "cut")
    blocks = [max(0, 100 - 10 * j) * max(0, 160 - 20 * j) for j in range(9)]
    expected = [blocks[j] - blocks[j + 1] for j in range(8)]
    flat = depth_rows(result, 3)
    t.equal(depth_column(flat, "Bin_Pixels"), expected, "cut ROI, border not outline: bin pixels")
    t.equal(sum(expected), 100 * 160, "the construction covers the cut ROI")
    t.equal(set(depth_column(flat, "ROI_Touches_Image_Border")), set([1]), "the cut ROI is flagged")
    t.equal(set(depth_column(flat, "ROI_Border_Contact_Pixels")), set([160]), "160 of its pixels lie on the image border")
    t.equal(set(depth_column(flat, "Image_Border_Is_Surface")), set([0]), "the table records the border setting")
    t.equal(result["summary_rows"][0]["ROI_Touches_Image_Border"], 1, "the summary flags the cut ROI")
    t.close(result["summary_rows"][0]["ROI_Max_Depth"], 39.75, "cut ROI, border not outline: maximum depth")
    # (10,100): min(99 - 10, 80, 79) = ring 79 -> bin 8; (60,60): min(39, 40, 119) = ring 39 -> bin 4
    t.equal(depth_column(depth_rows(result, 1), "Object_Count"), [0, 0, 0, 1, 0, 0, 0, 1], "cut ROI, border not outline: objects")

    # Border counts as outline: an ordinary 100 x 160 rectangle, 50 rings.
    result = depth_run(E, imp, configs, [cut], tmp, "cut_surface", {"depth_image_border_is_surface": True})
    flat = depth_rows(result, 3)
    t.equal(depth_column(flat, "Bin_Pixels"), depth_rect_bins(100, 160, 10), "cut ROI, border is outline: bin pixels")
    t.equal(len(flat), 5, "cut ROI, border is outline: 5 bins")
    t.equal(set(depth_column(flat, "Image_Border_Is_Surface")), set([1]), "the table records the border setting")
    t.equal(set(depth_column(flat, "ROI_Touches_Image_Border")), set([1]), "the cut ROI is still flagged")
    t.close(result["summary_rows"][0]["ROI_Max_Depth"], 24.75, "cut ROI, border is outline: maximum depth")
    # (10,100): ring 10 -> bin 2; (60,60): min(60, 39, 40, 119) = ring 39 -> bin 4
    t.equal(depth_column(depth_rows(result, 1), "Object_Count"), [0, 1, 0, 1, 0], "cut ROI, border is outline: objects")
    imp.close()

    # The whole image has no visible outline: a note, no numbers, and the
    # standard results are still there.
    imp = depth_image("depth_whole.tif")
    whole = Roi(0, 0, WIDTH, HEIGHT)
    result = depth_run(E, imp, configs, [whole], tmp, "whole")
    table = result["extra_tables"]["depth_profile.csv"]["rows"]
    t.equal(len(table), 4, "whole image: one note row per channel")
    t.equal(set(depth_column(table, "Note")), set([E.DEPTH_NOTE_NO_SURFACE]), "whole image: the note explains why")
    t.equal(set(depth_column(table, "Bin_Index")), set([""]), "whole image: no bins")
    t.equal(set(depth_column(table, "Bin_Pixels")), set([""]), "whole image: no numbers")
    summary = result["summary_rows"][0]
    t.equal(summary.get("Depth_Profile_Note"), E.DEPTH_NOTE_NO_SURFACE, "whole image: the summary carries the note")
    t.equal(summary.get("ROI_Max_Depth", ""), "", "whole image: no maximum depth")
    t.equal(summary["C1_DAPI_Object_Count"], len(DEPTH_DISKS), "whole image: objects are still counted")
    t.equal(set([row.get("Object_Depth_From_Outline", "") for row in result["object_rows"]]), set([""]), "whole image: no object depth")
    text = read_utf8(os.path.join(result["run_dir"], "depth_profile.csv"))
    t.check("covers the whole image" in text, "whole image: the note is in the file")

    result = depth_run(E, imp, configs, [whole], tmp, "whole_surface", {"depth_image_border_is_surface": True})
    flat = depth_rows(result, 3)
    t.equal(depth_column(flat, "Bin_Pixels"), depth_rect_bins(WIDTH, HEIGHT, 10), "whole image, border is outline: rectangular frames")
    t.equal(flat[0]["Bin_Pixels"], 240 * 200 - 220 * 180, "whole image, border is outline: bin 1")
    t.close(sum(depth_column(flat, "Bin_Area")), result["summary_rows"][0]["ROI_Area"], "whole image, border is outline: areas add up")
    t.equal(sum(depth_column(depth_rows(result, 1), "Object_Count")), len(DEPTH_DISKS), "whole image, border is outline: every object is in a bin")
    imp.close()


def test_depth_corner_triangle(E, t, tmp):
    """A triangle in the corner of the image, cut by two image borders: the
    distance map and the bins are compared, pixel by pixel, with distances
    worked out from the coordinates, for both settings of the border."""
    polygon = Polygon()
    polygon.addPoint(0, 0)
    polygon.addPoint(40, 0)
    polygon.addPoint(0, 30)
    roi = PolygonRoi(polygon, Roi.POLYGON)
    roi.setName("corner")
    inside = E.roi_inside_mask(roi, WIDTH, HEIGHT)
    roi_pixels = inside.getHistogram()[255]
    t.check(550 <= roi_pixels <= 650, "the triangle holds about 600 pixels (got " + str(roi_pixels) + ")")
    configs = depth_configs(E)
    imp = depth_image("depth_corner.tif")
    for border_is_surface in [False, True]:
        label = "corner triangle, border " + ("is outline" if border_is_surface else "not outline")
        truth = depth_brute_force(inside, (-1, -1, 41, 31), border_is_surface)
        t.equal(len(truth), roi_pixels, label + ": the check covers every ROI pixel")

        geometry = E.depth_roi_geometry(inside, roi.getBounds(), WIDTH, HEIGHT)
        E.depth_distance_map(geometry, border_is_surface)
        t.equal((geometry["x"], geometry["y"]), (0, 0), label + ": no margin where the ROI is on the image border")
        worst = 0.0
        for key in truth:
            worst = max(worst, abs(geometry["edm"].getf(key[0] - geometry["x"], key[1] - geometry["y"]) - 0.5 - truth[key]))
        t.check(worst < 0.0001, label + ": every pixel depth equals the distance worked out from coordinates (worst " + repr(worst) + ")")
        t.close((geometry["max_edm"] - 0.5), max(truth.values()), label + ": largest depth", 0.0001)

        # bins of 2.5 um = 5 pixels
        result = depth_run(E, imp, configs, [roi], tmp, "corner_" + str(border_is_surface), {
            "depth_bin_width": "2.5", "depth_image_border_is_surface": border_is_surface
        })
        expected = {}
        for key in truth:
            index = int(math.floor(truth[key] / 5.0))
            expected[index] = expected.get(index, 0) + 1
        flat = depth_rows(result, 3)
        t.equal(depth_column(flat, "Bin_Pixels"), [expected.get(index, 0) for index in range(max(expected.keys()) + 1)], label + ": bin pixels")
        summary = result["summary_rows"][0]
        t.close(sum(depth_column(flat, "Bin_Area")), summary["ROI_Area"], label + ": bin areas add up to ROI_Area")
        t.equal(summary["ROI_Touches_Image_Border"], 1, label + ": flagged as touching the image border")
        t.close(summary["ROI_Max_Depth"], max(truth.values()) * PIXEL, label + ": ROI_Max_Depth", 0.0001)
        t.close(summary["ROI_Mean_Depth"], sum(truth.values()) * PIXEL / len(truth), label + ": ROI_Mean_Depth", 0.0001)
        # The disk at (8, 8) is the only one in the triangle.
        nuclei = depth_rows(result, 1)
        t.equal(sum(depth_column(nuclei, "Object_Count")), 1, label + ": one object")
        row = depth_object_row(result, 8, 8)
        t.check(row is not None, label + ": the disk at (8, 8) has an object row")
        if row is not None:
            t.close(row["Object_Depth_From_Outline"], truth[(8, 8)] * PIXEL, label + ": object depth", 0.0001)
            t.equal(row["Object_Depth_Bin"], int(math.floor(truth[(8, 8)] / 5.0)) + 1, label + ": object bin")
    imp.close()
    # With the border as outline the pixel (8, 8) is 8.5 pixels deep (nearest:
    # the positions just beyond the top or left border); without, only the
    # slanted side counts. That side is the line 30 x + 40 y = 1200; the
    # centre (8.5, 8.5) of the pixel is (1200 - 595) / 50 = 12.1 pixels from
    # it, and along a slanted outline depth is within 0.75 pixel of that.
    t.close(depth_brute_force(inside, (-1, -1, 41, 31), True)[(8, 8)], 8.5, "the construction: (8, 8) is 8.5 pixels from the border")
    slanted = depth_brute_force(inside, (-1, -1, 41, 31), False)[(8, 8)]
    t.check(abs(slanted - 12.1) <= 0.75, "the construction: (8, 8) is 12.1 pixels from the slanted side (got " + repr(slanted) + ")")


def test_depth_hole(E, t, tmp):
    """A ROI with a hole: the hole is outline too. Checked pixel by pixel
    against the distance computed from the coordinates."""
    hole = (100, 80, 40, 40)
    roi = getattr(ShapeRoi(Roi(*DEPTH_RECT)), "not")(ShapeRoi(Roi(*hole)))
    roi.setName("with_hole")
    imp = depth_image("depth_hole.tif")
    result = depth_run(E, imp, depth_configs(E), [roi], tmp, "hole")
    imp.close()

    # Distance from a pixel centre to the centre of the nearest pixel outside
    # the ROI: the nearest pixel beyond a side of the rectangle, or the
    # nearest pixel of the hole, whichever is closer. Minus half a pixel.
    expected = {}
    deepest = 0.0
    for y in range(20, 180):
        for x in range(20, 220):
            dx = max(100 - x, 0, x - 139)
            dy = max(80 - y, 0, y - 119)
            if dx == 0 and dy == 0:
                continue
            depth = min(x - 19, 220 - x, y - 19, 180 - y, math.sqrt(dx * dx + dy * dy)) - 0.5
            deepest = max(deepest, depth)
            index = int(math.floor(depth / 10.0))
            expected[index] = expected.get(index, 0) + 1
    flat = depth_rows(result, 3)
    t.equal(depth_column(flat, "Bin_Pixels"), [expected.get(index, 0) for index in range(max(expected.keys()) + 1)],
            "ROI with a hole: every bin matches the distances computed from coordinates")
    t.equal(sum(depth_column(flat, "Bin_Pixels")), 200 * 160 - 40 * 40, "ROI with a hole: the hole is not part of the ROI")
    t.close(sum(depth_column(flat, "Bin_Area")), result["summary_rows"][0]["ROI_Area"], "ROI with a hole: areas add up to ROI_Area")
    t.close(result["summary_rows"][0]["ROI_Max_Depth"], deepest * PIXEL, "ROI with a hole: maximum depth")
    # Around the hole, bin 1 holds four 40 x 10 strips plus the rounded
    # corners: pixels (a, b) from a hole corner with sqrt(a^2 + b^2) < 10.5.
    corner = len([1 for a in range(1, 11) for b in range(1, 11) if a * a + b * b < 10.5 * 10.5])
    t.equal(flat[0]["Bin_Pixels"], (200 * 160 - 180 * 140) + 4 * 400 + 4 * corner, "ROI with a hole: bin 1 includes a rounded band around the hole")
    t.equal(flat[0]["ROI_Touches_Image_Border"], 0, "ROI with a hole: not on the image border")
    # The disk at (120,32) is 12 rings from the top side; the disk at
    # (120,100) lies in the hole and is not part of the ROI at all.
    row = depth_object_row(result, 120, 32)
    t.check(row is not None and row.get("Object_Depth_Bin") == 2, "ROI with a hole: an object far from the hole keeps its bin")
    t.equal(sum(depth_column(depth_rows(result, 1), "Object_Count")), result["summary_rows"][0]["C1_DAPI_Object_Count"],
            "ROI with a hole: object counts add up")


def test_depth_bins_and_degenerate(E, t, tmp):
    """Maximum depth, very narrow and very wide bins, a one-pixel ROI, an
    empty mask, pixels that are not square, a 32-bit channel."""
    configs = depth_configs(E)
    rect = Roi(*DEPTH_RECT)
    imp = depth_image("depth_bins.tif")

    # Maximum depth 20 um = 4 bins, then one open bin from ring 40 inward.
    result = depth_run(E, imp, configs, [rect], tmp, "cap", {"depth_max_depth": "20"})
    flat = depth_rows(result, 3)
    regular = depth_rect_bins(200, 160, 10)
    t.equal(depth_column(flat, "Bin_Pixels"), regular[:4] + [120 * 80], "maximum depth: deeper pixels go into one last bin")
    t.equal(depth_column(flat, "Note"), [""] * 4 + [E.DEPTH_NOTE_OPEN_BIN], "maximum depth: the last bin is marked open-ended")
    t.close(flat[4]["Depth_Start"], 20.0, "maximum depth: the open bin starts at the limit")
    t.close(flat[4]["Depth_End"], 39.75, "maximum depth: the open bin ends at the maximum depth of the ROI")
    t.close(flat[4]["Relative_Depth_End"], 1.0, "maximum depth: the open bin ends at relative depth 1")
    t.close(sum(depth_column(flat, "Bin_Area")), result["summary_rows"][0]["ROI_Area"], "maximum depth: areas still add up to ROI_Area")
    # (60,60) at 20.25 um and (120,100) at 39.75 um both fall in the open bin
    t.equal(depth_column(depth_rows(result, 1), "Object_Count"), [1, 2, 1, 0, 2], "maximum depth: objects per bin")
    t.equal(depth_object_row(result, 120, 100)["Object_Depth_Bin"], 5, "maximum depth: the deepest object is in the open bin")
    t.close(depth_object_row(result, 120, 100)["Object_Depth_From_Outline"], 39.75, "maximum depth: its depth is still reported exactly")

    # A limit that is not a multiple of the bin width is rounded up: 12 -> 15.
    result = depth_run(E, imp, configs, [rect], tmp, "cap_round", {"depth_max_depth": "12"})
    flat = depth_rows(result, 3)
    t.equal(depth_column(flat, "Bin_Pixels"), regular[:3] + [140 * 100], "maximum depth is rounded up to a whole bin")
    t.close(flat[3]["Depth_Start"], 15.0, "the open bin then starts at 15")

    # A limit beyond the ROI changes nothing.
    result = depth_run(E, imp, configs, [rect], tmp, "cap_far", {"depth_max_depth": "500"})
    t.equal(depth_column(depth_rows(result, 3), "Bin_Pixels"), regular, "a maximum depth beyond the ROI changes nothing")
    t.equal(set(depth_column(depth_rows(result, 3), "Note")), set([""]), "and no bin is marked open-ended")

    # Bins of a quarter pixel: ring k (depth k + 0.5 pixel) is bin 4 k + 3,
    # 319 bins in all, most of them empty. More than 255 bins also exercises
    # the path that measures intensity in groups of bins.
    result = depth_run(E, imp, configs, [rect], tmp, "narrow", {"depth_bin_width": "0.125"})
    flat = depth_rows(result, 3)
    ramp = depth_rows(result, 4)
    band = depth_rows(result, 2)
    t.equal(len(flat), 4 * 79 + 3, "quarter-pixel bins: 319 bins")
    wrong = 0
    for index in range(len(flat)):
        k = (index - 2) // 4
        ring = depth_rect_ring(200, 160, k) if (index - 2) % 4 == 0 else 0
        ok = flat[index]["Bin_Pixels"] == ring and band[index]["Positive_Pixels"] == (ring if k < 20 else 0)
        if ring > 0:
            ok = ok and abs(flat[index]["Mean_Intensity"] - 500.0) < 1e-9 and abs(ramp[index]["Mean_Intensity"] - 23907.0) < 1e-6
        else:
            ok = ok and flat[index]["Mean_Intensity"] != flat[index]["Mean_Intensity"]
        if not ok:
            wrong += 1
    t.equal(wrong, 0, "quarter-pixel bins: every ring is alone in its bin with the right pixels, positives and means")
    t.equal(depth_object_row(result, 27, 100)["Object_Depth_Bin"], 4 * 7 + 3, "quarter-pixel bins: object bin")
    text = read_utf8(os.path.join(result["run_dir"], "depth_profile.csv"))
    t.check("nan" not in text.lower(), "an empty bin is written as an empty cell, not as NaN")

    # A bin width that would give tens of thousands of bins is refused. The
    # deepest ring is at 39.75 um = 39750 x 0.001 um, exactly on the lower
    # edge of bin 39751, and a depth on an edge belongs to the deeper bin.
    result = depth_run(E, imp, configs, [rect], tmp, "too_many", {"depth_bin_width": "0.001"})
    table = result["extra_tables"]["depth_profile.csv"]["rows"]
    t.equal(len(table), 4, "too many bins: one note row per channel")
    t.check("would give 39751 bins" in table[0]["Note"], "too many bins: the note names the number (" + table[0]["Note"] + ")")
    t.equal(result["summary_rows"][0]["C1_DAPI_Object_Count"], 6, "too many bins: the standard results are still there")

    # One pixel: one bin, depth half a pixel.
    result = depth_run(E, imp, configs, [Roi(150, 150, 1, 1)], tmp, "one_pixel")
    flat = depth_rows(result, 3)
    t.equal(depth_column(flat, "Bin_Pixels"), [1], "one-pixel ROI: one bin with one pixel")
    t.close(result["summary_rows"][0]["ROI_Max_Depth"], 0.25, "one-pixel ROI: depth is half a pixel")
    t.close(flat[0]["Relative_Depth_End"], 1.0, "one-pixel ROI: relative depth ends at 1")
    nuclei = depth_rows(result, 1)
    t.equal(depth_column(nuclei, "Object_Count"), [0], "one-pixel ROI without objects: count 0")
    t.close(nuclei[0]["Object_Density_Per_Area"], 0.0, "one-pixel ROI without objects: density 0")
    median = result["summary_rows"][0]["C1_DAPI_Median_Object_Depth"]
    t.check(median != median, "no objects: the median object depth is empty")

    # An empty mask: the bright band of the marker is 30 pixels wide, so 16
    # rounds of erosion remove all of it and nothing is positive anywhere.
    result = depth_run(E, imp, depth_configs(E, {"erode_iterations": 16}), [rect], tmp, "empty_mask")
    band = depth_rows(result, 2)
    t.equal(result["summary_rows"][0]["C2_MAP2_Positive_Area"], 0.0, "empty mask: roi_summary has no positive area either")
    t.equal(depth_column(band, "Positive_Pixels"), [0] * 8, "empty mask: no positive pixels in any bin")
    t.equal(depth_column(band, "Positive_Area_Fraction_Percent"), [0.0] * 8, "empty mask: fraction 0 in every bin")
    t.equal(depth_column(band, "Mean_Intensity"), [900.0, 900.0] + [100.0] * 6, "empty mask: the intensity profile is still measured")
    mean_depth = result["summary_rows"][0]["C2_MAP2_Positive_Area_Mean_Depth"]
    t.check(mean_depth != mean_depth, "empty mask: the mean depth of positive pixels is empty")

    # A ROI in the dark middle of the marker with no disk in it (x 130..189,
    # y 60..139, 60 x 80 pixels = 30 rings = 3 bins): nothing to count.
    result = depth_run(E, imp, configs, [Roi(130, 60, 60, 80)], tmp, "nothing_inside")
    t.equal(result["summary_rows"][0]["C1_DAPI_Object_Count"], 0, "ROI without objects: the standard count is 0")
    nuclei = depth_rows(result, 1)
    t.equal(depth_column(nuclei, "Bin_Pixels"), depth_rect_bins(60, 80, 10), "ROI without objects: three rectangular frames")
    t.equal(depth_column(nuclei, "Object_Count"), [0, 0, 0], "ROI without objects: count 0 in every bin")
    t.equal(depth_column(nuclei, "Object_Density_Per_mm2"), [0.0, 0.0, 0.0], "ROI without objects: density 0 in every bin")
    t.equal(depth_column(depth_rows(result, 2), "Positive_Pixels"), [0, 0, 0], "ROI without marker: no positive pixels in any bin")

    # Pixels that are not square: a note instead of numbers.
    result = depth_run(E, imp, configs, [rect], tmp, "not_square", {
        "override_calibration": True, "pixel_width": "0.5", "pixel_height": "0.6", "spatial_unit": "um"
    })
    table = result["extra_tables"]["depth_profile.csv"]["rows"]
    t.equal(set(depth_column(table, "Note")), set([E.DEPTH_NOTE_NOT_SQUARE]), "pixels that are not square: a note")
    t.equal(len(table), 4, "pixels that are not square: one note row per channel")
    t.close(result["summary_rows"][0]["ROI_Area"], 200 * 160 * 0.5 * 0.6, "pixels that are not square: the standard results are still there")
    imp.close()

    # No calibration: depth and bin width are in pixels (default width 25),
    # and there is no density per mm2. The ramp channel is the median of two
    # identical Z slices, which makes it 32-bit (an "Average Z" projection
    # stays 16-bit in this version, so it would not test this): its bins are
    # measured with ImageJ's own statistics and must give the same means.
    imp = depth_image("depth_pixels.tif", slices=2, calibrated=False)
    ramp_configs = depth_configs(E, ramp_projection="Median Z")
    probe = E.AnalysisEngine().prepare_channels(imp, [ramp_configs[3]], 1, None)
    t.equal(probe[4].getBitDepth(), 32, "the median projection of the ramp is 32-bit")
    probe[4].close()
    result = depth_run(E, imp, ramp_configs, [rect], tmp, "pixels", {"depth_bin_width": "25"})
    imp.close()
    flat = depth_rows(result, 3)
    ramp = depth_rows(result, 4)
    t.equal(depth_column(flat, "Bin_Pixels"), depth_rect_bins(200, 160, 25), "no calibration: bins of 25 pixels")
    t.equal(depth_column(flat, "Bin_Pixels"), [15500, 10500, 5500, 500], "no calibration: the four frames by hand")
    t.equal(set(depth_column(flat, "Depth_Unit")), set(["pixel"]), "no calibration: the unit is the pixel")
    t.close(result["summary_rows"][0]["ROI_Max_Depth"], 79.5, "no calibration: maximum depth in pixels")
    t.equal(set(depth_column(depth_rows(result, 1), "Object_Density_Per_mm2")), set([""]), "no calibration: no density per mm2")
    for row in ramp:
        t.close(row["Mean_Intensity"], 23907.0, "32-bit channel, bin " + str(row["Bin_Index"]))

    # Depths that lie exactly on a bin edge. Bins of 2.3 pixel: ring k (depth
    # k + 0.5 pixel) belongs to bin floor((k + 0.5) / 2.3) + 1, which in whole
    # numbers is (20 k + 10) // 46 + 1. Ring 34 has depth 34.5 = 15 x 2.3,
    # exactly the lower edge of bin 16, and must be counted there; so must
    # the disk whose centre pixel (54, 100) lies in ring 34.
    imp = depth_image("depth_edges.tif", [(54, 100)], calibrated=False)
    result = depth_run(E, imp, configs, [rect], tmp, "edges", {"depth_bin_width": "2.3"})
    imp.close()
    expected = [0] * ((20 * 79 + 10) // 46 + 1)
    for k in range(80):
        expected[(20 * k + 10) // 46] += depth_rect_ring(200, 160, k)
    flat = depth_rows(result, 3)
    t.equal(len(expected), 35, "bins of 2.3 pixel: the construction gives 35 bins")
    t.equal(depth_column(flat, "Bin_Pixels"), expected, "bins of 2.3 pixel: a depth exactly on a bin edge goes into the deeper bin")
    t.equal((20 * 34 + 10) % 46, 0, "ring 34 lies exactly on a bin edge")
    # Bin 16 covers depths 34.5 to 36.8 pixel: rings 34, 35 and 36.
    t.equal(flat[15]["Bin_Pixels"], sum([depth_rect_ring(200, 160, k) for k in [34, 35, 36]]), "bin 16 holds rings 34, 35 and 36")
    t.close(flat[15]["Depth_Start"], 34.5, "bin 16 starts at depth 34.5 pixel")
    row = depth_object_row_pixels(result, 54, 100)
    t.check(row is not None, "the disk on the bin edge has an object row")
    if row is not None:
        t.close(row["Object_Depth_From_Outline"], 34.5, "the disk on the bin edge has depth 34.5 pixel")
        t.equal(row["Object_Depth_Bin"], 16, "the disk on the bin edge is in the deeper bin")


def test_depth_crescent_object(E, t, tmp):
    """An object whose centroid lies outside its own outline (a C shape):
    its depth is the mean depth of its pixels, and the row says so."""
    plane = blank_plane(100)
    draw_disk(plane, 60, 100, 20, 2000)
    draw_disk(plane, 60, 100, 14, 100)
    draw_rect(plane, (60, 96, 25, 9), 100)
    imp = depth_image("depth_crescent.tif", first_plane=plane)
    result = depth_run(E, imp, depth_configs(E), [Roi(*DEPTH_RECT)], tmp, "crescent")
    imp.close()
    t.equal(len(result["object_rows"]), 1, "the C shape is one object")
    row = result["object_rows"][0]
    t.equal(row["Object_Depth_Measured_At"], E.DEPTH_AT_OBJECT_MEAN, "the row says the object mean was used")
    # Every pixel of the C is nearer to the left side of the ROI (x = 20)
    # than to any other side, so depth = column - 20 + 0.5 pixel: a straight
    # line in the column. The mean of a straight line over the object is its
    # value at the centroid: (Centroid_X in pixels) - 20, in um Centroid_X - 10.
    t.close(row["Object_Depth_From_Outline"], row["Centroid_X"] - 10.0, "mean depth of the C shape")
    t.equal(row["Object_Depth_Bin"], int(math.floor(row["Object_Depth_From_Outline"] / 5.0)) + 1, "bin of the C shape")
    nuclei = depth_rows(result, 1)
    t.equal(sum(depth_column(nuclei, "Object_Count")), 1, "the C shape is counted once")
    t.equal(nuclei[row["Object_Depth_Bin"] - 1]["Object_Count"], 1, "and in the bin its row names")


def test_depth_frames_and_batch(E, t, tmp):
    """Two frames with different objects, and a folder batch with one ROI
    file per image."""
    configs = depth_configs(E)
    rect = Roi(*DEPTH_RECT)
    rect.setName("section")
    imp = depth_image("depth_frames.tif", [[(27, 100), (45, 100)], [(120, 100)]], frames=2)
    result = depth_run(E, imp, configs, [rect, triangle_roi()], tmp, "frames", frames=[1, 2])
    imp.close()
    table = result["extra_tables"]["depth_profile.csv"]["rows"]
    t.equal(set(depth_column(table, "Frame")), set([1, 2]), "both frames are in the table")
    t.equal(set(depth_column(table, "Parent_ROI_Name")), set(["section", "triangle"]), "both parent ROIs are in the table")
    t.equal(depth_column(depth_rows(result, 1, 1, 1), "Object_Count"), [1, 0, 1, 0, 0, 0, 0, 0], "frame 1: objects per bin")
    t.equal(depth_column(depth_rows(result, 1, 1, 2), "Object_Count"), [0, 0, 0, 0, 0, 0, 0, 1], "frame 2: objects per bin")
    t.equal(depth_column(depth_rows(result, 3, 1, 1), "Bin_Pixels"), depth_column(depth_rows(result, 3, 1, 2), "Bin_Pixels"),
            "the bins of a ROI are the same in every frame")
    for index in range(len(result["summary_rows"])):
        summary = result["summary_rows"][index]
        rows = depth_rows(result, 1, summary["Parent_ROI_Index"], summary["Frame"])
        t.close(sum(depth_column(rows, "Bin_Area")), summary["ROI_Area"], "summary row " + str(index + 1) + ": bin areas add up to ROI_Area")
        t.equal(sum(depth_column(rows, "Object_Count")), summary["C1_DAPI_Object_Count"], "summary row " + str(index + 1) + ": object counts add up")
        band = depth_rows(result, 2, summary["Parent_ROI_Index"], summary["Frame"])
        t.close(sum(depth_column(band, "Positive_Pixels")) * PIXEL * PIXEL, summary["C2_MAP2_Positive_Area"],
                "summary row " + str(index + 1) + ": positive pixels add up to Positive_Area")
        weighted = sum([row["Mean_Intensity"] * row["Bin_Pixels"] for row in band if row["Bin_Pixels"] > 0])
        t.close(weighted / sum(depth_column(band, "Bin_Pixels")), summary["C2_MAP2_Mean"],
                "summary row " + str(index + 1) + ": bin means weighted by pixels give the ROI mean")

    folder = os.path.join(tmp, "depth_batch")
    os.makedirs(folder)
    for index in range(2):
        imp = depth_image("d" + str(index) + ".tif")
        FileSaver(imp).saveAsTiff(os.path.join(folder, "d" + str(index) + ".tif"))
        imp.close()
        E.save_roi_array([E.clone_roi(rect)], os.path.join(folder, "d" + str(index) + "_ROIs.zip"))
    options = E.default_global_options()
    options["depth_enabled"] = True
    options["depth_bin_width"] = "5"
    for key in DEPTH_SAVE_KEYS:
        options[key] = False
    batch = E.run_batch(folder, os.path.join(tmp, "depth_batch_out"), configs, configs, options,
                        {"extensions": "tif", "roi_policy": E.BATCH_ROI_FILE_SKIP})
    t.equal(batch["analyzed"], 2, "the batch analyzed both images")
    combined = os.path.join(batch["batch_dir"], "batch_depth_profile.csv")
    t.check(os.path.isfile(combined), "the batch writes batch_depth_profile.csv")
    if os.path.isfile(combined):
        lines = [line for line in read_utf8(combined).split("\n") if line.strip() != ""]
        t.equal(len(lines), 1 + 2 * 4 * 8, "the combined depth table has one row per image, channel and bin")
        t.check(lines[0].endswith('"Note","Batch_ID","Image_Relative_Path"'), "the combined depth table carries the batch columns")
        t.check('"d0.tif"' in lines[1] and '"d1.tif"' in lines[-1], "its rows name their image")
    header = read_utf8(os.path.join(batch["batch_dir"], "batch_roi_summary.csv")).split("\n")[0]
    t.check('"ROI_Max_Depth","ROI_Mean_Depth","ROI_Touches_Image_Border","Depth_Profile_Note"' in header,
            "the combined summary carries the depth columns")
    header = read_utf8(os.path.join(batch["batch_dir"], "batch_objects.csv")).split("\n")[0]
    t.check('"Object_Depth_Measured_At","Batch_ID","Image_Relative_Path"' in header, "the combined object table carries the depth columns")


# -----------------------------------------------------------------------------
# Known-answer cases added by the independent check of the module
# -----------------------------------------------------------------------------

DEPTH_BYTE_DISKS = [(90, 90), (20, 60), (50, 125), (8, 90), (150, 100)]


def depth_byte_image(title):
    """Two 8-bit channels at 0.25 micron per pixel, with the unit written
    "micron" and a stage origin of (7, 3) pixels, as Bio-Formats may set.
      C1 "DAPI"  disks of radius 6, value 200 on 20, at DEPTH_BYTE_DISKS.
      C2 "Ramp"  value = column (0 .. 239)."""
    nuclei = ByteProcessor(WIDTH, HEIGHT)
    nuclei.setValue(20.0)
    nuclei.fill()
    for disk in DEPTH_BYTE_DISKS:
        draw_disk(nuclei, disk[0], disk[1], DISK_RADIUS, 200)
    ramp = ByteProcessor(WIDTH, HEIGHT)
    for y in range(HEIGHT):
        for x in range(WIDTH):
            ramp.set(x, y, x)
    stack = ImageStack(WIDTH, HEIGHT)
    stack.addSlice("", nuclei)
    stack.addSlice("", ramp)
    imp = ImagePlus(title, stack)
    imp.setDimensions(2, 1, 1)
    imp = CompositeImage(imp, CompositeImage.GRAYSCALE)
    cal = imp.getCalibration()
    cal.pixelWidth = 0.25
    cal.pixelHeight = 0.25
    cal.setUnit("micron")
    cal.xOrigin = 7.0
    cal.yOrigin = 3.0
    return imp


def depth_byte_object(result, cx, cy):
    """Object row of the disk around pixel (cx, cy) of depth_byte_image:
    calibrated centroid = (pixel centre - origin) x 0.25."""
    for row in result["object_rows"]:
        if (abs(row["Centroid_X"] - (cx + 0.5 - 7.0) * 0.25) < 0.001 and
                abs(row["Centroid_Y"] - (cy + 0.5 - 3.0) * 0.25) < 0.001):
            return row
    return None


def test_depth_byte_micron_cut(E, t, tmp):
    """8-bit image, 0.25 micron pixels with the unit "micron" and a stage
    origin, ROI reaching past the left image border. Every number derived
    by hand below.

    ROI Roi(-30, 40, 130, 100) is clipped by the image to x 0..99, y 40..139
    (100 x 100 pixels, 625 square micron). The left side is the image border.
    Border NOT outline (default): the nearest outside pixel is beyond the
    right side (x = 100), the top (y = 39) or the bottom (y = 140), so
        depth = min(99.5 - x, y - 39.5, 139.5 - y) pixels.
    Pixels at least 10 j pixels deep form the block x 0..99-10j,
    y 40+10j..139-10j of (100 - 10 j)(100 - 20 j) pixels: 10000, 7200, 4800,
    2800, 1200, 0. Bins of 2.5 micron = 10 pixels are the differences:
    2800, 2400, 2000, 1600, 1200. The deepest pixels (x <= 50, y = 89 or 90)
    are 49.5 pixels = 12.375 micron deep.
    Ramp sum over a block = rows x (0 + ... + (99 - 10j)):
    495000, 320400, 189600, 96600, 35400 -> bin sums 174600, 130800, 93000,
    61200, 35400 -> means 62.357.., 54.5, 46.5, 38.25, 29.5.
    Positive (ramp >= 50) pixels of a block: rows x (50 - 10 j): 5000, 3200,
    1800, 800, 200 -> per bin 1800, 1400, 1000, 600, 200.
    Mean depth: with m = depth - 0.5, sum of m = sum over t = 1..49 of
    (100 - t)(100 - 2t) = 203350, so the mean is (203350 + 5000) / 10000 =
    20.835 pixels = 5.20875 micron. Over the positive pixels: sum over t of
    (50 - t)(100 - 2t) = 2 (1^2 + ... + 49^2) = 80850 -> (80850 + 2500) /
    5000 = 16.67 pixels = 4.1675 micron.
    Disks: (90,90) 9.5 px -> bin 1; (50,125) 14.5 -> bin 2; (20,60) 20.5 ->
    bin 3; (8,90) 49.5 -> bin 5; (150,100) lies outside the ROI. Median of
    2.375, 3.625, 5.125, 12.375 micron = 4.375.
    Border IS outline: a plain 100 x 100 square; ring k holds 396 - 8 k
    pixels, so bins of 10 rings hold 3600, 2800, 2000, 1200, 400; the disk at
    (8,90) is then 8.5 pixels deep (bin 1), the other three keep their
    depths: counts 2, 1, 1, 0, 0."""
    imp = depth_byte_image("depth_byte.tif")
    configs = [
        nuclei_config(E, {"manual_low": 100}),
        E.normalize_channel_config(2, {
            "enabled": True, "name": "Ramp", "analysis_mode": "Positive area",
            "threshold_mode": "Manual", "manual_low": 50, "manual_high": -1
        })
    ]
    t.equal(imp.getBitDepth(), 8, "the test image is 8-bit")
    roi = Roi(-30, 40, 130, 100)
    result = depth_run(E, imp, configs, [roi], tmp, "byte_cut", {"depth_bin_width": "2.5"})
    surface = depth_run(E, imp, configs, [roi], tmp, "byte_cut_surface", {
        "depth_bin_width": "2.5", "depth_image_border_is_surface": True
    })
    imp.close()
    summary = result["summary_rows"][0]

    t.close(summary["ROI_Area"], 625.0, "the ROI is clipped by the image to 100 x 100 pixels")
    ramp = depth_rows(result, 2)
    t.equal(depth_column(ramp, "Bin_Pixels"), [2800, 2400, 2000, 1600, 1200], "8-bit, cut ROI: bin pixels by hand")
    for index in range(5):
        t.close(ramp[index]["Bin_Area"], [175.0, 150.0, 125.0, 100.0, 75.0][index], "bin " + str(index + 1) + " area in square micron")
        t.close(ramp[index]["Depth_Start"], 2.5 * index, "bin " + str(index + 1) + " starts at a multiple of 2.5 micron")
        t.close(ramp[index]["Relative_Depth_Start"], 2.5 * index / 12.375, "bin " + str(index + 1) + " relative start")
    t.close(ramp[4]["Depth_End"], 12.5, "the last bin ends at 12.5 micron")
    t.close(ramp[4]["Relative_Depth_End"], 1.0, "and at relative depth 1")
    t.close(summary["ROI_Max_Depth"], 12.375, "maximum depth 49.5 pixels = 12.375 micron")
    t.close(summary["ROI_Mean_Depth"], 5.20875, "mean depth by hand")
    t.equal(summary["ROI_Touches_Image_Border"], 1, "the ROI is flagged as cut by the image border")
    t.equal(ramp[0]["ROI_Border_Contact_Pixels"], 100, "100 ROI pixels lie in image column 0")
    t.equal(ramp[0]["Depth_Unit"], "micron", "the unit is the image's own")
    sums = [174600.0, 130800.0, 93000.0, 61200.0, 35400.0]
    pixels = [2800, 2400, 2000, 1600, 1200]
    for index in range(5):
        t.close(ramp[index]["Mean_Intensity"], sums[index] / pixels[index], "8-bit ramp mean, bin " + str(index + 1))
    t.equal(depth_column(ramp, "Positive_Pixels"), [1800, 1400, 1000, 600, 200], "positive pixels per bin by hand")
    for index in range(5):
        t.close(ramp[index]["Positive_Area_Fraction_Percent"], 100.0 * [1800, 1400, 1000, 600, 200][index] / pixels[index],
                "positive fraction, bin " + str(index + 1))
    t.close(summary["C2_Ramp_Positive_Area"], 5000 * 0.0625, "roi_summary agrees: 5000 positive pixels")
    t.close(summary["C2_Ramp_Positive_Area_Mean_Depth"], 4.1675, "mean depth of the positive pixels by hand")

    nuclei = depth_rows(result, 1)
    t.equal(depth_column(nuclei, "Object_Count"), [1, 1, 1, 0, 1], "objects per bin by hand")
    t.equal(summary["C1_DAPI_Object_Count"], 4, "four disks lie in the ROI")
    t.close(nuclei[0]["Object_Density_Per_mm2"], 1.0e6 / 175.0, "density per mm2 with the unit written micron")
    expected = [((90, 90), 2.375, 1, 2.375, 1), ((50, 125), 3.625, 2, 3.625, 2),
                ((20, 60), 5.125, 3, 5.125, 3), ((8, 90), 12.375, 5, 2.125, 1)]
    for center, depth, bin_index, surface_depth, surface_bin in expected:
        row = depth_byte_object(result, center[0], center[1])
        other = depth_byte_object(surface, center[0], center[1])
        t.check(row is not None and other is not None, "the disk at " + str(center) + " is found by its calibrated centroid")
        if row is None or other is None:
            continue
        t.equal(row["Object_Depth_Measured_At"], E.DEPTH_AT_CENTROID, "the stage origin is undone: read at the centroid " + str(center))
        t.close(row["Object_Depth_From_Outline"], depth, "depth of the disk at " + str(center))
        t.equal(row["Object_Depth_Bin"], bin_index, "bin of the disk at " + str(center))
        t.close(row["Object_Relative_Depth"], depth / 12.375, "relative depth of the disk at " + str(center))
        t.close(other["Object_Depth_From_Outline"], surface_depth, "border as outline: depth of the disk at " + str(center))
        t.equal(other["Object_Depth_Bin"], surface_bin, "border as outline: bin of the disk at " + str(center))
        t.check(other["Object_Depth_From_Outline"] <= row["Object_Depth_From_Outline"],
                "the two border settings bracket the depth of the disk at " + str(center))
    t.close(summary["C1_DAPI_Median_Object_Depth"], 4.375, "median object depth by hand")

    ramp = depth_rows(surface, 2)
    t.equal(depth_column(ramp, "Bin_Pixels"), [3600, 2800, 2000, 1200, 400], "border as outline: square frames by hand")
    t.equal(depth_column(depth_rows(surface, 1), "Object_Count"), [2, 1, 1, 0, 0], "border as outline: objects per bin")
    t.close(surface["summary_rows"][0]["ROI_Max_Depth"], 12.375, "border as outline: maximum depth")


def test_depth_tiny_rois_two_frames(E, t, tmp):
    """Two parent ROIs, two frames, no objects at all in frame 2.

    ROI 1 "tiny" = Roi(150, 150, 5, 3), bins of 0.5 um = 1 pixel. Ring 0
    holds 5 x 3 - 3 x 1 = 12 pixels at depth 0.5 pixel (0.25 um), ring 1 the
    3 middle pixels at 1.5 pixels (0.75 um): bins 12 and 3, areas 3.0 and
    0.75 um2, mean depth (12 x 0.25 + 3 x 0.75) / 15 = 0.35 um, relative
    bin edges 0, 0.5 / 0.75, 1. Its columns 150..154 are symmetric about
    152, so the ramp mean is 200 x 152 + 7 = 30407 in both bins; it lies
    in the dark middle of the marker (100, nothing positive) and holds no
    disk. ROI 2 is the main rectangle: with 1-pixel bins the disk in ring r
    is in bin r + 1 (rings 7, 12, 14, 25, 40, 79 -> bins 8, 13, 15, 26, 41,
    80) in frame 1. In frame 2 the only disk is at (8, 8), outside both
    ROIs (it keeps the image maximum above the Manual threshold, without
    which the engine refuses the frame), so both ROIs have no object."""
    tiny = Roi(150, 150, 5, 3)
    tiny.setName("tiny")
    rect = Roi(*DEPTH_RECT)
    rect.setName("section")
    imp = depth_image("depth_tiny.tif", [DEPTH_DISKS, [(8, 8)]], frames=2)
    result = depth_run(E, imp, depth_configs(E), [tiny, rect], tmp, "tiny", {"depth_bin_width": "0.5"}, frames=[1, 2])
    imp.close()
    table = result["extra_tables"]["depth_profile.csv"]["rows"]
    t.equal(len(table), 2 * (2 + 80) * 4, "rows: 2 frames x (2 + 80 bins) x 4 channels")
    for frame in [1, 2]:
        label = "frame " + str(frame) + ", tiny ROI: "
        flat = depth_rows(result, 3, 1, frame)
        t.equal(depth_column(flat, "Bin_Pixels"), [12, 3], label + "bin pixels by hand")
        t.equal(depth_column(flat, "Bin_Area"), [3.0, 0.75], label + "bin areas")
        t.equal(depth_column(flat, "Parent_ROI_Name"), ["tiny", "tiny"], label + "named")
        t.close(flat[1]["Relative_Depth_Start"], 0.5 / 0.75, label + "relative start of bin 2")
        t.close(flat[0]["Relative_Depth_End"], 0.5 / 0.75, label + "relative end of bin 1")
        t.close(flat[1]["Relative_Depth_End"], 1.0, label + "relative end of bin 2")
        t.equal(depth_column(depth_rows(result, 4, 1, frame), "Mean_Intensity"), [30407.0, 30407.0], label + "ramp means")
        band = depth_rows(result, 2, 1, frame)
        t.equal(depth_column(band, "Positive_Pixels"), [0, 0], label + "nothing positive")
        t.equal(depth_column(band, "Mean_Intensity"), [100.0, 100.0], label + "marker means")
        nuclei = depth_rows(result, 1, 1, frame)
        t.equal(depth_column(nuclei, "Object_Count"), [0, 0], label + "no objects")
        t.equal(depth_column(nuclei, "Object_Density_Per_mm2"), [0.0, 0.0], label + "density 0")
    summaries = {}
    for summary in result["summary_rows"]:
        summaries[(summary["Frame"], summary["Parent_ROI_Index"])] = summary
    t.equal(sorted(summaries.keys()), [(1, 1), (1, 2), (2, 1), (2, 2)], "one summary row per frame and ROI")
    for frame in [1, 2]:
        t.close(summaries[(frame, 1)]["ROI_Max_Depth"], 0.75, "frame " + str(frame) + ": tiny ROI maximum depth")
        t.close(summaries[(frame, 1)]["ROI_Mean_Depth"], 0.35, "frame " + str(frame) + ": tiny ROI mean depth by hand")
        median = summaries[(frame, 1)]["C1_DAPI_Median_Object_Depth"]
        t.check(median != median, "frame " + str(frame) + ": tiny ROI has no median object depth")
        t.close(summaries[(frame, 2)]["ROI_Max_Depth"], 39.75, "frame " + str(frame) + ": rectangle maximum depth")
    counts = [0] * 80
    for ring in [7, 12, 14, 25, 40, 79]:
        counts[ring] += 1
    t.equal(depth_column(depth_rows(result, 1, 2, 1), "Object_Count"), counts, "frame 1: each disk in the bin of its ring")
    t.equal(depth_column(depth_rows(result, 1, 2, 2), "Object_Count"), [0] * 80, "frame 2: no objects in any bin")
    t.equal(summaries[(2, 2)]["C1_DAPI_Object_Count"], 0, "frame 2: roi_summary counts no object either")
    median = summaries[(2, 2)]["C1_DAPI_Median_Object_Depth"]
    t.check(median != median, "frame 2: no median object depth")
    t.close(summaries[(1, 2)]["C1_DAPI_Median_Object_Depth"], 10.0, "frame 1: median object depth as in the rectangle test")
    t.equal(depth_column(depth_rows(result, 3, 2, 1), "Bin_Pixels"), [depth_rect_ring(200, 160, k) for k in range(80)],
            "frame 1: one ring per 1-pixel bin")
    t.equal(depth_column(depth_rows(result, 3, 2, 2), "Bin_Pixels"), [depth_rect_ring(200, 160, k) for k in range(80)],
            "frame 2: the same rings")


def depth_star_roi(cx, cy):
    """A non-convex 7-pointed star (outer radius 45, inner 18, squashed to
    90 percent in y). Its reflex corners make the depth contours curved."""
    polygon = Polygon()
    for index in range(14):
        radius = 45 if index % 2 == 0 else 18
        angle = math.pi * 2 * index / 14.0 + 0.3
        polygon.addPoint(int(round(cx + radius * math.cos(angle))), int(round(cy + radius * 0.9 * math.sin(angle))))
    roi = PolygonRoi(polygon, Roi.POLYGON)
    roi.setName("star")
    return roi


def depth_nearest_outside(inside, box):
    """Depth in pixels (centre-to-centre distance to the nearest pixel that
    is not in the ROI, minus half a pixel) of every ROI pixel in box =
    (x0, y0, x1, y1), which must hold the ROI with a pixel to spare and lie
    inside the image. Only outside pixels with an ROI pixel beside them are
    tried: the nearest outside pixel is always one of them, because a step
    from it towards the ROI pixel along x or y is nearer still."""
    members = []
    frontier = []
    for y in range(box[1], box[3] + 1):
        for x in range(box[0], box[2] + 1):
            if inside.get(x, y) != 0:
                members.append((x, y))
            elif (inside.get(x - 1, y) != 0 or inside.get(x + 1, y) != 0 or
                  inside.get(x, y - 1) != 0 or inside.get(x, y + 1) != 0):
                frontier.append((x, y))
    depths = {}
    for member in members:
        best = None
        for other in frontier:
            dx = member[0] - other[0]
            dy = member[1] - other[1]
            squared = dx * dx + dy * dy
            if best is None or squared < best:
                best = squared
        depths[member] = math.sqrt(best) - 0.5
    return depths


def test_depth_star_32bit(E, t, tmp):
    """A non-convex star with a 32-bit channel: every bin, mean, positive
    count and object bin against depths found by searching every outside
    pixel, and the marker band worked out from how it is drawn.
    Bins of 1.5 um = 3 pixels: pixel depth d goes to bin floor(d / 3) + 1.
    No depth lies on a bin edge (d / 3 whole would need sqrt(n) = 3 k + 0.5).
    The marker is 900, i.e. positive, inside x 10..229, y 10..189 except
    x 40..199, y 40..159; the ramp is 200 x + 7."""
    roi = depth_star_roi(120, 62)
    inside = E.roi_inside_mask(roi, WIDTH, HEIGHT)
    bounds = roi.getBounds()
    t.check(bounds.x > 1 and bounds.y > 1 and bounds.x + bounds.width < WIDTH - 1 and bounds.y + bounds.height < HEIGHT - 1,
            "the star lies inside the image")
    truth = depth_nearest_outside(inside, (bounds.x - 1, bounds.y - 1, bounds.x + bounds.width, bounds.y + bounds.height))
    t.equal(len(truth), inside.getHistogram()[255], "the search covers every pixel of the star")

    imp = depth_image("depth_star.tif", slices=2)
    configs = depth_configs(E, ramp_projection="Median Z")
    result = depth_run(E, imp, configs, [roi], tmp, "star", {"depth_bin_width": "1.5"})
    imp.close()

    pixels = {}
    ramp_sums = {}
    positives = {}
    positive_depth = [0.0, 0]
    for key in truth:
        index = int(math.floor(truth[key] / 3.0))
        pixels[index] = pixels.get(index, 0) + 1
        ramp_sums[index] = ramp_sums.get(index, 0) + 200 * key[0] + 7
        band = 10 <= key[0] <= 229 and 10 <= key[1] <= 189 and not (40 <= key[0] <= 199 and 40 <= key[1] <= 159)
        if band:
            positives[index] = positives.get(index, 0) + 1
            positive_depth[0] += truth[key]
            positive_depth[1] += 1
    count = max(pixels.keys()) + 1
    t.check(count >= 4, "the star spans several bins (" + str(count) + ")")
    t.check(positive_depth[1] > 0 and positive_depth[1] < len(truth), "the marker band covers part of the star")
    flat = depth_rows(result, 3)
    ramp = depth_rows(result, 4)
    band = depth_rows(result, 2)
    t.equal(depth_column(flat, "Bin_Pixels"), [pixels.get(index, 0) for index in range(count)], "star: bin pixels")
    t.equal(depth_column(band, "Positive_Pixels"), [positives.get(index, 0) for index in range(count)], "star: positive pixels per bin")
    wrong = []
    for index in range(count):
        if pixels.get(index, 0) == 0:
            continue
        if abs(ramp[index]["Mean_Intensity"] - ramp_sums[index] / float(pixels[index])) > 1e-6 * ramp[index]["Mean_Intensity"]:
            wrong.append(index + 1)
        if abs(flat[index]["Mean_Intensity"] - 500.0) > 1e-9:
            wrong.append(index + 1)
    t.equal(wrong, [], "star: 32-bit ramp and 16-bit flat means in every bin")
    summary = result["summary_rows"][0]
    t.close(summary["ROI_Max_Depth"], max(truth.values()) * PIXEL, "star: maximum depth", 1e-5)
    t.close(summary["ROI_Mean_Depth"], sum(truth.values()) * PIXEL / len(truth), "star: mean depth", 1e-5)
    t.close(summary["C2_MAP2_Positive_Area_Mean_Depth"], positive_depth[0] * PIXEL / positive_depth[1],
            "star: mean depth of the positive pixels", 1e-5)
    t.equal(summary["ROI_Touches_Image_Border"], 0, "star: not on the image border")

    objects = [row for row in result["object_rows"] if row["Segmentation_Channel"] == 1]
    t.check(len(objects) > 0, "star: there are objects")
    expected = [0] * count
    mismatched = 0
    for row in objects:
        key = (int(math.floor(row["Centroid_X"] / PIXEL + 1e-7)), int(math.floor(row["Centroid_Y"] / PIXEL + 1e-7)))
        if row["Object_Depth_Measured_At"] == E.DEPTH_AT_CENTROID and key in truth:
            index = int(math.floor(truth[key] / 3.0))
            if row["Object_Depth_Bin"] != index + 1 or abs(row["Object_Depth_From_Outline"] - truth[key] * PIXEL) > 1e-5:
                mismatched += 1
        else:
            index = row["Object_Depth_Bin"] - 1
        expected[index] += 1
    t.equal(mismatched, 0, "star: every object read at its centroid has the depth and bin of that pixel")
    t.equal(depth_column(depth_rows(result, 1), "Object_Count"), expected, "star: objects per bin")


def test_depth_calibrated_bin_means(E, t, tmp):
    """Bin means of a 16-bit channel with a calibration function (signed
    16-bit: value = raw - 32768) must be calibrated, like roi_summary's
    Mean; an 8-bit channel goes through the exact histogram method. The
    bin codes are drawn by hand: column x of a 10 x 4 image is bin
    x // 4 + 1 (bins of 16, 16 and 8 pixels); raw value 40000 + x."""
    codes = ShortProcessor(10, 4)
    signed = ShortProcessor(10, 4)
    small = ByteProcessor(10, 4)
    for y in range(4):
        for x in range(10):
            codes.set(x, y, x // 4 + 1)
            signed.set(x, y, 40000 + x)
            small.set(x, y, 3 * x)
    signed_imp = ImagePlus("signed", signed)
    signed_imp.getCalibration().setSigned16BitCalibration()
    small_imp = ImagePlus("small", small)
    t.check(signed_imp.getCalibration().calibrated(), "the 16-bit channel has a calibration function")
    configs = [
        E.normalize_channel_config(1, {"enabled": True, "name": "Signed", "analysis_mode": "Intensity only"}),
        E.normalize_channel_config(2, {"enabled": True, "name": "Small", "analysis_mode": "Intensity only"})
    ]
    context = {"configs": configs, "prepared": {1: signed_imp, 2: small_imp}}
    geometry = {"x": 0, "y": 0, "width": 10, "height": 4}
    means = E.depth_bin_means(context, geometry, codes, [16, 16, 8])
    # columns 0-3 -> mean x 1.5; 4-7 -> 5.5; 8-9 -> 8.5
    t.equal(len(means[1]), 3, "one calibrated mean per bin")
    for index in range(3):
        t.close(means[1][index], [40001.5, 40005.5, 40008.5][index] - 32768, "calibrated 16-bit mean, bin " + str(index + 1), 1e-12)
        t.close(means[2][index], [4.5, 16.5, 25.5][index], "8-bit mean, bin " + str(index + 1), 1e-12)
    t.equal(signed.get(0, 0), 40000, "the channel image is left as it was")


def test_depth_profile_round_trip(E, t, tmp):
    """Every depth setting survives a save and load; the example profile
    shipped with the application (written before this module) loads with
    the depth defaults; true and false are read the way profiles write them."""
    configs = depth_configs(E)
    options = E.default_global_options()
    options["depth_enabled"] = True
    options["depth_bin_width"] = "12.5"
    options["depth_max_depth"] = "40"
    options["depth_image_border_is_surface"] = True
    path = os.path.join(tmp, "depth_round_trip.properties")
    E.write_settings_properties(path, options, configs)
    text = read_utf8(path)
    for line in ["depth.enabled=true", "depth.bin_width=12.5", "depth.max_depth=40", "depth.image_border_is_surface=true"]:
        t.check(line in text, "the profile holds " + line)
    loaded = E.load_profile_settings(path)[0]
    t.equal(E.depth_settings(loaded), {"enabled": True, "bin_width": 12.5, "max_depth": 40.0, "border_is_surface": True},
            "all depth settings load back")
    for key in ["depth_enabled", "depth_bin_width", "depth_max_depth", "depth_image_border_is_surface"]:
        t.equal(loaded[key], options[key], "option " + key + " loads back unchanged")
    example = os.path.join(E.CELLQUANT_PACKAGE_DIR, "examples", "Example_4Channel_SOX2_Casp3_MAP2_DAPI.properties")
    if os.path.isfile(example):
        t.check("depth." not in read_utf8(example), "the example profile has no depth keys")
        old = E.load_profile_settings(example)[0]
        t.equal(E.depth_settings(old), {"enabled": False, "bin_width": 25.0, "max_depth": 0.0, "border_is_surface": False},
                "the example profile loads with the depth defaults")
        t.check(not any([extension.name == "depth_profile" for extension in E.active_run_extensions(old, configs)]),
                "the example profile does not switch the depth profile on")
    else:
        t.check(False, "the example profile is missing: " + example)


TESTS.append(("depth profile: a rectangle gives exact frames, means, positives and object bins", test_depth_rectangle))
TESTS.append(("depth profile: switched off changes nothing, switched on only appends", test_depth_switched_off))
TESTS.append(("depth profile: circular ROI against ring areas and pixel geometry", test_depth_disk_roi))
TESTS.append(("depth profile: ROI cut by the image border and whole-image ROI", test_depth_image_border))
TESTS.append(("depth profile: triangle cut by two image borders, pixel by pixel", test_depth_corner_triangle))
TESTS.append(("depth profile: a hole in the ROI is outline too", test_depth_hole))
TESTS.append(("depth profile: maximum depth, bin widths and degenerate cases", test_depth_bins_and_degenerate))
TESTS.append(("depth profile: object with its centroid outside itself", test_depth_crescent_object))
TESTS.append(("depth profile: frames and folder batch", test_depth_frames_and_batch))
TESTS.append(("depth profile: 8-bit, micron with stage origin, ROI past the image border, by hand", test_depth_byte_micron_cut))
TESTS.append(("depth profile: tiny ROI by hand, two ROIs, two frames, no objects", test_depth_tiny_rois_two_frames))
TESTS.append(("depth profile: non-convex star with a 32-bit channel against a full search", test_depth_star_32bit))
TESTS.append(("depth profile: calibrated 16-bit and 8-bit bin means", test_depth_calibrated_bin_means))
TESTS.append(("depth profile: every setting survives a profile round trip", test_depth_profile_round_trip))

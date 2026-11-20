# Cell Quant Studio Universal - tests of the automatic tissue outline
# (cellquanttissue.py).
#
# Loaded by run_tests.py, which provides the helpers used here (make_image,
# nuclei_config, run_image, read_utf8, ...). Every image is synthetic and its
# correct outline is known by construction: the tissue is drawn with the same
# ImageJ ROI that is then filled into a reference mask, so "the right answer"
# is that mask, pixel for pixel.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij.process import Blitter, ImageStatistics, LUT, FloatProcessor
from java.awt import Color
from ij.measure import Measurements, Calibration
from ij.io import RoiEncoder, RoiDecoder
from java.io import FileOutputStream

TT_SIZE = 300


def tt_disk(cx, cy, radius):
    return OvalRoi(cx - radius, cy - radius, 2 * radius, 2 * radius)


def tt_plane(background, shapes, width=TT_SIZE, height=TT_SIZE, noise=0.0, seed=7):
    """16-bit plane of one background value; each (roi, value) of shapes is
    then filled in order, so a later shape can cut a hole into an earlier
    one. Gaussian noise of the given standard deviation is added last."""
    ip = ShortProcessor(width, height)
    ip.setValue(float(background))
    ip.fill()
    for shape in shapes:
        ip.setValue(float(shape[1]))
        ip.fill(shape[0])
    ip.resetRoi()
    if noise > 0:
        ImageProcessor.setRandomSeed(seed)
        ip.noise(noise)
    return ip


def tt_truth(shapes, width=TT_SIZE, height=TT_SIZE):
    """Reference mask: each (roi, is_tissue) of shapes filled in order with
    255 (tissue) or 0 (background)."""
    ip = ByteProcessor(width, height)
    for shape in shapes:
        ip.setValue(255.0 if shape[1] else 0.0)
        ip.fill(shape[0])
    ip.resetRoi()
    return ip


def tt_count(mask):
    mask.resetRoi()
    return mask.getHistogram()[255]


def tt_image(title, planes, pixel=None, unit="um", slices=1):
    """Image from planes in ImageJ stack order (channel fastest, then Z)."""
    stack = ImageStack(planes[0].getWidth(), planes[0].getHeight())
    for plane in planes:
        stack.addSlice("", plane)
    imp = ImagePlus(title, stack)
    channels = int(len(planes) / slices)
    imp.setDimensions(channels, slices, 1)
    if channels > 1:
        imp = CompositeImage(imp, CompositeImage.GRAYSCALE)
    elif slices > 1:
        imp.setOpenAsHyperStack(True)
    if pixel is not None:
        cal = imp.getCalibration()
        cal.pixelWidth = pixel
        cal.pixelHeight = pixel
        cal.setUnit(unit)
    return imp


def tt_options(E, values):
    """Global options with tissue settings given without the "tissue_" prefix."""
    options = E.default_global_options()
    for key in values:
        options["tissue_" + key] = values[key]
    return options


def tt_detect(E, imp, values, configs=None):
    return E.detect_tissue_outline(None, imp, configs or [], tt_options(E, values))


def tt_union(E, rois, width, height):
    mask = ByteProcessor(width, height)
    for roi in rois:
        mask.copyBits(E.roi_inside_mask(roi, width, height), 0, 0, Blitter.OR)
    return mask


def tt_mismatch(E, rois, truth):
    """Number of pixels where the detected ROIs and the reference differ."""
    mask = tt_union(E, rois, truth.getWidth(), truth.getHeight())
    mask.copyBits(truth, 0, 0, Blitter.DIFFERENCE)
    return tt_count(mask)


def tt_overlap(E, rois, truth):
    """Number of reference pixels that lie inside the detected ROIs."""
    mask = tt_union(E, rois, truth.getWidth(), truth.getHeight())
    mask.copyBits(truth, 0, 0, Blitter.AND)
    return tt_count(mask)


def tt_centroid(mask):
    ip = mask.duplicate()
    ip.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
    stats = ImageStatistics.getStatistics(ip, Measurements.CENTROID | Measurements.LIMIT, None)
    return (stats.xCentroid, stats.yCentroid)


def tt_inverse_normal(p):
    """x with Phi(x) = p (standard normal), by bisection on math.erf."""
    low = -10.0
    high = 10.0
    for step in range(200):
        middle = (low + high) / 2.0
        if 0.5 * (1.0 + math.erf(middle / math.sqrt(2.0))) < p:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def tt_read_table(path):
    """(headers, rows) of a CSV written by write_csv; fields hold no '","'."""
    lines = [line for line in read_utf8(path).split("\n") if line.strip() != ""]

    def split(line):
        return [part.replace('""', '"') for part in line.strip()[1:-1].split('","')]

    headers = split(lines[0])
    return (headers, [dict(zip(headers, split(line))) for line in lines[1:]])


# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

def test_tissue_disk(E, t, tmp):
    """A disk of known radius on a dark background."""
    # (a) Two-level image, no smoothing: any threshold between the two
    # levels gives exactly the drawn disk, so the outline must equal the
    # reference mask pixel for pixel.
    disk = tt_disk(150, 140, 80)
    truth = tt_truth([(disk, True)])
    pixels = tt_count(truth)
    t.close(pixels, math.pi * 80 * 80, "construction: the drawn disk has about pi r^2 pixels", 0.01)
    imp = tt_image("disk.tif", [tt_plane(100, [(disk, 1000)])])
    rois, info = tt_detect(E, imp, {"sigma": "0"})
    t.equal(len(rois), 1, "one tissue region")
    t.equal(E.to_text(rois[0].getName()), "Tissue_001", "the region is named Tissue_001")
    t.equal(info["area_pixels"], pixels, "detected area equals the drawn disk")
    t.equal(tt_mismatch(E, rois, truth), 0, "the outline matches the disk pixel for pixel")
    expected = tt_centroid(truth)
    t.close(info["regions"][0]["centroid_x"], expected[0], "centroid x is the disk centroid")
    t.close(info["regions"][0]["centroid_y"], expected[1], "centroid y is the disk centroid")
    t.check(abs(expected[0] - 150.0) < 0.5 and abs(expected[1] - 140.0) < 0.5, "construction: the disk is centred at (150, 140)")
    t.check(100 < info["threshold"] <= 1000, "the threshold lies between background and tissue")
    t.equal(info["scale"], "log", "the default threshold scale is logarithmic")
    t.equal(info["length_unit"], "pixel", "an uncalibrated image uses pixels")
    t.equal((info["regions_found"], info["regions_kept"]), (1, 1), "region counts")
    t.equal((info["touches_border"], info["border_sides"]), (False, ""), "the disk does not touch the border")
    t.close(info["area_percent"], 100.0 * pixels / (TT_SIZE * TT_SIZE), "area as a percentage of the image")
    t.check(E.to_text(info["description"]).startswith("Automatic tissue outline: 1 region(s)"), "the description says what was found")
    t.equal(imp.getRoi(), None, "the source image is left without a selection")
    imp.close()

    # (b) Noise (SD 40) and smoothing with sigma 1.5 um on 0.5 um pixels,
    # that is exactly 3 pixels. The blurred edge has the profile
    # I(d) = 100 + 900 * Phi(-d / 3) at outward distance d, so a threshold
    # L puts the outline at d = -3 * InvPhi((L - 100) / 900). The detected
    # equivalent radius must be the drawn one plus that d, within 1 pixel
    # (noise, edge curvature and pixel steps).
    # Which side d is on follows from the histogram Triangle sees: it works
    # from the tallest peak. Smoothing by 3 pixels leaves noise of about
    # 40 / (2 sqrt(pi) 3) = 3.8 in both regions. On the linear scale both
    # peaks are equally wide and the background (78 percent of the pixels)
    # is taller, so the threshold sits at its foot and d > 0 (outside). On
    # the log scale the background peak is 10 times wider (3.8 / 100
    # against 3.8 / 1000) and the tissue peak, with 22 percent of the
    # pixels, is taller: the threshold sits at the foot of the tissue peak
    # and d < 0 (inside). Either way the threshold lies inside the blurred
    # edge, so |d| is at most about 2.5 sigma.
    imp = tt_image("noisy.tif", [tt_plane(100, [(disk, 1000)], noise=40)], pixel=0.5)
    expected = tt_centroid(truth)
    for log_scale, side in [(True, -1), (False, 1)]:
        label = "noisy disk, " + ("log" if log_scale else "linear") + " scale: "
        rois, info = tt_detect(E, imp, {"sigma": "1.5", "log_scale": log_scale})
        t.equal(info["length_unit"], "micrometer", label + "a micrometer calibration is recognised")
        t.close(info["sigma_pixels"], 3.0, label + "1.5 um at 0.5 um per pixel is 3 pixels")
        t.equal(len(rois), 1, label + "one region")
        if len(rois) != 1:
            continue
        predicted = -3.0 * tt_inverse_normal((info["threshold"] - 100.0) / 900.0)
        shift = math.sqrt(info["area_pixels"] / math.pi) - math.sqrt(pixels / math.pi)
        t.check(abs(shift - predicted) <= 1.0, label + "outline offset matches the blurred-edge prediction (" +
                ("%.2f" % shift) + " vs " + ("%.2f" % predicted) + " pixels)")
        t.check(side * predicted > 0 and abs(predicted) <= 7.5, label + "offset on the expected side and within 2.5 sigma")
        inside = tt_overlap(E, rois, truth)
        if side > 0:
            t.check(inside >= 0.995 * pixels, label + "the outline encloses the whole disk")
        else:
            t.check(inside >= 0.995 * info["area_pixels"], label + "the outline lies within the disk")
        t.check(abs(info["regions"][0]["centroid_x"] - expected[0]) < 0.5 and
                abs(info["regions"][0]["centroid_y"] - expected[1]) < 0.5, label + "centroid within half a pixel")
    # A threshold half way between the levels (Otsu on the linear scale, for
    # well separated levels) puts the outline on the true edge: d = 0.
    rois, info = tt_detect(E, imp, {"sigma": "1.5", "method": "Otsu", "log_scale": False})
    t.close(info["area_pixels"], pixels, "a mid-level threshold gives the disk area within 2 percent", 0.02)
    t.close(info["area_pixels"], math.pi * 80 * 80, "a mid-level threshold gives pi r^2 within 3 percent", 0.03)
    t.equal(info["scale"], "linear", "the linear scale is reported")
    imp.close()


def test_tissue_two_blobs(E, t, tmp):
    """Two separate regions: all kept, the largest only, or a minimum area."""
    big = tt_disk(90, 150, 60)
    small = tt_disk(230, 150, 35)
    big_truth = tt_truth([(big, True)])
    small_truth = tt_truth([(small, True)])
    big_pixels = tt_count(big_truth)
    small_pixels = tt_count(small_truth)
    image_pixels = float(TT_SIZE * TT_SIZE)
    imp = tt_image("pair.tif", [tt_plane(100, [(big, 1000), (small, 600)])])

    rois, info = tt_detect(E, imp, {"sigma": "0"})
    t.equal(len(rois), 2, "both regions are kept by default")
    t.equal([E.to_text(roi.getName()) for roi in rois], ["Tissue_001", "Tissue_002"], "named by decreasing area")
    t.equal(tt_mismatch(E, [rois[0]], big_truth), 0, "Tissue_001 is the larger disk")
    t.equal(tt_mismatch(E, [rois[1]], small_truth), 0, "Tissue_002 is the smaller disk")
    t.equal((info["regions_found"], info["regions_kept"]), (2, 2), "two found, two kept")
    t.equal(info["area_pixels"], big_pixels + small_pixels, "total area is the sum of both disks")
    t.close(info["rejected_area_percent"], 0.0, "nothing rejected")

    rois, info = tt_detect(E, imp, {"sigma": "0", "largest_only": True})
    t.equal(len(rois), 1, "largest only: one region")
    t.equal(tt_mismatch(E, rois, big_truth), 0, "largest only: the larger disk")
    t.equal((info["regions_found"], info["regions_kept"]), (2, 1), "largest only: two found, one kept")
    t.close(info["rejected_area_percent"], 100.0 * small_pixels / image_pixels, "the smaller disk is reported as rejected area")
    t.check("largest of 2" in E.to_text(info["description"]), "the description mentions the other region")

    # The small disk is about 4.3 percent of the image, the large one 12.6.
    t.check(small_pixels / image_pixels < 0.06 < big_pixels / image_pixels, "construction: 6 percent lies between the two areas")
    rois, info = tt_detect(E, imp, {"sigma": "0", "min_area_percent": "6"})
    t.equal(len(rois), 1, "minimum area 6 percent: one region")
    t.equal(info["regions_found"], 1, "regions below the minimum are not counted as found")
    t.equal(tt_mismatch(E, rois, big_truth), 0, "minimum area: the larger disk is kept")
    t.close(info["rejected_area_percent"], 100.0 * small_pixels / image_pixels, "minimum area: rejected area")
    imp.close()


def test_tissue_border(E, t, tmp):
    """Tissue cut by the image border."""
    # (a) A ring cut by the left border: its inner disk is open to the
    # border, so it is background, not a hole, and fill holes must leave it.
    outer = tt_disk(10, 150, 90)
    inner = tt_disk(10, 150, 40)
    truth = tt_truth([(outer, True), (inner, False)])
    pixels = tt_count(truth)
    truth.setRoi(0, 0, 1, TT_SIZE)
    contact = truth.getHistogram()[255]
    truth.resetRoi()
    imp = tt_image("cut.tif", [tt_plane(100, [(outer, 1000), (inner, 100)])])
    rois, info = tt_detect(E, imp, {"sigma": "0", "fill_holes": True})
    t.equal(len(rois), 1, "cut ring: one region")
    t.equal(tt_mismatch(E, rois, truth), 0, "cut ring: outline equals the part inside the image")
    t.equal(info["filled_hole_pixels"], 0, "a bay open to the border is not filled")
    t.equal((info["touches_border"], info["border_sides"]), (True, "left"), "touches the left border")
    t.equal(info["regions"][0]["border_contact_pixels"], contact, "border contact counts the tissue pixels in column 0")
    t.check("touches the image border (left)" in E.to_text(info["description"]), "the description names the side")
    t.equal(len(E.validate_parent_rois(rois, TT_SIZE, TT_SIZE)), 1, "the outline is a valid parent ROI")
    t.equal(info["area_pixels"], pixels, "cut ring: area")
    imp.close()

    # (b) A disk in the lower right corner touches two sides. A corner
    # pixel belongs to both border lines and is counted once.
    corner = tt_disk(280, 285, 60)
    truth = tt_truth([(corner, True)])
    edge = TT_SIZE - 1
    truth.setRoi(edge, 0, 1, TT_SIZE)
    right = truth.getHistogram()[255]
    truth.setRoi(0, edge, TT_SIZE, 1)
    bottom = truth.getHistogram()[255]
    truth.resetRoi()
    both = 1 if truth.get(edge, edge) == 255 else 0
    imp = tt_image("corner.tif", [tt_plane(100, [(corner, 1000)])])
    rois, info = tt_detect(E, imp, {"sigma": "0"})
    t.equal(tt_mismatch(E, rois, truth), 0, "corner disk: outline equals the part inside the image")
    t.equal(info["border_sides"], "right+bottom", "corner disk touches right and bottom")
    t.equal(info["regions"][0]["border_contact_pixels"], right + bottom - both, "corner contact pixels")
    imp.close()


def test_tissue_empty(E, t, tmp):
    """Images without tissue: no outline and a clear reason, no error."""
    for value in [0, 500]:
        imp = tt_image("blank.tif", [tt_plane(value, [])])
        rois, info = tt_detect(E, imp, {})
        t.equal(rois, [], "uniform image (" + str(value) + "): no outline")
        t.equal(info["regions_kept"], 0, "uniform image: no region kept")
        t.check(math.isnan(info["threshold"]), "uniform image: no threshold")
        description = E.to_text(info["description"])
        t.check(description.startswith("Automatic tissue outline: no tissue found") and "uniform" in description,
                "uniform image: the description says why (" + description + ")")
        rois, note = E.tissue_batch_rois({"imp": imp, "path": "", "configs": [], "global_options": tt_options(E, {}), "log": None})
        t.equal(rois, None, "batch policy: an image without tissue is skipped")
        t.check("no tissue found" in E.to_text(note), "batch policy: the reason goes to the report")
        imp.close()
    # A speck of 3 x 3 pixels gives a threshold but no region reaches 1
    # percent of the image (900 pixels).
    imp = tt_image("speck.tif", [tt_plane(100, [(Roi(100, 100, 3, 3), 1000)])])
    rois, info = tt_detect(E, imp, {"sigma": "0"})
    t.equal(rois, [], "a speck is not tissue")
    t.check(not math.isnan(info["threshold"]), "the speck still gives a threshold")
    t.check("no connected region reaches" in E.to_text(info["description"]), "the reason names the minimum area")
    t.equal(info["regions_found"], 0, "no region found")
    imp.close()


def test_tissue_holes(E, t, tmp):
    """A blob with holes, with fill holes on and off."""
    outer = tt_disk(150, 150, 100)
    hole_a = tt_disk(120, 130, 20)
    hole_b = tt_disk(190, 170, 15)
    full = tt_truth([(outer, True)])
    ring = tt_truth([(outer, True), (hole_a, False), (hole_b, False)])
    hole_pixels = tt_count(full) - tt_count(ring)
    shapes = [(outer, 1000), (hole_a, 100), (hole_b, 100)]
    imp = tt_image("holes.tif", [tt_plane(100, shapes)])
    rois, info = tt_detect(E, imp, {"sigma": "0", "fill_holes": True})
    t.equal(len(rois), 1, "fill holes: one region")
    t.equal(tt_mismatch(E, rois, full), 0, "fill holes: the outline is the outer disk")
    t.equal(info["filled_hole_pixels"], hole_pixels, "fill holes: the filled pixels are the two holes")
    rois, info = tt_detect(E, imp, {"sigma": "0", "fill_holes": False})
    t.equal(len(rois), 1, "holes kept: one region")
    t.equal(tt_mismatch(E, rois, ring), 0, "holes kept: the region excludes both holes exactly")
    t.equal(info["area_pixels"], tt_count(ring), "holes kept: area excludes the holes")
    t.equal(info["filled_hole_pixels"], 0, "holes kept: nothing filled")
    imp.close()

    # Textured tissue (noise SD 250 on a level of 1000) smoothed by 2
    # pixels: the holes must be inside the outline when filled and outside
    # it when kept.
    imp = tt_image("textured.tif", [tt_plane(100, shapes, noise=250)])
    filled, info_filled = tt_detect(E, imp, {"sigma": "2", "fill_holes": True})
    kept, info_kept = tt_detect(E, imp, {"sigma": "2", "fill_holes": False})
    t.equal((len(filled), len(kept)), (1, 1), "textured: one region either way")
    t.check(filled[0].contains(120, 130) and filled[0].contains(190, 170), "textured, filled: the hole centres are inside")
    t.check(not kept[0].contains(120, 130) and not kept[0].contains(190, 170), "textured, kept: the hole centres are outside")
    difference = info_filled["area_pixels"] - info_kept["area_pixels"]
    t.check(0.5 * hole_pixels < difference < 1.5 * hole_pixels, "textured: filling adds about the hole area")
    imp.close()


def test_tissue_prefs_independent(E, t, tmp):
    """The same outline whatever the black-background and edge-padding
    options are, and the options are left as they were."""
    shapes = [(tt_disk(150, 150, 100), 1000), (tt_disk(120, 130, 20), 100), (Roi(250, 20, 40, 40), 1000)]
    imp = tt_image("prefs.tif", [tt_plane(100, shapes, noise=60)], pixel=0.5)
    values = {"sigma": "1", "close_radius": "1", "shrink": "1"}
    saved = (Prefs.blackBackground, Prefs.padEdges)
    results = []
    try:
        for setting in [(True, True), (False, False), (False, True)]:
            Prefs.blackBackground = setting[0]
            Prefs.padEdges = setting[1]
            rois, info = tt_detect(E, imp, values)
            t.equal((Prefs.blackBackground, Prefs.padEdges), setting, "detection leaves the options unchanged")
            results.append((rois, info))
    finally:
        Prefs.blackBackground = saved[0]
        Prefs.padEdges = saved[1]
    reference = tt_union(E, results[0][0], TT_SIZE, TT_SIZE)
    for index in range(1, len(results)):
        t.equal(len(results[index][0]), len(results[0][0]), "same number of regions")
        t.equal(results[index][1]["threshold"], results[0][1]["threshold"], "same threshold")
        t.equal(tt_mismatch(E, results[index][0], reference), 0, "identical outline")
    t.equal(len(results[0][0]), 2, "the disk and the square are found")
    imp.close()


def test_tissue_channels(E, t, tmp):
    """Which channel makes the tissue image, and its Z processing."""
    disk_a = tt_disk(90, 150, 60)
    disk_b = tt_disk(230, 150, 35)
    truth_a = tt_truth([(disk_a, True)])
    truth_b = tt_truth([(disk_b, True)])
    # C1 is dim and shows disk A; C2 is twenty times brighter and shows
    # disk B. Averaged after rescaling each to 0-1, both disks have the
    # same value 0.5, so channel 0 must find both.
    imp = tt_image("channels.tif", [tt_plane(100, [(disk_a, 300)]), tt_plane(200, [(disk_b, 6000)])])
    rois, info = tt_detect(E, imp, {"sigma": "0", "channel": "1"})
    t.equal((len(rois), info["channels"]), (1, "C1"), "channel 1: one region from C1")
    t.equal(tt_mismatch(E, rois, truth_a), 0, "channel 1 finds disk A")
    rois, info = tt_detect(E, imp, {"sigma": "0", "channel": "2"})
    t.equal(tt_mismatch(E, rois, truth_b), 0, "channel 2 finds disk B")
    rois, info = tt_detect(E, imp, {"sigma": "0", "channel": "0"})
    t.equal(info["channels"], "C1+C2", "channel 0 without settings uses every channel")
    t.equal(len(rois), 2, "channel 0 finds both disks")
    t.equal(tt_mismatch(E, rois, tt_truth([(disk_a, True), (disk_b, True)])), 0, "channel 0: both disks exactly")
    only_c2 = [E.normalize_channel_config(2, {"enabled": True, "name": "B", "analysis_mode": "Positive area"})]
    rois, info = tt_detect(E, imp, {"sigma": "0", "channel": "0"}, only_c2)
    t.equal(info["channels"], "C2", "channel 0 uses the analysed channels")
    failed = False
    try:
        tt_detect(E, imp, {"channel": "3"})
    except ValueError as exc:
        failed = "does not exist" in E.to_text(exc)
    t.check(failed, "a channel the image does not have is refused clearly")
    imp.close()

    # One channel, two Z slices: Z1 shows disk A, Z2 disk B. Not analysed,
    # the channel is Max-projected (both disks); analysed with Specific Z 2,
    # only disk B is seen.
    imp = tt_image("zstack.tif", [tt_plane(100, [(disk_a, 1000)]), tt_plane(100, [(disk_b, 1000)])], slices=2)
    rois, info = tt_detect(E, imp, {"sigma": "0", "channel": "1"})
    t.equal(len(rois), 2, "Max Z by default: both disks")
    z2 = [E.normalize_channel_config(1, {"enabled": True, "name": "N", "analysis_mode": "Positive area",
                                         "projection": "Specific Z", "z_slice": 2})]
    rois, info = tt_detect(E, imp, {"sigma": "0", "channel": "1"}, z2)
    t.equal(len(rois), 1, "the channel's own Z setting is used")
    t.equal(tt_mismatch(E, rois, truth_b), 0, "Specific Z 2: disk B")
    imp.close()


def test_tissue_settings(E, t, tmp):
    """Defaults, profile round trip, old profiles, bad values, lengths."""
    defaults = E.tissue_settings({})
    t.equal(defaults, {
        "channel": 0, "sigma": 5.0, "method": "Triangle", "log_scale": True, "min_area_percent": 1.0,
        "largest_only": False, "fill_holes": True, "close_radius": 0.0, "shrink": 0.0,
        "save_qc": True, "save_rois": True
    }, "documented defaults")
    t.equal(E.tissue_settings(E.default_global_options()), defaults, "default global options give the defaults")
    keys = [spec[0] for spec in E.GLOBAL_SETTING_SPECS]
    form_keys = []
    for form in E.SETTINGS_FORMS:
        if form["title"] == "Tissue outline":
            form_keys = [field[0] for field in form["fields"]]
    for spec in E.TISSUE_SETTING_SPECS:
        t.check(spec[0] in keys, spec[0] + " is saved in profiles")
        t.check(spec[1] in form_keys, spec[1] + " is shown in the settings form")
    t.check(E.TISSUE_POLICY in E.BATCH_ROI_POLICIES, "the batch policy is registered")
    t.equal(E.ENGINE_MODULE_VERSIONS.get("cellquanttissue.py"), E.APP_VERSION, "module version")

    options = tt_options(E, {
        "channel": "2", "sigma": "3.5", "method": "Li", "log_scale": False, "min_area_percent": "2.5",
        "largest_only": True, "fill_holes": False, "close_radius": "1.5", "shrink": "2",
        "save_qc": False, "save_rois": False
    })
    path = os.path.join(tmp, "tissue_profile.properties")
    E.write_settings_properties(path, options, [nuclei_config(E)])
    loaded = E.load_profile_settings(path)[0]
    t.equal(E.tissue_settings(loaded), E.tissue_settings(options), "settings survive a profile round trip")
    t.equal(E.tissue_settings(loaded)["method"], "Li", "round trip keeps the method")

    # A profile from before this module: the same file without tissue keys.
    props = E.read_properties_file(path)
    for spec in E.TISSUE_SETTING_SPECS:
        props.remove(spec[0])
    old_path = os.path.join(tmp, "tissue_old_profile.properties")
    stream = FileOutputStream(old_path)
    try:
        props.store(stream, "old profile")
    finally:
        stream.close()
    t.equal(E.tissue_settings(E.load_profile_settings(old_path)[0]), defaults, "an older profile loads with the defaults")

    for bad, word in [({"tissue_method": "Nonsense"}, "unknown threshold method"),
                      ({"tissue_min_area_percent": "150"}, "at most 100"),
                      ({"tissue_sigma": "-1"}, "at least"),
                      ({"tissue_channel": "x"}, "integer")]:
        raised = ""
        try:
            E.tissue_settings(bad)
        except ValueError as exc:
            raised = E.to_text(exc)
        t.check(word in raised, "bad setting refused with a clear message: " + repr(bad))

    for unit, size, expected in [("micron", 0.25, ("micrometer", 0.25, 0.25)),
                                 (MICRO + "m", 0.5, ("micrometer", 0.5, 0.5)),
                                 ("nm", 100.0, ("micrometer", 0.1, 0.1)),
                                 ("pixel", 1.0, ("pixel", 1.0, 1.0)),
                                 ("", 2.0, ("pixel", 1.0, 1.0))]:
        imp = tt_image("cal.tif", [tt_plane(0, [], width=8, height=8)], pixel=size, unit=unit)
        scale = E.tissue_length_scale(imp)
        t.equal(scale[0], expected[0], "length unit for '" + E.to_text(unit) + "'")
        t.close(scale[1], expected[1], "micrometers per pixel for '" + E.to_text(unit) + "'")
        imp.close()

    # The calibration override of the settings applies: 4 um at 2 um per
    # pixel is 2 pixels.
    imp = tt_image("override.tif", [tt_plane(100, [(tt_disk(150, 150, 80), 1000)])])
    options = tt_options(E, {"sigma": "4"})
    options.update({"override_calibration": True, "pixel_width": "2", "pixel_height": "2", "spatial_unit": "um"})
    rois, info = E.detect_tissue_outline(None, imp, [], options)
    t.close(info["sigma_pixels"], 2.0, "the calibration override sets the smoothing in pixels")
    imp.close()


def tt_plain_copy(E, roi):
    """A copy of an ROI without the detection record. (Roi.clone() shares
    its properties with the original, so the copy goes through a ROI file.)"""
    copy = RoiDecoder(RoiEncoder.saveAsByteArray(roi), "copy").getRoi()
    keys = [pair[0] for pair in E.TISSUE_RECORD] + ["Width", "Height", "Bounds"]
    for key in keys:
        copy.setProperty(E.TISSUE_TAG + key, "")
    return copy


def test_tissue_in_analysis(E, t, tmp):
    """A detected outline used as parent ROI: the run records how it was
    made, saves it, and measures exactly like the same ROI drawn by hand."""
    disk = tt_disk(150, 140, 80)
    pixels = tt_count(tt_truth([(disk, True)]))
    imp = tt_image("tissue_run.tif", [tt_plane(100, [(disk, 1000)])], pixel=0.5)
    rois, info = tt_detect(E, imp, {"sigma": "0"})
    configs = [nuclei_config(E, {"manual_low": 500})]
    override = {"tissue_sigma": "0"}
    tagged = run_image(E, imp, configs, rois, os.path.join(tmp, "tissue_run"), options_override=override)
    t.check("tissue_outline.csv" in tagged["extra_table_names"], "the run writes tissue_outline.csv")
    table = tagged["extra_tables"].get("tissue_outline.csv", {"headers": [], "rows": []})
    t.equal(table["headers"], E.TISSUE_TABLE_HEADERS, "tissue_outline.csv columns")
    t.equal(len(table["rows"]), 1, "one row per detected parent ROI")
    headers, rows = tt_read_table(os.path.join(tagged["run_dir"], "tissue_outline.csv"))
    t.equal(headers, E.TISSUE_TABLE_HEADERS, "the written file has the same columns")
    if len(rows) == 1:
        row = rows[0]
        t.equal(row["Outline_Status"], "as detected", "an unchanged outline is reported as detected")
        t.equal(int(row["ROI_Area_Pixels"]), pixels, "ROI area in pixels")
        t.close(float(row["ROI_Area"]), pixels * 0.25, "ROI area in square micrometers")
        t.equal(row["Detected_Area_Pixels"], str(pixels), "the detection record travels with the ROI")
        t.equal((row["Threshold_Method"], row["Threshold_Scale"]), ("Triangle", "log"), "method and scale recorded")
        t.equal(row["Touches_Image_Border"], "false", "border flag")
        t.equal(row["Parent_ROI_Name"], "Tissue_001", "parent ROI name")
        t.equal(row["Run_ID"], tagged["run_id"], "Run_ID")
        t.equal(row["Application_Version"], E.APP_VERSION, "Application_Version")
    t.close(tagged["summary_rows"][0]["ROI_Area"], pixels * 0.25, "roi_summary ROI_Area is the disk area")
    saved = os.path.join(tagged["run_dir"], "ROIs", "Tissue_Outline_ROIs.zip")
    t.check(os.path.isfile(saved), "the outline is saved inside the run folder")
    t.check(os.path.isfile(os.path.join(tagged["run_dir"], "Overlays", "tissue_run_Tissue_Outline_QC.tif")), "the check image is saved")
    reloaded = E.load_rois_from_file(saved) if os.path.isfile(saved) else []
    t.equal([E.to_text(roi.getName()) for roi in reloaded], ["Tissue_001"], "the saved outline reloads with its name")
    if len(reloaded) == 1:
        t.equal(E.tissue_tag(reloaded[0], "Method"), "Triangle", "the saved outline keeps its detection record")

        # Moved by hand after detection: reported as edited.
        moved = reloaded[0]
        bounds = moved.getBounds()
        moved.setLocation(bounds.x + 5, bounds.y)
        edited = run_image(E, imp, configs, [moved], os.path.join(tmp, "tissue_edit"), options_override=override)
        status = [r["Outline_Status"] for r in edited["extra_tables"]["tissue_outline.csv"]["rows"]]
        t.equal(status, ["edited or reused"], "a moved outline is reported as edited")

    # The same outline without the detection record measures identically
    # and writes no tissue files.
    plain = tt_plain_copy(E, rois[0])
    t.equal(E.tissue_tag(plain, "Method"), "", "construction: the copy has no detection record")
    t.equal(E.tissue_tag(rois[0], "Method"), "Triangle", "construction: the original keeps its record")
    untagged = run_image(E, imp, configs, [plain], os.path.join(tmp, "tissue_plain"), options_override=override)
    t.equal(untagged["summary_headers"], tagged["summary_headers"], "the outline adds no roi_summary column")
    for key in tagged["summary_headers"]:
        if key in ["Run_ID"]:
            continue
        t.check(E.to_text(tagged["summary_rows"][0].get(key)) == E.to_text(untagged["summary_rows"][0].get(key)),
                "detected and plain outline give the same " + key)
    t.equal(untagged["extra_table_names"].count("tissue_outline.csv"), 0, "a plain ROI writes no tissue table")
    t.check(not os.path.exists(os.path.join(untagged["run_dir"], "ROIs", "Tissue_Outline_ROIs.zip")), "a plain ROI saves no tissue ROI file")

    # Switched off, the check image and ROI file are not written.
    quiet = run_image(E, imp, configs, rois, os.path.join(tmp, "tissue_quiet"),
                      options_override={"tissue_save_qc": False, "tissue_save_rois": False})
    t.check("tissue_outline.csv" in quiet["extra_table_names"], "the record table is still written")
    t.check(not os.path.exists(os.path.join(quiet["run_dir"], "ROIs", "Tissue_Outline_ROIs.zip")), "save_rois off: no ROI file")
    t.check(not os.path.exists(os.path.join(quiet["run_dir"], "Overlays", "tissue_run_Tissue_Outline_QC.tif")), "save_qc off: no check image")
    imp.close()


def test_tissue_unused_changes_nothing(E, t, tmp):
    """With ordinary parent ROIs every standard table is byte for byte what
    it is without this module, whatever the tissue settings are."""
    tissue_extensions = [extension for extension in E.RUN_EXTENSIONS if extension.name == "tissue_outline"]
    t.equal(len(tissue_extensions), 1, "the tissue extension is registered")
    configs = [nuclei_config(E), marker_config(E)]

    def run(folder, override):
        imp = make_image("plain.tif")
        rois = [Roi(0, 0, WIDTH, HEIGHT), triangle_roi()]
        rois[0].setName("whole")
        options = E.default_global_options()
        options.update(override)
        run_dir = os.path.join(tmp, folder, "same_run")
        os.makedirs(run_dir)
        try:
            return E.AnalysisEngine().analyze_image(
                imp, "", configs, rois, E.engine_options_from_global(options), run_dir, [1],
                settings=(options, configs)
            )
        finally:
            imp.close()

    saved = list(E.RUN_EXTENSIONS)
    try:
        E.RUN_EXTENSIONS[:] = [extension for extension in saved if extension.name != "tissue_outline"]
        without = run("tissue_without", {})
    finally:
        E.RUN_EXTENSIONS[:] = saved
    with_default = run("tissue_with", {})
    with_changed = run("tissue_changed", {"tissue_sigma": "2", "tissue_method": "Li", "tissue_largest_only": True})
    for result in [with_default, with_changed]:
        t.equal(result["summary_headers"], without["summary_headers"], "roi_summary columns unchanged")
        t.equal(result["object_headers"], without["object_headers"], "objects columns unchanged")
        t.equal(result["extra_table_names"], without["extra_table_names"], "no extra table")
        for name in ["roi_summary.csv", "objects.csv", "thresholds.csv"]:
            t.check(read_utf8(os.path.join(result["run_dir"], name)) == read_utf8(os.path.join(without["run_dir"], name)),
                    name + " is byte for byte unchanged")
        t.check(not os.path.exists(os.path.join(result["run_dir"], "tissue_outline.csv")), "no tissue_outline.csv")
        t.check(not os.path.exists(os.path.join(result["run_dir"], "ROIs", "Tissue_Outline_ROIs.zip")), "no tissue ROI file")


def test_tissue_batch_policy(E, t, tmp):
    """The batch policy on a small folder: outlines found, an image without
    tissue skipped with its reason, files written only in the output, and
    the saved outlines reusable as ROI files per image."""
    folder = os.path.join(tmp, "tissue_batch_in")
    os.makedirs(folder)
    disk = tt_disk(150, 140, 80)
    big = tt_disk(90, 150, 60)
    small = tt_disk(230, 150, 35)
    planes = {
        "disk.tif": tt_plane(100, [(disk, 1000)]),
        "pair.tif": tt_plane(100, [(big, 1000), (small, 600)]),
        "blank.tif": tt_plane(100, [])
    }
    for name in planes:
        imp = tt_image(name, [planes[name]], pixel=0.5)
        FileSaver(imp).saveAsTiff(os.path.join(folder, name))
        imp.close()
    before = sorted(os.listdir(folder))
    disk_area = tt_count(tt_truth([(disk, True)])) * 0.25
    pair_areas = [tt_count(tt_truth([(big, True)])) * 0.25, tt_count(tt_truth([(small, True)])) * 0.25]

    configs = [nuclei_config(E, {"manual_low": 500})]
    options = E.default_global_options()
    for key in ["save_prepared", "save_masks", "save_overlays", "save_object_rois"]:
        options[key] = False
    options["tissue_sigma"] = "0"
    out = os.path.join(tmp, "tissue_batch_out")
    result = E.run_batch(folder, out, configs, configs, options, {"extensions": "tif", "roi_policy": E.TISSUE_POLICY})
    t.equal((result["analyzed"], result["skipped"], result["failed"]), (2, 1, 0), "two analyzed, the blank one skipped")
    t.equal(sorted(os.listdir(folder)), before, "nothing is written beside the source images")
    batch_dir = result["batch_dir"]

    headers, report = tt_read_table(os.path.join(batch_dir, "batch_report.csv"))
    by_file = dict([(row["Image_File"], row) for row in report])
    t.equal(by_file["blank.tif"]["Status"], "SKIPPED", "the blank image is skipped")
    t.check("no tissue found" in by_file["blank.tif"]["Message"], "the report says why it was skipped")
    t.check(by_file["disk.tif"]["ROI_Source"].startswith("Automatic tissue outline: 1 region(s)"), "the report says what was found")
    t.check(by_file["pair.tif"]["ROI_Source"].startswith("Automatic tissue outline: 2 region(s)"), "both regions of the pair")

    headers, summary = tt_read_table(os.path.join(batch_dir, "batch_roi_summary.csv"))
    areas = {}
    for row in summary:
        areas.setdefault(row["Source_Image"], []).append(float(row["ROI_Area"]))
    t.equal(sorted(areas.keys()), ["disk.tif", "pair.tif"], "summary rows for the two analyzed images")
    t.check(len(areas.get("disk.tif", [])) == 1 and abs(areas["disk.tif"][0] - disk_area) < 1e-6, "disk ROI area")
    t.check(sorted(areas.get("pair.tif", [])) == sorted(pair_areas), "pair ROI areas")

    combined = os.path.join(batch_dir, "batch_tissue_outline.csv")
    t.check(os.path.isfile(combined), "the batch combines the outline records")
    if os.path.isfile(combined):
        headers, rows = tt_read_table(combined)
        t.equal(headers, E.TISSUE_TABLE_HEADERS + ["Batch_ID", "Image_Relative_Path"], "combined record columns")
        t.equal(len(rows), 3, "one record per outline")
        t.equal(set([row["Outline_Status"] for row in rows]), set(["as detected"]), "all as detected")
    collection = os.path.join(batch_dir, E.TISSUE_COLLECTION_DIR)
    for stem in ["disk", "pair"]:
        t.check(os.path.isfile(os.path.join(collection, stem + "_ROIs.zip")), stem + ": outline collected for the batch")
        t.check(os.path.isfile(os.path.join(collection, stem + "_Tissue_Outline_QC.tif")), stem + ": check image collected")
    t.check(not os.path.exists(os.path.join(collection, "blank_ROIs.zip")), "no outline file for the skipped image")

    # Reuse the collected outlines with the policy "ROI file per image".
    again = E.run_batch(folder, out, configs, configs, options,
                        {"extensions": "tif", "roi_policy": E.BATCH_ROI_FILE_SKIP, "roi_dir": collection})
    t.equal((again["analyzed"], again["skipped"]), (2, 1), "reused outlines: same images analyzed")
    headers, summary = tt_read_table(os.path.join(again["batch_dir"], "batch_roi_summary.csv"))
    reused = {}
    for row in summary:
        reused.setdefault(row["Source_Image"], []).append(float(row["ROI_Area"]))
    t.equal(dict([(k, sorted(v)) for k, v in reused.items()]), dict([(k, sorted(v)) for k, v in areas.items()]),
            "reused outlines measure the same areas")
    combined = os.path.join(again["batch_dir"], "batch_tissue_outline.csv")
    if os.path.isfile(combined):
        headers, rows = tt_read_table(combined)
        t.equal(set([row["Outline_Status"] for row in rows]), set(["as detected"]), "reused outlines keep their record")
    else:
        t.check(False, "reused outlines are still recorded")


def test_tissue_check_image(E, t, tmp):
    """The check image: RGB, at most 1024 pixels on its longer side, the
    outline drawn in yellow, the input left unchanged."""
    plane = tt_plane(100, [(tt_disk(400, 250, 200), 1000)], width=2000, height=1000)
    # A colour table on the input (as Bio-Formats gives a CZI channel)
    # must not colour the check image: it is shown in gray.
    plane.setColorModel(LUT.createLutFromColor(Color.red))
    imp = tt_image("wide.tif", [plane])
    rois, info = tt_detect(E, imp, {"sigma": "0"}, None)
    before = ImageStatistics.getStatistics(plane, Measurements.MEAN, None).mean
    path = E.save_tissue_qc(os.path.join(tmp, "tissue_qc", "wide_QC.tif"), [plane], rois, "wide")
    t.check(os.path.isfile(path), "the check image is written (folder created)")
    picture = IJ.openImage(path)
    t.equal((picture.getBitDepth(), picture.getWidth(), picture.getHeight()), (24, 1024, 512), "RGB, reduced to 1024 x 512")
    pixels = picture.getProcessor().getPixels()
    yellow = 0
    colored = 0
    for value in pixels:
        rgb = value & 0xffffff
        if rgb == 0xffff00:
            yellow += 1
        elif (rgb >> 16) != (rgb & 0xff):
            colored += 1
    t.check(yellow > 500, "the outline is drawn in yellow (" + str(yellow) + " pixels)")
    t.equal(colored, 0, "everything but the outline is gray")
    t.equal(ImageStatistics.getStatistics(plane, Measurements.MEAN, None).mean, before, "the tissue image is not changed")
    picture.close()
    rois, info = E.detect_tissue_outline(None, imp, [], tt_options(E, {"sigma": "2"}), 1, True)
    kept = info.get("tissue_image")
    t.check(kept is not None and kept.getBitDepth() == 32, "keep_image returns the smoothed 32-bit tissue image")
    if kept is not None:
        kept.close()
    imp.close()


# -----------------------------------------------------------------------------
# Further known-answer cases added by the independent check
# -----------------------------------------------------------------------------

def tt_calibrated(title, ip, pixel_x, pixel_y, unit):
    imp = ImagePlus(title, ip)
    cal = imp.getCalibration()
    cal.pixelWidth = pixel_x
    cal.pixelHeight = pixel_y
    cal.setUnit(unit)
    return imp


def tt_frames(title, planes, pixel=0.5):
    """One channel, one Z, one frame per plane, calibrated in micrometers."""
    stack = ImageStack(planes[0].getWidth(), planes[0].getHeight())
    for plane in planes:
        stack.addSlice("", plane)
    imp = ImagePlus(title, stack)
    imp.setDimensions(1, 1, len(planes))
    imp.setOpenAsHyperStack(True)
    cal = imp.getCalibration()
    cal.pixelWidth = pixel
    cal.pixelHeight = pixel
    cal.setUnit("um")
    return imp


def test_tissue_anisotropic_edge(E, t, tmp):
    """Non-square pixels in 'micron': each direction is smoothed by its own
    sigma in pixels, which a straight edge shows."""
    # Tissue is the rectangle of pixel columns 0-99 and rows 50-249 (cut by
    # the left border), level 1000 on 100. Pixels are 0.5 x 0.25 micron, so
    # a sigma of 1 micron is 2 pixels along x and 4 along y. A pixel whose
    # centre lies at distance d outside a straight edge has, after the blur,
    # about 100 + 900 * Phi(-d / s), so the threshold L keeps the pixels
    # with d <= D = -s * InvPhi((L - 100) / 900). With pixel i covering
    # [i, i + 1): the right edge is at x = 100, so the last column kept is
    # floor(99.5 + D_x); the top edge is at y = 50, so the first row kept is
    # ceil(49.5 - D_y), and by symmetry the outline is 200 + 2 * (50 - first)
    # rows tall. The blur extends the border pixels outward, so the left
    # side stays at column 0. Tolerance 1 pixel: the sampled kernel is not
    # exactly Phi.
    rect = Roi(0, 50, 100, 200)
    imp = tt_calibrated("aniso.tif", tt_plane(100, [(rect, 1000)]), 0.5, 0.25, "micron")
    for log_scale in [True, False]:
        label = "anisotropic, " + ("log" if log_scale else "linear") + " scale: "
        rois, info = tt_detect(E, imp, {"sigma": "1", "log_scale": log_scale})
        t.equal(info["length_unit"], "micrometer", label + "'micron' is a micrometer calibration")
        t.close(info["sigma_pixels"], 2.0, label + "sigma along x is 1 / 0.5 = 2 pixels")
        t.close(info["sigma_pixels_y"], 4.0, label + "sigma along y is 1 / 0.25 = 4 pixels")
        t.equal(len(rois), 1, label + "one region")
        if len(rois) != 1:
            continue
        q = (info["threshold"] - 100.0) / 900.0
        t.check(0 < q < 1, label + "the threshold lies between the levels")
        d_x = -2.0 * tt_inverse_normal(q)
        d_y = -4.0 * tt_inverse_normal(q)
        t.check(abs(d_y - d_x) >= 3, "construction: the two directions differ by 3 pixels or more, so a swap would show")
        bounds = rois[0].getBounds()
        width = int(math.floor(99.5 + d_x)) + 1
        first = int(math.ceil(49.5 - d_y))
        t.equal(bounds.x, 0, label + "the outline stays on the left border")
        t.check(abs(bounds.width - width) <= 1, label + "right edge moved by sigma_x (" + str(bounds.width) + " vs " + str(width) + ")")
        t.check(abs(bounds.y - first) <= 1, label + "top edge moved by sigma_y (" + str(bounds.y) + " vs " + str(first) + ")")
        t.check(abs(bounds.height - (200 + 2 * (50 - first))) <= 2, label + "outline height")
        t.equal(info["border_sides"], "left", label + "touches the left border only")
        region = info["regions"][0]
        t.equal(region["area_unit"], "micrometer^2", label + "area unit")
        t.close(region["area"], region["area_pixels"] * 0.125, label + "area is pixels x 0.5 x 0.25 square micrometers")
        t.equal(E.tissue_tag(rois[0], "Sigma_Pixels_Y"), "4.0", label + "the record keeps sigma along y")
        t.equal(info["warnings"], [], label + "no warning without closing or shift")
    imp.close()


def test_tissue_shift_and_closing(E, t, tmp):
    """The inward shift and the gap closing, exact on rectangles."""
    # Inward shift. A rectangle cut by the left border (columns 0-99, rows
    # 50-249), pixels 0.5 x 0.25 micron, shift 0.7071 micron. The radius
    # uses the geometric mean pixel size sqrt(0.5 * 0.25) = 0.35355 micron,
    # so it is 2.0 pixels; ImageJ's Minimum filter of radius 2 reaches 2
    # pixels along x and y (its kernel holds dx^2 + dy^2 <= r^2 + 1). A
    # rectangle eroded by any kernel whose extent is +-2 is the rectangle
    # shrunk by 2 on each side; pixels beyond the border copy the border,
    # so the left side does not move: columns 0-97, rows 52-247.
    rect = Roi(0, 50, 100, 200)
    imp = tt_calibrated("shift.tif", tt_plane(100, [(rect, 1000)]), 0.5, 0.25, "micron")
    rois, info = tt_detect(E, imp, {"sigma": "0", "shrink": "0.7071"})
    t.equal(len(rois), 1, "shift: one region")
    truth = tt_truth([(Roi(0, 52, 98, 196), True)])
    t.equal(tt_mismatch(E, rois, truth), 0, "shift: the rectangle moved inward by 2 pixels, not from the border")
    t.equal(info["area_pixels"], 98 * 196, "shift: area")
    t.check(len([w for w in info["warnings"] if "not square" in w]) == 1, "shift: non-square pixels are warned about")
    t.check("not square" in E.to_text(info["description"]), "shift: the warning is in the description")
    imp.close()

    # Gap closing. Two rectangles (columns 20-119 and 121-220, rows 50-249)
    # separated by one empty column. Square 0.5 micron pixels; a closing
    # radius of 0.5 micron is 1 pixel, whose kernel is the full 3 x 3
    # square. Dilation by it fills the gap; the erosion then removes only
    # what the dilation added outside, so the result is exactly the one
    # rectangle of columns 20-220: 201 x 200 pixels.
    shapes = [(Roi(20, 50, 100, 200), 1000), (Roi(121, 50, 100, 200), 1000)]
    imp = tt_image("gap.tif", [tt_plane(100, shapes)], pixel=0.5)
    rois, info = tt_detect(E, imp, {"sigma": "0"})
    t.equal(len(rois), 2, "gap not closed: two regions (one empty column keeps them apart)")
    rois, info = tt_detect(E, imp, {"sigma": "0", "close_radius": "0.5"})
    t.equal(len(rois), 1, "gap closed: one region")
    t.equal(tt_mismatch(E, rois, tt_truth([(Roi(20, 50, 201, 200), True)])), 0, "gap closed: exactly the joined rectangle")
    t.equal(info["warnings"], [], "square pixels: no warning")
    imp.close()


def test_tissue_fills_field(E, t, tmp):
    """Tissue everywhere except an inner lumen: the outline is the whole
    image, touching all four sides."""
    # Level 1000 everywhere except a dark disk (a lumen, 100) inside. The
    # lumen is the only background and is not connected to the border, so
    # with fill holes the mask is the whole 300 x 300 image: border contact
    # is every border pixel once, 2 * 300 + 2 * 300 - 4 = 1196, the
    # centroid is (150, 150) and the filled pixels are the lumen. Without
    # fill holes the region is the image minus the lumen, and its centroid
    # is (90000 * 150 - n * c) / (90000 - n) for a lumen of n pixels with
    # centroid c.
    lumen = tt_disk(130, 120, 30)
    lumen_truth = tt_truth([(lumen, True)])
    n = tt_count(lumen_truth)
    c = tt_centroid(lumen_truth)
    imp = tt_image("field.tif", [tt_plane(1000, [(lumen, 100)])])
    rois, info = tt_detect(E, imp, {"sigma": "0"})
    t.equal(len(rois), 1, "filled: one region")
    t.equal(info["area_pixels"], TT_SIZE * TT_SIZE, "filled: the whole image")
    t.close(info["area_percent"], 100.0, "filled: 100 percent of the image")
    t.equal(info["filled_hole_pixels"], n, "filled: the lumen pixels were filled")
    t.equal(info["border_sides"], "left+top+right+bottom", "filled: all four sides")
    t.equal(info["regions"][0]["border_contact_pixels"], 4 * TT_SIZE - 4, "filled: each border pixel counted once")
    t.close(info["regions"][0]["centroid_x"], 150.0, "filled: centroid x")
    t.close(info["regions"][0]["centroid_y"], 150.0, "filled: centroid y")
    t.equal(len(E.validate_parent_rois(rois, TT_SIZE, TT_SIZE)), 1, "filled: a valid parent ROI")
    rois, info = tt_detect(E, imp, {"sigma": "0", "fill_holes": False})
    t.equal(len(rois), 1, "lumen kept: one region")
    t.equal(info["area_pixels"], TT_SIZE * TT_SIZE - n, "lumen kept: image minus lumen")
    t.equal(tt_mismatch(E, rois, tt_truth([(Roi(0, 0, TT_SIZE, TT_SIZE), True), (lumen, False)])), 0, "lumen kept: exact")
    total = float(TT_SIZE * TT_SIZE)
    t.close(info["regions"][0]["centroid_x"], (total * 150.0 - n * c[0]) / (total - n), "lumen kept: centroid x", 1e-9)
    t.close(info["regions"][0]["centroid_y"], (total * 150.0 - n * c[1]) / (total - n), "lumen kept: centroid y", 1e-9)
    t.equal(info["regions"][0]["border_contact_pixels"], 4 * TT_SIZE - 4, "lumen kept: still every border pixel")
    imp.close()


def test_tissue_bit_depths(E, t, tmp):
    """8-bit and 32-bit images, including negative values (background
    subtracted data) and values far below 1."""
    # Two levels and no smoothing: the outline is the drawn disk exactly,
    # and the means above and below the threshold are the two levels.
    disk = tt_disk(150, 140, 80)
    truth = tt_truth([(disk, True)])
    cases = [("8-bit", ByteProcessor(TT_SIZE, TT_SIZE), 10.0, 200.0),
             ("32-bit negative background", FloatProcessor(TT_SIZE, TT_SIZE), -50.0, 400.0),
             ("32-bit small values", FloatProcessor(TT_SIZE, TT_SIZE), 0.001, 0.02)]
    for name, ip, low, high in cases:
        ip.setValue(low)
        ip.fill()
        ip.setValue(high)
        ip.fill(disk)
        ip.resetRoi()
        imp = ImagePlus(name + ".tif", ip)
        for log_scale in [True, False]:
            label = name + ", " + ("log" if log_scale else "linear") + ": "
            rois, info = tt_detect(E, imp, {"sigma": "0", "log_scale": log_scale})
            t.equal(len(rois), 1, label + "one region")
            t.equal(tt_mismatch(E, rois, truth), 0, label + "the disk exactly")
            t.check(low < info["threshold"] <= high, label + "threshold between the levels")
            t.close(info["mean_above"], high, label + "mean above the threshold is the tissue level", 1e-6)
            t.close(info["mean_below"], low, label + "mean below the threshold is the background level", 1e-6)
        imp.close()


def test_tissue_several_parents_two_frames(E, t, tmp):
    """Two detected regions beside a hand-drawn ROI, analysed over two
    frames: one record row per detected ROI, with the right index."""
    big = tt_disk(90, 150, 60)
    small = tt_disk(230, 150, 35)
    big_truth = tt_truth([(big, True)])
    small_truth = tt_truth([(small, True)])
    plane = tt_plane(100, [(big, 1000), (small, 600)])
    imp = tt_frames("frames.tif", [plane, plane.duplicate()])
    rois, info = tt_detect(E, imp, {"sigma": "0"})
    t.equal(len(rois), 2, "construction: two regions")
    hand = Roi(10, 10, 50, 40)
    hand.setName("hand")
    configs = [nuclei_config(E, {"manual_low": 500})]
    result = run_image(E, imp, configs, [hand] + rois, os.path.join(tmp, "tissue_frames"), frames=[1, 2],
                       options_override={"tissue_sigma": "0"})
    t.equal(len(result["summary_rows"]), 6, "construction: 3 ROIs x 2 frames in roi_summary")
    rows = result["extra_tables"].get("tissue_outline.csv", {"rows": []})["rows"]
    t.equal(len(rows), 2, "one tissue row per detected ROI, not per frame")
    expected = [(2, "Tissue_001", big_truth), (3, "Tissue_002", small_truth)]
    for index in range(min(len(rows), 2)):
        row = rows[index]
        number, name, truth = expected[index]
        label = name + ": "
        t.equal(int(row["Parent_ROI_Index"]), number, label + "Parent_ROI_Index counts the hand-drawn ROI too")
        t.equal(E.to_text(row["Parent_ROI_Name"]), name, label + "name")
        pixels = tt_count(truth)
        t.equal(int(row["ROI_Area_Pixels"]), pixels, label + "ROI_Area_Pixels")
        t.close(float(row["ROI_Area"]), pixels * 0.25, label + "ROI_Area in square micrometers")
        t.close(float(row["Image_Area_Percent"]), 100.0 * pixels / (TT_SIZE * TT_SIZE), label + "Image_Area_Percent")
        centre = tt_centroid(truth)
        t.close(float(row["Centroid_X_Pixels"]), centre[0], label + "Centroid_X_Pixels")
        t.close(float(row["Centroid_Y_Pixels"]), centre[1], label + "Centroid_Y_Pixels")
        t.equal(E.to_text(row["Regions_Kept"]), "2", label + "Regions_Kept")
        t.equal(E.to_text(row["Outline_Status"]), "as detected", label + "status")
        t.equal(E.to_text(row["Touches_Image_Border"]), "false", label + "border flag")
    saved = os.path.join(result["run_dir"], "ROIs", "Tissue_Outline_ROIs.zip")
    names = [E.to_text(roi.getName()) for roi in E.load_rois_from_file(saved)] if os.path.isfile(saved) else []
    t.equal(names, ["Tissue_001", "Tissue_002"], "only the detected ROIs go into the tissue ROI file")
    imp.close()

    # Detection on frame 2 of a time series whose frames differ.
    plane_1 = tt_plane(100, [(big, 1000)])
    plane_2 = tt_plane(100, [(small, 1000)])
    imp = tt_frames("timelapse.tif", [plane_1, plane_2])
    rois, info = E.detect_tissue_outline(None, imp, [], tt_options(E, {"sigma": "0"}), 2)
    t.equal(tt_mismatch(E, rois, small_truth), 0, "frame 2 is used when asked for")
    t.equal((info["frame"], E.tissue_tag(rois[0], "Frame") if rois else ""), (2, "2"), "the frame is recorded")
    imp.close()


def test_tissue_global_calibration_and_units(E, t, tmp):
    """Pixel results do not follow Analyze > Set Scale > Global, and length
    units are recognised or warned about."""
    disk = tt_disk(150, 140, 80)
    expected = tt_centroid(tt_truth([(disk, True)]))
    probe = ImagePlus("probe", ByteProcessor(4, 4))
    global_cal = Calibration()
    global_cal.pixelWidth = 0.5
    global_cal.pixelHeight = 0.5
    global_cal.setUnit("um")
    saved = probe.getGlobalCalibration()
    try:
        probe.setGlobalCalibration(global_cal)
        imp = ImagePlus("global.tif", tt_plane(100, [(disk, 1000)]))
        rois, info = tt_detect(E, imp, {"sigma": "0"})
        t.close(info["regions"][0]["centroid_x"], expected[0], "global calibration: centroid x stays in pixels")
        t.close(info["regions"][0]["centroid_y"], expected[1], "global calibration: centroid y stays in pixels")
        t.equal(E.tissue_tag(rois[0], "Centroid_X"), E.to_text(repr(info["regions"][0]["centroid_x"])), "the record holds the same value")
        imp.close()
    finally:
        probe.setGlobalCalibration(saved)
    t.equal(probe.getGlobalCalibration(), saved, "the global calibration is restored")

    for unit, size, expected_scale in [("micrometre", 0.5, ("micrometer", 0.5)),
                                       ("nanometers", 200.0, ("micrometer", 0.2)),
                                       ("millimetre", 0.001, ("micrometer", 1.0))]:
        imp = tt_calibrated("unit.tif", tt_plane(0, [], width=8, height=8), size, size, unit)
        scale = E.tissue_length_scale(imp)
        t.equal(scale[0], expected_scale[0], "unit '" + unit + "' is a length unit")
        t.close(scale[1], expected_scale[1], "unit '" + unit + "' in micrometers per pixel")
        imp.close()
    imp = tt_calibrated("furlong.tif", tt_plane(100, [(disk, 1000)]), 0.5, 0.5, "furlong")
    rois, info = tt_detect(E, imp, {"sigma": "2"})
    t.equal(info["length_unit"], "pixel", "an unknown unit: lengths in pixels")
    t.close(info["sigma_pixels"], 2.0, "an unknown unit: sigma 2 is 2 pixels")
    t.check("not a known length unit" in E.to_text(info["description"]), "an unknown unit is warned about")
    imp.close()
    imp = tt_image("plain.tif", [tt_plane(100, [(disk, 1000)])])
    rois, info = tt_detect(E, imp, {"sigma": "2"})
    t.equal(info["warnings"], [], "an uncalibrated image gives no unit warning")
    imp.close()


def test_tissue_batch_same_names(E, t, tmp):
    """Two images with the same name in different subfolders: the second
    outline must not overwrite the first in Tissue_Outlines."""
    folder = os.path.join(tmp, "tissue_same_in")
    for sub, shapes in [("a", [(tt_disk(150, 140, 80), 1000)]),
                        ("b", [(tt_disk(90, 150, 60), 1000), (tt_disk(230, 150, 35), 600)])]:
        os.makedirs(os.path.join(folder, sub))
        imp = tt_image("x.tif", [tt_plane(100, shapes)], pixel=0.5)
        FileSaver(imp).saveAsTiff(os.path.join(folder, sub, "x.tif"))
        imp.close()
    options = E.default_global_options()
    for key in ["save_prepared", "save_masks", "save_overlays", "save_object_rois"]:
        options[key] = False
    options["tissue_sigma"] = "0"
    configs = [nuclei_config(E, {"manual_low": 500})]
    result = E.run_batch(folder, os.path.join(tmp, "tissue_same_out"), configs, configs, options,
                         {"extensions": "tif", "roi_policy": E.TISSUE_POLICY, "recursive": True})
    t.equal(result["analyzed"], 2, "both images analyzed")
    collection = os.path.join(result["batch_dir"], E.TISSUE_COLLECTION_DIR)
    files = sorted(os.listdir(collection)) if os.path.isdir(collection) else []
    t.equal(files, ["x_2_ROIs.zip", "x_2_Tissue_Outline_QC.tif", "x_ROIs.zip", "x_Tissue_Outline_QC.tif"],
            "the second image's outline is collected under a new name")
    counts = sorted([len(E.load_rois_from_file(os.path.join(collection, name))) for name in files if name.endswith(".zip")])
    t.equal(counts, [1, 2], "both outlines are kept (1 and 2 regions)")


TESTS.append(("tissue outline: disk of known radius", test_tissue_disk))
TESTS.append(("tissue outline: two regions, largest only, minimum area", test_tissue_two_blobs))
TESTS.append(("tissue outline: tissue cut by the image border", test_tissue_border))
TESTS.append(("tissue outline: images without tissue", test_tissue_empty))
TESTS.append(("tissue outline: holes filled or kept", test_tissue_holes))
TESTS.append(("tissue outline: independent of black background and edge padding", test_tissue_prefs_independent))
TESTS.append(("tissue outline: channel choice and Z processing", test_tissue_channels))
TESTS.append(("tissue outline: settings and profiles", test_tissue_settings))
TESTS.append(("tissue outline: used as parent ROI in a run", test_tissue_in_analysis))
TESTS.append(("tissue outline: nothing changes when unused", test_tissue_unused_changes_nothing))
TESTS.append(("tissue outline: batch ROI policy", test_tissue_batch_policy))
TESTS.append(("tissue outline: check image", test_tissue_check_image))
TESTS.append(("tissue outline: non-square micron pixels, straight edge", test_tissue_anisotropic_edge))
TESTS.append(("tissue outline: inward shift and gap closing on rectangles", test_tissue_shift_and_closing))
TESTS.append(("tissue outline: tissue fills the field around a lumen", test_tissue_fills_field))
TESTS.append(("tissue outline: 8-bit and 32-bit images", test_tissue_bit_depths))
TESTS.append(("tissue outline: several parent ROIs over two frames", test_tissue_several_parents_two_frames))
TESTS.append(("tissue outline: global calibration and length units", test_tissue_global_calibration_and_units))
TESTS.append(("tissue outline: batch images with the same name", test_tissue_batch_same_names))

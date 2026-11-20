# Cell Quant Studio Universal - tests of the robustness module
# (cellquantrobust.py: a group difference across a grid of settings)
#
# Loaded by run_tests.py, which provides the helpers used here.
# The source is kept plain ASCII with no encoding declaration on purpose.

import re
from ij.process import FloatProcessor


def robust_close_each(t, actual, expected, message):
    """t.close for each pair of two equally long sequences."""
    t.equal(len(actual), len(expected), message + ": number of values")
    for index in range(min(len(actual), len(expected))):
        t.close(actual[index], expected[index], message + " [" + str(index + 1) + "]")

#
# Every image is WIDTH x HEIGHT = 240 x 200 pixels at 0.5 micrometer per
# pixel, so the whole image is 48000 pixels = 12000 square micrometers.
# Objects are disks of radius 6 at the positions in DISKS, which do not
# touch each other or the image border; the first n positions give n disks.

#
# ImageJ stores the unit "um" as the micro sign followed by "m", so that is
# what the tables hold.

ROBUST_AREA = WIDTH * HEIGHT * PIXEL * PIXEL
ROBUST_UM = MICRO + "m"


def robust_plane(disk_count, value, salt=None, extra_disks=None, background=100):
    ip = blank_plane(background)
    for disk in DISKS[:disk_count] + (extra_disks or []):
        draw_disk(ip, disk[0], disk[1], DISK_RADIUS, value)
    for point in (salt or []):
        ip.set(point[0], point[1], int(value))
    return ip


def robust_save(folder, name, ip):
    imp = ImagePlus(name, ip)
    cal = imp.getCalibration()
    cal.pixelWidth = PIXEL
    cal.pixelHeight = PIXEL
    cal.setUnit("um")
    FileSaver(imp).saveAsTiff(os.path.join(folder, name))
    imp.close()


def robust_folder(tmp, name, images, background=100):
    """images: list of (file name, number of disks, disk value)."""
    folder = os.path.join(tmp, name)
    os.makedirs(folder)
    for image in images:
        robust_save(folder, image[0], robust_plane(image[1], image[2], background=background))
    return folder


def robust_without_run_names(text):
    """A batch table with the time stamps of the batch and run folders
    (CellQuant_Batch_<date>_<time>, <image>_CellQuant_<date>_<time>, with
    _2, _3 ... when a name was taken) replaced by a fixed word."""
    return re.sub("CellQuant_(Batch_)?[0-9]{8}_[0-9]{6}(_[0-9]+)?", "CellQuant_ID", text)


def robust_table(E, folder, name):
    """A written table as a list of dictionaries (every cell is text)."""
    rows = E.robust_parse_table(read_utf8(os.path.join(folder, name)))
    return [dict(zip(rows[0], row)) for row in rows[1:]]


def robust_run(E, folder, tmp, grid, groups=None, config=None, **keywords):
    configs = [config if config is not None else nuclei_config(E)]
    options = E.default_global_options()
    batch_options = keywords.pop("batch_options", {"extensions": "tif"})
    return E.run_robustness(
        folder, os.path.join(tmp, "robust_out"), configs, configs, options, batch_options,
        configs[0]["channel"], grid, groups if groups is not None else {"control_group": "ctrl"}, **keywords
    )


def robust_refused(function, *arguments):
    """The message of the ValueError a call raises, or None."""
    try:
        function(*arguments)
    except ValueError as exc:
        # The message itself, not str() of it: it may hold non-ASCII text.
        return exc.args[0] if len(exc.args) > 0 else u""
    return None


def robust_t4_two_sided_p(t):
    """Two-sided p of Student's t with 4 degrees of freedom, from the closed
    form of its distribution function:
    F(t) = 1/2 + 3/8 * (t / sqrt(1 + t*t/4)) * (1 - t*t / (12 * (1 + t*t/4)))."""
    q = 1.0 + t * t / 4.0
    cdf = 0.5 + 0.375 * (t / math.sqrt(q)) * (1.0 - t * t / (12.0 * q))
    return 2.0 * (1.0 - cdf)


# Quantile 0.975 of Student's t with 4 degrees of freedom (any t table).
ROBUST_T4_975 = 2.7764451


def test_robust_statistics(E, t, tmp):
    """The group statistics against numbers derived by hand."""
    # A = 10, 12, 14: mean 12, variance (4 + 0 + 4) / 2 = 4, SD 2.
    # B = 17, 19, 21: mean 19, variance 4, SD 2.
    # difference = 7; ratio = 19 / 12.
    # pooled SD = sqrt((2*4 + 2*4) / 4) = 2; d = 3.5; J = 1 - 3/(4*4 - 1) = 0.8; g = 2.8.
    # Welch: standard error = sqrt(4/3 + 4/3) = sqrt(8/3); t = 7 / sqrt(8/3) = 4.2866070;
    # df = (8/3)^2 / ((4/3)^2/2 + (4/3)^2/2) = (64/9) / (16/9) = 4.
    # 95 percent interval = 7 +- 2.7764451 * sqrt(8/3) = 7 +- 4.5339152.
    # Mann-Whitney: every B exceeds every A, so U of B = 9 and the rank sum of
    # B is as extreme as it can be; of the C(6,3) = 20 ways to divide six
    # ranks, 2 are that extreme (all low or all high): exact p = 0.1.
    s = E.robust_two_group_statistics([10, 12, 14], [17, 19, 21])
    t.equal((s["n_a"], s["n_b"]), (3, 3), "three units per group")
    t.close(s["mean_a"], 12.0, "mean of A")
    t.close(s["sd_b"], 2.0, "SD of B")
    t.close(s["difference"], 7.0, "difference is B minus A")
    t.close(s["ratio"], 19.0 / 12.0, "ratio is B over A")
    t.close(s["hedges_g"], 2.8, "Hedges g with the small-sample correction")
    t.close(s["welch_t"], 7.0 / math.sqrt(8.0 / 3.0), "Welch t")
    t.close(s["welch_df"], 4.0, "Welch degrees of freedom")
    t.equal(s["sign"], 1, "B larger than A is a positive sign")
    t.close(s["mw_u"], 9.0, "Mann-Whitney U of B")
    t.close(s["mw_p"], 0.1, "exact Mann-Whitney p for complete separation of 3 and 3")
    t.equal(s["mw_method"], "exact permutation", "small samples use the exact test")
    if E.ROBUST_MATH["available"]:
        t.close(s["welch_p"], robust_t4_two_sided_p(7.0 / math.sqrt(8.0 / 3.0)), "Welch p equals the closed form for 4 df")
        t.close(s["welch_p"], 0.012780, "Welch p, value worked out by hand", 1e-5)
        t.close(s["ci_low"], 7.0 - ROBUST_T4_975 * math.sqrt(8.0 / 3.0), "lower end of the 95 percent interval")
        t.close(s["ci_high"], 7.0 + ROBUST_T4_975 * math.sqrt(8.0 / 3.0), "upper end of the 95 percent interval")

    # Swapping the groups reverses the signs and keeps the p values.
    r = E.robust_two_group_statistics([17, 19, 21], [10, 12, 14])
    t.close(r["hedges_g"], -2.8, "g changes sign when the groups are swapped")
    t.equal(r["sign"], -1, "B smaller than A is a negative sign")
    t.close(r["mw_p"], 0.1, "Mann-Whitney p does not depend on the order of the groups")
    t.close(r["mw_u"], 0.0, "U of B is zero when every B is below every A")

    # Unequal variances and sizes: A = 10, 12, 14 (variance 4), B = 15, 19, 23, 27
    # (mean 21, variance (36 + 4 + 4 + 36) / 3 = 80/3).
    # Welch: se^2 = 4/3 + (80/3)/4 = 8; t = 9 / sqrt(8) = 3.1819805;
    # df = 8^2 / ((4/3)^2/2 + (20/3)^2/3) = 64 / (8/9 + 400/27) = 64 / (424/27) = 4.0754717.
    # pooled variance = (2*4 + 3*80/3) / 5 = 17.6; J = 1 - 3/(4*5 - 1) = 16/19;
    # g = (16/19) * 9 / sqrt(17.6) = 1.8065556.
    u = E.robust_two_group_statistics([10, 12, 14], [15, 19, 23, 27])
    t.close(u["welch_t"], 9.0 / math.sqrt(8.0), "Welch t with unequal variances")
    t.close(u["welch_df"], 64.0 * 27.0 / 424.0, "Welch-Satterthwaite degrees of freedom")
    t.close(u["hedges_g"], 16.0 / 19.0 * 9.0 / math.sqrt(17.6), "Hedges g with unequal group sizes")
    # every B above every A: 2 of C(7,3) = 35 divisions are as extreme
    t.close(u["mw_p"], 2.0 / 35.0, "exact Mann-Whitney p for complete separation of 3 and 4")

    # Ties: A = 1, 2, 2 and B = 2, 3, 4, 5. Pooled ranks: 1 -> 1; the three
    # 2s share (2+3+4)/3 = 3; 3 -> 5; 4 -> 6; 5 -> 7. Rank sum of A = 7,
    # expected 3 * 8 / 2 = 12, distance 5. Of the 35 ways to choose three of
    # the ranks {1, 3, 3, 3, 5, 6, 7}, those with a sum of at most 7 are
    # {1, 3, 3} (3 ways) and those with at least 17 are {5, 6, 7} (1 way):
    # p = 4 / 35. U of A = 7 - 6 = 1, so U of B = 12 - 1 = 11.
    w = E.robust_two_group_statistics([1, 2, 2], [2, 3, 4, 5])
    t.equal(E.robust_midranks([1, 2, 2, 2, 3, 4, 5]), [1.0, 3.0, 3.0, 3.0, 5.0, 6.0, 7.0], "tied values share their mean rank")
    t.close(w["mw_u"], 11.0, "U of B with ties")
    t.close(w["mw_p"], 4.0 / 35.0, "exact Mann-Whitney p with ties")

    # Fewer than three units in a group: means only, and the reason.
    few = E.robust_two_group_statistics([1.0, 3.0], [5.0, 6.0, 7.0])
    t.close(few["difference"], 4.0, "with two units the difference of the means is still reported")
    t.check(few["hedges_g"] is None and few["welch_p"] is None and few["mw_p"] is None and few["ci_low"] is None,
            "with two units in a group there is no effect size, interval, or test")
    t.check("Fewer than 3 units" in few["note"], "and the note says why")

    # Degenerate inputs.
    flat = E.robust_two_group_statistics([2.0, 2.0, 2.0], [5.0, 5.0, 5.0])
    t.check(flat["hedges_g"] is None and flat["welch_t"] is None, "no variation between units: g and t are undefined")
    t.check("No variation" in flat["note"], "and the note says so")
    t.close(flat["mw_p"], 0.1, "the rank test is still defined without variation inside the groups")
    same = E.robust_two_group_statistics([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    t.equal(same["sign"], 0, "equal means have sign zero")
    t.check(same["ratio"] is None, "a ratio over a zero mean is left empty")
    t.close(same["mw_p"], 1.0, "identical groups give a rank-test p of 1")
    none = E.robust_two_group_statistics([], [1.0, 2.0, 3.0])
    t.check(none["difference"] is None and none["sign"] is None, "a group without units gives no difference")

    # Larger samples (more than 20000 ways to divide the ranks): the normal
    # approximation with continuity correction.
    # A = 1, 2, ..., 11 and B = 4.5, 5.5, ..., 15.5 (12 values), no ties.
    # U of B = pairs with b > a: for b = k + 0.5 that is min(k, 11) values
    # of A: 4+5+6+7+8+9+10+11 + 4*11 = 104. Expected 11*12/2 = 66;
    # variance = 11*12*(23+1)/12 = 264; z = (104 - 66 - 0.5) / sqrt(264) = 2.30797;
    # two-sided p = erfc(z / sqrt(2)) = 0.0210.
    big_a = [float(v) for v in range(1, 12)]
    big_b = [float(v) + 0.5 for v in range(4, 16)]
    large = E.robust_two_group_statistics(big_a, big_b)
    t.equal(large["mw_method"], "normal approximation (tie and continuity corrected)", "large samples use the normal approximation")
    t.close(large["mw_u"], 104.0, "U of B for the large samples")
    t.close(large["mw_p"], math.erfc(37.5 / math.sqrt(264.0) / math.sqrt(2.0)), "normal approximation with continuity correction")
    t.close(large["mw_p"], 0.0210, "normal approximation, value worked out by hand", 2e-4)
    t.close(E.robust_two_group_statistics(big_b, big_a)["mw_p"], large["mw_p"], "the approximation does not depend on the order of the groups")
    # With ties: A = five 1s and five 2s, B = five 2s and five 3s (C(20,10) =
    # 184756 ways). Mid-ranks: the 1s share 3, the ten 2s share 10.5, the 3s
    # share 18. Rank sum of B = 5*10.5 + 5*18 = 142.5; U of B = 142.5 - 55 = 87.5;
    # expected 50. Tie term = (5^3-5) + (10^3-10) + (5^3-5) = 1230;
    # variance = 100/12 * (21 - 1230/(20*19)) = 148.026316;
    # z = (37.5 - 0.5) / sqrt(148.026316) = 3.04111; p = erfc(z / sqrt(2)) = 0.00236.
    tied = E.robust_mann_whitney([1.0] * 5 + [2.0] * 5, [2.0] * 5 + [3.0] * 5)
    t.close(tied[0], 87.5, "U of B with many ties")
    t.close(tied[1], math.erfc(37.0 / math.sqrt(100.0 / 12.0 * (21.0 - 1230.0 / 380.0)) / math.sqrt(2.0)), "normal approximation with the correction for ties")
    t.close(tied[1], 0.00236, "tie-corrected approximation, value worked out by hand", 1e-5)
    t.close(E.robust_mann_whitney([4.0] * 10, [4.0] * 10)[1], 1.0, "all values equal: p is 1")
    # The exact test and the approximation agree roughly where both apply:
    # 8 against 8 completely separated, exact p = 2 / C(16,8) = 2 / 12870.
    t.close(E.robust_mann_whitney([float(v) for v in range(8)], [float(v) for v in range(8, 16)])[1], 2.0 / 12870.0,
            "exact p for complete separation of 8 and 8")

    # Without Commons Math: everything that needs no distribution is still there.
    saved = E.ROBUST_MATH["available"]
    E.ROBUST_MATH["available"] = False
    try:
        bare = E.robust_two_group_statistics([10, 12, 14], [17, 19, 21])
        t.close(bare["hedges_g"], 2.8, "without Commons Math: Hedges g is unchanged")
        t.check(bare["welch_p"] is None and bare["ci_low"] is None, "without Commons Math: p and interval are empty")
        t.check("Commons Math is missing" in bare["note"], "without Commons Math: the note says so")
        t.close(bare["mw_p"], 0.1, "without Commons Math: the exact rank test still works")
        t.close(E.robust_two_group_statistics(big_a, big_b)["mw_p"], large["mw_p"], "without Commons Math: so does its approximation")
    finally:
        E.ROBUST_MATH["available"] = saved

    t.close(E.robust_median([3.0, 1.0, 2.0]), 2.0, "median of an odd number of values")
    t.close(E.robust_median([4.0, 1.0, 2.0, 3.0]), 2.5, "median of an even number of values")


def test_robust_grid(E, t, tmp):
    """The grid of specifications: full factorial, reference included,
    duplicates removed, and the refusals."""
    config = nuclei_config(E)
    density = E.ROBUST_METRIC_DENSITY
    t.equal(E.parse_robust_grid("threshold_factor=0.8, 1.2; watershed=true,false"),
            {"threshold_factor": ["0.8", "1.2"], "watershed": ["true", "false"]}, "the grid text is parsed")

    plan = E.robust_build_specifications(config, {"threshold_factor": [1.2, 0.8], "watershed": [True]}, density)
    specs = plan["specs"]
    t.equal(len(specs), 6, "3 threshold factors (1.0 added) x 2 watershed values (false added) = 6 specifications")
    t.equal([spec["id"] for spec in specs], ["S001", "S002", "S003", "S004", "S005", "S006"], "specifications are numbered")
    references = [spec for spec in specs if spec["is_reference"]]
    t.equal(len(references), 1, "exactly one specification is the reference")
    t.equal(references[0]["label"], "threshold_factor=1.0; watershed=false", "the reference is the configured setting")
    t.check(references[0]["config"] == config, "the reference specification is the configuration itself")
    t.equal(E.robust_grid_text(plan["axes"]), "threshold_factor=0.8,1.0,1.2; watershed=false,true",
            "the recorded grid lists every value, the configured one included")
    text_plan = E.robust_build_specifications(config, "threshold_factor=0.8,1.2; watershed=true", density)
    t.equal([spec["label"] for spec in text_plan["specs"]], [spec["label"] for spec in specs], "a grid given as text is the same grid")

    # A setting that has no effect in part of the grid does not multiply the
    # specifications: with an automatic threshold the 3 factors count (the
    # local radius does not), with the local threshold the 2 radii count
    # (the factor does not): 3 + 2 = 5, not 2 x 3 x 2 = 12.
    auto = nuclei_config(E, {"threshold_mode": "Auto"})
    mixed = E.robust_build_specifications(
        auto, {"threshold_mode": ["Auto", "Local Phansalkar"], "threshold_factor": [0.8, 1.2], "local_radius": [20]}, density
    )
    t.equal(len(mixed["specs"]), 5, "combinations that segment identically are one specification")
    labels = [spec["label"] for spec in mixed["specs"]]
    t.check("threshold_factor=0.8; threshold_mode=Auto" in labels, "a label lists only the settings that apply (automatic)")
    t.check("threshold_mode=Local Phansalkar; local_radius=20" in labels, "a label lists only the settings that apply (local)")

    message = robust_refused(E.robust_build_specifications, config, {"local_radius": [20, 50]}, density)
    t.check(message is not None and "no effect" in message, "a local-threshold setting on a manual threshold is refused")
    message = robust_refused(E.robust_build_specifications, config, {"watershed": [True]}, E.ROBUST_METRIC_AREA_FRACTION)
    t.check(message is not None and "no effect" in message, "an object setting is refused for the positive-area fraction")
    message = robust_refused(E.robust_build_specifications, config, {"lut": ["Red"]}, density)
    t.check(message is not None and "not a setting that can be varied" in message, "an unknown grid setting is refused")
    message = robust_refused(E.robust_build_specifications, config, {"threshold_factor": [1.0]}, density)
    t.check(message is not None and "only the configured setting" in message, "a grid without any alternative is refused")
    message = robust_refused(E.robust_build_specifications, config, {"watershed": ["perhaps"]}, density)
    t.check(message is not None and "true or false" in message, "a yes/no setting does not accept other text")
    message = robust_refused(E.robust_build_specifications, config, {"threshold_mode": ["Display minimum"]}, density)
    t.check(message is not None and "Display-minimum" in message, "a Display-minimum specification is refused")
    message = robust_refused(E.robust_build_specifications, nuclei_config(E, {"min_size": 5, "max_size": 50}),
                             {"min_size": [100]}, density)
    t.check(message is not None and "not usable" in message, "a specification the engine would reject is refused before any file is read")

    # The cap: 4 x 4 x 4 = 64 specifications against the default limit of 60.
    big = {"threshold_factor": [0.7, 0.85, 1.15], "min_size": [5, 10, 20], "max_circularity": [0.7, 0.8, 0.9]}
    message = robust_refused(E.robust_build_specifications, config, big, density)
    t.check(message is not None and "64 specifications" in message and "60" in message,
            "more specifications than the limit are refused, not cut")
    t.equal(len(E.robust_build_specifications(config, big, density, 64)["specs"]), 64, "the limit can be raised on purpose")

    # Requests that are refused before anything is read or written.
    configs = [config]
    plan_call = lambda **changes: E.robust_plan(
        changes.get("configs", configs), changes.get("channel", 1), {"threshold_factor": [0.8]}, {},
        changes.get("metric", density), changes.get("aggregation", E.ROBUST_UNIT_MEAN), 60
    )
    display = [nuclei_config(E, {"threshold_mode": "Display minimum"})]
    message = robust_refused(lambda: plan_call(configs=display))
    t.check(message is not None and "Display-minimum" in message, "a Display-minimum threshold is refused as in a batch")
    message = robust_refused(lambda: plan_call(metric=E.ROBUST_METRIC_INTENSITY))
    t.check(message is not None and "does not depend on the segmentation" in message, "mean intensity is refused with the reason")
    message = robust_refused(lambda: plan_call(metric=E.ROBUST_METRIC_AREA_FRACTION, configs=[nuclei_config(E, {"analysis_mode": "Objects only"})]))
    t.check(message is not None and "Positive area" in message, "the area fraction needs a positive-area channel")
    message = robust_refused(lambda: plan_call(channel=3))
    t.check(message is not None and "not enabled" in message, "a channel that is not analysed is refused")

    request = E.robust_request_from_options(E.default_global_options())
    t.equal(request["grid"], {"threshold_factor": ["0.8", "0.9", "1.0", "1.1", "1.2"]}, "the default settings hold a threshold sweep")
    t.equal(request["metric"], density, "the default metric is the object density")
    t.equal(request["max_specifications"], 60, "the default limit is 60 specifications")


def test_robust_groups(E, t, tmp):
    """Group and unit assignment by file-name rule, expression, and table."""
    names = ["ctrl_o1_a.tif", "ctrl_o1_b.tif", "mut_o2_a.tif", "Blank.tif", os.path.join("sub", "mut_o3_a.tif")]

    rule = E.robust_group_rule({})
    rows = E.robust_assign_groups(names, rule)
    t.equal([row["group"] for row in rows], ["ctrl", "ctrl", "mut", "", "mut"], "default rule: text before the first underscore")
    t.equal(rows[0]["unit"], "ctrl_o1_a.tif", "default unit: the image file")
    t.equal(rows[4]["unit"], "sub/mut_o3_a.tif", "the unit of a file in a subfolder is its relative path")
    t.equal(rows[3]["status"], "NOT ASSIGNED", "a file name without the separator gets no group")
    t.check("no text before" in rows[3]["message"], "and the reason is given")

    rule = E.robust_group_rule({
        "source": E.ROBUST_GROUP_REGEX, "group_regex": "(?:^|/)([a-z]+)_", "unit_regex": "([a-z]+_o[0-9]+)_"
    })
    rows = E.robust_assign_groups(names, rule)
    t.equal([row["group"] for row in rows], ["ctrl", "ctrl", "mut", "", "mut"], "group by regular expression")
    t.equal([row["unit"] for row in rows], ["ctrl_o1", "ctrl_o1", "mut_o2", "", "mut_o3"], "unit by regular expression: the organoid")
    t.equal(rows[3]["status"], "NOT ASSIGNED", "a path the expression does not match gets no group")
    message = robust_refused(E.robust_group_rule, {"source": E.ROBUST_GROUP_REGEX, "group_regex": "[a-z]+_"})
    t.check(message is not None and "one pair of parentheses" in message, "an expression without a capture group is refused")
    message = robust_refused(E.robust_group_rule, {"source": E.ROBUST_GROUP_REGEX, "group_regex": "([a-z"})
    t.check(message is not None and "not a valid regular expression" in message, "a malformed expression is refused")
    message = robust_refused(E.robust_group_rule, {"source": E.ROBUST_GROUP_CSV, "csv_path": os.path.join(tmp, "missing.csv")})
    t.check(message is not None and "does not exist" in message, "a missing group table is refused")

    # A table as a spreadsheet would save it: semicolons, a quoted cell,
    # a row for an image that does not exist, and no row for Blank.tif.
    table = os.path.join(tmp, "groups.csv")
    handle = codecs.open(table, "w", "utf-8")
    handle.write(u"Image;Group;Unit\r\n")
    handle.write(u"ctrl_o1_a.tif;control;org 1\r\n")
    handle.write(u"ctrl_o1_b;control;org 1\r\n")
    handle.write(u"\"mut_o2_a.tif\";\"mutant; line 2\";\r\n")
    handle.write(u"sub/mut_o3_a.tif;mutant; line 2;org 3\r\n".replace(u"mutant; line 2", u"\"mutant; line 2\""))
    handle.write(u"gone.tif;control;org 9\r\n")
    handle.close()
    rule = E.robust_group_rule({"source": E.ROBUST_GROUP_CSV, "csv_path": table})
    rows = E.robust_assign_groups(names, rule)
    t.equal(len(rows), 6, "one row per image and one per table row without an image")
    t.equal([row["group"] for row in rows[:5]], ["control", "control", "mutant; line 2", "", "mutant; line 2"], "group from the table")
    t.equal([row["unit"] for row in rows[:5]], ["org 1", "org 1", "mut_o2_a.tif", "", "org 3"],
            "unit from the table, or the image file where the cell is empty")
    t.equal(rows[3]["status"], "NOT ASSIGNED", "an image the table does not list gets no group")
    t.equal((rows[5]["status"], rows[5]["relative"]), ("NO IMAGE FOR TABLE ROW", "gone.tif"), "a table row without an image is reported")

    t.equal(E.robust_group_order(["mut", "ctrl", "Alpha"], ""), ["Alpha", "ctrl", "mut"], "groups in alphabetical order")
    t.equal(E.robust_group_order(["mut", "ctrl", "Alpha"], "mut"), ["mut", "Alpha", "ctrl"], "the control group comes first")
    message = robust_refused(E.robust_group_order, ["mut", "ctrl"], "wt")
    t.check(message is not None and "not among the groups" in message, "an unknown control group is refused")


def test_robust_more_objects(E, t, tmp):
    """Group B has more disks than group A at every setting: the direction
    holds in every specification, with values known at the reference."""
    # Three images (= three units) per group, disks of value 2000 on 100:
    #   ctrl: 3, 4, 5 disks -> densities 3, 4, 5 per 12000 um^2; mean 4, SD 1
    #   mut:  7, 8, 9 disks -> mean 8, SD 1 (all per 12000 um^2)
    # Thresholds 500 to 1500 lie between background and disks, and a disk
    # (more than 100 pixels = 25 um^2) passes a minimum size of 10 um^2, so
    # every one of the 5 x 2 specifications counts the same disks.
    # difference = 4/12000; ratio = 2; pooled SD = 1/12000, d = 4, J = 0.8,
    # g = 3.2; Welch t = 4 / sqrt(1/3 + 1/3) = 4.8989795 with 4 df;
    # Mann-Whitney exact p = 2/20 = 0.1 (complete separation, 3 against 3).
    images = [("ctrl_a.tif", 3, 2000), ("ctrl_b.tif", 4, 2000), ("ctrl_c.tif", 5, 2000),
              ("mut_a.tif", 7, 2000), ("mut_b.tif", 8, 2000), ("mut_c.tif", 9, 2000)]
    folder = robust_folder(tmp, "robust_more", images)
    handle = open(os.path.join(folder, "mut_broken.tif"), "w")
    handle.write("not an image")
    handle.close()
    robust_save(folder, "Blank.tif", robust_plane(2, 2000))
    before = sorted(os.listdir(folder))

    configs = [nuclei_config(E)]
    options = E.default_global_options()
    for key in ["save_prepared", "save_masks", "save_overlays", "save_object_rois", "save_parent_rois"]:
        options[key] = False
    batch_before = E.run_batch(folder, os.path.join(tmp, "robust_batch_1"), configs, configs, options, {"extensions": "tif"})

    progress_calls = []
    result = robust_run(
        E, folder, tmp, {"threshold_factor": [0.5, 0.8, 1.2, 1.5], "min_size": [10]},
        progress=lambda done, total, name: progress_calls.append((done, total))
    )
    out = result["run_dir"]
    t.check(os.path.basename(out).startswith("CellQuant_Robustness_"), "the result folder is named CellQuant_Robustness_<time>")
    t.equal((result["files"], result["assigned"], result["not_assigned"]), (8, 7, 1), "seven of eight files get a group")
    t.equal((result["analyzed"], result["failed"], result["skipped"]), (6, 1, 0), "six images analysed, the unreadable one reported")
    t.equal(result["specifications"], 10, "5 threshold factors x 2 minimum sizes")
    t.equal(result["groups"], ["ctrl", "mut"], "the control group is group A")
    t.equal(progress_calls[0], (0, 7), "progress starts at the first file")
    t.equal(progress_calls[-1], (7, 7), "progress ends at the last file")
    for name in ["robustness_measurements.csv", "robustness_units.csv", "robustness_groups.csv",
                 "robustness_specifications.csv", "robustness_specification_parameters.csv",
                 "robustness_parameter_values.csv", "robustness_summary.csv", "group_assignment.csv",
                 "robustness_image_report.csv", "settings.properties", "robustness_log.txt"]:
        t.check(os.path.isfile(os.path.join(out, name)), name + " is written")
    t.equal(sorted(os.listdir(folder)), before, "nothing is written into the input folder")

    assignment = robust_table(E, out, "group_assignment.csv")
    by_file = dict([(row["Image_File"], row) for row in assignment])
    t.equal(by_file["ctrl_a.tif"]["Group"], "ctrl", "the assignment table lists the group of each image")
    t.equal(by_file["Blank.tif"]["Status"], "NOT ASSIGNED", "the assignment table lists the image without a group")
    t.check("first separator" in by_file["ctrl_a.tif"]["Rule"], "the assignment table states the rule")
    report = robust_table(E, out, "robustness_image_report.csv")
    statuses = dict([(row["Image_File"], row["Status"]) for row in report])
    t.equal(statuses["mut_broken.tif"], "FAILED", "the unreadable file is reported as failed")
    t.check("Blank.tif" not in statuses, "an image without a group is not opened")

    measurements = robust_table(E, out, "robustness_measurements.csv")
    t.equal(len(measurements), 60, "one measurement per specification, image, and parent ROI")
    t.check(measurements[0]["Run_ID"] == os.path.basename(out) and measurements[0]["Application_Version"] == E.APP_VERSION,
            "every row carries the run and the application version")
    reference = [row for row in measurements if row["Is_Reference"] == "True"]
    t.equal(len(reference), 6, "one reference measurement per image")
    counts = dict([(row["Image_File"], int(row["Object_Count"])) for row in reference])
    t.equal([counts[image[0]] for image in images], [3, 4, 5, 7, 8, 9], "the reference counts are the disks drawn")
    t.close(float(reference[0]["ROI_Area"]), ROBUST_AREA, "the whole image is 12000 square micrometers")
    t.close(float(reference[0]["Value"]) * ROBUST_AREA, 3.0, "density = count / ROI area")
    t.equal(reference[0]["Metric_Unit"], "objects per " + ROBUST_UM + "^2", "the unit of the metric is stated")
    t.equal(reference[0]["ROI_Area_Unit"], ROBUST_UM + "^2", "the unit of the ROI area is stated")
    t.close(float(reference[0]["Threshold_Low"]), 1000.0, "the reference threshold is the configured one")
    low = [row for row in measurements if row["Specification"] == "threshold_factor=0.5; min_size=0.0"]
    t.close(float(low[0]["Threshold_Low"]), 500.0, "factor 0.5 halves the manual threshold")

    # The reference specification is the standard analysis: its values are
    # exactly those in the roi_summary of a batch with the same settings.
    batch_rows = robust_table(E, batch_before["batch_dir"], "batch_roi_summary.csv")
    batch_density = dict([(row["Source_Image"], row["C1_DAPI_Object_Density_Per_Area"]) for row in batch_rows])
    for row in reference:
        t.equal(row["Value"], batch_density[row["Image_File"]], "reference value equals the batch value for " + row["Image_File"])

    units = robust_table(E, out, "robustness_units.csv")
    t.equal(len(units), 60, "one unit row per specification and unit")
    groups = robust_table(E, out, "robustness_groups.csv")
    t.equal(len(groups), 20, "one group row per specification and group")
    t.close(float(groups[0]["Mean"]) * ROBUST_AREA, 4.0, "group mean of the control group")

    specifications = robust_table(E, out, "robustness_specifications.csv")
    t.equal(len(specifications), 10, "one row per specification")
    t.equal([row["Curve_Rank"] for row in specifications], [str(i + 1) for i in range(10)], "the rows are ranked for a specification curve")
    for row in specifications:
        label = row["Specification_ID"]
        t.equal((row["Group_A"], row["Group_B"], row["N_Units_A"], row["N_Units_B"]), ("ctrl", "mut", "3", "3"), label + ": three units per group")
        t.close(float(row["Mean_A"]) * ROBUST_AREA, 4.0, label + ": mean of A")
        t.close(float(row["Mean_B"]) * ROBUST_AREA, 8.0, label + ": mean of B")
        t.close(float(row["SD_A"]) * ROBUST_AREA, 1.0, label + ": SD of A")
        t.close(float(row["Difference_B_Minus_A"]) * ROBUST_AREA, 4.0, label + ": difference")
        t.close(float(row["Ratio_B_Over_A"]), 2.0, label + ": ratio")
        t.close(float(row["Hedges_g"]), 3.2, label + ": Hedges g")
        t.close(float(row["Welch_t"]), 4.0 / math.sqrt(2.0 / 3.0), label + ": Welch t")
        t.close(float(row["Welch_df"]), 4.0, label + ": Welch df")
        t.close(float(row["Mann_Whitney_P"]), 0.1, label + ": exact Mann-Whitney p")
        t.equal((row["Sign"], row["Same_Sign_As_Reference"], row["Complete"]), ("positive", "True", "True"), label + ": direction")
        if E.ROBUST_MATH["available"]:
            t.close(float(row["Welch_P"]), robust_t4_two_sided_p(4.0 / math.sqrt(2.0 / 3.0)), label + ": Welch p")
            half = ROBUST_T4_975 * math.sqrt(2.0 / 3.0)
            t.close(float(row["Difference_CI95_Low"]) * ROBUST_AREA, 4.0 - half, label + ": interval, lower end")
            t.close(float(row["Difference_CI95_High"]) * ROBUST_AREA, 4.0 + half, label + ": interval, upper end")

    summary = robust_table(E, out, "robustness_summary.csv")
    t.equal(len(summary), 1, "the summary is one row")
    summary = summary[0]
    t.equal((summary["N_Specifications"], summary["N_Specifications_Compared"], summary["N_Specifications_Incomplete"]),
            ("10", "10", "0"), "ten specifications, all compared")
    t.close(float(summary["Fraction_Same_Sign_As_Reference"]), 1.0, "the direction holds in every specification")
    t.equal((summary["N_Same_Sign_As_Reference"], summary["N_Opposite_Sign"]), ("10", "0"), "ten of ten, none reversed")
    t.equal(summary["Sign_Flipping_Values"], "none", "no value of any setting reverses the direction")
    t.close(float(summary["Hedges_g_Minimum"]), 3.2, "smallest effect size")
    t.close(float(summary["Hedges_g_Maximum"]), 3.2, "largest effect size")
    t.equal((summary["N_Images_Analyzed"], summary["N_Images_Failed"], summary["N_Images_Not_Assigned"]), ("6", "1", "1"),
            "the summary counts the images")
    t.equal(summary["Descriptive_N_Mann_Whitney_P_Below_0.05_Same_Sign"], "0", "an exact rank test of 3 against 3 cannot reach 0.05")
    if E.ROBUST_MATH["available"]:
        t.close(float(summary["Descriptive_Fraction_Welch_P_Below_0.05_Same_Sign"]), 1.0, "share of specifications with Welch p below 0.05")
    t.check("not significance tests" in summary["Statistics_Note"], "the summary says the fractions are descriptive")
    t.equal(summary["Grid"], "threshold_factor=0.5,0.8,1.0,1.2,1.5; min_size=0.0,10.0", "the summary records the grid")

    parameters = robust_table(E, out, "robustness_specification_parameters.csv")
    t.equal(len(parameters), 20, "one parameter row per specification and grid setting")
    values = robust_table(E, out, "robustness_parameter_values.csv")
    t.equal(len(values), 7, "one row per value of each grid setting")
    t.equal([row["N_Specifications"] for row in values], ["2", "2", "2", "2", "2", "5", "5"], "specifications per value")
    settings = read_utf8(os.path.join(out, "settings.properties"))
    t.check("robust.grid=threshold_factor\\=0.5,0.8,1.0,1.2,1.5; min_size\\=0.0,10.0" in settings, "the settings record holds the grid")
    t.check("channel.1.manual_low=1000.0" in settings, "the settings record holds the channel settings")
    log_text = read_utf8(os.path.join(out, "robustness_log.txt"))
    t.check("10 of 10 specifications" in log_text and "Descriptive only" in log_text, "the log states the outcome and its limits")

    # Feature off: the module adds no extension and no column, and a batch
    # after the robustness run writes exactly what it wrote before it.
    t.equal([extension.name for extension in E.RUN_EXTENSIONS if "robust" in extension.name.lower()], [], "the module adds nothing to a normal run")
    engine = E.AnalysisEngine()
    t.equal(engine.build_summary_headers(configs, E.default_global_options()), engine.build_summary_headers(configs),
            "the summary columns of a normal run are unchanged")
    batch_after = E.run_batch(folder, os.path.join(tmp, "robust_batch_2"), configs, configs, options, {"extensions": "tif"})
    for name in ["batch_roi_summary.csv", "batch_objects.csv", "batch_thresholds.csv"]:
        first = robust_without_run_names(read_utf8(os.path.join(batch_before["batch_dir"], name)))
        second = robust_without_run_names(read_utf8(os.path.join(batch_after["batch_dir"], name)))
        t.check(len(first.split("\n")) > 7, name + " of the batch holds the rows of the seven images")
        t.check(first == second, name + " of a batch is the same before and after a robustness run")



def test_robust_automatic_threshold(E, t, tmp):
    """With an automatic threshold the factor multiplies the level the
    method finds on each image, after the configured multiplier."""
    # Disks of value 2000 on a background of 0: ctrl 3, 4, 5 and mut 7, 8, 9
    # disks. An image of two grey levels leaves a threshold method no choice
    # but to separate them, so the level L it reports satisfies 0 < L <= 2000
    # (ImageJ reports the lowest value of the bright class).
    #   factor 0.5 and 0.9: 0 < 0.5 L < 0.9 L < 2000, so every disk is counted
    #   and nothing else: the counts, and with them the direction, stay.
    # The second part needs L > 200 (checked first; Otsu's level is nowhere
    # near the background), so that ten times the level is above the image:
    #   multiplier 10, factor 1.0:  10 L > 2000 -> nothing is positive (EMPTY)
    #   multiplier 10, factor 0.05: 0.5 L, between 0 and 2000 -> every disk
    images = [("ctrl_a.tif", 3, 2000), ("ctrl_b.tif", 4, 2000), ("ctrl_c.tif", 5, 2000),
              ("mut_a.tif", 7, 2000), ("mut_b.tif", 8, 2000), ("mut_c.tif", 9, 2000)]
    folder = robust_folder(tmp, "robust_auto", images, background=0)
    auto = nuclei_config(E, {"threshold_mode": "Auto", "auto_method": "Otsu"})
    result = robust_run(E, folder, tmp, {"threshold_factor": [0.5, 0.9]}, config=auto)
    rows = robust_table(E, result["run_dir"], "robustness_measurements.csv")
    one = dict([(row["Specification"], row) for row in rows if row["Image_File"] == "ctrl_a.tif"])
    base_level = float(one["threshold_factor=1.0"]["Threshold_Low"])
    t.check(200.0 < base_level <= 2000.0, "the automatic threshold separates background and disks")
    t.close(float(one["threshold_factor=0.5"]["Threshold_Low"]), 0.5 * base_level, "factor 0.5 of the automatic level")
    t.close(float(one["threshold_factor=0.9"]["Threshold_Low"]), 0.9 * base_level, "factor 0.9 of the automatic level")
    t.equal([row["Object_Count"] for row in rows if row["Image_File"] == "mut_c.tif"], ["9", "9", "9"], "the counts do not change")
    t.equal(set([row["Status"] for row in rows]), set(["OK"]), "every measurement succeeded")
    t.close(float(result["summary"]["Fraction_Same_Sign_As_Reference"]), 1.0, "the direction holds with an automatic threshold")
    t.close(result["summary"]["Reference_Difference_B_Minus_A"] * ROBUST_AREA, 4.0, "reference difference: 8 - 4 disks per image")

    # The reference specification of an automatic threshold is the standard
    # analysis: the same level and count as the engine's own run of the image.
    imp = IJ.openImage(os.path.join(folder, "ctrl_a.tif"))
    engine = E.AnalysisEngine()
    data = engine.build_masks(imp, auto, "check")
    t.close(base_level, data["threshold_low"], "the reference uses the level of a normal run")
    for key in ["raw", "final"]:
        E.safe_close_image(data[key])
    imp.close()

    # The multiplier and the factor together; a level above the image is an
    # empty result, and a factor can bring it back into the image.
    result = robust_run(E, folder, tmp, {"auto_scale": [10], "threshold_factor": [0.05]}, config=auto)
    rows = robust_table(E, result["run_dir"], "robustness_measurements.csv")
    t.equal(len(rows), 24, "2 multipliers x 2 factors x 6 images")
    by_spec = {}
    for row in rows:
        by_spec.setdefault(row["Specification"], []).append(row)
    high = by_spec["threshold_factor=1.0; auto_scale=10.0"]
    t.equal(set([row["Status"] for row in high]), set(["EMPTY"]), "ten times the automatic level lies above the image: an empty result")
    t.equal(set([row["Object_Count"] for row in high]), set(["0"]), "an empty result counts zero objects")
    back = dict([(row["Image_File"], row) for row in by_spec["threshold_factor=0.05; auto_scale=10.0"]])
    t.close(float(back["ctrl_a.tif"]["Threshold_Low"]), 0.5 * base_level, "level x multiplier 10 x factor 0.05 = half the level")
    t.equal([back[image[0]]["Object_Count"] for image in images], ["3", "4", "5", "7", "8", "9"],
            "with the level back inside the image every disk is counted")
    t.equal(set([row["Status"] for row in back.values()]), set(["OK"]), "and the measurement is a normal one")
    specifications = dict([(row["Specification"], row) for row in result["tables"]["specifications"]])
    t.equal(specifications["threshold_factor=1.0; auto_scale=10.0"]["Sign"], "zero", "both groups are empty: no difference")
    t.equal((result["summary"]["N_Same_Sign_As_Reference"], result["summary"]["N_Zero_Or_Undefined"],
             result["summary"]["N_Opposite_Sign"]), (3, 1, 0), "three specifications keep the direction, one is zero")
    t.equal(result["summary"]["Sign_Flipping_Values"], "none", "a zero difference is not a reversal")


def test_robust_sign_flip(E, t, tmp):
    """Group A has many dim objects, group B few bright ones: the direction
    of the difference reverses above a threshold known by construction."""
    #   ctrl: 9, 10, 11 disks of value 1200 -> mean 10, SD 1 (per 12000 um^2)
    #   mut:  4, 5, 6 disks of value 3000   -> mean 5, SD 1
    # Manual threshold 1000 x factor: 500, 800, 1000, 1100 keep every disk
    # (difference 5 - 10 = -5, negative); 1300, 1500, 2000 lie above the dim
    # disks (and above every pixel of a control image, which is then empty)
    # but below the bright ones (difference 5 - 0 = +5, positive).
    # So 4 of 7 specifications have the sign of the reference, and the
    # factors 1.3, 1.5 and 2.0 reverse it.
    # Hedges g: below the flip d = -5 / 1, g = -4; above it the control SD is
    # 0, pooled SD = sqrt((0 + 2) / 4) = 0.7071, d = 7.0711, g = 5.6569.
    images = [("ctrl_a.tif", 9, 1200), ("ctrl_b.tif", 10, 1200), ("ctrl_c.tif", 11, 1200),
              ("mut_a.tif", 4, 3000), ("mut_b.tif", 5, 3000), ("mut_c.tif", 6, 3000)]
    folder = robust_folder(tmp, "robust_flip", images)
    factors = [0.5, 0.8, 1.1, 1.3, 1.5, 2.0]
    result = robust_run(E, folder, tmp, {"threshold_factor": factors})
    out = result["run_dir"]
    rows = dict([(row["Specification"], row) for row in robust_table(E, out, "robustness_specifications.csv")])
    t.equal(len(rows), 7, "six factors and the configured one")
    for factor in [0.5, 0.8, 1.0, 1.1]:
        row = rows["threshold_factor=" + repr(factor)]
        t.close(float(row["Difference_B_Minus_A"]) * ROBUST_AREA, -5.0, "factor " + repr(factor) + ": B has 5 fewer")
        t.close(float(row["Hedges_g"]), -4.0, "factor " + repr(factor) + ": g")
        t.equal((row["Sign"], row["Same_Sign_As_Reference"]), ("negative", "True"), "factor " + repr(factor) + ": direction of the reference")
    for factor in [1.3, 1.5, 2.0]:
        row = rows["threshold_factor=" + repr(factor)]
        t.close(float(row["Mean_A"]), 0.0, "factor " + repr(factor) + ": no dim disk is left")
        t.close(float(row["Difference_B_Minus_A"]) * ROBUST_AREA, 5.0, "factor " + repr(factor) + ": B has 5 more")
        t.close(float(row["Hedges_g"]), 0.8 * 5.0 / math.sqrt(0.5), "factor " + repr(factor) + ": g")
        t.equal((row["Sign"], row["Same_Sign_As_Reference"]), ("positive", "False"), "factor " + repr(factor) + ": reversed")
    ordered = robust_table(E, out, "robustness_specifications.csv")
    t.equal([row["Sign"] for row in ordered], ["negative"] * 4 + ["positive"] * 3, "the curve is sorted by effect size")
    t.equal(rows["threshold_factor=1.0"]["Is_Reference"], "True", "factor 1.0 is the reference")

    summary = result["summary"]
    t.equal((summary["N_Same_Sign_As_Reference"], summary["N_Opposite_Sign"], summary["N_Specifications_Compared"]), (4, 3, 7),
            "four specifications keep the direction, three reverse it")
    t.close(summary["Fraction_Same_Sign_As_Reference"], 4.0 / 7.0, "fraction with the sign of the reference")
    t.equal(summary["Reference_Sign"], "negative", "the reference direction")
    t.equal(summary["Sign_Flipping_Values"],
            "threshold_factor=1.3 (1 of 1); threshold_factor=1.5 (1 of 1); threshold_factor=2.0 (1 of 1)",
            "the summary names the values that reverse the direction")
    t.close(summary["Hedges_g_Minimum"], -4.0, "smallest g")
    t.close(summary["Hedges_g_Median"], -4.0, "median g (four values of -4, three of 5.66)")
    t.close(summary["Hedges_g_Maximum"], 0.8 * 5.0 / math.sqrt(0.5), "largest g")
    values = dict([(row["Value"], row) for row in robust_table(E, out, "robustness_parameter_values.csv")])
    t.equal((values["1.1"]["N_Opposite_Sign"], values["1.3"]["N_Opposite_Sign"]), ("0", "1"), "the flip lies between factor 1.1 and 1.3")
    t.equal(values["1.0"]["Is_Reference_Value"], "True", "the configured value is marked")

    measurements = robust_table(E, out, "robustness_measurements.csv")
    empty = [row for row in measurements if row["Status"] == "EMPTY"]
    t.equal(len(empty), 9, "a threshold above every pixel of an image is an empty result (3 factors x 3 control images)")
    t.equal(set([row["Object_Count"] for row in empty]), set(["0"]), "an empty result counts zero objects")
    t.check("nothing is positive" in empty[0]["Note"], "and says why")

    # The engine reports an empty threshold range either with a fixed error
    # text or, in a later engine, with an empty mask marked "empty_range".
    # The module recognises both; keep the engine within one of them.
    probe = ImagePlus("probe", robust_plane(3, 1200))
    message = None
    flagged = False
    try:
        data = E.AnalysisEngine().build_masks(probe, nuclei_config(E, {"manual_low": 5000}), "probe")
        flagged = data.get("empty_range", False) is True
        for key in ["raw", "final"]:
            E.safe_close_image(data[key])
    except ValueError as exc:
        message = E.to_text(exc)
    t.check(flagged or (message is not None and E.ROBUST_EMPTY_RANGE_TEXT in message),
            "the engine reports an empty threshold range in a way the module recognises")
    scaled = nuclei_config(E, {"threshold_mode": "Auto", "auto_method": "Otsu"})
    result = robust_run(E, folder, tmp, {"auto_scale": [10]}, config=scaled)
    rows = robust_table(E, result["run_dir"], "robustness_measurements.csv")
    high = [row for row in rows if row["Specification"] == "auto_scale=10.0"]
    t.equal(set([row["Status"] for row in high]), set(["EMPTY"]), "an automatic threshold scaled above the image is an empty result too")

    # Mean object area: every disk has the same area, so the groups do not
    # differ at the reference (sign zero); above the flip the control group
    # has no object and therefore no value.
    disk_area = disk_pixel_count() * PIXEL * PIXEL
    result = robust_run(E, folder, tmp, {"threshold_factor": [1.5]}, metric=E.ROBUST_METRIC_OBJECT_AREA)
    rows = dict([(row["Specification"], row) for row in result["tables"]["specifications"]])
    t.close(rows["threshold_factor=1.0"]["Mean_A"], disk_area, "mean object area is the area of a disk")
    t.equal(rows["threshold_factor=1.0"]["Sign"], "zero", "equal object areas give no difference")
    t.equal(rows["threshold_factor=1.5"]["N_Units_A"], 0, "without objects a unit has no mean object area")
    t.equal(rows["threshold_factor=1.5"]["Difference_B_Minus_A"], "", "and the specification has no difference")
    t.check("zero at the reference" in result["summary"]["Note"], "a zero reference difference is explained")
    t.equal(result["summary"]["Metric_Unit"], ROBUST_UM + "^2", "the unit of a mean object area")
    undefined = [row for row in robust_table(E, result["run_dir"], "robustness_measurements.csv") if row["Status"] == "UNDEFINED"]
    t.equal(len(undefined), 3, "the three control images without objects have an undefined mean object area")


def robust_l_shape():
    """An L-shaped parent ROI that touches the image border: the top band
    y < 100 over the full width plus the left half below it.
    Area = 240 * 100 + 120 * 100 = 36000 pixels = 9000 square micrometers."""
    polygon = Polygon()
    for point in [(0, 0), (240, 0), (240, 100), (120, 100), (120, 200), (0, 200)]:
        polygon.addPoint(point[0], point[1])
    roi = PolygonRoi(polygon, Roi.POLYGON)
    roi.setName("L_shape")
    return roi


def test_robust_units_and_rois(E, t, tmp):
    """Several images of one organoid count as one unit; parent ROIs come
    from ROI files, including a non-rectangular one on the image border."""
    # ctrl_o1_a: 3 disks, ROI = left half (120 x 200 px = 6000 um^2), which holds all 3 -> 3/6000
    # ctrl_o1_b: 3 disks, whole image                                                    -> 3/12000
    # ctrl_o2_a: 5 disks -> 5/12000;  ctrl_o3_a: 7 disks -> 7/12000
    # mut_o1_a: the first 14 disks plus one centred at (3, 190), which the image
    #           border cuts. ROI = the L shape (9000 um^2): it holds the 10 disks
    #           of the two top rows, (50,140) and (100,140) of the third row, and
    #           the cut disk = 13 objects; (150,140) and (200,140) lie outside -> 13/9000
    # mut_o1_b: 9 disks -> 9/12000;  mut_o2_a: 10 disks -> 10/12000
    # mut_o3_a: no disk -> 0;  mut_o3_b: no disk, ROI of 2 x 2 pixels -> 0/1
    # Units (organoids) by mean of their measurements:
    #   ctrl: o1 = (3/6000 + 3/12000)/2 = 4.5/12000, o2 = 5/12000, o3 = 7/12000 -> mean 5.5/12000
    #   mut:  o1 = (13/9000 + 9/12000)/2, o2 = 10/12000, o3 = 0
    # Pooled instead: ctrl o1 = 6/18000, mut o1 = 22/21000, mut o3 = 0/12001.
    folder = os.path.join(tmp, "robust_units")
    os.makedirs(folder)
    for name, count in [("ctrl_o1_a.tif", 3), ("ctrl_o1_b.tif", 3), ("ctrl_o2_a.tif", 5), ("ctrl_o3_a.tif", 7),
                        ("mut_o1_b.tif", 9), ("mut_o2_a.tif", 10), ("mut_o3_a.tif", 0), ("mut_o3_b.tif", 0)]:
        robust_save(folder, name, robust_plane(count, 2000))
    robust_save(folder, "mut_o1_a.tif", robust_plane(14, 2000, extra_disks=[(3, 190)]))
    E.save_roi_array([Roi(0, 0, 120, 200)], os.path.join(folder, "ctrl_o1_a_ROIs.zip"))
    E.save_roi_array([robust_l_shape()], os.path.join(folder, "mut_o1_a_ROIs.zip"))
    E.save_roi_array([Roi(0, 0, 2, 2)], os.path.join(folder, "mut_o3_b_ROIs.zip"))

    groups = {"source": E.ROBUST_GROUP_REGEX, "group_regex": "^([a-z]+)_", "unit_regex": "^([a-z]+_o[0-9]+)_"}
    batch_options = {"extensions": "tif", "roi_policy": E.BATCH_ROI_FILE_WHOLE}
    result = robust_run(E, folder, tmp, {"min_size": [10]}, groups=groups, batch_options=batch_options)
    t.equal(result["analyzed"], 9, "nine images analysed")
    measurements = robust_table(E, result["run_dir"], "robustness_measurements.csv")
    reference = dict([(row["Image_File"], row) for row in measurements if row["Is_Reference"] == "True"])
    t.close(float(reference["mut_o1_a.tif"]["ROI_Area"]), 9000.0, "the L-shaped ROI has its exact area")
    t.equal(reference["mut_o1_a.tif"]["Object_Count"], "13", "objects are counted inside the true outline of the ROI")
    t.equal(reference["mut_o1_a.tif"]["Parent_ROI_Name"], "L_shape", "the ROI of the ROI file is used")
    t.close(float(reference["ctrl_o1_a.tif"]["ROI_Area"]), 6000.0, "the rectangular ROI has its area")
    t.close(float(reference["mut_o3_b.tif"]["ROI_Area"]), 1.0, "a 2 x 2 pixel ROI is one square micrometer")
    t.equal((reference["mut_o3_b.tif"]["Object_Count"], float(reference["mut_o3_b.tif"]["Value"])), ("0", 0.0), "a tiny empty ROI has density zero")
    # mut_o3_a holds no pixel above the threshold. A normal run stops on such
    # an image; here it is a measurement of zero, marked EMPTY, because
    # leaving out exactly the images without signal would bias the groups.
    t.equal((reference["mut_o3_a.tif"]["Status"], reference["mut_o3_a.tif"]["Object_Count"], float(reference["mut_o3_a.tif"]["Value"])),
            ("EMPTY", "0", 0.0), "an image without any pixel above the threshold is a measurement of zero")

    units = [row for row in result["tables"]["units"] if row["Is_Reference"]]
    by_unit = dict([(row["Unit"], row) for row in units])
    t.equal(sorted(by_unit.keys()), ["ctrl_o1", "ctrl_o2", "ctrl_o3", "mut_o1", "mut_o2", "mut_o3"], "the units are the organoids")
    t.equal((by_unit["ctrl_o1"]["N_Images"], by_unit["ctrl_o1"]["N_Measurements"]), (2, 2), "two images of one organoid are one unit")
    t.close(by_unit["ctrl_o1"]["Value"] * 12000.0, 4.5, "a unit is the mean of its measurements")
    t.close(by_unit["mut_o1"]["Value"] * 36000.0, (13.0 / 9000.0 + 9.0 / 12000.0) / 2.0 * 36000.0, "unit with the L-shaped ROI")
    t.close(by_unit["mut_o3"]["Value"], 0.0, "a unit without objects has density zero")
    rows = dict([(row["Specification"], row) for row in result["tables"]["specifications"]])
    row = rows["min_size=0.0"]
    t.equal((row["N_Units_A"], row["N_Units_B"]), (3, 3), "three units per group, not four and five images")
    t.close(row["Mean_A"] * 12000.0, 5.5, "group mean over units (over images it would be 5.25)")
    expected_b = ((13.0 / 9000.0 + 9.0 / 12000.0) / 2.0 + 10.0 / 12000.0 + 0.0) / 3.0
    t.close(row["Mean_B"] * 36000.0, expected_b * 36000.0, "group mean of B over units")
    t.equal(rows["min_size=10.0"]["Mean_B"], row["Mean_B"], "a 10 um^2 minimum keeps every disk, also the one cut by the border")

    pooled = robust_run(E, folder, tmp, {"min_size": [10]}, groups=groups, batch_options=batch_options,
                        unit_aggregation=E.ROBUST_UNIT_POOLED)
    by_unit = dict([(row["Unit"], row) for row in pooled["tables"]["units"] if row["Is_Reference"]])
    t.close(by_unit["ctrl_o1"]["Value"] * 18000.0, 6.0, "pooled unit: 6 objects in 18000 square micrometers")
    t.close(by_unit["mut_o1"]["Value"] * 21000.0, 22.0, "pooled unit: 22 objects in 21000 square micrometers")
    t.close(by_unit["mut_o3"]["Value"], 0.0, "pooled unit without objects")
    t.equal(pooled["summary"]["Unit_Aggregation"], E.ROBUST_UNIT_POOLED, "the summary states how units were combined")

    # The same assignment from a table gives the same result; an image the
    # table does not list is left out, and so is a row without an image.
    table = os.path.join(tmp, "robust_units.csv")
    handle = codecs.open(table, "w", "utf-8")
    handle.write(u"Image,Group,Unit\n")
    for name in sorted(reference.keys()):
        if name != "mut_o3_b.tif":
            handle.write(name + u"," + name.split("_")[0] + u"," + name[:-6] + u"\n")
    handle.write(u"missing.tif,ctrl,ctrl_o9\n")
    handle.close()
    from_table = robust_run(E, folder, tmp, {"min_size": [10]}, batch_options=batch_options,
                            groups={"source": E.ROBUST_GROUP_CSV, "csv_path": table, "control_group": "ctrl"})
    t.equal((from_table["assigned"], from_table["not_assigned"], from_table["analyzed"]), (8, 1, 8), "the unlisted image is not analysed")
    table_rows = dict([(row["Specification"], row) for row in from_table["tables"]["specifications"]])
    t.equal(table_rows["min_size=0.0"]["Mean_A"], row["Mean_A"], "group A from the table equals group A from the expression")
    # without mut_o3_b the unit mut_o3 is the single empty image: still 0
    t.equal(table_rows["min_size=0.0"]["Mean_B"], row["Mean_B"], "group B from the table equals group B from the expression")
    assignment = robust_table(E, from_table["run_dir"], "group_assignment.csv")
    statuses = dict([(row["Image_Relative_Path"], row["Status"]) for row in assignment])
    t.equal(statuses["mut_o3_b.tif"], "NOT ASSIGNED", "the assignment table shows the unlisted image")
    t.equal(statuses["missing.tif"], "NO IMAGE FOR TABLE ROW", "the assignment table shows the row without an image")

    # ROI files required: images without one are skipped and reported.
    strict = dict(batch_options)
    strict["roi_policy"] = E.BATCH_ROI_FILE_SKIP
    skipped = robust_run(E, folder, tmp, {"min_size": [10]}, groups=groups, batch_options=strict)
    t.equal((skipped["analyzed"], skipped["skipped"]), (3, 6), "images without a ROI file are skipped as in a batch")
    t.equal(skipped["summary"]["N_Images_Skipped"], 6, "and counted in the summary")
    t.equal(skipped["tables"]["specifications"][0]["N_Units_A"], 1, "one control unit is left")
    t.check("Fewer than 3 units" in skipped["tables"]["specifications"][0]["Note"], "too few units: means only, with the reason")


def test_robust_small_designs(E, t, tmp):
    """Fewer than three units, more than two groups, one group, cancel."""
    images = [("ctrl_a.tif", 2, 2000), ("ctrl_b.tif", 4, 2000), ("mut_a.tif", 6, 2000), ("mut_b.tif", 8, 2000),
              ("het_a.tif", 5, 2000), ("het_b.tif", 5, 2000)]
    folder = robust_folder(tmp, "robust_small", images)
    two_groups = [os.path.join(folder, image[0]) for image in images[:4]]

    # Two units per group: means 3 and 7 per 12000 um^2, difference 4, ratio 7/3.
    result = robust_run(E, folder, tmp, {"threshold_factor": [1.2]}, file_list=two_groups)
    row = result["tables"]["specifications"][0]
    t.equal((row["N_Units_A"], row["N_Units_B"]), (2, 2), "two units per group")
    t.close(row["Mean_A"] * ROBUST_AREA, 3.0, "mean of A is reported")
    t.close(row["Difference_B_Minus_A"] * ROBUST_AREA, 4.0, "the difference of the means is reported")
    t.close(row["Ratio_B_Over_A"], 7.0 / 3.0, "the ratio of the means is reported")
    t.equal((row["Hedges_g"], row["Welch_P"], row["Mann_Whitney_P"], row["Difference_CI95_Low"]), ("", "", "", ""),
            "no effect size, test, or interval with two units")
    t.check("Fewer than 3 units" in row["Note"], "the note says why")
    t.close(result["summary"]["Fraction_Same_Sign_As_Reference"], 1.0, "the direction is still compared")
    t.equal(result["summary"]["Hedges_g_Minimum"], "", "the summary has no effect-size range")
    # No test was run, so the descriptive shares of p below 0.05 are empty,
    # not 0 (which would read as "tested and never below 0.05").
    for name in ["Welch", "Mann_Whitney"]:
        for prefix in ["Descriptive_N_", "Descriptive_Fraction_"]:
            t.equal(result["summary"].get(prefix + name + "_P_Below_0.05_Same_Sign", ""), "",
                    prefix + name + " is empty when no specification has a p value")
    written = robust_table(E, result["run_dir"], "robustness_summary.csv")[0]
    t.equal(written["Descriptive_Fraction_Welch_P_Below_0.05_Same_Sign"], "", "and the written summary cell is empty")
    t.equal(result["summary"]["Curve_Sorted_By"], "Difference_B_Minus_A", "without effect sizes the curve is sorted by the difference")

    # Three groups: group means only.
    result = robust_run(E, folder, tmp, {"threshold_factor": [1.2]}, groups={})
    t.equal(result["groups"], ["ctrl", "het", "mut"], "three groups are found")
    means = dict([(row["Group"], row["Mean"]) for row in result["tables"]["groups"] if row["Is_Reference"]])
    t.close(means["het"] * ROBUST_AREA, 5.0, "group means are reported for every group")
    row = result["tables"]["specifications"][0]
    t.equal((row["N_Groups"], row.get("Difference_B_Minus_A", "")), (3, ""), "no comparison with three groups")
    t.check("More than two groups" in row["Note"], "the specification row says why")
    t.check("exactly two" in result["summary"]["Note"], "the summary says why")

    # One group only: refused before anything is written.
    message = robust_refused(lambda: robust_run(E, folder, tmp, {"threshold_factor": [1.2]}, file_list=two_groups[:2]))
    t.check(message is not None and "at least two groups" in message, "a single group is refused")
    message = robust_refused(lambda: robust_run(E, folder, tmp, {"threshold_factor": [1.2]}, groups={"control_group": "wt"}))
    t.check(message is not None and "not among the groups" in message, "an unknown control group is refused")
    message = robust_refused(lambda: robust_run(E, os.path.join(folder, "nowhere"), tmp, {"threshold_factor": [1.2]}))
    t.check(message is not None and "does not exist" in message, "a missing input folder is refused")

    # Cancel after the first image: its measurements stay, no statistics.
    calls = {"n": 0}

    def cancel_after_first():
        calls["n"] += 1
        return calls["n"] > 3

    result = robust_run(E, folder, tmp, {"threshold_factor": [1.2]}, file_list=two_groups, is_cancelled=cancel_after_first)
    t.check(result["cancelled"], "the run reports that it was cancelled")
    t.check(result["summary"] is None, "a cancelled run returns no summary")
    t.check(not os.path.exists(os.path.join(result["run_dir"], "robustness_summary.csv")), "a cancelled run writes no summary")
    t.check(not os.path.exists(os.path.join(result["run_dir"], "robustness_specifications.csv")), "a cancelled run writes no statistics")
    t.check(os.path.isfile(os.path.join(result["run_dir"], "group_assignment.csv")), "the group assignment is written before the first image")
    t.check("Cancelled by user" in read_utf8(os.path.join(result["run_dir"], "robustness_log.txt")), "the log records the cancellation")
    t.check(result["analyzed"] < 4, "not every image was analysed")

    # Settings stored as text (a profile or the window) run the same analysis.
    options = E.default_global_options()
    options.update({"robust_grid": "threshold_factor=1.2", "robust_control_group": "ctrl"})
    configs = [nuclei_config(E)]
    result = E.run_robustness_from_options(folder, os.path.join(tmp, "robust_out"), configs, configs, options,
                                           {"extensions": "tif"}, file_list=two_groups)
    t.close(result["tables"]["specifications"][0]["Mean_B"] * ROBUST_AREA, 7.0, "run_robustness_from_options reads its arguments from the settings")


def test_robust_area_and_preparation(E, t, tmp):
    """Positive-area fraction with a known flip, and a preparation setting
    (median filter) in the grid."""
    # Positive-area fraction. Each image holds bright rectangles on 100:
    #   ctrl: a w x 80 rectangle of value 900 with a 60 x 40 core of value 1500, w = 96, 100, 104
    #   mut:  a w x 60 rectangle of value 1500, w = 76, 80, 84
    # Manual threshold 800 x factor over 48000 pixels:
    #   400, 800: ctrl = 7680, 8000, 8320 px (mean 16.6667 percent); mut = 4560, 4800, 5040 px (mean 10) -> negative
    #   1200:     ctrl = 2400 px each (5 percent, SD 0);             mut mean 10 percent                  -> positive
    #   2000:     above every pixel: 0 for both                                                         -> zero
    # So of 4 compared specifications 2 share the reference sign, 1 is reversed, 1 is zero.
    folder = os.path.join(tmp, "robust_area")
    os.makedirs(folder)
    for index, width in enumerate([96, 100, 104]):
        ip = blank_plane(100)
        draw_rect(ip, (20, 20, width, 80), 900)
        draw_rect(ip, (30, 30, 60, 40), 1500)
        robust_save(folder, "ctrl_" + str(index) + ".tif", ip)
    for index, width in enumerate([76, 80, 84]):
        ip = blank_plane(100)
        draw_rect(ip, (20, 20, width, 60), 1500)
        robust_save(folder, "mut_" + str(index) + ".tif", ip)
    marker = E.normalize_channel_config(1, {
        "enabled": True, "name": "MAP2", "analysis_mode": "Positive area",
        "threshold_mode": "Manual", "manual_low": 800, "manual_high": -1
    })
    result = robust_run(E, folder, tmp, {"threshold_factor": [0.5, 1.5, 2.5]}, config=marker,
                        metric=E.ROBUST_METRIC_AREA_FRACTION)
    rows = dict([(row["Specification"], row) for row in result["tables"]["specifications"]])
    for factor in ["0.5", "1.0"]:
        t.close(rows["threshold_factor=" + factor]["Mean_A"], 8000.0 / 480.0, "factor " + factor + ": the control mean is 16.67 percent")
        t.close(rows["threshold_factor=" + factor]["Mean_B"], 10.0, "factor " + factor + ": the mutant mean is 10 percent")
        t.equal(rows["threshold_factor=" + factor]["Sign"], "negative", "factor " + factor + ": B below A")
    t.close(rows["threshold_factor=1.5"]["Mean_A"], 5.0, "above 900 only the core is positive: 5 percent")
    t.close(rows["threshold_factor=1.5"]["SD_A"], 0.0, "the core has the same size in every control image")
    t.equal(rows["threshold_factor=1.5"]["Sign"], "positive", "above 900 the direction is reversed")
    t.equal(rows["threshold_factor=2.5"]["Sign"], "zero", "above every pixel both groups are zero")
    summary = result["summary"]
    t.equal((summary["N_Same_Sign_As_Reference"], summary["N_Opposite_Sign"], summary["N_Zero_Or_Undefined"]), (2, 1, 1),
            "two keep the direction, one reverses it, one is zero")
    t.close(summary["Fraction_Same_Sign_As_Reference"], 0.5, "fraction with the sign of the reference")
    t.equal(summary["Sign_Flipping_Values"], "threshold_factor=1.5 (1 of 1)", "the reversing value is named")
    t.equal(summary["Metric_Unit"], "percent", "the unit of the area fraction")
    measurement = [row for row in robust_table(E, result["run_dir"], "robustness_measurements.csv")
                   if row["Is_Reference"] == "True" and row["Image_File"] == "ctrl_1.tif"][0]
    t.close(float(measurement["Positive_Area"]), 8000 * PIXEL * PIXEL, "positive area of the reference")
    t.close(float(measurement["Value"]), 8000.0 / 480.0, "positive-area fraction in percent")

    # A preparation setting. ctrl images hold 3, 4, 5 disks plus 6 single
    # bright pixels; mut images hold 6, 7, 8 disks. Unfiltered, each single
    # pixel is an object: ctrl = 9, 10, 11 against mut = 6, 7, 8 (negative).
    # A median filter of radius 1 (a 3 x 3 neighbourhood) removes a single
    # pixel and keeps a disk of radius 6: ctrl = 3, 4, 5 (positive).
    folder = os.path.join(tmp, "robust_median")
    os.makedirs(folder)
    salt = [(20, 110), (60, 110), (100, 110), (140, 110), (180, 110), (220, 110)]
    for index, count in enumerate([3, 4, 5]):
        robust_save(folder, "ctrl_" + str(index) + ".tif", robust_plane(count, 2000, salt=salt))
    for index, count in enumerate([6, 7, 8]):
        robust_save(folder, "mut_" + str(index) + ".tif", robust_plane(count, 2000))
    result = robust_run(E, folder, tmp, {"median_radius": [1], "threshold_factor": [0.8]})
    rows = dict([(row["Specification"], row) for row in result["tables"]["specifications"]])
    t.equal(len(rows), 4, "2 median radii x 2 threshold factors")
    for factor in ["0.8", "1.0"]:
        raw = rows["threshold_factor=" + factor + "; median_radius=0.0"]
        filtered = rows["threshold_factor=" + factor + "; median_radius=1.0"]
        t.close(raw["Mean_A"] * ROBUST_AREA, 10.0, "factor " + factor + ", unfiltered: single pixels are counted")
        t.close(filtered["Mean_A"] * ROBUST_AREA, 4.0, "factor " + factor + ", median filter: only the disks are left")
        t.close(filtered["Mean_B"] * ROBUST_AREA, 7.0, "factor " + factor + ", median filter: the disks of B are kept")
        t.equal((raw["Sign"], filtered["Sign"]), ("negative", "positive"), "factor " + factor + ": the filter reverses the direction")
    t.equal(result["summary"]["Sign_Flipping_Values"], "median_radius=1.0 (2 of 2)", "the preparation setting is named as the cause")


def test_robust_fast_count(E, t, tmp):
    """The quick object count used here gives the count and total area of
    analyze_objects, for every option that changes what is counted."""
    # Disks of radius 6; (40,40) and (50,40) overlap and form one object
    # that a watershed splits; (5,100) is cut by the image border.
    ip = ByteProcessor(WIDTH, HEIGHT)
    ip.setValue(255.0)
    for center in DISKS + [(40, 40), (5, 100)]:
        ip.fill(OvalRoi(center[0] - DISK_RADIUS, center[1] - DISK_RADIUS, 2 * DISK_RADIUS + 1, 2 * DISK_RADIUS + 1))
    ip.resetRoi()
    ip.setValue(0.0)
    ip.fill(Roi(109, 79, 3, 3))    # a hole in the disk at (110, 80)
    ip.resetRoi()
    mask = ImagePlus("mask", ip)
    cal = mask.getCalibration()
    cal.pixelWidth = PIXEL
    cal.pixelHeight = PIXEL
    cal.setUnit("um")
    engine = E.AnalysisEngine()
    rois = [Roi(0, 0, WIDTH, HEIGHT), triangle_roi(), robust_l_shape(), Roi(0, 0, 2, 2)]
    variants = [{}, {"min_size": 20}, {"max_size": 30}, {"min_circularity": 0.9}, {"include_holes": True},
                {"four_connected": True}, {"size_units": "Pixel^2", "min_size": 120}, {"exclude_edges": True}]
    for overrides in variants:
        config = nuclei_config(E, overrides)
        for roi in rois:
            masks = E.roi_masks(roi, WIDTH, HEIGHT)
            objects, info = engine.analyze_objects(mask, roi, config, masks)
            expected = (len(objects), sum([obj["area"] for obj in objects]))
            got = E.robust_count_objects(engine, mask, roi, config, lambda: masks)
            label = repr(sorted(overrides.items())) + " in " + E.to_text(roi.getName() or roi.getBounds())
            t.equal(got[0], expected[0], "count as analyze_objects: " + label)
            t.close(got[1], expected[1], "total area as analyze_objects: " + label, 1e-9)
    whole = E.robust_count_objects(engine, mask, rois[0], nuclei_config(E), lambda: E.roi_masks(rois[0], WIDTH, HEIGHT))
    t.equal(whole[0], len(DISKS) + 1, "15 separate disks, the overlapping pair as one object, the cut disk")
    empty = ImagePlus("empty", ByteProcessor(WIDTH, HEIGHT))
    t.equal(E.robust_count_objects(engine, empty, rois[1], nuclei_config(E), lambda: E.roi_masks(rois[1], WIDTH, HEIGHT)),
            (0, 0.0), "an empty mask has no objects and no area")
    t.equal(E.robust_metric_value(E.ROBUST_METRIC_OBJECT_AREA, 100.0, 0, 0.0, None)[0], None, "no objects: the mean object area is undefined")
    t.close(E.robust_metric_value(E.ROBUST_METRIC_DENSITY, 100.0, 5, 0.0, None)[0], 0.05, "density = count / area")
    t.close(E.robust_metric_value(E.ROBUST_METRIC_AREA_FRACTION, 200.0, None, None, 50.0)[0], 25.0, "fraction = positive / area x 100")


# -----------------------------------------------------------------------------
# Known-answer cases added by the independent check
# -----------------------------------------------------------------------------

def robust_log_gamma(x):
    """log Gamma(x) for x > 0 (Lanczos, g = 7, nine terms; about 15 digits)."""
    coefficients = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
                    -176.61502916214059, 12.507343278686905, -0.13857109526572012,
                    9.9843695780195716e-6, 1.5056327351493116e-7]
    if x < 0.5:
        return math.log(math.pi / abs(math.sin(math.pi * x))) - robust_log_gamma(1.0 - x)
    x -= 1.0
    total = coefficients[0]
    for index in range(1, 9):
        total += coefficients[index] / (x + index)
    shifted = x + 7.5
    return 0.5 * math.log(2.0 * math.pi) + (x + 0.5) * math.log(shifted) - shifted + math.log(total)


def robust_incomplete_beta(a, b, x):
    """Regularized incomplete beta function I_x(a, b) by its continued
    fraction (modified Lentz). Written here so that the Welch p value and
    interval are checked against code that shares nothing with Commons Math."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    if x > (a + 1.0) / (a + b + 2.0):
        return 1.0 - robust_incomplete_beta(b, a, 1.0 - x)
    front = math.exp(robust_log_gamma(a + b) - robust_log_gamma(a) - robust_log_gamma(b) +
                     a * math.log(x) + b * math.log(1.0 - x)) / a
    tiny = 1e-300
    c = 1.0
    d = 1.0 - (a + b) * x / (a + 1.0)
    d = 1.0 / (d if abs(d) > tiny else tiny)
    value = d
    for m in range(1, 300):
        for numerator in [m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m)),
                          -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0))]:
            d = 1.0 + numerator * d
            d = 1.0 / (d if abs(d) > tiny else tiny)
            c = 1.0 + numerator / c
            c = c if abs(c) > tiny else tiny
            value *= c * d
        if abs(c * d - 1.0) < 1e-15:
            break
    return front * value


def robust_t_two_sided(t_value, df):
    """Two-sided p of Student's t: I_{df/(df+t^2)}(df/2, 1/2)."""
    return robust_incomplete_beta(df / 2.0, 0.5, df / (df + t_value * t_value))


def robust_t_quantile_975(df):
    """The t with two-sided p = 0.05, by bisection."""
    low, high = 0.0, 1000.0
    for step in range(200):
        middle = (low + high) / 2.0
        if robust_t_two_sided(middle, df) > 0.05:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def robust_byte_plane(disks, background, value):
    ip = ByteProcessor(WIDTH, HEIGHT)
    ip.setValue(float(background))
    ip.fill()
    for disk in disks:
        draw_disk(ip, disk[0], disk[1], DISK_RADIUS, value)
    return ip


def robust_notched_roi():
    """Rectilinear, non-rectangular parent ROI on the right and bottom image
    border: the rectangle x 120..240, y 110..200 without its lower left
    corner x 120..160, y 160..200. 120*90 - 40*40 = 9200 pixels."""
    polygon = Polygon()
    for point in [(120, 110), (240, 110), (240, 200), (160, 200), (160, 160), (120, 160)]:
        polygon.addPoint(point[0], point[1])
    roi = PolygonRoi(polygon, Roi.POLYGON)
    roi.setName("notched")
    return roi


def test_robust_micron_byte_rois(E, t, tmp):
    """8-bit images, a calibration override in "micron", two parent ROIs from
    the ROI Manager set (one non-rectangular), edge exclusion, and an empty
    threshold, with every value derived by hand."""
    # Calibration override 0.25 micron per pixel: a pixel is 0.0625 micron^2.
    # R1 = rectangle x 0..120, y 0..110 (top-left corner of the image):
    #      13200 px = 825 micron^2.
    # R2 = the notched polygon of robust_notched_roi: 9200 px = 575 micron^2.
    # Disks of radius 6, value 200 on 10, threshold Manual 100:
    #   k disks from P1, each wholly inside R1;  m disks from P2, wholly inside R2;
    #   in every image also: E3 = (3, 50) cut by the left image border, inside R1;
    #   E1 = (236, 130) cut by the right image border, inside R2; E2 = (190, 108)
    #   crossing the top outline of R2 (only its rows y >= 110 are in R2); and
    #   (60, 160), which lies in neither ROI.
    # With exclude_edges (configured): R1 counts k, R2 counts m.
    # Without it: R1 counts k + 1 (E3), R2 counts m + 2 (E1 and the part of E2).
    # Images (one unit each): ctrl (k, m) = (1, 1), (2, 1), (3, 2);
    #                         mut  (k, m) = (4, 2), (5, 3), (6, 3).
    # Mean of the two ROIs per unit: (k/825 + m/575) / 2, so
    #   Mean_A = (6/825 + 4/575) / 6, Mean_B = (15/825 + 8/575) / 6, and without
    #   edge exclusion (9/825 + 10/575) / 6 and (18/825 + 14/575) / 6: the
    #   difference (9/825 + 4/575) / 6 is the same in both.
    # Pooled per unit: (k + m) / 1400 (edges excluded), (k + m + 3) / 1400 (not):
    #   Mean_A = 10/4200, Mean_B = 23/4200 (excluded).
    # Threshold factor 0.6 (level 60) segments the same pixels; factor 2.5
    # (level 250) is above every pixel: every measurement is EMPTY, zero.
    p1 = [(30, 30), (70, 30), (110, 30), (30, 80), (70, 80), (110, 80)]
    p2 = [(150, 140), (200, 140), (200, 180)]
    always = [(3, 50), (236, 130), (190, 108), (60, 160)]
    design = [("ctrl_a.tif", 1, 1), ("ctrl_b.tif", 2, 1), ("ctrl_c.tif", 3, 2),
              ("mut_a.tif", 4, 2), ("mut_b.tif", 5, 3), ("mut_c.tif", 6, 3)]
    folder = os.path.join(tmp, "robust_micron")
    os.makedirs(folder)
    for name, k, m in design:
        imp = ImagePlus(name, robust_byte_plane(p1[:k] + p2[:m] + always, 10, 200))
        FileSaver(imp).saveAsTiff(os.path.join(folder, name))
        imp.close()
    r1 = Roi(0, 0, 120, 110)
    r1.setName("corner")
    rois = [r1, robust_notched_roi()]
    config = nuclei_config(E, {"manual_low": 100, "exclude_edges": True})
    options = E.default_global_options()
    options.update({"override_calibration": True, "pixel_width": "0.25", "pixel_height": "0.25",
                    "spatial_unit": "micron"})
    batch_options = {"extensions": "tif", "roi_policy": E.BATCH_ROI_MANAGER}
    grid = {"threshold_factor": [0.6, 2.5], "exclude_edges": ["false"]}
    result = E.run_robustness(folder, os.path.join(tmp, "robust_out"), [config], [config], options, batch_options,
                              1, grid, {"control_group": "ctrl"}, manager_rois=rois)
    t.equal(result["specifications"], 6, "3 threshold factors x 2 edge settings")
    rows = robust_table(E, result["run_dir"], "robustness_measurements.csv")
    t.equal(len(rows), 72, "6 specifications x 6 images x 2 parent ROIs")
    reference = [row for row in rows if row["Is_Reference"] == "True"]
    t.equal(len(reference), 12, "two reference measurements per image")
    areas = dict([(row["Parent_ROI_Name"], float(row["ROI_Area"])) for row in reference])
    t.close(areas["corner"], 825.0, "rectangle on the image corner: 13200 px x 0.0625 micron^2")
    t.close(areas["notched"], 575.0, "notched polygon on the image border: 9200 px x 0.0625 micron^2")
    t.equal(set([row["ROI_Area_Unit"] for row in rows]), set(["micron^2"]), "the area unit is the override unit")
    t.equal(set([row["Metric_Unit"] for row in rows]), set(["objects per micron^2"]), "the metric unit follows it")
    by_key = dict([((row["Specification"], row["Image_File"], row["Parent_ROI_Name"]), row) for row in rows])
    for name, k, m in design:
        for factor in ["0.6", "1.0"]:
            on = by_key[("threshold_factor=" + factor + "; exclude_edges=true", name, "corner")]
            t.equal(on["Object_Count"], str(k), name + " factor " + factor + ": objects wholly inside R1")
            t.close(float(on["Value"]), k / 825.0, name + " factor " + factor + ": density in R1")
            on = by_key[("threshold_factor=" + factor + "; exclude_edges=true", name, "notched")]
            t.equal(on["Object_Count"], str(m), name + " factor " + factor + ": objects wholly inside R2")
            off = by_key[("threshold_factor=" + factor + "; exclude_edges=false", name, "corner")]
            t.equal(off["Object_Count"], str(k + 1), name + " factor " + factor + ": R1 with the border-cut disk")
            off = by_key[("threshold_factor=" + factor + "; exclude_edges=false", name, "notched")]
            t.equal(off["Object_Count"], str(m + 2), name + " factor " + factor + ": R2 with the cut and crossing disks")
        t.close(float(by_key[("threshold_factor=1.0; exclude_edges=true", name, "corner")]["Threshold_Low"]), 100.0,
                name + ": the configured threshold")
        t.close(float(by_key[("threshold_factor=0.6; exclude_edges=true", name, "corner")]["Threshold_Low"]), 60.0,
                name + ": 0.6 x 100")
        for edges in ["true", "false"]:
            high = by_key[("threshold_factor=2.5; exclude_edges=" + edges, name, "notched")]
            t.equal((high["Status"], high["Object_Count"], float(high["Value"])), ("EMPTY", "0", 0.0),
                    name + ": a level of 250 lies above every 8-bit pixel of 200")
    specs = dict([(row["Specification"], row) for row in result["tables"]["specifications"]])
    difference = (9.0 / 825.0 + 4.0 / 575.0) / 6.0
    for factor in ["0.6", "1.0"]:
        on = specs["threshold_factor=" + factor + "; exclude_edges=true"]
        off = specs["threshold_factor=" + factor + "; exclude_edges=false"]
        t.close(on["Mean_A"], (6.0 / 825.0 + 4.0 / 575.0) / 6.0, "factor " + factor + ": mean of A over units (edges excluded)")
        t.close(on["Mean_B"], (15.0 / 825.0 + 8.0 / 575.0) / 6.0, "factor " + factor + ": mean of B over units (edges excluded)")
        t.close(off["Mean_A"], (9.0 / 825.0 + 10.0 / 575.0) / 6.0, "factor " + factor + ": mean of A with edge objects")
        t.close(off["Mean_B"], (18.0 / 825.0 + 14.0 / 575.0) / 6.0, "factor " + factor + ": mean of B with edge objects")
        for row in [on, off]:
            t.close(row["Difference_B_Minus_A"], difference, "factor " + factor + ": edge objects add the same to both groups")
            t.equal(row["Sign"], "positive", "factor " + factor + ": B above A")
    for edges in ["true", "false"]:
        t.equal(specs["threshold_factor=2.5; exclude_edges=" + edges]["Sign"], "zero", "an empty threshold gives no difference")
    summary = result["summary"]
    t.equal((summary["N_Same_Sign_As_Reference"], summary["N_Zero_Or_Undefined"], summary["N_Opposite_Sign"]), (4, 2, 0),
            "four specifications keep the direction, the two empty ones are zero")
    t.equal(summary["Sign_Flipping_Values"], "none", "nothing reverses the direction")
    t.equal(summary["Metric_Unit"], "objects per micron^2", "the summary states the unit")

    pooled = E.run_robustness(folder, os.path.join(tmp, "robust_out"), [config], [config], options, batch_options,
                              1, grid, {"control_group": "ctrl"}, unit_aggregation=E.ROBUST_UNIT_POOLED,
                              manager_rois=rois)
    specs = dict([(row["Specification"], row) for row in pooled["tables"]["specifications"]])
    t.close(specs["threshold_factor=1.0; exclude_edges=true"]["Mean_A"], 10.0 / 4200.0, "pooled: (2 + 3 + 5) objects / 1400 micron^2 / 3 units")
    t.close(specs["threshold_factor=1.0; exclude_edges=true"]["Mean_B"], 23.0 / 4200.0, "pooled: (6 + 8 + 9) / 1400 / 3")
    t.close(specs["threshold_factor=1.0; exclude_edges=false"]["Mean_A"], 19.0 / 4200.0, "pooled with edge objects: (10 + 9) / 1400 / 3")

    # The settings record of the run reproduces it when loaded as a profile.
    loaded_options, enabled, all_configs = E.load_profile_settings(os.path.join(result["run_dir"], "settings.properties"))
    again = E.run_robustness_from_options(folder, os.path.join(tmp, "robust_out"), enabled, all_configs, loaded_options,
                                          {"extensions": "tif", "roi_policy": loaded_options["batch_roi_policy"]},
                                          manager_rois=rois)
    first = dict([(row["Specification_ID"], row) for row in result["tables"]["specifications"]])
    second = dict([(row["Specification_ID"], row) for row in again["tables"]["specifications"]])
    t.equal(sorted(second.keys()), sorted(first.keys()), "the reloaded settings give the same specifications")
    for key in sorted(first.keys()):
        for column in ["Specification", "Mean_A", "Mean_B", "Hedges_g", "Welch_P", "Mann_Whitney_P", "Sign"]:
            t.equal(second[key][column], first[key][column], "reloaded settings, " + key + " " + column)


def robust_two_frame_image(name, counts):
    """A one-channel, two-frame 16-bit image at 0.5 um per pixel; frame f
    holds the first counts[f] disks of value 2000 on 100."""
    stack = ImageStack(WIDTH, HEIGHT)
    for count in counts:
        stack.addSlice("", robust_plane(count, 2000))
    imp = ImagePlus(name, stack)
    imp.setDimensions(1, 1, len(counts))
    imp.setOpenAsHyperStack(True)
    cal = imp.getCalibration()
    cal.pixelWidth = PIXEL
    cal.pixelHeight = PIXEL
    cal.setUnit("um")
    return imp


def test_robust_two_frames(E, t, tmp):
    """Two frames per image, two images of one organoid, a frame without
    objects, and Welch statistics with non-integer degrees of freedom checked
    against an independent implementation of Student's t."""
    # Disks per frame (frame 1, frame 2), whole image of 12000 um^2:
    #   ctrl_o1_a (2, 4), ctrl_o1_b (0, 6), ctrl_o2_a (3, 3), ctrl_o3_a (5, 3)
    #   mut_o1_a (8, 6), mut_o2_a (7, 9), mut_o3_a (10, 10)
    # Frame 1 of ctrl_o1_b holds no pixel above 1000: EMPTY, density 0.
    # Units (organoids), mean of all frames of all their images, per 12000 um^2:
    #   A = o1 (2+4+0+6)/4 = 3, o2 3, o3 4;   B = 7, 8, 10.
    # mean A = 10/3, var A = 1/3; mean B = 25/3, var B = 7/3; difference 5.
    # pooled var = (2/3 + 14/3) / 4 = 4/3; d = 5 / sqrt(4/3); J = 1 - 3/15 = 0.8;
    #   g = 4 / sqrt(4/3) = 2 sqrt(3) = 3.4641016.
    # Welch: se^2 = (1/3)/3 + (7/3)/3 = 8/9; t = 5 / sqrt(8/9) = 15 / (2 sqrt 2) = 5.3033009;
    #   df = (8/9)^2 / ((1/9)^2/2 + (7/9)^2/2) = (64/81) / (50/162) = 2.56.
    # Mann-Whitney: complete separation of 3 and 3 (the tie inside A does not
    #   matter): U of B = 9, exact p = 2/20 = 0.1.
    # Thresholds 500 and 1500 count the same disks.
    design = [("ctrl_o1_a.tif", (2, 4)), ("ctrl_o1_b.tif", (0, 6)), ("ctrl_o2_a.tif", (3, 3)), ("ctrl_o3_a.tif", (5, 3)),
              ("mut_o1_a.tif", (8, 6)), ("mut_o2_a.tif", (7, 9)), ("mut_o3_a.tif", (10, 10))]
    folder = os.path.join(tmp, "robust_frames")
    os.makedirs(folder)
    for name, counts in design:
        imp = robust_two_frame_image(name, counts)
        FileSaver(imp).saveAsTiffStack(os.path.join(folder, name))
        imp.close()
    configs = [nuclei_config(E)]
    options = E.default_global_options()
    options["frame_mode"] = "All frames"
    groups = {"source": E.ROBUST_GROUP_REGEX, "group_regex": "^([a-z]+)_", "unit_regex": "^([a-z]+_o[0-9]+)_",
              "control_group": "ctrl"}
    result = E.run_robustness(folder, os.path.join(tmp, "robust_out"), configs, configs, options, {"extensions": "tif"},
                              1, {"threshold_factor": [0.5, 1.5]}, groups)
    rows = robust_table(E, result["run_dir"], "robustness_measurements.csv")
    t.equal(len(rows), 42, "7 images x 2 frames x 3 specifications")
    reference = dict([((row["Image_File"], row["Frame"]), row) for row in rows if row["Is_Reference"] == "True"])
    for name, counts in design:
        for frame in [1, 2]:
            t.equal(reference[(name, str(frame))]["Object_Count"], str(counts[frame - 1]), name + " frame " + str(frame))
    t.equal(reference[("ctrl_o1_b.tif", "1")]["Status"], "EMPTY", "a frame without signal is an empty measurement")
    t.equal(float(reference[("ctrl_o1_b.tif", "1")]["Value"]), 0.0, "and counts as density zero")
    report = robust_table(E, result["run_dir"], "robustness_image_report.csv")
    t.equal(set([row["Frames"] for row in report]), set(["2"]), "both frames of every image are measured")
    units = dict([(row["Unit"], row) for row in result["tables"]["units"] if row["Is_Reference"]])
    t.equal((units["ctrl_o1"]["N_Images"], units["ctrl_o1"]["N_Measurements"], units["ctrl_o1"]["N_Measurements_Used"]),
            (2, 4, 4), "two images with two frames each are one unit of four measurements")
    for unit, value in [("ctrl_o1", 3.0), ("ctrl_o2", 3.0), ("ctrl_o3", 4.0), ("mut_o1", 7.0), ("mut_o2", 8.0), ("mut_o3", 10.0)]:
        t.close(units[unit]["Value"] * ROBUST_AREA, value, unit + ": mean over its frames and images")
    df = 2.56
    t_value = 15.0 / (2.0 * math.sqrt(2.0))
    se = math.sqrt(8.0 / 9.0)
    t.close(robust_t_two_sided(4.2866070, 4.0), robust_t4_two_sided_p(4.2866070), "the independent t distribution agrees with the closed form for 4 df")
    t.close(robust_t_quantile_975(4.0), ROBUST_T4_975, "and gives the tabulated quantile for 4 df")
    for row in result["tables"]["specifications"]:
        label = row["Specification"]
        t.equal((row["N_Units_A"], row["N_Units_B"]), (3, 3), label + ": three organoids per group")
        t.close(row["Mean_A"] * ROBUST_AREA, 10.0 / 3.0, label + ": mean of A")
        t.close(row["SD_A"] * ROBUST_AREA, math.sqrt(1.0 / 3.0), label + ": SD of A")
        t.close(row["Difference_B_Minus_A"] * ROBUST_AREA, 5.0, label + ": difference")
        t.close(row["Hedges_g"], 2.0 * math.sqrt(3.0), label + ": Hedges g")
        t.close(row["Welch_t"], t_value, label + ": Welch t")
        t.close(row["Welch_df"], df, label + ": Welch-Satterthwaite df is not an integer")
        t.close(row["Mann_Whitney_U"], 9.0, label + ": U of B")
        t.close(row["Mann_Whitney_P"], 0.1, label + ": exact Mann-Whitney p")
        if E.ROBUST_MATH["available"]:
            t.close(row["Welch_P"], robust_t_two_sided(t_value, df), label + ": Welch p from an independent t distribution", 1e-8)
            half = robust_t_quantile_975(df) * se
            t.close(row["Difference_CI95_Low"] * ROBUST_AREA, 5.0 - half, label + ": lower end of the Welch interval", 1e-8)
            t.close(row["Difference_CI95_High"] * ROBUST_AREA, 5.0 + half, label + ": upper end of the Welch interval", 1e-8)
    t.close(result["summary"]["Fraction_Same_Sign_As_Reference"], 1.0, "the direction holds")


def test_robust_light_background(E, t, tmp):
    """Dark objects on a light background with an automatic threshold: the
    factor moves the upper threshold level. Images without calibration
    report pixel units."""
    # Disks of value 1000 on 3000, no calibration (48000 pixels):
    #   ctrl 2, 3, 4 disks (mean 3); mut 5, 6, 7 disks (mean 6) per 48000 px.
    # Otsu with a light background thresholds [1000, L]: the dark disks, with
    # 1000 <= L < 3000 (ImageJ reports a level just below the background).
    # The test needs 2500 <= L < 3000, checked first. Then, with the factor f
    # on the upper level f L (the lower level stays at the image minimum 1000):
    #   f = 0.5, 1.0: 1000 <= f L < 3000: exactly the disks -> difference +3 per 48000 px
    #   f = 0.3: 0.3 L < 1000 = darkest pixel: nothing is positive (EMPTY), zero
    #   f = 1.2: 1.2 L >= 3000: every pixel is positive, one object per image, zero difference
    images = [("ctrl_a.tif", 2), ("ctrl_b.tif", 3), ("ctrl_c.tif", 4), ("mut_a.tif", 5), ("mut_b.tif", 6), ("mut_c.tif", 7)]
    folder = os.path.join(tmp, "robust_light")
    os.makedirs(folder)
    for name, count in images:
        imp = ImagePlus(name, robust_plane(count, 1000, background=3000))
        FileSaver(imp).saveAsTiff(os.path.join(folder, name))
        imp.close()
    config = nuclei_config(E, {"threshold_mode": "Auto", "auto_method": "Otsu", "dark_background": False})
    result = robust_run(E, folder, tmp, {"threshold_factor": [0.3, 0.5, 1.2]}, config=config)
    rows = robust_table(E, result["run_dir"], "robustness_measurements.csv")
    by_key = dict([((row["Specification"], row["Image_File"]), row) for row in rows])
    for name, count in images:
        reference = by_key[("threshold_factor=1.0", name)]
        level = float(reference["Threshold_High"])
        t.check(2500.0 <= level < 3000.0, name + ": the light-background level lies between the disks and the background")
        t.close(float(reference["Threshold_Low"]), 1000.0, name + ": the lower level is the darkest pixel")
        t.equal(reference["Object_Count"], str(count), name + ": the reference counts the dark disks")
        half = by_key[("threshold_factor=0.5", name)]
        t.close(float(half["Threshold_High"]), 0.5 * level, name + ": factor 0.5 halves the upper level")
        t.close(float(half["Threshold_Low"]), 1000.0, name + ": and leaves the lower level")
        t.equal(half["Object_Count"], str(count), name + ": the same disks at half the level")
        low = by_key[("threshold_factor=0.3", name)]
        t.equal((low["Status"], low["Object_Count"]), ("EMPTY", "0"), name + ": 0.3 L is below every pixel")
        high = by_key[("threshold_factor=1.2", name)]
        t.equal((high["Status"], high["Object_Count"]), ("OK", "1"), name + ": 1.2 L takes in the background: one object")
        t.close(float(high["Threshold_High"]), 3000.0, name + ": the upper level is clipped at the brightest pixel")
        t.equal(reference["ROI_Area_Unit"], "pixel^2", name + ": an uncalibrated image is measured in pixels")
        t.close(float(reference["ROI_Area"]), 48000.0, name + ": the whole image in pixels")
    specs = dict([(row["Specification"], row) for row in result["tables"]["specifications"]])
    for factor in ["0.5", "1.0"]:
        t.close(specs["threshold_factor=" + factor]["Difference_B_Minus_A"] * 48000.0, 3.0, "factor " + factor + ": 6 - 3 disks")
    for factor in ["0.3", "1.2"]:
        t.equal(specs["threshold_factor=" + factor]["Sign"], "zero", "factor " + factor + ": no difference")
    summary = result["summary"]
    t.equal((summary["N_Same_Sign_As_Reference"], summary["N_Zero_Or_Undefined"], summary["N_Opposite_Sign"]), (2, 2, 0),
            "two keep the direction, two have none")
    t.equal(summary["Metric_Unit"], "objects per pixel^2", "the metric unit in pixels")


def test_robust_settings_and_isolation(E, t, tmp):
    """Profiles: old ones load with this module's defaults, its settings make
    a round trip, and loading the module changes no table of a normal run."""
    # A profile without any robust.* key (the shipped example) gets the defaults.
    example = os.path.join(E.CELLQUANT_PACKAGE_DIR, "examples", "Example_4Channel_SOX2_Casp3_MAP2_DAPI.properties")
    if os.path.isfile(example):
        t.check("robust." not in read_utf8(example), "the example profile predates this module")
        options, enabled, all_configs = E.load_profile_settings(example)
        for spec in E.ROBUST_SETTING_SPECS:
            t.equal(options[spec[1]], spec[3], "an old profile gets the default of " + spec[0])
        request = E.robust_request_from_options(options)
        t.equal((request["channel"], request["metric"], request["max_specifications"]), (1, E.ROBUST_METRIC_DENSITY, 60),
                "the defaults give a usable request")
    # Non-default values make a round trip through a profile.
    configs = [nuclei_config(E)]
    options = E.default_global_options()
    chosen = {
        "robust_channel": "1", "robust_metric": E.ROBUST_METRIC_AREA_FRACTION,
        "robust_grid": "threshold_factor=0.7,1.3; fill_holes=true", "robust_max_specifications": "12",
        "robust_group_source": E.ROBUST_GROUP_REGEX, "robust_separator": "-", "robust_group_regex": "^(wt|ko)",
        "robust_unit_regex": "(org[0-9]+)", "robust_group_csv": os.path.join(tmp, "groups.csv"),
        "robust_control_group": "wt", "robust_unit_aggregation": E.ROBUST_UNIT_POOLED
    }
    options.update(chosen)
    path = os.path.join(tmp, "robust_profile.properties")
    E.write_settings_properties(path, options, configs)
    loaded, enabled, all_configs = E.load_profile_settings(path)
    for key in sorted(chosen.keys()):
        t.equal(loaded[key], chosen[key], "profile round trip of " + key)
    request = E.robust_request_from_options(loaded)
    t.equal(request["grid"], {"threshold_factor": ["0.7", "1.3"], "fill_holes": ["true"]}, "the grid is read back")
    t.equal(request["groups"]["unit_regex"], "(org[0-9]+)", "the unit expression is read back")
    t.equal(request["max_specifications"], 12, "the limit is read back")

    # Loading the module changes nothing a normal batch writes: an engine
    # alone and an engine with this module write the same tables.
    folder = robust_folder(tmp, "robust_isolation", [("a_1.tif", 4, 2000), ("b_1.tif", 7, 2000)])
    engine_path = os.path.join(E.CELLQUANT_PACKAGE_DIR, "cellquantengine.py")
    module_path = os.path.join(E.CELLQUANT_PACKAGE_DIR, "cellquantrobust.py")
    plain = {"CELLQUANT_NO_GUI": True}
    E.run_source_file(engine_path, plain)
    loaded_ns = {"CELLQUANT_NO_GUI": True}
    E.run_source_file(engine_path, loaded_ns)
    E.run_source_file(module_path, loaded_ns)
    t.check("run_robustness" in loaded_ns and "run_robustness" not in plain, "one engine has the module, one has not")
    texts = []
    for namespace, name in [(plain, "robust_plain"), (loaded_ns, "robust_loaded")]:
        channel = [namespace["normalize_channel_config"](1, {
            "enabled": True, "name": "DAPI", "analysis_mode": "Objects + area", "threshold_mode": "Manual", "manual_low": 1000
        })]
        batch = namespace["run_batch"](folder, os.path.join(tmp, name), channel, channel,
                                       namespace["default_global_options"](), {"extensions": "tif"})
        tables = {}
        for table in ["batch_roi_summary.csv", "batch_objects.csv", "batch_thresholds.csv"]:
            tables[table] = robust_without_run_names(read_utf8(os.path.join(batch["batch_dir"], table)))
        props = namespace["read_properties_file"](os.path.join(batch["batch_dir"], "settings.properties"))
        tables["settings"] = dict([(E.to_text(key), E.to_text(props.getProperty(key))) for key in props.stringPropertyNames()])
        texts.append(tables)
    for table in ["batch_roi_summary.csv", "batch_objects.csv", "batch_thresholds.csv"]:
        t.check(len(texts[0][table].split("\n")) > 2, table + " holds rows")
        t.check(texts[0][table] == texts[1][table], table + " is identical with and without the module")
    extra = sorted([key for key in texts[1]["settings"] if key not in texts[0]["settings"]])
    t.equal(extra, sorted([spec[0] for spec in E.ROBUST_SETTING_SPECS]), "the settings record only gains the robust.* keys")
    t.check(len([key for key in texts[0]["settings"] if texts[1]["settings"].get(key) != texts[0]["settings"][key]]) == 0,
            "every other setting is recorded with the same value")


def test_robust_table_names(E, t, tmp):
    """A group table can name a top-level image whose file name also occurs
    in a subfolder, and a repeated row without a unit keeps the unit."""
    names = ["a.tif", os.path.join("sub", "a.tif"), "b.tif", "c.tif"]
    table = os.path.join(tmp, "robust_names.csv")
    handle = codecs.open(table, "w", "utf-8")
    handle.write(u"Image,Group,Unit\n")
    handle.write(u"a.tif,ctrl,org1\n")         # the top-level a.tif exactly
    handle.write(u"sub/a.tif,mut,org2\n")      # the one in the subfolder
    handle.write(u"b,ctrl,org3\n")             # without its extension
    handle.write(u"b.tif,ctrl,\n")             # repeated, no unit: keeps org3
    handle.write(u"c.tif,ctrl,org4\n")
    handle.write(u"c.tif,mut,org4\n")          # repeated with another group: conflict
    handle.close()
    rows = E.robust_assign_groups(names, E.robust_group_rule({"source": E.ROBUST_GROUP_CSV, "csv_path": table}))
    t.equal([(row["group"], row["unit"], row["status"]) for row in rows[:3]],
            [("ctrl", "org1", "ASSIGNED"), ("mut", "org2", "ASSIGNED"), ("ctrl", "org3", "ASSIGNED")],
            "exact relative paths are matched before file names; a repeated row keeps its unit")
    t.equal(rows[3]["status"], "NOT ASSIGNED", "an image listed with two groups is not assigned")
    t.check("more than once" in rows[3]["message"], "and the reason is given")
    t.equal(len(rows), 4, "no table row is left without an image")


def robust_float_image(name, n_small, n_big):
    """A 32-bit image, 0.5 um per pixel, background 0.0:
    n_small squares of 4 x 4 px (16 px = 4 um^2) of value 3.0 in the row y = 80,
    n_big squares of 8 x 8 px (64 px = 16 um^2) of value 1.0 in the row y = 110,
    and two decoy 4 x 4 squares of value 3.0 at (8, 8) and (226, 186), which
    lie inside the bounding box of the oval parent ROI but outside the oval."""
    ip = FloatProcessor(WIDTH, HEIGHT)
    for x in [40, 70, 100, 130, 160, 190][:n_small]:
        draw_rect(ip, (x, 80, 4, 4), 3.0)
    for x in [50, 90, 130, 170][:n_big]:
        draw_rect(ip, (x, 110, 8, 8), 1.0)
    for corner in [(8, 8), (226, 186)]:
        draw_rect(ip, (corner[0], corner[1], 4, 4), 3.0)
    imp = ImagePlus(name, ip)
    cal = imp.getCalibration()
    cal.pixelWidth = PIXEL
    cal.pixelHeight = PIXEL
    cal.setUnit("um")
    return imp


def test_robust_object_area_float(E, t, tmp):
    """Mean object area on 32-bit images inside an oval parent ROI, with
    units whose value is undefined and specifications without any unit."""
    # Parent ROI: OvalRoi(5, 5, 230, 190), centre (120, 100), half-axes 115 and 95.
    # Every square of the two rows lies inside it (the farthest corner, (40, 80),
    # gives (80/115)^2 + (20/95)^2 = 0.53 < 1); the decoys give 1.74 > 1.
    # Manual threshold 0.5 (32-bit: compared exactly); grid threshold_factor
    # 4, 8 and min_size 5 um^2 (20 px, removes the 16 px squares):
    #   (n_small, n_big): ctrl (6,2) (4,2) (2,2); mut (2,4) (1,4) (0,4).
    #   f=1, min 0 (reference): mean object area = (4 n_small + 16 n_big) / (n_small + n_big)
    #       ctrl 56/8 = 7, 48/6 = 8, 40/4 = 10; mut 72/6 = 12, 68/5 = 13.6, 64/4 = 16.
    #       A: mean 25/3, var ((4/3)^2 + (1/3)^2 + (5/3)^2)/2 = 7/3.
    #       B: mean 208/15, deviations -28/15, -4/15, 32/15: var 1824/225/2 = 912/225.
    #       difference 83/15; ratio 208/125; pooled var (2*7/3 + 2*912/225)/4 = 1437/450;
    #       g = (1 - 3/15) * (83/15) / sqrt(1437/450).
    #       Welch se^2 = 7/9 + 304/225 = 479/225, t = 83 / sqrt(479),
    #       df = (479/225)^2 / ((7/9)^2/2 + (304/225)^2/2). Mann-Whitney: separation, p = 0.1.
    #   f=1, min 5: only the big squares: 16 everywhere, difference 0 (sign zero).
    #   f=4 (level 2.0), min 0: only the small squares: 4 everywhere; mut (0,4) has no
    #       object (UNDEFINED), leaving 2 units in B: means only, difference 0.
    #   f=4, min 5: no object anywhere: no unit has a value, no difference.
    #   f=8 (level 4.0 above the maximum 3.0): nothing positive, no unit has a value.
    # So 3 specifications are compared: the reference (positive) and two zeros.
    design = [("ctrl_a.tif", 6, 2), ("ctrl_b.tif", 4, 2), ("ctrl_c.tif", 2, 2),
              ("mut_a.tif", 2, 4), ("mut_b.tif", 1, 4), ("mut_c.tif", 0, 4)]
    folder = os.path.join(tmp, "robust_float")
    os.makedirs(folder)
    for name, n_small, n_big in design:
        imp = robust_float_image(name, n_small, n_big)
        FileSaver(imp).saveAsTiff(os.path.join(folder, name))
        imp.close()
    oval = OvalRoi(5, 5, 230, 190)
    oval.setName("oval")
    config = nuclei_config(E, {"manual_low": 0.5})
    result = E.run_robustness(folder, os.path.join(tmp, "robust_out"), [config], [config], E.default_global_options(),
                              {"extensions": "tif", "roi_policy": E.BATCH_ROI_MANAGER}, 1,
                              {"threshold_factor": [4, 8], "min_size": [5]}, {"control_group": "ctrl"},
                              metric=E.ROBUST_METRIC_OBJECT_AREA, manager_rois=[oval])
    t.equal(result["specifications"], 6, "3 threshold factors x 2 minimum sizes")
    rows = robust_table(E, result["run_dir"], "robustness_measurements.csv")
    t.equal(len(rows), 36, "6 specifications x 6 images x 1 parent ROI")
    by_key = dict([((row["Specification"], row["Image_File"]), row) for row in rows])
    for name, n_small, n_big in design:
        reference = by_key[("threshold_factor=1.0; min_size=0.0", name)]
        t.equal(reference["Object_Count"], str(n_small + n_big), name + ": the decoys outside the oval are not counted")
        t.close(float(reference["Total_Object_Area"]), 4.0 * n_small + 16.0 * n_big, name + ": total object area in um^2")
        t.close(float(reference["Value"]), (4.0 * n_small + 16.0 * n_big) / (n_small + n_big), name + ": mean object area")
        t.close(float(by_key[("threshold_factor=4.0; min_size=0.0", name)]["Threshold_Low"]), 2.0, name + ": 4 x 0.5")
        high = by_key[("threshold_factor=8.0; min_size=0.0", name)]
        t.equal((high["Status"], high["Object_Count"], high["Value"]), ("UNDEFINED", "0", ""),
                name + ": a level above every pixel leaves no object and no mean object area")
        t.check("nothing is positive" in high["Note"] and "undefined" in high["Note"], name + ": and the note gives both reasons")
    t.equal(reference["Metric_Unit"], ROBUST_UM + "^2", "a mean object area is in square micrometers")
    t.equal(by_key[("threshold_factor=4.0; min_size=0.0", "mut_c.tif")]["Status"], "UNDEFINED",
            "an image without small squares has no object at level 2.0")
    t.equal(len([row for row in rows if row["Status"] == "UNDEFINED"]), 19, "1 + 6 + 12 measurements without a value")

    specs = dict([(row["Specification"], row) for row in result["tables"]["specifications"]])
    ref = specs["threshold_factor=1.0; min_size=0.0"]
    t.close(ref["Mean_A"], 25.0 / 3.0, "reference: mean of A")
    t.close(ref["Mean_B"], 208.0 / 15.0, "reference: mean of B")
    t.close(ref["SD_A"], math.sqrt(7.0 / 3.0), "reference: SD of A")
    t.close(ref["SD_B"], math.sqrt(912.0 / 225.0), "reference: SD of B")
    t.close(ref["Difference_B_Minus_A"], 83.0 / 15.0, "reference: difference")
    t.close(ref["Ratio_B_Over_A"], 208.0 / 125.0, "reference: ratio")
    t.close(ref["Hedges_g"], 0.8 * (83.0 / 15.0) / math.sqrt(1437.0 / 450.0), "reference: Hedges g")
    t.close(ref["Welch_t"], 83.0 / math.sqrt(479.0), "reference: Welch t")
    df = (479.0 / 225.0) ** 2 / ((7.0 / 9.0) ** 2 / 2.0 + (304.0 / 225.0) ** 2 / 2.0)
    t.close(ref["Welch_df"], df, "reference: Welch-Satterthwaite df")
    t.close(ref["Mann_Whitney_P"], 0.1, "reference: exact Mann-Whitney p")
    if E.ROBUST_MATH["available"]:
        t.close(ref["Welch_P"], robust_t_two_sided(83.0 / math.sqrt(479.0), df), "reference: Welch p", 1e-8)
        half = robust_t_quantile_975(df) * math.sqrt(479.0 / 225.0)
        t.close(ref["Difference_CI95_Low"], 83.0 / 15.0 - half, "reference: interval, lower end", 1e-8)
        t.close(ref["Difference_CI95_High"], 83.0 / 15.0 + half, "reference: interval, upper end", 1e-8)
    big_only = specs["threshold_factor=1.0; min_size=5.0"]
    robust_close_each(t, (big_only["Mean_A"], big_only["Mean_B"]), (16.0, 16.0), "min_size 5: only the 16 um^2 squares are left")
    t.equal((big_only["Sign"], big_only["Hedges_g"]), ("zero", ""), "equal groups without variation: zero, no g")
    t.close(big_only["Mann_Whitney_P"], 1.0, "all six values tied: rank-test p is 1")
    small_only = specs["threshold_factor=4.0; min_size=0.0"]
    t.equal((small_only["N_Units_A"], small_only["N_Units_B"], small_only["Sign"]), (3, 2, "zero"),
            "level 2.0: the unit without objects drops out, the rest are 4 um^2")
    t.check("Fewer than 3 units" in small_only["Note"], "and with two units only the means are given")
    for label in ["threshold_factor=4.0; min_size=5.0", "threshold_factor=8.0; min_size=0.0", "threshold_factor=8.0; min_size=5.0"]:
        t.equal((specs[label]["N_Units_A"], specs[label]["Sign"]), (0, ""), label + ": no unit has a value")
    summary = result["summary"]
    t.equal((summary["N_Specifications_Compared"], summary["N_Same_Sign_As_Reference"], summary["N_Zero_Or_Undefined"],
             summary["N_Opposite_Sign"]), (3, 1, 2, 0), "three compared: the reference and two zero differences")
    t.close(summary["Fraction_Same_Sign_As_Reference"], 1.0 / 3.0, "the zero differences count against the share")
    t.equal(summary["Sign_Flipping_Values"], "none", "a zero difference is not a reversal")
    robust_close_each(t, (summary["Difference_Minimum"], summary["Difference_Median"], summary["Difference_Maximum"]),
            (0.0, 0.0, 83.0 / 15.0), "range of the difference over the compared specifications")
    robust_close_each(t, (summary["Ratio_Minimum"], summary["Ratio_Maximum"]), (1.0, 208.0 / 125.0), "range of the ratio")
    robust_close_each(t, (summary["Hedges_g_Minimum"], summary["Hedges_g_Maximum"]), (ref["Hedges_g"], ref["Hedges_g"]),
            "only the reference has an effect size")


def robust_record(spec, group, unit, image, numerator, denominator, status="OK"):
    value = None if status == "FAILED" else float(numerator) / float(denominator)
    return {"spec": spec, "group": group, "unit": unit, "image": image, "value": value,
            "numerator": None if status == "FAILED" else numerator,
            "denominator": None if status == "FAILED" else denominator, "status": status}


def test_robust_tables_from_records(E, t, tmp):
    """Units, groups, the specification curve, and the summary from hand-made
    measurements, plus a small Welch and Mann-Whitney example worked out by hand."""
    # Specifications of threshold_factor 0.5, 1.0, 2.0: S001, S002 (reference), S003.
    # Group ctrl: unit s1 = two images (1 object / 1 area, 6 / 2), units a2 = 3/1, a3 = 4/1.
    # Group mut:  unit s1 (the same name, another group: another unit) and b2, b3.
    #   S002: mut = 5, 6, 7. Unit means: A = 2 (mean of 1 and 3), 3, 4; B = 5, 6, 7:
    #         difference 3, SDs 1, g = 0.8 * 3 = 2.4, Welch t = 3 / sqrt(2/3), df 4.
    #         Pooled over a unit instead: A's s1 = (1 + 6) / (1 + 2) = 7/3.
    #   S001: mut = 0, 1, 2: difference -2, g = -1.6 (reversed).
    #   S003: like S002, but the measurement of a2 failed: incomplete, left out of
    #         the summary; with 2 units left in A it has means only (difference 3).
    density = E.ROBUST_METRIC_DENSITY
    config = nuclei_config(E)
    plan = E.robust_plan([config], 1, {"threshold_factor": [0.5, 2.0]}, {}, density, E.ROBUST_UNIT_MEAN, 60)
    t.equal([spec["label"] for spec in plan["specs"]],
            ["threshold_factor=0.5", "threshold_factor=1.0", "threshold_factor=2.0"], "three specifications")
    records = []
    for spec in range(3):
        records += [robust_record(spec, "ctrl", "s1", "i1", 1, 1.0), robust_record(spec, "ctrl", "s1", "i2", 6, 2.0),
                    robust_record(spec, "ctrl", "a2", "i3", 3, 1.0, "FAILED" if spec == 2 else "OK"),
                    robust_record(spec, "ctrl", "a3", "i4", 4, 1.0)]
        offset = 0 if spec == 0 else 5
        records += [robust_record(spec, "mut", "s1", "i5", offset, 1.0), robust_record(spec, "mut", "b2", "i6", offset + 1, 1.0),
                    robust_record(spec, "mut", "b3", "i7", offset + 2, 1.0)]
    base = {"Run_ID": "R", "Application_Version": E.APP_VERSION}
    tables = E.robust_statistics_tables(plan, records, ["ctrl", "mut"], base, "objects per pixel^2")
    units = [row for row in tables["units"] if row["Specification_ID"] == "S002"]
    t.equal([(row["Group"], row["Unit"]) for row in units],
            [("ctrl", "a2"), ("ctrl", "a3"), ("ctrl", "s1"), ("mut", "b2"), ("mut", "b3"), ("mut", "s1")],
            "a unit name used in two groups is two units")
    s1 = [row for row in units if row["Unit"] == "s1" and row["Group"] == "ctrl"][0]
    t.equal((s1["N_Images"], s1["N_Measurements"], s1["N_Measurements_Used"]), (2, 2, 2), "unit s1 of ctrl holds two images")
    t.close(s1["Value"], 2.0, "mean of the measurements of a unit")
    failed_unit = [row for row in tables["units"] if row["Specification_ID"] == "S003" and row["Unit"] == "a2"][0]
    t.equal((failed_unit["N_Measurements_Used"], failed_unit["Value"]), (0, ""), "a failed measurement gives the unit no value")
    specs = dict([(row["Specification_ID"], row) for row in tables["specifications"]])
    robust_close_each(t, (specs["S002"]["Difference_B_Minus_A"], specs["S002"]["Hedges_g"], specs["S002"]["Welch_df"]), (3.0, 2.4, 4.0),
            "reference: difference, g, df")
    t.close(specs["S002"]["Welch_t"], 3.0 / math.sqrt(2.0 / 3.0), "reference: Welch t")
    robust_close_each(t, (specs["S001"]["Difference_B_Minus_A"], specs["S001"]["Hedges_g"]), (-2.0, -1.6), "factor 0.5: reversed")
    t.equal((specs["S001"]["Same_Sign_As_Reference"], specs["S002"]["Same_Sign_As_Reference"]), (False, True), "directions")
    t.equal((specs["S003"]["Complete"], specs["S003"]["N_Failed_Measurements"], specs["S003"]["N_Measurements"]), (False, 1, 7),
            "a failed measurement makes the specification incomplete")
    t.equal(specs["S003"]["Same_Sign_As_Reference"], "", "an incomplete specification is not compared")
    t.check("A measurement failed" in specs["S003"]["Note"] and "Fewer than 3 units" in specs["S003"]["Note"], "and says why")
    t.close(specs["S003"]["Difference_B_Minus_A"], 3.0, "its means are still reported: 6 - 3")
    t.equal([row["Specification_ID"] for row in tables["specifications"]], ["S001", "S002", "S003"],
            "without an effect size for every specification the curve is sorted by the difference (ties by number)")
    summary = tables["summary"]
    t.equal((summary["N_Specifications"], summary["N_Specifications_Incomplete"], summary["N_Specifications_Compared"]),
            (3, 1, 2), "three specifications, one incomplete, two compared")
    t.equal((summary["N_Same_Sign_As_Reference"], summary["N_Opposite_Sign"]), (1, 1), "one keeps, one reverses")
    t.close(summary["Fraction_Same_Sign_As_Reference"], 0.5, "share of compared specifications with the reference sign")
    robust_close_each(t, (summary["Hedges_g_Minimum"], summary["Hedges_g_Median"], summary["Hedges_g_Maximum"]), (-1.6, 0.4, 2.4),
            "range of g over the compared specifications only")
    t.equal(summary["Sign_Flipping_Values"], "threshold_factor=0.5 (1 of 1)", "the reversing value is named")
    values = dict([(row["Value"], row) for row in tables["values"]])
    t.equal((values["2.0"]["N_Specifications"], values["2.0"]["Fraction_Opposite_Sign"]), (0, ""),
            "the incomplete specification is not counted for its value")
    if E.ROBUST_MATH["available"]:
        # S002: t = 3.6742346 with 4 df: p = 0.0213 < 0.05 in the reference direction; S001 is reversed.
        t.close(specs["S002"]["Welch_P"], robust_t4_two_sided_p(3.0 / math.sqrt(2.0 / 3.0)), "reference: Welch p")
        t.equal((summary["Descriptive_N_Welch_P_Below_0.05_Same_Sign"], summary["Descriptive_N_Mann_Whitney_P_Below_0.05_Same_Sign"]),
                (1, 0), "descriptive counts of p below 0.05 in the reference direction")
    plan["aggregation"] = E.ROBUST_UNIT_POOLED
    pooled = E.robust_statistics_tables(plan, records, ["ctrl", "mut"], base, "objects per pixel^2")
    pooled_specs = dict([(row["Specification_ID"], row) for row in pooled["specifications"]])
    t.close(pooled_specs["S002"]["Mean_A"], (7.0 / 3.0 + 3.0 + 4.0) / 3.0, "pooled: s1 = 7 objects / 3 area units")

    # A = 1, 4, 6 and B = 2, 3, 5, 7 (no ties, ranks are the values).
    # Rank sum of A = 11, expected 3 * 8 / 2 = 12. Of the 35 three-element
    # subsets of 1..7 exactly 5 have the sum 12 ({1,4,7} {1,5,6} {2,3,7} {2,4,6}
    # {3,4,5}); all others lie at least 1 away: exact p = 30/35 = 6/7.
    # U of B = pairs with b > a = 1 + 1 + 2 + 3 = 7.
    # mean A = 11/3, var A = 19/3; mean B = 17/4, var B = 59/12; difference 7/12.
    # Welch se^2 = 19/9 + 59/48 = 481/144, t = 7 / sqrt(481);
    # df = (481/144)^2 / ((19/9)^2/2 + (59/48)^2/3).
    # pooled var = (2 * 19/3 + 3 * 59/12) / 5 = 329/60; J = 1 - 3/19 = 16/19.
    s = E.robust_two_group_statistics([1, 4, 6], [2, 3, 5, 7])
    t.close(s["mw_u"], 7.0, "U of B, counted pair by pair")
    t.close(s["mw_p"], 6.0 / 7.0, "exact Mann-Whitney p, enumerated by hand")
    t.close(s["difference"], 7.0 / 12.0, "difference of the means")
    t.close(s["welch_t"], 7.0 / math.sqrt(481.0), "Welch t")
    df = (481.0 / 144.0) ** 2 / ((19.0 / 9.0) ** 2 / 2.0 + (59.0 / 48.0) ** 2 / 3.0)
    t.close(s["welch_df"], df, "Welch-Satterthwaite df")
    t.close(s["hedges_g"], 16.0 / 19.0 * (7.0 / 12.0) / math.sqrt(329.0 / 60.0), "Hedges g")
    if E.ROBUST_MATH["available"]:
        t.close(s["welch_p"], robust_t_two_sided(7.0 / math.sqrt(481.0), df), "Welch p from an independent t distribution", 1e-8)


def test_robust_manual_zero_level(E, t, tmp):
    """A threshold factor multiplies the Manual lower level; a level of 0
    cannot be scaled, and a negative one would move the wrong way."""
    density = E.ROBUST_METRIC_DENSITY
    # Dark objects on a light background thresholded as the band [0, 800]:
    # 0.5 x 0 = 2 x 0 = 0, so the factor changes nothing and is refused,
    # instead of counting three identical segmentations as three results.
    band = nuclei_config(E, {"manual_low": 0, "manual_high": 800})
    message = robust_refused(E.robust_build_specifications, band, {"threshold_factor": [0.5, 2.0]}, density)
    t.check(message is not None and "no effect" in message, "a factor on a lower level of 0 is refused")
    # Together with a manual_low of 100 the factor applies there only:
    # manual_low 0 gives one specification, manual_low 100 gives three.
    plan = E.robust_build_specifications(band, {"threshold_factor": [0.5, 2.0], "manual_low": [100]}, density)
    t.equal([spec["label"] for spec in plan["specs"]],
            ["manual_low=0.0", "threshold_factor=0.5; manual_low=100.0", "threshold_factor=1.0; manual_low=100.0",
             "threshold_factor=2.0; manual_low=100.0"], "the factor is counted only where the level is not 0")
    t.equal([spec["is_reference"] for spec in plan["specs"]], [True, False, False, False], "the configured setting is the reference")
    negative = nuclei_config(E, {"manual_low": -50})
    message = robust_refused(E.robust_build_specifications, negative, {"threshold_factor": [0.5]}, density)
    t.check(message is not None and "negative" in message, "a factor on a negative lower level is refused")


TESTS.append(("robustness: group statistics against hand-derived numbers", test_robust_statistics))
TESTS.append(("robustness: grid of specifications and refusals", test_robust_grid))
TESTS.append(("robustness: group and unit assignment", test_robust_groups))
TESTS.append(("robustness: direction holds in every specification", test_robust_more_objects))
TESTS.append(("robustness: automatic threshold, multiplier, and factor", test_robust_automatic_threshold))
TESTS.append(("robustness: direction reverses above a known threshold", test_robust_sign_flip))
TESTS.append(("robustness: units, ROI files, and a non-rectangular ROI", test_robust_units_and_rois))
TESTS.append(("robustness: small designs, several groups, cancel", test_robust_small_designs))
TESTS.append(("robustness: area fraction and a preparation setting", test_robust_area_and_preparation))
TESTS.append(("robustness: quick object count equals analyze_objects", test_robust_fast_count))
TESTS.append(("robustness: 8-bit, micron override, two ROIs on the border", test_robust_micron_byte_rois))
TESTS.append(("robustness: two frames and Welch with fractional df", test_robust_two_frames))
TESTS.append(("robustness: light background moves the upper level", test_robust_light_background))
TESTS.append(("robustness: profiles and isolation from a normal run", test_robust_settings_and_isolation))
TESTS.append(("robustness: group table names and repeated rows", test_robust_table_names))
TESTS.append(("robustness: 32-bit mean object area in an oval ROI", test_robust_object_area_float))
TESTS.append(("robustness: tables from hand-made measurements", test_robust_tables_from_records))
TESTS.append(("robustness: threshold factor on a Manual level of 0", test_robust_manual_zero_level))

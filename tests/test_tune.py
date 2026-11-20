# Cell Quant Studio Universal - tests of tuning settings from point marks
#
# Loaded by run_tests.py, which provides the helpers used here.
# The source is kept plain ASCII with no encoding declaration on purpose.
#
# The synthetic fields (240 x 200 pixels, background 100):
#   "nuclei"       bright disks of radius 6 (value 2000), one mark at each centre
#   "distractors"  dim disks of the same size (value 600), never marked
#   "specks"       2 x 2 pixel squares (value 2000, 4 pixels), never marked
# With a Manual threshold L (a pixel is foreground when its value is >= L)
# and a minimum object area M in pixels, by construction:
#   nuclei are detected       when L <= 2000 and M <= pixels of a disk
#   distractors are detected  when L <= 600  and M <= pixels of a disk
#   specks are detected       when L <= 2000 and M <= 4
# so every detection is a nucleus (F1 = 1) exactly when 600 < L <= 2000 and
# 4 < M <= pixels of a disk. Outside that window the number of false
# positives is the number of distractors and/or specks, counted below.

from ij.io import RoiEncoder


TUNE_FIELD_A = {
    "nuclei": [(30, 30), (70, 30), (110, 30), (150, 30), (190, 30), (30, 80)],
    "distractors": [(70, 80), (110, 80)],
    "specks": [(150, 80), (190, 80), (50, 140)]
}
# Field B is scored inside the triangle (100,30) (170,170) (30,170), whose
# half-width at height y is (y - 30) / 2. Five nuclei, one distractor and two
# specks lie wholly inside it; three nuclei and one distractor wholly outside.
TUNE_FIELD_B = {
    "nuclei": [(100, 70), (100, 100), (100, 140), (70, 150), (130, 150), (30, 40), (200, 60), (210, 150)],
    "distractors": [(100, 120), (200, 120)],
    "specks": [(85, 160), (115, 160)]
}
# Field C has three nuclei cut by the image border (left, right, bottom).
# A cut disk keeps more than 60 of its pixels, so it is still an object.
TUNE_FIELD_C = {
    "nuclei": [(3, 50), (236, 100), (120, 196), (60, 60), (120, 100), (180, 150)],
    "distractors": [(60, 150)],
    "specks": [(200, 30), (20, 180)]
}
# Field D holds no nucleus and no mark.
TUNE_FIELD_D = {
    "nuclei": [],
    "distractors": [(50, 50), (150, 100)],
    "specks": [(100, 150), (200, 40)]
}


def tune_plane(layout, nucleus_value=2000, distractor_value=600):
    ip = blank_plane(100)
    for center in layout["distractors"]:
        draw_disk(ip, center[0], center[1], DISK_RADIUS, distractor_value)
    for center in layout["nuclei"]:
        draw_disk(ip, center[0], center[1], DISK_RADIUS, nucleus_value)
    for corner in layout["specks"]:
        draw_rect(ip, (corner[0], corner[1], 2, 2), 2000)
    return ip


def tune_image(title, plane):
    """Two channels like make_image: C1 the tuning field, C2 the marker."""
    stack = ImageStack(WIDTH, HEIGHT)
    stack.addSlice("", plane)
    stack.addSlice("", marker_plane())
    imp = ImagePlus(title, stack)
    imp.setDimensions(2, 1, 1)
    imp = CompositeImage(imp, CompositeImage.GRAYSCALE)
    cal = imp.getCalibration()
    cal.pixelWidth = PIXEL
    cal.pixelHeight = PIXEL
    cal.setUnit("um")
    return imp


def tune_config(E, overrides=None):
    values = {"analysis_mode": "Objects only", "size_units": "Pixel^2", "manual_low": 300, "min_size": 0}
    if overrides:
        values.update(overrides)
    return nuclei_config(E, values)


def tune_point_roi(points):
    roi = PointRoi(points[0][0], points[0][1])
    for point in points[1:]:
        roi.addPoint(point[0], point[1])
    return roi


def tune_standard_fields():
    """The four fields above as dictionaries for run_tuning. The caller
    closes the images."""
    fields = []
    layouts = [("A", TUNE_FIELD_A), ("B", TUNE_FIELD_B), ("C", TUNE_FIELD_C), ("D", TUNE_FIELD_D)]
    for entry in layouts:
        fields.append({
            "name": entry[0], "source": tune_image(entry[0] + ".tif", tune_plane(entry[1])),
            "points": list(entry[1]["nuclei"]), "rois": [triangle_roi()] if entry[0] == "B" else None
        })
    return fields


def tune_close_fields(fields):
    for field in fields:
        if field.get("source") is not None:
            field["source"].close()


def tune_run(E, fields, config, space, out_root, options=None, engine=None, all_configs=None):
    if all_configs is None:
        all_configs = [config, marker_config(E)]
    if engine is None:
        engine = E.AnalysisEngine()
    return E.run_tuning(engine, fields, config, all_configs, E.default_global_options(), space, out_root, options)


def tune_csv_rows(path):
    """Read one of the tuning tables back as a list of dictionaries."""
    lines = [line for line in read_utf8(path).split("\n") if line.strip() != ""]
    headers = [cell.strip('"') for cell in lines[0].split('","')]
    rows = []
    for line in lines[1:]:
        cells = [cell.strip('"') for cell in line.split('","')]
        rows.append(dict(zip(headers, cells)))
    return rows


TUNE_STANDARD_SPACE = {"manual_low": [300, 1000, 1500], "min_size": [0, 20, 60]}


def test_tune_known_window(E, t, tmp):
    """Four fields with one shared window of good settings. Expected counts
    (true positives, false positives) per field, from the layouts above:
        settings                 A        B (triangle)   C        D
        L=300,  M=0              6, 5     5, 3           6, 3     0, 4
        L=300,  M=20 or 60       6, 2     5, 1           6, 1     0, 2
        L>=1000, M=0             6, 3     5, 2           6, 2     0, 2
        L>=1000, M=20 or 60      6, 0     5, 0           6, 0     0, 0
    Marks inside a parent ROI: 6 + 5 + 6 + 0 = 17 (three marks of B lie
    outside the triangle and are ignored)."""
    fields = tune_standard_fields()
    try:
        config = tune_config(E)
        result = tune_run(E, fields, config, TUNE_STANDARD_SPACE, os.path.join(tmp, "tune_window"))
        t.equal(result["cancelled"], False, "the run finished")
        t.check(os.path.basename(result["run_dir"]).startswith("CellQuant_Tuning_"), "the run folder is named CellQuant_Tuning_<timestamp>")
        for name in ["tuning_grid.csv", "tuning_field_scores.csv", "tuning_folds.csv", "tuning_stability.csv",
                     "tuning_summary.csv", "tuning_fields.csv", "tuned_profile.properties",
                     "settings.properties", "tuning_log.txt"]:
            t.check(os.path.isfile(os.path.join(result["run_dir"], name)), name + " is written")

        # Four combinations have F1 = 1; all change both settings. The tie
        # goes to the smallest relative change from L=300, M=0: L=1000
        # (700/300 against 1200/300) and M=20 (20/60 against 60/60).
        t.equal(result["tuned_values"], {"manual_low": 1000.0, "min_size": 20.0}, "the tuned values lie inside the window, nearest the start")
        t.equal(result["design"], "leave-one-field-out cross-validation", "four fields give leave-one-field-out")
        summary = result["summary_row"]
        t.equal(summary["Fields"], 4, "four fields are used")
        t.equal(summary["Combinations"], 9, "the grid has 3 x 3 combinations")
        t.equal(summary["Marks_Inside_ROIs"], 17, "marks inside a parent ROI")
        t.equal(summary["Marks_Outside_ROIs"], 3, "marks outside every parent ROI are reported")
        t.equal(summary["Folds_Scored"], 4, "every fold is scored")
        t.equal(summary["HeldOut_True_Positive_Objects"], 17, "held-out true positives")
        t.equal(summary["HeldOut_False_Positive_Objects"], 0, "held-out false positives")
        t.equal(summary["HeldOut_Missed_Annotations"], 0, "held-out missed marks")
        t.close(summary["HeldOut_F1"], 1.0, "pooled held-out F1 is 1 by construction")
        t.close(summary["Resubstitution_Optimistic_F1"], 1.0, "resubstitution F1 is 1 by construction")
        # Starting settings L=300, M=0: TP 17, FP 5 + 3 + 3 + 4 = 15.
        t.equal(summary["AllFields_Baseline_False_Positive_Objects"], 15, "baseline false positives")
        t.close(summary["AllFields_Baseline_Precision"], 17.0 / 32.0, "baseline precision")
        t.close(summary["AllFields_Baseline_F1"], 34.0 / 49.0, "baseline F1 = 2 TP / (TP + FP + marks)")
        t.close(summary["HeldOut_Baseline_F1"], 34.0 / 49.0, "the baseline on the held-out fields is the same number")
        t.equal(summary["HeldOut_Fields"], 4, "each field is held out once")
        t.close(summary["HeldOut_F1_Field_Min"], 1.0, "lowest per-field held-out F1 (field D has none)")

        expected = {
            (300.0, 0.0): (17, 15), (300.0, 20.0): (17, 6), (300.0, 60.0): (17, 6),
            (1000.0, 0.0): (17, 9), (1000.0, 20.0): (17, 0), (1000.0, 60.0): (17, 0),
            (1500.0, 0.0): (17, 9), (1500.0, 20.0): (17, 0), (1500.0, 60.0): (17, 0)
        }
        t.equal(len(result["grid_rows"]), 9, "one grid row per combination")
        for row in result["grid_rows"]:
            key = (row["P_manual_low"], row["P_min_size"])
            t.equal((row["True_Positive_Objects"], row["False_Positive_Objects"]), expected[key],
                    "pooled counts of L=" + str(key[0]) + " M=" + str(key[1]))
            t.equal(row["Status"], "OK", "combination " + str(row["Combination"]) + " is usable")
        ranks = sorted([row["Rank"] for row in result["grid_rows"]])
        t.equal(ranks, list(range(1, 10)), "every usable combination has a distinct rank")
        selected = [row for row in result["grid_rows"] if row["Selected_Final"] == "true"]
        t.equal(len(selected), 1, "exactly one combination is selected")
        t.equal(selected[0]["Rank"], 1, "the selected combination has rank 1")
        t.equal(selected[0]["Folds_Selected"], 4, "the same combination wins every fold")
        start = [row for row in result["grid_rows"] if row["Is_Start"] == "true"]
        t.equal([(row["P_manual_low"], row["P_min_size"]) for row in start], [(300.0, 0.0)], "the starting combination is flagged")

        folds = result["fold_rows"]
        t.equal([row["Held_Out_Field"] for row in folds], ["A", "B", "C", "D"], "one fold per field, in order")
        t.equal([row["Manual_Count"] for row in folds], [6, 5, 6, 0], "held-out marks per field")
        t.equal([row["False_Positive_Objects"] for row in folds], [0, 0, 0, 0], "held-out false positives per field")
        t.equal([row["Baseline_False_Positive_Objects"] for row in folds], [5, 3, 3, 4], "baseline false positives per field")
        t.close(folds[0]["F1"], 1.0, "held-out F1 of field A")
        t.check(folds[3]["F1"] != folds[3]["F1"], "a field without marks or detections has no F1")
        t.equal(folds[1]["Training_Marks"], 12, "training marks when B is held out")
        stable = {}
        for row in result["stability_rows"]:
            stable[(row["Parameter"], row["Value"])] = (row["Folds_Selected"], row["Folds_Scored"], row["Selected_Final"])
        t.equal(stable[("manual_low", "1000")], (4, 4, "true"), "the threshold 1000 is selected in every fold")
        t.equal(stable[("min_size", "20")], (4, 4, "true"), "the minimum area 20 is selected in every fold")
        t.equal(stable[("manual_low", "300")], (0, 4, "false"), "the starting threshold is never selected")

        # Field report and file contents.
        by_name = {}
        for row in result["field_rows"]:
            by_name[row["Field"]] = row
        t.equal(by_name["B"]["Marks_Total"], 8, "field B has eight marks")
        t.equal(by_name["B"]["Marks_Inside_ROIs"], 5, "five of them inside the triangle")
        t.equal(by_name["B"]["Marks_Outside_ROIs"], 3, "three outside every ROI")
        t.equal(by_name["D"]["Marks_Total"], 0, "field D has no mark")
        saved = tune_csv_rows(os.path.join(result["run_dir"], "tuning_summary.csv"))
        t.equal(len(saved), 1, "the summary file has one row")
        t.equal(saved[0]["Run_ID"], result["run_id"], "tables carry the run identifier")
        t.equal(saved[0]["Application_Version"], E.APP_VERSION, "tables carry the application version")
        t.equal(saved[0]["Validation_Design"], "leave-one-field-out cross-validation", "the design is written to the summary")
        t.equal(len(tune_csv_rows(os.path.join(result["run_dir"], "tuning_field_scores.csv"))), 36, "one score row per field and combination")
        log_text = read_utf8(os.path.join(result["run_dir"], "tuning_log.txt"))
        t.check("HELD-OUT result" in log_text and "optimistic" in log_text, "the log separates the held-out from the optimistic score")

        # Determinism: a second run gives the same tables apart from Run_ID and timing.
        again = tune_run(E, fields, config, TUNE_STANDARD_SPACE, os.path.join(tmp, "tune_window"))
        t.check(again["run_dir"] != result["run_dir"], "a second run gets its own folder")
        for index in range(9):
            for key in E.TUNE_GRID_HEADERS:
                if key != "Run_ID":
                    first = result["grid_rows"][index].get(key)
                    second = again["grid_rows"][index].get(key)
                    t.check(first == second or (first != first and second != second), "repeat run: grid column " + key + " is identical")

        # Tie-breaking: a start inside the window is kept unchanged, and
        # fewer changed settings beat a smaller relative change.
        inside = tune_run(E, fields, tune_config(E, {"manual_low": 1500, "min_size": 60}), TUNE_STANDARD_SPACE, os.path.join(tmp, "tune_window"))
        t.equal(inside["tuned_values"], {"manual_low": 1500.0, "min_size": 60.0}, "a starting configuration that is already best is kept")
        fewer = tune_run(E, fields, tune_config(E, {"manual_low": 1500, "min_size": 0}), TUNE_STANDARD_SPACE, os.path.join(tmp, "tune_window"))
        t.equal(fewer["tuned_values"], {"manual_low": 1500.0, "min_size": 20.0}, "ties go to the fewest changed settings")

        # The saved profile reproduces the reported score through the
        # ordinary validation tool, and leaves the other channel alone.
        options, enabled, everything = E.load_profile_settings(result["profile_path"])
        t.equal(len(everything), 2, "the tuned profile holds both channels")
        for key in result["tuned_config"]:
            t.equal(enabled[0][key], result["tuned_config"][key], "tuned profile C1 '" + key + "'")
        marker = marker_config(E)
        for key in marker:
            t.equal(everything[1][key], marker[key], "tuned profile leaves C2 '" + key + "' unchanged")
        t.equal(enabled[0]["manual_low"], 1000.0, "the tuned threshold is in the profile")
        engine = E.AnalysisEngine()
        totals = {"True_Positive_Objects": 0, "False_Positive_Objects": 0, "Missed_Annotations": 0, "Manual_Count": 0}
        for field in fields:
            prepared = engine.prepare_single_channel(field["source"], enabled[0], 1)
            rois = [E.whole_image_roi(prepared)] if field["rois"] is None else field["rois"]
            check = E.validate_objects_against_points(
                engine, prepared, enabled[0], rois, field["points"], E.unique_dir(tmp, "tune_recheck"), field["name"], 1
            )
            for key in totals:
                totals[key] += check["rows"][0][key]
            prepared.close()
        for key in totals:
            t.equal(totals[key], summary["Resubstitution_Optimistic_" + key], "the reloaded profile reproduces " + key)
    finally:
        tune_close_fields(fields)


def test_tune_held_out_is_honest(E, t, tmp):
    """One field needs a different threshold, so settings tuned without it
    fail on it. Fields A2 and B2: 8 nuclei (2000) and 8 distractors (600).
    Field C2: 8 dim nuclei (500), no distractor. Specks (2000) in every
    field keep the threshold 1000 valid and are removed by M = 20.
    Grid L in {400, 1000}; counts (TP, FP, missed) by construction:
                 A2 or B2      C2
        L=400    8, 8, 0       8, 0, 0
        L=1000   8, 0, 0       0, 0, 8
    All three fields: L=1000 gives F1 = 2*16/(16+0+24) = 0.8, L=400 gives
    2*24/(24+16+24) = 0.75, so the final setting is 1000 with an optimistic
    F1 of 0.8. Leaving A2 out: on B2+C2, L=400 scores 2*16/(16+8+16) = 0.8
    against 2*8/(8+0+16) = 0.667, so 400 is chosen and A2 gets TP 8, FP 8;
    the same for B2. Leaving C2 out: L=1000 is perfect on A2+B2 and misses
    all 8 marks of C2. Pooled held-out: TP 16, FP 16, marks 24, so
    precision 1/2, recall 2/3, F1 = 2*16/(16+16+24) = 4/7."""
    nuclei = DISKS[0:8]
    distractors = DISKS[8:15] + [(30, 175)]
    specks = [(200, 180), (10, 110)]
    fields = []
    for name in ["A2", "B2"]:
        plane = tune_plane({"nuclei": nuclei, "distractors": distractors, "specks": specks})
        fields.append({"name": name, "source": tune_image(name + ".tif", plane), "points": list(nuclei)})
    plane = tune_plane({"nuclei": nuclei, "distractors": [], "specks": specks}, nucleus_value=500)
    fields.append({"name": "C2", "source": tune_image("C2.tif", plane), "points": list(nuclei)})
    try:
        config = tune_config(E, {"manual_low": 1000, "min_size": 20})
        result = tune_run(E, fields, config, {"manual_low": [400, 1000]}, os.path.join(tmp, "tune_honest"))
        summary = result["summary_row"]
        t.equal(result["tuned_values"], {"manual_low": 1000.0}, "the final setting is the best compromise on all fields")
        t.close(summary["Resubstitution_Optimistic_F1"], 0.8, "optimistic score of the final setting")
        t.equal(summary["HeldOut_True_Positive_Objects"], 16, "held-out true positives")
        t.equal(summary["HeldOut_False_Positive_Objects"], 16, "held-out false positives")
        t.equal(summary["HeldOut_Missed_Annotations"], 8, "held-out missed marks")
        t.close(summary["HeldOut_Precision"], 0.5, "held-out precision")
        t.close(summary["HeldOut_Recall"], 2.0 / 3.0, "held-out recall")
        t.close(summary["HeldOut_F1"], 4.0 / 7.0, "held-out F1")
        t.check(summary["HeldOut_F1"] < summary["Resubstitution_Optimistic_F1"], "the held-out score is lower than the optimistic one")
        t.close(summary["HeldOut_Baseline_F1"], 0.8, "the untuned start on the same fields")
        t.close(summary["HeldOut_F1_Field_Min"], 2.0 / 3.0, "worst scored field: TP 8, FP 8 gives 2/3 (C2 has no F1)")
        t.equal([row["P_manual_low"] for row in result["fold_rows"]], [400.0, 400.0, 1000.0], "the setting chosen in each fold")
        t.equal([row["Missed_Annotations"] for row in result["fold_rows"]], [0, 0, 8], "missed marks per held-out field")
        t.close(result["fold_rows"][0]["Training_F1_Optimistic"], 0.8, "training score of the first fold")
        stable = {}
        for row in result["stability_rows"]:
            stable[row["Value"]] = (row["Folds_Selected"], row["Selected_Final"])
        t.equal(stable["400"], (2, "false"), "threshold 400 wins two of three folds but is not the final choice")
        t.equal(stable["1000"], (1, "true"), "threshold 1000 wins one fold and is the final choice")
        t.close([row for row in result["stability_rows"] if row["Value"] == "400"][0]["Fraction_Of_Folds"], 2.0 / 3.0, "stability as a fraction")

        # The count objective: sum over fields of |detected - marks|.
        # L=400: 8 + 8 + 0 = 16. L=1000: 0 + 0 + 8 = 8, so 1000 is chosen.
        counted = tune_run(E, fields, config, {"manual_low": [400, 1000]}, os.path.join(tmp, "tune_honest"),
                           {"objective": E.TUNE_OBJECTIVE_COUNT})
        errors = dict([(row["P_manual_low"], row["Abs_Count_Error"]) for row in counted["grid_rows"]])
        t.equal(errors, {400.0: 16, 1000.0: 8}, "absolute count error per combination")
        t.equal(counted["tuned_values"], {"manual_low": 1000.0}, "the count objective picks the smaller count error")
        t.close(counted["grid_rows"][0]["Abs_Count_Error_Percent"], 100.0 * 16 / 24, "count error as a percentage of the marks")

        # Two fields: a two-way split. One field: no held-out score at all.
        two = tune_run(E, fields[0:2], config, {"manual_low": [400, 1000]}, os.path.join(tmp, "tune_honest"))
        t.equal(two["design"], "two-way split (each field held out once)", "two fields give a two-way split")
        t.close(two["summary_row"]["HeldOut_F1"], 1.0, "two like fields: held-out F1 is 1")
        one = tune_run(E, fields[0:1], config, {"manual_low": [400, 1000]}, os.path.join(tmp, "tune_honest"))
        t.equal(one["design"], "none: tuned and scored on the same field (optimistic)", "one field is labelled as optimistic")
        t.check(one["held_out"] is None, "one field gives no held-out result")
        t.check(one["summary_row"].get("HeldOut_F1", "") == "", "one field leaves the held-out columns empty")
        t.equal(one["summary_row"]["Folds"], 0, "one field has no fold")
        t.close(one["summary_row"]["Resubstitution_Optimistic_F1"], 1.0, "its optimistic score is still reported")
        t.check("No held-out estimate" in one["summary_row"]["Note"], "the summary says there is no held-out estimate")
        t.check(os.path.isfile(one["profile_path"]), "a profile is still written")

        # Fields that belong together are held out together.
        for index in range(3):
            fields[index]["group"] = "organoid_1" if index < 2 else "organoid_2"
        grouped = tune_run(E, fields, config, {"manual_low": [400, 1000]}, os.path.join(tmp, "tune_honest"))
        t.equal(grouped["design"], "leave-one-group-out cross-validation (2 groups)", "groups give leave-one-group-out")
        t.equal(grouped["summary_row"]["Folds"], 2, "one fold per group")
        # Leaving organoid_1 out: tuned on C2 alone (400), scored on A2+B2: TP 16, FP 16.
        # Leaving organoid_2 out: tuned on A2+B2 (1000), scored on C2: 8 missed.
        t.equal(grouped["summary_row"]["HeldOut_False_Positive_Objects"], 16, "grouped held-out false positives")
        t.equal(grouped["summary_row"]["HeldOut_Missed_Annotations"], 8, "grouped held-out missed marks")
    finally:
        tune_close_fields(fields)


def test_tune_annotated_folder(E, t, tmp):
    """The folder convention: <image>_points.zip or .roi beside the image,
    parent ROIs from the usual ROI file, and a report of what is missing."""
    folder = os.path.join(tmp, "tune_folder")
    os.makedirs(folder)
    layouts = {"a": TUNE_FIELD_A, "b": TUNE_FIELD_B, "c": TUNE_FIELD_C, "d": TUNE_FIELD_D, "e": TUNE_FIELD_A, "f": TUNE_FIELD_C}
    for name in sorted(layouts.keys()):
        imp = tune_image(name + ".tif", tune_plane(layouts[name]))
        FileSaver(imp).saveAsTiff(os.path.join(folder, name + ".tif"))
        imp.close()
    # a: one multi-point ROI saved as .roi, no parent ROI file
    RoiEncoder.save(tune_point_roi(TUNE_FIELD_A["nuclei"]), os.path.join(folder, "a_points.roi"))
    # b: triangle as parent ROI file, marks as a ROI Manager set of two point ROIs
    E.save_roi_array([triangle_roi()], os.path.join(folder, "b_ROIs.zip"))
    E.save_roi_array([tune_point_roi(TUNE_FIELD_B["nuclei"][0:3]), tune_point_roi(TUNE_FIELD_B["nuclei"][3:])],
                     os.path.join(folder, "b_points.zip"))
    # c: no points file, but the ROI file holds a rectangle and the marks
    whole = Roi(0, 0, WIDTH, HEIGHT)
    whole.setName("all")
    E.save_roi_array([whole, tune_point_roi(TUNE_FIELD_C["nuclei"])], os.path.join(folder, "c_ROIs.zip"))
    # d: a points file without any point ROI means "annotated, zero marks"
    E.save_roi_array([Roi(0, 0, 10, 10)], os.path.join(folder, "d_points.zip"))
    # e: no annotation at all.  f: a damaged points file.
    handle = open(os.path.join(folder, "f_points.zip"), "wb")
    try:
        handle.write(b"not a ROI file")
    finally:
        handle.close()

    t.equal(os.path.basename(E.find_points_file(os.path.join(folder, "a.tif"))), "a_points.roi", "the points file is found beside the image")
    t.check(E.find_points_file(os.path.join(folder, "e.tif")) is None, "no points file is reported as none")
    t.check(E.find_roi_file(os.path.join(folder, "a.tif")) is None, "a points file is not mistaken for a parent ROI file")
    found = E.discover_annotated_fields(folder, "tif")
    status = dict([(row["Field"], row["Status"]) for row in found["report"]])
    t.equal(status, {"a.tif": "ANNOTATED", "b.tif": "ANNOTATED", "c.tif": "ANNOTATED", "d.tif": "ANNOTATED",
                     "e.tif": "NO_ANNOTATION", "f.tif": "FAILED"}, "every image is reported with its status")
    t.equal([os.path.basename(path) for path in found["missing"]], ["e.tif"], "the image without marks is listed as missing")
    t.equal([field["name"] for field in found["fields"]], ["a.tif", "b.tif", "c.tif", "d.tif"], "annotated fields, sorted by name")
    t.equal([len(field["points"]) for field in found["fields"]], [6, 8, 6, 0], "marks read per field")
    t.equal(sorted(found["fields"][1]["points"]), sorted(TUNE_FIELD_B["nuclei"]), "marks of a ROI Manager set keep their positions")
    t.check(found["fields"][0]["rois"] is None, "no ROI file means the whole image")
    t.equal(len(found["fields"][1]["rois"]), 1, "the parent ROI file is used")
    t.equal(len(found["fields"][2]["rois"]), 1, "point ROIs in the ROI file are not parent ROIs")
    t.check("c_ROIs.zip" in found["fields"][2]["points_source"], "the source of the marks is recorded")
    grouped = E.discover_annotated_fields(folder, "tif", False, None, None, E.TUNE_GROUP_FOLDER)
    t.equal(set([field["group"] for field in grouped["fields"]]), set(["tune_folder"]), "grouping by subfolder uses the folder name")

    # The whole route from a folder and the settings of a profile gives the
    # numbers of test_tune_known_window (the same four fields).
    config = tune_config(E)
    options = E.default_global_options()
    t.equal(options["tune_objective"], E.TUNE_OBJECTIVE_F1, "the tuning settings have defaults")
    options["tune_search_space"] = "manual_low = 300, 1000, 1500; min_size = 0, 20, 60"
    options["batch_extensions"] = "tif"
    messages = []
    engine = E.AnalysisEngine(log=messages.append)
    progress = []
    result = E.run_tuning_folder(
        engine, folder, os.path.join(tmp, "tune_folder_out"), [config, marker_config(E)], [config, marker_config(E)],
        options, None, lambda done, total, name: progress.append((done, total))
    )
    summary = result["summary_row"]
    t.equal(result["tuned_values"], {"manual_low": 1000.0, "min_size": 20.0}, "the folder route finds the same settings")
    t.equal(summary["Fields"], 4, "four annotated fields are used")
    t.equal(summary["Images_Without_Annotation"], 2, "images that could not be used are counted")
    t.equal(summary["Marks_Inside_ROIs"], 17, "marks inside parent ROIs")
    t.equal(summary["Marks_Outside_ROIs"], 3, "marks outside the triangle are ignored and counted")
    t.close(summary["HeldOut_F1"], 1.0, "held-out F1 from files")
    t.close(summary["AllFields_Baseline_F1"], 34.0 / 49.0, "baseline F1 from files")
    t.equal(summary["Search_Space"], "manual_low=300,1000,1500; min_size=0,20,60", "the search space is recorded")
    rows = tune_csv_rows(os.path.join(result["run_dir"], "tuning_fields.csv"))
    t.equal(sorted([row["Status"] for row in rows]), ["FAILED", "NO_ANNOTATION", "USED", "USED", "USED", "USED"], "the field table lists unused images too")
    t.equal(progress[-1], (40, 40), "progress ends at fields x (combinations + baseline)")
    t.check(len([m for m in messages if "Not used for tuning: e.tif" in E.to_text(m)]) == 1, "the missing annotation is logged")
    tuned_options, enabled, everything = E.load_profile_settings(result["profile_path"])
    t.equal(tuned_options["tune_search_space"], options["tune_search_space"], "the tuning settings are kept in the profile")

    empty = os.path.join(tmp, "tune_folder_empty")
    os.makedirs(empty)
    imp = tune_image("x.tif", tune_plane(TUNE_FIELD_A))
    FileSaver(imp).saveAsTiff(os.path.join(empty, "x.tif"))
    imp.close()
    refused = ""
    try:
        E.run_tuning_folder(engine, empty, os.path.join(tmp, "tune_folder_out"), [config], [config], options)
    except ValueError as exc:
        refused = E.to_text(exc)
    t.check("_points.zip" in refused, "a folder without annotation is refused with the naming rule")


def test_tune_search_space_rules(E, t, tmp):
    config = tune_config(E)
    space = E.normalize_tune_search_space(config, {"min_size": ["20", 20, 40.0], "manual_low": [800], "watershed": ["true", False]})
    t.equal(space, [("manual_low", [800.0]), ("watershed", [True, False]), ("min_size", [20.0, 40.0])],
            "values are converted, repeated values dropped, settings put in a fixed order")
    t.equal(E.parse_tune_search_space("manual_low = 800, 1000;\n watershed=true,false"),
            {"manual_low": ["800", "1000"], "watershed": ["true", "false"]}, "the text form is parsed")
    t.equal(E.parse_tune_search_space(E.describe_tune_search_space(space)),
            {"manual_low": ["800"], "watershed": ["true", "false"], "min_size": ["20", "40"]}, "the description reads back")

    combos = E.tune_build_combinations(config, E.normalize_tune_search_space(config, {"manual_low": [300, 900], "min_size": [0, 30, 60]}))
    t.equal(len(combos), 7, "2 x 3 combinations and the baseline")
    t.equal([(c["values"]["manual_low"], c["values"]["min_size"]) for c in combos[0:6]],
            [(300.0, 0.0), (300.0, 30.0), (300.0, 60.0), (900.0, 0.0), (900.0, 30.0), (900.0, 60.0)], "grid order: first setting slowest")
    t.equal([c["changes"] for c in combos[0:6]], [0, 1, 1, 1, 2, 2], "number of settings changed from the start")
    t.close(combos[4]["distance"], 2.0 + 0.5, "relative change: 600/300 plus 30/60 (a zero start is scaled by the largest value)")
    t.check(combos[6]["baseline"] and combos[6]["changes"] == 0, "the baseline is the starting configuration")

    def refused(search_space, start=None, limit=None):
        try:
            E.normalize_tune_search_space(start if start is not None else config, search_space, limit)
        except ValueError as exc:
            return E.to_text(exc)
        return ""

    big = {"manual_low": list(range(1, 21)), "min_size": list(range(1, 21)), "max_size": list(range(100, 120))}
    message = refused(big)
    t.check("8000" in message and "2000" in message and "Shrink" in message, "a grid above the limit is refused with both numbers")
    t.check("limit is 5" in refused({"manual_low": [1, 2, 3], "min_size": [1, 2]}, None, 5), "the limit can be set")
    t.check(refused({"manual_low": [1, 2, 3], "min_size": [1, 2]}, None, 6) == "", "a grid at the limit is accepted")
    t.check("cannot be tuned" in refused({"lut": ["Red"]}), "an unknown setting is refused")
    t.check("Auto" in refused({"auto_scale": [0.8, 1.2]}), "a setting of another threshold mode is refused")
    t.check("Manual" in refused({"manual_low": [500]}, tune_config(E, {"threshold_mode": "Auto"})), "manual_low is refused in Auto mode")
    t.check("Display-minimum" in refused({"min_size": [5]}, tune_config(E, {"threshold_mode": "Display minimum"})), "a display-minimum threshold is refused")
    t.check("Objects" in refused({"min_size": [5]}, tune_config(E, {"analysis_mode": "Positive area"})), "a channel without objects is refused")
    t.check("method" in refused({"auto_method": ["NoSuchMethod"]}, tune_config(E, {"threshold_mode": "Auto"})), "an unknown threshold method is refused")
    t.check("numeric" in refused({"min_size": ["abc"]}), "a value that is not a number is refused")
    t.check("empty" in refused({}), "an empty search space is refused")
    t.check("CLAHE is off" in refused({"clahe_max_slope": [2, 3]}, tune_config(E, {"threshold_mode": "Local Phansalkar"})),
            "CLAHE slope without CLAHE is refused")
    t.check(refused({"clahe_max_slope": [2, 3], "clahe": [False, True]}, tune_config(E, {"threshold_mode": "Local Phansalkar"})) == "",
            "CLAHE slope with CLAHE in the search is accepted")

    default = E.default_tune_search_space(tune_config(E, {"manual_low": 1000, "min_size": 40}))
    t.equal(default["manual_low"], [700.0, 850.0, 1000.0, 1150.0, 1300.0], "default grid: threshold times the factors")
    t.equal(default["min_size"], [28.0, 34.0, 40.0, 46.0, 52.0], "default grid: minimum area times the factors")
    t.equal(default["watershed"], [False, True], "default grid: watershed off and on, start first")
    t.check("auto_scale" not in default and "local_radius" not in default, "default grid: only the threshold setting of the mode in use")
    t.check("min_size" not in E.default_tune_search_space(config), "default grid: a setting at 0 is left out")
    local = E.default_tune_search_space(tune_config(E, {"threshold_mode": "Local Phansalkar", "local_radius": 3}), "0.9, 1.1")
    t.equal(local["local_radius"], [3], "default grid: radii are whole numbers without repeats")
    t.equal(E.tune_config_from_options({"tune_channel": "0"}, [marker_config(E), config])["channel"], 1, "channel 0 means the first object channel")
    try:
        E.tune_config_from_options({"tune_channel": "3"}, [config])
        t.check(False, "a channel that is not enabled is refused")
    except ValueError:
        t.check(True, "a channel that is not enabled is refused")


def test_tune_size_shortcut_and_failures(E, t, tmp):
    """The size-limit shortcut gives exactly what a direct analysis gives,
    also at the limits themselves, with calibrated areas, a non-rectangular
    ROI and objects cut by the image border; and combinations that cannot
    be run are reported, not selected."""
    disk = disk_pixel_count()
    area = PIXEL * PIXEL
    engine = E.AnalysisEngine()
    # Calibrated units: one pixel is 0.25 um^2, so a limit of n * 0.25 is n pixels exactly.
    config = tune_config(E, {"manual_low": 1000, "size_units": "Calibrated area"})
    sizes = {"min_size": [0, 4 * area, 5 * area, disk * area, (disk + 1) * area],
             "max_size": [-1, disk * area, (disk - 1) * area, 3 * area]}
    space = E.normalize_tune_search_space(config, sizes)
    combos = E.tune_build_combinations(config, space)
    layouts = [("A", TUNE_FIELD_A, None), ("B", TUNE_FIELD_B, [triangle_roi()]), ("C", TUNE_FIELD_C, None)]
    for entry in layouts:
        source = tune_image(entry[0] + ".tif", tune_plane(entry[1]))
        field = {"name": entry[0], "source": source, "points": list(entry[1]["nuclei"]), "rois": entry[2]}
        work = E.tune_open_field(engine, field, config, False, None, False, None)
        try:
            fast = E.tune_score_field(engine, work, combos, True)
            direct = E.tune_score_field(engine, work, combos, False)
        finally:
            E.tune_close_field(work)
            source.close()
        t.equal(fast, direct, "field " + entry[0] + ": shortcut and direct analysis agree for all " + str(len(combos)) + " combinations")
        if entry[0] == "A":
            # 6 whole disks and 3 specks pass the threshold 1000.
            by_values = {}
            for index in range(len(combos) - 1):
                by_values[(combos[index]["values"]["min_size"], combos[index]["values"]["max_size"])] = fast[index]
            t.equal(by_values[(0.0, -1.0)]["detected"], 9, "no limits: disks and specks")
            t.equal(by_values[(4 * area, -1.0)]["detected"], 9, "a speck of exactly the minimum area is kept")
            t.equal(by_values[(5 * area, -1.0)]["detected"], 6, "one pixel more removes the specks")
            t.equal(by_values[(disk * area, -1.0)]["detected"], 6, "a disk of exactly the minimum area is kept")
            t.equal(by_values[((disk + 1) * area, -1.0)]["detected"], 0, "one pixel more removes the disks")
            t.equal(by_values[(0.0, disk * area)]["detected"], 9, "a disk of exactly the maximum area is kept")
            t.equal(by_values[(0.0, (disk - 1) * area)]["detected"], 3, "one pixel less leaves only the specks")
            t.equal(by_values[(5 * area, 3 * area)].get("error", "")[0:7], "INVALID", "a maximum below the minimum is marked invalid")
            t.equal(by_values[(0.0, 3 * area)]["detected"], 0, "a maximum below every object leaves an empty list")

    # A threshold above the brightest pixel (2000) of every field holds no
    # pixel value: the true answer is "nothing detected", so all 17 marks
    # inside parent ROIs are missed (TP 0, FP 0). It is scored, ranked
    # below the working threshold, and never selected.
    fields = tune_standard_fields()
    try:
        result = tune_run(E, fields, tune_config(E), {"manual_low": [1000, 5000], "min_size": [20]}, os.path.join(tmp, "tune_fail"))
        rows = result["grid_rows"]
        t.equal(rows[0]["Status"], "OK", "the valid threshold is usable")
        t.check(rows[1]["Status"].startswith("OK: nothing detected on 4 field(s)"), "a threshold above the image maximum detects nothing and says so")
        t.equal((rows[1]["Detected_Count"], rows[1]["True_Positive_Objects"], rows[1]["Missed_Annotations"], rows[1]["Manual_Count"]),
                (0, 0, 17, 17), "nothing detected: every mark inside a parent ROI is missed")
        t.equal(rows[1]["Objective_Value"], 0.0, "nothing detected scores F1 objective 0")
        t.equal(rows[1]["Rank"], 2, "the empty threshold range ranks below the working one")
        t.equal(result["summary_row"]["Combinations_Usable"], 2, "both combinations are usable")
        t.equal(result["tuned_values"], {"manual_low": 1000.0, "min_size": 20.0}, "the empty threshold range is not selected")
        empty_rows = [row for row in result["score_rows"] if row["Combination"] == 2]
        t.equal([row["Status"] for row in empty_rows], [E.TUNE_EMPTY_RANGE_STATUS] * 4, "the field scores say why nothing was detected")
        # A grid in which no combination can run at all (every maximum area
        # below its minimum area) stops with a clear message.
        stopped = ""
        try:
            tune_run(E, fields, tune_config(E), {"min_size": [50], "max_size": [10]}, os.path.join(tmp, "tune_fail"))
        except ValueError as exc:
            stopped = E.to_text(exc)
        t.check("No combination could be run on every field" in stopped and "INVALID" in stopped,
                "a grid with no usable combination stops with a clear message")

        # Degenerate inputs.
        no_marks = [{"name": "D", "source": fields[3]["source"], "points": []}]
        stopped = ""
        try:
            tune_run(E, no_marks, tune_config(E), {"min_size": [0, 20]}, os.path.join(tmp, "tune_fail"))
        except ValueError as exc:
            stopped = E.to_text(exc)
        t.check("no mark" in stopped, "tuning without any mark is refused")
        # A tiny parent ROI that holds no object and no mark, next to a field with marks:
        # with A as the only training field the fold of the tiny field is scored
        # (0 detections), the fold of A cannot be tuned (no training marks).
        tiny = Roi(0, 190, 4, 4)
        pair = [{"name": "A", "source": fields[0]["source"], "points": list(TUNE_FIELD_A["nuclei"])},
                {"name": "tiny", "source": fields[0]["source"], "points": [(150, 30)], "rois": [tiny]}]
        result = tune_run(E, pair, tune_config(E), {"manual_low": [300, 1000], "min_size": [0, 20]}, os.path.join(tmp, "tune_fail"))
        t.equal(result["summary_row"]["Marks_Outside_ROIs"], 1, "a mark outside the tiny ROI is ignored and counted")
        t.equal(result["summary_row"]["Folds"], 2, "two folds")
        t.equal(result["summary_row"]["Folds_Scored"], 1, "the fold whose training fields hold no marks is not scored")
        t.check("no marks" in result["fold_rows"][0]["Note"], "the unscored fold says why")
        t.equal(result["fold_rows"][1]["Detected_Count"], 0, "nothing is detected in the tiny ROI")
        t.equal(result["summary_row"]["HeldOut_Fields"], 1, "only scored folds enter the held-out result")
        # A field that cannot be used is reported and the others are still tuned.
        outside = Roi(200, 150, 100, 100)
        mixed = [{"name": "A", "source": fields[0]["source"], "points": list(TUNE_FIELD_A["nuclei"])},
                 {"name": "bad", "source": fields[2]["source"], "points": [(1, 1)], "rois": [outside]}]
        result = tune_run(E, mixed, tune_config(E), {"min_size": [0, 20]}, os.path.join(tmp, "tune_fail"))
        t.equal([row["Status"] for row in result["field_rows"]], ["USED", "FAILED"], "a field with a ROI outside the image is reported as failed")
        t.equal(result["summary_row"]["Fields_Failed"], 1, "failed fields are counted in the summary")
        t.equal(result["summary_row"]["Fields"], 1, "the remaining field is used")
    finally:
        tune_close_fields(fields)


def test_tune_preparation_filters_and_cancel(E, t, tmp):
    fields = tune_standard_fields()
    try:
        engine = E.AnalysisEngine()
        # Re-preparation: an image prepared without the filters, then
        # filtered here, equals a full preparation by the engine.
        full_config = tune_config(E, {"median_radius": 2, "gaussian_sigma": 1.5, "subtract_background": True, "background_radius": 20})
        full = engine.prepare_single_channel(fields[2]["source"], full_config, 1)
        base_config = dict(full_config)
        base_config["median_radius"] = 0
        base_config["gaussian_sigma"] = 0
        base = engine.prepare_single_channel(fields[2]["source"], base_config, 1)
        again = E.tune_apply_filters(base, full_config)
        a = full.getProcessor().getPixels()
        b = again.getProcessor().getPixels()
        differences = len([i for i in range(len(a)) if a[i] != b[i]])
        t.equal(differences, 0, "filtering a prepared image equals the engine's full preparation, pixel for pixel")
        t.check(len([i for i in range(len(a)) if a[i] != base.getProcessor().getPixels()[i]]) > 0, "the filters do change the image")
        t.equal(E.to_text(again.getTitle()), E.to_text(full.getTitle()), "the re-prepared image keeps its title")
        t.equal(again.getCalibration().pixelWidth, full.getCalibration().pixelWidth, "the re-prepared image keeps its calibration")
        full.close()
        base.close()
        again.close()

        # Searching a preparation filter from source images works. A Gaussian
        # blur of sigma 2 spreads a 2 x 2 speck (1900 above background) so
        # that its brightest pixel keeps less than (0.39)^2 = 0.16 of it:
        # below 100 + 300 = 400, under the threshold 1000. The centre of a
        # disk of radius 6 keeps 1 - exp(-36 / 8) = 0.989 of its 1900: it
        # stays above 1000 and still holds its mark. So on fields A, B, C
        # (3 + 2 inside the triangle + 2 = 7 specks, 6 + 5 + 6 = 17 nuclei)
        # only the blurred setting removes the speck detections.
        config = tune_config(E, {"manual_low": 1000, "min_size": 0})
        result = tune_run(E, fields[0:3], config, {"gaussian_sigma": [0, 2]}, os.path.join(tmp, "tune_prep"))
        rows = result["grid_rows"]
        t.equal(rows[0]["False_Positive_Objects"], 7, "no blur: the seven specks inside parent ROIs are detected")
        t.equal(rows[0]["True_Positive_Objects"], 17, "no blur: every nucleus is found")
        t.equal(rows[1]["False_Positive_Objects"], 0, "blur: the specks fall below the threshold")
        t.equal(rows[1]["True_Positive_Objects"], 17, "blur: every nucleus is still found")
        t.equal(result["tuned_values"], {"gaussian_sigma": 2.0}, "the blur is selected")
        t.close(result["summary_row"]["HeldOut_F1"], 1.0, "blur: held-out F1 is 1 by construction")
        options, enabled, everything = E.load_profile_settings(result["profile_path"])
        t.equal(enabled[0]["gaussian_sigma"], 2.0, "the tuned filter is in the profile")

        # Field D holds only distractors (600) and specks: after the blur
        # its brightest pixel (about 100 + 0.989 * 500 = 594) is below the
        # Manual threshold 1000, so the blur detects nothing on D, which is
        # the right answer for D (no nucleus, no mark). The blur therefore
        # still has F1 = 1 on all four fields and is selected; without it
        # the nine specks inside parent ROIs (7 + 2 of D) are false positives.
        messages = []
        logging = E.AnalysisEngine(log=messages.append)
        result = tune_run(E, fields, config, {"gaussian_sigma": [0, 2]}, os.path.join(tmp, "tune_prep"), None, logging)
        rows = result["grid_rows"]
        t.equal(rows[0]["False_Positive_Objects"], 9, "four fields, no blur: nine specks inside parent ROIs")
        t.check(rows[1]["Status"].startswith("OK: nothing detected on 1 field(s), first D"), "the blur detects nothing on field D and says so")
        t.equal((rows[1]["True_Positive_Objects"], rows[1]["False_Positive_Objects"]), (17, 0), "blur on four fields: every nucleus, no speck")
        t.equal(rows[1]["Rank"], 1, "nothing detected on a field without nuclei is a correct result, not a failure")
        t.equal(result["tuned_values"], {"gaussian_sigma": 2.0}, "the blur is selected on four fields too")
        t.close(result["summary_row"]["HeldOut_F1"], 1.0, "four fields with blur: held-out F1 is 1 by construction")
        t.check(len([m for m in messages if "could not be run on every field" in E.to_text(m)]) == 0, "no combination is reported as left out")

        # ... and is refused clearly for a field given as a prepared image.
        prepared = engine.prepare_single_channel(fields[0]["source"], config, 1)
        given = [{"name": "prepared A", "prepared": prepared, "points": list(TUNE_FIELD_A["nuclei"])}]
        result = None
        stopped = ""
        try:
            tune_run(E, given, config, {"gaussian_sigma": [0, 2]}, os.path.join(tmp, "tune_prep"))
        except ValueError as exc:
            stopped = E.to_text(exc)
        t.check("already prepared image" in stopped, "a preparation filter is refused for a prepared image")
        # Other settings work on a prepared image, which stays open.
        result = tune_run(E, given, config, {"min_size": [0, 20]}, os.path.join(tmp, "tune_prep"))
        t.equal(result["tuned_values"], {"min_size": 20.0}, "a prepared image can be tuned")
        t.equal(result["grid_rows"][0]["False_Positive_Objects"], 3, "prepared image: three specks without a size limit")
        t.check(prepared.getProcessor() is not None and prepared.getWidth() == WIDTH, "the caller's prepared image is left open")
        prepared.close()

        # Cancel: the run stops, says so, and writes no profile.
        calls = {"count": 0}

        def cancel_after_some():
            calls["count"] += 1
            return calls["count"] > 6

        cancelling = E.AnalysisEngine(is_cancelled=cancel_after_some)
        result = tune_run(E, fields, tune_config(E), TUNE_STANDARD_SPACE, os.path.join(tmp, "tune_cancel"), None, cancelling)
        t.equal(result["cancelled"], True, "a cancelled run reports it")
        t.check(result["profile_path"] is None, "a cancelled run returns no profile")
        t.check(not os.path.exists(os.path.join(result["run_dir"], "tuned_profile.properties")), "a cancelled run writes no profile")
        t.check("cancelled" in read_utf8(os.path.join(result["run_dir"], "tuning_log.txt")), "the log records the cancellation")
    finally:
        tune_close_fields(fields)


def test_tune_leaves_standard_runs_unchanged(E, t, tmp):
    """Tuning is a separate tool: it adds nothing to an ordinary run and
    leaves nothing behind that could change one. The same analysis before
    and after a tuning run writes identical files."""
    t.equal(E.ENGINE_MODULE_VERSIONS.get("cellquanttune.py"), E.APP_VERSION, "the module registers with the application version")
    t.equal([extension.name for extension in E.RUN_EXTENSIONS if "tun" in extension.name.lower()], [], "the module adds no per-run measurement")
    configs = [nuclei_config(E), marker_config(E)]
    engine = E.AnalysisEngine()
    t.equal(engine.build_summary_headers(configs, E.default_global_options()), engine.build_summary_headers(configs),
            "the default settings add no summary column")

    def standard_run(where):
        imp = make_image("standard.tif")
        run_dir = os.path.join(tmp, where, "RUN")
        os.makedirs(run_dir)
        options = E.default_global_options()
        try:
            engine.analyze_image(imp, "", configs, [Roi(0, 0, WIDTH, HEIGHT), triangle_roi()],
                                 E.engine_options_from_global(options), run_dir, [1], settings=(options, configs))
        finally:
            imp.close()
        return run_dir

    before = standard_run("tune_before")
    fields = tune_standard_fields()
    try:
        tune_run(E, fields, tune_config(E), {"manual_low": [300, 1000], "watershed": [False, True], "fill_holes": [False, True]},
                 os.path.join(tmp, "tune_between"))
    finally:
        tune_close_fields(fields)
    after = standard_run("tune_after")
    for name in ["roi_summary.csv", "objects.csv", "thresholds.csv"]:
        t.check(read_utf8(os.path.join(before, name)) == read_utf8(os.path.join(after, name)), name + " is byte-for-byte the same before and after tuning")
    t.check(not os.path.exists(os.path.join(after, "tuning_grid.csv")), "an ordinary run writes no tuning table")
    t.equal(sorted(os.listdir(before)), sorted(os.listdir(after)), "an ordinary run writes the same set of files")


# -----------------------------------------------------------------------------
# Known-answer cases added by the independent check
# -----------------------------------------------------------------------------

# Two frames of one 8-bit image (200 x 120, background 10, objects 200),
# calibrated at 0.4 micron per pixel (one pixel is 0.16 micron^2). Disks
# have radius 6 (well over 7 pixels), specks are 2 x 2 = 4 pixels = 0.64
# micron^2, and a "pair" is two overlapping disks 8 pixels apart that the
# threshold joins into one object and the watershed splits in two.
# Parent ROIs: A, an ellipse OvalRoi(0, 0, 100, 120) that touches the left,
# top and bottom image border; B, the rectangle (120, 20, 70, 80).
TUNE_MICRON_FRAMES = [
    {   # frame 1: A holds 3 disks, a pair, a speck; B holds 2 disks, a speck,
        # a mark on background (always missed); one mark lies outside A and B.
        "disks": [(30, 60), (50, 40), (50, 80), (140, 40), (170, 70)],
        "pairs": [((70, 60), (78, 60))],
        "specks": [(60, 95), (150, 85)],
        "extra_marks": [(160, 50), (110, 10)]
    },
    {   # frame 2: A holds 2 disks and a speck; B holds a disk, a pair, a speck.
        "disks": [(30, 60), (50, 40), (170, 80)],
        "pairs": [((140, 50), (148, 50))],
        "specks": [(60, 95), (150, 85)],
        "extra_marks": []
    }
]


def tune_micron_image():
    stack = ImageStack(200, 120)
    for layout in TUNE_MICRON_FRAMES:
        ip = ByteProcessor(200, 120)
        ip.setValue(10.0)
        ip.fill()
        centres = list(layout["disks"])
        for pair in layout["pairs"]:
            centres += [pair[0], pair[1]]
        for centre in centres:
            draw_disk(ip, centre[0], centre[1], DISK_RADIUS, 200)
        for corner in layout["specks"]:
            draw_rect(ip, (corner[0], corner[1], 2, 2), 200)
        stack.addSlice("", ip)
    imp = ImagePlus("micron.tif", stack)
    imp.setDimensions(1, 1, 2)
    imp.setOpenAsHyperStack(True)
    cal = imp.getCalibration()
    cal.pixelWidth = 0.4
    cal.pixelHeight = 0.4
    cal.setUnit("micron")
    return imp


def tune_micron_marks(layout):
    marks = list(layout["disks"])
    for pair in layout["pairs"]:
        marks += [pair[0], pair[1]]
    return marks + list(layout["extra_marks"])


def test_tune_micron_frames_and_rois(E, t, tmp):
    """8-bit, calibrated in micron, two frames of one image as two fields,
    two parent ROIs (one an ellipse on the image border), merged marks and
    a mark outside every ROI. Search: watershed off/on x minimum area 0 or
    1.0 micron^2 (= 6.25 pixels: removes the 4-pixel specks, keeps disks).

    Counts per frame (TP, FP, missed, merged; marks inside ROIs 8 and 5):
                          frame 1        frame 2
        off, 0            6, 2, 1, 1     4, 2, 0, 1
        on,  0            7, 2, 1, 0     5, 2, 0, 0
        off, 1.0          6, 0, 1, 1     4, 0, 0, 1
        on,  1.0          7, 0, 1, 0     5, 0, 0, 0
    Pooled F1 = 2 TP / (TP + FP + marks): 20/27, 24/29, 20/23, 24/25, so
    watershed on with 1.0 micron^2 is tuned (24/25 optimistic). Each fold
    (two-way split) picks the same combination: frame 1 held out scores
    14/15 (one mark on background), frame 2 scores 1; pooled held-out TP
    12, FP 0, missed 1: F1 = 24/25, precision 1, recall 12/13.
    |detected - marks| per parent ROI (A, B), frame 1: off,0 (0, 0); on,0
    (1, 0); off,1 (1, 1); on,1 (0, 1). Frame 2: off,0 (1, 0); on,0 (1, 1);
    off,1 (0, 1); on,1 (0, 0). Totals 1, 3, 3, 1: under the count objective
    the start (off, 0) ties with (on, 1.0) and is kept because it changes
    nothing; held out, frame 1 picks (on, 1.0) from frame 2 (error 1) and
    frame 2 picks (off, 0) from frame 1 (error 1): held-out error 2."""
    imp = tune_micron_image()
    try:
        rois = [OvalRoi(0, 0, 100, 120), Roi(120, 20, 70, 80)]
        fields = []
        for index in range(2):
            fields.append({"name": "frame " + str(index + 1), "source": imp, "frame": index + 1,
                           "rois": rois, "points": tune_micron_marks(TUNE_MICRON_FRAMES[index])})
        config = E.normalize_channel_config(1, {
            "enabled": True, "name": "Nuclei", "analysis_mode": "Objects only", "threshold_mode": "Manual",
            "manual_low": 100, "size_units": "Calibrated area", "min_size": 0
        })
        space = {"watershed": [False, True], "min_size": [0, 1.0]}
        result = tune_run(E, fields, config, space, os.path.join(tmp, "tune_micron"), None, None, [config])
        expected = {
            (False, 0.0): (10, 4, 1, 2, 14, 1), (True, 0.0): (12, 4, 1, 0, 16, 3),
            (False, 1.0): (10, 0, 1, 2, 10, 3), (True, 1.0): (12, 0, 1, 0, 12, 1)
        }
        for row in result["grid_rows"]:
            key = (row["P_watershed"] == "true", row["P_min_size"])
            t.equal((row["True_Positive_Objects"], row["False_Positive_Objects"], row["Missed_Annotations"],
                     row["Merged_Annotations"], row["Detected_Count"], row["Abs_Count_Error"]), expected[key],
                    "micron frames: pooled counts of watershed " + row["P_watershed"] + ", min " + str(key[1]))
            t.equal(row["Manual_Count"], 13, "micron frames: 13 marks inside parent ROIs")
        f1 = dict([((row["P_watershed"], row["P_min_size"]), row["F1"]) for row in result["grid_rows"]])
        t.close(f1[("false", 0.0)], 20.0 / 27.0, "micron frames: F1 of the start")
        t.close(f1[("true", 0.0)], 24.0 / 29.0, "micron frames: F1 with watershed")
        t.close(f1[("false", 1.0)], 20.0 / 23.0, "micron frames: F1 with the area limit")
        t.close(f1[("true", 1.0)], 24.0 / 25.0, "micron frames: F1 with both")
        t.equal(result["tuned_values"], {"watershed": True, "min_size": 1.0}, "micron frames: tuned settings")
        t.equal(result["design"], "two-way split (each field held out once)", "two frames are two fields")
        summary = result["summary_row"]
        t.equal((summary["Marks_Inside_ROIs"], summary["Marks_Outside_ROIs"]), (13, 1), "micron frames: marks inside and outside")
        t.equal((summary["HeldOut_True_Positive_Objects"], summary["HeldOut_False_Positive_Objects"],
                 summary["HeldOut_Missed_Annotations"], summary["HeldOut_Merged_Annotations"]), (12, 0, 1, 0),
                "micron frames: pooled held-out counts")
        t.close(summary["HeldOut_F1"], 24.0 / 25.0, "micron frames: held-out F1")
        t.close(summary["HeldOut_Recall"], 12.0 / 13.0, "micron frames: held-out recall")
        t.close(summary["HeldOut_Precision"], 1.0, "micron frames: held-out precision")
        t.close(summary["HeldOut_Baseline_F1"], 20.0 / 27.0, "micron frames: the start on the same folds")
        t.close(summary["HeldOut_F1_Field_Min"], 14.0 / 15.0, "micron frames: worst held-out field")
        t.equal(summary["HeldOut_Count_Difference"], -1, "micron frames: held-out detected minus marks")
        folds = result["fold_rows"]
        t.equal([(row["True_Positive_Objects"], row["Missed_Annotations"], row["Manual_Count"]) for row in folds],
                [(7, 1, 8), (5, 0, 5)], "micron frames: held-out counts per frame")
        t.equal([(row["Baseline_True_Positive_Objects"], row["Baseline_False_Positive_Objects"], row["Baseline_Merged_Annotations"])
                 for row in folds], [(6, 2, 1), (4, 2, 1)], "micron frames: start counts per frame")
        t.close(folds[0]["Training_F1_Optimistic"], 1.0, "micron frames: tuned on frame 2 alone")
        t.close(folds[1]["Training_F1_Optimistic"], 14.0 / 15.0, "micron frames: tuned on frame 1 alone")
        options, enabled, everything = E.load_profile_settings(result["profile_path"])
        t.equal((enabled[0]["watershed"], enabled[0]["min_size"], enabled[0]["size_units"]), (True, 1.0, "Calibrated area"),
                "micron frames: the tuned profile holds the calibrated area limit")

        counted = tune_run(E, fields, config, space, os.path.join(tmp, "tune_micron"),
                           {"objective": E.TUNE_OBJECTIVE_COUNT}, None, [config])
        t.equal(counted["tuned_values"], {"watershed": False, "min_size": 0.0}, "count objective: a tie keeps the start")
        t.equal(counted["summary_row"]["Resubstitution_Optimistic_Abs_Count_Error"], 1, "count objective: optimistic error")
        t.equal(counted["summary_row"]["HeldOut_Abs_Count_Error"], 2, "count objective: held-out error")
        t.equal([row["P_watershed"] for row in counted["fold_rows"]], ["true", "false"], "count objective: choice per fold")
        t.close(counted["summary_row"]["HeldOut_Abs_Count_Error_Percent"], 200.0 / 13.0, "count objective: held-out error in percent")
    finally:
        imp.close()


def tune_oval_pixels(diameter):
    ip = ByteProcessor(40, 40)
    ip.setValue(255.0)
    ip.fill(OvalRoi(2, 2, diameter, diameter))
    return ip.getHistogram()[255]


def test_tune_size_limits_with_holes(E, t, tmp):
    """A ring (outer oval 19 px wide, hole 9 px wide: R pixels, F = R + hole
    when filled) and a disk of radius 6 (P pixels, P < R). The particle
    analysis filters by the area it measures: the ring alone (R) without
    'include holes', the filled outline (F) with it. Size grid (pixels):
    minimum 0 or R + 1, maximum none or R. Expected detected objects:
                          holes excluded   holes included
        min 0,     no max      2               2
        min R + 1, no max      0               1  (ring only)
        min 0,     max R       2               1  (disk only)
    The size shortcut must give exactly the direct analysis, including the
    zero-object case."""
    ring = tune_oval_pixels(19) - tune_oval_pixels(9)
    engine = E.AnalysisEngine()
    ip = ByteProcessor(100, 60)
    ip.setValue(10.0)
    ip.fill()
    draw_disk(ip, 30, 30, 9, 200)
    draw_disk(ip, 30, 30, 4, 10)
    draw_disk(ip, 75, 30, DISK_RADIUS, 200)
    imp = ImagePlus("ring.tif", ip)
    cal = imp.getCalibration()
    cal.pixelWidth = 0.4
    cal.pixelHeight = 0.4
    cal.setUnit("micron")
    try:
        t.check(disk_pixel_count() < ring, "the disk is smaller than the ring")
        expected = {False: {(0.0, -1.0): 2, (ring + 1.0, -1.0): 0, (0.0, float(ring)): 2},
                    True: {(0.0, -1.0): 2, (ring + 1.0, -1.0): 1, (0.0, float(ring)): 1}}
        for holes in [False, True]:
            config = E.normalize_channel_config(1, {
                "enabled": True, "name": "Ring", "analysis_mode": "Objects only", "threshold_mode": "Manual",
                "manual_low": 100, "size_units": "Pixel^2", "include_holes": holes
            })
            combos = E.tune_build_combinations(config, E.normalize_tune_search_space(
                config, {"min_size": [0, ring + 1], "max_size": [-1, ring]}))
            field = {"name": "ring", "source": imp, "points": [(30, 23), (75, 30)]}
            work = E.tune_open_field(engine, field, config, False, None, False, None)
            try:
                fast = E.tune_score_field(engine, work, combos, True)
                direct = E.tune_score_field(engine, work, combos, False)
            finally:
                E.tune_close_field(work)
            t.equal(fast, direct, "holes " + str(holes) + ": the size shortcut equals the direct analysis")
            for index in range(len(combos) - 1):
                key = (combos[index]["values"]["min_size"], combos[index]["values"]["max_size"])
                if key in expected[holes]:
                    t.equal(fast[index]["detected"], expected[holes][key],
                            "holes " + str(holes) + ", min " + str(key[0]) + ", max " + str(key[1]) + ": objects")
                    t.equal(fast[index]["missed"], 2 - fast[index]["true_positive"], "holes " + str(holes) + ": missed marks")
                else:
                    t.check("error" in fast[index], "a maximum below the minimum is not run")
    finally:
        imp.close()


def tune_counts(manual, true_positive, false_positive):
    """Counts of one field with one parent ROI and no merged marks."""
    detected = true_positive + false_positive
    return {"manual": manual, "detected": detected, "true_positive": true_positive,
            "false_positive": false_positive, "missed": manual - true_positive, "merged": 0,
            "abs_count_error": abs(detected - manual)}


def test_tune_selection_arithmetic(E, t, tmp):
    """Selection and cross-validation on a hand-made table of counts (no
    image). Start manual_low = 1000; grid 900, 1000, 1100 (900 and 1100 both
    change one setting by 0.1). (TP, FP) per field, marks 10, 10, 5, 4:
                900      1000     1100
        F0     8, 1     8, 0     7, 0
        F1     9, 1     error    9, 1
        F2     4, 0     5, 0     4, 0
        F3     3, 0     4, 0     3, 0
    F1 objective 2 TP / (TP + FP + marks). 1000 cannot run on F1, so it is
    ranked only in the fold that leaves F1 out.
      All fields: 900 = 48/55, 1100 = 46/53: 900 is tuned.
      Leave F0 out: 900 = 32/36 = 1100 (tie, same change): grid order, 900.
      Leave F1 out: 900 = 30/35, 1000 = 34/36, 1100 = 28/33: 1000, which
        cannot run on F1: every mark of F1 is missed.
      Leave F2 out: 900 = 40/46, 1100 = 38/44: 900.
      Leave F3 out: 900 = 42/48, 1100 = 40/46: 900.
    Held out: F0 (8, 1), F1 (0, 0), F2 (4, 0), F3 (3, 0): TP 15, FP 1,
    marks 29, F1 = 30/45 = 2/3. The start on the same fields (F1 again all
    missed): TP 17, FP 0: 34/46 = 17/23, better than the tuned settings,
    which the tables must show. Count objective, |detected - marks|:
    900 = 1, 0, 1, 1; 1000 = 2, -, 0, 0; 1100 = 3, 0, 1, 1. All fields:
    900 (3) against 1100 (5). Folds: F0 out tie 2 = 2, grid order 900
    (error 1); F1 out 1000 (2) against 3 and 5 (error 10: all missed);
    F2 out 900 (error 1); F3 out 900 (error 1): held-out error 13."""
    config = tune_config(E, {"manual_low": 1000})
    space = E.normalize_tune_search_space(config, {"manual_low": [900, 1000, 1100]})
    combos = E.tune_build_combinations(config, space)
    t.equal([c["values"]["manual_low"] for c in combos], [900.0, 1000.0, 1100.0, 1000.0], "grid and baseline")
    t.equal(combos[0]["distance"], combos[2]["distance"], "900 and 1100 are equally far from the start")
    error = {"error": "test failure"}
    scores = [
        [tune_counts(10, 8, 1), tune_counts(10, 8, 0), tune_counts(10, 7, 0)],
        [tune_counts(10, 9, 1), error, tune_counts(10, 9, 1)],
        [tune_counts(5, 4, 0), tune_counts(5, 5, 0), tune_counts(5, 4, 0)],
        [tune_counts(4, 3, 0), tune_counts(4, 4, 0), tune_counts(4, 3, 0)]
    ]
    for row in scores:
        row.append(row[1])  # the baseline is the start, 1000
    names = ["F0", "F1", "F2", "F3"]
    marks = [10, 10, 5, 4]
    folds = E.tune_cross_validate(combos, scores, names, marks, E.TUNE_OBJECTIVE_F1)
    t.equal([fold["selected"] for fold in folds], [0, 1, 0, 0], "selected combination per fold (900, 1000, 900, 900)")
    t.check("could not be run" in folds[1]["note"] and "starting settings" in folds[1]["note"],
            "the fold says that neither the selected nor the start could run on F1")
    t.close(E.tune_metrics(folds[1]["train_totals"])["f1"], 34.0 / 36.0, "training score of the fold without F1")
    tables = E.tune_tables("RUN", config, space, combos, names, names, marks, scores, folds, E.TUNE_OBJECTIVE_F1)
    t.equal(tables["final"], 0, "900 is tuned on all fields")
    t.equal([row["Rank"] for row in tables["grid_rows"]], [1, "", 2], "ranks; the combination that fails on F1 has none")
    t.check(tables["grid_rows"][1]["Status"].startswith("FAILED on 1 field(s), first F1"), "the failure is reported")
    t.equal([row["Folds_Selected"] for row in tables["grid_rows"]], [3, 1, 0], "folds won per combination")
    summary = E.tune_summary_row(tables, combos, scores, folds, marks, E.TUNE_OBJECTIVE_F1)
    row = summary["row"]
    t.equal((row["HeldOut_True_Positive_Objects"], row["HeldOut_False_Positive_Objects"], row["HeldOut_Missed_Annotations"],
             row["HeldOut_Manual_Count"]), (15, 1, 14, 29), "pooled held-out counts")
    t.close(row["HeldOut_F1"], 2.0 / 3.0, "pooled held-out F1")
    t.close(row["HeldOut_Precision"], 15.0 / 16.0, "pooled held-out precision")
    t.close(row["HeldOut_Recall"], 15.0 / 29.0, "pooled held-out recall")
    t.close(row["HeldOut_Baseline_F1"], 17.0 / 23.0, "the start on the same folds")
    t.close(row["AllFields_Baseline_F1"], 17.0 / 23.0, "the start on all fields (F1 counted as all missed)")
    t.close(row["Resubstitution_Optimistic_F1"], 48.0 / 55.0, "optimistic score of the tuned settings")
    t.check(row["HeldOut_F1"] < row["HeldOut_Baseline_F1"], "a tuned result worse than the start is reported as such")
    # Per held-out field: F0 2*8/(9+10) = 16/19, F2 2*4/(4+5) = 8/9,
    # F3 2*3/(3+4) = 6/7; F1 (nothing detected) has no F1 and is left out.
    t.close(row["HeldOut_F1_Field_Max"], 8.0 / 9.0, "best held-out field: F2")
    t.close(row["HeldOut_F1_Field_Min"], 16.0 / 19.0, "worst scored held-out field: F0")
    stability = dict([(r["Value"], (r["Folds_Selected"], r["Folds_Scored"], r["Selected_Final"])) for r in tables["stability_rows"]])
    t.equal(stability, {"900": (3, 4, "true"), "1000": (1, 4, "false"), "1100": (0, 4, "false")}, "stability of the threshold")

    counted = E.tune_cross_validate(combos, scores, names, marks, E.TUNE_OBJECTIVE_COUNT)
    t.equal([fold["selected"] for fold in counted], [0, 1, 0, 0], "count objective: selected per fold")
    count_tables = E.tune_tables("RUN", config, space, combos, names, names, marks, scores, counted, E.TUNE_OBJECTIVE_COUNT)
    t.equal(count_tables["final"], 0, "count objective: 900 on all fields")
    count_row = E.tune_summary_row(count_tables, combos, scores, counted, marks, E.TUNE_OBJECTIVE_COUNT)["row"]
    t.equal(count_row["HeldOut_Abs_Count_Error"], 13, "count objective: held-out error")
    t.equal(count_row["Resubstitution_Optimistic_Abs_Count_Error"], 3, "count objective: optimistic error")


def test_tune_profile_compatibility(E, t, tmp):
    """Profiles without tuning keys (2.2 and early 2.3) load with the
    defaults; tuning settings survive a save and reload; and loading this
    module changes nothing in an ordinary run: the same analysis with the
    engine alone and with the engine plus cellquanttune.py writes
    byte-identical tables, and the run record only gains the tune.* keys."""
    from java.util import Properties
    from java.io import FileOutputStream
    defaults = dict([(spec[1], spec[3]) for spec in E.TUNE_SETTING_SPECS])
    t.equal(defaults["tune_max_combinations"], "2000", "default combination limit")
    props = E.settings_to_properties(E.default_global_options(), [nuclei_config(E)])
    old = Properties()
    for key in props.stringPropertyNames():
        if not E.to_text(key).startswith("tune."):
            old.setProperty(key, props.getProperty(key))
    old_path = os.path.join(tmp, "tune_old_profile.properties")
    stream = FileOutputStream(old_path)
    try:
        old.store(stream, "profile without tuning keys")
    finally:
        stream.close()
    options, enabled, everything = E.load_profile_settings(old_path)
    for key in defaults:
        t.equal(options[key], defaults[key], "an old profile gets the default " + key)
    t.equal(enabled[0]["manual_low"], 1000.0, "an old profile keeps its channel settings")

    changed = E.default_global_options()
    changed.update({"tune_channel": "1", "tune_objective": E.TUNE_OBJECTIVE_COUNT,
                    "tune_search_space": "manual_low = 800, 1000; watershed = true, false",
                    "tune_factors": "0.5, 1, 2", "tune_max_combinations": "3", "tune_group_by": E.TUNE_GROUP_FOLDER})
    path = os.path.join(tmp, "tune_round_trip.properties")
    E.write_settings_properties(path, changed, [nuclei_config(E)])
    loaded, enabled, everything = E.load_profile_settings(path)
    for key in defaults:
        t.equal(loaded[key], changed[key], "round trip of " + key)
    config = tune_config(E)
    space = E.normalize_tune_search_space(config, E.parse_tune_search_space(loaded["tune_search_space"]))
    t.equal(space, [("manual_low", [800.0, 1000.0]), ("watershed", [True, False])], "the saved search space reads back")
    refused = ""
    try:
        E.normalize_tune_search_space(config, E.parse_tune_search_space(loaded["tune_search_space"]),
                                      int(loaded["tune_max_combinations"]))
    except ValueError as exc:
        refused = E.to_text(exc)
    t.check("limit is 3" in refused, "the saved combination limit applies")

    # The engine alone against the engine plus this module, same input.
    engine_path = os.path.join(E.CELLQUANT_PACKAGE_DIR, "cellquantengine.py")
    tune_path = os.path.join(E.CELLQUANT_PACKAGE_DIR, "cellquanttune.py")
    outputs = {}
    for label in ["engine_only", "engine_and_tune"]:
        namespace = {"CELLQUANT_NO_GUI": True}
        E.run_source_file(engine_path, namespace)
        if label == "engine_and_tune":
            E.run_source_file(tune_path, namespace)
        N = Namespace(namespace)
        t.equal("run_tuning" in namespace, label == "engine_and_tune", label + ": the module is loaded only when asked")
        configs = [nuclei_config(N), marker_config(N)]
        options = N.default_global_options()
        run_dir = os.path.join(tmp, "tune_compat_" + label, "RUN")
        os.makedirs(run_dir)
        imp = make_image("standard.tif")
        try:
            N.AnalysisEngine().analyze_image(imp, "", configs, [Roi(0, 0, WIDTH, HEIGHT), triangle_roi()],
                                             N.engine_options_from_global(options), run_dir, [1], settings=(options, configs))
        finally:
            imp.close()
        outputs[label] = run_dir
    for name in ["roi_summary.csv", "objects.csv", "thresholds.csv"]:
        t.check(read_utf8(os.path.join(outputs["engine_only"], name)) == read_utf8(os.path.join(outputs["engine_and_tune"], name)),
                name + " is byte-for-byte the same with and without the tuning module")
    t.equal(sorted(os.listdir(outputs["engine_only"])), sorted(os.listdir(outputs["engine_and_tune"])),
            "the same files are written with and without the tuning module")
    plain = E.read_properties_file(os.path.join(outputs["engine_only"], "settings.properties"))
    tuned = E.read_properties_file(os.path.join(outputs["engine_and_tune"], "settings.properties"))
    plain_keys = set([E.to_text(key) for key in plain.stringPropertyNames()])
    tuned_keys = set([E.to_text(key) for key in tuned.stringPropertyNames()])
    t.equal(sorted(tuned_keys - plain_keys), sorted([spec[0] for spec in E.TUNE_SETTING_SPECS]), "the run record only gains the tune.* keys")
    t.equal(sorted(plain_keys - tuned_keys), [], "the run record loses no key")
    t.check(len([key for key in plain_keys if plain.getProperty(key) != tuned.getProperty(key)]) == 0,
            "every shared key of the run record has the same value")


def test_tune_hold_out_by_subfolder(E, t, tmp):
    """Fields of one organoid are held out together when each organoid has
    its own subfolder: org1 holds fields A and C, org2 a second copy of A.
    Subfolders are searched even when the batch setting 'recursive' is off.
    Grid L in {300, 1000} x M in {0, 20}; only L = 1000, M = 20 has F1 = 1
    on A and C (see test_tune_known_window), so both folds select it and
    the pooled held-out F1 is 1. With every field in one subfolder there is
    no held-out estimate, and the summary says why."""
    def save(folder, name, layout):
        if not os.path.isdir(folder):
            os.makedirs(folder)
        imp = tune_image(name + ".tif", tune_plane(layout))
        FileSaver(imp).saveAsTiff(os.path.join(folder, name + ".tif"))
        imp.close()
        RoiEncoder.save(tune_point_roi(layout["nuclei"]), os.path.join(folder, name + "_points.roi"))

    root = os.path.join(tmp, "tune_by_organoid")
    save(os.path.join(root, "org1"), "a", TUNE_FIELD_A)
    save(os.path.join(root, "org1"), "c", TUNE_FIELD_C)
    save(os.path.join(root, "org2"), "a2", TUNE_FIELD_A)
    options = E.default_global_options()
    options.update({"tune_search_space": "manual_low = 300, 1000; min_size = 0, 20", "tune_group_by": E.TUNE_GROUP_FOLDER,
                    "batch_extensions": "tif", "batch_recursive": False})
    config = tune_config(E)
    messages = []
    engine = E.AnalysisEngine(log=messages.append)
    result = E.run_tuning_folder(engine, root, os.path.join(tmp, "tune_by_organoid_out"), [config], [config], options)
    t.equal(result["design"], "leave-one-group-out cross-validation (2 groups)", "fields are held out by organoid")
    t.equal([E.to_text(row["Held_Out_Group"]) for row in result["fold_rows"]], ["org1", "org1", "org2"], "fold groups")
    t.equal([row["Fold"] for row in result["fold_rows"]], [1, 1, 2], "the two fields of org1 form one fold")
    t.equal(result["fold_rows"][0]["Training_Fields"], 1, "org1 held out: tuned on org2's one field")
    t.equal(result["tuned_values"], {"manual_low": 1000.0, "min_size": 20.0}, "tuned settings")
    t.close(result["summary_row"]["HeldOut_F1"], 1.0, "held-out F1 by organoid")
    t.equal(result["summary_row"]["HeldOut_True_Positive_Objects"], 18, "held-out marks found: 6 + 6 + 6")
    t.check(len([m for m in messages if "subfolders are searched" in E.to_text(m)]) == 1, "the log says subfolders were searched")

    single = os.path.join(tmp, "tune_one_organoid")
    save(os.path.join(single, "org1"), "a", TUNE_FIELD_A)
    save(os.path.join(single, "org1"), "c", TUNE_FIELD_C)
    one = E.run_tuning_folder(engine, single, os.path.join(tmp, "tune_by_organoid_out"), [config], [config], options)
    t.equal(one["design"], "none: tuned and scored on the same fields (optimistic)", "one organoid: no held-out design")
    t.check(one["held_out"] is None, "one organoid: no held-out result")
    t.check("belong to one group" in one["summary_row"]["Note"], "one organoid: the summary says why there is no held-out score")


TESTS.append(("tuning from point marks finds the known window and scores it on held-out fields", test_tune_known_window))
TESTS.append(("tuning: a field-specific optimum lowers the held-out score, not the optimistic one", test_tune_held_out_is_honest))
TESTS.append(("tuning: annotated folder convention and missing annotation report", test_tune_annotated_folder))
TESTS.append(("tuning: search space rules, grid order and the combination limit", test_tune_search_space_rules))
TESTS.append(("tuning: size-limit shortcut equals direct analysis; failed combinations and degenerate fields", test_tune_size_shortcut_and_failures))
TESTS.append(("tuning: preparation filters are re-applied exactly; cancel writes no profile", test_tune_preparation_filters_and_cancel))
TESTS.append(("tuning leaves ordinary runs byte-for-byte unchanged", test_tune_leaves_standard_runs_unchanged))
TESTS.append(("tuning: 8-bit micron frames, ellipse ROI on the border, merged marks (known counts)", test_tune_micron_frames_and_rois))
TESTS.append(("tuning: size shortcut with holes and zero objects (known counts)", test_tune_size_limits_with_holes))
TESTS.append(("tuning: selection, ties and cross-validation on hand-computed counts", test_tune_selection_arithmetic))
TESTS.append(("tuning: old profiles, round trip, and module-off runs byte-identical", test_tune_profile_compatibility))
TESTS.append(("tuning: fields of one organoid are held out together (subfolders)", test_tune_hold_out_by_subfolder))


def test_tune_java_error_is_reported(E, t, tmp):
    """An error raised inside ImageJ or a plugin (a Java error) for one
    combination is reported for that combination and the others are still
    scored. Before the fix only Python errors were caught, so the whole
    tuning run stopped. A cancel still stops it."""
    from java.lang import IllegalStateException

    class FailingEngine(E.AnalysisEngine):
        failure = None

        def build_masks(self, prepared_imp, config, title_prefix):
            if config["manual_low"] > 1200:
                raise self.failure
            return E.AnalysisEngine.build_masks(self, prepared_imp, config, title_prefix)

    engine = FailingEngine()
    engine.failure = IllegalStateException("plugin failure")
    config = tune_config(E, {"manual_low": 1000})
    space = E.normalize_tune_search_space(config, {"manual_low": [1000, 1500], "min_size": [20]})
    combos = E.tune_build_combinations(config, space)
    source = tune_image("J.tif", tune_plane(TUNE_FIELD_A))
    field = {"name": "J", "source": source, "points": list(TUNE_FIELD_A["nuclei"]), "rois": None}
    work = E.tune_open_field(engine, field, config, False, None, False, None)
    try:
        scores = E.tune_score_field(engine, work, combos, True)
        by_threshold = {}
        for index in range(len(combos)):
            if combos[index]["status"] == "":
                by_threshold[combos[index]["config"]["manual_low"]] = scores[index]
        t.check("plugin failure" in by_threshold[1500.0].get("error", ""), "the failing combination carries the error text")
        t.check("error" not in by_threshold[1000.0] and by_threshold[1000.0]["detected"] > 0,
                "the working combination is still scored")

        engine.failure = E.AnalysisCancelled("Cancelled by user.")
        cancelled = False
        try:
            E.tune_score_field(engine, work, combos, True)
        except E.AnalysisCancelled:
            cancelled = True
        t.check(cancelled, "a cancel still ends the scoring")
    finally:
        E.tune_close_field(work)
        source.close()


TESTS.append(("tuning: a Java error in one combination is reported, not fatal", test_tune_java_error_is_reported))

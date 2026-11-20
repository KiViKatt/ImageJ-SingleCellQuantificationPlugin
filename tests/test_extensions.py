# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach
# Copyright (C) 2026 Maksym Taran

# Cell Quant Studio Universal - tests of the optional-module mechanism
#
# Loaded by run_tests.py, which provides the helpers used here.
# The source is kept plain ASCII with no encoding declaration on purpose.


def test_run_extension_mechanism(E, t, tmp):
    """A registered extension adds columns after every existing column,
    writes its own table, and changes nothing when it is switched off."""

    class CountingExtension(E.RunExtension):
        name = "test_counting"

        def is_enabled(self, global_options, configs):
            return E.parse_bool(global_options.get("test_counting_enabled"), False)

        def summary_headers(self, configs, global_options):
            return ["Test_Frame_Objects"]

        def object_headers(self, configs, global_options):
            return ["Test_Object_Tag"]

        def begin_run(self, context):
            context["state"][self.name]["rows"] = []

        def frame_done(self, context):
            total = 0
            for key in context["objects"]:
                for obj in context["objects"][key]:
                    obj["row"]["Test_Object_Tag"] = "C" + str(key[0]) + "R" + str(key[1] + 1)
                    total += 1
            for index in range(len(context["summary_rows"])):
                context["summary_rows"][index]["Test_Frame_Objects"] = total
                inside = context["get_roi_masks"](index)[0]
                context["state"][self.name]["rows"].append({
                    "Frame": context["frame"],
                    "Parent_ROI_Index": index + 1,
                    "ROI_Pixels": inside.getHistogram()[255],
                    "Has_Prepared": len(context["prepared"]),
                    "Has_Masks": len(context["masks"])
                })

        def end_run(self, context):
            headers = ["Frame", "Parent_ROI_Index", "ROI_Pixels", "Has_Prepared", "Has_Masks"]
            return [("test_counting.csv", headers, context["state"][self.name]["rows"])]

    saved_extensions = list(E.RUN_EXTENSIONS)
    saved_specs = list(E.GLOBAL_SETTING_SPECS)
    try:
        E.register_run_extension(CountingExtension())
        E.register_run_extension(CountingExtension())
        t.equal(len(E.RUN_EXTENSIONS), len(saved_extensions) + 1, "registering the same name twice keeps one extension")
        E.register_global_settings([("test.counting", "test_counting_enabled", "bool", False)])
        E.register_global_settings([("test.counting", "test_counting_enabled", "bool", True)])
        t.equal(len(E.GLOBAL_SETTING_SPECS), len(saved_specs) + 1, "a setting is registered once")
        t.equal(E.default_global_options()["test_counting_enabled"], False, "the first registered default is kept")

        configs = [nuclei_config(E), marker_config(E)]
        rois = [Roi(0, 0, WIDTH, HEIGHT), triangle_roi()]
        rois[0].setName("whole")
        engine = E.AnalysisEngine()
        base_headers = engine.build_summary_headers(configs)
        base_object_headers = engine.build_object_headers(configs)

        def run(enabled, out_name):
            imp = make_image("ext.tif", frames=2)
            options = E.default_global_options()
            options["test_counting_enabled"] = enabled
            engine_options = E.engine_options_from_global(options)
            for key in ["save_prepared", "save_masks", "save_overlays", "save_object_rois", "save_parent_rois"]:
                engine_options[key] = False
            run_dir = os.path.join(tmp, out_name)
            os.makedirs(run_dir)
            try:
                return engine.analyze_image(
                    imp, "", configs, [E.clone_roi(roi) for roi in rois], engine_options, run_dir, [1, 2],
                    settings=(options, configs)
                )
            finally:
                imp.close()

        off = run(False, "ext_off")
        on = run(True, "ext_on")

        t.equal(off["summary_headers"], base_headers, "switched off: summary columns are unchanged")
        t.equal(off["object_headers"], base_object_headers, "switched off: object columns are unchanged")
        t.equal(off["extra_table_names"], [], "switched off: no extra table")
        t.check(not os.path.exists(os.path.join(off["run_dir"], "test_counting.csv")), "switched off: no extra file")

        t.equal(on["summary_headers"], base_headers + ["Test_Frame_Objects"], "summary columns are appended after all others")
        t.equal(on["object_headers"], base_object_headers + ["Test_Object_Tag"], "object columns are appended after all others")
        t.equal(on["extra_table_names"], ["test_counting.csv"], "the extra table is reported")
        t.check(os.path.isfile(os.path.join(on["run_dir"], "test_counting.csv")), "the extra table is written")
        t.equal(len(on["extra_tables"]["test_counting.csv"]["rows"]), 4, "one extra row per frame and parent ROI")
        t.equal(on["extra_tables"]["test_counting.csv"]["rows"][0]["ROI_Pixels"], WIDTH * HEIGHT, "the ROI mask is available to the extension")
        t.equal(on["extra_tables"]["test_counting.csv"]["rows"][0]["Has_Prepared"], 2, "prepared images are available to the extension")
        t.equal(on["extra_tables"]["test_counting.csv"]["rows"][0]["Has_Masks"], 2, "masks are available to the extension")
        t.equal(on["summary_rows"][0]["Test_Frame_Objects"], len(DISKS) + on["summary_rows"][1]["C1_DAPI_Object_Count"],
                "the extension sees the objects of every parent ROI of the frame")
        tags = set([row.get("Test_Object_Tag") for row in on["object_rows"]])
        t.equal(tags, set(["C1R1", "C1R2"]), "every object row is reachable from its object")
        for key in base_headers:
            if key == "Run_ID":
                continue
            t.check(off["summary_rows"][0].get(key) == on["summary_rows"][0].get(key), "an extension does not change column " + key)
        header_line = read_utf8(os.path.join(on["run_dir"], "roi_summary.csv")).split("\n")[0]
        t.check(header_line.endswith('"Test_Frame_Objects"'), "the appended column is the last one in the file")

        # the same through a batch: combined table batch_<name>
        folder = os.path.join(tmp, "ext_batch")
        os.makedirs(folder)
        for index in range(2):
            imp = make_image("b" + str(index) + ".tif")
            FileSaver(imp).saveAsTiff(os.path.join(folder, "b" + str(index) + ".tif"))
            imp.close()
        options = E.default_global_options()
        options["test_counting_enabled"] = True
        for key in ["save_prepared", "save_masks", "save_overlays", "save_object_rois", "save_parent_rois"]:
            options[key] = False
        result = E.run_batch(folder, os.path.join(tmp, "ext_batch_out"), configs, configs, options,
                             {"extensions": "tif", "roi_policy": E.BATCH_ROI_WHOLE})
        t.equal(result["analyzed"], 2, "the batch analyzed both images")
        combined = os.path.join(result["batch_dir"], "batch_test_counting.csv")
        t.check(os.path.isfile(combined), "the batch writes the combined extra table")
        if os.path.isfile(combined):
            lines = [line for line in read_utf8(combined).split("\n") if line.strip() != ""]
            t.equal(len(lines), 3, "the combined extra table has one row per image and ROI")
            t.check(lines[0].endswith('"Batch_ID","Image_Relative_Path"'), "the combined extra table carries the batch columns")
        summary_header = read_utf8(os.path.join(result["batch_dir"], "batch_roi_summary.csv")).split("\n")[0]
        t.check('"Test_Frame_Objects","Batch_ID","Image_Relative_Path"' in summary_header, "the combined summary carries the appended column")
    finally:
        E.RUN_EXTENSIONS[:] = saved_extensions
        E.GLOBAL_SETTING_SPECS[:] = saved_specs


TESTS.append(("optional analysis modules plug in without changing existing tables", test_run_extension_mechanism))

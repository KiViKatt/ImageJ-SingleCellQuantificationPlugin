#@ File (label="Input folder", style="directory") input_dir
#@ File (label="Settings profile (.properties)", style="file") profile
#@ File (label="Output folder", style="directory") output_dir
#@ String (label="Regions", choices={"Whole image", "ROI file per image (skip image if missing)", "ROI file per image (whole image if missing)", "Automatic tissue outline"}, value="Whole image", required=false, persist=false) roi_policy
#@ String (label="File extensions", value="tif, tiff, czi, lif, nd2, lsm, oib, oir, ics, dv, vsi") extensions
#@ Boolean (label="Include subfolders", value=false) recursive
#@ Boolean (label="Always open with Bio-Formats", value=false) force_bioformats
#@ Boolean (label="Analyze every series of a multi-series file", value=false) all_series

# Cell Quant Studio Universal - batch without the window
#
# Analyzes every image of a folder with a settings profile saved from the
# application ("Save settings profile"), using the same engine as the window.
# The optional analysis modules that lie beside the engine (depth profile,
# threshold spectrum, tissue outline, ...) are loaded too, so a profile that
# switches one on gives its tables here exactly as in the window.
#
# From the Script Editor: open this file, choose Python (Jython), click Run,
# and fill in the dialog.
#
# Without a display (the launcher name depends on the Fiji build):
#   ImageJ-linux64 --ij2 --headless --console --run cellquant_batch_headless.py \
#     'input_dir="/data/images",profile="/data/profile.properties",output_dir="/data/results",roi_policy="Whole image",extensions="tif, czi",recursive=false,force_bioformats=false,all_series=false'
#
# "Regions" (roi_policy) may be left out; it is then "Whole image". The
# policy "Automatic tissue outline" needs cellquanttissue.py beside the
# engine. "Current ROI Manager set for every image" is not offered here
# because a run without the window has no ROI Manager.
#
# Give every other parameter on the command line: Fiji fills a parameter
# that is left out with the value used the last time this script ran on
# this computer, not with the value written above.
#
# Fiji ends with exit status 0 even when the script stops with an error.
# Read the last printed line or batch_report.csv to know how the batch went.
#
# Keep this file beside cellquantengine.py. No window is opened and
# cellquantcore.py is not needed.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij import IJ
import os


def cellquant_path_text(value):
    try:
        return value.getAbsolutePath()
    except:
        return value


def cellquant_find_engine():
    folders = []
    try:
        folders.append(os.path.dirname(os.path.abspath(__file__)))
    except:
        pass
    try:
        imagej_root = IJ.getDirectory("imagej")
        if imagej_root:
            folders.append(os.path.join(imagej_root, "plugins", "Scripts", "Cell Quant Studio"))
    except:
        pass
    for folder in folders:
        if folder and os.path.isfile(os.path.join(folder, "cellquantengine.py")):
            return os.path.join(folder, "cellquantengine.py")
    raise IOError(
        "Cell Quant Studio could not find cellquantengine.py. Keep cellquant_batch_headless.py beside it, "
        "or install the application in Fiji.app/plugins/Scripts/Cell Quant Studio/."
    )


def cellquant_print(message):
    # A console that cannot show a character must not stop the batch.
    try:
        print(to_text(message))
    except:
        print(repr(message))


def cellquant_choose_policy(value):
    """The batch ROI policy asked for, "Whole image" when none is given.
    Policies registered by optional modules are accepted too."""
    policy = ""
    if value is not None:
        policy = to_text(value).strip()
    if policy == "":
        return BATCH_ROI_WHOLE
    usable = [name for name in BATCH_ROI_POLICIES if name != BATCH_ROI_MANAGER]
    if policy not in usable:
        if policy == BATCH_ROI_MANAGER:
            reason = "A batch without the window has no ROI Manager, so '" + policy + "' cannot be used here."
        else:
            reason = "Unknown region choice '" + policy + "'."
        raise ValueError(reason + " Choose one of: " + "; ".join(usable) + ".")
    return policy


def cellquant_check_objects_without_display(configs):
    """Without a display (Fiji started with --headless) some ImageJ classes
    cannot be created. Count objects once on a tiny test image first, so
    that such a problem stops the batch with one clear message instead of
    failing every image."""
    from java.awt import GraphicsEnvironment
    from ij import ImagePlus
    from ij.process import ByteProcessor
    if not GraphicsEnvironment.isHeadless():
        return
    object_configs = [config for config in configs if needs_objects(config)]
    if len(object_configs) == 0:
        return
    probe = ByteProcessor(16, 16)
    probe.setValue(255.0)
    probe.fill(Roi(4, 4, 6, 6))
    probe.resetRoi()
    probe_imp = ImagePlus("CellQuant_probe", probe)
    try:
        AnalysisEngine().analyze_objects(probe_imp, whole_image_roi(probe_imp), object_configs[0])
    except ANY_ERROR as exc:
        raise RuntimeError(
            "Object counting does not work in this Fiji without a display (started with --headless): " +
            to_text(exc) + ". Run the batch from the Script Editor, or start Fiji without --headless "
            "(on a Linux server, for example through xvfb-run). Channels that only measure positive "
            "area or intensity are not affected."
        )
    finally:
        probe_imp.close()


# The parameters are taken before the engine and its modules are loaded into
# this namespace, so that no name defined there can replace one of them.
cellquant_inputs = {
    "input_dir": cellquant_path_text(input_dir), "profile": cellquant_path_text(profile),
    "output_dir": cellquant_path_text(output_dir), "roi_policy": globals().get("roi_policy"),
    "extensions": extensions, "recursive": recursive,
    "force_bioformats": force_bioformats, "all_series": all_series
}

CELLQUANT_NO_GUI = True
cellquant_engine_path = cellquant_find_engine()
execfile(cellquant_engine_path, globals())
cellquant_modules = load_engine_modules(globals(), os.path.dirname(os.path.abspath(cellquant_engine_path)))
if len(cellquant_modules) > 0:
    cellquant_print("Optional modules loaded: " + ", ".join(cellquant_modules))

cellquant_global_options, cellquant_configs, cellquant_all_configs = load_profile_settings(
    to_text(cellquant_inputs["profile"])
)
cellquant_batch_options = {
    "extensions": to_text(cellquant_inputs["extensions"]).strip() or DEFAULT_BATCH_EXTENSIONS,
    "recursive": parse_bool(cellquant_inputs["recursive"]),
    "roi_policy": cellquant_choose_policy(cellquant_inputs["roi_policy"]),
    "roi_dir": "",
    "force_bioformats": parse_bool(cellquant_inputs["force_bioformats"]),
    "all_series": parse_bool(cellquant_inputs["all_series"])
}
cellquant_check_objects_without_display(cellquant_configs)
cellquant_result = run_batch(
    to_text(cellquant_inputs["input_dir"]), to_text(cellquant_inputs["output_dir"]),
    cellquant_configs, cellquant_all_configs, cellquant_global_options, cellquant_batch_options,
    log=cellquant_print
)
cellquant_print(
    "Cell Quant Studio batch: " + str(cellquant_result["analyzed"]) + " analyzed, " +
    str(cellquant_result["skipped"]) + " skipped, " + str(cellquant_result["failed"]) +
    " failed. Results: " + to_text(cellquant_result["batch_dir"])
)
if cellquant_result["failed"] > 0:
    raise RuntimeError(
        "Cell Quant Studio batch finished with " + str(cellquant_result["failed"]) +
        " failed file(s). See batch_report.csv in " + to_text(cellquant_result["batch_dir"])
    )

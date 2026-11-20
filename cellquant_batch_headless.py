#@ File (label="Input folder", style="directory") input_dir
#@ File (label="Settings profile (.properties)", style="file") profile
#@ File (label="Output folder", style="directory") output_dir
#@ String (label="Regions", choices={"Whole image", "ROI file per image (skip image if missing)", "ROI file per image (whole image if missing)"}, value="Whole image") roi_policy
#@ String (label="File extensions", value="tif, tiff, czi, lif, nd2, lsm, oib, oir, ics, dv, vsi") extensions
#@ Boolean (label="Include subfolders", value=false) recursive
#@ Boolean (label="Always open with Bio-Formats", value=false) force_bioformats
#@ Boolean (label="Analyze every series of a multi-series file", value=false) all_series

# Cell Quant Studio Universal - batch without the window
#
# Analyzes every image of a folder with a settings profile saved from the
# application ("Save settings profile"), using the same engine as the window.
#
# From the Script Editor: open this file, choose Python (Jython), click Run,
# and fill in the dialog.
#
# Without a display (the launcher name depends on the Fiji build):
#   ImageJ-linux64 --ij2 --headless --console --run cellquant_batch_headless.py \
#     'input_dir="/data/images",profile="/data/profile.properties",output_dir="/data/results",roi_policy="Whole image",extensions="tif, czi",recursive=false,force_bioformats=false,all_series=false'
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


CELLQUANT_NO_GUI = True
execfile(cellquant_find_engine(), globals())

cellquant_global_options, cellquant_configs, cellquant_all_configs = load_profile_settings(
    to_text(cellquant_path_text(profile))
)
cellquant_batch_options = {
    "extensions": to_text(extensions),
    "recursive": parse_bool(recursive),
    "roi_policy": to_text(roi_policy),
    "roi_dir": "",
    "force_bioformats": parse_bool(force_bioformats),
    "all_series": parse_bool(all_series)
}
cellquant_result = run_batch(
    to_text(cellquant_path_text(input_dir)), to_text(cellquant_path_text(output_dir)),
    cellquant_configs, cellquant_all_configs, cellquant_global_options, cellquant_batch_options,
    log=cellquant_print
)
print(
    "Cell Quant Studio batch: " + str(cellquant_result["analyzed"]) + " analyzed, " +
    str(cellquant_result["skipped"]) + " skipped, " + str(cellquant_result["failed"]) +
    " failed. Results: " + to_text(cellquant_result["batch_dir"])
)
if cellquant_result["failed"] > 0:
    raise RuntimeError(
        "Cell Quant Studio batch finished with " + str(cellquant_result["failed"]) +
        " failed file(s). See batch_report.csv in " + to_text(cellquant_result["batch_dir"])
    )

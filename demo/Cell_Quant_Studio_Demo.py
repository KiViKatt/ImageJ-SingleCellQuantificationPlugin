# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach
# Copyright (C) 2026 Maksym Taran

# Cell Quant Studio Universal - one-click demo
#
# Opens the application, loads a demo image and a demo profile, outlines the
# tissue, adds three depth zones that do not overlap (the outer 0-100 um, 100-
# 200 um, and the core deeper than 200 um, measured from the tissue edge
# only, not from the edge of the picture), runs the analysis, and ends on the
# Results tab. Each step pauses briefly so an audience can follow it.
#
# Install: copy this file beside Cell_Quant_Studio_final.py (for example in
# Fiji.app/plugins/Scripts/Cell Quant Studio/) and run it from the Plugins
# menu or the Script Editor. The first run asks for the demo image and the
# demo profile and remembers them (Fiji preferences). To choose them again,
# hold Shift while starting the demo.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

import os
import threading
import time

from ij import IJ, ImagePlus, Prefs
from ij.gui import Roi
from ij.io import OpenDialog
from ij.plugin.filter import EDM, ThresholdToSelection
from ij.plugin.frame import RoiManager
from ij.process import ByteProcessor, ImageProcessor
from java.awt import Window
from java.lang import Runnable
from java.lang import Throwable as JavaThrowable
from javax.swing import JDialog, SwingUtilities

DEMO_ERRORS = (Exception, JavaThrowable)
DEMO_PAUSE = 1.5
DEMO_IMAGE_KEY = "cellquant.demo.image"
DEMO_PROFILE_KEY = "cellquant.demo.profile"
DEMO_OUTPUT = "C:/CQ_demo" if os.name == "nt" or os.sep == "\\" else os.path.join(os.path.expanduser("~"), "CQ_demo")
MICRO = b"\xc2\xb5".decode("utf-8")


class DemoTask(Runnable):
    def __init__(self, function):
        self.function = function
        self.done = threading.Event()
        self.error = None
        self.value = None

    def run(self):
        try:
            self.value = self.function()
        except DEMO_ERRORS as exc:
            self.error = exc
        self.done.set()


def demo_ui(function):
    """Run function on the event thread and wait for it."""
    task = DemoTask(function)
    SwingUtilities.invokeLater(task)
    task.done.wait(120)
    if task.error is not None:
        raise task.error
    return task.value


def demo_package_dir():
    folders = []
    try:
        folders.append(os.path.dirname(os.path.abspath(__file__)))
    except DEMO_ERRORS:
        pass
    root = IJ.getDirectory("imagej")
    if root:
        folders.append(os.path.join(root, "plugins", "Scripts", "Cell Quant Studio"))
    for folder in folders:
        if os.path.isfile(os.path.join(folder, "Cell_Quant_Studio_final.py")):
            return folder
    raise IOError("Keep Cell_Quant_Studio_Demo.py beside Cell_Quant_Studio_final.py.")


def demo_window():
    for window in Window.getWindows():
        try:
            if window.isDisplayable() and window.getTitle().startswith("Cell Quant Studio") and hasattr(window, "run_analysis"):
                return window
        except DEMO_ERRORS:
            pass
    return None


def demo_settings_file():
    """Paths from Cell_Quant_Studio_Demo.txt beside this script, if present:
    lines such as image=C:/data/demo.czi and profile=C:/data/demo.properties."""
    values = {}
    try:
        path = os.path.join(demo_package_dir(), "Cell_Quant_Studio_Demo.txt")
        if os.path.isfile(path):
            handle = open(path, "r")
            try:
                for line in handle.read().splitlines():
                    if "=" in line and not line.strip().startswith("#"):
                        key, value = line.split("=", 1)
                        values[key.strip()] = value.strip()
            finally:
                handle.close()
    except DEMO_ERRORS:
        pass
    return values


def demo_choose(key, title):
    preset = demo_settings_file().get("image" if key == DEMO_IMAGE_KEY else "profile", "")
    if preset and os.path.isfile(preset) and not IJ.shiftKeyDown():
        return preset
    path = Prefs.get(key, "")
    if path and os.path.isfile(path) and not IJ.shiftKeyDown():
        return path
    dialog = OpenDialog(title)
    if dialog.getFileName() is None:
        raise RuntimeError("Demo cancelled.")
    path = os.path.join(dialog.getDirectory(), dialog.getFileName())
    Prefs.set(key, path)
    Prefs.savePreferences()
    return path


def demo_zones(tissue, width, height, pixel):
    """Three zones that do not overlap: 0-100 um from the tissue edge,
    100-200 um, and deeper than 200 um. Distances are measured to the
    background inside the picture only; the edge of the picture is not a
    tissue edge."""
    mask = ByteProcessor(width, height)
    mask.setValue(255)
    mask.fill(tissue)
    distance = EDM().makeFloatEDM(mask, 0, False)
    first = 100.0 / pixel
    second = 200.0 / pixel
    zones = []
    for name, low, high in [
        ("0-100 " + MICRO + "m", 0.5, first),
        ("100-200 " + MICRO + "m", first + 1e-3, second),
        ("Over 200 " + MICRO + "m", second + 1e-3, 1e9),
    ]:
        band = distance.duplicate()
        band.setThreshold(low, high, ImageProcessor.NO_LUT_UPDATE)
        roi = ThresholdToSelection.run(ImagePlus("zone", band))
        if roi is not None and roi.getBounds().width > 3:
            roi.setName(name)
            zones.append(roi)
    return zones


def demo_main():
    package = demo_package_dir()
    image = demo_choose(DEMO_IMAGE_KEY, "Demo image")
    profile = demo_choose(DEMO_PROFILE_KEY, "Demo profile (.properties)")

    frame = demo_window()
    if frame is None:
        launcher = os.path.join(package, "Cell_Quant_Studio_final.py")
        namespace = {"__file__": launcher, "__name__": "__main__"}
        execfile(launcher, namespace)
        for attempt in range(200):
            frame = demo_window()
            if frame is not None and frame.isShowing():
                break
            time.sleep(0.1)
    if frame is None:
        raise RuntimeError("The Cell Quant Studio window did not open.")

    def say(text):
        demo_ui(lambda: frame.append_status("Demo: " + text))

    look = demo_settings_file().get("look", "Dark")
    demo_set_look(frame, look)
    # The demo shows the organoid map (the Pitch view of the Results tab)
    # from the moment the analysis ends, without changing the view the
    # user chose for their own work.
    demo_pitch_view(frame)
    demo_ui(lambda: frame.main_tabs.setSelectedIndex(0))
    say("opening the image")
    demo_ui(lambda: frame.open_image_path(image, True))
    imp = frame.source_imp
    if imp is None:
        raise RuntimeError("The demo image did not open.")
    time.sleep(DEMO_PAUSE)

    say("loading the settings profile")
    demo_ui(lambda: frame.load_profile_from_path(profile))

    def outputs():
        if not os.path.isdir(DEMO_OUTPUT):
            os.makedirs(DEMO_OUTPUT)
        frame.output_dir_field.setText(DEMO_OUTPUT)
        frame.save_prepared_cb.setSelected(True)
        for box in [frame.show_summary_cb, frame.keep_prepared_open_cb, frame.keep_masks_open_cb]:
            box.setSelected(False)
    demo_ui(outputs)
    demo_ui(lambda: frame.main_tabs.setSelectedIndex(0))
    time.sleep(DEMO_PAUSE)

    say("finding the tissue outline")
    globals_of_app = demo_app_globals(frame)
    engine = globals_of_app["AnalysisEngine"]()
    rois, info = globals_of_app["detect_tissue_outline"](
        engine, imp, frame.collect_channel_configs(), frame.collect_global_options(), frame=1
    )
    if len(rois) == 0:
        raise RuntimeError("No tissue was found in the demo image.")
    pixel = imp.getCalibration().pixelWidth or 1.0
    zones = demo_zones(rois[0], imp.getWidth(), imp.getHeight(), pixel)
    manager = demo_ui(lambda: RoiManager.getRoiManager())
    demo_ui(lambda: manager.reset())
    for zone in zones:
        demo_ui(lambda: manager.addRoi(zone))
    demo_ui(lambda: frame.roi_source_combo.setSelectedItem("ROI Manager"))
    say("3 depth zones from the tissue edge")
    time.sleep(DEMO_PAUSE)

    say("analysing")
    demo_ui(frame.run_analysis)
    end = time.time() + 600
    time.sleep(1.0)
    while time.time() < end and (frame.background_active or demo_close_finished_dialog()):
        time.sleep(0.5)
        demo_close_finished_dialog()
    time.sleep(0.5)

    def to_results():
        for index in range(frame.main_tabs.getTabCount()):
            if frame.main_tabs.getTitleAt(index) == "Results":
                frame.main_tabs.setSelectedIndex(index)
    demo_set_look(frame, look)
    demo_pitch_view(frame)
    demo_ui(to_results)
    say("done")
    try:
        handle = open(os.path.join(DEMO_OUTPUT, "demo_done.txt"), "w")
        handle.write(time.strftime("%Y-%m-%d %H:%M:%S") + "\n")
        handle.close()
    except DEMO_ERRORS:
        pass


def demo_close_finished_dialog():
    """Close the "Analysis complete" message, so the demo goes straight to
    the Results tab. Returns True when one was closed. Uses invokeLater
    without waiting, because the message is modal."""
    closed = []
    for window in Window.getWindows():
        try:
            if isinstance(window, JDialog) and window.isShowing() and \
                    window.getTitle().startswith("Cell Quant Studio") and \
                    demo_dialog_text(window).find("Analysis complete") >= 0:
                SwingUtilities.invokeLater(DemoTask(window.dispose))
                closed.append(window)
        except DEMO_ERRORS:
            pass
    return len(closed) > 0


def demo_dialog_text(container):
    words = []
    for child in container.getComponents():
        try:
            if hasattr(child, "getText") and child.getText() is not None:
                words.append(child.getText())
            elif hasattr(child, "getMessage") and child.getMessage() is not None:
                words.append(str(child.getMessage()))
        except DEMO_ERRORS:
            pass
        if hasattr(child, "getComponents"):
            words.append(demo_dialog_text(child))
    return " ".join(words)


def demo_pitch_view(frame):
    """Show the Results tab as the organoid map for this demo. The view the
    user chose (Researcher or Pitch) stays the remembered one."""
    if not hasattr(frame, "set_results_view"):
        return
    try:
        demo_ui(lambda: frame.set_results_view("Pitch", False))
    except DEMO_ERRORS as exc:
        IJ.log("Cell Quant Studio demo: could not show the Pitch view (" + str(exc) + ")")


def demo_set_look(frame, look):
    """Show the window in the demo's look, and check that it took."""
    if not hasattr(frame, "set_look"):
        return
    for attempt in range(3):
        try:
            demo_ui(lambda: frame.set_look(look))
        except DEMO_ERRORS as exc:
            IJ.log("Cell Quant Studio demo: could not set the look (" + str(exc) + ")")
        if Prefs.get("cellquant.look", "") == look:
            return
        time.sleep(0.5)


def demo_app_globals(frame):
    """The namespace the application was loaded into (engine and modules)."""
    method = frame.__class__.run_analysis
    function = getattr(method, "im_func", method)
    names = getattr(function, "func_globals", None) or getattr(function, "__globals__", None)
    if names is None or "detect_tissue_outline" not in names:
        raise RuntimeError("The tissue-outline module is not loaded; keep cellquanttissue.py beside the launcher.")
    return names


try:
    demo_main()
except DEMO_ERRORS as demo_error:
    IJ.log("Cell Quant Studio demo: " + str(demo_error))
    IJ.showMessage("Cell Quant Studio demo", "The demo stopped: " + str(demo_error))

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach

# Cell Quant Studio Universal - pictures of every tab of the window
#
# Opens the application the way a user does, shows each tab in turn, and saves
# a PNG of the window for each tab, in each look. Used to review the interface
# design without clicking through it by hand.
#
#   bash tests/run_tests.sh <Fiji folder> "<absolute path>/tests/ui_snapshots.py"
#
# Environment variables:
#   CELLQUANT_SNAPSHOT_DIR     folder for the pictures (default: a "ui_snapshots"
#                              folder in the system temp folder)
#   CELLQUANT_SNAPSHOT_LOOKS   the looks, comma-separated ("Light,Dark" by
#                              default; the names the window offers)
#   CELLQUANT_SNAPSHOT_WIDTH, CELLQUANT_SNAPSHOT_HEIGHT   window size (1440 x 900)
#   CELLQUANT_SNAPSHOT_DEMO    "1": first analyse a synthetic organoid section,
#                              so that the Session and Results tabs show content
#   CELLQUANT_SNAPSHOT_IMAGE   instead, an image file to open (with Bio-Formats)
#   CELLQUANT_SNAPSHOT_PROFILE and a settings profile to load for it: the
#                              tissue outline and two zones inside it (100 and
#                              200 micrometers in) become the parent ROIs, and
#                              the image is analysed with the profile. Output
#                              goes to a temporary folder that is removed.
#
# The window and everything the demo opened are closed afterwards; the look
# remembered in Fiji's preferences is put back.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

import math
import os
import random
import shutil
import tempfile
import threading
import time

from ij import IJ, ImagePlus, ImageStack, CompositeImage, WindowManager, Prefs
from ij.gui import OvalRoi
from ij.plugin.frame import RoiManager
from ij.process import ShortProcessor, LUT
from java.awt import Window, Color
from java.awt.image import BufferedImage
from java.io import File
from java.lang import Runnable
from java.lang import Throwable as JavaThrowable
from javax.imageio import ImageIO
from javax.swing import SwingUtilities, JDialog

ANY_ERROR_SNAP = (Exception, JavaThrowable)


class SnapTask(Runnable):
    def __init__(self, function):
        self.function = function
        self.done = threading.Event()
        self.error = None
        self.value = None

    def run(self):
        try:
            self.value = self.function()
        except ANY_ERROR_SNAP as exc:
            self.error = exc
        self.done.set()


def snap_on_ui(function):
    task = SnapTask(function)
    SwingUtilities.invokeLater(task)
    task.done.wait(60)
    if task.error is not None:
        raise task.error
    return task.value


def snap_find_launcher():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(os.path.dirname(here), "Cell_Quant_Studio_final.py")


def snap_find_frame(title_start):
    for window in Window.getWindows():
        try:
            if window.isDisplayable() and window.getTitle().startswith(title_start):
                return window
        except ANY_ERROR_SNAP:
            pass
    return None


def snap_settle(frame):
    """Lay the window out twice (wrapped texts take their height in the
    second pass) and let the event queue empty."""
    for attempt in range(3):
        snap_on_ui(frame.validate)
        time.sleep(0.15)


def snap_paint(frame, path):
    def paint():
        frame.validate()
        root = frame.getRootPane()
        # Painted at the scale of the screen (for example 150 percent), so
        # that text is laid out exactly as the window shows it.
        transform = frame.getGraphicsConfiguration().getDefaultTransform()
        scale_x = transform.getScaleX()
        scale_y = transform.getScaleY()
        image = BufferedImage(int(math.ceil(root.getWidth() * scale_x)), int(math.ceil(root.getHeight() * scale_y)),
                              BufferedImage.TYPE_INT_RGB)
        graphics = image.createGraphics()
        graphics.scale(scale_x, scale_y)
        root.paint(graphics)
        graphics.dispose()
        ImageIO.write(image, "png", File(path))
        return (root.getWidth(), root.getHeight())
    return snap_on_ui(paint)


def snap_size(frame, width, height):
    if frame.getExtendedState() != 0:
        frame.setExtendedState(0)
    if frame.getWidth() != width or frame.getHeight() != height:
        frame.setSize(width, height)
    frame.validate()


def snap_slug(text):
    keep = []
    for char in text.lower():
        keep.append(char if char.isalnum() else "_")
    return "".join(keep).strip("_")


def snap_close_dialogs():
    def close():
        for window in Window.getWindows():
            if isinstance(window, JDialog) and window.isShowing():
                window.dispose()
    snap_on_ui(close)


def snap_demo_image():
    """A synthetic organoid section: nuclei (C1) in a round tissue, and a
    marker (C2) in patches. 0.65 micrometer per pixel."""
    chance = random.Random(11)
    size = 640
    middle = size / 2.0
    nuclei = ShortProcessor(size, size)
    nuclei.setValue(180.0)
    nuclei.fill()
    marker = ShortProcessor(size, size)
    marker.setValue(160.0)
    marker.fill()
    nuclei.setValue(2600.0)
    placed = 0
    while placed < 320:
        x = chance.uniform(40, size - 40)
        y = chance.uniform(40, size - 40)
        if (x - middle) ** 2 + (y - middle) ** 2 > 265 ** 2:
            continue
        radius = chance.uniform(4.5, 7.0)
        nuclei.fill(OvalRoi(x - radius, y - radius, 2 * radius, 2 * radius))
        placed += 1
    marker.setValue(1500.0)
    for patch in range(40):
        x = chance.uniform(90, size - 90)
        y = chance.uniform(90, size - 90)
        if (x - middle) ** 2 + (y - middle) ** 2 > 230 ** 2:
            continue
        radius = chance.uniform(10, 30)
        marker.fill(OvalRoi(x - radius, y - radius, 2 * radius, 2 * radius))
    nuclei.resetRoi()
    marker.resetRoi()
    nuclei.blurGaussian(1.2)
    marker.blurGaussian(4.0)
    stack = ImageStack(size, size)
    stack.addSlice("", nuclei)
    stack.addSlice("", marker)
    imp = ImagePlus("Organoid section (synthetic).tif", stack)
    imp.setDimensions(2, 1, 1)
    imp = CompositeImage(imp, CompositeImage.COMPOSITE)
    imp.setChannelLut(LUT.createLutFromColor(Color(70, 130, 255)), 1)
    imp.setChannelLut(LUT.createLutFromColor(Color(60, 230, 110)), 2)
    cal = imp.getCalibration()
    cal.pixelWidth = 0.65
    cal.pixelHeight = 0.65
    cal.setUnit("um")
    return imp


def snap_demo(frame, output):
    """Analyse the synthetic section in three parent ROIs, as a user would."""
    imp = snap_demo_image()
    snap_on_ui(imp.show)
    snap_on_ui(lambda: frame.use_image_as_source(imp))

    def setup():
        frame.channel_panels[0].set_values({
            "enabled": True, "name": "DAPI", "analysis_mode": "Objects + area", "lut": "Blue",
            "threshold_mode": "Manual", "manual_low": 1000, "manual_high": -1, "projection": "Max Z",
            "clahe": False, "watershed": True, "fill_holes": True, "subtract_background": False
        })
        frame.channel_panels[1].set_values({
            "enabled": True, "name": "MAP2", "analysis_mode": "Positive area", "lut": "Green",
            "threshold_mode": "Manual", "manual_low": 700, "manual_high": -1, "projection": "Max Z",
            "clahe": False, "subtract_background": False
        })
        for panel in frame.channel_panels[2:]:
            panel.set_values({"enabled": False})
        for box in [frame.show_summary_cb, frame.keep_prepared_open_cb, frame.keep_masks_open_cb]:
            box.setSelected(False)
        frame.output_dir_field.setText(output)
        manager = RoiManager.getRoiManager()
        manager.reset()
        for name, inset in [("Outer zone", 50), ("Middle zone", 130), ("Core", 210)]:
            roi = OvalRoi(inset, inset, 640 - 2 * inset, 640 - 2 * inset)
            roi.setName(name)
            manager.addRoi(roi)
        frame.roi_source_combo.setSelectedItem("ROI Manager")
    snap_on_ui(setup)
    snap_on_ui(frame.run_analysis)
    end = time.time() + 120
    while time.time() < end:
        time.sleep(0.3)
        snap_close_dialogs()
        if not frame.background_active:
            break
    time.sleep(0.5)
    snap_close_dialogs()
    snap_on_ui(lambda: frame.main_tabs.setSelectedIndex(0))
    return imp


def snap_real(frame, namespace, image, profile, output):
    """Open a real image, apply a profile, outline the tissue and two zones
    inside it as parent ROIs, and analyse it, as a user would."""
    from ij.plugin import RoiEnlarger
    snap_on_ui(lambda: frame.open_image_path(image, True))
    imp = frame.source_imp
    if imp is None:
        raise RuntimeError("The image did not open: " + image)
    if profile:
        snap_on_ui(lambda: frame.load_profile_from_path(profile))

    def outputs():
        # The picture of the Results tab is made from the prepared images;
        # nothing else needs to stay open.
        frame.save_prepared_cb.setSelected(True)
        for box in [frame.show_summary_cb, frame.keep_prepared_open_cb, frame.keep_masks_open_cb]:
            box.setSelected(False)
        frame.output_dir_field.setText(output)
    snap_on_ui(outputs)
    engine = namespace["AnalysisEngine"]()
    rois, info = namespace["detect_tissue_outline"](
        engine, imp, frame.collect_channel_configs(), frame.collect_global_options(), frame=1
    )
    if len(rois) == 0:
        raise RuntimeError("No tissue outline was found.")
    tissue = rois[0]
    pixel = imp.getCalibration().pixelWidth or 1.0
    manager = snap_on_ui(lambda: RoiManager.getRoiManager())
    snap_on_ui(lambda: manager.reset())
    for name, inset in [("Section", 0.0), ("100 um in", 100.0), ("200 um in", 200.0)]:
        zone = tissue if inset == 0 else RoiEnlarger.enlarge(tissue, -int(round(inset / pixel)))
        if zone is None or zone.getBounds().width < 4:
            continue
        zone.setName(name)
        snap_on_ui(lambda: manager.addRoi(zone))
    snap_on_ui(lambda: frame.roi_source_combo.setSelectedItem("ROI Manager"))
    started = time.time()
    snap_on_ui(frame.run_analysis)
    end = time.time() + 600
    while time.time() < end:
        time.sleep(0.5)
        snap_close_dialogs()
        if not frame.background_active:
            break
    print("analysis took %.0f s" % (time.time() - started))
    time.sleep(1.0)
    snap_close_dialogs()
    snap_on_ui(lambda: frame.main_tabs.setSelectedIndex(0))
    return imp


def snap_extras(frame, demo):
    """More states worth a picture: (name, set up, put back). The run log
    open on the Session tab; with the demo, a ring pointed at and the
    picture of the run enlarged on the Results tab."""
    extras = []
    if getattr(frame, "toggle_run_log", None) is not None:
        def log_open():
            frame.main_tabs.setSelectedIndex(0)
            frame.toggle_run_log()
        extras.append(("09_session_run_log", log_open, lambda: frame.toggle_run_log()))
    canvas = getattr(frame, "results_canvas", None)
    if demo and canvas is not None and len(canvas.pictures) > 0:
        def point():
            frame.main_tabs.setSelectedIndex(frame.main_tabs.indexOfTab("Results"))
            frame.validate()
            shape = canvas.geometry()
            if len(shape["rings"]) > 1:
                canvas.point_at(shape["x"] + shape["rings"][1][0], shape["y"])

        def unpoint():
            canvas.point_at(-1000, -1000)

        def zoom():
            frame.main_tabs.setSelectedIndex(frame.main_tabs.indexOfTab("Results"))
            canvas.zoom = True
            canvas.repaint()

        def unzoom():
            canvas.zoom = False
            canvas.repaint()
        extras.append(("10_results_hover", point, unpoint))
        extras.append(("11_results_picture", zoom, unzoom))
    panels = getattr(frame, "channel_panels", [])
    if demo and getattr(frame, "refresh_channel_list", None) is not None and len(panels) > 1:
        # The check of channel names against the file: C2 called DAPI.
        kept = {}

        def misname():
            frame.main_tabs.setSelectedIndex(0)
            kept["name"] = panels[1].name_field.getText()
            panels[1].name_field.setText("DAPI")
            frame.refresh_channel_list()

        def rename():
            panels[1].name_field.setText(kept["name"])
            frame.refresh_channel_list()
        extras.append(("12_session_name_check", misname, rename))
    return extras


def snap_cleanup():
    def close():
        ids = WindowManager.getIDList()
        if ids is not None:
            for image_id in list(ids):
                image = WindowManager.getImage(image_id)
                if image is not None:
                    image.changes = False
                    image.close()
        manager = RoiManager.getInstance()
        if manager is not None:
            manager.reset()
            manager.close()
    snap_on_ui(close)


def snap_main():
    out = os.environ.get("CELLQUANT_SNAPSHOT_DIR") or os.path.join(tempfile.gettempdir(), "ui_snapshots")
    if not os.path.isdir(out):
        os.makedirs(out)
    looks = [name.strip() for name in (os.environ.get("CELLQUANT_SNAPSHOT_LOOKS") or "Light,Dark").split(",") if name.strip()]
    demo = (os.environ.get("CELLQUANT_SNAPSHOT_DEMO") or "") == "1"
    image = os.environ.get("CELLQUANT_SNAPSHOT_IMAGE") or ""
    profile = os.environ.get("CELLQUANT_SNAPSHOT_PROFILE") or ""
    remembered = Prefs.get("cellquant.look", "")
    launcher = snap_find_launcher()
    namespace = {"__file__": launcher, "__name__": "__main__"}
    execfile(launcher, namespace)
    frame = None
    for attempt in range(100):
        frame = snap_find_frame("Cell Quant Studio")
        if frame is not None and frame.isShowing():
            break
        time.sleep(0.1)
    if frame is None:
        raise RuntimeError("The application window did not open.")
    width = int(os.environ.get("CELLQUANT_SNAPSHOT_WIDTH") or "1440")
    height = int(os.environ.get("CELLQUANT_SNAPSHOT_HEIGHT") or "900")
    snap_on_ui(lambda: frame.setSize(width, height))
    written = []
    work = tempfile.mkdtemp(prefix="cellquant_snap_")
    try:
        if image:
            snap_real(frame, namespace, image, profile, work)
            demo = True
        elif demo:
            snap_demo(frame, work)
        for look in looks:
            chooser = getattr(frame, "set_look", None)
            combo = getattr(frame, "look_combo", None)
            try:
                if chooser is not None:
                    snap_on_ui(lambda: chooser(look))
                elif combo is not None:
                    snap_on_ui(lambda: combo.setSelectedItem(look))
            except ANY_ERROR_SNAP as exc:
                print("Look " + look + " not available: " + str(exc))
                continue
            time.sleep(0.5)
            tabs = frame.main_tabs
            count = snap_on_ui(lambda: tabs.getTabCount())
            for index in range(count):
                # The asked size, also when something maximized the window meanwhile.
                snap_on_ui(lambda: snap_size(frame, width, height))
                title = snap_on_ui(lambda: tabs.getTitleAt(index))
                snap_on_ui(lambda: tabs.setSelectedIndex(index))
                time.sleep(0.3)
                snap_settle(frame)
                path = os.path.join(out, "%s_%02d_%s.png" % (snap_slug(look), index + 1, snap_slug(title)))
                size = snap_paint(frame, path)
                written.append(path)
                print("snapshot %s (%dx%d)" % (path, size[0], size[1]))
            for name, before, after in snap_extras(frame, demo):
                snap_on_ui(lambda: snap_size(frame, width, height))
                snap_on_ui(before)
                snap_settle(frame)
                path = os.path.join(out, "%s_%s.png" % (snap_slug(look), name))
                snap_paint(frame, path)
                snap_on_ui(after)
                written.append(path)
                print("snapshot %s" % path)
    finally:
        snap_on_ui(lambda: frame.dispose())
        if demo:
            snap_cleanup()
        Prefs.set("cellquant.look", remembered)
        shutil.rmtree(work, True)
    print("%d snapshots in %s" % (len(written), out))


snap_main()

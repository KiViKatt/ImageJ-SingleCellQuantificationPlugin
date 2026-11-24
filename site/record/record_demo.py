# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach
# Copyright (C) 2026 Maksym Taran

# Cell Quant Studio Universal - pictures of the one-click demo for the landing page
#
# Takes the steps of demo/Cell_Quant_Studio_Demo.py in the same order (open
# the image, load the profile, find the tissue outline, three depth zones,
# analyse, Results) and saves a picture of the window at each one, in the
# Light and in the Dark look. See site/README.md.
#
#   bash tests/run_tests.sh <Fiji folder> "<absolute path>/site/record/record_demo.py"
#
# Environment variables:
#   CQ_DEMO_IMAGE    the image to open (required; make_demo_image.py builds one)
#   CQ_RESULTS       the output folder (required; its path is shown in the window)
#   CQ_SHOTS         where the pictures go (required)
#   CQ_DEMO_PROFILE  the settings profile (default: the demo profile in examples/)
#   CQ_REPO          the folder with Cell_Quant_Studio_final.py (default: two
#                    folders above this file)
#
# The window is on the screen while this runs. Do not click in it: a real
# click changes what the next picture shows.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

import os
import threading
import time

from ij import IJ, ImagePlus, Prefs, WindowManager
from ij.gui import Overlay, TextRoi
from ij.plugin.filter import EDM, ThresholdToSelection
from ij.plugin.frame import RoiManager
from ij.process import ByteProcessor, ImageProcessor
from java.awt import Window, Color, Font
from java.awt.image import BufferedImage
from java.io import File
from java.lang import Runnable
from java.lang import Throwable as JavaThrowable
from javax.imageio import ImageIO
from javax.swing import SwingUtilities, JDialog

ERRORS = (Exception, JavaThrowable)
REPO = os.environ.get("CQ_REPO") or os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMAGE = os.environ["CQ_DEMO_IMAGE"]
PROFILE = os.environ.get("CQ_DEMO_PROFILE") or os.path.join(REPO, "examples", "Demo_4Channel_HOPX_PAX6_DAPI_TBR2.properties")
RESULTS = os.environ["CQ_RESULTS"]
SHOTS = os.environ["CQ_SHOTS"]
LOOKS = ["Light", "Dark"]
MICRO = b"\xc2\xb5".decode("utf-8")
for folder in [RESULTS, SHOTS]:
    if not os.path.isdir(folder):
        os.makedirs(folder)


class Task(Runnable):
    def __init__(self, function):
        self.function = function
        self.done = threading.Event()
        self.error = None
        self.value = None

    def run(self):
        try:
            self.value = self.function()
        except ERRORS as exc:
            self.error = exc
        self.done.set()


def ui(function):
    task = Task(function)
    SwingUtilities.invokeLater(task)
    task.done.wait(180)
    if task.error is not None:
        raise task.error
    return task.value


def find_frame():
    for window in Window.getWindows():
        try:
            if window.isDisplayable() and window.getTitle().startswith("Cell Quant Studio") and hasattr(window, "run_analysis"):
                return window
        except ERRORS:
            pass
    return None


def settle(frame):
    for attempt in range(3):
        ui(frame.validate)
        time.sleep(0.2)


def paint(frame, name, before=None):
    def work():
        # The state is set and the picture taken in one go on the event
        # thread, so that no real mouse event can come between the two.
        if before is not None:
            before()
        frame.validate()
        root = frame.getRootPane()
        transform = frame.getGraphicsConfiguration().getDefaultTransform()
        scale = max(2.0, transform.getScaleX())
        image = BufferedImage(int(round(root.getWidth() * scale)), int(round(root.getHeight() * scale)), BufferedImage.TYPE_INT_RGB)
        graphics = image.createGraphics()
        graphics.scale(scale, scale)
        root.paint(graphics)
        graphics.dispose()
        ImageIO.write(image, "png", File(os.path.join(SHOTS, name + ".png")))
        return (image.getWidth(), image.getHeight())
    size = ui(work)
    print("shot %s %dx%d" % (name, size[0], size[1]))


def in_both_looks(frame, name, before=None):
    """Paint the window in each look; before() sets up the state again
    after the look changed (a change of look lays the tabs out anew)."""
    for look in LOOKS:
        ui(lambda: frame.set_look(look))
        time.sleep(0.6)
        ui(lambda: frame.setSize(1440, 900))
        if before is not None:
            ui(before)
        settle(frame)
        paint(frame, look.lower() + "_" + name, before)


def close_dialogs():
    closed = 0
    for window in Window.getWindows():
        try:
            if isinstance(window, JDialog) and window.isShowing() and window.getTitle().startswith("Cell Quant Studio"):
                SwingUtilities.invokeLater(Task(window.dispose))
                closed += 1
        except ERRORS:
            pass
    return closed


def zones_of(tissue, width, height, pixel):
    """As demo_zones() of the one-click demo: three zones that do not overlap."""
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


def zones_picture(imp, tissue, zones, pixel):
    """The image with the tissue outline and the zones drawn on it: what the
    audience sees on the image window at this step of the demo."""
    copy = imp.duplicate()
    try:
        for channel in range(1, copy.getNChannels() + 1):
            copy.setC(channel)
            copy.setDisplayRange(150.0, 3600.0)
        copy.setDisplayMode(IJ.COMPOSITE)
        flat = copy.flatten()
    finally:
        copy.close()
    overlay = Overlay()
    for index in range(len(zones)):
        outline = zones[index].clone()
        outline.setStrokeColor(Color.white)
        outline.setStrokeWidth(9.0 if index == 0 else 6.0)
        overlay.add(outline)
    bounds = tissue.getBounds()
    middle = bounds.x + bounds.width / 2.0
    font = Font("SansSerif", Font.BOLD, 76)
    depths = [50.0, 150.0, 262.0]
    for index in range(len(zones)):
        label = TextRoi(0, 0, " " + zones[index].getName() + " ", font)
        label.setStrokeColor(Color.white)
        label.setFillColor(Color(0, 0, 0, 190))
        box = label.getBounds()
        label.setLocation(middle - box.width / 2.0, bounds.y + depths[index] / pixel - box.height / 2.0)
        overlay.add(label)
    flat.setOverlay(overlay)
    drawn = flat.flatten()
    small = drawn.getProcessor().resize(1744, 1744, True)
    ImageIO.write(small.getBufferedImage(), "png", File(os.path.join(SHOTS, "zones.png")))
    flat.close()
    drawn.close()
    print("shot zones 1744x1744")


def main():
    remembered = Prefs.get("cellquant.look", "")
    launcher = os.path.join(REPO, "Cell_Quant_Studio_final.py")
    namespace = {"__file__": launcher, "__name__": "__main__"}
    execfile(launcher, namespace)
    frame = None
    for attempt in range(200):
        frame = find_frame()
        if frame is not None and frame.isShowing():
            break
        time.sleep(0.1)
    if frame is None:
        raise RuntimeError("The window did not open.")
    try:
        ui(lambda: frame.setSize(1440, 900))
        ui(lambda: frame.main_tabs.setSelectedIndex(0))

        # 1. the image
        ui(lambda: frame.open_image_path(IMAGE, True))
        imp = frame.source_imp
        if imp is None:
            raise RuntimeError("The image did not open.")
        time.sleep(5.0)
        in_both_looks(frame, "1_open")

        # 2. the profile
        ui(lambda: frame.load_profile_from_path(PROFILE))

        def outputs():
            frame.output_dir_field.setText(RESULTS)
            frame.save_prepared_cb.setSelected(True)
            for box in [frame.show_summary_cb, frame.keep_prepared_open_cb, frame.keep_masks_open_cb]:
                box.setSelected(False)
            frame.main_tabs.setSelectedIndex(0)
        ui(outputs)
        time.sleep(2.0)
        ui(lambda: frame.append_status("Demo: image and settings profile loaded"))
        time.sleep(0.5)
        in_both_looks(frame, "2_profile")

        # 3. the tissue outline and three depth zones
        engine = namespace["AnalysisEngine"]()
        rois, info = namespace["detect_tissue_outline"](
            engine, imp, frame.collect_channel_configs(), frame.collect_global_options(), frame=1
        )
        if len(rois) == 0:
            raise RuntimeError("No tissue was found.")
        print("tissue outlines found: %d" % len(rois))
        pixel = imp.getCalibration().pixelWidth or 1.0
        zones = zones_of(rois[0], imp.getWidth(), imp.getHeight(), pixel)
        print("zones: " + ", ".join(["%s %d px wide" % (zone.getName().encode("ascii", "replace"), zone.getBounds().width) for zone in zones]))
        zones_picture(imp, rois[0], zones, pixel)
        manager = ui(lambda: RoiManager.getRoiManager())
        ui(lambda: manager.reset())
        for zone in zones:
            ui(lambda: manager.addRoi(zone))
        ui(lambda: frame.roi_source_combo.setSelectedItem("ROI Manager"))

        # 4. the analysis
        started = time.time()
        ui(frame.run_analysis)
        time.sleep(1.0)
        end = time.time() + 900
        while time.time() < end and (frame.background_active or close_dialogs() > 0):
            time.sleep(0.5)
            close_dialogs()
        print("analysis took %.0f s" % (time.time() - started))
        time.sleep(1.0)
        close_dialogs()

        canvas = frame.results_canvas
        results_index = ui(lambda: frame.main_tabs.indexOfTab("Results"))

        def to_results():
            # The pictures of the page show the organoid map (the Pitch view).
            if hasattr(frame, "set_results_view"):
                frame.set_results_view("Pitch", False)
                if getattr(frame, "results_animation", None) is not None:
                    frame.results_animation.timer.stop()
                canvas.progress = 1.0
            frame.main_tabs.setSelectedIndex(results_index)
            canvas.zoom = False
            canvas.point_at(-1000, -1000)
            frame.validate()
        names = [metric["name"] for metric in canvas.metrics()]
        print("metrics: " + " | ".join([name.encode("ascii", "replace") for name in names]))
        print("tiles: " + " | ".join([("%s=%s" % (tile["title"], tile["value"])).encode("ascii", "replace") for tile in canvas.data["tiles"]]))
        for metric in canvas.metrics():
            print("  %s: %s" % (metric["name"].encode("ascii", "replace"), ", ".join(["%.4g" % value if value is not None else "-" for value in metric["values"]])))
        first = 0
        for index in range(len(names)):
            if names[index].startswith("DAPI") and canvas.metrics()[index]["unit"] == "count":
                first = index

        def tables():
            # What the Results tab shows first: the Researcher view, plain tables.
            frame.set_results_view("Researcher", False)
            frame.main_tabs.setSelectedIndex(results_index)
            frame.validate()
        if hasattr(frame, "set_results_view"):
            in_both_looks(frame, "3_research", tables)

        def counts():
            to_results()
            canvas.choose_metric(first)
        in_both_looks(frame, "3_results", counts)

        def pointed():
            counts()
            shape = canvas.geometry()
            canvas.point_at(shape["x"] + shape["rings"][0][0], shape["y"])
        in_both_looks(frame, "4_point", pointed)

        marker = None
        for index in range(len(names)):
            if names[index].startswith("TBR2"):
                marker = index

        def other_marker():
            to_results()
            canvas.choose_metric(marker)
        if marker is not None:
            in_both_looks(frame, "5_marker", other_marker)

        def zoomed():
            to_results()
            canvas.choose_metric(first)
            canvas.zoom = True
            canvas.repaint()
        in_both_looks(frame, "6_picture", zoomed)
        print("RECORD-DONE")
    finally:
        try:
            ui(lambda: frame.dispose())
        except ERRORS:
            pass
        for image_id in (WindowManager.getIDList() or []):
            image = WindowManager.getImage(image_id)
            if image is not None:
                image.changes = False
                image.close()
        manager = RoiManager.getInstance()
        if manager is not None:
            manager.close()
        Prefs.set("cellquant.look", remembered)


main()

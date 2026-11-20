# Cell Quant Studio Universal - end-to-end test of the application window
#
# What this checks: the window (cellquantcore.py), started the way a user
# starts it (by running Cell_Quant_Studio_final.py), on small synthetic images
# whose correct answers are known. The test presses the real buttons, lets the
# real code show its message boxes, closes them, and compares what the window
# wrote with a direct call of the analysis engine on the same image.
#
# How to run it from a terminal:
#   tests/run_tests.sh /path/to/Fiji tests/gui_smoke.py
# In Fiji: open this file in the Script Editor, Language > Python (Jython), Run.
#
# Windows really open on the screen for about half a minute. Every window and
# image the test opens is closed again, results are written to a temporary
# folder that is removed afterwards, and no file chooser is ever opened.
# Do not use the mouse or keyboard in Fiji while it runs.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij import IJ, ImagePlus, ImageStack, CompositeImage, WindowManager, Prefs, ImageJ
from ij.process import ShortProcessor, ImageConverter
from ij.gui import Roi, OvalRoi, PolygonRoi, PointRoi
from ij.io import FileSaver
from ij.plugin.frame import RoiManager
from java.awt import Polygon, Window, Dialog, Component
from java.awt.event import ActionListener
from java.io import FileOutputStream
from java.lang import Runnable, Thread as JavaThread, Throwable as JavaThrowable
from javax.swing import UIManager
from javax.swing import SwingUtilities, JOptionPane, JFileChooser, JButton, JFrame, JPanel, JCheckBox, JComboBox, JTextField
from javax.swing import Timer as SwingTimer
from javax.swing.border import TitledBorder

import os
import math
import time
import shutil
import tempfile
import traceback
import codecs
import threading


# Under Jython, "except Exception" alone does not catch Java exceptions.
ANY_ERROR = (Exception, JavaThrowable)

MICRO = b"\xc2\xb5".decode("utf-8")
WIDTH = 240
HEIGHT = 200
PIXEL = 0.5
DISK_RADIUS = 6
DISKS = [
    (30, 30), (70, 30), (110, 30), (150, 30), (190, 30),
    (30, 80), (70, 80), (110, 80), (150, 80), (190, 80),
    (50, 140), (100, 140), (150, 140), (200, 140), (120, 175)
]
MARKER_OUTER = (40, 50, 100, 80)
MARKER_INNER = (60, 70, 60, 40)
RECT = (10, 10, 220, 180)
TAB_TITLES = ["Session", "Results", "Channels", "ROIs and output", "Batch", "Validation and QC", "Analyses", "Help"]
TAB_TITLES_NO_MODULES = ["Session", "Channels", "ROIs and output", "Batch", "Validation and QC", "Help"]
ACTION_TIMEOUT = 90.0

C1_VALUES = {
    "enabled": True, "name": "DAPI", "analysis_mode": "Objects + area", "lut": "Blue",
    "threshold_mode": "Manual", "manual_low": 1000, "manual_high": -1, "cross_measure": True
}
C2_VALUES = {
    "enabled": True, "name": "MAP2", "analysis_mode": "Positive area", "lut": "Green",
    "threshold_mode": "Manual", "manual_low": 800, "manual_high": -1
}
C3_VALUES = {"enabled": True, "name": "Channel_3", "analysis_mode": "Intensity only"}


class Hang(Exception):
    """An action did not return in time: the window is blocked."""
    pass


class Namespace(object):
    def __init__(self, values):
        self.__dict__.update(values)


class Context(object):
    pass


class Recorder(object):
    def __init__(self, name):
        self.name = name
        self.checks = 0
        self.failures = []

    def check(self, condition, message):
        self.checks += 1
        if not condition:
            self.failures.append(message)

    def equal(self, actual, expected, message):
        self.check(actual == expected, message + " (expected " + repr(expected) + ", got " + repr(actual) + ")")

    def close(self, actual, expected, message, tolerance=1e-6):
        ok = False
        try:
            ok = abs(float(actual) - float(expected)) <= tolerance * max(1.0, abs(float(expected)))
        except ANY_ERROR:
            ok = False
        self.check(ok, message + " (expected " + repr(expected) + ", got " + repr(actual) + ")")


def as_text(value):
    if value is None:
        return u""
    if isinstance(value, unicode):
        return value
    try:
        return unicode(value)
    except ANY_ERROR:
        return unicode(repr(value))


# -----------------------------------------------------------------------------
# Threads: running an action on the event thread or on a script thread
# -----------------------------------------------------------------------------

class EdtTask(Runnable):
    def __init__(self, function):
        self.function = function

    def run(self):
        self.function()


def thread_dump():
    lines = ["Thread dump:"]
    traces = JavaThread.getAllStackTraces()
    for thread in traces.keySet():
        frames = traces.get(thread)
        if frames is None or len(frames) == 0:
            continue
        lines.append("  " + as_text(thread.getName()) + " [" + as_text(thread.getState()) + "]")
        for index in range(min(len(frames), 30)):
            lines.append("      at " + as_text(frames[index]))
    return "\n".join(lines)


def run_guarded(function, on_event_thread, timeout, what):
    """Run function on the Swing event thread or on a fresh script thread and
    wait for it. A Python or Java error inside it is raised here; if it does
    not return in time the threads are printed and Hang is raised."""
    box = {"value": None, "error": None}
    done = threading.Event()

    def body():
        try:
            box["value"] = function()
        except ANY_ERROR:
            box["error"] = traceback.format_exc()
        finally:
            done.set()

    if on_event_thread:
        SwingUtilities.invokeLater(EdtTask(body))
    else:
        worker = threading.Thread(target=body, name="cellquant-script-thread")
        worker.setDaemon(True)
        worker.start()
    done.wait(timeout)
    if not done.isSet():
        print(thread_dump())
        raise Hang(what + " did not return within " + str(timeout) + " s")
    if box["error"] is not None:
        raise RuntimeError(what + " raised:\n" + box["error"])
    return box["value"]


def on_edt(function, what="event-thread action", timeout=ACTION_TIMEOUT):
    return run_guarded(function, True, timeout, what)


def off_edt(function, what="script-thread action", timeout=ACTION_TIMEOUT):
    return run_guarded(function, False, timeout, what)


def nothing():
    return None


def wait_until(condition, timeout, what):
    end = time.time() + timeout
    while time.time() < end:
        if condition():
            return
        time.sleep(0.01)
    print(thread_dump())
    raise Hang("Timed out waiting for " + what)


def find_components(container, kind, found=None):
    if found is None:
        found = []
    try:
        children = container.getComponents()
    except ANY_ERROR:
        children = []
    for child in children:
        if isinstance(child, kind):
            found.append(child)
        find_components(child, kind, found)
    return found


# -----------------------------------------------------------------------------
# Modal dialogs: closed by a background thread, which records what they said
# -----------------------------------------------------------------------------

class DialogWatcher(object):
    LEVELS = {
        JOptionPane.ERROR_MESSAGE: "error", JOptionPane.INFORMATION_MESSAGE: "info",
        JOptionPane.WARNING_MESSAGE: "warning", JOptionPane.QUESTION_MESSAGE: "question",
        JOptionPane.PLAIN_MESSAGE: "plain"
    }

    def __init__(self):
        self.records = []
        self.answers = []
        self.seen = []
        self.stopped = False
        self.thread = threading.Thread(target=self.loop, name="cellquant-dialog-watcher")
        self.thread.setDaemon(True)

    def start(self):
        self.thread.start()

    def stop(self):
        self.stopped = True

    def answer_next(self, answer):
        """Queue 'yes' or 'no' for the next confirmation. The default is yes."""
        self.answers.append(answer)

    def take(self):
        taken = self.records
        self.records = []
        return taken

    def loop(self):
        while not self.stopped:
            try:
                self.scan()
            except ANY_ERROR:
                self.records.append({"kind": "watcher error", "level": "error", "title": "", "text": as_text(traceback.format_exc())})
            time.sleep(0.03)

    def showing_modal_dialogs(self):
        result = []
        for window in Window.getWindows():
            if isinstance(window, Dialog) and window.isShowing() and window.isModal():
                result.append(window)
        return result

    def scan(self):
        # Forget closed dialogs first: a dialog that closes between these two
        # lines is then neither forgotten while still listed as showing nor
        # handled twice.
        self.seen = [window for window in self.seen if window.isShowing()]
        showing = self.showing_modal_dialogs()
        for dialog in showing:
            known = False
            for seen in self.seen:
                if seen is dialog:
                    known = True
            if not known:
                self.seen.append(dialog)
                self.handle(dialog)

    def handle(self, dialog):
        record = {"kind": "other", "level": "plain", "title": as_text(dialog.getTitle()), "text": ""}
        choosers = find_components(dialog, JFileChooser)
        panes = find_components(dialog, JOptionPane)
        if len(choosers) > 0:
            record["kind"] = "file chooser"
            closer = choosers[0].cancelSelection
        elif len(panes) > 0:
            pane = panes[0]
            record["text"] = as_text(pane.getMessage())
            record["level"] = self.LEVELS.get(pane.getMessageType(), "plain")
            value = JOptionPane.OK_OPTION
            if pane.getOptionType() == JOptionPane.DEFAULT_OPTION:
                record["kind"] = "message"
            else:
                record["kind"] = "confirm"
                answer = "yes"
                if len(self.answers) > 0:
                    answer = self.answers.pop(0)
                record["answer"] = answer
                if answer != "yes":
                    if pane.getOptionType() == JOptionPane.OK_CANCEL_OPTION:
                        value = JOptionPane.CANCEL_OPTION
                    else:
                        value = JOptionPane.NO_OPTION

            def closer(pane=pane, value=value):
                pane.setValue(value)
        else:
            closer = dialog.dispose
        self.records.append(record)
        SwingUtilities.invokeLater(EdtTask(closer))


def settle(ctx, timeout=20.0):
    """Wait until the event queue is empty and no modal dialog is showing."""
    end = time.time() + timeout
    while True:
        on_edt(nothing, "event queue flush", timeout)
        if len(ctx.watcher.showing_modal_dialogs()) == 0:
            on_edt(nothing, "event queue flush", timeout)
            if len(ctx.watcher.showing_modal_dialogs()) == 0:
                return
        if time.time() > end:
            print(thread_dump())
            raise Hang("A modal dialog stayed open")
        time.sleep(0.03)


def wait_for_job(ctx, timeout=120.0):
    """Wait until a background job of the window (an analysis, validation,
    sensitivity or Analyses-tab run) has finished, its closing message box
    included. The batch has its own wait_for_batch."""
    if ctx.app is not None:
        wait_until(lambda: not ctx.app.background_active, timeout, "the background job to finish")


def expect_dialogs(ctx, t, expected, label):
    """expected: list of (level, text fragment), in the order the dialogs appeared."""
    wait_for_job(ctx)
    settle(ctx)
    records = ctx.watcher.take()
    described = [record["level"] + ": " + record["text"].replace("\n", " | ") for record in records]
    t.equal(len(records), len(expected), label + ": number of message boxes " + repr(described))
    for index in range(min(len(records), len(expected))):
        level, fragment = expected[index]
        record = records[index]
        t.check(
            record["level"] == level and fragment in record["text"],
            label + ": message box " + str(index + 1) + " should be " + level + " containing '" + fragment +
            "' but was " + described[index]
        )
    return records


def click(ctx, text, container=None, current_image=None, timeout=ACTION_TIMEOUT):
    """Press the one button with this label, on the event thread."""
    root = container if container is not None else ctx.app
    buttons = [button for button in find_components(root, JButton) if as_text(button.getText()) == text]
    if len(buttons) != 1:
        raise RuntimeError("Expected exactly one button labelled '" + text + "', found " + str(len(buttons)))
    button = buttons[0]

    def press():
        if not button.isEnabled():
            raise RuntimeError("Button '" + text + "' is disabled")
        if current_image is not None:
            WindowManager.setTempCurrentImage(current_image)
        try:
            button.doClick(0)
        finally:
            if current_image is not None:
                WindowManager.setTempCurrentImage(None)

    on_edt(press, "button '" + text + "'", timeout)
    settle(ctx)


def script_call(ctx, function, what, timeout=ACTION_TIMEOUT):
    """Call a window method from a script thread, as a Fiji script would."""
    value = off_edt(function, what, timeout)
    settle(ctx)
    return value


def set_widgets(panel, values):
    """Type values into the widgets of one channel panel (event thread)."""
    checkboxes = panel.checkbox_widgets()
    combos = panel.combo_widgets()
    fields = panel.text_widgets()

    def apply():
        for key in values:
            if key in checkboxes:
                checkboxes[key].setSelected(bool(values[key]))
            elif key in combos:
                combos[key].setSelectedItem(values[key])
                if as_text(combos[key].getSelectedItem()) != as_text(values[key]):
                    raise RuntimeError("Choice '" + as_text(values[key]) + "' is not offered for " + key)
            elif key in fields:
                fields[key].setText(as_text(values[key]))
            else:
                raise RuntimeError("No widget for channel setting " + key)

    on_edt(apply, "typing channel settings")


# -----------------------------------------------------------------------------
# Synthetic images with known answers (same construction as run_tests.py)
# -----------------------------------------------------------------------------

def blank_plane(value):
    ip = ShortProcessor(WIDTH, HEIGHT)
    ip.setValue(float(value))
    ip.fill()
    return ip


def nuclei_plane(disks):
    ip = blank_plane(100)
    ip.setValue(2000.0)
    for disk in disks:
        ip.fill(OvalRoi(disk[0] - DISK_RADIUS, disk[1] - DISK_RADIUS, 2 * DISK_RADIUS + 1, 2 * DISK_RADIUS + 1))
    ip.resetRoi()
    return ip


def marker_plane():
    ip = blank_plane(100)
    ip.setValue(900.0)
    ip.fill(Roi(MARKER_OUTER[0], MARKER_OUTER[1], MARKER_OUTER[2], MARKER_OUTER[3]))
    ip.setValue(1500.0)
    ip.fill(Roi(MARKER_INNER[0], MARKER_INNER[1], MARKER_INNER[2], MARKER_INNER[3]))
    ip.resetRoi()
    return ip


def make_image(title, disks=None, slices=2, frames=1, disks_per_frame=None):
    """Three channels: C1 bright disks ('nuclei'), C2 a two-level rectangle
    ('marker'), C3 a flat plane whose value is 200 + 100 * Z. Calibrated at
    0.5 micrometer per pixel."""
    stack = ImageStack(WIDTH, HEIGHT)
    for t in range(frames):
        frame_disks = DISKS if disks is None else disks
        if disks_per_frame is not None:
            frame_disks = disks_per_frame[t]
        for z in range(slices):
            stack.addSlice("", nuclei_plane(frame_disks))
            stack.addSlice("", marker_plane())
            stack.addSlice("", blank_plane(200 + 100 * (z + 1)))
    imp = ImagePlus(title, stack)
    imp.setDimensions(3, slices, frames)
    imp = CompositeImage(imp, CompositeImage.GRAYSCALE)
    cal = imp.getCalibration()
    cal.pixelWidth = PIXEL
    cal.pixelHeight = PIXEL
    cal.setUnit("um")
    return imp


def triangle_roi():
    polygon = Polygon()
    polygon.addPoint(100, 30)
    polygon.addPoint(170, 170)
    polygon.addPoint(30, 170)
    roi = PolygonRoi(polygon, Roi.POLYGON)
    roi.setName("triangle")
    return roi


def rect_roi():
    roi = Roi(RECT[0], RECT[1], RECT[2], RECT[3])
    roi.setName("rect")
    return roi


def save_image(imp, path):
    saver = FileSaver(imp)
    if imp.getStackSize() > 1:
        return saver.saveAsTiffStack(path)
    return saver.saveAsTiff(path)


def read_csv(path):
    """Read a table written by the application: every cell is quoted."""
    handle = codecs.open(path, "r", "utf-8")
    try:
        text = handle.read()
    finally:
        handle.close()
    table = []
    for line in text.split("\n"):
        if line.strip() == "":
            continue
        cells = []
        current = []
        quoted = False
        index = 0
        while index < len(line):
            char = line[index]
            if quoted:
                if char == u'"' and index + 1 < len(line) and line[index + 1] == u'"':
                    current.append(u'"')
                    index += 1
                elif char == u'"':
                    quoted = False
                else:
                    current.append(char)
            elif char == u'"':
                quoted = True
            elif char == u",":
                cells.append(u"".join(current))
                current = []
            elif char != u"\r":
                current.append(char)
            index += 1
        cells.append(u"".join(current))
        table.append(cells)
    headers = table[0]
    rows = []
    for cells in table[1:]:
        row = {}
        for index in range(len(headers)):
            row[headers[index]] = cells[index] if index < len(cells) else u""
        rows.append(row)
    return headers, rows


def compare_tables(t, label, actual, expected, skip):
    """Both arguments are (headers, rows) from read_csv."""
    t.equal(actual[0], expected[0], label + ": same columns in the same order")
    t.equal(len(actual[1]), len(expected[1]), label + ": same number of rows")
    different = []
    for index in range(min(len(actual[1]), len(expected[1]))):
        for header in expected[0]:
            if header in skip:
                continue
            if actual[1][index].get(header) != expected[1][index].get(header):
                different.append(
                    "row " + str(index + 1) + " " + header + ": window " + repr(actual[1][index].get(header)) +
                    ", engine " + repr(expected[1][index].get(header))
                )
    t.check(len(different) == 0, label + ": window and engine numbers differ: " + "; ".join(different[:8]))


def newest_dir(root, marker):
    best = None
    for name in os.listdir(root):
        path = os.path.join(root, name)
        if os.path.isdir(path) and marker in name:
            if best is None or os.path.getmtime(path) >= os.path.getmtime(best):
                best = path
    return best


def expected_configs(ctx, overrides=None):
    E = ctx.E
    values = [dict(C1_VALUES), dict(C2_VALUES), dict(C3_VALUES)]
    if overrides:
        for channel in overrides:
            values[channel - 1].update(overrides[channel])
    return [E.normalize_channel_config(index + 1, values[index]) for index in range(3)]


def engine_reference(ctx, imp, source_path, configs, rois, frames, global_overrides=None):
    """Run the engine directly, with no window, and return its tables."""
    E = ctx.E
    options = E.default_global_options()
    if global_overrides:
        options.update(global_overrides)
    run_dir = E.unique_run_dir(os.path.join(ctx.tmp, "reference"), E.strip_extension(imp.getTitle()))
    E.AnalysisEngine().analyze_image(
        imp, source_path, configs, E.validate_parent_rois(rois, imp.getWidth(), imp.getHeight()),
        E.engine_options_from_global(options), run_dir, frames, settings=(options, configs)
    )
    return {
        "summary": read_csv(os.path.join(run_dir, "roi_summary.csv")),
        "objects": read_csv(os.path.join(run_dir, "objects.csv")),
        "thresholds": read_csv(os.path.join(run_dir, "thresholds.csv"))
    }


def window_tables(run_dir):
    return {
        "summary": read_csv(os.path.join(run_dir, "roi_summary.csv")),
        "objects": read_csv(os.path.join(run_dir, "objects.csv")),
        "thresholds": read_csv(os.path.join(run_dir, "thresholds.csv"))
    }


def compare_run(t, label, window, engine):
    compare_tables(t, label + " roi_summary.csv", window["summary"], engine["summary"], ["Run_ID"])
    compare_tables(t, label + " objects.csv", window["objects"], engine["objects"], ["Run_ID"])
    compare_tables(t, label + " thresholds.csv", window["thresholds"], engine["thresholds"], ["Run_ID"])


def row_named(rows, name):
    for row in rows:
        if row.get("Parent_ROI_Name") == name:
            return row
    return {}


def use_as_source(ctx, imp):
    """Show an image and make it the source through 'Use active image'."""
    if imp.getWindow() is None:
        imp.show()
    click(ctx, "Use active image", current_image=imp)


def find_app_windows():
    found = []
    for window in Window.getWindows():
        if isinstance(window, JFrame) and window.isDisplayable():
            if as_text(window.getTitle()).startswith("Cell Quant Studio Universal"):
                found.append(window)
    return found


# -----------------------------------------------------------------------------
# Steps
# -----------------------------------------------------------------------------

def step_launch(ctx, t):
    namespace = {"__file__": ctx.launcher, "__name__": "__main__"}
    off_edt(lambda: execfile(ctx.launcher, namespace), "running the launcher")
    ctx.ns = namespace
    ctx.E = Namespace(namespace)
    windows = find_app_windows()
    t.equal(len(windows), 1, "the launcher opens exactly one application window")
    ctx.app = windows[0]
    settle(ctx)
    t.check(ctx.app.isShowing(), "the window is visible")
    t.equal(as_text(ctx.app.getTitle()), ctx.E.APP_NAME + " " + ctx.E.APP_VERSION, "window title")
    tabs = ctx.app.main_tabs
    t.equal([as_text(tabs.getTitleAt(i)) for i in range(tabs.getTabCount())], TAB_TITLES,
            "the tabs, with Analyses because modules are loaded")
    t.equal(ctx.app.channel_tabs.getTabCount(), ctx.E.MAX_CHANNELS, "one channel tab per supported channel")
    log = as_text(ctx.app.status_area.getText())
    t.check("Ready." in log and "Applied preset" in log, "start-up messages reached the status log")
    t.equal(ctx.E.CORE_VERSION, ctx.E.APP_VERSION, "window and engine versions match")
    again = off_edt(namespace["launch_cell_quant_studio"], "launching a second time")
    t.check(again is ctx.app, "launching again returns the open window")
    t.equal(len(find_app_windows()), 1, "launching again opens no second window")

    # The optional analysis modules beside the engine are loaded by the launcher.
    E = ctx.E
    present = [name for name in E.ENGINE_MODULE_FILES if os.path.isfile(os.path.join(ctx.repo, name))]
    t.equal(namespace.get("CELLQUANT_LOADED_MODULES"), present, "the launcher loads every analysis module found beside it")
    for name in present:
        t.equal(E.ENGINE_MODULE_VERSIONS.get(name), E.APP_VERSION, name + " has the version of the engine")
    combo = ctx.app.batch_roi_combo
    t.equal([as_text(combo.getItemAt(i)) for i in range(combo.getItemCount())], [as_text(p) for p in E.BATCH_ROI_POLICIES],
            "the batch tab offers every parent-ROI policy, including those added by modules")
    expect_dialogs(ctx, t, [], "launch")


def step_module_version(ctx, t):
    """A module from another release stops the launcher before any window opens."""
    folder = os.path.join(ctx.tmp, "mixed release")
    os.makedirs(folder)
    for name in ["Cell_Quant_Studio_final.py", "cellquantengine.py", "cellquantcore.py"]:
        shutil.copy(os.path.join(ctx.repo, name), os.path.join(folder, name))
    handle = open(os.path.join(folder, "cellquantdepth.py"), "w")
    try:
        handle.write('register_engine_module("cellquantdepth.py", "0.0.1")\n')
    finally:
        handle.close()
    launcher = os.path.join(folder, "Cell_Quant_Studio_final.py")
    namespace = {"__file__": launcher, "__name__": "__main__"}
    message = ""
    try:
        off_edt(lambda: execfile(launcher, namespace), "running a launcher with a module of another version")
    except Hang:
        raise
    except ANY_ERROR as exc:
        message = as_text(exc)
    t.check("cellquantdepth.py is version 0.0.1" in message, "the refusal names the module and its version: " + message[-300:])
    t.check("Re-extract the release" in message, "the refusal says what to do")
    t.check("launch_cell_quant_studio" not in namespace, "the window file is not loaded")
    t.equal(len(find_app_windows()), 1, "no second window opens")
    expect_dialogs(ctx, t, [], "module of another version")


def step_status_messages(ctx, t):
    app = ctx.app
    name = "caf" + b"\xc3\xa9".decode("utf-8") + " " + MICRO + "m"
    on_edt(lambda: app.append_status("event thread " + name), "status message on the event thread")
    off_edt(lambda: app.append_status("script thread " + name), "status message on a script thread")
    settle(ctx)
    log = as_text(app.status_area.getText())
    t.check(("event thread " + name) in log, "a message sent on the event thread is shown, non-ASCII included")
    t.check(("script thread " + name) in log, "a message sent from a script thread is shown, non-ASCII included")

    off_edt(lambda: app.report_exception("Script thread check", ValueError("bad " + name)), "error report from a script thread")
    expect_dialogs(ctx, t, [("error", "Script thread check: bad " + name)], "error reported from a script thread")
    on_edt(lambda: app.report_exception("Event thread check", ValueError("bad " + name)), "error report on the event thread")
    expect_dialogs(ctx, t, [("error", "Event thread check: bad " + name)], "error reported on the event thread")

    click(ctx, "1. Prepare enabled channels")
    expect_dialogs(ctx, t, [("error", "No source image is selected")], "preparing with no image")
    script_call(ctx, app.run_analysis_action, "run analysis with no image")
    expect_dialogs(ctx, t, [("error", "No source image is selected")], "analysis with no image, from a script thread")


def step_source_image(ctx, t):
    app = ctx.app
    ctx.imp1 = make_image("smoke one.tif")
    use_as_source(ctx, ctx.imp1)
    t.check(app.source_imp is ctx.imp1, "'Use active image' takes the active image as the source")
    info = as_text(app.image_info_area.getText())
    t.check("Channels: 3" in info and "Z slices: 2" in info, "the source description lists channels and slices")
    t.equal([app.channel_tabs.isEnabledAt(i) for i in range(8)], [True, True, True, False, False, False, False, False],
            "channel tabs beyond the image's channels are disabled")
    t.equal(as_text(app.pixel_width_field.getText()), "0.5", "pixel width is read from the image")
    # ImageJ itself turns the unit "um" into the micro sign followed by "m".
    t.equal(as_text(app.spatial_unit_field.getText()), MICRO + "m", "unit is read from the image")
    expect_dialogs(ctx, t, [], "selecting the source")


def step_configure(ctx, t):
    app = ctx.app
    E = ctx.E
    on_edt(lambda: app.preset_combo.setSelectedItem("Intensity-only quantification"), "choosing a preset")
    click(ctx, "Apply preset")
    configs = app.collect_channel_configs()
    t.equal([config["channel"] for config in configs], [1, 2, 3], "the intensity preset enables the image's three channels")
    t.equal([config["analysis_mode"] for config in configs], ["Intensity only"] * 3, "the preset sets the analysis mode")

    set_widgets(app.channel_panels[0], C1_VALUES)
    set_widgets(app.channel_panels[1], C2_VALUES)
    on_edt(lambda: app.output_dir_field.setText(ctx.out_dir), "typing the output folder")
    configs = app.collect_channel_configs()
    expected = expected_configs(ctx)
    t.equal(len(configs), 3, "three channels are enabled")
    for index in range(min(len(configs), 3)):
        for key in expected[index]:
            t.equal(configs[index][key], expected[index][key], "C" + str(index + 1) + " widget value '" + key + "'")

    set_widgets(app.channel_panels[0], {"manual_low": "abc"})
    click(ctx, "1. Prepare enabled channels")
    expect_dialogs(ctx, t, [("error", "must be numeric")], "a non-numeric threshold")
    set_widgets(app.channel_panels[0], {"manual_low": 1000})

    manager = RoiManager.getRoiManager()
    manager.reset()
    manager.addRoi(triangle_roi())
    manager.addRoi(rect_roi())
    ctx.manager = manager
    t.equal(manager.getCount(), 2, "two parent ROIs are in ROI Manager")
    on_edt(lambda: app.roi_source_combo.setSelectedItem("ROI Manager"), "choosing the ROI source")
    rois = app.get_parent_rois()
    t.equal([as_text(roi.getName()) for roi in rois], ["triangle", "rect"], "the window reads both parent ROIs")
    expect_dialogs(ctx, t, [], "configuration")


def step_prepare_and_preview(ctx, t):
    app = ctx.app
    click(ctx, "1. Prepare enabled channels")
    expect_dialogs(ctx, t, [], "preparing the channels")
    t.equal(sorted(app.prepared_images.keys()), [1, 2, 3], "three prepared images")
    for channel in sorted(app.prepared_images.keys()):
        imp = app.prepared_images[channel]
        t.check(WindowManager.getImage(imp.getID()) is not None, "prepared C" + str(channel) + " is shown")
        t.check(as_text(imp.getTitle()).endswith("_PREP"), "prepared C" + str(channel) + " is named _PREP")
        t.equal(imp.getStackSize(), 1, "prepared C" + str(channel) + " is one projected plane")
    t.check(app.prepared_cache_valid(app.collect_channel_configs(), 1), "the prepared images are current")
    t.equal(app.prepared_images[3].getStatistics().max, 400.0, "the maximum projection of C3 holds the brighter slice")
    ctx.prepared_ids = [app.prepared_images[channel].getID() for channel in [1, 2, 3]]

    click(ctx, "2. Preview all masks")
    expect_dialogs(ctx, t, [], "previewing the masks")
    t.equal(sorted([channel for channel in app.mask_preview_ids if len(app.mask_preview_ids[channel]) > 0]), [1, 2], "masks for the two segmented channels")
    for channel in [1, 2]:
        for image_id in app.mask_preview_ids.get(channel, []):
            t.check(WindowManager.getImage(image_id) is not None, "mask preview of C" + str(channel) + " is shown")
    t.equal([app.prepared_images[channel].getID() for channel in [1, 2, 3]], ctx.prepared_ids, "previewing reuses the prepared images")

    click(ctx, "Preview mask", container=app.channel_panels[1])
    expect_dialogs(ctx, t, [], "previewing one mask")
    click(ctx, "Preview mask", container=app.channel_panels[2])
    expect_dialogs(ctx, t, [("error", "Intensity only")], "previewing a mask of an intensity-only channel")
    click(ctx, "Prepare / refresh", container=app.channel_panels[0])
    expect_dialogs(ctx, t, [], "refreshing one channel")
    t.equal([app.prepared_images[channel].getID() for channel in [1, 2, 3]], ctx.prepared_ids, "refreshing reuses current prepared images")


def step_run_analysis(ctx, t):
    app = ctx.app
    E = ctx.E
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("info", "Analysis complete.")], "run analysis")
    run_dir = as_text(app.last_run_dir)
    t.check(os.path.isdir(run_dir) and run_dir.startswith(ctx.out_dir), "the run folder is inside the chosen output folder")
    t.check("smoke_one_CellQuant_" in os.path.basename(run_dir), "the run folder is named after the image")
    for name in ["roi_summary.csv", "objects.csv", "thresholds.csv", "settings.properties", "run_log.txt"]:
        t.check(os.path.isfile(os.path.join(run_dir, name)), name + " is written")
    t.check(os.path.isfile(os.path.join(run_dir, "ROIs", "Parent_ROIs_Input.zip")), "the parent ROIs are saved")
    t.equal(len(os.listdir(os.path.join(run_dir, "Prepared"))), 3, "three prepared images are saved")
    t.equal(len(os.listdir(os.path.join(run_dir, "Overlays"))), 1, "one object overlay is saved")
    log_handle = codecs.open(os.path.join(run_dir, "run_log.txt"), "r", "utf-8")
    try:
        log_text = log_handle.read()
    finally:
        log_handle.close()
    t.check("Started analysis run" in log_text and "Analysis completed" in log_text, "run_log.txt holds the messages of this run")

    window = window_tables(run_dir)
    reference_image = make_image("smoke one.tif")
    engine = engine_reference(ctx, reference_image, "", expected_configs(ctx), [triangle_roi(), rect_roi()], [1])
    reference_image.close()
    compare_run(t, "first image", window, engine)
    rect = row_named(window["summary"][1], "rect")
    t.equal(rect.get("C1_DAPI_Object_Count"), "15", "all 15 disks are counted in the rectangle")
    t.close(rect.get("C2_MAP2_Positive_Area"), MARKER_OUTER[2] * MARKER_OUTER[3] * PIXEL * PIXEL, "marker area in the rectangle")
    t.close(rect.get("ROI_Area"), RECT[2] * RECT[3] * PIXEL * PIXEL, "rectangle area")
    t.close(rect.get("C3_Channel_3_Mean"), 400.0, "C3 mean is the maximum projection")
    triangle = row_named(window["summary"][1], "triangle")
    t.check(0 < int(triangle.get("C1_DAPI_Object_Count", "0")) < 15, "the triangle holds fewer objects than its bounding box")

    table = WindowManager.getWindow("CellQuant ROI Summary")
    t.check(table is not None, "the summary table window is shown")
    if table is not None:
        t.equal(table.getTextPanel().getLineCount(), 2, "the summary table has one line per parent ROI")
    t.check(RoiManager.getInstance2() is ctx.manager, "ROI Manager is still connected after the analysis")
    t.equal(ctx.manager.getCount(), 2, "ROI Manager still holds the two parent ROIs")
    t.equal([app.prepared_images[channel].getID() for channel in [1, 2, 3]], ctx.prepared_ids, "the analysis used the prepared images")
    for image_id in ctx.prepared_ids:
        t.check(WindowManager.getImage(image_id) is not None, "prepared images stay open after the analysis")

    # The Results tab shows the same run: its numbers come from the summary table.
    tabs = app.main_tabs
    t.equal(as_text(tabs.getTitleAt(tabs.getSelectedIndex())), "Results", "the Results tab comes to the front after an analysis")
    shown = app.results_canvas.data
    t.equal(shown["regions"], 2, "the Results tab covers both parent ROIs")
    counted = sum([int(row.get("C1_DAPI_Object_Count", "0")) for row in window["summary"][1]])
    count_tiles = [tile for tile in shown["tiles"] if tile["kind"] == "count"]
    t.equal([(tile["title"], tile["value"]) for tile in count_tiles], [("DAPI objects", float(counted))],
            "the object tile is the sum of the summary table")
    positive = sum([float(row.get("C2_MAP2_Positive_Area")) for row in window["summary"][1]])
    measured = sum([float(row.get("ROI_Area")) for row in window["summary"][1]])
    marker_tiles = [tile for tile in shown["tiles"] if tile["title"] == "MAP2 positive area"]
    t.equal(len(marker_tiles), 1, "one positive-area tile for the marker")
    if len(marker_tiles) == 1:
        t.close(marker_tiles[0]["value"], positive / measured * 100.0, "the marker tile is positive area over measured area")
    t.equal(shown["bars"]["labels"], [row.get("Parent_ROI_Name") for row in window["summary"][1]], "one group of bars per parent ROI")
    captions = [as_text(picture[0]) for picture in app.results_canvas.pictures]
    t.equal(len(captions), 1, "the merged picture of the run is shown")
    t.check(len(captions) == 1 and captions[0].startswith("All channels") and (str(counted) + " detected objects") in captions[0],
            "the merged picture outlines every counted object")
    canvas = app.results_canvas
    t.equal(len(canvas.pictures[0][2]), 2, "both parent ROIs are kept as outlines for the picture")
    shape = canvas.geometry()
    t.equal(len(shape["rings"]), 2, "one ring per parent ROI")
    for ring_index in range(len(shape["rings"])):
        radius = shape["rings"][ring_index][0]
        t.equal(canvas.ring_at(shape["x"] + radius, shape["y"]), ring_index, "pointing at ring " + str(ring_index + 1) + " finds it")
    t.equal(canvas.ring_at(shape["x"], shape["y"]), None, "the picture in the middle is no ring")
    radii = [ring[0] for ring in shape["rings"]]
    t.check(radii == sorted(radii, reverse=True) and shape["picture"] < radii[-1], "rings lie around the picture, without overlap")
    on_edt(lambda: canvas.point_at(shape["x"] + radii[0], shape["y"]), "pointing at the first ring")
    t.equal(canvas.hover, 0, "the pointed ring is highlighted")
    on_edt(lambda: canvas.point_at(-1000, -1000), "moving the pointer away")
    t.equal(canvas.hover, None, "and released when the pointer leaves")
    t.check(app.source_preview.picture is not None, "the Session tab shows a small picture of the source image")
    t.equal(as_text(tabs.getTitleAt(0)), "Session", "the Session tab is still the first one")
    on_edt(lambda: tabs.setSelectedIndex(0), "back to the Session tab")


def step_second_image(ctx, t):
    """A second image, analyzed from a script thread, without touching ROI Manager."""
    app = ctx.app
    old_ids = list(app.generated_ids)
    ctx.imp2 = make_image("smoke two.tif", disks=DISKS[:7])
    ctx.imp2.show()

    def select():
        WindowManager.setTempCurrentImage(ctx.imp2)
        try:
            app.use_active_image()
        finally:
            WindowManager.setTempCurrentImage(None)

    script_call(ctx, select, "'Use active image' from a script thread")
    t.check(app.source_imp is ctx.imp2, "the second image is the source")
    t.check(len(old_ids) > 0, "the first image had generated images")
    t.equal([image_id for image_id in old_ids if WindowManager.getImage(image_id) is not None], [],
            "images generated for the first image are closed")
    t.equal(len(app.prepared_images), 0, "no prepared image of the first image is kept")

    script_call(ctx, app.run_analysis_action, "run analysis from a script thread")
    expect_dialogs(ctx, t, [("info", "Analysis complete.")], "second image")
    t.check(RoiManager.getInstance2() is ctx.manager and ctx.manager.getCount() == 2, "ROI Manager is untouched")
    run_dir = as_text(app.last_run_dir)
    t.check("smoke_two_CellQuant_" in os.path.basename(run_dir), "a new run folder for the second image")
    window = window_tables(run_dir)
    reference_image = make_image("smoke two.tif", disks=DISKS[:7])
    engine = engine_reference(ctx, reference_image, "", expected_configs(ctx), [triangle_roi(), rect_roi()], [1])
    reference_image.close()
    compare_run(t, "second image", window, engine)
    t.equal(row_named(window["summary"][1], "rect").get("C1_DAPI_Object_Count"), "7", "seven disks in the second image")

    script_call(ctx, app.quick_run_action, "quick run from a script thread")
    expect_dialogs(ctx, t, [("info", "Analysis complete.")], "quick run")
    t.check(as_text(app.last_run_dir) != run_dir, "quick run writes its own run folder")
    compare_run(t, "quick run", window_tables(as_text(app.last_run_dir)), engine)


def step_prepared_cache(ctx, t):
    app = ctx.app
    use_as_source(ctx, ctx.imp1)
    click(ctx, "1. Prepare enabled channels")
    first = [app.prepared_images[channel].getID() for channel in [1, 2, 3]]
    set_widgets(app.channel_panels[0], {"gaussian_sigma": 1})
    t.check(not app.prepared_cache_valid(app.collect_channel_configs(), 1), "changing a preprocessing setting makes the prepared images stale")
    click(ctx, "2. Preview all masks")
    second = [app.prepared_images[channel].getID() for channel in [1, 2, 3]]
    t.check(first[0] != second[0], "stale prepared images are rebuilt before a mask preview")
    t.equal([image_id for image_id in first if WindowManager.getImage(image_id) is not None], [], "the stale prepared images are closed")
    set_widgets(app.channel_panels[0], {"gaussian_sigma": 0})

    set_widgets(app.channel_panels[2], {"projection": "Current Z"})
    on_edt(lambda: ctx.imp1.setPosition(1, 1, 1), "moving to Z1")
    click(ctx, "1. Prepare enabled channels")
    t.equal(app.prepared_images[3].getStatistics().max, 300.0, "'Current Z' prepares the displayed slice")
    on_edt(lambda: ctx.imp1.setPosition(1, 2, 1), "moving to Z2")
    t.check(not app.prepared_cache_valid(app.collect_channel_configs(), 1), "moving the Z slider makes a 'Current Z' prepared image stale")
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("info", "Analysis complete.")], "analysis after moving the Z slider")
    summary = read_csv(os.path.join(as_text(app.last_run_dir), "roi_summary.csv"))
    t.close(row_named(summary[1], "rect").get("C3_Channel_3_Mean"), 400.0, "the analysis measures the slice now displayed, not a stale prepared image")
    set_widgets(app.channel_panels[2], {"projection": "Max Z"})
    on_edt(lambda: ctx.imp1.setPosition(1, 1, 1), "moving back to Z1")

    click(ctx, "1. Prepare enabled channels")
    prepared = app.prepared_images[2]
    click(ctx, "Use active image", current_image=prepared)
    records = expect_dialogs(ctx, t, [("warning", "was made by this application")], "choosing a prepared image as the source")
    t.check(app.source_imp is ctx.imp1, "a prepared image is not taken as the source")
    t.check(WindowManager.getImage(prepared.getID()) is not None, "the prepared image stays open")
    expect_dialogs(ctx, t, [], "prepared cache")


def step_display_minimum(ctx, t):
    app = ctx.app
    set_widgets(app.channel_panels[1], {"threshold_mode": "Display minimum"})
    click(ctx, "Quick run")
    expect_dialogs(ctx, t, [("error", "Quick run cannot")], "quick run with a display-minimum threshold")
    click(ctx, "1. Prepare enabled channels")
    prepared = app.prepared_images[2]
    on_edt(lambda: prepared.setDisplayRange(1200, 3000), "setting the display range")
    click(ctx, "Preview mask", container=app.channel_panels[1])
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("info", "Analysis complete.")], "analysis with a display-minimum threshold")
    run_dir = as_text(app.last_run_dir)
    thresholds = read_csv(os.path.join(run_dir, "thresholds.csv"))
    c2 = [row for row in thresholds[1] if row.get("Channel") == "2"]
    t.equal(len(c2), 1, "one threshold row for C2")
    if len(c2) == 1:
        t.equal(c2[0].get("Threshold_Mode"), "Display minimum", "the threshold mode is recorded")
        t.close(c2[0].get("Threshold_Low"), 1200.0, "the display minimum is the threshold")
    summary = read_csv(os.path.join(run_dir, "roi_summary.csv"))
    t.close(row_named(summary[1], "rect").get("C2_MAP2_Positive_Area"), MARKER_INNER[2] * MARKER_INNER[3] * PIXEL * PIXEL,
            "only the brighter inner rectangle is above the display minimum")

    set_widgets(app.channel_panels[0], {"gaussian_sigma": 1})
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("error", "Display-minimum thresholding requires the current prepared images")],
                   "display minimum with stale prepared images")
    # Validation, sensitivity and the capture button must not rebuild the
    # prepared image either: its display range would be a default one.
    on_edt(lambda: app.qc_channel_combo.setSelectedIndex(1), "choosing the validated channel")
    click(ctx, "Run sensitivity analysis")
    expect_dialogs(ctx, t, [("error", "Display-minimum thresholding requires the current prepared images")],
                   "sensitivity on a display-minimum channel with stale prepared images")
    click(ctx, "Run validation")
    expect_dialogs(ctx, t, [("error", "Display-minimum thresholding requires the current prepared images")],
                   "validation of a display-minimum channel with stale prepared images")
    click(ctx, "Capture display minimum", container=app.channel_panels[1])
    expect_dialogs(ctx, t, [("error", "There is no current prepared image to capture from")],
                   "capturing with stale prepared images")
    t.equal(app.channel_panels[1].get_config()["threshold_mode"], "Display minimum",
            "a refused capture leaves the threshold mode alone")
    set_widgets(app.channel_panels[0], {"gaussian_sigma": 0})

    click(ctx, "1. Prepare enabled channels")
    prepared = app.prepared_images[2]
    on_edt(lambda: prepared.setDisplayRange(850, 3000), "setting the display range")
    click(ctx, "Brightness / Contrast", container=app.channel_panels[1])
    adjusters = [window for window in Window.getWindows() if window.isShowing() and "ContrastAdjuster" in as_text(type(window))]
    t.equal(len(adjusters), 1, "the Brightness/Contrast window opens")
    for window in adjusters:
        on_edt(window.close, "closing Brightness/Contrast")
    t.equal(prepared.getDisplayRangeMin(), 850.0, "opening Brightness/Contrast keeps the display range")
    click(ctx, "Capture display minimum", container=app.channel_panels[1])
    expect_dialogs(ctx, t, [], "capturing the display minimum")
    config = app.channel_panels[1].get_config()
    t.equal(config["threshold_mode"], "Manual", "capturing switches the channel to a Manual threshold")
    t.equal(config["manual_low"], 850.0, "the captured lower threshold")
    t.equal(config["manual_high"], 1500.0, "the captured upper threshold is the image maximum")
    set_widgets(app.channel_panels[1], C2_VALUES)


def read_text(path):
    handle = codecs.open(path, "r", "utf-8")
    try:
        return handle.read()
    finally:
        handle.close()


def step_edited_source(ctx, t):
    """Pixels edited by hand after opening. One plane of the source set to 0
    inside the parent ROI drops out of a maximum projection without any
    error: one Backspace key press (Edit > Clear) on C1 Z1 reproduced the
    2.2 run of 2026-11-18 exactly. The window now asks first."""
    app = ctx.app
    ctx.imp4 = make_image("smoke edited.tif")
    use_as_source(ctx, ctx.imp4)
    t.check(not ctx.imp4.changes, "a new image has no unsaved changes")

    def edit_source():
        ctx.imp4.setPosition(3, 2, 1)
        ctx.imp4.setRoi(Roi(RECT[0], RECT[1], RECT[2], RECT[3]))
        IJ.run(ctx.imp4, "Set...", "value=0 slice")
        ctx.imp4.deleteRoi()
        ctx.imp4.setPosition(1, 1, 1)

    on_edt(edit_source, "setting C3 Z2 to 0 inside the rectangle")
    t.check(ctx.imp4.changes, "editing a plane marks the source as changed")
    before = as_text(app.last_run_dir)
    ctx.watcher.answer_next("no")
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("warning", "has been changed since it was opened")], "analysis of an edited source, declined")
    t.equal(as_text(app.last_run_dir), before, "declining writes no run folder")
    t.check("Analysis cancelled: the source image has unsaved changes." in as_text(app.status_area.getText()), "declining is logged")

    click(ctx, "1. Prepare enabled channels")
    expect_dialogs(ctx, t, [("warning", "Backspace or Delete key")], "preparing an edited source, accepted")
    t.equal(sorted(app.prepared_images.keys()), [1, 2, 3], "after accepting, the channels are prepared")
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("info", "Analysis complete.")], "the question is asked once per source image")
    run_dir = as_text(app.last_run_dir)
    summary = read_csv(os.path.join(run_dir, "roi_summary.csv"))
    t.close(row_named(summary[1], "rect").get("C3_Channel_3_Mean"), 300.0,
            "the edited plane no longer contributes to the maximum projection (400 before the edit)")
    t.check("unsaved changes in memory" in read_text(os.path.join(run_dir, "run_log.txt")),
            "run_log.txt records that the source was edited")

    # A prepared image edited by hand is not measured: the channel is prepared again.
    configs = app.collect_channel_configs()
    t.check(app.prepared_cache_valid(configs, 1), "the run leaves the prepared images current")
    prepared = app.prepared_images[3]
    t.check(not prepared.changes, "the window's own work does not mark a prepared image as edited")
    on_edt(lambda: IJ.run(prepared, "Add...", "value=1000"), "editing a prepared image by hand")
    t.check(prepared.changes, "editing marks the prepared image as changed")
    t.check(not app.prepared_cache_valid(configs, 1), "a prepared image edited by hand is not current")
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("info", "Analysis complete.")], "analysis after editing a prepared image")
    summary = read_csv(os.path.join(as_text(app.last_run_dir), "roi_summary.csv"))
    t.close(row_named(summary[1], "rect").get("C3_Channel_3_Mean"), 300.0, "the edited prepared image is not measured")
    t.check("was edited by hand" in as_text(app.status_area.getText()), "the user is told why the channels were prepared again")

    set_widgets(app.channel_panels[1], {"threshold_mode": "Display minimum"})
    click(ctx, "1. Prepare enabled channels")
    marker = app.prepared_images[2]
    on_edt(lambda: IJ.run(marker, "Add...", "value=10"), "editing the prepared marker image")
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("error", "do not press 'Apply'")], "display minimum on an edited prepared image")
    set_widgets(app.channel_panels[1], C2_VALUES)

    def close_edited():
        ctx.imp4.changes = False
        ctx.imp4.close()

    on_edt(close_edited, "closing the edited image")
    expect_dialogs(ctx, t, [], "edited source")


def step_module_settings(ctx, t):
    """Settings of the analysis modules (shown on the Analyses tab): a
    profile's values are kept, saved again, recorded, and used by a run."""
    app = ctx.app
    E = ctx.E
    if "depth_enabled" not in E.default_global_options():
        t.check(False, "the depth-profile module (cellquantdepth.py) is loaded")
        return
    base = os.path.join(ctx.tmp, "modules base.properties")
    on_edt(lambda: app.save_settings_to_path(base), "saving a profile")
    props = E.read_properties_file(base)
    t.equal(as_text(props.getProperty("depth.enabled")), "false", "a profile saved by the window holds the module settings")
    props.setProperty("depth.enabled", "true")
    props.setProperty("depth.bin_width", "10")
    path = os.path.join(ctx.tmp, "modules depth.properties")
    stream = FileOutputStream(path)
    try:
        props.store(stream, "depth on")
    finally:
        stream.close()

    on_edt(lambda: app.load_profile_from_path(path), "loading a profile that switches the depth profile on")
    expect_dialogs(ctx, t, [], "loading a profile with module settings")
    options = app.collect_global_options()
    t.equal(options.get("depth_enabled"), True, "the module setting is taken from the profile")
    t.equal(as_text(options.get("depth_bin_width")), "10", "a text module setting is taken from the profile")
    log = as_text(app.status_area.getText())
    t.check("depth.enabled=true" in log and "depth.bin_width=10" in log, "the status log names the module settings the profile sets")
    t.check(app.module_widgets["depth_enabled"][1].isSelected(), "the Analyses tab shows the module switch from the profile")
    t.equal(as_text(app.module_widgets["depth_bin_width"][1].getText()), "10", "the Analyses tab shows the module value from the profile")
    again = os.path.join(ctx.tmp, "modules again.properties")
    on_edt(lambda: app.save_settings_to_path(again), "saving the profile again")
    saved = E.load_profile_settings(again)[0]
    t.equal(saved["depth_enabled"], True, "saving again keeps the module setting")
    t.equal(as_text(saved["depth_bin_width"]), "10", "saving again keeps the text module setting")

    use_as_source(ctx, ctx.imp1)
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("info", "Analysis complete.")], "analysis with the depth profile on")
    run_dir = as_text(app.last_run_dir)
    t.check(os.path.isfile(os.path.join(run_dir, "depth_profile.csv")), "the run writes depth_profile.csv")
    t.equal(E.load_profile_settings(os.path.join(run_dir, "settings.properties"))[0]["depth_enabled"], True,
            "the run record holds the module setting")
    t.check("Analysis-module options in use: depth.enabled=true" in read_text(os.path.join(run_dir, "run_log.txt")),
            "run_log.txt names the module settings in use")
    window = window_tables(run_dir)
    reference_image = make_image("smoke one.tif")
    engine = engine_reference(ctx, reference_image, "", expected_configs(ctx), [triangle_roi(), rect_roi()], [1],
                              {"depth_enabled": True, "depth_bin_width": "10"})
    reference_image.close()
    compare_run(t, "depth profile on", window, engine)
    t.check("ROI_Max_Depth" in window["summary"][0], "the depth columns are appended to roi_summary.csv")

    # A profile without module keys (written before the modules existed): defaults.
    legacy = os.path.join(ctx.tmp, "modules legacy.properties")
    props = E.read_properties_file(path)
    for key in list(props.stringPropertyNames()):
        if as_text(key).startswith("depth."):
            props.remove(key)
    stream = FileOutputStream(legacy)
    try:
        props.store(stream, "no module keys")
    finally:
        stream.close()
    on_edt(lambda: app.load_profile_from_path(path), "loading the depth profile again")
    on_edt(lambda: app.load_profile_from_path(legacy), "loading a profile without module keys")
    options = app.collect_global_options()
    t.equal(options.get("depth_enabled"), False, "a module switch missing from a profile is off")
    t.equal(as_text(options.get("depth_bin_width")), "25", "a module value missing from a profile gets its default")
    on_edt(lambda: app.load_profile_from_path(base), "loading the first profile again")
    expect_dialogs(ctx, t, [], "module settings")


def step_all_frames(ctx, t):
    app = ctx.app
    per_frame = [DISKS[:5], DISKS[:9]]
    ctx.imp3 = make_image("smoke frames.tif", frames=2, disks_per_frame=per_frame)
    use_as_source(ctx, ctx.imp3)
    on_edt(lambda: ctx.imp3.setPosition(1, 1, 2), "moving to frame 2")
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("info", "Analysis complete.")], "current frame only")
    summary = read_csv(os.path.join(as_text(app.last_run_dir), "roi_summary.csv"))
    t.equal([row.get("Frame") for row in summary[1]], ["2", "2"], "'Current frame only' analyzes the displayed frame")
    t.equal(row_named(summary[1], "rect").get("C1_DAPI_Object_Count"), "9", "frame 2 holds nine disks")

    on_edt(lambda: app.frame_mode_combo.setSelectedItem("All frames"), "choosing all frames")
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("info", "Analysis complete.")], "all frames")
    window = window_tables(as_text(app.last_run_dir))
    reference_image = make_image("smoke frames.tif", frames=2, disks_per_frame=per_frame)
    engine = engine_reference(ctx, reference_image, "", expected_configs(ctx), [triangle_roi(), rect_roi()], [1, 2],
                              {"frame_mode": "All frames"})
    reference_image.close()
    compare_run(t, "all frames", window, engine)
    counts = [row.get("C1_DAPI_Object_Count") for row in window["summary"][1] if row.get("Parent_ROI_Name") == "rect"]
    t.equal(counts, ["5", "9"], "one row per frame with that frame's count")
    on_edt(lambda: app.frame_mode_combo.setSelectedItem("Current frame only"), "back to the current frame")


def step_profile(ctx, t):
    app = ctx.app
    E = ctx.E
    path = os.path.join(ctx.tmp, "profile " + MICRO + ".properties")
    on_edt(lambda: app.spatial_unit_field.setText(MICRO + "m"), "typing a unit")
    on_edt(lambda: app.save_settings_to_path(path), "saving a profile")
    t.check(os.path.isfile(path), "the profile file is written")
    before_options = app.collect_global_options()
    before_configs = app.collect_all_channel_configs()

    set_widgets(app.channel_panels[0], {"manual_low": 1234, "morph_count": 3, "fill_holes": True})
    set_widgets(app.channel_panels[1], {"enabled": False, "threshold_mode": "Auto"})
    set_widgets(app.channel_panels[4], {"name": "Changed"})

    def change_globals():
        app.frame_mode_combo.setSelectedItem("All frames")
        app.roi_source_combo.setSelectedItem("Whole image")
        app.save_masks_cb.setSelected(False)
        app.batch_recursive_cb.setSelected(True)
        app.spatial_unit_field.setText("inch")

    on_edt(change_globals, "changing global settings")
    t.check(app.collect_global_options() != before_options, "the settings really changed")
    click(ctx, "1. Prepare enabled channels")
    on_edt(lambda: app.load_profile_from_path(path), "loading a profile")
    expect_dialogs(ctx, t, [], "loading a profile")
    after_options = app.collect_global_options()
    for key in sorted(before_options.keys()):
        t.equal(after_options.get(key), before_options[key], "global setting '" + key + "' is restored from the profile")
    after_configs = app.collect_all_channel_configs()
    for index in range(len(before_configs)):
        for key in sorted(before_configs[index].keys()):
            t.equal(after_configs[index][key], before_configs[index][key],
                    "C" + str(index + 1) + " setting '" + key + "' is restored from the profile")
    t.check(not app.prepared_cache_valid(app.collect_channel_configs(), app.source_imp.getFrame()),
            "loading a profile makes the prepared images stale")
    loaded_options, enabled, everything = E.load_profile_settings(path)
    window_configs = app.collect_channel_configs()
    t.equal(len(enabled), len(window_configs), "the engine reads the same enabled channels from the profile")
    for index in range(min(len(enabled), len(window_configs))):
        t.equal(enabled[index], window_configs[index], "the engine reads C" + str(window_configs[index]["channel"]) + " as the window shows it")

    # A profile written by 2.2.x: no neighbor count, no batch keys.
    legacy_path = os.path.join(ctx.tmp, "legacy.properties")
    props = E.read_properties_file(path)
    for key in list(props.stringPropertyNames()):
        if as_text(key).startswith("batch.") or as_text(key).endswith(".morph_count"):
            props.remove(key)
    stream = FileOutputStream(legacy_path)
    try:
        props.store(stream, "legacy")
    finally:
        stream.close()
    set_widgets(app.channel_panels[0], {"morph_count": 5})
    on_edt(lambda: app.batch_recursive_cb.setSelected(True), "changing a batch option")
    on_edt(lambda: app.load_profile_from_path(legacy_path), "loading a 2.2 profile")
    expect_dialogs(ctx, t, [], "loading a 2.2 profile")
    t.equal(app.channel_panels[0].get_config()["morph_count"], 1, "a key missing from an old profile gets its default")
    t.equal(app.batch_recursive_cb.isSelected(), True, "batch options absent from an old profile are left alone")
    t.equal(app.channel_panels[0].get_config()["manual_low"], 1000.0, "the old profile's own values are applied")
    on_edt(lambda: app.batch_recursive_cb.setSelected(False), "restoring a batch option")

    example = os.path.join(ctx.repo, "examples", "Example_4Channel_SOX2_Casp3_MAP2_DAPI.properties")
    if os.path.isfile(example):
        on_edt(lambda: app.load_profile_from_path(example), "loading the example profile")
        expect_dialogs(ctx, t, [], "loading the example profile")
        example_options, example_enabled, example_all = E.load_profile_settings(example)
        shown = [config for config in app.collect_all_channel_configs() if config["enabled"]]
        t.equal([config["name"] for config in shown], [config["name"] for config in example_enabled],
                "the example profile enables the same channels in the window as in the engine")

    on_edt(lambda: app.load_profile_from_path(os.path.join(ctx.tmp, "missing.properties")), "loading a missing profile")
    expect_dialogs(ctx, t, [("error", "Could not load settings profile")], "loading a missing profile")
    on_edt(lambda: app.load_profile_from_path(path), "loading the saved profile again")
    on_edt(lambda: app.spatial_unit_field.setText("um"), "restoring the unit")
    expect_dialogs(ctx, t, [], "profile")


def step_validation(ctx, t):
    app = ctx.app
    E = ctx.E
    use_as_source(ctx, ctx.imp1)
    points = [DISKS[i] for i in range(13)]
    points.append((DISKS[0][0] + 2, DISKS[0][1] + 1))
    points.append((20, 110))
    marks = PointRoi(points[0][0], points[0][1])
    for point in points[1:]:
        marks.addPoint(point[0], point[1])

    def setup():
        ctx.imp1.setRoi(marks)
        app.qc_channel_combo.setSelectedIndex(0)
        app.annotation_source_combo.setSelectedItem(E.ANNOTATION_POINTS_IMAGE)
        app.roi_source_combo.setSelectedItem("ROI Manager")

    on_edt(setup, "marking points")
    state = press_and_watch(ctx, "Run validation")
    t.check(state.get("active") and as_text(state["thread"].getName()).startswith("cellquant-validation"),
            "the validation runs on a worker thread")
    expect_dialogs(ctx, t, [("info", "Validation complete.")], "validation")
    run_dir = as_text(app.last_run_dir)
    t.check("_CellQuant_Validation_" in os.path.basename(run_dir), "a validation folder is written")
    for name in ["validation_summary.csv", "validation_objects.csv", "validation_points.csv", "thresholds.csv", "settings.properties", "Parent_ROIs_Input.zip"]:
        t.check(os.path.isfile(os.path.join(run_dir, name)), name + " is written")
    summary = read_csv(os.path.join(run_dir, "validation_summary.csv"))
    rect = row_named(summary[1], "rect")
    t.equal(
        [rect.get(key) for key in ["Manual_Count", "Detected_Count", "True_Positive_Objects", "False_Positive_Objects", "Missed_Annotations", "Merged_Annotations"]],
        ["15", "15", "13", "2", "1", "1"], "validation counts in the rectangle"
    )
    t.close(rect.get("Precision"), 13.0 / 15.0, "precision")
    t.close(rect.get("Recall"), 13.0 / 15.0, "recall")
    t.equal([row.get("Parent_ROI_Name") for row in summary[1]], ["triangle", "rect", "ALL_ROIS_POOLED"], "one row per parent ROI and a pooled row")

    reference_image = make_image("smoke one.tif")
    engine = E.AnalysisEngine()
    config = expected_configs(ctx)[0]
    prepared = engine.prepare_single_channel(reference_image, config, 1)
    reference_dir = E.unique_dir(os.path.join(ctx.tmp, "reference"), "validation")
    E.validate_objects_against_points(
        engine, prepared, config, E.validate_parent_rois([triangle_roi(), rect_roi()], WIDTH, HEIGHT),
        points, reference_dir, "smoke one.tif", 1
    )
    prepared.close()
    reference_image.close()
    for name in ["validation_summary.csv", "validation_objects.csv", "validation_points.csv"]:
        compare_tables(t, name, read_csv(os.path.join(run_dir, name)), read_csv(os.path.join(reference_dir, name)), ["Run_ID"])
    t.check(WindowManager.getWindow("CellQuant Object Validation") is not None, "the validation table is shown")
    t.check(ctx.imp1.getRoi() is not None and ctx.imp1.getRoi().getType() == Roi.POINT, "the point marks stay on the image")

    on_edt(lambda: app.roi_source_combo.setSelectedItem("Current image ROI"), "choosing the image ROI")
    click(ctx, "Run validation")
    expect_dialogs(ctx, t, [("error", "cannot be both the parent ROI and the point marks")], "marks and parent ROI in one selection")
    on_edt(lambda: app.roi_source_combo.setSelectedItem("ROI Manager"), "back to ROI Manager")
    on_edt(ctx.imp1.deleteRoi, "removing the marks")
    click(ctx, "Run validation")
    expect_dialogs(ctx, t, [("error", "No point marks were found")], "validation without marks")

    marks_path = os.path.join(ctx.tmp, "marks.roi")
    E.save_roi_array([marks], marks_path)

    def from_file():
        app.annotation_source_combo.setSelectedItem(E.ANNOTATION_POINTS_FILE)
        app.annotation_file_field.setText(marks_path)

    on_edt(from_file, "choosing a point ROI file")
    script_call(ctx, app.run_validation_action, "validation from a script thread")
    expect_dialogs(ctx, t, [("info", "Validation complete.")], "validation from a point ROI file")
    compare_tables(t, "validation from a file", read_csv(os.path.join(as_text(app.last_run_dir), "validation_summary.csv")),
                   read_csv(os.path.join(reference_dir, "validation_summary.csv")), ["Run_ID"])

    mask_path = os.path.join(ctx.tmp, "reference_mask.tif")
    mask = ImagePlus("reference_mask", marker_plane())
    mask.getProcessor().setThreshold(1000, 65535, 2)
    mask = ImagePlus("reference_mask", mask.getProcessor().createMask())
    FileSaver(mask).saveAsTiff(mask_path)
    mask.close()

    def from_mask():
        app.annotation_source_combo.setSelectedItem(E.ANNOTATION_MASK_FILE)
        app.annotation_file_field.setText(mask_path)
        app.qc_channel_combo.setSelectedIndex(1)

    on_edt(from_mask, "choosing a reference mask")
    click(ctx, "Run validation")
    expect_dialogs(ctx, t, [("info", "Validation complete.")], "validation against a reference mask")
    mask_summary = read_csv(os.path.join(as_text(app.last_run_dir), "mask_validation_summary.csv"))
    rect = row_named(mask_summary[1], "rect")
    t.equal(rect.get("Test_Positive_Pixels"), str(MARKER_OUTER[2] * MARKER_OUTER[3]), "pixels of the window's mask")
    t.equal(rect.get("Reference_Positive_Pixels"), str(MARKER_INNER[2] * MARKER_INNER[3]), "pixels of the reference mask")
    t.close(rect.get("Recall"), 1.0, "the reference lies inside the mask")
    on_edt(lambda: app.annotation_source_combo.setSelectedItem(E.ANNOTATION_POINTS_IMAGE), "restoring the annotation source")


def step_sensitivity(ctx, t):
    app = ctx.app
    E = ctx.E

    def setup():
        app.qc_channel_combo.setSelectedIndex(1)
        app.sensitivity_parameter_combo.setSelectedItem(E.SENS_THRESHOLD)
        app.sensitivity_factors_field.setText("0.9, 1.5")

    on_edt(setup, "choosing the sweep")
    state = press_and_watch(ctx, "Run sensitivity analysis")
    t.check(state.get("active") and as_text(state["thread"].getName()).startswith("cellquant-sensitivity"),
            "the sensitivity analysis runs on a worker thread")
    expect_dialogs(ctx, t, [("info", "Sensitivity analysis complete.")], "sensitivity")
    run_dir = as_text(app.last_run_dir)
    t.check("_CellQuant_Sensitivity_" in os.path.basename(run_dir), "a sensitivity folder is written")
    table = read_csv(os.path.join(run_dir, "sensitivity.csv"))
    rect = [row for row in table[1] if row.get("Parent_ROI_Name") == "rect"]
    t.equal([row.get("Factor") for row in rect], ["0.9", "1.0", "1.5"], "one row per factor, the baseline included")
    outer = MARKER_OUTER[2] * MARKER_OUTER[3] * PIXEL * PIXEL
    inner = MARKER_INNER[2] * MARKER_INNER[3] * PIXEL * PIXEL
    if len(rect) == 3:
        t.close(rect[0].get("Positive_Area"), outer, "a lower threshold keeps the whole marker")
        t.close(rect[1].get("Positive_Area"), outer, "the configured threshold")
        t.close(rect[2].get("Positive_Area"), inner, "a higher threshold keeps the inner rectangle only")
        t.close(rect[2].get("Positive_Area_Change_Percent"), (inner - outer) / outer * 100.0, "percent change from the configured value")
    t.check(WindowManager.getWindow("CellQuant Sensitivity") is not None, "the sensitivity table is shown")

    on_edt(lambda: app.sensitivity_factors_field.setText("0.9, x"), "typing a bad factor")
    click(ctx, "Run sensitivity analysis")
    expect_dialogs(ctx, t, [("error", "Sensitivity analysis failed")], "a non-numeric factor")
    on_edt(lambda: app.sensitivity_factors_field.setText(E.DEFAULT_SENSITIVITY_FACTORS), "restoring the factors")


def step_open_from_file(ctx, t):
    app = ctx.app
    folder = os.path.join(ctx.tmp, "files " + MICRO)
    os.makedirs(folder)
    path = os.path.join(folder, "smoke file.tif")
    image = make_image("smoke file.tif", disks=DISKS[:4])
    save_image(image, path)
    image.close()

    opened = []

    def remember():
        if app.source_imp is not None:
            opened.append(app.source_imp)

    script_call(ctx, lambda: app.open_image_path(path), "opening a file from a script thread")
    expect_dialogs(ctx, t, [], "opening a file")
    remember()
    t.check(app.source_imp is not None and as_text(app.source_path) == os.path.abspath(path), "the opened file is the source")
    t.equal(as_text(app.image_path_field.getText()), os.path.abspath(path), "the path field shows the file")
    if app.source_imp is not None:
        t.equal((app.source_imp.getNChannels(), app.source_imp.getNSlices()), (3, 2), "channels and slices are read from the file")
        t.check(WindowManager.getImage(app.source_imp.getID()) is not None, "the opened image is shown")

    click(ctx, "Open selected image")
    expect_dialogs(ctx, t, [], "'Open selected image'")
    remember()
    click(ctx, "Open via Bio-Formats")
    expect_dialogs(ctx, t, [], "'Open via Bio-Formats'")
    remember()
    if app.source_imp is not None:
        t.equal((app.source_imp.getNChannels(), app.source_imp.getNSlices()), (3, 2), "Bio-Formats reads the same dimensions")
        click(ctx, "3. Run analysis")
        expect_dialogs(ctx, t, [("info", "Analysis complete.")], "analysis of a file opened with Bio-Formats")
        summary = read_csv(os.path.join(as_text(app.last_run_dir), "roi_summary.csv"))
        t.equal(row_named(summary[1], "rect").get("C1_DAPI_Object_Count"), "4", "four disks in the file")
        t.equal(row_named(summary[1], "rect").get("Source_Path"), os.path.abspath(path), "the source path is recorded, non-ASCII included")
        t.check(as_text(os.path.basename(as_text(app.last_run_dir))).startswith("smoke_file_CellQuant_"), "the run folder is named after the file")

    missing = os.path.join(folder, "missing " + MICRO + ".tif")
    on_edt(lambda: app.open_image_path(missing), "opening a missing file")
    expect_dialogs(ctx, t, [("error", "Image file does not exist: " + os.path.abspath(missing))], "a missing file")
    bad = os.path.join(folder, "not an image.tif")
    handle = open(bad, "w")
    handle.write("not an image")
    handle.close()
    script_call(ctx, lambda: app.open_image_path(bad), "opening an unreadable file")
    records = ctx.watcher.take()
    t.check(len(records) >= 1 and records[-1]["level"] == "error" and "Could not open image" in records[-1]["text"],
            "an unreadable file is reported: " + repr([record["text"] for record in records]))
    for imp in opened:
        if imp is not app.source_imp:
            imp.changes = False
            imp.close()
    ctx.file_path = path


def write_batch_folder(ctx, name, disk_counts):
    folder = os.path.join(ctx.tmp, name)
    os.makedirs(folder)
    for index in range(len(disk_counts)):
        title = "field_" + str(index + 1).zfill(2) + ".tif"
        image = make_image(title, disks=DISKS[:disk_counts[index]])
        save_image(image, os.path.join(folder, title))
        image.close()
    return folder


def wait_for_batch(ctx, timeout=120.0):
    wait_until(lambda: not ctx.app.batch_running, timeout, "the batch to finish")
    settle(ctx)


def step_batch(ctx, t):
    app = ctx.app
    E = ctx.E
    counts = [15, 7, 3]
    images_before = WindowManager.getImageCount()
    folder = write_batch_folder(ctx, "batch_in", counts)
    handle = open(os.path.join(folder, "notes.txt"), "w")
    handle.write("not an image")
    handle.close()

    def setup():
        app.batch_input_field.setText(folder)
        app.batch_extensions_field.setText("tif")
        app.batch_roi_combo.setSelectedItem(E.BATCH_ROI_WHOLE)
        app.output_dir_field.setText(ctx.out_dir)
        app.main_tabs.setSelectedIndex(TAB_TITLES.index("Batch"))

    on_edt(setup, "filling the Batch tab")
    click(ctx, "Scan folder")
    expect_dialogs(ctx, t, [], "scanning the folder")
    listing = as_text(app.batch_list_area.getText())
    t.check("3 image file(s) found" in listing and "field_02.tif" in listing and "notes.txt" not in listing, "the scan lists the three images")

    set_widgets(app.channel_panels[1], {"threshold_mode": "Display minimum"})
    click(ctx, "Run batch")
    expect_dialogs(ctx, t, [("error", "Display-minimum")], "a batch with a display-minimum threshold")
    t.check(not app.batch_running, "a refused batch does not start")
    set_widgets(app.channel_panels[1], C2_VALUES)

    ctx.watcher.answer_next("no")
    click(ctx, "Run batch")
    records = expect_dialogs(ctx, t, [("question", "Analyze 3 image file(s)")], "declining the batch")
    t.check(not app.batch_running, "declining the confirmation starts nothing")

    state = {}
    buttons = [button for button in find_components(app, JButton) if as_text(button.getText()) == "Run batch"]

    def start():
        buttons[0].doClick(0)
        state["running"] = app.batch_running
        state["run_enabled"] = app.batch_run_button.isEnabled()
        state["cancel_enabled"] = app.batch_cancel_button.isEnabled()
        state["actions"] = [button.isEnabled() for button in app.analysis_action_buttons]

    on_edt(start, "starting the batch")
    t.equal((state.get("running"), state.get("run_enabled"), state.get("cancel_enabled")), (True, False, True),
            "while a batch runs, Run is disabled and Cancel is enabled")
    on_edt(nothing, "the event thread during a batch", 5.0)
    wait_for_batch(ctx)
    expect_dialogs(ctx, t, [("question", "Analyze 3 image file(s)"), ("info", "Batch complete.")], "batch")
    t.equal((app.batch_run_button.isEnabled(), app.batch_cancel_button.isEnabled()), (True, False), "the buttons are restored after the batch")
    t.equal(as_text(app.batch_progress.getString()), "Finished", "the progress bar reports the end")
    t.equal((app.batch_progress.getValue(), app.batch_progress.getMaximum()), (3, 3), "the progress bar is full")
    t.check(len(state.get("actions", [])) > 0 and True not in state["actions"], "the Analyses buttons are disabled while a batch runs")
    t.check(False not in [button.isEnabled() for button in app.analysis_action_buttons], "the Analyses buttons are enabled again after the batch")
    batch_dir = as_text(app.last_run_dir)
    ctx.batch_dir = batch_dir
    t.check(os.path.basename(batch_dir).startswith("CellQuant_Batch_") and batch_dir.startswith(ctx.out_dir), "the batch folder is in the output folder")
    report = read_csv(os.path.join(batch_dir, "batch_report.csv"))
    t.equal([row.get("Status") for row in report[1]], ["OK", "OK", "OK"], "every image is reported OK")
    summary = read_csv(os.path.join(batch_dir, "batch_roi_summary.csv"))
    t.equal([row.get("C1_DAPI_Object_Count") for row in summary[1]], [str(count) for count in counts], "object counts of the three images")
    t.equal([row.get("Source_Image") for row in summary[1]], ["field_01.tif", "field_02.tif", "field_03.tif"], "images are named after their files")

    reference_image = make_image("field_02.tif", disks=DISKS[:7])
    whole = E.whole_image_roi(reference_image)
    engine = engine_reference(ctx, reference_image, os.path.join(folder, "field_02.tif"), expected_configs(ctx), [whole], [1],
                              {"batch_extensions": "tif"})
    reference_image.close()
    run_folder = os.path.join(batch_dir, report[1][1].get("Run_Folder"))
    compare_run(t, "batch image 2", window_tables(run_folder), engine)
    log = as_text(app.status_area.getText())
    t.check("Batch finished" in log, "batch messages reach the status log")
    t.equal(WindowManager.getImageCount(), images_before, "the batch leaves no image open")


def step_batch_cancel(ctx, t):
    app = ctx.app
    E = ctx.E
    total = 16
    folder = write_batch_folder(ctx, "batch_cancel", [5] * total)

    def setup():
        app.batch_input_field.setText(folder)
        app.batch_roi_combo.setSelectedItem(E.BATCH_ROI_MANAGER)

    on_edt(setup, "filling the Batch tab")
    buttons = [button for button in find_components(app, JButton) if as_text(button.getText()) == "Run batch"]
    on_edt(lambda: buttons[0].doClick(0), "starting the batch")
    t.check(app.batch_running, "the batch is running")
    wait_until(lambda: app.batch_progress.getValue() >= 1 or not app.batch_running, 60.0, "the first batch image")
    progress_text = as_text(app.batch_progress.getString())
    click(ctx, "3. Run analysis")
    click(ctx, "Cancel batch")
    wait_for_batch(ctx)
    records = ctx.watcher.take()
    texts = [record["text"] for record in records]
    t.check(len([text for text in texts if "A batch is running" in text]) == 1, "other actions are refused while a batch runs: " + repr(texts))
    t.check(len([text for text in texts if text.startswith("Batch cancelled.")]) == 1, "the batch reports that it was cancelled: " + repr(texts))
    t.check(" / " + str(total) in progress_text, "the progress bar counted images while running: " + repr(progress_text))
    t.equal(as_text(app.batch_progress.getString()), "Cancelled", "the progress bar reports the cancellation")
    t.equal((app.batch_run_button.isEnabled(), app.batch_cancel_button.isEnabled(), app.batch_running), (True, False, False),
            "the window is ready again after a cancelled batch")
    batch_dir = as_text(app.last_run_dir)
    report = read_csv(os.path.join(batch_dir, "batch_report.csv"))
    finished = [row for row in report[1] if row.get("Status") == "OK"]
    t.check(1 <= len(finished) < total, "finished images are kept and the rest are not analyzed (" + str(len(finished)) + " of " + str(total) + ")")
    summary = read_csv(os.path.join(batch_dir, "batch_roi_summary.csv"))
    t.equal(len(summary[1]), 2 * len(finished), "the combined table holds two parent ROIs per finished image")
    t.check(RoiManager.getInstance2() is ctx.manager and ctx.manager.getCount() == 2, "ROI Manager is untouched by the batch")
    log = as_text(app.status_area.getText())
    t.check("Batch cancel requested" in log and "Batch cancelled by user." in log, "the cancellation is logged")
    on_edt(lambda: app.main_tabs.setSelectedIndex(0), "back to the Session tab")


# -----------------------------------------------------------------------------
# Background jobs and the Analyses tab
# -----------------------------------------------------------------------------

ANALYSIS_ACTIONS = [
    ("detect_tissue_outline", "Detect tissue outline into ROI Manager"),
    ("run_tuning_folder", "Run tuning..."),
    ("run_count_correction", "Run count correction..."),
    ("run_robustness_from_options", "Run robustness...")
]


class TimerProbe(ActionListener):
    """Counts the ticks of a Swing timer, which fire only while the event
    thread is free, that come while the worker thread is alive."""
    def __init__(self):
        self.worker = None
        self.busy = 0

    def actionPerformed(self, event):
        if self.worker is not None and self.worker.isAlive():
            self.busy += 1


def button_labelled(container, text):
    buttons = [button for button in find_components(container, JButton) if as_text(button.getText()) == text]
    if len(buttons) != 1:
        raise RuntimeError("Expected exactly one button labelled '" + text + "', found " + str(len(buttons)))
    return buttons[0]


def press_and_watch(ctx, text, container=None, then=None):
    """Press a button and, in the same event-thread task, record the state
    of the window the moment the press has returned; then() runs there too.
    A job the button starts is then still running."""
    app = ctx.app
    button = button_labelled(container if container is not None else app, text)
    state = {}

    def press():
        if not button.isEnabled():
            raise RuntimeError("Button '" + text + "' is disabled")
        button.doClick(0)
        state["active"] = app.background_active
        state["running"] = app.batch_running
        state["thread"] = app.batch_thread
        state["edt"] = JavaThread.currentThread().getName()
        state["cancel"] = app.job_cancel_button.isEnabled()
        state["actions"] = [item.isEnabled() for item in app.analysis_action_buttons]
        state["batch_run"] = app.batch_run_button.isEnabled()
        state["batch_cancel"] = app.batch_cancel_button.isEnabled()
        if then is not None:
            then()

    on_edt(press, "button '" + text + "'")
    return state


def take_one(ctx, t, fragment, label, levels=("info", "warning")):
    """Wait for the job, then expect exactly one message box containing
    fragment, at one of the given levels."""
    wait_for_job(ctx)
    settle(ctx)
    records = ctx.watcher.take()
    texts = [record["level"] + ": " + record["text"].replace("\n", " | ") for record in records]
    t.check(len(records) == 1 and records[0]["level"] in levels and fragment in records[0]["text"],
            label + ": one message box with '" + fragment + "': " + repr(texts))
    return records


def compare_files(t, label, actual_dir, expected_dir, name, skip_words):
    """compare_tables on one file of two folders, leaving out every column
    whose name contains one of skip_words (run names, times, paths)."""
    actual_path = os.path.join(actual_dir, name)
    expected_path = os.path.join(expected_dir, name)
    if not (os.path.isfile(actual_path) and os.path.isfile(expected_path)):
        t.check(False, label + ": " + name + " is missing in " + actual_dir + " or " + expected_dir)
        return
    actual = read_csv(actual_path)
    expected = read_csv(expected_path)
    skip = [header for header in expected[0] if len([word for word in skip_words if word in header]) > 0]
    compare_tables(t, label + " " + name, actual, expected, skip)


def make_tissue_image(title, value):
    """Three channels with one disk of radius 70 pixels and the given value
    on a background of 100 (value 100: a uniform image, no tissue)."""
    stack = ImageStack(WIDTH, HEIGHT)
    for channel in range(3):
        ip = blank_plane(100)
        ip.setValue(float(value))
        ip.fill(OvalRoi(120 - 70, 100 - 70, 141, 141))
        ip.resetRoi()
        stack.addSlice("", ip)
    imp = ImagePlus(title, stack)
    imp.setDimensions(3, 1, 1)
    imp = CompositeImage(imp, CompositeImage.GRAYSCALE)
    cal = imp.getCalibration()
    cal.pixelWidth = PIXEL
    cal.pixelHeight = PIXEL
    cal.setUnit("um")
    return imp


def set_module_values(ctx, values):
    """Type values into Analyses-tab controls (event thread)."""
    app = ctx.app

    def apply():
        for key in values:
            kind, widget = app.module_widgets[key]
            if kind == "bool":
                widget.setSelected(bool(values[key]))
            elif kind == "choice":
                widget.setSelectedItem(values[key])
                if as_text(widget.getSelectedItem()) != as_text(values[key]):
                    raise RuntimeError("Choice '" + as_text(values[key]) + "' is not offered for " + key)
            else:
                widget.setText(as_text(values[key]))

    on_edt(apply, "typing Analyses settings")


def restore_parent_rois(ctx):
    def restore():
        ctx.manager.reset()
        ctx.manager.addRoi(triangle_roi())
        ctx.manager.addRoi(rect_roi())
        ctx.app.roi_source_combo.setSelectedItem("ROI Manager")

    on_edt(restore, "restoring ROI Manager")


def step_background_analysis(ctx, t):
    """'3. Run analysis' works on a worker thread: the window keeps
    answering while it runs, Cancel stops it, and its numbers equal a direct
    engine call."""
    app = ctx.app
    frames = 40
    ctx.imp5 = make_image("smoke long.tif", frames=frames)
    use_as_source(ctx, ctx.imp5)
    keep_masks = app.keep_masks_open_cb.isSelected()

    def setup():
        app.frame_mode_combo.setSelectedItem("All frames")
        # Showing 160 mask windows afterwards is not what this step is about.
        app.keep_masks_open_cb.setSelected(False)

    on_edt(setup, "choosing all frames")

    probe = TimerProbe()
    timer = SwingTimer(10, probe)
    on_edt(timer.start, "starting a Swing timer")
    state = press_and_watch(ctx, "3. Run analysis", then=lambda: setattr(probe, "worker", app.batch_thread))
    t.check(state.get("active") and state.get("running"), "the analysis is still running when the button press has returned")
    thread = state.get("thread")
    t.check(thread is not None and as_text(thread.getName()).startswith("cellquant-analysis"),
            "the analysis runs on its own thread: " + repr(thread))
    t.check(thread is not None and as_text(thread.getName()) != as_text(state.get("edt")), "that thread is not the event thread")
    t.check(state.get("cancel"), "Cancel is enabled while the analysis runs")
    t.check(len(state.get("actions", [])) > 0 and True not in state["actions"], "the Analyses buttons are disabled while it runs")
    # Questions to the event thread while the worker runs; an answer counts
    # when the worker was still running at the moment it was given.
    answers = []
    started = time.time()
    while thread is not None and thread.isAlive() and len(answers) < 25:
        sent = time.time()
        alive = on_edt(lambda: thread.isAlive(), "the event thread during an analysis", 10.0)
        if alive:
            answers.append((round(time.time() - sent, 3), round(sent - started, 2)))
        time.sleep(0.02)
    wait_for_job(ctx)
    settle(ctx)
    records = ctx.watcher.take()
    on_edt(timer.stop, "stopping the Swing timer")
    texts = [record["level"] + ": " + record["text"].replace("\n", " | ") for record in records]
    t.check(len(records) == 1 and records[0]["level"] == "info" and "Analysis complete." in records[0]["text"],
            "the analysis reports its end in one message box: " + repr(texts))
    t.check(probe.busy >= 3, "a Swing timer kept firing while the analysis ran (" + str(probe.busy) + " ticks)")
    t.check(len(answers) >= 3 and max([answer[0] for answer in answers]) < 1.0,
            "the event thread answered within 1 s while the analysis ran (seconds, at): " + repr(answers))
    t.check(not app.job_cancel_button.isEnabled() and False not in [b.isEnabled() for b in app.analysis_action_buttons],
            "Cancel is disabled and the Analyses buttons enabled again afterwards")
    t.equal(as_text(app.job_progress.getString()), "Finished", "the job bar reports the end")
    window = window_tables(as_text(app.last_run_dir))
    reference_image = make_image("smoke long.tif", frames=frames)
    engine = engine_reference(ctx, reference_image, "", expected_configs(ctx), [triangle_roi(), rect_roi()],
                              list(range(1, frames + 1)), {"frame_mode": "All frames"})
    reference_image.close()
    compare_run(t, "background analysis of " + str(frames) + " frames", window, engine)
    t.equal(len(window["summary"][1]), 2 * frames, "one row per parent ROI and frame")

    # Cancel, pressed straight after the start.
    images_before = WindowManager.getImageCount()
    state = press_and_watch(ctx, "3. Run analysis", then=lambda: app.job_cancel_button.doClick(0))
    t.check(state.get("active"), "the second analysis started")
    take_one(ctx, t, "Analysis cancelled.", "a cancelled analysis", ("info",))
    run_dir = as_text(app.last_run_dir)
    t.check(os.path.isdir(run_dir) and not os.path.isfile(os.path.join(run_dir, "roi_summary.csv")),
            "a cancelled analysis writes no result tables")
    t.check(os.path.isfile(os.path.join(run_dir, "run_log.txt")) and "cancelled by user" in read_text(os.path.join(run_dir, "run_log.txt")),
            "run_log.txt records the cancel")
    t.check("Analysis cancel requested" in as_text(app.status_area.getText()), "the cancel is logged")
    t.equal(as_text(app.job_progress.getString()), "Cancelled", "the job bar reports the cancel")
    t.equal((app.background_active, app.batch_running), (False, False), "the window is ready again after the cancel")
    t.equal(WindowManager.getImageCount(), images_before, "a cancelled analysis leaves no image open")

    def restore():
        app.frame_mode_combo.setSelectedItem("Current frame only")
        app.keep_masks_open_cb.setSelected(keep_masks)

    on_edt(restore, "back to the current frame")

    # The window stays usable during a run, so it checks afterwards that the
    # images it measured were left alone, and says so when they were not.
    def edit_source():
        ctx.imp5.changes = True

    press_and_watch(ctx, "3. Run analysis", then=edit_source)
    records = take_one(ctx, t, "Analysis complete.", "an analysis whose source was edited while it ran", ("warning",))
    t.check(len(records) == 1 and "the pixels of the source image were edited" in records[0]["text"],
            "the closing message names the change: " + repr([r["text"] for r in records]))
    t.check("the pixels of the source image were edited" in read_text(os.path.join(as_text(app.last_run_dir), "run_log.txt")),
            "run_log.txt records the change")
    ctx.imp5.changes = False
    # The other checks, on the check itself (a selection drawn during a real
    # run could race with the worker).
    configs = app.collect_channel_configs()
    configs[2]["projection"] = "Current Z"

    def watched(change):
        job = {"label": "Analysis", "values": {"watch": app.image_watch(ctx.imp5, configs, {})}}
        change()
        return app.changed_during_job(job)

    t.equal(on_edt(lambda: watched(lambda: None), "checking an unchanged image"), "", "no warning when nothing changed")
    text = on_edt(lambda: watched(lambda: ctx.imp5.setRoi(Roi(5, 5, 20, 20))), "checking a new selection")
    t.check("the selection on the source image was changed" in text, "a new selection is reported: " + text)
    text = on_edt(lambda: watched(lambda: ctx.imp5.setPosition(1, 2, 1)), "checking a moved Z slider")
    t.check("the Z slider of the source image was moved" in text, "a moved Z slider is reported for 'Current Z': " + text)

    def restore_position():
        ctx.imp5.deleteRoi()
        ctx.imp5.setPosition(1, 1, 1)

    on_edt(restore_position, "removing the selection")
    use_as_source(ctx, ctx.imp1)
    on_edt(lambda: ctx.imp5.close(), "closing the long image")
    expect_dialogs(ctx, t, [], "background analysis")


def step_analyses_tab(ctx, t):
    """The Analyses tab: one control per field of every settings form the
    modules register, saved and loaded with profiles."""
    app = ctx.app
    E = ctx.E
    forms = E.SETTINGS_FORMS
    t.check(len(forms) > 0, "the loaded modules register settings forms")
    tab = app.main_tabs.getComponentAt(TAB_TITLES.index("Analyses"))
    inside = find_components(tab, Component)
    defaults = E.default_global_options()
    registered = [spec[1] for spec in E.GLOBAL_SETTING_SPECS]
    keys = []
    for form in forms:
        for field in form["fields"]:
            key = field[0]
            kind = field[2]
            keys.append(key)
            entry = app.module_widgets.get(key)
            t.check(entry is not None, "a control for " + key)
            t.check(key in registered, key + " is a registered setting")
            if entry is None:
                continue
            got_kind, widget = entry
            t.equal(got_kind, kind, key + ": kind of control")
            t.check(len([item for item in inside if item is widget]) == 1, key + ": the control is on the Analyses tab")
            if kind == "bool":
                t.check(isinstance(widget, JCheckBox), key + " is a check box")
                t.equal(as_text(widget.getText()), as_text(field[1]), key + ": label")
                t.equal(widget.isSelected(), E.parse_bool(defaults.get(key), False), key + ": shows its default")
            elif kind == "choice":
                t.check(isinstance(widget, JComboBox), key + " is a list")
                t.equal([as_text(widget.getItemAt(i)) for i in range(widget.getItemCount())],
                        [as_text(choice) for choice in field[3]], key + ": the allowed values")
                t.equal(as_text(widget.getSelectedItem()), as_text(defaults.get(key)), key + ": shows its default")
            else:
                t.check(isinstance(widget, JTextField), key + " is a text field")
                t.equal(as_text(widget.getText()), as_text(defaults.get(key)), key + ": shows its default")
    t.equal(sorted(app.module_widgets.keys()), sorted(keys), "one control per field and no other")
    borders = [as_text(panel.getBorder().getTitle()) for panel in find_components(tab, JPanel)
               if isinstance(panel.getBorder(), TitledBorder)]
    for form in forms:
        t.check(as_text(form["title"]) in borders, "a section titled '" + as_text(form["title"]) + "'")
    for function, label in ANALYSIS_ACTIONS:
        if function in ctx.ns:
            buttons = [b for b in find_components(tab, JButton) if as_text(b.getText()) == label]
            t.check(len(buttons) == 1 and buttons[0].isEnabled(), "the Analyses tab has the button '" + label + "', enabled")
    t.equal(len(app.analysis_action_buttons), len([a for a in ANALYSIS_ACTIONS if a[0] in ctx.ns]), "one action button per module run")

    # Typing values, then a profile round trip.
    base = os.path.join(ctx.tmp, "analyses base.properties")
    on_edt(lambda: app.save_settings_to_path(base), "saving the profile before typing")
    typed = {}

    def type_values():
        for key in keys:
            kind, widget = app.module_widgets[key]
            if kind == "bool":
                widget.setSelected(not widget.isSelected())
                typed[key] = widget.isSelected()
            elif kind == "choice":
                widget.setSelectedIndex((widget.getSelectedIndex() + 1) % widget.getItemCount())
                typed[key] = as_text(widget.getSelectedItem())
            else:
                widget.setText("typed " + key)
                typed[key] = u"typed " + key

    on_edt(type_values, "typing every Analyses setting")
    options = app.collect_global_options()
    for key in keys:
        t.equal(options.get(key), typed[key], key + ": the typed value is a window setting")
    path = os.path.join(ctx.tmp, "analyses typed.properties")
    on_edt(lambda: app.save_settings_to_path(path), "saving the typed profile")
    props = E.read_properties_file(path)
    for spec in E.GLOBAL_SETTING_SPECS:
        if spec[1] in typed:
            expected = E.bool_text(typed[spec[1]]) if spec[2] == "bool" else typed[spec[1]]
            t.equal(as_text(props.getProperty(spec[0])), as_text(expected), spec[0] + " is saved in the profile")
    on_edt(lambda: app.load_profile_from_path(base), "loading the profile from before typing")
    after = app.module_widget_options()
    for key in keys:
        if app.module_widgets[key][0] == "bool":
            expected = E.parse_bool(defaults[key], False)
        else:
            expected = as_text(defaults[key])
        t.equal(after[key], expected, key + ": loading the earlier profile restores the default")
    on_edt(lambda: app.load_profile_from_path(path), "loading the typed profile")
    after = app.module_widget_options()
    for key in keys:
        t.equal(after[key], typed[key], key + ": loading the profile shows the typed value")
    engine_options = E.load_profile_settings(path)[0]
    for key in keys:
        value = engine_options.get(key)
        if app.module_widgets[key][0] == "bool":
            value = E.parse_bool(value, None)
        t.equal(value, typed[key], key + ": the engine reads the same value from the profile")
    if "tissue_sigma" in keys:
        t.check("tissue.sigma=typed tissue_sigma" in as_text(app.status_area.getText()),
                "the status log names the module settings of the loaded profile")
    on_edt(lambda: app.load_profile_from_path(base), "loading the earlier profile again")

    # A module reports its own bad values before any job starts.
    if "tissue_sigma" in keys and "detect_tissue_outline" in ctx.ns:
        use_as_source(ctx, ctx.imp1)
        set_module_values(ctx, {"tissue_sigma": "abc"})
        click(ctx, "Detect tissue outline into ROI Manager")
        expect_dialogs(ctx, t, [("error", "Tissue outline detection failed")], "a non-numeric tissue setting")
        t.equal(app.background_active, False, "a refused job does not start")
        set_module_values(ctx, {"tissue_sigma": defaults["tissue_sigma"]})
    expect_dialogs(ctx, t, [], "Analyses tab")


def step_tissue(ctx, t):
    app = ctx.app
    E = ctx.E
    if "detect_tissue_outline" not in ctx.ns:
        t.check(False, "the tissue module (cellquanttissue.py) is loaded")
        return
    ctx.tissue_imp = make_tissue_image("smoke tissue.tif", 1200)
    use_as_source(ctx, ctx.tissue_imp)
    set_module_values(ctx, {"tissue_channel": "1"})
    before = ctx.manager.getCount()
    state = press_and_watch(ctx, "Detect tissue outline into ROI Manager")
    t.check(state.get("active"), "the detection runs in the background")
    expect_dialogs(ctx, t, [], "tissue outline of a disk")
    t.equal(ctx.manager.getCount(), before + 1, "one tissue outline is added to ROI Manager")
    t.check(RoiManager.getInstance2() is ctx.manager, "ROI Manager is still the same window")
    roi = ctx.manager.getRoi(ctx.manager.getCount() - 1) if ctx.manager.getCount() > before else None
    t.check(roi is not None and as_text(roi.getName()).startswith("Tissue"), "the outline is named Tissue_...: " + repr(roi))
    area = len(roi.getContainedPoints()) if roi is not None else 0
    rim = [(int(120 + 68 * math.cos(a / 10.0)), int(100 + 68 * math.sin(a / 10.0))) for a in range(63)]
    t.check(roi is not None and False not in [roi.contains(x, y) for x, y in rim + [(120, 100)]],
            "the outline holds the whole disk (" + str(area) + " pixels; the disk has " + str(int(math.pi * 70 * 70)) + ")")
    t.equal(as_text(app.roi_source_combo.getSelectedItem()), "ROI Manager", "the ROI source becomes ROI Manager")
    log = as_text(app.status_area.getText())
    t.check(E.TISSUE_POLICY + ": 1 region(s)" in log, "the description of the detection is in the status log")
    copy = make_tissue_image("smoke tissue.tif", 1200)
    rois, info = E.detect_tissue_outline(E.AnalysisEngine(), copy, app.collect_channel_configs(), app.collect_global_options(), frame=1)
    copy.close()
    t.equal(len(rois), 1, "a direct call finds one outline too")
    if len(rois) == 1 and roi is not None:
        t.equal(len(rois[0].getContainedPoints()), area, "the window's outline equals the direct call's, pixel for pixel")
        t.equal(as_text(rois[0].getBounds()), as_text(roi.getBounds()), "same bounds")

    ctx.blank_imp = make_tissue_image("smoke blank.tif", 100)
    use_as_source(ctx, ctx.blank_imp)
    click(ctx, "Detect tissue outline into ROI Manager")
    expect_dialogs(ctx, t, [("warning", "No tissue outline was found")], "a uniform image")
    t.equal(ctx.manager.getCount(), before + 1, "nothing is added for a uniform image")

    restore_parent_rois(ctx)
    set_module_values(ctx, {"tissue_channel": E.default_global_options()["tissue_channel"]})
    use_as_source(ctx, ctx.imp1)
    on_edt(lambda: ctx.tissue_imp.close(), "closing the tissue image")
    on_edt(lambda: ctx.blank_imp.close(), "closing the blank image")
    expect_dialogs(ctx, t, [], "tissue outline")


def write_annotated_field(folder, title, count):
    """An image with count disks and its point marks: one per disk except
    the last, and one on the background."""
    image = make_image(title, disks=DISKS[:count])
    save_image(image, os.path.join(folder, title))
    image.close()
    points = [DISKS[i] for i in range(count - 1)] + [(20, 110)]
    marks = PointRoi(points[0][0], points[0][1])
    for point in points[1:]:
        marks.addPoint(point[0], point[1])
    return marks


def step_tuning(ctx, t):
    app = ctx.app
    E = ctx.E
    if "run_tuning_folder" not in ctx.ns:
        t.check(False, "the tuning module (cellquanttune.py) is loaded")
        return
    folder = os.path.join(ctx.tmp, "tune_in")
    os.makedirs(folder)
    for title, count in [("field_01.tif", 15), ("field_02.tif", 7)]:
        marks = write_annotated_field(folder, title, count)
        E.save_roi_array([marks], os.path.join(folder, E.strip_extension(title) + "_points.zip"))

    def setup():
        app.tune_folder_field.setText(folder)
        app.batch_extensions_field.setText("tif")
        app.batch_recursive_cb.setSelected(False)
        app.batch_roi_dir_field.setText("")
        app.output_dir_field.setText(ctx.out_dir)

    on_edt(setup, "choosing the annotated folder")
    set_module_values(ctx, {"tune_search_space": "manual_low = 800, 1000, 3000", "tune_channel": "0"})
    ctx.watcher.answer_next("no")
    state = press_and_watch(ctx, "Run tuning...")
    t.check(state.get("active"), "tuning runs in the background")
    records = expect_dialogs(ctx, t, [("question", "Tuning complete.")], "tuning")
    t.check(len(records) == 1 and "Load the tuned profile" in records[0]["text"] and "Held-out F1" in records[0]["text"],
            "the window offers the tuned profile and shows the held-out score: " + repr([r["text"] for r in records]))
    run_dir = as_text(app.last_run_dir)
    t.check(os.path.basename(run_dir).startswith("CellQuant_Tuning_") and run_dir.startswith(ctx.out_dir),
            "the tuning folder is in the output folder: " + run_dir)
    for name in ["tuned_profile.properties", "tuning_summary.csv", "tuning_grid.csv", "tuning_log.txt", "settings.properties"]:
        t.check(os.path.isfile(os.path.join(run_dir, name)), name + " is written")
    if os.path.isfile(os.path.join(run_dir, "tuning_summary.csv")):
        summary = read_csv(os.path.join(run_dir, "tuning_summary.csv"))
        t.equal(summary[1][0].get("Fields") if len(summary[1]) > 0 else None, "2", "both annotated fields are used")
    t.equal(app.channel_panels[0].get_config()["manual_low"], 1000.0, "answering No leaves the window's settings alone")
    t.equal(as_text(app.job_progress.getString()), "Finished", "the job bar reports the end")

    reference = os.path.join(ctx.tmp, "reference", "tuning")
    os.makedirs(reference)
    result = E.run_tuning_folder(
        E.AnalysisEngine(), folder, reference, app.collect_channel_configs(), app.collect_all_channel_configs(),
        app.collect_global_options(), {"extensions": "tif", "recursive": False, "roi_dir": None, "force_bioformats": False}
    )
    for name in ["tuning_summary.csv", "tuning_grid.csv", "tuning_field_scores.csv"]:
        compare_files(t, "tuning", run_dir, as_text(result["run_dir"]), name, ["Run_ID", "Seconds", "Per_Second", "Image_Path"])

    # An error inside the worker reaches a message box.
    empty = os.path.join(ctx.tmp, "tune_no_marks")
    os.makedirs(empty)
    image = make_image("field_01.tif", disks=DISKS[:5])
    save_image(image, os.path.join(empty, "field_01.tif"))
    image.close()
    on_edt(lambda: app.tune_folder_field.setText(empty), "choosing a folder without marks")
    state = press_and_watch(ctx, "Run tuning...")
    t.check(state.get("active"), "the job starts; the module finds the problem on the worker thread")
    expect_dialogs(ctx, t, [("error", "Tuning failed: No annotated image was found")], "tuning without marks")
    t.equal(as_text(app.job_progress.getString()), "Failed", "the job bar reports the failure")
    t.equal((app.background_active, app.batch_running), (False, False), "the window is ready again after a failed job")
    t.check(False not in [b.isEnabled() for b in app.analysis_action_buttons], "the Analyses buttons are enabled after a failed job")
    on_edt(lambda: app.tune_folder_field.setText(""), "clearing the tuning folder")
    defaults = E.default_global_options()
    set_module_values(ctx, {"tune_search_space": defaults["tune_search_space"], "tune_channel": defaults["tune_channel"]})


def step_correction(ctx, t):
    app = ctx.app
    E = ctx.E
    if "run_count_correction" not in ctx.ns:
        t.check(False, "the count-correction module (cellquantcorrect.py) is loaded")
        return
    table = os.path.join(ctx.batch_dir, "batch_roi_summary.csv")
    original = read_text(table)

    def setup():
        app.correct_source_field.setText(ctx.out_dir)
        app.correct_table_field.setText(table)
        app.correct_channel_field.setText("1")
        app.output_dir_field.setText(ctx.out_dir)

    on_edt(setup, "choosing the validation runs and the table")
    state = press_and_watch(ctx, "Run count correction...")
    t.check(state.get("active"), "the count correction runs in the background")
    take_one(ctx, t, "Count correction complete.", "count correction")
    out_dir = as_text(app.last_run_dir)
    t.check("_CellQuant_CountCorrection_" in os.path.basename(out_dir) and out_dir.startswith(ctx.out_dir),
            "a new correction folder in the output folder: " + out_dir)
    for name in ["count_calibration.csv", "count_calibration_fit.csv", "batch_roi_summary_corrected.csv"]:
        t.check(os.path.isfile(os.path.join(out_dir, name)), name + " is written")
    t.equal(read_text(table), original, "the table that was corrected is left unchanged")
    if os.path.isfile(os.path.join(out_dir, "count_calibration.csv")):
        pairs = read_csv(os.path.join(out_dir, "count_calibration.csv"))
        t.equal(len(pairs[1]), 4, "the two object validations of the validation step give four calibration fields")
    t.check("Count calibration" in as_text(app.status_area.getText()), "the module's messages reach the status log")

    reference = os.path.join(ctx.tmp, "reference", "correction")
    os.makedirs(reference)
    result = E.run_count_correction(
        calibration_folder=ctx.out_dir, table_path=table, channel="1", output_root=reference,
        global_options=app.collect_global_options()
    )
    for name in ["batch_roi_summary_corrected.csv", "count_calibration.csv", "count_calibration_fit.csv"]:
        compare_files(t, "correction", out_dir, as_text(result["out_dir"]), name, ["Run_ID", "Source_File", "Path"])

    def clear():
        app.correct_source_field.setText("")
        app.correct_table_field.setText("")
        app.correct_channel_field.setText("")

    on_edt(clear, "clearing the correction inputs")
    expect_dialogs(ctx, t, [], "count correction")


def write_group_folder(ctx, name, titles_and_counts):
    folder = os.path.join(ctx.tmp, name)
    os.makedirs(folder)
    for title, count in titles_and_counts:
        image = make_image(title, disks=DISKS[:count])
        save_image(image, os.path.join(folder, title))
        image.close()
    return folder


def step_robustness(ctx, t):
    app = ctx.app
    E = ctx.E
    if "run_robustness_from_options" not in ctx.ns:
        t.check(False, "the robustness module (cellquantrobust.py) is loaded")
        return
    folder = write_group_folder(ctx, "robust_in", [("ctrl_01.tif", 15), ("ctrl_02.tif", 13), ("treat_01.tif", 5), ("treat_02.tif", 7)])

    def setup(path):
        app.batch_input_field.setText(path)
        app.batch_extensions_field.setText("tif")
        app.batch_roi_combo.setSelectedItem(E.BATCH_ROI_WHOLE)
        app.batch_roi_dir_field.setText("")
        app.batch_recursive_cb.setSelected(False)
        app.output_dir_field.setText(ctx.out_dir)

    on_edt(lambda: setup(folder), "filling the Batch tab")
    set_module_values(ctx, {"robust_channel": "1", "robust_grid": "threshold_factor=0.9,1.1"})
    state = press_and_watch(ctx, "Run robustness...")
    t.check(state.get("active"), "the robustness analysis runs in the background")
    t.equal((state.get("batch_run"), state.get("batch_cancel")), (False, True),
            "it uses the Batch tab: Run batch is disabled and its Cancel enabled")
    expect_dialogs(ctx, t, [("question", "Run the robustness analysis on 4 image file(s)"), ("info", "Robustness analysis complete.")],
                   "robustness")
    run_dir = as_text(app.last_run_dir)
    t.check(os.path.basename(run_dir).startswith("CellQuant_Robustness_") and run_dir.startswith(ctx.out_dir),
            "the robustness folder is in the output folder: " + run_dir)
    for name in ["robustness_summary.csv", "robustness_specifications.csv", "robustness_groups.csv", "robustness_log.txt"]:
        t.check(os.path.isfile(os.path.join(run_dir, name)), name + " is written")
    t.equal(as_text(app.batch_progress.getString()), "Finished", "the Batch tab's progress bar reports the end")
    t.equal((app.batch_run_button.isEnabled(), app.batch_cancel_button.isEnabled()), (True, False), "the Batch tab buttons are restored")
    reference = os.path.join(ctx.tmp, "reference", "robustness")
    os.makedirs(reference)
    result = E.run_robustness_from_options(
        folder, reference, app.collect_channel_configs(), app.collect_all_channel_configs(), app.collect_global_options(),
        app.collect_batch_options()
    )
    for name in ["robustness_summary.csv", "robustness_specifications.csv", "robustness_groups.csv", "robustness_units.csv"]:
        compare_files(t, "robustness", run_dir, as_text(result["run_dir"]), name, ["Run_ID", "Seconds", "Duration", "Path"])

    # Cancel, pressed straight after the start.
    many = write_group_folder(ctx, "robust_cancel", [
        (("ctrl_" if index % 2 == 0 else "treat_") + str(index + 1).zfill(2) + ".tif", 5 + index) for index in range(8)
    ])
    on_edt(lambda: setup(many), "choosing a larger folder")
    state = press_and_watch(ctx, "Run robustness...", then=lambda: app.batch_cancel_button.doClick(0))
    t.check(state.get("active"), "the second robustness analysis started")
    expect_dialogs(ctx, t, [("question", "8 image file(s)"), ("info", "Robustness analysis cancelled")], "a cancelled robustness analysis")
    run_dir = as_text(app.last_run_dir)
    t.check(os.path.isdir(run_dir) and not os.path.isfile(os.path.join(run_dir, "robustness_summary.csv")),
            "a cancelled robustness analysis writes no statistics")
    t.equal(as_text(app.batch_progress.getString()), "Cancelled", "the Batch tab's progress bar reports the cancel")
    t.check("Robustness analysis cancel requested" in as_text(app.status_area.getText()), "the cancel is logged")

    # One group only: the module refuses on the worker thread; the window says why.
    one = write_group_folder(ctx, "robust_one", [("ctrl_01.tif", 5), ("ctrl_02.tif", 6)])
    on_edt(lambda: setup(one), "choosing a folder of one group")
    click(ctx, "Run robustness...")
    expect_dialogs(ctx, t, [("question", "2 image file(s)"),
                            ("error", "Robustness analysis failed: A comparison needs images of at least two groups")],
                   "robustness with one group")
    t.equal((app.background_active, app.batch_running, app.batch_run_button.isEnabled()), (False, False, True),
            "the window is ready again after the refusal")

    def clear():
        app.batch_input_field.setText("")
        app.batch_roi_combo.setSelectedItem(E.BATCH_ROI_WHOLE)

    on_edt(clear, "clearing the Batch tab")
    defaults = E.default_global_options()
    set_module_values(ctx, {"robust_channel": defaults["robust_channel"], "robust_grid": defaults["robust_grid"]})
    expect_dialogs(ctx, t, [], "robustness")


def step_other_actions(ctx, t):
    app = ctx.app
    use_as_source(ctx, ctx.imp1)
    on_edt(lambda: ctx.imp1.setRoi(Roi(20, 20, 50, 40)), "drawing a rectangle")
    click(ctx, "Add current ROI")
    t.equal(ctx.manager.getCount(), 3, "'Add current ROI' adds the image selection to ROI Manager")
    on_edt(ctx.imp1.deleteRoi, "removing the selection")
    click(ctx, "Add current ROI")
    expect_dialogs(ctx, t, [("warning", "no active ROI")], "'Add current ROI' without a selection")
    ctx.watcher.answer_next("no")
    click(ctx, "Clear ROI Manager")
    expect_dialogs(ctx, t, [("question", "Clear all ROIs")], "declining 'Clear ROI Manager'")
    t.equal(ctx.manager.getCount(), 3, "declining keeps the ROIs")
    click(ctx, "Clear ROI Manager")
    expect_dialogs(ctx, t, [("question", "Clear all ROIs")], "'Clear ROI Manager'")
    t.equal(ctx.manager.getCount(), 0, "the ROIs are cleared")
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("error", "ROI Manager contains no parent ROIs")], "analysis with an empty ROI Manager")
    on_edt(lambda: app.roi_source_combo.setSelectedItem("Whole image"), "choosing the whole image")
    click(ctx, "3. Run analysis")
    expect_dialogs(ctx, t, [("info", "Analysis complete.")], "whole-image analysis")
    summary = read_csv(os.path.join(as_text(app.last_run_dir), "roi_summary.csv"))
    t.equal([row.get("Parent_ROI_Name") for row in summary[1]], ["Whole_Image"], "the whole image is one parent ROI")

    click(ctx, "1. Prepare enabled channels")
    generated = list(app.generated_ids)
    t.check(len(generated) >= 3, "generated images are tracked")
    click(ctx, "Close generated images")
    t.equal([image_id for image_id in generated if WindowManager.getImage(image_id) is not None], [], "'Close generated images' closes them")
    t.check(WindowManager.getImage(ctx.imp1.getID()) is not None, "the source image stays open")

    ctx.watcher.answer_next("no")
    click(ctx, "Reset for next image")
    expect_dialogs(ctx, t, [("question", "Close the current source image")], "'Reset for next image'")
    t.check(app.source_imp is None, "the reset forgets the source image")
    t.check(WindowManager.getImage(ctx.imp1.getID()) is not None, "answering No keeps the image open")
    t.equal(as_text(app.image_info_area.getText()), "No source image selected.", "the source description is cleared")

    use_as_source(ctx, ctx.imp1)
    click(ctx, "Close all ImageJ images")
    expect_dialogs(ctx, t, [("warning", "Close every open image")], "'Close all ImageJ images'")
    t.equal(WindowManager.getImageCount(), 0, "every image is closed")
    t.check(app.source_imp is None, "the source is forgotten")
    click(ctx, "2. Preview all masks")
    expect_dialogs(ctx, t, [("error", "No source image is selected")], "previewing with no image")


def step_relaunch(ctx, t):
    launch = ctx.ns["launch_cell_quant_studio"]
    on_edt(ctx.app.dispose, "closing the window")
    settle(ctx)
    t.equal(len(find_app_windows()), 0, "the window closes")
    for attempt in range(5):
        window = off_edt(launch, "launch " + str(attempt + 1) + " from a script thread", 60.0)
        settle(ctx)
        t.check(window is not ctx.app and window.isShowing(), "launch " + str(attempt + 1) + " opens a new window")
        t.check("Applied preset" in as_text(window.status_area.getText()), "launch " + str(attempt + 1) + " shows its start-up messages")
        on_edt(window.dispose, "closing the window")
        settle(ctx)
    window = on_edt(launch, "launch on the event thread", 60.0)
    settle(ctx)
    t.check(window.isShowing(), "launching on the event thread works")
    on_edt(window.dispose, "closing the window")
    settle(ctx)
    t.equal(len(find_app_windows()), 0, "every window is closed again")
    expect_dialogs(ctx, t, [], "relaunch")


def step_fewer_channels(ctx, t):
    """An image with fewer channels must not erase the settings of the
    channels it lacks."""
    launch = ctx.ns["launch_cell_quant_studio"]
    app = on_edt(launch, "opening a window", 60.0)
    settle(ctx)
    ctx.app = app
    three = make_image("smoke three channels.tif")
    single = ImagePlus("smoke one channel.tif", nuclei_plane(DISKS))
    try:
        use_as_source(ctx, three)
        set_widgets(app.channel_panels[1], C2_VALUES)
        set_widgets(app.channel_panels[2], {"enabled": True, "name": "Third", "analysis_mode": "Intensity only"})
        use_as_source(ctx, single)
        t.equal([app.channel_tabs.isEnabledAt(i) for i in range(3)], [True, False, False],
                "the tabs of missing channels are disabled")
        second = app.channel_panels[1].get_config()
        third = app.channel_panels[2].get_config()
        t.equal((second["enabled"], second["name"], second["analysis_mode"], second["manual_low"]),
                (True, "MAP2", "Positive area", 800.0), "C2 keeps its settings behind a one-channel image")
        t.equal((third["enabled"], third["name"], third["analysis_mode"]), (True, "Third", "Intensity only"),
                "C3 keeps its settings behind a one-channel image")
        t.check(max([0] + [config["channel"] for config in app.collect_channel_configs()]) <= 1,
                "only the channels of the open image are analyzed")
        folder_channels = [config["channel"] for config in app.collect_channel_configs(limit_to_source=False)]
        t.check(2 in folder_channels and 3 in folder_channels, "a folder job still uses every enabled channel")
        use_as_source(ctx, three)
        back = [config["channel"] for config in app.collect_channel_configs()]
        t.check(2 in back and 3 in back, "the channels are used again with the next three-channel image")
        expect_dialogs(ctx, t, [], "switching between images")
    finally:
        for imp in [three, single]:
            imp.changes = False
            on_edt(imp.close, "closing a test image")
        on_edt(app.dispose, "closing the window")
        settle(ctx)


def step_look(ctx, t):
    """The dark look applies to the whole of Fiji, so it must be gone again
    once the window has closed; and the small helpers of the Results tab."""
    ns = ctx.ns
    t.equal(ns["format_count"](1234567.4), "1,234,567", "counts are written with thousands separators")
    t.equal(ns["format_count"](12), "12", "a small count is left alone")
    headers = ["Parent_ROI_Name", "ROI_Area", "C1_A_Object_Count", "C1_A_Edge_Touching_Object_Count",
               "C2_B_Positive_Area", "C2_B_Positive_Area_Fraction_Percent"]
    rows = [
        {"Parent_ROI_Name": "one", "ROI_Area": 100.0, "C1_A_Object_Count": 3, "C1_A_Edge_Touching_Object_Count": 1,
         "C2_B_Positive_Area": 10.0, "C2_B_Positive_Area_Fraction_Percent": 10.0},
        {"Parent_ROI_Name": "two", "ROI_Area": 300.0, "C1_A_Object_Count": 5, "C1_A_Edge_Touching_Object_Count": 0,
         "C2_B_Positive_Area": 150.0, "C2_B_Positive_Area_Fraction_Percent": 50.0}
    ]
    shown = ns["results_from_summary"](headers, rows, [])
    t.equal([(tile["title"], tile["kind"], tile["value"]) for tile in shown["tiles"]],
            [("A objects", "count", 8.0), ("B positive area", "percent", 40.0)],
            "counts are added up, fractions are weighted by area, edge counts get no tile")
    t.equal(ns["results_from_summary"]([], [], [])["tiles"], [], "an empty table gives an empty Results tab")
    t.equal(ns["results_pictures"](os.path.join(ctx.out_dir, "no such run")), [], "a run without overlays shows no picture")

    before = UIManager.getLookAndFeel()
    state = ns["LOOK_STATE"]
    t.check(not state["active"], "no look is switched on while no window is open")
    launch = ns["launch_cell_quant_studio"]
    window = on_edt(launch, "opening a window", 60.0)
    settle(ctx)
    wanted = ns["look_preference"]()
    if state["active"]:
        t.check(UIManager.getLookAndFeel() is not before, "the chosen look is on while the window is open")
        t.equal(ns["look_is_dark"](), wanted == ns["LOOK_DARK"], "it is dark exactly when Dark is chosen")
        for kind in ["ui", "ui_bold", "figure"]:
            t.check(kind in ns["LOOK_FONTS"], "the bundled font '" + kind + "' is loaded")
        t.check(ns["look_font"]("figure", 20).getSize() == 20, "look_font gives the asked size")
        t.equal(as_text(window.look_combo.getSelectedItem()), wanted, "the title row shows the look in use")
        t.equal(state["name"], wanted, "the look in use is the remembered choice")
        # Switching the look while the window is open: Fiji gets its own look
        # back at once, the tabs are shaped again for the other font, and
        # the labels keep their texts.
        some_label = [item for item in find_components(window, Component)
                      if as_text(type(item).__name__) == "JLabel" and item.getClientProperty("cellquant.text") is not None][0]
        label_text = as_text(some_label.getClientProperty("cellquant.text"))
        on_edt(lambda: window.look_combo.setSelectedItem(ns["LOOK_FIJI"]), "choosing Fiji's look")
        settle(ctx)
        t.check(not state["active"] and UIManager.getLookAndFeel() is before, "'Fiji default' gives Fiji its look back at once")
        on_edt(lambda: window.look_combo.setSelectedItem(wanted), "choosing the remembered look again")
        settle(ctx)
        t.check(state["active"] and state["name"] == wanted, "and the window's look comes back")
        t.equal(as_text(some_label.getClientProperty("cellquant.text")), label_text, "labels keep their text through both changes")
        t.equal(ns["look_preference"](), wanted, "the remembered choice is as it was")
    else:
        t.check(UIManager.getLookAndFeel() is before, "without a look of its own nothing about Fiji changes")
    # Text that does not fit is shortened, never drawn past its room.
    from java.awt.image import BufferedImage
    pen = BufferedImage(10, 10, BufferedImage.TYPE_INT_RGB).createGraphics()
    wide = pen.getFontMetrics().stringWidth("a rather long caption for a small card")
    t.check(0 < ns["paint_fitted"](pen, "a rather long caption for a small card", 0, 9, wide // 2) <= wide // 2,
            "a caption that is too long is shortened to its room")
    t.equal(ns["paint_fitted"](pen, "short", 0, 9, 500), pen.getFontMetrics().stringWidth("short"), "a caption that fits is drawn whole")
    t.equal(ns["paint_fitted"](pen, "anything", 0, 9, 0), 0, "nothing is drawn where there is no room")
    pen.dispose()
    on_edt(window.dispose, "closing the window")
    settle(ctx)
    wait_until(lambda: not state["active"], 10.0, "the look to be put back")
    t.check(UIManager.getLookAndFeel() is before, "Fiji has its own look again after the window closed")
    expect_dialogs(ctx, t, [], "switching the look")


def step_source_rules(ctx, t):
    for name in ["cellquantcore.py", "Cell_Quant_Studio_final.py"]:
        handle = open(os.path.join(ctx.repo, name), "rb")
        try:
            data = handle.read()
        finally:
            handle.close()
        text = data.decode("latin-1")
        t.check(len([char for char in text if ord(char) > 127]) == 0, name + " is plain ASCII")
        t.check("invokeAnd" + "Wait" not in text, name + " never blocks a worker thread on the event thread")
        t.check("except Exception" not in text, name + " catches Java errors too (no 'except Exception')")
        t.check("coding:" not in text.split("\n")[0] and "coding:" not in text.split("\n")[1], name + " has no encoding line")


STEPS = [
    ("launch the application", step_launch),
    ("a module of another version is refused", step_module_version),
    ("status messages and error reports from both threads", step_status_messages),
    ("select a source image", step_source_image),
    ("configure channels, output folder, and parent ROIs", step_configure),
    ("prepare channels and preview masks", step_prepare_and_preview),
    ("run analysis and compare with the engine", step_run_analysis),
    ("second image without touching ROI Manager", step_second_image),
    ("prepared-image cache", step_prepared_cache),
    ("display-minimum threshold workflow", step_display_minimum),
    ("source or prepared image edited by hand", step_edited_source),
    ("current frame and all frames", step_all_frames),
    ("save and load a profile", step_profile),
    ("analysis-module settings in profiles and runs", step_module_settings),
    ("validation against marks and a reference mask", step_validation),
    ("parameter sensitivity", step_sensitivity),
    ("open an image from a file", step_open_from_file),
    ("folder batch", step_batch),
    ("cancel a batch", step_batch_cancel),
    ("analysis in the background: responsive, cancel, same numbers", step_background_analysis),
    ("Analyses tab: one control per module setting, profiles", step_analyses_tab),
    ("tissue outline into ROI Manager", step_tissue),
    ("tuning from point marks", step_tuning),
    ("count correction from the validation runs", step_correction),
    ("robustness over a grid of settings, cancel, refusal", step_robustness),
    ("ROI, cleanup, and reset actions", step_other_actions),
    ("close and relaunch", step_relaunch),
    ("an image with fewer channels keeps the other channels' settings", step_fewer_channels),
    ("the dark look is put back; Results tab helpers", step_look),
    ("source rules", step_source_rules)
]


# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------

def find_repo():
    candidates = []
    if os.environ.get("CELLQUANT_PACKAGE_DIR"):
        candidates.append(os.environ.get("CELLQUANT_PACKAGE_DIR"))
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        candidates.append(os.path.dirname(here))
        candidates.append(here)
    except ANY_ERROR:
        pass
    try:
        root = IJ.getDirectory("imagej")
        if root:
            candidates.append(os.path.join(root, "plugins", "Scripts", "Cell Quant Studio"))
    except ANY_ERROR:
        pass
    for candidate in candidates:
        if candidate and os.path.isfile(os.path.join(candidate, "Cell_Quant_Studio_final.py")):
            return candidate
    raise IOError("Could not find Cell_Quant_Studio_final.py. Keep the tests folder inside the Cell Quant Studio folder.")


def close_everything(ctx):
    """Close every image and window this test may have opened."""
    def close_all():
        ids = WindowManager.getIDList()
        if ids is not None:
            for image_id in list(ids):
                imp = WindowManager.getImage(image_id)
                if imp is not None:
                    imp.changes = False
                    imp.close()
        for window in Window.getWindows():
            if window in ctx.windows_before:
                continue
            try:
                window.dispose()
            except ANY_ERROR:
                pass

    try:
        on_edt(close_all, "closing every window", 30.0)
    except ANY_ERROR:
        print("Could not close every window:\n" + traceback.format_exc())


def run_all():
    ctx = Context()
    ctx.repo = find_repo()
    ctx.launcher = os.path.join(ctx.repo, "Cell_Quant_Studio_final.py")
    ctx.windows_before = list(Window.getWindows())
    ctx.started_imagej = False
    if IJ.getInstance() is None and not os.environ.get("CELLQUANT_GUI_NO_IMAGEJ"):
        # A user always has Fiji's main window; without one, ImageJ prints
        # its messages instead of showing them.
        ImageJ(ImageJ.NO_SHOW)
        ctx.started_imagej = True
    ctx.tmp = tempfile.mkdtemp(prefix="cellquant_gui_")
    ctx.out_dir = os.path.join(ctx.tmp, "output")
    os.makedirs(ctx.out_dir)
    os.makedirs(os.path.join(ctx.tmp, "reference"))
    ctx.app = None
    ctx.manager = None
    ctx.watcher = DialogWatcher()
    ctx.watcher.start()
    saved = (Prefs.blackBackground, Prefs.padEdges, ImageConverter.getDoScaling())
    lines = []
    passed = 0
    failed = 0
    total_checks = 0
    aborted = False
    hung = False
    start = time.time()
    try:
        for name, function in STEPS:
            recorder = Recorder(name)
            if aborted:
                lines.append("NOT RUN  " + name)
                continue
            step_start = time.time()
            try:
                function(ctx, recorder)
                settle(ctx)
            except Hang:
                recorder.failures.append("the window stopped responding:\n" + traceback.format_exc())
                aborted = True
                hung = True
            except ANY_ERROR:
                recorder.failures.append("unexpected error:\n" + traceback.format_exc())
            if ctx.app is None:
                # Without the window no later step can say anything useful.
                aborted = True
            leftover = ctx.watcher.take()
            for record in leftover:
                recorder.failures.append(
                    "unexpected " + record["kind"] + " (" + record["level"] + "): " + record["title"] + " | " + record["text"]
                )
            total_checks += recorder.checks
            status = "PASS"
            if len(recorder.failures) > 0:
                status = "FAIL"
                failed += 1
            else:
                passed += 1
            line = status.ljust(8) + " " + name + "  [" + str(recorder.checks) + " checks, " + ("%.1f" % (time.time() - step_start)) + " s]"
            lines.append(line)
            print(line)
            for failure in recorder.failures:
                lines.append("         - " + failure)
                print("         - " + failure)
    finally:
        ctx.watcher.stop()
        if failed > 0:
            try:
                log_text = IJ.getLog()
                if log_text:
                    print("Last lines of the ImageJ Log window:")
                    print("\n".join(as_text(log_text).split("\n")[-60:]))
            except ANY_ERROR:
                pass
        if not hung:
            close_everything(ctx)
        shutil.rmtree(ctx.tmp, True)
    if (Prefs.blackBackground, Prefs.padEdges, ImageConverter.getDoScaling()) != saved:
        failed += 1
        lines.append("FAIL     Fiji options were not restored after the test")
    if os.path.isdir(ctx.tmp):
        lines.append("NOTE     the temporary folder could not be removed: " + ctx.tmp)
    lines.append("")
    lines.append(
        "Window test: " + str(passed) + " passed, " + str(failed) + " failed of " + str(len(STEPS)) + " steps; " +
        str(total_checks) + " checks; " + ("%.0f" % (time.time() - start)) + " s."
    )
    print("")
    print("\n".join(lines))
    return failed


CELLQUANT_GUI_FAILURES = run_all()
if CELLQUANT_GUI_FAILURES > 0:
    raise RuntimeError("Cell Quant Studio window test failed: " + str(CELLQUANT_GUI_FAILURES) + " step(s). See the report above.")

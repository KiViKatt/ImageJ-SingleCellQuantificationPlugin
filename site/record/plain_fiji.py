# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach
# Copyright (C) 2026 Maksym Taran

# Cell Quant Studio Universal - the same job in plain Fiji, for the landing page
#
# Opens Fiji's own tools for the job, with no Cell Quant Studio: the image,
# the Channels Tool, the Threshold window, ROI Manager and the Analyze
# Particles dialog (stage 1), then the Results and Summary tables (stage 2).
# The windows are only opened and placed; the pictures of them are taken from
# outside, window by window, and arranged for the "before" side of the page.
# See site/README.md.
#
#   <Fiji launcher> --run "<absolute path>/site/record/plain_fiji.py"
#
# Environment variables:
#   CQ_DEMO_IMAGE   the image to open (make_demo_image.py builds one)
#   CQ_FLAGS        a folder for the hand-over: this script writes
#                   stage1.ready and stage2.ready there, and waits for the
#                   files stage2.go and stop.go to appear
#
# The source is kept plain ASCII with no encoding declaration on purpose.

import os
import threading
import time

from ij import IJ, WindowManager
from ij.gui import OvalRoi
from ij.plugin.frame import RoiManager
from java.awt import Window, Dialog
from java.lang import Throwable as JavaThrowable

ERRORS = (Exception, JavaThrowable)
IMAGE = os.environ["CQ_DEMO_IMAGE"]
FLAGS = os.environ["CQ_FLAGS"]
LEFT = 60
TOP = 60


def flag(name):
    handle = open(os.path.join(FLAGS, name), "w")
    handle.write("ok\n")
    handle.close()


def wait_for(name, seconds=900):
    end = time.time() + seconds
    while time.time() < end and not os.path.isfile(os.path.join(FLAGS, name)):
        time.sleep(0.3)


def window_titled(start):
    for window in Window.getWindows():
        try:
            if window.isShowing() and window.getTitle() is not None and window.getTitle().startswith(start):
                return window
        except ERRORS:
            pass
    return None


def place(start, x, y):
    for attempt in range(40):
        window = window_titled(start)
        if window is not None:
            window.setLocation(LEFT + x, TOP + y)
            return window
        time.sleep(0.25)
    print("window not found: " + start)
    return None


IJ.run("Bio-Formats Importer", "open=[" + IMAGE + "] autoscale color_mode=Composite view=Hyperstack stack_order=XYCZT")
imp = IJ.getImage()
IJ.getInstance().setLocation(LEFT, TOP)
IJ.run(imp, "Set... ", "zoom=12.5 x=2042 y=2042")
imp.getWindow().setLocation(LEFT, TOP + 110)
imp.setC(3)
# one channel in gray, as you look at it when you set a threshold by eye
imp.setDisplayMode(IJ.GRAYSCALE)
imp.updateAndDraw()

IJ.run("Channels Tool...")
place("Channels", 560, 110)
IJ.run("Threshold...")
IJ.setAutoThreshold(imp, "Otsu dark")
place("Threshold", 560, 330)

manager = RoiManager.getRoiManager()
manager.reset()
center = imp.getWidth() / 2.0
for name, reach in [("zone 0-100 um", 1750), ("zone 100-200 um", 1270), ("zone over 200 um", 790)]:
    roi = OvalRoi(center - reach, center - reach, 2 * reach, 2 * reach)
    roi.setName(name)
    manager.addRoi(roi)
manager.setLocation(LEFT + 880, TOP + 110)
manager.select(0)


def particles_dialog():
    try:
        IJ.run(imp, "Analyze Particles...", "")
    except ERRORS:
        pass


threading.Thread(target=particles_dialog).start()
place("Analyze Particles", 880, 400)
time.sleep(1.0)
flag("stage1.ready")
print("STAGE1-READY")

wait_for("stage2.go")
dialog = window_titled("Analyze Particles")
if dialog is not None:
    dialog.dispose()
time.sleep(0.8)
for title in ["Threshold", "Channels"]:
    extra = window_titled(title)
    if extra is not None:
        extra.dispose()
imp.setC(3)
IJ.setAutoThreshold(imp, "Otsu dark")
manager.select(0)
IJ.run(imp, "Analyze Particles...", "size=8-200 display clear summarize")
place("Results", 560, 110)
place("Summary", 560, 560)
time.sleep(1.0)
flag("stage2.ready")
print("STAGE2-READY")

wait_for("stop.go")
print("STOPPING")

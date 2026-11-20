# Cell Quant Studio Universal - application window
#
# Swing user interface for the analysis engine in cellquantengine.py:
# channel configuration, preparation and mask preview, single-image analysis,
# folder batches, validation against manual annotation, and settings profiles.
#
# The launcher (Cell_Quant_Studio_final.py) loads cellquantengine.py and then
# this file into one namespace, so everything the engine defines is available
# here. This file is not meant to be started on its own.
#
# Compatibility target: Fiji with ImageJ 1.53+ and Jython 2.7.
# The source is kept plain ASCII with no encoding declaration on purpose.

CORE_VERSION = "2.3.0"

if "AnalysisEngine" not in globals():
    raise RuntimeError(
        "Cell Quant Studio: cellquantengine.py has not been loaded. Start the application with "
        "Cell_Quant_Studio_final.py from the same package as cellquantcore.py."
    )
if APP_VERSION != CORE_VERSION:
    raise RuntimeError(
        "Cell Quant Studio: cellquantcore.py is version " + CORE_VERSION + " but cellquantengine.py is version " +
        APP_VERSION + ". Re-extract the release so that all files come from one package."
    )

from ij import IJ, WindowManager, Prefs
from ij.process import ImageProcessor
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable, Measurements

from javax.swing import (
    JFrame, JPanel, JTabbedPane, JScrollPane, JTextArea, JButton, JLabel,
    JFileChooser, JOptionPane, JSplitPane, JTextField, JCheckBox, JComboBox,
    JProgressBar, BoxLayout, BorderFactory, SwingUtilities, UIManager, Box
)
from javax.swing import Timer as SwingTimer
from javax.swing import JViewport, SwingConstants
from javax.swing.plaf.basic import BasicHTML
from javax.swing.text import View
from java.awt import LayoutManager, GraphicsEnvironment, RadialGradientPaint, AlphaComposite
from java.awt.event import MouseAdapter
from java.awt.geom import Ellipse2D, Path2D
from java.awt import (
    BorderLayout, GridBagLayout, GridBagConstraints, Insets, Dimension, Font,
    Color, FlowLayout, Cursor, Window, RenderingHints, GradientPaint, BasicStroke
)
from java.awt.event import WindowAdapter, ActionListener
from java.awt.geom import RoundRectangle2D, Arc2D
from java.awt.image import BufferedImage
from ij.process import ImageStatistics, FloatProcessor, ColorProcessor, Blitter
from java.util import HashMap
from java.awt.dnd import DropTarget, DnDConstants, DropTargetAdapter
from java.awt.datatransfer import DataFlavor
from java.io import File
from java.lang import Runnable

import os
import math
import time
import traceback
import threading


# -----------------------------------------------------------------------------
# User-interface helpers
# -----------------------------------------------------------------------------

def text_value(component):
    return to_text(component.getText()).strip()


def set_combo(combo, value):
    model = combo.getModel()
    for i in range(model.getSize()):
        item = to_text(model.getElementAt(i))
        if item == to_text(value):
            combo.setSelectedIndex(i)
            return


def selection_text(roi):
    """A short description of a selection, to tell whether it was changed."""
    if roi is None:
        return ""
    try:
        bounds = roi.getBounds()
        return to_text(roi.getTypeAsString()) + repr((bounds.x, bounds.y, bounds.width, bounds.height, roi.size()))
    except ANY_ERROR:
        return "unknown"


def parse_float_text(component, label, minimum=None, allow_negative_one=False):
    return parse_float_value(text_value(component), label, minimum, allow_negative_one)


def parse_int_text(component, label, minimum=None):
    return parse_int_value(text_value(component), label, minimum)


class UiTask(Runnable):
    def __init__(self, function):
        self.function = function

    def run(self):
        self.function()


def run_on_ui(function):
    """Run function on the Swing event thread. Background work (the batch)
    reports progress through this."""
    if SwingUtilities.isEventDispatchThread():
        function()
    else:
        SwingUtilities.invokeLater(UiTask(function))


def call_on_ui(function):
    """Run function on the Swing event thread and return what it returns.

    A caller that is already on the event thread runs it directly. A caller on
    another thread (a script, or the launcher) hands it to the event thread and
    waits until it has finished; an error raised by function is raised again in
    the caller. The batch worker must never use this, because the event thread
    may be showing a message box for as long as the user likes: the batch
    reports through run_on_ui, which does not wait."""
    if SwingUtilities.isEventDispatchThread():
        return function()
    box = {}
    done = threading.Event()

    def task():
        try:
            box["value"] = function()
        except ANY_ERROR as exc:
            box["error"] = exc
        finally:
            done.set()

    SwingUtilities.invokeLater(UiTask(task))
    while not done.isSet():
        done.wait(0.5)
    if "error" in box:
        raise box["error"]
    return box.get("value")


def ui_action(method):
    """Mark a window method as one that changes Swing components.

    Swing components may only be changed on the event thread. A button click
    already arrives there. When the same method is called from a script
    thread, it is moved to the event thread and the caller waits for it,
    instead of two threads changing the window at once (which can block both
    for good)."""
    def on_event_thread(self, *args, **kwargs):
        if SwingUtilities.isEventDispatchThread():
            return method(self, *args, **kwargs)
        return call_on_ui(lambda: method(self, *args, **kwargs))

    on_event_thread.__name__ = method.__name__
    on_event_thread.__doc__ = method.__doc__
    return on_event_thread


def make_note(text):
    note = JTextArea(text)
    note.setEditable(False)
    note.setLineWrap(True)
    note.setWrapStyleWord(True)
    note.setOpaque(False)
    note.setBorder(BorderFactory.createEmptyBorder(4, 6, 4, 6))
    return note


# -----------------------------------------------------------------------------
# Drag and drop
# -----------------------------------------------------------------------------

class FileDropTarget(DropTargetAdapter):
    def __init__(self, app):
        self.app = app

    def drop(self, event):
        try:
            event.acceptDrop(DnDConstants.ACTION_COPY)
            transferable = event.getTransferable()
            if transferable.isDataFlavorSupported(DataFlavor.javaFileListFlavor):
                files = transferable.getTransferData(DataFlavor.javaFileListFlavor)
                if files is not None and len(files) > 0:
                    self.app.open_image_path(files[0].getAbsolutePath())
                    event.dropComplete(True)
                    return
            event.dropComplete(False)
        except ANY_ERROR as exc:
            self.app.report_exception("File drop failed", exc)
            event.dropComplete(False)


# -----------------------------------------------------------------------------
# Per-channel settings panel
# -----------------------------------------------------------------------------

class ChannelSettingsPanel(JPanel):
    def __init__(self, app, channel_number):
        JPanel.__init__(self, BorderLayout(6, 6))
        self.app = app
        self.channel_number = channel_number

        self.enabled_cb = JCheckBox("Enable channel", channel_number <= 4)
        self.name_field = JTextField("Channel_" + str(channel_number), 18)
        self.analysis_combo = JComboBox(ANALYSIS_MODES)
        self.analysis_combo.setSelectedItem("Intensity only")

        self.projection_combo = JComboBox(PROJECTION_MODES)
        self.projection_combo.setSelectedItem("Max Z")
        self.z_slice_field = JTextField("1", 8)
        self.subtract_bg_cb = JCheckBox("Subtract rolling-ball background", False)
        self.bg_radius_field = JTextField("30", 8)
        self.light_background_cb = JCheckBox("Light background", False)
        self.median_field = JTextField("0", 8)
        self.gaussian_field = JTextField("0", 8)
        self.lut_combo = JComboBox(LUT_NAMES)
        self.lut_combo.setSelectedItem("Grays")

        self.threshold_combo = JComboBox(THRESHOLD_MODES)
        self.threshold_combo.setSelectedItem("Auto")
        self.auto_method_combo = JComboBox(AUTO_METHODS)
        set_combo(self.auto_method_combo, "Otsu")
        self.dark_background_cb = JCheckBox("Dark background / bright foreground", True)
        self.auto_scale_field = JTextField("1.0", 8)
        self.manual_low_field = JTextField("500", 8)
        self.manual_high_field = JTextField("-1", 8)

        # Local adaptive thresholding / segmentation-only CLAHE.
        # These settings affect only a duplicate used to make the mask.
        self.local_radius_field = JTextField("35", 8)
        self.local_parameter_1_field = JTextField("0.25", 8)
        self.local_parameter_2_field = JTextField("0.5", 8)
        self.clahe_cb = JCheckBox("Use CLAHE on segmentation copy", False)
        self.clahe_block_size_field = JTextField("127", 8)
        self.clahe_histogram_bins_field = JTextField("256", 8)
        self.clahe_max_slope_field = JTextField("3.0", 8)

        self.open_iterations_field = JTextField("0", 6)
        self.close_iterations_field = JTextField("0", 6)
        self.erode_iterations_field = JTextField("0", 6)
        self.dilate_iterations_field = JTextField("0", 6)
        self.morph_count_field = JTextField("1", 6)
        self.fill_holes_cb = JCheckBox("Fill holes", False)
        self.watershed_cb = JCheckBox("Watershed final mask", False)

        self.size_units_combo = JComboBox(SIZE_UNITS)
        self.size_units_combo.setSelectedItem("Calibrated area")
        self.min_size_field = JTextField("0", 8)
        self.max_size_field = JTextField("-1", 8)
        self.min_circularity_field = JTextField("0.0", 8)
        self.max_circularity_field = JTextField("1.0", 8)
        self.exclude_edges_cb = JCheckBox("Exclude objects touching ROI/image edges", False)
        self.include_holes_cb = JCheckBox("Include holes", False)
        self.composite_rois_cb = JCheckBox("Composite object ROIs", False)
        self.four_connected_cb = JCheckBox("Use 4-connected tracing", False)
        self.cross_measure_cb = JCheckBox("Measure these object ROIs on all enabled channels", False)

        self.build_panel()

    def add_form_row(self, panel, row, label_text, component, extra=None):
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = row
        gbc.anchor = GridBagConstraints.WEST
        gbc.insets = Insets(4, 6, 4, 6)
        panel.add(JLabel(label_text), gbc)

        gbc = GridBagConstraints()
        gbc.gridx = 1
        gbc.gridy = row
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(4, 6, 4, 6)
        panel.add(component, gbc)

        if extra is not None:
            gbc = GridBagConstraints()
            gbc.gridx = 2
            gbc.gridy = row
            gbc.anchor = GridBagConstraints.WEST
            gbc.insets = Insets(4, 6, 4, 6)
            panel.add(extra, gbc)

    def build_panel(self):
        header = JPanel(GridBagLayout())
        header.setBorder(BorderFactory.createTitledBorder("Channel identity and analysis role"))
        self.add_form_row(header, 0, "Physical channel", JLabel("C" + str(self.channel_number)))
        self.add_form_row(header, 1, "Use", self.enabled_cb)
        self.add_form_row(header, 2, "Channel name", self.name_field)
        self.add_form_row(header, 3, "Analysis mode", self.analysis_combo)

        tabs = JTabbedPane()
        tabs.add("Preprocessing", JScrollPane(self.build_preprocessing_panel()))
        tabs.add("Segmentation", JScrollPane(self.build_segmentation_panel()))
        tabs.add("Object counting", JScrollPane(self.build_object_panel()))

        actions = JPanel(FlowLayout(FlowLayout.LEFT))
        actions.setBorder(BorderFactory.createTitledBorder("Channel actions"))
        actions.add(JButton(
            "Prepare / refresh",
            actionPerformed=lambda event, ch=self.channel_number: self.app.preview_prepared_channel(ch)
        ))
        actions.add(JButton(
            "Brightness / Contrast",
            actionPerformed=lambda event, ch=self.channel_number: self.app.open_brightness_contrast(ch)
        ))
        actions.add(JButton(
            "Capture display minimum",
            actionPerformed=lambda event, ch=self.channel_number: self.app.capture_display_threshold(ch)
        ))
        actions.add(JButton(
            "Preview mask",
            actionPerformed=lambda event, ch=self.channel_number: self.app.preview_mask(ch)
        ))

        self.add(header, BorderLayout.NORTH)
        self.add(tabs, BorderLayout.CENTER)
        self.add(actions, BorderLayout.SOUTH)

    def build_preprocessing_panel(self):
        panel = JPanel(GridBagLayout())
        panel.setBorder(BorderFactory.createEmptyBorder(6, 6, 6, 6))
        self.add_form_row(panel, 0, "Z processing", self.projection_combo)
        self.add_form_row(panel, 1, "Specific Z slice", self.z_slice_field, JLabel("Used only for Specific Z"))
        self.add_form_row(panel, 2, "Background subtraction", self.subtract_bg_cb)
        self.add_form_row(panel, 3, "Rolling-ball radius", self.bg_radius_field)
        self.add_form_row(panel, 4, "Background polarity", self.light_background_cb)
        self.add_form_row(panel, 5, "Median radius", self.median_field, JLabel("0 disables"))
        self.add_form_row(panel, 6, "Gaussian sigma", self.gaussian_field, JLabel("0 disables"))
        self.add_form_row(panel, 7, "Display LUT", self.lut_combo)
        return panel

    def build_segmentation_panel(self):
        panel = JPanel(GridBagLayout())
        panel.setBorder(BorderFactory.createEmptyBorder(6, 6, 6, 6))
        self.add_form_row(panel, 0, "Threshold mode", self.threshold_combo)
        self.add_form_row(panel, 1, "Automatic method", self.auto_method_combo)
        self.add_form_row(panel, 2, "Automatic polarity", self.dark_background_cb)
        self.add_form_row(panel, 3, "Automatic threshold multiplier", self.auto_scale_field)
        self.add_form_row(panel, 4, "Manual lower threshold", self.manual_low_field)
        self.add_form_row(panel, 5, "Manual upper threshold", self.manual_high_field, JLabel("-1 = image maximum"))
        self.add_form_row(panel, 6, "Local threshold radius", self.local_radius_field, JLabel("pixels; Local Phansalkar only"))
        self.add_form_row(panel, 7, "Phansalkar parameter 1", self.local_parameter_1_field)
        self.add_form_row(panel, 8, "Phansalkar parameter 2", self.local_parameter_2_field)
        self.add_form_row(panel, 9, "Segmentation-only CLAHE", self.clahe_cb)
        self.add_form_row(panel, 10, "CLAHE block size", self.clahe_block_size_field)
        self.add_form_row(panel, 11, "CLAHE histogram bins", self.clahe_histogram_bins_field)
        self.add_form_row(panel, 12, "CLAHE maximum slope", self.clahe_max_slope_field)
        self.add_form_row(panel, 13, "Open iterations", self.open_iterations_field)
        self.add_form_row(panel, 14, "Close iterations", self.close_iterations_field)
        self.add_form_row(panel, 15, "Erode iterations", self.erode_iterations_field)
        self.add_form_row(panel, 16, "Dilate iterations", self.dilate_iterations_field)
        self.add_form_row(panel, 17, "Binary neighbor count", self.morph_count_field, JLabel("1-8; 1 = ImageJ default"))
        self.add_form_row(panel, 18, "Hole filling", self.fill_holes_cb)
        self.add_form_row(panel, 19, "Object separation", self.watershed_cb)

        note = JTextArea(
            "Display minimum uses the lower Brightness/Contrast display limit without changing pixel values. "
            "Local Phansalkar calculates a spatially adaptive threshold and is useful when dim and bright objects coexist. "
            "CLAHE is applied only to the temporary segmentation copy, never to the prepared image used for intensity measurements. "
            "Open with N iterations is N erosions followed by N dilations; Close is the reverse. "
            "The iteration counts and the neighbor count set here are used as given and are stored in the profile; "
            "Fiji's Process > Binary > Options do not affect them. "
            "The raw mask is saved after cleanup but before watershed. The final mask is a separate copy after watershed."
        )
        note.setEditable(False)
        note.setLineWrap(True)
        note.setWrapStyleWord(True)
        note.setOpaque(False)
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = 20
        gbc.gridwidth = 3
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(10, 6, 6, 6)
        panel.add(note, gbc)
        return panel

    def build_object_panel(self):
        panel = JPanel(GridBagLayout())
        panel.setBorder(BorderFactory.createEmptyBorder(6, 6, 6, 6))
        self.add_form_row(panel, 0, "Particle-size units", self.size_units_combo)
        self.add_form_row(panel, 1, "Minimum object area", self.min_size_field)
        self.add_form_row(panel, 2, "Maximum object area", self.max_size_field, JLabel("-1 = no maximum"))
        self.add_form_row(panel, 3, "Minimum circularity", self.min_circularity_field)
        self.add_form_row(panel, 4, "Maximum circularity", self.max_circularity_field)
        self.add_form_row(panel, 5, "Edge exclusion", self.exclude_edges_cb)
        self.add_form_row(panel, 6, "Hole behavior", self.include_holes_cb)
        self.add_form_row(panel, 7, "ROI topology", self.composite_rois_cb)
        self.add_form_row(panel, 8, "Connectivity", self.four_connected_cb)
        self.add_form_row(panel, 9, "Cross-channel measurement", self.cross_measure_cb)

        note = JTextArea(
            "Use Objects only or Objects + area to count segmented objects. Objects are detected inside the true outline of each parent ROI; "
            "an object crossing the outline is clipped there. Edge exclusion removes every object with a pixel on the ROI outline or the image border. "
            "Cross-channel measurement is optional. "
            "Leave it disabled for markers such as MAP2 when nuclear object ROIs are not biologically meaningful."
        )
        note.setEditable(False)
        note.setLineWrap(True)
        note.setWrapStyleWord(True)
        note.setOpaque(False)
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = 10
        gbc.gridwidth = 3
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(10, 6, 6, 6)
        panel.add(note, gbc)
        return panel

    def checkbox_widgets(self):
        return {
            "enabled": self.enabled_cb,
            "subtract_background": self.subtract_bg_cb,
            "light_background": self.light_background_cb,
            "dark_background": self.dark_background_cb,
            "clahe": self.clahe_cb,
            "fill_holes": self.fill_holes_cb,
            "watershed": self.watershed_cb,
            "exclude_edges": self.exclude_edges_cb,
            "include_holes": self.include_holes_cb,
            "composite_rois": self.composite_rois_cb,
            "four_connected": self.four_connected_cb,
            "cross_measure": self.cross_measure_cb
        }

    def combo_widgets(self):
        return {
            "analysis_mode": self.analysis_combo,
            "projection": self.projection_combo,
            "lut": self.lut_combo,
            "threshold_mode": self.threshold_combo,
            "auto_method": self.auto_method_combo,
            "size_units": self.size_units_combo
        }

    def text_widgets(self):
        return {
            "name": self.name_field,
            "z_slice": self.z_slice_field,
            "background_radius": self.bg_radius_field,
            "median_radius": self.median_field,
            "gaussian_sigma": self.gaussian_field,
            "auto_scale": self.auto_scale_field,
            "manual_low": self.manual_low_field,
            "manual_high": self.manual_high_field,
            "local_radius": self.local_radius_field,
            "local_parameter_1": self.local_parameter_1_field,
            "local_parameter_2": self.local_parameter_2_field,
            "clahe_block_size": self.clahe_block_size_field,
            "clahe_histogram_bins": self.clahe_histogram_bins_field,
            "clahe_max_slope": self.clahe_max_slope_field,
            "open_iterations": self.open_iterations_field,
            "close_iterations": self.close_iterations_field,
            "erode_iterations": self.erode_iterations_field,
            "dilate_iterations": self.dilate_iterations_field,
            "morph_count": self.morph_count_field,
            "min_size": self.min_size_field,
            "max_size": self.max_size_field,
            "min_circularity": self.min_circularity_field,
            "max_circularity": self.max_circularity_field
        }

    def get_config(self):
        values = {}
        checkboxes = self.checkbox_widgets()
        for key in checkboxes:
            values[key] = checkboxes[key].isSelected()
        combos = self.combo_widgets()
        for key in combos:
            values[key] = to_text(combos[key].getSelectedItem())
        fields = self.text_widgets()
        for key in fields:
            values[key] = text_value(fields[key])
        return normalize_channel_config(self.channel_number, values)

    def set_values(self, values):
        checkboxes = self.checkbox_widgets()
        combos = self.combo_widgets()
        fields = self.text_widgets()
        for key in values:
            if key in checkboxes:
                checkboxes[key].setSelected(parse_bool(values[key]))
            elif key in combos:
                set_combo(combos[key], values[key])
            elif key in fields:
                fields[key].setText(to_text(values[key]))

    def write_properties(self, props):
        prefix = "channel." + str(self.channel_number) + "."
        config = self.get_config()
        for key in config:
            if key == "channel":
                continue
            value = config[key]
            props.setProperty(prefix + key, bool_text(value) if isinstance(value, bool) else to_text(value))

    def read_properties(self, props):
        # Start from the defaults, so that a setting the profile does not hold
        # (one added after the profile was saved) gets its default value
        # instead of whatever the window happened to show. This is what
        # load_profile_settings does for a batch that runs without the window.
        values = dict(CHANNEL_DEFAULTS)
        values["name"] = "Channel_" + str(self.channel_number)
        loaded = channel_values_from_properties(props, self.channel_number)
        for key in loaded:
            values[key] = loaded[key]
        self.set_values(values)


# -----------------------------------------------------------------------------
# Look of the window, and the Results tab
# -----------------------------------------------------------------------------
#
# The look is the window's own: deep-water blues with one glowing accent,
# Manrope for the interface and Instrument Serif italic for the figures
# (both under the SIL Open Font License, in the fonts folder). It is built
# on FlatLaf's dark theme, which ships with Fiji. A Swing look applies to
# every window of Fiji, so it is switched on when the window opens and the
# previous look is put back when the window closes. "Fiji default" (top
# right) leaves Fiji alone; the choice is remembered.

LOOK_DARK = "Dark"
LOOK_FIJI = "Fiji default"
LOOK_NAMES = [LOOK_DARK, LOOK_FIJI]
LOOK_PREF_KEY = "cellquant.look"
LOOK_STATE = {"active": False, "name": None, "previous": None, "extras": None}
LOOK_FONTS = {"loaded": False}
LOOK_FONT_FILES = [
    ("ui", "Manrope-Regular.ttf"), ("ui_semibold", "Manrope-SemiBold.ttf"), ("ui_bold", "Manrope-Bold.ttf"),
    ("figure", "InstrumentSerif-Italic.ttf"), ("figure_upright", "InstrumentSerif-Regular.ttf")
]
# The palette.
LOOK_DEEP = Color(6, 18, 31)
LOOK_PANEL = Color(11, 34, 51)
LOOK_LINE = Color(22, 54, 77)
LOOK_TEXT = Color(221, 243, 245)
LOOK_FAINT = Color(127, 163, 181)
LOOK_ACCENT = Color(94, 242, 224)
LOOK_ACCENT_HEX = "#5ef2e0"
LOOK_CORAL = Color(255, 122, 107)
LOOK_VIOLET = Color(155, 140, 255)
LOOK_EXTRAS = [
    ("@background", "#081a29"),
    ("@foreground", "#ddf3f5"),
    ("@accentColor", LOOK_ACCENT_HEX),
    ("Label.disabledForeground", "#7fa3b5"),
    ("Component.borderColor", "#1b4460"),
    ("Component.arc", "14"),
    ("TextComponent.arc", "14"),
    ("Button.arc", "999"),
    ("Button.background", "#0f2b40"),
    ("Button.hoverBackground", "#143750"),
    ("Button.pressedBackground", "#0b2233"),
    ("Button.borderColor", "#1b4460"),
    ("TextField.background", "#0b2233"),
    ("TextArea.background", "#0b2233"),
    ("ComboBox.background", "#0b2233"),
    ("ComboBox.buttonBackground", "#0b2233"),
    ("TabbedPane.underlineColor", LOOK_ACCENT_HEX),
    ("TabbedPane.tabHeight", "36"),
    ("TabbedPane.showTabSeparators", "false"),
    ("ProgressBar.arc", "999"),
    ("ProgressBar.foreground", LOOK_ACCENT_HEX),
    ("ProgressBar.background", "#16364d"),
    ("ScrollBar.width", "12"),
    ("ScrollBar.thumbArc", "999"),
    ("ScrollBar.thumbInsets", "2,2,2,2"),
    ("SplitPane.dividerSize", "10"),
    ("TitledBorder.titleColor", "#7fa3b5"),
    ("TitledBorder.position", "ABOVE_TOP"),
    ("TitledBorder.border", "8,10,10,10,#16364d,1,16")
]
CHANNEL_SWATCHES = {
    "Blue": Color(120, 180, 255), "Green": Color(94, 242, 160), "Red": Color(255, 122, 107),
    "Magenta": Color(232, 140, 255), "Cyan": Color(94, 242, 224), "Yellow": Color(255, 214, 110)
}
CHANNEL_SERIES = [Color(94, 242, 224), Color(155, 140, 255), Color(255, 122, 107), Color(255, 214, 110)]
# One color per ring of the organoid map, from the outside in.
RING_COLORS = [
    Color(94, 242, 224), Color(155, 140, 255), Color(255, 122, 107),
    Color(255, 214, 110), Color(120, 180, 255), Color(232, 140, 255)
]
RING_LIMIT = 8
RESULTS_TAB_TITLE = "Results"
RESULTS_PICTURE_SIZE = 900
CARD_ARC = 16.0


def look_preference():
    try:
        look_value = to_text(Prefs.get(LOOK_PREF_KEY, LOOK_DARK))
    except ANY_ERROR:
        return LOOK_DARK
    return look_value if look_value in LOOK_NAMES else LOOK_DARK


def look_fonts_load():
    """Read the fonts in the fonts folder beside this file, once. A font
    that is missing is simply not used."""
    if LOOK_FONTS["loaded"]:
        return
    LOOK_FONTS["loaded"] = True
    folders = []
    if globals().get("package_dir"):
        folders.append(os.path.join(to_text(globals().get("package_dir")), "fonts"))
    try:
        folders.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts"))
    except ANY_ERROR:
        pass
    for folder in folders:
        if not os.path.isdir(folder):
            continue
        for pair in LOOK_FONT_FILES:
            path = os.path.join(folder, pair[1])
            if pair[0] in LOOK_FONTS or not os.path.isfile(path):
                continue
            try:
                font = Font.createFont(Font.TRUETYPE_FONT, File(path))
                GraphicsEnvironment.getLocalGraphicsEnvironment().registerFont(font)
                LOOK_FONTS[pair[0]] = font
            except ANY_ERROR:
                pass


def look_font(kind, size):
    """A font of the look: "ui", "ui_semibold", "ui_bold", "figure" (the
    serif italic of the large numbers) or "figure_upright". Without the
    bundled font, the nearest system font."""
    look_fonts_load()
    font = LOOK_FONTS.get(kind)
    if font is not None:
        return font.deriveFont(float(size))
    if kind == "figure":
        return Font("Serif", Font.ITALIC, int(size))
    if kind == "figure_upright":
        return Font("Serif", Font.PLAIN, int(size))
    base = UIManager.getFont("Label.font")
    if base is None:
        base = Font("SansSerif", Font.PLAIN, 13)
    return base.deriveFont(Font.PLAIN if kind == "ui" else Font.BOLD, float(size))


def look_install(name=LOOK_DARK):
    """Switch Fiji to the window's look (event thread). True when it is on."""
    if LOOK_STATE["active"] and LOOK_STATE["name"] == name:
        return True
    try:
        from com.formdev.flatlaf import FlatLaf, FlatDarkLaf
        look_fonts_load()
        if not LOOK_STATE["active"]:
            # The look object itself is kept, to be handed back by look_restore().
            LOOK_STATE["previous"] = UIManager.getLookAndFeel()
            LOOK_STATE["extras"] = FlatLaf.getGlobalExtraDefaults()
        extras = HashMap()
        for pair in LOOK_EXTRAS:
            extras.put(pair[0], pair[1])
        if "ui" in LOOK_FONTS:
            extras.put("defaultFont", "13 Manrope")
        FlatLaf.setGlobalExtraDefaults(extras)
        FlatDarkLaf.setup()
        FlatLaf.updateUI()
        LOOK_STATE["active"] = True
        LOOK_STATE["name"] = name
    except ANY_ERROR as exc:
        # No FlatLaf in this Fiji: the window keeps Fiji's look.
        try:
            IJ.log(APP_NAME + ": the window's own look is not available (" + to_text(exc) + ").")
        except ANY_ERROR:
            pass
        return False
    return True


def look_restore():
    """Put back the look Fiji had before look_install() (event thread)."""
    if not LOOK_STATE["active"]:
        return
    LOOK_STATE["active"] = False
    LOOK_STATE["name"] = None
    try:
        from com.formdev.flatlaf import FlatLaf
        FlatLaf.setGlobalExtraDefaults(LOOK_STATE["extras"])
        if LOOK_STATE["previous"] is not None:
            UIManager.setLookAndFeel(LOOK_STATE["previous"])
        FlatLaf.updateUI()
    except ANY_ERROR:
        pass


def look_color(key, fallback):
    try:
        found = UIManager.getColor(key)
        if found is not None:
            return found
    except ANY_ERROR:
        pass
    return fallback


def look_is_dark():
    back = look_color("Panel.background", Color.white)
    return (back.getRed() + back.getGreen() + back.getBlue()) < 384


def look_mix(first, second, share):
    """A color between first and second: share 0 is first, 1 is second."""
    return Color(
        int(round(first.getRed() + (second.getRed() - first.getRed()) * share)),
        int(round(first.getGreen() + (second.getGreen() - first.getGreen()) * share)),
        int(round(first.getBlue() + (second.getBlue() - first.getBlue()) * share))
    )


def look_smooth(graphics):
    graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON)
    graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON)
    graphics.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR)


def style_button(button, kind="plain"):
    """One filled, glowing button per section ("primary"), the rest plain.
    Looks other than FlatLaf ignore this."""
    try:
        if kind == "primary":
            button.putClientProperty(
                "FlatLaf.style",
                "background: #5ef2e0; hoverBackground: #8af7ea; pressedBackground: #3fd6c4; "
                "foreground: #06121f; borderWidth: 0; font: bold; margin: 5,14,5,14"
            )
        else:
            button.putClientProperty("FlatLaf.style", "margin: 5,14,5,14")
    except ANY_ERROR:
        pass
    return button


def paint_fitted(graphics, text, left, baseline, room):
    """Draw text that never runs past room pixels: it is shortened with
    three dots when it is too long. Returns the width drawn."""
    text = to_text(text)
    metrics = graphics.getFontMetrics()
    if room <= 0:
        return 0
    if metrics.stringWidth(text) > room:
        while len(text) > 1 and metrics.stringWidth(text + "...") > room:
            text = text[:-1]
        text = text + "..."
        if metrics.stringWidth(text) > room:
            return 0
    graphics.drawString(text, int(left), int(baseline))
    return metrics.stringWidth(text)


def channel_swatch(config, position):
    """Color that stands for a channel in the Results tab: its lookup table,
    or one of a fixed series when it has none."""
    if config is not None and config.get("lut") in CHANNEL_SWATCHES:
        return CHANNEL_SWATCHES[config["lut"]]
    return CHANNEL_SERIES[position % len(CHANNEL_SERIES)]


def format_count(value):
    """1234567 -> '1,234,567'."""
    digits = str(int(round(value)))
    sign = ""
    if digits.startswith("-"):
        sign = "-"
        digits = digits[1:]
    groups = []
    while len(digits) > 3:
        groups.insert(0, digits[-3:])
        digits = digits[:-3]
    groups.insert(0, digits)
    return sign + ",".join(groups)


def result_number(row, key):
    """A table value as a float, or None when it is missing or not a number."""
    try:
        number = float(row.get(key, ""))
    except (TypeError, ValueError):
        return None
    if number != number:
        return None
    return number


# The superscript two of "per mm2", built without a non-ASCII character in the source.
SQUARED = b"\xc2\xb2".decode("utf-8")
METRIC_ORDER = {"count": 0, "density": 1, "percent": 2}


def format_metric(value, unit):
    """A value of the Results tab as text: a count, a density per square
    millimeter, or a percentage."""
    if value is None:
        return "-"
    if unit == "percent":
        return "%.1f%%" % value
    if unit == "density":
        return format_count(value) + " /mm" + SQUARED
    return format_count(value)


def results_from_summary(headers, rows, configs):
    """What the Results tab shows, worked out from the ROI summary table.

    Returns a dictionary:
      "tiles"    the headline figures: {"title", "value", "kind" ("count" or
                 "percent"), "note", "color"}. Counts are added up over the
                 regions; a positive-area fraction is the positive area of
                 all regions divided by their area.
      "metrics"  what the rings can show, one value per region each:
                 {"name", "short", "unit" ("count", "density" or "percent"),
                 "color", "values"}; object counts first, then object
                 densities per square millimeter, then positive-area shares.
      "labels"   the name of every region (summary row).
      "areas"    the area of every region as text, with its unit.
      "bars"     the first metric in the older shape {"labels", "series",
                 "unit"}, or None.
      "regions"  the number of summary rows."""
    by_channel = {}
    for config in configs or []:
        by_channel[config["channel"]] = config
    tiles = []
    metrics = []
    position = 0
    for header in headers:
        kind = None
        if header.endswith("_Object_Count") and not header.endswith("_Edge_Touching_Object_Count") \
                and not header.endswith("_Excluded_Edge_Object_Count"):
            kind = "count"
            prefix = header[:-len("_Object_Count")]
        elif header.endswith("_Positive_Area_Fraction_Percent"):
            kind = "percent"
            prefix = header[:-len("_Positive_Area_Fraction_Percent")]
        if kind is None:
            continue
        config = None
        name = prefix
        if prefix.startswith("C") and "_" in prefix:
            try:
                config = by_channel.get(int(prefix[1:prefix.index("_")]))
            except ValueError:
                config = None
            name = prefix[prefix.index("_") + 1:].replace("_", " ")
        if config is not None:
            name = config["name"]
        color = channel_swatch(config, position)
        position += 1
        per_region = []
        for row in rows:
            per_region.append(result_number(row, header))
        present = [number for number in per_region if number is not None]
        if len(present) == 0:
            continue
        if kind == "count":
            note = "counted in " + str(len(present)) + (" region" if len(present) == 1 else " regions")
            densities = [result_number(row, prefix + "_Object_Density_Per_mm2") for row in rows]
            # The area of each region in square millimeters follows from its
            # count and its density; the overall density is all objects over
            # all of that area.
            area_mm2 = 0.0
            complete = len(rows) > 0
            for index in range(len(rows)):
                if per_region[index] is None or densities[index] is None:
                    complete = False
                elif per_region[index] > 0 and densities[index] > 0:
                    area_mm2 += per_region[index] / densities[index]
                elif per_region[index] > 0:
                    complete = False
            extra = ""
            if complete and area_mm2 > 0:
                extra = format_metric(sum(present) / area_mm2, "density") + " overall"
            tiles.append({
                "title": name + " objects", "value": sum(present), "kind": "count", "color": color,
                "note": note, "extra": extra
            })
            metrics.append({"name": name + " objects", "short": name, "unit": "count", "color": color, "values": per_region})
            if len([number for number in densities if number is not None]) > 0:
                metrics.append({
                    "name": name + " density", "short": name, "unit": "density", "color": color, "values": densities
                })
        else:
            positive = 0.0
            area = 0.0
            for row in rows:
                part = result_number(row, prefix + "_Positive_Area")
                whole = result_number(row, "ROI_Area")
                if part is not None and whole is not None and whole > 0:
                    positive += part
                    area += whole
            value = positive / area * 100.0 if area > 0 else sum(present) / float(len(present))
            tiles.append({
                "title": name + " positive area", "value": value, "kind": "percent", "color": color,
                "note": "of the measured area"
            })
            metrics.append({"name": name + " positive area", "short": name, "unit": "percent", "color": color, "values": per_region})
    # Counts first: they are what a viewer looks for.
    tiles = [tile for tile in tiles if tile["kind"] == "count"] + [tile for tile in tiles if tile["kind"] != "count"]
    metrics = sorted(metrics, key=lambda metric: METRIC_ORDER[metric["unit"]])
    labels = []
    areas = []
    several_frames = len(set([to_text(other.get("Frame", "")) for other in rows])) > 1
    for row in rows:
        label = to_text(row.get("Parent_ROI_Name", ""))
        if several_frames:
            label = label + " (frame " + to_text(row.get("Frame", "")) + ")"
        labels.append(label)
        size = result_number(row, "ROI_Area")
        areas.append("" if size is None else format_count(size) + " " + to_text(row.get("ROI_Area_Unit", "")))
    bars = None
    if len(metrics) > 0 and len(rows) > 0:
        first = metrics[0]
        bars = {
            "labels": labels, "unit": first["unit"],
            "series": [{"name": first["short"], "color": first["color"], "values": first["values"]}]
        }
    return {"tiles": tiles, "bars": bars, "regions": len(rows), "metrics": metrics, "labels": labels, "areas": areas}


def results_pictures(run_dir, configs=None):
    """The picture of a run for the Results tab, as a list holding at most
    one (caption, java.awt.Image, parents) entry: all channels merged in
    color with the detected objects outlined. parents are the parent ROIs
    as (name, x points, y points) in picture coordinates. Never raises."""
    merged = results_merge_picture(run_dir, configs)
    if merged is None:
        return []
    return [merged]


MERGE_COLORS = {
    "Blue": (70, 120, 255), "Green": (40, 255, 110), "Red": (255, 60, 60), "Magenta": (255, 60, 255),
    "Cyan": (40, 230, 255), "Yellow": (255, 225, 40), "Grays": (255, 255, 255), "Original": (255, 255, 255)
}
MERGE_SERIES = [(70, 120, 255), (40, 255, 110), (255, 60, 60), (255, 60, 255)]
RESULTS_MAX_OUTLINES = 40000


def picture_small_size(width, height, longest):
    factor = min(1.0, longest / float(max(width, height)))
    return (max(1, int(round(width * factor))), max(1, int(round(height * factor))), factor)


def picture_unit_range(processor):
    """A float copy of an image plane scaled so that its minimum is 0 and its
    99.8th percentile is 1 (brighter pixels stay above 1 and are clipped later)."""
    plane = processor.convertToFloat()
    stats = ImageStatistics.getStatistics(plane, Measurements.MIN_MAX, None)
    high = stats.max
    if stats.max > stats.min:
        wanted = 0.998 * stats.pixelCount
        seen = 0
        counts = stats.getHistogram()
        for index in range(len(counts)):
            seen += int(counts[index])
            if seen >= wanted:
                high = stats.histMin + (index + 1) * stats.binSize
                break
    span = high - stats.min
    if span <= 0:
        span = 1.0
    plane.subtract(stats.min)
    plane.multiply(1.0 / span)
    return plane


def results_merge_picture(run_dir, configs):
    """All prepared channels of a run in one color picture, each in the color
    of its lookup table, with the parent ROIs and the detected objects
    outlined. Returns (caption, java.awt.Image) or None. Never raises, and
    reads only files the run has written."""
    try:
        run_dir = to_text(run_dir)
        folder = os.path.join(run_dir, "Prepared")
        if not os.path.isdir(folder):
            return None
        by_channel = {}
        for config in configs or []:
            by_channel[config["channel"]] = config
        sums = None
        size = None
        full = None
        used = 0
        parent_outlines = []
        for name in sorted(os.listdir(folder)):
            text = to_text(name)
            if not text.endswith("_PREP.tif"):
                continue
            channel = None
            for number in range(1, MAX_CHANNELS + 1):
                if ("_C" + str(number) + "_") in text:
                    channel = number
            imp = IJ.openImage(os.path.join(folder, text))
            if imp is None:
                continue
            try:
                processor = imp.getProcessor()
                if processor.getBitDepth() == 24:
                    continue
                if size is None:
                    full = (processor.getWidth(), processor.getHeight())
                    size = picture_small_size(full[0], full[1], RESULTS_PICTURE_SIZE)
                    sums = [FloatProcessor(size[0], size[1]) for part in range(3)]
                elif (processor.getWidth(), processor.getHeight()) != full:
                    continue
                if (size[0], size[1]) != full:
                    processor.setInterpolationMethod(ImageProcessor.BILINEAR)
                    processor = processor.resize(size[0], size[1], True)
                plane = picture_unit_range(processor)
                config = by_channel.get(channel)
                color = MERGE_SERIES[used % len(MERGE_SERIES)]
                if config is not None and config.get("lut") in MERGE_COLORS:
                    color = MERGE_COLORS[config["lut"]]
                for part in range(3):
                    if color[part] > 0:
                        tinted = plane.duplicate()
                        tinted.multiply(float(color[part]))
                        sums[part].copyBits(tinted, 0, 0, Blitter.ADD)
                used += 1
            finally:
                safe_close_image(imp)
        if used == 0:
            return None
        parts = []
        for part in range(3):
            sums[part].setMinAndMax(0.0, 255.0)
            parts.append(sums[part].convertToByte(True).getPixels())
        merged = ColorProcessor(size[0], size[1])
        merged.setRGB(parts[0], parts[1], parts[2])
        image = merged.getBufferedImage()
        graphics = image.createGraphics()
        try:
            look_smooth(graphics)
            factor = size[2]
            roi_folder = os.path.join(run_dir, "ROIs")
            drawn = 0
            if os.path.isdir(roi_folder):
                graphics.setStroke(BasicStroke(1.0))
                graphics.setColor(Color(255, 255, 255, 200))
                for name in sorted(os.listdir(roi_folder)):
                    if not to_text(name).endswith("_Objects.zip") or drawn >= RESULTS_MAX_OUTLINES:
                        continue
                    for roi in load_rois_from_file(os.path.join(roi_folder, name)):
                        if drawn >= RESULTS_MAX_OUTLINES:
                            break
                        outline = roi.getPolygon()
                        graphics.drawPolygon(
                            [int(round(value * factor)) for value in outline.xpoints[:outline.npoints]],
                            [int(round(value * factor)) for value in outline.ypoints[:outline.npoints]],
                            outline.npoints
                        )
                        drawn += 1
                parents = os.path.join(roi_folder, "Parent_ROIs_Input.zip")
                if os.path.isfile(parents):
                    for roi in load_rois_from_file(parents):
                        outline = roi.getPolygon()
                        parent_outlines.append((
                            roi_name(roi, ""),
                            [value * factor for value in outline.xpoints[:outline.npoints]],
                            [value * factor for value in outline.ypoints[:outline.npoints]]
                        ))
        finally:
            graphics.dispose()
        caption = "All channels" + (" - " + format_count(drawn) + " detected objects outlined" if drawn > 0 else "")
        return (caption, image, parent_outlines)
    except ANY_ERROR:
        return None


class PreviewPanel(JPanel):
    """A small picture of the source image on the Session tab."""

    def __init__(self):
        JPanel.__init__(self)
        self.picture = None
        self.setPreferredSize(Dimension(260, 200))

    def show_image(self, imp):
        """Event thread: take a small copy of what the image window shows."""
        self.picture = None
        try:
            if imp is not None:
                shown = imp.getImage()
                size = picture_small_size(imp.getWidth(), imp.getHeight(), 640)
                small = BufferedImage(size[0], size[1], BufferedImage.TYPE_INT_RGB)
                graphics = small.createGraphics()
                try:
                    look_smooth(graphics)
                    graphics.drawImage(shown, 0, 0, size[0], size[1], None)
                finally:
                    graphics.dispose()
                self.picture = small
        except ANY_ERROR:
            self.picture = None
        self.repaint()

    def paintComponent(self, graphics):
        self.super__paintComponent(graphics)
        look_smooth(graphics)
        width = self.getWidth()
        height = self.getHeight()
        back = look_color("Panel.background", Color(245, 245, 245))
        text = look_color("Label.foreground", Color.black)
        if self.picture is None:
            graphics.setColor(look_mix(text, back, 0.6))
            graphics.setFont(self.getFont().deriveFont(Font.PLAIN, 13.0))
            words = "The image appears here"
            graphics.drawString(words, (width - graphics.getFontMetrics().stringWidth(words)) // 2, height // 2)
            return
        room_width = width - 16
        room_height = height - 16
        if room_width < 10 or room_height < 10:
            return
        scale = min(room_width / float(self.picture.getWidth()), room_height / float(self.picture.getHeight()))
        shown_width = int(self.picture.getWidth() * scale)
        shown_height = int(self.picture.getHeight() * scale)
        left = (width - shown_width) // 2
        top = (height - shown_height) // 2
        old_clip = graphics.getClip()
        graphics.clip(RoundRectangle2D.Float(left, top, shown_width, shown_height, 14.0, 14.0))
        graphics.drawImage(self.picture, left, top, shown_width, shown_height, None)
        graphics.setClip(old_clip)


class DropZonePanel(JPanel):
    """The place to drop an image file on the Session tab."""

    STEPS = [
        "Open an image and apply a preset",
        "Prepare the enabled channels",
        "Set thresholds; preview the masks",
        "Load or draw the regions to measure",
        "Run the analysis"
    ]

    def paintComponent(self, graphics):
        self.super__paintComponent(graphics)
        look_smooth(graphics)
        width = self.getWidth()
        height = self.getHeight()
        back = look_color("Panel.background", Color(245, 245, 245))
        text = look_color("Label.foreground", Color.black)
        faint = look_mix(text, back, 0.45)
        graphics.setColor(look_mix(back, text, 0.04))
        graphics.fill(RoundRectangle2D.Float(4.0, 8.0, width - 12.0, height - 12.0, CARD_ARC, CARD_ARC))
        graphics.setColor(look_mix(back, text, 0.28))
        graphics.setStroke(BasicStroke(1.5, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 1.0, [6.0, 6.0], 0.0))
        graphics.draw(RoundRectangle2D.Float(4.0, 8.0, width - 12.0, height - 12.0, CARD_ARC, CARD_ARC))
        graphics.clipRect(14, 14, max(0, width - 32), max(0, height - 24))
        center = width // 2
        base = self.getFont()
        if height < 170:
            # A low zone: the headline alone, in the middle. Never a cut-off part.
            graphics.setFont(base.deriveFont(Font.BOLD, 15.0))
            graphics.setColor(text)
            headline = "Drop a microscopy image here"
            shown_width = min(graphics.getFontMetrics().stringWidth(headline), width - 40)
            paint_fitted(graphics, headline, center - shown_width // 2, height // 2 + 6, width - 40)
            return
        block = 136 + 22 * len(self.STEPS)
        if height < block + 44:
            block = 110
        top = max(22, (height - block) // 2)
        # An arrow into a tray.
        graphics.setStroke(BasicStroke(2.5, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND))
        graphics.setColor(LOOK_ACCENT)
        graphics.drawLine(center, top, center, top + 34)
        graphics.drawLine(center - 12, top + 22, center, top + 34)
        graphics.drawLine(center + 12, top + 22, center, top + 34)
        graphics.drawPolyline([center - 26, center - 26, center + 26, center + 26], [top + 36, top + 48, top + 48, top + 36], 4)
        base = self.getFont()
        graphics.setFont(base.deriveFont(Font.BOLD, 17.0))
        graphics.setColor(text)
        headline = "Drop a microscopy image here"
        shown_width = min(graphics.getFontMetrics().stringWidth(headline), width - 40)
        paint_fitted(graphics, headline, center - shown_width // 2, top + 84, width - 40)
        graphics.setFont(base.deriveFont(Font.PLAIN, 13.0))
        graphics.setColor(faint)
        hint = "or use Browse, or the image that is open in Fiji"
        shown_width = min(graphics.getFontMetrics().stringWidth(hint), width - 40)
        paint_fitted(graphics, hint, center - shown_width // 2, top + 106, width - 40)
        metrics = graphics.getFontMetrics()
        widest = 0
        for step in self.STEPS:
            widest = max(widest, metrics.stringWidth(step))
        widest = min(widest, width - 80)
        left = center - (widest + 30) // 2
        for index in range(len(self.STEPS)):
            line = top + 144 + 22 * index
            if top + 144 + 22 * (len(self.STEPS) - 1) > height - 22:
                # Not room for the whole list: none of it, never a part.
                break
            graphics.setColor(look_mix(back, text, 0.12))
            graphics.fillOval(left, line - 14, 20, 20)
            graphics.setColor(text)
            number = str(index + 1)
            graphics.drawString(number, left + 10 - metrics.stringWidth(number) // 2, line)
            graphics.setColor(faint)
            paint_fitted(graphics, self.STEPS[index], left + 30, line, widest)


class ResultsCanvas(JPanel):
    """The organoid map: the picture of the run in a round window, one ring
    per parent ROI around it (its arc is that region's share of the largest
    value), the headline figures on the left and the regions on the right.

    Pointing at a ring or at a region lights it up, with its outline on the
    picture. The row of names above the map chooses what the rings show.
    A click on the picture shows all of it, uncropped; a click closes it."""

    def __init__(self):
        JPanel.__init__(self)
        self.data = None
        self.pictures = []
        self.progress = 1.0
        self.hover = None
        self.metric = 0
        self.zoom = False
        self.setOpaque(True)
        watcher = RingWatcher(self)
        self.addMouseMotionListener(watcher)
        self.addMouseListener(watcher)

    # ----- what is shown -----

    def metrics(self):
        if self.data is None:
            return []
        return self.data.get("metrics", [])

    def current_metric(self):
        found = self.metrics()
        if len(found) == 0:
            return None
        if self.metric >= len(found):
            self.metric = 0
        return found[self.metric]

    def ring_values(self):
        """(labels, values, unit) of the rings, at most RING_LIMIT of them."""
        chosen = self.current_metric()
        if chosen is None:
            return ([], [], "count")
        return (self.data.get("labels", [])[:RING_LIMIT], chosen["values"][:RING_LIMIT], chosen["unit"])

    # ----- geometry, shared by painting and pointing -----

    def geometry(self):
        width = self.getWidth()
        height = self.getHeight()
        labels, values, unit = self.ring_values()
        count = len(labels)
        side = 270 if width >= 1100 else 210
        show_right = width >= 820
        map_left = side
        map_right = width - (side if show_right else 0)
        # The names that choose what the rings show, in one row above the map.
        pills = []
        found = self.metrics()
        pill_space = 0
        if len(found) > 1:
            measure = self.getFontMetrics(look_font("ui_semibold", 12))
            widths = [measure.stringWidth(metric["name"]) + 26 for metric in found]
            total = sum(widths) + 6 * (len(widths) - 1)
            if total <= map_right - map_left - 20:
                start = (map_left + map_right) / 2.0 - total / 2.0
                for index in range(len(widths)):
                    pills.append((start, 14.0, float(widths[index]), 28.0))
                    start += widths[index] + 6
                pill_space = 46
        hint_space = 26
        center_x = (map_left + map_right) / 2.0
        center_y = pill_space + (height - pill_space - hint_space) / 2.0
        outer = max(40.0, min(map_right - map_left, height - pill_space - hint_space) / 2.0 - 22.0)
        band = 16.0
        gap = 10.0
        if outer - count * (band + gap) < 80:
            band = 9.0
            gap = 6.0
        rings = []
        for index in range(count):
            rings.append((outer - index * (band + gap) - band / 2.0, band))
        picture = max(20.0, outer - count * (band + gap) - 6.0)
        rows = []
        if show_right:
            row_high = 60
            row_top = max(64.0, height / 2.0 - (count * row_high) / 2.0)
            for index in range(count):
                if row_top + (index + 1) * row_high > height - 30:
                    break
                rows.append((map_right + 6, row_top + index * row_high, side - 30, row_high))
        return {
            "width": width, "height": height, "side": side, "show_right": show_right,
            "map_left": map_left, "map_right": map_right, "x": center_x, "y": center_y,
            "rings": rings, "band": band, "picture": picture, "rows": rows, "pills": pills,
            "outer": outer
        }

    def ring_at(self, x, y):
        """Index of the ring or region row under a point, or None."""
        shape = self.geometry()
        distance = math.sqrt((x - shape["x"]) ** 2 + (y - shape["y"]) ** 2)
        for index in range(len(shape["rings"])):
            radius, band = shape["rings"][index]
            if abs(distance - radius) <= band / 2.0 + 4:
                return index
        for index in range(len(shape["rows"])):
            left, top, wide, high = shape["rows"][index]
            if left <= x <= left + wide and top <= y <= top + high:
                return index
        return None

    def pill_at(self, x, y):
        shape = self.geometry()
        for index in range(len(shape["pills"])):
            left, top, wide, high = shape["pills"][index]
            if left <= x <= left + wide and top <= y <= top + high:
                return index
        return None

    def on_picture(self, x, y):
        shape = self.geometry()
        return len(self.pictures) > 0 and math.sqrt((x - shape["x"]) ** 2 + (y - shape["y"]) ** 2) <= shape["picture"]

    def point_at(self, x, y):
        found = None if self.zoom else self.ring_at(x, y)
        clickable = self.zoom or self.pill_at(x, y) is not None or self.on_picture(x, y)
        try:
            self.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR if clickable else Cursor.DEFAULT_CURSOR))
        except ANY_ERROR:
            pass
        if found != self.hover:
            self.hover = found
            self.repaint()

    def click_at(self, x, y):
        if self.zoom:
            self.zoom = False
        else:
            pill = self.pill_at(x, y)
            if pill is not None:
                self.metric = pill
            elif self.on_picture(x, y):
                self.zoom = True
                self.hover = None
        self.repaint()

    # ----- painting -----

    def colors(self):
        own = LOOK_STATE["active"]
        back = LOOK_DEEP if own else look_color("Panel.background", Color(245, 245, 245))
        text = LOOK_TEXT if own else look_color("Label.foreground", Color.black)
        faint = LOOK_FAINT if own else look_mix(text, back, 0.45)
        line = LOOK_LINE if own else look_mix(back, text, 0.18)
        return (own, back, text, faint, line)

    def paintComponent(self, graphics):
        self.super__paintComponent(graphics)
        look_smooth(graphics)
        shape = self.geometry()
        width = shape["width"]
        height = shape["height"]
        own, back, text, faint, line = self.colors()
        if own:
            # Deep water: lighter around the map, dark at the edges.
            reach = float(max(width, height)) * 0.75
            graphics.setPaint(RadialGradientPaint(
                float(shape["x"]), float(shape["y"]), max(1.0, reach), [0.0, 1.0], [LOOK_PANEL, LOOK_DEEP]
            ))
            graphics.fillRect(0, 0, width, height)
        if self.data is None or (len(self.data["tiles"]) == 0 and len(self.pictures) == 0):
            self.paint_empty(graphics, width, height, text, faint)
            return
        self.paint_figures(graphics, shape, text, faint)
        self.paint_map(graphics, shape, back, text, faint, line)
        self.paint_pills(graphics, shape, back, text, faint, line)
        if shape["show_right"]:
            self.paint_regions(graphics, shape, text, faint, line)
        if len(self.pictures) > 0 and not self.zoom:
            graphics.setFont(look_font("ui", 12))
            graphics.setColor(faint)
            hint = "Click the picture to see all of it."
            room = shape["map_right"] - shape["map_left"] - 20
            shown_width = min(graphics.getFontMetrics().stringWidth(hint), room)
            paint_fitted(graphics, hint, shape["x"] - shown_width / 2.0, height - 10, room)
        if self.zoom:
            self.paint_zoom(graphics, width, height, text, faint)

    def paint_empty(self, graphics, width, height, text, faint):
        graphics.setFont(look_font("figure", 34))
        graphics.setColor(text)
        headline = "No results yet"
        graphics.drawString(headline, (width - graphics.getFontMetrics().stringWidth(headline)) // 2, height // 2 - 6)
        graphics.setFont(look_font("ui", 13))
        graphics.setColor(faint)
        hint = "Open an image on the Session tab and press '3. Run analysis'."
        shown_width = min(graphics.getFontMetrics().stringWidth(hint), width - 40)
        paint_fitted(graphics, hint, (width - shown_width) // 2, height // 2 + 22, width - 40)

    def paint_figures(self, graphics, shape, text, faint):
        """Left: one large figure per measurement."""
        room = shape["side"] - 44
        tiles = self.data["tiles"]
        step = 124
        fit = max(1, int((shape["height"] - 40) // step))
        tiles = tiles[:fit]
        top = max(28.0, shape["height"] / 2.0 - len(tiles) * step / 2.0)
        for index in range(len(tiles)):
            tile = tiles[index]
            base = top + index * step
            graphics.setColor(tile["color"])
            graphics.fillOval(26, int(base + 3), 8, 8)
            graphics.setFont(look_font("ui_semibold", 11))
            graphics.setColor(faint)
            paint_fitted(graphics, tile["title"].upper(), 42, base + 12, room - 16)
            shown = tile["value"] * self.progress
            unit = "percent" if tile["kind"] == "percent" else "count"
            figure = format_metric(shown, unit)
            final = format_metric(tile["value"], unit)
            size = 54
            graphics.setFont(look_font("figure", size))
            while size > 24 and graphics.getFontMetrics().stringWidth(final) > room:
                size -= 4
                graphics.setFont(look_font("figure", size))
            graphics.setColor(text)
            graphics.drawString(figure, 26, int(base + 62))
            graphics.setFont(look_font("ui", 12))
            graphics.setColor(faint)
            paint_fitted(graphics, tile["note"], 26, base + 84, room)
            if tile.get("extra", "") != "":
                paint_fitted(graphics, tile["extra"], 26, base + 101, room)

    def paint_pills(self, graphics, shape, back, text, faint, line):
        """Above the map: what the rings show. The chosen one is filled."""
        found = self.metrics()
        graphics.setFont(look_font("ui_semibold", 12))
        metrics = graphics.getFontMetrics()
        for index in range(len(shape["pills"])):
            left, top, wide, high = shape["pills"][index]
            outline = RoundRectangle2D.Float(left, top, wide, high, high, high)
            if index == self.metric:
                graphics.setColor(LOOK_ACCENT)
                graphics.fill(outline)
                graphics.setColor(LOOK_DEEP)
            else:
                graphics.setStroke(BasicStroke(1.0))
                graphics.setColor(line)
                graphics.draw(outline)
                graphics.setColor(faint)
            name = found[index]["name"]
            graphics.drawString(name, int(left + (wide - metrics.stringWidth(name)) / 2.0), int(top + 19))

    def paint_map(self, graphics, shape, back, text, faint, line):
        labels, values, unit = self.ring_values()
        highest = 0.0
        for number in values:
            if number is not None and number > highest:
                highest = number
        if unit == "percent":
            highest = 100.0
        if highest <= 0:
            highest = 1.0
        center_x = shape["x"]
        center_y = shape["y"]
        for index in range(len(shape["rings"])):
            radius, band = shape["rings"][index]
            color = RING_COLORS[index % len(RING_COLORS)]
            dimmed = self.hover is not None and self.hover != index
            graphics.setStroke(BasicStroke(float(band)))
            graphics.setColor(line)
            graphics.draw(Ellipse2D.Float(center_x - radius, center_y - radius, 2 * radius, 2 * radius))
            number = values[index]
            if number is None or number <= 0:
                continue
            sweep = -360.0 * min(1.0, number / highest) * self.progress
            arc = Arc2D.Float(center_x - radius, center_y - radius, 2 * radius, 2 * radius, 90.0, sweep, Arc2D.OPEN)
            if not dimmed:
                # The glow: the same arc, wider and faint, under the arc itself.
                for halo in [(band + 14.0, 26), (band + 7.0, 46)]:
                    graphics.setStroke(BasicStroke(halo[0], BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND))
                    graphics.setColor(Color(color.getRed(), color.getGreen(), color.getBlue(), halo[1]))
                    graphics.draw(arc)
            graphics.setStroke(BasicStroke(band - 4.0, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND))
            graphics.setColor(look_mix(color, back, 0.62) if dimmed else color)
            graphics.draw(arc)
        self.paint_picture(graphics, shape, text, faint, line)

    def paint_outlines(self, graphics, parents, left, top, scale, named):
        """The parent ROIs on a picture drawn at (left, top) with a scale."""
        for index in range(len(parents)):
            name, xs, ys = parents[index]
            if len(xs) < 3:
                continue
            color = RING_COLORS[index % len(RING_COLORS)]
            path = Path2D.Float()
            path.moveTo(left + xs[0] * scale, top + ys[0] * scale)
            for corner in range(1, len(xs)):
                path.lineTo(left + xs[corner] * scale, top + ys[corner] * scale)
            path.closePath()
            if named:
                graphics.setStroke(BasicStroke(2.0))
                graphics.setColor(color)
                graphics.draw(path)
                graphics.setFont(look_font("ui_semibold", 13))
                paint_fitted(graphics, name, left + min(xs) * scale + 8, top + min(ys) * scale + 20,
                             (max(xs) - min(xs)) * scale - 16)
            elif self.hover == index:
                graphics.setColor(Color(color.getRed(), color.getGreen(), color.getBlue(), 46))
                graphics.fill(path)
                graphics.setStroke(BasicStroke(2.5))
                graphics.setColor(color)
                graphics.draw(path)
            elif self.hover is None:
                graphics.setStroke(BasicStroke(1.0))
                graphics.setColor(Color(221, 243, 245, 90))
                graphics.draw(path)

    def paint_picture(self, graphics, shape, text, faint, line):
        """The picture of the run in the round window at the center."""
        center_x = shape["x"]
        center_y = shape["y"]
        radius = shape["picture"]
        window = Ellipse2D.Float(center_x - radius, center_y - radius, 2 * radius, 2 * radius)
        if len(self.pictures) == 0:
            graphics.setColor(line)
            graphics.setStroke(BasicStroke(1.5))
            graphics.draw(window)
            return
        picture = self.pictures[0]
        image = picture[1]
        parents = picture[2] if len(picture) > 2 else []
        # Cover the round window: the shorter side of the picture fills it.
        scale = 2 * radius / float(min(image.getWidth(), image.getHeight()))
        left = center_x - image.getWidth() * scale / 2.0
        top = center_y - image.getHeight() * scale / 2.0
        old_clip = graphics.getClip()
        graphics.clip(window)
        graphics.drawImage(image, int(left), int(top), int(image.getWidth() * scale), int(image.getHeight() * scale), None)
        self.paint_outlines(graphics, parents, left, top, scale, False)
        graphics.setClip(old_clip)
        graphics.setStroke(BasicStroke(2.0))
        graphics.setColor(line)
        graphics.draw(window)

    def paint_zoom(self, graphics, width, height, text, faint):
        """All of the picture, uncropped, over the map; a click closes it."""
        picture = self.pictures[0]
        image = picture[1]
        parents = picture[2] if len(picture) > 2 else []
        graphics.setColor(LOOK_DEEP)
        graphics.fillRect(0, 0, width, height)
        room_width = width - 56
        room_height = height - 76
        if room_width < 40 or room_height < 40:
            return
        scale = min(room_width / float(image.getWidth()), room_height / float(image.getHeight()))
        shown_width = int(image.getWidth() * scale)
        shown_height = int(image.getHeight() * scale)
        left = (width - shown_width) // 2
        top = 22 + (room_height - shown_height) // 2
        old_clip = graphics.getClip()
        graphics.clip(RoundRectangle2D.Float(left, top, shown_width, shown_height, CARD_ARC, CARD_ARC))
        graphics.drawImage(image, left, top, shown_width, shown_height, None)
        self.paint_outlines(graphics, parents, left, top, scale, True)
        graphics.setClip(old_clip)
        graphics.setFont(look_font("ui", 12))
        graphics.setColor(faint)
        caption = to_text(picture[0]) + ".  Click anywhere to go back to the map."
        shown = min(graphics.getFontMetrics().stringWidth(caption), width - 40)
        paint_fitted(graphics, caption, (width - shown) // 2, height - 20, width - 40)

    def paint_regions(self, graphics, shape, text, faint, line):
        """Right: the regions, one row per ring. The large figure is what
        the rings show; the line under the name holds the other values."""
        labels, values, unit = self.ring_values()
        rows = shape["rows"]
        if len(rows) == 0:
            return
        found = self.metrics()
        chosen = self.current_metric()
        areas = self.data.get("areas", [])
        left = rows[0][0]
        wide = rows[0][2]
        graphics.setFont(look_font("ui_semibold", 11))
        graphics.setColor(faint)
        paint_fitted(graphics, (chosen["name"] + " by region").upper(), left + 8, rows[0][1] - 16, wide - 8)
        for index in range(len(rows)):
            row_left, row_top, row_wide, row_high = rows[index]
            color = RING_COLORS[index % len(RING_COLORS)]
            if self.hover == index:
                graphics.setColor(Color(color.getRed(), color.getGreen(), color.getBlue(), 34))
                graphics.fill(RoundRectangle2D.Float(row_left, row_top + 2, row_wide, row_high - 4, 12.0, 12.0))
            graphics.setColor(line)
            graphics.setStroke(BasicStroke(1.0))
            graphics.drawLine(int(row_left + 8), int(row_top + row_high), int(row_left + row_wide - 8), int(row_top + row_high))
            number = values[index]
            figure = format_metric(None if number is None else number * self.progress, unit)
            size = 26
            graphics.setFont(look_font("figure", size))
            while size > 16 and graphics.getFontMetrics().stringWidth(format_metric(number, unit)) > row_wide * 0.5:
                size -= 2
                graphics.setFont(look_font("figure", size))
            figure_wide = graphics.getFontMetrics().stringWidth(format_metric(number, unit))
            graphics.setColor(text)
            graphics.drawString(figure, int(row_left + row_wide - 10 - figure_wide), int(row_top + 28))
            graphics.setColor(color)
            graphics.fillOval(int(row_left + 10), int(row_top + 15), 9, 9)
            graphics.setFont(look_font("ui", 13))
            graphics.setColor(text if self.hover in (None, index) else faint)
            paint_fitted(graphics, labels[index], row_left + 28, row_top + 25, row_wide - 28 - figure_wide - 22)
            # The other values of this region, and its area.
            others = []
            for metric_index in range(len(found)):
                if metric_index == self.metric or index >= len(found[metric_index]["values"]):
                    continue
                other = found[metric_index]
                word = other["short"] + (" area" if other["unit"] == "percent" else "")
                if other["unit"] == "density":
                    word = ""
                others.append((word + " " + format_metric(other["values"][index], other["unit"])).strip())
            if index < len(areas) and areas[index] != "":
                others.append(areas[index])
            # As many of them as fit whole; none is cut short.
            graphics.setFont(look_font("ui", 11))
            graphics.setColor(faint)
            measure = graphics.getFontMetrics()
            line_text = ""
            for other_text in others:
                longer = other_text if line_text == "" else line_text + ",  " + other_text
                if measure.stringWidth(longer) > row_wide - 38:
                    break
                line_text = longer
            graphics.drawString(line_text, int(row_left + 28), int(row_top + 46))
        note_top = rows[-1][1] + rows[-1][3] + 24
        if note_top < shape["height"] - 14:
            graphics.setFont(look_font("ui", 12))
            graphics.setColor(faint)
            more = len(self.data.get("labels", [])) - len(rows)
            note = "Point at a ring to see its region." if more <= 0 else str(more) + " more regions are in the tables."
            paint_fitted(graphics, note, left + 8, note_top, wide - 8)


class RingWatcher(MouseAdapter):
    """Tells the organoid map where the pointer is."""

    def __init__(self, canvas):
        self.canvas = canvas

    def mouseMoved(self, event):
        self.canvas.point_at(event.getX(), event.getY())

    def mouseExited(self, event):
        self.canvas.point_at(-1000, -1000)

    def mouseClicked(self, event):
        self.canvas.click_at(event.getX(), event.getY())


class ResultsAnimation(ActionListener):
    """Lets the numbers and bars of the Results tab grow for half a second."""

    def __init__(self, canvas):
        self.canvas = canvas
        self.started = time.time()
        self.timer = SwingTimer(16, self)

    def actionPerformed(self, event):
        elapsed = (time.time() - self.started) / 0.6
        if elapsed >= 1.0:
            self.canvas.progress = 1.0
            self.timer.stop()
        else:
            # Ease out: fast at first, slow at the end.
            self.canvas.progress = 1.0 - (1.0 - elapsed) ** 3
        self.canvas.repaint()


class LookWindowListener(WindowAdapter):
    """Puts Fiji's own look back when the window has closed."""

    def windowClosed(self, event):
        look_restore()


# -----------------------------------------------------------------------------
# One layout for every tab
# -----------------------------------------------------------------------------
#
# The tabs are built row by row with their own grids. tidy_tab() then gives
# all of them the same shape: the content sits in one centered column of a
# fixed largest width, every label column of a tab has one width with its
# labels right-aligned (long labels wrap), buttons in a row have one size,
# and a form in a scrolling area starts at the top.

CONTENT_WIDTH = 1040
FORM_LABEL_CAP = 240
FORM_HINT_CAP = 250
BUTTON_ROW_LIMIT = 980


class ColumnLayout(LayoutManager):
    """Gives the one child the full height and at most `limit` pixels of
    width, centered."""

    def __init__(self, limit):
        self.limit = limit

    def addLayoutComponent(self, name, component):
        pass

    def removeLayoutComponent(self, component):
        pass

    def preferredLayoutSize(self, parent):
        if parent.getComponentCount() == 0:
            return Dimension(0, 0)
        wanted = parent.getComponent(0).getPreferredSize()
        return Dimension(min(wanted.width, self.limit), wanted.height)

    def minimumLayoutSize(self, parent):
        if parent.getComponentCount() == 0:
            return Dimension(0, 0)
        wanted = parent.getComponent(0).getMinimumSize()
        return Dimension(min(wanted.width, self.limit), wanted.height)

    def layoutContainer(self, parent):
        if parent.getComponentCount() == 0:
            return
        width = min(parent.getWidth(), self.limit)
        parent.getComponent(0).setBounds((parent.getWidth() - width) // 2, 0, width, parent.getHeight())


def in_column(component):
    """The component in a centered column of at most CONTENT_WIDTH pixels."""
    holder = JPanel(ColumnLayout(CONTENT_WIDTH))
    holder.add(component)
    return holder


def html_text(text):
    text = to_text(text)
    for pair in [("&", "&amp;"), ("<", "&lt;"), (">", "&gt;")]:
        text = text.replace(pair[0], pair[1])
    return text


def tidy_collect(container, grids, rows, notes):
    for child in container.getComponents():
        if isinstance(child, JPanel):
            layout = child.getLayout()
            if isinstance(layout, GridBagLayout):
                grids.append(child)
            elif isinstance(layout, FlowLayout):
                rows.append(child)
        elif isinstance(child, JTextArea) and not child.isEditable() and not child.isOpaque():
            notes.append(child)
        if hasattr(child, "getComponents"):
            tidy_collect(child, grids, rows, notes)


def tidy_plain(label):
    """Put a label back to its own text and size, so that tidy_tab can run
    again after the look (and with it the font) has changed."""
    kept = label.getClientProperty("cellquant.text")
    if kept is None:
        label.putClientProperty("cellquant.text", to_text(label.getText()))
    else:
        label.setText(kept)
    label.setPreferredSize(None)
    label.setMinimumSize(None)


def tidy_tab(tab):
    """Give the forms and button rows of one tab one consistent shape.
    It can be called again, for example after the look has changed."""
    grids = []
    rows = []
    notes = []
    tidy_collect(tab, grids, rows, notes)

    # Label columns: one width for the whole tab, right-aligned.
    labels = []
    widest = 0
    for panel in grids:
        layout = panel.getLayout()
        for child in panel.getComponents():
            constraints = layout.getConstraints(child)
            if isinstance(child, JLabel) and constraints.gridx == 0 and constraints.gridwidth == 1:
                tidy_plain(child)
                labels.append((layout, child, constraints))
                widest = max(widest, child.getPreferredSize().width)
            elif isinstance(child, JLabel) and constraints.gridx >= 2:
                # A hint behind a field: secondary text, wrapped when long so
                # that it never forces the form wider than its column.
                child.putClientProperty("FlatLaf.style", "foreground: $Label.disabledForeground")
                tidy_plain(child)
                hint_size = child.getPreferredSize()
                if hint_size.width > FORM_HINT_CAP:
                    child.setText("<html>" + html_text(child.getText()) + "</html>")
                    hint_view = child.getClientProperty(BasicHTML.propertyKey)
                    if hint_view is not None:
                        hint_view.setSize(float(FORM_HINT_CAP), 0.0)
                        child.setPreferredSize(Dimension(
                            FORM_HINT_CAP, int(math.ceil(hint_view.getPreferredSpan(View.Y_AXIS)))
                        ))
            elif isinstance(child, JButton) and constraints.gridx >= 2:
                # Buttons behind fields (Browse, Apply preset): one width per form.
                constraints.fill = GridBagConstraints.HORIZONTAL
                layout.setConstraints(child, constraints)
    # A few pixels to spare: text is measured a little differently from one
    # screen to the next, and a label must never lose its last letter.
    width = min(widest + 6, FORM_LABEL_CAP)
    for entry in labels:
        layout, label, constraints = entry
        wanted = label.getPreferredSize()
        height = wanted.height
        if wanted.width + 6 > width:
            # Too long for the column: wrap it instead of widening the column.
            label.setText("<html><div style='text-align: right'>" + html_text(label.getText()) + "</div></html>")
            view = label.getClientProperty(BasicHTML.propertyKey)
            if view is not None:
                view.setSize(float(width), 0.0)
                height = int(math.ceil(view.getPreferredSpan(View.Y_AXIS)))
        label.setPreferredSize(Dimension(width, height))
        label.setMinimumSize(Dimension(width, height))
        label.setHorizontalAlignment(SwingConstants.RIGHT)
        constraints.anchor = GridBagConstraints.EAST
        constraints.fill = GridBagConstraints.NONE
        constraints.insets = Insets(5, 8, 5, 4)
        layout.setConstraints(label, constraints)

    # A form inside a scrolling area starts at the top instead of floating
    # in the middle of the free space.
    for panel in grids:
        if isinstance(panel.getParent(), JViewport) and panel.getClientProperty("cellquant.pinned") is None:
            panel.putClientProperty("cellquant.pinned", "yes")
            filler = GridBagConstraints()
            filler.gridx = 0
            filler.gridy = 1000
            filler.weighty = 1.0
            filler.fill = GridBagConstraints.VERTICAL
            panel.add(Box.createVerticalGlue(), filler)

    # Button rows: every button of a row has the size of the largest one,
    # unless the row would then no longer fit.
    for panel in rows:
        buttons = [child for child in panel.getComponents() if isinstance(child, JButton)]
        if len(buttons) == 0:
            continue
        layout = panel.getLayout()
        layout.setHgap(8)
        for button in buttons:
            button.setPreferredSize(None)
        high = max([button.getPreferredSize().height for button in buttons])
        wide = max([button.getPreferredSize().width for button in buttons])
        even = len(buttons) >= 2 and (wide + 8) * len(buttons) <= BUTTON_ROW_LIMIT
        for button in buttons:
            button.setPreferredSize(Dimension((wide if even else button.getPreferredSize().width) + 4, high))

    # Explanations: secondary text.
    for note in notes:
        note.putClientProperty("FlatLaf.style", "foreground: $Label.disabledForeground")


# -----------------------------------------------------------------------------
# Main application
# -----------------------------------------------------------------------------

class CellQuantStudioUniversal(JFrame):
    def __init__(self):
        if IJ.versionLessThan("1.53"):
            raise RuntimeError(APP_NAME + " requires ImageJ 1.53 or newer.")
        JFrame.__init__(self, APP_NAME + " " + APP_VERSION)
        self.setSize(1440, 940)
        self.setMinimumSize(Dimension(1120, 760))
        self.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE)

        self.source_imp = None
        self.source_path = ""
        self.prepared_images = {}
        self.prepared_signature = None
        self.prepared_source_id = None
        self.prepared_frame = None
        self.mask_preview_ids = {}
        self.generated_ids = set()
        self.last_run_dir = ""
        self.run_log_messages = []
        self.channel_panels = []
        self.engine = AnalysisEngine(log=self.append_status)
        self.batch_running = False
        self.batch_cancel_requested = False
        self.batch_thread = None
        # Settings of the optional analysis modules that have no control on
        # the Analyses tab: kept as a loaded profile set them (see
        # collect_global_options).
        self.module_options = {}
        # Controls of the Analyses tab: option key -> (kind, widget).
        self.module_widgets = {}
        # Buttons of the Analyses tab, disabled while any background job runs.
        self.analysis_action_buttons = []
        # A background job other than the batch (see start_background).
        # batch_running is also True while one runs, so that every action
        # that refuses to start during a batch also refuses during it.
        self.background_active = False
        self.background_label = ""
        self.background_uses_images = False
        # ID of the source image whose unsaved changes the user has accepted.
        self.accepted_changes_id = None

        # What the Results tab shows: set by show_results().
        self.results_animation = None

        if look_preference() != LOOK_FIJI:
            look_install(look_preference())
        self.addWindowListener(LookWindowListener())
        self.build_ui()
        self.apply_preset_by_name(PRESET_NAMES[0])
        self.setVisible(True)

    # ----- UI construction ----------------------------------------------------

    def build_ui(self):
        root = JPanel(BorderLayout(8, 8))
        root.setBorder(BorderFactory.createEmptyBorder(12, 16, 10, 16))
        root.add(self.build_header(), BorderLayout.NORTH)

        self.main_tabs = JTabbedPane()
        self.main_tabs.add("Session", self.build_session_tab())
        self.main_tabs.add(RESULTS_TAB_TITLE, self.build_results_tab())
        self.main_tabs.add("Channels", in_column(self.build_channels_tab()))
        self.main_tabs.add("ROIs and output", in_column(self.build_roi_output_tab()))
        self.main_tabs.add("Batch", in_column(self.build_batch_tab()))
        self.main_tabs.add("Validation and QC", in_column(self.build_validation_tab()))
        if len(SETTINGS_FORMS) > 0:
            self.main_tabs.add("Analyses", in_column(self.build_analyses_tab()))
        self.main_tabs.add("Help", in_column(self.build_help_tab()))
        root.add(self.main_tabs, BorderLayout.CENTER)
        root.add(self.build_job_bar(), BorderLayout.SOUTH)

        self.getContentPane().add(root)
        self.style_buttons(root)
        for tab_index in range(self.main_tabs.getTabCount()):
            tidy_tab(self.main_tabs.getComponentAt(tab_index))

    def build_header(self):
        """The title row: name, one line on what the tool does, and the look."""
        header = JPanel(BorderLayout(12, 0))
        header.setBorder(BorderFactory.createEmptyBorder(6, 8, 4, 4))
        words = JPanel()
        words.setLayout(BoxLayout(words, BoxLayout.Y_AXIS))
        base = UIManager.getFont("Label.font")
        if base is None:
            base = Font("SansSerif", Font.PLAIN, 13)
        title = JLabel(APP_NAME)
        title.setFont(look_font("ui_bold", 24))
        subtitle = JLabel(
            "Count cells and measure markers in fluorescence microscopy images, the same way every time.  Version " + APP_VERSION
        )
        subtitle.setFont(base.deriveFont(Font.PLAIN, 13.0))
        # The secondary text color of the look, also after the look changes.
        subtitle.putClientProperty("FlatLaf.style", "foreground: $Label.disabledForeground")
        words.add(title)
        words.add(Box.createVerticalStrut(2))
        words.add(subtitle)
        side = JPanel(FlowLayout(FlowLayout.RIGHT, 8, 10))
        look_label = JLabel("Appearance")
        look_label.putClientProperty("FlatLaf.style", "foreground: $Label.disabledForeground")
        self.look_combo = JComboBox(LOOK_NAMES)
        self.look_combo.setSelectedItem(look_preference())
        self.look_combo.setToolTipText(
            "Dark gives the whole of Fiji this window's look while it is open; Fiji default leaves Fiji alone. The choice is remembered."
        )
        self.look_combo.addActionListener(self.look_changed)
        side.add(look_label)
        side.add(self.look_combo)
        header.add(words, BorderLayout.CENTER)
        header.add(side, BorderLayout.EAST)
        return header

    def look_changed(self, event=None):
        chosen = to_text(self.look_combo.getSelectedItem())
        try:
            Prefs.set(LOOK_PREF_KEY, chosen)
        except ANY_ERROR:
            pass
        if chosen == LOOK_FIJI:
            look_restore()
        else:
            look_install(chosen)
        # Another look means another font: shape the tabs again.
        for tab_index in range(self.main_tabs.getTabCount()):
            tidy_tab(self.main_tabs.getComponentAt(tab_index))
        self.revalidate()
        self.repaint()

    def style_buttons(self, container):
        """Round every button that has no style of its own yet."""
        for child in container.getComponents():
            if isinstance(child, JButton):
                if child.getClientProperty("FlatLaf.style") is None:
                    style_button(child)
            elif hasattr(child, "getComponents"):
                self.style_buttons(child)

    def build_results_tab(self):
        """The latest analysis as a picture: headline numbers, a chart of the
        object counts per region, and the overlay images. The tables and
        files of the run are unchanged; this only shows them."""
        panel = JPanel(BorderLayout(8, 8))
        self.results_canvas = ResultsCanvas()
        self.results_caption = JLabel("The latest analysis is shown here.")
        bar = JPanel(BorderLayout())
        bar.setBorder(BorderFactory.createEmptyBorder(8, 18, 0, 12))
        bar.add(self.results_caption, BorderLayout.CENTER)
        buttons = JPanel(FlowLayout(FlowLayout.RIGHT))
        buttons.add(JButton("Show result files", actionPerformed=self.open_last_output_folder))
        bar.add(buttons, BorderLayout.EAST)
        panel.add(bar, BorderLayout.NORTH)
        panel.add(self.results_canvas, BorderLayout.CENTER)
        return panel

    def show_results(self, result, configs, run_dir, source_title=""):
        """Event thread: fill the Results tab from a finished analysis and
        bring it to the front. A failure here never fails the analysis."""
        try:
            self.results_canvas.data = results_from_summary(result["summary_headers"], result["summary_rows"], configs)
            self.results_canvas.pictures = results_pictures(run_dir, configs)
            self.results_canvas.metric = 0
            self.results_canvas.zoom = False
            regions = self.results_canvas.data["regions"]
            self.results_caption.setText(
                to_text(source_title) + "   -   " + str(regions) + (" region" if regions == 1 else " regions") +
                "   -   " + time.strftime("%H:%M:%S")
            )
            if self.results_animation is not None:
                self.results_animation.timer.stop()
            self.results_canvas.progress = 0.0
            self.results_animation = ResultsAnimation(self.results_canvas)
            self.results_animation.timer.start()
            self.main_tabs.setSelectedIndex(self.main_tabs.indexOfTab(RESULTS_TAB_TITLE))
        except ANY_ERROR as exc:
            self.append_status("The Results tab could not be filled: " + to_text(exc))

    def build_job_bar(self):
        """The strip under the tabs that shows a background job (an analysis,
        a validation, a sensitivity, tuning, tissue, count correction or
        robustness run) and lets the user cancel it."""
        bar = JPanel(FlowLayout(FlowLayout.LEFT))
        self.job_label = JLabel("No background job is running.")
        self.job_progress = JProgressBar(0, 100)
        # A thin bar; what it stands for is written beside it (job_label).
        # The string is still set, for scripts and tests that read the state.
        self.job_progress.setStringPainted(False)
        self.job_progress.setString("Idle")
        self.job_progress.setPreferredSize(Dimension(220, 6))
        self.job_cancel_button = JButton("Cancel", actionPerformed=self.cancel_background_action)
        self.job_cancel_button.setEnabled(False)
        bar.add(self.job_progress)
        bar.add(self.job_cancel_button)
        bar.add(self.job_label)
        return bar

    def add_grid_row(self, panel, row, label_text, component, button=None):
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = row
        gbc.anchor = GridBagConstraints.WEST
        gbc.insets = Insets(4, 6, 4, 6)
        panel.add(JLabel(label_text), gbc)

        gbc = GridBagConstraints()
        gbc.gridx = 1
        gbc.gridy = row
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(4, 6, 4, 6)
        panel.add(component, gbc)

        if button is not None:
            gbc = GridBagConstraints()
            gbc.gridx = 2
            gbc.gridy = row
            gbc.insets = Insets(4, 6, 4, 6)
            panel.add(button, gbc)

    def build_session_tab(self):
        panel = JPanel(BorderLayout(8, 8))

        inputs = JPanel(GridBagLayout())
        inputs.setBorder(BorderFactory.createTitledBorder("Input and workflow"))
        self.image_path_field = JTextField("", 60)
        self.output_dir_field = JTextField("", 60)
        self.preset_combo = JComboBox(PRESET_NAMES)
        self.frame_mode_combo = JComboBox(FRAME_MODES)

        self.add_grid_row(
            inputs, 0, "Image file", self.image_path_field,
            JButton("Browse", actionPerformed=self.choose_image_file)
        )
        self.add_grid_row(
            inputs, 1, "Output folder", self.output_dir_field,
            JButton("Browse", actionPerformed=self.choose_output_dir)
        )
        self.add_grid_row(
            inputs, 2, "Analysis preset", self.preset_combo,
            JButton("Apply preset", actionPerformed=self.apply_selected_preset)
        )
        self.add_grid_row(inputs, 3, "Frames", self.frame_mode_combo)

        # Two rows: getting an image, then the three steps of an analysis.
        image_row = JPanel(FlowLayout(FlowLayout.LEFT, 8, 4))
        image_row.add(JButton("Open selected image", actionPerformed=self.open_selected_image))
        image_row.add(JButton("Open via Bio-Formats", actionPerformed=self.open_selected_with_bioformats))
        image_row.add(JButton("Use active image", actionPerformed=self.use_active_image))
        image_row.add(JButton("Reset for next image", actionPerformed=self.reset_for_next_image))
        steps_row = JPanel(FlowLayout(FlowLayout.LEFT, 8, 4))
        steps_row.add(style_button(JButton("1. Prepare enabled channels", actionPerformed=self.prepare_enabled_channels)))
        steps_row.add(style_button(JButton("2. Preview all masks", actionPerformed=self.preview_all_masks)))
        steps_row.add(style_button(JButton("3. Run analysis", actionPerformed=self.run_analysis_action), "primary"))
        steps_row.add(JButton("Quick run", actionPerformed=self.quick_run_action))
        actions = JPanel()
        actions.setLayout(BoxLayout(actions, BoxLayout.Y_AXIS))
        actions.setBorder(BorderFactory.createTitledBorder("Workflow actions"))
        actions.add(image_row)
        actions.add(steps_row)

        self.drop_area = DropZonePanel()
        self.drop_area.setPreferredSize(Dimension(520, 230))
        self.drop_area.setMinimumSize(Dimension(320, 120))
        self.drop_area.setToolTipText(
            "Recommended manual-threshold workflow: open an image and apply a preset; prepare the enabled channels; "
            "use each channel's Brightness/Contrast button; preview masks and adjust thresholds/morphology; "
            "load or draw parent ROIs; run the analysis."
        )
        DropTarget(self.drop_area, FileDropTarget(self))

        self.image_info_area = JTextArea()
        self.image_info_area.setEditable(False)
        self.image_info_area.setLineWrap(True)
        self.image_info_area.setWrapStyleWord(True)
        self.image_info_area.setText("No source image selected.")
        self.image_info_area.setBorder(BorderFactory.createEmptyBorder(8, 10, 8, 10))
        self.source_preview = PreviewPanel()
        info_scroll = JScrollPane(self.image_info_area)
        info_scroll.setBorder(BorderFactory.createEmptyBorder())
        info_scroll.setPreferredSize(Dimension(340, 200))
        source_panel = JPanel(BorderLayout(6, 6))
        source_panel.setBorder(BorderFactory.createTitledBorder("Source image"))
        source_panel.add(self.source_preview, BorderLayout.CENTER)
        source_panel.add(info_scroll, BorderLayout.EAST)

        center = JSplitPane(JSplitPane.HORIZONTAL_SPLIT, self.drop_area, source_panel)
        center.setResizeWeight(0.4)

        self.status_area = JTextArea()
        self.status_area.setEditable(False)
        self.status_area.setLineWrap(True)
        self.status_area.setWrapStyleWord(True)
        self.status_area.setFont(Font("Monospaced", Font.PLAIN, 12))
        self.status_area.setBorder(BorderFactory.createTitledBorder("Status and run log"))
        self.append_status("Ready.")

        top = JPanel(BorderLayout(6, 6))
        top.add(inputs, BorderLayout.NORTH)
        top.add(actions, BorderLayout.SOUTH)

        # The log keeps a few lines even in the smallest window.
        status_scroll = JScrollPane(self.status_area)
        status_scroll.setMinimumSize(Dimension(200, 96))
        status_scroll.setPreferredSize(Dimension(600, 150))
        source_panel.setMinimumSize(Dimension(320, 120))
        center.setMinimumSize(Dimension(400, 130))
        split = JSplitPane(JSplitPane.VERTICAL_SPLIT, center, status_scroll)
        split.setResizeWeight(0.7)

        panel.add(top, BorderLayout.NORTH)
        panel.add(split, BorderLayout.CENTER)
        return panel

    def build_channels_tab(self):
        panel = JPanel(BorderLayout(6, 6))
        note = JTextArea(
            "Configure each physical channel independently. Positive-area and object measurements are always based on that channel's own mask. "

            "Object ROIs are transferred to other channels only when the explicit cross-channel option is enabled."
        )
        note.setEditable(False)
        note.setLineWrap(True)
        note.setWrapStyleWord(True)
        note.setOpaque(False)
        note.setBorder(BorderFactory.createEmptyBorder(4, 6, 4, 6))

        self.channel_tabs = JTabbedPane()
        for channel_number in range(1, MAX_CHANNELS + 1):
            channel_panel = ChannelSettingsPanel(self, channel_number)
            self.channel_panels.append(channel_panel)
            self.channel_tabs.add("C" + str(channel_number), JScrollPane(channel_panel))

        panel.add(note, BorderLayout.NORTH)
        panel.add(self.channel_tabs, BorderLayout.CENTER)
        return panel

    def build_roi_output_tab(self):
        outer = JPanel()
        outer.setLayout(BoxLayout(outer, BoxLayout.Y_AXIS))

        roi_panel = JPanel(GridBagLayout())
        roi_panel.setBorder(BorderFactory.createTitledBorder("Parent ROI source"))
        self.roi_source_combo = JComboBox(ROI_SOURCES)
        self.roi_source_combo.setSelectedItem("ROI Manager")
        self.add_grid_row(roi_panel, 0, "Quantification region", self.roi_source_combo)

        roi_buttons = JPanel(FlowLayout(FlowLayout.LEFT))
        roi_buttons.add(JButton("Open ROI Manager", actionPerformed=self.open_roi_manager))
        roi_buttons.add(JButton("Load ROI set", actionPerformed=self.load_roi_set))
        roi_buttons.add(JButton("Add current ROI", actionPerformed=self.add_current_roi))
        roi_buttons.add(JButton("Clear ROI Manager", actionPerformed=self.clear_roi_manager))
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = 1
        gbc.gridwidth = 3
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(4, 6, 4, 6)
        roi_panel.add(roi_buttons, gbc)

        calibration_panel = JPanel(GridBagLayout())
        calibration_panel.setBorder(BorderFactory.createTitledBorder("Spatial calibration"))
        self.override_calibration_cb = JCheckBox("Override source calibration", False)
        self.pixel_width_field = JTextField("1.0", 10)
        self.pixel_height_field = JTextField("1.0", 10)
        self.spatial_unit_field = JTextField("pixel", 12)
        self.add_grid_row(calibration_panel, 0, "Use override", self.override_calibration_cb)
        self.add_grid_row(calibration_panel, 1, "Pixel width", self.pixel_width_field)
        self.add_grid_row(calibration_panel, 2, "Pixel height", self.pixel_height_field)
        self.add_grid_row(calibration_panel, 3, "Spatial unit", self.spatial_unit_field)

        output_panel = JPanel(GridBagLayout())
        output_panel.setBorder(BorderFactory.createTitledBorder("Saved outputs"))
        self.save_prepared_cb = JCheckBox("Save prepared measurement images", True)
        self.save_masks_cb = JCheckBox("Save raw and final masks", True)
        self.save_parent_rois_cb = JCheckBox("Save input parent ROIs", True)
        self.save_object_rois_cb = JCheckBox("Save detected object ROIs", True)
        self.save_overlays_cb = JCheckBox("Save object-outline overlays on their segmentation channels", True)
        self.show_summary_cb = JCheckBox("Show summary table after analysis", True)
        self.add_grid_row(output_panel, 0, "Prepared images", self.save_prepared_cb)
        self.add_grid_row(output_panel, 1, "Masks", self.save_masks_cb)
        self.add_grid_row(output_panel, 2, "Parent ROIs", self.save_parent_rois_cb)
        self.add_grid_row(output_panel, 3, "Object ROIs", self.save_object_rois_cb)
        self.add_grid_row(output_panel, 4, "Overlays", self.save_overlays_cb)
        self.add_grid_row(output_panel, 5, "Results display", self.show_summary_cb)

        lifecycle_panel = JPanel(GridBagLayout())
        lifecycle_panel.setBorder(BorderFactory.createTitledBorder("Image lifecycle"))
        self.keep_prepared_open_cb = JCheckBox("Keep prepared images open after analysis", True)
        self.keep_masks_open_cb = JCheckBox("Keep analysis masks open after analysis", True)
        self.close_source_cb = JCheckBox("Close source image after successful analysis", False)
        self.close_generated_after_cb = JCheckBox("Close all app-generated images after successful analysis", False)
        self.tile_prepared_cb = JCheckBox("Tile prepared images after preparation", False)
        self.add_grid_row(lifecycle_panel, 0, "Prepared windows", self.keep_prepared_open_cb)
        self.add_grid_row(lifecycle_panel, 1, "Mask windows", self.keep_masks_open_cb)
        self.add_grid_row(lifecycle_panel, 2, "Source", self.close_source_cb)
        self.add_grid_row(lifecycle_panel, 3, "Complete cleanup", self.close_generated_after_cb)
        self.add_grid_row(lifecycle_panel, 4, "Window layout", self.tile_prepared_cb)

        profile_panel = JPanel(FlowLayout(FlowLayout.LEFT))
        profile_panel.setBorder(BorderFactory.createTitledBorder("Profiles and output access"))
        profile_panel.add(JButton("Save settings profile", actionPerformed=self.save_profile))
        profile_panel.add(JButton("Load settings profile", actionPerformed=self.load_profile))
        profile_panel.add(JButton("Open last output folder", actionPerformed=self.open_last_output_folder))
        profile_panel.add(JButton("Close generated images", actionPerformed=self.close_generated_images_action))
        profile_panel.add(JButton("Close all ImageJ images", actionPerformed=self.close_all_images_action))

        outer.add(roi_panel)
        outer.add(calibration_panel)
        outer.add(output_panel)
        outer.add(lifecycle_panel)
        outer.add(profile_panel)

        return JScrollPane(outer)

    def build_batch_tab(self):
        panel = JPanel(BorderLayout(6, 6))
        controls = JPanel()
        controls.setLayout(BoxLayout(controls, BoxLayout.Y_AXIS))

        input_panel = JPanel(GridBagLayout())
        input_panel.setBorder(BorderFactory.createTitledBorder("Batch input"))
        self.batch_input_field = JTextField("", 60)
        self.batch_extensions_field = JTextField(DEFAULT_BATCH_EXTENSIONS, 60)
        self.batch_recursive_cb = JCheckBox("Include subfolders", False)
        self.batch_force_bioformats_cb = JCheckBox("Always open with Bio-Formats", False)
        self.batch_all_series_cb = JCheckBox("Analyze every series of a multi-series file (otherwise the first series)", False)
        self.add_grid_row(
            input_panel, 0, "Input folder", self.batch_input_field,
            JButton("Browse", actionPerformed=self.choose_batch_input_dir)
        )
        self.add_grid_row(input_panel, 1, "File extensions", self.batch_extensions_field)
        self.add_grid_row(input_panel, 2, "Subfolders", self.batch_recursive_cb)
        self.add_grid_row(input_panel, 3, "File opening", self.batch_force_bioformats_cb)
        self.add_grid_row(input_panel, 4, "Series", self.batch_all_series_cb)

        roi_panel = JPanel(GridBagLayout())
        roi_panel.setBorder(BorderFactory.createTitledBorder("Batch regions"))
        self.batch_roi_combo = JComboBox(BATCH_ROI_POLICIES)
        self.batch_roi_dir_field = JTextField("", 60)
        self.add_grid_row(roi_panel, 0, "Quantification region", self.batch_roi_combo)
        self.add_grid_row(
            roi_panel, 1, "ROI folder", self.batch_roi_dir_field,
            JButton("Browse", actionPerformed=self.choose_batch_roi_dir)
        )
        roi_buttons = JPanel(FlowLayout(FlowLayout.LEFT))
        roi_buttons.add(JButton("Save ROI Manager set for the current image", actionPerformed=self.save_rois_for_current_image))
        roi_buttons.add(JLabel("Leave the ROI folder empty to look beside each image."))
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = 2
        gbc.gridwidth = 3
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(4, 6, 4, 6)
        roi_panel.add(roi_buttons, gbc)

        actions = JPanel(FlowLayout(FlowLayout.LEFT))
        actions.setBorder(BorderFactory.createTitledBorder("Batch actions"))
        self.batch_scan_button = JButton("Scan folder", actionPerformed=self.scan_batch_folder)
        self.batch_run_button = JButton("Run batch", actionPerformed=self.run_batch_action)
        self.batch_cancel_button = JButton("Cancel batch", actionPerformed=self.cancel_batch_action)
        self.batch_cancel_button.setEnabled(False)
        self.batch_progress = JProgressBar(0, 100)
        self.batch_progress.setStringPainted(True)
        self.batch_progress.setString("Idle")
        self.batch_progress.setPreferredSize(Dimension(300, 18))
        actions.add(self.batch_scan_button)
        actions.add(self.batch_run_button)
        actions.add(self.batch_cancel_button)
        actions.add(self.batch_progress)

        note = make_note(
            "A batch applies the channel settings, calibration, and saved-output choices of this window to every image in the folder. "
            "Results go to a new CellQuant_Batch folder inside the output folder chosen on the Session tab (or inside the input folder if none is set). "
            "Each image gets its own run folder; batch_roi_summary.csv, batch_objects.csv, and batch_thresholds.csv combine all images, "
            "and batch_report.csv lists every file with its status. "
            "An ROI file is a ROI Manager set saved as <image name>_ROIs.zip (or <image name>.zip / .roi). "
            "Display-minimum thresholds cannot be used in a batch: capture them as Manual thresholds first. "
            "With 'Current frame only' a batch analyzes the first time point of each image."
        )

        controls.add(input_panel)
        controls.add(roi_panel)
        controls.add(actions)
        controls.add(note)

        self.batch_list_area = JTextArea()
        self.batch_list_area.setEditable(False)
        self.batch_list_area.setFont(Font("Monospaced", Font.PLAIN, 11))
        self.batch_list_area.setText("Choose an input folder and click Scan folder to list the images a batch would analyze.")
        scan_scroll = JScrollPane(self.batch_list_area)
        scan_scroll.setBorder(BorderFactory.createTitledBorder("Folder scan"))

        panel.add(controls, BorderLayout.NORTH)
        panel.add(scan_scroll, BorderLayout.CENTER)
        return panel

    def build_validation_tab(self):
        outer = JPanel()
        outer.setLayout(BoxLayout(outer, BoxLayout.Y_AXIS))

        channel_panel = JPanel(GridBagLayout())
        channel_panel.setBorder(BorderFactory.createTitledBorder("Channel under test"))
        self.qc_channel_combo = JComboBox(["C" + str(i) for i in range(1, MAX_CHANNELS + 1)])
        self.add_grid_row(channel_panel, 0, "Channel", self.qc_channel_combo, JLabel("Tested with its current settings, on the current frame, inside the parent ROIs."))

        validation_panel = JPanel(GridBagLayout())
        validation_panel.setBorder(BorderFactory.createTitledBorder("Validation against manual annotation"))
        self.annotation_source_combo = JComboBox(ANNOTATION_SOURCES)
        self.annotation_file_field = JTextField("", 60)
        self.add_grid_row(validation_panel, 0, "Annotation", self.annotation_source_combo)
        self.add_grid_row(
            validation_panel, 1, "Annotation file", self.annotation_file_field,
            JButton("Browse", actionPerformed=self.choose_annotation_file)
        )
        validation_buttons = JPanel(FlowLayout(FlowLayout.LEFT))
        validation_buttons.add(JButton("Run validation", actionPerformed=self.run_validation_action))
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = 2
        gbc.gridwidth = 3
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(4, 6, 4, 6)
        validation_panel.add(validation_buttons, gbc)
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = 3
        gbc.gridwidth = 3
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(4, 6, 4, 6)
        validation_panel.add(make_note(
            "Object detection: mark one point per object with the Multi-point tool (on the source or a prepared image), or load saved point ROIs. "
            "Each detected object is then scored inside every parent ROI: an object containing a mark is a true positive, an object without one is a false positive, "
            "a mark outside every object is missed, and extra marks inside one object are counted as merged. Precision, recall, F1, and count error are reported. "
            "Positive area: load a reference binary mask of the same size (nonzero = positive) to get Dice, Jaccard, precision, recall, and area difference. "
            "Results are saved in a CellQuant_Validation folder with a color-coded overlay (green = matched, orange = merged, red = false positive, magenta = missed mark)."
        ), gbc)

        sensitivity_panel = JPanel(GridBagLayout())
        sensitivity_panel.setBorder(BorderFactory.createTitledBorder("Parameter sensitivity"))
        self.sensitivity_parameter_combo = JComboBox(SENSITIVITY_PARAMETERS)
        self.sensitivity_factors_field = JTextField(DEFAULT_SENSITIVITY_FACTORS, 40)
        self.add_grid_row(sensitivity_panel, 0, "Parameter to vary", self.sensitivity_parameter_combo)
        self.add_grid_row(sensitivity_panel, 1, "Multiply by", self.sensitivity_factors_field, JLabel("comma-separated; 1.0 = the configured value"))
        sensitivity_buttons = JPanel(FlowLayout(FlowLayout.LEFT))
        sensitivity_buttons.add(JButton("Run sensitivity analysis", actionPerformed=self.run_sensitivity_action))
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = 2
        gbc.gridwidth = 3
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(4, 6, 4, 6)
        sensitivity_panel.add(sensitivity_buttons, gbc)
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = 3
        gbc.gridwidth = 3
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(4, 6, 4, 6)
        sensitivity_panel.add(make_note(
            "Re-segments the channel while one setting is scaled by each factor and records positive area and object count per parent ROI, "
            "with the percent change from the configured value. This shows how strongly a reported number depends on a threshold or filter choice. "
            "Results are saved as sensitivity.csv in a CellQuant_Sensitivity folder."
        ), gbc)

        outer.add(channel_panel)
        outer.add(validation_panel)
        outer.add(sensitivity_panel)
        # Keep the three sections at their natural height at the top of the tab.
        wrapper = JPanel(BorderLayout())
        wrapper.add(outer, BorderLayout.NORTH)
        return JScrollPane(wrapper)

    # ----- Analyses tab (settings of the optional analysis modules) ----------

    def build_analyses_tab(self):
        """One section per settings form the analysis modules registered
        (SETTINGS_FORMS), in registration order. The controls hold profile
        settings like every other control of the window; the window does no
        checking of its own, the modules' readers report bad values."""
        outer = JPanel()
        outer.setLayout(BoxLayout(outer, BoxLayout.Y_AXIS))
        outer.add(make_note(
            "Optional analyses loaded beside the engine. These settings are saved in profiles and in the "
            "settings.properties of every run, like all other settings. Analyses that are switched on here "
            "(for example the depth profile or the threshold spectrum) run inside '3. Run analysis' and batches. "
            "The buttons start a separate run in the background; the strip at the bottom of the window shows "
            "its progress and can cancel it."
        ))
        for form in SETTINGS_FORMS:
            outer.add(self.build_settings_form(form))
        wrapper = JPanel(BorderLayout())
        wrapper.add(outer, BorderLayout.NORTH)
        return JScrollPane(wrapper)

    def build_settings_form(self, form):
        panel = JPanel(GridBagLayout())
        panel.setBorder(BorderFactory.createTitledBorder(to_text(form.get("title", ""))))
        defaults = default_global_options()
        row = 0
        if to_text(form.get("note", "")).strip() != "":
            self.add_wide_row(panel, row, make_note(to_text(form["note"])))
            row += 1
        keys = []
        for field in form.get("fields", []):
            key = field[0]
            label = to_text(field[1]) if len(field) > 1 else to_text(key)
            kind = field[2] if len(field) > 2 else "text"
            choices = field[3] if len(field) > 3 else None
            if kind == "bool":
                widget = JCheckBox(label, parse_bool(defaults.get(key), False))
                self.add_grid_row(panel, row, "", widget)
            elif kind == "choice" and choices:
                widget = JComboBox([to_text(choice) for choice in choices])
                self.set_choice(widget, defaults.get(key, ""))
                self.add_grid_row(panel, row, label, widget)
            else:
                kind = "text"
                widget = JTextField(to_text(defaults.get(key, "")), 20)
                self.add_grid_row(panel, row, label, widget)
            self.module_widgets[key] = (kind, widget)
            keys.append(key)
            row += 1
        actions = self.module_action_panel(keys)
        if actions is not None:
            self.add_wide_row(panel, row, actions)
        return panel

    def add_wide_row(self, panel, row, component):
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = row
        gbc.gridwidth = 3
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(4, 6, 4, 6)
        panel.add(component, gbc)

    def set_choice(self, combo, value):
        """Select value; a value that is not offered (for example a typing
        error in a profile) is added, so that saving keeps it unchanged and
        the module reports it when it is used."""
        value = to_text(value)
        model = combo.getModel()
        for i in range(model.getSize()):
            if to_text(model.getElementAt(i)) == value:
                combo.setSelectedIndex(i)
                return
        combo.addItem(value)
        combo.setSelectedIndex(combo.getItemCount() - 1)

    def analysis_action_button(self, label, handler):
        button = JButton(label, actionPerformed=handler)
        self.analysis_action_buttons.append(button)
        return button

    def module_action_panel(self, keys):
        """The button row (and inputs that are not settings) of a form whose
        module offers a separate run, or None."""
        def has(prefix, function):
            return function in globals() and len([key for key in keys if key.startswith(prefix)]) > 0

        if has("tissue_", "detect_tissue_outline"):
            row = JPanel(FlowLayout(FlowLayout.LEFT))
            row.add(self.analysis_action_button("Detect tissue outline into ROI Manager", self.detect_tissue_action))
            row.add(JLabel("Uses the source image, its current frame, and the enabled channels."))
            return row
        if has("tune_", "run_tuning_folder"):
            box = JPanel(GridBagLayout())
            self.tune_folder_field = JTextField("", 40)
            self.add_grid_row(
                box, 0, "Annotated folder (empty = Batch input folder)", self.tune_folder_field,
                JButton("Browse", actionPerformed=lambda event: self.choose_into(self.tune_folder_field, True))
            )
            buttons = JPanel(FlowLayout(FlowLayout.LEFT))
            buttons.add(self.analysis_action_button("Run tuning...", self.run_tuning_action))
            buttons.add(JLabel("File extensions, subfolders, ROI folder and Bio-Formats come from the Batch tab."))
            self.add_wide_row(box, 1, buttons)
            return box
        if has("correct_", "run_count_correction"):
            box = JPanel(GridBagLayout())
            self.correct_source_field = JTextField("", 40)
            self.correct_table_field = JTextField("", 40)
            self.correct_channel_field = JTextField("", 8)
            self.add_grid_row(
                box, 0, "Validation runs (folder) or a table of counts (CSV)", self.correct_source_field,
                JButton("Browse", actionPerformed=lambda event: self.choose_into(self.correct_source_field, None))
            )
            self.add_grid_row(
                box, 1, "Result table to correct (optional)", self.correct_table_field,
                JButton("Browse", actionPerformed=lambda event: self.choose_into(self.correct_table_field, False))
            )
            self.add_grid_row(box, 2, "Channel (number, C-label or name)", self.correct_channel_field)
            buttons = JPanel(FlowLayout(FlowLayout.LEFT))
            buttons.add(self.analysis_action_button("Run count correction...", self.run_correction_action))
            buttons.add(JLabel("Nothing existing is changed; results go to a new folder."))
            self.add_wide_row(box, 3, buttons)
            return box
        if has("robust_", "run_robustness_from_options"):
            row = JPanel(FlowLayout(FlowLayout.LEFT))
            row.add(self.analysis_action_button("Run robustness...", self.run_robustness_action))
            row.add(JLabel("Uses the input folder, regions and file options of the Batch tab; progress and Cancel are on the Batch tab and below."))
            return row
        return None

    @ui_action
    def choose_into(self, field, directories):
        """Browse for a folder (True), a file (False) or either (None)."""
        chooser = JFileChooser()
        if directories is None:
            chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES)
        elif directories:
            chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
        else:
            chooser.setFileSelectionMode(JFileChooser.FILES_ONLY)
        if chooser.showOpenDialog(self) == JFileChooser.APPROVE_OPTION:
            field.setText(chooser.getSelectedFile().getAbsolutePath())

    def module_widget_options(self):
        values = {}
        for key in self.module_widgets:
            kind, widget = self.module_widgets[key]
            if kind == "bool":
                values[key] = widget.isSelected()
            elif kind == "choice":
                values[key] = to_text(widget.getSelectedItem())
            else:
                values[key] = text_value(widget)
        return values

    def apply_module_widget_options(self, options):
        for key in self.module_widgets:
            if key not in options:
                continue
            kind, widget = self.module_widgets[key]
            if kind == "bool":
                widget.setSelected(parse_bool(options[key], False))
            elif kind == "choice":
                self.set_choice(widget, options[key])
            else:
                widget.setText(to_text(options[key]))

    def build_help_tab(self):
        text = JTextArea()
        text.setEditable(False)
        text.setLineWrap(True)
        text.setWrapStyleWord(True)
        text.setFont(Font("SansSerif", Font.PLAIN, 13))
        text.setText(
            "CELL QUANT STUDIO UNIVERSAL " + APP_VERSION + "\n\n"
            "Analysis modes\n"
            "- Intensity only: reports area, mean, standard deviation, minimum, maximum, median, integrated density, and raw integrated density inside each parent ROI.\n"
            "- Positive area: independently thresholds the channel and reports positive area and positive-area fraction/density.\n"
            "- Objects only: thresholds the channel, applies morphology and optional watershed, then counts objects using size and circularity filters.\n"
            "- Objects + area: reports both positive-area and object-count metrics.\n\n"
            "MAP2 + DAPI preset\n"
            "The MAP2/DAPI presets configure independent marker area and nuclear object analysis. The CLAHE preset uses C3 MAP2 and C2 DAPI with local Phansalkar thresholding, segmentation-only CLAHE, hole filling, and watershed. "
            "DAPI object ROIs are not transferred onto MAP2 unless you explicitly enable cross-channel object measurement. This keeps MAP2 density independent of nuclear segmentation.\n\n"
            "Mask definitions\n"
            "The raw mask is the thresholded mask after optional open, close, erosion, dilation, and hole filling, but before watershed. "
            "The final mask is a separate copy after optional watershed. Positive-area measurements use the raw mask; object counting uses the final mask. "
            "Open with N iterations means N erosions then N dilations (Close is the reverse). Iteration and neighbor counts come from this window and the profile only.\n\n"
            "Objects and parent ROIs\n"
            "Objects are detected inside the true outline of each parent ROI, whatever its shape. An object crossing the outline is clipped there and still counted. "
            "'Exclude objects touching ROI/image edges' removes every object with a pixel on the ROI outline or the image border. "
            "objects.csv flags such objects in Touches_ROI_Edge, and roi_summary.csv reports how many were found and how many were excluded.\n\n"
            "Manual thresholds\n"
            "Prepare the channel, open Brightness/Contrast, and adjust the lower display limit. Display minimum uses that lower limit without rescaling image pixels. "
            "For a reproducible numeric threshold, click Capture display minimum; the value is copied into the Manual threshold fields. "
            "Do not click Apply in the Brightness/Contrast panel unless you intentionally want to rescale the pixel values; that operation changes the measurement image.\n\n"
            "ROI behavior\n"
            "Whole image uses the entire XY field. Current image ROI uses the source image's active area selection. ROI Manager uses every loaded area ROI as a separate parent region. "
            "Line, point, and angle ROIs are rejected because densities require an area.\n\n"
            "Frames\n"
            "Current frame supports manual display thresholds. All frames processes each time point independently using the same parent ROIs and requires Auto, Manual, or Local Phansalkar thresholds.\n\n"
            "File opening\n"
            "Use Open via Bio-Formats for proprietary microscopy formats, multiple series, or files that do not open with the expected channel, Z, or time dimensions.\n\n"
            "Batch\n"
            "The Batch tab analyzes every image of a folder with the settings of this window. Regions can be the whole image, an ROI file saved for each image "
            "(<image name>_ROIs.zip), or the current ROI Manager set. Each image gets a normal run folder; combined tables and a per-file report are written to the batch folder. "
            "A batch runs in the background and can be cancelled; finished images are kept. The same batch can run without the window through cellquant_batch_headless.py.\n\n"
            "Validation and QC\n"
            "The Validation and QC tab scores object detection against manual point marks (precision, recall, F1), compares a positive-area mask with a reference mask (Dice, Jaccard), "
            "and measures how much positive area and object count change when one setting is varied.\n\n"
            "Analyses (optional modules)\n"
            "When analysis modules sit beside the engine, the Analyses tab shows their settings, which are saved in profiles and run records. "
            "Its buttons detect the tissue outline into ROI Manager, tune settings against point marks, correct counts from validated fields, "
            "and test whether a group difference holds across a grid of settings.\n\n"
            "Background runs\n"
            "Run analysis, Quick run, validation, sensitivity and the Analyses buttons work in the background: the window stays usable, "
            "the strip at the bottom shows progress, and its Cancel button stops the run after the current step. "
            "While such a run uses the prepared images, preparing again and closing generated images wait until it ends.\n\n"
            "Outputs\n"
            "Each run creates a timestamped folder containing roi_summary.csv, objects.csv, thresholds.csv, settings.properties, run_log.txt, and optional prepared images, masks, ROI sets, and overlays. "
            "Tables are UTF-8 text. thresholds.csv records every threshold and segmentation setting actually used, together with the application and ImageJ versions."
        )
        return JScrollPane(text)

    # ----- status and error handling -----------------------------------------

    def append_status(self, message):
        text = "[" + time.strftime("%H:%M:%S") + "] " + to_text(message)
        self.run_log_messages.append(text)

        def update():
            try:
                self.status_area.append(text + "\n")
                self.status_area.setCaretPosition(self.status_area.getDocument().getLength())
                if self.background_active:
                    # The latest message of a background job, beside its progress bar.
                    shown = to_text(message)
                    self.job_label.setText(self.background_label + ": " + (shown if len(shown) <= 110 else shown[:107] + "..."))
            except ANY_ERROR:
                pass

        # While the window is still being built (not yet displayable) the text
        # area is updated directly, exactly as before. Handing these early
        # messages to the event thread would let them run while this thread is
        # showing the window, and the two can block each other.
        if self.isDisplayable():
            run_on_ui(update)
        else:
            update()
        try:
            IJ.log(APP_NAME + ": " + to_text(message))
        except ANY_ERROR:
            pass

    def report_exception(self, context, exc):
        message = context + ": " + to_text(exc)
        self.append_status("ERROR - " + message)
        try:
            # Outside an except block there is no traceback; format_exc()
            # then returns the text "None", which is not worth a log line.
            details = traceback.format_exc()
            if details is not None and "Traceback" in details:
                IJ.log(details)
        except ANY_ERROR:
            pass

        def show():
            JOptionPane.showMessageDialog(self, message, APP_NAME, JOptionPane.ERROR_MESSAGE)

        run_on_ui(show)

    def set_busy(self, busy):
        try:
            if busy:
                self.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR))
            else:
                self.setCursor(Cursor.getDefaultCursor())
        except ANY_ERROR:
            pass

    def blocked_by_batch(self):
        """True (with a message) while a batch or another background job is running."""
        if not self.batch_running:
            return False
        if self.background_active:
            text = (
                self.background_label + " is running in the background. "
                "Wait for it to finish, or click Cancel at the bottom of the window."
            )
        else:
            text = "A batch is running. Wait for it to finish or cancel it on the Batch tab."
        JOptionPane.showMessageDialog(self, text, APP_NAME, JOptionPane.WARNING_MESSAGE)
        return True

    def blocked_by_image_job(self):
        """True (with a message) while a background job reads the source or
        the prepared images, which preparing again or closing the generated
        images would change under it. During a batch these stay allowed."""
        if not (self.background_active and self.background_uses_images):
            return False
        JOptionPane.showMessageDialog(
            self,
            self.background_label + " is using the prepared images. "
            "Wait for it to finish, or click Cancel at the bottom of the window.",
            APP_NAME,
            JOptionPane.WARNING_MESSAGE
        )
        return True

    # ----- background jobs ----------------------------------------------------

    def start_background(self, label, collect, work, finished, error_context,
                         uses_images=False, batch_bar=False):
        """Run one long job without blocking the window.

        On the event thread: collect() reads every control it needs and
        returns a dictionary of plain values (or None to stop quietly; it may
        ask questions and raise errors, which are shown under error_context).
        Then work(job) runs on a daemon thread with job["values"] and
        job["engine"] (an AnalysisEngine whose cancel check is the Cancel
        button). It must not touch Swing: it reports through append_status
        and job["progress"](done, total, text), which hand over to the event
        thread with invokeLater and never wait. Finally finished(job) runs
        on the event thread with job["result"], job["error"] and
        job["cancelled"]; an error is then reported in a message box.
        Returns True when the job was started."""
        if self.blocked_by_batch():
            return False
        self.set_busy(True)
        try:
            values = collect()
        except ANY_ERROR as exc:
            self.set_busy(False)
            self.report_exception(error_context, exc)
            return False
        self.set_busy(False)
        if values is None:
            return False
        job = {
            "label": label, "values": values, "work": work, "finished": finished,
            "error_context": error_context, "batch_bar": batch_bar,
            "result": None, "error": None, "cancelled": False
        }
        job["engine"] = AnalysisEngine(log=self.append_status, is_cancelled=self.batch_is_cancelled)
        job["progress"] = lambda done, total, text: self.job_progress_update(job, done, total, text)
        self.batch_cancel_requested = False
        self.batch_running = True
        self.background_active = True
        self.background_label = label
        self.background_uses_images = uses_images
        self.set_job_controls(job, True)
        worker = threading.Thread(target=self.background_worker, args=(job,), name="cellquant-" + label.lower().replace(" ", "-"))
        worker.setDaemon(True)
        self.batch_thread = worker
        worker.start()
        return True

    def background_worker(self, job):
        try:
            job["result"] = job["work"](job)
        except AnalysisCancelled:
            job["cancelled"] = True
        except ANY_ERROR as exc:
            job["error"] = exc
            try:
                IJ.log(traceback.format_exc())
            except ANY_ERROR:
                pass
        finally:
            def finish():
                self.background_finished(job)

            run_on_ui(finish)

    @ui_action
    def background_finished(self, job):
        """Event thread, after the worker has ended. Other actions stay
        refused until finished(job) has returned (its message box included),
        so that no second job can start while this one is still being shown."""
        self.set_job_controls(job, False)
        IJ.showStatus("Cell Quant Studio ready")
        try:
            if job["cancelled"]:
                self.append_status(job["label"] + " cancelled by user.")
            try:
                job["finished"](job)
            except ANY_ERROR as exc:
                if job["error"] is None:
                    job["error"] = exc
            if job["error"] is not None:
                self.report_exception(job["error_context"], job["error"])
        finally:
            self.batch_running = False
            self.batch_thread = None
            self.background_active = False
            self.background_uses_images = False

    def image_watch(self, source, configs, prepared):
        """Event thread, when a job that reads the open images starts: what it
        takes for granted about them. The window stays usable during the job,
        so changed_during_job() checks afterwards that nobody changed them."""
        display = [config["channel"] for config in configs
                   if needs_mask(config) and config["threshold_mode"] == "Display minimum"]
        watch = {
            "source": source, "source_changes": bool(source.changes), "z": source.getZ(),
            "source_open": WindowManager.getImage(source.getID()) is not None,
            "selection": selection_text(source.getRoi()),
            "current_z": len([config for config in configs if config["projection"] == "Current Z"]) > 0,
            "prepared": []
        }
        for channel in sorted((prepared or {}).keys()):
            imp = prepared[channel]
            if WindowManager.getImage(imp.getID()) is not None:
                watch["prepared"].append((channel, imp, bool(imp.changes), imp.getDisplayRangeMin(), channel in display))
        return watch

    def changed_during_job(self, job):
        """Event thread, when the job has ended: a warning text, or "" when
        the images it read were left alone while it ran."""
        watch = job["values"].get("watch")
        if watch is None:
            return ""
        notes = []
        source = watch["source"]
        if WindowManager.getImage(source.getID()) is None:
            if watch["source_open"]:
                notes.append("the source image was closed")
        else:
            if source.changes and not watch["source_changes"]:
                notes.append("the pixels of the source image were edited")
            if watch["current_z"] and source.getZ() != watch["z"]:
                notes.append("the Z slider of the source image was moved (a 'Current Z' channel uses the slice shown)")
            if selection_text(source.getRoi()) != watch["selection"]:
                notes.append("the selection on the source image was changed")
        for channel, imp, changes, low, display in watch["prepared"]:
            if WindowManager.getImage(imp.getID()) is None:
                notes.append("the prepared image of C" + str(channel) + " was closed")
            elif imp.changes and not changes:
                notes.append("the pixels of the prepared image of C" + str(channel) + " were edited")
            elif display and imp.getDisplayRangeMin() != low:
                notes.append("the display range of the prepared image of C" + str(channel) + " was moved")
        if len(notes) == 0:
            return ""
        text = (
            "While the " + job["label"].lower() + " ran, " + "; ".join(notes) + ". "
            "Its results may come from the images before or after that change: run it again without touching them."
        )
        self.append_status("Warning: " + text)
        return text

    def mark_job_cancelled(self, job):
        """For a module run that stops by itself on Cancel and returns a
        result instead of raising AnalysisCancelled."""
        job["cancelled"] = True
        self.job_label.setText(job["label"] + " cancelled.")
        self.job_progress.setValue(0)
        self.job_progress.setString("Cancelled")
        if job["batch_bar"]:
            self.batch_progress.setString("Cancelled")
        self.append_status(job["label"] + " cancelled by user.")

    def set_job_controls(self, job, running):
        if running:
            self.job_label.setText(job["label"] + " is running.")
            self.job_progress.setIndeterminate(True)
            self.job_progress.setString(job["label"])
        else:
            state = "cancelled" if job["cancelled"] else ("failed" if job["error"] is not None else "finished")
            self.job_label.setText(job["label"] + " " + state + ".")
            self.job_progress.setIndeterminate(False)
            self.job_progress.setValue(self.job_progress.getMaximum() if state == "finished" else 0)
            self.job_progress.setString(state[0].upper() + state[1:])
        self.job_cancel_button.setEnabled(running)
        for button in self.analysis_action_buttons:
            button.setEnabled(not running)
        if job["batch_bar"]:
            self.batch_run_button.setEnabled(not running)
            self.batch_cancel_button.setEnabled(running)
            if running:
                self.batch_progress_update(0, 1, "starting " + job["label"].lower())
            else:
                self.batch_progress.setString(state[0].upper() + state[1:])

    def job_progress_update(self, job, done, total, text):
        """Progress of a background job; safe to call from the worker."""
        def update():
            try:
                self.job_progress.setIndeterminate(False)
                self.job_progress.setMaximum(max(1, int(total)))
                self.job_progress.setValue(int(done))
                self.job_progress.setString(str(int(done)) + " / " + str(int(total)))
                if text:
                    self.job_label.setText(job["label"] + ": " + to_text(text))
            except ANY_ERROR:
                pass

        run_on_ui(update)
        if job["batch_bar"]:
            self.batch_progress_update(done, total, text)

    @ui_action
    def cancel_background_action(self, event=None):
        if not self.background_active:
            return
        self.cancel_batch_action()

    # ----- presets ------------------------------------------------------------

    def reset_channel_defaults(self):
        for panel in self.channel_panels:
            values = dict(CHANNEL_DEFAULTS)
            values["name"] = "Channel_" + str(panel.channel_number)
            panel.set_values(values)

    @ui_action
    def apply_selected_preset(self, event=None):
        if self.blocked_by_batch():
            return
        self.apply_preset_by_name(to_text(self.preset_combo.getSelectedItem()))

    @ui_action
    def apply_preset_by_name(self, name):
        set_combo(self.preset_combo, name)
        if name == "Custom":
            self.append_status("Custom settings retained.")
            return

        self.close_generated_images()
        self.reset_channel_defaults()

        if name == "MAP2 + DAPI (C3 marker, C2 nuclei, CLAHE)":
            self.channel_panels[2].set_values({
                "enabled": True,
                "name": "MAP2",
                "analysis_mode": "Positive area",
                "projection": "Max Z",
                "lut": "Green",
                "threshold_mode": "Display minimum",
                "fill_holes": False,
                "watershed": False,
                "cross_measure": False
            })
            self.channel_panels[1].set_values({
                "enabled": True,
                "name": "DAPI",
                "analysis_mode": "Objects + area",
                "projection": "Max Z",
                "lut": "Blue",
                "threshold_mode": "Local Phansalkar",
                "local_radius": 35,
                "local_parameter_1": 0.25,
                "local_parameter_2": 0.5,
                "clahe": True,
                "clahe_block_size": 127,
                "clahe_histogram_bins": 256,
                "clahe_max_slope": 3.0,
                "open_iterations": 0,
                "fill_holes": True,
                "watershed": True,
                "size_units": "Calibrated area",
                "min_size": 8,
                "max_size": 200,
                "cross_measure": False
            })
            self.roi_source_combo.setSelectedItem("ROI Manager")

        elif name == "MAP2 + DAPI (C3 marker, C4 nuclei)":
            self.channel_panels[2].set_values({
                "enabled": True,
                "name": "MAP2",
                "analysis_mode": "Positive area",
                "projection": "Max Z",
                "lut": "Magenta",
                "threshold_mode": "Display minimum",
                "fill_holes": False,
                "watershed": False,
                "cross_measure": False
            })
            self.channel_panels[3].set_values({
                "enabled": True,
                "name": "DAPI",
                "analysis_mode": "Objects + area",
                "projection": "Max Z",
                "lut": "Blue",
                "threshold_mode": "Display minimum",
                "fill_holes": True,
                "watershed": True,
                "size_units": "Calibrated area",
                "min_size": 8,
                "max_size": 200,
                "cross_measure": False
            })
            self.roi_source_combo.setSelectedItem("ROI Manager")

        elif name == "Nuclear objects + independent marker":
            self.channel_panels[0].set_values({
                "enabled": True,
                "name": "Marker",
                "analysis_mode": "Positive area",
                "projection": "Max Z",
                "lut": "Magenta",
                "threshold_mode": "Auto",
                "auto_method": "Otsu"
            })
            self.channel_panels[1].set_values({
                "enabled": True,
                "name": "DAPI",
                "analysis_mode": "Objects + area",
                "projection": "Max Z",
                "lut": "Blue",
                "threshold_mode": "Auto",
                "auto_method": "Otsu",
                "fill_holes": True,
                "watershed": True,
                "min_size": 8,
                "max_size": 200,
                "cross_measure": False
            })

        elif name == "Four-channel cell-based":
            colors = ["Green", "Red", "Magenta", "Blue"]
            names = ["Channel_1", "Channel_2", "Cell_Segmentation", "Channel_4"]
            for i in range(4):
                mode = "Objects + area" if i == 2 else "Intensity only"
                self.channel_panels[i].set_values({
                    "enabled": True,
                    "name": names[i],
                    "analysis_mode": mode,
                    "projection": "Max Z",
                    "lut": colors[i],
                    "threshold_mode": "Auto",
                    "auto_method": "Otsu",
                    "fill_holes": i == 2,
                    "watershed": i == 2,
                    "cross_measure": i == 2
                })

        elif name == "Positive-area quantification":
            count = self.source_imp.getNChannels() if self.source_imp is not None else 4
            count = min(count, MAX_CHANNELS)
            colors = ["Green", "Red", "Magenta", "Blue", "Cyan", "Yellow", "Grays", "Grays"]
            for i in range(count):
                self.channel_panels[i].set_values({
                    "enabled": True,
                    "name": "Marker_" + str(i + 1),
                    "analysis_mode": "Positive area",
                    "projection": "Max Z",
                    "lut": colors[i],
                    "threshold_mode": "Auto",
                    "auto_method": "Otsu"
                })

        elif name == "Intensity-only quantification":
            count = self.source_imp.getNChannels() if self.source_imp is not None else 4
            count = min(count, MAX_CHANNELS)
            for i in range(count):
                self.channel_panels[i].set_values({
                    "enabled": True,
                    "name": "Channel_" + str(i + 1),
                    "analysis_mode": "Intensity only",
                    "projection": "Max Z",
                    "lut": "Grays"
                })

        self.invalidate_prepared_cache()
        self.append_status("Applied preset: " + name)

    # ----- image selection and opening ---------------------------------------

    @ui_action
    def choose_image_file(self, event=None):
        chooser = JFileChooser()
        chooser.setFileSelectionMode(JFileChooser.FILES_ONLY)
        if chooser.showOpenDialog(self) == JFileChooser.APPROVE_OPTION:
            path = chooser.getSelectedFile().getAbsolutePath()
            self.image_path_field.setText(path)
            self.open_image_path(path)

    @ui_action
    def choose_output_dir(self, event=None):
        chooser = JFileChooser()
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
        if chooser.showOpenDialog(self) == JFileChooser.APPROVE_OPTION:
            self.output_dir_field.setText(chooser.getSelectedFile().getAbsolutePath())

    @ui_action
    def open_selected_image(self, event=None):
        path = text_value(self.image_path_field)
        if path == "":
            JOptionPane.showMessageDialog(self, "Choose an image file first.", APP_NAME, JOptionPane.WARNING_MESSAGE)
            return
        self.open_image_path(path)

    @ui_action
    def open_selected_with_bioformats(self, event=None):
        path = text_value(self.image_path_field)
        if path == "":
            JOptionPane.showMessageDialog(self, "Choose an image file first.", APP_NAME, JOptionPane.WARNING_MESSAGE)
            return
        self.open_image_path(path, force_bioformats=True)

    @ui_action
    def open_image_path(self, path, force_bioformats=False):
        if self.blocked_by_batch():
            return
        self.set_busy(True)
        try:
            if path is None or to_text(path).strip() == "":
                raise ValueError("No image path was supplied.")
            path = os.path.abspath(to_text(path))
            if not os.path.isfile(path):
                raise ValueError("Image file does not exist: " + path)

            if force_bioformats:
                imp = self.open_with_bioformats(path)
            else:
                imp = IJ.openImage(path)
                if imp is None:
                    imp = self.open_with_bioformats(path)
            if imp is None:
                raise RuntimeError("Fiji could not open the selected image.")

            imp.show()
            self.image_path_field.setText(path)
            self.set_source_image(imp, path)
            if text_value(self.output_dir_field) == "":
                self.output_dir_field.setText(os.path.dirname(path))
            self.append_status("Opened source image: " + path)
        except ANY_ERROR as exc:
            self.report_exception("Could not open image", exc)
        finally:
            self.set_busy(False)

    def open_with_bioformats(self, path):
        try:
            from loci.plugins import BF
            images = BF.openImagePlus(path)
            if images is None or len(images) == 0:
                return None
            if len(images) == 1:
                return images[0]

            combo = JComboBox([images[i].getTitle() for i in range(len(images))])
            result = JOptionPane.showConfirmDialog(
                self, combo, "Select Bio-Formats series", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE
            )
            if result != JOptionPane.OK_OPTION:
                for imp in images:
                    safe_close_image(imp)
                return None
            selected = combo.getSelectedIndex()
            chosen = images[selected]
            for i in range(len(images)):
                if i != selected:
                    safe_close_image(images[i])
            self.append_status("Bio-Formats returned " + str(len(images)) + " series; selected series " + str(selected + 1) + ".")
            return chosen
        except ANY_ERROR as exc:
            IJ.log("Bio-Formats fallback failed: " + to_text(exc))
            return None

    def use_active_image(self, event=None):
        # ImageJ keeps "the current image" per thread (a script can have its
        # own), so it is looked up on the calling thread before the rest runs
        # on the event thread.
        self.use_image_as_source(WindowManager.getCurrentImage())

    @ui_action
    def use_image_as_source(self, imp):
        if self.blocked_by_batch():
            return
        if imp is None:
            JOptionPane.showMessageDialog(self, "No active image is open.", APP_NAME, JOptionPane.WARNING_MESSAGE)
            return
        if self.is_generated_image(imp):
            # After "Prepare" the active image is usually a prepared image or
            # a mask. Taking it as the source would close it (generated images
            # are closed when the source changes) and leave no source at all.
            JOptionPane.showMessageDialog(
                self,
                "The active image (" + to_text(imp.getTitle()) + ") was made by this application "
                "(a prepared image or a mask), so it cannot be the source image.\n"
                "Click the window of the original image first, then click 'Use active image'.",
                APP_NAME,
                JOptionPane.WARNING_MESSAGE
            )
            return
        path = ""
        try:
            info = imp.getOriginalFileInfo()
            if info is not None and info.directory is not None and info.fileName is not None:
                path = os.path.join(to_text(info.directory), to_text(info.fileName))
        except ANY_ERROR:
            path = ""
        self.set_source_image(imp, path)
        if path != "":
            self.image_path_field.setText(path)
            if text_value(self.output_dir_field) == "":
                self.output_dir_field.setText(os.path.dirname(path))
        self.append_status("Using active image as source: " + to_text(imp.getTitle()))

    @ui_action
    def set_source_image(self, imp, path):
        self.close_generated_images()
        self.source_imp = imp
        self.source_path = path if path is not None else ""
        self.invalidate_prepared_cache()
        self.update_source_info()
        self.update_channel_availability()
        self.populate_calibration_fields()

    def update_source_info(self):
        self.source_preview.show_image(self.source_imp)
        if self.source_imp is None:
            self.image_info_area.setText("No source image selected.")
            return
        imp = self.source_imp
        cal = imp.getCalibration()
        unit = cal.getUnit() if cal is not None else "pixel"
        lines = [
            "Title: " + to_text(imp.getTitle()),
            "Path: " + (self.source_path if self.source_path else "<active image; no file path>"),
            "Dimensions: " + str(imp.getWidth()) + " x " + str(imp.getHeight()),
            "Channels: " + str(imp.getNChannels()),
            "Z slices: " + str(imp.getNSlices()),
            "Time frames: " + str(imp.getNFrames()),
            "Bit depth: " + str(imp.getBitDepth()),
            "Current position: C" + str(imp.getC()) + ", Z" + str(imp.getZ()) + ", T" + str(imp.getT()),
            "Calibration: " + str(cal.pixelWidth if cal is not None else 1) + " x " + str(cal.pixelHeight if cal is not None else 1) + " " + to_text(unit)
        ]
        if imp.getNChannels() > MAX_CHANNELS:
            lines.append("Warning: this application exposes the first " + str(MAX_CHANNELS) + " channels.")
        if imp.getBitDepth() == 24:
            lines.append("Warning: RGB images are not supported for threshold/object analysis. Open channels as grayscale with Bio-Formats.")
        self.image_info_area.setText("\n".join(lines))

    def update_channel_availability(self):
        count = self.source_imp.getNChannels() if self.source_imp is not None else MAX_CHANNELS
        for i in range(MAX_CHANNELS):
            enabled = i < count
            # Only the tab is switched off. The channel's settings stay as
            # they are, so that an image with fewer channels does not erase
            # them; collect_channel_configs leaves such channels out.
            self.channel_tabs.setEnabledAt(i, enabled)

    def populate_calibration_fields(self):
        if self.source_imp is None:
            return
        cal = self.source_imp.getCalibration()
        if cal is None:
            return
        self.pixel_width_field.setText(str(cal.pixelWidth))
        self.pixel_height_field.setText(str(cal.pixelHeight))
        unit = cal.getUnit()
        self.spatial_unit_field.setText(to_text(unit if unit else "pixel"))

    # ----- ROI Manager actions ------------------------------------------------

    def open_roi_manager(self, event=None):
        RoiManager.getRoiManager()

    @ui_action
    def load_roi_set(self, event=None):
        chooser = JFileChooser()
        chooser.setFileSelectionMode(JFileChooser.FILES_ONLY)
        if chooser.showOpenDialog(self) != JFileChooser.APPROVE_OPTION:
            return
        path = chooser.getSelectedFile().getAbsolutePath()
        manager = RoiManager.getRoiManager()
        if manager.open(path):
            self.roi_source_combo.setSelectedItem("ROI Manager")
            self.append_status("Loaded ROI set: " + path)
        else:
            JOptionPane.showMessageDialog(self, "Could not load ROI set:\n" + path, APP_NAME, JOptionPane.ERROR_MESSAGE)

    @ui_action
    def add_current_roi(self, event=None):
        if self.source_imp is None:
            JOptionPane.showMessageDialog(self, "Select a source image first.", APP_NAME, JOptionPane.WARNING_MESSAGE)
            return
        roi = self.source_imp.getRoi()
        if roi is None:
            JOptionPane.showMessageDialog(self, "The source image has no active ROI.", APP_NAME, JOptionPane.WARNING_MESSAGE)
            return
        manager = RoiManager.getRoiManager()
        manager.addRoi(clone_roi(roi))
        self.roi_source_combo.setSelectedItem("ROI Manager")
        self.append_status("Added current source ROI to ROI Manager.")

    @ui_action
    def clear_roi_manager(self, event=None):
        manager = RoiManager.getInstance2()
        if manager is None or manager.getCount() == 0:
            return
        result = JOptionPane.showConfirmDialog(
            self, "Clear all ROIs from ROI Manager?", APP_NAME, JOptionPane.YES_NO_OPTION
        )
        if result == JOptionPane.YES_OPTION:
            manager.reset()
            self.append_status("Cleared ROI Manager.")

    # ----- configuration collection and validation ---------------------------

    def collect_channel_configs(self, limit_to_source=True):
        # limit_to_source leaves out the channels that the open image does
        # not have. Jobs that read a folder of other images pass False.
        available = MAX_CHANNELS
        if limit_to_source and self.source_imp is not None:
            available = self.source_imp.getNChannels()
        configs = []
        for panel in self.channel_panels:
            analysis_mode = to_text(panel.analysis_combo.getSelectedItem())
            if not panel.enabled_cb.isSelected() or analysis_mode == "Ignore":
                continue
            if panel.channel_number > available:
                continue
            config = panel.get_config()
            configs.append(config)
        return configs

    def collect_all_channel_configs(self):
        return [panel.get_config() for panel in self.channel_panels]

    def collect_global_options(self):
        """Every registered setting: the window's controls (the Analyses tab
        included), plus any module setting without a control, which keeps
        the value the last loaded profile gave it (or its default). So saving
        a profile does not lose it, and a run from the window uses the same
        settings as a batch run without the window."""
        options = default_global_options()
        options.update(self.module_options)
        options.update(self.collect_widget_options())
        return options

    def module_option_notes(self):
        """'property=value' for each analysis-module setting (every setting
        that is not one of the main controls of the window) that differs
        from its default."""
        notes = []
        main_keys = self.main_widget_options()
        current = self.collect_global_options()
        for spec in GLOBAL_SETTING_SPECS:
            if spec[1] in main_keys or spec[1] not in current:
                continue
            value = current[spec[1]]
            if spec[2] == "bool":
                # Written as in the profile file (true / false).
                value = parse_bool(value, spec[3])
                differs = value != parse_bool(spec[3], False)
                text = bool_text(value)
            else:
                differs = to_text(value) != to_text(spec[3])
                text = to_text(value)
            if differs:
                notes.append(spec[0] + "=" + text)
        return notes

    def collect_widget_options(self):
        """The values of every control: the main controls and the Analyses tab."""
        options = self.main_widget_options()
        options.update(self.module_widget_options())
        return options

    def main_widget_options(self):
        return {
            "frame_mode": to_text(self.frame_mode_combo.getSelectedItem()),
            "roi_source": to_text(self.roi_source_combo.getSelectedItem()),
            "preset": to_text(self.preset_combo.getSelectedItem()),
            "override_calibration": self.override_calibration_cb.isSelected(),
            "pixel_width": text_value(self.pixel_width_field),
            "pixel_height": text_value(self.pixel_height_field),
            "spatial_unit": text_value(self.spatial_unit_field),
            "save_prepared": self.save_prepared_cb.isSelected(),
            "save_masks": self.save_masks_cb.isSelected(),
            "save_parent_rois": self.save_parent_rois_cb.isSelected(),
            "save_object_rois": self.save_object_rois_cb.isSelected(),
            "save_overlays": self.save_overlays_cb.isSelected(),
            "show_summary": self.show_summary_cb.isSelected(),
            "keep_prepared_open": self.keep_prepared_open_cb.isSelected(),
            "keep_masks_open": self.keep_masks_open_cb.isSelected(),
            "close_source": self.close_source_cb.isSelected(),
            "close_generated_after": self.close_generated_after_cb.isSelected(),
            "tile_prepared": self.tile_prepared_cb.isSelected(),
            "batch_extensions": text_value(self.batch_extensions_field),
            "batch_recursive": self.batch_recursive_cb.isSelected(),
            "batch_roi_policy": to_text(self.batch_roi_combo.getSelectedItem()),
            "batch_force_bioformats": self.batch_force_bioformats_cb.isSelected(),
            "batch_all_series": self.batch_all_series_cb.isSelected()
        }

    @ui_action
    def apply_global_options(self, options):
        set_combo(self.frame_mode_combo, options["frame_mode"])
        set_combo(self.roi_source_combo, options["roi_source"])
        set_combo(self.preset_combo, options["preset"])
        self.override_calibration_cb.setSelected(parse_bool(options["override_calibration"]))
        self.pixel_width_field.setText(to_text(options["pixel_width"]))
        self.pixel_height_field.setText(to_text(options["pixel_height"]))
        self.spatial_unit_field.setText(to_text(options["spatial_unit"]))
        self.save_prepared_cb.setSelected(parse_bool(options["save_prepared"]))
        self.save_masks_cb.setSelected(parse_bool(options["save_masks"]))
        self.save_parent_rois_cb.setSelected(parse_bool(options["save_parent_rois"]))
        self.save_object_rois_cb.setSelected(parse_bool(options["save_object_rois"]))
        self.save_overlays_cb.setSelected(parse_bool(options["save_overlays"]))
        self.show_summary_cb.setSelected(parse_bool(options["show_summary"]))
        self.keep_prepared_open_cb.setSelected(parse_bool(options["keep_prepared_open"]))
        self.keep_masks_open_cb.setSelected(parse_bool(options["keep_masks_open"]))
        self.close_source_cb.setSelected(parse_bool(options["close_source"]))
        self.close_generated_after_cb.setSelected(parse_bool(options["close_generated_after"]))
        self.tile_prepared_cb.setSelected(parse_bool(options["tile_prepared"]))
        self.batch_extensions_field.setText(to_text(options["batch_extensions"]))
        self.batch_recursive_cb.setSelected(parse_bool(options["batch_recursive"]))
        set_combo(self.batch_roi_combo, options["batch_roi_policy"])
        self.batch_force_bioformats_cb.setSelected(parse_bool(options["batch_force_bioformats"]))
        self.batch_all_series_cb.setSelected(parse_bool(options["batch_all_series"]))
        self.apply_module_widget_options(options)
        widget_keys = self.collect_widget_options()
        for spec in GLOBAL_SETTING_SPECS:
            if spec[1] not in widget_keys and spec[1] in options:
                self.module_options[spec[1]] = options[spec[1]]

    def confirm_source_changes(self, action):
        """Ask before measuring a source image whose pixels were changed in
        memory after it was opened. One Backspace or Delete key press while a
        selection is on the image runs Edit > Clear on the plane shown; one
        such key press reproduces exactly the 2.2 run of 2026-11-18 that lost
        Z1 of channel 1 inside its ROI.
        Returns False when the user declines. Asked once per source image."""
        imp = self.source_imp
        if imp is None or not imp.changes or self.accepted_changes_id == imp.getID():
            return True
        answer = JOptionPane.showConfirmDialog(
            self,
            "The source image '" + to_text(imp.getTitle()) + "' has been changed since it was opened, "
            "and the changes are not saved to its file.\n\n"
            "Common causes: the Backspace or Delete key while a selection is on the image (this runs "
            "Edit > Clear and fills the selection on the plane shown), or 'Apply' in Brightness/Contrast.\n"
            "The measurement will use the pixels as they are now.\n\n"
            "Continue anyway? (Answer No, then open the file again, to measure the original data.)",
            APP_NAME,
            JOptionPane.YES_NO_OPTION,
            JOptionPane.WARNING_MESSAGE
        )
        if answer != JOptionPane.YES_OPTION:
            self.append_status(action + " cancelled: the source image has unsaved changes.")
            return False
        self.accepted_changes_id = imp.getID()
        self.append_status("Warning: continuing with a source image that has unsaved changes: " + to_text(imp.getTitle()))
        return True

    def validate_configs(self, configs, for_analysis=True):
        if self.source_imp is None:
            raise ValueError("No source image is selected.")
        if WindowManager.getImage(self.source_imp.getID()) is None:
            raise ValueError("The selected source image is no longer open. Open it again or choose another active image.")
        validate_channel_configs(self.source_imp, configs)

        if for_analysis and to_text(self.frame_mode_combo.getSelectedItem()) == "All frames":
            for config in configs:
                if self.needs_mask(config) and config["threshold_mode"] == "Display minimum":
                    raise ValueError(
                        "Display-minimum thresholds are not reproducible across all frames. "
                        "Use Auto or Manual for C" + str(config["channel"]) + ", or process the current frame only."
                    )

    def needs_mask(self, config):
        return needs_mask(config)

    def needs_objects(self, config):
        return needs_objects(config)

    def needs_positive_area(self, config):
        return needs_positive_area(config)

    def get_calibration_override(self):
        if not self.override_calibration_cb.isSelected():
            return None
        pixel_width = parse_float_text(self.pixel_width_field, "Pixel width", 0.0000000001)
        pixel_height = parse_float_text(self.pixel_height_field, "Pixel height", 0.0000000001)
        unit = text_value(self.spatial_unit_field)
        if unit == "":
            unit = "pixel"
        return (pixel_width, pixel_height, unit)

    def get_output_root(self):
        root = text_value(self.output_dir_field)
        if root == "":
            if self.source_path:
                root = os.path.dirname(self.source_path)
            else:
                root = IJ.getDirectory("home")
        if root is None or to_text(root).strip() == "":
            raise ValueError("Choose an output folder.")
        root = os.path.abspath(to_text(root))
        ensure_dir(root)
        self.output_dir_field.setText(root)
        return root

    # ----- prepared-image cache and preprocessing -----------------------------

    def invalidate_prepared_cache(self):
        self.prepared_signature = None
        self.prepared_source_id = None
        self.prepared_frame = None

    def prep_signature(self, configs, frame):
        values = [self.source_imp.getID(), frame]
        calibration = self.get_calibration_override()
        values.append(calibration)
        # A prepared image carries the calibration the source had when it was
        # made, so a later Set Scale on the source makes it out of date.
        source_calibration = self.source_imp.getCalibration()
        if source_calibration is not None:
            values.append((
                source_calibration.pixelWidth, source_calibration.pixelHeight,
                to_text(source_calibration.getUnit())
            ))
        for config in configs:
            # "Current Z" prepares the slice the source image is showing, so
            # the prepared image is out of date once the Z slider is moved.
            shown_z = self.source_imp.getZ() if config["projection"] == "Current Z" else 0
            values.append((
                config["channel"], config["name"], config["projection"], config["z_slice"],
                config["subtract_background"], config["background_radius"], config["light_background"],
                config["median_radius"], config["gaussian_sigma"], config["lut"], shown_z
            ))
        return repr(values)

    def prepared_cache_valid(self, configs, frame):
        signature = self.prep_signature(configs, frame)
        if self.prepared_signature != signature:
            return False
        if self.prepared_source_id != self.source_imp.getID() or self.prepared_frame != frame:
            return False
        for config in configs:
            imp = self.prepared_images.get(config["channel"])
            if imp is None or WindowManager.getImage(imp.getID()) is None:
                return False
            # Pixels edited by hand (Edit > Clear, Brightness/Contrast
            # 'Apply', a filter run on the window) no longer follow from the
            # settings, so the image is not used for a measurement.
            if imp.changes:
                return False
        return True

    def edited_prepared_channels(self, configs):
        """Channels whose open prepared image has pixels edited by hand."""
        edited = []
        for config in configs:
            imp = self.prepared_images.get(config["channel"])
            if imp is not None and WindowManager.getImage(imp.getID()) is not None and imp.changes:
                edited.append(config["channel"])
        return edited

    def register_generated_image(self, imp):
        if imp is not None:
            self.generated_ids.add(imp.getID())

    def is_generated_image(self, imp):
        """True for a prepared image or a mask this window made."""
        if imp is None:
            return False
        if imp.getID() in self.generated_ids:
            return True
        for prepared in self.prepared_images.values():
            if prepared is imp:
                return True
        return False

    def close_prepared_images(self):
        for imp in list(self.prepared_images.values()):
            if imp is not None:
                self.generated_ids.discard(imp.getID())
                safe_close_image(imp)
        self.prepared_images = {}
        self.invalidate_prepared_cache()

    @ui_action
    def prepare_enabled_channels(self, event=None):
        if self.blocked_by_image_job():
            return
        self.set_busy(True)
        try:
            configs = self.collect_channel_configs()
            self.validate_configs(configs, for_analysis=False)
            if not self.confirm_source_changes("Preparation"):
                return
            frame = self.source_imp.getFrame()
            self.prepare_channels_for_frame(configs, frame, show_images=True, replace_cache=True)
            self.append_status("Prepared " + str(len(configs)) + " enabled channels for frame " + str(frame) + ".")
            if self.tile_prepared_cb.isSelected():
                IJ.run("Tile")
        except ANY_ERROR as exc:
            self.report_exception("Channel preparation failed", exc)
        finally:
            self.set_busy(False)

    @ui_action
    def preview_prepared_channel(self, channel_number):
        if self.blocked_by_image_job():
            return
        self.set_busy(True)
        try:
            configs = self.collect_channel_configs()
            self.validate_configs(configs, for_analysis=False)
            target = None
            for config in configs:
                if config["channel"] == channel_number:
                    target = config
                    break
            if target is None:
                raise ValueError("C" + str(channel_number) + " is not enabled.")
            frame = self.source_imp.getFrame()
            if not self.prepared_cache_valid(configs, frame):
                self.prepare_channels_for_frame(configs, frame, show_images=True, replace_cache=True)
            imp = self.prepared_images[channel_number]
            imp.show()
            WindowManager.setCurrentWindow(imp.getWindow())
            self.append_status("Prepared/refreshed C" + str(channel_number) + " (" + target["name"] + ").")
        except ANY_ERROR as exc:
            self.report_exception("Channel preview failed", exc)
        finally:
            self.set_busy(False)

    def prepare_channels_for_frame(self, configs, frame, show_images, replace_cache=False):
        if replace_cache:
            self.close_prepared_images()
        prepared = {}
        for config in configs:
            IJ.showStatus("Preparing C" + str(config["channel"]) + ": " + config["name"])
            imp = self.prepare_single_channel(self.source_imp, config, frame)
            # A new image holds no edit by hand; see prepared_cache_valid.
            imp.changes = False
            prepared[config["channel"]] = imp
            self.register_generated_image(imp)
            if show_images:
                imp.show()
        if replace_cache:
            self.prepared_images = prepared
            self.prepared_signature = self.prep_signature(configs, frame)
            self.prepared_source_id = self.source_imp.getID()
            self.prepared_frame = frame
        return prepared

    def prepare_single_channel(self, source, config, frame):
        return self.engine.prepare_single_channel(source, config, frame, self.get_calibration_override())

    def apply_lut(self, imp, lut_name):
        apply_lut(imp, lut_name)

    def require_display_range_current(self, configs, target):
        # A display-minimum threshold is the Brightness/Contrast minimum set
        # by eye on one prepared image. Once that image has to be rebuilt the
        # setting is gone, so stop instead of thresholding at a default.
        if not (self.needs_mask(target) and target["threshold_mode"] == "Display minimum"):
            return
        if self.prepared_cache_valid(configs, self.source_imp.getFrame()):
            return
        raise ValueError(
            "Display-minimum thresholding requires the current prepared images. "
            "Click '1. Prepare enabled channels', adjust Brightness/Contrast for C" +
            str(target["channel"]) + " (" + target["name"] + "), and try again without changing preprocessing settings."
        )

    def ensure_current_prepared(self, configs):
        frame = self.source_imp.getFrame()
        if not self.prepared_cache_valid(configs, frame):
            self.prepare_channels_for_frame(configs, frame, show_images=True, replace_cache=True)
        return self.prepared_images

    @ui_action
    def open_brightness_contrast(self, channel_number):
        if self.blocked_by_image_job():
            return
        try:
            configs = self.collect_channel_configs()
            self.validate_configs(configs, for_analysis=False)
            prepared = self.ensure_current_prepared(configs)
            imp = prepared.get(channel_number)
            if imp is None:
                raise ValueError("C" + str(channel_number) + " is not enabled.")
            imp.show()
            WindowManager.setCurrentWindow(imp.getWindow())
            IJ.run("Brightness/Contrast...")
            self.append_status(
                "Opened Brightness/Contrast for C" + str(channel_number) +
                ". Adjust the display range only; do not click Apply unless pixel rescaling is intentional."
            )
        except ANY_ERROR as exc:
            self.report_exception("Could not open Brightness/Contrast", exc)

    @ui_action
    def capture_display_threshold(self, channel_number):
        if self.blocked_by_image_job():
            return
        try:
            configs = self.collect_channel_configs()
            self.validate_configs(configs, for_analysis=False)
            if not self.prepared_cache_valid(configs, self.source_imp.getFrame()):
                # Preparing now would capture a default display range, not
                # the one that was set by eye.
                raise ValueError(
                    "There is no current prepared image to capture from. Click '1. Prepare enabled channels', "
                    "adjust Brightness/Contrast for C" + str(channel_number) + ", and capture again."
                )
            prepared = self.ensure_current_prepared(configs)
            imp = prepared.get(channel_number)
            if imp is None:
                raise ValueError("C" + str(channel_number) + " is not enabled.")
            low = imp.getDisplayRangeMin()
            stats = imp.getStatistics(Measurements.MIN_MAX)
            high = stats.max
            panel = self.channel_panels[channel_number - 1]
            panel.manual_low_field.setText(str(low))
            panel.manual_high_field.setText(str(high))
            panel.threshold_combo.setSelectedItem("Manual")
            self.append_status(
                "Captured C" + str(channel_number) + " display minimum as manual threshold: " + str(low) + " to " + str(high)
            )
        except ANY_ERROR as exc:
            self.report_exception("Could not capture display threshold", exc)

    # ----- mask construction and preview -------------------------------------

    @ui_action
    def preview_all_masks(self, event=None):
        if self.blocked_by_batch():
            return
        self.set_busy(True)
        try:
            configs = self.collect_channel_configs()
            self.validate_configs(configs, for_analysis=False)
            prepared = self.ensure_current_prepared(configs)
            found = 0
            for config in configs:
                if self.needs_mask(config):
                    self.preview_mask_internal(config, prepared[config["channel"]])
                    found += 1
            if found == 0:
                raise ValueError("No enabled channel uses Positive area or Objects analysis.")
            self.append_status("Previewed masks for " + str(found) + " channel(s).")
        except ANY_ERROR as exc:
            self.report_exception("Mask preview failed", exc)
        finally:
            self.set_busy(False)

    @ui_action
    def preview_mask(self, channel_number):
        if self.blocked_by_batch():
            return
        self.set_busy(True)
        try:
            configs = self.collect_channel_configs()
            self.validate_configs(configs, for_analysis=False)
            target = None
            for config in configs:
                if config["channel"] == channel_number:
                    target = config
                    break
            if target is None:
                raise ValueError("C" + str(channel_number) + " is not enabled.")
            if not self.needs_mask(target):
                raise ValueError("C" + str(channel_number) + " is set to Intensity only. Select Positive area or an Objects mode to create a mask.")
            prepared = self.ensure_current_prepared(configs)
            self.preview_mask_internal(target, prepared[channel_number])
        except ANY_ERROR as exc:
            self.report_exception("Mask preview failed", exc)
        finally:
            self.set_busy(False)

    def preview_mask_internal(self, config, prepared_imp):
        channel = config["channel"]
        self.close_mask_preview(channel)
        prefix = to_text(prepared_imp.getTitle()).replace("_PREP", "") + "_PREVIEW"
        mask_data = self.build_masks(prepared_imp, config, prefix)
        ids = []
        raw = mask_data["raw"]
        final = mask_data["final"]
        raw.show()
        self.register_generated_image(raw)
        ids.append(raw.getID())
        if final is not raw:
            final.show()
            self.register_generated_image(final)
            ids.append(final.getID())
        self.mask_preview_ids[channel] = ids
        if mask_data["threshold_mode"] == "Local Phansalkar":
            self.append_status(
                "C" + str(channel) + " mask preview: Local Phansalkar radius=" +
                str(config["local_radius"]) +
                (" with CLAHE." if config["clahe"] else ".")
            )
        else:
            self.append_status(
                "C" + str(channel) + " mask preview: threshold " + str(mask_data["threshold_low"]) +
                " to " + str(mask_data["threshold_high"])
            )

    def close_mask_preview(self, channel):
        ids = self.mask_preview_ids.get(channel, [])
        for image_id in ids:
            imp = WindowManager.getImage(image_id)
            if imp is not None:
                self.generated_ids.discard(image_id)
                safe_close_image(imp)
        self.mask_preview_ids[channel] = []

    def build_masks(self, prepared_imp, config, title_prefix):
        # Segmentation lives in AnalysisEngine so the window, a batch, and the
        # tests all build masks with exactly the same code.
        return self.engine.build_masks(prepared_imp, config, title_prefix)

    # ----- measurement helpers ----------------------------------------------

    def get_manager_rois(self):
        manager = RoiManager.getInstance2()
        if manager is None or manager.getCount() == 0:
            raise ValueError("ROI Manager contains no parent ROIs.")
        return [clone_roi(roi) for roi in manager.getRoisAsArray()]

    def get_parent_rois(self):
        if self.source_imp is None:
            raise ValueError("No source image selected.")
        source = to_text(self.roi_source_combo.getSelectedItem())
        rois = []

        if source == "Whole image":
            rois = [whole_image_roi(self.source_imp)]

        elif source == "Current image ROI":
            roi = self.source_imp.getRoi()
            if roi is None:
                raise ValueError("The source image has no active ROI.")
            rois = [clone_roi(roi)]

        elif source == "ROI Manager":
            rois = self.get_manager_rois()
        else:
            raise ValueError("Unsupported ROI source: " + source)

        return validate_parent_rois(rois, self.source_imp.getWidth(), self.source_imp.getHeight())

    def measure_intensity(self, imp, roi):
        return self.engine.measure_intensity(imp, roi)

    def measure_positive_area(self, mask_imp, roi):
        return self.engine.measure_positive_area(mask_imp, roi)

    def particle_size_pixels(self, config, calibration):
        return self.engine.particle_size_pixels(config, calibration)

    def analyze_objects(self, mask_imp, parent_roi, config):
        return self.engine.analyze_objects(mask_imp, parent_roi, config)

    # ----- full analysis ------------------------------------------------------

    @ui_action
    def quick_run_action(self, event=None):
        if self.blocked_by_batch():
            return
        try:
            configs = self.collect_channel_configs()
            self.validate_configs(configs, for_analysis=True)
            for config in configs:
                if self.needs_mask(config) and config["threshold_mode"] == "Display minimum":
                    raise ValueError(
                        "Quick run cannot pause for manual Brightness/Contrast. Change Display minimum to Auto or Manual, "
                        "or use Prepare followed by Run analysis."
                    )
            self.run_analysis_action()
        except ANY_ERROR as exc:
            self.report_exception("Quick run cannot start", exc)

    @ui_action
    def run_analysis_action(self, event=None):
        self.run_analysis()

    @ui_action
    def run_analysis(self):
        """Start the analysis of the source image in the background.
        Returns True when it started; the result is shown when it ends."""
        return self.start_background(
            "Analysis", self.collect_analysis_job, self.analysis_work, self.analysis_finished,
            "Analysis failed", uses_images=True
        )

    def create_run_folder(self, output_root, source_title, kind, configs=None, frames=None):
        """Name and create the folder of one single-image run: kind "analysis"
        (3. Run analysis and Quick run), "Validation" or "Sensitivity".
        This is the only place the window names such a folder."""
        # The engine names it: <image>_CellQuant_[<kind>_]<time>, with the
        # image name shortened only when a result path would otherwise be
        # too long for Windows programs (plan_run_folder).
        base = sanitize_name(strip_extension(source_title))
        multi_frame = self.source_imp is None or self.source_imp.getNFrames() > 1
        options = self.collect_global_options() if kind == "analysis" else None
        return unique_run_dir(output_root, base, configs, "" if kind == "analysis" else kind,
                              multi_frame=multi_frame, options=options)

    def collect_analysis_job(self):
        """Event thread: everything the analysis reads from the window, the
        checks that can refuse it, and its run folder."""
        configs = self.collect_channel_configs()
        self.validate_configs(configs, for_analysis=True)
        if not self.confirm_source_changes("Analysis"):
            return None
        source_changed = self.source_imp.changes
        parent_rois = self.get_parent_rois()
        output_root = self.get_output_root()
        global_options = self.collect_global_options()
        all_configs = self.collect_all_channel_configs()
        engine_options = engine_options_from_global(global_options)

        frame_mode = to_text(self.frame_mode_combo.getSelectedItem())
        if frame_mode == "All frames":
            frames = range(1, self.source_imp.getNFrames() + 1)
        else:
            frames = [self.source_imp.getFrame()]

        # A display-minimum threshold is meaningful only on the exact prepared
        # image whose Brightness/Contrast range the user adjusted. Do not
        # silently regenerate that image and substitute a default display range.
        display_threshold_channels = [
            config for config in configs
            if self.needs_mask(config) and config["threshold_mode"] == "Display minimum"
        ]
        if len(display_threshold_channels) > 0:
            current_frame = frames[0]
            if frame_mode != "Current frame only" or not self.prepared_cache_valid(configs, current_frame):
                labels = ", ".join([
                    "C" + str(config["channel"]) + " (" + config["name"] + ")"
                    for config in display_threshold_channels
                ])
                edited = self.edited_prepared_channels(configs)
                raise ValueError(
                    "Display-minimum thresholding requires the current prepared images. "
                    "Click '1. Prepare enabled channels', adjust Brightness/Contrast for " + labels +
                    ", preview the masks, and then run the analysis without changing preprocessing settings." +
                    ("" if len(edited) == 0 else
                     "\nThe pixels of the prepared image of " + ", ".join(["C" + str(channel) for channel in edited]) +
                     " were edited (for example with 'Apply' in Brightness/Contrast, or the Backspace key while a "
                     "selection was on it). Move only the display range sliders; do not press 'Apply'.")
                )

        run_dir = self.create_run_folder(output_root, self.source_imp.getTitle(), "analysis", configs, frames)
        self.last_run_dir = run_dir
        self.run_log_messages = []
        # Recorded in run_log.txt of this run.
        if source_changed:
            self.append_status(
                "Note: the source image had unsaved changes in memory (it was edited after it was opened); "
                "this run measures the edited pixels."
            )
        module_notes = self.module_option_notes()
        if len(module_notes) > 0:
            self.append_status("Analysis-module options in use: " + ", ".join(module_notes))
        edited = self.edited_prepared_channels(configs)
        if len(edited) > 0:
            self.append_status(
                "The prepared image of " + ", ".join(["C" + str(number) for number in edited]) +
                " was edited by hand, so the channels are prepared again from the source image for this run."
            )

        use_cache = (
            frame_mode == "Current frame only" and
            self.prepared_cache_valid(configs, frames[0])
        )
        close_generated = self.close_generated_after_cb.isSelected()
        # A copy: the worker keeps reading these images whatever the window does.
        prepared_cache = dict(self.prepared_images) if use_cache else None
        return {
            "source": self.source_imp, "source_path": self.source_path, "configs": configs,
            "parent_rois": parent_rois, "engine_options": engine_options, "run_dir": run_dir,
            "frames": list(frames),
            "prepared_cache": prepared_cache,
            "watch": self.image_watch(self.source_imp, configs, prepared_cache),
            "settings": (global_options, all_configs),
            "retain_masks": self.keep_masks_open_cb.isSelected() and not close_generated,
            "close_generated": close_generated,
            "show_summary": self.show_summary_cb.isSelected(),
            "keep_prepared_open": self.keep_prepared_open_cb.isSelected(),
            "close_source": self.close_source_cb.isSelected()
        }

    def analysis_work(self, job):
        """Worker thread: the measurement itself. No Swing here."""
        values = job["values"]
        run_dir = values["run_dir"]
        try:
            result = job["engine"].analyze_image(
                values["source"], values["source_path"], values["configs"], values["parent_rois"],
                values["engine_options"], run_dir, values["frames"],
                prepared_cache=values["prepared_cache"],
                settings=values["settings"],
                retain_masks=values["retain_masks"]
            )
        except AnalysisCancelled:
            self.append_status("Analysis cancelled by user. The run folder keeps only what was written before: " + to_text(run_dir))
            self.save_run_log(os.path.join(run_dir, "run_log.txt"))
            raise
        self.save_run_log(os.path.join(run_dir, "run_log.txt"))
        return result

    def analysis_finished(self, job):
        """Event thread: show what the worker made, as before."""
        values = job["values"]
        result = job["result"]
        if result is None:
            if job["cancelled"]:
                JOptionPane.showMessageDialog(
                    self,
                    "Analysis cancelled.\n\nThe run folder holds only what was written before the cancel:\n" +
                    to_text(values["run_dir"]),
                    APP_NAME,
                    JOptionPane.INFORMATION_MESSAGE
                )
            return
        warning = self.changed_during_job(job)
        if warning != "":
            # run_log.txt was written by the worker; add the warning to it.
            self.save_run_log(os.path.join(values["run_dir"], "run_log.txt"))

        if values["show_summary"]:
            self.show_summary_table(result["summary_headers"], result["summary_rows"])

        if values["close_generated"]:
            self.close_generated_images()
        else:
            for imp in result["mask_images"]:
                self.register_generated_image(imp)
                if WindowManager.getImage(imp.getID()) is None:
                    imp.show()
            if not values["keep_prepared_open"]:
                self.close_prepared_images()

        if values["close_source"] and self.source_imp is not None and self.source_imp is values["source"]:
            safe_close_image(self.source_imp)
            self.source_imp = None
            self.source_path = ""
            self.update_source_info()

        self.show_results(result, values["configs"], values["run_dir"], values["source"].getTitle())
        self.show_job_result("Analysis complete.", values["run_dir"], warning)

    def show_job_result(self, headline, run_dir, warning=""):
        """The closing message box of a background job; a warning from
        changed_during_job() turns it into a warning box."""
        text = headline + "\n\nSaved to:\n" + to_text(run_dir)
        if warning != "":
            text += "\n\nWarning: " + warning
        JOptionPane.showMessageDialog(
            self, text, APP_NAME, JOptionPane.WARNING_MESSAGE if warning != "" else JOptionPane.INFORMATION_MESSAGE
        )

    # ----- saving and display -------------------------------------------------

    def save_tiff(self, imp, path):
        save_tiff(imp, path)

    def save_roi_array(self, rois, path):
        save_roi_array(rois, path)

    def show_table(self, title, headers, rows):
        table = ResultsTable()
        for row in rows:
            table.incrementCounter()
            for header in headers:
                value = row.get(header, "")
                try:
                    if isinstance(value, bool):
                        table.addValue(header, 1 if value else 0)
                    elif is_number(value):
                        table.addValue(header, float(value))
                    else:
                        table.addValue(header, to_text(value))
                except ANY_ERROR:
                    pass
        table.show(title)

    def show_summary_table(self, headers, rows):
        self.show_table("CellQuant ROI Summary", headers, rows)

    def save_run_log(self, path):
        write_text_file(path, [APP_NAME + " " + APP_VERSION] + self.run_log_messages)

    # ----- profiles -----------------------------------------------------------

    @ui_action
    def save_profile(self, event=None):
        chooser = JFileChooser()
        chooser.setSelectedFile(File("CellQuant_Profile.properties"))
        if chooser.showSaveDialog(self) != JFileChooser.APPROVE_OPTION:
            return
        path = chooser.getSelectedFile().getAbsolutePath()
        if not path.lower().endswith(".properties"):
            path += ".properties"
        try:
            self.save_settings_to_path(path)
            self.append_status("Saved settings profile: " + path)
        except ANY_ERROR as exc:
            self.report_exception("Could not save settings profile", exc)

    def save_settings_to_path(self, path):
        write_settings_properties(path, self.collect_global_options(), self.collect_all_channel_configs())

    @ui_action
    def load_profile(self, event=None):
        chooser = JFileChooser()
        if chooser.showOpenDialog(self) != JFileChooser.APPROVE_OPTION:
            return
        self.load_profile_from_path(chooser.getSelectedFile().getAbsolutePath())

    @ui_action
    def load_profile_from_path(self, path):
        try:
            props = read_properties_file(path)
            loaded = global_options_from_properties(props)
            options = self.collect_global_options()
            main_keys = self.main_widget_options()
            for spec in GLOBAL_SETTING_SPECS:
                present = props.getProperty(spec[0]) is not None
                # A profile from 2.2.x has no batch keys: leave those as they are.
                # An analysis-module setting the profile does not hold gets its
                # default, as in a batch run without the window, whether or
                # not it has a control on the Analyses tab.
                if present or (spec[2] == "bool" and not spec[0].startswith("batch.")) or spec[1] not in main_keys:
                    options[spec[1]] = loaded[spec[1]]
            self.apply_global_options(options)

            for panel in self.channel_panels:
                panel.read_properties(props)
            self.invalidate_prepared_cache()
            saved_version = props.getProperty("application.version")
            self.append_status(
                "Loaded settings profile: " + to_text(path) +
                (" (saved by version " + to_text(saved_version) + ")" if saved_version is not None else "")
            )
            notes = self.module_option_notes()
            if len(notes) > 0:
                self.append_status(
                    "The profile also sets analysis-module options (shown on the Analyses tab); "
                    "they are used by runs and batches: " + ", ".join(notes)
                )
        except ANY_ERROR as exc:
            self.report_exception("Could not load settings profile", exc)

    # ----- batch --------------------------------------------------------------

    @ui_action
    def choose_batch_input_dir(self, event=None):
        chooser = JFileChooser()
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
        if chooser.showOpenDialog(self) == JFileChooser.APPROVE_OPTION:
            self.batch_input_field.setText(chooser.getSelectedFile().getAbsolutePath())

    @ui_action
    def choose_batch_roi_dir(self, event=None):
        chooser = JFileChooser()
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
        if chooser.showOpenDialog(self) == JFileChooser.APPROVE_OPTION:
            self.batch_roi_dir_field.setText(chooser.getSelectedFile().getAbsolutePath())

    def collect_batch_options(self):
        return {
            "extensions": text_value(self.batch_extensions_field),
            "recursive": self.batch_recursive_cb.isSelected(),
            "roi_policy": to_text(self.batch_roi_combo.getSelectedItem()),
            "roi_dir": text_value(self.batch_roi_dir_field),
            "force_bioformats": self.batch_force_bioformats_cb.isSelected(),
            "all_series": self.batch_all_series_cb.isSelected()
        }

    def get_batch_input_dir(self):
        input_dir = text_value(self.batch_input_field)
        if input_dir == "":
            raise ValueError("Choose a batch input folder.")
        input_dir = os.path.abspath(input_dir)
        if not os.path.isdir(input_dir):
            raise ValueError("Batch input folder does not exist: " + input_dir)
        return input_dir

    @ui_action
    def scan_batch_folder(self, event=None):
        try:
            options = self.collect_batch_options()
            input_dir = self.get_batch_input_dir()
            extensions = parse_extension_list(options["extensions"])
            if len(extensions) == 0:
                raise ValueError("List at least one file extension for the batch.")
            files = discover_batch_images(input_dir, extensions, options["recursive"])
            uses_roi_files = options["roi_policy"] in [BATCH_ROI_FILE_SKIP, BATCH_ROI_FILE_WHOLE]
            roi_dir = options["roi_dir"] if options["roi_dir"] != "" else None
            lines = []
            without_roi = 0
            for path in files:
                relative = to_text(os.path.relpath(path, input_dir))
                if uses_roi_files:
                    roi_path = find_roi_file(path, roi_dir)
                    if roi_path is None:
                        without_roi += 1
                        lines.append(relative + "    [no ROI file]")
                    else:
                        lines.append(relative + "    [" + to_text(os.path.basename(roi_path)) + "]")
                else:
                    lines.append(relative)
            header = [str(len(files)) + " image file(s) found in " + to_text(input_dir)]
            if uses_roi_files:
                header.append(str(len(files) - without_roi) + " with an ROI file, " + str(without_roi) + " without.")
            header.append("")
            self.batch_list_area.setText("\n".join(header + lines))
            self.batch_list_area.setCaretPosition(0)
            self.append_status("Batch scan: " + str(len(files)) + " image file(s) in " + to_text(input_dir))
        except ANY_ERROR as exc:
            self.report_exception("Folder scan failed", exc)

    @ui_action
    def save_rois_for_current_image(self, event=None):
        try:
            if self.source_imp is None or self.source_path == "":
                raise ValueError("Open the image from a file first, so the ROI set can be named after it.")
            rois = self.get_manager_rois()
            roi_dir = text_value(self.batch_roi_dir_field)
            folder = roi_dir if roi_dir != "" else os.path.dirname(self.source_path)
            path = os.path.join(folder, strip_extension(self.source_path) + "_ROIs.zip")
            if os.path.exists(path):
                result = JOptionPane.showConfirmDialog(
                    self, "Replace the existing ROI file?\n" + to_text(path), APP_NAME, JOptionPane.YES_NO_OPTION
                )
                if result != JOptionPane.YES_OPTION:
                    return
            ensure_dir(folder)
            save_roi_array(rois, path)
            self.append_status("Saved " + str(len(rois)) + " ROI(s) for batch use: " + to_text(path))
        except ANY_ERROR as exc:
            self.report_exception("Could not save the ROI set", exc)

    @ui_action
    def run_batch_action(self, event=None):
        if self.blocked_by_batch():
            return
        try:
            configs = self.collect_channel_configs(limit_to_source=False)
            if len(configs) == 0:
                raise ValueError("No channels are enabled for analysis.")
            validate_unattended_configs(configs, "Batch processing")
            all_configs = self.collect_all_channel_configs()
            global_options = self.collect_global_options()
            engine_options_from_global(global_options)
            batch_options = self.collect_batch_options()
            input_dir = self.get_batch_input_dir()
            output_root = text_value(self.output_dir_field)
            if output_root == "":
                output_root = input_dir
            output_root = os.path.abspath(output_root)
            ensure_dir(output_root)
            extensions = parse_extension_list(batch_options["extensions"])
            if len(extensions) == 0:
                raise ValueError("List at least one file extension for the batch.")
            files = discover_batch_images(input_dir, extensions, batch_options["recursive"])
            if len(files) == 0:
                raise ValueError("No image files with the listed extensions were found in " + input_dir)
            manager_rois = None
            if batch_options["roi_policy"] == BATCH_ROI_MANAGER:
                manager_rois = self.get_manager_rois()
            result = JOptionPane.showConfirmDialog(
                self,
                "Analyze " + str(len(files)) + " image file(s) with the current settings?\n\n" +
                "Regions: " + batch_options["roi_policy"] + "\nOutput folder: " + to_text(output_root),
                APP_NAME,
                JOptionPane.YES_NO_OPTION
            )
            if result != JOptionPane.YES_OPTION:
                return
        except ANY_ERROR as exc:
            self.report_exception("Batch cannot start", exc)
            return

        self.batch_cancel_requested = False
        self.batch_running = True
        self.batch_run_button.setEnabled(False)
        self.batch_cancel_button.setEnabled(True)
        for button in self.analysis_action_buttons:
            button.setEnabled(False)
        self.batch_progress_update(0, len(files), "starting")
        worker = threading.Thread(
            target=self.batch_worker,
            args=(input_dir, output_root, configs, all_configs, global_options, batch_options, manager_rois, files)
        )
        worker.setDaemon(True)
        self.batch_thread = worker
        worker.start()

    def batch_is_cancelled(self):
        return self.batch_cancel_requested

    def batch_worker(self, input_dir, output_root, configs, all_configs, global_options,
                     batch_options, manager_rois, files):
        outcome = {"result": None, "error": None}
        try:
            outcome["result"] = run_batch(
                input_dir, output_root, configs, all_configs, global_options, batch_options,
                manager_rois=manager_rois, log=self.append_status, progress=self.batch_progress_update,
                is_cancelled=self.batch_is_cancelled, file_list=files
            )
        except ANY_ERROR as exc:
            outcome["error"] = exc
            try:
                IJ.log(traceback.format_exc())
            except ANY_ERROR:
                pass
        finally:
            def finish():
                self.batch_finished(outcome)

            run_on_ui(finish)

    def batch_progress_update(self, done, total, label):
        def update():
            try:
                self.batch_progress.setMaximum(max(1, total))
                self.batch_progress.setValue(done)
                text = str(done) + " / " + str(total)
                if label:
                    text += "   " + to_text(label)
                self.batch_progress.setString(text)
            except ANY_ERROR:
                pass

        run_on_ui(update)

    @ui_action
    def cancel_batch_action(self, event=None):
        if not self.batch_running:
            return
        self.batch_cancel_requested = True
        self.batch_cancel_button.setEnabled(False)
        if self.background_active:
            # A robustness run uses the Batch tab's Cancel button too.
            self.job_cancel_button.setEnabled(False)
            self.append_status(self.background_label + " cancel requested; stopping after the current step.")
        else:
            self.append_status("Batch cancel requested; stopping after the current step.")

    @ui_action
    def batch_finished(self, outcome):
        self.batch_running = False
        self.batch_thread = None
        self.batch_run_button.setEnabled(True)
        self.batch_cancel_button.setEnabled(False)
        for button in self.analysis_action_buttons:
            button.setEnabled(True)
        IJ.showStatus("Cell Quant Studio ready")
        result = outcome["result"]
        if result is None:
            self.batch_progress.setString("Stopped")
            error = outcome["error"]
            self.report_exception("Batch failed", error if error is not None else "The batch stopped unexpectedly; see the Log window.")
            return
        self.last_run_dir = result["batch_dir"]
        self.batch_progress.setString("Cancelled" if result["cancelled"] else "Finished")
        message = (
            ("Batch cancelled.\n\n" if result["cancelled"] else "Batch complete.\n\n") +
            str(result["analyzed"]) + " analyzed, " + str(result["skipped"]) + " skipped, " +
            str(result["failed"]) + " failed of " + str(result["files"]) + " file(s).\n" +
            "See batch_report.csv for the status of every file.\n\nSaved to:\n" + to_text(result["batch_dir"])
        )
        JOptionPane.showMessageDialog(
            self, message, APP_NAME,
            JOptionPane.WARNING_MESSAGE if result["failed"] > 0 else JOptionPane.INFORMATION_MESSAGE
        )

    # ----- validation and sensitivity ----------------------------------------

    @ui_action
    def choose_annotation_file(self, event=None):
        chooser = JFileChooser()
        chooser.setFileSelectionMode(JFileChooser.FILES_ONLY)
        if chooser.showOpenDialog(self) == JFileChooser.APPROVE_OPTION:
            self.annotation_file_field.setText(chooser.getSelectedFile().getAbsolutePath())

    def get_annotation_file(self):
        path = text_value(self.annotation_file_field)
        if path == "":
            raise ValueError("Choose the annotation file.")
        if not os.path.isfile(path):
            raise ValueError("Annotation file does not exist: " + path)
        return path

    def qc_context(self, prepare_now=True):
        """Everything a validation or sensitivity run needs: the channel under
        test, its prepared image for the current frame, and the parent ROIs.
        With prepare_now False nothing is prepared here: "prepared_imp" is
        then None when the prepared images are not current, and
        prepare_qc_channels() makes them on the worker thread."""
        configs = self.collect_channel_configs()
        self.validate_configs(configs, for_analysis=False)
        channel_number = self.qc_channel_combo.getSelectedIndex() + 1
        target = None
        for config in configs:
            if config["channel"] == channel_number:
                target = config
                break
        if target is None:
            raise ValueError("C" + str(channel_number) + " is not enabled.")
        # Read point marks before anything else touches the image selections.
        selections = [self.source_imp.getRoi()]
        for imp in self.prepared_images.values():
            if imp is not None:
                selections.append(imp.getRoi())
        self.require_display_range_current(configs, target)
        frame = self.source_imp.getFrame()
        if prepare_now:
            prepared = self.ensure_current_prepared(configs)
        elif self.prepared_cache_valid(configs, frame):
            prepared = dict(self.prepared_images)
        else:
            prepared = None
        return {
            "configs": configs,
            "config": target,
            "prepared": prepared,
            "prepared_imp": prepared[channel_number] if prepared is not None else None,
            "selections": selections,
            "frame": frame,
            "source_title": to_text(self.source_imp.getTitle()),
            # What prepare_qc_channels needs, read here on the event thread.
            "source": self.source_imp,
            "calibration": self.get_calibration_override(),
            "signature": self.prep_signature(configs, frame)
        }

    def prepare_qc_channels(self, job):
        """Worker thread: prepare the enabled channels when the prepared
        images were not current. They are shown and become the window's
        prepared images when the job ends (install_qc_prepared), exactly as
        preparing them on the event thread did before."""
        context = job["values"]["context"]
        if context["prepared_imp"] is not None:
            return context["prepared_imp"]
        prepared = {}
        job["values"]["new_prepared"] = prepared
        for config in context["configs"]:
            IJ.showStatus("Preparing C" + str(config["channel"]) + ": " + config["name"])
            imp = job["engine"].prepare_single_channel(context["source"], config, context["frame"], context["calibration"])
            imp.changes = False
            prepared[config["channel"]] = imp
        return prepared[context["config"]["channel"]]

    def install_qc_prepared(self, job):
        """Event thread: show the images prepare_qc_channels made."""
        values = job["values"]
        prepared = values.get("new_prepared")
        if not prepared:
            return
        context = values["context"]
        if self.source_imp is not context["source"] or len(prepared) != len(context["configs"]):
            for imp in prepared.values():
                safe_close_image(imp)
            return
        self.close_prepared_images()
        for channel in sorted(prepared.keys()):
            self.register_generated_image(prepared[channel])
            prepared[channel].show()
        self.prepared_images = prepared
        self.prepared_signature = context["signature"]
        self.prepared_source_id = context["source"].getID()
        self.prepared_frame = context["frame"]

    @ui_action
    def run_validation_action(self, event=None):
        self.start_background(
            "Validation", self.collect_validation_job, self.validation_work, self.validation_finished,
            "Validation failed", uses_images=True
        )

    def collect_validation_job(self):
        source_kind = to_text(self.annotation_source_combo.getSelectedItem())
        if source_kind == ANNOTATION_POINTS_IMAGE and to_text(self.roi_source_combo.getSelectedItem()) == "Current image ROI":
            raise ValueError(
                "The image selection cannot be both the parent ROI and the point marks. "
                "Add the parent ROI to ROI Manager (or use Whole image), or load the marks from a point ROI file."
            )
        context = self.qc_context(prepare_now=False)
        points = None
        reference = None
        if source_kind == ANNOTATION_MASK_FILE:
            reference = IJ.openImage(self.get_annotation_file())
            if reference is None:
                raise ValueError("Could not open the reference mask image.")
        elif source_kind == ANNOTATION_POINTS_FILE:
            points = annotation_points_from_rois(load_rois_from_file(self.get_annotation_file()))
            if len(points) == 0:
                raise ValueError("The annotation file contains no point ROIs.")
        else:
            points = annotation_points_from_rois(context["selections"])
            if len(points) == 0:
                raise ValueError(
                    "No point marks were found. Mark one point per object with the Multi-point tool "
                    "on the source image, then run the validation while that selection is active."
                )
        started = False
        try:
            if not self.confirm_source_changes("Validation"):
                return None
            parent_rois = self.get_parent_rois()
            output_root = self.get_output_root()
            run_dir = self.create_run_folder(output_root, context["source_title"], "Validation", context["configs"], [context["frame"]])
            self.last_run_dir = run_dir
            self.save_settings_to_path(os.path.join(run_dir, "settings.properties"))
            save_roi_array(parent_rois, os.path.join(run_dir, "Parent_ROIs_Input.zip"))
            started = True
        finally:
            if not started:
                safe_close_image(reference)
        return {
            "context": context, "points": points, "reference": reference,
            "parent_rois": parent_rois, "run_dir": run_dir,
            "watch": self.image_watch(context["source"], context["configs"], context["prepared"])
        }

    def validation_work(self, job):
        values = job["values"]
        context = values["context"]
        try:
            prepared_imp = self.prepare_qc_channels(job)
            if values["reference"] is not None:
                return validate_mask_against_reference(
                    job["engine"], prepared_imp, context["config"], values["parent_rois"], values["reference"],
                    values["run_dir"], context["source_title"], context["frame"]
                )
            return validate_objects_against_points(
                job["engine"], prepared_imp, context["config"], values["parent_rois"], values["points"],
                values["run_dir"], context["source_title"], context["frame"]
            )
        finally:
            safe_close_image(values["reference"])

    def validation_finished(self, job):
        self.install_qc_prepared(job)
        result = job["result"]
        if result is None:
            return
        warning = self.changed_during_job(job)
        values = job["values"]
        config = values["context"]["config"]
        run_dir = values["run_dir"]
        if values["reference"] is not None:
            table_title = "CellQuant Mask Validation"
        else:
            table_title = "CellQuant Object Validation"
        self.show_table(table_title, result["headers"], result["rows"])
        self.append_status("Validation of C" + str(config["channel"]) + " (" + config["name"] + ") saved: " + to_text(run_dir))
        self.show_job_result("Validation complete.", run_dir, warning)

    @ui_action
    def run_sensitivity_action(self, event=None):
        self.start_background(
            "Sensitivity analysis", self.collect_sensitivity_job, self.sensitivity_work, self.sensitivity_finished,
            "Sensitivity analysis failed", uses_images=True
        )

    def collect_sensitivity_job(self):
        parameter = to_text(self.sensitivity_parameter_combo.getSelectedItem())
        factors = parse_factor_list(text_value(self.sensitivity_factors_field))
        context = self.qc_context(prepare_now=False)
        if not self.confirm_source_changes("Sensitivity analysis"):
            return None
        parent_rois = self.get_parent_rois()
        output_root = self.get_output_root()
        run_dir = self.create_run_folder(output_root, context["source_title"], "Sensitivity", context["configs"], [context["frame"]])
        self.last_run_dir = run_dir
        self.save_settings_to_path(os.path.join(run_dir, "settings.properties"))
        save_roi_array(parent_rois, os.path.join(run_dir, "Parent_ROIs_Input.zip"))
        return {
            "context": context, "parameter": parameter, "factors": factors,
            "parent_rois": parent_rois, "run_dir": run_dir,
            "watch": self.image_watch(context["source"], context["configs"], context["prepared"])
        }

    def sensitivity_work(self, job):
        values = job["values"]
        context = values["context"]
        prepared_imp = self.prepare_qc_channels(job)
        return run_sensitivity(
            job["engine"], prepared_imp, context["config"], values["parent_rois"], values["parameter"],
            values["factors"], values["run_dir"], context["source_title"], context["frame"]
        )

    def sensitivity_finished(self, job):
        self.install_qc_prepared(job)
        result = job["result"]
        if result is None:
            return
        warning = self.changed_during_job(job)
        values = job["values"]
        config = values["context"]["config"]
        run_dir = values["run_dir"]
        self.show_table("CellQuant Sensitivity", result["headers"], result["rows"])
        self.append_status(
            "Sensitivity of C" + str(config["channel"]) + " (" + config["name"] + ") to '" + values["parameter"] +
            "' saved: " + to_text(run_dir)
        )
        self.show_job_result("Sensitivity analysis complete.", run_dir, warning)

    # ----- separate runs of the analysis modules (Analyses tab) --------------

    @ui_action
    def detect_tissue_action(self, event=None):
        self.start_background(
            "Tissue outline detection", self.collect_tissue_job, self.tissue_work, self.tissue_finished,
            "Tissue outline detection failed", uses_images=True
        )

    def collect_tissue_job(self):
        configs = self.collect_channel_configs()
        self.validate_configs(configs, for_analysis=False)
        options = self.collect_global_options()
        # The module reads and checks its own settings; a bad value is
        # reported now, before the job starts.
        if "tissue_settings" in globals():
            tissue_settings(options)
        if not self.confirm_source_changes("Tissue outline detection"):
            return None
        return {
            "source": self.source_imp, "configs": configs, "options": options, "frame": self.source_imp.getFrame(),
            "watch": self.image_watch(self.source_imp, configs, None)
        }

    def tissue_work(self, job):
        values = job["values"]
        return detect_tissue_outline(job["engine"], values["source"], values["configs"], values["options"], frame=values["frame"])

    def tissue_finished(self, job):
        if job["result"] is None:
            return
        rois, info = job["result"][0], job["result"][1]
        description = to_text(info.get("description", ""))
        # The module logs its description through the engine; say it here
        # when it did not.
        if description != "" and len([line for line in self.run_log_messages[-50:] if description in line]) == 0:
            self.append_status(description)
        for warning in info.get("warnings", []):
            self.append_status("Tissue outline warning: " + to_text(warning))
        changed = self.changed_during_job(job)
        if len(rois) == 0:
            JOptionPane.showMessageDialog(
                self, "No tissue outline was found.\n\n" + description, APP_NAME, JOptionPane.WARNING_MESSAGE
            )
            return
        if changed != "":
            JOptionPane.showMessageDialog(self, changed + "\n\nThe outline was not added.", APP_NAME, JOptionPane.WARNING_MESSAGE)
            return
        manager = RoiManager.getRoiManager()
        for roi in rois:
            manager.addRoi(roi)
        self.roi_source_combo.setSelectedItem("ROI Manager")
        self.append_status(
            "Added " + str(len(rois)) + " tissue outline(s) to ROI Manager, which now holds " +
            str(manager.getCount()) + " ROI(s). The ROI source is now 'ROI Manager'."
        )

    def analysis_input_folder(self, field):
        folder = text_value(field)
        if folder == "":
            return self.get_batch_input_dir()
        folder = os.path.abspath(folder)
        if not os.path.isdir(folder):
            raise ValueError("Folder does not exist: " + folder)
        return folder

    def analysis_output_root(self, fallback):
        root = text_value(self.output_dir_field)
        if root == "":
            root = fallback
        root = os.path.abspath(to_text(root))
        ensure_dir(root)
        return root

    @ui_action
    def run_tuning_action(self, event=None):
        self.start_background("Tuning", self.collect_tuning_job, self.tuning_work, self.tuning_finished, "Tuning failed")

    def collect_tuning_job(self):
        configs = self.collect_channel_configs(limit_to_source=False)
        if len(configs) == 0:
            raise ValueError("No channels are enabled for analysis.")
        options = self.collect_global_options()
        engine_options_from_global(options)
        # Only the channel that is tuned must not use a by-eye threshold (the
        # module checks it again); a bad tune_channel is reported now.
        if "tune_config_from_options" in globals():
            validate_unattended_configs([tune_config_from_options(options, configs)], "Tuning")
        batch_options = self.collect_batch_options()
        folder = self.analysis_input_folder(self.tune_folder_field)
        return {
            "folder": folder, "output_root": self.analysis_output_root(folder), "configs": configs,
            "all_configs": self.collect_all_channel_configs(), "options": options,
            "folder_options": {
                "extensions": batch_options["extensions"], "recursive": batch_options["recursive"],
                "roi_dir": batch_options["roi_dir"] if batch_options["roi_dir"] != "" else None,
                "force_bioformats": batch_options["force_bioformats"]
            }
        }

    def tuning_work(self, job):
        values = job["values"]
        return run_tuning_folder(
            job["engine"], values["folder"], values["output_root"], values["configs"], values["all_configs"],
            values["options"], values["folder_options"], job["progress"]
        )

    def tuning_finished(self, job):
        result = job["result"]
        if result is None:
            return
        run_dir = to_text(result.get("run_dir", ""))
        if run_dir != "":
            self.last_run_dir = run_dir
        if result.get("cancelled"):
            self.mark_job_cancelled(job)
            JOptionPane.showMessageDialog(
                self, "Tuning cancelled. No tuned profile was written.\n\nSaved to:\n" + run_dir,
                APP_NAME, JOptionPane.INFORMATION_MESSAGE
            )
            return
        profile = result.get("profile_path")
        text = "Tuning complete.\n\n"
        if result.get("tuned_values") is not None:
            text += "Tuned settings: " + ", ".join([
                to_text(key) + "=" + to_text(result["tuned_values"][key]) for key in sorted(result["tuned_values"].keys())
            ]) + "\n"
        if result.get("design"):
            text += "Held-out design: " + to_text(result["design"]) + "\n"
        try:
            if result.get("held_out") is not None:
                scores = tune_metrics(result["held_out"])
                text += "Held-out F1 %.3f (precision %.3f, recall %.3f)\n" % (scores["f1"], scores["precision"], scores["recall"])
            else:
                text += "No held-out estimate: annotate at least two fields for a validation result.\n"
        except ANY_ERROR:
            pass
        text += "See tuning_log.txt and tuning_summary.csv for the scores.\n\nSaved to:\n" + run_dir
        if profile is None:
            JOptionPane.showMessageDialog(self, text, APP_NAME, JOptionPane.INFORMATION_MESSAGE)
            return
        answer = JOptionPane.showConfirmDialog(
            self, text + "\n\nLoad the tuned profile into the window now?\n" + to_text(profile),
            APP_NAME, JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE
        )
        if answer == JOptionPane.YES_OPTION:
            self.load_profile_from_path(profile)

    @ui_action
    def run_correction_action(self, event=None):
        self.start_background(
            "Count correction", self.collect_correction_job, self.correction_work, self.correction_finished,
            "Count correction failed"
        )

    def collect_correction_job(self):
        source = text_value(self.correct_source_field)
        if source == "":
            raise ValueError("Choose the folder of the validation runs, or a table of manual and detected counts.")
        source = os.path.abspath(source)
        if not os.path.exists(source):
            raise ValueError("Does not exist: " + source)
        table = text_value(self.correct_table_field)
        output_root = text_value(self.output_dir_field)
        return {
            "folder": source if os.path.isdir(source) else None,
            "csv": None if os.path.isdir(source) else source,
            "table": os.path.abspath(table) if table != "" else None,
            "channel": text_value(self.correct_channel_field),
            "output_root": os.path.abspath(output_root) if output_root != "" else None,
            "options": self.collect_global_options()
        }

    def correction_work(self, job):
        values = job["values"]
        return run_count_correction(
            calibration_folder=values["folder"], calibration_csv=values["csv"], table_path=values["table"],
            channel=values["channel"], output_root=values["output_root"], global_options=values["options"],
            log=self.append_status
        )

    def correction_finished(self, job):
        result = job["result"]
        if result is None:
            return
        out_dir = to_text(result["out_dir"])
        self.last_run_dir = out_dir
        warnings = result.get("warnings", [])
        text = "Count correction complete.\n\n"
        if result.get("corrected_path"):
            text += "Corrected table: " + to_text(os.path.basename(result["corrected_path"])) + "\n"
        if len(warnings) > 0:
            text += str(len(warnings)) + " warning(s); see the status log.\n"
        text += "\nSaved to:\n" + out_dir
        JOptionPane.showMessageDialog(
            self, text, APP_NAME, JOptionPane.WARNING_MESSAGE if len(warnings) > 0 else JOptionPane.INFORMATION_MESSAGE
        )

    @ui_action
    def run_robustness_action(self, event=None):
        self.start_background(
            "Robustness analysis", self.collect_robustness_job, self.robustness_work, self.robustness_finished,
            "Robustness analysis failed", batch_bar=True
        )

    def collect_robustness_job(self):
        """The same inputs as a batch (Batch tab), plus the robust_* settings."""
        configs = self.collect_channel_configs(limit_to_source=False)
        if len(configs) == 0:
            raise ValueError("No channels are enabled for analysis.")
        validate_unattended_configs(configs, "Robustness analysis")
        global_options = self.collect_global_options()
        engine_options_from_global(global_options)
        if "robust_request_from_options" in globals():
            robust_request_from_options(global_options)
        batch_options = self.collect_batch_options()
        input_dir = self.get_batch_input_dir()
        output_root = self.analysis_output_root(input_dir)
        extensions = parse_extension_list(batch_options["extensions"])
        if len(extensions) == 0:
            raise ValueError("List at least one file extension for the batch.")
        files = discover_batch_images(input_dir, extensions, batch_options["recursive"])
        if len(files) == 0:
            raise ValueError("No image files with the listed extensions were found in " + input_dir)
        manager_rois = None
        if batch_options["roi_policy"] == BATCH_ROI_MANAGER:
            manager_rois = self.get_manager_rois()
        answer = JOptionPane.showConfirmDialog(
            self,
            "Run the robustness analysis on " + str(len(files)) + " image file(s)?\n\n" +
            "Regions: " + batch_options["roi_policy"] + "\nOutput folder: " + to_text(output_root),
            APP_NAME,
            JOptionPane.YES_NO_OPTION
        )
        if answer != JOptionPane.YES_OPTION:
            return None
        return {
            "input_dir": input_dir, "output_root": output_root, "configs": configs,
            "all_configs": self.collect_all_channel_configs(), "options": global_options,
            "batch_options": batch_options, "manager_rois": manager_rois, "files": files
        }

    def robustness_work(self, job):
        values = job["values"]
        return run_robustness_from_options(
            values["input_dir"], values["output_root"], values["configs"], values["all_configs"],
            values["options"], values["batch_options"], manager_rois=values["manager_rois"],
            log=self.append_status, progress=job["progress"], is_cancelled=self.batch_is_cancelled,
            file_list=values["files"]
        )

    def robustness_finished(self, job):
        result = job["result"]
        if result is None:
            return
        if result.get("cancelled"):
            self.mark_job_cancelled(job)
        self.last_run_dir = to_text(result["run_dir"])
        text = (
            ("Robustness analysis cancelled; no statistics were written.\n\n" if result.get("cancelled")
             else "Robustness analysis complete.\n\n") +
            str(result.get("analyzed", 0)) + " analyzed, " + str(result.get("skipped", 0)) + " skipped, " +
            str(result.get("failed", 0)) + " failed of " + str(result.get("files", 0)) + " file(s); " +
            str(result.get("specifications", 0)) + " specifications.\n\nSaved to:\n" + to_text(result["run_dir"])
        )
        JOptionPane.showMessageDialog(
            self, text, APP_NAME,
            JOptionPane.WARNING_MESSAGE if result.get("failed", 0) > 0 else JOptionPane.INFORMATION_MESSAGE
        )

    # ----- cleanup and navigation --------------------------------------------

    @ui_action
    def close_generated_images(self):
        for image_id in list(self.generated_ids):
            imp = WindowManager.getImage(image_id)
            if imp is not None:
                safe_close_image(imp)
        self.generated_ids = set()
        self.prepared_images = {}
        self.mask_preview_ids = {}
        self.invalidate_prepared_cache()

    @ui_action
    def close_generated_images_action(self, event=None):
        if self.blocked_by_image_job():
            return
        self.close_generated_images()
        self.append_status("Closed all images generated by this application.")

    @ui_action
    def close_all_images_action(self, event=None):
        if self.blocked_by_batch():
            return
        result = JOptionPane.showConfirmDialog(
            self,
            "Close every open image in ImageJ, including unrelated images?\nUnsaved image changes will be discarded.",
            APP_NAME,
            JOptionPane.YES_NO_OPTION,
            JOptionPane.WARNING_MESSAGE
        )
        if result != JOptionPane.YES_OPTION:
            return
        ids = WindowManager.getIDList()
        if ids is not None:
            for image_id in list(ids):
                imp = WindowManager.getImage(image_id)
                safe_close_image(imp)
        self.generated_ids = set()
        self.prepared_images = {}
        self.mask_preview_ids = {}
        self.source_imp = None
        self.source_path = ""
        self.invalidate_prepared_cache()
        self.update_source_info()
        self.append_status("Closed all ImageJ image windows.")

    @ui_action
    def reset_for_next_image(self, event=None):
        if self.blocked_by_batch():
            return
        self.close_generated_images()
        if self.source_imp is not None:
            result = JOptionPane.showConfirmDialog(
                self,
                "Close the current source image as part of the reset?",
                APP_NAME,
                JOptionPane.YES_NO_OPTION
            )
            if result == JOptionPane.YES_OPTION:
                safe_close_image(self.source_imp)
        self.source_imp = None
        self.source_path = ""
        self.image_path_field.setText("")
        self.update_source_info()
        self.append_status("Reset complete. Ready for the next image.")

    @ui_action
    def open_last_output_folder(self, event=None):
        if self.last_run_dir == "" or not os.path.isdir(self.last_run_dir):
            JOptionPane.showMessageDialog(self, "No completed output folder is available yet.", APP_NAME, JOptionPane.WARNING_MESSAGE)
            return
        try:
            from java.awt import Desktop
            Desktop.getDesktop().open(File(self.last_run_dir))
        except ANY_ERROR as exc:
            self.report_exception("Could not open output folder", exc)


def launch_cell_quant_studio():
    """Open the window, or bring an already open one to the front.

    The window is built and shown on the Swing event thread, also when a
    script thread (the launcher) asks for it: building it on one thread while
    the event thread paints it can block both."""
    if not SwingUtilities.isEventDispatchThread():
        return call_on_ui(launch_cell_quant_studio)
    existing = None
    for window in Window.getWindows():
        try:
            if window.isDisplayable() and window.getTitle() == APP_NAME + " " + APP_VERSION:
                existing = window
                break
        except ANY_ERROR:
            pass

    if existing is not None:
        existing.setVisible(True)
        existing.toFront()
        existing.requestFocus()
        return existing
    return CellQuantStudioUniversal()


# The launcher loads this file to open the window. The headless batch script
# and the tests define CELLQUANT_NO_GUI = True first and use the functions only.
if not globals().get("CELLQUANT_NO_GUI", False):
    launch_cell_quant_studio()

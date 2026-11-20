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

from ij import IJ, WindowManager
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable, Measurements

from javax.swing import (
    JFrame, JPanel, JTabbedPane, JScrollPane, JTextArea, JButton, JLabel,
    JFileChooser, JOptionPane, JSplitPane, JTextField, JCheckBox, JComboBox,
    JProgressBar, BoxLayout, BorderFactory, SwingUtilities
)
from java.awt import (
    BorderLayout, GridBagLayout, GridBagConstraints, Insets, Dimension, Font,
    Color, FlowLayout, Cursor, Window
)
from java.awt.dnd import DropTarget, DnDConstants, DropTargetAdapter
from java.awt.datatransfer import DataFlavor
from java.io import File
from java.lang import Runnable

import os
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
        self.set_values(channel_values_from_properties(props, self.channel_number))


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

        self.build_ui()
        self.apply_preset_by_name(PRESET_NAMES[0])
        self.setVisible(True)

    # ----- UI construction ----------------------------------------------------

    def build_ui(self):
        root = JPanel(BorderLayout(8, 8))
        root.setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8))

        header = JPanel(BorderLayout())
        title = JLabel(APP_NAME + "  " + APP_VERSION)
        title.setFont(Font("SansSerif", Font.BOLD, 18))
        subtitle = JLabel("One Fiji application for channel preparation, independent masks, ROI measurements, object quantification, batches, and validation")
        subtitle.setFont(Font("SansSerif", Font.PLAIN, 12))
        header.add(title, BorderLayout.NORTH)
        header.add(subtitle, BorderLayout.SOUTH)
        root.add(header, BorderLayout.NORTH)

        self.main_tabs = JTabbedPane()
        self.main_tabs.add("Session", self.build_session_tab())
        self.main_tabs.add("Channels", self.build_channels_tab())
        self.main_tabs.add("ROIs and output", self.build_roi_output_tab())
        self.main_tabs.add("Batch", self.build_batch_tab())
        self.main_tabs.add("Validation and QC", self.build_validation_tab())
        self.main_tabs.add("Help", self.build_help_tab())
        root.add(self.main_tabs, BorderLayout.CENTER)

        self.getContentPane().add(root)

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

        actions = JPanel(FlowLayout(FlowLayout.LEFT))
        actions.setBorder(BorderFactory.createTitledBorder("Workflow actions"))
        actions.add(JButton("Open selected image", actionPerformed=self.open_selected_image))
        actions.add(JButton("Open via Bio-Formats", actionPerformed=self.open_selected_with_bioformats))
        actions.add(JButton("Use active image", actionPerformed=self.use_active_image))
        actions.add(JButton("1. Prepare enabled channels", actionPerformed=self.prepare_enabled_channels))
        actions.add(JButton("2. Preview all masks", actionPerformed=self.preview_all_masks))
        actions.add(JButton("3. Run analysis", actionPerformed=self.run_analysis_action))
        actions.add(JButton("Quick run", actionPerformed=self.quick_run_action))
        actions.add(JButton("Reset for next image", actionPerformed=self.reset_for_next_image))

        self.drop_area = JTextArea()
        self.drop_area.setEditable(False)
        self.drop_area.setLineWrap(True)
        self.drop_area.setWrapStyleWord(True)
        self.drop_area.setFont(Font("SansSerif", Font.PLAIN, 15))
        self.drop_area.setText(
            "Drop one microscopy image here.\n\n"
            "Recommended manual-threshold workflow:\n"
            "1. Open an image and apply a preset.\n"
            "2. Prepare the enabled channels.\n"
            "3. Use each channel's Brightness/Contrast button.\n"
            "4. Preview masks and adjust thresholds/morphology.\n"
            "5. Load or draw parent ROIs.\n"
            "6. Run analysis."
        )
        self.drop_area.setBorder(BorderFactory.createTitledBorder("Drag and drop"))
        DropTarget(self.drop_area, FileDropTarget(self))

        self.image_info_area = JTextArea()
        self.image_info_area.setEditable(False)
        self.image_info_area.setLineWrap(True)
        self.image_info_area.setWrapStyleWord(True)
        self.image_info_area.setText("No source image selected.")
        self.image_info_area.setBorder(BorderFactory.createTitledBorder("Source image"))

        center = JSplitPane(
            JSplitPane.HORIZONTAL_SPLIT,
            JScrollPane(self.drop_area),
            JScrollPane(self.image_info_area)
        )
        center.setResizeWeight(0.45)

        self.status_area = JTextArea()
        self.status_area.setEditable(False)
        self.status_area.setLineWrap(True)
        self.status_area.setWrapStyleWord(True)
        self.status_area.setFont(Font("Monospaced", Font.PLAIN, 11))
        self.status_area.setBorder(BorderFactory.createTitledBorder("Status and run log"))
        self.append_status("Ready.")

        top = JPanel(BorderLayout(6, 6))
        top.add(inputs, BorderLayout.NORTH)
        top.add(actions, BorderLayout.SOUTH)

        split = JSplitPane(JSplitPane.VERTICAL_SPLIT, center, JScrollPane(self.status_area))
        split.setResizeWeight(0.57)

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
        self.batch_progress.setPreferredSize(Dimension(420, 22))
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
            IJ.log(traceback.format_exc())
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
        """True (with a message) while a batch is running in the background."""
        if not self.batch_running:
            return False
        JOptionPane.showMessageDialog(
            self,
            "A batch is running. Wait for it to finish or cancel it on the Batch tab.",
            APP_NAME,
            JOptionPane.WARNING_MESSAGE
        )
        return True

    # ----- presets ------------------------------------------------------------

    def reset_channel_defaults(self):
        for panel in self.channel_panels:
            values = dict(CHANNEL_DEFAULTS)
            values["name"] = "Channel_" + str(panel.channel_number)
            panel.set_values(values)

    def apply_selected_preset(self, event=None):
        if self.blocked_by_batch():
            return
        self.apply_preset_by_name(to_text(self.preset_combo.getSelectedItem()))

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

    def choose_image_file(self, event=None):
        chooser = JFileChooser()
        chooser.setFileSelectionMode(JFileChooser.FILES_ONLY)
        if chooser.showOpenDialog(self) == JFileChooser.APPROVE_OPTION:
            path = chooser.getSelectedFile().getAbsolutePath()
            self.image_path_field.setText(path)
            self.open_image_path(path)

    def choose_output_dir(self, event=None):
        chooser = JFileChooser()
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
        if chooser.showOpenDialog(self) == JFileChooser.APPROVE_OPTION:
            self.output_dir_field.setText(chooser.getSelectedFile().getAbsolutePath())

    def open_selected_image(self, event=None):
        path = text_value(self.image_path_field)
        if path == "":
            JOptionPane.showMessageDialog(self, "Choose an image file first.", APP_NAME, JOptionPane.WARNING_MESSAGE)
            return
        self.open_image_path(path)

    def open_selected_with_bioformats(self, event=None):
        path = text_value(self.image_path_field)
        if path == "":
            JOptionPane.showMessageDialog(self, "Choose an image file first.", APP_NAME, JOptionPane.WARNING_MESSAGE)
            return
        self.open_image_path(path, force_bioformats=True)

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
        if self.blocked_by_batch():
            return
        imp = WindowManager.getCurrentImage()
        if imp is None:
            JOptionPane.showMessageDialog(self, "No active image is open.", APP_NAME, JOptionPane.WARNING_MESSAGE)
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

    def set_source_image(self, imp, path):
        self.close_generated_images()
        self.source_imp = imp
        self.source_path = path if path is not None else ""
        self.invalidate_prepared_cache()
        self.update_source_info()
        self.update_channel_availability()
        self.populate_calibration_fields()

    def update_source_info(self):
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
            self.channel_tabs.setEnabledAt(i, enabled)
            if not enabled:
                self.channel_panels[i].enabled_cb.setSelected(False)
                self.channel_panels[i].analysis_combo.setSelectedItem("Ignore")

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

    def collect_channel_configs(self):
        configs = []
        for panel in self.channel_panels:
            analysis_mode = to_text(panel.analysis_combo.getSelectedItem())
            if not panel.enabled_cb.isSelected() or analysis_mode == "Ignore":
                continue
            config = panel.get_config()
            configs.append(config)
        return configs

    def collect_all_channel_configs(self):
        return [panel.get_config() for panel in self.channel_panels]

    def collect_global_options(self):
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
        for config in configs:
            values.append((
                config["channel"], config["name"], config["projection"], config["z_slice"],
                config["subtract_background"], config["background_radius"], config["light_background"],
                config["median_radius"], config["gaussian_sigma"], config["lut"]
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
        return True

    def register_generated_image(self, imp):
        if imp is not None:
            self.generated_ids.add(imp.getID())

    def close_prepared_images(self):
        for imp in list(self.prepared_images.values()):
            if imp is not None:
                self.generated_ids.discard(imp.getID())
                safe_close_image(imp)
        self.prepared_images = {}
        self.invalidate_prepared_cache()

    def prepare_enabled_channels(self, event=None):
        self.set_busy(True)
        try:
            configs = self.collect_channel_configs()
            self.validate_configs(configs, for_analysis=False)
            frame = self.source_imp.getFrame()
            self.prepare_channels_for_frame(configs, frame, show_images=True, replace_cache=True)
            self.append_status("Prepared " + str(len(configs)) + " enabled channels for frame " + str(frame) + ".")
            if self.tile_prepared_cb.isSelected():
                IJ.run("Tile")
        except ANY_ERROR as exc:
            self.report_exception("Channel preparation failed", exc)
        finally:
            self.set_busy(False)

    def preview_prepared_channel(self, channel_number):
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

    def ensure_current_prepared(self, configs):
        frame = self.source_imp.getFrame()
        if not self.prepared_cache_valid(configs, frame):
            self.prepare_channels_for_frame(configs, frame, show_images=True, replace_cache=True)
        return self.prepared_images

    def open_brightness_contrast(self, channel_number):
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

    def capture_display_threshold(self, channel_number):
        try:
            configs = self.collect_channel_configs()
            self.validate_configs(configs, for_analysis=False)
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

    def run_analysis_action(self, event=None):
        if self.blocked_by_batch():
            return
        self.set_busy(True)
        try:
            self.run_analysis()
        except ANY_ERROR as exc:
            self.report_exception("Analysis failed", exc)
        finally:
            self.set_busy(False)
            IJ.showStatus("Cell Quant Studio ready")

    def run_analysis(self):
        configs = self.collect_channel_configs()
        self.validate_configs(configs, for_analysis=True)
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
                raise ValueError(
                    "Display-minimum thresholding requires the current prepared images. "
                    "Click '1. Prepare enabled channels', adjust Brightness/Contrast for " + labels +
                    ", preview the masks, and then run the analysis without changing preprocessing settings."
                )

        base = sanitize_name(strip_extension(self.source_imp.getTitle()))
        run_dir = unique_run_dir(output_root, base)
        self.last_run_dir = run_dir
        self.run_log_messages = []

        use_cache = (
            frame_mode == "Current frame only" and
            self.prepared_cache_valid(configs, frames[0])
        )
        close_generated = self.close_generated_after_cb.isSelected()
        retain_masks = self.keep_masks_open_cb.isSelected() and not close_generated

        result = self.engine.analyze_image(
            self.source_imp, self.source_path, configs, parent_rois, engine_options, run_dir, frames,
            prepared_cache=self.prepared_images if use_cache else None,
            settings=(global_options, all_configs),
            retain_masks=retain_masks
        )
        write_text_file(os.path.join(run_dir, "run_log.txt"), [APP_NAME + " " + APP_VERSION] + self.run_log_messages)

        if self.show_summary_cb.isSelected():
            self.show_summary_table(result["summary_headers"], result["summary_rows"])

        if close_generated:
            self.close_generated_images()
        else:
            for imp in result["mask_images"]:
                self.register_generated_image(imp)
                if WindowManager.getImage(imp.getID()) is None:
                    imp.show()
            if not self.keep_prepared_open_cb.isSelected():
                self.close_prepared_images()

        if self.close_source_cb.isSelected() and self.source_imp is not None:
            safe_close_image(self.source_imp)
            self.source_imp = None
            self.source_path = ""
            self.update_source_info()

        JOptionPane.showMessageDialog(
            self,
            "Analysis complete.\n\nSaved to:\n" + to_text(run_dir),
            APP_NAME,
            JOptionPane.INFORMATION_MESSAGE
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

    def load_profile(self, event=None):
        chooser = JFileChooser()
        if chooser.showOpenDialog(self) != JFileChooser.APPROVE_OPTION:
            return
        self.load_profile_from_path(chooser.getSelectedFile().getAbsolutePath())

    def load_profile_from_path(self, path):
        try:
            props = read_properties_file(path)
            loaded = global_options_from_properties(props)
            options = self.collect_global_options()
            for spec in GLOBAL_SETTING_SPECS:
                present = props.getProperty(spec[0]) is not None
                # A profile from 2.2.x has no batch keys: leave those as they are.
                if present or (spec[2] == "bool" and not spec[0].startswith("batch.")):
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
        except ANY_ERROR as exc:
            self.report_exception("Could not load settings profile", exc)

    # ----- batch --------------------------------------------------------------

    def choose_batch_input_dir(self, event=None):
        chooser = JFileChooser()
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
        if chooser.showOpenDialog(self) == JFileChooser.APPROVE_OPTION:
            self.batch_input_field.setText(chooser.getSelectedFile().getAbsolutePath())

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

    def run_batch_action(self, event=None):
        if self.blocked_by_batch():
            return
        try:
            configs = self.collect_channel_configs()
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

    def cancel_batch_action(self, event=None):
        if not self.batch_running:
            return
        self.batch_cancel_requested = True
        self.batch_cancel_button.setEnabled(False)
        self.append_status("Batch cancel requested; stopping after the current step.")

    def batch_finished(self, outcome):
        self.batch_running = False
        self.batch_thread = None
        self.batch_run_button.setEnabled(True)
        self.batch_cancel_button.setEnabled(False)
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

    def qc_context(self):
        """Everything a validation or sensitivity run needs: the channel under
        test, its prepared image for the current frame, and the parent ROIs."""
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
        prepared = self.ensure_current_prepared(configs)
        return {
            "configs": configs,
            "config": target,
            "prepared": prepared,
            "prepared_imp": prepared[channel_number],
            "selections": selections,
            "frame": self.source_imp.getFrame(),
            "source_title": to_text(self.source_imp.getTitle())
        }

    def run_validation_action(self, event=None):
        if self.blocked_by_batch():
            return
        self.set_busy(True)
        try:
            source_kind = to_text(self.annotation_source_combo.getSelectedItem())
            if source_kind == ANNOTATION_POINTS_IMAGE and to_text(self.roi_source_combo.getSelectedItem()) == "Current image ROI":
                raise ValueError(
                    "The image selection cannot be both the parent ROI and the point marks. "
                    "Add the parent ROI to ROI Manager (or use Whole image), or load the marks from a point ROI file."
                )
            context = self.qc_context()
            config = context["config"]
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
            parent_rois = self.get_parent_rois()
            output_root = self.get_output_root()
            base = sanitize_name(strip_extension(self.source_imp.getTitle()))
            run_dir = unique_dir(output_root, base + "_CellQuant_Validation_" + timestamp_text())
            self.last_run_dir = run_dir
            self.save_settings_to_path(os.path.join(run_dir, "settings.properties"))
            save_roi_array(parent_rois, os.path.join(run_dir, "Parent_ROIs_Input.zip"))
            if reference is not None:
                try:
                    result = validate_mask_against_reference(
                        self.engine, context["prepared_imp"], config, parent_rois, reference,
                        run_dir, context["source_title"], context["frame"]
                    )
                finally:
                    safe_close_image(reference)
                table_title = "CellQuant Mask Validation"
            else:
                result = validate_objects_against_points(
                    self.engine, context["prepared_imp"], config, parent_rois, points,
                    run_dir, context["source_title"], context["frame"]
                )
                table_title = "CellQuant Object Validation"
            self.show_table(table_title, result["headers"], result["rows"])
            self.append_status("Validation of C" + str(config["channel"]) + " (" + config["name"] + ") saved: " + to_text(run_dir))
            JOptionPane.showMessageDialog(
                self, "Validation complete.\n\nSaved to:\n" + to_text(run_dir), APP_NAME, JOptionPane.INFORMATION_MESSAGE
            )
        except ANY_ERROR as exc:
            self.report_exception("Validation failed", exc)
        finally:
            self.set_busy(False)
            IJ.showStatus("Cell Quant Studio ready")

    def run_sensitivity_action(self, event=None):
        if self.blocked_by_batch():
            return
        self.set_busy(True)
        try:
            parameter = to_text(self.sensitivity_parameter_combo.getSelectedItem())
            factors = parse_factor_list(text_value(self.sensitivity_factors_field))
            context = self.qc_context()
            config = context["config"]
            parent_rois = self.get_parent_rois()
            output_root = self.get_output_root()
            base = sanitize_name(strip_extension(self.source_imp.getTitle()))
            run_dir = unique_dir(output_root, base + "_CellQuant_Sensitivity_" + timestamp_text())
            self.last_run_dir = run_dir
            self.save_settings_to_path(os.path.join(run_dir, "settings.properties"))
            save_roi_array(parent_rois, os.path.join(run_dir, "Parent_ROIs_Input.zip"))
            result = run_sensitivity(
                self.engine, context["prepared_imp"], config, parent_rois, parameter, factors,
                run_dir, context["source_title"], context["frame"]
            )
            self.show_table("CellQuant Sensitivity", result["headers"], result["rows"])
            self.append_status(
                "Sensitivity of C" + str(config["channel"]) + " (" + config["name"] + ") to '" + parameter + "' saved: " + to_text(run_dir)
            )
            JOptionPane.showMessageDialog(
                self, "Sensitivity analysis complete.\n\nSaved to:\n" + to_text(run_dir), APP_NAME, JOptionPane.INFORMATION_MESSAGE
            )
        except ANY_ERROR as exc:
            self.report_exception("Sensitivity analysis failed", exc)
        finally:
            self.set_busy(False)
            IJ.showStatus("Cell Quant Studio ready")

    # ----- cleanup and navigation --------------------------------------------

    def close_generated_images(self):
        for image_id in list(self.generated_ids):
            imp = WindowManager.getImage(image_id)
            if imp is not None:
                safe_close_image(imp)
        self.generated_ids = set()
        self.prepared_images = {}
        self.mask_preview_ids = {}
        self.invalidate_prepared_cache()

    def close_generated_images_action(self, event=None):
        self.close_generated_images()
        self.append_status("Closed all images generated by this application.")

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
    """Open the window, or bring an already open one to the front."""
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

#@ String(label="Cell Quant Studio", value="Launch from Plugins > Cell Quant Studio") banner

from ij import IJ, WindowManager
from javax.swing import JFrame, JPanel, JTabbedPane, JScrollPane, JTextArea, JButton, JLabel, JFileChooser, JOptionPane, JSplitPane, JTextField
from javax.swing import BorderFactory
from java.awt import BorderLayout, GridBagLayout, GridBagConstraints, Insets, Dimension, Font
from java.awt.dnd import DropTarget, DnDConstants, DropTargetAdapter
from java.awt.datatransfer import DataFlavor
from java.io import File
import os

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
                    path = files[0].getAbsolutePath()
                    self.app.load_image_path(path)
                    event.dropComplete(True)
                    return
            event.dropComplete(False)
        except Exception as e:
            self.app.append_status("Drop failed: " + str(e))
            event.dropComplete(False)

class CellQuantStudio(JFrame):
    def __init__(self):
        JFrame.__init__(self, "Cell Quant Studio")
        self.setSize(1320, 900)
        self.setMinimumSize(Dimension(1080, 760))
        self.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE)
        self.fields = {}
        self.editors = {}
        self.status_area = None
        self.image_info_area = None
        self.drop_area = None
        self.detected_image_path = ""
        self.build_ui()
        self.autofill_macro_paths()
        self.load_all_editors_from_paths()

    def append_status(self, msg):
        if self.status_area is not None:
            self.status_area.append(msg + "\n")
            self.status_area.setCaretPosition(self.status_area.getDocument().getLength())

    def add_labeled_field(self, panel, row, key, label_text, default_value, browse_mode=None, load_editor_key=None):
        gbc = GridBagConstraints()
        gbc.gridx = 0
        gbc.gridy = row
        gbc.anchor = GridBagConstraints.WEST
        gbc.insets = Insets(4, 6, 4, 6)
        panel.add(JLabel(label_text), gbc)

        field = JTextField(str(default_value), 50)
        self.fields[key] = field

        gbc = GridBagConstraints()
        gbc.gridx = 1
        gbc.gridy = row
        gbc.weightx = 1.0
        gbc.fill = GridBagConstraints.HORIZONTAL
        gbc.insets = Insets(4, 6, 4, 6)
        panel.add(field, gbc)

        if browse_mode is not None:
            btn = JButton("Browse", actionPerformed=lambda e, mode=browse_mode, k=key, editor_key=load_editor_key: self.browse_into_field(k, mode, editor_key))
            gbc = GridBagConstraints()
            gbc.gridx = 2
            gbc.gridy = row
            gbc.insets = Insets(4, 6, 4, 6)
            panel.add(btn, gbc)

        if load_editor_key is not None:
            btn2 = JButton("Load", actionPerformed=lambda e, editor_key=load_editor_key: self.load_editor_from_field(editor_key))
            gbc = GridBagConstraints()
            gbc.gridx = 3
            gbc.gridy = row
            gbc.insets = Insets(4, 6, 4, 6)
            panel.add(btn2, gbc)

    def create_readonly_text(self, title, text, rows=6):
        area = JTextArea(rows, 40)
        area.setEditable(False)
        area.setLineWrap(True)
        area.setWrapStyleWord(True)
        area.setText(text)
        area.setBorder(BorderFactory.createTitledBorder(title))
        return area

    def build_ui(self):
        root = JPanel(BorderLayout(8, 8))
        root.setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8))

        header = self.create_readonly_text(
            "Overview",
            "Cell Quant Studio\n\n"
            "Use the Session tab for the workflow: drag in an image, detect the active image, choose an output folder, and run Prep, Quant, or the full sequence.\n\n"
            "Use the settings tabs to edit the actual macro files directly inside the UI."
        )
        root.add(header, BorderLayout.NORTH)

        tabs = JTabbedPane()
        tabs.add("Session", self.build_session_tab())
        tabs.add("Macro paths", self.build_macro_tab())
        tabs.add("Prep settings", self.build_settings_tab("prep_macro_path", "Prep macro editor", "Edit the full Prep macro here."))
        tabs.add("Quant settings", self.build_settings_tab("run_macro_path", "Quant macro editor", "Edit the full Quant macro here."))
        root.add(tabs, BorderLayout.CENTER)

        self.getContentPane().add(root)

    def build_session_tab(self):
        panel = JPanel(BorderLayout(8, 8))

        top = JPanel(GridBagLayout())
        top.setBorder(BorderFactory.createTitledBorder("Session input and output"))
        self.add_labeled_field(top, 0, "image_path", "Image file", "", "file")
        self.add_labeled_field(top, 1, "output_dir", "Output folder", "", "dir")

        actions = JPanel()
        actions.setBorder(BorderFactory.createTitledBorder("Session actions"))
        actions.add(JButton("Open selected image", actionPerformed=self.open_selected_image))
        actions.add(JButton("Detect active image", actionPerformed=self.detect_active_image))
        actions.add(JButton("Run Prep", actionPerformed=self.run_prep))
        actions.add(JButton("Run Quant", actionPerformed=self.run_run))
        actions.add(JButton("Run All", actionPerformed=self.run_all))
        actions.add(JButton("Save session", actionPerformed=self.save_session))
        actions.add(JButton("Load session", actionPerformed=self.load_session))

        self.drop_area = JTextArea()
        self.drop_area.setEditable(False)
        self.drop_area.setLineWrap(True)
        self.drop_area.setWrapStyleWord(True)
        self.drop_area.setText(
            "Drop a microscopy file here.\n\n"
            "Recommended order:\n"
            "1. Drag in one file\n"
            "2. Verify image detection on the right\n"
            "3. Choose an output folder\n"
            "4. Adjust macro settings in the Prep and Quant tabs if needed\n"
            "5. Run Prep, Quant, or Run All"
        )
        self.drop_area.setFont(Font("SansSerif", Font.PLAIN, 15))
        self.drop_area.setBorder(BorderFactory.createTitledBorder("Drag and drop zone"))
        DropTarget(self.drop_area, FileDropTarget(self))

        self.image_info_area = JTextArea()
        self.image_info_area.setEditable(False)
        self.image_info_area.setLineWrap(True)
        self.image_info_area.setWrapStyleWord(True)
        self.image_info_area.setText("No image detected yet.")
        self.image_info_area.setBorder(BorderFactory.createTitledBorder("Detected image information"))

        center_split = JSplitPane(JSplitPane.HORIZONTAL_SPLIT, JScrollPane(self.drop_area), JScrollPane(self.image_info_area))
        center_split.setResizeWeight(0.45)

        self.status_area = JTextArea()
        self.status_area.setEditable(False)
        self.status_area.setLineWrap(True)
        self.status_area.setWrapStyleWord(True)
        self.status_area.setBorder(BorderFactory.createTitledBorder("Session status"))
        self.append_status("Ready.")

        top_wrapper = JPanel(BorderLayout(8, 8))
        top_wrapper.add(top, BorderLayout.NORTH)
        top_wrapper.add(actions, BorderLayout.CENTER)

        main_split = JSplitPane(JSplitPane.VERTICAL_SPLIT, center_split, JScrollPane(self.status_area))
        main_split.setResizeWeight(0.62)

        panel.add(top_wrapper, BorderLayout.NORTH)
        panel.add(main_split, BorderLayout.CENTER)
        return panel

    def build_macro_tab(self):
        panel = JPanel(BorderLayout(8, 8))

        fields_panel = JPanel(GridBagLayout())
        fields_panel.setBorder(BorderFactory.createTitledBorder("Macro files"))
        self.add_labeled_field(fields_panel, 0, "prep_macro_path", "Prep macro", "", "file", "prep_macro_path")
        self.add_labeled_field(fields_panel, 1, "run_macro_path", "Quant macro", "", "file", "run_macro_path")

        buttons = JPanel()
        buttons.setBorder(BorderFactory.createTitledBorder("Macro file actions"))
        buttons.add(JButton("Load all editors from paths", actionPerformed=self.load_all_editors))
        buttons.add(JButton("Save all editors to paths", actionPerformed=self.save_all_editors))
        buttons.add(JButton("Autofill macro paths", actionPerformed=self.re_autofill_paths))

        note = self.create_readonly_text(
            "Notes",
            "The settings tabs are direct text editors for the macro files, so every setting exposed in the macro source can be changed from the UI."
        )

        panel.add(fields_panel, BorderLayout.NORTH)
        panel.add(buttons, BorderLayout.CENTER)
        panel.add(note, BorderLayout.SOUTH)
        return panel

    def build_settings_tab(self, macro_key, title, help_text):
        panel = JPanel(BorderLayout(8, 8))

        help_area = self.create_readonly_text("How to use this tab", help_text, 4)

        editor = JTextArea()
        editor.setLineWrap(False)
        editor.setWrapStyleWord(False)
        editor.setFont(Font("Monospaced", Font.PLAIN, 12))
        editor.setText("Open or load a macro file to edit its settings here.")
        editor.setBorder(BorderFactory.createTitledBorder(title))
        self.editors[macro_key] = editor

        btns = JPanel()
        btns.setBorder(BorderFactory.createTitledBorder("Editor actions"))
        btns.add(JButton("Load from path", actionPerformed=lambda e, k=macro_key: self.load_editor_from_field(k)))
        btns.add(JButton("Save to path", actionPerformed=lambda e, k=macro_key: self.save_editor_to_field(k)))
        btns.add(JButton("Save as...", actionPerformed=lambda e, k=macro_key: self.save_editor_as(k)))
        btns.add(JButton("Reload from disk", actionPerformed=lambda e, k=macro_key: self.reload_editor_from_disk(k)))

        panel.add(help_area, BorderLayout.NORTH)
        panel.add(JScrollPane(editor), BorderLayout.CENTER)
        panel.add(btns, BorderLayout.SOUTH)
        return panel

    def browse_into_field(self, key, mode, editor_key=None):
        chooser = JFileChooser()
        if mode == "dir":
            chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
        result = chooser.showOpenDialog(self)
        if result != JFileChooser.APPROVE_OPTION:
            return
        path = chooser.getSelectedFile().getAbsolutePath()
        self.fields[key].setText(path)
        if key == "image_path":
            self.load_image_path(path)
        elif editor_key is not None:
            self.load_editor_from_field(editor_key)

    def autofill_macro_paths(self):
        try:
            here = os.path.dirname(os.path.abspath(__file__))
        except Exception:
            here = ""
        if here == "":
            return

        candidates = {
            "prep_macro_path": ["CellQuant_01_Prep.ijm", "Prep.ijm"],
            "run_macro_path": ["CellQuant_02_Quant.ijm", "Run.ijm"]
        }
        for key in candidates:
            for name in candidates[key]:
                p = os.path.join(here, name)
                if os.path.exists(p):
                    self.fields[key].setText(p)
                    break

    def re_autofill_paths(self, event=None):
        self.autofill_macro_paths()
        self.append_status("Re-autofilled macro paths.")

    def load_editor_from_field(self, macro_key):
        path = self.fields[macro_key].getText().strip()
        if path == "":
            self.append_status("No macro path set for: " + macro_key)
            return
        if not os.path.exists(path):
            self.append_status("Macro file does not exist: " + path)
            return
        try:
            text = open(path, 'r').read()
            self.editors[macro_key].setText(text)
            self.editors[macro_key].setCaretPosition(0)
            self.append_status("Loaded editor from: " + path)
        except Exception as e:
            self.append_status("Could not load macro editor: " + str(e))
            JOptionPane.showMessageDialog(self, "Could not load macro file:\n" + str(e), "Cell Quant Studio", JOptionPane.ERROR_MESSAGE)

    def save_editor_to_field(self, macro_key):
        path = self.fields[macro_key].getText().strip()
        if path == "":
            self.append_status("No macro path set for: " + macro_key)
            return
        try:
            text = self.editors[macro_key].getText()
            f = open(path, 'w')
            f.write(text)
            f.close()
            self.append_status("Saved editor to: " + path)
        except Exception as e:
            self.append_status("Could not save macro editor: " + str(e))
            JOptionPane.showMessageDialog(self, "Could not save macro file:\n" + str(e), "Cell Quant Studio", JOptionPane.ERROR_MESSAGE)

    def save_editor_as(self, macro_key):
        chooser = JFileChooser()
        if chooser.showSaveDialog(self) != JFileChooser.APPROVE_OPTION:
            return
        path = chooser.getSelectedFile().getAbsolutePath()
        self.fields[macro_key].setText(path)
        self.save_editor_to_field(macro_key)

    def reload_editor_from_disk(self, macro_key):
        self.load_editor_from_field(macro_key)

    def load_all_editors(self, event=None):
        self.load_all_editors_from_paths()

    def load_all_editors_from_paths(self):
        for key in ["prep_macro_path", "run_macro_path"]:
            try:
                self.load_editor_from_field(key)
            except Exception:
                pass

    def save_all_editors(self, event=None):
        for key in ["prep_macro_path", "run_macro_path"]:
            self.save_editor_to_field(key)
        self.append_status("Saved all macro editors to their current paths.")

    def load_image_path(self, path):
        if path is None or path.strip() == "":
            return
        path = path.strip()
        self.fields["image_path"].setText(path)
        imp = IJ.openImage(path)
        if imp is None:
            self.append_status("Could not open image: " + path)
            JOptionPane.showMessageDialog(self, "Could not open image:\n" + path, "Cell Quant Studio", JOptionPane.ERROR_MESSAGE)
            return
        imp.show()
        self.detected_image_path = path
        self.append_status("Opened image: " + path)
        self.update_image_info(imp, path)

    def update_image_info(self, imp, path):
        msg = []
        msg.append("Path: " + path)
        msg.append("Title: " + imp.getTitle())
        msg.append("Width x Height: %d x %d" % (imp.getWidth(), imp.getHeight()))
        msg.append("Channels: %d" % imp.getNChannels())
        msg.append("Slices: %d" % imp.getNSlices())
        msg.append("Frames: %d" % imp.getNFrames())
        msg.append("Bit depth: %d" % imp.getBitDepth())
        try:
            cal = imp.getCalibration()
            if cal is not None:
                msg.append("Pixel size: %s x %s %s" % (str(cal.pixelWidth), str(cal.pixelHeight), cal.getUnit()))
                msg.append("Voxel depth: %s %s" % (str(cal.pixelDepth), cal.getUnit()))
        except Exception:
            pass
        if imp.getNChannels() == 4:
            msg.append("Status: 4-channel image detected.")
        else:
            msg.append("Status: channel count is not 4. Prep expects 4 channels.")
        self.image_info_area.setText("\n".join(msg))

    def open_selected_image(self, event=None):
        path = self.fields["image_path"].getText().strip()
        if path == "":
            self.append_status("No image path selected.")
            return
        self.load_image_path(path)

    def detect_active_image(self, event=None):
        imp = WindowManager.getCurrentImage()
        if imp is None:
            self.append_status("No active image found.")
            JOptionPane.showMessageDialog(self, "No active image found.", "Cell Quant Studio", JOptionPane.WARNING_MESSAGE)
            return
        info = imp.getOriginalFileInfo()
        path = ""
        try:
            if info is not None and info.directory is not None and info.fileName is not None:
                path = os.path.join(info.directory, info.fileName)
        except Exception:
            path = ""
        if path != "":
            self.fields["image_path"].setText(path)
            self.detected_image_path = path
        self.update_image_info(imp, path if path != "" else "<active image has no file path>")
        self.append_status("Detected active image: " + imp.getTitle())

    def ensure_output_dir_default(self):
        output_dir = self.fields["output_dir"].getText().strip()
        image_path = self.fields["image_path"].getText().strip()
        if output_dir == "" and image_path != "":
            try:
                self.fields["output_dir"].setText(os.path.dirname(image_path))
            except Exception:
                pass

    def build_prep_arg(self):
        self.ensure_output_dir_default()
        image_path = self.fields["image_path"].getText().strip().replace("\\", "/")
        output_dir = self.fields["output_dir"].getText().strip().replace("\\", "/")
        parts = []
        if image_path != "":
            parts.append("imagePath=" + image_path)
        if output_dir != "":
            parts.append("outputDir=" + output_dir)
        return "\n".join(parts)

    def build_output_arg(self):
        self.ensure_output_dir_default()
        output_dir = self.fields["output_dir"].getText().strip().replace("\\", "/")
        if output_dir == "":
            return ""
        return "outputDir=" + output_dir

    def run_macro_at(self, key, arg_text=None):
        path = self.fields[key].getText().strip()
        if path == "":
            self.append_status("No macro path set for: " + key)
            JOptionPane.showMessageDialog(self, "No macro path set for: " + key, "Cell Quant Studio", JOptionPane.WARNING_MESSAGE)
            return False
        if not File(path).exists():
            self.append_status("Macro file does not exist: " + path)
            JOptionPane.showMessageDialog(self, "Macro file does not exist:\n" + path, "Cell Quant Studio", JOptionPane.ERROR_MESSAGE)
            return False
        try:
            if arg_text is None:
                IJ.runMacroFile(path)
            else:
                IJ.runMacroFile(path, arg_text)
            self.append_status("Launched macro: " + path)
            return True
        except Exception as e:
            self.append_status("Macro launch failed: " + str(e))
            JOptionPane.showMessageDialog(self, "Macro launch failed:\n" + str(e), "Cell Quant Studio", JOptionPane.ERROR_MESSAGE)
            return False

    def run_prep(self, event=None):
        self.save_editor_to_field("prep_macro_path")
        arg_text = self.build_prep_arg()
        self.append_status("Running Prep with args:\n" + arg_text)
        self.run_macro_at("prep_macro_path", arg_text)

    def run_run(self, event=None):
        self.save_editor_to_field("run_macro_path")
        self.activate_first_prep_window()
        arg_text = self.build_output_arg()
        self.run_macro_at("run_macro_path", arg_text if arg_text != "" else None)

    def run_all(self, event=None):
        self.save_all_editors()
        prep_ok = self.run_macro_at("prep_macro_path", self.build_prep_arg())
        if not prep_ok:
            return
        self.activate_first_prep_window()
        run_arg = self.build_output_arg()
        run_ok = self.run_macro_at("run_macro_path", run_arg if run_arg != "" else None)
        if run_ok:
            self.append_status("Run All completed: Prep -> Quant")

    def activate_first_prep_window(self):
        try:
            ids = WindowManager.getIDList()
            if ids is None:
                return False
            for iid in ids:
                imp = WindowManager.getImage(iid)
                if imp is None:
                    continue
                title = imp.getTitle()
                if (title.endswith("_C1_AvgZ_Med") or title.endswith("_C2_AvgZ_Med") or title.endswith("_C3_AvgZ_Med") or title.endswith("_C4_AvgZ_Med") or
                    title.endswith("_C1_PREP") or title.endswith("_C2_PREP") or title.endswith("_C3_PREP") or title.endswith("_C4_PREP") or
                    "_C1_AvgZ_Med-" in title or "_C2_AvgZ_Med-" in title or "_C3_AvgZ_Med-" in title or "_C4_AvgZ_Med-" in title):
                    imp.show()
                    WindowManager.setCurrentWindow(imp.getWindow())
                    self.append_status("Activated prep window: " + title)
                    return True
            return False
        except Exception as e:
            self.append_status("Could not activate prep window: " + str(e))
            return False

    def save_session(self, event=None):
        chooser = JFileChooser()
        if chooser.showSaveDialog(self) != JFileChooser.APPROVE_OPTION:
            return
        path = chooser.getSelectedFile().getAbsolutePath()
        items = []
        for k in sorted(self.fields.keys()):
            items.append(k + "=" + self.fields[k].getText().strip())
        f = open(path, 'w')
        f.write("\n".join(items))
        f.close()
        self.append_status("Saved session: " + path)

    def load_session(self, event=None):
        chooser = JFileChooser()
        if chooser.showOpenDialog(self) != JFileChooser.APPROVE_OPTION:
            return
        path = chooser.getSelectedFile().getAbsolutePath()
        lines = open(path, 'r').read().splitlines()
        for line in lines:
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            if k in self.fields:
                self.fields[k].setText(v)
        self.load_all_editors_from_paths()
        self.append_status("Loaded session: " + path)

app = CellQuantStudio()
app.setVisible(True)

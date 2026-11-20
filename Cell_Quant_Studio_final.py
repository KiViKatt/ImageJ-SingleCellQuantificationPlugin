# Cell Quant Studio Universal - Fiji launcher
# Keep cellquantengine.py and cellquantcore.py in the same folder as this file.

from ij import IJ
import os

LAUNCHER_VERSION = "2.3.0"
PACKAGE_FILES = ["cellquantengine.py", "cellquantcore.py"]

print("Cell Quant Studio Universal " + LAUNCHER_VERSION + " launcher started")

folders = []

try:
    folders.append(os.path.dirname(os.path.abspath(__file__)))
except:
    pass

try:
    imagej_root = IJ.getDirectory("imagej")
    if imagej_root:
        folders.append(
            os.path.join(
                imagej_root,
                "plugins",
                "Scripts",
                "Cell Quant Studio"
            )
        )
except:
    pass

package_dir = None
for folder in folders:
    if not folder:
        continue
    complete = True
    for name in PACKAGE_FILES:
        if not os.path.isfile(os.path.join(folder, name)):
            complete = False
    if complete:
        package_dir = folder
        break

if package_dir is None:
    raise IOError(
        "Cell Quant Studio could not find cellquantengine.py and cellquantcore.py. "
        "Keep Cell_Quant_Studio_final.py, cellquantengine.py, and cellquantcore.py together in "
        "Fiji.app/plugins/Scripts/Cell Quant Studio/."
    )

# Important: execfile loads each source file directly rather than sending
# the whole application through Fiji's JSR-223 eval path. The application is
# split into an engine file and a window file so that each one stays well
# inside the size Jython can compile.
CELLQUANT_NO_GUI = False

engine_path = os.path.join(package_dir, PACKAGE_FILES[0])
try:
    print("Loading Cell Quant Studio engine: " + engine_path)
except:
    pass
execfile(engine_path, globals())

if APP_VERSION != LAUNCHER_VERSION:
    raise RuntimeError(
        "Cell Quant Studio: the launcher is version " + LAUNCHER_VERSION +
        " but cellquantengine.py is version " + APP_VERSION +
        ". Re-extract the release so that all three files come from one package."
    )

core_path = os.path.join(package_dir, PACKAGE_FILES[1])
try:
    print("Loading Cell Quant Studio window: " + core_path)
except:
    pass
execfile(core_path, globals())

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach
# Copyright (C) 2026 Maksym Taran

# Cell Quant Studio Universal - a synthetic test image for the landing-page demo
#
# Builds a 4-channel "organoid section" as an OME-TIFF with channel metadata.
# It is NOT real tissue: nuclei are drawn as blurred ovals that do not
# overlap, and a share of them, larger near the surface, is lit in each
# marker channel. Its numbers say nothing about biology.
#   C1 HOPX (ex 633), C2 PAX6 (ex 561), C3 DAPI (ex 405), C4 TBR2 (ex 488):
#   the channel order of examples/Demo_4Channel_HOPX_PAX6_DAPI_TBR2.properties.
# Scale as the files that profile was made for: 4084 x 4084, 0.2076 um/pixel.
#
#   bash tests/run_tests.sh <Fiji folder> "<absolute path>/site/record/make_demo_image.py"
#
# Environment variables:
#   CQ_DEMO_IMAGE    the file to write (required), for example
#                    "/data/organoids/Organoid section (synthetic).ome.tif"
#   CQ_DEMO_NUCLEI   how many nuclei to draw (default 2600)
#   CQ_DEMO_SIZE     width and height in pixels (default 4084)
# Give Java enough memory: CELLQUANT_JAVA_OPTS=-Xmx8g.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

import math
import os

from java.util import Random
from ij import ImagePlus
from ij.gui import OvalRoi, PolygonRoi, Roi
from ij.plugin.filter import EDM
from ij.process import FloatProcessor, ByteProcessor

SIZE = int(os.environ.get("CQ_DEMO_SIZE", "4084"))
PIXEL = 0.2076 * 4084.0 / SIZE
OUT = os.environ["CQ_DEMO_IMAGE"]
rnd = Random(20261120)

center = SIZE / 2.0
radius = 355.0 / PIXEL
xs = []
ys = []
for step in range(720):
    angle = 2 * math.pi * step / 720.0
    reach = radius * (1 + 0.05 * math.sin(3 * angle + 0.7) + 0.035 * math.sin(5 * angle + 2.1) + 0.02 * math.sin(9 * angle + 4.0))
    xs.append(int(round(center + reach * math.cos(angle))))
    ys.append(int(round(center + reach * math.sin(angle))))
tissue = PolygonRoi(xs, ys, len(xs), Roi.POLYGON)
mask = ByteProcessor(SIZE, SIZE)
mask.setValue(255)
mask.fill(tissue)
depth_map = EDM().makeFloatEDM(mask, 0, False)


def between(depth, points):
    """Linear interpolation through (depth in um, value) points."""
    if depth <= points[0][0]:
        return points[0][1]
    for index in range(1, len(points)):
        if depth <= points[index][0]:
            low = points[index - 1]
            high = points[index]
            return low[1] + (high[1] - low[1]) * (depth - low[0]) / float(high[0] - low[0])
    return points[-1][1]


# How likely a place is to hold a nucleus, and a nucleus to carry a marker,
# by depth from the tissue edge (um). Progenitor markers near the surface.
DENSITY = [(0, 1.0), (110, 1.0), (190, 0.72), (260, 0.5), (400, 0.42)]
MARKERS = {
    "HOPX": [(0, 0.16), (70, 0.30), (150, 0.26), (230, 0.08), (400, 0.04)],
    "PAX6": [(0, 0.56), (80, 0.50), (160, 0.26), (240, 0.08), (400, 0.05)],
    "TBR2": [(0, 0.30), (90, 0.36), (170, 0.18), (250, 0.05), (400, 0.03)]
}

# Nuclei that do not overlap, found by trying random places.
cell = int(60.0 * 0.2076 / PIXEL) + 1
grid = {}
nuclei = []
tries = 0
target = int(os.environ.get("CQ_DEMO_NUCLEI", "2600"))
while len(nuclei) < target and tries < 900000:
    tries += 1
    x = rnd.nextDouble() * SIZE
    y = rnd.nextDouble() * SIZE
    size = (3.5 + rnd.nextDouble() * 1.5) / PIXEL
    edge = depth_map.getf(int(x), int(y))
    if edge < size + 3:
        continue
    depth = edge * PIXEL
    if rnd.nextDouble() > between(depth, DENSITY):
        continue
    gx = int(x / cell)
    gy = int(y / cell)
    free = True
    for nx in range(gx - 1, gx + 2):
        for ny in range(gy - 1, gy + 2):
            for other in grid.get((nx, ny), ()):
                if (other[0] - x) ** 2 + (other[1] - y) ** 2 < (other[2] + size + 3.0) ** 2:
                    free = False
                    break
            if not free:
                break
        if not free:
            break
    if not free:
        continue
    squash = 0.82 + rnd.nextDouble() * 0.36
    entry = (x, y, size, squash, depth)
    nuclei.append(entry)
    grid.setdefault((gx, gy), []).append(entry)
print("nuclei placed: %d after %d tries" % (len(nuclei), tries))


def plane(kind):
    ip = FloatProcessor(SIZE, SIZE)
    # A faint haze inside the tissue, as autofluorescence gives.
    haze = 260.0 if kind == "DAPI" else 150.0
    ip.setValue(haze)
    ip.fill(tissue)
    ip.blurGaussian(14.0 * 0.2076 / PIXEL)
    lit = 0
    for x, y, size, squash, depth in nuclei:
        if kind == "DAPI":
            level = 2600.0 + rnd.nextDouble() * 2200.0
        else:
            if rnd.nextDouble() > between(depth, MARKERS[kind]):
                continue
            level = 1900.0 + rnd.nextDouble() * 1700.0
        lit += 1
        wide = 2 * size
        high = 2 * size * squash
        ip.setValue(level)
        ip.fill(OvalRoi(x - wide / 2.0, y - high / 2.0, wide, high))
    ip.resetRoi()
    ip.blurGaussian(2.4 * 0.2076 / PIXEL)
    ip.noise(22.0)
    ip.add(118.0)
    ip.min(0.0)
    print("%s: %d nuclei lit" % (kind, lit))
    return ip.convertToShort(False)


ORDER = ["HOPX", "PAX6", "DAPI", "TBR2"]
FACTS = {
    "HOPX": ("Alexa Fluor 647", 633.0, 668.0, (255, 0, 255)),
    "PAX6": ("Alexa Fluor 568", 561.0, 603.0, (255, 0, 0)),
    "DAPI": ("DAPI", 405.0, 461.0, (0, 0, 255)),
    "TBR2": ("Alexa Fluor 488", 488.0, 520.0, (0, 255, 0))
}
planes = [plane(name) for name in ORDER]

from loci.common import DataTools
from loci.formats import MetadataTools
from loci.formats.out import OMETiffWriter
from ome.units import UNITS
from ome.units.quantity import Length
from ome.xml.model.primitives import Color as OmeColor

meta = MetadataTools.createOMEXMLMetadata()
MetadataTools.populateMetadata(meta, 0, "Organoid section (synthetic)", False, "XYCZT", "uint16", SIZE, SIZE, 1, len(ORDER), 1, 1)
meta.setPixelsPhysicalSizeX(Length(PIXEL, UNITS.MICROMETER), 0)
meta.setPixelsPhysicalSizeY(Length(PIXEL, UNITS.MICROMETER), 0)
for index in range(len(ORDER)):
    fluor, excitation, emission, shade = FACTS[ORDER[index]]
    meta.setChannelName(fluor, 0, index)
    meta.setChannelFluor(fluor, 0, index)
    meta.setChannelExcitationWavelength(Length(excitation, UNITS.NANOMETER), 0, index)
    meta.setChannelEmissionWavelength(Length(emission, UNITS.NANOMETER), 0, index)
    meta.setChannelColor(OmeColor(shade[0], shade[1], shade[2], 255), 0, index)
if os.path.isfile(OUT):
    os.remove(OUT)
writer = OMETiffWriter()
writer.setMetadataRetrieve(meta)
writer.setCompression("LZW")
writer.setId(OUT)
for index in range(len(planes)):
    writer.saveBytes(index, DataTools.shortsToBytes(planes[index].getPixels(), False))
writer.close()
print("IMAGE-WRITTEN %s (%d x %d, %.4f um/pixel, %d nuclei)" % (OUT, SIZE, SIZE, PIXEL, len(nuclei)))

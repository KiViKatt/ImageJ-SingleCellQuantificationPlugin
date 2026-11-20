# Cell Quant Studio Universal - tests of opening files through Bio-Formats
#
# Loaded by run_tests.py, which provides the helpers used here. Both tests
# are reported as SKIPPED when Bio-Formats is not installed.
#
# No microscope file is needed. Two kinds of input are made on the spot:
#  - "fake" files: Bio-Formats synthesizes an image from an empty file whose
#    name ends in .fake and carries the image parameters, for example
#    "multi&series=3&sizeC=2&sizeX=256&sizeY=200&pixelType=uint16.fake".
#    A pixel of such an image holds its X position, apart from the first few
#    pixels of each row.
#  - a small multi-series OME-TIFF written with the Bio-Formats writer, whose
#    series differ in size and hold the disks of the other tests, so the
#    object counts are known.
#
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij import WindowManager
from java.awt import Window as AwtWindow

import time


def bioformats_or_skip():
    try:
        from loci.plugins import BF
        from loci.formats import ImageWriter
    except ANY_ERROR:
        raise SkipTest("Bio-Formats is not installed")


def fake_file(folder, name):
    path = os.path.join(folder, name)
    handle = open(path, "wb")
    handle.close()
    return path


def text_file(folder, name, text):
    path = os.path.join(folder, name)
    handle = open(path, "wb")
    try:
        handle.write(text)
    finally:
        handle.close()
    return path


def showing_window_count():
    return len([window for window in AwtWindow.getWindows() if window.isShowing()])


def can_be_renamed(path):
    """False while a reader still holds the file open (matters on Windows).
    A file that ImageJ could not read is released by Java's garbage
    collector, a moment after the engine asks for it; allow one second."""
    for attempt in range(6):
        try:
            os.rename(path, path + ".moved")
            os.rename(path + ".moved", path)
            return True
        except ANY_ERROR:
            time.sleep(0.2)
    return False


def bioformats_log_level():
    try:
        from org.slf4j import LoggerFactory
        logger = LoggerFactory.getLogger("loci")
        return (to_plain(logger.getLevel()), to_plain(logger.getEffectiveLevel()))
    except ANY_ERROR:
        return None


def to_plain(value):
    return None if value is None else value.toString()


def csv_table(path):
    """Read a CSV file written by the engine into a list of dictionaries."""
    rows = csv_rows(path)
    return [dict(zip(rows[0], values)) for values in rows[1:]]


def csv_headers(path):
    return csv_rows(path)[0]


def csv_rows(path):
    """Read a CSV file written by the engine into a list of lists."""
    text = read_utf8(path)
    rows = []
    row = []
    field = []
    quoted = False
    index = 0
    while index < len(text):
        ch = text[index]
        if quoted:
            if ch == '"':
                if index + 1 < len(text) and text[index + 1] == '"':
                    field.append('"')
                    index += 1
                else:
                    quoted = False
            else:
                field.append(ch)
        elif ch == '"':
            quoted = True
        elif ch == ",":
            row.append("".join(field))
            field = []
        elif ch == "\n":
            row.append("".join(field))
            field = []
            rows.append(row)
            row = []
        elif ch != "\r":
            field.append(ch)
        index += 1
    if len(field) > 0 or len(row) > 0:
        row.append("".join(field))
        rows.append(row)
    return rows


def small_plane(width, height, disks, background=100):
    ip = ShortProcessor(width, height)
    ip.setValue(float(background))
    ip.fill()
    for disk in disks:
        draw_disk(ip, disk[0], disk[1], DISK_RADIUS, 2000)
    return ip


def write_ome_tiff(path, series):
    """Write a multi-series OME-TIFF. series is a list of dictionaries with
    name, pixel (micrometers) and planes (one 16-bit processor per channel)."""
    from loci.formats import ImageWriter, MetadataTools
    from loci.common import DataTools
    from loci.common.services import ServiceFactory
    from loci.formats.services import OMEXMLService
    from ome.units import UNITS
    from ome.units.quantity import Length

    metadata = ServiceFactory().getInstance(OMEXMLService).createOMEXMLMetadata()
    for index in range(len(series)):
        planes = series[index]["planes"]
        MetadataTools.populateMetadata(
            metadata, index, series[index]["name"], True, "XYCZT", "uint16",
            planes[0].getWidth(), planes[0].getHeight(), 1, len(planes), 1, 1
        )
        size = Length(series[index]["pixel"], UNITS.MICROMETER)
        metadata.setPixelsPhysicalSizeX(size, index)
        metadata.setPixelsPhysicalSizeY(size, index)
    writer = ImageWriter()
    writer.setMetadataRetrieve(metadata)
    try:
        writer.setId(path)
        for index in range(len(series)):
            writer.setSeries(index)
            planes = series[index]["planes"]
            for plane in range(len(planes)):
                writer.saveBytes(plane, DataTools.shortsToBytes(planes[plane].getPixels(), True))
    finally:
        writer.close()


SMALL_DISKS = [(30, 30), (70, 30), (60, 70)]
FAKE_MULTI = "multi&series=3&sizeC=2&sizeZ=3&sizeT=1&sizeX=256&sizeY=200&pixelType=uint16.fake"


def make_slide(folder, name="slide.ome.tif"):
    """Three series: a full field, a smaller field at another pixel size,
    and the photograph of the slide label that slide scanners add."""
    path = os.path.join(folder, name)
    write_ome_tiff(path, [
        {"name": "field A", "pixel": PIXEL, "planes": [nuclei_plane(DISKS), marker_plane()]},
        {"name": "field B", "pixel": 0.25, "planes": [small_plane(120, 100, SMALL_DISKS), small_plane(120, 100, [(60, 50)])]},
        {"name": "label image", "pixel": 10.0, "planes": [small_plane(64, 48, [])]}
    ])
    return path


def test_bioformats_open(E, t, tmp):
    bioformats_or_skip()
    folder = os.path.join(tmp, "bioformats_open")
    os.makedirs(folder)
    notes = []

    def log(message):
        notes.append(E.to_text(message))

    def close_all(series):
        for label, imp in series:
            imp.close()
            imp.flush()

    def dimensions(imp):
        return (imp.getWidth(), imp.getHeight(), imp.getNChannels(), imp.getNSlices(), imp.getNFrames())

    windows_before = showing_window_count()
    images_before = WindowManager.getImageCount()
    level_before = bioformats_log_level()

    t.check(E.bioformats_importer_class("ImporterOptions") is not None,
            "the importer options class (package loci.plugins.in) can be reached from Jython")

    # One series.
    single = fake_file(folder, "single&sizeC=2&sizeZ=3&sizeT=1&sizeX=256&sizeY=200&pixelType=uint16.fake")
    series = E.open_image_series(single, False, False, log)
    t.equal(len(series), 1, "a single-series file opens as one image")
    t.equal(series[0][0], "", "a single image has no series label")
    t.equal(dimensions(series[0][1]), (256, 200, 2, 3, 1), "its dimensions are read")
    t.equal(series[0][1].getBitDepth(), 16, "its bit depth is read")
    t.check(series[0][1].getWindow() is None, "it is opened without a window")
    t.check(not series[0][1].getStack().isVirtual(), "its pixels are in memory, not a virtual stack")
    t.equal(series[0][1].getStack().getProcessor(1).get(200, 100), 200, "its pixel values are read")
    close_all(series)
    series = E.open_image_series(single, False, True, log)
    t.equal((len(series), series[0][0]), (1, ""), "'every series' on a single-series file still gives one unlabeled image")
    close_all(series)
    t.equal(len(notes), 0, "nothing is logged for a file that opens normally")

    # Several series.
    multi = fake_file(folder, FAKE_MULTI)
    series = E.open_image_series(multi, False, False, log)
    t.equal((len(series), series[0][0]), (1, ""), "without 'every series' only the first series is opened")
    t.check(len(notes) == 1 and "3 series" in notes[0], "and the log says that the file holds more series")
    close_all(series)
    del notes[:]
    series = E.open_image_series(multi, False, True, log)
    t.equal([label for label, imp in series], ["series_001", "series_002", "series_003"],
            "'every series' opens each series with a numbered label")
    t.equal([dimensions(imp) for label, imp in series], [(256, 200, 2, 3, 1)] * 3, "each series has its dimensions")
    t.equal(len(set([E.to_text(imp.getTitle()) for label, imp in series])), 3, "each series has its own title")
    t.equal(len(notes), 0, "'every series' needs no fallback and logs nothing")
    close_all(series)
    t.equal(len(E.open_with_bioformats_api(multi, True)), 3, "open_with_bioformats_api returns every series")
    t.equal(len(E.open_with_bioformats_api(multi)), 1, "open_with_bioformats_api returns the first series by default")

    # The series of a file are not read until their turn comes.
    planned = E.plan_image_series(multi, False, True, log)
    t.equal(len(planned), 3, "a multi-series file is listed without reading pixels")
    t.check(len([entry for entry in planned if entry.image is not None]) == 0, "no series is in memory before it is needed")
    t.equal([entry.name for entry in planned], ["multi", "multi 2", "multi 3"], "the series names of the file are kept")
    second = planned[1].open()
    t.equal(dimensions(second), (256, 200, 2, 3, 1), "one series can be opened by itself")
    t.check(planned[0].image is None and planned[2].image is None, "opening one series leaves the others unread")
    planned[1].release()
    t.check(planned[1].image is None, "a released series is no longer held")

    # Calibration, channels, Z and time.
    calibrated = fake_file(folder, "cal&series=2&sizeC=3&sizeZ=2&sizeT=2&sizeX=128&sizeY=96&pixelType=uint8"
                                   "&physicalSizeX=0.25&physicalSizeY=0.25&physicalSizeZ=1.5.fake")
    series = E.open_image_series(calibrated, False, True, log)
    t.equal(len(series), 2, "both calibrated series open")
    imp = series[1][1]
    t.equal(dimensions(imp), (128, 96, 3, 2, 2), "channels, Z slices and time points are told apart")
    t.equal(imp.getBitDepth(), 8, "an 8-bit series stays 8-bit")
    cal = imp.getCalibration()
    t.close(cal.pixelWidth, 0.25, "pixel width comes from the file")
    t.close(cal.pixelHeight, 0.25, "pixel height comes from the file")
    t.close(cal.pixelDepth, 1.5, "Z step comes from the file")
    t.check(E.is_micrometer_unit(cal.getUnit()), "the unit Bio-Formats reports is recognised as micrometers")
    close_all(series)
    plain = fake_file(folder, "plain&sizeC=1&sizeX=64&sizeY=64.fake")
    series = E.open_image_series(plain, False, False, log)
    t.check(not E.is_micrometer_unit(series[0][1].getCalibration().getUnit()), "an uncalibrated file is not taken for micrometers")
    close_all(series)

    # A native TIFF, by ImageJ and forced through Bio-Formats.
    native = os.path.join(folder, "native.tif")
    source = make_image("native.tif", slices=2)
    FileSaver(source).saveAsTiff(native)
    source.close()
    by_imagej = E.open_image_series(native, False, False, log)
    by_bioformats = E.open_image_series(native, True, False, log)
    t.equal((len(by_imagej), len(by_bioformats)), (1, 1), "a TIFF opens either way")
    a = by_imagej[0][1]
    b = by_bioformats[0][1]
    t.equal(dimensions(b), dimensions(a), "Bio-Formats and ImageJ agree on the dimensions of a TIFF")
    t.close(b.getCalibration().pixelWidth, a.getCalibration().pixelWidth, "and on its pixel size")
    t.check(E.is_micrometer_unit(a.getCalibration().getUnit()) and E.is_micrometer_unit(b.getCalibration().getUnit()),
            "and both units are recognised as micrometers")
    same = True
    for plane in range(1, a.getStackSize() + 1):
        pa = a.getStack().getProcessor(plane).getPixels()
        pb = b.getStack().getProcessor(plane).getPixels()
        if list(pa) != list(pb):
            same = False
    t.check(same, "and on every pixel")
    close_all(by_imagej)
    close_all(by_bioformats)
    t.equal(len(notes), 0, "none of this needed a log message")

    # Several series of different size, one of them a slide label.
    slide = make_slide(folder)
    planned = E.plan_image_series(slide, False, True, log)
    t.equal([entry.label for entry in planned], ["series_001", "series_002", "series_003"], "an OME-TIFF lists its series")
    t.equal([entry.skip_reason != "" for entry in planned], [False, False, True], "the label image is marked as not to be analyzed")
    t.check("label image" in planned[2].skip_reason, "and the reason names it")
    series = E.open_image_series(slide, False, True, log)
    t.equal([label for label, imp in series], ["series_001", "series_002"], "open_image_series leaves the label image out")
    t.equal([dimensions(imp) for label, imp in series], [(WIDTH, HEIGHT, 2, 1, 1), (120, 100, 2, 1, 1)],
            "series of different size are each opened at their own size")
    t.close(series[1][1].getCalibration().pixelWidth, 0.25, "and with their own pixel size")
    close_all(series)
    everything = E.plan_image_series(slide, False, True, log, False)
    t.equal(len([entry for entry in everything if entry.skip_reason != ""]), 0, "label images can be kept on request")
    series = E.open_image_series(slide, False, False, log)
    t.equal(len(series), 1, "an OME-TIFF without 'every series' opens its first series")
    t.equal(dimensions(series[0][1]), (WIDTH, HEIGHT, 2, 1, 1), "at the right size")
    t.check(len(notes) == 1 and "3 series" in notes[0], "and the log mentions the other series")
    close_all(series)
    del notes[:]
    flagged = fake_file(folder, "thumb&series=2&thumbnail=true&sizeC=2&sizeX=64&sizeY=48.fake")
    t.equal(len([entry for entry in E.plan_image_series(flagged, False, True, log) if entry.skip_reason != ""]), 2,
            "series that the file itself marks as thumbnails are not analyzed")

    # Files that cannot be read: nothing is returned, nothing is raised, and
    # the file is not left open.
    broken_czi = text_file(folder, "broken.czi", b"this is not an image")
    broken_tif = text_file(folder, "broken.tif", b"this is not an image")
    missing = os.path.join(folder, "missing.czi")
    del notes[:]
    t.equal(len(E.open_image_series(broken_czi, False, False, log)), 0, "a damaged proprietary file opens as nothing")
    t.check(len(notes) > 0 and "broken.czi" in notes[-1], "and the log gives the reason")
    t.equal(len(E.open_image_series(broken_czi, False, True, log)), 0, "also with 'every series'")
    t.equal(len(E.open_image_series(broken_tif, True, False, log)), 0, "a damaged TIFF forced through Bio-Formats opens as nothing")
    t.equal(len(E.open_image_series(broken_tif, False, False, log)), 0, "a damaged TIFF opens as nothing")
    t.equal(len(E.open_image_series(missing, False, False, log)), 0, "a file that does not exist opens as nothing")
    t.equal(len(E.open_image_series(broken_czi)), 0, "a damaged file needs no log function")
    for path in [broken_czi, broken_tif, single, multi, slide, native]:
        t.check(can_be_renamed(path), "no reader is left holding " + os.path.basename(path)[:12])

    t.equal(showing_window_count(), windows_before, "no window or dialog was opened")
    t.equal(WindowManager.getImageCount(), images_before, "no image was left open in ImageJ")
    t.equal(bioformats_log_level(), level_before, "the Bio-Formats log level is put back")


def test_bioformats_batch(E, t, tmp):
    bioformats_or_skip()
    folder = os.path.join(tmp, "bioformats_batch")
    os.makedirs(folder)
    out = os.path.join(tmp, "bioformats_batch_out")

    slide = make_slide(folder)
    multi = fake_file(folder, FAKE_MULTI)
    dotted = fake_file(folder, "exp.1&series=2&sizeC=2&sizeX=256&sizeY=200&pixelType=uint16.fake")
    one_channel = fake_file(folder, "onechannel&sizeC=1&sizeX=256&sizeY=200&pixelType=uint16.fake")
    three_channels = fake_file(folder, "three&sizeC=3&sizeX=256&sizeY=200&pixelType=uint16.fake")
    movie = fake_file(folder, "movie&sizeC=2&sizeT=3&sizeX=256&sizeY=200&pixelType=uint16.fake")
    text_file(folder, "broken.czi", b"this is not an image")

    # Thresholds that suit both kinds of input: the disks (2000) and the
    # marker rectangle (900 and 1500) of the OME-TIFF lie above them, and so
    # does the right-hand part of a fake image, whose values run from 0 to 255.
    configs = [nuclei_config(E, {"manual_low": 150}), marker_config(E, {"manual_low": 200})]
    options = E.default_global_options()
    for key in ["save_prepared", "save_masks", "save_overlays", "save_object_rois"]:
        options[key] = False
    batch_options = {"extensions": "tif, fake, czi", "recursive": False, "roi_policy": E.BATCH_ROI_WHOLE, "all_series": True}
    windows_before = showing_window_count()
    images_before = WindowManager.getImageCount()
    lines = []

    result = E.run_batch(folder, out, configs, configs, options, batch_options, log=lambda message: lines.append(E.to_text(message)))
    # slide 2 + multi 3 + exp.1 2 + three 1 + movie 1 analyzed; the label image skipped; one channel and broken failed.
    t.equal((result["analyzed"], result["skipped"], result["failed"]), (9, 1, 2),
            "every series of every file is analyzed, the label image skipped, two files fail")
    report = csv_table(os.path.join(result["batch_dir"], "batch_report.csv"))
    t.equal(len(report), 12, "the report has one row per series and per failed file")

    def rows_of(name, table=None):
        return [row for row in (report if table is None else table) if row["Image_File"] == name]

    slide_rows = rows_of("slide.ome.tif")
    t.equal([row["Series"] for row in slide_rows], ["series_001", "series_002", "series_003"], "the series of a file are reported in order")
    t.equal([row["Status"] for row in slide_rows], ["OK", "OK", "SKIPPED"], "the label image is reported as skipped")
    t.equal([row["Series_Name"] for row in slide_rows], ["field A", "field B", "label image"], "the report keeps the series names")
    t.equal([(row["Image_Width"], row["Image_Height"], row["Image_Channels"]) for row in slide_rows[:2]],
            [(str(WIDTH), str(HEIGHT), "2"), ("120", "100", "2")], "the report gives the size of each series")
    t.equal([float(row["Pixel_Width"]) for row in slide_rows[:2]], [PIXEL, 0.25], "and its pixel size")
    t.check(E.is_micrometer_unit(slide_rows[0]["Spatial_Unit"]), "and its unit")
    t.equal([row["Object_Rows"] for row in slide_rows[:2]], [str(len(DISKS)), str(len(SMALL_DISKS))],
            "the objects of each series are counted")
    t.check(slide_rows[0]["Run_Folder"].startswith("slide.ome_series_001_CellQuant_"), "the run folder is named after file and series")
    t.check(os.path.isfile(os.path.join(result["batch_dir"], slide_rows[1]["Run_Folder"], "roi_summary.csv")),
            "each series has its own run folder")

    summary = csv_table(os.path.join(result["batch_dir"], "batch_roi_summary.csv"))
    t.equal(len(summary), 9, "the combined summary has one row per analyzed series")
    sources = [row["Source_Image"] for row in summary]
    t.equal(len(set(sources)), 9, "every series has its own Source_Image")
    t.check("slide.ome_series_001.tif" in sources and "slide.ome_series_002.tif" in sources,
            "Source_Image names the file and the series")
    by_source = dict([(row["Source_Image"], row) for row in summary])
    field_b = by_source["slide.ome_series_002.tif"]
    t.equal(int(float(field_b["C1_DAPI_Object_Count"])), len(SMALL_DISKS), "the smaller series is measured at its own size")
    t.close(float(field_b["ROI_Area"]), 120 * 100 * 0.25 * 0.25, "and with its own calibration")
    t.close(float(by_source["slide.ome_series_001.tif"]["C2_MAP2_Positive_Area"]),
            MARKER_OUTER[2] * MARKER_OUTER[3] * PIXEL * PIXEL, "the marker area of the first series is right")
    t.close(float(field_b["C2_MAP2_Positive_Area"]), disk_pixel_count() * 0.25 * 0.25, "and that of the second")
    t.equal(len(set([row["Run_ID"] for row in summary])), 9, "and every series has its own run")

    dotted_rows = rows_of(os.path.basename(dotted))
    t.equal(len(dotted_rows), 2, "a file with a dot in its name gives one row per series")
    # The full name is used, or, where the result paths below the test's
    # temporary folder would pass the Windows limit, its shortened form,
    # whose code is computed from the full name (including what follows the dot).
    dotted_stem = "exp.1&series=2&sizeC=2&sizeX=256&sizeY=200&pixelType=uint16_series_002"
    dotted_folder = dotted_rows[1]["Run_Folder"]
    t.check(dotted_folder.startswith(dotted_stem + "_CellQuant_") or
            (dotted_folder.startswith(dotted_stem[:20]) and ("_" + E.name_code(dotted_stem) + "_CellQuant_") in dotted_folder),
            "and the dot does not cut its run folder name short")

    failed = rows_of(os.path.basename(one_channel))
    t.equal((len(failed), failed[0]["Status"]), (1, "FAILED"), "a file with fewer channels than the settings need fails")
    t.check("C2" in failed[0]["Message"], "and the message names the missing channel")
    t.equal(failed[0]["Image_Channels"], "1", "and the report shows how many channels it has")
    more = rows_of(os.path.basename(three_channels))
    t.equal(more[0]["Status"], "OK", "a file with more channels than the settings use is analyzed")
    t.check("3 channel" in more[0]["Message"], "with a note that its channel count differs from the first image")
    t.check(len([line for line in lines if "3 channel" in line]) == 1, "the note is also in the log")
    broken = rows_of("broken.czi")
    t.equal(broken[0]["Status"], "FAILED", "an unreadable file is reported as failed")
    t.check("Could not open" in broken[0]["Message"] and "Bio-Formats" in broken[0]["Message"], "with the reason Bio-Formats gave")
    t.equal(rows_of(os.path.basename(movie))[0]["Image_Frames"], "3", "the report shows the time points of a time series")

    t.equal(showing_window_count(), windows_before, "the batch opened no window")
    t.equal(WindowManager.getImageCount(), images_before, "the batch left no image open")
    for path in [slide, multi, movie]:
        t.check(can_be_renamed(path), "the batch left no reader holding " + os.path.basename(path)[:12])

    # Without "every series" each file gives one image.
    first_only = dict(batch_options)
    first_only["all_series"] = False
    del lines[:]
    result = E.run_batch(folder, out, configs, configs, options, first_only, file_list=[slide, multi],
                         log=lambda message: lines.append(E.to_text(message)))
    t.equal((result["analyzed"], result["skipped"], result["failed"]), (2, 0, 0), "without 'every series' one image per file is analyzed")
    summary = csv_table(os.path.join(result["batch_dir"], "batch_roi_summary.csv"))
    t.equal(sorted([row["Source_Image"] for row in summary]), sorted(["slide.ome.tif", FAKE_MULTI]),
            "and Source_Image is the plain file name")
    t.equal(len([line for line in lines if "only the first is analyzed" in line]), 2, "the log says that series were left out")

    # A time series with "All frames".
    all_frames = dict(options)
    all_frames["frame_mode"] = "All frames"
    result = E.run_batch(folder, out, configs, configs, all_frames, batch_options, file_list=[movie])
    summary = csv_table(os.path.join(result["batch_dir"], "batch_roi_summary.csv"))
    t.equal([row["Frame"] for row in summary], ["1", "2", "3"], "'All frames' analyzes every time point of a Bio-Formats image")

    # A known answer from a fake image: with a threshold of 100 the pixels
    # from X = 100 to 255 form one object.
    gradient = [nuclei_config(E, {"manual_low": 100, "watershed": False})]
    result = E.run_batch(folder, out, gradient, gradient, options, batch_options, file_list=[multi])
    summary = csv_table(os.path.join(result["batch_dir"], "batch_roi_summary.csv"))
    t.equal([int(float(row["C1_DAPI_Object_Count"])) for row in summary], [1, 1, 1], "each series of the fake file holds one object")
    t.equal([float(row["C1_DAPI_Positive_Area"]) for row in summary], [156.0 * 200.0] * 3, "of exactly the expected area")

    # ROI files: one for the file, and one for a single series.
    top_row = Roi(10, 10, 205, 45)
    top_row.setName("top_row")
    E.save_roi_array([top_row], os.path.join(folder, "slide_ROIs.zip"))
    one_disk = Roi(10, 10, 40, 40)
    one_disk.setName("one_disk")
    E.save_roi_array([one_disk], os.path.join(folder, "slide.ome_series_002_ROIs.zip"))
    with_rois = dict(batch_options)
    with_rois["roi_policy"] = E.BATCH_ROI_FILE_SKIP
    result = E.run_batch(folder, out, configs, configs, options, with_rois, file_list=[slide, multi])
    report = csv_table(os.path.join(result["batch_dir"], "batch_report.csv"))
    slide_rows = rows_of("slide.ome.tif", report)
    t.equal([row["Object_Rows"] for row in slide_rows[:2]], ["5", "1"],
            "a series uses its own ROI file when there is one, and the file's otherwise")
    t.check("slide.ome_series_002_ROIs.zip" in slide_rows[1]["ROI_Source"], "the report names the ROI file of the series")
    t.equal((result["analyzed"], result["skipped"]), (2, 4), "series without any ROI file are skipped, like the label image")


class OutputCollector(object):
    """Stands in for sys.stdout while the headless script prints its log."""

    def __init__(self):
        self.parts = []

    def write(self, text):
        self.parts.append(text)

    def flush(self):
        pass

    def text(self):
        return "".join([to_plain_text(part) for part in self.parts])


def to_plain_text(value):
    try:
        return value if isinstance(value, type(u"")) else value.decode("utf-8")
    except ANY_ERROR:
        return repr(value)


def run_headless_script(E, values):
    """Run cellquant_batch_headless.py the way Fiji's script service does:
    every #@ parameter is a global of the script, folders and files arrive
    as java.io.File. Returns (namespace, printed text, error or None)."""
    import sys
    from java.io import File
    path = os.path.join(E.CELLQUANT_PACKAGE_DIR, "cellquant_batch_headless.py")
    namespace = {"__file__": path, "__name__": "__main__"}
    for key in values:
        value = values[key]
        if key in ["input_dir", "profile", "output_dir"]:
            value = File(value)
        namespace[key] = value
    collector = OutputCollector()
    previous = sys.stdout
    sys.stdout = collector
    error = None
    try:
        try:
            E.run_source_file(path, namespace)
        except ANY_ERROR as exc:
            error = exc
    finally:
        sys.stdout = previous
    return (namespace, collector.text(), error)


def test_headless_script_modules(E, t, tmp):
    """The headless script loads the optional modules: a profile that
    switches the depth profile on gives the same depth tables as a direct
    engine run, and the batch ROI policies of modules can be chosen."""
    if not hasattr(E, "DEPTH_TABLE_NAME"):
        raise SkipTest("cellquantdepth.py is not beside the engine")
    folder = os.path.join(tmp, "headless_in")
    out = os.path.join(tmp, "headless_out")
    os.makedirs(folder)
    os.makedirs(out)
    images = {
        "a.tif": make_image("a.tif", slices=1),
        "b.tif": make_image("b.tif", slices=1, disks_per_frame=[DISKS[:6]])
    }
    for name in images:
        FileSaver(images[name]).saveAsTiff(os.path.join(folder, name))
    configs = [nuclei_config(E), marker_config(E)]
    options = E.default_global_options()
    for key in ["save_prepared", "save_masks", "save_overlays", "save_object_rois", "save_parent_rois"]:
        options[key] = False
    t.equal(options.get("depth_enabled"), False, "the depth profile is off by default")
    plain_profile = os.path.join(tmp, "plain.properties")
    E.write_settings_properties(plain_profile, options, configs)
    depth_options = dict(options)
    depth_options["depth_enabled"] = True
    depth_options["depth_bin_width"] = "10"
    # The parent ROI is the whole image, so the depth is measured from the
    # image border; otherwise there is no outline and only a note is written.
    depth_options["depth_image_border_is_surface"] = True
    depth_profile = os.path.join(tmp, "depth.properties")
    E.write_settings_properties(depth_profile, depth_options, configs)
    t.check("depth.enabled=true" in read_utf8(depth_profile).replace(" ", "").lower(),
            "the profile records that the depth profile is on")
    windows_before = showing_window_count()
    images_before = WindowManager.getImageCount()
    values = {
        "input_dir": folder, "profile": depth_profile, "output_dir": out,
        "extensions": "tif", "recursive": False, "force_bioformats": False, "all_series": False
    }

    # Depth profile on; no region choice given, so the whole image is used.
    namespace, printed, error = run_headless_script(E, values)
    t.check(error is None, "the headless script runs (" + E.to_text(error) + ")")
    if error is not None:
        return
    t.check("cellquantdepth.py" in printed and "Optional modules loaded" in printed,
            "the script loads the optional modules and says so")
    result = namespace["cellquant_result"]
    t.equal((result["analyzed"], result["skipped"], result["failed"]), (2, 0, 0), "both images are analyzed")
    batch_dir = result["batch_dir"]
    report = csv_table(os.path.join(batch_dir, "batch_report.csv"))
    t.equal([row["ROI_Source"] for row in report], ["Whole image", "Whole image"],
            "without a region choice the whole image is the parent ROI")
    combined_path = os.path.join(batch_dir, "batch_" + E.DEPTH_TABLE_NAME)
    t.check(os.path.isfile(combined_path), "the batch has a combined depth table")
    if not os.path.isfile(combined_path):
        return
    combined = csv_table(combined_path)
    t.equal(sorted(set([row["Source_Image"] for row in combined])), ["a.tif", "b.tif"],
            "the combined depth table holds both images")
    binned = [row for row in combined if row["Bin_Index"].strip() != ""]
    t.equal(sorted(set([row["Source_Image"] for row in binned])), ["a.tif", "b.tif"],
            "with depth bins for both images, not only notes")
    t.check(len(set([row["Bin_Index"] for row in binned])) >= 3, "the 120 x 100 micrometer images give several 10 micrometer bins")
    for row in report:
        t.check(os.path.isfile(os.path.join(batch_dir, row["Run_Folder"], E.DEPTH_TABLE_NAME)),
                "each run folder has its own depth table")
    summary_headers = csv_headers(os.path.join(batch_dir, "batch_roi_summary.csv"))
    plain_headers = E.AnalysisEngine().build_summary_headers(configs, options)
    depth_headers = E.depth_summary_headers(configs)
    t.equal(summary_headers[:len(plain_headers)], plain_headers, "the summary keeps its usual columns first")
    t.check(len(depth_headers) > 0 and len([h for h in depth_headers if h in summary_headers]) == len(depth_headers),
            "and has the depth columns added")

    # The same image analyzed directly by the engine gives the same depth rows.
    direct_root = os.path.join(tmp, "headless_direct")
    os.makedirs(direct_root)
    direct = run_image(E, images["a.tif"], configs, [E.whole_image_roi(images["a.tif"])], direct_root,
                       options_override=depth_options)
    direct_rows = csv_table(os.path.join(direct["run_dir"], E.DEPTH_TABLE_NAME))
    batch_rows = [row for row in combined if row["Source_Image"] == "a.tif"]
    ignore = ["Run_ID", "Batch_ID", "Image_Relative_Path"]

    def comparable(rows):
        return [sorted([(key, row[key]) for key in row if key not in ignore]) for row in rows]

    t.check(len(direct_rows) > 0, "the direct run has depth rows")
    t.equal(comparable(batch_rows), comparable(direct_rows),
            "the headless batch gives the same depth rows as a direct engine run")

    # Depth profile off: no depth table and no depth columns.
    plain_values = dict(values)
    plain_values["profile"] = plain_profile
    plain_values["roi_policy"] = "Whole image"
    namespace, printed, error = run_headless_script(E, plain_values)
    t.check(error is None, "the script runs with the depth profile off")
    if error is None:
        plain_dir = namespace["cellquant_result"]["batch_dir"]
        t.check(not os.path.exists(os.path.join(plain_dir, "batch_" + E.DEPTH_TABLE_NAME)),
                "with the depth profile off there is no depth table")
        t.equal(csv_headers(os.path.join(plain_dir, "batch_roi_summary.csv")),
                plain_headers + ["Batch_ID", "Image_Relative_Path"], "and the summary has only its usual columns")

    # Region choices, including those registered by modules.
    choose = namespace["cellquant_choose_policy"]
    t.equal(choose(None), E.BATCH_ROI_WHOLE, "a missing region choice means the whole image")
    t.equal(choose("  "), E.BATCH_ROI_WHOLE, "so does an empty one")
    t.equal(choose(E.BATCH_ROI_FILE_SKIP), E.BATCH_ROI_FILE_SKIP, "the ROI file policies are accepted")
    for bad, words in [(E.BATCH_ROI_MANAGER, "no ROI Manager"), ("Nonsense", "Unknown region choice")]:
        try:
            choose(bad)
            message = ""
        except ValueError as exc:
            message = E.to_text(exc)
        t.check(words in message and E.BATCH_ROI_WHOLE in message,
                "'" + bad + "' is refused with a message that lists the choices")
    before = sorted(os.listdir(out))
    bad_values = dict(values)
    bad_values["roi_policy"] = "Nonsense"
    namespace, printed, error = run_headless_script(E, bad_values)
    t.check(error is not None and "Unknown region choice" in E.to_text(error), "the script stops on an unknown region choice")
    t.equal(sorted(os.listdir(out)), before, "and writes nothing")
    if hasattr(E, "TISSUE_POLICY"):
        t.equal(choose(E.TISSUE_POLICY), E.TISSUE_POLICY, "the tissue module's policy can be chosen")
        tissue_values = dict(values)
        tissue_values["roi_policy"] = E.TISSUE_POLICY
        namespace, printed, error = run_headless_script(E, tissue_values)
        report_path = os.path.join(E.to_text(namespace.get("cellquant_result", {}).get("batch_dir", "")), "batch_report.csv")
        t.check(os.path.isfile(report_path), "a batch with the tissue policy runs (" + E.to_text(error) + ")")
        if os.path.isfile(report_path):
            rows = csv_table(report_path)
            t.equal(len(rows), 2, "and reports both images")
            t.check(len([row for row in rows if row["ROI_Source"].strip() != ""]) == 2,
                    "with the source of their parent ROIs")

    t.equal(showing_window_count(), windows_before, "the headless script opened no window")
    t.equal(WindowManager.getImageCount(), images_before, "and left no image open")
    for name in images:
        images[name].close()


TESTS.append(("Bio-Formats opens single and multi-series files", test_bioformats_open))
TESTS.append(("Bio-Formats folder batch with every series", test_bioformats_batch))
TESTS.append(("headless batch script loads the optional modules (depth profile)", test_headless_script_modules))

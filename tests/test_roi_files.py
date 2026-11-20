# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach
# Copyright (C) 2026 Maksym Taran

# Cell Quant Studio Universal - tests of ROI file reading
#
# Loaded by run_tests.py, which provides the helpers used here.
# The source is kept plain ASCII with no encoding declaration on purpose.

from ij.io import RoiEncoder


def test_single_roi_saved_as_zip(E, t, tmp):
    """ROI Manager writes a plain .roi file when a set holds one ROI, even
    when the file is named .zip. Version 2.2 saved its parent ROIs that way,
    so such a file must load as one ROI."""
    folder = os.path.join(tmp, "roi_files")
    os.makedirs(folder)
    roi = triangle_roi()

    disguised = os.path.join(folder, "Parent_ROIs_Input.zip")
    RoiEncoder.save(roi, disguised)
    t.equal(E.roi_file_kind(disguised), "roi", "a .roi file named .zip is recognised by its content")
    loaded = E.load_rois_from_file(disguised)
    t.equal(len(loaded), 1, "one ROI is read from a .roi file named .zip")
    t.equal(loaded[0].getBounds(), roi.getBounds(), "its outline is unchanged")
    t.check(loaded[0].getName() is not None and not E.to_text(loaded[0].getName()).lower().endswith(".zip"),
            "its name does not end in .zip")

    real_zip = os.path.join(folder, "set_ROIs.zip")
    E.save_roi_array([roi, Roi(5, 5, 20, 20)], real_zip)
    t.equal(E.roi_file_kind(real_zip), "zip", "a real archive is recognised")
    t.equal(len(E.load_rois_from_file(real_zip)), 2, "a real archive still loads every ROI")

    plain = os.path.join(folder, "one.roi")
    RoiEncoder.save(roi, plain)
    t.equal(len(E.load_rois_from_file(plain)), 1, "a plain .roi file still loads")

    garbage = os.path.join(folder, "broken_ROIs.zip")
    handle = open(garbage, "wb")
    try:
        handle.write(b"this is not a ROI file")
    finally:
        handle.close()
    raised = False
    try:
        E.load_rois_from_file(garbage)
    except ANY_ERROR:
        raised = True
    t.check(raised, "a damaged file raises instead of returning nothing")


TESTS.append(("a one-ROI set saved under a .zip name loads", test_single_roi_saved_as_zip))

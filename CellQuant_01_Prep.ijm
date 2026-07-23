// ==========================================================
// CellQuant_01_Prep.ijm
//
// Creates:
//   base_C1_MEAS
//   base_C2_MEAS
//   base_C3_MEAS
//   base_C4_MEAS
//   base_C3_SEG
// ==========================================================


// =========================
// USER SETTINGS
// =========================
closeSourceAfterPrep = 0;

// Channel 1
c1PrepMode = "Max";
c1Slice = 1;
c1DoSubtractBackground = 0;
c1BgRadius = 30;
c1MedianRadius = 1;
c1GaussianSigma = 0;

// Channel 2
c2PrepMode = "Max";
c2Slice = 1;
c2DoSubtractBackground = 0;
c2BgRadius = 30;
c2MedianRadius = 1;
c2GaussianSigma = 0;

// Channel 3 measurement image
c3PrepMode = "Slice";
c3Slice = 2;
c3DoSubtractBackground = 0;
c3BgRadius = 30;
c3MedianRadius = 1;
c3GaussianSigma = 0;

// Channel 3 segmentation settings
c3ThresholdMethod = "Otsu";
c3ThresholdMode = "dark";
c3ThresholdScale = 0.5;
c3UseManualThreshold = 0;
c3ManualLow = 800;
c3ManualHigh = 65535;

c3DoOpen = 0;
c3DoClose = 0;
c3DoErode = 0;
c3DoDilate = 0;
c3FillHoles = 0;
c3DoWatershed = 1;

// Channel 4
c4PrepMode = "Max";
c4Slice = 1;
c4DoSubtractBackground = 0;
c4BgRadius = 30;
c4MedianRadius = 1;
c4GaussianSigma = 0;


// =========================
// MAIN
// =========================
requires("1.53");

allArgs = getArgument();
imagePath = getArgValue(allArgs, "imagePath", "");

openedHere = 0;
sourceTitle = "";

if (imagePath != "") {
    open(imagePath);
    sourceTitle = getTitle();
    openedHere = 1;
} else {
    if (nImages == 0) exit("No image is open.");
    sourceTitle = getTitle();
}

selectWindow(sourceTitle);
Stack.getDimensions(imgW, imgH, nCh, nSl, nFr);

if (nCh != 4) {
    exit("This Prep macro expects a 4-channel image. Detected channels: " + nCh);
}

base = buildBaseName(sourceTitle);

// Build grayscale measurement outputs
c1MeasID = buildMeasurementImage(sourceTitle, base, 1, c1PrepMode, c1Slice, c1DoSubtractBackground, c1BgRadius, c1MedianRadius, c1GaussianSigma);
c2MeasID = buildMeasurementImage(sourceTitle, base, 2, c2PrepMode, c2Slice, c2DoSubtractBackground, c2BgRadius, c2MedianRadius, c2GaussianSigma);
c3MeasID = buildMeasurementImage(sourceTitle, base, 3, c3PrepMode, c3Slice, c3DoSubtractBackground, c3BgRadius, c3MedianRadius, c3GaussianSigma);
c4MeasID = buildMeasurementImage(sourceTitle, base, 4, c4PrepMode, c4Slice, c4DoSubtractBackground, c4BgRadius, c4MedianRadius, c4GaussianSigma);

// Build channel 3 segmentation from C3_MEAS
c3MeasTitle = base + "_C3_MEAS";
if (!isOpen(c3MeasID)) exit("Could not find " + c3MeasTitle);

// Remove only a stale internal work window from an interrupted prior run.
closeInternalWindow("__TMP_C3_SEG_WORK");

selectImage(c3MeasID);
run("Duplicate...", "title=__TMP_C3_SEG_WORK");
selectWindow("__TMP_C3_SEG_WORK");
segID = getImageID();

// Threshold
if (c3UseManualThreshold == 1) {
    setThreshold(c3ManualLow, c3ManualHigh);
} else {
    if (startsWith(toLowerCase(c3ThresholdMode), "dark"))
        setAutoThreshold(c3ThresholdMethod + " dark");
    else
        setAutoThreshold(c3ThresholdMethod);

    getThreshold(tLow, tHigh);
    tLow = tLow * c3ThresholdScale;
    if (tLow < 0) tLow = 0;
    setThreshold(tLow, tHigh);
}

// Convert to true binary mask
run("Convert to Mask");

// Force exact final segmentation title
rename(base + "_C3_SEG");

// Re-select exact segmentation image by ID in case Fiji changed focus
selectImage(segID);

// Optional morphology
if (c3DoOpen == 1) {
    run("Open");
    rename(base + "_C3_SEG");
    selectImage(segID);
}
if (c3DoClose == 1) {
    run("Close");
    rename(base + "_C3_SEG");
    selectImage(segID);
}
if (c3DoErode == 1) {
    run("Erode");
    rename(base + "_C3_SEG");
    selectImage(segID);
}
if (c3DoDilate == 1) {
    run("Dilate");
    rename(base + "_C3_SEG");
    selectImage(segID);
}
if (c3FillHoles == 1) {
    run("Fill Holes");
    rename(base + "_C3_SEG");
    selectImage(segID);
}
if (c3DoWatershed == 1) {
    run("Watershed");
    rename(base + "_C3_SEG");
    selectImage(segID);
}

if (closeSourceAfterPrep == 1 && openedHere == 1) {
    if (isOpen(sourceTitle)) {
        selectWindow(sourceTitle);
        close();
    }
}

showMessage("Prep complete",
    "Created:\n" +
    base + "_C1_MEAS\n" +
    base + "_C2_MEAS\n" +
    base + "_C3_MEAS\n" +
    base + "_C4_MEAS\n" +
    base + "_C3_SEG"
);


// =========================
// HELPERS
// =========================
function buildMeasurementImage(srcTitle, base, ch, prepMode, sliceNum, doBG, bgRadius, medRadius, gaussSigma) {
    selectWindow(srcTitle);

    tmpTitle = "__TMP_C" + ch;
    sliceTmpTitle = "__TMP_SLICE_C" + ch;

    closeInternalWindow(tmpTitle);
    closeInternalWindow(sliceTmpTitle);
    selectWindow(srcTitle);

    run("Duplicate...", "title=" + tmpTitle + " duplicate channels=" + ch);
    selectWindow(tmpTitle);

    Stack.getDimensions(w, h, c, z, t);

    if (prepMode == "Average") {
        if (z > 1) {
            oldID = getImageID();

            if (t > 1)
                run("Z Project...", "projection=[Average Intensity] all");
            else
                run("Z Project...", "projection=[Average Intensity]");

            newID = getImageID();
            if (newID == oldID)
                exit("Average Intensity Z projection failed for channel " + ch + ".");

            if (isOpen(oldID)) {
                selectImage(oldID);
                close();
            }
            selectImage(newID);
        }
    } else if (prepMode == "Max") {
        if (z > 1) {
            oldID = getImageID();

            if (t > 1)
                run("Z Project...", "projection=[Max Intensity] all");
            else
                run("Z Project...", "projection=[Max Intensity]");

            newID = getImageID();
            if (newID == oldID)
                exit("Max Intensity Z projection failed for channel " + ch + ".");

            if (isOpen(oldID)) {
                selectImage(oldID);
                close();
            }
            selectImage(newID);
        }
    } else {
        // Preserve the original behavior: every mode other than Average or Max
        // is treated as a single-Z-slice preparation.
        if (z > 1) {
            if (sliceNum < 1) sliceNum = 1;
            if (sliceNum > z) sliceNum = z;

            // setSlice() is correct for an ordinary Z stack. For a Z-T
            // hyperstack, explicitly select the requested Z at the current T.
            if (t > 1) {
                Stack.getPosition(currentC, currentZ, currentT);
                Stack.setPosition(1, sliceNum, currentT);
            } else {
                setSlice(sliceNum);
            }

            oldID = getImageID();
            run("Duplicate...", "title=" + sliceTmpTitle);
            newID = getImageID();

            if (newID == oldID)
                exit("Slice duplication failed for channel " + ch + ".");

            if (isOpen(oldID)) {
                selectImage(oldID);
                close();
            }
            selectImage(newID);
        }
    }

    finalTitle = base + "_C" + ch + "_MEAS";
    rename(finalTitle);
    finalID = getImageID();
    selectImage(finalID);

    if (doBG == 1)
        run("Subtract Background...", "rolling=" + bgRadius);

    if (medRadius > 0)
        run("Median...", "radius=" + medRadius);

    if (gaussSigma > 0)
        run("Gaussian Blur...", "sigma=" + gaussSigma);

    rename(finalTitle);
    selectImage(finalID);
    return finalID;
}

function closeInternalWindow(title) {
    if (isOpen(title)) {
        selectWindow(title);
        close();
    }
}

function buildBaseName(title) {
    s = replace(title, "\\", "/");
    lastSlash = lastIndexOfString(s, "/");
    if (lastSlash >= 0)
        s = substring(s, lastSlash + 1);

    dot = lastIndexOfString(s, ".");
    if (dot > 0)
        s = substring(s, 0, dot);

    return s;
}

function getArgValue(allArgs, key, defaultValue) {
    if (allArgs == "") return defaultValue;
    lines = split(allArgs, "\n");
    prefix = key + "=";
    for (i = 0; i < lines.length; i++) {
        if (startsWith(lines[i], prefix))
            return substring(lines[i], lengthOf(prefix));
    }
    return defaultValue;
}

function lastIndexOfString(s, token) {
    pos = -1;
    start = 0;
    while (true) {
        idx = indexOf(s, token, start);
        if (idx < 0) return pos;
        pos = idx;
        start = idx + 1;
    }
}

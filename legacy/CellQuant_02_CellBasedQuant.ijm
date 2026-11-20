// ==========================================================
// Cell-based quantification from channel 3 segmentation
//
// Required open images:
//   base_C1_MEAS
//   base_C2_MEAS
//   base_C3_MEAS
//   base_C4_MEAS
//   base_C3_SEG
//
// Required:
// - parent ROIs already loaded in ROI Manager
//
// Output:
// - one new self-contained run folder per execution
// ==========================================================


// =========================
// USER SETTINGS
// =========================
c1MeasTitle = "";
c2MeasTitle = "";
c3SegTitle  = "";
c3MeasTitle = "";
c4MeasTitle = "";

c3ParticleSizeMin = 8;
c3ParticleSizeMax = 200;
c3CircularityMin  = 0.00;
c3CircularityMax  = 1.00;
c3ExcludeEdgeParticles = 0;
c3ApplyWatershedToCrop = 0;

saveParentMaskImages = 1;

decimalPlaces = 3;
defaultOutputDir = "";


// =========================
// MAIN
// =========================
requires("1.53");

allArgs = getArgument();
rootOutputDir = getArgValue(allArgs, "outputDir", defaultOutputDir);

if (rootOutputDir == "") {
    rootOutputDir = getDirectory("Choose output folder");
    if (rootOutputDir == "") exit("No output folder selected.");
}
if (!endsWith(rootOutputDir, File.separator))
    rootOutputDir = rootOutputDir + File.separator;

if (!isOpen("ROI Manager")) exit("ROI Manager is not open.");
originalParentCount = roiManager("count");
if (originalParentCount < 1) exit("No parent ROIs found in ROI Manager.");

activeTitle = getTitle();
base = deriveBaseName(activeTitle);
safeBase = sanitizeFileName(base);
runTimestamp = "" + getTime();
runName = safeBase + "_" + runTimestamp;

// Create a fully self-contained run folder
outputDir = rootOutputDir + runName + File.separator;
File.makeDirectory(outputDir);
File.makeDirectory(outputDir + "C3_CellMasks");

// Save original parent ROI set immediately
parentRoiZip = outputDir + "Parent_ROIs_Input.zip";
roiManager("Save", parentRoiZip);

// Resolve open windows
if (c1MeasTitle == "") c1MeasTitle = resolveExactOrBase(base + "_C1_MEAS", base, "_C1_MEAS");
if (c2MeasTitle == "") c2MeasTitle = resolveExactOrBase(base + "_C2_MEAS", base, "_C2_MEAS");
if (c3SegTitle  == "") c3SegTitle  = resolveExactOrBase(base + "_C3_SEG",  base, "_C3_SEG");
if (c3MeasTitle == "") c3MeasTitle = resolveExactOrBase(base + "_C3_MEAS", base, "_C3_MEAS");
if (c4MeasTitle == "") c4MeasTitle = resolveExactOrBase(base + "_C4_MEAS", base, "_C4_MEAS");

if (c1MeasTitle == "") exit("Could not resolve channel 1 measurement window.");
if (c2MeasTitle == "") exit("Could not resolve channel 2 measurement window.");
if (c3SegTitle  == "") exit("Could not resolve channel 3 segmentation window.");
if (c3MeasTitle == "") exit("Could not resolve channel 3 measurement window.");
if (c4MeasTitle == "") exit("Could not resolve channel 4 measurement window.");

selectWindow(c1MeasTitle);
imgWidth = getWidth();
imgHeight = getHeight();
getVoxelSize(pixelWidth, pixelHeight, voxelDepth, pixelUnit);

// Results columns used later
run("Set Measurements...", "area mean integrated centroid perimeter shape feret redirect=None decimal=" + decimalPlaces);

// Output files
parentCsv = outputDir + "parent_roi_summary.csv";
cellCsv   = outputDir + "per_cell_summary.csv";

parentHeader =
"Run_Name,Base,Run_Timestamp,Parent_ROI_Index,Parent_ROI_Label,Parent_ROI_X,Parent_ROI_Y,Parent_ROI_Width,Parent_ROI_Height," +
"Image_Width,Image_Height,Pixel_Width,Pixel_Height,Pixel_Unit," +
"C1_Meas_Title,C2_Meas_Title,C3_Seg_Title,C3_Meas_Title,C4_Meas_Title," +
"C3_ParticleSizeMin,C3_ParticleSizeMax,C3_CircularityMin,C3_CircularityMax,C3_ExcludeEdgeParticles,C3_ApplyWatershedToCrop," +
"Detected_Cell_Count,C3_TotalParticleArea,C3_Mask_Path\n";

cellHeader =
"Run_Name,Base,Run_Timestamp,Parent_ROI_Index,Parent_ROI_Label,Cell_Index_In_Parent,Global_Cell_ID," +
"Parent_ROI_X,Parent_ROI_Y,Parent_ROI_Width,Parent_ROI_Height," +
"Cell_Bounds_X,Cell_Bounds_Y,Cell_Bounds_Width,Cell_Bounds_Height," +
"Centroid_X,Centroid_Y,Area,Perimeter,Circularity,Feret," +
"C1_Mean,C1_IntDen,C2_Mean,C2_IntDen,C3_Mean,C3_IntDen,C4_Mean,C4_IntDen\n";

File.saveString(parentHeader, parentCsv);
File.saveString(cellHeader, cellCsv);

globalCellID = 0;

// -------------------------
// Main loop over all saved parent ROIs
// -------------------------
for (r = 0; r < originalParentCount; r++) {

    // Always restore the original parent ROI set fresh for this iteration
    roiManager("Reset");
    roiManager("Open", parentRoiZip);

    if (roiManager("count") != originalParentCount) {
        exit("Failed to restore the original parent ROI set correctly.");
    }

    parentLabel = "ROI_" + d2(r + 1);

    // -------------------------
    // Get parent ROI bounds and coordinates on full image
    // -------------------------
    selectWindow(c3SegTitle);
    roiManager("Select", r);
    getSelectionBounds(parentX, parentY, parentW, parentH);
    getSelectionCoordinates(parentXP, parentYP);

    // -------------------------
    // Duplicate bounding box crop from C3 segmentation
    // -------------------------
    cropTitle = "__C3SEG_PARENT_" + (r + 1);
    if (isOpen(cropTitle)) {
        selectWindow(cropTitle);
        close();
    }

    selectWindow(c3SegTitle);
    makeRectangle(parentX, parentY, parentW, parentH);
    run("Duplicate...", "title=" + cropTitle);
    selectWindow(cropTitle);

    // Force proper binary state
    run("8-bit");
    setThreshold(1, 255);
    run("Convert to Mask");

    // Rebuild the actual parent ROI shape inside the crop and clear outside
    for (k = 0; k < parentXP.length; k++) {
        parentXP[k] = parentXP[k] - parentX;
        parentYP[k] = parentYP[k] - parentY;
    }
    makeSelection("polygon", parentXP, parentYP);
    run("Clear Outside");

    if (c3ApplyWatershedToCrop == 1)
        run("Watershed");

    // -------------------------
    // Analyze particles only inside the actual parent ROI
    // -------------------------
    rmCountBefore = roiManager("count");   // should equal originalParentCount here
    run("Clear Results");

    particleCmd = "size=" + c3ParticleSizeMin + "-" + c3ParticleSizeMax +
                  " circularity=" + c3CircularityMin + "-" + c3CircularityMax +
                  " show=Masks display clear add";

    if (c3ExcludeEdgeParticles == 1)
        particleCmd = particleCmd + " exclude";

    run("Analyze Particles...", particleCmd);

    detectedCellCount = nResults;
    c3TotalParticleArea = 0;
    for (i = 0; i < nResults; i++) {
        c3TotalParticleArea = c3TotalParticleArea + getResult("Area", i);
    }

    // Save parent mask if generated
    maskPath = "";
    maskTitle = "Mask of " + cropTitle;
    if (saveParentMaskImages == 1 && isOpen(maskTitle)) {
        selectWindow(maskTitle);
        maskPath = outputDir + "C3_CellMasks" + File.separator + safeBase + "_" + parentLabel + "_C3_mask.tif";
        saveAs("Tiff", maskPath);
    }

    // -------------------------
    // Convert local child ROIs back to full-image ROIs
    // -------------------------
    localStart = rmCountBefore;
    localEnd = roiManager("count") - 1;
    childCount = 0;

    for (localIdx = localStart; localIdx <= localEnd; localIdx++) {

        selectWindow(cropTitle);
        roiManager("Select", localIdx);
        getSelectionCoordinates(xp, yp);

        for (k = 0; k < xp.length; k++) {
            xp[k] = xp[k] + parentX;
            yp[k] = yp[k] + parentY;
        }

        // Add full-image child ROI
        selectWindow(c1MeasTitle);
        makeSelection("polygon", xp, yp);
        roiManager("Add");
        globalIdx = roiManager("count") - 1;

        childCount = childCount + 1;
        globalCellID = globalCellID + 1;

        // -------------------------
        // Shape from C3 segmentation
        // -------------------------
        run("Clear Results");
        selectWindow(c3SegTitle);
        roiManager("Select", globalIdx);
        run("Measure");

        cellArea = getResult("Area", 0);
        cellPerimeter = getResult("Perim.", 0);
        cellCircularity = getResult("Circ.", 0);
        cellFeret = getResult("Feret", 0);
        cellCentroidX = getResult("XM", 0);
        cellCentroidY = getResult("YM", 0);

        // Bounds from full-image ROI
        selectWindow(c1MeasTitle);
        roiManager("Select", globalIdx);
        getSelectionBounds(cellX, cellY, cellW, cellH);

        // -------------------------
        // C1 intensity
        // -------------------------
        run("Clear Results");
        selectWindow(c1MeasTitle);
        roiManager("Select", globalIdx);
        run("Measure");
        c1Mean = getResult("Mean", 0);
        c1IntDen = getResult("IntDen", 0);

        // -------------------------
        // C2 intensity
        // -------------------------
        run("Clear Results");
        selectWindow(c2MeasTitle);
        roiManager("Select", globalIdx);
        run("Measure");
        c2Mean = getResult("Mean", 0);
        c2IntDen = getResult("IntDen", 0);

        // -------------------------
        // C3 intensity
        // -------------------------
        run("Clear Results");
        selectWindow(c3MeasTitle);
        roiManager("Select", globalIdx);
        run("Measure");
        c3Mean = getResult("Mean", 0);
        c3IntDen = getResult("IntDen", 0);

        // -------------------------
        // C4 intensity
        // -------------------------
        run("Clear Results");
        selectWindow(c4MeasTitle);
        roiManager("Select", globalIdx);
        run("Measure");
        c4Mean = getResult("Mean", 0);
        c4IntDen = getResult("IntDen", 0);

        cellRow = "";
        cellRow = cellRow + q(runName) + ",";
        cellRow = cellRow + q(base) + ",";
        cellRow = cellRow + q(runTimestamp) + ",";
        cellRow = cellRow + q("" + (r + 1)) + ",";
        cellRow = cellRow + q(parentLabel) + ",";
        cellRow = cellRow + q("" + childCount) + ",";
        cellRow = cellRow + q("" + globalCellID) + ",";
        cellRow = cellRow + q("" + parentX) + ",";
        cellRow = cellRow + q("" + parentY) + ",";
        cellRow = cellRow + q("" + parentW) + ",";
        cellRow = cellRow + q("" + parentH) + ",";
        cellRow = cellRow + q("" + cellX) + ",";
        cellRow = cellRow + q("" + cellY) + ",";
        cellRow = cellRow + q("" + cellW) + ",";
        cellRow = cellRow + q("" + cellH) + ",";
        cellRow = cellRow + q("" + cellCentroidX) + ",";
        cellRow = cellRow + q("" + cellCentroidY) + ",";
        cellRow = cellRow + q("" + cellArea) + ",";
        cellRow = cellRow + q("" + cellPerimeter) + ",";
        cellRow = cellRow + q("" + cellCircularity) + ",";
        cellRow = cellRow + q("" + cellFeret) + ",";
        cellRow = cellRow + q("" + c1Mean) + ",";
        cellRow = cellRow + q("" + c1IntDen) + ",";
        cellRow = cellRow + q("" + c2Mean) + ",";
        cellRow = cellRow + q("" + c2IntDen) + ",";
        cellRow = cellRow + q("" + c3Mean) + ",";
        cellRow = cellRow + q("" + c3IntDen) + ",";
        cellRow = cellRow + q("" + c4Mean) + ",";
        cellRow = cellRow + q("" + c4IntDen) + "\n";

        File.append(cellRow, cellCsv);
    }

    parentRow = "";
    parentRow = parentRow + q(runName) + ",";
    parentRow = parentRow + q(base) + ",";
    parentRow = parentRow + q(runTimestamp) + ",";
    parentRow = parentRow + q("" + (r + 1)) + ",";
    parentRow = parentRow + q(parentLabel) + ",";
    parentRow = parentRow + q("" + parentX) + ",";
    parentRow = parentRow + q("" + parentY) + ",";
    parentRow = parentRow + q("" + parentW) + ",";
    parentRow = parentRow + q("" + parentH) + ",";
    parentRow = parentRow + q("" + imgWidth) + ",";
    parentRow = parentRow + q("" + imgHeight) + ",";
    parentRow = parentRow + q("" + pixelWidth) + ",";
    parentRow = parentRow + q("" + pixelHeight) + ",";
    parentRow = parentRow + q(pixelUnit) + ",";
    parentRow = parentRow + q(c1MeasTitle) + ",";
    parentRow = parentRow + q(c2MeasTitle) + ",";
    parentRow = parentRow + q(c3SegTitle) + ",";
    parentRow = parentRow + q(c3MeasTitle) + ",";
    parentRow = parentRow + q(c4MeasTitle) + ",";
    parentRow = parentRow + q("" + c3ParticleSizeMin) + ",";
    parentRow = parentRow + q("" + c3ParticleSizeMax) + ",";
    parentRow = parentRow + q("" + c3CircularityMin) + ",";
    parentRow = parentRow + q("" + c3CircularityMax) + ",";
    parentRow = parentRow + q("" + c3ExcludeEdgeParticles) + ",";
    parentRow = parentRow + q("" + c3ApplyWatershedToCrop) + ",";
    parentRow = parentRow + q("" + detectedCellCount) + ",";
    parentRow = parentRow + q("" + c3TotalParticleArea) + ",";
    parentRow = parentRow + q(maskPath) + "\n";

    File.append(parentRow, parentCsv);

    // Cleanup temp windows
    if (isOpen(cropTitle)) {
        selectWindow(cropTitle);
        close();
    }
    if (isOpen(maskTitle)) {
        selectWindow(maskTitle);
        close();
    }

    run("Clear Results");
}

// Restore the user's original parent ROI set at the end
roiManager("Reset");
roiManager("Open", parentRoiZip);

showMessage("Done",
    "Saved in:\n" + outputDir + "\n\n" +
    "Files:\n" +
    "Parent_ROIs_Input.zip\n" +
    "parent_roi_summary.csv\n" +
    "per_cell_summary.csv"
);


// =========================
// HELPERS
// =========================
function resolveExactOrBase(exactTitle, base, suffixTag) {
    titles = getList("image.titles");

    for (i = 0; i < titles.length; i++) {
        if (titles[i] == exactTitle) return titles[i];
    }

    for (i = 0; i < titles.length; i++) {
        if (startsWith(titles[i], base) && endsWith(titles[i], suffixTag)) return titles[i];
    }

    for (i = 0; i < titles.length; i++) {
        if (endsWith(titles[i], suffixTag)) return titles[i];
    }

    return "";
}

function deriveBaseName(title) {
    tags = newArray("_C1_MEAS", "_C2_MEAS", "_C3_MEAS", "_C4_MEAS", "_C3_SEG");
    cut = -1;
    for (i = 0; i < tags.length; i++) {
        p = indexOf(title, tags[i]);
        if (p >= 0) {
            if (cut < 0 || p < cut) cut = p;
        }
    }
    if (cut < 0) return title;
    return substring(title, 0, cut);
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

function sanitizeFileName(s) {
    s = replace(s, "\\", "_");
    s = replace(s, "/", "_");
    s = replace(s, ":", "_");
    s = replace(s, "*", "_");
    s = replace(s, "?", "_");
    s = replace(s, "\"", "_");
    s = replace(s, "<", "_");
    s = replace(s, ">", "_");
    s = replace(s, "|", "_");
    s = replace(s, " ", "_");
    return s;
}

function q(s) {
    s = "" + s;
    s = replace(s, "\"", "\"\"");
    return "\"" + s + "\"";
}

function d2(n) {
    if (n < 10) return "0" + n;
    return "" + n;
}

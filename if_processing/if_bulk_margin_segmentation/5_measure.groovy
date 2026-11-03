// ============================================================================
// Script 5: per-region measurements (area, perimeter, #cells) + export
// ----------------------------------------------------------------------------
// Reports, for each segmentation region:
//     - area       (mm^2)
//     - perimeter  (mm)
//     - number of cells (detections whose CENTROID falls inside the region)
//     - tumour+ cells and % (cells classified as the positive class)
//
// Cells are counted by CENTROID-IN-REGION, not by hierarchy nesting, so the
// count is independent of resolveHierarchy() quirks and the INFILTRATIVE
// MARGIN's cell count EXCLUDES the tumour core by construction: the margin ROI
// is the ring (envelope - core), so core centroids are geometrically outside it.
// Necrosis is a SUBREGION of the Tumour core (it sits inside it), so cells in
// necrosis also count toward Tumour; it is reported separately for reference.
// ============================================================================
resolveHierarchy()

def cal = getCurrentServer().getPixelCalibration()
double pw = cal.getPixelWidthMicrons()
double ph = cal.getPixelHeightMicrons()

// Positive (tumour) class used for the tumour+ column. Must match the class
// used in Script 2 and Script 3 ('Tumour' by default).
def posClass = getPathClass('Tumour')

// All detections, cached once with their centroids in image (pixel) coords.
def dets = getDetectionObjects()
def cellData = dets.collect { d ->
    def r = d.getROI()
    [x: r.getCentroidX(), y: r.getCentroidY(), pos: (d.getPathClass() == posClass)]
}
println "Total detections: " + dets.size()

// Order matters only for the printout; areas/perimeters/counts are independent.
def regionOrder = ['Tumour', 'Infiltrative margin', 'Non-tumour', 'Necrosis']

def rows = []
regionOrder.each { clsName ->
    def region = getAnnotationObjects().find { it.getPathClass() == getPathClass(clsName) }
    if (region == null) { println String.format("%-20s (none)", clsName); return }
    def roi = region.getROI()

    // Calibrated area (mm^2) and perimeter (mm).
    double areaMm2 = roi.getScaledArea(pw, ph) / 1e6
    double periMm  = roi.getScaledLength(pw, ph) / 1e3

    // Cells whose centroid falls inside this region ROI.
    int nCells = 0, nPos = 0
    for (c in cellData) {
        if (roi.contains(c.x, c.y)) {
            nCells++
            if (c.pos) nPos++
        }
    }
    double pctPos = nCells > 0 ? 100.0 * nPos / nCells : 0.0

    rows << [clsName, areaMm2, periMm, nCells, nPos, pctPos]
    println String.format("%-20s area=%.4f mm^2 | perimeter=%.3f mm | cells=%d | tumour+=%d (%.1f%%)",
        clsName, areaMm2, periMm, nCells, nPos, pctPos)
}

// Necrotic fraction of the tumour (necrosis area / core area).
def tumReg = getAnnotationObjects().find { it.getPathClass() == getPathClass('Tumour') }
def necReg = getAnnotationObjects().find { it.getPathClass() == getPathClass('Necrosis') }
if (tumReg != null && necReg != null) {
    double necA = necReg.getROI().getScaledArea(pw, ph)
    double tumA = tumReg.getROI().getScaledArea(pw, ph)   // filled core, INCLUDES the necrosis inside it
    double frac = tumA > 0 ? 100.0 * necA / tumA : 0
    println String.format("Necrotic fraction of tumour core: %.1f%% (necrosis %.4f / core %.4f mm^2)",
        frac, necA/1e6, tumA/1e6)
}

// ---- Export --------------------------------------------------------------
def name = getProjectEntry().getImageName().replaceAll('[^A-Za-z0-9._-]', '_')
def outDir = buildFilePath(PROJECT_BASE_DIR, 'exports'); mkdirs(outDir)

// Per-region summary CSV (area, perimeter, cell counts).
def csv = new StringBuilder()
csv.append("image,region,area_mm2,perimeter_mm,n_cells,n_tumour_pos,pct_tumour_pos\n")
rows.each { r ->
    csv.append(String.format("%s,%s,%.6f,%.6f,%d,%d,%.2f\n",
        name, r[0], r[1], r[2], r[3], r[4], r[5]))
}
def summaryPath = buildFilePath(outDir, name + '_region_summary.csv')
new File(summaryPath).text = csv.toString()

exportObjectsToGeoJson(getAnnotationObjects(),
    buildFilePath(outDir, name + '_regions.geojson'), 'FEATURE_COLLECTION')
saveDetectionMeasurements(buildFilePath(outDir, name + '_cells.csv'))

println "Exported region summary: " + summaryPath
println "Exported regions GeoJSON + per-cell measurements to: " + outDir
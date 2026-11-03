// ---- Script 4: density thresholds -> 3 regions ---------------------------
import qupath.lib.objects.PathObjects
import qupath.lib.roi.RoiTools
import qupath.lib.roi.RoiTools.CombineOp
import qupath.lib.roi.GeometryTools
import qupath.lib.regions.ImagePlane
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.operation.union.UnaryUnionOp
import qupath.lib.regions.RegionRequest
import qupath.lib.analysis.images.ContourTracing
import qupath.lib.analysis.images.SimpleImages
import java.util.Arrays

// --- Density-map lookup ----------------------------------------------------
// This MUST match Step 3:
//   densitymap_<current image name>_gw<Gaussian radius>
//
// Step 3 uses radiusMicrons = 100.0, so keep this at 100.0 unless you also
// change Step 3 and rebuild the density map.
double densityRadiusMicrons = 100.0
def densityMapName = 'densitymap_' +
        getCurrentImageNameWithoutExtension() +
        '_gw' +
        ((int) densityRadiusMicrons)

int densityChannel = 0          // single-class Step 3 map should use channel [0]

println "Current image: '${getCurrentImageName()}'"
println "Expected density map: '${densityMapName}'"

// --- Threshold mode -------------------------------------------------------
// 'absolute'   : T_high / T_low are raw density values (the original behaviour).
//                Use when staining/imaging is tightly controlled and you WANT to
//                detect genuine density differences between slides.
// 'percentile' : T_high / T_low are computed per slide as percentiles of THIS
//                slide's own tumour-density distribution, so "core" is always the
//                densest pctHigh% of the tumour regardless of the slide's overall
//                density. This is the cross-slide-consistent option: a biologically
//                equivalent tumour is classified the same way on a dense and a
//                sparse section. The percentiles are taken over density pixels
//                ABOVE minSampleDensity (see below) so empty background does not
//                drag them down.
def thresholdMode = 'percentile'   // 'percentile' | 'absolute'
double pctHigh = 20.0    // core   = top (100-pctHigh)% of tumour density  -> 80 = densest 20%
double pctLow  = 1.0    // envelope= top (100-pctLow)%                     -> 40 = densest 60%
double minSampleDensity = 50   // percentile pool = density pixels strictly above this
                                // (drop near-zero background so percentiles reflect tumour, not glass)

// Fallback / 'absolute'-mode thresholds (also used if percentile sampling fails):
double T_low  = 2.0
double T_high = 8.0
double minRegionArea = 5000.0   // um^2

// Artefact filter: keep only tumour foci that actually contain a dense CORE.
// Scattered low-density envelope specks (stray positive / autofluorescent
// cells) never build a core - dropping the coreless ones removes those islands
// automatically, without hand-tuning a size cutoff. Fully data-driven.
boolean requireCoreInFocus = true

// Minimum area of a whole tumour FOCUS (one envelope connected component,
// measured on the solid tumour footprint) to keep. This is the size filter for
// the small scattered patches: a focus below this is dropped entirely - its core,
// margin and necrosis go with it, so the region algebra stays coherent. Set 0 to
// keep every focus. NOTE: this is distinct from minRegionArea, which only drops
// tiny density-map speckle fragments up front; minFocusArea culls whole foci.
double minFocusArea = 200000.0   // um^2  (0.2 mm^2). Raise to keep only large tumours.

// Minimum area of a dense CORE (T_high connected component) to count as a real
// core. A tiny core speck otherwise anchors a whole envelope focus (via
// requireCoreInFocus), and since margin = envelope - core, almost the entire
// focus then becomes Infiltrative margin. Dropping sub-threshold core specks
// removes those spurious margin-only patches: with requireCoreInFocus on, a
// focus whose only core was a speck now has no qualifying core and is dropped.
double minCoreArea = 50000.0   // um^2. Set 0 to keep every core piece.

// --- Tissue mask mode -----------------------------------------------------
// 'image' : tissue = whole image. FAST, no cell geometry, low memory. Use this
//           when tissue fills most of the slide (little empty glass). Non-tumour
//           will then include any large empty background - usually fine, and you
//           can crop it out with a manual tissue annotation if needed (see below).
// 'cells' : tissue = footprint derived from detected cells (excludes glass).
//           Heavier - buffers every cell centroid; can exhaust memory on slides
//           with hundreds of thousands of cells. Only use when glass matters.
def tissueMode = 'image'
double tissueBuffer = 40.0      // um - ('cells' mode) disk radius around each centroid
boolean fillTissueHoles = true  // ('cells' mode) fill interior holes in the tissue mask

// --- Necrosis handling ----------------------------------------------------
// A necrotic core has few viable cells -> low tumour-cell density -> the
// density contours would punch a hole in the middle of the tumour. Filling
// enclosed holes in the core keeps the tumour solid. Optionally label those
// filled holes as their own 'Necrosis' region (only holes >= minNecrosisArea,
// so small gaps between tumour nests are NOT called necrosis).
boolean fillNecrosis     = true    // absorb enclosed low-density holes into the tumour core
boolean classifyNecrosis = true    // additionally label large filled holes as 'Necrosis'
double  minNecrosisArea  = 20000.0 // um^2 - minimum enclosed-hole area to count as necrosis
// The core (T_high) ring is often BROKEN (an open C, not a closed O) - e.g. when
// T_high is high, or the dense rim thins on one side. A broken ring has no
// enclosed hole, so the necrotic centre is missed (no 'Necrosis' label). Before
// extracting the hole we morphologically CLOSE the core - bridging rim gaps up
// to ~necrosisBridge wide - so the C becomes an O and the centre is detected.
// Raise it if the necrosis is still not labelled (gap wider than this); lower it
// if separate tumour nests get wrongly merged. Keep it well below the necrosis
// radius (closing also fills genuine holes narrower than ~necrosisBridge).
double  necrosisBridge   = 1000.0   // um - max core-ring gap width to bridge (0 = off)

def imageData = getCurrentImageData()
def plane = ImagePlane.getDefaultPlane()
def cal = getCurrentServer().getPixelCalibration()
double pxW = cal.getPixelWidthMicrons()
double px2 = pxW * cal.getPixelHeightMicrons()

// --- Load and validate the saved density map ------------------------------
// Density maps are project resources. Stop immediately if the expected map
// cannot be found; the rest of this script cannot run without it.
if (getProject() == null) {
    println "ERROR: No QuPath project is open."
    println "Add the image to a project, run Step 3, then run this script again."
    return
}

def densityMapBuilder
try {
    densityMapBuilder = loadDensityMap(densityMapName)
    def dserver = densityMapBuilder.buildServer(imageData)

    println "Successfully loaded density map '${densityMapName}'."
    println "Density-map channels:"
    dserver.getMetadata().getChannels().eachWithIndex { ch, i ->
        println "  [${i}] ${ch.getName()}"
    }

    int nChannels = dserver.getMetadata().getChannels().size()
    if (densityChannel < 0 || densityChannel >= nChannels) {
        println "ERROR: densityChannel=${densityChannel}, but the map has ${nChannels} channel(s)."
        return
    }
} catch (Exception e) {
    println ""
    println "ERROR: Unable to load density map '${densityMapName}'."
    println "Reason: ${e.getMessage()}"
    println ""
    println "Run Step 3 successfully for this same image in this same project."
    println "Step 3 should print:"
    println "  Saved density map '${densityMapName}'"
    println ""
    println "If Step 3 used a different Gaussian radius, change"
    println "densityRadiusMicrons at the top of this script to match it."
    return
}

// Density-map thresholding should run across the whole image here.
// Script 1 ends with selectDetections(), which can leave thousands of detections
// selected. Clear them so QuPath does not try to use detections as parent regions.
resetSelection()
println "Cleared object selection before density-map thresholding."

// Warn about an extremely broad envelope threshold.
if (thresholdMode == 'percentile' && pctLow <= 5.0) {
    println String.format(
        "WARNING: pctLow=%.1f keeps approximately %.1f%% of sampled density pixels in the envelope; this may create very large/fragmented geometry and run slowly.",
        pctLow, 100.0 - pctLow
    )
}

// --- STEP A: remove prior region annotations FIRST ------------------------
// (Must happen before we pick the tissue parent, or a leftover Tumour region
//  from a previous run gets mistaken for 'tissue' and Non-tumour comes out empty.)
def regionClasses = [getPathClass('Tumour'), getPathClass('Infiltrative margin'),
                     getPathClass('Non-tumour'), getPathClass('Necrosis')]
removeObjects(getAnnotationObjects().findAll { it.getPathClass() in regionClasses }, true)

// --- STEP B: define the TISSUE mask ---------------------------------------
def cells = getDetectionObjects().findAll { it.getROI() != null }
if (cells.isEmpty()) { println "No detections - run Script 1 first."; return }
def gf = cells[0].getROI().getGeometry().getFactory()
def tissueRoi

if (tissueMode == 'image') {
    // FAST path: whole image is the tissue. Or, if you have drawn a manual
    // tissue annotation and classified it 'Tissue', that is used instead.
    def manual = getAnnotationObjects().find { it.getPathClass() == getPathClass('Tissue') && it.getROI() != null }
    if (manual != null) {
        tissueRoi = manual.getROI()
        println "Tissue mask: manual 'Tissue' annotation."
    } else {
        def srv = getCurrentServer()
        tissueRoi = qupath.lib.roi.ROIs.createRectangleROI(0, 0, srv.getWidth(), srv.getHeight(), plane)
        println "Tissue mask: whole image."
    }
} else {
    // 'cells' path: footprint from cell CENTROIDS (buffer a MultiPoint once).
    // Heavier and memory-hungry on very large cell counts.
    def coords = cells.collect {
        new org.locationtech.jts.geom.Coordinate(it.getROI().getCentroidX(), it.getROI().getCentroidY())
    } as org.locationtech.jts.geom.Coordinate[]
    def mp = gf.createMultiPointFromCoords(coords)
    double bufPx = tissueBuffer / pxW
    int quadSeg = 2   // very coarse circle -> fewer vertices, less memory
    def tissueGeom = mp.buffer(bufPx, quadSeg).buffer(-bufPx * 0.4, quadSeg)
    if (fillTissueHoles) {
        def polys = []
        for (int i = 0; i < tissueGeom.getNumGeometries(); i++) {
            def g = tissueGeom.getGeometryN(i)
            if (g instanceof Polygon) polys << gf.createPolygon(((Polygon) g).getExteriorRing().getCoordinates())
        }
        if (!polys.isEmpty()) tissueGeom = UnaryUnionOp.union(gf.createGeometryCollection(polys as Geometry[]))
    }
    tissueRoi = GeometryTools.geometryToROI(tissueGeom, plane)
}
println String.format("Tissue mask area: %.3f mm^2", tissueRoi.getArea() * px2 / 1e6)

// --- STEP C: envelope / core contours from the density map ----------------
// NOTE the parentheses around the key: [(var): value] uses the variable's
// value as the key. The key must be the Integer channel index.
//
// Performance: a low threshold produces MANY speckle polygons. We (1) drop
// fragments below minRegionArea up front, (2) union what's left in ONE
// UnaryUnionOp call (not pairwise combineROIs, which is ~quadratic), and
// (3) simplify the jagged contour so downstream boolean ops are cheap.
double simplifyPx = 2.0   // Douglas-Peucker tolerance in pixels for contour smoothing
double smoothRadius = 30.0 // um - morphological round-off of the blocky density-grid edges (0 = off)
def envelopeRoiAt = { double thr, String label ->
    long t0 = System.currentTimeMillis()
    println String.format("Creating %s contour at threshold %.3f...", label, thr)

    // Be defensive: generated annotations or previous scripts may change selection.
    resetSelection()
    createAnnotationsFromDensityMap(densityMapName, [(densityChannel): thr], 'tmpDensity')

    def tmp = getAnnotationObjects().findAll {
        it.getPathClass() == getPathClass('tmpDensity')
    }
    println "  ${label}: density threshold produced ${tmp.size()} temporary annotation(s)."

    def geoms = tmp.collect { it.getROI()?.getGeometry() }
                   .findAll { it != null && it.getArea() * px2 >= minRegionArea }

    removeObjects(tmp, true)

    if (geoms.isEmpty()) {
        println "  ${label}: no components remained after the area filter."
        return null
    }

    println "  ${label}: unioning ${geoms.size()} retained component(s)..."
    def u = UnaryUnionOp.union(gf.createGeometryCollection(geoms as Geometry[]))
    u = org.locationtech.jts.simplify.DouglasPeuckerSimplifier.simplify(u, simplifyPx)

    println String.format(
        "  %s contour complete in %.1f s.",
        label,
        (System.currentTimeMillis() - t0) / 1000.0
    )
    return GeometryTools.geometryToROI(u, plane)
}
// --- STEP C0: percentile thresholds (per-slide, cross-slide consistent) ---
// Reads the saved density raster for densityChannel, pools all pixels strictly
// above minSampleDensity, and sets T_high / T_low to the requested percentiles
// of THAT pool. So the cut tracks each slide's own density spread.
if (thresholdMode == 'percentile') {
    try {
        def dsrv = densityMapBuilder.buildServer(imageData)
        int dw = dsrv.getWidth(), dh = dsrv.getHeight()
        // Whole-image request straight from the server: it maps the density grid
        // back to the original image's coordinate space for us. Building the
        // request by hand from the grid's own w/h gives "Invalid coordinates",
        // because RegionRequest x/y/w/h are in the ORIGINAL image's pixel space.
        def req = RegionRequest.createInstance(dsrv)
        def img = dsrv.readRegion(req)
        def raster = img.getRaster()
        int nch = raster.getNumBands()
        if (densityChannel >= nch) { println "densityChannel " + densityChannel + " >= " + nch + " channels - check densityChannel."; return }
        int rw = raster.getWidth(), rh = raster.getHeight()
        float[] samp = raster.getSamples(0, 0, rw, rh, densityChannel, (float[]) null)
        def pool = new ArrayList<Double>(samp.length)
        for (float f : samp) { if (Double.isFinite(f) && f > minSampleDensity) pool.add((double) f) }
        if (pool.size() < 100) {
            println "Percentile pool too small (" + pool.size() + " px > " + minSampleDensity + ") - falling back to absolute T_high/T_low."
        } else {
            double[] arr = pool as double[]
            Arrays.sort(arr)
            def pctile = { double p -> arr[(int) Math.min(arr.length - 1, Math.round(p / 100.0 * (arr.length - 1)))] }
            T_high = pctile(pctHigh)
            T_low  = pctile(pctLow)
            println String.format("Percentile thresholds over %d density px (> %.2f): T_high(p%.0f)=%.3f  T_low(p%.0f)=%.3f",
                arr.length, minSampleDensity, pctHigh, T_high, pctLow, T_low)
            println String.format("  density pool range: min=%.3f  median=%.3f  max=%.3f", arr[0], pctile(50), arr[arr.length-1])
            // VERIFY the cut follows the percentile: count pool pixels above each
            // threshold and confirm the fraction matches (100-pct)%. If these two
            // lines disagree with the target, the percentile did NOT take.
            int nHigh = 0, nLow = 0
            for (double d : arr) { if (d >= T_high) nHigh++; if (d >= T_low) nLow++ }
            println String.format("  CHECK core:   %.1f%% of density pool >= T_high  (target %.1f%%)",
                100.0 * nHigh / arr.length, 100.0 - pctHigh)
            println String.format("  CHECK envelope: %.1f%% of density pool >= T_low   (target %.1f%%)",
                100.0 * nLow / arr.length, 100.0 - pctLow)
        }
    } catch (Exception e) {
        println "Percentile sampling failed (" + e.getMessage() + ") - falling back to absolute T_high/T_low."
    }
}
// Safety: T_low must be strictly below T_high or the ring algebra is degenerate.
if (T_low >= T_high) { println String.format("T_low (%.3f) >= T_high (%.3f) - widen pctLow/pctHigh or the absolute values.", T_low, T_high); return }

def coreRoi     = envelopeRoiAt(T_high, 'core')
def envelopeRoi = envelopeRoiAt(T_low,  'envelope')
if (envelopeRoi == null) { println "Envelope empty - lower T_low."; return }
if (coreRoi == null)     { println "Core empty - lower T_high."; return }

// Clip contours to the tissue footprint so nothing spills onto glass
coreRoi     = RoiTools.combineROIs(coreRoi,     tissueRoi, CombineOp.INTERSECT)
envelopeRoi = RoiTools.combineROIs(envelopeRoi, tissueRoi, CombineOp.INTERSECT)

// Smooth the blocky density-grid edges: morphological OPEN then CLOSE
// (buffer -r,+r rounds outward corners; +r,-r rounds inward corners).
// Runs on the region polygon only, so it is cheap regardless of cell count.
def smooth = { roi ->
    if (roi == null || smoothRadius <= 0) return roi
    double r = smoothRadius / pxW
    // Simplify BEFORE buffering and use coarse arcs (qs=2): a low T_low makes a huge
    // jagged contour, and buffering that with JTS's default 8 quadrant-segments
    // explodes the vertex count -> OutOfMemoryError. Both cuts slash peak memory.
    int qs = 2
    def g0 = org.locationtech.jts.simplify.DouglasPeuckerSimplifier.simplify(roi.getGeometry(), simplifyPx)
    def g = g0.buffer(-r, qs).buffer(2*r, qs).buffer(-r, qs)  // open+close in one chain
    return g.isEmpty() ? roi : GeometryTools.geometryToROI(g, plane)
}
coreRoi     = smooth(coreRoi)
envelopeRoi = smooth(envelopeRoi)

// --- STEP C1a: drop tiny CORE specks --------------------------------------
// Keep only core components >= minCoreArea. A speck-sized core would otherwise
// keep its whole envelope focus alive (requireCoreInFocus) while contributing
// almost nothing to Tumour, so the focus reads as nearly all Infiltrative margin.
if (minCoreArea > 0 && coreRoi != null) {
    def cg = coreRoi.getGeometry()
    def bigCores = []
    for (int i = 0; i < cg.getNumGeometries(); i++) {
        def c = cg.getGeometryN(i)
        if (c.getArea() * px2 >= minCoreArea) bigCores << c
    }
    int coreDropped = cg.getNumGeometries() - bigCores.size()
    if (bigCores.isEmpty()) { println "All " + cg.getNumGeometries() + " core piece(s) below minCoreArea - lower minCoreArea or T_high."; return }
    coreRoi = GeometryTools.geometryToROI(
        UnaryUnionOp.union(gf.createGeometryCollection(bigCores as Geometry[])), plane)
    println String.format("Cores: kept %d, dropped %d speck(s) < %.0f um^2.",
        bigCores.size(), coreDropped, minCoreArea)
}

// --- STEP C1b: drop tumour "islands" with no dense core -------------------
// A genuine tumour focus is an envelope component that encloses a core. Stray
// low-density specks have no core inside them -> discard them. This removes the
// scattered artefact islands without any manual size threshold.
// Each envelope connected component is one candidate tumour focus. A focus is
// kept only if it (a) encloses a dense core (when requireCoreInFocus) AND
// (b) is at least minFocusArea. Because the WHOLE component is dropped, the
// focus's core/margin/necrosis all disappear together - no orphan margin rings.
if ((requireCoreInFocus || minFocusArea > 0) && coreRoi != null && envelopeRoi != null) {
    def coreGeomTest = coreRoi.getGeometry()
    def envGeom = envelopeRoi.getGeometry()
    def kept = []
    int dropCoreless = 0, dropSmall = 0
    for (int i = 0; i < envGeom.getNumGeometries(); i++) {
        def comp = envGeom.getGeometryN(i)
        if (requireCoreInFocus && !comp.intersects(coreGeomTest)) { dropCoreless++; continue }
        if (minFocusArea > 0 && comp.getArea() * px2 < minFocusArea) { dropSmall++; continue }
        kept << comp
    }
    if (!kept.isEmpty()) {
        envelopeRoi = GeometryTools.geometryToROI(
            UnaryUnionOp.union(gf.createGeometryCollection(kept as Geometry[])), plane)
    } else {
        println "All " + envGeom.getNumGeometries() + " foci were dropped - lower minFocusArea or T_low."
        return
    }
    println String.format("Tumour foci: kept %d (dropped %d coreless, %d below %.0f um^2).",
        kept.size(), dropCoreless, dropSmall, minFocusArea)
    // Restrict the core to the surviving envelope, or the core of a dropped focus
    // would reappear as Tumour (mk(coreRoi,'Tumour') below adds cores directly).
    coreRoi = RoiTools.combineROIs(coreRoi, envelopeRoi, CombineOp.INTERSECT)
    if (coreRoi == null || coreRoi.isEmpty()) { println "Core empty after focus filter - lower T_high or minFocusArea."; return }
}

// --- STEP C2: fill tumour-enclosed holes; classify large ones as necrosis --
// A low-density pocket FULLY ENCLOSED by tumour is internal necrosis/lumen, not
// non-tumour. We fill enclosed holes at the ENVELOPE level (the outer tumour
// boundary) so Non-tumour can never reach into them - this holds even if the
// higher-threshold CORE rim is broken (e.g. T_high set high), which used to let
// the necrotic centre fall through to Non-tumour.
def fillHoles = { roi ->
    def geom = roi.getGeometry()
    def outer = []
    for (int i = 0; i < geom.getNumGeometries(); i++) {
        def g = geom.getGeometryN(i)
        if (g instanceof Polygon)
            outer << gf.createPolygon(((Polygon) g).getExteriorRing().getCoordinates())
    }
    if (outer.isEmpty()) return [geom, null]
    def filled = UnaryUnionOp.union(gf.createGeometryCollection(outer as Geometry[]))
    return [filled, filled.difference(geom)]   // [filled geometry, enclosed holes]
}

def necrosisRoi = null
if (fillNecrosis) {
    def envRes  = fillHoles(envelopeRoi)
    envelopeRoi = GeometryTools.geometryToROI(envRes[0], plane)   // solid tumour territory
    // Extract the necrotic hole from the CORE (T_high) ring. In a ring /
    // peri-necrotic tumour the core is an annulus and its inner hole reaches out
    // to the dense-rim boundary (the purple T_high line), so necrosis defined
    // this way extends to the rim and absorbs the peri-necrotic shell.
    // BUT the core ring is often broken (an open C, not a closed O): then it has
    // no enclosed hole, the necrotic centre is missed AND the centre falls into
    // the margin. So first morphologically CLOSE the core - buffer(+d).buffer(-d),
    // d = necrosisBridge/2 - which bridges rim gaps up to ~necrosisBridge wide,
    // turning the C into an O. Then fillHoles on the CLOSED core gives both the
    // solid core (centre absorbed into Tumour) and the enclosed necrotic hole.
    def coreGeomRaw = coreRoi.getGeometry()          // T_high contour, BEFORE closing
    if (necrosisBridge > 0) {
        double dPx = (necrosisBridge / pxW) / 2.0
        int qs = 2   // coarse arcs: a large-radius close with default arcs can exhaust the heap
        coreGeomRaw = coreGeomRaw.buffer(dPx, qs).buffer(-dPx, qs)   // close: bridge gaps <= necrosisBridge
    }
    def coreRes = fillHoles(GeometryTools.geometryToROI(coreGeomRaw, plane))
    def coreHoles = coreRes[1]
    coreRoi     = GeometryTools.geometryToROI(coreRes[0], plane)  // solid dense core (centre filled)
    if (classifyNecrosis && coreHoles != null && !coreHoles.isEmpty()) {
        def bigHoles = []
        for (int i = 0; i < coreHoles.getNumGeometries(); i++) {
            def h = coreHoles.getGeometryN(i)
            if (h.getArea() * px2 >= minNecrosisArea) bigHoles << h
        }
        if (!bigHoles.isEmpty())
            necrosisRoi = GeometryTools.geometryToROI(
                UnaryUnionOp.union(gf.createGeometryCollection(bigHoles as Geometry[])), plane)
        double biggestUm2 = 0
        for (int i = 0; i < coreHoles.getNumGeometries(); i++)
            biggestUm2 = Math.max(biggestUm2, coreHoles.getGeometryN(i).getArea() * px2)
        println String.format("Necrosis: %d hole(s) >= %.0f um^2 kept; largest core hole = %.0f um^2.",
            bigHoles.size(), minNecrosisArea, biggestUm2)
        if (bigHoles.isEmpty())
            println "  -> none labelled. If you expect necrosis: raise necrosisBridge (broken core ring) or lower minNecrosisArea."
    } else if (classifyNecrosis) {
        println "Necrosis: no enclosed hole in the core. If you expect one, the core ring is likely open - raise necrosisBridge (currently " + necrosisBridge + " um)."
    }
}

// --- STEP D: region algebra -----------------------------------------------
// After hole-filling, envelopeRoi is the SOLID tumour footprint (core is a
// subset of it). Non-tumour is tissue minus that footprint, so it can never
// reach a tumour-enclosed pocket. tumourFootprint keeps this correct even if
// fillNecrosis is off (then envelope keeps its holes and core is added back).
def tumourFootprint = RoiTools.combineROIs(envelopeRoi, coreRoi, CombineOp.ADD)     // envelope + core
def marginRoi    = RoiTools.combineROIs(envelopeRoi, coreRoi,        CombineOp.SUBTRACT) // envelope - core
def nonTumourRoi = RoiTools.combineROIs(tissueRoi,   tumourFootprint, CombineOp.SUBTRACT) // tissue - footprint

// Keep only polygonal parts of a ROI. JTS boolean ops (SUBTRACT/ADD on the
// jagged multi-focus contours) can leave a GeometryCollection with stray
// LineStrings/Points where two boundaries touch at a line or a point. QuPath's
// resolveHierarchy() -> MaximumInscribedCircle then throws "Input geometry must
// be a Polygon or MultiPolygon". Rebuilding from the Polygon parts only removes
// those slivers without changing the area.
def polygonsOnly = { roi ->
    def g = roi.getGeometry()
    if (g == null || g.isEmpty()) return roi
    def polys = []
    for (int i = 0; i < g.getNumGeometries(); i++) {
        def gi = g.getGeometryN(i)
        if (gi instanceof Polygon && !gi.isEmpty()) polys << gi
    }
    if (polys.isEmpty()) return null
    def clean = (polys.size() == 1) ? polys[0]
                : UnaryUnionOp.union(gf.createGeometryCollection(polys as Geometry[]))
    return GeometryTools.geometryToROI(clean, plane)
}
def mk = { roi, clsName ->
    roi = (roi == null) ? null : polygonsOnly(roi)
    if (roi == null || roi.isEmpty()) return
    addObject(PathObjects.createAnnotationObject(roi, getPathClass(clsName)))
}
mk(coreRoi,      'Tumour')
mk(marginRoi,    'Infiltrative margin')
mk(nonTumourRoi, 'Non-tumour')
// Necrosis sits on top of (inside) the Tumour core; add last so resolveHierarchy
// nests it as a child region. Its cells stay counted under Tumour unless you
// prefer otherwise - see the note in the README.
if (necrosisRoi != null) mk(necrosisRoi, 'Necrosis')

// Drop tiny fragments
getAnnotationObjects().findAll { it.getPathClass() in regionClasses }.each {
    def r = it.getROI()
    if (r != null && r.getArea() * px2 < minRegionArea) removeObject(it, true)
}
resolveHierarchy()
fireHierarchyUpdate()
println "Regions: " + getAnnotationObjects().findAll{ it.getPathClass()!=null }.collect{ it.getPathClass().toString() }
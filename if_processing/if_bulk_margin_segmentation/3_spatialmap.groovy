// ============================================================================
// Step 3 (Groovy): build and SAVE the tumour-cell density map
// ----------------------------------------------------------------------------
// Programmatic replacement for Analyze > Density maps > Create density map, so
// Step 3 is scriptable / batchable instead of a GUI click. It builds a
// GAUSSIAN-WEIGHTED density of the POSITIVE class, bakes in the Jet
// colormap for display, and saves it into the PROJECT under a name DERIVED
// FROM THE IMAGE NAME, so Script 4 can load it by name
// (loadDensityMap / createAnnotationsFromDensityMap).
//
// Run AFTER classification (Script 2), which classifies detections as
// Tumour / Stroma. Set mainClassName to whichever class marks positive cells;
// 'Tumour' matches Script 2's output.
// The image MUST be in a QuPath PROJECT (density maps are saved as project
// resources; loadDensityMap can only find them there).
// ============================================================================
import qupath.lib.analysis.heatmaps.DensityMaps
import qupath.lib.analysis.heatmaps.DensityMaps.DensityMapBuilder
import qupath.lib.analysis.heatmaps.ColorModels
import qupath.lib.objects.PathObjectFilter
import qupath.lib.regions.RegionRequest
import static qupath.lib.objects.PathObjectPredicates.*

// ---- EDIT: positive class + density radius + colormap --------------------
def mainClassName    = 'Tumour'         // positive class from Script 2
double radiusMicrons = 100.0            // Gaussian radius in um (sigma of the weighting).
def colorMapName     = 'Jet'            // display palette baked into the saved map.
                                        // (Run once; if it fails, the catch below
                                        //  prints the exact names your build offers.)
// --------------------------------------------------------------------------

// Name derived from the image: densitymap_<image>_gw<radius>  (gw = Gaussian-weighted).
// Use the SAME expression for densityMapName in Script 4 so they always match,
// then the whole pipeline batches across a project with no hand-editing.
def imgName        = getCurrentImageNameWithoutExtension()
def densityMapName = 'densitymap_' + imgName + '_gw' + ((int) radiusMicrons)

// Objects that count toward the density: DETECTIONS of the positive class.
// IMPORTANT: the predicate must be JSON-SERIALIZABLE, because the density map is
// saved as a project resource (Gson). A closure coerced with 'as
// PathObjectPredicate' becomes a java.lang.reflect.Proxy, which Gson cannot
// serialize on modern Java (-> "does not opens java.lang.reflect" at save time).
// The PathObjectPredicates STATIC HELPERS build a serializable predicate:
//   filter(PathObjectFilter.DETECTIONS_ALL) -> all detections (serializable)
//   exactClassification(cls)                -> object carries this exact class
def cls = getPathClass(mainClassName)
def mainFilter = filter(PathObjectFilter.DETECTIONS_ALL).and(exactClassification(cls))

int nPos = getDetectionObjects().count { it.getPathClass() == cls }
if (nPos == 0) {
    println "No detections classified '" + mainClassName + "'. Run classification first, " +
            "or set mainClassName to your positive class."
    return
}

// GAUSSIAN-weighted density (DensityMapType.GAUSSIAN): equivalent to applying a
// Gaussian filter to the positive-cell count per pixel, area-normalised. This is
// the smooth surface Script 4 percentiles.
def builder = DensityMaps.builder(mainFilter)
                .type(DensityMaps.DensityMapType.GAUSSIAN)
                .radius(radiusMicrons)

def imageData = getCurrentImageData()

// --- Determine a sensible display max, then bake in the Jet colormap.
// The colormap is DISPLAY-ONLY: it changes how the map renders, NOT the values
// Script 4 thresholds. We read the (low-res) density raster once to get the peak
// so Jet spans 0..peak; the Jet rendering bug needs maxDisplay > minDisplay.
double dispMax = 1.0
try {
    def server = builder.buildServer(imageData)
    def req    = RegionRequest.createInstance(server)      // whole-image convenience overload
    def img    = server.readRegion(req)
    def raster = img.getRaster()
    int rw = raster.getWidth(), rh = raster.getHeight()
    float[] samp = raster.getSamples(0, 0, rw, rh, 0, (float[]) null)
    double mx = 0.0
    for (float f : samp) { if (Double.isFinite(f) && f > mx) mx = f }
    if (mx > 0) dispMax = mx
    println String.format("Density peak (channel 0) = %.4f -> Jet spans 0..%.4f", mx, dispMax)
} catch (Exception e) {
    println "  (could not read raster to scale colormap: " + e.getMessage() + " - using 0..1)"
}

// Attach the Jet colormap to the builder so it is saved with the map.
try {
    def band = ColorModels.createBand(colorMapName, 0, 0.0, dispMax)   // (name, band=0, min, max)
    def cmBuilder = ColorModels.createColorModelBuilder(band, null)    // main channel, no alpha
    builder.colorModel(cmBuilder)
    println "Colormap '" + colorMapName + "' attached to the saved density map."
} catch (Exception e) {
    println "  (could not attach colormap '" + colorMapName + "': " + e.getMessage() +
            " - map saved without a baked palette; set it in the GUI viewer instead.)"
    try {
        println "  available colormap names in this build: " +
                qupath.lib.color.ColorMaps.getColorMaps().keySet()
    } catch (Exception e2) { /* ignore */ }
}

// Save into the project so Script 4 can load it by name.
def project = getProject()
if (project == null) {
    println "No project open. Add this image to a QuPath PROJECT so the density map " +
            "can be saved as a resource and loaded by name in Script 4."
    return
}
project.getResources(DensityMaps.PROJECT_LOCATION, DensityMapBuilder.class, "json")
       .put(densityMapName, builder)
println String.format("Saved density map '%s'  (%d '%s' detections, Gaussian radius %d um).",
        densityMapName, nPos, mainClassName, (int) radiusMicrons)

// Confirm: reload and print channels, so you can copy the exact name into Script 4
// and verify densityChannel (should be 0 for this single-class map).
try {
    def b2 = loadDensityMap(densityMapName)
    def dserver = b2.buildServer(getCurrentImageData())
    println "  channels:"
    dserver.getMetadata().getChannels().eachWithIndex { ch, i -> println "    [" + i + "] " + ch.getName() }
    println "  -> in Script 4 set: densityMapName = '" + densityMapName + "'  and  densityChannel = 0"
} catch (Exception e) {
    println "  (saved, but could not reload for channel check: " + e.getMessage() + ")"
}

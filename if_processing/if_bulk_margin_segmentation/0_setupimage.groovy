// ---- Script 0: setup (run once per image) --------------------------------
import qupath.lib.images.servers.ImageServerMetadata

setImageType('FLUORESCENCE')

// EDIT: name your channels to match your acquisition.
// Order must match the channel order in the image (channel 0 first).
setChannelNames(
    'tdtomato',
    'cd31',
	'dapi',
	'iba1'     
)

// --- Verify pixel calibration (microns) -----------------------------------
def cal = getCurrentServer().getPixelCalibration()
println "Pixel width  (um): " + cal.getPixelWidthMicrons()
println "Pixel height (um): " + cal.getPixelHeightMicrons()
if (!cal.hasPixelSizeMicrons()) {
    println "!! No pixel size set. Set it via Image tab > double-click 'Pixel size', "
    println "   or uncomment the block below with your known value."
}
// // EDIT if calibration is missing (value in microns/pixel):
// def md = getCurrentServer().getMetadata()
// def newMd = new ImageServerMetadata.Builder(md)
//         .pixelSizeMicrons(1.1, 1.1)   // <-- your value
//         .build()
// getCurrentImageData().updateServerMetadata(newMd)
println "Setup done."
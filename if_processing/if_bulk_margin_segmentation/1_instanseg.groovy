// ---- Script 1: InstanSeg cell detection -----------------------------------
setImageType('FLUORESCENCE')

// InstanSeg (like every QuPath detection plugin) runs INSIDE selected
// annotations, so there must be a parent region to detect in; without one,
// detectObjects() has nothing to run on. Create a whole-image annotation when
// the image has none, then select every annotation.
if (getAnnotationObjects().isEmpty())
    createFullImageAnnotation(true)
selectAnnotations()

// EDIT: absolute path to the downloaded InstanSeg model directory. In QuPath this
// is normally <QuPath user directory>/instanseg/downloaded/<model name>. Download the
// model via Extensions > InstanSeg first. Set the INSTANSEG_MODEL environment
// variable to override without editing this file.
def modelPath = System.getenv('INSTANSEG_MODEL') ?:
        '/path/to/QuPath/instanseg/downloaded/fluorescence_nuclei_and_cells-0.1.1'
if (!new File(modelPath).exists()) {
    println "InstanSeg model not found: ${modelPath}"
    println "Edit modelPath in this script, or set INSTANSEG_MODEL."
    return
}

qupath.ext.instanseg.core.InstanSeg.builder()
    .modelPath(modelPath)
    .device("cpu")
    .inputChannels([ColorTransforms.createChannelExtractor("tdtomato"), ColorTransforms.createChannelExtractor("cd31"), ColorTransforms.createChannelExtractor("dapi"), ColorTransforms.createChannelExtractor("iba1")])
    .outputChannels()
    .tileDims(512)
    .interTilePadding(32)
    .nThreads(4)
    .makeMeasurements(false)
    .randomColors(true)
    .outputType("Default")
    .build()
    .detectObjects()
selectDetections();
runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons":1.0,"region":"ROI","tileSizeMicrons":25.0,"channel1":true,"channel2":true,"channel3":true,"channel4":true,"doMean":true,"doStdDev":true,"doMinMax":true,"doMedian":true,"doHaralick":false,"haralickMin":NaN,"haralickMax":NaN,"haralickDistance":1,"haralickBins":32}')
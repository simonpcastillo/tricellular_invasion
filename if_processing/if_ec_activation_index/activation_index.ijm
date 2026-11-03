// ============================================================================
// Endothelial activation index (ImageJ / Fiji macro)
// ----------------------------------------------------------------------------
// Measures VCAM1 and CCL2 signal in the tdTomato-NEGATIVE compartment of a
// multichannel IF image, normalised to the number of nuclei in that compartment.
//
// Steps:
//   1. Threshold the tdTomato channel and INVERT it, giving the region outside
//      the tdTomato+ tumour ("outside_tdTomato", ROI index 0).
//   2. Within that region, segment DAPI nuclei (threshold + watershed) and count
//      them with Analyze Particles.
//   3. Measure mean VCAM1 and mean CCL2 intensity over the same region and
//      divide each by the nucleus count.
//
// Output: a Results table row with NucleiCount, VCAM_Mean, CCL2_Mean and the
// two per-nucleus values. Run on one open image at a time.
//
// EDIT the channel indices and thresholding settings below to match your
// acquisition; they are 1-based and must follow the channel order in the image.
// ============================================================================

// CONFIG
origID = getImageID();
origTitle = getTitle();
roiManager("Reset");
run("Clear Results");
vcamCh = 1; tdtCh = 2; ccl2Ch = 3; dapiCh = 4;
tdtMethod = "Triangle"; dapiMethod = "Otsu";
tdtSigma = 1.0; dapiSigma = 1.0;
minNucArea = 30; maxNucArea = 1e9;

// CREATE OUTSIDE tdTom ROI
selectImage(origID);
Stack.setChannel(tdtCh);
run("Duplicate...", "title=tdt");
run("Gaussian Blur...", "sigma=" + tdtSigma);
setAutoThreshold(tdtMethod + " dark");
run("Convert to Mask");
run("Invert");
run("Create Selection");
roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "outside_tdTomato");
close();

// COUNT DAPI NUCLEI
selectImage(origID);
Stack.setChannel(dapiCh);
run("Duplicate...", "title=DAPI");
roiManager("Select", 0);
run("Clear Outside");
run("Gaussian Blur...", "sigma=" + dapiSigma);
setAutoThreshold(dapiMethod + " dark");
run("Convert to Mask");
run("Watershed");
run("Set Measurements...", "area mean standard modal min display redirect=None decimal=3");
run("Analyze Particles...", "size=" + minNucArea + "-" + maxNucArea + " show=Nothing clear");
nNuclei = nResults;
print("Nuclei outside tdTomato = " + nNuclei);
close();

// MEASURE VCAM & CCL2
run("Set Measurements...", "mean redirect=None decimal=3");
selectImage(origID);
roiManager("Select", 0);
Stack.setChannel(vcamCh);
getStatistics(area, meanV);
Stack.setChannel(ccl2Ch);
getStatistics(area, meanC);

// NORMALIZE & OUTPUT
vcamPerNuc = meanV / nNuclei;
ccl2PerNuc = meanC / nNuclei;
setResult("Region", 0, "outside_tdTomato");
setResult("NucleiCount", 0, nNuclei);
setResult("VCAM_Mean", 0, meanV);
setResult("VCAM_Mean_perNucleus", 0, vcamPerNuc);
setResult("CCL2_Mean", 0, meanC);
setResult("CCL2_Mean_perNucleus", 0, ccl2PerNuc);
updateResults();
print("VCAM mean = " + meanV + ", CCL2 mean = " + meanC);
selectImage(origID);
run("Select None");
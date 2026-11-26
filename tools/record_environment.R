#!/usr/bin/env Rscript
# Regenerate sessionInfo.txt: the R version, platform, and the versions of every
# package this repository uses. Run from the repository root:
#   Rscript tools/record_environment.R
# Packages that are not installed are listed as such rather than omitted, so the
# record shows which modules a given machine can actually run.

pk <- c("Seurat","SeuratObject","Matrix","data.table","jsonlite","png","ggplot2","patchwork",
        "mclust","optparse","homologene","ggrepel","tibble","dplyr","stringr","RANN",
        "spatstat.geom","scales","scattermore","svglite","ggh4x","cluster","tidyr","readr",
        "RcppML","liana","future","purrr","limma","edgeR","ggnewscale","fgsea",
        "AnnotationDbi","GO.db","org.Mm.eg.db")

out <- "sessionInfo.txt"
con <- file(out, "wt")
writeLines(c(
  "Environment record for github.com/simonpcastillo/tricellular_invasion",
  "Regenerate with: Rscript tools/record_environment.R",
  paste("Recorded:", format(Sys.Date())),
  "", strrep("-", 72), "R", strrep("-", 72), ""), con)
writeLines(c(R.version.string,
             paste("Platform:", R.version$platform),
             paste("Running under:", sessionInfo()$running), ""), con)
writeLines(c("Package versions used by this repository:", ""), con)
for (p in sort(pk)) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
  writeLines(sprintf("  %-18s %s", p, ifelse(is.na(v), "not installed", v)), con)
}
writeLines(c("", strrep("-", 72), "Other software", strrep("-", 72), "",
  "  QuPath             0.7   (InstanSeg extension; model fluorescence_nuclei_and_cells-0.1.1)",
  "  Fiji/ImageJ        TODO  record the version used",
  "",
  "Python environment the code was verified against (all four entry points import and",
  "show --help). This is a known-good set, not necessarily the set used for the paper:",
  "",
  "  Python             3.11.18",
  "  numpy              2.4.6",
  "  pandas             2.3.3      (scanpy constrains pandas to 2.x)",
  "  scipy              1.17.1",
  "  scikit-learn       1.9.1",
  "  matplotlib         3.11.2",
  "  scikit-image       0.26.0",
  "  shapely            2.2.0",
  "  tifffile           2026.3.3",
  "  imagecodecs        2026.3.6",
  "  scanpy             1.11.5",
  "  anndata            0.12.19",
  "  bin2cell           0.3.5"), con)
close(con)
cat("Wrote", out, "\n")

#!/usr/bin/env Rscript
# Stage 1 -- export the small tables the rescue model needs from a Seurat object,
# so stage 2 never loads the expression matrix. Sample-agnostic: everything comes
# from the command line.
#
# Usage:
#   Rscript 2a2_export_from_seurat.R <seurat.rds> <outdir> <reporter_gene_or_AUTO> \
#                                   [signature_genes.txt] [region_col]
#
# CRITICAL -- the SCT assay must have been built with depth regressed out per gene:
#   SCTransform(seu, assay = "Spatial",
#               vars.to.regress = c("nCount_Spatial", "nFeature_Spatial",
#                                   "percent.mt", "bin_count"),
#               return.only.var.genes = FALSE)
# Do NOT substitute NormalizeData()/LogNormalize (CP10K). Under bin2cell
# segmentation library size is partly cell size, which is biology; dividing it out
# removes the tumour signal along with the technical component.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: 2a2_export_from_seurat.R <seurat.rds> <outdir> <reporter_gene_or_AUTO> ",
       "[signature_genes.txt] [region_col]")
}
rds_path <- args[1]
outdir   <- args[2]
reporter <- args[3]
sig_file <- if (length(args) >= 4 && nzchar(args[4]) && args[4] != "NA") args[4] else NA
region_col <- if (length(args) >= 5 && nzchar(args[5]) && args[5] != "NA") args[5] else NA
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("loading", rds_path, "\n"); flush.console()
seu <- readRDS(rds_path)
cat("cells:", ncol(seu), " assays:", paste(names(seu@assays), collapse = ","), "\n")
stopifnot("SCT" %in% names(seu@assays))

# ---- reporter counts, from the RAW counts assay, never the normalised one.
# Normalisation zeroes genuine low-count detections (697 cells on the pilot
# sample). Detection is "any signal > 0", not "rounded > 0".
raw_assay <- if ("Spatial" %in% names(seu@assays)) "Spatial" else DefaultAssay(seu)
cnt <- GetAssayData(seu, assay = raw_assay, layer = "counts")
if (toupper(reporter) == "AUTO") {
  matches <- grep("tdtomato", rownames(cnt), ignore.case = TRUE, value = TRUE)
  if (length(matches) != 1L) {
    stop("AUTO reporter detection expected exactly one TdTomato feature, found ",
         length(matches), ": ", paste(matches, collapse = ", "),
         ". Pass the exact feature name instead.")
  }
  reporter <- matches[[1]]
  cat("AUTO-detected reporter: ", reporter, "\n", sep = "")
}
if (!reporter %in% rownames(cnt)) {
  stop(sprintf("reporter '%s' not in %s assay rownames", reporter, raw_assay))
}
rep_counts <- as.numeric(cnt[reporter, ])
cat(sprintf("reporter %s: %d cells with any signal (%.2f%%)\n",
            reporter, sum(rep_counts > 0), 100 * mean(rep_counts > 0)))
if (sum(rep_counts > 0) < 200) {
  cat("WARNING: fewer than 200 reporter-positive cells. The PU model needs a\n",
      "  few hundred positives to be worth fitting; check the reporter name.\n")
}

# ---- PCA on SCT with the reporter EXCLUDED.
# If the reporter stays in, the model relearns the probe gate and every
# validation number becomes meaningless.
DefaultAssay(seu) <- "SCT"
vf <- VariableFeatures(seu)
if (length(vf) == 0) {
  seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)
  vf  <- VariableFeatures(seu)
}
vf <- setdiff(vf, reporter)
cat("variable features used:", length(vf), "(reporter excluded)\n")
stopifnot(!reporter %in% vf)

seu <- ScaleData(seu, features = vf, verbose = FALSE)
seu <- RunPCA(seu, features = vf, npcs = 50, verbose = FALSE)
emb <- Embeddings(seu, "pca")
cat("embedding:", nrow(emb), "x", ncol(emb), "\n")

# ---- confounder signature score (optional but strongly recommended).
# This is what catches a contaminating non-tumour population being rescued. On
# sample 7003 a 214-cell neuronal cluster was rescued at 64x uplift; the score
# below is how stage 2 detects that automatically instead of by eye.
# signature_genes.txt: one gene symbol per line, '#' comments allowed.
sig_score <- rep(NA_real_, ncol(seu))
n_sig <- 0L
if (!is.na(sig_file)) {
  if (!file.exists(sig_file)) stop("signature file not found: ", sig_file)
  raw_lines <- readLines(sig_file, warn = FALSE)
  g <- trimws(raw_lines); g <- g[nzchar(g) & !startsWith(g, "#")]
  sct <- GetAssayData(seu, assay = "SCT", layer = "data")
  g <- intersect(g, rownames(sct))
  g <- setdiff(g, reporter)          # never score on the reporter itself
  n_sig <- length(g)
  cat("confounder signature:", n_sig, "genes matched in SCT\n")
  if (n_sig >= 5) {
    sig_score <- Matrix::colMeans(sct[g, , drop = FALSE])
  } else {
    cat("WARNING: <5 signature genes matched; skipping the score. Check symbol\n",
        "  case and species (Snap25 vs SNAP25).\n")
    n_sig <- 0L
  }
}

# ---- metadata: depth, clusters, coordinates, optional region membership
meta_cols <- intersect(c("nCount_Spatial", "nFeature_Spatial", "nCount_RNA",
                         "nFeature_RNA", "nCount_SCT", "nFeature_SCT",
                         "percent.mt", "percent_mt", "bin_count",
                         "seurat_clusters"),
                       colnames(seu@meta.data))
meta <- data.frame(cell = colnames(seu), seu@meta.data[, meta_cols, drop = FALSE],
                   check.names = FALSE)

# spatial coordinates -- QC and audit only, NEVER a model feature
if ("spatial" %in% names(seu@reductions)) {
  xy <- Embeddings(seu, "spatial")
  meta$x <- xy[, 1]; meta$y <- xy[, 2]
} else if (all(c("x", "y") %in% colnames(seu@meta.data))) {
  meta$x <- seu$x; meta$y <- seu$y
} else {
  cat("NOTE: no spatial coordinates found; spatial QC will be skipped.\n")
}

# optional histology region membership (e.g. in_tumour_region from a QuPath
# annotation). Used by the audit as an independent check, never as a feature.
if (!is.na(region_col)) {
  if (!region_col %in% colnames(seu@meta.data)) {
    stop("region_col '", region_col, "' not in meta.data")
  }
  meta$in_region <- as.logical(seu@meta.data[[region_col]])
  cat(sprintf("region column %s: %.1f%% of cells inside\n",
              region_col, 100 * mean(meta$in_region, na.rm = TRUE)))
}

if (n_sig > 0) meta$confounder_score <- as.numeric(sig_score)

write.csv(data.frame(cell = colnames(seu), emb, check.names = FALSE),
          gzfile(file.path(outdir, "embedding.csv.gz")), row.names = FALSE)
write.csv(meta, gzfile(file.path(outdir, "meta.csv.gz")), row.names = FALSE)
write.csv(data.frame(cell = colnames(seu), reporter_count = rep_counts),
          gzfile(file.path(outdir, "reporter.csv.gz")), row.names = FALSE)

# provenance, so a run can be traced back to its inputs
writeLines(c(
  paste0("rds=", normalizePath(rds_path)),
  paste0("reporter=", reporter),
  paste0("raw_assay=", raw_assay),
  paste0("n_cells=", ncol(seu)),
  paste0("n_reporter_pos=", sum(rep_counts > 0)),
  paste0("n_variable_features=", length(vf)),
  paste0("n_signature_genes=", n_sig),
  paste0("region_col=", ifelse(is.na(region_col), "none", region_col)),
  paste0("seurat_version=", as.character(packageVersion("Seurat"))),
  paste0("exported=", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
), file.path(outdir, "export_provenance.txt"))

cat("wrote embedding.csv.gz, meta.csv.gz, reporter.csv.gz to", outdir, "\n")

#!/usr/bin/env Rscript
# =============================================================================
# Step 1 | QC plots + cell filtering of Step 0 Seurat RDS files
# =============================================================================
# Run after Step 0: Rscript 1_qc.R
#
# Default paths (relative to this script):
#   output/00_seurat/<sample_id>.rds         (input from Step 0)
#   output/01_qc/<sample_id>_filtered.rds    (output)
#   output/01_qc/<sample_id>_QC.pdf          (output)
#   output/01_qc/QC_summary.csv              (output)
#
# Override the shared output root with PIPELINE_OUTPUT_DIR.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# ---- Paths and QC thresholds -------------------------------------------------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
  normalizePath(getwd())
}

OUTPUT_DIR <- Sys.getenv("PIPELINE_OUTPUT_DIR", unset = file.path(SCRIPT_DIR, "output"))
SEURAT_DIR <- file.path(OUTPUT_DIR, "00_seurat")
QC_DIR     <- file.path(OUTPUT_DIR, "01_qc")

QC_MIN_FEATURES <- 50L           # retain cells with >50 genes
QC_MIN_COUNTS   <- 50L           # retain cells with >50 UMIs
QC_MAX_MT_PCT   <- 20            # retain cells with <20% mitochondrial UMIs

QC_METRICS <- c("nFeature_Spatial", "nCount_Spatial", "percent.mt")
QC_SUMMARY <- file.path(QC_DIR, "QC_summary.csv")

if (!dir.exists(SEURAT_DIR)) {
  stop("Step 0 output folder does not exist: ", SEURAT_DIR,
       "\nRun Rscript 0_b2c_seurat.R first.")
}
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ----------------------------------------------------------------
make_qc_violins <- function(object, title, show_thresholds = FALSE) {
  cutoffs <- c(nFeature_Spatial = QC_MIN_FEATURES,
               nCount_Spatial = QC_MIN_COUNTS,
               percent.mt = QC_MAX_MT_PCT)

  plots <- lapply(QC_METRICS, function(metric) {
    p <- Seurat::VlnPlot(
      object, features = metric, pt.size = 0,
      group.by = "sample_id",
      log = metric %in% c("nFeature_Spatial", "nCount_Spatial")
    ) +
      Seurat::NoLegend() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        plot.title = ggplot2::element_text(size = 9)
      )
    if (show_thresholds) {
      p <- p + ggplot2::geom_hline(
        yintercept = cutoffs[[metric]], linetype = "dashed",
        colour = "firebrick", linewidth = 0.5
      )
    }
    p
  })

  patchwork::wrap_plots(plots, nrow = 1) +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"))
    )
}

recalculate_mt_fraction <- function(object) {
  # Recompute consistently even for RDS files created by an older script.
  mito_genes <- grep("^mt-", rownames(object), ignore.case = TRUE, value = TRUE)
  if (!length(mito_genes)) {
    stop("No mitochondrial genes found (^MT-/^mt-). ",
         "Check features.tsv gene symbols and rerun Step 0.")
  }
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, features = mito_genes)
  object
}

upsert_qc_summary <- function(new_row) {
  if (file.exists(QC_SUMMARY)) {
    previous <- utils::read.csv(QC_SUMMARY, stringsAsFactors = FALSE)
    if (!"sample_id" %in% names(previous)) {
      stop("Incompatible older QC_summary.csv; move it aside and rerun Step 1")
    }
    previous <- previous[previous$sample_id != new_row$sample_id, , drop = FALSE]
    rows <- rbind(previous, new_row)
  } else {
    rows <- new_row
  }
  rows <- rows[order(rows$sample_id), , drop = FALSE]
  utils::write.csv(rows, QC_SUMMARY, row.names = FALSE)
}

# ---- Main -------------------------------------------------------------------
rds_files <- sort(list.files(SEURAT_DIR, pattern = "\\.rds$",
                             full.names = TRUE, ignore.case = TRUE))
if (!length(rds_files)) {
  stop("No RDS files from Step 0 were found in: ", SEURAT_DIR)
}

for (in_path in rds_files) {
  sample_id <- tools::file_path_sans_ext(basename(in_path))
  out_rds <- file.path(QC_DIR, paste0(sample_id, "_filtered.rds"))
  out_pdf <- file.path(QC_DIR, paste0(sample_id, "_QC.pdf"))

  summary_has_sample <- FALSE
  if (file.exists(QC_SUMMARY)) {
    previous_ids <- utils::read.csv(QC_SUMMARY, stringsAsFactors = FALSE)$sample_id
    summary_has_sample <- sample_id %in% previous_ids
  }
  if (file.exists(out_rds) && file.exists(out_pdf) && summary_has_sample) {
    message("[SKIP] QC outputs and summary exist: ", sample_id)
    next
  }

  message("\n[Step 1] QC: ", sample_id)
  obj <- readRDS(in_path)
  DefaultAssay(obj) <- "Spatial"
  obj$sample_id <- sample_id
  obj <- recalculate_mt_fraction(obj)

  missing_metrics <- setdiff(QC_METRICS, colnames(obj@meta.data))
  if (length(missing_metrics)) {
    stop("Missing QC metrics for ", sample_id, ": ", paste(missing_metrics, collapse = ", "))
  }

  n_before <- ncol(obj)
  if (n_before == 0L) stop("No cells in Seurat object: ", sample_id)
  keep <- with(obj@meta.data,
               nFeature_Spatial > QC_MIN_FEATURES &
               nCount_Spatial > QC_MIN_COUNTS &
               percent.mt < QC_MAX_MT_PCT)
  keep[is.na(keep)] <- FALSE
  kept_cells <- rownames(obj@meta.data)[keep]
  n_after <- length(kept_cells)
  if (n_after == 0L) {
    stop(sample_id, ": no cells pass QC. Review the thresholds and input data.")
  }

  filtered <- subset(obj, cells = kept_cells)
  pct_removed <- round(100 * (n_before - n_after) / n_before, 1)
  message("  Cells: ", n_before, " -> ", n_after, " (", pct_removed, "% removed)")

  before_plot <- make_qc_violins(
    obj, paste0(sample_id, " | before QC (n=", n_before, ")"),
    show_thresholds = TRUE
  )
  after_plot <- make_qc_violins(
    filtered, paste0(sample_id, " | after QC (n=", n_after,
                     "; ", pct_removed, "% removed)")
  )

  grDevices::pdf(out_pdf, width = 10, height = 5)
  print(before_plot)
  print(after_plot)
  grDevices::dev.off()
  saveRDS(filtered, out_rds)

  qc_row <- data.frame(
    sample_id = sample_id,
    cells_before = n_before,
    cells_after = n_after,
    cells_removed = n_before - n_after,
    pct_removed = pct_removed,
    median_nFeature_before = median(obj$nFeature_Spatial),
    median_nCount_before = median(obj$nCount_Spatial),
    median_mt_before = median(obj$percent.mt),
    median_nFeature_after = median(filtered$nFeature_Spatial),
    median_nCount_after = median(filtered$nCount_Spatial),
    median_mt_after = median(filtered$percent.mt)
  )
  upsert_qc_summary(qc_row)
  message("  Saved: ", out_rds, " | ", out_pdf)
  rm(obj, filtered, before_plot, after_plot)
  invisible(gc(verbose = FALSE))
}

message("\n[Step 1] Finished. QC summary: ", QC_SUMMARY)

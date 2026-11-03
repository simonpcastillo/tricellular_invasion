#!/usr/bin/env Rscript
# =============================================================================
# Step 0 | bin2cell counts + metadata + image -> Seurat RDS
# =============================================================================
# Run from anywhere: Rscript 0_b2c_seurat.R
#
# Default paths (relative to this script):
#   input/bin2cell/<sample_id>/
#   input/spaceranger/<sample_id>/.../spatial/
#   output/00_seurat/<sample_id>.rds
#
# Override the two root folders with the environment variables
# PIPELINE_INPUT_DIR and PIPELINE_OUTPUT_DIR, if required.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(png)
})

# ---- Paths ------------------------------------------------------------------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
  normalizePath(getwd()) # convenient when sourced from the repository root
}

INPUT_DIR  <- Sys.getenv("PIPELINE_INPUT_DIR",  unset = file.path(SCRIPT_DIR, "input"))
OUTPUT_DIR <- Sys.getenv("PIPELINE_OUTPUT_DIR", unset = file.path(SCRIPT_DIR, "output"))

BIN2CELL_DIR    <- file.path(INPUT_DIR,  "bin2cell")
SPACERANGER_DIR <- file.path(INPUT_DIR,  "spaceranger")
SEURAT_DIR      <- file.path(OUTPUT_DIR, "00_seurat")

if (!dir.exists(BIN2CELL_DIR)) stop("Missing bin2cell input folder: ", BIN2CELL_DIR)
if (!dir.exists(SPACERANGER_DIR)) stop("Missing Space Ranger input folder: ", SPACERANGER_DIR)
dir.create(SEURAT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Helper functions -------------------------------------------------------
read_counts <- function(sample_dir) {
  required <- file.path(sample_dir, c("matrix.mtx", "features.tsv", "barcodes.tsv"))
  missing <- required[!file.exists(required)]
  if (length(missing)) stop("Missing count-matrix file(s): ", paste(missing, collapse = ", "))

  counts <- as(Matrix::readMM(required[[1L]]), "CsparseMatrix")
  features <- data.table::fread(required[[2L]], header = FALSE)
  barcodes <- as.character(data.table::fread(required[[3L]], header = FALSE)[[1L]])

  # Standard 10x features.tsv: first column = stable ID, second = gene symbol.
  # Use gene symbols so mitochondrial QC works for human and mouse data.
  gene_names <- as.character(features[[if (ncol(features) >= 2L) 2L else 1L]])
  if (anyNA(gene_names) || any(!nzchar(gene_names))) {
    stop("Empty or missing gene names in: ", required[[2L]])
  }
  if (nrow(counts) != length(gene_names) || ncol(counts) != length(barcodes)) {
    stop("matrix.mtx dimensions do not match features.tsv / barcodes.tsv in ", sample_dir)
  }
  if (anyDuplicated(barcodes)) stop("Duplicate cell barcodes in: ", required[[3L]])

  rownames(counts) <- make.unique(gene_names)
  colnames(counts) <- barcodes
  counts
}

read_coordinates <- function(path) {
  xy <- data.table::fread(path)
  if (!all(c("x_full", "y_full") %in% colnames(xy))) {
    stop("Expected x_full and y_full in: ", path)
  }
  cell_ids <- as.character(xy[[1L]])
  if (anyNA(cell_ids) || any(!nzchar(cell_ids)) || anyDuplicated(cell_ids)) {
    stop("The first column must contain unique, nonempty cell IDs: ", path)
  }
  if (any(!is.finite(xy$x_full)) || any(!is.finite(xy$y_full))) {
    stop("Missing or non-finite full-resolution coordinates in: ", path)
  }

  # Seurat's VisiumV1 coordinates are stored in FULL-RESOLUTION pixels.
  # The real Space Ranger hires scale is applied when the image is plotted.
  data.frame(
    tissue   = 1L,
    row      = as.integer(round(xy$y_full)),
    col      = as.integer(round(xy$x_full)),
    imagerow = as.numeric(xy$y_full),
    imagecol = as.numeric(xy$x_full),
    row.names = cell_ids
  )
}

find_spatial_files <- function(sample_id) {
  # All spatial files must belong to this exact sample directory.
  sample_root <- file.path(SPACERANGER_DIR, sample_id)
  if (!dir.exists(sample_root)) {
    stop("Missing matching Space Ranger sample directory: ", sample_root)
  }

  scale_files <- list.files(sample_root, "^scalefactors_json\\.json$",
                            full.names = TRUE, recursive = TRUE)
  image_files <- list.files(sample_root, "^tissue_hires_image\\.png$",
                            full.names = TRUE, recursive = TRUE)
  possible_dirs <- intersect(dirname(scale_files), dirname(image_files))
  if (!length(possible_dirs)) {
    stop("Could not find scalefactors_json.json and tissue_hires_image.png ",
         "together under: ", sample_root)
  }

  # When multiple Space Ranger bin sizes are present, prefer 2-um bins.
  preferred <- possible_dirs[grepl("(^|/)square_002um(/|$)", possible_dirs)]
  if (length(preferred)) possible_dirs <- preferred
  if (length(possible_dirs) != 1L) {
    stop("Multiple matching spatial directories for ", sample_id,
         ". Keep one or select the intended bin size: ",
         paste(possible_dirs, collapse = " | "))
  }

  list(
    scales = file.path(possible_dirs, "scalefactors_json.json"),
    image  = file.path(possible_dirs, "tissue_hires_image.png")
  )
}

make_visium_image <- function(image_path, sf, coordinates) {
  hires <- as.numeric(sf$tissue_hires_scalef)
  if (length(hires) != 1L || !is.finite(hires) || hires <= 0) {
    stop("Invalid/missing tissue_hires_scalef in Space Ranger JSON")
  }

  he <- png::readPNG(image_path)
  if (length(dim(he)) == 2L) {
    he <- array(rep(he, 3L), dim = c(dim(he), 3L))
  }

  # bin2cell coordinates represent cell centres, not Visium spots.
  # The 1-pixel spot diameter is a plotting placeholder, not a cell size.
  # Only the HIRES PNG is loaded. Give BOTH plot scale choices its scale factor
  # to keep default SpatialPlot(image.scale="lowres") overlays aligned.
  # The original Space Ranger JSON is retained in obj@misc below.
  image <- methods::new(
    "VisiumV1",
    assay = "Spatial",
    key = "slice1_",
    image = he,
    coordinates = coordinates,
    scale.factors = Seurat::scalefactors(
      spot = 1, fiducial = 1, hires = hires, lowres = hires
    )
  )
  image@spot.radius <- SeuratObject::Radius(image, scale = "hires")
  image
}

add_mitochondrial_qc <- function(object) {
  mito_genes <- grep("^mt-", rownames(object), ignore.case = TRUE, value = TRUE)
  if (!length(mito_genes)) {
    stop("No mitochondrial genes (^MT-/^mt-) were found. ",
         "Check that features.tsv contains gene symbols (usually column 2).")
  }
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, features = mito_genes)
  object
}

# ---- Main -------------------------------------------------------------------
sample_dirs <- list.files(BIN2CELL_DIR, full.names = TRUE)
sample_dirs <- sort(sample_dirs[dir.exists(sample_dirs)])
if (!length(sample_dirs)) stop("No sample directories found in: ", BIN2CELL_DIR)

for (sample_dir in sample_dirs) {
  sample_id <- basename(sample_dir)
  out_path <- file.path(SEURAT_DIR, paste0(sample_id, ".rds"))
  if (file.exists(out_path)) {
    message("[SKIP] Seurat RDS exists: ", sample_id)
    next
  }

  message("\n[Step 0] Processing ", sample_id)
  spatial <- find_spatial_files(sample_id)
  sf <- jsonlite::fromJSON(spatial$scales)
  counts <- read_counts(sample_dir)

  label_path <- file.path(sample_dir, "cell_labels.csv")
  coord_path <- file.path(sample_dir, "cdata_spatial_fullres.csv")
  if (!file.exists(label_path) || !file.exists(coord_path)) {
    stop("Missing cell_labels.csv or cdata_spatial_fullres.csv for ", sample_id)
  }
  labels <- utils::read.csv(label_path, row.names = 1L, check.names = FALSE)
  coords <- read_coordinates(coord_path)

  # Preserve the order of barcodes in the count matrix.
  common_cells <- colnames(counts)[
    colnames(counts) %in% rownames(labels) & colnames(counts) %in% rownames(coords)
  ]
  if (!length(common_cells)) stop("No shared cell IDs across count/labels/coords for ", sample_id)
  if (length(common_cells) < ncol(counts)) {
    warning(sample_id, ": keeping ", length(common_cells), "/", ncol(counts),
            " cells with BOTH annotations and spatial coordinates")
  }

  obj <- Seurat::CreateSeuratObject(counts = counts[, common_cells, drop = FALSE],
                                    assay = "Spatial", project = sample_id)
  obj <- Seurat::AddMetaData(obj, labels[common_cells, , drop = FALSE])
  obj$sample_id <- sample_id
  obj <- add_mitochondrial_qc(obj)

  coords <- coords[common_cells, , drop = FALSE]
  obj[["slice1"]] <- make_visium_image(spatial$image, sf, coords)
  obj@misc$spaceranger_scalefactors <- sf

  img_dim <- dim(obj@images$slice1@image)
  x_hires <- coords$imagecol * sf$tissue_hires_scalef
  y_hires <- coords$imagerow * sf$tissue_hires_scalef
  inside <- x_hires >= 0 & x_hires <= img_dim[2L] &
            y_hires >= 0 & y_hires <= img_dim[1L]
  message("  Cells: ", ncol(obj), " | genes: ", nrow(obj),
          " | hires image: ", img_dim[2L], "x", img_dim[1L],
          " | coordinates within image: ", round(100 * mean(inside), 1), "%")
  if (mean(inside) < 0.95) {
    warning(sample_id, ": >5% of coordinates fall outside the hires image; ",
            "verify image and coordinate alignment")
  }

  saveRDS(obj, out_path)
  message("  Saved: ", out_path)
  rm(obj, counts, coords, labels)
  invisible(gc(verbose = FALSE))
}

message("\n[Step 0] Finished.")

# =============================================================================
# Export lightweight tables from classified Seurat object
# =============================================================================
# Usage:
#   Rscript export_tables.R \
#     --rds     /path/to/sample_classified.rds \
#     --out_dir /path/to/output_dir \
#     --lib_path /path/to/R/libs   # optional
# =============================================================================

suppressPackageStartupMessages(library(optparse))

opt_list <- list(
  make_option("--rds",      type = "character", help = "Path to classified RDS"),
  make_option("--out_dir",  type = "character", help = "Output directory"),
  make_option("--lib_path", type = "character", default = NULL,
              help = "Optional R library path")
)

opt <- parse_args(OptionParser(option_list = opt_list))
for (arg in c("rds", "out_dir"))
  if (is.null(opt[[arg]])) stop("Missing required argument: --", arg)

if (!is.null(opt$lib_path)) .libPaths(opt$lib_path)

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

RDS_PATH <- opt$rds
OUT_DIR  <- opt$out_dir

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

case_id <- str_extract(basename(RDS_PATH), "[0-9]{4}_Scan[0-9]+")
message("Loading: ", RDS_PATH)
seu <- readRDS(RDS_PATH)
message("  Sample: ", case_id, "  Cells: ", ncol(seu))
message("  Assays: ", paste(Assays(seu), collapse = ", "))

# detect score columns present in this object
score_cols <- grep("^Sig_", colnames(seu@meta.data), value = TRUE)
message("  Score columns found: ", length(score_cols))

# base metadata columns — take what exists, skip what doesn't
base_meta_cols <- c("cell_class", "pval_min", "seurat_clusters",
                    "nCount_Spatial", "nFeature_Spatial", "percent.mt", 'tdtom_label')
base_meta_cols <- base_meta_cols[base_meta_cols %in% colnames(seu@meta.data)]

meta_df <- seu@meta.data %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, all_of(base_meta_cols), all_of(score_cols))

# =============================================================================
# 1. umap_df
# =============================================================================

message("Building umap_df ...")
umap_df <- as.data.frame(Embeddings(seu, "umap")) %>%
  rownames_to_column("cell_id") %>%
  left_join(meta_df, by = "cell_id")

umap_path <- file.path(OUT_DIR, paste0(case_id, "_umap_df.csv"))
write.csv(umap_df, umap_path, row.names = FALSE)
message("umap_df -> ", umap_path,
        "  (", nrow(umap_df), " cells x ", ncol(umap_df), " cols)")

# =============================================================================
# 2. plot_df
# =============================================================================

message("Building plot_df ...")

# handle both old (imagerow/imagecol) and new (x/y) coordinate naming
coords <- GetTissueCoordinates(seu, image = "slice1")

plot_df <- coords %>%
  rownames_to_column("cell_id") %>%
  left_join(meta_df, by = "cell_id")

plot_path <- file.path(OUT_DIR, paste0(case_id, "_plot_df.csv"))
write.csv(plot_df, plot_path, row.names = FALSE)
message("plot_df -> ", plot_path,
        "  (", nrow(plot_df), " cells x ", ncol(plot_df), " cols)")

# =============================================================================
# 3. Sparse expression matrix
# =============================================================================

message("Saving sparse expression matrix ...")
expr_sparse <- GetAssayData(seu, assay = "SCT", layer = "data")

expr_meta <- seu@meta.data %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, all_of(base_meta_cols))

expr_out  <- list(
  matrix   = expr_sparse,  # dgCMatrix: genes x cells
  metadata = expr_meta      # data.frame: one row per cell
)

expr_path <- file.path(OUT_DIR, paste0(case_id, "_expr_sparse.rds"))
saveRDS(expr_out, expr_path)
message("expr_sparse -> ", expr_path,
        "  (", nrow(expr_sparse), " genes x ", ncol(expr_sparse), " cells)",
        "  (", round(file.size(expr_path) / 1024^2, 1), " MB)")

message("\nDone. All outputs in: ", OUT_DIR)
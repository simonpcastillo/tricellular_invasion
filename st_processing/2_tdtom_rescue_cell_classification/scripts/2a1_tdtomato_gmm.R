# =============================================================================
# 2a1_tdtomato_gmm.R
# SCTransform + clustering + UMAP + GMM gating of TdTomato (G=2)
# No signature scoring, no permutation test, no cell type classification
# =============================================================================
# Usage:
#   Rscript 2a1_tdtomato_gmm.R \
#     --rds      /path/to/sample_filtered.rds \
#     --out_dir  /path/to/output_dir \
#     --lib_path /optional/path/to/R/libs
# =============================================================================

suppressPackageStartupMessages(library(optparse))

opt_list <- list(
  make_option("--rds",      type = "character", help = "Path to input RDS file"),
  make_option("--out_dir",  type = "character", help = "Output directory"),
  make_option("--lib_path", type = "character", default = NULL),
  make_option("--seed",     type = "integer",   default = 123L),
  make_option("--n_pcs",    type = "integer",   default = 30L),
  make_option("--res",      type = "double",    default = 0.4),
  make_option("--gmm_g",    type = "integer",   default = 2L,
              help = "Number of GMM components (default: 2)")
)

opt <- parse_args(OptionParser(option_list = opt_list))

for (arg in c("rds", "out_dir"))
  if (is.null(opt[[arg]])) stop("Missing required argument: --", arg)

if (!is.null(opt$lib_path)) .libPaths(opt$lib_path)

suppressPackageStartupMessages({
  library(Seurat)
  library(mclust)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(jsonlite)
})

options(future.globals.maxSize = Inf)

RDS_PATH    <- opt$rds
OUT_DIR     <- opt$out_dir
SEED        <- opt$seed
N_PCS       <- opt$n_pcs
CLUSTER_RES <- opt$res
GMM_G       <- opt$gmm_g

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

message("=== Parameters ===")
message("  RDS:         ", RDS_PATH)
message("  OUT_DIR:     ", OUT_DIR)
message("  N_PCS:       ", N_PCS)
message("  CLUSTER_RES: ", CLUSTER_RES)
message("  GMM_G:       ", GMM_G)

# =============================================================================
# 1. Load & preprocess
# =============================================================================

message("\n=== Step 1: Load & preprocess ===")

seu     <- readRDS(RDS_PATH)
seu     <- UpdateSeuratObject(seu)
case_id <- str_extract(basename(RDS_PATH), "[0-9]{4}_Scan[0-9]+")
message("  Sample: ", case_id, "  Cells: ", ncol(seu))

message("Computing percent.mt ...")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-|^MT-")

message("Running SCTransform ...")
seu <- SCTransform(seu,
                   assay                 = "Spatial",
                   vars.to.regress       = c("nCount_Spatial", "nFeature_Spatial",
                                             "percent.mt", "bin_count"),
                   verbose               = TRUE,
                   return.only.var.genes = FALSE)

message("Running PCA ...")
seu <- RunPCA(seu, assay = "SCT", npcs = N_PCS, verbose = FALSE)

message("Finding neighbours & clusters ...")
seu <- FindNeighbors(seu, dims = 1:N_PCS, verbose = FALSE)
seu <- FindClusters(seu,  resolution = CLUSTER_RES, verbose = FALSE)

message("Running UMAP ...")
seu <- RunUMAP(seu, dims = 1:N_PCS, verbose = FALSE)

DefaultAssay(seu) <- "SCT"
message("  Clusters: ", nlevels(seu$seurat_clusters))

# =============================================================================
# 2. GMM gating of TdTomato
# =============================================================================

message("\n=== Step 2: GMM gating (G=", GMM_G, ") ===")

tdtom_gene <- grep("tdtomato", rownames(seu), value = TRUE, ignore.case = TRUE)[1]
if (is.na(tdtom_gene)) stop("TdTomato gene not found in rownames(seu)")
message("  TdTomato gene: ", tdtom_gene)

tdtom_expr        <- as.numeric(GetAssayData(seu, assay = "SCT", layer = "data")[tdtom_gene, ])
names(tdtom_expr) <- colnames(seu)

# GMM fit — NA-safe
valid_idx  <- !is.na(tdtom_expr)
gmm_input  <- tdtom_expr[valid_idx]

set.seed(SEED)
gmm_fit <- mclust::Mclust(gmm_input, G = GMM_G, verbose = FALSE)

if (is.null(gmm_fit))
  stop("GMM failed to converge — try a different G or check TdTomato expression distribution")

message("  GMM model selected: ", gmm_fit$modelName)
message("  Component means:    ", paste(round(gmm_fit$parameters$mean, 3), collapse = ", "))
message("  Component props:    ", paste(round(gmm_fit$parameters$pro,  3), collapse = ", "))

# The positive component = highest mean
pos_component <- which.max(gmm_fit$parameters$mean)
message("  TdTom_pos component: ", pos_component,
        " (mean = ", round(gmm_fit$parameters$mean[pos_component], 3), ")")

# Gate threshold = midpoint between the two component means
comp_means <- sort(gmm_fit$parameters$mean)
gmm_threshold <- mean(comp_means)
message("  Gate threshold (midpoint of means): ", round(gmm_threshold, 4))

# Assign labels using posterior probabilities: P(pos_component | x) > 0.5
classifications          <- rep(NA_character_, length(tdtom_expr))
names(classifications)   <- names(tdtom_expr)

post_pos <- gmm_fit$z[, pos_component]   # posterior P(positive component)

classifications[valid_idx] <- ifelse(post_pos > 0.5, "TdTom_pos", "TdTom_neg")

seu$tdtom_expr      <- tdtom_expr
seu$tdtom_label     <- classifications
seu$tdtom_gmm_prob  <- NA_real_
seu$tdtom_gmm_prob[valid_idx] <- post_pos

message("  TdTom_pos: ", sum(seu$tdtom_label == "TdTom_pos", na.rm = TRUE))
message("  TdTom_neg: ", sum(seu$tdtom_label == "TdTom_neg", na.rm = TRUE))
if (any(is.na(seu$tdtom_label)))
  message("  NA (missing expression): ", sum(is.na(seu$tdtom_label)))

# =============================================================================
# 3. Plots
# =============================================================================

message("\n=== Step 3: Plots ===")

# GMM density plot
gmm_df <- tibble(
  expr  = gmm_input,
  label = classifications[valid_idx]
)

# Build component density curves from GMM parameters
x_seq <- seq(min(gmm_input), max(gmm_input), length.out = 500)
dens_df <- purrr::map_dfr(seq_len(GMM_G), function(k) {
  tibble(
    x    = x_seq,
    y    = gmm_fit$parameters$pro[k] *
           dnorm(x_seq,
                 mean = gmm_fit$parameters$mean[k],
                 sd   = sqrt(gmm_fit$parameters$variance$sigmasq[
                               min(k, length(gmm_fit$parameters$variance$sigmasq))])),
    comp = paste0("Component ", k,
                  " (mu=", round(gmm_fit$parameters$mean[k], 2), ")")
  )
})

p_gmm <- ggplot(gmm_df, aes(x = expr)) +
  geom_histogram(aes(y = after_stat(density), fill = label),
                 bins = 80, alpha = 0.5) +
  geom_line(data = dens_df, aes(x = x, y = y, colour = comp),
            linewidth = 1) +
  geom_vline(xintercept = gmm_threshold, linetype = "dashed",
             colour = "firebrick", linewidth = 0.8) +
  annotate("text", x = gmm_threshold, y = Inf,
           label = sprintf("threshold = %.3f", gmm_threshold),
           hjust = -0.1, vjust = 1.5, size = 3, colour = "firebrick") +
  scale_fill_manual(values = c("TdTom_pos" = "#e6550d", "TdTom_neg" = "#d9d9d9"),
                    na.value = "grey70") +
  labs(title = paste0(case_id, " — GMM gating of ", tdtom_gene),
       x     = "SCT normalised expression",
       y     = "Density",
       fill  = "Gate",
       colour = "GMM component") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Posterior probability distribution
p_post <- ggplot(tibble(prob = post_pos, label = classifications[valid_idx]),
                 aes(x = prob, fill = label)) +
  geom_histogram(bins = 60, alpha = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "firebrick") +
  scale_fill_manual(values = c("TdTom_pos" = "#e6550d", "TdTom_neg" = "#d9d9d9")) +
  labs(title = "Posterior P(TdTom_pos component)",
       x = "Posterior probability", y = "Count", fill = "Gate") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# UMAP plots
umap_df <- as.data.frame(Embeddings(seu, "umap")) %>%
  rownames_to_column("cell_id") %>%
  left_join(
    seu@meta.data %>%
      rownames_to_column("cell_id") %>%
      dplyr::select(cell_id, tdtom_label, tdtom_expr,
                    tdtom_gmm_prob, seurat_clusters),
    by = "cell_id"
  )

cluster_centres <- umap_df %>%
  group_by(seurat_clusters) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

p_label <- ggplot(umap_df, aes(umap_1, umap_2, colour = tdtom_label)) +
  geom_point(size = 0.3, alpha = 0.7) +
  scale_colour_manual(values = c("TdTom_pos" = "#e6550d", "TdTom_neg" = "#d9d9d9"),
                      na.value = "grey70") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  labs(title = "TdTomato GMM label", colour = NULL,
       x = "UMAP 1", y = "UMAP 2") +
  theme_classic() + theme(plot.title = element_text(face = "bold"))

p_prob <- ggplot(umap_df, aes(umap_1, umap_2, colour = tdtom_gmm_prob)) +
  geom_point(size = 0.3, alpha = 0.7) +
  scale_colour_viridis_c(option = "magma", name = "P(pos)") +
  labs(title = "Posterior P(TdTom_pos)", x = "UMAP 1", y = "UMAP 2") +
  theme_classic() + theme(plot.title = element_text(face = "bold"))

p_expr <- ggplot(umap_df, aes(umap_1, umap_2, colour = tdtom_expr)) +
  geom_point(size = 0.3, alpha = 0.7) +
  scale_colour_viridis_c(option = "magma", name = "SCT expr") +
  labs(title = paste0(tdtom_gene, " expression"), x = "UMAP 1", y = "UMAP 2") +
  theme_classic() + theme(plot.title = element_text(face = "bold"))

p_cl <- ggplot(umap_df, aes(umap_1, umap_2, colour = factor(seurat_clusters))) +
  geom_point(size = 0.3, alpha = 0.5) +
  geom_label_repel(data = cluster_centres, aes(label = seurat_clusters),
                   colour = "black", fill = "white", size = 3,
                   fontface = "bold", max.overlaps = Inf) +
  guides(colour = "none") +
  labs(title = "Clusters", x = "UMAP 1", y = "UMAP 2") +
  theme_classic() + theme(plot.title = element_text(face = "bold"))

pdf_path <- file.path(OUT_DIR, paste0(case_id, "_gmm_plots.pdf"))
pdf(pdf_path, width = 14, height = 7)
print(p_gmm | p_post)
print(p_label | p_prob)
print(p_expr  | p_cl)
dev.off()
message("  Plots saved -> ", pdf_path)

# =============================================================================
# 4. Save GMM summary
# =============================================================================

gmm_summary <- tibble(
  sample          = case_id,
  gene            = tdtom_gene,
  gmm_model       = gmm_fit$modelName,
  n_components    = GMM_G,
  pos_component   = pos_component,
  mean_neg        = comp_means[1],
  mean_pos        = comp_means[2],
  threshold       = gmm_threshold,
  n_pos           = sum(seu$tdtom_label == "TdTom_pos", na.rm = TRUE),
  n_neg           = sum(seu$tdtom_label == "TdTom_neg", na.rm = TRUE),
  pct_pos         = round(mean(seu$tdtom_label == "TdTom_pos", na.rm = TRUE) * 100, 2)
)

write.csv(gmm_summary,
          file.path(OUT_DIR, paste0(case_id, "_gmm_summary.csv")),
          row.names = FALSE)
message("  GMM summary saved")

# =============================================================================
# 5. GeoJSON export — TdTomato label
# =============================================================================

message("\n=== Step 5: GeoJSON export ===")

# Ensure imagerow / imagecol are available
if (!all(c("imagerow", "imagecol") %in% colnames(seu@meta.data))) {
  message("  imagerow/imagecol not in metadata — fetching from tissue coordinates")
  tissue_coords      <- GetTissueCoordinates(seu)
  seu@meta.data$imagerow <- tissue_coords[rownames(seu@meta.data), "imagerow"]
  seu@meta.data$imagecol <- tissue_coords[rownames(seu@meta.data), "imagecol"]
}

geo_df <- seu@meta.data %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, imagerow, imagecol,
                tdtom_label, tdtom_gmm_prob, tdtom_expr)

tdtom_colours <- list(TdTom_pos = "#e6550d", TdTom_neg = "#d9d9d9")

features <- lapply(seq_len(nrow(geo_df)), function(i) {
  cls <- geo_df$tdtom_label[i]
  list(
    type     = "Feature",
    geometry = list(
      type        = "Point",
      coordinates = list(geo_df$imagecol[i], geo_df$imagerow[i])
    ),
    properties = list(
      cell_id      = geo_df$cell_id[i],
      tdtom_label  = cls,
      tdtom_prob   = round(geo_df$tdtom_gmm_prob[i], 4),
      tdtom_expr   = round(geo_df$tdtom_expr[i], 4),
      color        = if (!is.na(cls) && cls %in% names(tdtom_colours))
                       tdtom_colours[[cls]] else "lightgrey"
    )
  )
})

gj <- list(type = "FeatureCollection", features = features)
gj_path <- file.path(OUT_DIR, paste0(case_id, "_tdtomato_gmm.geojson"))
write(toJSON(gj, auto_unbox = TRUE), gj_path)
message("  GeoJSON saved -> ", gj_path)

# =============================================================================
# 6. Save RDS
# =============================================================================

message("\n=== Step 6: Save RDS ===")

rds_out <- file.path(OUT_DIR, paste0(case_id, "_gmm_classified.rds"))
saveRDS(seu, rds_out)
message("  RDS saved -> ", rds_out)
message("\nDone.")
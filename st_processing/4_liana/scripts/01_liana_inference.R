################################################################################
# 01_liana_inference.R
# Project: Perivascular invasion paper
# Author: Simon Castillo
# Last update: 18/09/2026
################################################################################
##
##
# 0. preamble ----
##

suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(liana)
  library(dplyr);  library(readr);  library(ggplot2)
})

# Resolve paths independently of the current working directory.
.args_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(.args_file)) stop("Run this script with Rscript")
.script_path <- normalizePath(sub("^--file=", "", .args_file[[1]]), mustWork = TRUE)
PROJECT_ROOT <- dirname(dirname(.script_path))
.config_file <- Sys.getenv("LIANA_CONFIG", file.path(PROJECT_ROOT, "config", "config.R"))
if (!file.exists(.config_file)) stop("Configuration not found: ", .config_file)
source(.config_file, local = FALSE)

RUN_SAMPLES <- unname(SAMPLES)
_selected <- Sys.getenv("LIANA_SAMPLE", "")
if (nzchar(_selected)) {
  if (!(_selected %in% RUN_SAMPLES)) stop("Unknown LIANA_SAMPLE: ", _selected)
  RUN_SAMPLES <- _selected
}
if (!file.exists(RESOURCE_CSV)) stop("Missing curated mouse LR resource: ", RESOURCE_CSV)
curated <- read.csv(RESOURCE_CSV, stringsAsFactors = FALSE)
if (ncol(curated) != 2L) stop("Curated mouse resource must have exactly 2 columns: ", RESOURCE_CSV)
colnames(curated) <- c("source_genesymbol", "target_genesymbol")
OUT_CSV <- LIANA_RESULTS
OUT_FIG <- LIANA_FIGURES
OUT_LIANA <- LIANA_OBJECTS
N_PERM <- LIANA_N_PERM
SEED <- SEED_LIANA

for(SAMPLE in RUN_SAMPLES){ 
  OUT_LIANA_sample <- file.path(OUT_LIANA, SAMPLE)
  OUT_CSV_sample <- file.path(OUT_CSV, SAMPLE)
  OUT_FIG_sample <- file.path(OUT_FIG, SAMPLE)

  # ---- 1. load --------------------------------------------------------------
  obj  <- readRDS(file.path(IN_DIR, paste0(SAMPLE, "_expr_sparse.rds")))
  mat  <- obj$matrix
  meta <- as.data.frame(obj$metadata)
  stopifnot(inherits(mat, "dgCMatrix"), "cell_id" %in% names(meta))
  
  read_opt <- function(f) if (file.exists(f)) as.data.frame(read_csv(f, show_col_types = FALSE)) else NULL
  plot_df <- read_opt(file.path(IN_DIR, paste0(SAMPLE, "_plot_df.csv")))
  umap_df <- read_opt(file.path(IN_DIR, paste0(SAMPLE, "_umap_df.csv")))
  
  # ---- 2. align matrix columns to metadata rows -----------------------------
  # colnames(mat) are cell ids as characters; fall back to positional order.
  ids <- as.character(meta$cell_id)
  if (!is.null(colnames(mat)) && mean(colnames(mat) %in% ids) > 0.9) {
    keep <- intersect(colnames(mat), ids)
    mat  <- mat[, keep, drop = FALSE]
    meta <- meta[match(keep, ids), , drop = FALSE]
  } else {
    stopifnot(ncol(mat) == nrow(meta))
    message("colnames(mat) do not match cell_id - assuming column order == metadata row order.")
    colnames(mat) <- ids
  }
  if (anyDuplicated(ids)) stop("Duplicate cell IDs in expression metadata for ", SAMPLE)
  rownames(meta) <- colnames(mat)
  message(sprintf("%d genes x %d cells", nrow(mat), ncol(mat)))
  
  # merge extra per-cell tables (spatial coords, anything not already in meta) by cell id
  merge_by_id <- function(meta, df, tag) {
    if (is.null(df)) return(meta)
    idc <- intersect(c("cell_id", "cell", "cell_ID", "barcode", "id"), names(df))
    if (!length(idc)) { message("No cell-id column in ", tag, " - skipped."); return(meta) }
    i   <- match(rownames(meta), as.character(df[[idc[1]]]))
    add <- setdiff(names(df), c(names(meta), idc[1]))
    if (length(add)) meta[add] <- df[i, add, drop = FALSE]
    meta
  }
  meta <- merge_by_id(meta, plot_df, "plot_df")
  meta <- merge_by_id(meta, umap_df, "umap_df")
  
  # ---- 3. Seurat object (respect the pre-logged data) -----------------------
  if (DATA_IS_LOG) {
    # Reconstruct pseudo-counts solely to satisfy Seurat constructor.
    # SCT 'data' is NOT equivalent to log1p(raw counts): see README caveat.
    cnts <- mat; cnts@x <- round(pmax(expm1(cnts@x), 0))
  } else {
    cnts <- mat
  }
  so <- CreateSeuratObject(counts = cnts, assay = "Spatial", meta.data = meta,
                           project = SAMPLE, min.cells = 0, min.features = 0)
  if (DATA_IS_LOG) {
    if (packageVersion("Seurat") >= "5.0.0") {
      SeuratObject::LayerData(so, assay = "Spatial", layer = "data") <- mat
    } else {
      so <- SetAssayData(so, assay = "Spatial", slot = "data", new.data = mat)
    }
  } else {
    so <- NormalizeData(so, verbose = FALSE)
  }
  rm(cnts); invisible(gc())
  
  # ---- 4. define the groups -------------------------------------------------
  stopifnot(LABEL_COL %in% names(so@meta.data))
  so$liana_group <- as.character(so@meta.data[[LABEL_COL]])
  
  cells <- colnames(so)[!is.na(so$liana_group) & !(so$liana_group %in% DROP_LABELS)]
  if (!is.null(MAX_PVAL) && "pval_min" %in% names(so@meta.data)) {
    cells <- intersect(cells, colnames(so)[!is.na(so$pval_min) & so$pval_min <= MAX_PVAL])
  }
  so <- subset(so, cells = cells)
  
  tab  <- table(so$liana_group)
  small <- names(tab)[tab < MIN_CELLS]
  if (length(small)) {
    message("Dropping groups with < ", MIN_CELLS, " cells: ", paste(small, collapse = ", "))
    so <- subset(so, cells = colnames(so)[!(so$liana_group %in% small)])
  }
  so$liana_group <- factor(so$liana_group)
  Idents(so) <- "liana_group"
  print(table(so$liana_group))
  
  # ---- 5. ligand-receptor resource -----------------------------------------
  # Use the explicitly supplied curated mouse resource. Do not query or
  # regenerate another resource silently; the original LIANA call used "custom".
  message("Using custom ligand-receptor resource: ", RESOURCE_CSV)

  message('Starting with mode LIANA:')
  
  ## CCC modes ----
  modes <- MODES
  for (mode in modes) {
    print(mode)
    if(mode == 'E_T'){
      SUBSET_CLASS <- c( "AC_like",
                         'MES_like',
                         'OPC_like',
                         'NPC_like',
                         "Endothelial")
    }
    if(mode == 'M_T'){
      SUBSET_CLASS <- c( "AC_like",
                         'MES_like',
                         'OPC_like',
                         'NPC_like',
                         "Microglia_core", 
                         "Macrophage_BAM")
    }
    if(mode == 'M_E'){
      SUBSET_CLASS <- c( "Microglia_core", "Macrophage_BAM","Endothelial")
    }
    
    OUT_LIANA_mode <- file.path(OUT_LIANA_sample, mode)
    OUT_CSV_mode <- file.path(OUT_CSV_sample, mode)
    OUT_FIG_mode <- file.path(OUT_FIG_sample, mode)
    
    
    if (!dir.exists(OUT_LIANA_mode)) {
      dir.create(OUT_LIANA_mode, showWarnings = FALSE, recursive = TRUE)
    }
    
    if (!dir.exists(OUT_CSV_mode)) {
      dir.create(OUT_CSV_mode, showWarnings = FALSE, recursive = TRUE)
    }
    
    if (!dir.exists(OUT_FIG_mode)) {
      dir.create(OUT_FIG_mode, showWarnings = FALSE, recursive = TRUE)
    }
    
  
    set.seed(SEED)

    # ---- 6. run LIANA ---------------------------------------------------------
    # Respect LABEL_COL rather than hardcoding 'cell_class' in the mode subset.
    mode_cells <- colnames(so)[as.character(so$liana_group) %in% SUBSET_CLASS]
    if (length(mode_cells) < MIN_CELLS) {
      warning(sprintf("[%s %s] too few cells (%d); skipping mode", SAMPLE, mode, length(mode_cells)))
      next
    }
    so_subset <- subset(so, cells = mode_cells)
    if (length(unique(as.character(so_subset$liana_group))) < 2L) {
      warning(sprintf("[%s %s] <2 cell groups; skipping mode", SAMPLE, mode))
      next
    }                       
    so_subset$liana_group <- droplevels(so_subset$liana_group)
    Idents(so_subset) <- "liana_group"
    sce <- as.SingleCellExperiment(
      so_subset,
      assay = "Spatial"
    )
    if (!"logcounts" %in% SummarizedExperiment::assayNames(sce))
      stop("LIANA requires a logcounts assay. Check sparse export and data layer.")
    library(future)
    
    future::plan("sequential")
    
    liana_res <-liana_wrap(
      sce,
      method             = METHODS,
      #resource           = resource,
      #external_resource  = op_resource,
      resource          = "custom",
      external_resource = curated,
      idents_col         = "liana_group",
      assay              = "Spatial",
      assay.type         = "logcounts", 
      expr_prop          = EXPR_PROP,
      supp_columns = SUPP_COL,  
      permutation.params = list(nperms = N_PERM,
                                parallelize = N_WORKERS > 1,
                                workers = N_WORKERS),
      verbose            = TRUE
    )
    
    liana_agg <- liana_aggregate(liana_res, verbose = TRUE)
    top_25_quant <- quantile(liana_agg$aggregate_rank, c(0.25), na.rm = TRUE)
    top <- liana_agg |> filter(aggregate_rank <= top_25_quant, 
                               cellphonedb.pvalue< 0.05) |> 
      arrange(aggregate_rank)
    if(length(grep("Tlr*", top$receptor.complex, value = T))>0) {message("Yes, there's Tlr: ", unique(grep("Tlr*", top$receptor.complex, value = T)))}
    write_csv(liana_agg, file.path(OUT_CSV_mode, "liana_aggregate.csv"))
    write_csv(bind_rows(lapply(names(liana_res),
                               function(m) mutate(liana_res[[m]], method = m))),
              file.path(OUT_CSV_mode, "liana_all_methods.csv"))
    write_csv(top, file.path(OUT_CSV_mode, "liana_top_interactions.csv"))
    message(sprintf("%d interactions in top quartile passing CellPhoneDB p<0.05", nrow(top)))
    
    # ---- 7. figures -----------------------------------------------------------
    grp <- levels(droplevels(so_subset$liana_group))
    p_dot <- liana_dotplot(liana_agg,
                           source_groups = grp, target_groups = grp,
                           ntop = NTOP_DOT,
                           specificity = "natmi.edge_specificity",
                           magnitude   = "sca.LRscore") +
      ggtitle(sprintf("%s: top %d ranked interactions", SAMPLE, NTOP_DOT)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggsave(file.path(OUT_FIG_mode, "liana_dotplot.pdf"), plot = p_dot,
           width = 4 + 1.2 * length(grp), height = 10, limitsize = FALSE)
    
    if (nrow(top) > 0L) {
      pdf(file.path(OUT_FIG_mode, "liana_frequency_heatmap.pdf"), width = 7, height = 6)
      print(liana::heat_freq(top))
      dev.off()
    }
    
    ucols <- grep("^(umap|UMAP)_?[12]$", names(so@meta.data), value = TRUE)
    if (length(ucols) == 2) {
      um <- as.matrix(so@meta.data[, ucols]); colnames(um) <- c("UMAP_1", "UMAP_2")
      so[["umap"]] <- CreateDimReducObject(um, key = "UMAP_", assay = "Spatial")
      ggsave(file.path(OUT_FIG_mode, "liana_umap.pdf"),
             plot = DimPlot(so, reduction = "umap", group.by = "liana_group",
                            label = TRUE, repel = TRUE) + ggtitle(SAMPLE),
             width = 7.5, height = 6)
    }
    
    saveRDS(list(liana = liana_res, aggregate = liana_agg,
                 groups = table(so_subset$liana_group), resource = basename(RESOURCE_CSV),
                 config = list(label = LABEL_COL, expr_prop = EXPR_PROP,
                               nperms = N_PERM, organism = ORGANISM)),
            file.path(OUT_LIANA_mode, paste0(SAMPLE, "_liana.rds")))
    message("Done. Outputs in ", normalizePath(OUT_LIANA_mode))
    
    
    
    suppressPackageStartupMessages({
      library(dplyr); library(readr); library(ggplot2); library(liana)
    })
    
    AGG     <- file.path(OUT_CSV_mode, "liana_aggregate.csv")
    
    ALPHA        <- 0.01   # aggregate_rank cutoff (rank-aggregation p-value)
    USE_FDR      <- TRUE   # if TRUE, use BH-adjusted aggregate_rank <= FDR_ALPHA instead
    FDR_ALPHA    <- 0.05
    MIN_LRSCORE  <- 0.6    # sca.LRscore magnitude floor (0-1, ~0.5 = noise)
    CPDB_P       <- 0.05   # cellphonedb permutation p-value (NA-tolerant)
    TOP_PER_PAIR <- 10     # keep this many best interactions per source->target
    SELF_TALK    <- TRUE   # FALSE drops autocrine (source == target) edges
    
    liana_agg <- read_csv(AGG, show_col_types = FALSE)
    
    # ---- Extra 1. significance + magnitude -----------------------------------------
    res <- liana_agg |>
      mutate(qval = p.adjust(aggregate_rank, method = "BH"))
    
    sig <- res |>
      filter(if (USE_FDR) qval <= FDR_ALPHA else aggregate_rank <= ALPHA)
    
    if ("sca.LRscore" %in% names(sig))
      sig <- filter(sig, is.na(sca.LRscore) | sca.LRscore >= MIN_LRSCORE)
    if ("cellphonedb.pvalue" %in% names(sig))
      sig <- filter(sig, is.na(cellphonedb.pvalue) | cellphonedb.pvalue <= CPDB_P)
    if (!SELF_TALK)
      sig <- filter(sig, source != target)
    
    sig <- arrange(sig, aggregate_rank)
    
    # ---- Extra 2. depth-balanced shortlist per cell-type pair -----------------------
    top_pair <- sig |>
      group_by(source, target) |>
      #slice_min(aggregate_rank, n = TOP_PER_PAIR, with_ties = FALSE) |>
      ungroup() |>
      arrange(source, target, aggregate_rank)
    
    write_csv(sig,      file.path(OUT_CSV_mode, "liana_significant.csv"))
    write_csv(top_pair, file.path(OUT_CSV_mode, "liana_significant_top_per_pair.csv"))
    
    cat(sprintf("%d of %d interactions retained (%d source-target pairs)\n",
                nrow(sig), nrow(res), dplyr::n_distinct(paste(sig$source, sig$target))))
    print(sig |> count(source, target, name = "n_sig") |> arrange(desc(n_sig)))
    
    # ---- Extra 3. figures on the filtered set --------------------------------------
    if (nrow(top_pair) == 0L) {
      message("No retained interactions for ", SAMPLE, "/", mode, "; skipping filtered plots")
      next
    }
    p <- liana_dotplot(top_pair,
                       source_groups = unique(top_pair$source),
                       target_groups = unique(top_pair$target),
                       ntop = NULL,
                       specificity = "natmi.edge_specificity",
                       magnitude   = "sca.LRscore") +
      theme(axis.text.x = element_text(angle = 90, 
                                       hjust = 1, 
                                       vjust = 0.5)); print(p)
    
    ggsave(file.path(OUT_FIG_mode, "liana_significant_dotplot.pdf"), plot = p,
           width = 6 + 1.2 * dplyr::n_distinct(top_pair$target),
           height = 4 + 0.22 * nrow(top_pair), limitsize = FALSE)
    
    pdf(file.path(OUT_FIG_mode, "liana_significant_heatmap.pdf"), width = 7, height = 6)
    print(liana::heat_freq(sig))   # n significant interactions per source x target
    dev.off()
  }  
}


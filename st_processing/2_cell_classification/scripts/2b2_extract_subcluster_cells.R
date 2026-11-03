#!/usr/bin/env Rscript
# =============================================================================
# extract_subcluster_cells.R
#
# Recover the per-cell -> sub-cluster assignment from a completed
# tdtom_cluster_enrichment_v4.R run.
#
# WHY THIS EXISTS
#   The enrichment script builds `sub_cells` (cell, parent_cluster, sub_id) in
#   memory for Step 8b's Idents(), then discards it. Only the aggregated rate
#   table survives to disk. Anything downstream that needs to know WHICH CELLS
#   are in a sub-cluster - e.g. reclassifying TdTom_neg cells in tumour-like
#   sub-clusters - therefore has no join key.
#
#   This script replays ONLY the Step 8 sub-clustering loop, with the same seed
#   and the same parameters, and writes the mapping out. It does not rerun the
#   permutation tests, the enrichment, or any DE, so it costs a small fraction
#   of the full pipeline.
#
# WHY IT IS SAFE TO REPLAY
#   The loop calls set.seed(SEED) once per target cluster, immediately before
#   the PCA/neighbours/clustering. Each target is therefore independent of the
#   others and of iteration order, and with --subcluster_max_cells 0 the
#   sample() call never fires. Given identical inputs and package versions the
#   result is bit-identical.
#
#   That "given" is an assumption, not a guarantee - so this script VERIFIES it
#   rather than trusting it. Every replayed sub-cluster must match the run's own
#   _subcluster_tdtom.csv on both n_cells and n_tdtom_pos, exactly. If Seurat,
#   igraph or the BLAS changed under you, the counts diverge and the script
#   aborts instead of writing a mapping that silently belongs to a different
#   clustering than the table it will be joined to.
#
# USAGE
#   Rscript 2b2_extract_subcluster_cells.R \
#       --rds       <same --rds as the original run> \
#       --tdtomato  <same --tdtomato as the original run> \
#       --out_dir   <the run's output dir, containing _subcluster_tdtom.csv>
#
#   Every other flag defaults to the v4 default and only needs passing if the
#   original run overrode it.
#
# OUTPUT
#   <out_dir>/<case>_subcluster_cells.csv
#       cell_id, parent_cluster, sub_id, sub_id_full, tdtom_label
#   sub_id_full is "<parent>_<sub>" and matches sub_id_full in
#   _subcluster_markers.csv. Read it with colClasses="character": "7_1" and
#   "7_10" are distinct, and the underscore is what stops a reader from
#   collapsing them the way a dot would.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(stringr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--rds",       type = "character", default = NULL),
  make_option("--tdtomato",  type = "character", default = NULL),
  make_option("--out_dir",   type = "character", default = NULL),
  make_option("--use_predicted", type = "character", default = NULL),
  # These must match the original run. Defaults are v4's defaults.
  make_option("--cluster_col",          type = "character", default = "seurat_clusters"),
  make_option("--de_assay",             type = "character", default = "SCT"),
  make_option("--tdtom_feature_pattern", type = "character", default = "tdtom|td-tom|tdTomato"),
  make_option("--seed",                 type = "integer",   default = 123L),
  make_option("--subcluster_res",       type = "double",    default = 0.5),
  make_option("--subcluster_dims",      type = "integer",   default = 20L),
  make_option("--subcluster_min_cells", type = "integer",   default = 200L),
  make_option("--subcluster_max_cells", type = "integer",   default = 0L),
  make_option("--strict", type = "logical", default = TRUE,
              help = paste("Abort if the replayed sub-clusters do not reproduce",
                           "the run's _subcluster_tdtom.csv exactly [TRUE].",
                           "FALSE downgrades to a warning - only for inspecting",
                           "a mismatch, never for producing a mapping to use."))
)))

if (is.null(opt$rds) || is.null(opt$out_dir))
  stop("--rds and --out_dir are required")
if (is.null(opt$tdtomato) && is.null(opt$use_predicted))
  stop("Supply --tdtomato (rescued CSV) or --use_predicted, as in the original run")

SEED <- opt$seed
case_id <- str_extract(basename(opt$rds), "[0-9]{4}_Scan[0-9]+")
if (is.na(case_id)) stop("Could not parse <nnnn>_Scan<n> case id from --rds filename")
message("=== Replaying sub-clustering for ", case_id, " ===")

# --- the run's own rate table is the reference we validate against -----------
ref_path <- file.path(opt$out_dir, paste0(case_id, "_subcluster_tdtom.csv"))
if (!file.exists(ref_path))
  stop("Reference table not found: ", ref_path,
       "\n  The original run must have completed Step 8 before this can be replayed.")
ref <- read.csv(ref_path, colClasses = c(sub_id = "character",
                                         parent_cluster = "character"),
                stringsAsFactors = FALSE)
if (!nrow(ref)) stop("Reference table is empty: ", ref_path)

# Targets come from the run's OWN OUTPUT, not from recomputing the enrichment
# that chose them. Step 8 picks targets from the Step 5 enrichment results, so
# recomputing the choice would mean replaying the permutation tests as well -
# and any drift in those would silently change which clusters get sub-clustered.
# parent_cluster in the rate table is the authoritative record of what was
# actually sub-clustered in the run being recovered.
targets <- unique(ref$parent_cluster)
message("  Targets recovered from ", basename(ref_path), ": ",
        paste(targets, collapse = ", "))

if (any(ref$parent_downsampled %in% TRUE))
  warning("The original run downsampled at least one parent. sample() is seeded ",
          "and will replay identically, but the resulting sub-clusters cover only ",
          "part of the parent - cells outside the sample get no assignment.",
          call. = FALSE)

# --- Step 1 / 1b / 2 / 3: reproduce the object state exactly -----------------
message("  Loading ", opt$rds)
seu <- readRDS(opt$rds)
if ("SCT" %in% Assays(seu)) DefaultAssay(seu) <- "SCT"

if (!is.null(opt$use_predicted)) {
  pred_meta <- readRDS(opt$use_predicted)[[case_id]]
  if (is.null(pred_meta)) stop("Sample ", case_id, " not found in predicted meta RDS")
  shared <- intersect(colnames(seu), rownames(pred_meta))
  if (!length(shared)) stop("No shared barcodes with predicted meta")
  seu$tdtom_label <- NA_character_
  seu$tdtom_label[shared] <- pred_meta[shared, "tdtom_predicted"]
} else {
  td <- read.csv(opt$tdtomato, stringsAsFactors = FALSE)
  td$cell_id <- as.character(td$cell_id)
  if (anyDuplicated(td$cell_id)) stop("Duplicate cell_id values in TdTomato CSV")
  idx <- match(rownames(seu@meta.data), td$cell_id)
  if (anyNA(idx) || nrow(td) != ncol(seu))
    stop("TdTomato CSV does not cover exactly the Seurat cells: ",
         sum(is.na(idx)), " missing; expected ", ncol(seu),
         " rows, got ", nrow(td))
  for (col in setdiff(colnames(td), "cell_id"))
    seu@meta.data[[col]] <- td[[col]][idx]
  if (!"tdtom_label" %in% colnames(seu@meta.data))
    stop("tdtom_label not found after CSV merge")
}

seu$cluster_id <- if (opt$cluster_col %in% colnames(seu@meta.data))
  as.character(seu@meta.data[[opt$cluster_col]]) else as.character(Idents(seu))

# --- Step 8 loop, sub-clustering only ----------------------------------------
sub_cells <- tibble()
for (cl in targets) {
  cells <- colnames(seu)[as.character(seu$cluster_id) == cl & !is.na(seu$tdtom_label)]
  if (length(cells) < opt$subcluster_min_cells) {
    warning("cluster ", cl, ": only ", length(cells), " labelled cells now, below ",
            opt$subcluster_min_cells, " - it was sub-clustered in the original run, ",
            "so the inputs differ", call. = FALSE)
    next
  }
  set.seed(SEED)
  if (opt$subcluster_max_cells > 0 && length(cells) > opt$subcluster_max_cells)
    cells <- sample(cells, opt$subcluster_max_cells)

  sub <- subset(seu, cells = cells)
  DefaultAssay(sub) <- if (opt$de_assay %in% Assays(sub)) opt$de_assay else DefaultAssay(sub)
  if (DefaultAssay(sub) != "SCT") {
    sub <- NormalizeData(sub, verbose = FALSE)
    sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
    sub <- ScaleData(sub, verbose = FALSE)
  } else if (!length(VariableFeatures(sub))) {
    sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
  }
  vf <- setdiff(VariableFeatures(sub),
                grep(opt$tdtom_feature_pattern, rownames(sub),
                     ignore.case = TRUE, value = TRUE))
  npc <- min(opt$subcluster_dims, length(cells) - 1, length(vf) - 1)
  sub <- RunPCA(sub, features = vf, npcs = npc, verbose = FALSE)
  sub <- FindNeighbors(sub, dims = 1:npc, verbose = FALSE)
  sub <- FindClusters(sub, resolution = opt$subcluster_res, verbose = FALSE)

  sub_cells <- bind_rows(sub_cells, tibble(
    cell_id        = colnames(sub),
    parent_cluster = cl,
    sub_id         = as.character(Idents(sub)),
    tdtom_label    = sub$tdtom_label))
  message("    cluster ", cl, ": ", length(cells), " cells -> ",
          dplyr::n_distinct(Idents(sub)), " sub-clusters")
}

if (!nrow(sub_cells)) stop("No sub-clusters were produced - nothing to write")
sub_cells$sub_id_full <- paste0(sub_cells$parent_cluster, "_", sub_cells$sub_id)

# --- VALIDATION: does this replay describe the run we are recovering? --------
# Compare on the two quantities the rate table records per sub-cluster. Matching
# both, for every sub-cluster, means the same cells landed in the same groups.
rep_tbl <- sub_cells %>%
  group_by(parent_cluster, sub_id) %>%
  summarise(n_cells_replay = n(),
            n_pos_replay   = sum(tdtom_label == "TdTom_pos"), .groups = "drop")

cmp <- ref %>%
  transmute(parent_cluster, sub_id, n_cells_ref = n_cells, n_pos_ref = n_tdtom_pos) %>%
  full_join(rep_tbl, by = c("parent_cluster", "sub_id"))

bad <- cmp %>% filter(is.na(n_cells_ref) | is.na(n_cells_replay) |
                      n_cells_ref != n_cells_replay | n_pos_ref != n_pos_replay)

message("\n=== Validation against ", basename(ref_path), " ===")
message("  sub-clusters in run:    ", nrow(ref))
message("  sub-clusters replayed:  ", nrow(rep_tbl))
message("  exact matches:          ", nrow(cmp) - nrow(bad), " / ", nrow(cmp))

if (nrow(bad)) {
  print(as.data.frame(bad))
  msg <- paste0("Replay did NOT reproduce the original sub-clustering (",
                nrow(bad), " mismatched sub-cluster(s)). The mapping does not ",
                "correspond to the tables it would be joined to. Most likely a ",
                "Seurat/igraph/BLAS version change since the run, or different ",
                "--subcluster_* parameters than the original invocation.")
  if (isTRUE(opt$strict)) stop(msg) else warning(msg, call. = FALSE)
} else {
  message("  PASS - every sub-cluster reproduced exactly on size and TdTom+ count.")
}

# A mapping that failed validation must not be able to masquerade as a good one
# later: --strict FALSE is for inspecting a mismatch, so its output gets a name
# nothing downstream will pick up by accident.
suffix <- if (nrow(bad)) "_subcluster_cells_UNVALIDATED.csv" else "_subcluster_cells.csv"
out <- file.path(opt$out_dir, paste0(case_id, suffix))
if (nrow(bad))
  message("  Writing under an UNVALIDATED name — do not use this for rescue.")
write.csv(sub_cells %>%
            transmute(cell_id, parent_cluster, sub_id, sub_id_full, tdtom_label),
          out, row.names = FALSE)
message("\n=== Wrote ", out, " (", nrow(sub_cells), " cells) ===")
message("Read it with colClasses = c(sub_id_full = \"character\") - the ids are ",
        "not numbers.")

#!/usr/bin/env Rscript
# =============================================================================
# Cluster x TdTomato enrichment — which cluster(s) are enriched for TdTom_pos?
# =============================================================================
# This script reuses Steps 1 / 1b / 2 of the Neftel classification pipeline
# verbatim (load RDS -> merge rescued TdTomato CSV -> define tdtom_label), then
# replaces everything downstream with a single question:
#
#     is the TdTomato+ fraction of a cluster higher than expected given the
#     sample-wide TdTomato+ rate, and by how much?
#
# DELIBERATELY NOT INCLUDED (present in the full pipeline, removed here):
#   - signature scoring (Neftel / cell type / proliferation / EC / myeloid)
#   - expression x detection matched permutation nulls and their calibration
#   - hybrid state calling, winner-takes-all cell-class assignment
#   - every signature-derived GeoJSON layer (neftel_normal,
#     proliferation, intersection, state_label, ec_state). The two layers that
#     survive are the two that need no signature: tdtomato and seurat_cluster.
# Nothing here depends on a signature, so no null model is needed: the
# comparison is count-based (cluster membership vs TdTomato label), and the
# reference is the marginal TdTomato+ rate of the sample itself.
#
# -----------------------------------------------------------------------------
# WHY THE RUN DOES NOT STOP AT THE MATRIX
# -----------------------------------------------------------------------------
# The enrichment test locates TdTom+ cells. Three further questions decide
# whether an enriched cluster is usable as a tumour definition, and each is a
# step below:
#
#   Step 6, purity.  Enrichment is not purity. Given a background TdTomato+
#     rate fp (estimated from the strongly depleted clusters) and an assumed
#     detection sensitivity, the implied true tumour fraction of a cluster is
#     f = (obs - fp)/(sens - fp). Sensitivity is not estimable here, so it is
#     swept. Read f at sens = 1: if even a perfect gate implies f well below 1,
#     the cluster is a MIXTURE and no threshold on this label will purify it.
#     Step 6 also tabulates the purity/completeness trade-off over the candidate
#     cluster rules, because a rule that is pure usually leaves most TdTom+
#     cells outside it.
#
#   Step 7, DE.  Two modes, chosen with --de_mode:
#
#     all_vs_rest (default) — every tested cluster vs all other cells. These are
#       cluster markers: they say what each cluster IS. Read this first. An
#       enrichment call is only interpretable once you know the identity of the
#       cluster carrying it, and the marker table is also what settles an
#       ambiguous cluster: if a cluster sitting at the sample-wide TdTomato rate
#       turns out to be a large stromal or immune population, its share of the
#       TdTom+ pool is most likely background; if it looks like the enriched
#       clusters, the tumour programme is leaking outside them.
#
#     targeted — the contrasts that test the tumour question directly
#       (enriched vs reference, enriched vs enriched, within-cluster label).
#
#     both — runs the marker contrasts and the targeted ones.
#
#   The comparator differs between the two modes ON PURPOSE. all_vs_rest uses
#   "everything else", which is correct for marker discovery. The targeted
#   contrasts do NOT, and that is the point:
#
#   Two design choices there matter more than the test used:
#       * the reference is the DEPLETED clusters, not "everything not enriched".
#         An `ns` cluster sitting at the sample-wide rate can still hold a large
#         share of all TdTom+ cells; pooling it into the comparator puts tumour
#         on both sides and flattens the genes being looked for.
#       * enriched clusters are tested separately AND against each other. Each
#         vs reference asks "is this tumour?"; enriched-vs-enriched asks whether
#         two enriched clusters are one population — the contrast that decides
#         whether to merge them.
#     Within-cluster TdTom+ vs TdTom- is also run: if the label adds nothing
#     inside an enriched cluster, it is adding noise, not resolution. The
#     TdTomato feature is flagged and excluded from ranked plots throughout —
#     the labels derive from it, so its rank is circular.
#
#   Step 8, sub-cluster diagnostic.  If a cluster is a mixture, either the
#     resolution is too coarse (a sub-population carries the tumour) or the
#     label is lossy (TdTom+ cells are spread evenly). Re-clustering each target
#     in isolation — with the TdTomato feature removed from the variable
#     features, so the sub-structure is not defined by the label being tested —
#     distinguishes these.
#
# All three are switchable (--de, --subcluster) and none changes the enrichment
# result: the matrix and its calls are computed before any of them run.
#
# -----------------------------------------------------------------------------
# WHAT THE MATRIX IS
# -----------------------------------------------------------------------------
# A cluster x label contingency table (rows = clusters, columns = TdTom_pos /
# TdTom_neg) rendered three ways, because the three answer different questions:
#
#   (1) log2 observed/expected      "is this cell of the table over-occupied?"
#         E_ij = row_i * col_j / N  (independence expectation)
#         log2((O + 0.5) / (E + 0.5)) — the 0.5 pseudocount keeps empty cells
#         finite; it shrinks small counts toward 0, so it is conservative.
#   (2) standardised (adjusted) Pearson residual
#         (O - E) / sqrt(E * (1 - row_i/N) * (1 - col_j/N))
#         ~ N(0,1) under independence, so |r| > 2 is the usual eyeball cut.
#         This is the calibrated version of (1) — read effect size off (1),
#         read "is it real" off (2) or off the per-cluster test below.
#   (3) row-normalised composition (% of the cluster's cells)
#         the descriptive view; no test attached.
#
# -----------------------------------------------------------------------------
# THE PER-CLUSTER TEST
# -----------------------------------------------------------------------------
# For each cluster: 2x2 Fisher exact test of (this cluster vs all other cells)
# x (TdTom_pos vs TdTom_neg). Fisher rather than chi-square because small
# clusters and skewed TdTomato rates put expected counts below 5 routinely.
# Reported per cluster: two-sided p (the omnibus "different rate"), one-sided
# greater p (the directional "enriched" question), conditional odds ratio and
# its exact CI. BH correction is applied across the TESTED clusters only.
#
# Clusters with fewer than --min_cluster_cells cells are reported in the table
# and drawn in the heatmap but are NOT tested and NOT counted in the BH
# denominator (tested = FALSE). They have no power at any realistic TdTomato
# rate, so testing them only inflates the multiple-testing burden for the
# clusters that can actually answer the question. This filter uses cluster size
# only — it never looks at the TdTomato labels — so it does not bias the tests
# it keeps.
#
# A global chi-square test of the full table is also reported (Monte-Carlo p
# when any expected count < 5) as the "is there any cluster structure in
# TdTomato at all" sanity check that precedes the per-cluster calls.
#
# CALLS: q_two_sided < ALPHA and log2 O/E for TdTom_pos > 0  -> "enriched"
#        q_two_sided < ALPHA and log2 O/E for TdTom_pos < 0  -> "depleted"
#        otherwise                                           -> "ns"
#
# -----------------------------------------------------------------------------
# ONE CAVEAT WORTH KEEPING IN VIEW
# -----------------------------------------------------------------------------
# Clusters were derived from the expression matrix; if TdTomato transcript
# counts contributed to the HVGs/PCA used for clustering, a cluster can be
# enriched for TdTom_pos partly because TdTomato helped define it. That is a
# property of the upstream clustering, not of this test. --check_tdtom_feature
# (default TRUE) reports whether a TdTomato feature is present in the object
# and whether it is in VariableFeatures, so the caveat is visible rather than
# silent. With --use_predicted the labels come from a model instead, and the
# same circularity applies through whatever features that model used.
#
# =============================================================================
# Usage
# -----------------------------------------------------------------------------
#   Rscript tdtom_cluster_enrichment.R \
#     --rds       /path/to/1234_Scan1_classified.rds \
#     --tdtomato  /path/to/rescued_tdtomato.csv \
#     --out_dir   /path/to/out
#
#   Optional:
#     --use_predicted /path/to/t0_validation_meta_t2model.rds  # predicted labels
#     --cluster_col   seurat_clusters   # metadata column holding clusters
#     --min_cluster_cells 10            # below this: reported, not tested
#     --order_by      numeric|enrichment|size
#     --alpha         0.05
#     --save_rds      FALSE             # write object with tdtom_label attached
#
# Outputs (<case> = e.g. 1234_Scan1)
#   <case>_classification_plots.pdf          matrices, ranked, purity, sub-cluster,
#                                            DE overview + volcanoes, UMAP, spatial
#   <case>_cluster_tdtom_matrix.csv          wide contingency + O/E + residuals
#   <case>_cluster_tdtom_enrichment.csv      per-cluster Fisher stats + calls
#   <case>_cluster_tdtom_long.csv            tidy cluster x label cell values
#   <case>_cluster_purity.csv                implied tumour fraction per cluster
#                                            at several assumed sensitivities
#   <case>_candidate_rules.csv               purity vs completeness per rule
#   <case>_cluster_markers.csv               top up-genes per cluster vs rest
#   <case>_DE_summary.csv                    one row per contrast: n sig, top genes
#   <case>_DE_all_contrasts.csv              every gene x contrast (long)
#   <case>_DE_overlap.csv                    up-gene sharing across contrasts
#   <case>_DE_<contrast>.csv                 per-contrast files (--de_split_files)
#   <case>_subcluster_tdtom.csv              sub-cluster TdTom+ rates per target
#   <case>_subcluster_markers.csv            sub-cluster marker DE (Step 8b);
#                                            comparator column = section | parent
#
# NOTE ON SUB-CLUSTER IDS: sub_id_full is "<parent>_<sub>" (e.g. 7_10), NOT
# "<parent>.<sub>". With a dot, "7.1" and "7.10" are the same number, so any
# reader that type-infers the column (read.csv, pandas.read_csv, Excel) merges
# those two sub-clusters into one row and silently loses the other. Underscore
# cannot parse as numeric, so the id survives type inference everywhere.
#
# NOTE ON SUB-CLUSTER CELL COUNTS: --subcluster_max_cells defaults to 0 (no cap).
# If it is set, the n_cells of every sub-cluster under a capped parent is a SAMPLE
# of the real population, not its size. Those runs are marked in
# _subcluster_tdtom.csv and _subcluster_markers.csv by parent_downsampled = TRUE,
# with parent_n_labelled giving the true parent size; rates are unaffected,
# absolute counts are not comparable to an uncapped run.
#   <case>_tdtomato.geojson                  per-cell points coloured by label
#   <case>_seurat_cluster.geojson            per-cell points coloured by cluster,
#                                            each carrying that cluster's
#                                            log2 O/E, q and call as properties
#   <case>_run_parameters.csv                settings that produced the above
#
# Both GeoJSON layers use the same point convention as the full pipeline —
# coordinates are (imagecol, imagerow) and every feature carries a `color`
# property plus a `classification` property holding the value it is coloured by,
# so the two files behave identically in the viewer. They are written only if
# imagerow/imagecol are available (metadata or the image slot); the run reports
# a skip rather than failing if they are not.
# =============================================================================

suppressPackageStartupMessages(library(optparse))

opt_list <- list(
  make_option("--rds",      type = "character", help = "Path to classified RDS"),
  make_option("--tdtomato", type = "character", default = NULL,
              help = "Path to rescued TdTomato CSV (must contain cell_id)"),
  make_option("--out_dir",  type = "character", help = "Output directory"),
  make_option("--lib_path", type = "character", default = NULL),
  make_option("--use_predicted", type = "character", default = NULL,
              help = "Path to t0_validation_meta_t2model.rds — use predicted TdTomato labels"),
  make_option("--cluster_col", type = "character", default = "seurat_clusters",
              help = "Metadata column with cluster IDs; falls back to Idents() [seurat_clusters]"),
  make_option("--min_cluster_cells", type = "integer", default = 10L,
              help = "Clusters smaller than this are reported but not tested [10]"),
  make_option("--order_by", type = "character", default = "numeric",
              help = "Row order of the matrix: numeric | enrichment | size [numeric]"),
  make_option("--alpha", type = "double", default = 0.05),
  make_option("--seed",  type = "integer", default = 123L),
  make_option("--reduction", type = "character", default = NULL,
              help = "UMAP reduction name; auto-detected if omitted"),
  make_option("--check_tdtom_feature", type = "logical", default = TRUE,
              help = "Report whether a TdTomato feature is a variable feature [TRUE]"),
  make_option("--geojson", type = "logical", default = TRUE,
              help = "Write tdtomato.geojson and seurat_cluster.geojson [TRUE]"),
  make_option("--save_rds", type = "logical", default = FALSE,
              help = "Save the object with tdtom_label / cluster call attached [FALSE]"),

  # --- purity model ---------------------------------------------------------
  make_option("--fp_l2oe_cut", type = "double", default = -2.0,
              help = paste("Depleted clusters below this log2 O/E form the TdTomato",
                           "false-positive floor [-2.0]")),
  make_option("--tdtom_feature_pattern", type = "character", default = "tdtom|td-tom|tdTomato",
              help = "Regex matching the TdTomato feature name(s) [tdtom|td-tom|tdTomato]"),

  # --- differential expression ----------------------------------------------
  make_option("--de", type = "logical", default = TRUE,
              help = "Run the DE step [TRUE]"),
  make_option("--de_mode", type = "character", default = "all_vs_rest",
              help = paste("all_vs_rest | targeted | both.",
                           "all_vs_rest = every tested cluster vs all other cells",
                           "(cluster markers). targeted = enriched-vs-reference,",
                           "enriched-vs-enriched, within-cluster. [all_vs_rest]")),
  make_option("--de_only_pos", type = "logical", default = FALSE,
              help = "all_vs_rest: keep only genes up in the cluster [FALSE]"),
  make_option("--de_marker_top_n", type = "integer", default = 5L,
              help = "Markers per cluster in the marker heatmap and table [5]"),
  make_option("--de_reference", type = "character", default = "depleted",
              help = paste("Comparator: depleted | nonenriched | explicit.",
                           "'depleted' keeps ns clusters OUT of the reference",
                           "(they may contain tumour) [depleted]")),
  make_option("--de_ref_clusters", type = "character", default = NULL,
              help = "Comma-separated cluster IDs; required with --de_reference explicit"),
  make_option("--de_assay", type = "character", default = "SCT",
              help = "Assay for DE; falls back to DefaultAssay if absent [SCT]"),
  make_option("--de_test", type = "character", default = "wilcox",
              help = "FindMarkers test.use [wilcox]"),
  make_option("--de_min_pct", type = "double", default = 0.1,
              help = "FindMarkers min.pct [0.1]"),
  make_option("--de_logfc", type = "double", default = 0.25,
              help = "FindMarkers logfc.threshold [0.25]"),
  make_option("--de_min_cells", type = "integer", default = 30L,
              help = "Skip a contrast if either side has fewer cells than this [30]"),
  make_option("--de_max_cells", type = "integer", default = 3000L,
              help = paste("Downsample each side to at most this many cells;",
                           "0 = no cap. Effect sizes are unaffected [3000]")),
  make_option("--de_pairwise", type = "logical", default = TRUE,
              help = "Test enriched clusters against each other (merge question) [TRUE]"),
  make_option("--de_ambiguous", type = "logical", default = TRUE,
              help = "Also test ns clusters holding a large share of TdTom+ cells [TRUE]"),
  make_option("--de_within", type = "logical", default = TRUE,
              help = "Test TdTom+ vs TdTom- inside each target cluster [TRUE]"),
  make_option("--de_ns_pool_pct", type = "double", default = 10.0,
              help = "An ns cluster is tested if it holds >= this %% of all TdTom+ cells [10]"),
  make_option("--de_merge_hint", type = "integer", default = 50L,
              help = "Below this many DE genes, enriched clusters are flagged as mergeable [50]"),
  make_option("--de_split_files", type = "logical", default = FALSE,
              help = "Also write one CSV per contrast [FALSE]"),
  make_option("--de_top_n", type = "integer", default = 12L,
              help = "Genes labelled per direction in the DE plots [12]"),

  # --- sub-cluster diagnostic -----------------------------------------------
  make_option("--subcluster", type = "logical", default = TRUE,
              help = "Re-cluster target clusters to test the mixture explanation [TRUE]"),
  make_option("--subcluster_clusters", type = "character", default = NULL,
              help = "Comma-separated cluster IDs; default = enriched + large ns"),
  make_option("--subcluster_res", type = "double", default = 0.5,
              help = "FindClusters resolution for the sub-clustering [0.5]"),
  make_option("--subcluster_dims", type = "integer", default = 20L,
              help = "PCs used in the sub-clustering [20]"),
  make_option("--subcluster_min_cells", type = "integer", default = 200L,
              help = "Skip sub-clustering a cluster below this size [200]"),
  make_option("--subcluster_max_cells", type = "integer", default = 0L,
              help = paste("Downsample each target to this many cells before",
                           "sub-clustering; 0 = no cap [0]. A cap makes every",
                           "sub-cluster n_cells for that parent a sample of the",
                           "real population, so it defaults off; set it only if",
                           "a large cluster makes re-clustering too slow.")),
  make_option("--subcluster_gain", type = "double", default = 1.5,
              help = paste("Best sub-cluster must exceed parent rate by this factor",
                           "to call resolution limiting [1.5]")),
  make_option("--subcluster_de", type = "logical", default = TRUE,
              help = paste("Step 8b: marker DE with Idents set to sub-cluster identity,",
                           "so sub-clusters can be scored on transcriptional grounds",
                           "rather than on reporter rate alone [TRUE]")),
  make_option("--subcluster_de_mode", type = "character", default = "both",
              help = paste("Comparator for the sub-cluster marker DE:",
                           "'section' = sub vs all other cells in the section (the",
                           "analogue of cluster-vs-rest, and the one to score against",
                           "cell-type panels); 'parent' = sub vs the rest of its own",
                           "parent cluster (what makes this sub-cluster different",
                           "from its siblings); 'both' [both]")),
  make_option("--subcluster_de_min_cells", type = "integer", default = 50L,
              help = paste("Skip a sub-cluster below this size in Step 8b. Lower than",
                           "--de_min_cells on purpose: sub-clusters are small by",
                           "construction and a 200-cell floor would drop most of them [50]"))
)

opt <- parse_args(OptionParser(option_list = opt_list))
for (arg in c("rds", "out_dir"))
  if (is.null(opt[[arg]])) stop("Missing required argument: --", arg)
if (!is.null(opt$lib_path)) .libPaths(opt$lib_path)

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(jsonlite)
})

options(future.globals.maxSize = Inf)

RDS_PATH      <- opt$rds
TDTOMATO_PATH <- opt$tdtomato
OUT_DIR       <- opt$out_dir
USE_PREDICTED <- opt$use_predicted
CLUSTER_COL   <- opt$cluster_col
MIN_CL_CELLS  <- opt$min_cluster_cells
ORDER_BY      <- opt$order_by
ALPHA         <- opt$alpha
SEED          <- opt$seed
REDUCTION     <- opt$reduction
CHECK_FEATURE <- isTRUE(opt$check_tdtom_feature)
DO_GEOJSON    <- isTRUE(opt$geojson)
SAVE_RDS      <- isTRUE(opt$save_rds)

FP_L2OE_CUT   <- opt$fp_l2oe_cut
TDTOM_FEATURE_PATTERN <- opt$tdtom_feature_pattern

DO_DE            <- isTRUE(opt$de)
DE_MODE          <- opt$de_mode
DE_ONLY_POS      <- isTRUE(opt$de_only_pos)
DE_MARKER_TOP_N  <- opt$de_marker_top_n
DE_REFERENCE     <- opt$de_reference
DE_REF_EXPLICIT  <- if (is.null(opt$de_ref_clusters)) character(0) else
  trimws(strsplit(opt$de_ref_clusters, ",")[[1]])
DE_ASSAY         <- opt$de_assay
DE_TEST          <- opt$de_test
DE_MIN_PCT       <- opt$de_min_pct
DE_LOGFC         <- opt$de_logfc
DE_MIN_CELLS     <- opt$de_min_cells
DE_MAX_CELLS     <- opt$de_max_cells
DE_PAIRWISE      <- isTRUE(opt$de_pairwise)
DE_AMBIGUOUS     <- isTRUE(opt$de_ambiguous)
DE_WITHIN        <- isTRUE(opt$de_within)
DE_NS_POOL_PCT   <- opt$de_ns_pool_pct
DE_MERGE_HINT    <- opt$de_merge_hint
DE_SPLIT_FILES   <- isTRUE(opt$de_split_files)
DE_TOP_N         <- opt$de_top_n

DO_SUBCLUSTER        <- isTRUE(opt$subcluster)
SUBCLUSTER_CLUSTERS  <- if (is.null(opt$subcluster_clusters)) character(0) else
  trimws(strsplit(opt$subcluster_clusters, ",")[[1]])
SUBCLUSTER_RES       <- opt$subcluster_res
SUBCLUSTER_DIMS      <- opt$subcluster_dims
SUBCLUSTER_MIN_CELLS <- opt$subcluster_min_cells
SUBCLUSTER_MAX_CELLS <- opt$subcluster_max_cells
SUBCLUSTER_GAIN      <- opt$subcluster_gain
SUBCLUSTER_DE        <- isTRUE(opt$subcluster_de)
SUBCLUSTER_DE_MODE   <- opt$subcluster_de_mode
SUBCLUSTER_DE_MIN    <- opt$subcluster_de_min_cells

if (!ORDER_BY %in% c("numeric", "enrichment", "size"))
  stop("--order_by must be 'numeric', 'enrichment' or 'size' (got '", ORDER_BY, "')")
if (!DE_MODE %in% c("all_vs_rest", "targeted", "both"))
  stop("--de_mode must be 'all_vs_rest', 'targeted' or 'both' (got '", DE_MODE, "')")
if (!DE_REFERENCE %in% c("depleted", "nonenriched", "explicit"))
  stop("--de_reference must be 'depleted', 'nonenriched' or 'explicit' (got '",
       DE_REFERENCE, "')")
if (DE_REFERENCE == "explicit" && !length(DE_REF_EXPLICIT))
  stop("--de_reference explicit requires --de_ref_clusters")
if (!SUBCLUSTER_DE_MODE %in% c("section", "parent", "both"))
  stop("--subcluster_de_mode must be 'section', 'parent' or 'both' (got '",
       SUBCLUSTER_DE_MODE, "')")
# Step 8b consumes the sub-cluster assignments Step 8 produces, so it cannot run
# on its own. Fail loudly rather than silently writing no sub-cluster markers.
if (SUBCLUSTER_DE && !DO_SUBCLUSTER)
  stop("--subcluster_de TRUE requires --subcluster TRUE (Step 8b needs Step 8's ",
       "sub-cluster assignments)")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

TDTOM_COLOURS <- c(TdTom_pos = "#1f4fd8", TdTom_neg = "#c9c9c9")
CALL_COLOURS  <- c(enriched = "#b2182b", depleted = "#2166ac",
                   ns = "grey70", untested = "grey88")

message("=== Parameters ===")
message("  RDS:            ", RDS_PATH)
message("  OUT_DIR:        ", OUT_DIR)
message("  TdTomato CSV:   ", ifelse(is.null(TDTOMATO_PATH), "None", TDTOMATO_PATH))
message("  USE_PREDICTED:  ", ifelse(is.null(USE_PREDICTED), "No (gene expression / CSV)", USE_PREDICTED))
message("  Cluster column: ", CLUSTER_COL)
message("  Min cluster n:  ", MIN_CL_CELLS, " (smaller clusters reported, not tested)")
message("  Alpha:          ", ALPHA, "  (BH across tested clusters)")

# =============================================================================
# 1. Load Seurat object
# =============================================================================

message("\n=== Step 1: Loading Seurat object ===")
seu     <- readRDS(RDS_PATH)
seu     <- UpdateSeuratObject(seu)
case_id <- str_extract(basename(RDS_PATH), "[0-9]{4}_Scan[0-9]+")
if (is.na(case_id)) case_id <- tools::file_path_sans_ext(basename(RDS_PATH))
message("  Sample: ", case_id, "  Cells: ", ncol(seu))
# The full pipeline hard-set this to "SCT". Kept as the preference, but falling
# back rather than erroring so the script also runs on objects normalised
# another way (the enrichment step is count-based and assay-independent).
if ("SCT" %in% Assays(seu)) {
  DefaultAssay(seu) <- "SCT"
} else {
  message("  No SCT assay; using DefaultAssay '", DefaultAssay(seu), "'")
}

# =============================================================================
# 1b. Update Seurat object TdTomato gating
# =============================================================================

message("\n=== Step 1b: Updating TdTomato metadata ===")

if (!is.null(TDTOMATO_PATH)) {

  rescued_tdtomato_df <- read.csv(TDTOMATO_PATH, stringsAsFactors = FALSE)
  rescued_tdtomato_df$cell_id <- as.character(rescued_tdtomato_df$cell_id)

  seu_cell_ids <- rownames(seu@meta.data)

  if (anyDuplicated(rescued_tdtomato_df$cell_id))
    stop("Duplicate cell_id values found in TdTomato CSV")

  idx <- match(as.character(seu_cell_ids), rescued_tdtomato_df$cell_id)

  message("  Seurat cells: ", length(seu_cell_ids))
  message("  CSV cells:    ", nrow(rescued_tdtomato_df))
  message("  Matched:      ", sum(!is.na(idx)))
  message("  Unmatched:    ", sum(is.na(idx)))

  if (anyNA(idx) || nrow(rescued_tdtomato_df) != length(seu_cell_ids))
    stop("TdTomato CSV must cover exactly the Seurat cells: ", sum(is.na(idx)),
         " missing in CSV, ", nrow(rescued_tdtomato_df), " CSV rows vs ",
         length(seu_cell_ids), " Seurat cells. Check the sample and source.")
  if (!all(rescued_tdtomato_df$tdtom_label %in% c("TdTom_pos", "TdTom_neg")))
    stop("TdTomato CSV contains labels outside TdTom_pos/TdTom_neg")

  new_cols <- setdiff(colnames(rescued_tdtomato_df), "cell_id")
  for (col in new_cols) seu@meta.data[[col]] <- rescued_tdtomato_df[[col]][idx]

} else {

  if (is.null(USE_PREDICTED))
    stop("Supply --tdtomato (rescued CSV) or --use_predicted (model labels)")
  message("  No --tdtomato given; labels will come from --use_predicted")

}

# =============================================================================
# 2. Define TdTomato+ / TdTomato- cells
# =============================================================================

message("\n=== Step 2: Defining TdTomato labels ===")

if (!is.null(USE_PREDICTED)) {

  message("  Using predicted labels from: ", USE_PREDICTED)
  pred_list <- readRDS(USE_PREDICTED)
  pred_meta <- pred_list[[case_id]]
  if (is.null(pred_meta))
    stop("Sample ", case_id, " not found in predicted meta RDS.")

  shared <- intersect(colnames(seu), rownames(pred_meta))
  if (length(shared) == 0)
    stop("No shared cell barcodes between Seurat object and predicted meta")

  seu$tdtom_label <- NA_character_
  seu$tdtom_label[shared] <- pred_meta[shared, "tdtom_predicted"]
  seu$tdtom_source <- "predicted"

} else {

  if (!"tdtom_label" %in% colnames(seu@meta.data))
    stop("tdtom_label not found in Seurat metadata after CSV merge")
  seu$tdtom_source <- "rescued_csv"

}

tumour_cells    <- rownames(seu@meta.data)[which(seu$tdtom_label == "TdTom_pos")]
nontumour_cells <- rownames(seu@meta.data)[which(seu$tdtom_label == "TdTom_neg")]

message("  TdTom+ (tumour):     ", length(tumour_cells))
message("  TdTom- (non-tumour): ", length(nontumour_cells))

if (length(tumour_cells) == 0) stop("No TdTomato+ cells found")
if (length(nontumour_cells) == 0)
  stop("No TdTomato- cells found — an enrichment test needs both labels")

# =============================================================================
# 3. Cluster labels
# =============================================================================

message("\n=== Step 3: Resolving cluster labels ===")

if (CLUSTER_COL %in% colnames(seu@meta.data)) {
  seu$cluster_id <- as.character(seu@meta.data[[CLUSTER_COL]])
  message("  Using metadata column: ", CLUSTER_COL)
} else {
  seu$cluster_id <- as.character(Idents(seu))
  message("  Column '", CLUSTER_COL, "' not found — using Idents(seu) instead")
}

# Cells missing either coordinate of the table cannot contribute to it. They are
# counted here and excluded, rather than silently dropped by table().
keep <- !is.na(seu$cluster_id) &
        seu$tdtom_label %in% c("TdTom_pos", "TdTom_neg")
n_dropped <- sum(!keep)
message("  Cells usable: ", sum(keep), "   dropped (NA cluster or NA/other label): ", n_dropped)

meta <- seu@meta.data %>%
  rownames_to_column("cell_id") %>%
  filter(keep) %>%
  mutate(tdtom_label = factor(tdtom_label, levels = c("TdTom_pos", "TdTom_neg")))

# Cluster level order: numeric where the IDs are numeric, else lexical.
cl_levels <- unique(meta$cluster_id)
suppressWarnings(cl_num <- as.numeric(cl_levels))
cl_levels <- if (!any(is.na(cl_num))) cl_levels[order(cl_num)] else sort(cl_levels)
meta$cluster_id <- factor(meta$cluster_id, levels = cl_levels)

message("  Clusters: ", length(cl_levels), "  (", paste(head(cl_levels, 15), collapse = ", "),
        if (length(cl_levels) > 15) ", ..." else "", ")")

# One palette for clusters, defined once and reused by the UMAP panels, the
# spatial panels and the GeoJSON export, so cluster 7 is the same colour in the
# PDF and in the viewer. scales::hue_pal() reproduces the ggplot default hues.
cluster_colours <- setNames(scales::hue_pal()(length(cl_levels)), cl_levels)

# --- circularity check: did TdTomato itself help define the clusters? --------
tdtom_feature_note <- "not checked"
if (CHECK_FEATURE) {
  # rownames() of an assay can fail on unusual assay classes; the check is
  # diagnostic only, so a failure must not take the run down with it.
  all_feats <- unique(unlist(lapply(Assays(seu), function(a)
    tryCatch(rownames(seu[[a]]), error = function(e) character(0)))))
  tdtom_feats <- grep("tdtom|td-tom|tdtomato", all_feats, ignore.case = TRUE, value = TRUE)
  if (!length(tdtom_feats)) {
    tdtom_feature_note <- "no TdTomato feature in object"
  } else {
    vf <- tryCatch(VariableFeatures(seu), error = function(e) character(0))
    in_vf <- intersect(tdtom_feats, vf)
    tdtom_feature_note <- paste0(
      "features: ", paste(tdtom_feats, collapse = ","),
      "; in VariableFeatures: ",
      if (length(in_vf)) paste(in_vf, collapse = ",") else "none")
  }
  message("  TdTomato feature check — ", tdtom_feature_note)
  if (grepl("in VariableFeatures: (?!none)", tdtom_feature_note, perl = TRUE))
    message("  NOTE: TdTomato is a variable feature, so it may have contributed to ",
            "the clustering.\n        Cluster enrichment is then partly circular — ",
            "interpret alongside the\n        upstream clustering, not as independent evidence.")
}

# =============================================================================
# 4. The matrix: cluster x TdTomato label
# =============================================================================

message("\n=== Step 4: Building cluster x TdTomato matrix ===")

tab <- table(cluster = meta$cluster_id, label = meta$tdtom_label)
tab <- tab[, c("TdTom_pos", "TdTom_neg"), drop = FALSE]

N        <- sum(tab)
row_tot  <- rowSums(tab)
col_tot  <- colSums(tab)
exp_mat  <- outer(row_tot, col_tot) / N
pos_rate <- unname(col_tot[["TdTom_pos"]] / N)

# (1) log2 observed / expected, pseudocount 0.5 so empty cells stay finite
l2oe_mat <- log2((tab + 0.5) / (exp_mat + 0.5))

# (2) standardised (adjusted) Pearson residuals — approx N(0,1) under independence
denom    <- sqrt(exp_mat *
                 outer(1 - row_tot / N, 1 - col_tot / N))
resid_mat <- (tab - exp_mat) / denom
resid_mat[!is.finite(resid_mat)] <- 0

# (3) row-normalised composition (% of cluster)
pct_mat <- sweep(tab, 1, pmax(row_tot, 1), "/") * 100

message("  Table: ", nrow(tab), " clusters x ", ncol(tab), " labels;  N = ", N)
message("  Sample-wide TdTom+ rate: ", sprintf("%.2f%%", 100 * pos_rate))

# --- global test of independence --------------------------------------------
use_mc <- any(exp_mat < 5)
gl <- suppressWarnings(
  chisq.test(tab, simulate.p.value = use_mc,
             B = if (use_mc) 1e5 else 2000))
cramers_v <- sqrt(unname(gl$statistic) / (N * (min(dim(tab)) - 1)))
message("  Global chi-square: X2 = ", signif(unname(gl$statistic), 5),
        ", p = ", signif(gl$p.value, 4),
        if (use_mc) "  (Monte-Carlo, some expected < 5)" else "",
        ";  Cramer's V = ", signif(cramers_v, 3))

# =============================================================================
# 5. Per-cluster Fisher exact tests
# =============================================================================

message("\n=== Step 5: Per-cluster enrichment tests ===")

enr <- lapply(rownames(tab), function(cl) {
  in_pos  <- tab[cl, "TdTom_pos"]; in_neg  <- tab[cl, "TdTom_neg"]
  out_pos <- col_tot[["TdTom_pos"]] - in_pos
  out_neg <- col_tot[["TdTom_neg"]] - in_neg
  m <- matrix(c(in_pos, in_neg, out_pos, out_neg), nrow = 2,
              dimnames = list(c("pos", "neg"), c("cluster", "rest")))
  n_cl   <- in_pos + in_neg
  tested <- n_cl >= MIN_CL_CELLS

  if (tested) {
    ft2 <- fisher.test(m)                        # two-sided: different rate
    ftg <- fisher.test(m, alternative = "greater")# one-sided: enriched
    p2  <- ft2$p.value; pg <- ftg$p.value
    or  <- unname(ft2$estimate)
    lo  <- ft2$conf.int[1]; hi <- ft2$conf.int[2]
  } else {
    p2 <- NA_real_; pg <- NA_real_; or <- NA_real_; lo <- NA_real_; hi <- NA_real_
  }

  tibble(cluster      = cl,
         n_cells      = n_cl,
         n_tdtom_pos  = in_pos,
         n_tdtom_neg  = in_neg,
         pct_tdtom_pos = 100 * in_pos / max(n_cl, 1),
         exp_tdtom_pos = exp_mat[cl, "TdTom_pos"],
         log2_oe_pos  = l2oe_mat[cl, "TdTom_pos"],
         resid_pos    = resid_mat[cl, "TdTom_pos"],
         odds_ratio   = or,
         or_ci_low    = lo,
         or_ci_high   = hi,
         p_two_sided  = p2,
         p_enrich_one_sided = pg,
         tested       = tested)
}) %>% bind_rows()

# BH across tested clusters only; untested rows keep NA (see header rationale)
enr <- enr %>%
  mutate(q_two_sided = { q <- rep(NA_real_, n())
                         q[tested] <- p.adjust(p_two_sided[tested], "BH"); q },
         q_enrich    = { q <- rep(NA_real_, n())
                         q[tested] <- p.adjust(p_enrich_one_sided[tested], "BH"); q },
         call = case_when(
           !tested                                  ~ "untested",
           q_two_sided < ALPHA & log2_oe_pos > 0     ~ "enriched",
           q_two_sided < ALPHA & log2_oe_pos < 0     ~ "depleted",
           TRUE                                      ~ "ns"),
         stars = case_when(
           !tested                 ~ "",
           q_two_sided < 0.001     ~ "***",
           q_two_sided < 0.01      ~ "**",
           q_two_sided < ALPHA     ~ "*",
           TRUE                    ~ ""),
         sample = case_id) %>%
  mutate(cluster = factor(cluster, levels = levels(meta$cluster_id)))

n_tested <- sum(enr$tested)
message("  Clusters tested: ", n_tested, " / ", nrow(enr),
        " (", nrow(enr) - n_tested, " below --min_cluster_cells = ", MIN_CL_CELLS, ")")
message("  ENRICHED for TdTom_pos (q < ", ALPHA, "): ",
        { e <- enr %>% filter(call == "enriched") %>% arrange(desc(log2_oe_pos))
          if (!nrow(e)) "none" else paste0(e$cluster, " (", sprintf("%.1f%%", e$pct_tdtom_pos),
                                           ", log2 O/E ", sprintf("%+.2f", e$log2_oe_pos),
                                           e$stars, ")", collapse = "; ") })
message("  DEPLETED for TdTom_pos: ",
        { d <- enr %>% filter(call == "depleted") %>% arrange(log2_oe_pos)
          if (!nrow(d)) "none" else paste(d$cluster, collapse = ", ") })

# --- row ordering for the heatmaps ------------------------------------------
ord_levels <- switch(ORDER_BY,
  numeric    = levels(meta$cluster_id),
  enrichment = as.character(enr$cluster[order(-enr$log2_oe_pos)]),
  size       = as.character(enr$cluster[order(-enr$n_cells)]))

# y-axis labels carry n and significance stars so the matrix is self-contained
ylab_map <- setNames(
  paste0(enr$cluster, "  (n=", enr$n_cells, ")", ifelse(enr$stars == "", "", paste0(" ", enr$stars))),
  as.character(enr$cluster))

mat_to_long <- function(m, value_name) {
  as.data.frame(as.table(as.matrix(m))) %>%
    setNames(c("cluster", "label", value_name)) %>%
    mutate(cluster = as.character(cluster), label = as.character(label))
}

long <- mat_to_long(tab, "n") %>%
  left_join(mat_to_long(exp_mat,  "expected"),  by = c("cluster", "label")) %>%
  left_join(mat_to_long(l2oe_mat, "log2_oe"),   by = c("cluster", "label")) %>%
  left_join(mat_to_long(resid_mat,"std_resid"), by = c("cluster", "label")) %>%
  left_join(mat_to_long(pct_mat,  "pct_of_cluster"), by = c("cluster", "label")) %>%
  left_join(enr %>% mutate(cluster = as.character(cluster)) %>%
              dplyr::select(cluster, n_cells, q_two_sided, call, stars),
            by = "cluster") %>%
  mutate(cluster_f = factor(cluster, levels = rev(ord_levels)),
         label_f   = factor(label, levels = c("TdTom_pos", "TdTom_neg")),
         sample    = case_id)

# =============================================================================
# 6. Purity / mixture arithmetic
# =============================================================================
# The enrichment test says *where* TdTom+ cells concentrate. It cannot say
# whether an enriched cluster is PURE. That needs two extra numbers:
#
#   fp   = the TdTomato+ rate in cells that are almost certainly not tumour.
#          Taken from the strongly depleted clusters (log2 O/E < FP_L2OE_CUT),
#          i.e. the practical false-positive floor of the gate.
#   sens = the probability of detecting TdTomato in a cell that truly carries it.
#          NOT estimable from this data — so we sweep it.
#
# Then, for an observed rate p in a cluster:  p = sens*f + fp*(1-f)
#   =>  f = (p - fp) / (sens - fp)                      [clamped to 0..1]
#
# f is the implied fraction of that cluster that is truly tumour. Reading f at
# sens = 1 gives the OPTIMISTIC bound: if even a perfect gate implies f well
# below 1, the cluster is a mixture and no threshold on this label will purify
# it. That is a statement about the clustering resolution, not about the test.
# =============================================================================

message("\n=== Step 6: Purity / mixture arithmetic ===")

P_TOTAL  <- sum(enr$n_tdtom_pos)
fp_pool  <- enr %>% filter(tested, call == "depleted", log2_oe_pos < FP_L2OE_CUT)

if (!nrow(fp_pool)) {
  # No strongly depleted cluster: fall back to the single lowest-rate tested
  # cluster and say so, because the floor then rests on one population.
  fp_pool  <- enr %>% filter(tested) %>% slice_min(pct_tdtom_pos, n = 1)
  fp_basis <- paste0("fallback: lowest-rate tested cluster (", 
                     paste(fp_pool$cluster, collapse = ","), ")")
} else {
  fp_basis <- paste0("depleted clusters with log2 O/E < ", FP_L2OE_CUT, " (",
                     paste(fp_pool$cluster, collapse = ","), ")")
}
fp_rate <- sum(fp_pool$n_tdtom_pos) / sum(fp_pool$n_cells)

message("  Background TdTom+ floor: ", sprintf("%.2f%%", 100 * fp_rate),
        "  [", fp_basis, "]")

# Vectorised over BOTH arguments: the plot passes a grid of sensitivities and
# the tables pass a vector of observed rates, so a scalar `if` would be wrong.
implied_frac <- function(p_obs, sens, fp = fp_rate) {
  out <- (p_obs - fp) / (sens - fp)
  out[sens <= fp] <- NA_real_           # gate no better than its own background
  pmin(1, pmax(0, out))
}

# Per-cluster implied tumour fraction at a few reference sensitivities.
sens_ref <- c(0.60, 0.75, 0.90, 1.00)
purity_tbl <- enr %>%
  transmute(cluster, n_cells, pct_tdtom_pos, log2_oe_pos, call,
            obs_rate = pct_tdtom_pos / 100) %>%
  bind_cols(as.data.frame(setNames(
    lapply(sens_ref, function(s) round(100 * implied_frac(.$obs_rate, s), 2)),
    paste0("implied_tumour_pct_at_sens", round(100 * sens_ref))))) %>%
  mutate(fp_rate_pct = round(100 * fp_rate, 4),
         fp_basis    = fp_basis,
         sample      = case_id)

# A denser grid for the plot.
purity_curve <- expand_grid(cluster = enr$cluster,
                            sens    = seq(max(0.05, ceiling(fp_rate * 100) / 100 + 0.01),
                                          1, length.out = 80)) %>%
  left_join(enr %>% transmute(cluster, obs = pct_tdtom_pos / 100), by = "cluster") %>%
  mutate(implied_pct = 100 * implied_frac(obs, sens),
         sample      = case_id)

# --- candidate cluster rules: purity vs completeness -------------------------
# Every rule that calls a set of clusters "tumour" trades one against the other.
# Ranking enriched clusters by effect size and accumulating them shows the whole
# trade-off rather than asserting one cut.
enr_ranked <- enr %>% filter(call == "enriched") %>% arrange(desc(log2_oe_pos))
rule_sets  <- if (!nrow(enr_ranked)) list() else
  lapply(seq_len(nrow(enr_ranked)), function(k) as.character(enr_ranked$cluster[1:k]))
# plus the permissive rule that also absorbs the ambiguous (ns) clusters
ns_cl <- as.character(enr$cluster[enr$call == "ns"])
if (length(ns_cl) && nrow(enr_ranked))
  rule_sets <- c(rule_sets, list(c(as.character(enr_ranked$cluster), ns_cl)))

rules_tbl <- if (!length(rule_sets)) tibble() else
  lapply(rule_sets, function(cs) {
    m  <- as.character(enr$cluster) %in% cs
    nc <- sum(enr$n_cells[m]); tp <- sum(enr$n_tdtom_pos[m])
    tibble(rule            = paste(cs, collapse = " + "),
           n_clusters      = length(cs),
           n_cells_called  = nc,
           pct_of_section  = round(100 * nc / sum(enr$n_cells), 2),
           tdtom_pos_in    = tp,
           purity_pct      = round(100 * tp / nc, 2),
           captured_pct    = round(100 * tp / P_TOTAL, 2),
           pos_left_out    = P_TOTAL - tp,
           pos_left_out_pct= round(100 * (P_TOTAL - tp) / P_TOTAL, 2))
  }) %>% bind_rows() %>% mutate(sample = case_id)

if (nrow(rules_tbl)) {
  best <- rules_tbl %>% slice_max(captured_pct, n = 1)
  message("  Enriched-cluster rules: purity ",
          sprintf("%.1f-%.1f%%", min(rules_tbl$purity_pct), max(rules_tbl$purity_pct)),
          "; best capture ", sprintf("%.1f%%", best$captured_pct),
          " of all TdTom+ (", sprintf("%.1f%%", best$pos_left_out_pct), " left outside)")
  for (i in seq_len(nrow(enr_ranked))) {
    r_ <- enr_ranked[i, ]
    message("    cluster ", r_$cluster, ": implied tumour fraction ",
            sprintf("%.0f%%", 100 * implied_frac(r_$pct_tdtom_pos / 100, 1.0)),
            " even at perfect sensitivity",
            if (implied_frac(r_$pct_tdtom_pos / 100, 1.0) < 0.8) "  <- MIXTURE" else "")
  }
} else {
  message("  No enriched cluster: purity model has nothing to rank.")
}

# =============================================================================
# 7. Differential expression
# =============================================================================
# Two modes (--de_mode):
#   all_vs_rest (default) : every tested cluster vs all other cells = cluster
#                           markers, i.e. what each cluster IS. Read first.
#   targeted              : the contrasts below, which test the tumour question.
#   both                  : all of them.
#
# Design rationale for the TARGETED contrasts (the part that is easy to get
# wrong; all_vs_rest deliberately uses 'everything else', which is correct for
# marker discovery but wrong for 'is this tumour'):
#
#   * The reference is the DEPLETED clusters, not "everything not enriched".
#     An `ns` cluster sitting at the section-wide rate can still hold a large
#     share of all TdTom+ cells; pooling it into the reference puts tumour on
#     both sides of the contrast and flattens the genes you are looking for.
#     --de_reference nonenriched restores the naive behaviour if you want it.
#
#   * Enriched clusters are tested SEPARATELY, and against each other. Each vs
#     reference asks "is this tumour?"; enriched-vs-enriched asks "are these the
#     same tumour population, or two states?" — the contrast that decides
#     whether to merge them.
#
#   * Ambiguous (`ns`) clusters holding a large share of the TdTom+ pool are
#     tested too, because that is the direct test of whether the tumour
#     programme leaks outside the enriched clusters.
#
#   * Within-cluster TdTom+ vs TdTom- asks whether the LABEL carries signal
#     beyond cluster identity. If it returns nothing inside an enriched
#     cluster, the label is adding noise, not resolution.
#
#   * The TdTomato feature itself is flagged and excluded from the ranked plots.
#     The labels are derived from it, so it is guaranteed to top every list and
#     tells you only that the gate worked.
# =============================================================================

de_all <- tibble(); de_summary <- tibble(); de_status <- character()

if (DO_DE) {
  message("\n=== Step 7: Differential expression ===")

  ref_clusters <- switch(
    DE_REFERENCE,
    depleted    = as.character(enr$cluster[enr$tested & enr$call == "depleted"]),
    nonenriched = as.character(enr$cluster[enr$tested & enr$call %in% c("depleted", "ns")]),
    explicit    = DE_REF_EXPLICIT)
  ref_clusters <- intersect(ref_clusters, as.character(enr$cluster))
  n_ref <- sum(enr$n_cells[as.character(enr$cluster) %in% ref_clusters])

  message("  Reference (--de_reference ", DE_REFERENCE, "): ",
          if (!length(ref_clusters)) "EMPTY" else paste(ref_clusters, collapse = ", "),
          "  (", n_ref, " cells)")
  if (DE_REFERENCE == "depleted" && length(ns_cl))
    message("    note: ns cluster(s) ", paste(ns_cl, collapse = ", "),
            " held OUT of the reference (they may contain tumour)")

  # --- assay selection ------------------------------------------------------
  de_assay <- if (DE_ASSAY %in% Assays(seu)) DE_ASSAY else DefaultAssay(seu)
  if (de_assay != DE_ASSAY)
    message("    assay '", DE_ASSAY, "' absent; using '", de_assay, "'")
  if (de_assay == "SCT") {
    ok <- tryCatch({ seu <- PrepSCTFindMarkers(seu, verbose = FALSE); TRUE },
                   error = function(e) FALSE)
    if (!ok) message("    PrepSCTFindMarkers unavailable/failed; ",
                     "SCT counts used as-is (single sample, so corrected counts ",
                     "are already on one scale)")
  }

  tdtom_feats <- grep(TDTOM_FEATURE_PATTERN, rownames(seu[[de_assay]]),
                      ignore.case = TRUE, value = TRUE)
  if (length(tdtom_feats))
    message("    TdTomato feature(s) in assay: ", paste(tdtom_feats, collapse = ", "),
            " — flagged, excluded from ranked plots")

  # --- downsampling ---------------------------------------------------------
  # Wilcoxon on tens of thousands of cells is slow and its p-values are already
  # far past any decision threshold; capping cells per side keeps effect sizes
  # (avg_log2FC, pct.1/pct.2) essentially unchanged while bounding runtime.
  set.seed(SEED)
  pick <- function(cells) if (DE_MAX_CELLS > 0 && length(cells) > DE_MAX_CELLS)
    sample(cells, DE_MAX_CELLS) else cells

  cells_of_clusters <- function(cs)
    colnames(seu)[as.character(seu$cluster_id) %in% cs & !is.na(seu$tdtom_label)]
  cells_of_label <- function(cs, lab)
    colnames(seu)[as.character(seu$cluster_id) %in% cs & seu$tdtom_label == lab &
                    !is.na(seu$tdtom_label)]

  run_contrast <- function(c1, c2, contrast, kind, note = NA_character_) {
    if (length(c1) < DE_MIN_CELLS || length(c2) < DE_MIN_CELLS) {
      msg <- paste0("SKIP ", contrast, " (n=", length(c1), " vs ", length(c2),
                    "; --de_min_cells = ", DE_MIN_CELLS, ")")
      message("    ", msg); de_status <<- c(de_status, msg); return(NULL)
    }
    a <- pick(c1); b <- pick(c2)
    res <- tryCatch(
      FindMarkers(seu, ident.1 = a, ident.2 = b, assay = de_assay,
                  test.use = DE_TEST, min.pct = DE_MIN_PCT,
                  logfc.threshold = DE_LOGFC,
                  only.pos = DE_ONLY_POS && kind == "cluster_vs_rest",
                  verbose = FALSE),
      error = function(e) { msg <- paste0("FAILED ", contrast, ": ", conditionMessage(e))
                            message("    ", msg); de_status <<- c(de_status, msg); NULL })
    if (is.null(res) || !nrow(res)) return(NULL)
    res %>%
      mutate(gene = rownames(.), contrast = contrast, kind = kind,
             n_cells_1 = length(a), n_cells_2 = length(b),
             n_avail_1 = length(c1), n_avail_2 = length(c2),
             downsampled = length(a) < length(c1) || length(b) < length(c2),
             is_tdtom_feature = gene %in% tdtom_feats,
             direction = ifelse(avg_log2FC > 0, "up_in_1", "up_in_2"),
             note = note, sample = case_id) %>%
      as_tibble()
  }

  contrasts <- list()
  enr_cl <- as.character(enr_ranked$cluster)
  do_avr  <- DE_MODE %in% c("all_vs_rest", "both")
  do_tgt  <- DE_MODE %in% c("targeted", "both")

  # (0) EVERY tested cluster vs all remaining cells — the standard marker
  #     contrast. This is what identifies each cluster, and it is the contrast
  #     that answers the question the enrichment table cannot: WHAT is the
  #     ambiguous cluster made of? Note the comparator here is deliberately
  #     "everything else" — for marker discovery that is correct, because the
  #     question is what distinguishes this cluster from the rest of the
  #     section. It is NOT the right comparator for "is this tumour" (see the
  #     targeted contrasts below, where TdTom+ cells are kept out of the
  #     reference on purpose).
  if (do_avr) {
    avr_cl <- as.character(enr$cluster[enr$tested])
    # ordered so the interesting clusters run first if the job is cut short
    avr_cl <- c(intersect(enr_cl, avr_cl), setdiff(avr_cl, enr_cl))
    for (cl in avr_cl) {
      c1 <- cells_of_clusters(cl)
      c2 <- colnames(seu)[as.character(seu$cluster_id) != cl & !is.na(seu$tdtom_label) &
                            as.character(seu$cluster_id) %in% as.character(enr$cluster)]
      contrasts[[length(contrasts) + 1]] <- list(
        c1 = c1, c2 = c2,
        contrast = paste0("cl", cl, "_vs_rest"), kind = "cluster_vs_rest",
        note = paste0(enr$call[as.character(enr$cluster) == cl], " cluster; ",
                      sprintf("%.1f%%", enr$pct_tdtom_pos[as.character(enr$cluster) == cl]),
                      " TdTom+"))
    }
  }

  # (a) each enriched cluster vs the reference
  if (do_tgt && length(ref_clusters)) for (cl in enr_cl)
    contrasts[[length(contrasts) + 1]] <- list(
      c1 = cells_of_clusters(cl), c2 = cells_of_clusters(ref_clusters),
      contrast = paste0("cl", cl, "_vs_reference"), kind = "enriched_vs_reference",
      note = paste0("reference = ", paste(ref_clusters, collapse = "+")))

  # (b) enriched vs enriched — are they the same population?
  if (do_tgt && DE_PAIRWISE && length(enr_cl) > 1) {
    cmb <- combn(enr_cl, 2, simplify = FALSE)
    for (pr in cmb)
      contrasts[[length(contrasts) + 1]] <- list(
        c1 = cells_of_clusters(pr[2]), c2 = cells_of_clusters(pr[1]),
        contrast = paste0("cl", pr[2], "_vs_cl", pr[1]), kind = "enriched_vs_enriched",
        note = "few DE genes => candidates to merge")
  }

  # (c) ambiguous clusters that hold a lot of the TdTom+ pool
  amb <- enr %>% filter(call == "ns",
                        100 * n_tdtom_pos / P_TOTAL >= DE_NS_POOL_PCT) %>%
    arrange(desc(n_tdtom_pos))
  if (do_tgt && DE_AMBIGUOUS && nrow(amb) && length(ref_clusters)) for (cl in as.character(amb$cluster))
    contrasts[[length(contrasts) + 1]] <- list(
      c1 = cells_of_clusters(cl), c2 = cells_of_clusters(ref_clusters),
      contrast = paste0("cl", cl, "_vs_reference"), kind = "ambiguous_vs_reference",
      note = paste0("ns cluster holding ",
                    sprintf("%.1f%%", 100 * enr$n_tdtom_pos[as.character(enr$cluster) == cl] / P_TOTAL),
                    " of all TdTom+"))

  # (d) within-cluster TdTom+ vs TdTom- — does the label add anything?
  if (do_tgt && DE_WITHIN) for (cl in unique(c(enr_cl, as.character(amb$cluster))))
    contrasts[[length(contrasts) + 1]] <- list(
      c1 = cells_of_label(cl, "TdTom_pos"), c2 = cells_of_label(cl, "TdTom_neg"),
      contrast = paste0("cl", cl, "_TdTomPos_vs_Neg"), kind = "within_cluster_label",
      note = "tests whether the label adds signal beyond cluster identity")

  message("  Contrasts queued: ", length(contrasts))
  de_all <- lapply(contrasts, function(x)
    run_contrast(x$c1, x$c2, x$contrast, x$kind, x$note)) %>% bind_rows()

  if (nrow(de_all)) {
    de_summary <- de_all %>%
      group_by(contrast, kind, n_cells_1, n_cells_2, note) %>%
      summarise(n_tested_genes = n(),
                n_sig          = sum(p_val_adj < ALPHA, na.rm = TRUE),
                n_sig_up_in_1  = sum(p_val_adj < ALPHA & avg_log2FC > 0, na.rm = TRUE),
                n_sig_up_in_2  = sum(p_val_adj < ALPHA & avg_log2FC < 0, na.rm = TRUE),
                n_sig_excl_tdtom = sum(p_val_adj < ALPHA & !is_tdtom_feature, na.rm = TRUE),
                max_abs_log2FC = max(abs(avg_log2FC), na.rm = TRUE),
                top_up_in_1    = paste(head(gene[order(-avg_log2FC)][
                                        !is_tdtom_feature[order(-avg_log2FC)]], 8), collapse = ";"),
                top_up_in_2    = paste(head(gene[order(avg_log2FC)][
                                        !is_tdtom_feature[order(avg_log2FC)]], 8), collapse = ";"),
                .groups = "drop") %>%
      mutate(sample = case_id) %>%
      arrange(match(kind, c("cluster_vs_rest", "enriched_vs_reference",
                            "ambiguous_vs_reference", "enriched_vs_enriched",
                            "within_cluster_label")),
              desc(n_sig))

    # --- marker table: top genes per cluster from the all-vs-rest contrasts ---
    # This is the table to read first — it says what each cluster IS, which is
    # the context every enrichment call has to be interpreted in.
    if (any(de_all$kind == "cluster_vs_rest")) {
      markers_tbl <- de_all %>%
        filter(kind == "cluster_vs_rest", p_val_adj < ALPHA,
               avg_log2FC > 0, !is_tdtom_feature) %>%
        mutate(cluster = sub("^cl(.*)_vs_rest$", "\\1", contrast)) %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        mutate(marker_rank = row_number()) %>%
        ungroup() %>%
        left_join(enr %>% transmute(cluster = as.character(cluster),
                                    pct_tdtom_pos, log2_oe_pos, call),
                  by = "cluster") %>%
        dplyr::select(sample, cluster, call, pct_tdtom_pos, log2_oe_pos,
                      marker_rank, gene, avg_log2FC, pct.1, pct.2,
                      p_val, p_val_adj)
      write.csv(markers_tbl,
                file.path(OUT_DIR, paste0(case_id, "_cluster_markers.csv")),
                row.names = FALSE)

      message("  Top markers per cluster (up, q < ", ALPHA, ", TdTomato excluded):")
      for (cl in as.character(enr$cluster)) {
        g <- markers_tbl %>% filter(cluster == cl) %>%
          slice_head(n = DE_MARKER_TOP_N)
        if (!nrow(g)) next
        message("    cl ", format(cl, width = 4),
                format(enr$call[as.character(enr$cluster) == cl], width = 9),
                sprintf("%5.1f%% TdTom+  ", enr$pct_tdtom_pos[as.character(enr$cluster) == cl]),
                paste(g$gene, collapse = ", "))
      }
    }

    write.csv(de_all, file.path(OUT_DIR, paste0(case_id, "_DE_all_contrasts.csv")),
              row.names = FALSE)
    write.csv(de_summary, file.path(OUT_DIR, paste0(case_id, "_DE_summary.csv")),
              row.names = FALSE)
    if (DE_SPLIT_FILES) for (ct in unique(de_all$contrast))
      write.csv(de_all %>% filter(contrast == ct),
                file.path(OUT_DIR, paste0(case_id, "_DE_", ct, ".csv")), row.names = FALSE)

    for (i in seq_len(nrow(de_summary))) {
      s_ <- de_summary[i, ]
      message("    ", format(s_$contrast, width = 26), " ",
              s_$n_sig_excl_tdtom, " sig genes (q<", ALPHA, ", excl TdTomato)",
              "  up1=", s_$n_sig_up_in_1, " up2=", s_$n_sig_up_in_2)
    }

    # The merge question, stated explicitly rather than left to the reader.
    pw <- de_summary %>% filter(kind == "enriched_vs_enriched")
    if (nrow(pw)) for (i in seq_len(nrow(pw)))
      message("  ", pw$contrast[i], ": ", pw$n_sig_excl_tdtom[i],
              " DE genes between enriched clusters — ",
              if (pw$n_sig_excl_tdtom[i] < DE_MERGE_HINT) 
                paste0("below ", DE_MERGE_HINT, ", consider MERGING")
              else "distinct transcriptional states, keep separate")

    # Shared vs unique programme across the vs-reference contrasts.
    vr <- de_all %>% filter(kind %in% c("enriched_vs_reference", "ambiguous_vs_reference"),
                            p_val_adj < ALPHA, !is_tdtom_feature, avg_log2FC > 0)
    if (n_distinct(vr$contrast) > 1) {
      ov <- vr %>% count(gene, name = "n_contrasts")
      overlap_tbl <- vr %>% left_join(ov, by = "gene") %>%
        group_by(contrast) %>%
        summarise(n_up = n(),
                  n_shared_with_any = sum(n_contrasts > 1),
                  pct_shared = round(100 * sum(n_contrasts > 1) / n(), 1),
                  .groups = "drop") %>% mutate(sample = case_id)
      write.csv(overlap_tbl, file.path(OUT_DIR, paste0(case_id, "_DE_overlap.csv")),
                row.names = FALSE)
      message("  Up-gene sharing across vs-reference contrasts:")
      for (i in seq_len(nrow(overlap_tbl)))
        message("    ", format(overlap_tbl$contrast[i], width = 26),
                overlap_tbl$n_up[i], " up, ", overlap_tbl$pct_shared[i], "% shared")
    }
  } else {
    message("  No DE results (all contrasts skipped or empty).")
  }
} else {
  message("\n=== Step 7: Differential expression — SKIPPED (--de FALSE) ===")
}

# =============================================================================
# 8. Sub-cluster diagnostic  (optional; tests the mixture explanation)
# =============================================================================
# If an enriched cluster is a mixture (step 6), one of two things is true:
#   (i)  the resolution is too coarse — tumour cells sit inside it as a
#        sub-population, and re-clustering the cluster alone will expose a
#        sub-cluster with a much higher TdTom+ rate; or
#   (ii) the label is lossy — TdTom+ cells are spread evenly through the
#        cluster and no sub-structure tracks them.
# Re-clustering each target in isolation distinguishes these. A high
# heterogeneity statistic with one high-rate sub-cluster supports (i) and means
# sub-clustering will do more for you than any cluster-level tumour call.
# =============================================================================

sub_tbl <- tibble()
# Per-cell sub-cluster assignments, kept so Step 8b can set Idents() to them.
# Step 8 alone only needs the aggregated rates, but discarding the assignments
# would mean re-running the whole sub-clustering to get markers out of it.
sub_cells <- tibble()

if (DO_SUBCLUSTER) {
  message("\n=== Step 8: Sub-cluster diagnostic ===")
  targets <- if (length(SUBCLUSTER_CLUSTERS)) SUBCLUSTER_CLUSTERS else
    unique(c(as.character(enr_ranked$cluster),
             as.character((enr %>% filter(call == "ns",
                                          100 * n_tdtom_pos / P_TOTAL >= DE_NS_POOL_PCT))$cluster)))
  targets <- intersect(targets, as.character(enr$cluster))
  message("  Targets: ", if (!length(targets)) "none" else paste(targets, collapse = ", "))

  for (cl in targets) {
    cells <- colnames(seu)[as.character(seu$cluster_id) == cl & !is.na(seu$tdtom_label)]
    if (length(cells) < SUBCLUSTER_MIN_CELLS) {
      message("    cluster ", cl, ": SKIP (", length(cells), " cells < ",
              SUBCLUSTER_MIN_CELLS, ")"); next
    }
    # --subcluster_max_cells downsamples before re-clustering. When it bites, every
    # sub-cluster n_cells for this parent is a SAMPLE of the real population, not
    # its size — so the fact is recorded in the output table (parent_n_labelled /
    # parent_downsampled), not only in this log line. Rates are unaffected;
    # absolute counts are.
    set.seed(SEED)
    n_labelled <- length(cells)
    was_capped <- SUBCLUSTER_MAX_CELLS > 0 && n_labelled > SUBCLUSTER_MAX_CELLS
    if (was_capped) {
      cells <- sample(cells, SUBCLUSTER_MAX_CELLS)
      message("    cluster ", cl, ": DOWNSAMPLED ", n_labelled, " -> ",
              SUBCLUSTER_MAX_CELLS, " cells for sub-clustering",
              " (sub-cluster n_cells are a sample of the parent;",
              " use --subcluster_max_cells 0 to sub-cluster all cells)")
    }

    res <- tryCatch({
      sub <- subset(seu, cells = cells)
      DefaultAssay(sub) <- if (DE_ASSAY %in% Assays(sub)) DE_ASSAY else DefaultAssay(sub)
      if (DefaultAssay(sub) != "SCT") {
        sub <- NormalizeData(sub, verbose = FALSE)
        sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
        sub <- ScaleData(sub, verbose = FALSE)
      } else if (!length(VariableFeatures(sub))) {
        sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
      }
      # Drop the TdTomato feature from the sub-clustering so the sub-structure
      # is not defined by the very label we are about to test it against.
      vf <- setdiff(VariableFeatures(sub),
                    grep(TDTOM_FEATURE_PATTERN, rownames(sub), ignore.case = TRUE, value = TRUE))
      VariableFeatures(sub) <- vf
      npc <- min(SUBCLUSTER_DIMS, length(cells) - 1, length(vf) - 1)
      sub <- RunPCA(sub, features = vf, npcs = npc, verbose = FALSE)
      sub <- FindNeighbors(sub, dims = 1:npc, verbose = FALSE)
      sub <- FindClusters(sub, resolution = SUBCLUSTER_RES, verbose = FALSE)
      tibble(cell = colnames(sub), sub_id = as.character(Idents(sub)),
             lab = sub$tdtom_label)
    }, error = function(e) { message("    cluster ", cl, ": FAILED (",
                                     conditionMessage(e), ")"); NULL })
    if (is.null(res)) next

    # Namespace the sub-cluster id by its parent: FindClusters restarts numbering
    # from 0 in every target, so a bare "1" is ambiguous across parents.
    #
    # The separator is "_" and NOT ".". With a dot, "7.1" and "7.10" are the same
    # number, so any reader that lets the column be type-inferred (read.csv,
    # pandas.read_csv, Excel) silently merges those two sub-clusters into one row
    # and loses the other. Underscore cannot be parsed as numeric, so the id
    # survives type inference in every reader.
    sub_cells <- bind_rows(sub_cells,
                           res %>% transmute(cell, parent_cluster = cl,
                                             sub_id = paste0(cl, "_", sub_id)))

    parent_rate <- mean(res$lab == "TdTom_pos")
    st <- res %>% group_by(sub_id) %>%
      summarise(n_cells = n(),
                n_tdtom_pos = sum(lab == "TdTom_pos"),
                .groups = "drop") %>%
      mutate(pct_tdtom_pos = round(100 * n_tdtom_pos / n_cells, 3),
             parent_cluster = cl,
             # provenance of the counts: n_cells sums to parent_n_subclustered,
             # which is < parent_n_labelled exactly when parent_downsampled is TRUE
             parent_n_labelled = n_labelled,
             parent_n_subclustered = length(cells),
             parent_downsampled = was_capped,
             parent_pct_tdtom_pos = round(100 * parent_rate, 3),
             log2_oe_vs_parent = round(log2((n_tdtom_pos + 0.5) /
                                            (n_cells * parent_rate + 0.5)), 4),
             implied_tumour_pct_at_sens100 = round(100 * implied_frac(pct_tdtom_pos / 100, 1), 2),
             sample = case_id) %>%
      arrange(desc(pct_tdtom_pos))

    # Heterogeneity is only defined across >= 2 sub-clusters. On a 1-row table
    # chisq.test silently switches to a goodness-of-fit test against 50/50,
    # which is not the question — so it is not run at all in that case.
    ht <- if (nrow(st) < 2) NULL else
      tryCatch(suppressWarnings(
        chisq.test(cbind(st$n_tdtom_pos, st$n_cells - st$n_tdtom_pos))),
        error = function(e) NULL)
    st$het_chisq <- if (is.null(ht)) NA_real_ else signif(unname(ht$statistic), 6)
    st$het_p     <- if (is.null(ht)) NA_real_ else signif(ht$p.value, 6)
    st$het_note  <- if (nrow(st) < 2)
      "single sub-cluster: no sub-structure at this resolution" else NA_character_
    sub_tbl <- bind_rows(sub_tbl, st)

    top <- st[1, ]
    message("    cluster ", cl, ": ", nrow(st), " sub-cluster",
            if (nrow(st) == 1) "" else "s", ", parent ",
            sprintf("%.1f%%", 100 * parent_rate), " TdTom+",
            if (nrow(st) < 2) "" else paste0(
              "; best sub-cluster ", top$sub_id, " at ",
              sprintf("%.1f%%", top$pct_tdtom_pos), " (n=", top$n_cells,
              "); heterogeneity p = ",
              if (is.na(st$het_p[1])) "NA" else format(st$het_p[1], digits = 3)))
    message("      -> ", if (nrow(st) < 2)
              paste0("no sub-structure found at resolution ", SUBCLUSTER_RES,
                     " — raise --subcluster_res to probe further")
            else if (!is.na(st$het_p[1]) && st$het_p[1] < ALPHA &&
                     top$pct_tdtom_pos >= 100 * parent_rate * SUBCLUSTER_GAIN)
              paste0("sub-structure tracks TdTomato: RESOLUTION is limiting; ",
                     "sub-cluster ", cl, " before calling tumour")
            else "no sub-cluster tracks TdTomato: label is likely lossy, not resolution")
  }
  if (nrow(sub_tbl))
    write.csv(sub_tbl, file.path(OUT_DIR, paste0(case_id, "_subcluster_tdtom.csv")),
              row.names = FALSE)
} else {
  message("\n=== Step 8: Sub-cluster diagnostic — SKIPPED (--subcluster FALSE) ===")
}

# =============================================================================
# 8b. Sub-cluster marker DE — Idents() set to sub-cluster identity
#
# Step 8 ranks sub-clusters on reporter rate alone, which is the weaker axis:
# it cannot say what a sub-cluster IS. This step runs FindMarkers with the
# per-cell sub-cluster assignments as identities, so each sub-cluster gets a
# marker table of the same shape as the cluster-level one and can be scored
# against cell-type / programme panels on transcriptional grounds.
#
# Two comparators, because they answer different questions:
#   section — sub vs every other labelled cell in the section. The analogue of
#             cluster-vs-rest, and the table to score against marker panels.
#   parent  — sub vs the remaining cells of its own parent cluster. Says what
#             distinguishes this sub-cluster from its siblings; genes shared
#             with the parent cancel out, so this table is NOT suitable for
#             identity scoring — a tumour sub-cluster inside a tumour cluster
#             shows no tumour markers here.
# =============================================================================

sub_markers <- tibble(); sub_de_status <- character()

if (DO_SUBCLUSTER && SUBCLUSTER_DE && nrow(sub_cells)) {
  message("\n=== Step 8b: Sub-cluster marker DE (Idents <- sub_id) ===")

  # Step 7 defines de_assay / tdtom_feats / pick() / DE_* inside its own
  # conditional; when --de FALSE they were never created, so rebuild them here.
  # Cheap, and it keeps Step 8b runnable without a full DE pass.
  if (!exists("de_assay")) {
    de_assay <- if (DE_ASSAY %in% Assays(seu)) DE_ASSAY else DefaultAssay(seu)
    message("  DE assay: ", de_assay, " (Step 7 did not run)")
  }
  if (!exists("tdtom_feats_de"))
    tdtom_feats_de <- grep(TDTOM_FEATURE_PATTERN, rownames(seu[[de_assay]]),
                           ignore.case = TRUE, value = TRUE)
  sub_pick <- function(cells) {
    if (DE_MAX_CELLS > 0 && length(cells) > DE_MAX_CELLS)
      sample(cells, DE_MAX_CELLS) else cells
  }

  # Assignments are namespaced (parent.sub) and cover only the sub-clustered
  # targets — and only the cells that survived --subcluster_max_cells. Cells
  # outside that set are the comparator pool for mode 'section'.
  sub_cells <- sub_cells %>% filter(cell %in% colnames(seu))
  sub_map    <- setNames(sub_cells$sub_id, sub_cells$cell)
  parent_map <- setNames(sub_cells$parent_cluster, sub_cells$cell)
  labelled <- colnames(seu)[!is.na(seu$tdtom_label)]
  message("  Sub-clusters with assignments: ", dplyr::n_distinct(sub_map),
          " across ", dplyr::n_distinct(sub_cells$parent_cluster), " parents (",
          length(sub_map), " cells)")

  run_sub_contrast <- function(c1, c2, sid, parent, mode) {
    if (length(c1) < SUBCLUSTER_DE_MIN || length(c2) < SUBCLUSTER_DE_MIN) {
      msg <- paste0("SKIP sub ", sid, " (", mode, "): n=", length(c1), " vs ",
                    length(c2), "; --subcluster_de_min_cells = ", SUBCLUSTER_DE_MIN)
      message("    ", msg); sub_de_status <<- c(sub_de_status, msg); return(NULL)
    }
    set.seed(SEED)
    a <- sub_pick(c1); b <- sub_pick(c2)
    res <- tryCatch(
      FindMarkers(seu, ident.1 = a, ident.2 = b, assay = de_assay,
                  test.use = DE_TEST, min.pct = DE_MIN_PCT,
                  logfc.threshold = DE_LOGFC, only.pos = FALSE, verbose = FALSE),
      error = function(e) {
        msg <- paste0("FAILED sub ", sid, " (", mode, "): ", conditionMessage(e))
        message("    ", msg); sub_de_status <<- c(sub_de_status, msg); NULL })
    if (is.null(res) || !nrow(res)) return(NULL)
    res %>%
      mutate(gene = rownames(.), sub_id_full = sid, parent_cluster = parent,
             comparator = mode,
             n_cells_1 = length(a), n_cells_2 = length(b),
             n_avail_1 = length(c1), n_avail_2 = length(c2),
             downsampled = length(a) < length(c1) || length(b) < length(c2),
             is_tdtom_feature = gene %in% tdtom_feats_de) %>%
      as_tibble()
  }

  modes <- if (SUBCLUSTER_DE_MODE == "both") c("section", "parent")
           else SUBCLUSTER_DE_MODE
  out <- list()
  for (sid in sort(unique(sub_map))) {
    in_sub <- names(sub_map)[sub_map == sid]
    parent <- unname(parent_map[in_sub[1]])
    for (mode in modes) {
      comp <- if (mode == "section") setdiff(labelled, in_sub)
              else setdiff(names(parent_map)[parent_map == parent], in_sub)
      out[[paste(sid, mode)]] <- run_sub_contrast(in_sub, comp, sid, parent, mode)
    }
  }
  sub_markers <- bind_rows(out)

  if (nrow(sub_markers)) {
    # Recover the bare sub_id for the join to sub_tbl (which keys on the bare id
    # within parent). Done by stripping exactly nchar(parent) + 1 leading
    # characters rather than by regex: a regex on the separator would break the
    # moment a cluster label itself contained one, and this cannot.
    stopifnot(all(substring(sub_markers$sub_id_full, 1,
                            nchar(sub_markers$parent_cluster) + 1) ==
                    paste0(sub_markers$parent_cluster, "_")))
    sub_markers <- sub_markers %>%
      mutate(sample = case_id,
             sub_id = substring(sub_id_full, nchar(parent_cluster) + 2)) %>%
      left_join(sub_tbl %>%
                  transmute(parent_cluster = as.character(parent_cluster),
                            sub_id = as.character(sub_id),
                            sub_n_cells = n_cells,
                            sub_pct_tdtom_pos = pct_tdtom_pos,
                            parent_pct_tdtom_pos, log2_oe_vs_parent,
                            parent_n_labelled, parent_downsampled),
                by = c("parent_cluster", "sub_id")) %>%
      # Rank up-genes within each sub-cluster x comparator, TdTomato excluded,
      # so the table can be read the same way as _cluster_markers.csv.
      group_by(comparator, sub_id_full) %>%
      arrange(desc(avg_log2FC), .by_group = TRUE) %>%
      mutate(marker_rank = ifelse(p_val_adj < ALPHA & avg_log2FC > 0 &
                                    !is_tdtom_feature,
                                  cumsum(p_val_adj < ALPHA & avg_log2FC > 0 &
                                           !is_tdtom_feature),
                                  NA_integer_)) %>%
      ungroup() %>%
      dplyr::select(sample, comparator, parent_cluster, sub_id, sub_id_full,
                    sub_n_cells, sub_pct_tdtom_pos, parent_pct_tdtom_pos,
                    parent_n_labelled, parent_downsampled,
                    log2_oe_vs_parent, marker_rank, gene, avg_log2FC,
                    pct.1, pct.2, p_val, p_val_adj, n_cells_1, n_cells_2,
                    n_avail_1, n_avail_2, downsampled, is_tdtom_feature) %>%
      arrange(comparator, parent_cluster, sub_id, marker_rank)

    write.csv(sub_markers,
              file.path(OUT_DIR, paste0(case_id, "_subcluster_markers.csv")),
              row.names = FALSE)

    for (mode in modes) {
      message("  Top markers per sub-cluster [", mode, " comparator] (up, q < ",
              ALPHA, ", TdTomato excluded):")
      md <- sub_markers %>% filter(comparator == mode, !is.na(marker_rank))
      for (sid in sort(unique(md$sub_id_full))) {
        g <- md %>% filter(sub_id_full == sid) %>% slice_head(n = DE_MARKER_TOP_N)
        if (!nrow(g)) next
        message("    sub ", format(sid, width = 7),
                " (n=", format(g$sub_n_cells[1], width = 5), ", ",
                sprintf("%5.1f%%", g$sub_pct_tdtom_pos[1]), " TdTom+): ",
                paste(g$gene, collapse = ", "))
      }
    }
    n_scored <- sub_markers %>% filter(comparator == modes[1],
                                       !is.na(marker_rank)) %>%
      dplyr::pull(sub_id_full) %>% dplyr::n_distinct()
    message("  Sub-clusters with >=1 significant up-gene: ", n_scored,
            " of ", dplyr::n_distinct(sub_map))
  } else {
    message("  No sub-cluster markers returned — nothing written.")
  }
} else if (DO_SUBCLUSTER && SUBCLUSTER_DE) {
  message("\n=== Step 8b: Sub-cluster marker DE — SKIPPED ",
          "(no sub-cluster assignments from Step 8) ===")
} else {
  message("\n=== Step 8b: Sub-cluster marker DE — SKIPPED (--subcluster_de FALSE) ===")
}

# =============================================================================
# 9. Plots — heatmap matrices first, then the UMAP / spatial context
# =============================================================================

message("\n=== Step 9: Generating plots ===")
pdf_path <- file.path(OUT_DIR, paste0(case_id, "_classification_plots.pdf"))
# The base pdf() device mangles UTF-8 (em dashes come out as "..."), so use
# the cairo device when the R build has it and fall back otherwise.
if (isTRUE(capabilities("cairo"))) {
  cairo_pdf(pdf_path, width = 14, height = 8, onefile = TRUE)
} else {
  pdf(pdf_path, width = 14, height = 8)
}

scale_x_log2 <- function(...) {
  if ("transform" %in% names(formals(ggplot2::scale_x_continuous)))
    ggplot2::scale_x_continuous(transform = "log2", ...)
  else
    ggplot2::scale_x_continuous(trans = "log2", ...)
}

heat_base <- function(p) {
  p +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(labels = function(x) ylab_map[x]) +
    coord_equal(ratio = 0.6) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(family = "mono"),
          plot.title.position = "plot")
}

lim_oe <- max(abs(long$log2_oe), na.rm = TRUE)
p_h_oe <- heat_base(ggplot(long, aes(label_f, cluster_f, fill = log2_oe))) +
  geom_text(aes(label = sprintf("%+.2f", log2_oe),
                colour = abs(log2_oe) > 0.6 * lim_oe), size = 2.9, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_oe, lim_oe),
                       name = "log2 O/E") +
  labs(title = "log2 observed / expected",
       subtitle = "effect size; expected = row x col / N",
       x = NULL, y = NULL)

lim_r <- max(abs(long$std_resid), na.rm = TRUE)
p_h_res <- heat_base(ggplot(long, aes(label_f, cluster_f, fill = std_resid))) +
  geom_text(aes(label = sprintf("%+.1f", std_resid),
                colour = abs(std_resid) > 0.6 * lim_r), size = 2.9, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_r, lim_r),
                       name = "std resid") +
  labs(title = "Standardised residuals",
       subtitle = "~N(0,1) under independence; |r| > 2",
       x = NULL, y = NULL) +
  theme(axis.text.y = element_blank())

p_h_pct <- heat_base(ggplot(long, aes(label_f, cluster_f, fill = pct_of_cluster))) +
  geom_text(aes(label = paste0(n, "\n", sprintf("%.1f%%", pct_of_cluster)),
                colour = pct_of_cluster > 60), size = 2.7,
            lineheight = 0.95, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_fill_gradient(low = "white", high = "#4a1486", name = "% of cluster",
                      limits = c(0, 100)) +
  labs(title = "Counts and composition",
       subtitle = "n cells / % of that cluster",
       x = NULL, y = NULL) +
  theme(axis.text.y = element_blank())

# page 1: the matrix, three encodings side by side. The shared caption lives on
# the patchwork rather than on the first panel, otherwise a long subtitle runs
# underneath its neighbours' titles.
print((p_h_oe | p_h_res | p_h_pct) +
        plot_annotation(
          title = paste0(case_id, " - cluster x TdTomato enrichment matrix"),
          subtitle = paste0(
            "Sample-wide TdTom+ rate ", sprintf("%.1f%%", 100 * pos_rate),
            " (", col_tot[["TdTom_pos"]], "/", N, " cells).  ",
            "Stars = BH q from the two-sided Fisher test of cluster vs rest: ",
            "* q<", ALPHA, ", ** q<0.01, *** q<0.001.  ",
            "Global chi-square p = ", signif(gl$p.value, 3),
            ", Cramer's V = ", signif(cramers_v, 3), "."),
          theme = theme(plot.title = element_text(face = "bold", size = 13),
                        plot.subtitle = element_text(size = 9, colour = "grey25"))))

# page 2: ranked effect sizes with exact CIs + the composition bars
enr_plot <- enr %>%
  mutate(cluster_f = factor(as.character(cluster),
                            levels = as.character(cluster)[order(log2_oe_pos)]))

p_rank <- ggplot(enr_plot, aes(log2_oe_pos, cluster_f, colour = call)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_segment(aes(x = 0, xend = log2_oe_pos, yend = cluster_f), linewidth = 0.6) +
  geom_point(aes(size = n_cells)) +
  # Labels sit off the end of each lollipop — nudged in the direction the bar
  # points, so they never land on top of the marker.
  geom_text(aes(label = paste0(sprintf("%.1f%%", pct_tdtom_pos), stars),
                hjust = ifelse(log2_oe_pos >= 0, -0.35, 1.35)),
            size = 2.8, show.legend = FALSE) +
  scale_colour_manual(values = CALL_COLOURS, name = NULL) +
  scale_size_area(max_size = 6, name = "cells") +
  scale_x_continuous(expand = expansion(mult = c(0.22, 0.22))) +
  labs(title = paste0(case_id, " - TdTom+ enrichment per cluster"),
       subtitle = paste0("label = % TdTom+ in cluster; dashed line = sample rate ",
                         sprintf("%.1f%%", 100 * pos_rate), "\n",
                         nrow(enr) - n_tested, " cluster(s) untested (n < ",
                         MIN_CL_CELLS, ")"),
       x = "log2 observed/expected (TdTom_pos)", y = "cluster") +
  theme_bw(base_size = 10)

# A cluster that is 100% TdTom+ (or 0%) has an infinite conditional odds ratio
# and a half-open exact CI. Dropping those rows would silently remove the most
# extreme clusters — exactly the ones being looked for — so they are clamped to
# the panel limits and drawn as open triangles with the bound left unclosed.
or_plot <- enr_plot %>% filter(tested)
or_finite <- c(or_plot$odds_ratio, or_plot$or_ci_low, or_plot$or_ci_high)
or_finite <- or_finite[is.finite(or_finite) & or_finite > 0]
or_lo_lim <- if (length(or_finite)) min(or_finite) / 2 else 0.5
or_hi_lim <- if (length(or_finite)) max(or_finite) * 2 else 2
or_plot <- or_plot %>%
  mutate(or_clamped  = !is.finite(odds_ratio) | odds_ratio <= 0,
         or_plot_val = pmin(pmax(ifelse(or_clamped,
                                        ifelse(log2_oe_pos > 0, or_hi_lim, or_lo_lim),
                                        odds_ratio), or_lo_lim), or_hi_lim),
         ci_lo_plot  = pmax(ifelse(is.finite(or_ci_low)  & or_ci_low  > 0,
                                   or_ci_low,  or_lo_lim), or_lo_lim),
         ci_hi_plot  = pmin(ifelse(is.finite(or_ci_high) & or_ci_high > 0,
                                   or_ci_high, or_hi_lim), or_hi_lim))

p_or <- ggplot(or_plot, aes(or_plot_val, cluster_f, colour = call)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  # geom_errorbar + orientation="y": geom_errorbarh is deprecated in ggplot2 4.0
  geom_errorbar(aes(xmin = ci_lo_plot, xmax = ci_hi_plot),
                orientation = "y", width = 0.25, linewidth = 0.5) +
  geom_point(aes(shape = or_clamped), size = 2, fill = NA) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 24), guide = "none") +
  scale_colour_manual(values = CALL_COLOURS, name = NULL, guide = "none") +
  # `trans` was renamed to `transform` in ggplot2 3.5; pick the one this
  # installation accepts so the script runs on either side of that change.
  scale_x_log2() +
  labs(title = "Fisher conditional odds ratio (cluster vs rest)",
       subtitle = paste0("exact 95% CI; log2 x-axis",
                         if (any(or_plot$or_clamped))
                           "\ntriangle = infinite OR (cluster is 100% / 0% TdTom+), clamped to axis"
                         else ""),
       x = "odds ratio, TdTom+ in cluster vs rest", y = NULL) +
  theme_bw(base_size = 10)

p_stack <- ggplot(long, aes(y = cluster_f, x = n, fill = label_f)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = TDTOM_COLOURS, name = NULL) +
  scale_y_discrete(labels = function(x) ylab_map[x]) +
  labs(title = "Cluster composition (absolute)",
       subtitle = "cluster size is the power available to each test",
       x = "Cells", y = NULL) +
  theme_bw(base_size = 10) + theme(axis.text.y = element_text(family = "mono"))

print(p_rank | (p_or / p_stack))

# --- page: purity / mixture --------------------------------------------------
# The point of this page is that enrichment and purity are different questions.
# Left: implied tumour fraction against assumed gate sensitivity, one line per
# cluster — a line that stays low even at sens = 100% is a mixture. Right: the
# purity/completeness trade-off over the candidate cluster rules.

pc <- purity_curve %>%
  left_join(enr %>% transmute(cluster, call, n_cells), by = "cluster") %>%
  filter(call %in% c("enriched", "ns") | n_cells >= MIN_CL_CELLS)

curve_ends <- pc %>% group_by(cluster) %>% slice_max(sens, n = 1) %>% ungroup() %>%
  filter(implied_pct > 1)

p_purity <- ggplot(pc, aes(100 * sens, implied_pct,
                           group = cluster,
                           colour = factor(cluster, levels = cl_levels))) +
  geom_hline(yintercept = 100, linetype = "dotted", colour = "grey45") +
  geom_line(aes(linewidth = call == "enriched"), alpha = 0.9, show.legend = FALSE) +
  scale_linewidth_manual(values = c(`TRUE` = 1.1, `FALSE` = 0.4)) +
  ggrepel::geom_text_repel(data = curve_ends,
                           aes(label = paste0("cl ", cluster)),
                           size = 2.7, direction = "y", hjust = 0,
                           nudge_x = 2, segment.size = 0.2,
                           min.segment.length = 0, box.padding = 0.15) +
  scale_colour_manual(values = cluster_colours, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  coord_cartesian(ylim = c(0, 105)) +
  labs(title = "Implied tumour fraction vs assumed TdTomato sensitivity",
       subtitle = paste0("f = (obs - fp)/(sens - fp), fp = ",
                         sprintf("%.2f%%", 100 * fp_rate),
                         " from ", fp_basis,
                         "\nthick = enriched; a line still low at sens=100% is a MIXTURE, ",
                         "not a tumour cluster"),
       x = "assumed TdTomato detection sensitivity (%)",
       y = "implied true tumour fraction (%)") +
  theme_bw(base_size = 10)

if (nrow(rules_tbl)) {
  rl <- rules_tbl %>%
    mutate(rule_f = factor(rule, levels = rule)) %>%
    dplyr::select(rule_f, purity_pct, captured_pct) %>%
    pivot_longer(-rule_f, names_to = "metric", values_to = "pct") %>%
    mutate(metric = recode(metric,
                           purity_pct   = "TdTom+ within the call (purity)",
                           captured_pct = "of all TdTom+ captured (completeness)"))
  p_rules <- ggplot(rl, aes(rule_f, pct, fill = metric)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_hline(aes(yintercept = 100 * pos_rate, linetype = "section-wide TdTom+ rate"),
               colour = "grey35") +
    geom_text(aes(label = sprintf("%.0f", pct)),
              position = position_dodge(width = 0.75),
              vjust = -0.4, size = 2.6) +
    scale_fill_manual(values = c("TdTom+ within the call (purity)" = "#1f4fd8",
                                 "of all TdTom+ captured (completeness)" = "#f0a202"),
                      name = NULL) +
    scale_linetype_manual(values = c("section-wide TdTom+ rate" = "dashed"), name = NULL) +
    coord_cartesian(ylim = c(0, 108)) +
    labs(title = "Purity vs completeness across candidate tumour rules",
         subtitle = paste0("clusters added in order of effect size; ",
                           "no rule is usually both pure and complete"),
         x = "clusters called tumour", y = "per cent") +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", legend.box = "vertical",
          legend.margin = margin(b = -4),
          axis.text.x = element_text(angle = 20, hjust = 1))
  print(p_purity | p_rules)
} else {
  print(p_purity)
}

# --- page: sub-cluster diagnostic --------------------------------------------
if (nrow(sub_tbl)) {
  sub_plot <- sub_tbl %>%
    mutate(lab = paste0("cl ", parent_cluster, " / sub ", sub_id),
           lab = factor(lab, levels = rev(lab[order(parent_cluster, -pct_tdtom_pos)])))
  p_sub <- ggplot(sub_plot, aes(pct_tdtom_pos, lab)) +
    geom_segment(aes(x = parent_pct_tdtom_pos, xend = pct_tdtom_pos,
                     y = lab, yend = lab), colour = "grey75", linewidth = 0.4,
                 inherit.aes = FALSE) +
    geom_point(aes(colour = parent_cluster, size = n_cells)) +
    geom_vline(aes(xintercept = 100 * pos_rate), linetype = "dashed", colour = "grey35") +
    geom_text(aes(label = sprintf("%.1f%%", pct_tdtom_pos)),
              hjust = 0, nudge_x = 2.2, size = 2.5, show.legend = FALSE) +
    scale_colour_manual(values = cluster_colours, name = "parent") +
    scale_size_area(max_size = 6, name = "cells") +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    facet_grid(parent_cluster ~ ., scales = "free_y", space = "free_y") +
    labs(title = "Sub-cluster TdTomato rates within each target cluster",
         subtitle = paste0("grey link = parent cluster rate; dashed = section-wide. ",
                           "One sub-cluster far above its parent means\nresolution ",
                           "is limiting: sub-cluster before calling tumour."),
         x = "% TdTom+ within sub-cluster", y = NULL) +
    theme_bw(base_size = 10) + theme(strip.text.y = element_text(angle = 0))
  print(p_sub)
}

# --- pages: DE ---------------------------------------------------------------
if (nrow(de_all)) {

  # Overview: how many significant genes each contrast returns, by kind. This is
  # the page that answers "is this cluster tumour" and "should these merge".
  ds <- de_summary %>%
    mutate(contrast_f = factor(contrast, levels = rev(contrast))) %>%
    dplyr::select(contrast_f, kind, n_sig_up_in_1, n_sig_up_in_2) %>%
    pivot_longer(starts_with("n_sig_up"), names_to = "dir", values_to = "n") %>%
    mutate(dir = recode(dir, n_sig_up_in_1 = "up in group 1",
                             n_sig_up_in_2 = "up in group 2"))

  p_de_n <- ggplot(ds, aes(n, contrast_f, fill = dir)) +
    geom_col(width = 0.7) +
    geom_text(data = de_summary %>% mutate(contrast_f = factor(contrast, levels = rev(contrast))),
              aes(x = n_sig_up_in_1 + n_sig_up_in_2, y = contrast_f,
                  label = paste0(" ", n_cells_1, " vs ", n_cells_2, " cells")),
              inherit.aes = FALSE, hjust = 0, size = 2.5, colour = "grey30") +
    scale_fill_manual(values = c("up in group 1" = "#b2182b",
                                 "up in group 2" = "#2166ac"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.28))) +
    facet_grid(kind ~ ., scales = "free_y", space = "free_y") +
    labs(title = "Significant genes per contrast",
         subtitle = paste0("q < ", ALPHA, ", ", DE_TEST, " on ", de_assay,
                           "; group 1 is the first term in the contrast name.\n",
                           "enriched_vs_enriched near zero => the two clusters are ",
                           "one population (merge)."),
         x = paste0("genes with q < ", ALPHA), y = NULL) +
    theme_bw(base_size = 9) +
    theme(strip.text.y = element_text(angle = 0, size = 7))
  print(p_de_n)

  # Marker heatmap from the all-vs-rest contrasts, with the TdTomato rate and
  # enrichment call carried on the x axis so cluster identity and TdTomato
  # enrichment are read off the same figure.
  if (any(de_all$kind == "cluster_vs_rest")) {
    mk <- de_all %>%
      filter(kind == "cluster_vs_rest", p_val_adj < ALPHA,
             avg_log2FC > 0, !is_tdtom_feature) %>%
      mutate(cluster = sub("^cl(.*)_vs_rest$", "\\1", contrast)) %>%
      group_by(cluster) %>% slice_max(avg_log2FC, n = DE_MARKER_TOP_N) %>%
      ungroup()
    if (nrow(mk)) {
      gene_ord <- mk %>% arrange(match(cluster, ord_levels), desc(avg_log2FC)) %>%
        pull(gene) %>% unique()
      # Complete the grid: a gene/cluster pair FindMarkers never returned (below
      # min.pct or logfc.threshold in that contrast) has no row, so geom_tile
      # would draw nothing and the cell would read as panel background rather
      # than as missing. Filling it explicitly makes na.value apply.
      mk_full <- de_all %>%
        filter(kind == "cluster_vs_rest", gene %in% gene_ord) %>%
        mutate(cluster = sub("^cl(.*)_vs_rest$", "\\1", contrast)) %>%
        tidyr::complete(gene = gene_ord,
                        cluster = intersect(ord_levels, as.character(enr$cluster)))
      xlab_map <- setNames(
        paste0("cl ", enr$cluster, "\n", sprintf("%.1f%%", enr$pct_tdtom_pos),
               "\n", enr$call),
        as.character(enr$cluster))
      p_mk <- ggplot(mk_full,
                     aes(factor(cluster, levels = ord_levels),
                         factor(gene, levels = rev(gene_ord)),
                         fill = avg_log2FC)) +
        geom_tile(colour = "white", linewidth = 0.2) +
        scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                             midpoint = 0, na.value = "grey72",
                             name = "log2 FC\nvs rest") +
        scale_x_discrete(labels = function(x) xlab_map[x], position = "top") +
        labs(title = paste0("Cluster markers - top ", DE_MARKER_TOP_N,
                            " up genes per cluster vs all other cells"),
             subtitle = paste0("q < ", ALPHA, "; TdTomato excluded. ",
                               "grey = gene not returned for that contrast (below min.pct / logfc).\n",
                               "This is what each cluster IS - read it before the ",
                               "enrichment calls, because an enriched cluster's ",
                               "identity decides whether the call is plausible."),
             x = NULL, y = NULL) +
        theme_bw(base_size = 7) +
        theme(axis.text.x.top = element_text(size = 6),
              axis.text.y = element_text(size = 5.5),
              panel.grid = element_blank())
      print(p_mk)
    }
  }

  # Volcanoes, one panel per contrast. The TdTomato feature is drawn as an open
  # triangle and never labelled: the labels are derived from it, so its
  # appearance at the top of a list is circular, not a finding.
  vol_df <- de_all %>%
    mutate(neglog10q = -log10(pmax(p_val_adj, .Machine$double.xmin)),
           sig = p_val_adj < ALPHA,
           shown = case_when(is_tdtom_feature ~ "TdTomato (excluded)",
                             sig & avg_log2FC > 0 ~ "up in group 1",
                             sig & avg_log2FC < 0 ~ "up in group 2",
                             TRUE ~ "ns"))
  lab_df <- de_all %>%
    filter(p_val_adj < ALPHA, !is_tdtom_feature) %>%
    group_by(contrast) %>%
    slice_max(abs(avg_log2FC), n = DE_TOP_N) %>%
    ungroup() %>%
    mutate(neglog10q = -log10(pmax(p_val_adj, .Machine$double.xmin)))

  p_vol <- ggplot(vol_df, aes(avg_log2FC, neglog10q)) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed",
               colour = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_point(aes(colour = shown, shape = is_tdtom_feature),
               size = 0.8, alpha = 0.75) +
    ggrepel::geom_text_repel(data = lab_df, aes(label = gene), size = 2.1,
                             max.overlaps = Inf, segment.size = 0.15,
                             min.segment.length = 0, box.padding = 0.12) +
    scale_colour_manual(values = c("up in group 1" = "#b2182b",
                                   "up in group 2" = "#2166ac",
                                   "ns" = "grey78",
                                   "TdTomato (excluded)" = "#f0a202"), name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2), guide = "none") +
    facet_wrap(~ contrast, scales = "free", ncol = 3) +
    labs(title = "DE volcanoes",
         subtitle = paste0("dashed = q ", ALPHA,
                           "; top ", DE_TOP_N, " genes per direction labelled. ",
                           "TdTomato itself is drawn as an open triangle and never\n",
                           "labelled - the cell labels derive from it, so its rank is ",
                           "circular."),
         x = "avg log2 FC  (positive = up in group 1)",
         y = expression(-log[10]~q)) +
    theme_bw(base_size = 8) + theme(legend.position = "top")
  print(p_vol)

  # Shared vs unique up-genes: does one tumour programme run through every
  # enriched cluster, or does each have its own?
  vr2 <- de_all %>% filter(kind %in% c("enriched_vs_reference", "ambiguous_vs_reference"),
                           p_val_adj < ALPHA, !is_tdtom_feature, avg_log2FC > 0)
  if (n_distinct(vr2$contrast) > 1) {
    top_shared <- vr2 %>% count(gene, name = "n_contrasts") %>%
      filter(n_contrasts > 1) %>% slice_max(n_contrasts, n = 40)
    hm <- vr2 %>% filter(gene %in% top_shared$gene) %>%
      dplyr::select(contrast, gene, avg_log2FC) %>%
      mutate(gene = factor(gene, levels = rev(top_shared$gene[order(-top_shared$n_contrasts)])))
    if (nrow(hm)) {
      p_shared <- ggplot(hm, aes(contrast, gene, fill = avg_log2FC)) +
        geom_tile(colour = "white", linewidth = 0.3) +
        geom_text(aes(label = sprintf("%.1f", avg_log2FC)), size = 2.2) +
        scale_fill_gradient(low = "#fdd9c8", high = "#b2182b", name = "log2 FC") +
        labs(title = "Up-genes shared across vs-reference contrasts",
             subtitle = paste0("genes significant and up in >1 contrast; ",
                               "a dense column block = one shared programme,\n",
                               "sparse = each cluster has its own"),
             x = NULL, y = NULL) +
        theme_bw(base_size = 8) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
      print(p_shared)
    }
  }
}

# --- UMAP context ------------------------------------------------------------
red_name <- REDUCTION
if (is.null(red_name)) {
  cand <- grep("umap", Reductions(seu), ignore.case = TRUE, value = TRUE)
  red_name <- if (length(cand)) cand[1] else NA_character_
}

if (!is.na(red_name) && red_name %in% Reductions(seu)) {

  message("  UMAP reduction: ", red_name)
  emb <- as.data.frame(Embeddings(seu, red_name))[, 1:2]
  colnames(emb) <- c("umap_1", "umap_2")

  umap_df <- emb %>%
    rownames_to_column("cell_id") %>%
    left_join(meta %>% dplyr::select(cell_id, cluster_id, tdtom_label), by = "cell_id") %>%
    left_join(enr %>% mutate(cluster_id = as.character(cluster)) %>%
                dplyr::select(cluster_id, log2_oe_pos, call, pct_tdtom_pos),
              by = "cluster_id") %>%
    filter(!is.na(cluster_id))

  cluster_centres <- umap_df %>%
    group_by(cluster_id) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

  p_u_cl <- ggplot(umap_df, aes(umap_1, umap_2,
                               colour = factor(cluster_id, levels = cl_levels))) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_colour_manual(values = cluster_colours) +
    ggrepel::geom_label_repel(
               data = cluster_centres,
               aes(x = umap_1, y = umap_2, label = cluster_id),
               colour = "black", fill = "white", size = 2.8, fontface = "bold",
               label.padding = unit(0.12, "lines"), max.overlaps = Inf,
               inherit.aes = FALSE) +
    guides(colour = "none") + labs(title = "Clusters") + theme_classic()

  p_u_td <- ggplot(umap_df %>% arrange(tdtom_label == "TdTom_pos"),
                   aes(umap_1, umap_2, colour = tdtom_label)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_colour_manual(values = TDTOM_COLOURS) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    labs(title = "TdTomato label", colour = NULL) + theme_classic()

  p_u_oe <- ggplot(umap_df, aes(umap_1, umap_2, colour = log2_oe_pos)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b",
                           midpoint = 0, name = "log2 O/E") +
    labs(title = "Cluster TdTom+ enrichment (per-cluster value on cells)") +
    theme_classic()

  p_u_call <- ggplot(umap_df, aes(umap_1, umap_2, colour = call)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_colour_manual(values = CALL_COLOURS) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    labs(title = paste0("Call (BH q < ", ALPHA, ")"), colour = NULL) + theme_classic()

  print(((p_u_cl | p_u_td) / (p_u_oe | p_u_call)) +
          plot_annotation(title = paste0(case_id, " - UMAP context")))

} else {
  message("  No UMAP reduction found - skipping UMAP panels")
  umap_df <- NULL
}

# --- spatial context ---------------------------------------------------------
has_coords <- all(c("imagerow", "imagecol") %in% colnames(seu@meta.data))
if (!has_coords) {
  tissue_coords <- tryCatch(GetTissueCoordinates(seu), error = function(e) NULL)
  if (!is.null(tissue_coords) &&
      all(c("imagerow", "imagecol") %in% colnames(tissue_coords))) {
    seu@meta.data$imagerow <- tissue_coords[rownames(seu@meta.data), "imagerow"]
    seu@meta.data$imagecol <- tissue_coords[rownames(seu@meta.data), "imagecol"]
    has_coords <- TRUE
    message("  Tissue coordinates fetched from image slot")
  }
}

if (has_coords) {
  sp_df <- seu@meta.data %>%
    rownames_to_column("cell_id") %>%
    dplyr::select(cell_id, imagerow, imagecol) %>%
    inner_join(meta %>% dplyr::select(cell_id, cluster_id, tdtom_label), by = "cell_id") %>%
    left_join(enr %>% mutate(cluster_id = as.character(cluster)) %>%
                dplyr::select(cluster_id, log2_oe_pos, call), by = "cluster_id")

  p_sp_td <- ggplot(sp_df %>% arrange(tdtom_label == "TdTom_pos"),
                    aes(imagerow, imagecol, colour = tdtom_label)) +
    geom_point(size = 0.4, alpha = 0.85, shape = 16) +
    scale_colour_manual(values = TDTOM_COLOURS) +
    scale_x_reverse() + coord_fixed() + theme_void() +
    labs(title = paste0(case_id, " - TdTomato label"), colour = NULL)

  p_sp_cl <- ggplot(sp_df, aes(imagerow, imagecol,
                              colour = factor(cluster_id, levels = cl_levels))) +
    geom_point(size = 0.4, alpha = 0.85, shape = 16) +
    scale_colour_manual(values = cluster_colours) +
    scale_x_reverse() + coord_fixed() + theme_void() + guides(colour = "none") +
    labs(title = paste0(case_id, " - clusters"))

  p_sp_oe <- ggplot(sp_df, aes(imagerow, imagecol, colour = log2_oe_pos)) +
    geom_point(size = 0.4, alpha = 0.85, shape = 16) +
    scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b",
                           midpoint = 0, name = "log2 O/E") +
    scale_x_reverse() + coord_fixed() + theme_void() +
    labs(title = paste0(case_id, " - cluster TdTom+ enrichment in situ"))

  print(p_sp_td); print(p_sp_cl); print(p_sp_oe)
} else {
  message("  No imagerow/imagecol available - skipping spatial panels")
}

dev.off()
message("  Plots saved: ", pdf_path)

# =============================================================================
# 10. Tables
# =============================================================================

message("\n=== Step 10: Writing tables ===")

wide_tbl <- tibble(cluster = rownames(tab)) %>%
  mutate(n_cells       = as.integer(row_tot[cluster]),
         n_tdtom_pos   = as.integer(tab[cbind(cluster, "TdTom_pos")]),
         n_tdtom_neg   = as.integer(tab[cbind(cluster, "TdTom_neg")]),
         pct_tdtom_pos = round(100 * n_tdtom_pos / n_cells, 3),
         exp_tdtom_pos = round(exp_mat[cbind(cluster, "TdTom_pos")], 3),
         exp_tdtom_neg = round(exp_mat[cbind(cluster, "TdTom_neg")], 3),
         log2_oe_pos   = round(l2oe_mat[cbind(cluster, "TdTom_pos")], 4),
         log2_oe_neg   = round(l2oe_mat[cbind(cluster, "TdTom_neg")], 4),
         resid_pos     = round(resid_mat[cbind(cluster, "TdTom_pos")], 4),
         resid_neg     = round(resid_mat[cbind(cluster, "TdTom_neg")], 4),
         sample        = case_id)
write.csv(wide_tbl, file.path(OUT_DIR, paste0(case_id, "_cluster_tdtom_matrix.csv")),
          row.names = FALSE)

write.csv(enr %>% mutate(across(where(is.numeric), ~ signif(.x, 6))),
          file.path(OUT_DIR, paste0(case_id, "_cluster_tdtom_enrichment.csv")),
          row.names = FALSE)

write.csv(long %>% dplyr::select(sample, cluster, label, n, expected, log2_oe,
                                 std_resid, pct_of_cluster, n_cells,
                                 q_two_sided, call) %>%
            mutate(across(where(is.numeric), ~ signif(.x, 6))),
          file.path(OUT_DIR, paste0(case_id, "_cluster_tdtom_long.csv")),
          row.names = FALSE)

write.csv(purity_tbl %>% dplyr::select(-obs_rate) %>%
            mutate(across(where(is.numeric), ~ signif(.x, 6))),
          file.path(OUT_DIR, paste0(case_id, "_cluster_purity.csv")),
          row.names = FALSE)

if (nrow(rules_tbl))
  write.csv(rules_tbl, file.path(OUT_DIR, paste0(case_id, "_candidate_rules.csv")),
            row.names = FALSE)

params_tbl <- data.frame(
  parameter = c("sample", "rds", "tdtomato_csv", "label_source", "cluster_col",
                "n_clusters", "n_clusters_tested", "min_cluster_cells",
                "n_cells_used", "n_cells_dropped", "n_tdtom_pos", "n_tdtom_neg",
                "tdtom_pos_rate", "alpha", "multiple_testing", "seed",
                "global_chisq", "global_p", "global_p_method", "cramers_v",
                "n_enriched", "n_depleted", "enriched_clusters",
                "tdtom_feature_check", "order_by",
                "fp_background_rate", "fp_basis", "fp_l2oe_cut",
                "implied_tumour_pct_at_sens100_enriched",
                "best_rule", "best_rule_purity_pct", "best_rule_captured_pct",
                "de_run", "de_mode", "de_reference", "de_reference_clusters", "de_assay",
                "de_test", "de_min_pct", "de_logfc", "de_max_cells",
                "de_n_contrasts", "de_contrasts_skipped",
                "subcluster_run", "subcluster_targets", "subcluster_res",
                "subcluster_de_run", "subcluster_de_mode",
                "subcluster_de_min_cells", "subcluster_de_n_subclusters",
                "subcluster_de_n_scored", "subcluster_de_skipped"),
  value = c(case_id, RDS_PATH,
            ifelse(is.null(TDTOMATO_PATH), NA, TDTOMATO_PATH),
            unique(seu$tdtom_source)[1], CLUSTER_COL,
            nrow(enr), n_tested, MIN_CL_CELLS,
            N, n_dropped, sum(enr$n_tdtom_pos), sum(enr$n_tdtom_neg),
            signif(pos_rate, 6), ALPHA, "BH across tested clusters", SEED,
            signif(unname(gl$statistic), 6), signif(gl$p.value, 6),
            if (use_mc) "Monte-Carlo (1e5)" else "asymptotic",
            signif(cramers_v, 6),
            sum(enr$call == "enriched"), sum(enr$call == "depleted"),
            paste(enr$cluster[enr$call == "enriched"], collapse = ";"),
            tdtom_feature_note, ORDER_BY,
            signif(100 * fp_rate, 6), fp_basis, FP_L2OE_CUT,
            if (!nrow(enr_ranked)) NA else
              paste0(enr_ranked$cluster, "=",
                     sprintf("%.0f%%", 100 * implied_frac(enr_ranked$pct_tdtom_pos / 100, 1)),
                     collapse = ";"),
            if (nrow(rules_tbl)) rules_tbl$rule[which.max(rules_tbl$captured_pct)] else NA,
            if (nrow(rules_tbl)) rules_tbl$purity_pct[which.max(rules_tbl$captured_pct)] else NA,
            if (nrow(rules_tbl)) max(rules_tbl$captured_pct) else NA,
            DO_DE, if (DO_DE) DE_MODE else NA,
            if (DO_DE && DE_MODE != "all_vs_rest") DE_REFERENCE else NA,
            if (DO_DE && exists("ref_clusters")) paste(ref_clusters, collapse = ";") else NA,
            if (DO_DE && exists("de_assay")) de_assay else NA,
            if (DO_DE) DE_TEST else NA, if (DO_DE) DE_MIN_PCT else NA,
            if (DO_DE) DE_LOGFC else NA, if (DO_DE) DE_MAX_CELLS else NA,
            if (nrow(de_summary)) nrow(de_summary) else 0L,
            if (length(de_status)) length(de_status) else 0L,
            DO_SUBCLUSTER,
            if (nrow(sub_tbl)) paste(unique(sub_tbl$parent_cluster), collapse = ";") else NA,
            if (DO_SUBCLUSTER) SUBCLUSTER_RES else NA,
            DO_SUBCLUSTER && SUBCLUSTER_DE,
            if (DO_SUBCLUSTER && SUBCLUSTER_DE) SUBCLUSTER_DE_MODE else NA,
            if (DO_SUBCLUSTER && SUBCLUSTER_DE) SUBCLUSTER_DE_MIN else NA,
            if (nrow(sub_cells)) dplyr::n_distinct(sub_cells$sub_id) else 0L,
            if (nrow(sub_markers))
              sub_markers %>% filter(!is.na(marker_rank)) %>%
                dplyr::pull(sub_id_full) %>% dplyr::n_distinct() else 0L,
            if (length(sub_de_status)) length(sub_de_status) else 0L),
  stringsAsFactors = FALSE)
write.csv(params_tbl, file.path(OUT_DIR, paste0(case_id, "_run_parameters.csv")),
          row.names = FALSE)

message("  Tables saved.")
print(enr %>%
        dplyr::select(cluster, n_cells, n_tdtom_pos, pct_tdtom_pos,
                      log2_oe_pos, resid_pos, odds_ratio, q_two_sided, call) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 4))) %>%
        arrange(desc(log2_oe_pos)), n = 60)

# =============================================================================
# 11. GeoJSON export — tdtomato + seurat_cluster
# =============================================================================
# Same point convention as the full pipeline: coordinates are
# (imagecol, imagerow), one Feature per cell, `classification` holds the value
# the layer is coloured by and `color` the hex for it. Both layers carry the
# cluster enrichment result as properties, so a cell clicked in the viewer
# answers "which cluster, is that cluster TdTom+ enriched, and by how much".
#
# Built column-wise rather than with lapply(seq_len(nrow)): the original
# per-row closure allocates a fresh list per cell, which on a full section
# (>100k cells) dominates the runtime of the whole script.

if (DO_GEOJSON) {

  message("\n=== Step 8: GeoJSON export ===")

  if (!has_coords) {
    message("  No imagerow/imagecol available — skipping GeoJSON export")
  } else {

    geo_df <- seu@meta.data %>%
      rownames_to_column("cell_id") %>%
      dplyr::select(cell_id, imagerow, imagecol) %>%
      inner_join(meta %>% dplyr::select(cell_id, cluster_id, tdtom_label),
                 by = "cell_id") %>%
      left_join(enr %>%
                  mutate(cluster_id = as.character(cluster)) %>%
                  dplyr::select(cluster_id, n_cells, pct_tdtom_pos, log2_oe_pos,
                                resid_pos, odds_ratio, q_two_sided, call),
                by = "cluster_id") %>%
      filter(!is.na(imagerow), !is.na(imagecol))

    n_no_coord <- N - nrow(geo_df)
    if (n_no_coord > 0)
      message("  ", n_no_coord, " labelled cell(s) lack coordinates and are omitted")

    # Shared properties, one vector per field — recycled into every feature.
    base_props <- list(
      cell_id            = geo_df$cell_id,
      tdtom_label        = geo_df$tdtom_label,
      seurat_cluster     = geo_df$cluster_id,
      cluster_n_cells    = geo_df$n_cells,
      cluster_pct_tdtom_pos = round(geo_df$pct_tdtom_pos, 3),
      cluster_log2_oe_pos   = round(geo_df$log2_oe_pos, 4),
      cluster_resid_pos     = round(geo_df$resid_pos, 4),
      cluster_odds_ratio    = signif(geo_df$odds_ratio, 5),
      cluster_q_two_sided   = signif(geo_df$q_two_sided, 5),
      cluster_call          = geo_df$call,
      sample                = case_id)

    make_geojson_points <- function(df, class_vec, colour_map, props) {
      cls <- as.character(class_vec)
      col <- unname(colour_map[cls])
      features <- lapply(seq_len(nrow(df)), function(i) {
        p <- lapply(props, function(v) v[i])
        p$classification <- cls[i]
        if (!is.na(col[i])) p$color <- col[i]
        list(type       = "Feature",
             geometry   = list(type = "Point",
                               coordinates = list(df$imagecol[i], df$imagerow[i])),
             properties = p)
      })
      list(type = "FeatureCollection", features = features)
    }

    message("  Writing TdTomato binary GeoJSON ...")
    gj_tdtom <- make_geojson_points(geo_df, geo_df$tdtom_label,
                                    TDTOM_COLOURS, base_props)
    write(toJSON(gj_tdtom, auto_unbox = TRUE, na = "null", digits = 8),
          file.path(OUT_DIR, paste0(case_id, "_tdtomato.geojson")))

    message("  Writing Seurat cluster GeoJSON ...")
    gj_cluster <- make_geojson_points(geo_df, geo_df$cluster_id,
                                     cluster_colours, base_props)
    write(toJSON(gj_cluster, auto_unbox = TRUE, na = "null", digits = 8),
          file.path(OUT_DIR, paste0(case_id, "_seurat_cluster.geojson")))

    # The cluster -> colour mapping, so the viewer legend can be rebuilt and the
    # PDF panels and the GeoJSON can be shown to use the same palette.
    write.csv(data.frame(cluster = names(cluster_colours),
                         color   = unname(cluster_colours),
                         n_cells = as.integer(enr$n_cells[match(names(cluster_colours),
                                                                as.character(enr$cluster))]),
                         call    = enr$call[match(names(cluster_colours),
                                                  as.character(enr$cluster))],
                         sample  = case_id,
                         stringsAsFactors = FALSE),
              file.path(OUT_DIR, paste0(case_id, "_cluster_colours.csv")),
              row.names = FALSE)

    message("  GeoJSON files saved (", nrow(geo_df), " features each).")
  }
}

if (SAVE_RDS) {
  seu$cluster_tdtom_call    <- enr$call[match(seu$cluster_id, as.character(enr$cluster))]
  seu$cluster_tdtom_log2_oe <- enr$log2_oe_pos[match(seu$cluster_id, as.character(enr$cluster))]
  rds_out <- file.path(OUT_DIR, paste0(case_id, "_tdtom_cluster.rds"))
  saveRDS(seu, rds_out)
  message("  RDS saved: ", rds_out)
}

message("\n=== Done: ", case_id, " ===")
message("Read the matrix as: effect size from log2 O/E, calibration from the ",
        "standardised residual,\n  and the call from BH q on the per-cluster ",
        "Fisher test. A cluster with a large log2 O/E\n  but n below --min_cluster_cells ",
        "is untested, not negative.")

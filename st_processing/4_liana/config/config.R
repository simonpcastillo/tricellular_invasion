# LIANA (cell-type interaction -> replicate consensus -> PVN scoring)
# Relative paths are anchored to the repository root. Override directory roots via
# environment variables or copy this file and set LIANA_CONFIG to its absolute path.
if (!exists("PROJECT_ROOT")) stop("PROJECT_ROOT must be defined before sourcing config.R")
DATA_ROOT <- Sys.getenv("LIANA_DATA_ROOT", file.path(PROJECT_ROOT, "data"))
RESULTS_ROOT <- Sys.getenv("LIANA_RESULTS_ROOT", file.path(PROJECT_ROOT, "results"))
PVN_RESULTS_ROOT <- Sys.getenv("PVN_OUT_ROOT", file.path(PROJECT_ROOT, "..", "3_pvn", "results"))

SAMPLES <- c(mouse1 = "7000_Scan2", mouse2 = "7003_Scan2")
IN_DIR <- file.path(DATA_ROOT, "tables")
EXPR <- setNames(file.path(IN_DIR, paste0(unname(SAMPLES), "_expr_sparse.rds")), names(SAMPLES))
QUPATH <- setNames(file.path(DATA_ROOT, "cell_qupath", paste0(unname(SAMPLES), "_cell_qupath.txt")), names(SAMPLES))
RESOURCE_CSV <- file.path(IN_DIR, "mouseconsensus_v3.csv")
PVN_ASSIGN <- Sys.getenv("PVN_ASSIGN", file.path(PVN_RESULTS_ROOT, "PVN_results", "PVN_v2_assignments.csv"))

# First-stage outputs feed the second stage without manual copying.
LIANA_RESULTS <- file.path(RESULTS_ROOT, "LIANA_results")
LIANA_FIGURES <- file.path(RESULTS_ROOT, "LIANA_figures")
LIANA_OBJECTS <- file.path(RESULTS_ROOT, "LIANA_dat")
CONSENSUS_DIR <- file.path(LIANA_RESULTS, "consensus")
PVN_LR_OUT <- file.path(CONSENSUS_DIR, "outputs_lr_axes_consensus")

# Curation: gene-axis CSV is generated on first use; inspect/edit it then rerun.
# Per-complex overrides are optional; if present, they take precedence.
GENE_AXIS <- file.path(LIANA_RESULTS, "curation", "pvn_gene_axis_assignment_updated.csv")
COMBO_AXIS <- file.path(DATA_ROOT, "annotations", "pvn_lr_combo_axis.csv")
AXIS_FILE <- file.path(PVN_LR_OUT, "pvn_lr_axis_assignment.csv")

MODES <- c("M_T", "M_E", "E_T")
LABEL_COL <- "cell_class"
DROP_LABELS <- c("not_assigned", "NA", "unassigned")
MAX_PVAL <- NULL
MIN_CELLS <- 10L
ORGANISM <- 10090L  # annotation only; resource file must already use mouse symbols
DATA_IS_LOG <- TRUE  # sparse input is SCT 'data' as exported by 2b5; see README
METHODS <- c("natmi", "connectome", "logfc", "sca", "cellphonedb")
SUPP_COL <- c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat", "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR")
EXPR_PROP <- 0.0025
LIANA_N_PERM <- 1000L
N_WORKERS <- 4L
NTOP_DOT <- 25L
SEED_LIANA <- 1234L

# Cross-replicate consensus: historical code's *last* filter was "at least one
# CellPhoneDB p<0.05"; it overwrote a prior Fisher calculation. Default keeps
# that historical selection, but alternative modes are available.
CONSENSUS_FILTER <- "legacy_cpdb_any"  # legacy_cpdb_any | fisher_q | fisher_q_and_cpdb_any
FDR <- 0.05
CPDB_P <- 0.05
SHARED_ONLY <- FALSE
P_FLOOR <- 1e-6
CONSENSUS_MIN_LRSCORE <- NA_real_  # optional additional magnitude filter

# PVN-specific LR scoring; radius must match PVN discovery.
R_NB <- 30
EC_CLASS <- "Endothelial"
Q_MAX <- 0.05
REPLICATING_ONLY <- FALSE
MIN_CELLS_PVN <- 1L
SCALE <- "geomean_linear"  # original: expm1 then sqrt(ligand * receptor)
COMPLEX_FN <- min
MIN_PAIRS_PER_AXIS <- 2L
MIN_COMBOS_PER_AXIS <- 2L
PVN_N_PERM <- 1000L
PERM_WITHIN_REGION <- FALSE
SEED_PVN <- 12345L
EC_ONLY <- FALSE
FOCAL_AXIS <- NA_character_
PVN_LABELS <- character(0)
DROP_IF_LIGAND_EXCLUDED <- FALSE

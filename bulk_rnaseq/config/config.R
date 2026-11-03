# Bulk RNA-seq module configuration
# Relative paths are anchored to BULKRNASEQ_ROOT, which run.sh sets to this directory.
PROJECT_ROOT <- normalizePath(Sys.getenv("BULKRNASEQ_ROOT", unset = "."), mustWork = TRUE)
RDATA <- Sys.getenv("BULKRNASEQ_RDATA", unset = file.path(PROJECT_ROOT, "data", "bulk_counts.RData"))
ROOT_OUT <- Sys.getenv("BULKRNASEQ_OUT", unset = file.path(PROJECT_ROOT, "results"))
EXPORT_DIR <- NULL  # optional existing directory to receive selected summary files

# Input metadata conventions. RData must contain RAW_COUNTS2 and RAW_META2.
# RAW_META2: sample, type, mouse; counts rows are genes, columns are sample IDs.
EC_TYPES <- c("BEC", "CmEC", "ImEC", "IpEC")
CONTROL_TYPES <- c("Bt", "It")
MERGE <- list(mEC = c("CmEC", "ImEC"))
DROP_MOUSE <- NULL
PAIRS_IN <- NULL  # optional precomputed pairwise result with LFC, PADJ and genes

# Signature definition and internal validation
LFC_MIN <- 0.5
ALPHA <- 0.05
N_TOP <- 30
MIN_DET <- 10
NPERM <- 1000

# Figures; controls are displayed but excluded from gene selection
TOP_N <- 5
CONTROLS <- c("Bt", "It")
SHOW <- "samples"       # samples | regions
MARKER_SCALE <- "auto"  # auto | z | absolute
endothelial_markers <- c("Pecam1", "Cdh5", "Kdr", "Emcn", "Vwf")

# GO: GSEA on all tested genes, ORA on selected panel genes
RUN_GO <- TRUE
GO_ORGDB <- "org.Mm.eg.db"  # mouse (change for other species)
GO_ONT <- "BP"             # BP, MF, CC
GO_MINSIZE <- 15
GO_MAXSIZE <- 500
GO_RANK <- "mean_lfc"    # mean_lfc | min_lfc
GO_NPERM <- 10000
GO_TOP_PLOT <- 5
GO_JACCARD <- 0.5
GO_PADJ <- 0.05

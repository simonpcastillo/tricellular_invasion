# PVN analysis configuration — edit this file or set PVN_DATA_ROOT / PVN_OUT_ROOT.
# Paths are relative to the repository unless absolute environment variables are used.
# Two images are currently required by the cross-image PVN discovery method.

PROJECT_ROOT <- .project_root  # resolved from the executing script
DATA_ROOT <- Sys.getenv('PVN_DATA_ROOT', file.path(PROJECT_ROOT, 'data'))
ROOT_OUT  <- Sys.getenv('PVN_OUT_ROOT',  file.path(PROJECT_ROOT, 'results'))

SAMPLES <- c(mouse1 = '7000_Scan2', mouse2 = '7003_Scan2')
FILES <- setNames(file.path(DATA_ROOT, 'cell_qupath', paste0(SAMPLES, '_cell_qupath.txt')), names(SAMPLES))
EXPR <- setNames(file.path(DATA_ROOT, 'tables', paste0(SAMPLES, '_expr_sparse.rds')), names(SAMPLES))
REGION_FRACTIONS_CSV <- file.path(DATA_ROOT, 'tables', 'celltype_region_fractions_v7.csv')

OUT_CSV <- file.path(ROOT_OUT, 'PVN_results')
OUT_FIG <- file.path(ROOT_OUT, 'PVN_figures')
OUT_RES <- file.path(OUT_CSV, 'program_scoring')
OUT_LOG <- file.path(ROOT_OUT, 'logs')
PVN_ASSIGN <- file.path(OUT_CSV, 'PVN_v2_assignments.csv')
HELPERS_DIR <- file.path(PROJECT_ROOT, 'helpers')

# Set to a readable .gmt file to override/add endothelial gene sets; NA = built-ins only.
GMT <- NA_character_

# Analysis defaults retained from original files.
R_NB <- 30          # µm, assumes QuPath centroid coordinates are in µm
MIN_NB <- 2
K_CANDIDATES <- 3:10
MIN_CLUSTER_N <- 15
MIN_IMG_PROP_T <- 0.01
MAX_PCS <- 8
K_SPATIAL <- 5
N_PERM_MORAN <- 999
N_PERM_CENTROID <- 500
N_PERM_PIPELINE <- 200
REPRO_ALPHA <- 0.05
SPLIT_NEFTEL <- FALSE
NEFTEL <- c('MES_like','AC_like','OPC_like','NPC_like')
SEED <- 12345
# Original script overrides shuffled-centroid q95 with fixed 0.90.
# Keep the historical cutoff by default; set NA_real_ for the empirical cutoff.
REPLICATION_MIN_R <- 0.90

# EC program parameters
EC_CLASS <- 'Endothelial'
LABEL_COL <- 'cell_class'
ASSAY <- 'Spatial'
EXPR_IS_LOGNORM <- TRUE  # 2b5 sparse export contains SCT `data`, not raw counts
EC_CLASS_CHECK <- TRUE
REPLICATING_ONLY <- FALSE
NOT_ASSIGNED <- 'not_pvn_assigned'
EXCLUDE_PVN <- character(0)
REGION_COL <- 'region'
NR <- 300L
CTRL <- 100L
MIN_GENES <- 5L
MIN_CELLS <- 20L
Z_GATE <- 2
Q_GATE <- 0.05
DE_NOVO <- TRUE
N_DE_NOVO <- 10L
# Original EC program script uses SEED=1 independently of the PVN clustering seed.
EC_SEED <- 1L

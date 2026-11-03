# Rscript scripts/check_inputs.R [all|inference|consensus|pvn]
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[[1]] else "all"
.args_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
PROJECT_ROOT <- dirname(dirname(normalizePath(sub("^--file=", "", .args_file[[1]]))))
source(Sys.getenv("LIANA_CONFIG", file.path(PROJECT_ROOT,"config","config.R")))
check <- function(f, name) {
  if (!file.exists(f)) stop("MISSING ", name, ": ", f, call. = FALSE)
  cat("OK", name, ": ", f, "\n", sep = "")
}
if (phase %in% c("all", "inference")) {
  check(RESOURCE_CSV, "curated LR resource")
  for (s in unname(SAMPLES)) check(file.path(IN_DIR, paste0(s, "_expr_sparse.rds")), paste0("sparse expression: ", s))
}
if (phase %in% c("all", "consensus")) {
  for (s in unname(SAMPLES)) for (m in MODES)
    check(file.path(LIANA_RESULTS, s, m, "liana_aggregate.csv"), paste("aggregate",s,m))
}
if (phase %in% c("all", "pvn")) {
  check(PVN_ASSIGN, "PVN assignment")
  for (img in names(SAMPLES)) {
    check(QUPATH[[img]], paste("QuPath",img))
    check(EXPR[[img]], paste("expression",img))
  }
  for (m in MODES) check(file.path(CONSENSUS_DIR, m, "liana_consensus.csv"), paste("consensus",m))
  if (!file.exists(GENE_AXIS)) message("NOTE: gene-axis file will be seeded: ", GENE_AXIS)
  if (!file.exists(COMBO_AXIS)) message("NOTE: optional LR-complex axis overrides absent: ", COMBO_AXIS)
}
cat("Input checks completed (", phase, ")\n", sep="")

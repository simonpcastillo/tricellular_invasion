#!/usr/bin/env bash
# Bulk RNA-seq module runner.
#   bash run.sh check   preflight: R packages and input data
#   bash run.sh run     full analysis
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-run}"

export BULKRNASEQ_ROOT="${BULKRNASEQ_ROOT:-$ROOT}"
export BULKRNASEQ_CONFIG="${BULKRNASEQ_CONFIG:-$ROOT/config/config.R}"

if ! command -v Rscript >/dev/null 2>&1; then
  echo 'Rscript not found. Load or install R (>= 4.1) before running.' >&2; exit 1
fi
if [[ ! -f "$BULKRNASEQ_CONFIG" ]]; then
  echo "Configuration not found: $BULKRNASEQ_CONFIG" >&2; exit 1
fi

case "$MODE" in
  check)
    Rscript -e '
      cfg <- Sys.getenv("BULKRNASEQ_CONFIG"); e <- new.env(); sys.source(cfg, envir = e)
      req <- c("limma","edgeR","ggplot2","tidyr","ggnewscale","patchwork","dplyr","scales")
      opt <- c("svglite","fgsea","AnnotationDbi","GO.db", get("GO_ORGDB", envir = e))
      miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
      mopt <- opt[!vapply(opt, requireNamespace, logical(1), quietly = TRUE)]
      cat("ROOT:  ", Sys.getenv("BULKRNASEQ_ROOT"), "\n")
      cat("CONFIG:", cfg, "\n")
      rd <- get("RDATA", envir = e)
      cat(if (file.exists(rd)) "OK      " else "MISSING ", "input RData: ", rd, "\n", sep = "")
      if (length(miss)) cat("Missing required packages:", paste(miss, collapse = ", "), "\n")
      if (length(mopt)) cat("Missing optional packages:", paste(mopt, collapse = ", "),
                            "(SVG export / GO analysis)\n")
      if (getRversion() < "4.1.0") cat("R >= 4.1 required (native |> pipe)\n")
      if (length(miss) || !file.exists(rd))
        stop("Preflight failed.", call. = FALSE)
      message("Preflight passed.")' ;;
  run)
    Rscript "$ROOT/scripts/01_bulk_region_signatures.R" ;;
  *)
    echo 'Usage: bash run.sh [check|run]' >&2; exit 2 ;;
esac

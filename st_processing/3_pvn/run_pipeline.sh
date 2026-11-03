#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-all}"
if ! command -v Rscript >/dev/null 2>&1; then
  echo 'Rscript not found. Load an R module (e.g. on Myriad).' >&2; exit 1
fi
case "$MODE" in
  check)
    echo 'Validating R dependencies and required inputs...'
    Rscript "$ROOT/scripts/check_inputs.R" ;;
  pvn)
    Rscript "$ROOT/scripts/01_pvn_discovery.R" ;;
  ec)
    Rscript "$ROOT/scripts/02_ec_program_scoring.R" ;;
  regions)
    Rscript "$ROOT/scripts/03_region_composition.R"
    Rscript "$ROOT/scripts/04_whole_tissue_neftel.R" ;;
  all)
    Rscript "$ROOT/scripts/check_inputs.R"
    Rscript "$ROOT/scripts/03_region_composition.R"
    Rscript "$ROOT/scripts/04_whole_tissue_neftel.R"
    Rscript "$ROOT/scripts/01_pvn_discovery.R"
    Rscript "$ROOT/scripts/02_ec_program_scoring.R" ;;
  *) echo 'Usage: bash run_pipeline.sh [check|pvn|ec|regions|all]' >&2; exit 2;;
esac

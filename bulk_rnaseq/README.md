# Bulk RNA-seq | endothelial-region signature analysis

A standalone module for **bulk RNA-seq analysis of mouse endothelial-enriched regions**. This module characterises expression differences between region classes, identifies conservative gene signatures, evaluates their ability to distinguish regions across mice, visualises endothelial and control compartments, and performs GO enrichment. It follows the original `1_bulkRNAseq.R` provided with the project; the original remains in `original/` for comparison.

**Status:** Refactored for reusable file paths and configuration. *Not yet executed end-to-end on the original RData*, since the data and an R runtime are not available here. The original statistical code is mostly preserved; see `REVIEW_NOTES.md` for intentional changes and methodological caveats.

## Workflow

```text
Input: LOAD_MYRIAM.RData
  ├─ RAW_COUNTS2  [genes × samples, raw or estimated counts]
  └─ RAW_META2    [sample, type, mouse]
               │
               ▼
       01_bulk_region_signatures.R
               │
      ┌────────┴─────────────────────────────────────────┐
      │                                                  │
      ▼                                                  ▼
 (1) Select regions                             Controls for display
     BEC, CmEC, ImEC, IpEC                     Bt, It (not in DE)
      │                                                  │
      ▼                                                  │
 (2) Merge CmEC + ImEC → mEC (configurable)              │
      │                                                  │
      ▼                                                  │
 (3) limma–voom pairwise comparisons                     │
     mouse-adjusted where design permits                 │
      │                                                  │
      ▼                                                  │
 (4) Tier 1 / Tier 2 signature selection                 │
      │                                                  │
      ├── (5) Leave-one-mouse-out classification         │
      │                                                  │
      └───────────────┬──────────────────────────────────┘
                      ▼
 (6) Region signature / endothelial-marker heatmaps
                      │
                      ▼
 (7) GO GSEA + ORA (optional, RUN_GO = TRUE)
                      │
                      ▼
           results/{tables, figures}
```

## Repository structure

```text
bulkrnaseq_pipeline/
├── README.md
├── REVIEW_NOTES.md
├── run.sh                      # check, dry-run, run
├── config/
│   └── config.R                # project paths and tuning parameters
├── scripts/
│   └── 01_bulk_region_signatures.R
├── jobs/
│   ├── submit_myriad.sh        # submits one SGE job
│   └── bulk_signatures.sge.sh  # batch resources and R module
├── original/
│   └── 1_bulkRNAseq.R          # unchanged historical source
├── data/
│   └── README.md               # place input RData here or point elsewhere
└── results/
    ├── bulkRNAseq_results/    # generated tables and objects
    └── bulkRNAseq_figures/    # generated figures
```

## Script: `01_bulk_region_signatures.R`

### Stage 1 — Import and alignment

Reads `LOAD_MYRIAM.RData` containing `RAW_COUNTS2` (rows = gene symbols, columns = RNA-seq sample IDs) and `RAW_META2` (required columns `sample`, `type`, `mouse`). The script aligns the metadata to the count matrix, drops unmatched count columns, checks group labels and retains a record of the aligned metadata (`input_metadata_aligned.csv`). `mouse` is used for statistical adjustment and cross-validation grouping.

### Stage 2 — Define endothelial and control subsets

Uses the four region types `BEC`, `CmEC`, `ImEC`, and `IpEC` for *signature discovery*. The source script optionally merges `CmEC` and `ImEC` into `mEC`, producing three comparison groups: `BEC`, `mEC`, `IpEC`. The tumour-only regions `Bt` and `It` are loaded for **figure validation only** and excluded from differential expression and gene selection. All class names and the merge are configurable.

### Stage 3 — Pairwise differential expression

Runs `edgeR` library-size normalisation (TMM) and `limma`–`voom`. When the design allows, mouse is included as a blocking factor via `~ 0 + type + mouse`; otherwise a group-only design is used. All ordered region-versus-region contrasts are calculated using empirical Bayes moderation. The script retains log2 fold changes and BH-adjusted p-values for each comparison.

### Stage 4 — Conservative regional signatures

A **worst-case** gene-selection rule is applied: a gene is specific to region *g* only if it is sufficiently upregulated over **every** comparator, or all but one comparator under a relaxed tier.

- **Tier 1:** minimum pairwise log2 FC ≥ `LFC_MIN`, with maximum pairwise adjusted p-value < `ALPHA`.
- **Tier 2:** the second-worst comparison passes, but one named comparator (`hardest`) does not. Genes selected as tier 2 for multiple groups are excluded to avoid conflicting assignments.

The selected panel is capped at `N_TOP` genes per region. The original script writes long-format gene lists, a tier-labelled wide table, and an RDS of the score and classifier objects. **Tier 2 genes should not be treated as unique regional markers.**

### Stage 5 — Leave-one-mouse-out (LOMO) validation

For each held-out mouse, the script refits differential expression and reselects markers from training mice, then compares samples with training-region centroids. It reports a confusion matrix, accuracy by region, majority-class and uniform baselines, and a within-mouse permutation diagnostic. Training-only marker selection avoids direct feature-selection leakage, but see the limitations below regarding sample-combined normalisation and scaling.

### Stage 6 — Signatures and validation figures

Produces a regional gene heatmap, classification confusion matrix, and the final figure combining the selected regional genes with a separate endothelial marker block (`Pecam1`, `Cdh5`, `Kdr`, `Emcn`, `Vwf`). Expression for regional signature genes is scaled using **selection regions only**, then projected onto `Bt` and `It`. Endothelial-marker scaling is separately chosen to show depletion in tumour-only controls without exaggerating noise among EC samples.

### Stage 7 — GO enrichment and GSEA

If `RUN_GO = TRUE`, maps mouse gene symbols to Entrez identifiers using `org.Mm.eg.db`, builds GO gene sets restricted to genes tested in the analysis, and performs:

- **GSEA:** `fgsea` on the ranked, *full tested gene universe* (default `GO_RANK = 'mean_lfc'`); reports NES, adjusted p-values, leading-edge genes and GO plots.
- **ORA:** hypergeometric enrichment of **selected panel genes**, using the detected/mapped universe as background. ORA is descriptive of the selected panel and is not independent validation of it.

The two tests answer different questions. `GO_RANK = 'mean_lfc'` is chosen to preserve a signed, interpretable negative ranking tail. No new empirical results are generated unless this analysis is run on the data.

## Setup

**R dependencies:** `limma`, `edgeR`, `ggplot2`, `tidyr`, `ggnewscale`, `patchwork`, `dplyr`; `scales` is used in plotting via its namespace. Optional: `svglite` for SVG exports and `fgsea`, `AnnotationDbi`, `GO.db`, `org.Mm.eg.db` for GO analysis. The R script uses base R pipes (`|>`), so **R ≥ 4.1** is required.

1. Place your input `LOAD_MYRIAM.RData` in `data/`, or set `BULKRNASEQ_RDATA` to the existing file path.
2. Edit `config/config.R` for group names, signatures, merge policy, visualisation, and GO parameters. Do not edit analysis paths inside the R script.
3. Run the preflight check, then run the analysis.

```bash
bash run.sh check
bash run.sh dry-run
bash run.sh run
```

To keep input data in its existing location and send results elsewhere without editing the repository:

```bash
export BULKRNASEQ_RDATA=/path/to/LOAD_MYRIAM.RData
export BULKRNASEQ_OUT=/path/to/bulk_results
bash run.sh run
```

To disable GO (without deleting any code), edit `RUN_GO <- FALSE` in `config/config.R`.

### Myriad HPC (SGE)

```bash
bash jobs/submit_myriad.sh
```

The SGE template requests one CPU, 32 GB memory, 20 GB temporary space, and eight hours. Adjust in `jobs/bulk_signatures.sge.sh` if needed. By default it loads `r/r-4.4.2_bc-3.20` and uses `$HOME/Scratch/R/libs/spatial_transcriptomics`; replace with your actual R library if different.

## Outputs

The script creates the following directories under `ROOT_OUT`:

| Location | File | Contents |
| --- | --- | --- |
| `bulkRNAseq_results/` | `input_metadata_aligned.csv` | Sample metadata matched to counts |
| `bulkRNAseq_results/` | `region_panel.csv` | Long signature panel with tier, closest rival, and scores |
| `bulkRNAseq_results/` | `worst_case_scores.csv` | Candidate genes and worst/second-worst comparisons |
| `bulkRNAseq_results/` | `region_classifier.rds` | Panel, contrasts, accuracy and confusion matrix |
| `bulkRNAseq_results/` | `loo_predictions.csv` | Held-out per-sample predictions |
| `bulkRNAseq_results/` | `classification_accuracy.csv` | Region classification accuracy |
| `bulkRNAseq_results/` | `signature_wide.csv`, `signature_shapes.rds` | Wide signatures and tier-1-only lists (if a panel exists) |
| `bulkRNAseq_results/` | `figure_genes.csv`, `ec_signatures_20_wide.csv` | Genes used in the figure and wide gene table |
| `bulkRNAseq_figures/` | `region_panel_heatmap.png`, `confusion_matrix.png` | Signature QC and classifier plots |
| `bulkRNAseq_figures/` | `final_signature_heatmap.png`, `.svg` | Final region/marker validation heatmap |
| `bulkRNAseq_results/` | `go_gsea_all.csv`, `go_nes_matrix.csv` | GO GSEA results (optional) |
| `bulkRNAseq_results/` | `go_ora_panel.csv` | Panel GO ORA, when results exist (optional) |
| `bulkRNAseq_figures/` | `go_nes_heatmap.png`, `go_nes_dotplot.png`, `go_running_es.png` | GO plots (conditional) |

## Reproducibility and caveats

**No runtime validation has yet occurred.** The input RData and R/Bioconductor packages are required to verify that the outputs reproduce your previous results. The original script is preserved for diffing.

- **Statistical units:** mice are biological replicates; multiple dissected regions from one mouse are not independent mice. The mouse-adjusted model handles a fixed mouse effect where identifiable; inference should be interpreted in the context of the number of biological replicates.
- **LOMO caveat:** marker selection is nested within folds, but the historical classifier calculates logCPM and standardisation using combined training and held-out data. This is not strictly independent preprocessing and may bias performance. A nested, train-only normalisation/scaling implementation should be tested before reporting prospective accuracy. This was *not silently changed* here.
- **ORA caveat:** the historical ORA calculation applies BH to terms retained with ≥2 overlapping genes rather than to all tested GO terms. This may underestimate FDR. Treat ORA as exploratory until the multiple-testing family is corrected; this refactor preserves the existing statistics.
- **Input counts:** the original script allows non-integer estimated counts but does not automatically round before fitting; check that count input is appropriate for edgeR–voom.
- **Signature transfer:** the wide signature table discards or annotates tiers; downstream scoring should prefer `region_panel.csv` or `signature_shapes.rds$list_tier1` instead of treating all tier-2 genes as unique markers.

See [`REVIEW_NOTES.md`](REVIEW_NOTES.md) for detailed differences from the original source and future validation tasks.

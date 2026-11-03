# LIANA: cell-cell communication and perivascular niche signalling

A three-stage R workflow for **mouse glioblastoma spatial transcriptomics**:
(1) LIANA ligand-receptor inference by cell-type interaction mode, (2) cross-sample
interaction consensus, and (3) endothelial-centred **perivascular niche (PVN)**
ligand-receptor/functional-axis scoring. This pipeline follows the TdTomato/
Neftel classification and PVN discovery pipelines.


## Execution order

```text
Upstream: 2b5_export_tables_from_seurat.R
    results: <sample>_expr_sparse.rds, optionally _plot_df.csv / _umap_df.csv
                           │
        curated mouse LR resource (mouseconsensus_v3.csv)
                           │
                           ▼
           01_liana_inference.R [one job per sample]
                │               │                │
             E_T mode        M_E mode          M_T mode
                └───────────────┼────────────────┘
                                ▼
      LIANA_results/<sample>/<mode>/liana_aggregate.csv
                                │
                                ▼
             02_liana_consensus.R [across samples]
                                │
                                ▼
      LIANA_results/consensus/<mode>/liana_consensus.csv
                                │
       + upstream PVN_results/PVN_v2_assignments.csv
       + original QuPath exports for both sections
       + <sample>_expr_sparse.rds for both sections
       + optional curated gene-axis / LR-complex-axis tables
                                │
                                ▼
             03_liana_pvn_scoring.R
                                │
                                ▼
      PVN ligand-receptor scores, axis scores and reproducibility
```


## Step 1 — LIANA inference (`01_liana_inference.R`)

**Purpose:** infer cell-type-level ligand–receptor interactions separately within
three directed interaction contexts, using the curated mouse LR resource. Runs
independently for each sample.

- **Inputs:** `tables/<sample>_expr_sparse.rds` from the classification export
  (`matrix` genes × cells plus `metadata` with `cell_id`, `cell_class`), optional
  plot and UMAP CSVs, and `mouseconsensus_v3.csv` (exactly two ligand/receptor
  gene-symbol columns).
- **Cell groups:** `cell_class`, removing unassigned calls and groups with fewer
  than 10 cells (configurable). `E_T`: endothelial and four malignant Neftel
  states; `M_E`: microglia/BAM and endothelium; `M_T`: microglia/BAM and four
  malignant states. Modes specify **which groups are included**: LIANA estimates
  all available source→target directions within each subset; direction-specific
  filtering is done downstream.
- **Methods:** `liana_wrap` using NATMI, Connectome, logFC, SCA and CellPhoneDB;
  cell-group expression proportion threshold (`EXPR_PROP`), LIANA ranking and
  permutation settings. The custom mouse resource is read from disk; no
  internet-based resource conversion is attempted.
- **Outputs:** `LIANA_results/<sample>/<mode>/liana_aggregate.csv`,
  `liana_all_methods.csv`, `liana_top_interactions.csv`,
  `liana_significant.csv`, `liana_significant_top_per_pair.csv`; diagnostic PDF
  plots under `LIANA_figures/<sample>/<mode>/`; RDS under `LIANA_dat/<sample>/<mode>/`.

## Step 2 — cross-sample consensus (`02_liana_consensus.R`)

**Purpose:** combine interaction records from samples for each interaction mode.

- **Inputs:** per-sample, per-mode `liana_aggregate.csv` files from Step 1.
- **Method:** join on `source`, `target`, `ligand.complex`, `receptor.complex`;
  retain per-replicate ranks and CellPhoneDB p-values, compute mean LR magnitude
  and specificity, replicate rank concordance, and a historical Fisher-like
  combination of `aggregate_rank` values.
- **Default filter (`CONSENSUS_FILTER="legacy_cpdb_any"`):** reproduce the
  original script's **last** filter: CellPhoneDB p<0.05 in **any sample**.
  Interaction presence in just one sample is allowed; absent replicate p-values
  are not interpreted as failures. **Stage 3 additionally applies `Q_MAX`**
  to the consensus `q_fisher` column by default, even when this historical
  CellPhoneDB filter was used in Step 2. The earlier Fisher filter in the original
  script was overwritten and did not determine the final set.
- **Alternatives:** `fisher_q` or `fisher_q_and_cpdb_any`, applied via the config.
  **Important:** LIANA `aggregate_rank` is a rank-aggregation score, not
  necessarily a calibrated p-value. Consequently `p_fisher`/`q_fisher` are
  *historical heuristic summaries* and should not be reported as formally
  validated combined inferential p/q-values without further calibration.
- **Outputs:** `LIANA_results/consensus/<mode>/liana_replicate_wide.csv`,
  `liana_consensus.csv`, `liana_qc_concordance.csv`, and PDF concordance/dotplots.

## Step 3 — spatial PVN signalling (`03_liana_pvn_scoring.R`)

**Purpose:** identify consensus LR pairs expressed in endothelial-centred
neighbourhoods, quantify per-PVN signalling strength and its between-image
reproducibility.

- **Inputs:** Step 2 consensus tables; `PVN_results/PVN_v2_assignments.csv` from
  `01_pvn_discovery.R`; two QuPath cell exports, and expression matrices from
  Step 1. Cell IDs are matched between PVN assignment, QuPath and expression.
- **Geometry:** rebuild EC-centred edges to neighbours within `R_NB=30` µm.
  **Coordinates must be in micrometres**, matching the PVN discovery radius.
- **Biological axes:** assign LR pairs to functional programmes using curated
  per-gene/complex tables where present, or seed a candidate mapping using
  `AXIS_RULES` (vascular activation, ECM adhesion, immune/antigen, trophic/neural).
  Review ambiguous/unassigned axes **before trusting results**.
- **Score:** resolve complex subunits with the minimum-expression rule; mean
  ligand and receptor expression across sender/receiver cells within PVN,
  transformed to the linear scale (`expm1`) and combined as
  `sqrt(mean_ligand * mean_receptor)` (`SCALE="geomean_linear"`). Compare
  signalling axes between PVNs, optionally permuting PVN anchor assignments
  within anatomical region, and assess sign/FDR concordance across images.
- **Outputs:** `pvn_lr_pair_scores.csv`, `pvn_axis_scores.csv`,
  `pvn_axis_I_raw.csv`, `pvn_axis_coverage.csv`,
  `pvn_partner_cell_counts.csv`, `pvn_axis_reproducible.csv`,
  and `pvn_lr_axis_assignment.csv` under
  `LIANA_results/consensus/outputs_lr_axes_consensus/`.

**Interpretation:** these are *expression-supported candidate* LR associations in the perivascular space,
not evidence of physical contact or functional signalling. Differences in cell
abundance, low expression, sparse neighbourhood classes and curation decisions
can influence scores.

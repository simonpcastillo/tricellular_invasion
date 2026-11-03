# TdTomato rescue → cluster enrichment → GBM state classification

Nine scripts that take QC-filtered Seurat objects, establish which cells carry the
tdTomato lineage reporter, and assign malignant transcriptional states. Every input
and output path is passed on the command line — see the usage block at the top of
each script — so no paths are hardcoded in the source.

## Workflow

```text
INPUT: <sample>_filtered.rds
  │
  └─ 2a1_tdtomato_gmm.R → results/gmm/<sample>/<sample>_gmm_classified.rds
       │
       ├─ labels.mode = "pu" (new model-derived labels):
       │      2a2_export_from_seurat.R → embedding/meta/reporter.csv.gz
       │      2a3_rescue.py             → tumour_labels.csv.gz
       │      2a4_make_tdtomato_rescued_csv.py (--labels)
       │
       └─ labels.mode = "geojson" (re-export existing rescue-labelled GeoJSON):
              2a4_make_tdtomato_rescued_csv.py (--geojson)
                                  │
                                  ▼
results/csv_tomato_gated/<sample>_tdtomato_rescued.csv
                                  │
              2b1_tdtom_cluster_enrichment_v4.R
                        ├─ subcluster_tdtom.csv
                        └─ subcluster_markers.csv
                                  │
              2b2_extract_subcluster_cells.R (replays & STRICTLY validates)
                                  │
                       subcluster_cells.csv
                                  │
                cross-sample 2b4_tdtom_tumourlike_score.py
                                  │
                          tumourlike_units.csv
                                  │
                2b3_neftel_classify_hybrids_v7.R
                   (optional additional subcluster rescue)
                                  │
                     <sample>_classified.rds
                                  │
                2b5_export_tables_from_seurat.R
                                  │
                   UMAP/spatial/expression tables
```

## Step-by-step description

The pipeline has **two parts**: `2a*` establishes TdTomato labels and `2b*`
uses those labels for tumour-enrichment analysis and transcriptional
classification.

### Part A — tdTomato detection and rescue

#### `2a1_tdtomato_gmm.R` — Preprocess cells and establish the initial reporter gate

- **Input:** `<sample>_filtered.rds`, a QC-filtered Seurat object with raw
  expression counts (including a TdTomato feature) and spatial information.
- **Method:** Updates the Seurat object, calculates mitochondrial percentage,
  runs `SCTransform`, PCA, neighbours, clustering and UMAP, and fits a Gaussian
  mixture model (GMM; two components by default) to **SCT-normalised tdTomato
  expression**. The higher-expression component is designated positive, with
  cell labels assigned from the model's posterior probabilities.
- **Output:** `gmm/<sample>/<sample>_gmm_classified.rds`, containing
  `tdtom_label`, `tdtom_expr` and `tdtom_gmm_prob` in cell metadata. Also writes a
  GMM summary, diagnostic plots and a point GeoJSON for the initial gate.


#### `2a2_export_from_seurat.R` — Export features for positive–unlabelled learning

- **Input:** The GMM-classified RDS from `2a1`; a raw reporter feature name;
  optionally a list of confounder-signature genes and a histological-region
  column in the Seurat metadata.
- **Method:** Extracts **raw** TdTomato counts, where any count `> 0` means
  detected reporter. Recomputes PCA using SCT features **excluding the reporter
  gene**, so the rescue model cannot simply learn the reporter measurement.
  Exports library-depth and cluster metadata and, when supplied, confounder
  and spatial information for later quality control.
- **Output:** `pu_inputs/<sample>/embedding.csv.gz`, `meta.csv.gz`,
  `reporter.csv.gz` and an input-provenance record.

#### `2a3_rescue.py` — Infer additional putative reporter-positive cells

- **Input:** The three compressed tables produced by `2a2`.
- **Method:** Fits a **positive–unlabelled (PU)** classifier: cells with detected
  raw TdTomato counts provide positive examples, whereas cells with no detected
  reporter are **unlabelled**, not confirmed negative. Uses reporter-excluded
  expression PCs with sequencing-depth effects removed. Ranks unlabelled cells
  and selects a rescue set using a detection-rate-based prevalence floor.
- **Output:** `pu_rescue/<sample>/tumour_labels.csv.gz`, with each cell's
  `tumour_call` (`probe_positive`, `rescued` or `unlabeled`), plus rescue
  rankings, validation tables, audit files, a summary and provenance.

#### `2a4_make_tdtomato_rescued_csv.py` — Standardise the final tdTomato label file

- **Input:** **Either** the `tumour_labels.csv.gz` from `2a3` (`--labels`, PU
  mode) **or** a pre-existing rescue-labelled QuPath GeoJSON (`--geojson`,
  GeoJSON mode). These are *alternative input routes* to the same output
  contract, not steps to run sequentially.
- **Method:** Maps detailed provenance classes to the canonical labels
  `TdTom_pos` / `TdTom_neg`, normalises capitalisation, verifies unique cell
  IDs, and applies the configured handling for audit-flagged rescues
  (`positive`, `negative` or `exclude`). In GeoJSON mode, it prefers
  `properties.classification.name` over `properties.tdtom_label` by default
  because the latter may still contain the original GMM call.
- **Output:** `csv_tomato_gated/<sample>_tdtomato_rescued.csv` with **exactly**
  the two columns `cell_id,tdtom_label`, and a separate
  `<sample>_tdtomato_provenance.csv` recording `pos`, `pos_rescued`,
  `pos_rescued_suspect` or `neg`.

### Part B — Enrichment, cell-level rescue and state classification

#### `2b1_tdtom_cluster_enrichment_v4.R` — Identify TdTomato-enriched clusters

- **Input:** The `2a1` GMM-classified RDS and the **final** two-column CSV
  produced by `2a4`.
- **Method:** Joins rescued labels to cells and builds a cluster × TdTomato
  contingency table. Computes observed/expected enrichment, standardised
  residuals and per-cluster **Fisher exact tests** with multiple-testing
  correction. Evaluates possible cluster purity, performs differential
  expression (DE), and subclusters selected populations to examine whether
  tumour-associated cells are concentrated in finer groups. The reporter
  feature is excluded from subclustering features.
- **Output:** `enrichment/<sample>/` containing, among others,
  `<sample>_cluster_tdtom_enrichment.csv`, `<sample>_cluster_markers.csv`,
  `<sample>_subcluster_tdtom.csv` and
  `<sample>_subcluster_markers.csv`, alongside diagnostic plots and optional
  spatial GeoJSON layers.

#### `2b2_extract_subcluster_cells.R` — Recover and validate per-cell subcluster membership

- **Input:** The **same** GMM RDS and final rescued CSV used in `2b1`, plus
  `2b1`'s `<sample>_subcluster_tdtom.csv`.
- **Method:** Replays the subclustering stage with matching seed, assay,
  dimensions, resolution and cell limits. **Strictly checks** that each
  reconstructed subcluster has exactly the same total cells and TdTomato-positive
  count as the original `2b1` output; it fails instead of emitting a mismatched
  map.
- **Output:** `enrichment/<sample>/<sample>_subcluster_cells.csv` with cell IDs,
  parent cluster, subcluster and `sub_id_full`.

#### `2b4_tdtom_tumourlike_score.py` — Score tumour-like clusters and subclusters

- **Input:** Marker tables from `2b1` for **all configured samples**, and
  available subcluster reporter-rate/enrichment summaries. It does **not**
  require the per-cell mapping generated by `2b2`.
- **Method:** Evaluates enrichment for tumour-associated transcriptional
  programmes against normal-cell identity panels; uses marker-list-length-matched
  null comparisons and checks for misleading shared markers (for example,
  myeloid or interferon programmes). Assigns interpretable tumour-likeness
  tiers and reports agreement across samples where possible. Reporter rates
  are retained for comparison but **do not enter the tumour-likeness score**.
- **Output:** `tumourlike/tumourlike_units.csv` plus a panel-coverage audit,
  diagnostic plot and, where applicable, cross-sample reproducibility tables.

#### `2b3_neftel_classify_hybrids_v7.R` — Assign GBM states and other cell identities

- **Input:** The GMM RDS, final rescued-label CSV, required Neftel signature
  table, and optional endothelial/myeloid signatures. When **subcluster rescue**
  is enabled, also reads `2b2`'s cell mapping and `2b4`'s
  `tumourlike_units.csv`.
- **Method:** Optionally performs a **second, distinct rescue** by upgrading
  TdTomato-negative cells in selected transcriptionally tumour-like
  subclusters. Scores GBM cellular-state signatures (including hybrid states),
  proliferation and non-malignant cell-type programmes using matched
  permutation-based controls; optionally scores endothelial and myeloid
  substates. Produces cell assignments, diagnostics and spatial exports.
- **Output:** `classified/<sample>/<sample>_classified.rds`, plus per-cell and
  summary CSVs, calibration diagnostics and GeoJSON layers.


#### `2b5_export_tables_from_seurat.R` — Export lightweight analysis products

- **Input:** The fully classified Seurat RDS from `2b3`.
- **Method:** Exports cell metadata joined to UMAP coordinates and tissue
  coordinates, along with a sparse SCT expression matrix and cell metadata.
- **Output:** `tables/<sample>/<sample>_umap_df.csv`,
  `<sample>_plot_df.csv` and `<sample>_expr_sparse.rds`.

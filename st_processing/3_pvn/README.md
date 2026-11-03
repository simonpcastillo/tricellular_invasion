# Perivascular niche (PVN) discovery and endothelial programmes

R workflow for **EC-centred perivascular niche discovery**, anatomical region summaries and endothelial transcriptional programme analysis in two mouse GBM spatial datasets. The inputs are QuPath cell-level classifications/centroid coordinates and sparse expression exports from the upstream TdTomato/Neftel classification pipeline.

## Workflow overview

```text
Upstream classified cells (e.g. 2b5_export_tables_from_seurat.R)
  ├── tables/<sample>_expr_sparse.rds      [SCT expression + metadata]
  └── cell_qupath/<sample>_cell_qupath.txt  [QuPath positions, class, annotations]
                       │
                       ▼
             01_pvn_discovery.R
             • EC-centred neighbours (30 µm)
             • Hellinger-transformed cell-type proportions
             • cluster separately per image; match niche centroids
             • test replicability, enrichment, spatial structure
                       │
                       ▼
             PVN_results/PVN_v2_assignments.csv
                       │
                       ▼
             02_ec_program_scoring.R
             • merge per-EC PVN label with sparse expression
             • curated / optional de novo EC gene programmes
             • expression-matched nulls and diagnostics
                       │
                       ▼
             PVN_results/program_scoring/
               EC_PVN_programs.csv, calibration, confounds

Parallel, independent region-level analyses:
   tables/celltype_region_fractions_v7.csv ─▶ 03_region_composition.R
   cell_qupath/*_cell_qupath.txt ──────────▶ 04_whole_tissue_neftel.R
```


## Description of each step

### 01 — `01_pvn_discovery.R`: derive perivascular niches

**Purpose.** Identify reproducible cellular neighbourhoods around endothelial cells (ECs), rather than treating all vessels as a homogeneous compartment.

**Input:** the two `*_cell_qupath.txt` files, one for each mouse. They must contain `Classification`, a stable `cell_id`, centroid columns beginning `Centroid X` and `Centroid Y`, and signed-distance annotation fields required by the figure and region analyses. Centroid coordinates must be in **micrometres** for the default 30 µm radius to be meaningful.

**Method:**
1. Exclude ambiguous/unassigned classes; map `IpEC`, `mEC`, `BEC` to EC anchors; combine malignant Neftel states as `Tumour` by default, and merge `OPC`/`Oligodendro` into `Oligodendro_OPC`.
2. For each EC, count **non-EC** neighbours within **30 µm** (`R_NB=30`); retain ECs with at least `MIN_NB=2` phenotyped neighbours. Convert counts to relative cell-type proportions.
3. Apply the **Hellinger transformation** (`sqrt(proportion)`), then analyse neighbourhood composition with PCA and independent per-image Ward hierarchical clustering.
4. Evaluate candidate niche numbers `K=3..10` using silhouette, cluster size, per-image representation and cross-image similarity. Match centroids greedily across the two mice and assign shared `PVN_1`, `PVN_2`, etc. labels.
5. Examine niche composition, within-tumour AC/MES/NPC/OPC fractions, microglia enrichment and nulls, niche/region association, EC subtype effects, spatial autocorrelation and maps.


**Outputs** (under `results/PVN_results/`):
- `PVN_v2_assignments.csv` — the central handoff table: endothelial `cell_id`, `Image`, `PVN`, `replicating`, anatomical `region`, `subtype` and neighbour count.
- `PVN_v2_composition.csv` and `PVN_v2_composition_by_image.csv` — average neighbouring cell fractions.
- `PVN_v2_K_evaluation.csv`, `PVN_v2_replication.csv` — cluster selection and correspondence.
- `PVN_v2_microglia_enrichment.csv`, `PVN_v2_microglia_null.csv` — microglia enrichment checks.
- `PVN_v2_microglia_between_pvns.csv`, `PVN_v2_microglia_marginal_region_adjusted.csv` — permutation-based enrichment across PVNs (marginal and region-adjusted).
- `PVN_v2_microglia_glm.csv` — binomial-GLM cross-check.
- `PVN_v2_neftel_between_pvns.csv`, `PVN_v2_neftel_microglia_enriched_vs_rest.csv` — malignant-state composition contrasts between PVNs (group contrast omitted if no PVN is significantly microglia-enriched).
- `PVN_v2_microglia_vs_whole_tissue.csv`, `PVN_v2_neftel_vs_whole_tissue.csv` — perivascular vs whole-tissue tests with a **different permutation null**.
- `PVN_v2_morans_I.csv`, `PVN_v2_ec_retention.csv`, `PVN_v2_model.rds` — spatial structure, retained anchors and model objects.
- Plots (`SVG`, `PDF`) under `results/PVN_figures/`.

**Included helper scripts:**

- `helpers/pvn_enrichment_tests.R` — builds EC-level tables and tests microglia and Neftel-state enrichment **between PVNs**. Permutation tests can be restricted to image, or to image × anatomical region; the helper also provides a binomial GLM cross-check and a group contrast.
- `helpers/pvn_vs_global_tests.R` — compares microglia or malignant-state fractions near ECs with the **whole-tissue / whole-tumour** background by shuffling cell labels. This is a different null/reference from comparisons among PVNs.
- `helpers/pvn_composition_heatmap.R` — shows mean neighbouring cell-type proportions and across-PVN z-scores in composition heatmaps.

All three helper files are included in `helpers/`, and `01_pvn_discovery.R` sources them using its configured `HELPERS_DIR`. The preflight still verifies their presence.

### 02 — `02_ec_program_scoring.R`: EC programme associations with PVNs

**Purpose.** Identify which transcriptional programmes in vascular endothelium vary between the PVN groups identified in step 01.

**Inputs:** `PVN_v2_assignments.csv`; one `<sample>_expr_sparse.rds` per mouse from the upstream sparse-table export, each containing a list with `$matrix` (gene × cell) and `$metadata` (including `cell_id`, preferably `cell_class`). An optional GMT file may add or replace built-in EC gene sets.

**Method:**
1. Match EC anchor `cell_id`s to expression profiles, with an optional cross-check against `cell_class == "Endothelial"`; make IDs unique across images.
2. Restrict to genes shared across both images. Treat the exported **SCT `data` layer** as already normalised, rather than applying another normalisation.
3. Score curated programmes (e.g. arterial/capillary/venous, BBB transport, tip/stalk, inflammation, adhesion, immune recruitment, Notch, TLR-related signalling, EndMT and ECM). Optionally derive exploratory NMF or PCA programmes from a training half.
4. Construct expression/detection-matched random gene-set nulls and estimate PVN effects within each section. Compare across sections, report z-scores, normal and empirical p values, and BH-adjusted calls for curated programmes.
5. Export calibration/confounding diagnostics and a PVN × programme heatmap.

**Outputs** (under `results/PVN_results/program_scoring/`):
- `EC_PVN_programs.csv` — effect/z, adjusted q and call flags.
- `EC_PVN_program_calibration.csv` — null spread and precision checks.
- `EC_PVN_program_confounds.csv` — EC counts, detection and region composition.
- Heatmap: `results/PVN_figures/pvn_heatmap_signatures.svg`.

**Note:** source sparse matrices are SCT-normalised; the current script constructs a Seurat `counts` layer from that matrix for internal operations. These values shouldn't be interpreted as raw UMI counts. Assess the calibration table before interpreting significance.

### 03 — `03_region_composition.R`: anatomical region cell fractions (independent)

**Purpose.** Summarise cell-type proportions across **core**, **infiltrative margin** and **peritumoral** regions for the figure originally labelled *panel 2C*.

**Input:** `tables/celltype_region_fractions_v7.csv`, containing `sample`, `region` (`core`, `im`, `pt`) and columns `tumour`, `endothelial`, `mg`, `bam`, `macper`, `tcell`, `neurons`.

**Method:** Pivot to long format, average fractions across samples, and overlay each sample's fraction on a faceted bar chart.

**Outputs:** `region_fractions_summary.csv` plus `2C_region_fractions_faceted.svg`, `.png`, `.pdf`. The original script also mentioned QuPath-generated panels *2B* and *2D*, which are **not** generated by this R stage.

### 04 — `04_whole_tissue_neftel.R`: tumour-state composition by region (independent)

**Purpose.** Quantify the **whole-tissue** distribution of malignant Neftel states across anatomical parent regions, without restricting to EC neighbours or PVN assignments.

**Input:** the two QuPath cell-level tables with `Classification` (states `MES_like`, `AC_like`, `OPC_like`, `NPC_like`) and a `Parent` column taking `bulk`, `infiltrative`, or `peritumoral`.

**Method:** Count tumour cells by state × region × mouse; explicitly represent zero-count combinations; obtain within-region proportions and average mouse-specific proportions (equal weight per mouse).

**Outputs:** `whole_tissue_neftel_counts.csv`, `whole_tissue_neftel_composition.csv`, and the original comparison plot `neftel_stack_Parent_noPVN.svg` (currently bulk and infiltrative regions only, matching the source plot).


## Setup

Requires R and the packages **Seurat, SeuratObject, Matrix, data.table, ggplot2, patchwork, RANN, spatstat.geom, scales, scattermore, svglite, ggh4x, cluster, dplyr, tidyr, stringr, readr**. The optional NMF stage can use **RcppML**, otherwise the original script falls back to PCA.

1. Edit `config/config.R`: check `SAMPLES`, `R_NB`, `MIN_NB`, replication threshold, programme settings. Default sample names are `7000_Scan2` and `7003_Scan2`.
2. Choose **one** method of providing input/output roots:
   - Keep the repository-relative `data/` and `results/` defaults, or
   - Set `PVN_DATA_ROOT` and `PVN_OUT_ROOT` to absolute locations (e.g. Myriad / shared project storage).
3. Arrange inputs as follows (or change the paths in `config/config.R`):

```text
data/
├── cell_qupath/
│   ├── 7000_Scan2_cell_qupath.txt
│   └── 7003_Scan2_cell_qupath.txt
└── tables/
    ├── 7000_Scan2_expr_sparse.rds
    ├── 7003_Scan2_expr_sparse.rds
    └── celltype_region_fractions_v7.csv
```

4. The three supplied helper scripts are **already bundled** in `helpers/`. No additional copy is needed; see `helpers/README.md`.


## Reproducibility and caveats

- PVN assignment is **two-image-specific**, based on independently clustered mouse datasets and cross-image centroid matching. 
- EC programme calls depend on the null model, number of samples, and the calibration diagnostics; evaluate these before making inferential claims.
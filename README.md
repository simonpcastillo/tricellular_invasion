# Tri-cellular glioblastoma invasion

Code accompanying **"A tri-cellular inflammatory cross-talk drives perivascular invasion in Glioblastoma"**.

Analysis of mouse glioblastoma **Visium HD** spatial transcriptomics, bulk RNA-seq of
endothelium from tumour regions, and immunofluorescence of the invasive margin.

## Repository layout

```text
tricellular_invasion/
│
├── st_processing/                      Visium HD spatial transcriptomics
│   │
│   ├── 0_custom_bin2cell/              2 µm bins → cells, using externally supplied
│   │   ├── run_bin2cell_customseg.py   QuPath/InstanSeg segmentation instead of
│   │   ├── utilsfx/                    bin2cell's built-in StarDist step
│   │   │   ├── gjson_to_npz.py         GeoJSON polygons → sparse label mask
│   │   │   └── morphology.py           per-cell area, circularity, solidity
│   │   ├── requirements.txt
│   │   └── how_to_run.md
│   │
│   ├── 1_b2c_r_qc/                     bin2cell output → per-sample Seurat objects,
│   │   ├── scripts/                    cell-level QC and filtering
│   │   │   ├── 0_b2c_seurat.R          counts + metadata + H&E image → .rds
│   │   │   └── 1_qc.R                  nFeature / nCount / percent.mt filtering
│   │   └── README.md
│   │
│   ├── 2_tdtom_rescue_cell_classification/   tumour-cell identification and
│   │   ├── scripts/                          transcriptional state calling
│   │   │   ├── 2a1_tdtomato_gmm.R            GMM gate on tdTomato reporter
│   │   │   ├── 2a2_export_from_seurat.R      features for the rescue model
│   │   │   ├── 2a3_rescue.py                 positive–unlabelled rescue of
│   │   │   │                                 reporter-negative tumour cells
│   │   │   ├── 2a4_make_tdtomato_rescued_csv.py   canonical tdTomato label table
│   │   │   ├── 2b1_tdtom_cluster_enrichment_v4.R  cluster × tdTomato enrichment, DE
│   │   │   ├── 2b2_extract_subcluster_cells.R     per-cell subcluster map (validated)
│   │   │   ├── 2b3_neftel_classify_hybrids_v7.R   Neftel GBM states + cell types
│   │   │   ├── 2b4_tdtom_tumourlike_score.py      tumour-likeness of (sub)clusters
│   │   │   └── 2b5_export_tables_from_seurat.R    lightweight tables for steps 3–4
│   │   ├── requirements.txt
│   │   └── README.md
│   │
│   ├── 3_pvn/                          perivascular niche (PVN) discovery
│   │   ├── scripts/
│   │   │   ├── 01_pvn_discovery.R      EC-centred 30 µm neighbourhoods → niches
│   │   │   ├── 02_ec_program_scoring.R endothelial programmes across niches
│   │   │   ├── 03_region_composition.R cell fractions by anatomical region
│   │   │   ├── 04_whole_tissue_neftel.R  tumour-state composition by region
│   │   │   └── check_inputs.R          preflight on packages and input files
│   │   ├── helpers/                    enrichment tests, null models, heatmap
│   │   ├── config/config.R             samples, radii, thresholds, paths
│   │   ├── run_pipeline.sh             check | pvn | ec | regions | all
│   │   └── README.md
│   │
│   └── 4_liana/                        ligand–receptor communication
│       ├── scripts/
│       │   ├── 01_liana_inference.R    LIANA per sample, per interaction mode
│       │   ├── 02_liana_consensus.R    cross-sample consensus interactions
│       │   ├── 03_liana_pvn_scoring.R  signalling strength per PVN, by axis
│       │   └── check_inputs.R          preflight on packages and input files
│       ├── config/config.R             samples, modes, filters, paths
│       └── README.md
│
├── bulk_rnaseq/                        bulk RNA-seq of region-sorted endothelium
│   ├── scripts/
│   │   └── 01_bulk_region_signatures.R   region signatures, leave-one-mouse-out
│   │                                     classifier, GO/GSEA, paper figures
│   ├── config/config.R                 region classes, thresholds, GO settings
│   ├── run.sh                          check | run
│   └── README.md
│
├── if_processing/                      immunofluorescence quantification
│   ├── if_bulk_margin_segmentation/    QuPath: section → tumour core, infiltrative
│   │   ├── 0_setupimage.groovy         margin, non-tumour, necrosis
│   │   ├── 1_instanseg.groovy          InstanSeg cell detection + intensity features
│   │   ├── 2_classify.groovy           cells → Tumour / Stroma
│   │   ├── 3_spatialmap.groovy         Gaussian density map of tumour cells
│   │   ├── 4_segmentation_fixed_v2.groovy   density thresholds → the four regions
│   │   └── 5_measure.groovy            per-region area, perimeter, cell counts
│   ├── if_ec_activation_index/
│   │   └── activation_index.ijm        Fiji: VCAM1 / CCL2 per nucleus outside tdTomato
│   └── README.md
│
├── LICENSE                             Apache License 2.0
└── README.md
```

## How to read it

The spatial pipeline runs in order: `0 → 1 → 2 → 3`, with `4` taking input from `2` and `3`.
`bulk_rnaseq/` and `if_processing/` are independent of the spatial stages and of each other.
Within `if_bulk_margin_segmentation/` the Groovy scripts run in their numbered order inside a
QuPath project; `activation_index.ijm` is a standalone Fiji macro.

Each module directory carries its own README (or `how_to_run.md` for stage 0) with the inputs it
expects, the parameters it exposes, the outputs it writes, and the caveats that apply to them.
Stages `3_pvn`, `4_liana` and `bulk_rnaseq` are configured through `config/config.R`; run their
preflight (`check_inputs.R`, or `bash run.sh check`) before the analysis itself.

## System requirements

This is a collection of analysis scripts, not a packaged application. Each module runs
independently; install only what the stages you intend to run require.

**Operating systems.** All scripts are plain R, Python, Groovy and ImageJ macro code with
no platform-specific calls. The R environment recorded below is macOS 14 on Apple Silicon
(`aarch64-apple-darwin23`); the longer-running stages were also run on a Linux HPC cluster.
Not tested on Windows.

**No non-standard hardware is required** — no GPU, and nothing beyond a desktop or a
standard HPC node. Memory is the practical limit: the Visium HD stages
(`st_processing/1`–`4`) hold whole-section Seurat objects in memory and need substantially
more than the other modules. `bulk_rnaseq` is single-threaded and was run with 32 GB.

**R 4.6.0** (`aarch64-apple-darwin23`). Versions used:

| Package | Version | | Package | Version |
|---|---|---|---|---|
| Seurat | 5.5.1 | | liana | 0.1.14 |
| SeuratObject | 5.4.0 | | limma | 3.68.5 |
| Matrix | 1.7.6 | | edgeR | 4.10.5 |
| data.table | 1.18.6.1 | | fgsea | 1.38.0 |
| ggplot2 | 4.0.3 | | AnnotationDbi | 1.74.0 |
| patchwork | 1.3.2 | | GO.db | 3.23.1 |
| dplyr | 1.2.1 | | org.Mm.eg.db | 3.23.0 |
| tidyr | 1.3.2 | | mclust | 6.1.3 |
| readr | 2.2.0 | | homologene | 1.4.68.19.3.27 |
| stringr | 1.6.0 | | RANN | 2.6.3 |
| tibble | 3.3.1 | | spatstat.geom | 3.8.2 |
| scales | 1.4.0 | | scattermore | 1.2 |
| svglite | 2.2.2 | | ggh4x | 0.3.1 |
| cluster | 2.1.8.3 | | ggrepel | 0.9.8 |
| jsonlite | 2.0.0 | | png | 0.1.9 |
| future | 1.75.0 | | purrr | 1.2.2 |
| ggnewscale | 0.5.2 | | optparse | (module 2 only) |
| RcppML | 0.3.7.1 (optional; PCA fallback) | | | |

**Python 3** for `st_processing/0` and the two Python stages of `st_processing/2`.
Package names are in the per-module `requirements.txt`; see *Installation* below.

**QuPath v0.7** with the **InstanSeg** extension and the
`fluorescence_nuclei_and_cells-0.1.1` model, for `if_processing/if_bulk_margin_segmentation/`.
**Fiji/ImageJ** for `if_processing/if_ec_activation_index/`.

## Installation

No build or compilation step — clone the repository and install the interpreters and
packages listed above.

```bash
git clone https://github.com/simonpcastillo/tricellular_invasion.git
cd tricellular_invasion
```

R packages (CRAN and Bioconductor):

```r
install.packages(c("Seurat","Matrix","data.table","jsonlite","png","ggplot2","patchwork",
                   "dplyr","tidyr","readr","stringr","tibble","scales","svglite","cluster",
                   "RANN","spatstat.geom","scattermore","ggh4x","ggrepel","mclust",
                   "optparse","homologene","future","purrr","ggnewscale"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma","edgeR","fgsea","AnnotationDbi","GO.db","org.Mm.eg.db"))
remotes::install_github("saezlab/liana")          # not on CRAN
```

Python packages, per module:

```bash
pip install -r st_processing/0_custom_bin2cell/requirements.txt
pip install -r st_processing/2_tdtom_rescue_cell_classification/requirements.txt
```

QuPath and Fiji are installed from their own distributions; the InstanSeg model is
downloaded from inside QuPath via *Extensions > InstanSeg*.

Install time on a desktop depends almost entirely on whether CRAN/Bioconductor binaries
are available for your platform or the packages are built from source.

Before running a configured module, use its preflight to confirm the environment:

```bash
Rscript st_processing/3_pvn/scripts/check_inputs.R
Rscript st_processing/4_liana/scripts/check_inputs.R
bash bulk_rnaseq/run.sh check
```

## Demo

All data these scripts consume — both the spatial transcriptomics and the bulk
RNA-seq — are deposited on Zenodo:

> **TODO before submission — replace with the real DOI**
> Zenodo: https://doi.org/10.5281/zenodo.XXXXXXX

The deposit is larger than a demo needs, so cut a small subset first. The quickest
end-to-end check is `st_processing/1_b2c_r_qc`, which is single-sample and self-contained.

**1. Fetch and subset.** Download the bin2cell output for one sample from the Zenodo
record, then crop it to a few thousand cells:

```r
library(Matrix)
src <- "7000_Scan2"; dst <- "demo/input/bin2cell/demo_sample"
dir.create(dst, recursive = TRUE, showWarnings = FALSE)

m    <- readMM(file.path(src, "matrix.mtx"))
bc   <- readLines(file.path(src, "barcodes.tsv"))
keep <- seq_len(min(3000, ncol(m)))                 # first 3000 cells

writeMM(m[, keep], file.path(dst, "matrix.mtx"))
writeLines(bc[keep], file.path(dst, "barcodes.tsv"))
file.copy(file.path(src, "features.tsv"), dst)

for (f in c("cell_labels.csv", "cdata_spatial_fullres.csv")) {
  d <- read.csv(file.path(src, f), check.names = FALSE)
  write.csv(d[d[[1]] %in% bc[keep], ], file.path(dst, f), row.names = FALSE)
}
```

Both `cell_labels.csv` and `cdata_spatial_fullres.csv` are keyed by cell ID in their
first column and are matched against the barcodes, which is why filtering on that
column keeps the three files consistent.

Copy the matching Space Ranger `spatial/` folder (`scalefactors_json.json` and
`tissue_hires_image.png`) to `demo/input/spaceranger/demo_sample/spatial/`.

**2. Run.**

```bash
export PIPELINE_INPUT_DIR=$PWD/demo/input
export PIPELINE_OUTPUT_DIR=$PWD/demo/output
Rscript st_processing/1_b2c_r_qc/scripts/0_b2c_seurat.R
Rscript st_processing/1_b2c_r_qc/scripts/1_qc.R
```

**3. Expected output.** Under `demo/output/`:

| File | Contents |
|---|---|
| `00_seurat/demo_sample.rds` | Seurat object, 3000 cells, spatial coordinates in the `VisiumV1` image slot |
| `01_qc/demo_sample_filtered.rds` | The same object after QC filtering |
| `01_qc/demo_sample_QC.pdf` | Violin and scatter QC plots |
| `01_qc/QC_summary.csv` | One row: cells before and after filtering, and the retained fraction |

Cells are retained when `nFeature_Spatial > 50`, `nCount_Spatial > 50` and
`percent.mt < 20`, so `QC_summary.csv` should show a modest drop rather than a
near-total loss. Step 0 also prints the percentage of coordinates falling inside the
image, which should be close to 100%.

> **TODO before submission** — run the above once and record the demo run time on a
> desktop here.

**Bulk RNA-seq.** The bulk module needs no subsetting — the counts object from the same
Zenodo record is already small. Point the module at it and run:

```bash
export BULKRNASEQ_RDATA=/path/to/bulk_counts.RData
export BULKRNASEQ_OUT=$PWD/demo/bulk_results
bash bulk_rnaseq/run.sh check      # confirms packages and input
bash bulk_rnaseq/run.sh run
```

`check` prints the resolved paths and reports any missing package before the analysis
starts. A successful run writes `region_panel.csv`, `classification_accuracy.csv` and
`final_signature_heatmap.png` (plus the GO tables when `RUN_GO = TRUE`) under
`bulkRNAseq_results/` and `bulkRNAseq_figures/`.

**Beyond stage 1.** `3_pvn` and `4_liana` need **both** sections, because PVN discovery
clusters each image separately and matches niche centroids across them; a one-sample
demo cannot exercise them. Run those on the full deposit, using
`Rscript st_processing/3_pvn/scripts/check_inputs.R` to confirm the inputs first.

## Data availability

All data required to run this code are deposited in a single Zenodo record:

> **TODO before submission — replace with the real DOI, and add any raw-sequence
> repository accessions (e.g. GEO/ArrayExpress) if the raw reads are deposited separately**
> Zenodo: https://doi.org/10.5281/zenodo.XXXXXXX

The record contains:

- **Spatial transcriptomics** — bin2cell outputs per section (`matrix.mtx`,
  `features.tsv`, `barcodes.tsv`, `cell_labels.csv`, `cdata_spatial_fullres.csv`), the
  matching Space Ranger `spatial/` folders, QuPath cell exports
  (`<sample>_cell_qupath.txt`) and the sparse expression objects
  (`<sample>_expr_sparse.rds`) consumed by `st_processing/3_pvn` and `4_liana`.
- **Bulk RNA-seq** — the counts object read by `bulk_rnaseq/`, containing the raw count
  matrix `RAW_COUNTS2` (genes × samples) and the sample table `RAW_META2`
  (columns `sample`, `type`, `mouse`).

## Citation

*To be completed on publication.*

## Contact

Lead author: Simona Parrinello, UCL Cancer Institute, UCL — s.parrinello@ucl.ac.uk

Repository maintainer: Simon Castillo, UCL Cancer Institute, UCL — simon.castillo@ucl.ac.uk

## License

Released under the **Apache License 2.0**, an [OSI-approved](https://opensource.org/licenses/Apache-2.0)
licence — see [LICENSE](LICENSE). You may use, modify and redistribute this code, including
commercially, provided you retain the copyright and licence notices and state any changes.
The licence includes an express grant of patent rights from the contributors and provides
no warranty.

Source code: https://github.com/simonpcastillo/tricellular_invasion

## Generative AI Disclosure

Generative artificial intelligence (GenAI) was used to assist with code auditing, refactoring, documentation, and workflow organisation. This included reviewing code structure, identifying potential errors, improving readability, and standardising file paths and dependencies.

The underlying analytical approaches, scientific rationale, and interpretation of results remain the responsibility of the authors.


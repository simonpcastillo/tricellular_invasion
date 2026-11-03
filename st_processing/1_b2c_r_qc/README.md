# bin2cell → Seurat → quality control

Two sequential R scripts to convert bin2cell outputs to per-sample Seurat objects and perform cell-level QC. Paths are **relative to the repository root**, so the same code can run locally or on HPC without editing usernames or absolute paths.

## Repository layout

```text
1_b2c_r_qc/
├── scripts/
│   ├── 0_b2c_seurat.R
│   └── 1_qc.R
├── input/
│   ├── bin2cell/
│   │   └── <sample_id>/
│   │       ├── matrix.mtx
│   │       ├── features.tsv
│   │       ├── barcodes.tsv
│   │       ├── cell_labels.csv
│   │       └── cdata_spatial_fullres.csv
│   └── spaceranger/
│       └── <sample_id>/
│           └── .../spatial/
│               ├── scalefactors_json.json
│               └── tissue_hires_image.png
└── output/                         # generated; excluded from git
    ├── 00_seurat/
    │   └── <sample_id>.rds
    └── 01_qc/
        ├── <sample_id>_filtered.rds
        ├── <sample_id>_QC.pdf
        └── QC_summary.csv
```

`<sample_id>` can be any directory name (e.g. `7000_Scan2`), **but it must match** in `input/bin2cell/` and `input/spaceranger/`. The Space Ranger files may occur anywhere below that sample's directory (including `binned_outputs/square_002um/spatial/`). If multiple bins are found, the scripts prefer `square_002um`; an ambiguous match is rejected.

`features.tsv` is expected to contain gene symbols in its **second column** (standard 10x format). The scripts make duplicate symbols unique and use human or mouse mitochondrial gene prefixes (`MT-`/`mt-`) when calculating `percent.mt`. If no mitochondrial genes can be found, the pipeline stops rather than reporting a misleading QC value.

Create the input folders after cloning, then place each sample under its matching ID:

```bash
mkdir -p input/bin2cell input/spaceranger
```

## Install requirements

Install the R packages `Seurat`, `Matrix`, `data.table`, `jsonlite`, `png`, `ggplot2`, and `patchwork` in your normal R library. No machine-specific `.libPaths()` setting is needed.

## Run the pipeline

From the module root (`st_processing/1_b2c_r_qc/`):

```bash
Rscript scripts/0_b2c_seurat.R
Rscript scripts/1_qc.R
```

Each script resolves its default paths from its **own location**, so by default `input/` and `output/` are expected inside `scripts/`; the commands also work when launched from another working directory. Use the environment variables below to put the data anywhere else. To keep large data outside the Git repository, set one or both shared path environment variables before invoking either script:

```bash
export PIPELINE_INPUT_DIR=/path/to/input
export PIPELINE_OUTPUT_DIR=/path/to/output
Rscript scripts/0_b2c_seurat.R
Rscript scripts/1_qc.R
```

The output of Step 0 (`output/00_seurat/*.rds`) is exactly the input of Step 1. The scripts skip existing outputs; Step 1 also preserves and incrementally updates `QC_summary.csv`. Remove existing files if you deliberately want a full re-run with new settings.

## QC rules

Step 1 retains a cell when all conditions hold:

- `nFeature_Spatial > 50`
- `nCount_Spatial > 50`
- `percent.mt < 20`

Change the three `QC_*` constants near the beginning of `1_qc.R` as required. The thresholds are intentionally **strict inequalities** to match the original script.

## Spatial image conventions

The bin2cell `cdata_spatial_fullres.csv` file must contain unique cell IDs in its first column, plus `x_full` and `y_full` coordinates in full-resolution image pixels. Step 0 stores these unscaled positions in the Seurat `VisiumV1` image object, together with the **hires** Space Ranger scale factor. Because only the hires PNG is stored in the image slot, both the `hires` and `lowres` **plot scale selections** deliberately use the hires factor so that Seurat's default spatial plots align. The untouched original scale-factor JSON is also saved as `obj@misc$spaceranger_scalefactors`. The displayed image is `tissue_hires_image.png`; the marker diameter is a display placeholder, not a measured cell size. Check the reported percentage of coordinates inside the image and verify image registration if it is low.

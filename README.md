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

## Requirements

R (Seurat v5, LIANA, limma/edgeR and the packages named in each module README) and Python 3
(see the `requirements.txt` files). The immunofluorescence module needs QuPath with the
InstanSeg extension, and Fiji/ImageJ. Package lists are given per module; install only what
the stages you intend to run require.

## Data availability

*To be completed on publication — accession numbers for the Visium HD and bulk RNA-seq data.*

## Citation

*To be completed on publication.*

## Contact

Lead author: Simona Parrinello, UCL Cancer Institute, UCL — s.parrinello@ucl.ac.uk

Repository maintainer: Simon Castillo, UCL Cancer Institute, UCL — simon.castillo@ucl.ac.uk

## License

Apache License 2.0 — see [LICENSE](LICENSE).

## Generative AI Disclosure

Generative artificial intelligence (GenAI) was used to assist with code auditing, refactoring, documentation, and workflow organisation. This included reviewing code structure, identifying potential errors, improving readability, and standardising file paths and dependencies.

The underlying analytical approaches, scientific rationale, and interpretation of results remain the responsibility of the authors.


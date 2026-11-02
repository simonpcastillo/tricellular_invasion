## Running batch

Run the pipeline using:

```bash
python run_bin2cell_customseg.py \
    --spaceranger-root /st/spaceranger_outs \
    --geojson-root /st/segmentation_files \
    --outdir /st/bin2cell_updated2 \
    --skip-existing \
    --tile-size 300 \
    --tile-overlap 20
```

## Data structure

### Input directories

```text
spaceranger_outs/                           # --spaceranger-root
└── <sample>/
    ├── sr_wsi.tif                          # INPUT: Whole-slide image
    │                                       # Same image used for Space Ranger
    │
    └── binned_outputs/
        └── square_002um/
            ├── (Space Ranger count files)  # INPUT: Gene expression counts
            └── spatial/                    # INPUT: Spatial information

segmentation_files/                         # --geojson-root
└── <sample>.geojson                        # INPUT: Cell segmentations
```

### Output directory

```text
bin2cell_updated2/                          # --outdir
└── <sample>/
    │
    ├── rescaled_image/
    │   └── he.tiff                         # Rescaled H&E image (mpp-aligned)
    │
    ├── bin2cell_xform_full_to_crop.json    # Coordinate transform
    │                                       # (full → cropped)
    │
    ├── <sample>_cell_instanseg.npz          # Sparse integer label mask
    │
    ├── cell_labels.csv                     # Cell-level metadata
    ├── matrix.mtx                          # Cell × gene count matrix
    ├── barcodes.tsv                        # Cell barcodes
    ├── features.tsv                        # Gene names
    ├── spatial_coords.csv                  # Cropped-space coordinates (x, y)
    ├── cdata_spatial_fullres.csv           # Full-resolution coordinates
    │                                       # (x_full, y_full)
    │
    └── visualisations/
        └── tile_r{r0}-{r1}_c{c0}-{c1}.png   # One PNG per non-empty tile
```

### Output file descriptions

| File | Description |
|---|---|
| `he.tiff` | Rescaled H&E image aligned to the target microns-per-pixel resolution. |
| `bin2cell_xform_full_to_crop.json` | Coordinate transformation between full-resolution and cropped images. |
| `<sample>_cell_instanseg.npz` | Sparse integer label mask containing cell segmentation IDs. |
| `cell_labels.csv` | Cell-level metadata and morphological measurements. |
| `matrix.mtx` | Cell × gene expression count matrix in Matrix Market format. |
| `barcodes.tsv` | Cell identifiers. |
| `features.tsv` | Gene identifiers and names. |
| `spatial_coords.csv` | Cell coordinates in cropped-image space (`x`, `y`). |
| `cdata_spatial_fullres.csv` | Cell coordinates in full-resolution image space (`x_full`, `y_full`). |
| `visualisations/` | Tile-level quality-control images showing segmentation overlays. |

### Cell-level metadata

The `cell_labels.csv` file contains the following columns:

| Column | Description |
|---|---|
| `object_id` | Unique cell segmentation identifier. |
| `bin_count` | Number of spatial bins assigned to the cell. |
| `array_row` | Spatial array row. |
| `array_col` | Spatial array column. |
| `labels_joint_name` | Cell label identifier. |
| `area_um2` | Cell area in µm². |
| `circularity` | Cell circularity measurement. |
| `solidity` | Cell solidity measurement. |
| `aspect_ratio` | Cell aspect ratio. |

### Segmentation visualisations

Each non-empty tile produces a PNG file named:

```text
tile_r{r0}-{r1}_c{c0}-{c1}.png
```

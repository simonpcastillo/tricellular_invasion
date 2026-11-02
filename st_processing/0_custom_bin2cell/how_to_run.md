## Running batch
``python run_bin2cell_customseg.py \
    --spaceranger-root /st/spaceranger_outs \
    --geojson-root     /st/segmentation_files\
    --outdir           /st/bin2cell_updated2 \
    --skip-existing --tile-size 300 --tile-overlap 20
	``

## Data structure
spaceranger_outs/                          ← --spaceranger-root
└──  <sample>/                              
    ├── sr_wsi.tif                         INPUT  whole-slide image, this is the same used for running spaceranger
    └── binned_outputs/square_002um/
        ├── (Space Ranger count files)     INPUT
        └── spatial                        INPUT  Space Ranger spatial files

geojson_outputs/                           ← --geojson-root
└── <sample>.geojson                       INPUT  cell segmentations.

bin2cell_output/                           ← --outdir
└── <sample>/
    │
    ├── rescaled_image/
    │   └── he.tiff                        rescaled H&E image (mpp-aligned)
    │
    ├── bin2cell_xform_full_to_crop.json    coordinate transform (full → cropped)
    ├── <sample>_cell_instanseg.npz         sparse integer label mask
    │
    ├── cell_labels.csv                     cell-level metadata
    │                                            (object_id, bin_count, array_row,
    │                                       array_col, labels_joint_name,
    │                                       area_um2, circularity, solidity,
    │                                       aspect_ratio, class)
    ├── matrix.mtx                          cell × gene count matrix (MTX)
    ├── barcodes.tsv                        cell barcodes
    ├── features.tsv                        gene names
    ├── spatial_coords.csv                  cropped-space coordinates (x, y)
    ├── cdata_spatial_fullres.csv           full-res coordinates (x_full, y_full)
    │
    └── visualisations/
        └── tile_r{r0}-{r1}_c{c0}-{c1}.png   one PNG per non-empty tile
                                                   left panel  = H&E only
                                                   right panel = H&E + overlay
                                                   (red = PGCC, blue = non-PGCC)
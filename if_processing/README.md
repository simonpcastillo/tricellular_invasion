# Immunofluorescence processing

Two independent image-analysis workflows. Neither depends on the spatial
transcriptomics pipeline, and they do not depend on each other.

| Folder | Platform | Purpose |
|---|---|---|
| `if_bulk_margin_segmentation/` | QuPath (Groovy) | Segment a section into tumour core, infiltrative margin, non-tumour and necrosis, then measure each region |
| `if_ec_activation_index/` | ImageJ / Fiji (macro) | VCAM1 and CCL2 signal per nucleus in the tdTomato-negative compartment |

---

## `if_bulk_margin_segmentation/` — region segmentation in QuPath

Six scripts, run in order on an image **inside a QuPath project**. The project is
required: Script 3 saves the density map as a project resource and Script 4 loads
it back by name.

| Script | Does |
|---|---|
| `0_setupimage.groovy` | Sets the image type to fluorescence, names the channels, reports pixel calibration |
| `1_instanseg.groovy` | InstanSeg cell detection over a whole-image annotation, then intensity features per cell |
| `2_classify.groovy` | Splits cells into `Tumour` / `Stroma` on one intensity measurement |
| `3_spatialmap.groovy` | Gaussian-weighted density map of `Tumour` cells, saved to the project |
| `4_segmentation_fixed_v2.groovy` | Thresholds that density map into `Tumour`, `Infiltrative margin`, `Non-tumour`, `Necrosis` |
| `5_measure.groovy` | Per-region area, perimeter and cell counts; exports CSV and GeoJSON |

### Before you run

- **Channel names** (`0_setupimage.groovy`) must match your acquisition and the
  channel order in the image. The scripts as supplied expect
  `tdtomato`, `cd31`, `dapi`, `iba1`.
- **Pixel calibration in microns** must be set. Every area, radius and distance
  threshold downstream is in µm or µm², so an uncalibrated image gives
  meaningless regions. Script 0 reports the calibration and tells you how to set it.
- **InstanSeg model path** (`1_instanseg.groovy`): edit `modelPath`, or set the
  `INSTANSEG_MODEL` environment variable, to the model directory downloaded
  through Extensions > InstanSeg. The script stops with a message if the path
  does not exist.
- **Measurement name** (`2_classify.groovy`): `measName` must match the
  measurement column exactly, including the pixel-size prefix that
  `IntensityFeaturesPlugin` adds. Check one cell in the Measurements table first.

### How the regions are defined

Script 3 builds a Gaussian-weighted density of `Tumour`-classified cells
(default radius 100 µm). Script 4 cuts that surface at two levels:

- **core** (`T_high`) → annotated `Tumour`
- **envelope** (`T_low`) → the outer tumour territory
- **`Infiltrative margin`** = envelope − core, so the margin is a ring that
  excludes the core by construction
- **`Non-tumour`** = tissue − (envelope ∪ core)
- **`Necrosis`** = large low-density holes fully enclosed by the core

By default the thresholds are **percentiles of each slide's own density
distribution** (`thresholdMode = 'percentile'`), which keeps the definition
comparable across sections of differing overall brightness. Set
`thresholdMode = 'absolute'` to use fixed `T_high` / `T_low` values instead.

Script 4 then discards tumour foci that contain no dense core
(`requireCoreInFocus`) and foci below `minFocusArea`, which removes scattered
positive specks without a hand-tuned size cutoff.

### Counting convention

`5_measure.groovy` assigns a cell to a region when its **centroid** falls inside
that region's ROI, not by hierarchy nesting. Because the margin is a ring, its
count excludes the core geometrically.

`Necrosis` is a subregion *inside* the `Tumour` core, so cells in necrosis are
also counted under `Tumour`; necrosis is reported separately for reference. If
you need the two to be exclusive, subtract the necrosis row.

### Outputs

Written to `<project>/exports/`:

- `<image>_region_summary.csv` — area (mm²), perimeter (mm), cell count, tumour+ count and percentage per region
- `<image>_regions.geojson` — the region annotations
- `<image>_cells.csv` — per-cell measurements

---

## `if_ec_activation_index/` — endothelial activation index in Fiji

`activation_index.ijm` operates on one open multichannel image. It thresholds the
tdTomato channel and inverts it to obtain the tdTomato-negative compartment,
counts DAPI nuclei there, and measures mean VCAM1 and mean CCL2 intensity over
the same region, reporting each both raw and divided by the nucleus count.

Channel indices are 1-based and set at the top of the macro
(`vcamCh`, `tdtCh`, `ccl2Ch`, `dapiCh`), as are the autothreshold methods
(`Triangle` for tdTomato, `Otsu` for DAPI), the Gaussian sigmas and the nucleus
area bounds. Edit them to match your acquisition — the defaults assume a
different channel order from the QuPath scripts above.

Results appear in the Fiji Results table; the macro does not write a file, so
save or append the table yourself when batching.

## Requirements

- **QuPath** ≥ 0.5 with the **InstanSeg** extension and a downloaded fluorescence
  cell model. Script 4 uses JTS geometry classes bundled with QuPath.
- **Fiji / ImageJ** for the `.ijm` macro. No plugins beyond the standard
  distribution.

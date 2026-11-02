#!/usr/bin/env python3
"""
bin2cell segmentation replacement pipeline  (batch mode)
=========================================================
Loops over every sample directory found under ``--spaceranger-root`` and
runs the full bin2cell pipeline for each one.

Expected directory layout
--------------------------
    <spaceranger-root>/
        <sample-name>/
            sr_wsi.tif                          ← whole-slide image
            binned_outputs/square_002um/
                spatial/

    <geojson-root>/
        <sample-name>_pgcc_cellvit_multiclass_seg.geojson

The GeoJSON filename is: ``<sample-name> + --geojson-suffix``.

Usage
-----
    python run_bin2cell.py --help
    python run_bin2cell.py \\
        --spaceranger-root /path/to/spaceranger_outs \\
        --geojson-root     /path/to/cellvitpp_outputs \\
        --outdir           /path/to/bin2cell_updated2

Optional per-sample override (runs a single sample instead of all)
    python run_bin2cell.py ... --sample 241512-HD01
"""

import argparse
import json
import logging
import os
import sys
import traceback

import matplotlib
matplotlib.use("Agg")  # non-interactive backend — safe for batch/HPC
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scanpy as sc
import tifffile as tiff
from scipy import io, sparse

import bin2cell as b2c
from utilsfx.gjson_to_npz import geojson_to_sparse_labels_npz
from utilsfx.morphology import geojson_polygon_morphology, get_qupath_class

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="bin2cell segmentation-replacement pipeline (batch mode)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── required ──────────────────────────────────────────────────────────
    p.add_argument("--spaceranger-root", required=True, dest="spaceranger_root",
                   help="Root directory whose subdirectories are Space Ranger outputs")
    p.add_argument("--geojson-root", required=True, dest="geojson_root",
                   help="Root directory containing CellViT / QuPath GeoJSON files")
    p.add_argument("--outdir", required=True,
                   help="Root output directory; one subdirectory is created per sample")

    # ── optional ──────────────────────────────────────────────────────────
    p.add_argument("--sample", default=None,
                   help="Process only this sample (default: all samples found)")
    p.add_argument("--geojson-suffix",
                   default="_pgcc_cellvit_multiclass_seg.geojson",
                   dest="geojson_suffix",
                   help="Suffix appended to the sample name to form the GeoJSON filename")
    p.add_argument("--mpp", type=float, default=0.2628,
                   help="Microns per pixel for the rescaled H&E image")
    p.add_argument("--scale", type=float, default=1.0,
                   help="Scale factor applied to GeoJSON polygon coordinates")
    p.add_argument("--buffer", type=int, default=150,
                   help="Pixel buffer used by b2c.scaled_he_image")
    p.add_argument("--min-cells", type=int, default=3,
                   help="sc.pp.filter_genes min_cells threshold")
    p.add_argument("--min-counts", type=int, default=1,
                   help="sc.pp.filter_cells min_counts threshold")
    p.add_argument("--skip-existing", action="store_true", dest="skip_existing",
                   help="Skip a sample if its output directory already contains "
                        "cell_labels.csv (useful for resuming interrupted runs)")
    p.add_argument("--tile-size", type=int, default=200, dest="tile_size",
                   help="Side length (in array_row/col units) of each visualisation tile")
    p.add_argument("--tile-overlap", type=int, default=0, dest="tile_overlap",
                   help="Overlap between adjacent tiles in array_row/col units")
    p.add_argument("--no-viz", action="store_true", dest="no_viz",
                   help="Skip visualisation step (saves time for large datasets)")
    p.add_argument("--viz-only", action="store_true", dest="viz_only",
                   help="Skip all pipeline steps and only (re-)generate visualisations "
                        "from already-saved outputs (requires a completed prior run)")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Sample discovery
# ---------------------------------------------------------------------------

def discover_samples(spaceranger_root, geojson_root, geojson_suffix, target=None):
    """
    Walk ``spaceranger_root`` and return a list of dicts, one per valid sample.

    A directory is considered valid when it contains:
      - sr_wsi.tif
      - binned_outputs/square_002um/
      - a matching GeoJSON in geojson_root named <sample-name><geojson_suffix>
    """
    samples = []
    for name in sorted(os.listdir(spaceranger_root)):
        # honour --sample filter
        if target and name != target:
            continue

        sr_dir = os.path.join(spaceranger_root, name)
        if not os.path.isdir(sr_dir):
            continue

        binned = os.path.join(sr_dir, "binned_outputs", "square_002um")
        if not os.path.isdir(binned):
            log.debug("  %s - no binned_outputs/square_002um/, skipping", name)
            continue

        source_image = os.path.join(sr_dir, "sr_wsi.tif")
        if not os.path.isfile(source_image):
            log.warning("  %s - sr_wsi.tif not found, skipping", name)
            continue

        geojson = os.path.join(geojson_root, name + geojson_suffix)
        if not os.path.isfile(geojson):
            log.warning("  %s - GeoJSON not found: %s, skipping", name, geojson)
            continue

        samples.append(dict(
            sample=name,
            sr_dir=sr_dir,
            binned=binned,
            source_image=source_image,
            sr_image_dir=os.path.join(binned, "spatial"),
            geojson=geojson,
        ))

    return samples


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def build_paths(info, args):
    mpp    = args.mpp
    buf    = args.buffer
    outdir = os.path.join(args.outdir, info["sample"])

    return dict(
        outdir=outdir,
        he_tiff=os.path.join(outdir, "rescaled_image", "he.tiff"),
        npz_path=os.path.join(outdir, f"{info['sample']}_cell_instanseg.npz"),
        xform_json=os.path.join(outdir, "bin2cell_xform_full_to_crop.json"),
        spatial_key=f"spatial_cropped_{buf}_buffer",
        img_key=f"{mpp}_mpp_{buf}_buffer",
        **{k: info[k] for k in ("binned", "source_image", "sr_image_dir", "geojson")},
    )


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def step_load_visium(paths, args):
    log.info("  [1/8] Loading Visium data ...")
    adata = b2c.read_visium(
        paths["binned"],
        source_image_path=paths["source_image"],
        spaceranger_image_path=paths["sr_image_dir"],
    )
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    log.info("  %s", adata)
    return adata


def step_rescale_image(adata, paths, args):
    log.info("  [2/8] Rescaling H&E image ...")
    os.makedirs(os.path.dirname(paths["he_tiff"]), exist_ok=True)
    b2c.scaled_he_image(adata, mpp=args.mpp, save_path=paths["he_tiff"])

    log.info("  [3/8] Computing full->crop coordinate transform ...")
    src = np.asarray(adata.obsm["spatial"], float)
    dst = np.asarray(adata.obsm[paths["spatial_key"]], float)

    ax, b_offset = np.polyfit(src[:, 0], dst[:, 0], 1)
    ay, c_offset = np.polyfit(src[:, 1], dst[:, 1], 1)
    a = float((ax + ay) / 2)

    xf = {"a": a, "b": float(b_offset), "c": float(c_offset)}
    with open(paths["xform_json"], "w") as fh:
        json.dump(xf, fh, indent=2)
    log.info("  Transform saved -> %s", paths["xform_json"])
    return b_offset, c_offset


def step_destripe(adata):
    log.info("  [4/8] Destriping ...")
    b2c.destripe(adata)


def step_rasterise(paths, args, offset_x, offset_y):
    log.info("  [5/8] Rasterising GeoJSON polygons ...")
    with tiff.TiffFile(paths["he_tiff"]) as tf:
        page = tf.pages[0]
        image_shape_hw = (page.imagelength, page.imagewidth)

    label_to_name = geojson_to_sparse_labels_npz(
        paths["geojson"],
        paths["npz_path"],
        tuple(image_shape_hw),
        args.scale,
        offset_x=offset_x,
        offset_y=offset_y,
        return_name_map=True,
    )
    log.info("  %d cell labels rasterised", len(label_to_name))
    return label_to_name


def step_morphology(paths, args, label_to_name):
    log.info("  [6/8] Computing polygon morphology ...")
    morph_df = geojson_polygon_morphology(
        geojson_path=paths["geojson"],
        scale=args.scale,
        mpp=args.mpp,
        name_extractor=get_qupath_class,
    )
    morph_df["name"] = morph_df["label"].map(label_to_name)
    return morph_df.set_index("label")


def step_insert_and_expand_labels(adata, paths, args, label_to_name):
    log.info("  [7/8] Inserting and expanding labels ...")
    b2c.insert_labels(
        adata,
        labels_npz_path=paths["npz_path"],
        basis="spatial",
        spatial_key=paths["spatial_key"],
        mpp=args.mpp,
        labels_key="labels_he",
    )
    b2c.expand_labels(
        adata,
        labels_key="labels_he",
        expanded_labels_key="labels_he_expanded",
    )
    b2c.salvage_secondary_labels(
        adata,
        primary_label="labels_he_expanded",
        secondary_label="labels_he",
        labels_key="labels_joint",
    )
    adata.uns["labels_he_label_to_name"] = label_to_name
    adata.obs["labels_he_name"]    = adata.obs["labels_he"].astype(int).map(label_to_name)
    adata.obs["labels_joint_name"] = adata.obs["labels_joint"].astype(int).map(label_to_name)


def step_attach_bin_morphology(adata, morph_by_label):
    lab = adata.obs["labels_joint"].astype(int)
    for col in ("area_um2", "perimeter_um", "solidity", "circularity", "aspect_ratio"):
        if col in morph_by_label.columns:
            adata.obs[f"cell_{col}"] = lab.map(morph_by_label[col])
    adata.obs["cell_class"] = lab.map(morph_by_label["name"])


def step_bin_to_cell(adata, paths, args, label_to_name, morph_by_label):
    log.info("  [8/8] Aggregating bins -> cells ...")
    cdata = b2c.bin_to_cell(
        adata,
        labels_key="labels_he_expanded",
        spatial_keys=["spatial", paths["spatial_key"]],
    )
    cdata.uns["labels_he_label_to_name"] = label_to_name
    cdata.obs["labels_joint_name"] = (
        cdata.obs["object_id"].astype(int).map(label_to_name)
    )
    cell_ids = pd.to_numeric(pd.Index(cdata.obs_names), errors="coerce").astype("Int64")
    for col in ("area_um2", "circularity", "solidity", "aspect_ratio"):
        if col in morph_by_label.columns:
            cdata.obs[col] = cell_ids.map(morph_by_label[col])
    cdata.obs["class"] = cell_ids.map(morph_by_label["name"])
    log.info("  %s", cdata)
    return cdata



# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

# Colours for the two cell categories
_PGCC_COLOUR    = (1.0, 0.2, 0.2, 0.6)   # semi-transparent red
_NONPGCC_COLOUR = (0.2, 0.6, 1.0, 0.6)   # semi-transparent blue
_BG_COLOUR      = (0.0, 0.0, 0.0, 0.0)   # transparent (background label = 0)


def _is_pgcc(name: str) -> bool:
    """Return True if the cell-class name indicates a PGCC."""
    if not isinstance(name, str):
        return False
    return "pgcc" in name.lower()


def _label_mask_to_rgb(label_crop, label_to_name: dict) -> np.ndarray:
    """
    Colour a 2-D integer label crop as an RGBA image.

    PGCC cells  → red,  non-PGCC cells → blue,  background → transparent.
    """
    H, W = label_crop.shape
    rgba = np.zeros((H, W, 4), dtype=np.float32)

    unique_ids = np.unique(label_crop)
    for lid in unique_ids:
        if lid == 0:
            continue  # background
        name = label_to_name.get(int(lid), "")
        colour = _PGCC_COLOUR if _is_pgcc(name) else _NONPGCC_COLOUR
        mask = label_crop == lid
        rgba[mask] = colour

    return rgba


def _load_viz_inputs(paths):
    """Load H&E image and label mask from disk."""
    with tiff.TiffFile(paths["he_tiff"]) as tf:
        he = tf.pages[0].asarray()
    if he.ndim == 2:
        he = np.stack([he] * 3, axis=-1)   # greyscale -> RGB
    lab_sparse = sparse.load_npz(paths["npz_path"])
    lab = lab_sparse.toarray().astype(np.int32)
    return he, lab


def step_visualise(cdata, paths, args, label_to_name):
    """
    For every non-empty tile produce a side-by-side PNG:
      left panel  - H&E only
      right panel - H&E with PGCC/non-PGCC colour overlay

    Tiles are named: tile_r{r0}-{r1}_c{c0}-{c1}.png
    Only tiles that contain at least one labelled cell are written.
    """
    if args.no_viz:
        log.info("  [viz] Skipped (--no-viz)")
        return

    log.info("  [viz] Generating side-by-side H&E + segmentation tiles ...")

    viz_dir = os.path.join(paths["outdir"], "visualisations")
    os.makedirs(viz_dir, exist_ok=True)

    he, lab = _load_viz_inputs(paths)

    sk = paths["spatial_key"]
    coords_crop = np.asarray(cdata.obsm[sk], float)   # (n_cells, 2): x, y
    arr_row = cdata.obs["array_row"].values.astype(int)
    arr_col = cdata.obs["array_col"].values.astype(int)

    row_min, row_max = arr_row.min(), arr_row.max()
    col_min, col_max = arr_col.min(), arr_col.max()

    ts   = args.tile_size
    step = ts - args.tile_overlap

    legend_patches = [
        mpatches.Patch(color=_PGCC_COLOUR,    label="PGCC"),
        mpatches.Patch(color=_NONPGCC_COLOUR, label="Non-PGCC"),
    ]

    n_saved = 0
    r0 = row_min
    while r0 <= row_max:
        r1 = r0 + ts
        c0 = col_min
        while c0 <= col_max:
            c1 = c0 + ts

            in_tile = (
                (arr_row >= r0) & (arr_row < r1) &
                (arr_col >= c0) & (arr_col < c1)
            )
            if not in_tile.any():
                c0 += step
                continue

            tile_coords = coords_crop[in_tile]
            px_xmin = int(np.floor(tile_coords[:, 0].min()))
            px_xmax = int(np.ceil (tile_coords[:, 0].max()))
            px_ymin = int(np.floor(tile_coords[:, 1].min()))
            px_ymax = int(np.ceil (tile_coords[:, 1].max()))

            pad = 50
            px_xmin = max(0, px_xmin - pad)
            px_ymin = max(0, px_ymin - pad)
            px_xmax = min(he.shape[1], px_xmax + pad)
            px_ymax = min(he.shape[0], px_ymax + pad)

            if px_xmax <= px_xmin or px_ymax <= px_ymin:
                c0 += step
                continue

            he_crop  = he [px_ymin:px_ymax, px_xmin:px_xmax]
            lab_crop = lab[px_ymin:px_ymax, px_xmin:px_xmax]
            overlay  = _label_mask_to_rgb(lab_crop, label_to_name)

            title = f"row {r0}-{r1}  col {c0}-{c1}"

            # side-by-side: left = H&E only, right = H&E + overlay
            fig, (ax_he, ax_ov) = plt.subplots(
                1, 2,
                figsize=(16, 8),
                dpi=150,
                gridspec_kw={"wspace": 0.05},
            )

            ax_he.imshow(he_crop)
            ax_he.set_title("H&E", fontsize=9)
            ax_he.axis("off")

            ax_ov.imshow(he_crop)
            ax_ov.imshow(overlay, interpolation="nearest")
            ax_ov.set_title("H&E + segmentation", fontsize=9)
            ax_ov.axis("off")
            ax_ov.legend(
                handles=legend_patches,
                loc="upper right",
                fontsize=7,
                framealpha=0.7,
            )

            fig.suptitle(title, fontsize=10, y=1.01)

            fname = f"tile_r{r0}-{r1}_c{c0}-{c1}.png"
            fig.savefig(
                os.path.join(viz_dir, fname),
                bbox_inches="tight",
                dpi=150,
            )
            plt.close(fig)
            n_saved += 1

            c0 += step
        r0 += step

    log.info("  [viz] %d tile(s) saved to %s", n_saved, viz_dir)


def run_viz_only(info, args):
    """
    Re-generate visualisations from previously saved outputs without
    re-running any pipeline steps.

    Requires in <outdir>/<sample>/:
      - rescaled_image/he.tiff
      - <sample>_cell_instanseg.npz
      - cell_labels.csv
      - spatial_coords.csv
    """
    paths  = build_paths(info, args)
    outdir = paths["outdir"]

    required = [
        paths["he_tiff"],
        paths["npz_path"],
        os.path.join(outdir, "cell_labels.csv"),
        os.path.join(outdir, "spatial_coords.csv"),
    ]
    missing = [p for p in required if not os.path.isfile(p)]
    if missing:
        raise FileNotFoundError(
            "--viz-only requires prior pipeline outputs. Missing:\n" +
            "\n".join(f"  {m}" for m in missing)
        )

    log.info("  [viz-only] Loading saved outputs ...")

    cell_df = pd.read_csv(os.path.join(outdir, "cell_labels.csv"), index_col=0)
    coords  = pd.read_csv(os.path.join(outdir, "spatial_coords.csv"), index_col=0)

    if "labels_joint_name" not in cell_df.columns:
        raise ValueError("cell_labels.csv must contain a 'labels_joint_name' column")
    if "object_id" not in cell_df.columns:
        raise ValueError("cell_labels.csv must contain an 'object_id' column")

    label_to_name = dict(
        zip(cell_df["object_id"].astype(int),
            cell_df["labels_joint_name"].astype(str))
    )

    class _FakeCdata:
        pass

    cdata = _FakeCdata()
    cdata.obs  = cell_df.copy()
    cdata.obsm = {paths["spatial_key"]: coords[["x", "y"]].values}

    step_visualise(cdata, paths, args, label_to_name)


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

def export_outputs(cdata, outdir):
    log.info("  Exporting outputs ...")

    obs_cols = [
        "object_id", "bin_count", "array_row", "array_col",
        "labels_joint_name", "area_um2", "circularity",
        "solidity", "aspect_ratio", "class",
    ]
    available = [c for c in obs_cols if c in cdata.obs.columns]
    cdata.obs[available].to_csv(os.path.join(outdir, "cell_labels.csv"))

    counts = cdata.X.T if cdata.X.shape[0] == cdata.obs.shape[0] else cdata.X
    io.mmwrite(os.path.join(outdir, "matrix.mtx"), counts)

    pd.Series(cdata.obs_names).to_csv(
        os.path.join(outdir, "barcodes.tsv"), sep="\t", index=False, header=False
    )
    pd.Series(cdata.var_names).to_csv(
        os.path.join(outdir, "features.tsv"), sep="\t", index=False, header=False
    )

    if "spatial" in cdata.obsm:
        coords = pd.DataFrame(
            cdata.obsm["spatial"], index=cdata.obs_names, columns=["x", "y"]
        )
    elif "X_spatial" in cdata.obsm:
        coords = pd.DataFrame(
            cdata.obsm["X_spatial"], index=cdata.obs_names, columns=["x", "y"]
        )
    else:
        raise ValueError("No spatial coordinates found in cdata.obsm")
    coords.to_csv(os.path.join(outdir, "spatial_coords.csv"))

    xy = pd.DataFrame(
        cdata.obsm["spatial"], columns=["x_full", "y_full"], index=cdata.obs_names
    )
    xy.to_csv(os.path.join(outdir, "cdata_spatial_fullres.csv"))

    log.info("  Written to %s", outdir)


# ---------------------------------------------------------------------------
# Per-sample runner
# ---------------------------------------------------------------------------

def run_sample(info, args):
    paths  = build_paths(info, args)
    outdir = paths["outdir"]
    os.makedirs(outdir, exist_ok=True)

    # --viz-only: skip all pipeline steps, just regenerate tiles
    if args.viz_only:
        run_viz_only(info, args)
        return

    if args.skip_existing and os.path.isfile(os.path.join(outdir, "cell_labels.csv")):
        log.info("  Skipping - cell_labels.csv already exists")
        return

    adata = step_load_visium(paths, args)
    b_offset, c_offset = step_rescale_image(adata, paths, args)
    step_destripe(adata)
    label_to_name  = step_rasterise(paths, args, b_offset, c_offset)
    morph_by_label = step_morphology(paths, args, label_to_name)
    step_insert_and_expand_labels(adata, paths, args, label_to_name)
    step_attach_bin_morphology(adata, morph_by_label)
    cdata = step_bin_to_cell(adata, paths, args, label_to_name, morph_by_label)
    step_visualise(cdata, paths, args, label_to_name)
    export_outputs(cdata, outdir)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    samples = discover_samples(
        args.spaceranger_root,
        args.geojson_root,
        args.geojson_suffix,
        target=args.sample,
    )

    if not samples:
        log.error("No valid samples found - check --spaceranger-root and --geojson-root")
        sys.exit(1)

    log.info("Found %d sample(s) to process", len(samples))

    failed = []
    for i, info in enumerate(samples, 1):
        log.info("=" * 60)
        log.info("[%d/%d]  %s", i, len(samples), info["sample"])
        log.info("=" * 60)
        try:
            run_sample(info, args)
            log.info("Done: %s", info["sample"])
        except Exception:
            log.error("FAILED: %s\n%s", info["sample"], traceback.format_exc())
            failed.append(info["sample"])

    log.info("=" * 60)
    log.info("Finished: %d/%d succeeded", len(samples) - len(failed), len(samples))
    if failed:
        log.error("Failed samples: %s", failed)
        sys.exit(1)


if __name__ == "__main__":
    main()
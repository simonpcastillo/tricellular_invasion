"""
Utilities to convert QuPath GeoJSON polygon annotations into a bin2cell-compatible label mask.

bin2cell expects a 2D integer label image saved as a sparse .npz (scipy.sparse.save_npz),
where 0 = background and each object has a unique positive integer label.

This module additionally preserves QuPath object metadata (e.g. class name) by returning
a mapping from label_id -> class/name extracted from GeoJSON properties.
"""

from __future__ import annotations

import json
from typing import Dict, Optional, Tuple

import numpy as np
from shapely.geometry import shape
from skimage.draw import polygon as sk_polygon
from scipy import sparse


def _get_qupath_name(props: dict) -> Optional[str]:
    """
    Extract a human-meaningful label from QuPath GeoJSON feature properties.

    QuPath commonly stores the class label at:
      properties["classification"]["name"]

    Some exports also include:
      properties["class"]

    Returns None if no usable name is found.
    """
    if not isinstance(props, dict):
        return None

    # Preferred: QuPath classification.name
    cls = props.get("classification", {})
    if isinstance(cls, dict):
        name = cls.get("name", None)
        if name is not None and str(name).strip() != "":
            return str(name)

    # Fallback: properties["class"]
    name2 = props.get("class", None)
    if name2 is not None and str(name2).strip() != "":
        return str(name2)

    return None


def geojson_to_sparse_labels_npz(
    geojson_path: str,
    out_npz_path: str,
    image_shape_hw: Tuple[int, int],
    scale: float = 1.0,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    start_label: int = 1,
    return_name_map: bool = False,
) -> Optional[Dict[int, str]]:
    """
    Rasterize QuPath GeoJSON polygons into a 2D integer label mask and save as sparse .npz.

    Parameters
    ----------
    geojson_path:
        Path to a QuPath GeoJSON export (FeatureCollection preferred).
    out_npz_path:
        Output path for the sparse label mask (scipy sparse .npz).
    image_shape_hw:
        (H, W) of the *image used for segmentation / alignment* (e.g. your scaled H&E TIFF).
        This must match the image you later use with b2c.insert_labels().
    scale:
        Multiply GeoJSON coordinates by this factor before rasterization.
        Use 1.0 if GeoJSON coords already match your segmentation image pixel space.
        Example: if GeoJSON is in full-res pixels but you segment a 4x downsampled image, use 0.25.
    start_label:
        Starting integer label id (default 1).
    return_name_map:
        If True, return a dict mapping label_id -> name (QuPath class) extracted from properties.

    Returns
    -------
    Optional[Dict[int, str]]
        label_id -> name mapping if return_name_map=True, otherwise None.
    """
    H, W = image_shape_hw
    lab = np.zeros((H, W), dtype=np.int32)

    with open(geojson_path, "r", encoding="utf-8") as f:
        gj = json.load(f)

    features = gj["features"] if isinstance(gj, dict) and gj.get("type") == "FeatureCollection" else gj
    if not isinstance(features, list):
        raise ValueError("GeoJSON must be a FeatureCollection or a list of Feature objects.")

    current = int(start_label)
    label_to_name: Dict[int, str] = {}

    for feat in features:
        if not isinstance(feat, dict) or "geometry" not in feat:
            continue

        geom = shape(feat["geometry"])
        props = feat.get("properties", {}) or {}

        # Preserve QuPath class/name metadata
        name = _get_qupath_name(props)
        if name is None:
            name = f"obj_{current}"
        label_to_name[current] = str(name)

        # Polygon or MultiPolygon
        polys = [geom] if geom.geom_type == "Polygon" else list(getattr(geom, "geoms", []))
        for poly in polys:
            # exterior coords are (x, y) in pixel space
            xy = np.asarray(poly.exterior.coords, dtype=np.float64)
            #x = xy[:, 0] * scale
            #y = xy[:, 1] * scale
            x = xy[:, 0] * scale + offset_x
            y = xy[:, 1] * scale + offset_y

            rr, cc = sk_polygon(y, x, shape=lab.shape)  # row=y, col=x
            lab[rr, cc] = current

            # clear holes (interiors)
            for hole in poly.interiors:
                hxy = np.asarray(hole.coords, dtype=np.float64)
                hx = hxy[:, 0] * scale
                hy = hxy[:, 1] * scale
                hrr, hcc = sk_polygon(hy, hx, shape=lab.shape)
                lab[hrr, hcc] = 0

        current += 1

    sparse.save_npz(out_npz_path, sparse.csr_matrix(lab))

    if return_name_map:
        return label_to_name
    return None

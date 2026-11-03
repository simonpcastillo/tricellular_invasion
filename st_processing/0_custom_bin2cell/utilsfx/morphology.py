from typing import Optional
import json
import numpy as np
import pandas as pd
from shapely.geometry import shape

def get_qupath_class(props):
    cls = props.get("classification", {})
    if isinstance(cls, dict) and "name" in cls:
        return cls["name"]
    return props.get("class", None)

def geojson_polygon_morphology(
    geojson_path: str,
    scale: float = 1.0,
    mpp: Optional[float] = None,
    name_extractor=None,
):
    """
    Returns a dataframe with one row per GeoJSON Feature (cell instance),
    in the same order used when rasterizing with geojson_to_sparse_labels_npz.

    scale: same scale factor you used for rasterization (GeoJSON coords -> image coords)
    mpp: microns per pixel of the image you will align to (e.g. scaled H&E mpp)
    """
    with open(geojson_path, "r", encoding="utf-8") as f:
        gj = json.load(f)
    features = gj["features"] if isinstance(gj, dict) and gj.get("type") == "FeatureCollection" else gj

    rows = []
    label_id = 1

    for feat in features:
        geom = shape(feat["geometry"])
        props = feat.get("properties", {}) or {}

        # unify: Polygon or MultiPolygon -> treat as one object (union)
        g = geom
        if geom.geom_type == "MultiPolygon":
            # union is robust but sometimes slow; for morphology it's fine
            g = geom.buffer(0)

        # Apply scaling by scaling coordinates indirectly:
        # easiest: multiply derived measures by scale factors:
        # area scales by scale^2, perimeter by scale, bbox dims by scale, centroids by scale.
        area_px2 = g.area * (scale ** 2)
        perim_px = g.length * scale

        minx, miny, maxx, maxy = g.bounds
        bbox_w_px = (maxx - minx) * scale
        bbox_h_px = (maxy - miny) * scale
        aspect_ratio = (bbox_w_px / bbox_h_px) if bbox_h_px > 0 else np.nan

        # centroid
        cx_px = g.centroid.x * scale
        cy_px = g.centroid.y * scale

        # convex hull / solidity
        hull_area_px2 = g.convex_hull.area * (scale ** 2)
        solidity = (area_px2 / hull_area_px2) if hull_area_px2 > 0 else np.nan

        # circularity
        circularity = (4 * np.pi * area_px2 / (perim_px ** 2)) if perim_px > 0 else np.nan

        # optional class/name
        name = None
        if name_extractor is not None:
            name = name_extractor(props)

        row = {
            "label": label_id,
            "name": name,
            "area_px2": area_px2,
            "perimeter_px": perim_px,
            "circularity": circularity,
            "bbox_w_px": bbox_w_px,
            "bbox_h_px": bbox_h_px,
            "aspect_ratio": aspect_ratio,
            "solidity": solidity,
            "cx_px": cx_px,
            "cy_px": cy_px,
        }

        if mpp is not None:
            row.update({
                "area_um2": area_px2 * (mpp ** 2),
                "perimeter_um": perim_px * mpp,
                "bbox_w_um": bbox_w_px * mpp,
                "bbox_h_um": bbox_h_px * mpp,
                "cx_um": cx_px * mpp,
                "cy_um": cy_px * mpp,
            })

        rows.append(row)
        label_id += 1

    return pd.DataFrame(rows)
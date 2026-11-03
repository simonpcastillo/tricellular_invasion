#!/usr/bin/env python3
"""
Build <sample>_tdtomato_rescued.csv -- the post-rescue tdTomato call, in the
2-column form the Neftel pipeline reads from csv_tomato_gated/.

    python 2a4_make_tdtomato_rescued_csv.py --sample 7003_Scan2 \
        --geojson 7003_Scan2_..._rescued_flagged.geojson --outdir csv_tomato_gated

Input is EITHER the exported QuPath GeoJSON (--geojson) or the rescue
pipeline's label table (--labels tumour_labels.csv.gz). Either reproduces the
same CSV; the GeoJSON is what exists for already-processed samples, the label
table is what a fresh pipeline run writes.

WHY THIS IS NOT A ONE-LINER -- three things that silently corrupt the output:

 1. THE RESCUE LIVES IN DIFFERENT FIELDS IN DIFFERENT EXPORTS. On 7000_Scan2
    the 3-way call is in properties.tdtom_label. On 7003_Scan2 tdtom_label
    holds only the ORIGINAL 2-way GMM call and the rescue is in
    properties.classification.name. Reading tdtom_label on 7003 yields the
    un-rescued labels and looks entirely normal -- 13,043 positive instead of
    17,764. classification.name is therefore preferred and tdtom_label is only
    a fallback, with the field actually used reported on stderr.

 2. CASE DIFFERS BETWEEN EXPORTS. 7000 writes "tdtom_pos", 7003 writes
    "TdTom_pos". Any ``== "TdTom_pos"`` test silently fails on one of them.
    Labels are matched case-insensitively and emitted in one canonical casing.

 3. SUSPECT CELLS ARE A POLICY CHOICE, NOT A DETAIL. The specificity audit
    flagged 264 cells on 7003_Scan2 as a possible contaminating non-tumour
    population (a neuronal cluster, confirmed on H&E). The CSV currently in
    csv_tomato_gated/ counts them as POSITIVE. That is reproduced by default
    so this script is a drop-in, but it is an explicit --suspect flag rather
    than an accident, because everything downstream inherits it.

By default a third column `tdtom_class` is written carrying the un-collapsed
call (pos / pos_rescued / pos_rescued_suspect / neg). The 2-column contract is
unchanged -- readers that select cell_id and tdtom_label are unaffected -- but
it means a later question like "is this result driven by rescued cells?" can be
answered without rejoining against the GeoJSON. Pass --no-provenance for a
strictly 2-column file.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
from collections import Counter

# canonical output labels
POS, NEG = "TdTom_pos", "TdTom_neg"
# canonical provenance classes
C_POS, C_RESC, C_SUSP, C_NEG = "pos", "pos_rescued", "pos_rescued_suspect", "neg"

# every spelling seen across exports, lower-cased, -> provenance class
CLASS_MAP = {
    "tdtom_neg": C_NEG,
    "tdtom_pos": C_POS,
    "tdtom_pos_rescued": C_RESC,
    "tdtom_pos_rescued_suspect": C_SUSP,
}


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def iter_features(path: str):
    """Yield feature dicts without holding the parsed document in memory.

    json.load on a 57 MB export expands to well over a gigabyte of Python
    objects. raw_decode walks the features array one element at a time, so only
    the raw text plus a single feature are resident.
    """
    with open(path, "r", errors="replace") as fh:
        raw = fh.read()
    key = '"features"'
    i = raw.find(key)
    if i < 0:
        raise SystemExit("no \"features\" key in %s -- not a FeatureCollection?" % path)
    i = raw.index("[", i) + 1
    dec = json.JSONDecoder()
    n = len(raw)
    while i < n:
        while i < n and raw[i] in ", \t\r\n":
            i += 1
        if i >= n or raw[i] == "]":
            return
        obj, i = dec.raw_decode(raw, i)
        yield obj


def class_of(props: dict, prefer: str) -> tuple[str, str]:
    """Return (provenance_class, field_used) for one feature's properties."""
    cands = []
    cl = props.get("classification")
    if isinstance(cl, dict) and cl.get("name"):
        cands.append(("classification.name", str(cl["name"])))
    if props.get("tdtom_label") is not None:
        cands.append(("tdtom_label", str(props["tdtom_label"])))
    if prefer == "tdtom_label":
        cands.reverse()
    for field, val in cands:
        k = CLASS_MAP.get(val.strip().lower())
        if k is not None:
            return k, field
    seen = ", ".join("%s=%r" % (f, v) for f, v in cands) or "<no label fields>"
    raise SystemExit("unrecognised label (%s). Known: %s"
                     % (seen, ", ".join(sorted(CLASS_MAP))))


def from_geojson(path: str, prefer: str) -> list[tuple[str, str]]:
    rows, fields = [], Counter()
    for feat in iter_features(path):
        p = feat.get("properties") or {}
        cid = p.get("cell_id")
        if cid is None:
            raise SystemExit("feature without cell_id in %s" % path)
        k, field = class_of(p, prefer)
        fields[field] += 1
        rows.append((str(cid), k))
    log("  label field used: %s" % ", ".join("%s (%d)" % kv for kv in fields.items()))
    return rows


def from_labels(path: str) -> list[tuple[str, str]]:
    """Read the rescue pipeline's tumour_labels.csv.gz."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        rdr = csv.DictReader(fh)
        cols = set(rdr.fieldnames or [])
        need = {"cell", "tumour_call"}
        if not need <= cols:
            raise SystemExit("%s needs columns %s; has %s"
                             % (path, sorted(need), sorted(cols)))
        has_susp = "suspect" in cols
        rows = []
        for r in rdr:
            tc = (r["tumour_call"] or "").strip().lower()
            if tc == "probe_positive":
                k = C_POS
            elif tc == "rescued":
                k = C_SUSP if (has_susp and str(r["suspect"]).lower()
                               in ("true", "1")) else C_RESC
            elif tc == "unlabeled":
                k = C_NEG
            else:
                raise SystemExit(
                    f"unrecognised tumour_call {r['tumour_call']!r} "
                    f"for cell {r['cell']!r} in {path}"
                )
            if not str(r["cell"]).strip():
                raise SystemExit(f"blank cell ID in {path}")
            rows.append((str(r["cell"]), k))
    if not has_susp:
        log("  note: no `suspect` column -- no cell can be classed "
            "pos_rescued_suspect")
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--geojson", help="exported QuPath GeoJSON")
    src.add_argument("--labels", help="tumour_labels.csv.gz from the rescue pipeline")
    ap.add_argument("--sample", required=True, help="e.g. 7003_Scan2")
    ap.add_argument("--outdir", default=".", help="written as <outdir>/<sample>_tdtomato_rescued.csv")
    ap.add_argument("--suspect", choices=("positive", "negative", "exclude"),
                    default="positive",
                    help="audit-flagged rescued cells: count as positive "
                         "(default, reproduces the existing file), as negative, "
                         "or drop the rows entirely")
    ap.add_argument("--no-provenance", action="store_true",
                    help="omit the tdtom_class column (strict 2-column output)")
    ap.add_argument("--provenance-sidecar", action="store_true",
                    help="write a separate <sample>_tdtomato_provenance.csv with "
                         "cell_id and original provenance class, even when the "
                         "main CSV has only two columns")
    ap.add_argument("--prefer-field", choices=("classification", "tdtom_label"),
                    default="classification",
                    help="which GeoJSON field to trust first [classification]")
    ap.add_argument("--expect-total", type=int, default=None,
                    help="assert this many cells; fails loudly on a partial read")
    a = ap.parse_args(argv)

    log("sample %s" % a.sample)
    rows = (from_geojson(a.geojson, a.prefer_field) if a.geojson
            else from_labels(a.labels))

    counts = Counter(k for _, k in rows)
    log("  classes: " + " | ".join("%s %d" % (k, counts[k])
                                   for k in (C_POS, C_RESC, C_SUSP, C_NEG) if counts[k]))

    if a.expect_total is not None and len(rows) != a.expect_total:
        raise SystemExit("expected %d cells, read %d" % (a.expect_total, len(rows)))
    if not rows:
        raise SystemExit("no cells read")

    # ids must be unique or a downstream join silently multiplies rows
    ids = [c for c, _ in rows]
    if len(set(ids)) != len(ids):
        dup = [c for c, n in Counter(ids).items() if n > 1][:5]
        raise SystemExit("duplicate cell_id (e.g. %s)" % dup)

    out_rows = []
    dropped = 0
    for cid, k in rows:
        if k == C_SUSP and a.suspect == "exclude":
            dropped += 1
            continue
        if k == C_POS or k == C_RESC:
            lab = POS
        elif k == C_SUSP:
            lab = POS if a.suspect == "positive" else NEG
        else:
            lab = NEG
        out_rows.append((cid, lab, k))

    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, "%s_tdtomato_rescued.csv" % a.sample)
    header = ["cell_id", "tdtom_label"] + ([] if a.no_provenance else ["tdtom_class"])
    # lineterminator="\n": csv.writer defaults to \r\n, and tdtom_label is the
    # LAST column -- a trailing \r makes the label compare unequal to
    # "TdTom_pos" in any reader that does not strip it.
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        for cid, lab, k in out_rows:
            w.writerow([cid, lab] + ([] if a.no_provenance else [k]))

    if a.provenance_sidecar:
        prov = os.path.join(a.outdir, "%s_tdtomato_provenance.csv" % a.sample)
        with open(prov, "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["cell_id", "tdtom_class"])
            for cid, _lab, cls in out_rows:
                w.writerow([cid, cls])
        log("  provenance sidecar: %s" % prov)

    final = Counter(lab for _, lab, _ in out_rows)
    log("  suspect policy: %s%s" % (a.suspect,
                                    " (%d rows dropped)" % dropped if dropped else ""))
    log("  wrote %s" % out)
    log("  %s %d | %s %d | total %d"
        % (NEG, final[NEG], POS, final[POS], len(out_rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

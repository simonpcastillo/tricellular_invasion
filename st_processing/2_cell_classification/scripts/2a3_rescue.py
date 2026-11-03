#!/usr/bin/env python3
"""Rescue reporter-positive cells the probe failed to capture. Sample-agnostic.

Method
------
  1. Labels are positive-unlabeled, NOT positive-negative. A cell with reporter
     signal is a true positive; a cell without it is UNKNOWN, not negative.
  2. Features are SCT principal components with sequencing depth projected out
     linearly and quadratically. The reporter gene is excluded upstream (stage 1).
  3. A depth-detection floor sets HOW MANY cells to rescue: the detection rate in
     the deepest depth decile, applied sample-wide. This is a LOWER BOUND, not a
     point estimate -- the depth curve does not saturate, so no asymptote is
     identifiable.
  4. Validation is stratified on EXACT integer depth, never on decile bins.
     Approximate depth matching failed its own null control during development.
  5. A specificity audit flags rescued cells that look like a contaminating
     non-tumour population (confounder signature, region position, cluster
     uplift). See --help-audit.

Guards that exist because each of these went wrong at least once:
  - reporter must be absent from the feature matrix (asserted)
  - detection is "count > 0", never "round(count) > 0"
  - cell identifiers are compared as strings (an int64 join silently matched
    nothing and produced an all-negative rescue set)
  - the permuted-label null control must return ~0.5 or the run ABORTS
  - depth strata below --min-stratum cells are dropped from z-scoring
    (singleton strata have zero variance and silently produce NaN)

Usage:
  python 2a3_rescue.py --indir DIR --outdir DIR [options]

Outputs (all in --outdir):
  rescue_ranking.csv.gz     per-unlabeled-cell score, rank, rescued flag
  tumour_labels.csv.gz      per-cell tumour_call: probe_positive/rescued/unlabeled
  suspect_cells.csv         rescued cells flagged by the specificity audit
  holdout_validation.csv    precision/recall/stratified AUC per replicate
  auc_by_depth_decile.csv   expression AUC vs depth-alone AUC, per decile
  prevalence_by_cluster.csv observed and corrected prevalence per cluster
  spatial_coherence.csv     neighbour enrichment, prevalence-normalised
  summary.json              every headline metric, machine-readable
  run_provenance.txt        arguments, versions, input checksums
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

AUDIT_HELP = """
Specificity audit -- what it is and why it is on by default
-----------------------------------------------------------
A PU classifier trained on a sparse reporter can rescue a coherent non-tumour
population whose expression happens to resemble the positives. On sample 7003 a
214-cell neuronal cluster went from 0.9% observed prevalence to 59.8% after
rescue (64x uplift), sat 92% outside the histological tumour boundary, and was
confirmed neuronal on H&E. It was 2.7% of rescues but 47% of the rescues falling
outside the boundary.

Three independent flags, any of which marks a rescued cell 'suspect':
  1. confounder_score  -- a signature score from stage 1, z-scored WITHIN exact
                          depth strata, flagged above --suspect-z.
  2. in_region         -- rescued cells outside a histology region, when a
                          region column was exported.
  3. cluster uplift    -- membership of a cluster whose prevalence uplift exceeds
                          --max-cluster-uplift while its observed prevalence is
                          below --min-cluster-prev.

Suspect cells are FLAGGED, never dropped: they are written to suspect_cells.csv
and marked in tumour_labels.csv so you decide. The run WARNS (does not abort) if
the suspect fraction exceeds --max-suspect-frac.
"""


# ----------------------------------------------------------------- utilities
def log(msg: str) -> None:
    print(msg, flush=True)


def sha256(path: Path, blocks: int = 4096) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(blocks * 1024), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def strat_auc(score, target, depth):
    """AUC computed only between cells of IDENTICAL integer depth.

    Pools Mann-Whitney comparisons within each exact-depth stratum. Any
    depth-driven component of `score` cannot contribute, because every
    comparison holds depth fixed. Returns (auc, n_strata, n_pairs).
    """
    df = pd.DataFrame({"d": np.asarray(depth).astype(np.int64),
                       "s": np.asarray(score, float),
                       "t": np.asarray(target).astype(int)})
    num = den = 0.0
    strata = 0
    for _, g in df.groupby("d", sort=False):
        p = g.s.to_numpy()[g.t.to_numpy() == 1]
        n = g.s.to_numpy()[g.t.to_numpy() == 0]
        if len(p) == 0 or len(n) == 0:
            continue
        strata += 1
        diff = p[:, None] - n[None, :]
        num += (diff > 0).sum() + 0.5 * (diff == 0).sum()
        den += p.size * n.size
    return (num / den if den else np.nan), strata, int(den)


def project_out(X, covars):
    """Residualise X on covars (with intercept), then z-score."""
    D = np.column_stack([covars, np.ones(len(X))])
    beta, *_ = np.linalg.lstsq(D, X, rcond=None)
    R = X - D @ beta
    return ((R - R.mean(0)) / (R.std(0) + 1e-8)).astype(np.float32)


def depth_floor(detected, depth, n_bins=10):
    """Lower bound on true prevalence from the deepest depth bin.

    If detection still rises with depth at the deepest bin, that bin's rate is
    the best-observed capture and every shallower cell is under-sampled relative
    to it. Applying that rate sample-wide is therefore a FLOOR. Do not
    extrapolate a saturating asymptote -- the curve is unsaturated and the
    asymptote is not identifiable.
    """
    q = pd.qcut(depth, n_bins, labels=False, duplicates="drop")
    rate = np.array([detected[q == b].mean() for b in np.unique(q)])
    deepest = float(rate[-1])
    n_deep = int((q == np.unique(q)[-1]).sum())
    se = np.sqrt(deepest * (1 - deepest) / n_deep)
    monotone = bool(np.all(np.diff(rate[-4:]) > -0.005))
    return dict(floor=deepest, floor_lo=max(0.0, deepest - 1.96 * se),
                floor_hi=min(1.0, deepest + 1.96 * se),
                observed=float(detected.mean()), rate_by_bin=rate.tolist(),
                still_rising_at_deepest=monotone)


def _z_by_key(v, key, min_stratum):
    s = pd.DataFrame({"k": key, "v": v}).groupby("k").v
    mu, sd, n = (s.transform(f).to_numpy() for f in ("mean", "std", "size"))
    usable = (n >= min_stratum) & (sd > 1e-8) & np.isfinite(v)
    z = np.full(len(v), np.nan)
    z[usable] = (v[usable] - mu[usable]) / sd[usable]
    return z, usable


def z_within_depth(values, depth, min_stratum=20, min_coverage=0.80, n_bins=40):
    """z-score `values` within depth strata, adaptively.

    Exact integer depth is preferred: it is the strictest possible depth control.
    But exact strata only hold enough cells when the sample is large relative to
    the depth range -- on ~100k cells they cover >99% of cells, on a few thousand
    they cover almost none, because most depth values are near-unique. Strata
    below `min_stratum` cells or with zero variance MUST be dropped (a singleton
    stratum has zero SD and silently yields NaN, which propagates through every
    median downstream).

    So: try exact depth; if it covers less than `min_coverage` of cells, fall back
    to quantile bins of depth. The fallback is a weaker control but still removes
    the monotone depth trend, and it is reported in the returned mode so the
    summary records which was used.

    Returns (z, usable_mask, mode).
    """
    v = np.asarray(values, float)
    d = np.asarray(depth).astype(np.int64)
    z, usable = _z_by_key(v, d, min_stratum)
    if usable.mean() >= min_coverage:
        return z, usable, "exact_depth"
    q = pd.qcut(d.astype(float), n_bins, labels=False, duplicates="drop")
    z2, usable2 = _z_by_key(v, q, min_stratum)
    if usable2.mean() <= usable.mean():
        return z, usable, "exact_depth"
    return z2, usable2, f"depth_quantile_bins_{int(pd.Series(q).nunique())}"


# ---------------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # not marked required= so that --help-audit works on its own; checked below
    ap.add_argument("--indir", type=Path,
                    help="stage-1 output: embedding/meta/reporter .csv.gz")
    ap.add_argument("--outdir", type=Path)
    ap.add_argument("--sample", default=None,
                    help="sample label recorded in summary.json and prefixed to "
                         "output filenames; defaults to the indir parent name")
    ap.add_argument("--reporter-col", default="reporter_count")
    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--n-folds", type=int, default=5)
    ap.add_argument("--n-reps", type=int, default=5)
    ap.add_argument("--holdout-frac", type=float, default=0.2)
    ap.add_argument("--C", type=float, default=0.5,
                    help="inverse L2 strength for the logistic PU model")
    ap.add_argument("--seed", type=int, default=0)
    # --- specificity audit
    ap.add_argument("--suspect-z", type=float, default=2.0,
                    help="depth-controlled confounder z above which a rescued "
                         "cell is flagged suspect")
    ap.add_argument("--min-stratum", type=int, default=20,
                    help="minimum cells in a depth stratum for z-scoring")
    ap.add_argument("--max-cluster-uplift", type=float, default=5.0,
                    help="cluster prevalence uplift above which the cluster is "
                         "flagged, if its observed prevalence is also low")
    ap.add_argument("--min-cluster-prev", type=float, default=0.05,
                    help="observed prevalence below which a high-uplift cluster "
                         "is treated as suspect rather than genuinely enriched")
    ap.add_argument("--max-suspect-frac", type=float, default=0.15,
                    help="warn if more than this fraction of rescues are suspect")
    ap.add_argument("--null-tol", type=float, default=0.04,
                    help="permitted deviation of the null control from 0.5")
    ap.add_argument("--no-audit", action="store_true",
                    help="skip the specificity audit entirely (not recommended)")
    ap.add_argument("--help-audit", action="store_true",
                    help="explain the specificity audit and exit")
    a = ap.parse_args(argv)

    if a.help_audit:
        print(AUDIT_HELP)
        return 0
    for req in ("indir", "outdir"):
        if getattr(a, req) is None:
            ap.error(f"--{req} is required")

    a.outdir.mkdir(parents=True, exist_ok=True)
    sample = a.sample or a.indir.resolve().parent.name
    log(f"sample: {sample}")

    # -------------------------------------------------- load, join on strings
    paths = {n: a.indir / f"{n}.csv.gz" for n in ("embedding", "meta", "reporter")}
    for n, p in paths.items():
        if not p.exists():
            raise SystemExit(f"missing stage-1 output: {p}")
    emb = pd.read_csv(paths["embedding"])
    meta = pd.read_csv(paths["meta"])
    rep = pd.read_csv(paths["reporter"])
    for d in (emb, meta, rep):
        d["cell"] = d["cell"].astype(str)   # int64 join silently matches nothing

    df = meta.merge(rep, on="cell", validate="one_to_one")
    df = df.merge(emb, on="cell", validate="one_to_one")
    log(f"cells: {len(df):,}")
    if len(df) != len(meta):
        raise SystemExit(f"join lost cells: {len(meta):,} -> {len(df):,}")

    pcs = [c for c in emb.columns if c != "cell"][: a.n_pcs]
    X = df[pcs].to_numpy(np.float32)
    if X.shape[1] < a.n_pcs:
        log(f"warning: only {X.shape[1]} PCs available, requested {a.n_pcs}")

    # reporter must not be a feature
    assert not any(a.reporter_col.lower() in c.lower() for c in pcs), \
        "reporter column leaked into the feature matrix"

    # detection = ANY signal, never rounded
    if a.reporter_col not in df.columns:
        raise SystemExit(f"--reporter-col '{a.reporter_col}' not in the joined "
                         f"table; available: {sorted(df.columns)[:20]}")
    y = (df[a.reporter_col].fillna(0).to_numpy(float) > 0)
    log(f"reporter+ (any signal): {y.sum():,} ({y.mean():.2%})")
    if y.sum() < 100:
        raise SystemExit(f"only {int(y.sum())} reporter-positive cells; too few "
                         "to fit a PU model. Check the reporter gene name.")

    dep_col = next((c for c in ("nFeature_Spatial", "nFeature_SCT", "nFeature_RNA")
                    if c in df.columns), None)
    cnt_col = next((c for c in ("nCount_Spatial", "nCount_RNA") if c in df.columns), None)
    if dep_col is None:
        raise SystemExit("no nFeature_* column found; cannot control for depth")
    dep = df[dep_col].to_numpy(float)
    log(f"depth column: {dep_col} (median {np.median(dep):.0f})")

    # ------------------------------------------------------ prevalence floor
    fl = depth_floor(y, dep)
    n_total = len(df)
    fl["cells_floor"] = int(round(fl["floor"] * n_total))
    fl["cells_observed"] = int(y.sum())
    fl["n_to_rescue"] = max(0, fl["cells_floor"] - fl["cells_observed"])
    log(f"observed {fl['observed']:.2%} | depth floor {fl['floor']:.2%} "
        f"(95% CI {fl['floor_lo']:.2%}-{fl['floor_hi']:.2%}) "
        f"-> rescue {fl['n_to_rescue']:,} cells")
    if fl["n_to_rescue"] == 0:
        log("NOTE: the depth floor is at or below observed prevalence. There is "
            "nothing to rescue on this sample -- probe capture is already at the "
            "deepest-decile rate. This is a valid result, not an error.")
    if not fl["still_rising_at_deepest"]:
        log("NOTE: detection has flattened at the deepest bins. The floor is then "
            "close to a point estimate rather than a lower bound.")

    # -------------------------------------------------- depth-neutral features
    cov = [np.log1p(dep), np.log1p(dep) ** 2]
    if cnt_col is not None:
        c2 = df[cnt_col].to_numpy(float)
        cov += [np.log1p(c2), np.log1p(c2) ** 2]
    Xr = project_out(X, np.column_stack(cov))
    resid_r = abs(np.corrcoef(Xr[:, 0], np.log1p(dep))[0, 1])
    log(f"depth residual in PC1 after projection: r = {resid_r:.4f}")

    # ------------------------------------------- out-of-fold scores + per-decile
    def fit(Xm, ym):
        return LogisticRegression(max_iter=4000, C=a.C).fit(Xm, ym.astype(int))

    oof = np.zeros(len(y))
    for tr, te in StratifiedKFold(a.n_folds, shuffle=True,
                                  random_state=a.seed).split(Xr, y):
        oof[te] = fit(Xr[tr], y[tr]).predict_proba(Xr[te])[:, 1]
    auc_pool = roc_auc_score(y, oof)
    auc_depth_only = roc_auc_score(y, dep)
    log(f"out-of-fold AUC {auc_pool:.4f} | depth alone {auc_depth_only:.4f}")
    if auc_pool <= auc_depth_only:
        log("WARNING: pooled AUC does not exceed depth alone. At the pooled level "
            "the expression model adds nothing beyond sequencing depth. The "
            "per-decile and exact-depth-stratified numbers below are the ones "
            "that matter -- read them before using this rescue.")

    dq = pd.qcut(dep, 10, labels=False, duplicates="drop")
    dec = pd.DataFrame([
        dict(decile=int(b) + 1, n=int((dq == b).sum()),
             npos=int(y[dq == b].sum()),
             auc_expr=roc_auc_score(y[dq == b], oof[dq == b]),
             auc_depth=roc_auc_score(y[dq == b], dep[dq == b]),
             med_depth=float(np.median(dep[dq == b])))
        for b in np.unique(dq) if 0 < y[dq == b].sum() < (dq == b).sum()])
    dec.to_csv(a.outdir / "auc_by_depth_decile.csv", index=False)
    log(f"per-decile expression AUC: {dec.auc_expr.min():.3f}-{dec.auc_expr.max():.3f}"
        f" | depth alone: {dec.auc_depth.min():.3f}-{dec.auc_depth.max():.3f}")

    # ------------------------------------------- hold-out validation + null
    rng = np.random.default_rng(a.seed)
    pos = np.where(y)[0]
    ks = [fl["n_to_rescue"] // 4, fl["n_to_rescue"] // 2, fl["n_to_rescue"]]
    ks = sorted({k for k in ks if k > 0})
    rows, sa_list = [], []
    for rep_i in range(a.n_reps):
        hid = rng.choice(pos, int(len(pos) * a.holdout_frac), replace=False)
        ytr = y.copy()
        ytr[hid] = False
        s = fit(Xr, ytr).predict_proba(Xr)[:, 1]
        c = np.where(~ytr)[0]
        tgt = np.isin(c, hid)
        sa, nst, npr = strat_auc(s[c], tgt, dep[c])
        sa_list.append(sa)
        order = c[np.argsort(-s[c])]
        for k in ks:
            t = order[:k]
            rows.append(dict(rep=rep_i, k=k, precision=float(np.isin(t, hid).mean()),
                             recall=float(np.isin(hid, t).mean()),
                             baseline=float(tgt.mean()), strat_auc=sa,
                             n_strata=nst, n_pairs=npr))
    ho = pd.DataFrame(rows)
    ho.to_csv(a.outdir / "holdout_validation.csv", index=False)
    summ = ho.groupby("k")[["precision", "recall", "baseline"]].mean()
    summ["enrichment"] = summ.precision / summ.baseline
    log("hold-out (mean of reps):")
    for k, r in summ.iterrows():
        log(f"  top {int(k):>7,}: precision {r.precision:.3%} "
            f"recall {r.recall:.2%} enrichment {r.enrichment:.2f}x")
    log(f"depth-stratified AUC: {np.mean(sa_list):.4f} "
        f"(n={len(sa_list)} reps, {rows[-1]['n_strata']} strata)")

    # NULL CONTROL: permuted labels must give ~0.5 under the same machinery.
    # A stratification scheme that fails this is biased; abort rather than ship.
    yp = rng.permutation(y)
    sp = fit(Xr, yp).predict_proba(Xr)[:, 1]
    cp = np.where(~yp)[0]
    null_auc, _, _ = strat_auc(sp[cp], np.isin(cp, np.where(y)[0]), dep[cp])
    log(f"null control (permuted labels): {null_auc:.4f}")
    if not (0.5 - a.null_tol <= null_auc <= 0.5 + a.null_tol):
        raise SystemExit(
            f"ABORT: null control {null_auc:.4f} is not ~0.5 (tolerance "
            f"+/-{a.null_tol}). The depth stratification is biased; do not trust "
            "the reported AUCs.")

    # ------------------------------------------------------- final ranking
    final = fit(Xr, y).predict_proba(Xr)[:, 1]
    unl = np.where(~y)[0]
    ranks = np.argsort(np.argsort(-final[unl])) + 1
    out = pd.DataFrame(dict(cell=df.cell.to_numpy()[unl], score=final[unl],
                            rank=ranks, depth=dep[unl]))
    out["rescued"] = out["rank"] <= fl["n_to_rescue"]
    out = out.sort_values("rank")
    out.to_csv(a.outdir / "rescue_ranking.csv.gz", index=False)
    log(f"ranked {len(out):,} unlabeled cells; flagged {int(out.rescued.sum()):,}")

    # combined label column for downstream use
    lab = pd.DataFrame(dict(cell=df.cell, reporter_positive=y))
    lab = lab.merge(out[["cell", "score", "rank", "rescued"]], on="cell", how="left")
    lab["rescued"] = lab.rescued.astype(object).where(lab.rescued.notna(), False).astype(bool)
    lab["tumour_call"] = np.where(lab.reporter_positive, "probe_positive",
                                  np.where(lab.rescued, "rescued", "unlabeled"))
    resc_mask = lab.rescued.to_numpy()

    # ------------------------------------------- per-cluster prevalence + uplift
    pc, r_oc, r_oc_excl, flag_clusters = None, None, None, []
    if "seurat_clusters" in df.columns:
        cl = df.seurat_clusters.astype(str).to_numpy()
        corr = lab.tumour_call.isin(["probe_positive", "rescued"]).to_numpy()
        pc = pd.DataFrame({
            "n": pd.Series(cl).value_counts(),
            "observed": pd.Series(y).groupby(cl).mean(),
            "corrected": pd.Series(corr).groupby(cl).mean()})
        pc["uplift"] = pc.corrected / pc.observed.replace(0, np.nan)
        # A cluster with ZERO observed prevalence but non-zero corrected prevalence
        # has infinite uplift and is the MOST suspicious case there is -- the model
        # populated it from nothing. Treat NaN uplift with corrected > 0 as
        # infinite, or it silently escapes the comparison below.
        pc.loc[(pc.observed == 0) & (pc.corrected > 0), "uplift"] = np.inf
        pc["suspect_cluster"] = ((pc.uplift > a.max_cluster_uplift) &
                                 (pc.observed < a.min_cluster_prev)).fillna(False)
        flag_clusters = list(pc.index[pc.suspect_cluster])
        r_oc = float(np.corrcoef(pc.observed, pc.corrected)[0, 1])
        log(f"cluster prevalence corr(observed, corrected) = {r_oc:.3f}")
        if flag_clusters:
            keep = ~pc.index.isin(flag_clusters)
            if keep.sum() >= 3:
                r_oc_excl = float(np.corrcoef(pc.observed[keep], pc.corrected[keep])[0, 1])
            log(f"HIGH-UPLIFT CLUSTERS: {', '.join(map(str, flag_clusters))} "
                f"(uplift > {a.max_cluster_uplift}x with observed prevalence < "
                f"{a.min_cluster_prev:.0%})")
            if r_oc_excl is not None:
                log(f"  corr excluding them = {r_oc_excl:.3f} "
                    f"(was {r_oc:.3f}) -- a large jump means those clusters carry "
                    "the whole discrepancy")
        pc.sort_values("corrected", ascending=False).to_csv(
            a.outdir / "prevalence_by_cluster.csv")

    # ------------------------------------------------- specificity audit
    audit = dict(enabled=not a.no_audit)
    if not a.no_audit and resc_mask.sum() > 0:
        flags = pd.DataFrame(dict(cell=df.cell.to_numpy()))
        crit = []

        # 1. confounder signature, z-scored within exact depth strata
        if "confounder_score" in df.columns and df.confounder_score.notna().any():
            zc, usable, zmode = z_within_depth(df.confounder_score.to_numpy(), dep,
                                               a.min_stratum)
            flags["confounder_z"] = zc
            flags["flag_confounder"] = usable & (zc > a.suspect_z)
            crit.append("flag_confounder")
            m_r = resc_mask & usable
            audit["confounder_z_mode"] = zmode
            audit["confounder_z_coverage"] = float(usable.mean())
            log(f"confounder z (mode: {zmode}, {usable.sum():,} usable cells = "
                f"{usable.mean():.1%}):")
            if zmode != "exact_depth":
                log("  NOTE: exact-depth strata were too sparse on this sample, so "
                    "depth quantile bins were used. This is a weaker depth control; "
                    "it is recorded in summary.json.")
            for nm, m_ in [("reporter+", y & usable), ("rescued", m_r),
                           ("unlabeled", (~y) & (~resc_mask) & usable)]:
                if m_.sum():
                    log(f"  {nm:<12} n={m_.sum():>7,} median z {np.median(zc[m_]):+.3f} "
                        f"| z>{a.suspect_z:g}: {(zc[m_] > a.suspect_z).mean():.1%}")
            audit["confounder_median_z_rescued"] = float(np.median(zc[m_r])) if m_r.sum() else None
            audit["confounder_median_z_reporter_pos"] = (
                float(np.median(zc[y & usable])) if (y & usable).sum() else None)
        else:
            log("no confounder_score in meta -- signature flag skipped. Pass a "
                "signature file to stage 1 to enable it (see --help-audit).")

        # 2. outside a histology region
        if "in_region" in df.columns:
            inr = df.in_region.fillna(False).to_numpy().astype(bool)
            flags["flag_outside_region"] = ~inr
            crit.append("flag_outside_region")
            log(f"region: {inr.mean():.1%} of all cells inside | rescued inside "
                f"{inr[resc_mask].mean():.1%}")
            audit["frac_all_inside_region"] = float(inr.mean())
            audit["frac_rescued_inside_region"] = float(inr[resc_mask].mean())

        # 3. membership of a high-uplift, low-prevalence cluster
        if flag_clusters:
            cl = df.seurat_clusters.astype(str).to_numpy()
            flags["flag_suspect_cluster"] = np.isin(cl, flag_clusters)
            crit.append("flag_suspect_cluster")
            audit["suspect_clusters"] = [str(c) for c in flag_clusters]

        if crit:
            flags["suspect"] = flags[crit].any(axis=1) & resc_mask
            n_s = int(flags.suspect.sum())
            frac = n_s / max(1, int(resc_mask.sum()))
            log(f"SPECIFICITY AUDIT: {n_s:,} of {int(resc_mask.sum()):,} rescued "
                f"cells flagged suspect ({frac:.1%})")
            for c in crit:
                log(f"  {c:<22} {int((flags[c] & resc_mask).sum()):>6,}")
            if frac > a.max_suspect_frac:
                log(f"WARNING: suspect fraction {frac:.1%} exceeds "
                    f"{a.max_suspect_frac:.0%}. Inspect suspect_cells.csv against "
                    "histology before using the rescued set. Consider re-running "
                    "with the confounding population removed upstream.")
            sus = flags[flags.suspect].copy()
            for c in ("x", "y", "seurat_clusters"):
                if c in df.columns:
                    sus[c] = df.set_index("cell").loc[sus.cell, c].to_numpy()
            sus = sus.merge(out[["cell", "score", "rank"]], on="cell", how="left")
            sus.sort_values("rank").to_csv(a.outdir / "suspect_cells.csv", index=False)
            lab = lab.merge(flags[["cell", "suspect"] + crit], on="cell", how="left")
            lab["suspect"] = lab.suspect.fillna(False).astype(bool)
            lab["tumour_call_audited"] = np.where(
                lab.suspect, "rescued_suspect", lab.tumour_call)
            audit.update(n_suspect=n_s, frac_suspect=float(frac), criteria=crit)
        else:
            log("specificity audit: no criteria available (no signature, no "
                "region column, no high-uplift clusters). Nothing flagged.")
            audit["criteria"] = []

    lab.to_csv(a.outdir / "tumour_labels.csv.gz", index=False)

    # ------------------------------------------------- spatial QC (not a feature)
    sc = None
    if {"x", "y"}.issubset(df.columns):
        from sklearn.neighbors import NearestNeighbors
        xy = df[["x", "y"]].to_numpy(float)
        ok = np.isfinite(xy).all(1)
        k_nn = min(31, int(ok.sum()))
        _, idx = NearestNeighbors(n_neighbors=k_nn).fit(xy[ok]).kneighbors(xy[ok])
        idx = idx[:, 1:]

        def coherence(mask):
            """Neighbour enrichment over the group's OWN prevalence.

            Prevalence-normalised so groups of different size are comparable;
            a raw neighbour fraction is not.
            """
            m = mask[ok]
            if m.sum() == 0:
                return np.nan
            return float(m[idx].mean(1)[m].mean() / m.mean())

        half = np.zeros(len(df), bool)
        if fl["n_to_rescue"] > 1:
            half[unl[np.argsort(-final[unl])[: fl["n_to_rescue"] // 2]]] = True
        sc = pd.DataFrame([
            dict(group="probe_positive", n=int(y.sum()), enrichment=coherence(y)),
            dict(group="rescued_top_half", n=int(half.sum()), enrichment=coherence(half)),
            dict(group="rescued_all", n=int(resc_mask.sum()),
                 enrichment=coherence(resc_mask)),
            dict(group="random_matched", n=int(y.sum()),
                 enrichment=coherence(rng.permutation(y))),
        ])
        sc.to_csv(a.outdir / "spatial_coherence.csv", index=False)
        log("spatial coherence (x own prevalence): " +
            ", ".join(f"{r.group} {r.enrichment:.2f}" for r in sc.itertuples()))

    # ------------------------------------------------------------ manifest
    summary = dict(
        sample=sample, n_cells=int(n_total), reporter_col=a.reporter_col,
        depth_col=dep_col, n_pcs=int(X.shape[1]), prevalence=fl,
        auc_out_of_fold=float(auc_pool), auc_depth_alone=float(auc_depth_only),
        auc_per_decile_min=float(dec.auc_expr.min()),
        auc_per_decile_max=float(dec.auc_expr.max()),
        strat_auc_holdout=float(np.mean(sa_list)),
        null_control_auc=float(null_auc),
        depth_residual_pc1=float(resid_r),
        cluster_corr_observed_corrected=r_oc,
        cluster_corr_excl_suspect=r_oc_excl,
        n_rescued=int(resc_mask.sum()),
        audit=audit,
        holdout=summ.reset_index().to_dict("records"),
        seed=a.seed,
        params={k: (str(v) if isinstance(v, Path) else v)
                for k, v in vars(a).items()})
    (a.outdir / "summary.json").write_text(json.dumps(summary, indent=2))

    (a.outdir / "run_provenance.txt").write_text("\n".join([
        f"sample={sample}",
        f"indir={a.indir.resolve()}",
        *(f"sha256_{n}={sha256(p)}" for n, p in paths.items()),
        f"python={platform.python_version()}",
        f"numpy={np.__version__}", f"pandas={pd.__version__}",
        f"sklearn={sklearn.__version__}",
        f"argv={' '.join(sys.argv[1:])}",
    ]) + "\n")

    log(f"wrote outputs to {a.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

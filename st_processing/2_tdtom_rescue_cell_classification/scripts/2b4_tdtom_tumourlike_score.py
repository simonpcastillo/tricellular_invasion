#!/usr/bin/env python3
# =============================================================================
# 2b4_tdtom_tumourlike_score.py
#
# Score clusters AND sub-clusters on transcriptional tumour-likeness, from the
# marker tables written by 2b1_tdtom_cluster_enrichment_v4.R (Steps 7 and 8b).
#
# The question this answers: which units in a section carry a glioma programme,
# independently of their TdTomato reporter rate? Reporter rate is carried
# through so the two can be compared, but it never enters the score.
#
# WHY THIS IS NOT A ONE-LINER — four failure modes that silently corrupt naive
# panel scoring, all observed in this dataset:
#
#   1. MARKER-LIST LENGTH. A unit with 400 significant genes hits panel genes by
#      chance far more often than one with 13. Raw rank-weighted scores in
#      7000_Scan2 correlated with list length at rho = 0.57 — i.e. a third of the
#      apparent signal was list size. Fixed by z-scoring each observed score
#      against random gene sets of the SAME LENGTH drawn from that section's
#      tested universe (--null_draws).
#
#   2. PANEL ALIASING. Some panels cannot distinguish tumour from a normal cell
#      type because they share genes:
#        - AC-like vs normal astrocyte (Gfap/Aqp4/Clu/Slc1a3/Fabp7 ...)
#        - MES-like vs myeloid, entirely through Lgals3
#      In 7000_Scan2 the second-highest scoring sub-cluster in the whole section
#      (7.8, z = 7.4) was a macrophage whose only MES gene was Lgals3. AC-like is
#      therefore excluded from the tumour score outright, and a MES call carried
#      by <= MES_ALIAS_MAX_GENES panel genes alongside >= MYELOID_MIN_GENES
#      myeloid genes is reclassified as myeloid.
#
#   3. INTERFERON RESPONSE. Shared by tumour and infiltrate, so it is evidence of
#      neither. Units where interferon outscores every tumour panel are tiered
#      ambiguous rather than called.
#
#   4. sub_id_full COLLIDES AS A NUMBER. Older runs of the R script wrote ids as
#      "<parent>.<sub>", where "7.1" and "7.10" are the same float — reading with
#      default dtypes silently merges them and drops one sub-cluster. The R script
#      now emits "<parent>_<sub>"; this reader accepts BOTH, always reads ids as
#      strings, and reports how many would have collided.
#
# Panel genes are MOUSE symbols (Neftel-style states mapped to mouse). A panel
# gene absent from the tested universe contributes nothing and cannot be
# distinguished from a gene that was tested and not enriched — so coverage is
# audited and written out. That audit is how a corrupted gene name was found in
# an earlier version of these panels (a Cyrillic character in "Tuba1b"), which
# had silently scored zero.
#
# INPUT: one output directory per section, as written by the R script. Files are
# discovered by suffix; a missing optional file degrades one feature, not the run.
#   <case>_subcluster_markers.csv         REQUIRED for sub-cluster scoring (Step 8b)
#   <case>_cluster_markers.csv            REQUIRED for cluster scoring (Step 7)
#   <case>_subcluster_tdtom.csv           optional: cap detection, missing sub-clusters
#   <case>_cluster_tdtom_enrichment.csv   optional: parent sizes, enrichment calls
#
# USAGE
#   python 2b4_tdtom_tumourlike_score.py --in OUT_7000 OUT_7003 --out scored
#   python 2b4_tdtom_tumourlike_score.py --in OUT_7000 --out scored --top_n 150
#
# OUTPUT (in --out)
#   tumourlike_units.csv            every cluster and sub-cluster, all sections:
#                                   panel z-scores, tier, driver, alias diagnostics
#   tumourlike_reproducibility.csv  written when >= 2 sections: per-parent tier
#                                   composition, and whether the call reproduces
#   tumourlike_panel_audit.csv      panel gene coverage per section — READ THIS
#                                   before trusting any panel that scores low
#   tumourlike_scores.png           figure: score ranking, reporter independence,
#                                   within-parent spread
# =============================================================================

from __future__ import annotations
import argparse, glob, os, re, sys
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- marker panels
# Positive evidence: glioma cell states + stemness/cycling.
TUMOUR_PANELS = {
    "NPC-like":      ["Dll3","Ascl1","Dcx","Sox11","Sox4","Cd24a","Stmn2","Tubb3",
                      "Nnat","Igfbpl1","Rnd3"],
    "OPC-like":      ["Pdgfra","Cspg4","Olig1","Olig2","Sox10","Bcan","Ptprz1",
                      "Lhfpl3","Cnp","Apod"],
    "AC-like":       ["Gfap","Aqp4","Apoe","Clu","S100b","Slc1a3","Mlc1","Atp1a2",
                      "Fabp7","Sparcl1"],
    "MES-like":      ["Chi3l1","Cd44","Anxa1","Anxa2","Vim","Tnc","Lgals1","Lgals3",
                      "Timp1","Adm","Vegfa"],
    "stemness":      ["Sox2","Nes","Fabp7","Ptprz1","Hes5","Id4","Mdk","Ptn"],
    "proliferation": ["Mki67","Top2a","Ccnb1","Birc5","Ube2c","Cdk1","Tuba1b",
                      "Hmgb2","Cenpf","Rrm2"],
}
# Negative evidence: normal identities that a unit could be instead.
NORMAL_PANELS = {
    "mature oligo":      ["Plp1","Mbp","Mog","Mal","Cldn11","Mobp","Trf","Cnp"],
    "microglia":         ["Cx3cr1","P2ry12","Tmem119","Csf1r","Hexb","C1qa","C1qb",
                          "Ctss","Itgam"],
    "neuron":            ["Snap25","Syt1","Rbfox3","Meg3","Nefl","Nrgn","Camk2a",
                          "Slc17a7","Gad1","Gad2"],
    "endothelial/mural": ["Cldn5","Pecam1","Flt1","Rgs5","Des","Acta2","Vtn"],
    "lymphocyte":        ["Cd3e","Cd8a","Nkg7","Ms4a1","Cd79a"],
}
# Neither: shared by tumour and infiltrate, tracked to disqualify not to score.
CONTEXT_PANELS = {
    "interferon": ["Ifit1","Ifit3","Isg15","Irf7","Stat1","B2m","Bst2","Usp18"],
}
PANELS = {**TUMOUR_PANELS, **NORMAL_PANELS, **CONTEXT_PANELS}

# AC-like is excluded from the tumour score: it cannot separate AC-like tumour
# from a normal astrocyte, and in this dataset the normal reading was correct
# every time it was checked against the unit's own top markers.
TUMOUR_SCORED = ["NPC-like","OPC-like","MES-like","stemness","proliferation"]
NORMAL_SCORED = list(NORMAL_PANELS)

# Myeloid identity beyond the microglia panel — MHC-II, complement, Fc receptors.
# Used only to catch the Lgals3/MES alias, not as a scored panel.
MYELOID_EXTRA = ["Cd74","H2-Aa","H2-Ab1","H2-Eb1","H2-DMa","Lyz2","Ptprc","Irf8",
                 "Mpeg1","Fcer1g","Cd52","Tgfbi","Cfb","Ifi30","Cd68","Aif1"]

# THIRD ALIAS: the stemness panel also reads on normal astrocytes. These of its
# genes are expressed by normal astrocytes / adult NSCs, so a stemness call built
# only from them is not evidence of a tumour stem programme. Observed in
# 7000_Scan2 sub 4.0 (Daam2/Aldh1l1/Mlc1/Vcam1 astrocyte, 8 AC-like markers, all
# 4 stemness hits from this list) which the AC-like rule alone tiered B.
ASTRO_COMPATIBLE_STEM = ["Fabp7","Hes5","Id4","Ptn"]
STEM_ALIAS_MIN_AC_GENES = 5     # AC-like genes present alongside an aliased stem call

# ------------------------------------------------------------ tier thresholds
Z_CONFIDENT        = 6.0   # tumour z at/above this, with margin, = tier A
Z_PROBABLE         = 2.5   # tumour z at/above this, with margin, = tier B
MARGIN_CONFIDENT   = 3.0   # tumour z minus best normal z
MARGIN_PROBABLE    = 2.0
AC_ASTROCYTE_Z     = 8.0   # AC-like z above this with no other programme = astrocyte
MES_ALIAS_MAX_GENES= 2     # MES call resting on <= this many MES genes ...
MYELOID_MIN_GENES  = 3     # ... alongside >= this many myeloid genes = myeloid
MICROGLIA_Z        = 2.0   # normal microglia z above this, tumour z below Z_PROBABLE
N_PROG_HITS_Z      = 2.0   # panel z above this counts as a programme "hit"


# =============================================================================
# scoring
# =============================================================================
def rank_weighted(genes_ranked: list[str], panel: list[str]) -> float:
    """Rank-weighted panel score. Weight 1/sqrt(rank) so a panel gene at rank 3
    counts for more than one at rank 140; divided by sqrt(panel size) so panels
    of different length are on a comparable raw scale (the null model below
    handles the rest)."""
    w = {g: 1.0 / np.sqrt(i + 1) for i, g in enumerate(genes_ranked)}
    return sum(w.get(g, 0.0) for g in panel) / np.sqrt(len(panel))


def null_z(genes_ranked: list[str], universe: np.ndarray, panel_names: list[str],
           indicator: np.ndarray, rng: np.random.Generator, n_draws: int) -> dict:
    """Observed rank-weighted score per panel, z-scored against `n_draws` random
    gene sets of the same length drawn from `universe`.

    This is the length correction. Without it, score tracks how many significant
    genes the unit has rather than what those genes are.

    indicator: (len(universe), n_panels) float array, 1 where universe gene is in panel.
    """
    m = min(len(genes_ranked), len(universe))
    obs = np.array([rank_weighted(genes_ranked, PANELS[p]) for p in panel_names])
    if m < 3:                      # too few genes to say anything
        return dict(zip(panel_names, np.zeros(len(panel_names))))
    w = 1.0 / np.sqrt(np.arange(1, m + 1))                  # (m,)
    idx = np.array([rng.choice(len(universe), size=m, replace=False)
                    for _ in range(n_draws)])               # (B, m)
    # (B, m, P) -> weighted sum over m -> (B, P); then normalise per panel size
    draws = np.einsum("m,bmp->bp", w, indicator[idx])
    draws = draws / np.sqrt(np.array([len(PANELS[p]) for p in panel_names]))
    mu, sd = draws.mean(0), draws.std(0)
    return dict(zip(panel_names, (obs - mu) / np.where(sd > 1e-12, sd, np.inf)))


def build_indicator(universe: list[str], panel_names: list[str]) -> np.ndarray:
    pos = {g: i for i, g in enumerate(universe)}
    M = np.zeros((len(universe), len(panel_names)))
    for j, p in enumerate(panel_names):
        for g in PANELS[p]:
            if g in pos:
                M[pos[g], j] = 1.0
    return M


def audit_panels(universe: list[str], section: str) -> pd.DataFrame:
    """Which panel genes were testable at all. A gene absent from the universe
    scores zero and is indistinguishable from a tested-but-not-enriched gene, so
    a panel with low coverage cannot be interpreted as absence of that programme."""
    uset = set(universe)
    rows = []
    for p, genes in PANELS.items():
        present = [g for g in genes if g in uset]
        rows.append(dict(section=section, panel=p, n_genes=len(genes),
                         n_in_universe=len(present),
                         coverage=len(present) / len(genes),
                         missing=";".join(g for g in genes if g not in uset)))
    return pd.DataFrame(rows)


# =============================================================================
# tier assignment
# =============================================================================
def assign_tier(r: pd.Series) -> tuple[str, str]:
    """Tier + human-readable evidence string. Order matters: exclusions are
    tested before positive calls, because an aliased panel can outscore a real
    programme (observed: a macrophage sub-cluster at z = 7.4 on Lgals3 alone)."""
    if r["mes_is_myeloid"]:
        return ("E_non_tumour",
                f"myeloid: {r['n_myeloid_genes']} myeloid genes, MES call rests on "
                f"{r['n_mes_genes']} gene(s) (Lgals3 alias)")
    if r["normal_id"] == "microglia" and r["normal_z"] > MICROGLIA_Z \
            and r["tumour_z"] < Z_PROBABLE:
        return ("E_non_tumour", f"myeloid (microglia z = {r['normal_z']:.1f})")
    if r["AC-like"] > AC_ASTROCYTE_Z and r["n_prog_hits"] < 2:
        return ("E_non_tumour",
                f"astrocyte: AC-like z = {r['AC-like']:.1f} but no non-aliased "
                f"programme (AC-like excluded from score)")
    if r["stem_is_astro"]:
        return ("E_non_tumour",
                f"astrocyte: stemness call rests on {r['n_stem_genes']} "
                f"astrocyte-compatible gene(s) alongside {r['n_ac_genes']} AC-like "
                f"markers")
    if r["normal_z"] > r["tumour_z"] and r["normal_z"] > MICROGLIA_Z:
        return ("E_non_tumour", f"{r['normal_id']} (z = {r['normal_z']:.1f}) "
                                f"outscores tumour (z = {r['tumour_z']:.1f})")
    if r["ifn_dominant"] and r["tumour_z"] < Z_CONFIDENT:
        return ("C_ambiguous",
                f"interferon-response dominant (z = {r['interferon']:.1f}); shared "
                f"by tumour and infiltrate")
    if r["tumour_z"] >= Z_CONFIDENT and r["margin"] >= MARGIN_CONFIDENT:
        return ("A_confident_tumour",
                f"{r['driver']} z = {r['tumour_z']:.1f}, no competing normal identity")
    if r["tumour_z"] >= Z_PROBABLE and r["margin"] >= MARGIN_PROBABLE:
        return ("B_probable_tumour", f"{r['driver']} z = {r['tumour_z']:.1f}")
    return ("D_weak", f"no programme above noise (max z = {r['tumour_z']:.1f})")


# =============================================================================
# input discovery and loading
# =============================================================================
def find_one(d: str, suffix: str) -> str | None:
    hits = sorted(glob.glob(os.path.join(d, f"*{suffix}")))
    return hits[0] if hits else None


def case_of(path: str, suffix: str) -> str:
    return os.path.basename(path)[: -len(suffix)]


def load_section(d: str) -> dict:
    """Load one section's tables. sub_id / sub_id_full are read as STRINGS —
    see failure mode 4 in the header."""
    sm_p = find_one(d, "_subcluster_markers.csv")
    cm_p = find_one(d, "_cluster_markers.csv")
    if sm_p is None and cm_p is None:
        raise SystemExit(f"ERROR: no _subcluster_markers.csv or _cluster_markers.csv "
                         f"in {d}\n  (did Step 7 / Step 8b run? check "
                         f"*_run_parameters.csv for subcluster_de_run)")
    case = case_of(cm_p or sm_p, "_cluster_markers.csv"
                   if cm_p else "_subcluster_markers.csv")
    out = dict(case=case, dir=d, sm=None, cm=None, st=None, enr=None)

    if sm_p:
        sm = pd.read_csv(sm_p, dtype={"sub_id": str, "sub_id_full": str,
                                      "parent_cluster": str})
        # collision guard: if these ids were read as floats, how many would merge?
        as_float = pd.to_numeric(sm["sub_id_full"], errors="coerce")
        n_str, n_flt = sm["sub_id_full"].nunique(), as_float.dropna().nunique() + \
            sm.loc[as_float.isna(), "sub_id_full"].nunique()
        if n_flt < n_str:
            print(f"  NOTE {case}: {n_str - n_flt} sub-cluster id(s) would collide if "
                  f"this column were type-inferred (7.1 vs 7.10) — read as strings "
                  f"here; re-run the R script to get '_' separators")
        # parent_downsampled is authoritative when present (newer R runs); when the
        # column is absent, check_coverage() infers truncation from parent sizes.
        if "parent_downsampled" in sm.columns:
            capped = sm.loc[sm["parent_downsampled"].astype(str).str.upper()
                            .isin(["TRUE", "1"]), "parent_cluster"].unique()
            for p in sorted(capped):
                n_true = sm.loc[sm["parent_cluster"] == p, "parent_n_labelled"].iloc[0]
                print(f"  WARNING {case} cl{p}: parent_downsampled = TRUE "
                      f"(true size {n_true}) — sub-cluster n_cells are a SAMPLE; "
                      f"re-run with --subcluster_max_cells 0 for absolute counts")
        out["sm"] = sm

    if cm_p:
        out["cm"] = pd.read_csv(cm_p, dtype={"cluster": str})

    st_p = find_one(d, "_subcluster_tdtom.csv")
    if st_p:
        out["st"] = pd.read_csv(st_p, dtype={"sub_id": str, "parent_cluster": str})
    enr_p = find_one(d, "_cluster_tdtom_enrichment.csv")
    if enr_p:
        out["enr"] = pd.read_csv(enr_p, dtype={"cluster": str})
    return out


def check_coverage(sec: dict) -> pd.DataFrame:
    """Sub-cluster cells vs parent size. Step 8's --subcluster_max_cells caps the
    cells sent to re-clustering, so a capped parent's sub-clusters describe a
    SAMPLE of it. Also reports sub-clusters present in the rate table but absent
    from the marker table (skipped by --subcluster_de_min_cells)."""
    if sec["sm"] is None or sec["st"] is None:
        return pd.DataFrame()
    st, sm = sec["st"], sec["sm"]
    rows = []
    for parent, g in st.groupby("parent_cluster"):
        n_sub_cells = g["n_cells"].sum()
        parent_n = np.nan
        if sec["enr"] is not None:
            e = sec["enr"][sec["enr"]["cluster"] == parent]
            if len(e):
                parent_n = float(e["n_cells"].iloc[0])
        # match on the bare sub_id, so this works for '_' and legacy '.' ids alike
        ids_rate = set(g["sub_id"].astype(str))
        ids_mk = {i[len(parent) + 1:] for i in
                  sm.loc[sm["parent_cluster"] == parent, "sub_id_full"]}
        rows.append(dict(section=sec["case"], parent=parent,
                         n_subclusters_rate=len(g), n_subclusters_scored=len(ids_mk),
                         cells_in_subclusters=int(n_sub_cells), parent_n_cells=parent_n,
                         coverage=(n_sub_cells / parent_n) if parent_n == parent_n else np.nan,
                         not_scored=";".join(sorted(ids_rate - ids_mk))))
    return pd.DataFrame(rows)


# =============================================================================
# per-section scoring
# =============================================================================
def score_section(sec: dict, top_n: int, n_draws: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    panel_names = list(PANELS)
    units = []

    # --- sub-clusters: SECTION comparator only.
    # The 'parent' comparator asks what separates siblings; genes shared with the
    # parent cancel out, so a tumour sub-cluster inside a tumour cluster shows no
    # tumour markers there. It is the wrong table for identity scoring.
    if sec["sm"] is not None:
        sm = sec["sm"]
        if "comparator" in sm.columns:
            sm = sm[sm["comparator"] == "section"]
        sm = sm[sm["marker_rank"].notna()]
        universe = sorted(sec["sm"].loc[sec["sm"].get("comparator", "section") ==
                                        "section", "gene"].unique()) \
            if "comparator" in sec["sm"].columns else sorted(sec["sm"]["gene"].unique())
        ind = build_indicator(universe, panel_names)
        uarr = np.asarray(universe)
        for sid, g in sm.groupby("sub_id_full"):
            g = g.sort_values("marker_rank")
            genes = list(g["gene"].head(top_n))
            units.append(dict(
                section=sec["case"], unit=f"sub {sid}", kind="subcluster",
                parent=str(g["parent_cluster"].iloc[0]),
                n_cells=float(g["sub_n_cells"].iloc[0]),
                pct_tdtom_pos=float(g["sub_pct_tdtom_pos"].iloc[0]),
                parent_pct_tdtom_pos=float(g["parent_pct_tdtom_pos"].iloc[0]),
                n_markers=len(g), top5="; ".join(genes[:5]),
                **null_z(genes, uarr, panel_names, ind, rng, n_draws),
                **_alias_counts(genes)))

    # --- clusters, same method and same section, so the two are comparable
    if sec["cm"] is not None:
        cm = sec["cm"][sec["cm"]["marker_rank"].notna()]
        universe = sorted(sec["cm"]["gene"].unique())
        ind = build_indicator(universe, panel_names)
        uarr = np.asarray(universe)
        enr = sec["enr"].set_index("cluster") if sec["enr"] is not None else None
        for cl, g in cm.groupby("cluster"):
            g = g.sort_values("marker_rank")
            genes = list(g["gene"].head(top_n))
            n_cells = float(enr.loc[cl, "n_cells"]) if enr is not None and cl in enr.index else np.nan
            pct = float(g["pct_tdtom_pos"].iloc[0]) if "pct_tdtom_pos" in g.columns else np.nan
            units.append(dict(
                section=sec["case"], unit=f"cl {cl}", kind="cluster", parent=str(cl),
                n_cells=n_cells, pct_tdtom_pos=pct, parent_pct_tdtom_pos=np.nan,
                n_markers=len(g), top5="; ".join(genes[:5]),
                **null_z(genes, uarr, panel_names, ind, rng, n_draws),
                **_alias_counts(genes)))

    U = pd.DataFrame(units)
    if U.empty:
        return U

    # A panel with no genes in the tested universe scores exactly 0 for every unit.
    # Left in, idxmax picks it as the "best" normal identity whenever every real
    # normal panel is negative — which is how 'lymphocyte' became the reported
    # identity of 20 units in 7000_Scan2 at z = 0.00. Drop uncovered panels from
    # the argmax entirely; they carry no information either way.
    universe_all = set()
    for key in ("sm", "cm"):
        if sec[key] is not None:
            universe_all |= set(sec[key]["gene"])
    uncovered = [p for p in panel_names
                 if not any(g in universe_all for g in PANELS[p])]
    if uncovered:
        print(f"  NOTE excluded from best-identity calls (0 genes testable): "
              f"{', '.join(uncovered)}")
    t_use = [p for p in TUMOUR_SCORED if p not in uncovered]
    n_use = [p for p in NORMAL_SCORED if p not in uncovered]
    if not t_use:
        # Not fatal to the run — one section can have a gene namespace the panels
        # don't cover (wrong species symbols, or an object whose features are not
        # gene names) while its siblings are fine. Skip it and say so.
        print("  SKIPPING SECTION: no tumour panel has any testable gene here — "
              "check that feature names are mouse gene symbols")
        return pd.DataFrame()

    U["tumour_z"]     = U[t_use].max(axis=1)
    U["driver"]       = U[TUMOUR_SCORED].idxmax(axis=1)
    U["n_prog_hits"]  = (U[TUMOUR_SCORED] > N_PROG_HITS_Z).sum(axis=1)
    U["normal_z"]     = U[NORMAL_SCORED].max(axis=1)
    U["normal_id"]    = U[NORMAL_SCORED].idxmax(axis=1)
    U["margin"]       = U["tumour_z"] - U["normal_z"]
    U["ifn_dominant"] = U["interferon"] > U[TUMOUR_SCORED].max(axis=1)
    U["mes_is_myeloid"] = ((U["driver"] == "MES-like")
                           & (U["n_mes_genes"] <= MES_ALIAS_MAX_GENES)
                           & (U["n_myeloid_genes"] >= MYELOID_MIN_GENES))
    # a stemness call whose every hit is a gene normal astrocytes also express,
    # in a unit carrying many AC-like markers, is an astrocyte not a stem programme
    U["stem_is_astro"] = ((U["driver"] == "stemness") & U["stem_all_astro"]
                          & (U["n_ac_genes"] >= STEM_ALIAS_MIN_AC_GENES))
    U[["tier", "evidence"]] = U.apply(lambda r: pd.Series(assign_tier(r)), axis=1)
    U["reporter_delta"] = U["pct_tdtom_pos"] - U["parent_pct_tdtom_pos"]
    return U


def _alias_counts(genes: list[str]) -> dict:
    mes = set(TUMOUR_PANELS["MES-like"])
    mye = set(NORMAL_PANELS["microglia"]) | set(MYELOID_EXTRA)
    ifn = set(CONTEXT_PANELS["interferon"])
    ac = set(TUMOUR_PANELS["AC-like"])
    stem = set(TUMOUR_PANELS["stemness"])
    astro_ok = set(ASTRO_COMPATIBLE_STEM)
    stem_hits = [g for g in genes if g in stem]
    return dict(n_mes_genes=sum(g in mes for g in genes),
                n_myeloid_genes=sum(g in mye for g in genes),
                n_ifn_genes=sum(g in ifn for g in genes),
                n_ac_genes=sum(g in ac for g in genes),
                n_stem_genes=len(stem_hits),
                # every stemness gene found is one a normal astrocyte also expresses
                stem_all_astro=bool(stem_hits) and all(g in astro_ok for g in stem_hits))


# =============================================================================
# cross-section reproducibility
# =============================================================================
def reproducibility(U: pd.DataFrame) -> pd.DataFrame:
    """Per-parent-cluster tier composition in each section, and whether the call
    agrees. Cluster numbering is NOT stable across independent clusterings, so
    this compares by cluster id only as a starting point — the driver programme
    is the transferable quantity, and is reported alongside so a mismatch in id
    with a match in programme is visible."""
    secs = sorted(U["section"].unique())
    if len(secs) < 2:
        return pd.DataFrame()
    rows = []
    for parent in sorted(U["parent"].unique(), key=lambda x: (len(x), x)):
        rec = dict(parent=parent)
        calls = []
        for s in secs:
            g = U[(U["section"] == s) & (U["parent"] == parent)]
            if not len(g):
                rec[f"{s}_tier"] = "absent"; rec[f"{s}_driver"] = ""
                rec[f"{s}_tumour_z"] = np.nan; continue
            cl = g[g["kind"] == "cluster"]
            best = (cl if len(cl) else g).sort_values("tumour_z").iloc[-1]
            rec[f"{s}_tier"] = best["tier"]
            rec[f"{s}_driver"] = best["driver"]
            rec[f"{s}_tumour_z"] = round(float(best["tumour_z"]), 2)
            sub = g[g["kind"] == "subcluster"]
            rec[f"{s}_n_sub_AB"] = int(sub["tier"].str.startswith(("A", "B")).sum())
            rec[f"{s}_n_sub"] = int(len(sub))
            calls.append((best["tier"][0], best["driver"]))
        if len(calls) == len(secs):
            rec["tier_letter_agrees"] = len({c[0] for c in calls}) == 1
            rec["driver_agrees"] = len({c[1] for c in calls}) == 1
        rows.append(rec)
    return pd.DataFrame(rows)


# =============================================================================
# figure
# =============================================================================
TIER_COLOUR = {"A_confident_tumour": "#08519c", "B_probable_tumour": "#6baed6",
               "C_ambiguous": "#fdd0a2", "D_weak": "#bdbdbd",
               "E_non_tumour": "#a50f15"}
TIER_LABEL = {"A_confident_tumour": "A confident tumour",
              "B_probable_tumour": "B probable tumour",
              "C_ambiguous": "C ambiguous (interferon)",
              "D_weak": "D no programme", "E_non_tumour": "E not tumour"}


def make_figure(U: pd.DataFrame, path: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:                                    # house style if available
        from figure_style import apply_figure_style        # noqa
        apply_figure_style(sizes=(9, 8, 7))
    except Exception:
        plt.rcParams.update({"font.size": 8, "axes.spines.top": False,
                             "axes.spines.right": False, "figure.dpi": 110})

    sub = U[U["kind"] == "subcluster"]
    if sub.empty:
        sub = U
    secs = sorted(U["section"].unique())
    # Panel a needs a fixed vertical allowance PER ROW or the tick labels overlap
    # once several sections are pooled (62 rows collided at the old fixed height).
    # Give it 0.155 in per row and let the lower panels keep a constant height, so
    # the height ratio follows the row count instead of being hard-coded.
    # Margins are computed in INCHES then converted, not given as fractions: a
    # fractional bottom margin shrinks in absolute terms as the figure grows with
    # row count, which is how the shared legend ended up on top of the panel b/c
    # axis labels on a short figure.
    n_rows = len(sub)
    h_a, h_low = max(3.0, 0.155 * n_rows), 3.4
    # PAD_MID has to clear panel a's x-label AND the two-line titles of panels b/c
    PAD_TOP, PAD_BOT, PAD_MID = 0.75, 1.25, 1.30   # inches
    H = h_a + h_low + PAD_TOP + PAD_BOT + PAD_MID
    fig = plt.figure(figsize=(11.4, H))
    gs = fig.add_gridspec(2, 2, height_ratios=[h_a, h_low], wspace=.28,
                          hspace=PAD_MID / ((h_a + h_low) / 2),
                          left=.15, right=.97,
                          top=1 - PAD_TOP / H, bottom=PAD_BOT / H)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])

    # a  ranked tumour score, colour = tier (i.e. after alias correction)
    o = sub.sort_values("tumour_z")
    y = np.arange(len(o))
    axA.barh(y, o["tumour_z"], color=[TIER_COLOUR[t] for t in o["tier"]],
             height=.72, edgecolor="none")
    axA.set_yticks(y)
    axA.set_yticklabels([f"{r.section.split('_')[0]} {r.unit}" if len(secs) > 1
                         else r.unit for r in o.itertuples()])
    axA.axvline(0, color="0.3", lw=.7)
    axA.axvline(Z_PROBABLE, color="0.5", lw=.7, ls=":")
    axA.axvline(Z_CONFIDENT, color="0.5", lw=.7, ls="--")
    axA.set_xlabel("tumour-programme enrichment (z vs length-matched random gene sets)")
    axA.set_title("a  Sub-cluster tumour score — bar = raw score, colour = verdict after\n"
                  "    removing AC-like (astrocyte) and Lgals3-only (myeloid) aliases",
                  loc="left", pad=8)
    for xx, lab in [(Z_PROBABLE, "probable"), (Z_CONFIDENT, "confident")]:
        axA.text(xx, len(o) - .2, lab, fontsize=6.5, color="0.35", ha="center", va="bottom")

    # b  the reporter-independence result
    ok = sub[sub["pct_tdtom_pos"].notna()]
    for t, g in ok.groupby("tier"):
        axB.scatter(g["pct_tdtom_pos"], g["tumour_z"],
                    s=18 + g["n_cells"].fillna(200) / 90, c=TIER_COLOUR[t],
                    edgecolor="0.25", lw=.4, label=TIER_LABEL[t], zorder=3)
    rho = ok[["pct_tdtom_pos", "tumour_z"]].corr(method="spearman").iloc[0, 1] \
        if len(ok) > 2 else np.nan
    axB.axhline(Z_PROBABLE, color="0.5", lw=.7, ls=":")
    axB.set_xlabel("TdTomato+ rate of sub-cluster (%)")
    axB.set_ylabel("tumour-programme z")
    axB.set_title(f"b  Reporter rate vs tumour score\n"
                  f"    (Spearman $\\rho$ = {rho:+.2f}, n = {len(ok)} sub-clusters)",
                  loc="left", pad=6)

    # c  within-parent spread: what sub-clustering buys
    rows = []
    for (s, p), g in sub.groupby(["section", "parent"]):
        rows.append(dict(lab=(f"{s.split('_')[0]} cl{p}" if len(secs) > 1 else f"cl{p}"),
                         zmin=g["tumour_z"].min(), zmax=g["tumour_z"].max(),
                         n_sub=len(g),
                         n_ab=int(g["tier"].str.startswith(("A", "B")).sum()),
                         n_e=int((g["tier"] == "E_non_tumour").sum())))
    P = pd.DataFrame(rows).sort_values("zmax")
    yy = np.arange(len(P))
    axC.hlines(yy, P["zmin"], P["zmax"], color="0.55", lw=2.4, zorder=2)
    axC.scatter(P["zmax"], yy, s=34, c="#08519c", zorder=3, label="most tumour-like sub-cluster")
    axC.scatter(P["zmin"], yy, s=34, c="#a50f15", zorder=3, label="least")
    axC.set_yticks(yy)
    axC.set_yticklabels([f"{r.lab}  ({r.n_sub} sub)" for r in P.itertuples()])
    axC.axvline(Z_PROBABLE, color="0.5", lw=.7, ls=":")
    axC.set_xlabel("tumour-programme z")
    axC.set_title("c  Sub-clustering splits mixed parents\n    into tumour and non-tumour units",
                  loc="left", pad=6)
    for i, r in enumerate(P.itertuples()):
        axC.text(r.zmax + .4, i, f"{r.n_ab}A/B, {r.n_e}E", fontsize=6.5,
                 va="center", color="0.3")
    axC.legend(loc="lower right", frameon=False, fontsize=6.5)

    h, l = axB.get_legend_handles_labels()
    if h:
        fig.legend(h, l, loc="lower center", ncol=len(h), frameon=False,
                   bbox_to_anchor=(.55, 0.22 / H), fontsize=7,
                   handletextpad=.4, columnspacing=1.1)
    # Self-check: y tick labels in panel a must not overlap each other. This is the
    # defect that appeared only once two sections were pooled, so it is asserted
    # rather than eyeballed.
    fig.canvas.draw()
    bb = [t.get_window_extent() for t in axA.get_yticklabels() if t.get_text().strip()]
    clash = sum(1 for i in range(len(bb) - 1) if bb[i].overlaps(bb[i + 1]))
    if clash:
        print(f"  WARNING figure: {clash} overlapping row label(s) in panel a — "
              f"raise the per-row allowance in make_figure()")
    # ... panel a's x-label against the lower panels' titles ...
    xa = axA.xaxis.label.get_window_extent()
    if [a for a in (axB, axC) if xa.overlaps(a.title.get_window_extent())]:
        print("  WARNING figure: panel a x-label overlaps a lower panel title — "
              "raise PAD_MID in make_figure()")
    # ... and the shared legend against the lower panels' axis labels
    legs = [c for c in fig.legends]
    if legs:
        lb = legs[0].get_window_extent()
        hit = [a for a, nm in ((axB, "b"), (axC, "c"))
               if lb.overlaps(a.xaxis.label.get_window_extent())]
        if hit:
            print(f"  WARNING figure: legend overlaps {len(hit)} axis label(s) — "
                  f"raise PAD_BOT in make_figure()")
    fig.savefig(path, dpi=200)
    plt.close(fig)


# =============================================================================
def main() -> None:
    ap = argparse.ArgumentParser(
        description="Score clusters and sub-clusters on tumour-likeness from "
                    "2b1_tdtom_cluster_enrichment_v4.R marker tables.")
    ap.add_argument("--in", dest="indirs", nargs="+", required=True,
                    help="one output directory per section")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--top_n", type=int, default=150,
                    help="markers per unit to score, by rank [150]")
    ap.add_argument("--null_draws", type=int, default=2000,
                    help="random gene sets per unit for the length correction [2000]")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    all_units, all_audit, all_cov = [], [], []
    for d in args.indirs:
        print(f"=== {d}")
        sec = load_section(d)
        print(f"  case: {sec['case']}  "
              f"[{'sub' if sec['sm'] is not None else 'no sub'}, "
              f"{'cluster' if sec['cm'] is not None else 'no cluster'}]")
        univ = sorted(set(
            (sec["sm"]["gene"].tolist() if sec["sm"] is not None else []) +
            (sec["cm"]["gene"].tolist() if sec["cm"] is not None else [])))
        aud = audit_panels(univ, sec["case"])
        all_audit.append(aud)
        thin = aud[aud["coverage"] < 0.5]
        for r in thin.itertuples():
            print(f"  WARNING panel '{r.panel}' only {r.n_in_universe}/{r.n_genes} genes "
                  f"testable — a low score for it is NOT evidence of absence")
        cov = check_coverage(sec)
        if len(cov):
            all_cov.append(cov)
            for r in cov.itertuples():
                if r.coverage == r.coverage and r.coverage < 0.98:
                    print(f"  WARNING cl{r.parent}: sub-clusters cover "
                          f"{r.coverage:.0%} of {int(r.parent_n_cells)} cells — capped by "
                          f"--subcluster_max_cells; sub-cluster sizes are a SAMPLE")
                if r.not_scored:
                    print(f"  NOTE cl{r.parent}: not scored (below "
                          f"--subcluster_de_min_cells): {r.not_scored}")
        U = score_section(sec, args.top_n, args.null_draws, args.seed)
        if U.empty:
            print("  no scorable units"); continue
        s = U[U["kind"] == "subcluster"]
        if len(s) > 2:
            print(f"  length confound after correction: "
                  f"spearman(tumour_z, n_markers) = "
                  f"{s[['tumour_z','n_markers']].corr(method='spearman').iloc[0,1]:+.2f}")
            print(f"  reporter vs tumour score: spearman = "
                  f"{s[['pct_tdtom_pos','tumour_z']].corr(method='spearman').iloc[0,1]:+.2f}")
        print("  tiers: " + ", ".join(f"{k}={v}" for k, v in
                                      sorted(U["tier"].value_counts().items())))
        all_units.append(U)

    if not all_units:
        raise SystemExit("ERROR: nothing scored")
    U = pd.concat(all_units, ignore_index=True)

    lead = ["section", "unit", "kind", "parent", "n_cells", "pct_tdtom_pos",
            "parent_pct_tdtom_pos", "reporter_delta", "n_markers", "tier", "evidence",
            "tumour_z", "driver", "n_prog_hits", "normal_z", "normal_id", "margin",
            "interferon", "n_mes_genes", "n_myeloid_genes", "n_ifn_genes",
            "n_ac_genes", "n_stem_genes", "mes_is_myeloid", "stem_is_astro", "top5"]
    U = U[lead + [c for c in U.columns if c not in lead]]
    U = U.sort_values(["section", "tier", "tumour_z"], ascending=[True, True, False])
    p_units = os.path.join(args.out, "tumourlike_units.csv")
    U.to_csv(p_units, index=False)

    p_aud = os.path.join(args.out, "tumourlike_panel_audit.csv")
    pd.concat(all_audit, ignore_index=True).to_csv(p_aud, index=False)
    if all_cov:
        pd.concat(all_cov, ignore_index=True).to_csv(
            os.path.join(args.out, "tumourlike_coverage.csv"), index=False)
    R = reproducibility(U)
    if len(R):
        R.to_csv(os.path.join(args.out, "tumourlike_reproducibility.csv"), index=False)
        agree = R["driver_agrees"].sum() if "driver_agrees" in R else 0
        print(f"\n=== cross-section: driver programme agrees for {agree}/{len(R)} parents")
    p_fig = os.path.join(args.out, "tumourlike_scores.png")
    make_figure(U, p_fig)

    print(f"\n=== wrote {args.out}/")
    for f in sorted(os.listdir(args.out)):
        print(f"  {f}")
    pool = U[(U["kind"] == "subcluster") & U["tier"].str.startswith(("A", "B"))]
    if len(pool):
        print("\n=== DE-usable sub-clusters (tier A/B) ===")
        print(pool[["section", "unit", "n_cells", "pct_tdtom_pos", "tumour_z", "driver"]]
              .to_string(index=False, float_format=lambda x: f"{x:.1f}"))
        print(f"total cells {int(pool['n_cells'].sum())}; weighted reporter rate "
              f"{(pool['pct_tdtom_pos'] * pool['n_cells']).sum() / pool['n_cells'].sum():.1f}%")


if __name__ == "__main__":
    main()

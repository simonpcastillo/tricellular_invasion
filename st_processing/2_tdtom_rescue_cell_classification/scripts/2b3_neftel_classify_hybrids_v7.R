# =============================================================================
# Neftel GBM state + Cell type classification — expanded pipeline (v7)
# - TdTomato label from gene expression OR predicted (--use_predicted)
# - Neftel scoring on TdTom+ / cell-type scoring on TdTom-
# - Proliferation state (G1S / G2M) on ALL cells via permutation
# - EC signature scoring on Endothelial cells, unique genes, reported independently
# - MYELOID STATE sub-typing on microglia/macrophages (v3, Step 6c): two
#   published signature sets scored independently, permutation null + hybrids
# - DE analysis TdTom+ vs TdTom-
# - GeoJSON export: neftel_normal, proliferation, intersection, ec_state
# - Summary tables + null-calibration diagnostics
#   v7 (sparse-panel remedy). Five changes, all reachable from the CLI. Two are
#   correctness fixes and are ON by default; three are opt-in.
#
#   [FIX, on]  --borrow_strata TRUE. v6 line 1154 did pool_ex[empty] <- pool[empty],
#              which let the matched null RESAMPLE THE SIGNATURE GENES when a
#              stratum held nothing else. The null mean then climbs toward the
#              observed score and z shrinks, so the best-detected signatures were
#              silently penalised. v7 widens the mean-expression window until it
#              finds >= min_pool clean controls and never resamples the signature.
#              Median z of a true state in its own cells rose 0.47 -> 0.93 on a
#              simulated 300-plex panel. On whole-transcriptome input the empty-
#              stratum case it repairs is rare, so expect little or no change: the
#              run log prints how many strata actually borrowed, and if that count
#              is 0 this flag did nothing for you.
#
#   [FIX, on]  --auto_strata TRUE. Only ever REDUCES nbin/n_det_bin, and only
#              when the gene count cannot support them. On whole-transcriptome
#              input (~18000 genes) it is a NO-OP and returns 25x5 unchanged; on a
#              300-gene panel it coarsens to ~10x3. Kept on because it is safe in
#              both regimes, not because whole-transcriptome data needs it.
#
# ---------------------------------------------------------------------------
#   Measured on the two mouse sections (7000_Scan2 + 7003_Scan2, 208k cells,
#   whole transcriptome), by re-correcting the enrichment_stats.csv files.
#   Not simulated.
#
#   Your v6 gate contains TWO errors that push in OPPOSITE directions:
#
#     (i)  the denominator is too strict. BH across cells within a signature
#          makes m = your cell count (52,729), so a cell needs raw p < 9.5e-07,
#          i.e. z > 4.76, to be called.
#     (ii) the tail is too generous. With n_perm = 1000 the permutation p cannot
#          go below 1/1001 = 1.0e-03, so the gate operates ~1000x beyond what the
#          permutations resolve, and every call rests on pnorm(z) extrapolation.
#          Where both are measurable, that extrapolation is ANTI-CONSERVATIVE:
#          14x too small at p=0.004, 6x at 0.008, 3x at 0.015 — worsening outward.
#
#   Net effect: 24.4% / 17.1% of tumour cells assigned. Fixing only (i) inflates
#   false positives; fixing only (ii) collapses yield. Fix BOTH — --gate_unit
#   bh_cell --p_source emp — and assignment is essentially UNCHANGED (24.5% /
#   19.7%), because the two errors were cancelling. What changes is WHICH cells:
#
#       4,318 cells recovered that v6 left unassigned
#      23,880 v6 calls DROPPED — they were supported only by the bad tail
#
#   So the honest answer to "make it less stringent" is that v6 was not merely
#   too strict; it was mis-calibrated in both directions at once, and its
#   assignment rate was roughly right for the wrong reasons.
#
#   Composition is the part to check before anything else. Under v6 the
#   per-signature pass rate spanned 13.7-fold (MES1_like 10.45% vs NPC2_like
#   0.76%) because each signature gets its own adaptive BH threshold; under
#   bh_cell + emp that spread falls to 3.9-fold. Recovered proportions move by
#   up to 6.4 points (MES1_like 20.8% -> 14.4%), and MES1/AC lose share while
#   OPC/NPC gain in BOTH samples. If you have circulated Neftel state
#   proportions from v6, they are biased toward the signatures with the densest
#   p-distributions, and that is a systematic error, not noise.
#
#   Recommended for this dataset:
#     --gate_unit bh_cell --p_source emp        (and --n_perm 2000 if affordable)
# ---------------------------------------------------------------------------
#
#   WHICH INPUT DO YOU HAVE? The remedy differs, and the difference is large.
#
#   WHOLE TRANSCRIPTOME (~18000 genes, e.g. Visium / Stereo-seq / snPATHO-seq):
#     125 strata hold ~144 control genes each, so the matched null is well
#     estimated and BOTH default fixes above are essentially no-ops for you.
#     Your unassigned cells are NOT a strata problem. What binds instead:
#       1. --gate_unit. This is the whole remedy in your case. v6 corrects across
#          CELLS within a signature, so the BH denominator is your cell count.
#          At 50000 cells the effective per-cell cutoff is orders of magnitude
#          below 0.05 and the gate cannot fire regardless of sparsity. bh_cell
#          corrects across the K signatures within a cell — a denominator of ~6
#          instead of ~50000. Nothing about this depends on panel size.
#       2. Counts per unit. A 50-gene signature over ~2000 counts is a different
#          statistic than over ~200. If units are small bins, aggregating to
#          larger bins helps more than any gate change.
#       3. --pool_k, only if Moran's I says the biology is patchy.
#     Do NOT expect much from --auto_strata or --borrow_strata.
#
#   TARGETED PANEL (~300 genes, e.g. Xenium / MERFISH / CosMx):
#     the two fixes above matter, because most strata are then nearly empty.
#     Measured on a simulated 300-plex panel (4000 cells, 84% zeros, 4 states):
#
#         config                          assigned   acc   FP(no signal)
#         v6 as shipped                       4.2%  70.6%      2.30%
#         v7 defaults (both fixes)            1.8%  80.8%      0.32%
#         v7 defaults + --gate_unit bh_cell  13.3%  78.6%      7.60%
#
#     Note the fixes alone LOWER yield while raising accuracy — part of v6's
#     yield was miscalibration. --gate_unit buys the coverage back.
#
#   THOSE NUMBERS ARE PANEL NUMBERS. They are not a prediction for
#   whole-transcriptome input, and no equivalent table has been measured for it.
#   The --gate_unit argument above is structural (it follows from the BH
#   denominator, not from a simulation) and is the one to act on.
# ---------------------------------------------------------------------------
#
#              bh_cell caveat: on simulated no-signal panel data it assigned 7.6%
#              of cells, above the nominal 5%. Treat --alpha 0.05 under bh_cell as
#              approximate and use --alpha 0.01 if the gate must mean what it says.
#
#   [opt-in]   --gate_unit. THE STRINGENCY LEVER, and the one to think about.
#              v6 gated on padj_z = BH across CELLS within each signature, i.e.
#              every cell treated as an independent discovery. That is right if a
#              single positive cell is itself the claim; it is the wrong unit when
#              the per-cell label is an intermediate for composition, neighbourhood
#              or group comparison. bh_cell corrects across the K signatures WITHIN
#              a cell — the multiplicity of the decision actually being made.
#              Simulated: 1.8% -> 13.3% assigned at 78.6% accuracy (table above).
#              This is the single biggest lever on your unassigned count.
#
#              CONSEQUENCE WORTH READING. Because bh_global gives each signature
#              its own adaptive threshold, per-signature yield varied 8.4-fold in
#              simulation and a state truly 17% of cells came back as 80% of
#              assigned cells (total-variation distance from truth 0.63 vs 0.27
#              under bh_cell). If state PROPORTIONS have been reported from a
#              bh_global run, re-check them: accuracy among assigned looks fine
#              while composition is distorted, and nothing in the output flags it.
#
#   [opt-in]   --pool_k with --pool_coords. Average the log-normalised matrix over
#              each cell and its k nearest neighbours before scoring. On patchy
#              simulated tissue: 4.2% -> 11.6% assigned, accuracy 72 -> 87. BUT at
#              MATCHED coverage pooling was 3-6 points WORSE than no pooling — it
#              buys coverage by lowering the bar, not by ranking cells better. Use
#              real centroids, keep k at 6-10, check Moran's I first, and
#              never present a pooled label
#              as evidence that neighbouring cells share a state: the smoothing
#              assumes it. Report the headline result at k=0 as well.
#
#   [opt-in]   --skew_p. The permutation null of a signature mean is right-skewed
#              on sparse counts (median skewness 0.35 at 300 genes / 60 counts per
#              cell). pnorm() is accurate at p~0.05 (measured 0.051) but 3.5x
#              anti-conservative at p~0.001 — exactly where a global-BH gate
#              operates. Cornish-Fisher correction. Makes the gate STRICTER.
#              Measured: it pulls bh_cell's no-signal rate from 7.6% to 0.6% but
#              cuts yield 13.3% -> 2.2%, and with bh_global it drove assignment to
#              zero. Use it as an audit (does my headline result survive?), not as
#              a production setting, and never together with a low --alpha.
#
#   Also added: state_tier / tier_* metadata columns, which separate cells with no
#   programme evidence from cells whose evidence missed a multiplicity threshold.
#   v6 reported both as unassigned, which is what makes the output look empty.
#   Defaults keep cell_class identical to v6 apart from the two fixes above.
# =============================================================================
# WHAT CHANGED FROM v2 (Step 3, Step 6, new Step 6c, new plot 8g)
#
#   A. NEW Step 6c — myeloid state sub-typing. Cells whose cell_class is
#      Microglia_core / Macrophage_peripheral / Macrophage_BAM are sub-typed
#      against two published STATE signature sets (--miller_signatures,
#      --qian_signatures), each scored independently: AddModuleScore for
#      plotting, sig_enrich_perm() z + BH gate for inference, then
#      assign_states_with_hybrids() for single-state and hybrid calls. Identity
#      (Step 6) and state (Step 6c) stay separate. Full rationale, the power
#      table and the criterion-(d) evidence are in the Step 6c header.
#
#   B. Criterion (d) is DISABLED by default for myeloid hybrids
#      (--myeloid_hybrid_delta 0), a deliberate departure from the Neftel block.
#      The auto-calibrated delta is calibrated by shuffling score columns, so a
#      stronger signature widens its own null and the threshold rises with the
#      signal it should admit. Simulation with known truth: auto recovers 15% of
#      true hybrids, delta=0 recovers 43%, false positives 0 in both.
#      --myeloid_hybrid_delta auto restores v2 behaviour. The Neftel hybrid call
#      is UNCHANGED.
#
#   C. Order-free hybrid labels (canon_pair()). assign_states_with_hybrids()
#      names hybrids by score RANK, so MES_like|NPC_like and NPC_like|MES_like
#      were counted as two classes. New hybrid_class_canonical /
#      hybrid_sub_canonical columns fix this for the NEFTEL calls too; the
#      rank-ordered columns are kept because they say which state scored higher.
#
#   D. State sets have identity genes stripped by default
#      (--myeloid_strip_identity). The myeloid identity signatures are what
#      selects the compartment in Step 6, so a state gene that is also an
#      identity gene makes the state test partly re-test its own selection
#      criterion. Removed genes are logged either way, per state, to
#      *_myeloid_state_gene_audit.csv.
#
#   E. Qian MG-CT-Prolif + TR-Prolif-MG are collapsed to family MG-Prolif
#      (MYELOID_COLLAPSE). They share Jaccard 0.47 so cannot be ranked against
#      each other, and both are near-entirely cell-cycle genes, duplicating the
#      genome-wide proliferation test in Step 6 — Step 6c cross-tabulates the
#      family against prolif_state rather than reporting it as new.
#      collinearity_report equivalents run on the genes ACTUALLY tested, so the
#      collapse decision is re-checked on your object after orthology.
#
#   F. Macrophage_BAM colour was "#E41A1C" — byte-identical to MES_like, so BAM
#      cells and MES_like tumour cells were indistinguishable in every plot and
#      GeoJSON. Changed to "#22534E". See the note there about NPC2_like, which
#      is still the palette's closest pair.
#
#   G. read_myeloid_csv() strips a UTF-8 byte-order mark. qian_2026_microglia.csv
#      carries one; unstripped, its first state name silently fails every name
#      match (including the collapse spec) while still scoring under the mangled
#      name.
#
#   POWER CAVEAT, applies to every myeloid hybrid number. Hybrid recall is
#   roughly the SQUARE of single-state power, because criterion (b) requires both
#   programmes to clear the BH gate independently. ~80-gene sets: single-state
#   sensitivity ~0.78, hybrid recall ~0.54. ~25-gene sets: ~0.50 and ~0.22,
#   falling to 0.008 at a 4x effect. False positives never exceeded 0.53%
#   against a nominal 5% at any set size. So positives are trustworthy,
#   negatives are not: report small-set hybrid counts as a FLOOR and never
#   compare hybrid frequency between the two sets.
# =============================================================================
# WHAT CHANGED FROM v1 (all changes are in Step 5 / Step 6 / Step 6b; the
# loading, DE, plotting, table and GeoJSON sections are structurally unchanged)
#
#   1. modulescore_perm_pval() -> sig_enrich_perm().  The old null drew random
#      gene sets UNIFORMLY from all genes.  In a sparse matrix most genes are
#      near-undetected, so the null score had a small across-cell variance that
#      was IDENTICAL for every signature tested, and any well-detected marker
#      signature cleared it far too easily.  On simulated data with NO signature
#      structure the old null called 26.7% of cells at p<0.05 (nominal 5%);
#      matching the null on expression AND detection rate returns it to 3.4%.
#
#   2. Decisions are made on `z` (score standardised against its own matched
#      null), not on the empirical p.  The empirical p floors at 1/(n_perm+1),
#      so larger signatures saturate the floor more readily at equal effect and
#      which.min(pval) silently breaks ties by COLUMN ORDER.  z is continuous
#      and comparable across signatures of different sizes.
#
#   3. Significance is gated on `padj_z` (BH on the normal-tail p of z), NOT on
#      `padj` (BH on the empirical p).  BH at n cells needs the smallest p below
#      ALPHA/n; the empirical p cannot go below 1/(n_perm+1), so at n_perm=1000
#      and >20,000 cells NO cell can ever be called at any effect size.  You
#      would need n_perm > n_cells/ALPHA.  The z-tail p is continuous and BH
#      behaves sensibly on it.
#
#      The gate is a real trade-off and is exposed as --gate. On simulated data
#      (300 tumour cells, 4 true states, FC=3, 6 signatures, n_perm=300):
#          gate                  % assigned   accuracy among assigned
#          pval_z  (uncorrected)     42.3%          74.0%
#          padj_z  (BH, DEFAULT)     11.0%          93.9%
#          pval    (empirical)       34.0%          81.4%
#          padj    (BH on empirical)  0.0%             --   <- cannot ever fire
#      padj_z is the default because a state call that is wrong 26% of the time
#      is worse than no call. If your biological question needs coverage more
#      than precision, --gate pval_z is defensible, but report which you used:
#      with n cells the BH threshold is ALPHA/n at its strictest, so padj_z
#      necessarily gets stricter as the sample grows while pval_z does not.
#
#   4. expr_ref (the cells used to build the matching strata) is SELF-reference
#      everywhere: Neftel strata from TdTom+, cell-type strata from TdTom-.
#      Do NOT cross the populations.  A tumour-restricted gene is near-zero in
#      TdTom- cells, so a TdTom- reference files it in a bottom stratum, draws
#      near-zero matched controls, and the signature beats its own null in every
#      cell.  Simulated, on genes uniformly ON in every tumour cell with NO
#      cell-to-cell structure (right answer ~5%):
#          strata from TdTom+ (self)    1.0% of cells called
#          strata from all cells       82.8%
#          strata from TdTom- (cross) 100.0%
#      The usual no-signal calibration check does not catch this (3.9/5.6/5.1%)
#      because it only bites for POPULATION-RESTRICTED signatures — which is
#      exactly what tumour-state and cell-type signatures are.
#
#   5. Each signature gets its own RNG stream (seed + signature-name hash).  v1
#      passed the same `seed` to every call, so every signature was compared
#      against the SAME random gene sets; p-values across signatures were
#      perfectly dependent (cor = 1.0 between null means), which undermines
#      precisely the cross-signature comparisons the hybrid caller performs.
#
#   6. EC sub-states are tested on their UNIQUE genes and reported
#      INDEPENDENTLY (a cell may be positive for none, one, or several).
#      Winner-takes-all across signatures sharing a pan-EC core assigned only
#      22% of cells correctly in simulation when the truth was 100% one state —
#      worse than chance among three, because the shared core fires in every
#      cell and decides the ranking.
#
#   7. New Step 5c: null-calibration diagnostics written to CSV.  A gene-label
#      shuffle (destroys co-expression, preserves library size) should call
#      ~ALPHA of cells; a tumour-restricted-but-not-state-specific control
#      signature should NOT fire in nearly every TdTom+ cell.
#
#   8. The hybrid caller now receives z-scores.  Its criterion-(d) delta is
#      therefore dimensionless, which removes the "check delta against your
#      score scale" caveat that v1's comments correctly flagged.
#
# Usage:
#   Rscript 2b3_neftel_classify_hybrids_v7.R \
#     --rds           /path/to/sample_classified.rds \
#     --neftel        /path/to/neftel_209_full.csv \
#     --tdtomato      /path/to/rescued_tdtomato.csv \
#     --out_dir       /path/to/output_dir \
#     --ec_signatures /path/to/ec_signatures.csv \  # optional
#     --lib_path      /path/to/R/libs \             # optional
#     --use_predicted /path/to/t0_validation_meta_t2model.rds \  # optional
#     --n_perm        1000 \
#     --alpha         0.05 \
#     --seed          123
#
#   Invoked exactly as above, v3 behaves identically to v2: Step 6c is skipped
#   unless at least one myeloid signature file is supplied.
#
# Usage WITH myeloid state sub-typing (v3, Step 6c):
#   Rscript 2b3_neftel_classify_hybrids_v7.R \
#     --rds               /path/to/sample_classified.rds \
#     --neftel            /path/to/neftel_209_full.csv \
#     --tdtomato          /path/to/rescued_tdtomato.csv \
#     --out_dir           /path/to/output_dir \
#     --ec_signatures     /path/to/ec_signatures.csv \             # optional
#     --miller_signatures /path/to/miller_2025_microglia.csv \     # enables Step 6c
#     --qian_signatures   /path/to/qian_2026_microglia.csv \       # enables Step 6c
#     --n_perm            1000 \
#     --alpha             0.05 \
#     --seed              123
#
#   Either signature file alone enables Step 6c; supplying both additionally
#   produces the cross-set agreement table. Signature CSVs are one column per
#   state, gene symbols down the rows, ragged columns padded with blanks.
#   Miller columns are read as HUMAN symbols (converted via homologene);
#   Qian columns as MOUSE symbols (intersected only). A UTF-8 BOM is tolerated.
#
# Full option reference (defaults in brackets):
#   REQUIRED
#     --rds                     classified Seurat RDS
#     --neftel                  Neftel signature CSV
#     --out_dir                 output directory
#     --tdtomato                rescued TdTomato CSV (Step 1b fails without it)
#   OPTIONAL — inputs
#     --ec_signatures           EC sub-state CSV; enables Step 6b        [NULL]
#     --miller_signatures       Miller 2025 state CSV; enables Step 6c   [NULL]
#     --qian_signatures         Qian 2026 state CSV; enables Step 6c     [NULL]
#     --control_signature       tumour-restricted, non-state-specific genes;
#                               positive control for the expr_ref failure [NULL]
#     --use_predicted           predictions RDS for the TdTomato label   [NULL]
#     --lib_path                prepended to .libPaths()                 [NULL]
#   OPTIONAL — statistics
#     --n_perm                  permutations per signature               [1000]
#     --alpha                   significance level for every BH gate     [0.05]
#     --seed                    RNG seed                                 [123]
#     --nbin                    expression bins in the matched null      [25]
#     --n_det_bin               detection-rate bins in the matched null  [5]
#     --gate                    'padj_z' or 'pval_z'; 'padj' is refused  [padj_z]
#   OPTIONAL — Neftel hybrids
#     --hybrid_delta_q          quantile for the auto criterion-(d) delta [0.95]
#     --hybrid_pct_floor        min fraction of cells for a hybrid pair   [0.10]
#   OPTIONAL — EC sub-states (Step 6b)
#     --ec_unique_genes         score only genes unique to each sub-state [TRUE]
#     --ec_min_unique           below this, fall back to the full set     [5]
#   OPTIONAL — myeloid states (Step 6c)
#     --myeloid_min_cells       skip Step 6c below this many myeloid cells [30]
#     --myeloid_min_genes       skip a state with fewer testable genes     [3]
#     --myeloid_strip_identity  drop identity genes from state sets     [TRUE]
#     --myeloid_hybrid_delta    '0' disables criterion (d); 'auto' is
#                               Neftel-identical (see item B above)        ["0"]
#   OPTIONAL — diagnostics
#     --calibration_check       run the null-calibration diagnostics    [TRUE]
#     --calib_n_cells           cells sampled for calibration            [500]
#
# Recommended FIRST run on a new object — cheap, exercises every new code path:
#   Rscript 2b3_neftel_classify_hybrids_v7.R --rds subset.rds \
#     --neftel neftel_209_full.csv \
#     --tdtomato rescued_tdtomato.csv --out_dir results/smoketest \
#     --miller_signatures miller_2025_microglia.csv \
#     --qian_signatures qian_2026_microglia.csv \
#     --n_perm 100 --calibration_check FALSE
#   Then check the log for "Collapsed to family 'qian26_MG-Prolif'". If that line
#   is absent or warns, the state names in the CSV differ from MYELOID_COLLAPSE.
#
# Outputs (<case> = sample id parsed from the RDS filename):
#   <case>_classified.rds                     object with all new metadata columns
#   <case>_percell_states.csv                 per-cell calls + decision statistics
#   <case>_cell_summary.csv                   grouped counts
#   <case>_enrichment_stats.csv               per-signature enrichment
#   <case>_hybrid_pairs.csv                   Neftel hybrid pair counts
#   <case>_null_calibration.csv               calibration diagnostics
#   <case>_run_parameters.csv                 resolved parameters for this run
#   <case>_DE_tdtom_pos_vs_neg.csv            DE table, TdTom+ vs TdTom-
#   <case>_classification_plots.pdf           all figures, incl. section 8g
#   <case>_ec_substates.csv                   Step 6b, if --ec_signatures
#   <case>_myeloid_state_gene_audit.csv       Step 6c: published -> tested per state,
#                                             and which identity genes were stripped
#   <case>_myeloid_states_<set>.csv           Step 6c: per-cell scores/z/gate/p/calls
#   <case>_myeloid_state_summary.csv          Step 6c: positivity per state
#   <case>_myeloid_<set>_state_collinearity.csv  gene sharing on genes ACTUALLY tested
#   <case>_{neftel_normal,proliferation,intersection,state_label,tdtomato,
#           ec_state}.geojson                 spatial exports
#   NOTE: no GeoJSON is written for myeloid state — that plot is in the PDF only.
# =============================================================================

suppressPackageStartupMessages(library(optparse))

opt_list <- list(
  make_option("--rds",           type = "character", help = "Path to classified RDS"),
  make_option("--neftel",        type = "character", help = "Path to neftel_209_full.csv"),
  make_option("--tdtomato",      type = "character", help = "Path to new tdtomato"),
  make_option("--out_dir",       type = "character", help = "Output directory"),
  make_option("--ec_signatures", type = "character", default = NULL,
              help = "Path to EC signatures CSV (mouse genes; columns = signature names)"),
  make_option("--lib_path",      type = "character", default = NULL),
  make_option("--use_predicted", type = "character", default = NULL,
              help = "Path to t0_validation_meta_t2model.rds — use predicted TdTomato labels"),
  make_option("--n_perm",  type = "integer", default = 1000L),
  make_option("--alpha",   type = "double",  default = 0.05),
  make_option("--seed",    type = "integer", default = 123L),
  make_option("--nbin",       type = "integer", default = 25L,
              help = "Mean-expression strata for the matched null [25]"),
  make_option("--n_det_bin",  type = "integer", default = 5L,
              help = "Detection-rate strata for the matched null [5]"),
  make_option("--auto_strata", type = "logical", default = TRUE,
              help = paste("Derive nbin/n_det_bin from the number of genes present so each",
                           "stratum holds >=25 control genes. Targeted panels need this;",
                           "the 25x5 default is a whole-transcriptome setting. [TRUE]")),
  make_option("--strata_target", type = "integer", default = 10L,
              help = paste("Control genes per stratum targeted by --auto_strata. nbin is",
                           "capped at 25, so whole-transcriptome input reproduces the v6",
                           "default and only panels are coarsened. [10]")),
  make_option("--borrow_strata", type = "logical", default = TRUE,
              help = paste("Widen the matching window instead of ever resampling the",
                           "signature genes into their own null (v6 bug). [TRUE]")),
  make_option("--p_source", type = "character", default = "z",
              help = paste("Which p feeds the gate: 'z' (normal tail of the z score — the",
                           "default) or 'emp' (the permutation p directly). The normal tail is",
                           "anti-conservative: on 7000/7003_Scan2 it reports p~2e-4 where the",
                           "permutations measure ~0.017. 'emp' is honest but LOWERS yield to",
                           "roughly v6 levels, because part of v6's yield rested on that",
                           "overstatement -- it does NOT reduce the unassigned pile (that is",
                           "--gate_unit's job). Use 'emp' when a p-value itself is going into",
                           "a figure or a claim; ranking is unaffected either way. [z]")),
  make_option("--gate_unit", type = "character", default = "bh_cell",
              help = paste("Multiplicity unit for the significance gate:",
                           "bh_cell (BH across the competing SIGNATURES within each cell —",
                           "the default: the hypotheses in competition for one cell are its",
                           "own signatures, and this is what fixes the unassigned pile:",
                           "Neftel 24.4/17.1% -> 38.6/32.6% on 7000/7003_Scan2),",
                           "bh_global (BH across CELLS per signature — v6 behaviour; also",
                           "gives each signature its own adaptive threshold, which biases",
                           "recovered composition), raw (uncorrected). [bh_cell]")),
  make_option("--skew_p", type = "logical", default = FALSE,
              help = "Cornish-Fisher skew correction on the z-tail p. STRICTER. [FALSE]"),
  make_option("--pool_k", type = "integer", default = 0L,
              help = paste("Average expression over each cell + k nearest spatial",
                           "neighbours before scoring. 0 disables. 6-10 if used —",
                           "check Moran s I first. [0]")),
  make_option("--pool_coords", type = "character", default = NULL,
              help = paste("Comma-separated meta.data columns holding spatial",
                           "coordinates, e.g. x_centroid,y_centroid. Required by --pool_k.")),
  make_option("--tier_z_floor", type = "double", default = 1.0,
              help = "Top z above this, with no gate passed, is reported as tier weak [1.0]"),
  make_option("--hybrid_delta_q",   type = "double", default = 0.95,
              help = "Quantile of the shuffled-null gap_23 used as criterion-(d) delta; data-driven, dimensionless [0.95]"),
  make_option("--hybrid_pct_floor", type = "double", default = 0.10,
              help = "2nd program must exceed this quantile of its own top-cells (Neftel crit 2) [0.10]"),
  make_option("--gate", type = "character", default = "padj_z",
              help = paste("Significance gate: 'padj_z' (BH on z-tail p, FDR across cells,",
                           "conservative) or 'pval_z' (uncorrected z-tail p, more calls,",
                           "more false positives). 'padj' (BH on the empirical p) is",
                           "REJECTED — it cannot fire, see header note. [padj_z]")),
  make_option("--ec_unique_genes", type = "logical", default = TRUE,
              help = "Test EC sub-states on genes unique to each signature [TRUE]"),
  make_option("--ec_min_unique",   type = "integer", default = 5L,
              help = "Minimum unique genes required to test an EC sub-state [5]"),
  make_option("--calibration_check", type = "logical", default = TRUE,
              help = "Run the gene-label-shuffle null calibration diagnostics [TRUE]"),
  make_option("--calib_n_cells", type = "integer", default = 500L,
              help = "Cells subsampled for the calibration check [500]"),
  make_option("--control_signature", type = "character", default = NULL,
              help = "Optional CSV/txt of genes that are tumour-restricted but NOT state-specific (one column). Used as a positive control for the expr_ref failure mode."),

  # --- v3: myeloid (microglia / macrophage) STATE sub-typing -----------------
  make_option("--miller_signatures", type = "character", default = NULL,
              help = "CSV of Miller 2025 microglial state signatures (one column per state, HUMAN symbols). Enables Step 6c."),
  make_option("--qian_signatures", type = "character", default = NULL,
              help = "CSV of Qian 2026 microglial state signatures (one column per state, MOUSE symbols). Enables Step 6c."),
  make_option("--myeloid_min_cells", type = "integer", default = 30L,
              help = "Minimum myeloid cells required to run Step 6c [default 30]"),
  make_option("--myeloid_min_genes", type = "integer", default = 3L,
              help = "Minimum testable genes for a state to be scored [default 3]"),
  make_option("--myeloid_strip_identity", type = "logical", default = TRUE,
              help = paste("Remove genes that are also members of the myeloid IDENTITY signatures",
                           "(Microglia_core / Macrophage_peripheral / Macrophage_BAM) from the STATE",
                           "sets. Those genes selected the compartment, so leaving them in makes the",
                           "state test partly re-test its own selection criterion [default TRUE]")),
  make_option("--myeloid_hybrid_delta", type = "character", default = "0",
              help = paste("Criterion-(d) delta for myeloid hybrids: '0' (default, criterion disabled)",
                           "or 'auto' (Neftel-identical calibration). Simulation on this pipeline:",
                           "auto recovers 15% of true hybrids, 0 recovers 43%, false positives 0 in",
                           "both. See Step 6c header [default 0]")),

  # --- Step 2b: sub-cluster rescue (see the Step 2b header for the rationale) -
  make_option("--subcluster_rescue", type = "logical", default = FALSE,
              help = paste("Step 2b: relabel TdTom_neg cells as TdTom_pos when they sit in a",
                           "sub-cluster the marker scorer tiered as tumour-like. OFF by default:",
                           "with it off this script's output is unchanged [FALSE]")),
  make_option("--subcluster_units", type = "character", default = NULL,
              help = "tumourlike_units.csv from 2b4_tdtom_tumourlike_score.py (required if --subcluster_rescue)"),
  make_option("--subcluster_cells", type = "character", default = NULL,
              help = paste("<case>_subcluster_cells.csv (cell_id -> sub_id_full), from",
                           "2b2_extract_subcluster_cells.R or a v5+ enrichment run (required if",
                           "--subcluster_rescue)")),
  make_option("--rescue_tiers", type = "character",
              default = "A_confident_tumour,B_probable_tumour",
              help = paste("Comma-separated tiers eligible for rescue. Tier is assigned by",
                           "2b4_tdtom_tumourlike_score.py from marker panels, NOT from reporter rate",
                           "[A_confident_tumour,B_probable_tumour]")),
  make_option("--rescue_suspect_tiers", type = "character",
              default = "B_probable_tumour,C_ambiguous,D_weak",
              help = paste("Rescued tiers to report as 'suspect' in the counts file, to",
                           "show how much of the increment rests on weaker",
                           "evidence. Intersected with --rescue_tiers",
                           "[B_probable_tumour,C_ambiguous,D_weak]"))
)

opt <- parse_args(OptionParser(option_list = opt_list))
for (arg in c("rds", "neftel", "out_dir"))
  if (is.null(opt[[arg]])) stop("Missing required argument: --", arg)
if (!is.null(opt$lib_path)) .libPaths(opt$lib_path)

suppressPackageStartupMessages({
  library(Seurat)
  library(homologene)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(jsonlite)
})

options(future.globals.maxSize = Inf)

RDS_PATH      <- opt$rds
NEFTEL_CSV    <- opt$neftel
TDTOMATO_PATH <- opt$tdtomato
OUT_DIR       <- opt$out_dir
EC_SIG_CSV    <- opt$ec_signatures
USE_PREDICTED <- opt$use_predicted
N_PERM        <- opt$n_perm
ALPHA         <- opt$alpha
SEED          <- opt$seed
NBIN          <- opt$nbin
N_DET_BIN     <- opt$n_det_bin
GATE          <- opt$gate
AUTO_STRATA   <- isTRUE(opt$auto_strata)
STRATA_TARGET <- opt$strata_target
BORROW_STRATA <- isTRUE(opt$borrow_strata)
SKEW_P        <- isTRUE(opt$skew_p)
POOL_K        <- opt$pool_k
POOL_COORDS   <- if (is.null(opt$pool_coords)) NULL else
                   trimws(strsplit(opt$pool_coords, ",", fixed = TRUE)[[1]])
TIER_Z_FLOOR  <- opt$tier_z_floor
options(neftel.borrow_strata = BORROW_STRATA, neftel.skew_p = SKEW_P)
GATE_UNIT     <- opt$gate_unit
if (!GATE_UNIT %in% c("bh_global", "bh_cell", "raw"))
  stop("--gate_unit must be bh_global, bh_cell or raw (got '", GATE_UNIT, "')")
P_SOURCE <- opt$p_source
if (!P_SOURCE %in% c("z", "emp"))
  stop("--p_source must be z or emp (got '", P_SOURCE, "')")
# Resolution note. An empirical p cannot be smaller than 1/(n_perm+1), which
# limits how deep a gate can reach -- but NOT in the way an earlier version of
# this guard assumed. That version compared the floor against alpha/m and
# refused to run bh_global + emp. That was WRONG: alpha/m is the BH threshold
# for the SMALLEST p only. BH compares the i-th ordered p against i*alpha/m, so
# the usable threshold grows with the number of true signals. Measured on real
# whole-transcriptome data (2 samples, n_perm = 1000, m ~ 5e4 cells), bh_global
# + emp fired on 10.1% / 4.6% of cells for cell types and 2.1% / 0.0% for
# Neftel states -- degraded, not dead.
#
# So: report the arithmetic, warn when ties at the floor will dominate, and let
# the run proceed. The user decides.
if (P_SOURCE == "emp") {
  .p_floor <- 1 / (N_PERM + 1)
  message(sprintf("--p_source emp: permutation floor 1/(n_perm+1) = %.2e.", .p_floor))
  if (GATE_UNIT == "bh_cell") {
    # Denominator is the number of signatures tested per cell: small and known.
    message(sprintf(paste0("  bh_cell needs about %.2e (alpha/#sigs); floor gives ",
                           "%.1fx headroom."), ALPHA / 20, (ALPHA / 20) / .p_floor))
  } else if (GATE_UNIT == "bh_global") {
    warning("--p_source emp with --gate_unit bh_global is resolution-limited: many cells\n",
            "  tie at the ", sprintf("%.2e", .p_floor), " floor and BH breaks ties ",
            "conservatively, so\n  yield drops sharply (measured: 2-10% vs 17-44% under ",
            "--p_source z).\n  --gate_unit bh_cell is the better pairing for emp; raise ",
            "--n_perm if you\n  need bh_global specifically.", call. = FALSE)
  }
}
# --gate is retained for backward compatibility and mapped onto --gate_unit.
if (identical(GATE, "pval_z") && identical(GATE_UNIT, "bh_global")) {
  GATE_UNIT <- "raw"
  message("NOTE: --gate pval_z maps to --gate_unit raw.")
}
if (POOL_K > 0L && is.null(POOL_COORDS))
  stop("--pool_k requires --pool_coords (e.g. --pool_coords x_centroid,y_centroid)")
if (POOL_K > 0L)
  message("WARNING: neighbour pooling is ON (k=", POOL_K, "). Pooled labels must not ",
          "be used as evidence that neighbouring cells share a state, and the ",
          "headline result should also be reported at --pool_k 0.")
if (identical(GATE, "padj"))
  stop("--gate padj is not permitted: BH on the empirical permutation p cannot ",
       "fire at feasible n_perm (the floor 1/(n_perm+1) exceeds ALPHA/n_cells). ",
       "Use 'padj_z' (default) or 'pval_z'.")
if (!GATE %in% c("padj_z", "pval_z"))
  stop("--gate must be 'padj_z' or 'pval_z' (got '", GATE, "')")
EC_UNIQUE     <- isTRUE(opt$ec_unique_genes)
EC_MIN_UNIQUE <- opt$ec_min_unique
DO_CALIB      <- isTRUE(opt$calibration_check)
CALIB_N_CELLS <- opt$calib_n_cells
CONTROL_SIG   <- opt$control_signature

# --- v3: myeloid state sub-typing ------------------------------------------
MILLER_CSV       <- opt$miller_signatures
QIAN_CSV         <- opt$qian_signatures
USE_MYELOID      <- !is.null(MILLER_CSV) || !is.null(QIAN_CSV)
MYELOID_MIN_CELLS <- opt$myeloid_min_cells
MYELOID_MIN_GENES <- opt$myeloid_min_genes
MYELOID_STRIP    <- isTRUE(opt$myeloid_strip_identity)
MYELOID_DELTA    <- if (identical(tolower(opt$myeloid_hybrid_delta), "auto")) "auto" else
                    suppressWarnings(as.numeric(opt$myeloid_hybrid_delta))
if (!identical(MYELOID_DELTA, "auto") && (is.na(MYELOID_DELTA) || MYELOID_DELTA < 0))
  stop("--myeloid_hybrid_delta must be 'auto' or a non-negative number")

# --- Step 2b: sub-cluster rescue -------------------------------------------
SUBCLUSTER_RESCUE <- isTRUE(opt$subcluster_rescue)
SUBCLUSTER_UNITS  <- opt$subcluster_units
SUBCLUSTER_CELLS  <- opt$subcluster_cells
RESCUE_TIERS      <- trimws(strsplit(opt$rescue_tiers, ",")[[1]])
RESCUE_TIERS      <- RESCUE_TIERS[nzchar(RESCUE_TIERS)]
RESCUE_SUSPECT    <- trimws(strsplit(opt$rescue_suspect_tiers, ",")[[1]])
RESCUE_SUSPECT    <- intersect(RESCUE_SUSPECT[nzchar(RESCUE_SUSPECT)], RESCUE_TIERS)
if (SUBCLUSTER_RESCUE) {
  if (is.null(SUBCLUSTER_UNITS) || is.null(SUBCLUSTER_CELLS))
    stop("--subcluster_rescue needs both --subcluster_units and --subcluster_cells")
  for (p in c(SUBCLUSTER_UNITS, SUBCLUSTER_CELLS))
    if (!file.exists(p)) stop("File not found: ", p)
  if (!length(RESCUE_TIERS)) stop("--rescue_tiers is empty")
  # Fail on a typo rather than silently rescuing nothing: an unrecognised tier
  # name would otherwise match no rows and look like "no cells qualified".
  known_tiers <- c("A_confident_tumour", "B_probable_tumour", "C_ambiguous",
                   "D_weak", "E_non_tumour")
  if (length(setdiff(RESCUE_TIERS, known_tiers)))
    stop("Unknown tier(s) in --rescue_tiers: ",
         paste(setdiff(RESCUE_TIERS, known_tiers), collapse = ", "),
         "\n  Known tiers: ", paste(known_tiers, collapse = ", "))
}

# The myeloid compartment. Cells with these cell_class values are sub-typed in
# Step 6c. Identity is settled in Step 6; Step 6c only asks which STATE they are in.
MYELOID_CLASSES <- c("Microglia_core", "Macrophage_peripheral", "Macrophage_BAM")

# --- Hybrid-calling parameters (tumour cells only; see assign_states_with_hybrids) ---
# The criterion-(d) delta (min 2nd-vs-3rd score gap) is NOT hard-coded. It is
# calibrated from the data by shuffling each state's scores independently
# (Neftel's expected-hybrid null) and taking a high quantile of the resulting
# gap_23 distribution. In v2 the caller receives Z-SCORES, so delta is already
# dimensionless — the auto-calibration is still used (it adapts to how many
# states there are and how correlated they are), but the v1 caveat about score
# units no longer applies.
HYBRID_DELTA_Q   <- if (!is.null(opt$hybrid_delta_q))   opt$hybrid_delta_q   else 0.95
HYBRID_PCT_FLOOR <- if (!is.null(opt$hybrid_pct_floor)) opt$hybrid_pct_floor else 0.10

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

message("=== Parameters ===")
message("  RDS:           ", RDS_PATH)
message("  OUT_DIR:       ", OUT_DIR)
message("  EC_SIG_CSV:    ", ifelse(is.null(EC_SIG_CSV), "None", EC_SIG_CSV))
message("  USE_PREDICTED: ", ifelse(is.null(USE_PREDICTED), "No (gene expression)", USE_PREDICTED))
message("  N_PERM:        ", N_PERM)
message("  ALPHA:         ", ALPHA)
message("  Null:          expression x detection matched, ", NBIN, " x ", N_DET_BIN, " strata")
message("  Gate:          padj_z (BH on z-tail p)   Ranking: z")
message("  Myeloid states: ", if (!USE_MYELOID) "None" else
        paste0(paste(c(if (!is.null(MILLER_CSV)) "Miller", if (!is.null(QIAN_CSV)) "Qian"),
                     collapse = " + "),
               "  (strip_identity=", MYELOID_STRIP,
               ", hybrid delta=", if (identical(MYELOID_DELTA,"auto")) "auto" else MYELOID_DELTA, ")"))

# =============================================================================
# Signature definitions
# =============================================================================

neftel_csv       <- read.csv(NEFTEL_CSV)
neftel_sig_human <- lapply(colnames(neftel_csv), function(state) {
  genes <- neftel_csv[[state]]
  genes[!is.na(genes) & genes != ""]
})
names(neftel_sig_human) <- colnames(neftel_csv)

prolif_cols <- c("G1.S", "G2.M")
missing_prolif <- setdiff(prolif_cols, names(neftel_sig_human))
if (length(missing_prolif) > 0)
  stop("Proliferation columns not found in neftel CSV: ",
       paste(missing_prolif, collapse = ", "),
       "\nAvailable columns: ", paste(names(neftel_sig_human), collapse = ", "))

prolif_sig_human <- neftel_sig_human[prolif_cols]
neftel_sig_human <- neftel_sig_human[setdiff(names(neftel_sig_human), prolif_cols)]

message("Neftel states: ", paste(names(neftel_sig_human), collapse = ", "))
message("Prolif states: ", paste(prolif_cols, collapse = ", "))

# EC signatures — mouse genes, read directly if provided
USE_EC <- !is.null(EC_SIG_CSV)
if (USE_EC) {
  ec_csv <- read.csv(EC_SIG_CSV, stringsAsFactors = FALSE)
  ec_sig_mouse_raw <- lapply(colnames(ec_csv), function(sig) {
    genes <- ec_csv[[sig]]
    genes[!is.na(genes) & genes != ""]
  })
  names(ec_sig_mouse_raw) <- colnames(ec_csv)
  message("EC signatures loaded: ", paste(names(ec_sig_mouse_raw), collapse = ", "))
  message("EC signature sizes (raw): ", paste(sapply(ec_sig_mouse_raw, length), collapse = ", "))
} else {
  message("No EC signatures provided — skipping EC scoring.")
}

celltype_sig_mouse <- list(
  # Pan-neuronal core (merged with former Hypothalamic_neuron; region-specific
  # neuropeptide subtypes removed -- averaging Avp/Oxt/Npy/Agrp/Pomc/Sst/Pmch/Hcrt into
  # one mean diluted true single-subtype-positive cells, see specificity review)
  Neuron = c("Rbfox3", "Snap25", "Syt1", "Map2", "Tubb3", "Nefl", "Nefm", "Syp", "Camk2a", "Thy1", "Gabra1", "Grin2b"),
  # Six3 dropped (developmental patterning TF, not an established adult MSN marker)
  Striatal_MSN = c("Ppp1r1b", "Rgs9", "Gpr88", "Pde10a", "Penk", "Scn4b", "Rasd2", "Pde1b", "Arpp21", "Gng7", "Rasgrp2", "Adcy5", "Camkv", "Kcnab1", "Drd1", "Drd2", "Foxp1", "Foxp2", "Bcl11b"),
  # Gsn dropped (broadly expressed, not lineage-restricted)
  Oligodendro = c("Mbp", "Plp1", "Mog", "Cnp", "Mag", "Mobp", "Sox10", "Cldn11", "Opalin", "Fa2h", "Ermn"),
  # Vim dropped (canonical EMT/mesenchymal marker, elevated in MES-like GBM tumor
  # cells -- direct risk of misclassifying tumor cells as normal astrocytes).
  # Id3 dropped (weak, non-specific bHLH TF)
  Astrocyte = c("Gfap", "Aqp4", "Aldh1l1", "S100b", "Slc1a3", "Glul", "Ndrg2", "Cldn10", "Hopx", "Aldoc"),
  # unchanged -- no specificity flags
  Microglia_core = c("P2ry12", "Tmem119", "Fcrls", "Slc2a5", "Olfml3", "Gpr34", "Siglech", "Sall1", "Hexb"),
  # S100a4 dropped (canonical fibroblast/EMT marker, same MES-like risk as Vim)
  Macrophage_peripheral = c("Emilin2", "Gda", "Hp", "Sell", "Ly6c2", "Ccr2", "Plac8"),
  # unchanged -- no specificity flags
  Macrophage_BAM = c("Mrc1", "Ms4a7", "Cd163", "Lyve1", "Clec12a", "Ptgs1", "Cd36"),
  # State/subtype markers removed: Gzma, Gzmb, Nkg7 (cytotoxic-effector state),
  # Foxp3, Il2ra (Treg/activation state) -- these belong in a STATE signature, not
  # identity, mirroring the myeloid identity/state split. Ptprc (CD45) dropped --
  # pan-leukocyte, shared with the entire immune compartment, not T-cell-specific.
  # Left: core TCR complex + canonical coreceptors only
  Tcell = c("Cd3e", "Cd3d", "Cd3g", "Cd8a", "Cd4", "Trac", "Trbc1"),
  # Reduced to pan-EC lineage anchors only (see EC signature dilution analysis) --
  # BBB-transporter genes moved out to keep tumor-associated EC detectable
  Endothelial = c("Pecam1", 'Cldn5', 'Igfbp7', 'Itm2a', 'Ptprb', 'Egfl7',"Cdh5", "Kdr", "Flt1", "Erg", "Vwf", "Esam", "Icam2"),
  # S100a6 dropped (broadly expressed S100-family gene, not ependymal-restricted)
  Ependymal = c("Ccdc153", "Tmem212", "Foxj1", "Rsph1", "Aqp1", "Hydin", "Cfap43", "Dnai1", "Pifo"),
  # unchanged -- Olig1/Olig2 overlap with Neftel OPC-like tumor state is by design
  # (TdTom- vs TdTom+ compartment split)
  OPC = c("Pdgfra", "Cspg4", "Olig1", "Olig2", "Sox4", "Bcan", "Neu4", "Lhfpl3", "Epn2", "Cacng4"),
  # Acta2 dropped (smooth-muscle actin also marks myofibroblast/reactive stroma --
  # same MES-like tumor risk as Vim/S100a4)
  Pericyte = c("Pdgfrb", "Rgs5", "Des", "Kcnj8", "Abcc9", "Anpep", "Cd248", "Mcam", "Notch3")
)

# --- v3: myeloid STATE signatures (loaded from CSV) --------------------------
# These are STATE signatures, orthogonal to the identity signatures above.
# Identity ("is this cell a microglial cell?") is decided in Step 6 from
# celltype_sig_mouse; state ("which programme is it running?") in Step 6c.
#
# read_myeloid_csv() strips a UTF-8 byte-order mark from the header and the
# values. qian_2026_microglia.csv carries one: without stripping, the first
# state name arrives as "\ufeffMG-P2ry12" and silently fails every name match
# (including the collapse specification below) while still scoring under the
# mangled name.
read_myeloid_csv <- function(path) {
  csv <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  fileEncoding = "UTF-8-BOM")
  nms <- sub("^\ufeff", "", trimws(colnames(csv)))
  out <- lapply(colnames(csv), function(cn) {
    g <- sub("^\ufeff", "", trimws(as.character(csv[[cn]])))
    unique(g[!is.na(g) & nzchar(g) & g != "NA"])
  })
  names(out) <- nms
  out
}

# Collapse specification. Same principle as the Neftel MES1/MES2 -> MES_like
# collapse: states that share most of their genes cannot be ranked against each
# other, so they are reported as one family.
#
# Qian MG-CT-Prolif and TR-Prolif-MG share Jaccard 0.47 (9 and 5 unique genes).
# They are also almost entirely canonical cell-cycle genes, so they duplicate
# the genome-wide proliferation test in Step 6. Cross-tabulate the collapsed
# family against prolif_state rather than reading it as an independent finding.
#
# Members are matched PREFIX-TOLERANTLY by match_state_names() below: the CSV
# column names carry a study prefix ("qian26_MG-CT-Prolif"), so a spec written
# with bare state names would silently match nothing and the collapse would be
# skipped without error.
MYELOID_COLLAPSE <- list(
  qian = list(`MG-Prolif` = c("MG-CT-Prolif", "TR-Prolif-MG"))
)

# Match spec names against the actual state names: exact first, then ignoring a
# leading study prefix ("qian26_"), then case-insensitively. Returns the ACTUAL
# names. Warns on anything unmatched, because a silently-skipped collapse looks
# identical to a collapse that was never specified.
match_state_names <- function(wanted, actual) {
  bare <- sub("^[A-Za-z]+[0-9]*_", "", actual)
  hit  <- vapply(wanted, function(w) {
    i <- which(actual == w)
    if (!length(i)) i <- which(bare == w)
    if (!length(i)) i <- which(tolower(bare) == tolower(w))
    if (length(i)) actual[i[1]] else NA_character_
  }, character(1), USE.NAMES = FALSE)
  if (any(is.na(hit)))
    warning("Collapse spec names not found among the state names (collapse will be ",
            "incomplete): ", paste(wanted[is.na(hit)], collapse = ", "),
            " | available: ", paste(actual, collapse = ", "))
  hit[!is.na(hit)]
}

myeloid_sig_raw <- list()
if (!is.null(MILLER_CSV)) myeloid_sig_raw$miller <- read_myeloid_csv(MILLER_CSV)
if (!is.null(QIAN_CSV))   myeloid_sig_raw$qian   <- read_myeloid_csv(QIAN_CSV)
# Miller signatures are HUMAN symbols and are converted in Step 3.
# Qian signatures are already MOUSE and are only intersected.
MYELOID_SPECIES <- c(miller = "human", qian = "mouse")

class_colours <- c(
  # collapsed Neftel families (tumour cells are labelled with these) ...
  MES_like            = "#E41A1C", NPC_like            = "#984EA3",
  # ... and the sub-state names (kept for primary_sub / legacy references)
  MES1_like           = "#E41A1C", MES2_like           = "#FF7F00",
  AC_like             = "#4DAF4A", OPC_like            = "#377EB8",
  NPC1_like           = "#984EA3", NPC2_like           = "#A65628",
  Neuron              = "#FFD700", Striatal_MSN        = "#8FBC8F",
  Hypothalamic_neuron = "#CD853F", OPC                 = "#6495ED",
  Oligodendro         = "#00CED1", Astrocyte           = "#DA70D6",
  Microglia_core      = "#FF69B4", Macrophage_peripheral = "#8B4513",
  Tcell               = "#32CD32", Endothelial         = "#9400D3",
  Ependymal           = "#20B2AA", Pericyte            = "#708090",
  not_assigned        = "lightgrey",
  # v3: was "#E41A1C" — identical to MES_like, so BAM cells and MES_like tumour
  # cells were indistinguishable in every plot and GeoJSON. Replaced with a
  # colour chosen to maximise the minimum CIEDE2000 distance to every other
  # entry under normal vision plus simulated deuteranopia/protanopia/tritanopia
  # (worst case dE00 19.7; the warm region was already the most crowded).
  Macrophage_BAM      = "#22534E"
)

# v3: NPC2_like ("#A65628") sits at dE00 7.5 from Macrophage_peripheral
# ("#8B4513") — the closest pair in this palette, and closer still under
# deuteranopia. Not changed here because it would break continuity with earlier
# figures; set NPC2_like = "#7F7759" if you want it fixed.

myeloid_state_colours <- c(multiple = "#5E5E5E", none = "lightgrey",
                           undefined = "lightgrey", not_tested = "grey85")

prolif_colours <- c(G1S = "#E41A1C", G2M = "#377EB8", undefined = "lightgrey")

ec_colours <- c(
  BEC = "#1B7837",
  CmEC = "#762A83",
  ImEC = "#E08214",
  IpEC = "#0072B2",
  multiple       = "#5E5E5E",   # v2: cells positive for >1 sub-state
  none           = "lightgrey",
  undefined      = "lightgrey"
)

# =============================================================================
# 1. Load Seurat object
# =============================================================================

message("\n=== Step 1: Loading Seurat object ===")
seu     <- readRDS(RDS_PATH)
seu     <- UpdateSeuratObject(seu)
case_id <- str_extract(basename(RDS_PATH), "[0-9]{4}_Scan[0-9]+")
message("  Sample: ", case_id, "  Cells: ", ncol(seu))
DefaultAssay(seu) <- "SCT"

# =============================================================================
# 1b. Update Seurat object TdTomato gating
# =============================================================================

message("\n=== Step 1b: Updating TdTomato metadata ===")

rescued_tdtomato_df <- read.csv(TDTOMATO_PATH, stringsAsFactors = FALSE)
rescued_tdtomato_df$cell_id <- as.character(rescued_tdtomato_df$cell_id)

seu_cell_ids <- rownames(seu@meta.data)

if (anyDuplicated(rescued_tdtomato_df$cell_id))
  stop("Duplicate cell_id values found in TdTomato CSV")

idx <- match(as.character(seu_cell_ids), rescued_tdtomato_df$cell_id)

message("  Seurat cells: ", length(seu_cell_ids))
message("  CSV cells:    ", nrow(rescued_tdtomato_df))
message("  Matched:      ", sum(!is.na(idx)))
message("  Unmatched:    ", sum(is.na(idx)))

if (anyNA(idx) || nrow(rescued_tdtomato_df) != length(seu_cell_ids))
  stop("TdTomato CSV must cover exactly the Seurat cells: ", sum(is.na(idx)),
       " missing in CSV, ", nrow(rescued_tdtomato_df), " CSV rows vs ",
       length(seu_cell_ids), " Seurat cells. Check the sample and source.")
if (!all(rescued_tdtomato_df$tdtom_label %in% c("TdTom_pos", "TdTom_neg")))
  stop("TdTomato CSV contains labels outside TdTom_pos/TdTom_neg")

new_cols <- setdiff(colnames(rescued_tdtomato_df), "cell_id")
for (col in new_cols) seu@meta.data[[col]] <- rescued_tdtomato_df[[col]][idx]

# =============================================================================
# 2. Define TdTomato+ / TdTomato- cells
# =============================================================================

message("\n=== Step 2: Defining TdTomato labels ===")

if (!is.null(USE_PREDICTED)) {

  message("  Using predicted labels from: ", USE_PREDICTED)
  pred_list <- readRDS(USE_PREDICTED)
  pred_meta <- pred_list[[case_id]]
  if (is.null(pred_meta))
    stop("Sample ", case_id, " not found in predicted meta RDS.")

  shared <- intersect(colnames(seu), rownames(pred_meta))
  if (length(shared) == 0)
    stop("No shared cell barcodes between Seurat object and predicted meta")

  seu$tdtom_label <- NA_character_
  seu$tdtom_label[shared] <- pred_meta[shared, "tdtom_predicted"]
  seu$tdtom_source <- "predicted"

} else {

  if (!"tdtom_label" %in% colnames(seu@meta.data))
    stop("tdtom_label not found in Seurat metadata after CSV merge")
  seu$tdtom_source <- "rescued_csv"

}

tumour_cells    <- rownames(seu@meta.data)[which(seu$tdtom_label == "TdTom_pos")]
nontumour_cells <- rownames(seu@meta.data)[which(seu$tdtom_label == "TdTom_neg")]

message("  TdTom+ (tumour):     ", length(tumour_cells))
message("  TdTom- (non-tumour): ", length(nontumour_cells))

if (length(tumour_cells) == 0) stop("No TdTomato+ cells found")

# =============================================================================
# 2b. Sub-cluster rescue  (--subcluster_rescue; OFF by default)
# =============================================================================
#
# WHAT IT DOES
#   Relabels TdTom_neg -> TdTom_pos for cells sitting in a sub-cluster that
#   2b4_tdtom_tumourlike_score.py tiered as tumour-like. It runs HERE, before
#   tumour_cells / nontumour_cells are consumed by anything, so Steps 4, 5, 6,
#   6b, 6c, 7 and 8 all see the updated split with no further changes.
#
# WHY IT IS DEFENSIBLE
#   The tier comes from MARKER PANELS (proliferation, OPC-like, NPC-like, ...),
#   not from the reporter rate. So this is not "call it positive because it
#   looks positive" - it is transcriptional evidence overriding a reporter that
#   is known to be lossy. Reporter dropout is the whole reason the rescue CSV
#   exists upstream.
#
# WHY IT IS ALSO CIRCULAR, AND WHAT TO DO ABOUT IT
#   Step 6 assigns Neftel states from signatures that OVERLAP the panels used to
#   tier the sub-cluster. A sub-cluster rescued for being OPC-like will then be
#   called OPC-like by Step 6, largely because that is why it was rescued.
#   State PROPORTIONS over the rescued set are therefore inflated by
#   construction toward whatever programme drove the rescue.
#
#   This is not a reason to avoid the rescue - it is a reason not to read state
#   composition off the combined set without splitting it. Every rescued cell is
#   marked with tdtom_rescued_subcluster = TRUE and the sub-cluster that did it,
#   so any downstream comparison can hold rescued and probe-positive cells apart.
#   Treat "state composition of rescued cells" as descriptive, never as evidence
#   that the rescued population IS that state.
#
# WHAT IT NEVER DOES
#   It never turns a TdTom_pos cell negative, and never gives a label to a cell
#   whose tdtom_label is NA. It is strictly one-directional: neg -> pos.
# =============================================================================

seu$tdtom_label_original       <- seu$tdtom_label
seu$tdtom_rescued_subcluster   <- FALSE
seu$tdtom_rescue_source        <- NA_character_
rescue_stats <- NULL

if (SUBCLUSTER_RESCUE) {

  message("\n=== Step 2b: Sub-cluster rescue ===")

  # sub_id_full and sub_id are IDS, not numbers: "7_1" and "7_10" are different
  # sub-clusters. Forcing character stops any reader-side type inference from
  # merging or reformatting them.
  units <- read.csv(SUBCLUSTER_UNITS, stringsAsFactors = FALSE,
                    colClasses = c(unit = "character"))
  cmap  <- read.csv(SUBCLUSTER_CELLS, stringsAsFactors = FALSE,
                    colClasses = c(cell_id = "character", sub_id_full = "character"))

  for (col in c("kind", "tier", "unit", "section"))
    if (!col %in% names(units)) stop("--subcluster_units lacks column '", col, "'")
  for (col in c("cell_id", "sub_id_full"))
    if (!col %in% names(cmap)) stop("--subcluster_cells lacks column '", col, "'")

  # The units table can hold several sections. Taking rows for other samples
  # would rescue cells by an id that means something different there, so the
  # section must be matched to this run's case_id explicitly.
  if (!case_id %in% unique(units$section))
    stop("No rows for section '", case_id, "' in ", basename(SUBCLUSTER_UNITS),
         "\n  Sections present: ", paste(unique(units$section), collapse = ", "))
  units <- units[units$section == case_id, ]

  # unit is "sub 0.3" in the scorer's display form; the join key is the id.
  # Accept both separators: older tables used "." before the id-collision fix.
  units$sub_id_full <- gsub("^sub +", "", units$unit)
  units$sub_id_full <- gsub("\\.", "_", units$sub_id_full)

  sel <- units[units$kind == "subcluster" & units$tier %in% RESCUE_TIERS, ]
  message("  Tiers eligible: ", paste(RESCUE_TIERS, collapse = ", "))
  message("  Sub-clusters selected: ", nrow(sel),
          if (nrow(sel)) paste0(" (", paste(sel$sub_id_full, collapse = ", "), ")") else "")
  if (!nrow(sel))
    stop("No sub-clusters matched --rescue_tiers in section ", case_id,
         ". Nothing to rescue; re-run without --subcluster_rescue if that is intended.")

  # A units row with no cells in the map means the two files came from different
  # runs - the mapping would be wrong rather than merely incomplete.
  missing_units <- setdiff(sel$sub_id_full, unique(cmap$sub_id_full))
  if (length(missing_units))
    stop("Sub-cluster(s) selected for rescue are absent from --subcluster_cells: ",
         paste(missing_units, collapse = ", "),
         "\n  The units table and the cell map are from different runs.")

  cand <- cmap$cell_id[cmap$sub_id_full %in% sel$sub_id_full]
  in_obj <- intersect(cand, rownames(seu@meta.data))
  if (!length(in_obj))
    stop("No cells from the selected sub-clusters are present in this Seurat object.",
         "\n  Check that --subcluster_cells came from the same --rds as this run.")
  if (length(in_obj) < length(cand))
    message("  NOTE ", length(cand) - length(in_obj), " mapped cell(s) absent from ",
            "this object (expected if the objects differ in filtering)")

  # The one-directional rule, stated as code: only currently-negative cells flip.
  flip <- in_obj[which(seu@meta.data[in_obj, "tdtom_label"] == "TdTom_neg")]

  n_before_pos <- length(tumour_cells)
  n_before_neg <- length(nontumour_cells)

  seu@meta.data[flip, "tdtom_label"]             <- "TdTom_pos"
  seu@meta.data[flip, "tdtom_rescued_subcluster"] <- TRUE
  seu@meta.data[flip, "tdtom_rescue_source"]     <-
    cmap$sub_id_full[match(flip, cmap$cell_id)]

  # Recompute the two vectors every later step reads.
  tumour_cells    <- rownames(seu@meta.data)[which(seu$tdtom_label == "TdTom_pos")]
  nontumour_cells <- rownames(seu@meta.data)[which(seu$tdtom_label == "TdTom_neg")]

  stopifnot(length(tumour_cells) == n_before_pos + length(flip),
            length(nontumour_cells) == n_before_neg - length(flip))

  # Report the three fates separately. Lumping NA in with "already positive"
  # would hide cells that carry no reporter call at all.
  lab_in <- seu@meta.data[in_obj, "tdtom_label_original"]
  message("  Cells in selected sub-clusters: ", length(in_obj))
  message("  Already TdTom+:                 ", sum(lab_in == "TdTom_pos", na.rm = TRUE))
  message("  Unlabelled (NA, left as NA):    ", sum(is.na(lab_in)))
  message("  RESCUED (neg -> pos):           ", length(flip))
  message("  TdTom+ ", n_before_pos, " -> ", length(tumour_cells),
          "  |  TdTom- ", n_before_neg, " -> ", length(nontumour_cells))

  per_unit <- as.data.frame(table(seu@meta.data[flip, "tdtom_rescue_source"]),
                            stringsAsFactors = FALSE)
  if (nrow(per_unit)) {
    names(per_unit) <- c("sub_id_full", "n_rescued")
    per_unit$tier   <- sel$tier[match(per_unit$sub_id_full, sel$sub_id_full)]
    per_unit$driver <- if ("driver" %in% names(sel))
      sel$driver[match(per_unit$sub_id_full, sel$sub_id_full)] else NA_character_
    message("  Per sub-cluster:")
    print(per_unit[order(-per_unit$n_rescued), ], row.names = FALSE)
  }

  # Rescue-fraction counts file. Column names and the pre/post
  # convention are a fixed contract: n_tdtom_pos / n_tdtom_neg
  # are PRE-rescue and must sum to n_cells_total, which is why they are taken
  # from tdtom_label_original rather than from the updated labels.
  # n_rescued_suspect: the part of the increment that came from a weaker tier.
  # The plot script draws it as its own band, so a reader can see how much of
  # the gain rests on B/C/D calls rather than on the confident ones.
  suspect_units <- sel$sub_id_full[sel$tier %in% RESCUE_SUSPECT]
  n_suspect <- sum(seu@meta.data[flip, "tdtom_rescue_source"] %in% suspect_units)

  lab0 <- seu$tdtom_label_original
  rescue_stats <- data.frame(
    sample              = case_id,
    n_cells_total       = sum(!is.na(lab0)),
    n_tdtom_pos         = sum(lab0 == "TdTom_pos", na.rm = TRUE),
    n_tdtom_neg         = sum(lab0 == "TdTom_neg", na.rm = TRUE),
    n_tdtom_pos_rescued = length(flip),
    n_rescued_suspect   = n_suspect,
    rescue_tiers        = paste(RESCUE_TIERS, collapse = "|"),
    suspect_tiers       = paste(RESCUE_SUSPECT, collapse = "|"),
    n_subclusters_used  = nrow(sel),
    stringsAsFactors    = FALSE)
  stopifnot(rescue_stats$n_tdtom_pos + rescue_stats$n_tdtom_neg ==
              rescue_stats$n_cells_total,
            rescue_stats$n_tdtom_pos_rescued <= rescue_stats$n_tdtom_neg,
            rescue_stats$n_rescued_suspect <= rescue_stats$n_tdtom_pos_rescued)
  if (n_suspect)
    message("  Of which from weaker tiers (", paste(RESCUE_SUSPECT, collapse = ", "),
            "): ", n_suspect, " — reported as 'suspect' in the counts file")
  write.csv(rescue_stats,
            file.path(OUT_DIR, paste0(case_id, "_subcluster_rescue.csv")),
            row.names = FALSE)
  write.csv(data.frame(cell_id = flip,
                       sub_id_full = seu@meta.data[flip, "tdtom_rescue_source"],
                       tier = sel$tier[match(seu@meta.data[flip, "tdtom_rescue_source"],
                                             sel$sub_id_full)],
                       stringsAsFactors = FALSE),
            file.path(OUT_DIR, paste0(case_id, "_rescued_cells.csv")),
            row.names = FALSE)

  if (length(flip) == 0)
    warning("Sub-cluster rescue selected ", nrow(sel), " sub-cluster(s) but every ",
            "cell in them was already TdTom_pos - no labels changed.", call. = FALSE)

} else {
  message("\n=== Step 2b: Sub-cluster rescue — SKIPPED (--subcluster_rescue FALSE) ===")
}

# =============================================================================
# 3. Convert signatures human -> mouse; intersect EC sigs with features
# =============================================================================

message("\n=== Step 3: Converting/intersecting signatures ===")

neftel_sig_mouse <- lapply(neftel_sig_human, function(genes) {
  converted <- homologene(genes, inTax = 9606, outTax = 10090)
  intersect(converted$`10090`, rownames(seu))
})
names(neftel_sig_mouse) <- c("MES2_like","MES1_like","AC_like",
                             "OPC_like","NPC1_like","NPC2_like")

prolif_sig_mouse <- lapply(prolif_sig_human, function(genes) {
  converted <- homologene(genes, inTax = 9606, outTax = 10090)
  intersect(converted$`10090`, rownames(seu))
})
names(prolif_sig_mouse) <- c("G1S", "G2M")

celltype_sig_mouse <- lapply(celltype_sig_mouse, function(genes)
  intersect(genes, rownames(seu)))

if (USE_EC) {
  ec_sig_mouse <- lapply(ec_sig_mouse_raw, function(genes)
    intersect(genes, rownames(seu)))
  ec_names <- names(ec_sig_mouse)

  message("EC genes present in object:")
  print(sapply(ec_sig_mouse, length))

  # ---- v2: report the shared core explicitly ------------------------------
  # Winner-takes-all across signatures that share a core is decided by the core,
  # not by the genes that distinguish the sub-states. Quantify the overlap so
  # the log makes it obvious how much is shared.
  if (length(ec_sig_mouse) >= 2) {
    ec_core <- Reduce(intersect, ec_sig_mouse)
    message("  Shared EC core genes (in ALL sub-states): ", length(ec_core),
            if (length(ec_core)) paste0(" [", paste(head(ec_core, 12), collapse = ", "),
                                        if (length(ec_core) > 12) ", ..." else "", "]") else "")
    ec_sig_unique <- lapply(ec_names, function(s)
      setdiff(ec_sig_mouse[[s]], unlist(ec_sig_mouse[setdiff(ec_names, s)])))
    names(ec_sig_unique) <- ec_names
    message("  Unique genes per EC sub-state:")
    print(sapply(ec_sig_unique, length))
  } else {
    ec_sig_unique <- ec_sig_mouse
  }

  ec_too_sparse <- sapply(ec_sig_mouse, length) < 5
  if (any(ec_too_sparse))
    warning("EC signatures with <5 genes in object: ",
            paste(ec_names[ec_too_sparse], collapse = ", "),
            " — permutation results may be unreliable")
}

# --- v3: myeloid STATE signatures -------------------------------------------
# Three transformations, in this order, each logged:
#   (1) species  — Miller sets are human symbols, converted via homologene;
#                  Qian sets are already mouse.
#   (2) intersect with rownames(seu).
#   (3) identity strip — remove genes that are also in the myeloid IDENTITY
#       signatures. Those genes are what selects the compartment in Step 6, so
#       leaving them in the STATE sets means the state test partly re-tests its
#       own selection criterion. Disable with --myeloid_strip_identity FALSE to
#       score the published signatures verbatim.
myeloid_sig_mouse <- list()
myeloid_audit     <- NULL

if (USE_MYELOID) {
  identity_genes <- unique(unlist(celltype_sig_mouse[MYELOID_CLASSES]))
  message("\nMyeloid STATE signatures (identity genes available for stripping: ",
          length(identity_genes), ")")

  for (setnm in names(myeloid_sig_raw)) {
    raw <- myeloid_sig_raw[[setnm]]
    conv <- lapply(raw, function(genes) {
      if (identical(MYELOID_SPECIES[[setnm]], "human")) {
        cv <- homologene(genes, inTax = 9606, outTax = 10090)
        intersect(cv$`10090`, rownames(seu))
      } else {
        intersect(genes, rownames(seu))
      }
    })
    stripped <- lapply(conv, function(g)
      if (MYELOID_STRIP) setdiff(g, identity_genes) else g)

    myeloid_sig_mouse[[setnm]] <- stripped

    myeloid_audit <- rbind(myeloid_audit, data.frame(
      set        = setnm,
      species    = MYELOID_SPECIES[[setnm]],
      state      = names(raw),
      n_published = sapply(raw, length),
      n_in_object = sapply(conv, length),
      n_identity_removed = sapply(names(conv), function(k)
        length(intersect(conv[[k]], identity_genes))),
      identity_removed   = sapply(names(conv), function(k)
        paste(intersect(conv[[k]], identity_genes), collapse = ";")),
      n_tested   = sapply(stripped, length),
      row.names  = NULL, stringsAsFactors = FALSE))
  }

  message("Myeloid state genes: published -> in object -> tested (after identity strip)")
  print(myeloid_audit[, c("set", "state", "n_published", "n_in_object",
                          "n_identity_removed", "n_tested")], row.names = FALSE)

  drop_small <- myeloid_audit$n_tested < MYELOID_MIN_GENES
  if (any(drop_small))
    warning("Myeloid states with <", MYELOID_MIN_GENES, " testable genes (will be skipped): ",
            paste(paste0(myeloid_audit$set[drop_small], ":", myeloid_audit$state[drop_small]),
                  collapse = ", "))

  # Statistical power is a strong function of set size and is NOT comparable
  # between the two sets. Simulation on this pipeline (500-600 cells, 4-6x
  # effect): ~80 genes -> single-state sensitivity ~0.78, hybrid recall ~0.54;
  # ~25 genes -> ~0.50 and ~0.22, falling to 0.008 hybrid recall at 4x. False
  # positives never exceeded 0.53% against a nominal 5% at any size.
  for (setnm in names(myeloid_sig_mouse)) {
    med <- median(myeloid_audit$n_tested[myeloid_audit$set == setnm])
    if (med < 40)
      message("  NOTE [", setnm, "]: median ", med, " genes/state. Hybrid recall is ",
              "roughly the square of single-state power at this size — report hybrid ",
              "counts from this set as a FLOOR, not an estimate, and do not compare ",
              "hybrid frequencies between sets.")
  }
  write.csv(myeloid_audit,
            file.path(OUT_DIR, paste0(case_id, "_myeloid_state_gene_audit.csv")),
            row.names = FALSE)
}

message("Neftel genes after conversion:")
print(sapply(neftel_sig_mouse, length))
message("Proliferation genes after conversion:")
print(sapply(prolif_sig_mouse, length))
message("Cell-type genes present:")
print(sapply(celltype_sig_mouse, length))

# =============================================================================
# 4. Score signatures (AddModuleScore — kept for plotting / continuity)
#    Neftel        -> TdTom+
#    CellType      -> TdTom-
#    Proliferation -> ALL
#    EC            -> ALL
# -----------------------------------------------------------------------------
# NOTE: in v2 these AddModuleScore values are NOT what decides the calls. The
# calls come from sig_enrich_perm() z-scores in Step 5. AddModuleScore is
# retained because the score columns are used in plots and downstream tables,
# and because keeping both lets you compare them.
# =============================================================================

message("\n=== Step 4: Scoring signatures (AddModuleScore, descriptive) ===")
neftel_names <- names(neftel_sig_mouse)
ct_names     <- names(celltype_sig_mouse)
prolif_names <- names(prolif_sig_mouse)

existing_sig <- grep("^Sig_|^Prolif_|^EC_", colnames(seu@meta.data), value = TRUE)
if (length(existing_sig) > 0) {
  message("  Dropping pre-existing score columns: ", paste(existing_sig, collapse = ", "))
  seu@meta.data[, existing_sig] <- NULL
}

message("  Neftel on TdTom+ cells ...")
seu_t <- subset(seu, cells = tumour_cells)
seu_t@meta.data[, grep("^Sig_", colnames(seu_t@meta.data), value = TRUE)] <- NULL
seu_t <- AddModuleScore(seu_t, features = neftel_sig_mouse, name = "Sig_", seed = SEED)
seu_t@meta.data <- seu_t@meta.data %>%
  rename_with(~ paste0("Sig_", neftel_names),
              all_of(paste0("Sig_", seq_along(neftel_names))))

message("  Cell types on TdTom- cells ...")
seu_nt <- subset(seu, cells = nontumour_cells)
seu_nt@meta.data[, grep("^Sig_", colnames(seu_nt@meta.data), value = TRUE)] <- NULL
seu_nt <- AddModuleScore(seu_nt, features = celltype_sig_mouse, name = "Sig_", seed = SEED)
seu_nt@meta.data <- seu_nt@meta.data %>%
  rename_with(~ paste0("Sig_", ct_names),
              all_of(paste0("Sig_", seq_along(ct_names))))

message("  Proliferation on ALL cells ...")
seu <- AddModuleScore(seu, features = prolif_sig_mouse, name = "Prolif_", seed = SEED)
seu@meta.data <- seu@meta.data %>%
  rename_with(~ paste0("Prolif_", prolif_names),
              all_of(paste0("Prolif_", seq_along(prolif_names))))

if (USE_EC) {
  message("  EC signatures on ALL cells (scoring) ...")
  seu <- AddModuleScore(seu, features = ec_sig_mouse, name = "EC_", seed = SEED)
  seu@meta.data <- seu@meta.data %>%
    rename_with(~ paste0("EC_", ec_names),
                all_of(paste0("EC_", seq_along(ec_names))))
}

all_class_names <- c(neftel_names, ct_names)
score_cols      <- paste0("Sig_", all_class_names)
seu@meta.data[, score_cols] <- NA_real_
seu@meta.data[tumour_cells,    paste0("Sig_", neftel_names)] <-
  seu_t@meta.data[tumour_cells,    paste0("Sig_", neftel_names)]
seu@meta.data[nontumour_cells, paste0("Sig_", ct_names)] <-
  seu_nt@meta.data[nontumour_cells, paste0("Sig_", ct_names)]

# =============================================================================
# 5. Permutation tests — expression-matched null
# =============================================================================

message("\n=== Step 5: Permutation tests (", N_PERM, " permutations) ===")

# -----------------------------------------------------------------------------
# sig_enrich_perm() — per-cell signature enrichment with an expression-matched null
# Replaces modulescore_perm_pval(). See the v2 change list at the top of this
# file for the empirical justification of each change.
# -----------------------------------------------------------------------------
# RETURNS a data.frame, one row per cell:
#   score    signature mean log-expression in that cell
#   z        (score - null mean) / null sd     <- USE THIS TO RANK SIGNATURES
#   pval     empirical permutation p, floors at 1/(n_perm+1)
#   padj     BH on pval — DO NOT GATE ON THIS (see note below)
#   pval_z   normal-tail p on z, continuous, not floor-limited
#   padj_z   BH on pval_z                      <- USE THIS TO GATE CALLS
#   n_genes  signature genes actually present
#
# WHY NOT padj: BH at n cells needs the smallest p below ALPHA/n, but the
# empirical p cannot go below 1/(n_perm+1). At n_perm = 1000 and 1200 cells the
# floor is 1e-3 while BH needs 4.2e-5 — NO cell can be called at any effect
# size. You would need n_perm > n_cells/ALPHA. pval_z is continuous, so BH
# behaves sensibly on it; it assumes the null score is roughly normal across
# cells, which Step 5c checks on your data.
#
# expr_ref: the cells used to build the matching strata. DEFAULT (= expr) is
# self-reference and is the safe choice. Do NOT pass the opposite lineage
# population — a population-restricted gene lands in a bottom stratum of the
# other population, draws near-zero matched controls, and the signature then
# beats its own null in every cell (simulated: 100% false positives).
# -----------------------------------------------------------------------------
sig_enrich_perm <- function(expr, sig_genes,
                            expr_ref  = NULL,
                            n_perm    = 1000,
                            nbin      = 25,
                            n_det_bin = 5,
                            seed      = 1,
                            sig_name  = NULL,
                            min_genes = 5,
                            borrow_strata = getOption("neftel.borrow_strata", TRUE),
                            min_pool      = 5L,
                            skew_p        = getOption("neftel.skew_p", FALSE)) {
  stopifnot(is.matrix(expr) || inherits(expr, "Matrix"))
  if (is.null(expr_ref)) expr_ref <- expr
  stopifnot(all(rownames(expr) %in% rownames(expr_ref)))

  sig_genes <- intersect(unique(sig_genes), rownames(expr))
  n <- length(sig_genes)
  if (n < 3) stop("Fewer than 3 signature genes present in expr - cannot test.")
  if (n < min_genes)
    warning("Signature ", if (is.null(sig_name)) "" else paste0("'", sig_name, "' "),
            "has only ", n, " genes present; enrichment p-values are coarse and ",
            "z is unstable. Consider merging or dropping this signature.")

  # per-signature RNG stream: independent nulls across signatures
  seed_use <- if (is.null(sig_name)) seed else
    seed + sum(utf8ToInt(sig_name)) %% 100000L
  set.seed(seed_use)

  # --- 2D matching strata: mean expression x detection rate -------------------
  gmean  <- Matrix::rowMeans(expr_ref)
  gdet   <- Matrix::rowMeans(expr_ref > 0)
  gmean  <- gmean[rownames(expr)]
  gdet   <- gdet[rownames(expr)]
  b_mean <- cut(rank(gmean, ties.method = "random"), breaks = nbin,      labels = FALSE)
  b_det  <- cut(rank(gdet,  ties.method = "random"), breaks = n_det_bin, labels = FALSE)
  strat  <- paste(b_mean, b_det, sep = "_")
  names(strat) <- rownames(expr)

  pool    <- split(rownames(expr), strat)
  # exclude signature genes from their own pools so the null cannot resample them
  pool_ex <- lapply(pool, setdiff, sig_genes)

  # --- v7 FIX -----------------------------------------------------------------
  # v6 did:  empty <- lengths(pool_ex) == 0; pool_ex[empty] <- pool[empty]
  # i.e. when a stratum contained ONLY signature genes the null was allowed to
  # resample the signature genes themselves. mu then rises toward obs, z shrinks,
  # and the best-detected signature is penalised. On a targeted panel this fires
  # often (many near-empty strata), which is why it matters here and did not on
  # whole-transcriptome data. Fix: widen the mean-expression window until we find
  # >= min_pool clean controls; never admit a signature gene to the null.
  n_borrowed <- 0L
  if (isTRUE(borrow_strata)) {
    ordg <- names(sort(gmean))
    pos  <- setNames(seq_along(ordg), ordg)
    thin <- names(pool_ex)[vapply(pool_ex, length, integer(1)) < min_pool]
    for (k in thin) {
      centre <- round(stats::median(pos[pool[[k]]]))
      w <- max(min_pool, 10L)
      repeat {
        cand <- setdiff(ordg[max(1L, centre - w):min(length(ordg), centre + w)], sig_genes)
        if (length(cand) >= min_pool || w >= length(ordg)) break
        w <- w * 2L
      }
      if (length(cand)) { pool_ex[[k]] <- cand; n_borrowed <- n_borrowed + 1L }
      else              { pool_ex[[k]] <- setdiff(rownames(expr), sig_genes) }
    }
  } else {
    empty <- vapply(pool_ex, length, integer(1)) == 0L
    pool_ex[empty] <- pool[empty]          # v6 behaviour, retained for comparison
  }

  # --- statistic: plain signature mean ---------------------------------------
  # (v1 subtracted a fresh random control set on EVERY permutation, which
  #  injects control-set noise into the null and makes the observed and null
  #  statistics non-identical.)
  stat <- function(gs) Matrix::colMeans(expr[gs, , drop = FALSE])
  obs  <- stat(sig_genes)

  sig_strat <- strat[sig_genes]
  null_mat  <- vapply(seq_len(n_perm), function(p) {
    gs <- vapply(sig_strat, function(k) {
      pk <- pool_ex[[k]]
      pk[sample.int(length(pk), 1L)]
    }, character(1))
    stat(gs)
  }, numeric(ncol(expr)))

  mu  <- rowMeans(null_mat)
  sd0 <- apply(null_mat, 1, stats::sd)
  sd0[!is.finite(sd0) | sd0 == 0] <- NA_real_

  pval <- (1 + rowSums(null_mat >= obs)) / (n_perm + 1)
  z    <- as.numeric((obs - mu) / sd0)

  # A degenerate null (sd0 == 0: every control set scored identically, which
  # happens when a cell detects almost nothing) leaves z = NA but pval FINITE,
  # often at the 1/(n_perm+1) floor. That combination is a trap: a gate on the
  # empirical p ADMITS such a cell while every downstream ranking step uses z,
  # so which.max() gets an all-NA vector and returns integer(0) -- the
  # "values must be length 1, but FUN(X[[i]]) result is length 0" crash in
  # Step 6. v6 never hit it because pval_z is NA wherever z is, so its gate
  # excluded these cells implicitly.
  #
  # The gate and the ranking statistic must agree on which cells are evaluable.
  # A cell with no usable null is not evidence of enrichment -- it is absence of
  # information -- so NA the p as well and count it.
  .degen <- !is.finite(z)
  if (any(.degen)) {
    pval[.degen] <- NA_real_
    message(sprintf("    %d/%d cells (%.2f%%) have a degenerate null (sd0=0) -- p set NA",
                    sum(.degen), length(z), 100 * mean(.degen)))
  }

  # --- v7: optional Cornish-Fisher skew correction on the tail --------------
  # The null of a signature mean is right-skewed on sparse counts (median
  # skewness 0.35 at 300 genes / 60 counts per cell). pnorm() is accurate at
  # p ~ 0.05 but ~3.5x anti-conservative at p ~ 0.001, which is where a
  # BH-across-cells gate actually operates. Correcting makes the gate STRICTER.
  if (isTRUE(skew_p)) {
    g1 <- apply(null_mat, 1, function(x) {
      s <- stats::sd(x)
      if (!is.finite(s) || s == 0) NA_real_ else mean((x - mean(x))^3) / s^3
    })
    zc <- ifelse(is.na(g1) | abs(g1) < 1e-6, z,
                 (sqrt(pmax(1 + 2 * g1 * z, 1e-8)) - 1) / g1)
    pz <- stats::pnorm(zc, lower.tail = FALSE)
  } else {
    pz <- stats::pnorm(z, lower.tail = FALSE)
  }

  out <- data.frame(cell    = colnames(expr),
             score   = as.numeric(obs),
             z       = z,
             pval    = as.numeric(pval),
             padj    = stats::p.adjust(as.numeric(pval), method = "BH"),
             pval_z  = pz,
             padj_z  = stats::p.adjust(pz, method = "BH"),
             n_genes = n,
             row.names = NULL, stringsAsFactors = FALSE)
  attr(out, "n_strata_borrowed") <- n_borrowed
  attr(out, "n_strata_occupied") <- length(pool)
  attr(out, "median_pool")       <- stats::median(vapply(pool_ex, length, integer(1)))
  out
}

# -----------------------------------------------------------------------------
# auto_strata() — v7. Panel-aware matching strata.
# -----------------------------------------------------------------------------
# nbin=25 x n_det_bin=5 = 125 strata assumes a whole transcriptome. On a 300-gene
# panel only ~31 are occupied (mean and detection rate are strongly correlated),
# median 3 control genes each, ~40% holding <=2 — so sd0 is estimated off almost
# nothing and every z is noise-dominated. Target >=25 controls per stratum.
# ONLY EVER REDUCES, and only when the gene count cannot support the strata you
# asked for. It is passed the current CLI values and never exceeds them, so:
#
#   whole transcriptome (~18000 genes)  -> returns 25 x 5 unchanged: A NO-OP.
#   targeted panel (~300 genes)         -> coarsens to ~10 x 3.
#
# This matters because the two are opposite problems. With 18000 genes, 125
# strata hold ~144 control genes each and the matched null is well estimated —
# there is nothing to fix and coarsening would only make the matching worse.
# Validated on a simulated 300-plex panel, where 25x5 -> 10x3 raised accuracy
# among assigned cells from 71 to 79 while cutting the no-signal false-positive
# rate; that evidence does NOT transfer to whole-transcriptome input, which is
# why this function declines to act there.
auto_strata <- function(n_genes, target_per_stratum = 10L,
                        nbin_cur = 25L, n_det_bin_cur = 5L) {
  n_det_bin <- min(as.integer(n_det_bin_cur), 3L)
  cap_total <- floor(n_genes / target_per_stratum)
  if (cap_total >= as.integer(nbin_cur) * as.integer(n_det_bin_cur))
    return(list(nbin = as.integer(nbin_cur), n_det_bin = as.integer(n_det_bin_cur)))
  nbin <- min(as.integer(nbin_cur), max(4L, as.integer(floor(cap_total / n_det_bin))))
  list(nbin = as.integer(nbin), n_det_bin = as.integer(n_det_bin))
}

# -----------------------------------------------------------------------------
# gate_stat() — v7. Choose the MULTIPLICITY UNIT.
# -----------------------------------------------------------------------------
# Takes the per-signature pval_z columns for one population and returns the
# gate statistic under the requested unit.
#   bh_global  BH across CELLS within each signature (v6 padj_z).
#   bh_cell    BH across SIGNATURES within each cell.
#   raw        uncorrected pval_z.
# bh_global treats every cell as an independent discovery, which is the right
# unit only when a single positive cell is itself the claim. It also gives each
# signature its own adaptive threshold, so per-signature yield diverges and
# recovered COMPOSITION is distorted (see header).
gate_stat <- function(P, unit = c("bh_global", "bh_cell", "raw")) {
  unit <- match.arg(unit)
  stopifnot(is.matrix(P))
  Q <- switch(unit,
    bh_global = apply(P, 2, function(p) {
      out <- rep(NA_real_, length(p)); ok <- !is.na(p)
      out[ok] <- stats::p.adjust(p[ok], method = "BH"); out }),
    bh_cell   = t(apply(P, 1, function(p) {
      out <- rep(NA_real_, length(p)); ok <- !is.na(p)
      out[ok] <- stats::p.adjust(p[ok], method = "BH"); out })),
    raw       = P)
  if (!is.matrix(Q)) Q <- matrix(Q, nrow = nrow(P))
  dimnames(Q) <- dimnames(P)
  Q
}

# -----------------------------------------------------------------------------
# neighbour_pool() — v7, OPT-IN. Average expression over spatial neighbours.
# -----------------------------------------------------------------------------
# Raises per-cell evidence when the biology is spatially autocorrelated. It does
# NOT create information — it borrows from neighbours. On patchy simulated tissue
# it raised both coverage and accuracy; where the truth was NOT spatially
# structured, accuracy at MATCHED coverage was 3-6 points WORSE than no pooling.
# Use real centroids, k = 6-10, check Moran s I first, and never present a pooled
# label as evidence that neighbours share a state — the smoothing assumes it.
neighbour_pool <- function(expr, coords, k = 8L, max_dist = Inf) {
  if (!requireNamespace("RANN", quietly = TRUE))
    stop("neighbour_pool() needs the RANN package.")
  stopifnot(nrow(coords) == ncol(expr))
  keep <- stats::complete.cases(coords)
  if (!all(keep))
    warning("neighbour_pool(): ", sum(!keep), " cells have missing coordinates; ",
            "they are left unpooled.")
  nn  <- RANN::nn2(coords[keep, , drop = FALSE],
                   k = min(k + 1L, sum(keep)))
  idx <- nn$nn.idx; dst <- nn$nn.dists
  which_keep <- which(keep)
  out <- as.matrix(expr)
  for (j in seq_along(which_keep)) {
    use <- idx[j, dst[j, ] <= max_dist]
    if (!length(use)) next
    out[, which_keep[j]] <- Matrix::rowMeans(expr[, which_keep[use], drop = FALSE])
  }
  attr(out, "pool_k") <- k
  out
}

# =============================================================================
# assign_states_with_hybrids()  —  Neftel-faithful per-cell state + hybrid call
# -----------------------------------------------------------------------------
# Design: the enrichment tests decide WHICH programs are real; POSITIVE
# MAGNITUDE thresholds (Neftel's actual rule) decide whether a second real
# program is genuinely co-dominant. We deliberately do NOT use "difference
# within a null" to call hybrids — failing to reject equality is not evidence of
# equality, and it manufactures hybrids in the noisiest cells.
#
# A cell is HYBRID (A|B) iff ALL hold (A = top state, B = 2nd):
#   (0) collapse sub-state families first (MES1/MES2 -> MES, NPC1/NPC2 -> NPC)
#   (a) A passes enrichment (q_A < alpha)                       -> A is real
#   (b) B passes enrichment (q_B < alpha)                       -> B is real
#   (c) score_B > pct_floor-quantile of B among B-top cells     -> Neftel crit 2
#   (d) score_B - score_C >= delta  (C = 3rd)                   -> Neftel crit 3
# Otherwise single-state A (or not_assigned if A itself fails (a)).
#
# v2: score_mat is a matrix of Z-SCORES and pval_mat a matrix of BH-adjusted
# z-tail q-values. Ranking on z rather than on the empirical p removes the
# permutation-floor tie problem (with n_perm = 1000 the empirical p floors at
# 1/1001; in simulation 55% of true-positive cells sat exactly at the floor for
# a 40-gene signature, so which.min(p) broke ties by COLUMN ORDER, silently
# favouring whichever signature was listed first and, separately, whichever had
# more genes). z is continuous and comparable across signature sizes, and it
# makes delta dimensionless.
#
# Neftel's own criterion 1 ("score > 1") is scale-specific to their SmartSeq2
# relative-expression units and does not transfer; it is REPLACED by the
# adaptive percentile floor (c).
#
# Continuous outputs for plotting (no extra permutation needed):
#   gap_12 = score1 - score2   (small + both significant = strong hybrid)
#   gap_23 = score2 - score3   (the criterion-(d) quantity)
# =============================================================================

# -----------------------------------------------------------------------------
# .null_gap23_delta()  —  data-driven delta for criterion (d), no hard-coded value
# -----------------------------------------------------------------------------
# Repurposes Neftel's expected-hybrid shuffle: independently permute each
# collapsed state's scores across cells (destroys co-expression structure,
# preserves each state's marginal distribution), recompute gap_23 = score2 -
# score3 per shuffled cell, pool over shuffles, and return the delta_null_q
# quantile. A real cell must beat this to satisfy criterion (d).
# Interpretation: under independence, P(gap_23 >= delta) = 1 - delta_null_q per
# cell, i.e. delta_null_q is criterion (d)'s per-cell specificity.
.null_gap23_delta <- function(coll_score, delta_null_q = 0.95,
                              n_shuffle = 100, seed = 1, max_cells = 5000) {
  n <- nrow(coll_score); k <- ncol(coll_score)
  if (k < 3) return(0)                       # no 3rd state -> criterion (d) inert
  set.seed(seed)
  idx <- if (n > max_cells) sample.int(n, max_cells) else seq_len(n)
  sm  <- coll_score[idx, , drop = FALSE]
  gap23_of <- function(x) { xs <- sort(x, decreasing = TRUE); xs[2] - xs[3] }
  vals <- numeric(length(idx) * n_shuffle)
  for (s in seq_len(n_shuffle)) {
    shuf <- apply(sm, 2, sample)             # independent per-column permutation
    g    <- apply(shuf, 1, gap23_of)
    vals[((s - 1) * length(idx) + 1):(s * length(idx))] <- g
  }
  as.numeric(quantile(vals, probs = delta_null_q, names = FALSE, na.rm = TRUE))
}

assign_states_with_hybrids <- function(score_mat, pval_mat,
                                       collapse = list(
                                         MES_like = c("MES1_like","MES2_like"),
                                         NPC_like = c("NPC1_like","NPC2_like"),
                                         AC_like  = "AC_like",
                                         OPC_like = "OPC_like"),
                                       alpha        = 0.05,
                                       delta        = "auto",
                                       delta_null_q = 0.95,
                                       n_shuffle    = 100,
                                       pct_floor    = 0.10,
                                       min_top_n    = 10,
                                       seed         = 1) {
  stopifnot(identical(rownames(score_mat), rownames(pval_mat)))
  cells <- rownames(score_mat)
  fams  <- names(collapse)
  nf    <- length(fams)
  if (nf < 2) stop("Need >= 2 collapsed families to define hybrids.")

  # --- (0) collapse sub-state families: score = max member, pval = min member
  coll_score <- vapply(fams, function(f) {
    m <- intersect(collapse[[f]], colnames(score_mat))
    if (length(m) == 1) score_mat[, m] else apply(score_mat[, m, drop = FALSE], 1, max)
  }, numeric(length(cells)))
  coll_pval <- vapply(fams, function(f) {
    m <- intersect(collapse[[f]], colnames(pval_mat))
    if (length(m) == 1) pval_mat[, m] else apply(pval_mat[, m, drop = FALSE], 1, min)
  }, numeric(length(cells)))
  dimnames(coll_score) <- list(cells, fams)
  dimnames(coll_pval)  <- list(cells, fams)

  if (is.character(delta) && identical(delta, "auto")) {
    delta <- .null_gap23_delta(coll_score, delta_null_q = delta_null_q,
                               n_shuffle = n_shuffle, seed = seed)
    attr_delta_auto <- TRUE
  } else {
    stopifnot(is.numeric(delta), length(delta) == 1)
    attr_delta_auto <- FALSE
  }

  sub_winner <- vapply(fams, function(f) {
    m <- intersect(collapse[[f]], colnames(score_mat))
    if (length(m) == 1) rep(m, length(cells))
    else m[max.col(score_mat[, m, drop = FALSE], ties.method = "first")]
  }, character(length(cells)))
  dimnames(sub_winner) <- list(cells, fams)

  ord <- t(apply(coll_score, 1, order, decreasing = TRUE))
  ri  <- seq_len(nrow(coll_score))
  i1 <- ord[, 1]; i2 <- ord[, 2]
  i3 <- if (nf >= 3) ord[, 3] else rep(NA_integer_, length(ri))
  s1 <- coll_score[cbind(ri, i1)]
  s2 <- coll_score[cbind(ri, i2)]
  s3 <- if (nf >= 3) coll_score[cbind(ri, i3)] else rep(-Inf, length(ri))
  p1 <- coll_pval[cbind(ri, i1)]
  p2 <- coll_pval[cbind(ri, i2)]
  f1 <- fams[i1]; f2 <- fams[i2]

  # --- (c) per-family adaptive floor
  top_fam <- fams[i1]
  b_floor <- vapply(fams, function(f) {
    s <- coll_score[top_fam == f, f]
    if (length(s) < min_top_n) return(-Inf)
    as.numeric(quantile(s, probs = pct_floor, names = FALSE, na.rm = TRUE))
  }, numeric(1))
  names(b_floor) <- fams
  floor2 <- b_floor[f2]

  a_ok      <- !is.na(p1) & p1 < alpha
  primary   <- ifelse(a_ok, f1, "not_assigned")
  hybrid_ok <- a_ok & !is.na(p2) & p2 < alpha & (s2 > floor2) & ((s2 - s3) >= delta)
  hybrid    <- ifelse(hybrid_ok, paste0(f1, "|", f2), NA_character_)

  primary_sub <- ifelse(a_ok, sub_winner[cbind(ri, i1)], NA_character_)
  hybrid_sub  <- ifelse(hybrid_ok,
                        paste0(sub_winner[cbind(ri, i1)], "|", sub_winner[cbind(ri, i2)]),
                        NA_character_)

  out <- data.frame(
    cell         = cells,
    cell_class   = primary,
    primary_sub  = primary_sub,
    hybrid_class = hybrid,
    hybrid_sub   = hybrid_sub,
    top1_state = ifelse(a_ok, f1, NA_character_),
    top2_state = f2,
    top1_score = s1, top2_score = s2, top3_score = s3,
    gap_12     = s1 - s2,
    gap_23     = s2 - s3,
    p_top1     = p1, p_top2 = p2,
    is_hybrid  = hybrid_ok,
    row.names  = NULL, stringsAsFactors = FALSE
  )
  attr(out, "delta")        <- delta
  attr(out, "delta_auto")   <- attr_delta_auto
  attr(out, "delta_null_q") <- if (attr_delta_auto) delta_null_q else NA_real_
  out
}

# =============================================================================
# 5b. Run the enrichment tests
# -----------------------------------------------------------------------------
# expr_ref is left at its default (self-reference) for BOTH loops:
#   Neftel states   scored on TdTom+  with strata from TdTom+
#   Cell types      scored on TdTom-  with strata from TdTom-
# This is the conservative and correct pairing. Do NOT cross them.
# =============================================================================

expr_t   <- GetAssayData(seu_t,  assay = "SCT", layer = "data")
expr_nt  <- GetAssayData(seu_nt, assay = "SCT", layer = "data")
expr_all <- GetAssayData(seu,    assay = "SCT", layer = "data")

# --- v7: panel-aware matching strata ----------------------------------------
# NBIN/N_DET_BIN from the CLI are whole-transcriptome defaults. On a targeted
# panel they leave most strata with a handful of genes, so sd0 in the null is
# estimated off almost nothing. auto_strata() derives them from the genes
# actually present. Reported so the choice is visible in the log.
if (AUTO_STRATA) {
  .as <- auto_strata(nrow(expr_all), target_per_stratum = STRATA_TARGET,
                     nbin_cur = NBIN, n_det_bin_cur = N_DET_BIN)
  if (.as$nbin != NBIN || .as$n_det_bin != N_DET_BIN)
    message("  v7 auto_strata: ", nrow(expr_all), " genes -> nbin=", .as$nbin,
            " x n_det_bin=", .as$n_det_bin,
            " (CLI defaults were ", NBIN, " x ", N_DET_BIN, ")")
  NBIN <- .as$nbin; N_DET_BIN <- .as$n_det_bin
} else {
  message("  v7 auto_strata disabled; using nbin=", NBIN, " x n_det_bin=", N_DET_BIN)
}

# --- v7: optional neighbour pooling (OPT-IN, --pool_k) ----------------------
# Pooling borrows evidence from neighbours. It raises coverage, and raises
# accuracy only when the biology is genuinely patchy; at matched coverage it was
# WORSE than no pooling on unstructured simulated tissue. Never present a pooled
# label as evidence that neighbouring cells share a state.
if (POOL_K > 0L) {
  miss <- setdiff(POOL_COORDS, colnames(seu@meta.data))
  if (length(miss))
    stop("--pool_coords column(s) not in meta.data: ", paste(miss, collapse = ", "))
  .xy_all <- as.matrix(seu@meta.data[colnames(expr_all), POOL_COORDS, drop = FALSE])
  message("  v7 neighbour pooling: k=", POOL_K, " over ", paste(POOL_COORDS, collapse = ","))
  expr_t   <- neighbour_pool(expr_t,   .xy_all[colnames(expr_t),   , drop = FALSE], k = POOL_K)
  expr_nt  <- neighbour_pool(expr_nt,  .xy_all[colnames(expr_nt),  , drop = FALSE], k = POOL_K)
  expr_all <- neighbour_pool(expr_all, .xy_all,                                     k = POOL_K)
}

# -----------------------------------------------------------------------------
# canon_pair() — order-free hybrid label
# -----------------------------------------------------------------------------
# assign_states_with_hybrids() names hybrids by RANK: the same biology arrives as
# "A|B" in cells where A scores higher and "B|A" where B does, so one biological
# class is split across two labels and every table double-counts it.
#
# This also affects the NEFTEL calls in Step 6 — MES_like|NPC_like and
# NPC_like|MES_like are currently tallied as different classes. Both Step 6 and
# Step 6c now carry a *_canonical column alongside the rank-ordered one; the raw
# column is kept so nothing downstream that depends on rank order breaks.
canon_pair <- function(x) {
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    paste(sort(strsplit(s, "|", fixed = TRUE)[[1]]), collapse = "|")
  }, character(1), USE.NAMES = FALSE)
}

mk_mat <- function(cells, cols) matrix(NA_real_, nrow = length(cells), ncol = length(cols),
                                       dimnames = list(cells, cols))

z_mat    <- mk_mat(colnames(seu), all_class_names)   # z-scores (ranking)
q_mat    <- mk_mat(colnames(seu), all_class_names)   # the GATE statistic
pval_mat <- mk_mat(colnames(seu), all_class_names)   # raw empirical p (diagnostics only)
pz_mat   <- mk_mat(colnames(seu), all_class_names)   # v7: uncorrected z-tail p
message("  Gate statistic: ", GATE)

# --- Neftel states (scored on TdTom+, strata from TdTom+) ---
message("  -- Neftel states on TdTom+ (self-referenced null) --")
for (s in neftel_names) {
  message("  ", s, " (n=", length(neftel_sig_mouse[[s]]), ")")
  res <- sig_enrich_perm(expr_t, neftel_sig_mouse[[s]],
                         n_perm = N_PERM, nbin = NBIN, n_det_bin = N_DET_BIN,
                         seed = SEED, sig_name = s)
  m <- match(tumour_cells, res$cell)
  z_mat[tumour_cells, s]    <- res$z[m]
  q_mat[tumour_cells, s]    <- res[[GATE]][m]
  pval_mat[tumour_cells, s] <- res$pval[m]
  pz_mat[tumour_cells, s]   <- res$pval_z[m]
}

# --- Cell types (scored on TdTom-, strata from TdTom-) ---
message("  -- Cell types on TdTom- (self-referenced null) --")
for (s in ct_names) {
  message("  ", s, " (n=", length(celltype_sig_mouse[[s]]), ")")
  res <- sig_enrich_perm(expr_nt, celltype_sig_mouse[[s]],
                         n_perm = N_PERM, nbin = NBIN, n_det_bin = N_DET_BIN,
                         seed = SEED, sig_name = s)
  m <- match(nontumour_cells, res$cell)
  z_mat[nontumour_cells, s]    <- res$z[m]
  q_mat[nontumour_cells, s]    <- res[[GATE]][m]
  pval_mat[nontumour_cells, s] <- res$pval[m]
  pz_mat[nontumour_cells, s]   <- res$pval_z[m]
}

# --- Proliferation (scored on all cells, strata from all cells) ---
# Proliferation is legitimately a whole-tissue question, and G1S/G2M genes are
# not restricted to one lineage, so the all-cell reference is appropriate here.
z_prolif <- mk_mat(colnames(seu), prolif_names)
q_prolif <- mk_mat(colnames(seu), prolif_names)
pval_prolif <- mk_mat(colnames(seu), prolif_names)
pz_prolif   <- mk_mat(colnames(seu), prolif_names)   # v7

message("  -- Proliferation on ALL cells --")
for (s in prolif_names) {
  message("  ", s, " (n=", length(prolif_sig_mouse[[s]]), ")")
  res <- sig_enrich_perm(expr_all, prolif_sig_mouse[[s]],
                         n_perm = N_PERM, nbin = NBIN, n_det_bin = N_DET_BIN,
                         seed = SEED, sig_name = s)
  m <- match(colnames(seu), res$cell)
  z_prolif[, s]    <- res$z[m]
  q_prolif[, s]    <- res[[GATE]][m]
  pval_prolif[, s] <- res$pval[m]
  pz_prolif[, s]   <- res$pval_z[m]
}

if (USE_EC)
  message("  EC permutation deferred to post-classification (requires cell_class).")

# =============================================================================
# 5b. v7: apply the chosen MULTIPLICITY UNIT to the gate statistic
# -----------------------------------------------------------------------------
# q_mat currently holds whatever --gate named, computed inside sig_enrich_perm()
# over the cells of one population. --gate_unit re-derives it from the
# uncorrected z-tail p under the requested unit, per population, so the BH
# denominator is the set of decisions being made rather than the cell census.
#
# bh_global reproduces v6. bh_cell corrects across signatures within a cell.
# Read the header note on composition distortion before choosing.
if (!identical(GATE_UNIT, "bh_global") || !identical(GATE, "padj_z") ||
    !identical(P_SOURCE, "z")) {
  message("\n=== Step 5b: gate unit = ", GATE_UNIT, ", p source = ", P_SOURCE, " ===")
  # P_SOURCE 'emp' gates on the permutation p itself instead of the normal-tail
  # approximation to it. Measured on 7000/7003_Scan2: the normal tail is ~14x
  # anti-conservative at p=0.004 and worse further out, so under bh_global (which
  # needs p ~ 1e-6) the v6 calls rested entirely on that extrapolation.
  .P_of <- function(rows, cols)
    if (identical(P_SOURCE, "emp")) pval_mat[rows, cols, drop = FALSE]
    else                            pz_mat[rows, cols, drop = FALSE]
  .apply_unit <- function(rows, cols) {
    if (!length(rows) || !length(cols)) return(invisible(NULL))
    q_mat[rows, cols] <<- gate_stat(.P_of(rows, cols), GATE_UNIT)
  }
  .apply_unit(tumour_cells,    neftel_names)
  .apply_unit(nontumour_cells, ct_names)
  q_prolif[, ] <- gate_stat(if (identical(P_SOURCE, "emp")) pval_prolif else pz_prolif,
                            GATE_UNIT)

  for (nm in c("Neftel states (tumour)", "cell types (non-tumour)")) {
    rows <- if (grepl("tumour\\)", nm) && !grepl("non", nm)) tumour_cells else nontumour_cells
    cols <- if (grepl("Neftel", nm)) neftel_names else ct_names
    if (!length(rows)) next
    Q <- q_mat[rows, cols, drop = FALSE]
    message("  ", nm, ": ", round(100 * mean(rowSums(!is.na(Q) & Q < ALPHA) > 0), 1),
            "% of cells clear the gate")
  }
}

# =============================================================================
# 5c. Null-calibration diagnostics  (--calibration_check)
# -----------------------------------------------------------------------------
# Two checks, both on YOUR data. The simulation numbers in the header
# characterise the estimator, not this tissue.
#
#  (i) GENE-LABEL SHUFFLE. Permuting gene labels within each cell destroys
#      co-expression structure but preserves each cell's library size and each
#      gene's marginal. No signature should be enriched, so the fraction of
#      cells called should be ~ALPHA. Materially above ALPHA = the null is still
#      anti-conservative on this matrix; materially below = over-conservative.
#
# (ii) POPULATION-RESTRICTION CONTROL (--control_signature). Check (i) is BLIND
#      to the expr_ref failure mode: generic no-signal genes are well calibrated
#      under every reference choice (3.9 / 5.6 / 5.1% for self / all / cross in
#      simulation). The failure only appears for POPULATION-RESTRICTED gene
#      sets. So supply a signature you believe is tumour-restricted but NOT
#      state-specific and confirm it does not fire in nearly every TdTom+ cell.
#      If it does, the strata reference is wrong for your data.
# =============================================================================

if (DO_CALIB) {
  message("\n=== Step 5c: Null calibration on this dataset ===")
  set.seed(SEED)

  calib_rows <- list()

  shuffle_genes_within_cells <- function(m) {
    m2 <- as.matrix(m)
    for (j in seq_len(ncol(m2))) m2[, j] <- m2[sample.int(nrow(m2)), j]
    rownames(m2) <- rownames(m); colnames(m2) <- colnames(m)
    m2
  }

  run_calib <- function(expr_sub, sig_list, label, n_cells = CALIB_N_CELLS) {
    cs <- colnames(expr_sub)
    if (length(cs) > n_cells) cs <- sample(cs, n_cells)
    e  <- expr_sub[, cs, drop = FALSE]
    e  <- e[Matrix::rowSums(e) > 0, , drop = FALSE]
    es <- shuffle_genes_within_cells(e)
    for (s in names(sig_list)) {
      g <- intersect(sig_list[[s]], rownames(es))
      if (length(g) < 3) next
      r <- sig_enrich_perm(es, g, n_perm = min(N_PERM, 500),
                           nbin = NBIN, n_det_bin = N_DET_BIN,
                           seed = SEED, sig_name = s)
      calib_rows[[length(calib_rows) + 1]] <<- data.frame(
        check      = "gene_label_shuffle",
        population = label,
        signature  = s,
        n_genes    = length(g),
        n_cells    = ncol(es),
        pct_called_raw_p = round(100 * mean(r$pval   < ALPHA, na.rm = TRUE), 2),
        pct_called_gate  = round(100 * mean(r[[GATE]] < ALPHA, na.rm = TRUE), 2),
        median_z         = round(stats::median(r$z, na.rm = TRUE), 3),
        expected_pct     = 100 * ALPHA,
        stringsAsFactors = FALSE)
      message(sprintf("   [shuffle] %-22s %-20s raw %5.2f%%  %s %5.2f%%  (expect ~%.0f%% raw)",
                      label, s,
                      100 * mean(r$pval < ALPHA, na.rm = TRUE), GATE,
                      100 * mean(r[[GATE]] < ALPHA, na.rm = TRUE), 100 * ALPHA))
    }
  }

  run_calib(expr_t,  neftel_sig_mouse,   "TdTom_pos")
  run_calib(expr_nt, celltype_sig_mouse, "TdTom_neg")

  # ---- (ii) population-restriction control ---------------------------------
  if (!is.null(CONTROL_SIG) && file.exists(CONTROL_SIG)) {
    ctrl_genes <- unique(na.omit(read.csv(CONTROL_SIG, stringsAsFactors = FALSE)[[1]]))
    ctrl_genes <- intersect(ctrl_genes, rownames(expr_t))
    message("  Control signature genes present: ", length(ctrl_genes))
    if (length(ctrl_genes) >= 3) {
      rc_self <- sig_enrich_perm(expr_t, ctrl_genes, n_perm = N_PERM,
                                 nbin = NBIN, n_det_bin = N_DET_BIN,
                                 seed = SEED, sig_name = "control_self")
      rc_all  <- sig_enrich_perm(expr_t, ctrl_genes, expr_ref = expr_all,
                                 n_perm = N_PERM, nbin = NBIN, n_det_bin = N_DET_BIN,
                                 seed = SEED, sig_name = "control_all")
      rc_x    <- sig_enrich_perm(expr_t, ctrl_genes, expr_ref = expr_nt,
                                 n_perm = N_PERM, nbin = NBIN, n_det_bin = N_DET_BIN,
                                 seed = SEED, sig_name = "control_cross")
      for (nm in c("self_TdTomPos", "all_cells", "cross_TdTomNeg")) {
        r <- switch(nm, self_TdTomPos = rc_self, all_cells = rc_all, cross_TdTomNeg = rc_x)
        calib_rows[[length(calib_rows) + 1]] <- data.frame(
          check = "population_restriction_control", population = nm,
          signature = "control_signature", n_genes = length(ctrl_genes),
          n_cells = ncol(expr_t),
          pct_called_raw_p = round(100 * mean(r$pval   < ALPHA, na.rm = TRUE), 2),
          pct_called_gate  = round(100 * mean(r[[GATE]] < ALPHA, na.rm = TRUE), 2),
          median_z         = round(stats::median(r$z, na.rm = TRUE), 3),
          expected_pct     = NA_real_, stringsAsFactors = FALSE)
      }
      message("  [control] % of TdTom+ cells called, by strata reference:")
      message(sprintf("     self (used)  %5.2f%%   all cells %5.2f%%   cross TdTom- %5.2f%%",
                      100 * mean(rc_self[[GATE]] < ALPHA, na.rm = TRUE),
                      100 * mean(rc_all[[GATE]]  < ALPHA, na.rm = TRUE),
                      100 * mean(rc_x[[GATE]]    < ALPHA, na.rm = TRUE)))
      message("     A control signature that is tumour-restricted but not state-specific")
      message("     should NOT fire in nearly all cells. If 'self' does, investigate.")
    }
  } else if (!is.null(CONTROL_SIG)) {
    warning("--control_signature file not found: ", CONTROL_SIG)
  }

  if (length(calib_rows)) {
    calib_tbl <- do.call(rbind, calib_rows)
    calib_tbl$sample <- case_id
    write.csv(calib_tbl,
              file.path(OUT_DIR, paste0(case_id, "_null_calibration.csv")),
              row.names = FALSE)
    message("  Calibration table saved.")
  }
}

# =============================================================================
# 6. Assign final classes
# =============================================================================

message("\n=== Step 6: Assigning cell classes ===")

# --- Non-tumour cells: winner-takes-all cell-type assignment -----------------
# Gate on padj_z, rank on z. The cell-type signatures are largely disjoint, so
# winner-takes-all is defensible here (unlike the EC sub-states in Step 6b,
# which share a pan-EC core).
seu$cell_class <- NA_character_
seu$cell_class[nontumour_cells] <- vapply(nontumour_cells, function(cc) {
  q <- q_mat[cc, ct_names]; zz <- z_mat[cc, ct_names]
  # is.finite(zz) is load-bearing, not defensive padding: a cell can clear the
  # gate on an empirical p while having z = NA from a degenerate null, and
  # which.max() on an all-NA vector returns integer(0), which vapply rejects.
  ok <- !is.na(q) & q < ALPHA & is.finite(zz)
  if (!any(ok)) "not_assigned" else ct_names[ok][which.max(zz[ok])]
}, character(1))

# --- Tumour cells: Neftel-faithful state + hybrid assignment -----------------
# v2: the caller receives Z-SCORES (not AddModuleScore values) and BH-adjusted
# z-tail q-values, so ranking is floor-free and delta is dimensionless.
neftel_z_mat <- z_mat[tumour_cells, neftel_names, drop = FALSE]
neftel_q_mat <- q_mat[tumour_cells, neftel_names, drop = FALSE]

build_collapse <- function(nm) {
  fam_of <- function(x) sub("[12](_like)$", "\\1", x)   # MES1_like/MES2_like -> MES_like
  split(nm, fam_of(nm))
}
neftel_collapse <- build_collapse(neftel_names)
message("Hybrid collapse families: ",
        paste(sprintf("%s={%s}", names(neftel_collapse),
                      sapply(neftel_collapse, paste, collapse = ",")), collapse = "  "))

hyb <- assign_states_with_hybrids(
  score_mat    = neftel_z_mat,       # z-scores
  pval_mat     = neftel_q_mat,       # padj_z
  collapse     = neftel_collapse,
  alpha        = ALPHA,
  delta        = "auto",
  delta_null_q = HYBRID_DELTA_Q,
  n_shuffle    = 100,
  pct_floor    = HYBRID_PCT_FLOOR,
  seed         = SEED)
message("Hybrid criterion-(d) delta (auto, q=", HYBRID_DELTA_Q, "): ",
        signif(attr(hyb, "delta"), 4), " [z units — dimensionless]")

seu$cell_class[tumour_cells]   <- hyb$cell_class[match(tumour_cells, hyb$cell)]
seu$primary_sub                <- NA_character_
seu$primary_sub[tumour_cells]  <- hyb$primary_sub[match(tumour_cells, hyb$cell)]
seu$hybrid_class               <- NA_character_
seu$hybrid_class[tumour_cells] <- hyb$hybrid_class[match(tumour_cells, hyb$cell)]
seu$hybrid_sub                 <- NA_character_
seu$hybrid_sub[tumour_cells]   <- hyb$hybrid_sub[match(tumour_cells, hyb$cell)]
# v3: order-free labels. hybrid_class is RANK-ordered, so MES_like|NPC_like and
# NPC_like|MES_like are the same biology counted as two classes. Use the
# _canonical columns for tabulating and testing; the rank-ordered ones are kept
# because they carry which programme scored higher.
seu$hybrid_class_canonical               <- NA_character_
seu$hybrid_class_canonical[tumour_cells] <- canon_pair(seu$hybrid_class[tumour_cells])
seu$hybrid_sub_canonical                 <- NA_character_
seu$hybrid_sub_canonical[tumour_cells]   <- canon_pair(seu$hybrid_sub[tumour_cells])
seu$is_hybrid                  <- FALSE
seu$is_hybrid[tumour_cells]    <- hyb$is_hybrid[match(tumour_cells, hyb$cell)]
seu$state_gap12                <- NA_real_
seu$state_gap12[tumour_cells]  <- hyb$gap_12[match(tumour_cells, hyb$cell)]
seu$state_gap23                <- NA_real_
seu$state_gap23[tumour_cells]  <- hyb$gap_23[match(tumour_cells, hyb$cell)]

seu$state_label <- seu$cell_class
seu$state_label[!is.na(seu$hybrid_class)] <- seu$hybrid_class[!is.na(seu$hybrid_class)]

# carry the decision statistics into metadata
seu$z_max    <- apply(z_mat, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
seu$gate_min <- apply(q_mat, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
seu$pval_min <- apply(pval_mat, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))

message("\nTumour-cell state distribution (primary):")
print(table(seu$cell_class[tumour_cells], useNA = "ifany"))

# =============================================================================
# 6d. v7: graded confidence (state_tier)
# -----------------------------------------------------------------------------
# v6 emitted one binary outcome, so a cell with NO programme evidence and a cell
# whose best programme just missed the multiplicity threshold both arrived as
# unassigned. Those are different situations and the second is usually the bulk
# of them on a sparse panel. cell_class is NOT modified — these are additional
# columns, so nothing downstream changes unless you choose to use them.
#
#   called  gate passed at ALPHA under --gate_unit          -> publication-grade
#   weak    gate missed, but top z > --tier_z_floor and it
#           is the same programme the ranking prefers        -> usable for
#                                                               composition and
#                                                               neighbourhood
#                                                               statistics, not
#                                                               for per-cell claims
#   none    no programme stands out                          -> genuinely unassigned
#
# Report tier composition in any figure that uses "weak" cells, and run the
# headline test on tier == "called" alone as a sensitivity check.
seu$state_tier   <- NA_character_
seu$tier_top     <- NA_character_
seu$tier_top_z   <- NA_real_
seu$tier_top_q   <- NA_real_

.tier_one <- function(rows, cols) {
  if (!length(rows) || !length(cols)) return(invisible(NULL))
  Z <- z_mat[rows, cols, drop = FALSE]
  Q <- q_mat[rows, cols, drop = FALSE]
  top_i <- apply(Z, 1, function(v) if (all(is.na(v))) NA_integer_ else which.max(v))
  ok    <- !is.na(Q) & Q < ALPHA
  n_ok  <- rowSums(ok)
  top_z <- vapply(seq_along(top_i), function(i)
             if (is.na(top_i[i])) NA_real_ else Z[i, top_i[i]], numeric(1))
  top_q <- vapply(seq_along(top_i), function(i)
             if (is.na(top_i[i])) NA_real_ else Q[i, top_i[i]], numeric(1))
  tier <- ifelse(n_ok > 0, "called",
           ifelse(!is.na(top_z) & top_z > TIER_Z_FLOOR, "weak", "none"))
  seu$state_tier[rows] <<- tier
  seu$tier_top[rows]   <<- ifelse(is.na(top_i), NA_character_, cols[top_i])
  seu$tier_top_z[rows] <<- top_z
  seu$tier_top_q[rows] <<- top_q
  invisible(NULL)
}
.tier_one(tumour_cells,    neftel_names)
.tier_one(nontumour_cells, ct_names)

message("\n=== Step 6d: tier composition (v7) ===")
for (grp in list(list("TdTom+ (Neftel states)", tumour_cells),
                 list("TdTom- (cell types)",    nontumour_cells))) {
  rr <- grp[[2]]
  if (!length(rr)) next
  tt <- table(factor(seu$state_tier[rr], levels = c("called", "weak", "none")))
  message("  ", grp[[1]], ": ",
          paste(sprintf("%s %d (%.1f%%)", names(tt), as.integer(tt),
                        100 * as.integer(tt) / length(rr)), collapse = "  "))
}
message("  cell_class is unchanged; state_tier is additive. Cells in tier weak ",
        "are the ones v6 discarded silently.")

# recoverable-fraction diagnostic: how much of the loss is threshold, not signal
.rec <- table(factor(seu$state_tier[c(tumour_cells, nontumour_cells)],
                     levels = c("called", "weak", "none")))
if (sum(.rec) > 0 && .rec[["called"]] + .rec[["weak"]] > 0)
  message("  Of cells with ANY directional evidence, ",
          round(100 * .rec[["called"]] / (.rec[["called"]] + .rec[["weak"]]), 1),
          "% clear the gate. A low number here means the GATE is the binding ",
          "constraint, not the data.")
message("Hybrid pairs among tumour cells (rank-ordered — same pair may appear twice):")
print(sort(table(seu$hybrid_class[tumour_cells]), decreasing = TRUE))
message("Hybrid pairs among tumour cells (canonical — this is the one to report):")
print(sort(table(seu$hybrid_class_canonical[tumour_cells]), decreasing = TRUE))
message("Hybrid fraction of assigned tumour cells: ",
        round(mean(seu$is_hybrid[tumour_cells], na.rm = TRUE), 3))

# --- Proliferation: gate on padj_z, rank on z --------------------------------
seu$prolif_state <- vapply(seq_len(ncol(seu)), function(i) {
  q <- q_prolif[i, ]; zz <- z_prolif[i, ]
  ok <- !is.na(q) & q < ALPHA & is.finite(zz)   # see Step 6: guards which.max(all-NA)
  if (!any(ok)) "undefined" else prolif_names[ok][which.max(zz[ok])]
}, character(1))
seu$gate_prolif_min <- apply(q_prolif, 1, function(x)
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))

seu$combined_class <- paste0(seu$cell_class, "|", seu$prolif_state)

# =============================================================================
# 6b. EC sub-states — Endothelial cells only, UNIQUE genes, reported INDEPENDENTLY
# -----------------------------------------------------------------------------
# v1 used winner-takes-all across the three EC signatures. Those signatures
# share a pan-EC core, and in simulation the shared core made all three fire in
# every cell, so the "winner" was decided by the core rather than by the genes
# that distinguish the sub-states: with the truth set to 100% one state, only
# 22% of cells were assigned correctly — worse than chance among three.
#
# v2 therefore:
#   * tests each sub-state on its UNIQUE genes (--ec_unique_genes, default TRUE);
#   * reports each sub-state independently as ec_<name>_pos / _z / _gate, so a
#     cell may be positive for none, one, or several;
#   * keeps a single ec_state column for plotting/GeoJSON continuity, but it now
#     takes the values "none", the sub-state name (exactly one positive), or
#     "multiple" — it is a convenience summary, NOT an exclusive classification.
# Strata reference is expr_ec (self) — see the note in Step 5b. Widening the
# reference is the wrong fix for the shared-core problem; using unique genes is
# the right one.
# =============================================================================

USE_EC_RESULTS <- FALSE
seu$ec_state     <- NA_character_
seu$gate_ec_min <- NA_real_

if (USE_EC) {
  ec_cells <- rownames(seu@meta.data)[which(seu$cell_class == "Endothelial")]
  message("\n=== Step 6b: EC sub-states on Endothelial cells (n=",
          length(ec_cells), ") ===")

  if (length(ec_cells) < 10) {
    warning("Fewer than 10 Endothelial cells — skipping EC permutation.")
  } else {
    USE_EC_RESULTS <- TRUE
    expr_ec <- expr_all[, ec_cells, drop = FALSE]

    ec_test_sets <- if (EC_UNIQUE) ec_sig_unique else ec_sig_mouse
    n_unique     <- sapply(ec_test_sets, length)
    dropped      <- names(n_unique)[n_unique < EC_MIN_UNIQUE]
    if (EC_UNIQUE && length(dropped)) {
      warning("EC sub-states with <", EC_MIN_UNIQUE, " unique genes, falling back ",
              "to full signature (shared core will inflate them): ",
              paste(dropped, collapse = ", "))
      for (d in dropped) ec_test_sets[[d]] <- ec_sig_mouse[[d]]
    }
    message("  Genes used per sub-state (", if (EC_UNIQUE) "unique" else "full", "):")
    print(sapply(ec_test_sets, length))

    z_ec <- mk_mat(ec_cells, ec_names)
    q_ec <- mk_mat(ec_cells, ec_names)
    pval_ec <- mk_mat(ec_cells, ec_names)

    for (s in ec_names) {
      g <- ec_test_sets[[s]]
      if (length(g) < 3) {
        warning("EC sub-state ", s, " has <3 testable genes — skipped.")
        next
      }
      message("  ", s, " (n=", length(g), ") — Endothelial cells")
      res <- sig_enrich_perm(expr_ec, g, n_perm = N_PERM,
                             nbin = NBIN, n_det_bin = N_DET_BIN,
                             seed = SEED, sig_name = s)
      m <- match(ec_cells, res$cell)
      z_ec[, s]    <- res$z[m]
      q_ec[, s]    <- res[[GATE]][m]
      pval_ec[, s] <- res$pval[m]
    }

    # --- independent per-signature calls (the primary EC output) -------------
    for (s in ec_names) {
      seu@meta.data[[paste0("ec_", s, "_z")]]     <- NA_real_
      seu@meta.data[[paste0("ec_", s, "_gate")]]  <- NA_real_
      seu@meta.data[[paste0("ec_", s, "_pos")]]   <- NA
      seu@meta.data[ec_cells, paste0("ec_", s, "_z")]     <- z_ec[, s]
      seu@meta.data[ec_cells, paste0("ec_", s, "_gate")]  <- q_ec[, s]
      seu@meta.data[ec_cells, paste0("ec_", s, "_pos")]   <-
        !is.na(q_ec[, s]) & q_ec[, s] < ALPHA
    }

    n_pos <- rowSums(!is.na(q_ec) & q_ec < ALPHA)
    seu$ec_n_positive <- NA_integer_
    seu$ec_n_positive[ec_cells] <- n_pos

    seu$ec_state[ec_cells] <- ifelse(
      n_pos == 0, "none",
      ifelse(n_pos == 1,
             ec_names[max.col(matrix(!is.na(q_ec) & q_ec < ALPHA,
                                     nrow = length(ec_cells)), ties.method = "first")],
             "multiple"))
    seu$gate_ec_min[ec_cells] <- apply(q_ec, 1, function(x)
      if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))

    message("\nEC sub-state positivity (independent tests, Endothelial cells only):")
    print(data.frame(
      sub_state = ec_names,
      n_genes   = sapply(ec_names, function(s) length(ec_test_sets[[s]])),
      n_pos     = sapply(ec_names, function(s) sum(!is.na(q_ec[, s]) & q_ec[, s] < ALPHA)),
      pct_pos   = round(100 * sapply(ec_names, function(s)
        mean(!is.na(q_ec[, s]) & q_ec[, s] < ALPHA)), 2),
      row.names = NULL))
    message("Cells positive for 0 / 1 / >1 sub-states:")
    print(table(factor(ifelse(n_pos == 0, "0", ifelse(n_pos == 1, "1", ">1")),
                       levels = c("0", "1", ">1"))))
    message("NOTE: ec_state is a convenience summary. The primary output is the")
    message("      per-signature ec_<name>_pos columns — sub-states are NOT exclusive.")

    ec_tbl <- data.frame(cell_id = ec_cells, stringsAsFactors = FALSE)
    for (s in ec_names) {
      ec_tbl[[paste0(s, "_z")]]     <- z_ec[, s]
      ec_tbl[[paste0(s, "_gate")]]  <- q_ec[, s]
      ec_tbl[[paste0(s, "_pos")]]   <- !is.na(q_ec[, s]) & q_ec[, s] < ALPHA
    }
    ec_tbl$n_positive <- n_pos
    ec_tbl$sample     <- case_id
    write.csv(ec_tbl, file.path(OUT_DIR, paste0(case_id, "_ec_substates.csv")),
              row.names = FALSE)
    message("  EC per-cell sub-state table saved.")
  }
}

# =============================================================================
# 6c. MYELOID STATE SUB-TYPING — microglia / macrophages only
# -----------------------------------------------------------------------------
# Identity vs state. Step 6 decides identity (Microglia_core /
# Macrophage_peripheral / Macrophage_BAM) from celltype_sig_mouse. This step
# takes those cells as given and asks which transcriptional STATE each one is
# running, using two published state-signature sets scored independently:
#
#   Miller 2025 (human symbols, converted in Step 3)
#   Qian  2026 (mouse symbols)
#
# For each set, per state: AddModuleScore (descriptive, for plotting) +
# sig_enrich_perm() z and BH-adjusted tail p (the inference), then
# assign_states_with_hybrids() for single-state and hybrid calls. Structurally
# this is Step 6b applied to a different compartment, with three differences,
# each of which is a decision you can reverse:
#
# (1) NULL REFERENCE is expr_mye — the myeloid compartment itself, not all
#     cells. This is the Step 5b failure mode: a signature restricted to one
#     population beats an all-cell null in every cell of that population, so a
#     widened reference gives ~100% positives and measures nothing. Do not
#     change this.
#
# (2) CRITERION (d) IS DISABLED by default (--myeloid_hybrid_delta 0), unlike
#     the Neftel block. The auto-calibrated gap_23 delta is calibrated by
#     shuffling score columns, so a stronger signature widens its own null and
#     raises the threshold in step with the signal the threshold is meant to
#     admit. Simulation with known truth on this pipeline (600 cells, 4 disjoint
#     signatures, planted single states + planted hybrids + null cells):
#
#       criterion (d)                delta   true hybrids found   FP singles   FP null
#       "auto" (Neftel-identical)     3.38          0.15               0          0
#       scoped to gate-positive       3.83          0.08               0          0
#       0  (this default)             0             0.43               0          0
#
#     Specificity is carried by criterion (b) — both programmes must clear the
#     BH gate independently — not by (d). Recalibrating (d) on candidate cells
#     was tested and is worse. Set --myeloid_hybrid_delta auto for
#     Neftel-identical behaviour; the Neftel block itself is unchanged.
#
# (3) STATES ARE REPORTED INDEPENDENTLY as myeloid_<set>_<state>_pos, in
#     addition to the hybrid call. Same reasoning as Step 6b: a cell may run
#     none, one, or several programmes, and the summary column is a convenience.
#
# POWER. Hybrid recall is roughly the square of single-state power because of
# criterion (b). ~80-gene sets reach single-state sensitivity ~0.78 / hybrid
# recall ~0.54; ~25-gene sets reach ~0.50 / ~0.22, and only 0.008 hybrid recall
# at a 4x effect. False positives never exceeded 0.53% against a nominal 5% at
# any size, so positives are trustworthy and negatives are not. Report hybrid
# counts from small sets as a floor, and never compare hybrid frequency between
# the two sets.
#
# Myeloid cells misassigned as Neftel states among TdTom+ cells are NOT tested
# here — this step reads cell_class, which is a Neftel state for tumour cells.
# =============================================================================

USE_MYELOID_RESULTS <- FALSE
myeloid_tbl_list    <- list()
myeloid_summary     <- NULL

if (USE_MYELOID) {
  myeloid_cells <- rownames(seu@meta.data)[which(seu$cell_class %in% MYELOID_CLASSES)]
  message("\n=== Step 6c: Myeloid state sub-typing (n=", length(myeloid_cells),
          " cells in ", paste(MYELOID_CLASSES, collapse = "/"), ") ===")
  message("  Identity composition of the compartment being sub-typed:")
  print(table(seu$cell_class[myeloid_cells]))

  if (length(myeloid_cells) < MYELOID_MIN_CELLS) {
    warning("Fewer than ", MYELOID_MIN_CELLS, " myeloid cells (", length(myeloid_cells),
            ") — skipping Step 6c. The permutation null needs cells to be stable.")
  } else {
    USE_MYELOID_RESULTS <- TRUE
    expr_mye <- expr_all[, myeloid_cells, drop = FALSE]

    for (setnm in names(myeloid_sig_mouse)) {
      sets <- myeloid_sig_mouse[[setnm]]
      sets <- sets[sapply(sets, length) >= MYELOID_MIN_GENES]
      if (!length(sets)) {
        warning("Set '", setnm, "': no state has >=", MYELOID_MIN_GENES,
                " testable genes — skipped.")
        next
      }
      st_names <- names(sets)
      message("\n--- ", setnm, ": ", length(st_names), " states, genes tested per state ---")
      print(sapply(sets, length))

      # --- collinearity report on the genes ACTUALLY tested -----------------
      # The collapse specification below was decided on the published sets. After
      # orthology and identity stripping the overlaps change, so re-check them
      # here rather than trusting the published numbers.
      if (length(st_names) >= 2) {
        cmb <- combn(st_names, 2, simplify = FALSE)
        ov  <- do.call(rbind, lapply(cmb, function(k) {
          i <- intersect(sets[[k[1]]], sets[[k[2]]])
          u <- union(sets[[k[1]]], sets[[k[2]]])
          data.frame(state_a = k[1], state_b = k[2], n_shared = length(i),
                     jaccard = if (length(u)) round(length(i) / length(u), 3) else NA_real_,
                     stringsAsFactors = FALSE)
        }))
        ov <- ov[order(-ov$jaccard), ]
        flag <- ov[ov$jaccard > 0.25, ]
        if (nrow(flag)) {
          message("  Gene sharing above Jaccard 0.25 — these states cannot be RANKED",
                  " against each other:")
          print(flag, row.names = FALSE)
        } else {
          message("  No state pair shares more than Jaccard 0.25 (max ",
                  max(ov$jaccard, na.rm = TRUE), ") — ranking is meaningful.")
        }
        write.csv(ov, file.path(OUT_DIR, paste0(case_id, "_myeloid_", setnm,
                                               "_state_collinearity.csv")),
                  row.names = FALSE)
      }

      # --- AddModuleScore: descriptive only, kept for plotting --------------
      # Computed on the myeloid subset so the control-gene bins are drawn from
      # the same cells the permutation test uses.
      seu_mye <- subset(seu, cells = myeloid_cells)
      seu_mye <- AddModuleScore(seu_mye, features = sets, name = paste0("ms_", setnm, "_"),
                                nbin = min(24, NBIN), ctrl = 100, seed = SEED)
      ms_cols <- grep(paste0("^ms_", setnm, "_[0-9]+$"), colnames(seu_mye@meta.data),
                      value = TRUE)
      ms_cols <- ms_cols[order(as.integer(sub(".*_", "", ms_cols)))]

      # --- permutation test per state --------------------------------------
      z_my <- mk_mat(myeloid_cells, st_names)
      q_my <- mk_mat(myeloid_cells, st_names)
      p_my <- mk_mat(myeloid_cells, st_names)

      for (st in st_names) {
        message("  ", setnm, ":", st, " (n=", length(sets[[st]]), " genes)")
        res <- sig_enrich_perm(expr_mye, sets[[st]], n_perm = N_PERM,
                               nbin = NBIN, n_det_bin = N_DET_BIN,
                               seed = SEED, sig_name = paste0(setnm, ":", st))
        m <- match(myeloid_cells, res$cell)
        z_my[, st] <- res$z[m]
        q_my[, st] <- res[[GATE]][m]
        p_my[, st] <- res$pval[m]
      }

      # --- collapse near-duplicate states into families --------------------
      coll_spec <- MYELOID_COLLAPSE[[setnm]]
      collapse  <- as.list(setNames(st_names, st_names))
      if (!is.null(coll_spec)) {
        # Inherit the study prefix of the states being collapsed, so the family
        # name reads like its members ("qian26_MG-Prolif", not "MG-Prolif").
        st_prefix <- unique(regmatches(st_names, regexpr("^[A-Za-z]+[0-9]*_", st_names)))
        st_prefix <- if (length(st_prefix) == 1) st_prefix else ""
        for (fam in names(coll_spec)) {
          members  <- match_state_names(coll_spec[[fam]], st_names)
          fam_name <- paste0(st_prefix, fam)
          if (length(members) >= 2) {
            collapse[members]   <- NULL
            collapse[[fam_name]] <- members
            message("  Collapsed to family '", fam_name, "': ",
                    paste(members, collapse = " + "),
                    " (near-duplicate gene sets; see header)")
          } else if (length(members) == 1) {
            message("  Collapse family '", fam_name, "' has only one surviving member (",
                    members, ") — left uncollapsed.")
          } else {
            warning("Collapse family '", fam_name, "' matched no states — not collapsed.")
          }
        }
      }

      # --- hybrid / single-state calls --------------------------------------
      hyb_my <- assign_states_with_hybrids(
        score_mat  = z_my,
        pval_mat   = q_my,
        collapse   = collapse,
        alpha      = ALPHA,
        delta      = MYELOID_DELTA,
        pct_floor  = HYBRID_PCT_FLOOR,
        delta_null_q = HYBRID_DELTA_Q,
        seed       = SEED)
      hyb_my <- hyb_my[match(myeloid_cells, hyb_my$cell), ]
      message("  Criterion-(d) delta used: ",
              if (identical(MYELOID_DELTA, "auto"))
                paste0(signif(attr(hyb_my, "delta"), 4), " (auto)")
              else paste0(MYELOID_DELTA, " (criterion disabled — see header)"))

      # --- write back: independent per-state calls (primary output) ---------
      pref <- paste0("myeloid_", setnm, "_")
      for (i in seq_along(st_names)) {
        st <- st_names[i]
        seu@meta.data[[paste0(pref, st, "_score")]] <- NA_real_
        seu@meta.data[[paste0(pref, st, "_z")]]     <- NA_real_
        seu@meta.data[[paste0(pref, st, "_gate")]]  <- NA_real_
        seu@meta.data[[paste0(pref, st, "_pos")]]   <- NA
        seu@meta.data[myeloid_cells, paste0(pref, st, "_score")] <-
          seu_mye@meta.data[[ms_cols[i]]][match(myeloid_cells, colnames(seu_mye))]
        seu@meta.data[myeloid_cells, paste0(pref, st, "_z")]    <- z_my[, st]
        seu@meta.data[myeloid_cells, paste0(pref, st, "_gate")] <- q_my[, st]
        seu@meta.data[myeloid_cells, paste0(pref, st, "_pos")]  <-
          !is.na(q_my[, st]) & q_my[, st] < ALPHA
      }

      # --- write back: state / hybrid summary -------------------------------
      n_pos_my <- rowSums(!is.na(q_my) & q_my < ALPHA)
      seu@meta.data[[paste0(pref, "n_positive")]] <- NA_integer_
      seu@meta.data[myeloid_cells, paste0(pref, "n_positive")] <- n_pos_my

      for (cl in c("state", "state_canonical", "hybrid", "hybrid_canonical")) {
        seu@meta.data[[paste0(pref, cl)]] <- NA_character_
      }
      seu@meta.data[[paste0(pref, "is_hybrid")]] <- NA
      seu@meta.data[myeloid_cells, paste0(pref, "state")]  <- hyb_my$cell_class
      seu@meta.data[myeloid_cells, paste0(pref, "hybrid")] <- hyb_my$hybrid_class
      seu@meta.data[myeloid_cells, paste0(pref, "hybrid_canonical")] <-
        canon_pair(hyb_my$hybrid_class)
      seu@meta.data[myeloid_cells, paste0(pref, "is_hybrid")] <- hyb_my$is_hybrid
      # state_canonical = the label to plot: hybrid where called, else single state
      lab_my <- ifelse(!is.na(hyb_my$hybrid_class),
                       canon_pair(hyb_my$hybrid_class), hyb_my$cell_class)
      seu@meta.data[myeloid_cells, paste0(pref, "state_canonical")] <- lab_my

      # --- log ---------------------------------------------------------------
      st_summary <- data.frame(
        set     = setnm,
        state   = st_names,
        n_genes = sapply(st_names, function(x) length(sets[[x]])),
        n_pos   = sapply(st_names, function(x) sum(!is.na(q_my[, x]) & q_my[, x] < ALPHA)),
        pct_pos = round(100 * sapply(st_names, function(x)
                    mean(!is.na(q_my[, x]) & q_my[, x] < ALPHA)), 2),
        median_z = round(sapply(st_names, function(x) median(z_my[, x], na.rm = TRUE)), 3),
        row.names = NULL, stringsAsFactors = FALSE)
      message("\n  ", setnm, " — per-state positivity (independent tests):")
      print(st_summary, row.names = FALSE)
      myeloid_summary <- rbind(myeloid_summary, st_summary)

      message("  Cells positive for 0 / 1 / >1 states:")
      print(table(factor(ifelse(n_pos_my == 0, "0", ifelse(n_pos_my == 1, "1", ">1")),
                         levels = c("0", "1", ">1"))))
      message("  Hybrid rate: ", round(mean(hyb_my$is_hybrid, na.rm = TRUE), 3),
              "   (canonical hybrid pairs:)")
      print(sort(table(canon_pair(hyb_my$hybrid_class)), decreasing = TRUE))

      # State x identity. The state call should NOT simply recapitulate identity;
      # if it does, the state signature is carrying identity information.
      message("  State (canonical) x myeloid identity:")
      print(table(state = lab_my, identity = seu$cell_class[myeloid_cells]))

      # Proliferation cross-check: collapsed cell-cycle families duplicate the
      # genome-wide test in Step 6, so verify rather than double-report.
      prolif_fams <- grep("Prolif", unique(na.omit(lab_my)), value = TRUE)
      if (length(prolif_fams)) {
        message("  Proliferation family x genome-wide prolif_state ",
                "(these test the same biology — do not report both as findings):")
        print(table(myeloid_state = lab_my, prolif_state = seu$prolif_state[myeloid_cells]))
      }

      # --- per-cell table ---------------------------------------------------
      mt <- data.frame(cell_id = myeloid_cells,
                       identity = seu$cell_class[myeloid_cells],
                       stringsAsFactors = FALSE)
      for (i in seq_along(st_names)) {
        st <- st_names[i]
        mt[[paste0(st, "_score")]] <-
          seu_mye@meta.data[[ms_cols[i]]][match(myeloid_cells, colnames(seu_mye))]
        mt[[paste0(st, "_z")]]     <- z_my[, st]
        mt[[paste0(st, "_gate")]]  <- q_my[, st]
        mt[[paste0(st, "_pval")]]  <- p_my[, st]
        mt[[paste0(st, "_pos")]]   <- !is.na(q_my[, st]) & q_my[, st] < ALPHA
      }
      mt$n_positive       <- n_pos_my
      mt$state            <- hyb_my$cell_class
      mt$hybrid           <- hyb_my$hybrid_class
      mt$hybrid_canonical <- canon_pair(hyb_my$hybrid_class)
      mt$is_hybrid        <- hyb_my$is_hybrid
      mt$state_canonical  <- lab_my
      mt$gap_12           <- hyb_my$gap_12
      mt$gap_23           <- hyb_my$gap_23
      mt$prolif_state     <- seu$prolif_state[myeloid_cells]
      mt$tdtom_label      <- seu$tdtom_label[myeloid_cells]
      mt$sample           <- case_id
      myeloid_tbl_list[[setnm]] <- mt
      write.csv(mt, file.path(OUT_DIR, paste0(case_id, "_myeloid_states_", setnm, ".csv")),
                row.names = FALSE)
      message("  Per-cell table saved: ", case_id, "_myeloid_states_", setnm, ".csv")

      rm(seu_mye)
    }

    if (!is.null(myeloid_summary))
      write.csv(myeloid_summary,
                file.path(OUT_DIR, paste0(case_id, "_myeloid_state_summary.csv")),
                row.names = FALSE)

    # Cross-set agreement. The two sets are not comparable in power (see header),
    # so this is a qualitative check, not a concordance statistic.
    if (length(myeloid_tbl_list) == 2) {
      s1 <- names(myeloid_tbl_list)[1]; s2 <- names(myeloid_tbl_list)[2]
      message("\n  Cross-set state agreement (", s1, " rows x ", s2,
              " cols). Power differs between sets — read this qualitatively:")
      print(table(myeloid_tbl_list[[s1]]$state_canonical,
                  myeloid_tbl_list[[s2]]$state_canonical))
    }

    message("\nNOTE: myeloid_<set>_state_canonical is a convenience summary. The primary")
    message("      output is the per-state myeloid_<set>_<state>_pos columns — states")
    message("      are NOT mutually exclusive.")
  }
}

message("\nCell class distribution:")
print(table(seu$cell_class, useNA = "ifany"))
message("\nProliferation state distribution:")
print(table(seu$prolif_state))
message("\nTop combined classes:")
print(head(sort(table(seu$combined_class), decreasing = TRUE), 20))

# =============================================================================
# 7. DE analysis — TdTom+ vs TdTom-
# =============================================================================

message("\n=== Step 7: DE analysis TdTom+ vs TdTom- ===")

Idents(seu) <- seu$tdtom_label
de_results  <- FindMarkers(seu,
                           ident.1     = "TdTom_pos",
                           ident.2     = "TdTom_neg",
                           assay       = "SCT",
                           test.use    = "wilcox",
                           min.pct     = 0.1,
                           logfc.threshold = 0.25,
                           verbose     = FALSE)
de_results$gene      <- rownames(de_results)
de_results$direction <- ifelse(de_results$avg_log2FC > 0, "TdTom_pos", "TdTom_neg")

write.csv(de_results,
          file.path(OUT_DIR, paste0(case_id, "_DE_tdtom_pos_vs_neg.csv")),
          row.names = FALSE)
message("  DE genes found: ", nrow(de_results))
message("  TdTom+ enriched: ", sum(de_results$direction == "TdTom_pos" & de_results$p_val_adj < ALPHA))
message("  TdTom- enriched: ", sum(de_results$direction == "TdTom_neg" & de_results$p_val_adj < ALPHA))

# =============================================================================
# 8. Plots
# =============================================================================

message("\n=== Step 8: Generating plots ===")
pdf_path <- file.path(OUT_DIR, paste0(case_id, "_classification_plots.pdf"))
pdf(pdf_path, width = 14, height = 8)

# 8a. Decision-statistic distributions ---------------------------------------
# v2 plots BOTH the raw empirical p (to make the permutation floor visible) and
# the z-scores that actually drive the calls.
p_pval <- as.data.frame(pval_mat) %>%
  pivot_longer(everything(), names_to = "signature", values_to = "pval") %>%
  drop_na() %>%
  ggplot(aes(x = pval, fill = signature)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  geom_vline(xintercept = ALPHA, linetype = "dashed", colour = "firebrick") +
  geom_vline(xintercept = 1 / (N_PERM + 1), linetype = "dotted", colour = "grey30") +
  facet_wrap(~ signature, scales = "free_y", ncol = 4) +
  labs(title = paste0(case_id, " — empirical permutation p-values (Neftel + cell-type)"),
       subtitle = paste0("dotted line = permutation floor 1/(n_perm+1) = ",
                         signif(1 / (N_PERM + 1), 3),
                         "; calls are NOT made on this statistic"),
       x = "Empirical p-value", y = "Count") +
  theme_bw(base_size = 9) + theme(legend.position = "none")
print(p_pval)

p_z <- as.data.frame(z_mat) %>%
  pivot_longer(everything(), names_to = "signature", values_to = "z") %>%
  drop_na() %>%
  ggplot(aes(x = z, fill = signature)) +
  geom_histogram(bins = 60, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  facet_wrap(~ signature, scales = "free_y", ncol = 4) +
  labs(title = paste0(case_id, " — enrichment z vs matched null (decision statistic)"),
       subtitle = "z = (signature mean - null mean) / null sd; null matched on expression x detection",
       x = "z", y = "Count") +
  theme_bw(base_size = 9) + theme(legend.position = "none")
print(p_z)

p_pval_prolif <- as.data.frame(pval_prolif) %>%
  pivot_longer(everything(), names_to = "signature", values_to = "pval") %>%
  ggplot(aes(x = pval, fill = signature)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  geom_vline(xintercept = ALPHA, linetype = "dashed", colour = "firebrick") +
  facet_wrap(~ signature, scales = "free_y") +
  labs(title = paste0(case_id, " — permutation p-values (proliferation)"),
       x = "Empirical p-value", y = "Count") +
  theme_bw(base_size = 9) + theme(legend.position = "none")
print(p_pval_prolif)

# 8b. UMAP --------------------------------------------------------------------
umap_df <- as.data.frame(Embeddings(seu, "umap")) %>%
  rownames_to_column("cell_id") %>%
  left_join(
    seu@meta.data %>%
      rownames_to_column("cell_id") %>%
      dplyr::select(cell_id, cell_class, prolif_state, combined_class,
                    tdtom_label, seurat_clusters, state_label, is_hybrid,
                    any_of(if (USE_EC) c("ec_state", "ec_n_positive") else character(0))),
    by = "cell_id"
  )

cluster_centres <- umap_df %>%
  group_by(seurat_clusters) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

p_class <- ggplot(umap_df, aes(umap_1, umap_2, colour = cell_class)) +
  geom_point(size = 0.3, alpha = 0.7) +
  scale_colour_manual(values = class_colours) +
  guides(colour = guide_legend(override.aes = list(size = 4), ncol = 2)) +
  labs(title = "Cell class", colour = NULL) + theme_classic()

p_prolif <- ggplot(umap_df, aes(umap_1, umap_2, colour = prolif_state)) +
  geom_point(size = 0.3, alpha = 0.7) +
  scale_colour_manual(values = prolif_colours) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  labs(title = "Proliferation state", colour = NULL) + theme_classic()

p_tdtom <- ggplot(umap_df, aes(umap_1, umap_2, colour = tdtom_label)) +
  geom_point(size = 0.3, alpha = 0.7) +
  scale_colour_manual(values = c("TdTom_pos" = "blue", "TdTom_neg" = "lightgrey")) +
  labs(title = "TdTomato label", colour = NULL) + theme_classic()

p_cl <- ggplot(umap_df, aes(umap_1, umap_2, colour = factor(seurat_clusters))) +
  geom_point(size = 0.3, alpha = 0.5) +
  geom_label_repel(data = cluster_centres, aes(label = seurat_clusters),
                   colour = "black", fill = "white", size = 3,
                   fontface = "bold", max.overlaps = Inf) +
  guides(colour = "none") +
  labs(title = "Clusters") + theme_classic()

print((p_class | p_prolif) / (p_tdtom | p_cl))

# 8c. Spatial plots -----------------------------------------------------------
if (!all(c("imagerow", "imagecol") %in% colnames(seu@meta.data))) {
  message("  imagerow/imagecol not in metadata — fetching from tissue coordinates")
  tissue_coords <- GetTissueCoordinates(seu)
  seu@meta.data$imagerow <- tissue_coords[rownames(seu@meta.data), "imagerow"]
  seu@meta.data$imagecol <- tissue_coords[rownames(seu@meta.data), "imagecol"]
}

sp_df <- seu@meta.data %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, cell_class, prolif_state, combined_class,
                tdtom_label, imagerow, imagecol,
                any_of(if (USE_EC) "ec_state" else character(0)))

p_sp_class <- ggplot(sp_df, aes(x = imagerow, y = imagecol, colour = cell_class)) +
  geom_point(size = 0.4, alpha = 0.8, shape = 16) +
  scale_colour_manual(values = class_colours) +
  scale_x_reverse() + coord_fixed() + theme_void() +
  labs(title = paste0(case_id, " — cell class"), colour = NULL)

p_sp_prolif <- ggplot(sp_df, aes(x = imagerow, y = imagecol, colour = prolif_state)) +
  geom_point(size = 0.4, alpha = 0.8, shape = 16) +
  scale_colour_manual(values = prolif_colours) +
  scale_x_reverse() + coord_fixed() + theme_void() +
  labs(title = paste0(case_id, " — proliferation state"), colour = NULL)

p_sp_tdtom <- ggplot(sp_df, aes(x = imagerow, y = imagecol, colour = tdtom_label)) +
  geom_point(size = 0.4, alpha = 0.8, shape = 16) +
  scale_colour_manual(values = c("TdTom_pos" = "blue", "TdTom_neg" = "lightgrey")) +
  scale_x_reverse() + coord_fixed() + theme_void() +
  labs(title = paste0(case_id, " — TdTomato label"), colour = NULL)

print(p_sp_class)
print(p_sp_prolif)
print(p_sp_tdtom)

# 8d. EC sub-states — independent, not winner-takes-all -----------------------
if (USE_EC && USE_EC_RESULTS) {

  ec_long <- ec_tbl %>%
    dplyr::select(cell_id, ends_with("_z")) %>%
    pivot_longer(-cell_id, names_to = "signature", values_to = "z") %>%
    mutate(signature = sub("_z$", "", signature))

  p_ec_z <- ggplot(ec_long, aes(x = z, fill = signature)) +
    geom_histogram(bins = 50, alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
    facet_wrap(~ signature, scales = "free_y") +
    scale_fill_manual(values = ec_colours) +
    labs(title = paste0(case_id, " — EC sub-state z (",
                        if (EC_UNIQUE) "unique genes" else "full signatures", ")"),
         subtitle = "Sub-states are tested and reported INDEPENDENTLY; a cell may be positive for several",
         x = "z", y = "Count") +
    theme_bw(base_size = 9) + theme(legend.position = "none")
  print(p_ec_z)

  # positivity co-occurrence
  pos_mat <- ec_tbl %>% dplyr::select(ends_with("_pos"))
  colnames(pos_mat) <- sub("_pos$", "", colnames(pos_mat))
  combo <- apply(pos_mat, 1, function(r) {
    on <- colnames(pos_mat)[which(r)]
    if (!length(on)) "none" else paste(on, collapse = " + ")
  })
  p_ec_combo <- as.data.frame(table(combo), stringsAsFactors = FALSE) %>%
    setNames(c("combination", "n_cells")) %>%
    arrange(desc(n_cells)) %>%
    ggplot(aes(x = reorder(combination, n_cells), y = n_cells)) +
    geom_col(fill = "grey40", width = 0.65) +
    geom_text(aes(label = n_cells), hjust = -0.2, size = 3) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = paste0(case_id, " — EC sub-state co-positivity"),
         subtitle = "Cells positive for multiple sub-states are the reason winner-takes-all was dropped",
         x = NULL, y = "Cells") +
    theme_bw(base_size = 9)
  print(p_ec_combo)

  p_ec <- ggplot(umap_df, aes(umap_1, umap_2, colour = ec_state)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_colour_manual(values = ec_colours, na.value = "grey90") +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    labs(title = "EC state summary (Endothelial cells only)",
         subtitle = "'multiple' = positive for >1 sub-state",
         colour = NULL) +
    theme_classic()
  print(p_ec)

  p_sp_ec <- ggplot(sp_df, aes(x = imagerow, y = imagecol, colour = ec_state)) +
    geom_point(size = 0.4, alpha = 0.8, shape = 16) +
    scale_colour_manual(values = ec_colours, na.value = "grey90") +
    scale_x_reverse() + coord_fixed() + theme_void() +
    labs(title = paste0(case_id, " — EC state (Endothelial only)"), colour = NULL)
  print(p_sp_ec)
}

# 8g. Myeloid states — independent, not winner-takes-all ----------------------
if (USE_MYELOID && USE_MYELOID_RESULTS && length(myeloid_tbl_list)) {

  for (setnm in names(myeloid_tbl_list)) {
    mt <- myeloid_tbl_list[[setnm]]
    st_names <- sub("_z$", "", grep("_z$", colnames(mt), value = TRUE))

    # z distributions per state
    my_long <- mt %>%
      dplyr::select(cell_id, ends_with("_z")) %>%
      pivot_longer(-cell_id, names_to = "state", values_to = "z") %>%
      mutate(state = sub("_z$", "", state))

    p_my_z <- ggplot(my_long, aes(x = z)) +
      geom_histogram(bins = 40, fill = "grey35", alpha = 0.9) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
      facet_wrap(~ state, scales = "free_y") +
      labs(title = paste0(case_id, " — myeloid state z (", setnm, ")"),
           subtitle = paste("Microglia/macrophages only; null referenced on the myeloid",
                            "compartment. States tested INDEPENDENTLY."),
           x = "z (permutation)", y = "Cells") +
      theme_bw(base_size = 9)
    print(p_my_z)

    # module score vs z — shows what the descriptive score would have implied
    ms_long <- mt %>%
      dplyr::select(cell_id, ends_with("_score"), ends_with("_z")) %>%
      pivot_longer(-cell_id, names_to = "k", values_to = "v") %>%
      mutate(state = sub("_(score|z)$", "", k),
             what  = ifelse(grepl("_score$", k), "module_score", "z")) %>%
      dplyr::select(-k) %>%
      pivot_wider(names_from = what, values_from = v)

    p_ms_z <- ggplot(ms_long, aes(module_score, z)) +
      geom_point(size = 0.35, alpha = 0.5, colour = "grey25") +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      facet_wrap(~ state, scales = "free") +
      labs(title = paste0(case_id, " — module score vs permutation z (", setnm, ")"),
           subtitle = paste("Module score is descriptive and has no null; z carries the",
                            "inference. Cells at high score but z ~ 0 are the reason."),
           x = "AddModuleScore", y = "z") +
      theme_bw(base_size = 9)
    print(p_ms_z)

    # co-positivity
    pos_my <- mt %>% dplyr::select(ends_with("_pos"))
    colnames(pos_my) <- sub("_pos$", "", colnames(pos_my))
    combo_my <- apply(pos_my, 1, function(r) {
      on <- colnames(pos_my)[which(r)]
      if (!length(on)) "none" else paste(sort(on), collapse = " + ")
    })
    p_my_combo <- as.data.frame(table(combo_my), stringsAsFactors = FALSE) %>%
      setNames(c("combination", "n_cells")) %>%
      arrange(desc(n_cells)) %>% head(20) %>%
      ggplot(aes(x = reorder(combination, n_cells), y = n_cells)) +
      geom_col(fill = "grey40", width = 0.65) +
      geom_text(aes(label = n_cells), hjust = -0.2, size = 3) +
      coord_flip(clip = "off") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(title = paste0(case_id, " — myeloid state co-positivity (", setnm, ")"),
           subtitle = "Cells positive for several states are why the summary column is not a classification",
           x = NULL, y = "Cells") +
      theme_bw(base_size = 9)
    print(p_my_combo)

    # state x identity — the state call should not just recapitulate identity
    p_my_id <- mt %>%
      count(state_canonical, identity) %>%
      ggplot(aes(identity, state_canonical, fill = n)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = n), size = 2.8) +
      scale_fill_gradient(low = "grey93", high = "#1F6FB2") +
      labs(title = paste0(case_id, " — myeloid state x identity (", setnm, ")"),
           subtitle = "Identity is from Step 6; state from Step 6c. Strong diagonal structure would mean the state signature carries identity.",
           x = "identity (Step 6)", y = "state (Step 6c, canonical)", fill = "Cells") +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    print(p_my_id)

    # spatial + UMAP of the summary label
    my_lab_col <- paste0("myeloid_", setnm, "_state_canonical")
    if (my_lab_col %in% colnames(seu@meta.data)) {
      umap_df[[my_lab_col]] <- seu@meta.data[[my_lab_col]][match(rownames(umap_df),
                                                                 rownames(seu@meta.data))]
      p_my_umap <- ggplot(umap_df, aes(umap_1, umap_2)) +
        geom_point(data = umap_df[is.na(umap_df[[my_lab_col]]), ],
                   colour = "grey92", size = 0.25) +
        geom_point(data = umap_df[!is.na(umap_df[[my_lab_col]]), ],
                   aes(colour = .data[[my_lab_col]]), size = 0.6, alpha = 0.85) +
        guides(colour = guide_legend(override.aes = list(size = 4))) +
        labs(title = paste0(case_id, " — myeloid state (", setnm, ")"),
             subtitle = "Grey = not myeloid (not tested). Hybrid labels are canonical (order-free).",
             colour = NULL) +
        theme_classic()
      print(p_my_umap)

      sp_df[[my_lab_col]] <- seu@meta.data[[my_lab_col]][match(rownames(sp_df),
                                                               rownames(seu@meta.data))]
      p_my_sp <- ggplot(sp_df, aes(x = imagerow, y = imagecol)) +
        geom_point(data = sp_df[is.na(sp_df[[my_lab_col]]), ],
                   colour = "grey93", size = 0.3, shape = 16) +
        geom_point(data = sp_df[!is.na(sp_df[[my_lab_col]]), ],
                   aes(colour = .data[[my_lab_col]]), size = 0.55, shape = 16) +
        scale_x_reverse() + coord_fixed() + theme_void() +
        labs(title = paste0(case_id, " — myeloid state, spatial (", setnm, ")"),
             colour = NULL)
      print(p_my_sp)
    }
  }

  # gene-count / power context — the reason the two sets are not comparable
  if (!is.null(myeloid_summary)) {
    p_my_pow <- ggplot(myeloid_summary, aes(n_genes, pct_pos, colour = set)) +
      geom_point(size = 2) +
      ggrepel::geom_text_repel(aes(label = state), size = 2.5, max.overlaps = 20) +
      scale_colour_brewer(palette = "Dark2") +
      labs(title = paste0(case_id, " — positivity vs genes tested"),
           subtitle = paste("Sensitivity rises steeply with set size (simulation: ~25 genes",
                            "-> ~0.50, ~80 genes -> ~0.78). Sets with few genes give a FLOOR."),
           x = "genes tested", y = "% cells positive", colour = NULL) +
      theme_bw(base_size = 9)
    print(p_my_pow)
  }
}

# 8e. DE volcano --------------------------------------------------------------
sig_de <- de_results %>%
  filter(p_val_adj < ALPHA) %>%
  group_by(direction) %>%
  slice_max(abs(avg_log2FC), n = 10)

p_volcano <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300),
                                    colour = direction)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", colour = "grey50") +
  ggrepel::geom_text_repel(data = sig_de, aes(label = gene),
                           size = 3, max.overlaps = 20) +
  scale_colour_manual(values = c("TdTom_pos" = "blue", "TdTom_neg" = "#e6550d")) +
  labs(title = paste0(case_id, " — DE TdTom+ vs TdTom-"),
       x = "avg log2FC", y = "-log10(adj p-value)", colour = NULL) +
  theme_bw()
print(p_volcano)

# 8f. DE heatmap --------------------------------------------------------------
top_de_genes <- de_results %>%
  filter(p_val_adj < ALPHA) %>%
  group_by(direction) %>%
  slice_max(abs(avg_log2FC), n = 20) %>%
  pull(gene)

if (length(top_de_genes) > 0) {
  p_heat <- DoHeatmap(seu, features = top_de_genes, group.by = "tdtom_label",
                      assay = "SCT") +
    labs(title = paste0(case_id, " — Top DE genes")) +
    theme(axis.text.y = element_text(size = 7))
  print(p_heat)
}

dev.off()
message("  Plots saved: ", pdf_path)

# =============================================================================
# 9. Summary tables
# =============================================================================

message("\n=== Step 9: Summary tables ===")

summary_tbl <- seu@meta.data %>%
  as.data.frame() %>%
  {
    grp_vars <- c("cell_class", "prolif_state", "tdtom_label")
    if (USE_EC) grp_vars <- c(grp_vars, "ec_state")
    # v3: myeloid state summary label, one column per signature set
    if (USE_MYELOID && USE_MYELOID_RESULTS)
      grp_vars <- c(grp_vars, paste0("myeloid_", names(myeloid_tbl_list),
                                     "_state_canonical"))
    group_by(., across(all_of(grp_vars)))
  } %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(tdtom_label) %>%
  mutate(pct_within_tdtom = round(n_cells / sum(n_cells) * 100, 2)) %>%
  ungroup() %>%
  mutate(pct_total = round(n_cells / ncol(seu) * 100, 2),
         sample    = case_id)

write.csv(summary_tbl,
          file.path(OUT_DIR, paste0(case_id, "_cell_summary.csv")),
          row.names = FALSE)
message("  Summary saved.")
print(summary_tbl, n = 40)

# Per-cell table — carries hybrid labels, decision statistics and the EC
# independent calls, all of which the grouped summary collapses away.
percell_cols <- intersect(
  c("cell_class", "primary_sub", "hybrid_class", "hybrid_sub", "is_hybrid",
    "state_label", "state_gap12", "state_gap23", "prolif_state", "combined_class",
    "tdtom_label", "z_max", "gate_min", "pval_min", "gate_prolif_min",
    if (USE_EC) c("ec_state", "ec_n_positive", "gate_ec_min") else NULL,
    if (USE_EC && USE_EC_RESULTS)
      unlist(lapply(ec_names, function(s) paste0("ec_", s, c("_z", "_gate", "_pos"))))
    else NULL,
    # v3: order-free Neftel hybrid labels — report these, not the rank-ordered ones
    "hybrid_class_canonical", "hybrid_sub_canonical",
    # v3: myeloid state calls. Per-state *_pos columns are the primary output;
    # *_state_canonical is a convenience summary (states are not exclusive).
    if (USE_MYELOID && USE_MYELOID_RESULTS)
      unlist(lapply(names(myeloid_tbl_list), function(sn) {
        pf <- paste0("myeloid_", sn, "_")
        c(paste0(pf, c("state", "state_canonical", "hybrid", "hybrid_canonical",
                       "is_hybrid", "n_positive")),
          unlist(lapply(sub("_z$", "", grep("_z$", colnames(myeloid_tbl_list[[sn]]),
                                            value = TRUE)),
                        function(st) paste0(pf, st, c("_score", "_z", "_gate", "_pos")))))
      }))
    else NULL),
  colnames(seu@meta.data))

percell_tbl <- seu@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, all_of(percell_cols)) %>%
  mutate(sample = case_id)
write.csv(percell_tbl,
          file.path(OUT_DIR, paste0(case_id, "_percell_states.csv")),
          row.names = FALSE)
message("  Per-cell state table saved (", nrow(percell_tbl), " cells).")

# Full decision statistics, one row per cell x signature — lets you re-threshold
# without re-running the permutations.
stats_long <- as.data.frame(z_mat) %>% rownames_to_column("cell_id") %>%
  pivot_longer(-cell_id, names_to = "signature", values_to = "z") %>%
  left_join(as.data.frame(q_mat) %>% rownames_to_column("cell_id") %>%
              pivot_longer(-cell_id, names_to = "signature", values_to = "gate_stat"),
            by = c("cell_id", "signature")) %>%
  left_join(as.data.frame(pval_mat) %>% rownames_to_column("cell_id") %>%
              pivot_longer(-cell_id, names_to = "signature", values_to = "pval_emp"),
            by = c("cell_id", "signature")) %>%
  filter(!is.na(z)) %>%
  mutate(gate_statistic = GATE, sample = case_id)
write.csv(stats_long,
          file.path(OUT_DIR, paste0(case_id, "_enrichment_stats.csv")),
          row.names = FALSE)
message("  Enrichment statistics saved (", nrow(stats_long), " cell x signature rows).")

# Hybrid-pair frequency table among tumour cells
hyb_tbl <- as.data.frame(table(seu$hybrid_class[tumour_cells]),
                         stringsAsFactors = FALSE)
if (nrow(hyb_tbl)) {
  colnames(hyb_tbl) <- c("hybrid_pair", "n_cells")
  hyb_tbl$pct_of_tumour <- round(hyb_tbl$n_cells / length(tumour_cells) * 100, 2)
  hyb_tbl$sample <- case_id
  write.csv(hyb_tbl,
            file.path(OUT_DIR, paste0(case_id, "_hybrid_pairs.csv")),
            row.names = FALSE)
  message("  Hybrid-pair table saved.")
  print(hyb_tbl)
}

# Run parameters, so a result can be traced to the settings that produced it
params_tbl <- data.frame(
  parameter = c("sample", "n_perm", "alpha", "seed", "nbin", "n_det_bin",
                "null", "gate_statistic", "rank_statistic",
                "expr_ref_neftel", "expr_ref_celltype", "expr_ref_prolif",
                "expr_ref_ec", "ec_unique_genes",
                "hybrid_delta_q", "hybrid_pct_floor", "hybrid_delta_resolved",
                "n_tdtom_pos", "n_tdtom_neg", "calibration_check",
                "gate_unit", "p_source", "auto_strata", "strata_target",
                "borrow_strata", "skew_p", "pool_k", "tier_z_floor"),
  value = c(case_id, N_PERM, ALPHA, SEED, NBIN, N_DET_BIN,
            "expression x detection matched", GATE, "z",
            "self (TdTom+)", "self (TdTom-)", "all cells",
            if (USE_EC) "self (Endothelial)" else NA, EC_UNIQUE,
            HYBRID_DELTA_Q, HYBRID_PCT_FLOOR, signif(attr(hyb, "delta"), 6),
            length(tumour_cells), length(nontumour_cells), DO_CALIB,
            GATE_UNIT, P_SOURCE, AUTO_STRATA, STRATA_TARGET,
            BORROW_STRATA, SKEW_P, POOL_K, TIER_Z_FLOOR),
  stringsAsFactors = FALSE)
write.csv(params_tbl, file.path(OUT_DIR, paste0(case_id, "_run_parameters.csv")),
          row.names = FALSE)
message("  Run parameters saved.")

# =============================================================================
# 10. GeoJSON export
# =============================================================================

message("\n=== Step 10: GeoJSON export ===")

if (!all(c("imagerow", "imagecol") %in% colnames(seu@meta.data))) {
  tissue_coords <- GetTissueCoordinates(seu)
  seu@meta.data$imagerow <- tissue_coords[rownames(seu@meta.data), "imagerow"]
  seu@meta.data$imagecol <- tissue_coords[rownames(seu@meta.data), "imagecol"]
}

geo_df <- seu@meta.data %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, imagerow, imagecol,
                cell_class, prolif_state, combined_class, tdtom_label,
                state_label, is_hybrid,
                any_of(if (USE_EC) c("ec_state", "ec_n_positive") else character(0)))

make_geojson_points <- function(df, class_col, colour_map = NULL) {
  ec_col_present <- "ec_state" %in% colnames(df)
  features <- lapply(seq_len(nrow(df)), function(i) {
    cls   <- df[[class_col]][i]
    props <- list(
      cell_id        = df$cell_id[i],
      cell_class     = df$cell_class[i],
      state_tier     = df$state_tier[i],
      tier_top       = df$tier_top[i],
      tier_top_z     = df$tier_top_z[i],
      prolif_state   = df$prolif_state[i],
      combined_class = df$combined_class[i],
      tdtom_label    = df$tdtom_label[i],
      state_label    = df$state_label[i],
      is_hybrid      = df$is_hybrid[i],
      classification = cls
    )
    if (ec_col_present) {
      props$ec_state       <- df$ec_state[i]
      props$ec_n_positive  <- df$ec_n_positive[i]
    }
    if (!is.null(colour_map) && !is.na(cls) && cls %in% names(colour_map))
      props$color <- colour_map[[cls]]

    list(
      type     = "Feature",
      geometry = list(type = "Point",
                      coordinates = list(df$imagecol[i], df$imagerow[i])),
      properties = props
    )
  })
  list(type = "FeatureCollection", features = features)
}

message("  Writing neftel_normal GeoJSON ...")
gj_class <- make_geojson_points(geo_df, "cell_class", as.list(class_colours))
write(toJSON(gj_class, auto_unbox = TRUE),
      file.path(OUT_DIR, paste0(case_id, "_neftel_normal.geojson")))

message("  Writing proliferation GeoJSON ...")
gj_prolif <- make_geojson_points(geo_df, "prolif_state", as.list(prolif_colours))
write(toJSON(gj_prolif, auto_unbox = TRUE),
      file.path(OUT_DIR, paste0(case_id, "_proliferation.geojson")))

message("  Writing intersection GeoJSON ...")
combined_classes <- unique(geo_df$combined_class)
combined_colours <- setNames(
  sapply(combined_classes, function(cc) {
    base_class <- strsplit(cc, "\\|")[[1]][1]
    if (base_class %in% names(class_colours)) class_colours[[base_class]] else "lightgrey"
  }),
  combined_classes
)
gj_combined <- make_geojson_points(geo_df, "combined_class", as.list(combined_colours))
write(toJSON(gj_combined, auto_unbox = TRUE),
      file.path(OUT_DIR, paste0(case_id, "_intersection.geojson")))

message("  Writing TdTomato binary GeoJSON ...")
tdtom_colours <- list(TdTom_pos = "blue", TdTom_neg = "lightgrey")
gj_tdtom <- make_geojson_points(geo_df, "tdtom_label", tdtom_colours)
write(toJSON(gj_tdtom, auto_unbox = TRUE),
      file.path(OUT_DIR, paste0(case_id, "_tdtomato.geojson")))

# Hybrid state label (single states + "A|B" hybrids) — new in v2
message("  Writing state_label (incl. hybrids) GeoJSON ...")
state_lvls    <- unique(na.omit(geo_df$state_label))
state_colours <- setNames(
  sapply(state_lvls, function(sl) {
    base <- strsplit(sl, "\\|")[[1]][1]
    if (base %in% names(class_colours)) class_colours[[base]] else "lightgrey"
  }), state_lvls)
gj_state <- make_geojson_points(geo_df %>% filter(!is.na(state_label)),
                                "state_label", as.list(state_colours))
write(toJSON(gj_state, auto_unbox = TRUE),
      file.path(OUT_DIR, paste0(case_id, "_state_label.geojson")))

if (USE_EC && USE_EC_RESULTS) {
  message("  Writing EC state GeoJSON (Endothelial cells only) ...")
  geo_df_ec <- geo_df %>% filter(!is.na(ec_state))
  gj_ec <- make_geojson_points(geo_df_ec, "ec_state", as.list(ec_colours))
  write(toJSON(gj_ec, auto_unbox = TRUE),
        file.path(OUT_DIR, paste0(case_id, "_ec_state.geojson")))
}

message("  GeoJSON files saved.")

# =============================================================================
# 11. Save updated Seurat object
# =============================================================================

message("\n=== Step 11: Saving Seurat object ===")
rds_out <- file.path(OUT_DIR, paste0(case_id, "_classified.rds"))
saveRDS(seu, rds_out)
message("  RDS saved: ", rds_out)

message("\n=== Done: ", case_id, " ===")
message("Reminder: inspect ", case_id, "_null_calibration.csv before fixing thresholds.")
message("  (i)  gene-label shuffle should call ~", 100 * ALPHA, "% of cells on raw p;")
message("  (ii) a tumour-restricted but NOT state-specific control signature")
message("       should not fire in nearly every TdTom+ cell.")
